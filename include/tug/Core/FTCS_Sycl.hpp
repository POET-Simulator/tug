#pragma once

#include "../Boundary.hpp"
#include "../Grid.hpp"
#include "sycl/usm.hpp"

#include <cstddef>
#include <cstring>
#include <sycl/sycl.hpp>

namespace tug {

template <typename T>
static inline constexpr T boundApply(const T conc, const T conc_neighbor,
                                     const T alpha, const T alpha_neighbor,
                                     const T bcValue, const BC_TYPE bcType) {
  const T alpha_intercell = calcAlphaIntercell(alpha, alpha_neighbor);
  switch (bcType) {
  case BC_TYPE::BC_TYPE_CLOSED: {
    const T conc_change = conc_neighbor - conc;

    return alpha_intercell * conc_change;
  }
  case BC_TYPE::BC_TYPE_CONSTANT: {
    const T inner_term = alpha_intercell * conc_neighbor;
    const T mid_term = (alpha_intercell + 2 * alpha) * conc;
    const T outer_term = 2 * alpha * bcValue;

    return inner_term - mid_term + outer_term;
  }
  default:;
  }
}

template <typename T>
static inline constexpr T explicitChange(const T conc, const T conc_n1,
                                         const T conc_n2, const T alpha,
                                         const T alpha_n1, const T alpha_n2) {
  const T alpha_intercell_n1 = calcAlphaIntercell(alpha, alpha_n1);
  const T alpha_intercell_n2 = calcAlphaIntercell(alpha, alpha_n2);

  const T conc_term_n1 = alpha_intercell_n1 * conc_n1;
  const T conc_term_n2 = alpha_intercell_n2 * conc_n2;

  return conc_term_n1 + conc_term_n2 -
         (alpha_intercell_n1 + alpha_intercell_n2) * conc;
}

template <typename T>
static inline void
FTCS_SYCL_1D(sycl::queue &q, T *concentrations_inout, const T *alphas_in,
             const tug::BoundaryElement<T> &leftbound,
             const tug::BoundaryElement<T> &rightbound,
             const std::size_t colMax, const T dtOverDeltaColSquare,
             const std::size_t nIter) {

  const BC_TYPE bcType_left = leftbound.getType();
  const BC_TYPE bcType_right = rightbound.getType();

  const T bcValue_left = leftbound.getValue();
  const T bcValue_right = rightbound.getValue();

  {
    T *concentrations_device = sycl::malloc_device<T>(colMax, q);
    T *concentrations_t1_device = sycl::malloc_device<T>(colMax, q);
    T *alpha_device = sycl::malloc_device<T>(colMax, q);

    q.memcpy(concentrations_device, concentrations_inout, colMax * sizeof(T));
    q.memcpy(alpha_device, alphas_in, colMax * sizeof(T));

    q.wait();

    for (std::size_t i = 0; i < nIter; i++) {
      auto e_boundleft = q.submit([&](sycl::handler &h) {
        h.single_task([=]() {
          concentrations_t1_device[0] =
              concentrations_device[0] +
              dtOverDeltaColSquare *
                  boundApply(concentrations_device[0], concentrations_device[1],
                             alpha_device[0], alpha_device[1], bcValue_left,
                             bcType_left);
        });
      });

      auto e_boundright = q.submit([&](sycl::handler &h) {
        h.single_task([=]() {
          concentrations_t1_device[colMax - 1] =
              concentrations_device[colMax - 1] +
              dtOverDeltaColSquare *
                  boundApply(concentrations_device[colMax - 1],
                             concentrations_device[colMax - 2],
                             alpha_device[colMax - 1], alpha_device[colMax - 2],
                             bcValue_right, bcType_right);
        });
      });

      auto e_mid = q.submit([&](sycl::handler &h) {
        h.parallel_for(sycl::range<1>(colMax - 2), [=](sycl::id<1> idx) {
          const std::size_t i = idx[0] + 1;
          concentrations_t1_device[i] =
              concentrations_device[i] +
              dtOverDeltaColSquare *
                  explicitChange(concentrations_device[i],
                                 concentrations_device[i + 1],
                                 concentrations_device[i - 1], alpha_device[i],
                                 alpha_device[i + 1], alpha_device[i - 1]);
        });
      });

      q.memcpy(concentrations_device, concentrations_t1_device,
               colMax * sizeof(T), {e_boundleft, e_boundright, e_mid})
          .wait();
    }

    q.memcpy(concentrations_inout, concentrations_device, colMax * sizeof(T))
        .wait();
    sycl::free(concentrations_device, q);
    sycl::free(concentrations_t1_device, q);
    sycl::free(alpha_device, q);
  }
}

template <typename T>
void FTCS_Sycl_impl(sycl::queue &q, tug::Grid<T> &grid,
                    const tug::Boundary<T> &boundaries, const T dt,
                    const std::size_t nIter) {
  const std::size_t colMax = grid.getCol();
  const T deltaCol = grid.getDeltaCol();
  const T deltaColSquare = deltaCol * deltaCol;
  const T dtOverDeltaColSquare = dt / deltaColSquare;

  if (grid.getDim() == 1) {
    FTCS_SYCL_1D(q, grid.getConcentrations().data(), grid.getAlpha().data(),
                 boundaries.getBoundaryElement(BC_SIDE_LEFT, 0),
                 boundaries.getBoundaryElement(BC_SIDE_RIGHT, 0), colMax,
                 dtOverDeltaColSquare, nIter);
  }
}
} // namespace tug