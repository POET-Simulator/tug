#pragma once

#include "../Boundary.hpp"
#include "../Grid.hpp"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <sycl/sycl.hpp>

namespace tug {

template <typename T>
static inline constexpr T
evaluateBoundary(const T conc, const T conc_neighbor, const T alpha,
                 const T alpha_neighbor, const T bcValue,
                 const BC_TYPE bcType) {
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
static inline constexpr T
calculateExplicitChange(const std::size_t idx_base, const std::size_t idx_n1,
                        const std::size_t idx_n2, const T *conc,
                        const T *alpha) {
  const T alpha_intercell_n1 =
      calcAlphaIntercell(alpha[idx_base], alpha[idx_n1]);
  const T alpha_intercell_n2 =
      calcAlphaIntercell(alpha[idx_base], alpha[idx_n2]);

  const T conc_term_n1 = alpha_intercell_n1 * conc[idx_n1];
  const T conc_term_n2 = alpha_intercell_n2 * conc[idx_n2];

  return conc_term_n1 + conc_term_n2 -
         (alpha_intercell_n1 + alpha_intercell_n2) * conc[idx_base];
}

template <typename T>
static inline void
FTCS_SYCL_2D(sycl::queue &q, T *concentrations_inout, const T *alphas_x_in,
             const T *alphas_y_in, const tug::BoundaryElement<T> *boundaries,
             const std::size_t colMax, const std::size_t rowMax,
             const T dtOverDeltaColSquare, const T dtOverDeltaRowSquare,
             const std::size_t nIter) {

  T *concentrations_device = sycl::malloc_device<T>(colMax * rowMax, q);
  T *concentrations_t1_device = sycl::malloc_device<T>(colMax * rowMax, q);
  T *alpha_x_device = sycl::malloc_device<T>(colMax * rowMax, q);
  T *alpha_y_device = sycl::malloc_device<T>(colMax * rowMax, q);
  tug::BoundaryElement<T> *boundaries_device =
      sycl::malloc_device<tug::BoundaryElement<T>>(2 * colMax + 2 * rowMax, q);

  q.memcpy(concentrations_device, concentrations_inout,
           colMax * rowMax * sizeof(T));
  q.memcpy(alpha_x_device, alphas_x_in, colMax * rowMax * sizeof(T));
  q.memcpy(alpha_y_device, alphas_y_in, colMax * rowMax * sizeof(T));
  q.memcpy(boundaries_device, boundaries,
           (2 * colMax + 2 * rowMax) * sizeof(BoundaryElement<T>));

  const tug::BoundaryElement<T> *leftbound = &boundaries_device[0];
  const tug::BoundaryElement<T> *rightbound = &boundaries_device[rowMax];
  const tug::BoundaryElement<T> *topbound = &boundaries_device[2 * rowMax];
  const tug::BoundaryElement<T> *bottombound =
      &boundaries_device[2 * rowMax + colMax];

  const std::size_t wgMaxSize =
      q.get_device().get_info<sycl::info::device::max_work_group_size>();

  constexpr std::size_t CUDA_WARP_SIZE = 32;
  const std::size_t col_wg = std::min(colMax, CUDA_WARP_SIZE);
  const std::size_t row_wg = wgMaxSize / col_wg;

  auto round_to_multiple = [](std::size_t value, std::size_t multiple) {
    return ((value + multiple - 1) / multiple) * multiple;
  };

  sycl::range<2> wg_range(col_wg, row_wg);
  sycl::range<2> global_range = sycl::range<2>(
      round_to_multiple(colMax, col_wg), round_to_multiple(rowMax, row_wg));

  q.wait();

  for (std::size_t i = 0; i < nIter; i++) {
    q.submit([&](sycl::handler &h) {
      h.parallel_for(
          sycl::nd_range<2>(global_range, wg_range),
          [=](sycl::nd_item<2> item2DIndex) {
            const std::size_t col_i = item2DIndex.get_global_id(0);
            const std::size_t row_j = item2DIndex.get_global_id(1);

            if (col_i >= colMax || row_j >= rowMax) {
              return;
            }

            const std::size_t idx_self = col_i + row_j * colMax;
            const std::size_t idx_left = idx_self - 1;
            const std::size_t idx_right = idx_self + 1;
            const std::size_t idx_top = idx_self - colMax;
            const std::size_t idx_bottom = idx_self + colMax;

            T vertical_change;

            if (row_j == 0 || row_j == rowMax - 1) {
              const std::size_t neighbor_idx =
                  row_j == 0 ? idx_bottom : idx_top;
              const tug::BoundaryElement<T> &bound =
                  row_j == 0 ? topbound[col_i] : bottombound[col_i];

              vertical_change = evaluateBoundary(
                  concentrations_device[idx_self],
                  concentrations_device[neighbor_idx], alpha_y_device[idx_self],
                  alpha_y_device[neighbor_idx], bound.getValue(),
                  bound.getType());
            } else {
              vertical_change = calculateExplicitChange(
                  idx_self, idx_bottom, idx_top, concentrations_device,
                  alpha_y_device);
            }

            T horizontal_change;

            if (col_i == 0 || col_i == colMax - 1) {
              const std::size_t neighbor_idx =
                  col_i == 0 ? idx_right : idx_left;
              const tug::BoundaryElement<T> &bound =
                  col_i == 0 ? leftbound[row_j] : rightbound[row_j];

              horizontal_change = evaluateBoundary(
                  concentrations_device[idx_self],
                  concentrations_device[neighbor_idx], alpha_x_device[idx_self],
                  alpha_x_device[neighbor_idx], bound.getValue(),
                  bound.getType());
            } else {
              horizontal_change = calculateExplicitChange(
                  idx_self, idx_right, idx_left, concentrations_device,
                  alpha_x_device);
            }
            concentrations_t1_device[idx_self] =
                concentrations_device[idx_self] +
                dtOverDeltaColSquare * horizontal_change +
                dtOverDeltaRowSquare * vertical_change;
          });
    });

    q.wait();

    T *tmp = concentrations_device;
    concentrations_device = concentrations_t1_device;
    concentrations_t1_device = tmp;
  }

  q.memcpy(concentrations_inout, concentrations_device,
           colMax * rowMax * sizeof(T))
      .wait();

  sycl::free(concentrations_device, q);
  sycl::free(concentrations_t1_device, q);
  sycl::free(alpha_x_device, q);
  sycl::free(alpha_y_device, q);
  sycl::free(boundaries_device, q);
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

  const std::size_t max_wg_size =
      q.get_device().get_info<sycl::info::device::max_work_group_size>();

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
                  evaluateBoundary(concentrations_device[0],
                                   concentrations_device[1], alpha_device[0],
                                   alpha_device[1], bcValue_left, bcType_left);
        });
      });

      auto e_boundright = q.submit([&](sycl::handler &h) {
        h.single_task([=]() {
          concentrations_t1_device[colMax - 1] =
              concentrations_device[colMax - 1] +
              dtOverDeltaColSquare *
                  evaluateBoundary(concentrations_device[colMax - 1],
                                   concentrations_device[colMax - 2],
                                   alpha_device[colMax - 1],
                                   alpha_device[colMax - 2], bcValue_right,
                                   bcType_right);
        });
      });

      auto e_mid = q.submit([&](sycl::handler &h) {
        h.parallel_for(sycl::nd_range<1>(sycl::range<1>(colMax - 2),
                                         sycl::range<1>(max_wg_size)),
                       [=](sycl::nd_item<1> item) {
                         const std::size_t i = item.get_global_id(0) + 1;
                         concentrations_t1_device[i] =
                             concentrations_device[i] +
                             dtOverDeltaColSquare *
                                 calculateExplicitChange(i, i + 1, i - 1,
                                                         concentrations_device,
                                                         alpha_device);
                       });
      });

      q.wait();

      T *tmp = concentrations_device;
      concentrations_device = concentrations_t1_device;
      concentrations_t1_device = tmp;

      // q.memcpy(concentrations_device, concentrations_t1_device,
      //          colMax * sizeof(T), {e_boundleft, e_boundright, e_mid})
      //     .wait();
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
  } else if (grid.getDim() == 2) {
    const std::size_t rowMax = grid.getRow();
    const T deltaRow = grid.getDeltaRow();
    const T deltaRowSquare = deltaRow * deltaRow;
    const T dtOverDeltaRowSquare = dt / deltaRowSquare;

    const auto serialized_boundaries = boundaries.serialize();

    FTCS_SYCL_2D(q, grid.getConcentrations().data(), grid.getAlphaX().data(),
                 grid.getAlphaY().data(), serialized_boundaries.data(), colMax,
                 rowMax, dtOverDeltaColSquare, dtOverDeltaRowSquare, nIter);
  }
}
} // namespace tug