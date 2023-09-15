/**
 * @file BTCSv2.cpp
 * @brief Implementation of heterogenous BTCS (backward time-centered space)
 * solution of diffusion equation in 1D and 2D space. Internally the
 * alternating-direction implicit (ADI) method is used. Version 2, because
 * Version 1 was an implementation for the homogeneous BTCS solution.
 *
 */

#ifndef SCHEMES_H_
#define SCHEMES_H_

#include "TugUtils.hpp"

#include <tug/Boundary.hpp>
#include <tug/Grid.hpp>

namespace tug {

// entry point; differentiate between 1D and 2D grid
template <class T>
extern void FTCS(tug::Grid<T> &grid, tug::Boundary<T> &bc, T timestep,
                 int &numThreads);

// entry point for EigenLU solver; differentiate between 1D and 2D grid
template <class T>
extern void BTCS_LU(tug::Grid<T> &grid, tug::Boundary<T> &bc, T timestep,
                    int numThreads);

// entry point for Thomas algorithm solver; differentiate 1D and 2D grid
template <class T>
extern void BTCS_Thomas(tug::Grid<T> &grid, tug::Boundary<T> &bc, T timestep,
                        int numThreads);
} // namespace tug

#endif // SCHEMES_H_
