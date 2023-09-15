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

// entry point; differentiate between 1D and 2D grid
extern void FTCS(Grid &grid, Boundary &bc, double &timestep, int &numThreads);

// entry point for EigenLU solver; differentiate between 1D and 2D grid
extern void BTCS_LU(Grid &grid, Boundary &bc, double timestep, int numThreads);

// entry point for Thomas algorithm solver; differentiate 1D and 2D grid
extern void BTCS_Thomas(Grid &grid, Boundary &bc, double timestep,
                        int numThreads);

#endif // SCHEMES_H_
