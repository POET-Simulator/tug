#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include <Eigen/SparseCore>
#include <vector>

typedef Eigen::Triplet<double> T;

extern void BTCS2D(int x, int y, std::vector<double> &c, double alpha,
                   double timestep);
#endif // DIFFUSION_H_
