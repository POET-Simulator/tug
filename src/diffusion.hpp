#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include <Eigen/SparseCore>
#include <vector>

typedef Eigen::Triplet<double> T;

extern void BTCS1D(int x, std::vector<double> &c, std::vector<double> &alpha,
                   double timestep, std::vector<double> &bc);

extern void BTCS2D(int x, int y, std::vector<double> &c,
                   std::vector<double> &alpha, double timestep);
#endif // DIFFUSION_H_
