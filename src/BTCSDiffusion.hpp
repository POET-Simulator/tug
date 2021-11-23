#ifndef BTCSDIFFUSION_H_
#define BTCSDIFFUSION_H_

#include <vector>
#include <Eigen/SparseCore>

typedef int BCSide;
typedef Eigen::Triplet<double> T;

class BTCSDiffusion {

public:
  static const BCSide LEFT;
  static const BCSide RIGHT;

  BTCSDiffusion(int x);
  BTCSDiffusion(int x, int y);
  BTCSDiffusion(int x, int y, int z);

  void setBoundaryCondition(std::vector<double> input, BCSide side);
  void simulate(std::vector<double> &c, std::vector<double> &alpha,
                double timestep);

private:
  std::vector<double> bc;
  int grid_dim;
  int dim_x;
  int dim_y;
  int dim_z;
};

#endif // BTCSDIFFUSION_H_
