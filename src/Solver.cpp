#include <tug/Solver.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

auto tug::solver::EigenLU(const Eigen::SparseMatrix<double> &A_matrix,
                          const Eigen::VectorXd &b_vector) -> Eigen::VectorXd {

  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A_matrix);

  solver.factorize(A_matrix);

  Eigen::VectorXd x_vec = solver.solve(b_vector);

  return solver.solve(b_vector);
}

auto tug::solver::ThomasAlgorithm(const Eigen::SparseMatrix<double> &A_matrix,
                                  const Eigen::VectorXd &b_vector)
    -> Eigen::VectorXd {
  uint32_t n = b_vector.size();

  Eigen::VectorXd a_diag(n);
  Eigen::VectorXd b_diag(n);
  Eigen::VectorXd c_diag(n);
  Eigen::VectorXd x_vec = b_vector;

  // Fill diagonals vectors
  b_diag[0] = A_matrix.coeff(0, 0);
  c_diag[0] = A_matrix.coeff(0, 1);

  for (int i = 1; i < n - 1; i++) {
    a_diag[i] = A_matrix.coeff(i, i - 1);
    b_diag[i] = A_matrix.coeff(i, i);
    c_diag[i] = A_matrix.coeff(i, i + 1);
  }
  a_diag[n - 1] = A_matrix.coeff(n - 1, n - 2);
  b_diag[n - 1] = A_matrix.coeff(n - 1, n - 1);

  // start solving - c_diag and x_vec are overwritten
  n--;
  c_diag[0] /= b_diag[0];
  x_vec[0] /= b_diag[0];

  for (int i = 1; i < n; i++) {
    c_diag[i] /= b_diag[i] - a_diag[i] * c_diag[i - 1];
    x_vec[i] = (x_vec[i] - a_diag[i] * x_vec[i - 1]) /
               (b_diag[i] - a_diag[i] * c_diag[i - 1]);
  }

  x_vec[n] = (x_vec[n] - a_diag[n] * x_vec[n - 1]) /
             (b_diag[n] - a_diag[n] * c_diag[n - 1]);

  for (int i = n; i-- > 0;) {
    x_vec[i] -= c_diag[i] * x_vec[i + 1];
  }

  return x_vec;
}
