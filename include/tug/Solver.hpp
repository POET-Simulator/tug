#ifndef SOLVER_H_
#define SOLVER_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace tug {
namespace solver {

/**
 * Solving linear equation system with LU-decomposition implemented in Eigen
 * library.
 *
 * \param A_matrix The A matrix represented as a sparse matrix using Eigen
 * library.
 *
 * \param b_vector Right hand side vector of the linear equation system.
 *
 * \return Solution represented as vector.
 */
auto EigenLU(const Eigen::SparseMatrix<double> &A_matrix,
             const Eigen::VectorXd &b_vector) -> Eigen::VectorXd;

/**
 * Solving linear equation system with brutal implementation of the Thomas
 * algorithm.
 *
 * \param A_matrix The A matrix represented as a sparse matrix using Eigen
 * library.
 *
 * \param b_vector Right hand side vector of the linear equation system.
 *
 * \return Solution represented as vector.
 */
auto ThomasAlgorithm(const Eigen::SparseMatrix<double> &A_matrix,
                     const Eigen::VectorXd &b_vector) -> Eigen::VectorXd;

} // namespace solver
} // namespace tug

#endif // SOLVER_H_
