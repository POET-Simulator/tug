#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

namespace Diffusion {

/*!
 * Defines both types of boundary condition as a datatype.
 */
typedef int bctype;

/*!
 * Defines a closed/Neumann boundary condition.
 */
static const bctype BC_CLOSED = 0;

/*!
 * Defines a flux/Cauchy boundary condition.
 */
static const bctype BC_FLUX = 1;

/*!
 * Defines a constant/Dirichlet boundary condition.
 */
static const bctype BC_CONSTANT = 2;

typedef struct boundary_condition {
  bctype type;
  double value;
} boundary_condition;

} // namespace Diffusion

#endif // BOUNDARYCONDITION_H_
