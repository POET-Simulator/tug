
<a name="v0.1"></a>
## v0.1 (2022-08-09)

* Basic solving of diffusion problems with
  - 1D regular and rectangular grids using BTCS scheme and Eigen SparseLU solver
  - 2D regular and rectangular grids using 2D-ADI-BTCS scheme and Eigen SparseLU solver
* Definition of boundary conditions via class `BTCSBoundaryCondition` on ghost nodes 
* Boundaries types `CLOSED` and `CONSTANT` cells are provided for diffusion problem solving
* Software testing of both `BTCSDiffusion` and `BTCSBoundaryCondition` classes
* Description of both BTCS and 2D-ADI-BTCS schemes are provided in `doc`
* Example applications are attached in `app` subdirectory
