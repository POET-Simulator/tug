
<a name="v0.2"></a>
## [v0.2](https://git.gfz-potsdam.de/sec34/btcsdiffusion/compare/v0.1...v0.2) (2022-08-16)

### Build System

* fetch doctest during configuration

### Ci

* disable testing during static analyze
* add git as dependency

### Code Refactoring

* remove BTCSUtils header from include API

### Code Style

* fix various code style recommendations from clang
* Use enumerations for macros and use more useful function names

### Documentation

* update Roadmap and add Contributing section

### Features

* support for inner closed cells in diffusion module
* add setting of inner boundary conditions

### Housework

* configure git-chglog for new commit style

### Testing

* add new test case for diffusion module
* add tests for inner boundary conditions

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
  
