<a name="v0.3"></a>
## [v0.3](https://git.gfz-potsdam.de/sec34/tug/compare/v0.2...v0.3) (2022-09-08)

### Bug Fixes

* grid dimensions were stored and accessed incorrectly

### Build System

* remove BoundaryCondition as extra library

### Code Refactoring

* move includes into subdirectory 'tug'
* move BoundaryCondition header and source
* rename BoundaryCondition class
* rename and expand namespace

### Continious Integration

* linting needs to be triggered manually now

### Doc

* remove old stuff from ADI documentation

### Documentation

* Update and extending README

### Features

* allow undefined boundary conditions
* add helper functions to TugInput struct
* Remove class BTCSDiffusion

### Housework

* remove unneeded test file
* update Changelog link to new name
* Change URL of repo and and description for CI
* moved Comp*.R to scripts/

### Performance Improvements

* represent inner boundary conditions with a std::map

### Testing

* enable building of tests per default
* add target `check`

### BREAKING CHANGE

Functionality is now provided by function calls and
scheme generation is decoupled from LEqS solving.

<a name="v0.2"></a>
## [v0.2](https://git.gfz-potsdam.de/sec34/tug/compare/v0.1...v0.2) (2022-08-16)

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
  
