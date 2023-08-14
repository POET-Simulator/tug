<h1>src-Directory</h1>
This is the src-directory that holds the source code to the implementation of the TUG framework.


<pre>
src/  
  ├── CMakeFiles/
  ├── Boundary.cpp  
  ├── BoundaryCondition.cpp 
  ├── BTCS.cpp  
  ├── BTCSv2.cpp
  ├── CMakeLists.txt  
  ├── FTCS.cpp
  ├── Grid.cpp 
  ├── README.md
  ├── Simulation.cpp
  ├── Solver.cpp  
  └── TugUtils.hpp 
</pre>

**src/** Directory with the source code.

**CMakeFiles/** Automatically generated directory by CMake.

**Boundary.cpp** Implementation of Boundary class, that holds the boundary conditions.

**BoundaryCondition.cpp** <i>Outdated</i> implementation of boundary conditions.

**BTCS.cpp** <i>Outdated</i> implementation of BTCS solution to homogeneous diffusion.

**BTCSv2.cpp** Implementation of BTCS solution to heterogeneous diffusion in 1D and 2D. 

**CMakeLists.txt** CMakeLists for this directory.

**FTCS.cpp** Implementation of FTCS solution to heterogeneous diffusion in 1D and 2D. 

**Grid.cpp** Implementation of Grid class, that holds all of the concentrations alpha coefficient in x- and y-direction.

**README.md** <i>This</i> file.

**Simulation.cpp** Implementation of Simulation class, that holds all of the information for a specific simulation run, as well as the Boundary and Grid objects. 

**Solver.cpp** <i>Outdated</i> implementation of Eigen solvers.

**TugUtils.hpp** Helper functions for other source files. 

