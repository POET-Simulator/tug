Examples
========

At this point, some typical commented examples are presented to illustrate how Tug works.
In general, each simulation is divided into three blocks:
    - the initialization of the grid, which is to be simulated
    - the setting of the boundary conditions
    - the setting of the simulation parameters and the start of the simulation

Two dimensional grid with constant boundaries and FTCS method
-------------------------------------------------------------
**Initialization of the grid**

For example, the initalization of a grid with 20 by 20 cells using double values
and a domain size (physical extent of the grid) of also 20 by 20 length units
can be done as follows. The setting of the domain is optional here and is set to
the same size as the number of cells in the standard case. As seen in the code,
the cells of the grid are set to an initial value of 0 and only in the upper
left corner (0,0) the starting concentration is set to the value 20.

.. code-block:: cpp
    
    int row = 20
    int col = 20;
    Grid<double> grid(row,col);
    grid.setDomain(row, col);
    MatrixXd concentrations = MatrixXd::Constant(row,col,0);
    // or MatrixX<double> concentrations = MatrixX<double>::Constant(row,col,0);
    concentrations(0,0) = 20;
    grid.setConcentrations(concentrations);

**Setting of the boundary conditions**

First, a boundary class is created and then the corresponding boundary conditions are set. In this case, all four sides
of the grid are set as constant edges with a concentration of 0.

.. code-block:: cpp

    Boundary bc = Boundary(grid);
    bc.setBoundarySideConstant(BC_SIDE_LEFT, 0);
    bc.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
    bc.setBoundarySideConstant(BC_SIDE_TOP, 0);
    bc.setBoundarySideConstant(BC_SIDE_BOTTOM, 0);

**Setting of the simulation parameters and simulation start**

In the last block, a simulation class is created and the objects of the grid and
the boundary conditions are passed. The solution method is also specified
(either FCTS or BTCS). Furthermore, the desired time step and the number of
iterations are set. The penultimate parameter specifies the output of the
simulated results in a CSV file. In the present case, the result of each
iteration step is written one below the other into the corresponding CSV file.

.. code-block:: cpp
    
    Simulation<double, FTCS_APPROACH> simulation(grid, bc);
    simulation.setTimestep(0.1);
    simulation.setIterations(1000);
    simulation.setOutputCSV(CSV_OUTPUT_VERBOSE);
    simulation.run();






Setting special boundary conditions on individual cells
-------------------------------------------------------
