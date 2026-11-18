Examples
========

At this point, some typical commented examples are presented to illustrate how Tug works.
In general, each simulation is divided into three blocks:
    - the initialization of the grid, which is to be simulated
    - the setting of the boundary conditions
    - the setting of the simulation parameters and the start of the simulation

Two dimensional diffusion grid with constant boundaries and FTCS method
-----------------------------------------------------------------------
**Initialization of the grid**

For example, the initialization of a grid with 20 by 20 cells using double values
and a physical domain extent of 20 by 20 length units can be done as follows.
The setting of the domain is optional; by default each cell has dimensions of 1 by 1.
As seen in the code, the cells are initialized with a value of 0, and only in the
upper-left corner (0,0) the starting concentration is set to 20.

.. code-block:: cpp
    
    #include <tug/Diffusion/Diffusion.hpp>
    using namespace tug;

    int rows = 20;
    int cols = 20;
    RowMajMat<double> concentrations = RowMajMat<double>::Zero(rows, cols);
    concentrations(0, 0) = 20.0;

    Diffusion<double, FTCS_APPROACH> simulation(concentrations);
    simulation.setDomain(20.0, 20.0);

**Setting of the boundary conditions**

Boundary conditions are accessed directly from the simulation object. In this case,
all four outer boundaries of the grid are configured as constant Dirichlet boundaries
with a concentration of 0.

.. code-block:: cpp

    Boundary<double> &bc = simulation.getBoundary();
    bc.setBoundarySideConstant(BC_SIDE_LEFT, 0.0);
    bc.setBoundarySideConstant(BC_SIDE_RIGHT, 0.0);
    bc.setBoundarySideConstant(BC_SIDE_TOP, 0.0);
    bc.setBoundarySideConstant(BC_SIDE_BOTTOM, 0.0);

**Setting of the simulation parameters and simulation start**

In the last block, the simulation parameters are configured: the timestep and the
total number of iterations. Optionally, output of intermediate results to a CSV
file or the console can be enabled.

.. code-block:: cpp
    
    simulation.setTimestep(0.1);
    simulation.setIterations(1000);
    simulation.setOutputCSV(CSV_OUTPUT_VERBOSE);
    simulation.run();

Setting special boundary conditions on individual cells
-------------------------------------------------------

Individual cell boundaries can also be customized. For instance, inner sources or
impermeable obstacles (closed boundaries) can be placed at specific grid locations:

.. code-block:: cpp

    // Set an individual cell as a constant source
    bc.setBoundaryElementConstant(5, 5, 100.0);

    // Set an individual cell boundary as closed (no-flux)
    bc.setBoundaryElementClosed(10, 10);

Two dimensional advection simulation with Darcy velocities
----------------------------------------------------------

In this example, fluid flow through a porous medium is simulated by first determining
the hydraulic head distribution and staggered Darcy velocities, and then simulating
conservative advective solute transport using an upwind finite-volume scheme.

**Initialization of hydraulic heads, velocities, and concentrations**

A 20 by 20 cell grid is defined. First, a hydraulic head matrix is initialized with a
gradient from left (high head) to right (low head). A ``Velocities`` object is created
to solve for steady-state Darcy fluxes with a hydraulic conductivity of :math:`10^{-4}\,\text{m/s}`.
Then, a concentration matrix and an ``Advection`` object are created, referencing the
``Velocities`` object, with a porosity of :math:`\phi = 0.25`:

.. code-block:: cpp

    #include <tug/Advection/Velocities.hpp>
    #include <tug/Advection/Advection.hpp>
    using namespace tug;

    int rows = 20;
    int cols = 20;

    // 1. Hydraulic heads and Darcy velocities
    RowMajMat<double> heads = RowMajMat<double>::Constant(rows, cols, 10.0);
    Velocities<double, HYDRAULIC_MODE::STEADY_STATE, HYDRAULIC_RESOLVE::IMPLICIT> velocities(heads);
    velocities.setDomain(20.0, 20.0);

    RowMajMat<double> perm = RowMajMat<double>::Constant(rows, cols, 1e-4);
    velocities.setPermKX(perm);
    velocities.setPermKY(perm);

    // 2. Solute concentrations and advection
    RowMajMat<double> concentrations = RowMajMat<double>::Zero(rows, cols);
    Advection<double, HYDRAULIC_MODE::STEADY_STATE, HYDRAULIC_RESOLVE::IMPLICIT> advection(concentrations, velocities);
    advection.setPorosity(0.25);

**Setting of the boundary conditions**

Hydraulic boundary conditions establish the flow field (e.g., constant head on left and right,
closed no-flow on top and bottom). For solute transport, a continuous tracer source is
introduced at the inflow boundary:

.. code-block:: cpp

    // Hydraulic head boundaries: flow from left (15m) to right (10m)
    Boundary<double> &bcHeads = velocities.getBoundary();
    bcHeads.setBoundarySideConstant(BC_SIDE_LEFT, 15.0);
    bcHeads.setBoundarySideConstant(BC_SIDE_RIGHT, 10.0);
    bcHeads.setBoundarySideClosed(BC_SIDE_TOP);
    bcHeads.setBoundarySideClosed(BC_SIDE_BOTTOM);

    // Advection solute boundaries: constant inflow concentration of 1.0 mol/m³
    Boundary<double> &bcAdv = advection.getBoundary();
    bcAdv.setBoundarySideConstant(BC_SIDE_LEFT, 1.0);
    bcAdv.setBoundarySideClosed(BC_SIDE_TOP);
    bcAdv.setBoundarySideClosed(BC_SIDE_BOTTOM);

**Setting of the simulation parameters and simulation start**

The user specifies the desired macroscopic timestep and number of iterations.
If the chosen timestep exceeds the Courant-Friedrichs-Lewy (CFL) stability limit,
``Advection`` automatically subdivides each step into stable inner sub-steps.

.. code-block:: cpp

    advection.setTimestep(50.0);
    advection.setIterations(200);
    advection.setOutputConsole(CONSOLE_OUTPUT_VERBOSE);

    // Automatically solves steady Darcy velocities and advances advective transport
    advection.run();
