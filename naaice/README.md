This directory contains a concise benchmark designed for validating FPGA
offloading of the Thomas algorithm, primarily employed for solving linear
equation systems structured within a tridiagonal matrix.


# Benchmark Setup

The benchmark defines a domain measuring $1 \text{cm} \times 0.5 \text{cm}$ (easting $\times$ northing),
discretized in a $10 \times 5$ grid. Each grid cell initially
contains a specific concentration. The concentration in the left domain half is set to $6.92023 \times 10^{-7}$, while in the right half to
$2.02396 \times 10^{-8}$, creating an horizontal concentration discontinuity at
the center of the grid. These initial concentrations are read from headerless csv file [init_conc.csv](./init_conc.csv).

A diffusion time step is simulated with the
heterogeneous 2D-ADI approach detailed in the
[ADI_scheme.pdf](../doc/ADI_scheme.pdf) file. The x component of the
diffusion coefficients, read from headerless csv file [alphax.csv](./alphax.csv) ranges from $\alpha = 10^{-9}$ to $10^{-10}$ (distributed randomly), while the
y-component is held constant at $5 \times 10^{-10}$. Closed
boundary conditions are enforced at all domain boundaries, meaning that concentration cannot enter or exit
the system, or in other terms, that the sum of concentrations over the domain must stay constant. The benchmark simulates a single iteration with a
time step ($\Delta t$) of 360 seconds.


# Usage

To generate new makefiles using the `-DTUG_NAAICE_EXAMPLE=ON` option in CMake,
compile the executable, and run it to generate the benchmark output, follow
these steps:

1.  Navigate to your project's build directory.
2.  Run the following CMake command with the `-DTUG_NAAICE_EXAMPLE=ON` option to
    generate the makefiles:
    
        cmake -DTUG_NAAICE_EXAMPLE=ON ..

3.  After CMake configuration is complete, build the `naaice` executable by running `make`:
    
        make naaice

4.  Once the compilation is successful, navigate to the build directory by `cd
       <build_dir>/naaice`

5.  Finally, run the `naaice` executable to generate the benchmark output:
    
        ./naaice


## Output Files


### `Thomas_<n>.csv`

These files contain the values of the tridiagonal coefficient matrix $A$, where:

-   $Aa$ represents the leftmost value,
-   $Ab$ represents the middle value, and
-   $Ac$ represents the rightmost value of one row of the matrix.

Additionally, the corresponding values of the right-hand-side vector $b$ are
provided.

Since the 2D-ADI BTCS scheme processes each row first and then proceeds
column-wise through the grid, each iteration is saved separately in
consecutively numbered files.


### `BTCS_5_10_1.csv`

The result of the simulation, **separated by whitespaces**!

