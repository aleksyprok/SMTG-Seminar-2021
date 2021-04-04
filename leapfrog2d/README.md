# Leapfrog code readme
This Fortran 90 code solves the cold ideal linearised MHD equations using the leapfrog trapezoidal algorithm (see [1]) for more information on how the code works, see Manual/manual.pdf.

The makefile is set up for the gfortran compiler. To compile, type make, this will output an executable called leapfrog.exe which can be run using ./leapfrog.exe on the terminal.

Note that the command make clean will remove all .o, .mod, .exe files and empty the contents of the Data and Figures directories.

control.f90 contains the modules most likely to be edited by the user, such as set_initial_conditions, set_boundary_conditions and control_variables.

diagnotiscs.f90 contains modules associated with outputting the data.

leapfrog.f90 is the main file where the modules are run.

setup.f90 initialises the array and calculates the maximum allowed time step where the CFL condition is satisfied.

shared_data.f90 declares the global variables which are used across the code.

update_variables.f90 carries out a single step in time of the leapfrog algorithm.

Python/surf.py shows how a surface plot video of a variable can be made using ouputted simulation data.

[1] Rickard, G. J., and A. N. Wright, Alfven resonance excitation and fast wave propagation in magnetospheric waveguides, J. Geophys. Res., 99, 13455, 1994.
