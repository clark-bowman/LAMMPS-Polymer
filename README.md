# LAMMPS-Polymer

## polymer_fluid.cpp
C++ script which generates a LAMMPS script `polymer_fluid.in` and data file `polymer.dat`. Parameters are read in from `params_fluid.txt`.

## params_fluid.txt
Variable parameters for DPD simulation. Each parameter is the only entry on a line. For the order and description of parameters, see the commented parameter input section of `polymer_fluid.cpp`.

## polymer_fluid.in
LAMMPS script modeling polymer in a nanofluidic channel with DPD. `polymer.dat` must be in the same directory to run. This script runs with the version of LAMMPS compiled Feb 1, 2014.

_Outputs_:
* `xpos.out` - Average polymer x position over time
* `xmin.out` - Minimum polymer x position (i.e., position of the leftmost monomer) over time
* `xmax.out` - Maximum polymer x position (i.e., position of the rightmost monomer) over time
* `xsq.out` - Polymer [radius of gyration](http://lammps.sandia.gov/doc/compute_gyration.html) over time
* `xvel.out` - Average polymer x velocity over time
* `xfor.out` - Average polymer x force over time (i.e., average force applied to a monomer)
* `~.lammpstrj` - LAMMPS trajectory file for visualization. Loadable in, e.g., [VMD](http://www.ks.uiuc.edu/Research/vmd/)

Each output has two commented header lines, followed by a series of data in the form `timestep data` on each line. Timesteps are *not* the same as time units (convert with the timestep set in the script, default 0.001).

## polymer_nofluid.cpp
C++ script which generates a LAMMPS script `polymer_nofluid.in` and data file `polymer.dat`. Parameters are read in from `params_nofluid.txt`.

## params_nofluid.txt
Variable parameters for BD simulation. Each parameter is the only entry on a line. For the order and description of parameters, see the commented parameter input section of `polymer_nofluid.cpp`.

## polymer_nofluid.in
LAMMPS script modeling polymer in a nanofluidic channel with BD. `polymer.dat` must be in the same directory to run. This script runs with the version of LAMMPS compiled Feb 25, 2015.

_Outputs_:
* `xpos.out` - Average polymer x position over time
* `xmin.out` - Minimum polymer x position (i.e., position of the leftmost monomer) over time
* `xmax.out` - Maximum polymer x position (i.e., position of the rightmost monomer) over time
* `xsq.out` - Polymer [radius of gyration](http://lammps.sandia.gov/doc/compute_gyration.html) over time
* `~.lammpstrj` - LAMMPS trajectory file for visualization. Loadable in, e.g., [VMD](http://www.ks.uiuc.edu/Research/vmd/)

Each output has two commented header lines, followed by a series of data in the form `timestep data` on each line. Timesteps are *not* the same as time units (convert with the timestep set in the script, default 0.001).

## polymer.dat
[Data file](http://lammps.sandia.gov/doc/read_data.html) storing initial configuration and basic data such as box size, atom types, etc. for the polymer. Loaded into both LAMMPS scripts.
