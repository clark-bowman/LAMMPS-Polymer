/*
LAMMPS Script Maker - Polymer in Nanofluidic Channel, DPD Fluid (clark_bowman@brown.edu)
For LAMMPS version Feb. 1, 2014
*/

#include <fstream>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265359

double randDouble()
{
    return double(rand()) / double(RAND_MAX);
}

int main()
{
    using namespace std;
    srand (79212);

    int natoms, npits;
    double xlo, xhi, ylo, yhi, zlo, zhi, spacing, flow_velocity, border, micrometer, pitspacing, pitwidth, channeldepth, pitdepth, temperature;

    int v_run_length, v_dump_interval, v_out_interval, v_dump_interval_tracking;
    double v_bond_strength, v_angle_strength, v_dpd_conservative_fluid, v_dpd_cutoff_polymer, v_dpd_conservative_polymer, v_dpd_dissipative, v_dpd_cutoff_fluid, v_wall_coefficient, v_wall_radius, v_lj_cutoff, v_lattice_density, v_angle_neutral, v_timestep;

    // Read parameters from text file
    {
        ifstream file_in ("params_fluid.txt");
        file_in >> micrometer; // One micrometer is ___ lj units
        file_in >> npits; // Number of pits
        file_in >> pitspacing; // Micrometers between pits
        file_in >> pitwidth; // Micrometer pit width
        file_in >> channeldepth; // Micrometer channel depth
        file_in >> pitdepth; // Micrometer pit depth
        file_in >> border; // Micrometer border orthogonal to flow
        file_in >> flow_velocity; // Flow velocity (micrometers per lj unit time)
        file_in >> natoms; // Number of atoms in polymer
        file_in >> spacing; // Micrometer bond length in polymer
        file_in >> temperature; // LJ temperature
        file_in >> v_run_length; // Length of run in steps
        file_in >> v_dump_interval; // Whole simulation dump interval in steps (if used)
        file_in >> v_bond_strength; // Harmonic bond strength (LJ units)
        file_in >> v_angle_strength; // Harmonic angle strength (LJ units)
        file_in >> v_angle_neutral; // Angle neutral position (180.0 = straight)
        file_in >> v_wall_coefficient; // Wall LJ epsilon
        file_in >> v_wall_radius; // Wall LJ sigma
        file_in >> v_lj_cutoff; // Wall LJ cutoff
        file_in >> v_timestep; // Simulation timestep (LJ units)
        file_in >> v_out_interval; // Interval for outputting information (multiple of 10)
        file_in >> v_dpd_conservative_polymer; // DPD conservative coefficient (polymer)
        file_in >> v_dpd_cutoff_polymer; //  DPD cutoff radius (polymer)
        file_in >> v_dpd_conservative_fluid; // DPD conservative coefficient (fluid)
        file_in >> v_dpd_dissipative; //  DPD dissipative coefficient (fluid & polymer)
        file_in >> v_dpd_cutoff_fluid; //  DPD cutoff radius (fluid)
        file_in >> v_lattice_density; //  Lattice density (LJ units) for initialization of fluid in channel
        file_in >> v_dump_interval_tracking; // Interval at which to dump trajectory data on tracked particles
    }

    // Compute simulation bounds
    xlo = 0.;
    xhi = micrometer * (pitspacing + pitwidth) * double(npits);
    ylo = 0.;
    yhi = micrometer * (pitwidth + 2. * border);
    zlo = -pitdepth * micrometer;
    zhi = channeldepth * micrometer;

    // Generate polymer initial configuration
    {
        ofstream fileOut ("polymer.dat");
        fileOut.precision(5);
        fileOut << "# Polymer data file" << endl << endl;

        // Header with number of atoms, bonds, angles, dihedrals, and simulation information
        fileOut << natoms << " atoms" << endl;
        fileOut << natoms - 1 << " bonds" << endl;
        fileOut << natoms - 2 << " angles" << endl;
        fileOut << natoms - 3 << " dihedrals" << endl << endl;
        fileOut << "3 atom types" << endl << "1 bond types" << endl;
        fileOut << "1 angle types" << endl << "1 dihedral types" << endl << endl;
        fileOut << xlo << " " << xhi << " xlo xhi" << endl;
        fileOut << ylo << " " << yhi << " ylo yhi" << endl;
        fileOut << zlo - 1 << " " << zhi + 1 << " zlo zhi" << endl << endl;

        // Masses for each particle type
        fileOut << "Masses" << endl << endl;
        fileOut << "1 1.0" << endl << "2 1.0" << endl << "3 1.0" << endl << endl;

        // Loop to create atoms
        fileOut << "Atoms" << endl << endl;

        // Arrays to store locations of atoms in polymer (generated in order)
        double* curx = new double[natoms];
        double* cury = new double[natoms];
        double* curz = new double[natoms];
        double theta, phi;

        // Place the first atom at a random location in the first pit
        curx[0] = micrometer * 0.5 * pitspacing + 0.1 + randDouble() * (micrometer * pitwidth - 0.2);
        cury[0] = border * micrometer + 0.1 + randDouble() * (micrometer * pitwidth - 0.2);
        curz[0] = zlo + 0.1 + randDouble() * (zhi - zlo - 0.2);

        // For each successive atom...
        for (int i = 1; i < natoms; i++)
        {
            // Try to place the next atom at a random location one bond length away
            while (1)
            {
                theta = randDouble() * 2. * PI;
                phi = acos(2. * randDouble() - 1.);
                curx[i] = curx[i - 1] + micrometer * spacing * cos(theta) * sin(phi);
                cury[i] = cury[i - 1] + micrometer * spacing * sin(theta) * sin(phi);
                curz[i] = curz[i - 1] + micrometer * spacing * cos(phi);

                // If this is a valid location (not outside the channel), OK
                if (curx[i] > micrometer * 0.5 * pitspacing + 0.1 and curx[i] < (pitwidth + 0.5 * pitspacing) * micrometer - 0.1 and cury[i] > border * micrometer + 0.1 and cury[i] < (border + pitwidth) * micrometer - 0.1 and curz[i] > zlo + 0.1 and curz[i] < zhi - 0.1)
                    break;
            }
        }

        // Print out atom locations, bonds, angles, dihedrals
        for (int i = 0; i < natoms; i++)
            fileOut << i + 1 << " 1 1 " << curx[i] << " " << cury[i] << " " << curz[i] << endl;
        fileOut << endl << "Bonds" << endl << endl;
        for (int i = 0; i < natoms - 1; i++)
            fileOut << i + 1 << " 1 " << i + 1 << " " << i + 2 << endl;
        fileOut << endl << "Angles" << endl << endl;
        for (int i = 0; i < natoms - 2; i++)
            fileOut << i + 1 << " 1 " << i + 1 << " " << i + 2 << " " << i + 3 << endl;
        fileOut << endl << "Dihedrals" << endl << endl;
        for (int i = 0; i < natoms - 3; i++)
            fileOut << i + 1 << " 1 " << i + 1 << " " << i + 2 << " " << i + 3 << " " << i + 4 << endl;
    }

    // Create LAMMPS script
    {
        ofstream fileOut ("polymer_fluid.in");
        fileOut.precision(5);

        fileOut << "# LAMMPS Script - Polymer in Nanofluidic Channel, DPD Fluid (clark_bowman@brown.edu)" << endl;
        fileOut << "# LAMMPS version Feb. 1, 2014" << endl << endl << endl << endl << endl << endl;

        fileOut << "# SEC: This section defines the physical properties of the polymer. All units are LJ unless specified." << endl << endl;
        fileOut << "# 	VAR: harmonic bond coefficient between monomers in-chain" << endl;
        fileOut << "variable bond_strength equal " << v_bond_strength << endl << endl;
        fileOut << "# 	VAR: length of monomer-monomer bonds" << endl;
        fileOut << "variable bond_length equal " << spacing * micrometer << endl << endl;
        fileOut << "# 	VAR: harmonic angle coefficient in-chain" << endl;
        fileOut << "variable angle_strength equal " << v_angle_strength << endl << endl;
        fileOut << "# 	VAR: angle potentials are centered at this position in degrees (180.0 = straight)" << endl;
        fileOut << "variable angle_neutral equal " << v_angle_neutral << endl << endl;
        fileOut << "# 	VAR: DPD conservative coefficient for monomers" << endl;
        fileOut << "variable dpd_conservative_polymer equal " << v_dpd_conservative_polymer << endl << endl;
        fileOut << "# 	VAR: DPD cutoff radius for monomers" << endl;
        fileOut << "variable dpd_cutoff_polymer equal " << v_dpd_cutoff_polymer << endl << endl << endl;

        fileOut << "# SEC: This section defines properties of the LJ reflecting walls." << endl << endl;
        fileOut << "# 	VAR: coefficient epsilon for wall potential (energy scaling)" << endl;
        fileOut << "variable wall_coefficient equal " << v_wall_coefficient << endl << endl;
        fileOut << "# 	VAR: coefficient sigma for wall potential (distance scaling)" << endl;
        fileOut << "variable wall_radius equal " << v_wall_radius << endl << endl;
        fileOut << "# 	VAR: cutoff radius for wall potential" << endl;
        fileOut << "variable lj_cutoff equal " << v_lj_cutoff << endl << endl << endl;

        fileOut << "# SEC: This section defines properties of the fluid." << endl << endl;
        fileOut << "# 	VAR: DPD conservative coefficient for fluid" << endl;
        fileOut << "variable dpd_conservative_fluid equal " << v_dpd_conservative_fluid << endl << endl;
        fileOut << "# 	VAR: DPD dissipative coefficient (used for both fluid and polymer)" << endl;
        fileOut << "variable dpd_dissipative equal " << v_dpd_dissipative << endl << endl;
        fileOut << "# 	VAR: DPD cutoff radius for fluid particles" << endl;
        fileOut << "variable dpd_cutoff_fluid equal " << v_dpd_cutoff_fluid << endl << endl;
        fileOut << "# 	VAR: lattice density of fluid initialization in particles per cubic unit" << endl;
        fileOut << "variable lattice_density equal " << v_lattice_density << endl << endl;
        fileOut << "# 	VAR: enforced fluid flow velocity (critical velocity was around 0.05)" << endl;
        fileOut << "variable flow_velocity equal " << flow_velocity * micrometer << endl << endl << endl;

        fileOut << "# SEC: This section defines parameters of the simulation proper." << endl << endl;
        fileOut << "# 	VAR: thermostat temperature" << endl;
        fileOut << "variable temperature equal " << temperature << endl << endl;
        fileOut << "# 	VAR: total number of steps to run (use arbitrarily large number to run until process is killed)" << endl;
        fileOut << "variable run_length equal " << v_run_length << endl << endl;
        fileOut << "# 	VAR: interval between position dumps of tracked subset of particles" << endl;
        fileOut << "variable dump_interval_tracking equal " << v_dump_interval_tracking << endl << endl;
        fileOut << "# 	VAR: if dump is used, interval between position dumps of entire system" << endl;
        fileOut << "variable dump_interval equal " << v_dump_interval << endl << endl << endl;

        fileOut << "# SEC: This section initializes the simulation." << endl << endl;
        fileOut << "# 	Define units, atom style, log path, and neighbor settings; and read configuration data for the polymer." << endl;
        fileOut << "# 	Configuration data for the polymer in polymer.dat must be generated separately." << endl;
        fileOut << "# 	Final line communicates ghost data and is necessary for DPD parallelizing. On some versions of LAMMPS, use instead `comm_modify vel yes'" << endl << endl;
        fileOut << "units lj" << endl;
        fileOut << "atom_style molecular" << endl;
        fileOut << "log polymer_fluid.log" << endl;
        fileOut << "read_data polymer.dat" << endl;
        fileOut << "neighbor 0.3 bin" << endl;
        fileOut << "neigh_modify delay 5" << endl;
        fileOut << "communicate single vel yes" << endl << endl << endl;

        fileOut << "# 	Define bond, angle, pairwise interactions according to the values of the variables above." << endl;
        fileOut << "# 	For pairwise interactions, type 1 is polymer, type 2 is fluid, and type 3 is tracked fluid." << endl;
        fileOut << "#   35498 is a temperature seed; change for different randomness." << endl << endl;
        fileOut << "bond_style harmonic" << endl;
        fileOut << "bond_coeff 1 ${bond_strength} ${bond_length}" << endl;
        fileOut << "angle_style harmonic" << endl;
        fileOut << "angle_coeff 1 ${angle_strength} ${angle_neutral}" << endl;
        fileOut << "dihedral_style none" << endl;
        fileOut << "pair_style dpd ${temperature} ${dpd_cutoff_fluid} 35498" << endl;
        fileOut << "pair_coeff 1 * ${dpd_conservative_polymer} ${dpd_dissipative} ${dpd_cutoff_polymer}" << endl;
        fileOut << "pair_coeff 2*3 2*3 ${dpd_conservative_fluid} ${dpd_dissipative} ${dpd_cutoff_fluid}" << endl << endl << endl;

        fileOut << "# SEC: This section initializes the geometry." << endl << endl;
        fileOut << "#   Define safe regions where fluid may be placed. This is the union of safe1 (main channel) with safe2 - safe4 (three pits)." << endl;
        fileOut << "# 	The simulation box is " << xhi << " x " << yhi << " x " << zhi - zlo + 2. << ". The safe regions occur just inside the locations of walls to prevent singularities at the boundaries." << endl << endl;
        fileOut << "region safe1 block " << xlo + 0.1 << " " << xhi << " " << ylo + 0.1 << " " << yhi << " 0.1 " << zhi - 0.1 << " units box" << endl;
        for (int i = 0; i < npits; i++)
            fileOut << "region safe" << i + 2 << " block " << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) + 0.1 << " " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) - 0.1 << " " << ylo + micrometer * border + 0.1 << " " << yhi - micrometer * border - 0.1 << " " << zlo + 0.1 << " " << zhi - 0.2 << " units box" << endl;
        fileOut << "region safe union " << npits + 1;
        for (int i = -1; i < npits; i++)
            fileOut << " safe" << i + 2;
        fileOut << endl << endl;

        fileOut << "# 	Define unsafe regions where LJ walls will be imposed. This must be done as a series of boxes because of limitations of fix wall/region," << endl;
        fileOut << "# 	which does not handle internal corners well. `plates' is the upper and lower boundary, and the other blocks shape out the channel." << endl;
        fileOut << "# 	Final region defines where the force from the fluid will be applied." << endl << endl;
        fileOut << "region plates block " << xlo - 1. << " " << xhi + 1. << " " << ylo - 1. << " " << yhi + 1. << " " << zlo - 0.1 << " " << zhi + 0.1 << " units box side in" << endl;
        fileOut << "region block1 block " << xlo - 1. << " " << xhi + 1. << " " << ylo - 1. << " " << ylo + micrometer * border << " " << zlo - 1. << " 0 units box side out" << endl;
        fileOut << "region block2 block " << xlo - 1. << " " << xhi + 1. << " " << yhi - micrometer * border << " " << yhi + 1. << " " << zlo - 1. << " 0 units box side out" << endl;
        fileOut << "region block3 block " << xlo - 1. << " " << xlo + 0.5 * micrometer * pitspacing << " " << ylo + micrometer * border << " " << yhi - micrometer * border << " " << zlo - 1. << " 0 units box side out" << endl;
        fileOut << "region block4 block " << xhi - 0.5 * micrometer * pitspacing << " " << xhi + 1. << " " << ylo + micrometer * border << " " << yhi - micrometer * border << " " << zlo - 1. << " 0 units box side out" << endl;
        for (int i = 0; i < npits - 1; i++)
            fileOut << "region block" << i + 5 << " block " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) << " " << xlo + micrometer * (1.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) << " " << ylo + micrometer * border << " " << yhi - micrometer * border << " " << zlo - 1. << " 0 units box side out" << endl;
        fileOut << "region flow block " << xlo - 1. << " " << xhi + 1. << " " << ylo - 1. << " " << yhi + 1. << " 0 " << zhi << " units box" << endl << endl << endl;

        fileOut << "# SEC: This section initializes the fluid." << endl << endl;
        fileOut << "# 	At a lattice of the specified density, create the fluid atoms." << endl;
        fileOut << "# 	Define regions and groups to track. Tracked region is a section of fluid." << endl << endl;
        fileOut << "lattice fcc ${lattice_density}" << endl;
        fileOut << "create_atoms 2 region safe" << endl;
        fileOut << "region tracked block 60.3 60.5 -1 " << yhi + 1 << " " << zlo - 1 << " " << zhi + 1 << " units box" << endl;
        fileOut << "group topflow region flow" << endl;
        fileOut << "group polymer type 1" << endl;
        fileOut << "group fluid type 2" << endl;
        fileOut << "group stuff type 1 2" << endl << endl << endl;

        fileOut << "# SEC: This section introduces partial-slip boundary conditions by the method of Pivkin & Karniadakis (2004)." << endl << endl;
        fileOut << "# 	Each section here creates one section of frozen fluid particles just outside the wall." << endl;
        fileOut << "# 	Variable `v1' is the inter-particle spacing necessary to achieve the desired lattice density." << endl;
        fileOut << "# 	Particles are placed in a small region just outside the wall on the vertices of an fcc lattice." << endl;
        fileOut << "# 	The lattice, which by default has a particle on the origin, is shifted by variables `v2' onward so that a layer is centered in the target region." << endl << endl;

        double lo_outside = pow(4./v_lattice_density,1./3) * 0.5 / sqrt(2) - 0.1;
        double hi_outside = pow(4./v_lattice_density,1./3) * 0.5 / sqrt(2) + 0.1;

        fileOut << "# 	Top wall particles." << endl;
        fileOut << "variable v1 equal '(4.0/v_lattice_density)^(1/3)'" << endl;
        fileOut << "variable v2 equal '" << zhi << "/v_v1+0.5/sqrt(2)-floor(" << zhi << "/v_v1+0.5/sqrt(2))'" << endl;
        fileOut << "lattice fcc ${lattice_density} origin 0 0 ${v2}" << endl;
        fileOut << "region frozen1 block -1 " << xhi + 1 << " -1 " << yhi + 1 << " " << zhi + lo_outside << " " << zhi + hi_outside << " units box side in" << endl;
        fileOut << "create_atoms 2 region frozen1" << endl << endl;

        fileOut << "# 	Bottom wall particles (pit)." << endl;
        fileOut << "variable v3 equal '" << zlo << "/v_v1-0.5/sqrt(2)-floor(" << zlo << "/v_v1-0.5/sqrt(2))'" << endl;
        fileOut << "lattice fcc ${lattice_density} origin 0 0 ${v3}" << endl;
        for (int i = 0; i < npits; i++)
            fileOut << "region frozen2" << i << " block " << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) - 1 << " " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) + 1 << " " << ylo + micrometer * border - 1 << " " << yhi - micrometer * border + 1 << " " << zlo - hi_outside << " " << zlo - lo_outside << " units box side in" << endl;
        fileOut << "region frozen2 union " << npits;
        for (int i = 0; i < npits; i++)
            fileOut << " frozen2" << i;
        fileOut << endl;
        fileOut << "create_atoms 2 region frozen2" << endl << endl;

        fileOut << "# 	Bottom wall particles (channel)." << endl;
        fileOut << "variable v4 equal '-0.5/sqrt(2)+1'" << endl;
        fileOut << "lattice fcc ${lattice_density} origin 0 0 ${v4}" << endl;
        fileOut << "region frozen3a block -1 " << xlo + micrometer * (0.5 * pitspacing) - 1 << " " << ylo + micrometer * border - 2 << " " << yhi - micrometer * border + 2 << " " << -hi_outside << " " << -lo_outside << " units box side in" << endl;
        fileOut << "region frozen3b block " << xhi - micrometer * (0.5 * pitspacing) + 1 << " " << xhi + 1 << " " << ylo + micrometer * border - 2 << " " << yhi - micrometer * border + 2 << " " << -hi_outside << " " << -lo_outside << " units box side in" << endl;
        fileOut << "region frozen3c block -1 " << xhi + 1 << " -1 " << ylo + micrometer * border - 1 << " " << -hi_outside << " " << -lo_outside << " units box side in" << endl;
        fileOut << "region frozen3d block -1 " << xhi + 1 << " " << yhi - micrometer * border + 1 << " " << yhi + 1 << " " << -hi_outside << " " << -lo_outside << " units box side in" << endl;
        for (int i = 0; i < npits - 1; i++)
            fileOut << "region frozen3" << i << " block " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) + 1 << " " << xlo + micrometer * (1.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) - 1 << " " << ylo + micrometer * border - 2 << " " << yhi - micrometer * border + 2 << " " << -hi_outside << " " << -lo_outside << " units box side in" << endl;
        fileOut << "region frozen3 union " << npits + 3 << " frozen3a frozen3b frozen3c frozen3d";
        for (int i = 0; i < npits - 1; i++)
            fileOut << " frozen3" << i;
        fileOut << endl;
        fileOut << "create_atoms 2 region frozen3" << endl << endl;

        fileOut << "# 	Particles in left walls of pits." << endl;
        for (int i = 0; i < npits; i++)
        {
            fileOut << "variable v5" << i << " equal '" << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) << "/v_v1-0.5/sqrt(2)-floor(" << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) << "/v_v1-0.5/sqrt(2))'" << endl;
            fileOut << "lattice fcc ${lattice_density} origin ${v5" << i << "} 0 0" << endl;
            fileOut << "region frozen4" << i << " block " << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) - hi_outside << " " << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) - lo_outside << " " << ylo + micrometer * border - 1 << " " << yhi - micrometer * border + 1 << " " << zlo << " -1 units box side in" << endl;
            fileOut << "create_atoms 2 region frozen4" << i << endl;
        }
        fileOut << endl;

        fileOut << "# 	Particles in right walls of pits." << endl;
        for (int i = 0; i < npits; i++)
        {
            fileOut << "variable v6" << i << " equal '" << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) << "/v_v1+0.5/sqrt(2)-floor(" << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) << "/v_v1+0.5/sqrt(2))'" << endl;
            fileOut << "lattice fcc ${lattice_density} origin ${v6" << i << "} 0 0" << endl;
            fileOut << "region frozen5" << i << " block " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) + lo_outside << " " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) + hi_outside << " " << ylo + micrometer * border - 1 << " " << yhi - micrometer * border + 1 << " " << zlo << " -1 units box side in" << endl;
            fileOut << "create_atoms 2 region frozen5" << i << endl;
        }
        fileOut << endl;

        fileOut << "# 	Particles in front walls of pits." << endl;
        fileOut << "variable v7 equal '" << ylo + micrometer * border << "/v_v1-0.5/sqrt(2)-floor(" << ylo + micrometer * border << "/v_v1-0.5/sqrt(2))'" << endl;
        fileOut << "lattice fcc ${lattice_density} origin 0 ${v7} 0" << endl;
        for (int i = 0; i < npits; i++)
        {
            fileOut << "region frozen6" << i << " block " << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) << " " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) << " " << ylo + micrometer * border - hi_outside << " " << ylo + micrometer * border - lo_outside << " " << zlo << " -1 units box side in" << endl;
            fileOut << "create_atoms 2 region frozen6" << i << endl;
        }
        fileOut << endl;

        fileOut << "# 	Particles in back walls of pits." << endl;
        fileOut << "variable v8 equal '" << yhi - micrometer * border << "/v_v1+0.5/sqrt(2)-floor(" << yhi - micrometer * border << "/v_v1+0.5/sqrt(2))'" << endl;
        fileOut << "lattice fcc ${lattice_density} origin 0 ${v8} 0" << endl;
        for (int i = 0; i < npits; i++)
        {
            fileOut << "region frozen7" << i << " block " << xlo + micrometer * (0.5 * pitspacing + i * (pitspacing + pitwidth)) << " " << xlo + micrometer * (0.5 * pitspacing + pitwidth + i * (pitspacing + pitwidth)) << " " << yhi - micrometer * border + lo_outside << " " << yhi - micrometer * border + hi_outside << " " << zlo << " -1 units box side in" << endl;
            fileOut << "create_atoms 2 region frozen7" << i << endl;
        }
        fileOut << endl;


        fileOut << "# SEC: This section defines LAMMPS computes that will be used in the simulation." << endl << endl;
        fileOut << "# 	This compute finds in each step the average x flow rate of all particles in the upper channel." << endl;
        fileOut << "# 	By subtracting off the target velocity, a restoring force will be imposed." << endl << endl;
        fileOut << "compute av_vx fluid reduce/region flow ave vx" << endl;
        fileOut << "variable diff equal '10.0*(v_flow_velocity - c_av_vx)'" << endl << endl;
        fileOut << "# 	Intermediate computes used in the following computes." << endl << endl;
        fileOut << "compute atom_x polymer property/atom xu" << endl;
        fileOut << "compute atom_vx polymer property/atom vx" << endl;
        fileOut << "compute atom_fx polymer property/atom fx" << endl << endl;
        fileOut << "# 	Computes for gathering data: average polymer x position, radius of gyration, min/max x value, average x velocity, average x force." << endl << endl;
        fileOut << "compute xpos polymer reduce ave c_atom_x" << endl;
        fileOut << "compute xsq polymer gyration" << endl;
        fileOut << "compute xmin polymer reduce min c_atom_x" << endl;
        fileOut << "compute xmax polymer reduce max c_atom_x" << endl;
        fileOut << "compute xvel polymer reduce ave c_atom_vx" << endl;
        fileOut << "compute xfor polymer reduce ave c_atom_fx" << endl << endl << endl;

        fileOut << "# SEC: This section initializes the particle velocities for the fluid and upper channel." << endl << endl;
        fileOut << "velocity fluid create ${temperature} 12314" << endl;
        fileOut << "velocity topflow set ${flow_velocity} 0 0 sum yes" << endl << endl << endl;

        fileOut << "# SEC: This section defines fixes to impose forces in the simulation." << endl << endl;
        fileOut << "# 	NVE integration for fluid, with limit for initialization." << endl;
        fileOut << "# 	Fix 3 is a restoring force to maintain fluid speed in the channel." << endl;
        fileOut << "# 	Fixes 2, 4-" << npits + 6 << " are wall LJ potentials." << endl << endl;
        fileOut << "fix 1 fluid nve/limit 0.05" << endl;
        fileOut << "fix 2 stuff wall/region plates lj126 ${wall_coefficient} ${wall_radius} ${lj_cutoff}" << endl;
        fileOut << "fix 3 fluid addforce v_diff 0 0 region flow" << endl;
        for (int i = 0; i < npits + 3; i++)
            fileOut << "fix " << i + 4 << " stuff wall/region block" << i + 1 << " lj126 ${wall_coefficient} ${wall_radius} ${lj_cutoff}" << endl;
        fileOut << endl << endl;

        fileOut << "# SEC: This section runs the simulation." << endl << endl;
        fileOut << "# 	Simulation timestep (LJ time units)." << endl;
        fileOut << "timestep " << v_timestep << endl << endl;
        fileOut << "# 	How often to output thermo data and first run phase." << endl;
        fileOut << "thermo ${dump_interval_tracking}" << endl;
        fileOut << "run 25000" << endl << endl;
        fileOut << "# 	Release limit on fluid integrator and begin NVE integration on the polymer." << endl;
        fileOut << "unfix 1" << endl;
        fileOut << "fix 1 fluid nve" << endl;
        fileOut << "fix 1a polymer nve/limit 0.05" << endl;
        fileOut << "run 25000" << endl << endl;

        fileOut << "# 	At this time, collect a strip of fluid particles near the polymer into a group to be tracked." << endl;
        fileOut << "# 	Convert fluid particles in the strip to type 3 (tracked fluid)." << endl << endl;
        fileOut << "group grp_tracked_pre region tracked" << endl;
        fileOut << "group grp_tracked intersect fluid grp_tracked_pre" << endl;
        fileOut << "set group grp_tracked type 3" << endl;
        fileOut << "group grp_tracked type 1" << endl << endl;

        fileOut << "# 	Fixes for dumping information. Fixes " << npits + 7 << "-" << npits + 12 << " dump the computes at specified short intervals." << endl;
        fileOut << "# 	Dump 1 is a full simulation dump at the specified interval (leave commented)." << endl;
        fileOut << "# 	Dump 2 is a tracked particle dump at the specified interval." << endl << endl;
        fileOut << "fix " << npits + 7 << " polymer ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xpos file xpos.out" << endl;
        fileOut << "fix " << npits + 8 << " polymer ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xsq[0] file xsq.out" << endl;
        fileOut << "fix " << npits + 9 << " polymer ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xmin file xmin.out" << endl;
        fileOut << "fix " << npits + 10 << " polymer ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xmax file xmax.out" << endl;
        fileOut << "fix " << npits + 11 << " polymer ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xvel file xvel.out" << endl;
        fileOut << "fix " << npits + 12 << " polymer ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xfor file xfor.out" << endl;
        fileOut << "#dump 1 all atom ${dump_interval} polymer.lammpstrj" << endl;
        fileOut << "dump 2 grp_tracked atom ${dump_interval_tracking} tracked_polymer.lammpstrj" << endl << endl;
        fileOut << "# 	Run data-collecting simulation." << endl;
        fileOut << "run ${run_length}" << endl;
    }
    return 0;
}
