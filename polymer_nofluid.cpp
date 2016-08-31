/*
LAMMPS Script Maker - Polymer in Nanofluidic Channel, No Fluid (clark_bowman@brown.edu)
For LAMMPS version Feb. 25, 2014
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

    int v_run_length, v_dump_interval, v_out_interval;
    double v_bond_strength, v_angle_strength, v_lj_epsilon_polymer, v_lj_sigma_polymer, v_lj_cutoff_polymer, v_wall_coefficient, v_wall_radius, v_lj_cutoff, v_lattice_density, v_angle_neutral, v_timestep;

    // Read parameters from text file
    {
        ifstream file_in ("params_nofluid.txt");
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
        file_in >> v_dump_interval; // Visualization dump interval in steps
        file_in >> v_bond_strength; // Harmonic bond strength (LJ units)
        file_in >> v_angle_strength; // Harmonic angle strength (LJ units)
        file_in >> v_angle_neutral; // Angle neutral position (180.0 = straight)
        file_in >> v_lj_epsilon_polymer; // Polymer LJ epsilon
        file_in >> v_lj_sigma_polymer; // Polymer LJ sigma
        file_in >> v_lj_cutoff_polymer; // Polymer LJ cutoff
        file_in >> v_wall_coefficient; // Wall LJ epsilon
        file_in >> v_wall_radius; // Wall LJ sigma
        file_in >> v_lj_cutoff; // Wall LJ cutoff
        file_in >> v_timestep; // Simulation timestep (LJ units)
        file_in >> v_out_interval; // Interval for outputting information (multiple of 10)
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
        ofstream fileOut ("polymer_nofluid.in");
        fileOut.precision(5);

        fileOut << "# LAMMPS Script - Polymer in Nanofluidic Channel, No Fluid (clark_bowman@brown.edu)" << endl;
        fileOut << "# LAMMPS version Feb. 25, 2014" << endl << endl << endl << endl << endl << endl;

        fileOut << "# SEC: This section defines the physical properties of the polymer. All units are LJ unless specified." << endl << endl;
        fileOut << "# 	VAR: harmonic bond coefficient between monomers in-chain" << endl;
        fileOut << "variable bond_strength equal " << v_bond_strength << endl << endl;
        fileOut << "# 	VAR: length of monomer-monomer bonds" << endl;
        fileOut << "variable bond_length equal " << spacing * micrometer << endl << endl;
        fileOut << "# 	VAR: harmonic angle coefficient in-chain" << endl;
        fileOut << "variable angle_strength equal " << v_angle_strength << endl << endl;
        fileOut << "# 	VAR: angle potentials are centered at this position in degrees (180.0 = straight)" << endl;
        fileOut << "variable angle_neutral equal " << v_angle_neutral << endl << endl;
        fileOut << "# 	VAR: LJ epsilon coefficient for monomers (energy scaling)" << endl;
        fileOut << "variable lj_epsilon_polymer equal " << v_lj_epsilon_polymer << endl << endl;
        fileOut << "# 	VAR: LJ sigma coefficient for monomers (distance scaling)" << endl;
        fileOut << "variable lj_sigma_polymer equal " << v_lj_sigma_polymer << endl << endl;
        fileOut << "# 	VAR: cutoff radius for monomers" << endl;
        fileOut << "variable lj_cutoff_polymer equal " << v_lj_cutoff_polymer << endl << endl << endl;

        fileOut << "# SEC: This section defines properties of the LJ reflecting walls." << endl << endl;
        fileOut << "# 	VAR: coefficient epsilon for wall potential (energy scaling)" << endl;
        fileOut << "variable wall_coefficient equal " << v_wall_coefficient << endl << endl;
        fileOut << "# 	VAR: coefficient sigma for wall potential (distance scaling)" << endl;
        fileOut << "variable wall_radius equal " << v_wall_radius << endl << endl;
        fileOut << "# 	VAR: cutoff radius for wall potential" << endl;
        fileOut << "variable lj_cutoff equal " << v_lj_cutoff << endl << endl << endl;

        fileOut << "# SEC: This section defines properties of the fluid." << endl << endl;
        fileOut << "# 	VAR: enforced fluid flow velocity (critical velocity is NOT the same for BD simulation)" << endl;
        fileOut << "variable flow_velocity equal " << flow_velocity * micrometer << endl << endl << endl;

        fileOut << "# SEC: This section defines parameters of the simulation proper." << endl << endl;
        fileOut << "# 	VAR: thermostat temperature" << endl;
        fileOut << "variable temperature equal " << temperature << endl << endl;
        fileOut << "# 	VAR: total number of steps to run (use arbitrarily large number to run until process is killed)" << endl;
        fileOut << "variable run_length equal " << v_run_length << endl << endl;
        fileOut << "# 	VAR: if dump is used, interval between position dumps of entire system" << endl;
        fileOut << "variable dump_interval equal " << v_dump_interval << endl << endl << endl;

        fileOut << "# SEC: This section initializes the simulation." << endl << endl;
        fileOut << "# 	Define units, atom style, log path, and neighbor settings; and read configuration data for the polymer." << endl;
        fileOut << "# 	Configuration data for the polymer in polymer.dat must be generated separately." << endl << endl;
        fileOut << "units lj" << endl;
        fileOut << "atom_style molecular" << endl;
        fileOut << "log polymer_nofluid.log" << endl;
        fileOut << "read_data polymer.dat" << endl;
        fileOut << "neighbor 0.3 bin" << endl;
        fileOut << "neigh_modify delay 5" << endl << endl << endl;

        fileOut << "# 	Define bond, angle, pairwise interactions according to the values of the variables above." << endl;
        fileOut << "# 	For pairwise interactions, type 1 is polymer, type 2 is fluid, and type 3 is tracked fluid." << endl;
        fileOut << "#   For the no-fluid simulation, types 2 and 3 are not used." << endl << endl;
        fileOut << "bond_style harmonic" << endl;
        fileOut << "bond_coeff 1 ${bond_strength} ${bond_length}" << endl;
        fileOut << "angle_style harmonic" << endl;
        fileOut << "angle_coeff 1 ${angle_strength} ${angle_neutral}" << endl;
        fileOut << "dihedral_style none" << endl;
        fileOut << "pair_style lj/cut ${lj_cutoff_polymer}" << endl;
        fileOut << "pair_coeff * * ${lj_epsilon_polymer} ${lj_sigma_polymer}" << endl << endl << endl;

        fileOut << "# SEC: This section initializes the geometry." << endl << endl;
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
        fileOut << "region flow block " << xlo - 1. << " " << xhi + 1. << " " << ylo - 1. << " " << yhi + 1. << " 0 " << zhi << endl << endl << endl;

        fileOut << "# SEC: This section defines LAMMPS computes that will be used in the simulation." << endl << endl;
        fileOut << "# 	This compute finds in each step the average x flow rate of all particles." << endl;
        fileOut << "# 	By subtracting off the target velocity, a restoring force will be imposed." << endl << endl;
        fileOut << "compute atom_vx all property/atom vx" << endl;
        fileOut << "variable diff atom '(v_flow_velocity-c_atom_vx)*0.1'" << endl << endl;
        fileOut << "# 	Intermediate compute used in the following computes." << endl << endl;
        fileOut << "compute atom_x all property/atom xu" << endl << endl;
        fileOut << "# 	Computes for gathering data: average polymer x position, radius of gyration, min/max x value." << endl << endl;
        fileOut << "compute xpos all reduce ave c_atom_x" << endl;
        fileOut << "compute xsq all gyration" << endl;
        fileOut << "compute xmin all reduce min c_atom_x" << endl;
        fileOut << "compute xmax all reduce max c_atom_x" << endl << endl << endl;

        fileOut << "# SEC: This section defines fixes to impose forces in the simulation." << endl << endl;
        fileOut << "# 	NVE integration, with limit for initialization. BD fix imposes Langevin dynamics." << endl;
        fileOut << "# 	Fix 3 is a restoring force to maintain fluid speed in the channel, applied only in the upper channel." << endl;
        fileOut << "# 	Fixes 2, 4-" << npits + 6 << " are wall LJ potentials." << endl << endl;
        fileOut << "fix 1 all nve/limit 0.05" << endl;
        fileOut << "fix bd all langevin ${temperature} ${temperature} 10.0 4578" << endl;
        fileOut << "fix 2 all wall/region plates lj126 ${wall_coefficient} ${wall_radius} ${lj_cutoff}" << endl;
        fileOut << "fix 3 all addforce v_diff 0 0 region flow" << endl;
        for (int i = 0; i < npits + 3; i++)
            fileOut << "fix " << i + 4 << " all wall/region block" << i + 1 << " lj126 ${wall_coefficient} ${wall_radius} ${lj_cutoff}" << endl;
        fileOut << endl << endl;

        fileOut << "# SEC: This section runs the simulation." << endl << endl;
        fileOut << "# 	Simulation timestep (LJ time units)." << endl;
        fileOut << "timestep " << v_timestep << endl << endl;
        fileOut << "# 	How often to output thermo data and first run phase." << endl;
        fileOut << "thermo ${dump_interval}" << endl;
        fileOut << "run 5000" << endl << endl;
        fileOut << "# 	Release limit on integrator and run full simulation." << endl;
        fileOut << "# 	Fixes " << npits + 7 << "-" << npits + 10 << " dump the computes at specified short intervals." << endl;
        fileOut << "# 	Dump 1 is a full simulation dump at the specified interval." << endl << endl;
        fileOut << "unfix 1" << endl;
        fileOut << "fix 1 all nve" << endl;
        fileOut << "fix " << npits + 7 << " all ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xpos file xpos.out" << endl;
        fileOut << "fix " << npits + 8 << " all ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xsq[0] file xsq.out" << endl;
        fileOut << "fix " << npits + 9 << " all ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xmin file xmin.out" << endl;
        fileOut << "fix " << npits + 10 << " all ave/time " << v_out_interval / 10 << " " << v_out_interval / 5 << " " << v_out_interval << " c_xmax file xmax.out" << endl;
        fileOut << "dump 1 all atom ${dump_interval} polymer.lammpstrj" << endl;
        fileOut << "run ${run_length}" << endl;
    }
    return 0;
}
