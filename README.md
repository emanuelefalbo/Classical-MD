# Description 
This repository contains a small code to perform classical MD using the Velocity Verlet algorithm.
It performs MD through NVE ensemble by using periodic boundary condition (PBC). 
Only Lennard-Jones and Harmonic potentials are used to respectively describe 
non-bonded and bonded interactions on atoms. All bonds are treated equally at moment, which 
means that their equilibrium position is the same. Thus, the program works fairly well for small molecular systems 
such as small molecules or polymer chains. However, in some cases, these structures might need to be minimized
before running the MD to avoid unphysical contatcs between atoms. 

# Compilation

The program can be simply compiled by any compatible fortran compiler, despite only gfortran has been tested: 

`gfortran -o md config_io_module.f95 md_module_bodies.f95 md_nve_bodies.f95`

However, gfortran-x compiler with version x greater than 4.8 are advisable.  

# Use 
The code is written in FORTRAN95, the main code is the program `md_nve_bodies.f95` which uses the module
`md_module_bodies.f95`, and `config_io_modulef.95`. Three files need to be present in working directory: the `md_input.txt` file containing the settings for the dynamics, the `ff_parm.txt` containing the force field parameters for the system (hamronic forces + lennard-jones values), and a file named `connection.txt` containing the bound atoms. The `md_input.txt` file modules the criteria for the dynamics: dt = the time step,
nstep = total number of steps, nout = prints nout times the coordinates and thermodynamic quantities, decll = the lenght of the cubic cell,
and rattle = 1 (0) turns on (off) the contraint on the atoms in the `connection.txt`. This last file needs to be a two columns file, containing bonds informations. For example, if bonds are between the atoms indexed as 1 and 3 , and 2 and 4, the first column can be written as:

2 <br/>
1 3 <br/>
2 4


Before running the program, the`connection.txt` must be created. This is is done by running the program `$ ./build file.xyz`.
which convert the file.xyz into file.mol2 with babel and read the conncetion information from this latter. 
The proram can the be run as 

`$ ./md file.xyz`

where file.xyz is the system for which the dynamics is performed. The trajectory.xyz and md_output.dat files are returned 
containing respectively the dynamics along the time (xyz-formatted coordinates), and thermochemistry of the system. 

# Working Example

The `chain.xyz` is an atomic chain of four carbon atoms. Its values from the OPLS force field are in the ff_parm.txt file

`$ ./build chain.xyz` <br/>
`$ ./md chain.xyz`

