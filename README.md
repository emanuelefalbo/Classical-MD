# Description 
This repository contains a small code to perform classical MD using the Velocity Verlet algorithm.
It performs MD through NVE ensemble by using periodic boundary condition (PBC). 
Only Lennard-Jones and Harmonic potentials are used to respectively describe 
non-bonded and bonded interactions on atoms. All bonds are treated equally at moment, which 
means that their equilibrium position is the same. Thus, the program works fairly well for small molecular systems 
such as small molecules or polymer chains. However, in some cases, these structures might need to be minimized
before running the MD to avoid unphysical contatcs between atoms. 


# Compile & Use 
The code is written in FORTRAN95, the main code is the program `md_nve_bodies.f95` which uses the module
`md_module_bodies.f95`. Three files need to be present in working directory the `md.in` file, the coordinates of the systems
in an xyz-formatted file, and a file named `conncetion.txt` containing the bound atoms. The `md.in` file modules the criteria for the dynamics: dt is the time step,
nstep the total number of steps, nout print nout times the coordinates and thermodynamic quantities, decll is the lenght of the cubic cell,
 and rattle = 1 turns on 

takes in a xyz formatted file containing the coordinates, and it returns
the trajectory.xyz and result.dat files. The former contains the dynamics along the time, and thus the xyz-formatted coordinates, 
while result.dat contains in the respecitvely order: the current step, time, temperature, kinetic energy, potential energy, and total energy. 


can  be simply compiled by any compatible fortran compile, for example, one way is: 

`gfortran -o md md_module_bodies.f95 md_nve_bodies.f95`

The
