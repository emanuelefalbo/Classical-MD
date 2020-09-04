# Description 
This repository contains a small code to perform classical MD using the Velocity Verlet algorithm.
It performs MD through NVE ensemble by using periodic boundary condition (PBC). 
Only Lennard-Jones and Harmonic potentials are used to respectively describe 
non-bonded and bonded interactions on atoms. All bonds are treated equally at moment, which 
means that their equilibrium position is the same. Thus, the program works fairly well for small molecular systems 
such as small molecules or polymer chains. However, in some cases, these structures might need to be minimized
before running the MD to avoid unphysical contatcs between atoms. 


# Compile & Use 
The code is written in FORTRAN95 , and it can simply compiled by any fortran compile, for example, one way is 

