!*************************************************************************************************
!
!               Molecular Dynamics 
!
!  Program:      MD
!
!  Programmer:   Emanuele Falbo
!                School of Natural and Environmental Science
!                Newcastle University 
!                Newcastle upon Tyne NE1 7RU
!
!  Date:         October 2020 
!
!  Language:     Fortran-95
!
!  Description:  The program performs a molecular dynamics (MD) from coordinate file xyz type
!                The Verlet algorithm is used throughout. The potential only accounts for 
!                Lennard-Jones non-bonded potential and harmonic potential for bonded atoms 
!
!**********************************************************************************************

!**********************************************************************************************
!  Main program
!**********************************************************************************************


program main

 use iso_fortran_env
 use config_io_module
 use md_module
  
  implicit none
    
   integer, parameter :: dp = selected_real_kind(15, 307)
   character(len=*), parameter :: input_file="md_input.txt"
   character(len=100),dimension(:),allocatable :: label
   real(dp),dimension(:,:),allocatable :: r
   real(dp),dimension(:,:),allocatable :: rold
   real(dp),dimension(:,:),allocatable :: f
   real(dp),dimension(:,:),allocatable :: g
   real(dp),dimension(:,:),allocatable :: v
   integer,dimension(:),allocatable :: vec1,vec2
   real(dp) :: t,mass
   real(dp) :: dt
   real(dp) :: ekin,epot,etot,etot2,rms
   real(dp) :: temp,pressure
   real(dp) :: epotLJ,epot_pair
   real(dp) :: v2
   real(dp) :: Rcut
   real(dp) :: dcell
   real(dp) :: Req, keq                                    ! Req for Bond
   real(dp) :: lj_eps, lj_sigma                            ! LJ paramters
   real :: finish,start,time
   integer:: i,j,iatm,Natm,Nstep,nout
   integer:: rattle, pbc
   integer:: ios
   integer:: fu=20, fu1=21, fu2=31
   integer:: No

   namelist /ctrl/ dt, Nstep, dcell, Rcut, nout, rattle, pbc

   ! Print out  Initial Configuration
   write(*,fmt='(a)') 'md_nve_lj'
   write(*,fmt='(a)') 'Molecular dynamics, constant-NVE ensemble'
   write(*,fmt='(a)') 'Particle mass=1 throughout'   
   write(*,fmt='(a)') ' '

   ! Default Values
   dt=0.002
   Nstep=1000
   dcell=20.0
   Rcut=8.0
   nout=1
   pbc=1
   rattle=0

   call read_nml(input_file)           ! Read in input data

   ! Assign variable
   t=0.0
   mass=1.0
   call read_ff(Req,keq,lj_sigma,lj_eps)                            ! Read in FF parameters
   call read_natoms(Natm)
   call allocate_arrays(Natm,label,r,rold,v,f,g,dcell,Rcut) 
   call read_coords(Natm,label,r)
   call read_connection(vec1,vec2,No)

    ! MD Routine
    open(fu1,file="trajectory.xyz", status="replace")    
    open(fu2,file="md_output.txt", status="replace")    
    
    call cpu_time(start)  ! Time MD engine 

    do i=1,Nstep
      
       !Initialize force and energies
       call forces_LJ(r,f,Natm,dcell,epotLJ,Rcut,pbc,lj_sigma,lj_eps)
       call forces_spring(r,g,Natm,dcell,vec1,vec2,No,epot_pair,pbc,Req,keq)
       epot=epot_pair+epotLJ
       ekin=0.d0
      
       ! Verlet loop 
       do iatm=1,Natm
          rold(:,iatm)=r(:,iatm)  ! Store the nth-1 step coordinates
       end do

       do iatm=1,Natm
          v(:,iatm)=v(:,iatm)+0.5d0*dt*(f(:,iatm) + g(:,iatm))   !  Half-Kick
          r(:,iatm)=r(:,iatm)+v(:,iatm)*dt                       ! Drift rth
       end do  
    
       ! Turn on PBC 
       if (pbc==1) then 
           call set_pbc(r,Natm,dcell)
       end if
       ! Turn on RATTLE
       if (rattle==1) then
           call  rattle_a(r,rold,v,Natm,vec1,vec2,No,dcell)
       end if

       call forces_LJ(r,f,Natm,dcell,epotLJ,Rcut,pbc,lj_sigma,lj_eps)
       call forces_spring(r,g,Natm,dcell,vec1,vec2,No,epot_pair,pbc,Req,keq)
          
       do iatm=1,Natm
         v(:,iatm)=v(:,iatm)+0.5d0*dt*(f(:,iatm) + g(:,iatm))   ! Half-kick
          !call rattle_b(r,v,Natm)
          ! Compute energies
          v2=sum(v(:,iatm)**2) 
          ekin=ekin+(0.5d0*v2)
       end do
    
       ! Compute properties & Print every nout
       call thermo_quantities(dcell,Natm,ekin,epot,etot,temp,pressure)

       if (mod(i,nout)==0) then
           call print_output(Natm,label,r,fu1,fu2,nout,i,t,etot,epot,temp,ekin,epot_pair,epotLJ)
       end if
    
       t=t+dt

    end do

    close(fu1)
    close(fu2)
   
    call cpu_time(finish)
    time=finish-start
    print*, "Time = ", time, " seconds"      
    
end program main
