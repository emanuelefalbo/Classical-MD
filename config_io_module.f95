module config_io_module

 use iso_fortran_env

 contains 

 subroutine read_natoms(Natm)
  implicit none
  character(len=100) :: filename
  integer :: unit_out=19
  integer :: ios
  integer :: Natm

    !call get_command_argument(1,filename)
    !open(unit_out,file=trim(filename))
    open(unit_out,file='coordinates.xyz') 
    read(unit_out,*,iostat=ios) Natm
    write(*,*) "Number of atoms = ",  Natm
    close(unit_out)

  end subroutine read_natoms 

  subroutine read_coords(Natm,label,r)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    character(len=100), dimension(Natm) :: label
    character(len=100) :: filename
    real(dp), dimension(3,Natm) :: r
    integer :: unit_out=18
    integer :: ios
    integer :: iatm
    integer :: Natm

    !call get_command_argument(1,filename)
    open(unit_out,file='coordinates.xyz')
    !open(unit_out,file=trim(filename))
    read(unit_out,*,iostat=ios)
    read(unit_out,*,iostat=ios)
    do iatm=1,Natm
       read(unit_out,*,iostat=ios) label(iatm), r(:,iatm)
    end do
    close(unit_out)

  end subroutine 
   
  subroutine read_connection(vec1,vec2,Nbond)
  implicit none
  integer, dimension(:), allocatable :: vec1,vec2
  integer :: unit_out=55
  integer :: ios
  integer :: i
  integer :: Nbond

   open(unit_out,file="connection.txt",status='old')
   read(unit_out,*) Nbond
   write(*,*) " Number of bonds = ", Nbond
   allocate(vec1(Nbond),vec2(Nbond))
   do i=1,Nbond
     read(unit_out,*,iostat=ios) vec1(i), vec2(i) !connect(:,i)
   end do
   close(unit_out)

  end subroutine

  !subroutine read_nml(dt,Nstep,nout,dcell,rattle,Rcut,pbc)
  subroutine read_nml(filename)
  
  implicit none
  
   integer, parameter :: dp = selected_real_kind(15, 307)
   character(len=*), intent(in) :: filename
   real(dp) :: dt
   real(dp) :: dcell
   real(dp) :: Rcut
   integer :: Nstep, nout 
   integer :: rattle, pbc ! Rattle turns on Constraint, PBC = periodic boundary conditions
   integer :: ios
   integer :: fu=22
   
   namelist /ctrl/ dt, Nstep, dcell, Rcut, nout, rattle, pbc
 
   open(fu,file=filename, iostat=ios)
   read (fu, nml=ctrl, iostat=ios)
   print*, "dt = ", dt
   print*, "nstep = ", Nstep
   print*, "dcell = ", dcell
   print*, "Rcut = ", Rcut
   print*, "nout = ", nout
   print*, "rattle = ", rattle
   print*, "PBC = ", pbc
 
   if (ios /= 0) then
       print*, 'Invalid Namelist format.'
   end if
   close(fu)
 
   end subroutine
 
 
   subroutine print_output(Natm,label,r,fu1,fu2,nout,i,t,etot,epot,temp,ekin,epot_pair,epotLJ)
   
   implicit none
  
    integer, parameter :: dp = selected_real_kind(15, 307)
    character(len=100), dimension(Natm):: label
    real(dp),dimension(3,Natm):: r
    real(dp) :: t
    real(dp) :: etot,epot
    real(dp) :: temp,ekin
    real(dp) :: epotLJ, epot_pair
    integer :: i, nout 
    integer :: fu1, fu2
    integer :: iatm, Natm
     
    ! Print output every nout steps
       write(fu1,*) Natm
       write(fu1,*) "Cartesian"
  
       do iatm=1,Natm
          write(fu1,*) trim(adjustl(label(iatm))),r(:,(iatm)) !ekin,epot,etot
       end do
  
       write(fu2,"(a)") repeat("=",25) 
       write(fu2,"(a)") ' ' 
       write(fu2,"(a,t10,i10)") 'NSTEP = ', i
       write(fu2,"(a,t10,f15.4)") 'TIME (PS)  = ', t
       write(fu2,"(a,t10,f15.6)") 'ETOT = ', etot
       write(fu2,"(a,t10,f15.6)") 'EPOT = ', epot
       write(fu2,"(a)") ' '
       write(fu2,"(a,t10,f15.6)") 'TEMP = ', temp
       write(fu2,"(a,t10,f15.6)") 'EKIN = ', ekin
       write(fu2,"(a,t10,f15.6)") 'EBOND = ', epot_pair
       write(fu2,"(a,t10,f15.6)") 'EVWD = ', epotLJ
       write(fu2,"(a)") ' ' 
    
    end subroutine 
 
  subroutine read_ff(r_eq,keq,lj_sigma,lj_eps)
 
   implicit none
 
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp) :: r_eq, keq, theta_eq
    real(dp) :: lj_eps, lj_sigma
    integer :: fu=19
    integer :: ios
    integer :: nline
    integer :: var1, var4
    character(len=100) :: label
    character(len=100) :: match1, match2
 
    open(fu, file="ff_parms.txt", status="old")
 
    match1="BOND"
    match2="LJ"
    do
      read(fu,*,iostat=ios) label
      if (ios == iostat_end) exit
      if (label.eq.match1) then
         read(fu,*,iostat=ios) keq, r_eq
      else if (label.eq.match2) then
            read(fu,*,iostat=ios) lj_sigma, lj_eps
      end if
    end do
 
    close(fu)
 
  end subroutine

end module
