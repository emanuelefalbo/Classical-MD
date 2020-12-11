module md_module

 implicit none
 
 contains

  subroutine allocate_arrays(Natm,label,r,rold,v,f,g,dcell,Rcut)
   implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
   character(len=100), dimension(:), allocatable :: label
   real(dp), dimension(:,:), allocatable :: r, rold          ! nth and nth-1 position
   real(dp), dimension(:,:), allocatable :: v                ! velocity
   real(dp), dimension(:,:), allocatable :: f, g             ! These are the forces
   real(dp) :: dcell
   real(dp) :: Rcut
   real(dp) :: ratio
   integer :: Natm
  
   allocate(label(Natm))
   allocate(r(3,Natm))
   allocate(rold(3,Natm))
   allocate(v(3,Natm))
   allocate(f(3,Natm))
   allocate(g(3,Natm))
 
   ratio=Rcut/dcell
   if (ratio>0.5) then
      write(*,*) 'Rcut/box is too large', ratio
      write(*,*) 'Reduce Rcut or increase box size'
      stop 'Error in allocate_arrays'
   end if
  
  end subroutine allocate_arrays

  subroutine rattle_a(r,rold,v,Natm,vec1,vec2,No,dcell)
    ! RATTLE:
    !    
    !-------------------------------------------------------------
    implicit none
    
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp),dimension(3,Natm) :: r,rold
    real(dp),dimension(3,Natm) :: v
    real(dp),dimension(3) :: dr,dr_old,dr_adj
!   integer,dimension(2,No)::            connect
    integer,dimension(No) :: vec1,vec2
    real(dp) :: dr2
    real(dp) :: diffsq,g
    real(dp) :: dot
    real(dp) :: tol,tol2,dot_tol
    real(dp) :: bond                 ! It is  Req
    real(dp) :: dt
    real(dp) :: dcell
    integer ::  Natm,iatm,jatm,cstep
    integer ::  No                   ! No. of bonds found
    integer ::  i    
   
    bond=0.96 
    dt=0.002
    tol=1.0e-9
    tol2=2.0*tol 
    dot_tol=1.0e-9 
 
    do cstep=1,500 ! Until done 
    
       !do iatm=1,Natm-1,2  ! Skip every two atoms
       !   jatm=iatm+1
     
       do i=1,No
          iatm=vec1(i)
          jatm=vec2(i)

             dr=r(:,iatm)-r(:,jatm)
             dr=dr-dcell*ANINT(dr/dcell)   ! PBC 
             dr2=sum(dr**2)

             diffsq=bond**2-sum(dr**2)  ! amount of displacement 
   
             if (abs(diffsq) > tol2*bond**2) then  ! Test whther convergence is reached

                dr_old=rold(:,iatm)-rold(:,jatm) 
!                dr_old=dr_old-dcell*ANINT(dr_old/dcell)   ! PBC 
                dot=dot_product(dr_old,dr)

!                if ( dot < dot_tol*bond**2 ) then 
!                   write(*,*) "Constraint failure", dot, dot_tol, bond**2
!                   STOP
!                end if 
  
                   g=diffsq/(4.0*dot)
                   dr_adj=dr_old*g
                   r(:,iatm)=r(:,iatm)+dr_adj
                   r(:,jatm)=r(:,jatm)-dr_adj
                   v(:,iatm)=v(:,iatm)+dr_adj/dt
                   v(:,jatm)=v(:,jatm)-dr_adj/dt
           
             end if
       end do

       if (cstep>500) then
          write(*,*) "Erorr RATTLE_A: too many iteration steps"
       end if
 
   end do 

  end subroutine
 
!  subroutine rattle_b(r,v,Natm)
!    ! RATTLE:
!    !    
!    !-------------------------------------------------------------
!    implicit none
!    real(dp),dimension(3,Natm)::        r
!    real(dp),dimension(3,Natm)::        v
!    real(dp),dimension(3)::             dr,dv,dv_adj
!    real(dp)::                          g
!    real(dp)::                          dot
!    real(dp)::                          tol,tol2,dot_tol
!    real(dp)::                          bond
!    real(dp)::                          dt
!    integer::                                   Natm,iatm,jatm,cstep
!    
!    dt=0.002
!    bond=1.208 !1.50
!    tol=1.0e-9
!    tol2=2.0*tol 
!    dot_tol=1.0e-9 
!   
! 
!    do cstep=1,500 !50 ! Check each constraint  
!
!       do iatm=1,Natm-1
!          jatm=iatm+1
!       
!             dv=v(:,iatm)-v(:,jatm)
!             dr=r(:,iatm)-r(:,jatm)
!
!             dot=dot_product(dr,dv)
!             g=-dot/(2.0*bond**2)
!             if (abs(dot) > tol) then      ! test whether conistraint is satisfied 
!                   dv_adj=dr*g
!                   v(:,iatm)=v(:,iatm)+dv
!                   v(:,jatm)=v(:,jatm)-dv
!             end if
!  
!       end do
!   
!       if (cstep>500) then
!          write(*,*) "Erorr RATTLE_B: too many iteration steps"
!       end if
!
!   end do 
!
!  end subroutine
 
  subroutine forces_LJ(r,f,Natm,dcell,epotLJ,Rcut,pbc,sigma,eps)
  ! Calculate Intermolecular Force and Potential energy
  ! from Lennard-Jones Potential   
  !-------------------------------------------------------------
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp),dimension(3,Natm) :: r
  real(dp),dimension(3,Natm) :: f
  real(dp),dimension(3) :: dr,fij
  real(dp) :: r_sq,r2,fpr,fr2,fr6,flj
  real(dp) :: Rcut,Rcut2
  real(dp) :: h,dcell
  real(dp) :: sigma, eps, sigma2
  real(dp) :: epotLJ
  integer :: Natm,iatm,jatm
  integer :: pbc
  
  ! Assign LJ parameters
  sigma2=sigma**2
  Rcut2=Rcut**2
  ! Initialize 
  do iatm=1,Natm
     f(:,iatm)=0.0d0 
  end do
  epotLJ=0.d0
  
  do iatm=1,Natm-2 ! Natm-2 if not considering direct neighbouring atoms
     do jatm=iatm+2,Natm
  
        dr=r(:,iatm)-r(:,jatm)
 
        ! minimum image criterion
        if (pbc==1) then 
            dr=dr-dcell*ANINT(dr/dcell)
        end if

        r2=sum(dr**2)
        r_sq=sqrt(r2)
        
        ! Compute LJ Forces
        if (r_sq<Rcut) then

           fr2=sigma2/r2
           fr6=fr2**3
           flj=48.d0*eps*fr6*(fr6-0.5d0)/r2
           !flj=48.d0*eps*fr6*(fr6-0.5d0)/r2
           !flj=4*eps*(12*fr6**2-fr6)/r_sq

           fij=flj*dr ! force factor
  
           f(:,iatm)=f(:,iatm)+fij
           f(:,jatm)=f(:,jatm)-fij

           epotLJ=epotLJ+4.d0*eps*fr6*(fr6-1.d0) 

        end if

     end do
  end do
  
  end subroutine
 
  subroutine forces_spring(r,g,Natm,dcell,vec1,vec2,No,epot_pair,pbc,Req,keq)
  ! Calculate Intramolecular Force and Potential energy 
  ! from Harmonic Forces   
  !-------------------------------------------------------------
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp),dimension(3,Natm) :: r
  real(dp),dimension(3,Natm) :: g
  real(dp),dimension(3) :: dr,gij
  real(dp) :: r_sq,r2,g_pair
  real(dp) :: dbond
  real(dp), intent(in) :: Req
  real(dp) :: Req2
  real(dp), intent(in) :: keq 
  real(dp) :: k3, k4
  real(dp), intent(in) :: dcell
  real(dp) :: epot_pair
  integer,dimension(No) :: vec1,vec2
  integer :: No                   ! No. of bonds found
  integer :: Natm,iatm,jatm
  integer :: i
  integer, intent(in) :: pbc
  
  k3=keq !1.0d+02 ! spring constant cubic
  k4=k3 !1.0d+02 ! spring constant quadratic
  Req2=Req**2
  ! Initialize 
  do iatm=1,Natm
     g(:,iatm)=0.0d0 
  end do
  epot_pair=0.d0
  
       do i=1,No           ! Run over only bound atoms
          iatm=vec1(i)
          jatm=vec2(i)

          dr=r(:,iatm)-r(:,jatm)
 
          ! minimum image criterion
          if (pbc==1) then 
             dr=dr-dcell*ANINT(dr/dcell)
          end if 

          r2=sum(dr**2)
          r_sq=sqrt(r2)
          dbond=(r_sq-Req)        

         ! Compute Harmonic or Quadratic Forces 
          g_pair=-keq*dbond                                                     ! Harmonic form
!          g_pair=-0.5*dbond*(2*k+3*k3*(dbond)**2)                        ! Cubic form
!          g_pair=-0.5*(dbond)*(2*k+3*k3*(dbond)+4*k4*(dbond)**2)    ! Quadratic form
     
          gij=g_pair*dr ! force factor
   
          g(:,iatm)=g(:,iatm)+gij
          g(:,jatm)=g(:,jatm)-gij
     
          epot_pair=epot_pair+0.5*keq*(r_sq-Req)**2                                        ! Harmonic Potential  
!          epot_pair=epot_pair+0.5*k*(dbond)**2+k3*(dbond)**3                      ! Cubic Potential
!          epot_pair=epot_pair+0.5*(k*(dbond**2)+k3*(dbond**3)+k4*(dbond**4))    ! Quadratic Potential 
  
      end do
 ! end do
  
  end subroutine

 
  subroutine set_pbc(r,Natm,dcell)
  ! Adjusts the atomic coordinates for periodic boundary conditions
  ! according to the minimum image criterion
  !----------------------------------------------------------------------------
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp),dimension(3,Natm) :: r
  real(dp) :: dcell
  integer :: iatm,Natm
  
  do iatm=1,Natm
     r(:,iatm)=r(:,iatm)-dcell*ANINT(r(:,iatm)/dcell)
  end do
  
  end subroutine
  
  subroutine thermo_quantities(dcell,Natm,ekin,epot,etot,temp,pressure)
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp) :: etot,etot2,ekin,epot,rms_ene
  real(dp) :: temp,pressure
  real(dp) :: dcell,volume
  integer :: dgf,Natm,kB
  
  !Assign parameters
  dgf=3*Natm-3
  kB=1 
  temp=(2.d0/3.d0)*((kB*ekin)/dgf) !temp=temp+2*ekin/(3*dgf*kb)
  volume=dcell**3
  pressure=0.d0
  etot=(ekin+epot)
  etot2=etot**2
  
  end subroutine 


!  subroutine initialize(x,y,z,vx,vy,vz,Natm)
!  implicit none
!  real(dp),dimension(Natm)::x,y,z
!  real(dp),dimension(Natm)::vx,vy,vz
!  real(dp)::sumv,sumv2,fs
!  real(dp)::dt,r,temp
!  real(dp)::rn
!  integer::iatm,Natm
!  integer::i,nseed,clock,N=100
!  integer,allocatable,dimension(:)::seed
!
!
!  call random_seed(size=nseed) ! get the size of the seed array
!  allocate(seed(nseed))  ! allocate it into an array
!  call system_clock(count=clock)   ! get the current time
!  seed=clock ! copy the current time to the seed so that it'll always be different
!  call random_seed(put=seed)  ! use this to set the seed)
!  call random_number(r)
!
!   
!  dt=0.005
!  temp=1.d0
!  sumv=0.d0
!  sumv2=0.d0
!  do iatm=1,Natm
!     vx(iatm)=r-0.5
!     vy(iatm)=r-0.5
!     vz(iatm)=r-0.5
!     sumv=sumv+vx(iatm)+vy(iatm)+vz(iatm)
!     sumv2=sumv2+vx(iatm)**2+vy(iatm)**2+vz(iatm)**2
!  end do
!
!  sumv=sumv/Natm
!  sumv2=sumv2/Natm
!  fs=sqrt(3*temp/sumv2)
!  do iatm=1,Natm
!     vx(iatm)=(vx(iatm)-sumv)*fs
!     vy(iatm)=(vy(iatm)-sumv)*fs
!     vz(iatm)=(vz(iatm)-sumv)*fs
!     x(iatm)=x(iatm)-vx(iatm)*dt     
!     y(iatm)=y(iatm)-vy(iatm)*dt     
!     z(iatm)=z(iatm)-vz(iatm)*dt     
!  end do
!
!  end subroutine

end module md_module
