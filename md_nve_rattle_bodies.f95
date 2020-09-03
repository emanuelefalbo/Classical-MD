program main

 use md_module

  ! Takes coordinates in and carry out MD using velocit Verlet algorithm 
  ! Reduced units for Lennard-Jones:  mass=1, kB=1, time=2.17*10^-12 s
  ! velocity 1.57*10^2
  
  implicit none

    character(len=100)::                             arg1
    character(len=100),dimension(:),allocatable::    label
    double precision,dimension(:,:),allocatable::    r
    double precision,dimension(:,:),allocatable::    rold
    double precision,dimension(:,:),allocatable::    f
    double precision,dimension(:,:),allocatable::    g
    double precision,dimension(:,:),allocatable::    v
    integer,dimension(:,:),allocatable::             connect
    integer,dimension(:),allocatable::               vec1,vec2
    double precision::                               t,tmax,tmin,mass,alpha,k
    double precision::                               dt
    double precision::                               ekin,epot,etot,etot2,rms
    double precision::                               temp,pressure
    double precision::                               pi=4.0d0*atan(1.d0),dcell
    double precision::                               h1,h2
    double precision::                               epotLJ,epot_pair
    double precision::                               v2
    real ::                                          finish,start,time
    integer::                                        i,j,iatm,Natm,Mstep,nout
    integer::                                        rattle
    integer::                                        ios
    integer::                                        No,nlines
    
    ! Assign variable
    call control(dt,Mstep,nout,dcell,rattle)
    t=0.0
    mass=1.0
    
    ! Read in connections
    
       call execute_command_line("wc -l connection.txt| cut -c -2 > nlines.txt")
       open(50,file="nlines.txt")     
       read(50,*) No 
       close(50)
       allocate(vec1(No),vec2(No))
       open(40,file="connection.txt")    
       do i=1,No
         read(40,*,iostat=ios) vec1(i), vec2(i) !connect(:,i)
       end do
       
    ! Read in coordinates
    
    call getarg(1,arg1)
    open(18,file=TRIM(ADJUSTL(arg1)))    
    read(18,*) Natm
    ! Allocate vectors
    allocate(r(3,Natm),rold(3,Natm))
    allocate(v(3,Natm))
    allocate(f(3,Natm))
    allocate(g(3,Natm))
    allocate(label(Natm))
    ! End allocate
    read(18,*)
    do iatm=1,Natm
       read(18,*,iostat=ios) label(iatm), r(:,iatm) 
    end do
    close(18)
    
    ! MD Routine
    
    open(20,file="trajectory.xyz", status="replace")    
    open(30,file="result.dat", status="replace")    
    
    call cpu_time(start)  ! Time MD engine 

    do i=1,Mstep
      
     ! Initialize force and energies
       call forces_LJ(r,f,Natm,dcell,epotLJ)
       call forces_spring(r,g,Natm,dcell,vec1,vec2,No,epot_pair)
       epot=epot_pair+epotLJ
       !epot=epotLJ
       ekin=0.d0
      
       ! Verlet loop 
       do iatm=1,Natm
          rold(:,iatm)=r(:,iatm)  ! Store the nth-1 step coordinates
       end do
     
       do iatm=1,Natm
          !v(:,iatm)=v(:,iatm)+0.5d0*f(:,iatm)*dt
          v(:,iatm)=v(:,iatm)+0.5d0*dt*(f(:,iatm) + g(:,iatm)) ! LJ + Stretch
          r(:,iatm)=r(:,iatm)+v(:,iatm)*dt
       end do  
      
       call pbc(r,Natm,dcell)
!       call rattle_a(r,rold,v,Natm,vec1,vec2,No,dcell)
       call forces_LJ(r,f,Natm,dcell,epotLJ)
       call forces_spring(r,g,Natm,dcell,vec1,vec2,No,epot_pair)
          
       do iatm=1,Natm
!         v(:,iatm)=v(:,iatm)+0.5d0*f(:,iatm)*dt
         v(:,iatm)=v(:,iatm)+0.5d0*dt*(f(:,iatm) + g(:,iatm)) ! LJ + Stretch
          !call rattle_b(r,v,Natm)
     
          ! Compute energies
          v2=sum(v(:,iatm)**2) 
          ekin=ekin+(0.5d0*v2)
       end do
    
       ! Compute properties
       call thermo_quantities(dcell,Natm,ekin,epot,etot,temp,pressure)
       t=t+dt
     
       ! Print output every nout steps
       if (DBLE(i)/nout.eq.DBLE(FLOOR(DBLE(i)/nout))) then
          write(20,*) Natm
          write(20,*) "Cartesian"
    
          do iatm=1,Natm
             write(20,*) TRIM(ADJUSTL(label(iatm))),r(:,(iatm)) !ekin,epot,etot
          end do
    
          write(30,*) i,t,temp,ekin,epot_pair,epot,etot
       end if
    end do
    
    call cpu_time(finish)
    time=finish-start
    write(*,*) "Time = ", time, " seconds"      
    close(20)
    close(30)
    
end program main
