program main
!subroutine  build_and_read(var2,var3,No)

 use iso_fortran_env

  implicit none

   integer :: unitin=19
   integer :: unitout=22
   integer :: Natm, Nbond                 ! No atoms,and bonds
   integer :: Nsubst, Nfeat, Nset         ! No substrates, features, and something else   
   integer :: ios
   integer :: nline
   integer :: var1, var4 
   integer, dimension(:), allocatable :: var2, var3
   character(len=100) :: label
   character(len=100) :: match
   character(len=100) :: filename

   call get_command_argument(1,filename)
   call execute_command_line("babel -ixyz "//trim(filename)//" "//"-omol2 coord.mol2")
   open(unitin, file="coord.mol2", status="old")
   open(unitout, file="connection.txt", status="replace")

   match="@<TRIPOS>MOLECULE"
   do 
     read(unitin,*,iostat=ios) label
     if (ios == iostat_end) exit
     if (label.eq.match) then
        read(unitin,*,iostat=ios) label
        read(unitin,*,iostat=ios) Natm, Nbond, Nsubst, Nfeat, Nset       
        write(unitout,*) Nbond     
        allocate(var2(Nbond),var3(Nbond))
     else if (label.eq.'@<TRIPOS>BOND') then
        nline=1
        do while (nline.le.Nbond.and.ios.eq.0)       
           read(unitin,*,iostat=ios) var1, var2(nline), var3(nline), var4  
           write(unitout,*)  var2(nline), var3(nline)  
           nline=nline+1
        end do
     end if
   end do     

   close(unitout)
   close(unitin)


end program
