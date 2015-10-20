! ------------------------------------------------------ !
! Quantulaba - carbon_nanotube.f90 - Krzysztof Kolasinski 2015
!
! Example of hamiltonian creation for simple carbon
! nanotube.
! ------------------------------------------------------ !

program transporter

 use modunits ! unit conversion tools
 use modsys   ! eigen values
 use modlead  ! bandgap structure
 use modshape
 implicit none
 character(*),parameter :: output_folder = "carbon_nanotube_output/"
 type(qsys)                 :: tmpsystem,qsystem
 type(qshape)               :: lead_shape
 type(qlead)                :: lead
 doubleprecision            :: posA(3),posB(3),deltaR(3),Emin,Emax
 integer                    :: flagA,flagB
 integer :: i,j,k
 doubleprecision            :: lead_translation_vec(3)
 character(300) :: line

 ! Initalize system
 call qsystem%init()
 call tmpsystem%init()

!  ----------------------------------------------------------
!  1. Create mesh - read unit cell from file
!  ----------------------------------------------------------
open(unit=3,file=output_folder//"nt-17-0-1.xyz")
read(3,*) line
read(3,*) line
! ----------------------------------
! Reading unit cell
! ----------------------------------
do i=1,17*2
    read(3,*) line, (posA(j), j=1, 3)
    read(3,*) line, (posB(j), j=1, 3)
    call tmpsystem%qatom%init(posA,flag=1)
    call tmpsystem%add_atom(tmpsystem%qatom)
    call tmpsystem%qatom%init(posB,flag=2)
    call tmpsystem%add_atom(tmpsystem%qatom)
enddo
close(3)

! ----------------------------------
! Translating unit cell X-times
! ----------------------------------
do k = 1 , 10
    deltaR = (/0.0,0.0,(k-1)*4.26/)
    do i=1,tmpsystem%no_atoms
        posA  = tmpsystem%atoms(i)%atom_pos
        flagA = tmpsystem%atoms(i)%flag
        call qsystem%qatom%init(posA+deltaR,flag=flagA)
        call qsystem%add_atom(qsystem%qatom)
    enddo
enddo


!  ----------------------------------------------------------
!  2. Set nearest neightbour search radius and make Hmiltonian
!  ----------------------------------------------------------
qsystem%qnnbparam%distance   = 2.0
qsystem%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qsystem%make_lattice(qsystem%qnnbparam,c_simple=connect_simple)
call qsystem%save_lattice(output_folder//"lattice.xml")

!  ----------------------------------------------------------
!  3. Set up unit cell
!  ----------------------------------------------------------

lead_translation_vec = (/  0.0D0 , 0.0D0 , 4.26D0 /)
posA = (/0.0,0.0,-0.5/)
posB = lead_translation_vec ! direction and range
call lead_shape%init_range_3d(posA,posB)
call lead%init_lead(lead_shape,lead_translation_vec,qsystem%atoms)
call lead%save_lead(output_folder//"lead.xml")

!  ----------------------------------------------------------
!  4. Calculate band structure
!  ----------------------------------------------------------
Emin = -3.0
Emax =  3.0
call lead%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/60.0,Emin,Emax)
call lead%destroy()

call qsystem%destroy()
call tmpsystem%destroy()

print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Use Viewer program to see the structure and crated lead."
 contains

! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA) with spin s1 has hoping
! to atom B with spin s2, and what is the value of the coupling.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
logical function connect_simple(atomA,atomB,coupling_val)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_val ! you must overwrite this variable

    ! default return value
    connect_simple = .false.
    coupling_val   = 0.0
    if( atomA%flag /= atomB%flag ) then
        connect_simple = .true.
        coupling_val = 1.0
    endif
end function connect_simple
end program transporter
