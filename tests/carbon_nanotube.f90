! ------------------------------------------------------ !
! Quantulaba - carbon_nanotube.f90 - Krzysztof Kolasinski 2015
!
! Example of hamiltonian creation for simple carbon
! nanotube.
! ------------------------------------------------------ !

program cnt

use modunits     ! unit conversion tools
use modscatter   ! eigen problem and scattering problem
use modsys       ! eigen problem & atom container
use modlead      ! bandgap structure
use modshape
implicit none
character(*),parameter     :: output_folder = "carbon_nanotube_output/"
type(qsys)                 :: tmpsystem
type(qscatter)             :: qt

type(qshape)               :: lead_shape
doubleprecision            :: posA(3),posB(3),deltaR(3),Emin,Emax
integer                    :: flagA,flagB
integer :: i,j,k
doubleprecision            :: lead_translation_vec(3),Ef
character(300) :: line

! Initalize system
call qt%init_system()
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
        call qt%qatom%init(posA+deltaR,flag=flagA)
        if(i==1) call qt%qatom%init(posA+deltaR,flag=flagA,flag_aux0=1)
        call qt%qsystem%add_atom(qt%qatom)
    enddo
enddo


!  ----------------------------------------------------------
!  2. Set nearest neightbour search radius and make Hmiltonian
!  ----------------------------------------------------------
qt%qnnbparam%distance   = 2.0
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=coupling)


!  ----------------------------------------------------------
!  3. Set up unit cell
!  ----------------------------------------------------------

lead_translation_vec = (/  0.0D0 , 0.0D0 , 4.26D0 /)
posA = (/0.0,0.0,-0.5/)
posB = lead_translation_vec ! direction and range
call lead_shape%init_range_3d(posA,posB)
call qt%add_lead(lead_shape,lead_translation_vec)

!  ----------------------------------------------------------
!  4. Plot band structure
!  ----------------------------------------------------------
Emin = -3.1
Emax =  3.1
call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/60.0,Emin,Emax)


posA = (/0.0,0.0,42.0/)
posB = -lead_translation_vec ! direction and range
call lead_shape%init_range_3d(posA,posB)
call qt%add_lead(lead_shape,-lead_translation_vec)
call qt%save_system(output_folder//"system.xml")

Ef = 0.2
QSYS_DEBUG_LEVEL = 1
call qt%calculate_modes(Ef)
call qt%solve(1,Ef)
! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)


print*,"Performing energy scan..."
open(unit=111,file=output_folder//"T.dat")
!QSYS_DEBUG_LEVEL = 1 ! show more info
do Ef = -3.0 + 0.01 , 3.025 , 0.051
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_simple=coupling)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)

    print*,"Energy:",Ef
    write(111,"(2f20.6)"),Ef,sum(qt%Tn(:))
enddo
close(111)



call qt%destroy_system()
call tmpsystem%destroy()

print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")
print*,"Use Viewer program to see the structure and created lead."
 contains

! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA) with spin s1 has hoping
! to atom B with spin s2, and what is the value of the coupling.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
logical function coupling(atomA,atomB,coupling_val)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_val ! you must overwrite this variable

    ! default return value
    coupling = .false.
    coupling_val   = 0.0
    if( atomA%flag /= atomB%flag ) then
        coupling = .true.
        coupling_val = 1.0
!    else ! Add small perturbation to the lattice to remove periodic boundary degeneracy
!         ! this should not change !significantly! the resutls
!        coupling = .true.
!        coupling_val = atomB%flag_aux0*0.0000001
    endif
end function coupling
end program cnt
