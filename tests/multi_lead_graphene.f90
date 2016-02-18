! ------------------------------------------------------ !
! Bubel - multi_lead_graphene.f90 - Krzysztof Kolasinski 2016
! Example on adding multiple leads into system and s-matrix
! access. Demo is done on graphene flake.
! ------------------------------------------------------ !

program multi_lead
use modscatter
use modsys
use modshape
use modunits
implicit none
type(qscatter) :: qt
type(qshape) :: l_shape,lower_rect,middle_rect,upper_rect

character(*),parameter    :: output_folder = "multi_lead_graphene_output/"
character(*),parameter    :: names(4) = (/"1","2","3","4"/)
doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI
doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)
doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
doubleprecision           :: atom_pos(3),lead_translation_vec(3),range_base(3),range_dir(3),Ef
integer,parameter         :: atom_A = 1 , atom_B = 2
integer           :: atom,i,j,lead
double precision  :: xmin,xmax,ymin,ymax,l_start,l_width,g_delta
doubleprecision   :: Bz
logical :: test

print*,"-----------------------------------------------"
print*,"        Starting::multi_lead_graphene.f90"
print*,"-----------------------------------------------"

! Define physical dimensions of the system (coordinates)
xmin    =   0.0D0
xmax    = 100.5D0
ymin    =   0.0D0
ymax    = 100.2D0
l_start =   5.0D0   ! lenght of the input leads
l_width =  30.2D0   ! width of the leads, here we assume constant width
g_delta =   0.4     ! parameter which corrects position of lead shapes in order to fit them into lattice
call qt%init_system()
QSYS_DEBUG_LEVEL = 1

! We will use shape type to define device shape, we want to simulate such system:
!
!             (UL) = ====== = (UR)
!                    ======
!                    ======
!             (BL) = ====== = (BR)
!
! where U/B- stays for upper/bottom and L/R - left/right

call lower_rect%init_rect(SHAPE_RECTANGLE_XY,&
                        xmin-l_start,&
                        xmax+l_start,&
                        ymin,ymin+l_width)

call middle_rect%init_rect(SHAPE_RECTANGLE_XY,&
                        xmin,&
                        xmax,&
                        ymin,ymax)

call upper_rect%init_rect(SHAPE_RECTANGLE_XY,&
                        xmin-l_start,&
                        xmax+l_start,&
                        ymax-l_width,ymax)
! --------------------------------------------------------------------------
! 1. Generate atoms positions:
! --------------------------------------------------------------------------
do i = -500 , 500 ! some large values which will fill all the space necessary to make required flake
do j = -100 , 500
    do atom = atom_A , atom_B
    ! set atom position in space
        atom_pos(1:2)   = atoms_armchair(:,atom) + &
                   (i-1) * vecs_armchair(:,1) +  &
                   (j-1) * vecs_armchair(:,2)

        call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)

        ! check if atom is in the device shape
        test = lower_rect%is_inside(qt%qatom%atom_pos)
        test = test .or. middle_rect%is_inside(qt%qatom%atom_pos)
        test = test .or. upper_rect%is_inside(qt%qatom%atom_pos)
        if(test)then
            call qt%qsystem%add_atom(qt%qatom)
        endif
    enddo ! end of atom loop
enddo
enddo
! --------------------------------------------------------------------------
! 2. Initiate couplings between atoms: We use c_simple connection type
!    i.e. without including spin degree of freedom per atom
! --------------------------------------------------------------------------
qt%qnnbparam%distance   = 0.6 ! prefiltering - it will reject all the atoms which distance is bigger than 0.6
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

! Remove single bonds if necessary, loop through all atoms and
! check if the number of bonds is 2. (1st - for on site bond, 2nd for neightbour)
do atom = 1 ,qt%qsystem%no_atoms
    if(qt%qsystem%atoms(atom)%no_bonds == 2)   qt%qsystem%atoms(atom)%bActive  = .false.
enddo
! Calculate coupling again with make_lattice
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

! --------------------------------------------------------------------------
! 3. Add translational lead to the system. We use init_rect object to define
!    area where the lead is located in the system.
! --------------------------------------------------------------------------
! (1) bottom left lead
call l_shape%init_rect(SHAPE_RECTANGLE_XY,&
            xmin-l_start,xmin-l_start+1+g_delta,&
            ymin,ymin+l_width)
call qt%add_lead(l_shape,(/1.0D0,0.0D0,0.0D0/))

! (2) top left lead
call l_shape%init_rect(SHAPE_RECTANGLE_XY,&
            xmin-l_start,xmin-l_start+1+g_delta,&
            ymax-l_width,ymax)
call qt%add_lead(l_shape,(/1.0D0,0.0D0,0.0D0/))


! (3) bottom right lead
call l_shape%init_rect(SHAPE_RECTANGLE_XY,&
            xmax+l_start-1+g_delta,xmax+l_start+g_delta,&
            ymin,ymin+l_width)
call qt%add_lead(l_shape,(/-1.0D0,0.0D0,0.0D0/))

! (4) upper right lead
call l_shape%init_rect(SHAPE_RECTANGLE_XY,&
            xmax+l_start-1+g_delta,xmax+l_start+g_delta,&
            ymax-l_width,ymax)
call qt%add_lead(l_shape,(/-1.0D0,0.0D0,0.0D0/))


! debug system
call qt%save_system(output_folder//"system.xml")

! Set physical parameters
Bz = -0.03 ! magnetic field
Ef =   0.8 ! Fermi energy
! Update lattice, since we changed magnetic field, you need to recalculate Hamiltonian
call qt%qsystem%update_lattice(c_simple=connect)
call qt%leads(1)%update_lead(qt%qsystem%atoms) ! you can also use: call qt%calculate_modes(Ef), but this is faster
call qt%leads(1)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)

! Calculate modes in all the leads, this clean s-matrix also
call qt%calculate_modes(Ef)
! Solve scattering problem for electron incoming from all leads:
do lead = 1 , 4
    call qt%solve(lead,Ef)
    do i = 1 , size(qt%qauxvec)
        qt%qauxvec(i) = sum(qt%densities(:,i))
    enddo
    call qt%qsystem%save_data(output_folder//"densities"//names(lead)//".xml",array2d=qt%densities,array1d=qt%qauxvec)
enddo

! Scattering matrix can be accessed. In the loop above we have
! calculated full S-matrix -- qt%smatrix object. S-matrix has
! following form:
!   qt%smatrix(l_in,l_out)%Tnm(m_in,m_out)
! where:
!   l_in  - defines source lead (from)
!   l_out - defines output lead (to)
!   m_in  - defines incoming  mode in the lead l_in
!   m_out - defines outcoming mode in the lead l_out
! all values for given l_in are calculated in: call qt%solve(l_in,Ef)
! S-matrix is cleaned after: qt%calculate_modes(Ef)
! Note, that if you are only interested in transmission from the 1st lead
! you still can use qt%smatrix(1,l_out)%Tnm, since those matrix elements will
! be filled after: call qt%solve(1,Ef)
! Let us pring all transmisions probablities for each lead
do i = 1 , 4
do j = 1 , 4
    print"(A,i3,A,i3,A,f10.6)","T(",i,",",j,")=",sum(abs(qt%smatrix(i,j)%Tnm)**2)
enddo
enddo
! Or you can print total tranmission and reflection using auxiliary matrices
! but this remember only the last call of qt%solve(lead,Ef) and it contains
! summed tranmission probabilies for all leads
print*,"T     =",sum(qt%Tn(:))
print*,"R     =",sum(qt%Rn(:))



print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; python plot_bands.py")
call qt%destroy_system()
contains
! --------------------------------------------------------
! Calculate hoping between atoms, here we use Peltier phase
! to simulate magnetic field with gauge: A = (-Bz * y,0,0)
! i.e. B = (0,0,Bz)
! --------------------------------------------------------
logical function connect(atomA,atomB,coupling_val)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16 :: coupling_val
    doubleprecision :: xj,xi,yj,yi,phi_ij
    logical :: test
    connect = .false.

    xi = atomA%atom_pos(1)
    yi = atomA%atom_pos(2)
    xj = atomB%atom_pos(1)
    yj = atomB%atom_pos(2)
    phi_ij = 0.5*Bz*(yj+yi)*(xj-xi)

    ! In clean graphene atom from sublattice A always couples to atom from lattice B
    test = .not.(atomA%flag == atomB%flag)
    if(test) then
        connect      = .true.
        coupling_val = 1.0D0*exp(II*phi_ij)
    endif
end function connect
end program multi_lead
