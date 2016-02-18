! ------------------------------------------------------ !
! Quantulaba - simple_graphene.f90 - Krzysztof Kolasinski 2015
!
! ------------------------------------------------------ !

program graphene2
use modscatter
use modsys
use modshape
use modunits
implicit none
type(qscatter) :: qt
type(qshape) :: l_shape,lower_rect,middle_rect,upper_rect

character(*),parameter    :: output_folder = "plots/"
doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI
doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)
doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
doubleprecision           :: atom_pos(3),lead_translation_vec(3),range_base(3),range_dir(3),Ef
integer,parameter         :: atom_A = 1 , atom_B = 2
integer           :: atom,i,j
integer,parameter :: nx = 50  , ny = 10
double precision  :: xmin,xmax,ymin,ymax,l_start,l_width,g_delta
doubleprecision   :: Bz
logical :: test

xmin = 0.0D0
xmax = 200.5D0
ymin = 0.0D0
ymax = 100.2D0
l_start = 5.0D0
l_width = 30.2D0
g_delta = 0.4
call qt%init_system()

QSYS_DEBUG_LEVEL = 2

! Wykorzystam obiekt shape aby nadac ukladowi ksztalt litery H
! ale odwrocony
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
! Graphene test
! 1. Generate atoms positions:
! --------------------------------------------------------------------------
do i = -500 , 500
do j = -100 , 500
    do atom = atom_A , atom_B
    ! set atom position in space
        atom_pos(1:2)   = atoms_armchair(:,atom) + &
                   (i-1) * vecs_armchair(:,1) +  &
                   (j-1) * vecs_armchair(:,2)

        call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)

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
qt%qnnbparam%distance   = 0.6
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

! remove single bonds if necessary, loop through all amtoms and
! check if the number of bonds is 2. (1st - for on site bond, 2nd for neightbour)
do atom = 1 ,qt%qsystem%no_atoms
    if(qt%qsystem%atoms(atom)%no_bonds == 2)   qt%qsystem%atoms(atom)%bActive  = .false.
enddo
! 3. Calculate coupling again.
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

! --------------------------------------------------------------------------
! 3. Add translational lead to the system. We use range_3d object to define
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



call qt%save_system(output_folder//"system.xml")

Bz = -0.09
Ef = 0.8
call qt%qsystem%update_lattice(c_simple=connect)
call qt%leads(1)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)

call qt%calculate_modes(Ef)
call qt%solve(1,Ef)

! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)

! print transmission through each lead using scattering matrix
print*,"T(1,1)  =",sum(abs(qt%smatrix(1,1)%Tnm)**2)
print*,"T(1,2)  =",sum(abs(qt%smatrix(1,2)%Tnm)**2)
print*,"T(1,3)  =",sum(abs(qt%smatrix(1,3)%Tnm)**2)
print*,"T(1,4)  =",sum(abs(qt%smatrix(1,4)%Tnm)**2)
! print total tranmission and reflection using auxiliary matrices
print*,"T     =",sum(qt%Tn(:))
print*,"R     =",sum(qt%Rn(:))


!
!print*,"Performing energy scan..."
!open(unit=111,file=output_folder//"T.dat")
!!QSYS_DEBUG_LEVEL = 1 ! show more info
!do Ef = -3.0 , 3.025 , 0.025
!    ! Update hamiltonian elemenents value
!    call qt%qsystem%update_lattice(c_simple=connect)
!    call qt%calculate_modes(Ef)
!    call qt%solve(1,Ef)
!
!    print*,"Energy:",Ef
!    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
!enddo
!close(111)

print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
!print*,"Plotting Transmission..."
!call system("cd "//output_folder//"; ./plot_T.py")
!print*,"Use Viewer program to see the structure and created leads."

call qt%destroy_system()
contains

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
end program graphene2
