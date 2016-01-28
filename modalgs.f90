! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2016
! ---------------------------------------------------------------------------------------

module modalgs
use modutils
use modcommons
use lapack95
use blas95
use f95_precision
!use ifport

implicit none
public


contains


! ------------------------------------------------------------------------
!                             Helper functions
! ------------------------------------------------------------------------


! Conversion from row-col-value sparse matrix storage to
! HB compressed format.
subroutine convert_to_HB(no_vals,rows_cols,matA,out_rows)
      integer,intent(in)                  :: no_vals
      integer   ,intent(inout),dimension(:,:) :: rows_cols
      complex*16,intent(inout),dimension(:) :: matA
      integer,intent(inout),dimension(:)   :: out_rows
      integer :: iterator, irow ,  from , to
      integer :: i, n , j

      n        = no_vals
      iterator = 0
      irow     = 0
      do i = 1 , n
          if( rows_cols(i,1) /= irow ) then
            iterator = iterator + 1
            out_rows(iterator) = i
            irow = rows_cols(i,1)
          endif
      enddo
      out_rows(iterator+1) = n + 1


!DEC$ IF DEFINED  (USE_UMF_PACK)
  irow = size(out_rows)-1
    ! sortowanie  kolumn
  do i = 1 , irow
  from = out_rows(i)
  to   = out_rows(i+1)-1
      call sort_col_vals(rows_cols(from:to,2),matA(from:to))
  enddo

  ! przesuwanie indeksow do zera
  out_rows       = out_rows -1
  rows_cols(:,2) = rows_cols(:,2) -1
!DEC$ ENDIF

!DEC$ IF DEFINED  (USE_PARDISO)
  irow = size(out_rows)-1
  do i = 1 , irow
      from = out_rows(i)
      to   = out_rows(i+1)-1
      call sort_col_vals(rows_cols(from:to,2),matA(from:to))
  enddo
!DEC$ ENDIF
end subroutine convert_to_HB


! Conversion from row-col-value sparse matrix storage to
! HB compressed format. For double precision
subroutine dalg_convert2HB(no_vals,rows_cols,matA,out_rows)
      integer,intent(in)                  :: no_vals
      integer   ,intent(inout),dimension(:,:) :: rows_cols
      doubleprecision,intent(inout),dimension(:) :: matA
      integer,intent(inout),dimension(:)   :: out_rows
      integer :: iterator, irow ,  from , to
      integer :: i, n , j

      n        = no_vals
      iterator = 0
      irow     = 0
      do i = 1 , n
          if( rows_cols(i,1) /= irow ) then
            iterator = iterator + 1
            out_rows(iterator) = i
            irow = rows_cols(i,1)
          endif
      enddo
      out_rows(iterator+1) = n + 1

!DEC$ IF DEFINED  (USE_UMF_PACK)
      irow = size(out_rows)-1
        ! sortowanie  kolumn
      do i = 1 , irow
      from = out_rows(i)
      to   = out_rows(i+1)-1
          call dalg_SortColValues(rows_cols(from:to,2),matA(from:to))
      enddo

      ! przesuwanie indeksow do zera
      out_rows       = out_rows -1
      rows_cols(:,2) = rows_cols(:,2) -1
!DEC$ ENDIF

!DEC$ IF DEFINED  (USE_PARDISO)
      irow = size(out_rows)-1
      do i = 1 , irow
          from = out_rows(i)
          to   = out_rows(i+1)-1
          call dalg_SortColValues(rows_cols(from:to,2),matA(from:to))
      enddo
!DEC$ ENDIF
end subroutine dalg_convert2HB

subroutine dalg_SortColValues(cols,vals)
        integer,intent(inout),dimension(:)    :: cols
        doubleprecision,intent(inout),dimension(:) :: vals
        integer         :: tmp_col
        doubleprecision :: tmp_val
        integer :: i  , n
        logical :: test
        n = size(cols)

        test = .true.

        ! sortowanie bombelkowe
        do while(test)
          test = .false.
          do i = 1 , n-1
            if( cols(i) > cols(i+1)  ) then
            tmp_col   = cols(i)
            cols(i)   = cols(i+1)
            cols(i+1) = tmp_col

            tmp_val   = vals(i)
            vals(i)   = vals(i+1)
            vals(i+1) = tmp_val

            test = .true.
            exit
            endif
          enddo
        enddo
end subroutine dalg_SortColValues

subroutine sort_col_vals(cols,vals)
        integer,intent(inout),dimension(:)    :: cols
        complex*16,intent(inout),dimension(:) :: vals
        integer :: tmp_col
        complex*16 :: tmp_val
        integer :: i  , n
        logical :: test
        n = size(cols)

        test = .true.

        ! sortowanie bombelkowe
        do while(test)
          test = .false.
          do i = 1 , n-1
            if( cols(i) > cols(i+1)  ) then
            tmp_col   = cols(i)
            cols(i)   = cols(i+1)
            cols(i+1) = tmp_col

            tmp_val   = vals(i)
            vals(i)   = vals(i+1)
            vals(i+1) = tmp_val

            test = .true.
            exit
            endif
          enddo
        enddo
end subroutine sort_col_vals




  subroutine solve_SSOLEQ(no_rows,no_vals,colptr,rowind,values,b,iopt,mtype)
        use modutils
        implicit none
        integer,intent(in)                 :: no_rows
        integer,intent(in)                 :: no_vals
        integer,intent(in),dimension(:)    :: colptr,rowind
        complex*16,intent(in),dimension(:) :: values
        complex*16,intent(inout),dimension(:) :: b
        integer :: iopt,mtype
        integer n, nnz, nrhs, ldb

        integer, save    ::  info = 0
        integer*8 , save :: factors = 0

        doubleprecision,save :: total_time
!DEC$ IF DEFINED  (USE_UMF_PACK)
        ! UMFPACK constants
        type(c_ptr),save :: symbolic,numeric
        ! zero-based arrays
        real(8),save :: control(0:UMFPACK_CONTROL-1),umf_info(0:UMFPACK_INFO-1)
        complex*16,allocatable,dimension(:),save :: b_sol

!DEC$ ENDIF

!DEC$ IF DEFINED  (USE_PARDISO)

        INTEGER*8,save  :: pt(64)
        INTEGER,save    :: phase
        INTEGER,save    :: maxfct, mnum , error, msglvl
        INTEGER,save    :: iparm(64)
        complex*16,allocatable,dimension(:),save :: b_sol

        INTEGER    ,save::  idum(1)
        COMPLEX*16 ,save::  ddum(1)

!DEC$ ENDIF

        n    = no_rows
        nnz  = no_vals
        ldb  = n
        nrhs = 1


!DEC$ IF DEFINED  (USE_UMF_PACK)
      selectcase (iopt)
      case (1)
            total_time = get_clock();
            allocate(b_sol(size(b)))

            call umf4cdef (control)
            call umf4csym (n,n, rowind, colptr, values, symbolic, control, umf_info)
            call umf4cnum (rowind, colptr, values, symbolic, numeric, control, umf_info)
            call umf4cfsym (symbolic)
            !total_time =  umf_info(UMFPACK_NUMERIC_TIME)+umf_info(UMFPACK_SYMBOLIC_TIME)
            if (umf_info(UMFPACK_STATUS) .eq. 0) then
                if(TRANS_DEBUG) then
                     write (*,*) 'Factorization succeeded. Mem needed:', umf_info(UMFPACK_PEAK_MEMORY)/8.0/1024/1024 , "[MB]"
                endif
            else
                 write(*,*) 'SYS::UMF ERROR:: INFO from factorization step:', umf_info(UMFPACK_STATUS)
            endif

      case(2)
            b_sol = 0
            call umf4csolr (UMFPACK_Aat, rowind, colptr, values, b_sol, b, numeric, control, umf_info)
            b  = b_sol;

            if (umf_info(UMFPACK_STATUS) .eq. 0) then
                if(TRANS_DEBUG) then
                    write (*,*) 'Solve succeeded. Time needed:',umf_info(UMFPACK_SOLVE_WALLTIME)
                endif
            else
                 write(*,*) 'SYS::UMF ERROR: INFO from solve:', umf_info(UMFPACK_STATUS)
            endif

      case(3)
            print*,"UMFPACK Solved:"
            print*,"Solve time needed:",get_clock()-total_time,"[s]"
            call umf4cfnum (numeric)
            deallocate(b_sol)
      endselect

!DEC$ ELSE IF DEFINED  (USE_PARDISO)


 selectcase (iopt)
      case (1)
      allocate(b_sol(size(b)))
          total_time = get_clock();
          maxfct = 1 ! in many application this is 1
          mnum   = 1 ! same here
          iparm = 0
          iparm(1) = 1 ! no solver default
          iparm(2) = 2 ! fill-in reordering from METIS
          iparm(3) = 1 ! numbers of processors, value of OMP_NUM_THREADS
          iparm(4) = 0 ! 0 - no iterative-direct algorithm, if 1 multirecursive iterative algorithm 61, 31 para me
          iparm(5) = 0 ! no user fill-in reducing permutation
          iparm(6) = 0 ! =0 solution on the first n compoments of x
          iparm(7) = 0 ! not in use
          iparm(8) = 2 ! numbers of iterative refinement steps
          iparm(9) = 0 ! not in use
          iparm(10) = 10 ! perturbe the pivot elements with 1E-13
          iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
          iparm(12) = 0 ! not in use
          iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric).
          iparm(14) = 0 ! Output: number of perturbed pivots
          iparm(15) = 0 ! not in use
          iparm(16) = 0 ! not in use
          iparm(17) = 0 ! not in use
          iparm(18) = -1 ! Output: number of nonzeros in the factor LU
          iparm(19) = -1 ! Output: Mflops for LU factorization
          iparm(20) = 0 ! Output: Numbers of CG Iterations
          iparm(27) = 0 ! perform matrix check
          iparm(32) = 0 ! if 1 use multirecursive iterative algorithm

           error = 0 ! initialize error flag
          msglvl = 0 ! print statistical information
!          mtype     = 13     ! complex unsymmetric matrix
          phase     = 11      ! only reordering and symbolic factorization



          CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)

!          WRITE(*,*) 'Reordering completed ... '

          IF (error .NE. 0) THEN
            WRITE(*,*) 'SYS::PARDISO::The following ERROR was detected during the reordeing step:', error
            STOP 1
          END IF


!          WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

    !C.. Factorization.
          phase     = 22  ! only factorization
          CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)

!          WRITE(*,*) 'Factorization completed ...  '
          IF (error .NE. 0) THEN
             WRITE(*,*) 'SYS::PARDISO::The following ERROR was detected during the factorization step:', error
            STOP 1
          ENDIF
          if(QSYS_DEBUG_LEVEL > 0) then
          WRITE(*,*) 'Peak memory usage   = ',max (IPARM(15), IPARM(16)+IPARM(17))/1024.0,"[MB]"
          endif
      case(2)
          b_sol = 0
    !C.. Back substitution and iterative refinement
          phase     = 33  ! only factorization
          iparm(8)  = 3   ! max numbers of iterative refinement steps
          CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, b, b_sol, error)

          b  = b_sol;
          !WRITE(*,*) 'Solve completed ... '
          IF (error .NE. 0) THEN
             WRITE(*,*) 'SYS::PARDISO::The following ERROR was detected: ', error
          ENDIF

      case(3)
    !C.. Termination and release of memory
            phase     = -1           ! release internal memory
            CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)
            if(QSYS_DEBUG_LEVEL > 0) then
            print*,"SYS::PARDISO::Solve time needed:",get_clock()-total_time,"[s]"
            endif
            deallocate(b_sol)
      endselect


!DEC$ ELSE
      selectcase (iopt)
      case (1)
      total_time = get_clock();
! First, factorize the matrix. The factors are stored in *factors* handle.
      !iopt = 1
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr , rowind , b, ldb,factors, info )
!
      if (info .eq. 0) then
!         write (*,*) 'Factorization succeeded'
      else
         write(*,*) 'SYS::c_fortran_zgssv::SuperLU INFO from factorization step:', info
      endif
      case(2)
! Second, solve the system using the existing factors.
!      iopt = 2
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr,rowind ,  b, ldb,factors, info )
!
      if (info .eq. 0) then
!         write (*,*) 'Solve succeeded'
!         write (*,*) (b(i), i=1, n)
      else
         write(*,*) 'SYS::c_fortran_zgssv::SuperLU ERROR from triangular solver:', info
      endif
      case(3)
! Last, free the storage allocated inside SuperLU
!      iopt = 3
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr,rowind, b, ldb,factors, info )
      if(QSYS_DEBUG_LEVEL > 0) then
      print*,"SYS::SuperLU::Computations time:",get_clock()-total_time,"[s]"
      endif
      endselect
!DEC$ ENDIF

      endsubroutine solve_SSOLEQ


  subroutine dalg_SSOLEQ(no_rows,no_vals,colptr,rowind,values,b,iopt,pardiso_mtype)
        use modutils
        implicit none
        integer,intent(in)                 :: no_rows
        integer,intent(in)                 :: no_vals
        integer,intent(in),dimension(:)    :: colptr,rowind
        doubleprecision,intent(in),dimension(:) :: values
        doubleprecision,intent(inout),dimension(:) :: b
        integer :: iopt
        integer :: pardiso_mtype
        integer n, nnz, nrhs, ldb

        integer, save    ::  info = 0
        integer*8 , save :: factors = 0

        doubleprecision,save :: total_time
!DEC$ IF DEFINED  (USE_UMF_PACK)
        ! UMFPACK constants
        type(c_ptr),save :: symbolic,numeric
        ! zero-based arrays
        real(8),save :: control(0:UMFPACK_CONTROL-1),umf_info(0:UMFPACK_INFO-1)
        doubleprecision,allocatable,dimension(:),save :: b_sol

!DEC$ ENDIF

!DEC$ IF DEFINED  (USE_PARDISO)

        INTEGER*8,save  :: pt(64)
        INTEGER,save    :: phase
        INTEGER,save    :: maxfct, mnum, error, msglvl
        INTEGER,save    :: iparm(64)
        doubleprecision,allocatable,dimension(:),save :: b_sol

        INTEGER    ,save::  idum(1)
        doubleprecision ,save::  ddum(1)

!DEC$ ENDIF

        n    = no_rows
        nnz  = no_vals
        ldb  = n
        nrhs = 1


!DEC$ IF DEFINED  (USE_UMF_PACK)
      selectcase (iopt)
      case (1)
            total_time = get_clock();
            allocate(b_sol(size(b)))

            call umf4def (control)
            call umf4sym (n,n, rowind, colptr, values, symbolic, control, umf_info)
            call umf4num (rowind, colptr, values, symbolic, numeric, control, umf_info)
            call umf4fsym (symbolic)
            !total_time =  umf_info(UMFPACK_NUMERIC_TIME)+umf_info(UMFPACK_SYMBOLIC_TIME)
            if (umf_info(UMFPACK_STATUS) .eq. 0) then
                if(TRANS_DEBUG) then
                     write (*,*) 'Factorization succeeded. Mem needed:', umf_info(UMFPACK_PEAK_MEMORY)/8.0/1024/1024 , "[MB]"
                endif
            else
                 write(*,*) 'SYS::UMF ERROR:: INFO from factorization step:', umf_info(UMFPACK_STATUS)
            endif

      case(2)
            b_sol = 0
            call umf4solr (UMFPACK_Aat, rowind, colptr, values, b_sol, b, numeric, control, umf_info)
            b  = b_sol;

            if (umf_info(UMFPACK_STATUS) .eq. 0) then
                if(TRANS_DEBUG) then
                    write (*,*) 'Solve succeeded. Time needed:',umf_info(UMFPACK_SOLVE_WALLTIME)
                endif
            else
                 write(*,*) 'SYS::UMF ERROR: INFO from solve:', umf_info(UMFPACK_STATUS)
            endif

      case(3)
            print*,"UMFPACK Solved:"
            print*,"Solve time needed:",get_clock()-total_time,"[s]"
            call umf4fnum (numeric)
            deallocate(b_sol)
      endselect

!DEC$ ELSE IF DEFINED  (USE_PARDISO)


 selectcase (iopt)
      case (1)
      allocate(b_sol(size(b)))
          total_time = get_clock();
          maxfct = 1 ! in many application this is 1
          mnum   = 1 ! same here
          iparm = 0
          iparm(1) = 1 ! no solver default
          iparm(2) = 2 ! fill-in reordering from METIS
          iparm(3) = 1 ! numbers of processors, value of OMP_NUM_THREADS
          iparm(4) = 0 ! 0 - no iterative-direct algorithm, if 1 multirecursive iterative algorithm 61, 31 para me
          iparm(5) = 0 ! no user fill-in reducing permutation
          iparm(6) = 0 ! =0 solution on the first n compoments of x
          iparm(7) = 0 ! not in use
          iparm(8) = 2 ! numbers of iterative refinement steps
          iparm(9) = 0 ! not in use
          iparm(10) = 10 ! perturbe the pivot elements with 1E-13
          iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
          iparm(12) = 0 ! not in use
          iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric).
          iparm(14) = 0 ! Output: number of perturbed pivots
          iparm(15) = 0 ! not in use
          iparm(16) = 0 ! not in use
          iparm(17) = 0 ! not in use
          iparm(18) = -1 ! Output: number of nonzeros in the factor LU
          iparm(19) = -1 ! Output: Mflops for LU factorization
          iparm(20) = 0 ! Output: Numbers of CG Iterations
          iparm(27) = 0 ! perform matrix check
          iparm(32) = 0 ! if 1 use multirecursive iterative algorithm

           error = 0 ! initialize error flag
          msglvl = 0 ! print statistical information
!          pardiso_mtype     = 11     ! complex unsymmetric matrix
          phase     = 11      ! only reordering and symbolic factorization



          CALL pardiso (pt, maxfct, mnum, pardiso_mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)

!          WRITE(*,*) 'Reordering completed ... '

          IF (error .NE. 0) THEN
            WRITE(*,*) 'SYS::PARDISO::The following ERROR was detected during the reordeing step:', error
            STOP 1
          END IF


!          WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

    !C.. Factorization.
          phase     = 22  ! only factorization
          CALL pardiso (pt, maxfct, mnum, pardiso_mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)

!          WRITE(*,*) 'Factorization completed ...  '
          IF (error .NE. 0) THEN
             WRITE(*,*) 'SYS::PARDISO::The following ERROR was detected during the factorization step:', error
            STOP 1
          ENDIF
          if(QSYS_DEBUG_LEVEL > 0) then
          WRITE(*,*) 'Peak memory usage   = ',max (IPARM(15), IPARM(16)+IPARM(17))/1024.0,"[MB]"
          endif
      case(2)
          b_sol = 0
    !C.. Back substitution and iterative refinement
          phase     = 33  ! only factorization
          iparm(8)  = 3   ! max numbers of iterative refinement steps
          CALL pardiso (pt, maxfct, mnum, pardiso_mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, b, b_sol, error)

          b  = b_sol;
          !WRITE(*,*) 'Solve completed ... '
          IF (error .NE. 0) THEN
             WRITE(*,*) 'SYS::PARDISO::The following ERROR was detected: ', error
          ENDIF

      case(3)
    !C.. Termination and release of memory
            phase     = -1           ! release internal memory
            CALL pardiso (pt, maxfct, mnum, pardiso_mtype, phase, n, ddum, idum, idum,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)
            if(QSYS_DEBUG_LEVEL > 0) then
            print*,"SYS::PARDISO::Solve time needed:",get_clock()-total_time,"[s]"
            endif
            deallocate(b_sol)
      endselect


!DEC$ ELSE
      selectcase (iopt)
      case (1)
      total_time = get_clock();
! First, factorize the matrix. The factors are stored in *factors* handle.
      !iopt = 1
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, colptr , rowind , b, ldb,factors, info )
!
      if (info .eq. 0) then
!         write (*,*) 'Factorization succeeded'
      else
         write(*,*) 'SYS::c_fortran_zgssv::SuperLU INFO from factorization step:', info
      endif
      case(2)
! Second, solve the system using the existing factors.
!      iopt = 2
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, colptr,rowind ,  b, ldb,factors, info )
!
      if (info .eq. 0) then
!         write (*,*) 'Solve succeeded'
!         write (*,*) (b(i), i=1, n)
      else
         write(*,*) 'SYS::c_fortran_zgssv::SuperLU ERROR from triangular solver:', info
      endif
      case(3)
! Last, free the storage allocated inside SuperLU
!      iopt = 3
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, colptr,rowind, b, ldb,factors, info )
      if(QSYS_DEBUG_LEVEL > 0) then
      print*,"SYS::SuperLU::Computations time:",get_clock()-total_time,"[s]"
      endif
      endselect
!DEC$ ENDIF

      endsubroutine dalg_SSOLEQ

subroutine inverse_matrix(N,A)
  integer :: N
  complex*16,dimension(:,:):: A
  complex*16,allocatable,dimension(:)  :: WORK
  integer,allocatable,dimension (:)    :: IPIV
  integer info,error


  B_SINGULAR_MATRIX = .false.
  allocate(WORK(N),IPIV(N),stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:not enough memory"
    stop
  end if
  call ZGETRF(N,N,A,N,IPIV,info)
  if(info .eq. 0) then
!    write(*,*)"succeded"
  else
    write(*,*)"SYS::LEAD::ZGETRF::failed with info:",info
    write(*,*)"           It seems your matrix is singular check your code"
    B_SINGULAR_MATRIX = .true.
  end if
  call ZGETRI(N,A,N,IPIV,WORK,N,info)
  if(info .eq. 0) then
!    write(*,*)"succeded"
  else
   write(*,*)"SYS::LEAD::ZGETRI::failed with info:",info
  end if
  deallocate(IPIV,WORK,stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:fail to release"
    stop
  end if
end subroutine inverse_matrix

subroutine alg_ZHEEV(N,Mat,Vecs,Lambdas)
    integer :: N
    complex*16 :: Mat(:,:),Vecs(:,:),Lambdas(:)

    ! -------------------------------------------------
    !                    LAPACK
    ! -------------------------------------------------
    !     .. Local Scalars ..
    INTEGER          INFO, LWORK, LWMAX

    !     .. Local Arrays ..

    DOUBLE PRECISION,allocatable :: RWORK( : )
    COMPLEX*16,allocatable        :: WORK( : )


    ! Initialize lapack and allocate arrays
    LWMAX = N*50

    allocate(RWORK ( 3*N  ))
    allocate(WORK  ( LWMAX ))

    !
    !     Query the optimal workspace.
    !
    LWORK  = -1
    call zheev("N", "L", N, Mat, N, Lambdas, work, lwork, rwork, info)


    if( INFO /= 0 ) then
        print*,"  alg_ZHEEV: Error during solving with info:",INFO
        stop
    endif
    LWORK  = MIN( LWMAX, INT( WORK( 1 ) ) )

    deallocate( WORK)
    allocate(WORK  ( LWORK ))

    Vecs = Mat
    call zheev("V", "L", N, Vecs, N, Lambdas, work, lwork, rwork, info)

    deallocate(RWORK)
    deallocate(WORK)

end subroutine alg_ZHEEV

doubleprecision function cond_number(N,A) result(rval)
  integer :: N
  complex*16,dimension(:,:):: A
  complex*16,allocatable,dimension(:)       :: WORK
  doubleprecision,allocatable,dimension(:)  :: RWORK
  integer,allocatable,dimension (:)    :: IPIV
  complex*16,dimension(:,:),allocatable:: tmpA
  integer info,error
  doubleprecision :: anorm,norm
  external zlange
  doubleprecision zlange
  B_SINGULAR_MATRIX = .false.

  allocate(WORK(2*N),IPIV(N),RWORK(2*N),tmpA(N,N),stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:not enough memory"
    stop
  end if

  tmpA = A

  anorm = ZLANGE( '1', N, N, tmpA, N, WORK )

  call ZGETRF(N,N,tmpA,N,IPIV,info)
  call zgecon( '1', N, tmpA, N, anorm, norm, work, rwork, info )
  rval = norm


  deallocate(IPIV,WORK,RWORK,tmpA,stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:fail to release"
    stop
  end if
end function cond_number

doubleprecision function cond_SVD(N,A) result(rval)
    integer :: N
    complex*16,dimension(:,:):: A
    complex*16 , allocatable :: tmpA(:,:) , VT(:,:),tU(:,:),tVT(:,:)
    doubleprecision , allocatable :: S(:) , tS(:)
    integer :: i,j,M
    allocate(tmpA(N,N))
    tmpA = A
    call ZSVD(N,tmpA,tU,tS,tVT)
    M = 0
    do i = 1 , N
    if(abs(tS(i)) > 1.0d-20 ) M = M+1
    enddo
    rval = abs(tS(1)/tS(N))
    deallocate(tmpA,tU,tS,tVT)
end function cond_SVD


subroutine ZSVD(N,A,U,S,VT)
      integer :: N
      complex*16 , dimension(:,:) :: A
      complex*16 ,allocatable , dimension(:,:) :: U,VT
      doubleprecision,allocatable,dimension(:) :: S

!*     .. Parameters ..
      INTEGER          LDA, LDU, LDVT
      INTEGER          LWMAX
      INTEGER          INFO, LWORK
      complex*16,allocatable,dimension(:) :: WORK
      doubleprecision,allocatable,dimension(:) :: RWORK

      LDA   = N
      LDU   = N
      LDVT  = N
      LWORK = -1
      LWMAX = 40*N
      if(allocated(U)) deallocate(U)
      if(allocated(S)) deallocate(S)
      if(allocated(VT)) deallocate(VT)
      allocate(U( LDU, N ), VT( LDVT, N ), S( N ),WORK(LWMAX) , RWORK(5*N))


      CALL ZGESVD( 'All', 'All', N, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK , RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      CALL ZGESVD( 'All', 'All', N, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK , RWORK, INFO )

      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF


end subroutine ZSVD

subroutine inverse_svd(N,A)
    complex*16 :: A(:,:)
    integer :: N,i,j
    complex*16 ,allocatable , dimension(:,:) :: U,VT
    doubleprecision,allocatable,dimension(:) :: S


    call ZSVD(N,A,U,S,VT)

    do i = 1, N
    do j = 1, N
        A(i,j) = sum( conjg(Vt(:,i))*(S(:)**(-1))*conjg(U(j,:)) )
    enddo
    enddo

    deallocate(U,S,VT)
end subroutine inverse_svd


subroutine solve_GGEV(H0,Tau)
    complex*16 :: H0(:,:) , Tau(:,:)
    integer :: N,M
    complex*16 , allocatable  :: U(:,:),VT(:,:) , V(:,:) , tU(:,:),tVT(:,:)
    doubleprecision , allocatable  :: S(:) , tS(:)
    doubleprecision :: rcond,eps
    integer :: i,j
    logical :: bStabilize
    N = size(H0,1)
    print*,"Perform SVD",N

    call ZSVD(N,Tau,tU,tS,tVT)
    M = 0
    do i = 1 , N
       if(abs(tS(i)) > 1.0d-20 ) M = M+1
    enddo


    allocate(S(M),U(N,M),Vt(M,N),V(N,M))
    S = tS(1:M)
    print*,S,M
    U  = tU(1:N,1:M)
    Vt = tVt(1:M,1:N)
    V  = conjg(transpose(Vt))

    do i = 1 , N
    do j = 1 , M
        U(i,j) =   U(i,j) * sqrt( S(j) )
        V(i,j) =   V(i,j) * sqrt( S(j) )
    enddo
    enddo
    rcond = cond_number(N,tau)

    eps        = 1.0D-16
    bStabilize = .false.

    if( rcond < eps  ) bStabilize = .true.

    if(bStabilize) then
        do i = 1 , N
        do j = 1 , N
            tU(i,j) = sum(U(i,:)*conjg(U(j,:))) + sum(V(i,:)*conjg(V(j,:)))
        enddo
        enddo
        tVt = H0 + II * tU
        rcond = cond_number(N,tVt)
        if( rcond < eps ) then
            print*,"SYS::ERROR::Hoping matrix baldy defined"
            stop -1
        endif

    endif

end subroutine solve_GGEV

end module modalgs
