!!!
!!! File:   modutils.f
!!! Author: mkk12259
!!!
!!! Created on 5 styczeń 2012, 19:36
!!!

        MODULE modutils
        implicit none

        INTEGER,private :: clock1


        type,public :: cnum_val
            integer :: num
            double precision :: val
        endtype cnum_val


        public:: read_filed
        public:: read_filei
        public:: write_to_file
        public:: write_to_file_num
        public:: write_to_file_fnum
        public:: write_to_file_scale
        public:: write_to_filei_scale
        public:: save_all
        public:: printDate
        public:: vector_pearson_correlation
        public:: matrix_pearson_correlation
        contains

        double precision function vector_pearson_correlation(vecA,vecB) result(rval)
            double precision, dimension(:), intent(in) :: vecA, vecB

            integer :: n,i
            doubleprecision :: minmaxA(2) , minmaxB(2) , aveA,aveB
            double precision, dimension(:), allocatable :: vA , vB
            n = size(vecA)

            allocate(vA(n))
            allocate(vB(n))
            minmaxA(:) = vecA(1)
            minmaxB(:) = vecB(1)

            do i = 1 , n
                if( minmaxA(1) > vecA(i) ) minmaxA(1) = vecA(i)
                if( minmaxA(2) < vecA(i) ) minmaxA(2) = vecA(i)

                if( minmaxB(1) > vecB(i) ) minmaxB(1) = vecB(i)
                if( minmaxB(2) < vecB(i) ) minmaxB(2) = vecB(i)
            enddo

            print*,"Min:Max A:",minmaxA
            print*,"Min:Max B:",minmaxB
            vA = vecA
            vB = vecB
            vA = (vA - minmaxA(1))/(minmaxA(2)-minmaxA(1))
            vB = (vB - minmaxB(1))/(minmaxB(2)-minmaxB(1))
            aveA = sum(vA)/n
            aveB = sum(vB)/n


            rval = sum( (vA-aveA) * (vB-aveB)   )
            rval = rval / sqrt(sum( (vA-aveA)**2  ))
            rval = rval / sqrt(sum( (vB-aveB)**2  ))

            deallocate(va)
            deallocate(vb)
        end function vector_pearson_correlation

        double precision function matrix_pearson_correlation(matA,matB) result(rval)
            double precision, dimension(:,:), intent(in) :: matA, matB
            double precision, dimension(:), allocatable :: vA , vB

            integer :: n1, n2 , n , i , j
            n1 = size(matA,1)
            n2 = size(matA,2)
            n  = n1 * n2
            allocate(vA(n))
            allocate(vB(n))
            do i = 1 , n1
            do j = 1 , n2
                vA(i + (j-1)*n1) = matA(i,j)
                vB(i + (j-1)*n1) = matB(i,j)
            enddo
            enddo
            rval = vector_pearson_correlation(va,vb);
            deallocate(va)
            deallocate(vb)
        end function matrix_pearson_correlation

        subroutine write_to_file(indx,file_name,F,NX,NY)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              double precision,dimension(:,:),intent(in)::F
              integer,intent(in)::NX,NY
              integer       ::i,j

              !print*," <file> File:",file_name
              open(unit=indx, file=file_name)

              do i = 1, NX
              do j = 1, NY

             write(indx,"(f15.8,A,f15.8,A,e15.8)"), REAL(i),"  ",REAL(j),"   ",  F(i,j)

              end do
              write(indx,*)," "
              end do
             close(indx)
          end subroutine write_to_file
          ! complex
      subroutine cwrite_to_file(indx,file_name,F,NX,NY)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              complex*16,dimension(:,:),intent(in)::F
              integer,intent(in)::NX,NY
              integer       ::i,j

              !print*," <file> File:",file_name
              open(unit=indx, file=file_name)

              do i = 1, NX
              do j = 1, NY

             write(indx,"(f15.8,A,f15.8,A,e15.8)"),REAL(i),"  ",REAL(j),"   ", abs(F(i,j))**2

              end do
              write(indx,*)," "
              end do
             close(indx)
          end subroutine cwrite_to_file

 subroutine read_filed (UnitNum, FileName, NumRows, NumCols, Array )

      integer, intent (in) :: UnitNum
      character (len=*), intent (in) :: FileName
      integer, intent (in) :: NumRows, NumCols
      double precision, dimension (1:NumRows, 1:NumCols), intent (out) :: Array

      integer :: i, j

      open (unit=UnitNum, file=FileName, status='old', action='read' )



      do i=1, NumRows
         read (UnitNum, *) (Array (i, j), j=1,NumCols)
      end do

      close (UnitNum)

      return

   end subroutine read_filed

 subroutine read_filei (UnitNum, FileName, NumRows, NumCols, Array )

      integer, intent (in) :: UnitNum
      character (len=*), intent (in) :: FileName
      integer, intent (in) :: NumRows, NumCols
      integer, dimension (1:NumRows, 1:NumCols), intent (out) :: Array

      integer :: i, j

      open (unit=UnitNum, file=FileName, status='old', action='read' )



      do i=1, NumRows
         read (UnitNum, *) (Array (i, j), j=1,NumCols)
      end do

      close (UnitNum)

      return

   end subroutine read_filei


        subroutine write_to_file_scale(indx,file_name,F,NX,NY,scaleX,scaleY)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              double precision,dimension(:,:),intent(in)::F
              integer,intent(in)::NX,NY
              double precision,intent(in)::scaleX,scaleY
              integer       ::i,j

              !print*," <file> File:",file_name
              open(unit=indx, file=file_name)

              do i = 1, NX
              do j = 1, NY

             write(indx,"(f15.8,A,f15.8,A,e15.8)"),REAL(i)*scaleX,"  ",REAL(j)*scaleY,"   ",  F(i,j)

              end do
              write(indx,*)," "
              end do
             close(indx)
        end subroutine write_to_file_scale
        subroutine write_to_filei_scale(indx,file_name,F,NX,NY,scaleX,scaleY)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              integer,dimension(:,:),intent(in)::F
              integer,intent(in)::NX,NY
              integer,intent(in)::scaleX,scaleY
              integer       ::i,j

              !print*," <file> File:",file_name
              open(unit=indx, file=file_name)

              do i = 1, NX
              do j = 1, NY

             write(indx,"(f15.8,A,f15.8,A,e15.8)"), REAL(i)*scaleX,"  ",REAL(j)*scaleY,"   ",  DBLE(F(i,j))

              end do
              write(indx,*)," "
              end do
             close(indx)
        end subroutine write_to_filei_scale

      subroutine iwrite_to_file(indx,file_name,F,NX,NY)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              integer,dimension(:,:),intent(in)::F
              integer,intent(in)::NX,NY
              integer       ::i,j

              !print*," <file> File:",file_name
              open(unit=indx, file=file_name)

              do i = 1, NX
              do j = 1, NY

             write(indx,"(f15.8,A,f15.8,A,e15.8)"), REAL(i),"  ",REAL(j),"   ", REAL(F(i,j))

              end do
              write(indx,*)," "
              end do
             close(indx)
          end subroutine iwrite_to_file

         subroutine write_to_file_fnum(indx,file_name,val,F,NX,NY)
            integer,intent(in)::indx
            character(len=*),intent(in) :: file_name
            double precision, intent(in):: val
            double precision,dimension(:,:),intent(in)::F

            character(20) :: filename
!            character(*) :: name
            character(10) :: tmp
            integer,intent(in)::NX,NY
            integer       ::i


            i = len(file_name)
            filename = CHAR(0)! file_name
            filename(1:i) = file_name
            write(tmp,"(f10.6)"),val


            filename(i+1:i+8) = tmp
            filename(i+9:i+12) = ".txt"

            !allocate(name(i+9))
!            name = filename(1:i+9)
            print*,"name=",filename(1:i+12)
!            print*,"name=",name
            !deallocate(name)


            call write_to_file(indx,filename(1:i+12),F,NX,NY)


        end subroutine write_to_file_fnum

        ! numerowany zapis do pliku
        subroutine write_to_file_num(indx,file_name,iter,F,NX,NY)
            integer,intent(in)::indx
            character(len=*),intent(in) :: file_name
            integer, intent(in) :: iter
            double precision,dimension(:,:),intent(in)::F

            character(20) :: filename
!            character(*) :: name
            character(10) :: tmp
            integer,intent(in)::NX,NY
            integer       ::i


            i = len(file_name)
            filename = CHAR(0)! file_name
            filename(1:i) = file_name
            write(tmp,"(i5)"),iter+10000
            filename(i+1:i+5) = tmp
            filename(i+6:i+9) = ".txt"

            !allocate(name(i+9))
!            name = filename(1:i+9)
            print*,"name=",filename(1:i+9)
!            print*,"name=",name
            !deallocate(name)


            call write_to_file(indx,filename(1:i+9),F,NX,NY)


        end subroutine write_to_file_num

        subroutine write_cevecs(indx,file_name,F,NX,NY,nvals,ijkToL)
            integer,intent(in)::indx
            character(len=*),intent(in) :: file_name
            complex*16,dimension(:,:),intent(in)::F
            integer,intent(in)::NX,NY
            integer,intent(in)::nvals
            interface
            integer function ijkToL(ix,jx,kx) result(L)
            integer,intent(in)::ix,jx,kx

            end function ijkToL
            end interface

            integer       ::i,j,k
          print*," <file> File:",file_name
          open(unit=indx,file=file_name)
            do i=1,NX
            do j=1,NY

            write(indx,"(i10 A i10)",advance='no'),i,"    ",j
            do k = 1,nvals
                write(indx,"(A,f15.6)",advance='no'),"   ", abs(F(ijkToL(i,NY/2,j),k))
                enddo
                write(indx,*),""
            enddo
                write(indx,*),""
            enddo
        close(indx)


        end subroutine write_cevecs

        subroutine write_to_3dfile(indx,file_name,F,NX,NY,NZ)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              double precision,dimension(:,:,:),intent(in)::F
              integer,intent(in)::NX,NY,NZ
              integer       ::i,j,k

              print*," <file> Zapis w formacie 3D:",file_name
              open(unit=indx, file=file_name)
              write(indx,"(3i15)"),NX,NY,NZ
              do i = 1, NX
              do j = 1, NY
              do k = 1, NZ
                 write(indx,"(e15.8)"),F(i,j,k)
              end do
              end do
              end do


              close(indx)
          end subroutine write_to_3dfile

          subroutine write_to_filei(indx,file_name,F,NX,NY)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              integer,dimension(:,:),intent(in)::F
              integer,intent(in)::NX,NY
              integer       ::i,j

              print*," <file> File:",file_name
              open(unit=indx, file=file_name)

              do i = 1, NX
              do j = 1, NY

             write(indx,"(i15,A,i15,A,i15)"),  i,"  ",j,"   ",  F(i,j)

              end do
              write(indx,*)," "
              end do
              close(indx)
          end subroutine write_to_filei

        subroutine write_mapped_file(indx,file_name,INDICES,VEC)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              integer,dimension(:,:),intent(in)::INDICES
              double precision,dimension(:),intent(in)::VEC
              integer       ::i,j,nx,ny

              !print*," <file> File:",file_name
              open(unit=indx, file=file_name)
              nx = size(indices,1)
              ny = size(indices,2)

              do i = 1, NX
              do j = 1, NY
               if(indices(i,j) > 0) then
               write(indx,"(f15.8,A,f15.8,A,e15.8)"),REAL(i),"  ",REAL(j),"   ",  VEC(indices(i,j))
               else
               write(indx,"(f15.8,A,f15.8,A,e15.8)"),REAL(i),"  ",REAL(j),"   ",  DBLE(0.0)
               endif
              end do
              write(indx,*)," "
              end do
             close(indx)
          end subroutine write_mapped_file

          subroutine save_all(indx,file_name,F,NX,NY,NZ)
              integer,intent(in)::indx
              character(len=*),intent(in) :: file_name
              double precision,dimension(:,:,:),intent(in)::F
              integer,intent(in)::NX,NY,NZ
              integer       ::i,j,k

              print*," <file> File:",file_name
              open(unit=indx, file=file_name)

              do k = 1, NZ

              do i = 1, NX
              do j = 1, NY


             write(indx,"(e15.8,A)",advance='no'), F(i,j,k)  ,"    "

              end do
              write(indx,*),""
              end do


              enddo
              close(indx)
          end subroutine save_all


          real function clock() result(t)
            INTEGER :: clock_rate,c_time

            CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate)
            CALL SYSTEM_CLOCK(COUNT=c_time)
            t = REAL(c_time)/clock_rate

          end function clock

          doubleprecision function my_clock() result(t)
            INTEGER :: clock_rate,c_time

            CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate)
            CALL SYSTEM_CLOCK(COUNT=c_time)
            t = dble(c_time)/clock_rate

          end function my_clock

          subroutine reset_clock()
            CALL SYSTEM_CLOCK(COUNT=clock1)
          end subroutine reset_clock

          real function get_clock() result(c)
            INTEGER :: clock_rate,c_time
            CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate)
            CALL SYSTEM_CLOCK(COUNT=c_time)
            c = (real(c_time) - clock1)/clock_rate
          end function get_clock


          subroutine  printDate
              character(8)  :: date
              character(10) :: time
              character(5)  :: zone
              integer,dimension(8) :: values

              ! using keyword arguments
              call date_and_time(date,time,zone,values)
              call date_and_time(DATE=date,ZONE=zone)
              call date_and_time(TIME=time)
              call date_and_time(VALUES=values)

              print '(a,2x,a,A,2x,a)'," DATA:", date,"   TIME:", time


          end subroutine printDate


          ! ===============================================
          ! Tworzenie prostych plikow VTK
          ! ===============================================
          subroutine modutils_3darray2VTK(inArray,filename)
              double precision,intent(in),dimension(:,:,:)::inArray
              character(len=*),intent(in) :: filename

              character(len=40) :: text
              integer   :: tlen,nx,ny,nz,i,j,k

              nx = size(inArray,1)
              ny = size(inArray,2)
              nz = size(inArray,3)

              text = ''
              tlen = len(filename)+1
              text(1:tlen) = filename
              text(tlen:tlen+4) = ".vtk"
              print*,"modutils: VTK file=",text
              open(unit=111,file=text)
              write(111,"(A)"), "# vtk DataFile Version 2.0"
              write(111,"(A)"), filename
              write(111,"(A)"), "ASCII"
              write(111,"(A)"), "DATASET STRUCTURED_POINTS"
              write(111,"(A,I10,A,I10,A,I10)"),    "DIMENSIONS ",nx," ",ny," ",nz
              write(111,"(A)"), "ASPECT_RATIO 1 1 1"
              write(111,"(A)"), "ORIGIN 0 0 0"
              write(111,"(A)"), "ORIGIN 0 0 0"
              write(111,"(A,I10)"), "POINT_DATA ",nx*ny*nz
              write(111,"(A)"), "SCALARS volume_scalars float 1"
              write(111,"(A)"), "LOOKUP_TABLE default"


                do k=1,nz
                do j =1,ny
                do i =1,nx
                  write(111,"(e16.6,A)",advance='no'),inArray(i,j,k)," "
                enddo
                enddo
                    write(111,*)," "
                enddo
              close(111)

          endsubroutine modutils_3darray2VTK


 SUBROUTINE zgaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np
      complex*16 a(np,np),b(np,mp)
      integer,PARAMETER :: NMAX=5000
      INTEGER :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      complex*16 :: big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge. abs(big))then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END subroutine zgaussj

        END MODULE modutils

