! PARSER - parse an ini file
!     V1.0 6 AUG 03
!     Written by Douglas S. Elder elderdo@yahoo.com

!** common block that will contain the iniFile
      BLOCK DATA iniFile
      LOGICAL initialized
      INTEGER CNT
      CHARACTER LINES(100)*256
      CHARACTER iniFilename*256
      COMMON /Options/ CNT, LINES, initialized, iniFilename
      DATA initialized, CNT, iniFilename / .FALSE., 0, 'Options.ini' /
      END

      module modinip
      contains



!** allow you to use a different iniFile or to switch to another
      SUBROUTINE setIniFilename(value)
      CHARACTER value*(*)
      LOGICAL initialized
      INTEGER CNT
      CHARACTER LINES(100)*256
      CHARACTER iniFilename*256
      COMMON /Options/ CNT, LINES, initialized, iniFilename
      if (initialized .EQV. .TRUE.) then
        if (value .NE. iniFilename) then
      		iniFilename = value
!****** switching to a different ini file
                call loadOptions

        end if
      else
!****** overriding the default ini file
        iniFilename = value
      end if
      END SUBROUTINE setIniFilename

!***** search for the Section and keyword, if found return its value
!***** otherwise return the default
      SUBROUTINE getValue(section, kwd, value, default)
      IMPLICIT NONE
      CHARACTER section*(*)
      CHARACTER kwd*(*)
      CHARACTER value*(*)
      CHARACTER default*(*)

!**      WRITE(*,*) 'getValue', section, kwd, value, default
      value = ''
      call getValueX(section, kwd, value)
!**      WRITE(*,*) 'Got value ', value
      if (value .EQ. '') then
	 value = default
      end if
      END SUBROUTINE getValue

!**** read in the iniFile into the common Options block
      SUBROUTINE loadOptions
      IMPLICIT NONE
      INTEGER CNT
      LOGICAL initialized
      CHARACTER LINES(100)*256
      CHARACTER iniFilename*256
      COMMON /Options/ CNT, LINES, initialized, iniFilename
      CHARACTER LINE*256
!**      WRITE(*,*) 'loadOptions'
      CNT = 0
      OPEN(UNIT=33, FILE=iniFilename)
1     READ(33,'(A)', END=10) LINE
      CNT = CNT + 1
      if (CNT .GT. 100) then
!**         WRITE(0,*) 'Options.ini file > 100 lines.'
	 STOP 16
      else
	      LINES(CNT) = LINE
      end if
      GOTO 1
10    CONTINUE
!**      WRITE(*,*) 'CNT = ',CNT
      CLOSE (UNIT=33)
      initialized = .TRUE.
      END SUBROUTINE loadOptions


!*** try to find the Section and keyword, if found return its value
!*** otherwise return an empty string
      SUBROUTINE getValueX(section, kwd, value)
      IMPLICIT NONE
      CHARACTER section*(*)
      CHARACTER kwd*(*)
      CHARACTER value*(*)
      INTEGER I, J, STARTVAL
      INTEGER CNT
      LOGICAL initialized
      CHARACTER LINES(100)*256
      CHARACTER iniFilename*256
      COMMON /Options/ CNT, LINES, initialized, iniFilename
      INTEGER MAXLINE
      PARAMETER (MAXLINE = 256)

      LOGICAL foundSection, foundKwd
      if (initialized .EQV. .FALSE.) then
            call loadOptions
      end if
      foundSection = .FALSE.
      foundKwd = .FALSE.
      value = ''
!**      WRITE(*,*) 'Looking for ',section
      DO I=1, CNT
      	if (LINES(I)(1:1) .EQ. '[') then
		DO J=2, MAXLINE
		  if (LINES(I)(J:J) .EQ. ']') then
			if (section .EQ. LINES(I)(2:J-1)) then
				foundSection = .TRUE.
!**				Write(*,*) 'foundSection=', foundSection
				GOTO 15
			end if
		  end if
		ENDDO
	else
		if (foundSection .EQV. .TRUE.) then
                        STARTVAL = 0
!**			WRITE(*,*) 'Looking for keyword ', kwd
			DO J=1, MAXLINE
				if (LINES(I)(J:J) .EQ. '=') then
					if (kwd .EQ. LINES(I)(1:J-1)) then
						STARTVAL = J+1
						foundKwd = .TRUE.
!**						Write(*,*) 'found keyword ', kwd
						GOTO 13
					end if
				end if
			END DO
13			CONTINUE
			if (foundKwd .EQV. .TRUE.) then
				DO J=MAXLINE, STARTVAL, -1
!**					Write(*,*) 'Looking for value at ', STARTVAL, J
					if (LINES(I)(J:J) .NE. '') then
						value = LINES(I)(STARTVAL:J)
!**						Write(*,*) 'Got value at line ', I, ' pos ', STARTVAL, J
						GOTO 20
					end if
				END DO
			end if
		end if
	end if
15	CONTINUE
      ENDDO
20    CONTINUE
      END SUBROUTINE getValueX
      !> \brief
      !!  Pobiera z pliku wartosc double
       subroutine getDoubleValue(section, kwd, var)
      IMPLICIT NONE
      CHARACTER section*(*)
      CHARACTER kwd*(*)
      double precision :: var
      character*256     :: value
      character*15     :: text

      call getValueX(section ,kwd ,value)

      text(:) = " "
      text = kwd

      read(value,"(e16.6)"),var;
      print"(5A,f16.6)"," [",section,":",text,"]=",var

      end subroutine getDoubleValue
       subroutine getSDoubleValue(section, kwd, var)
      IMPLICIT NONE
      CHARACTER section*(*)
      CHARACTER kwd*(*)
      double precision :: var
      character*256     :: value
      character*15     :: text

      call getValueX(section ,kwd ,value)

      text(:) = " "
      text = kwd

      read(value,"(e16.6)"),var;
      print"(5A,e16.6)"," [",section,":",text,"]=",var

      end subroutine getSDoubleValue

       subroutine getIntValue(section, kwd, var)
      IMPLICIT NONE
      CHARACTER section*(*)
      CHARACTER kwd*(*)
      integer :: var
      character*256     :: value
      character*15     :: text

      call getValueX(section ,kwd ,value)

      text(:) = " "
      text = kwd

      read(value,"(I16)"),var;
      print"(5A,I16)"," [",section,":",text,"]=",var

      end subroutine getIntValue

      end module modinip

