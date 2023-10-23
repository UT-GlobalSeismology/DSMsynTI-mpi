module FileIO
  implicit none

  interface writeSPCFile
    module procedure write1dble, write2dble, write3dble
    module procedure write3int
    module procedure write1int2dble
  end interface

contains

  subroutine openSPCFile(fileName, unitNum, spcFormat, append)
    implicit none
    character(len=80), intent(in) :: fileName  ! Name of output spc file.
    integer, intent(in) :: unitNum  ! Unit number of output spc file.
    integer, intent(in) :: spcFormat  ! Format of output spc file (0:binary, 1:ascii).
    integer, intent(in) :: append  ! Whether to append to existing file (0:create, 1:append).

    if (spcFormat == 0 .and. append == 0) then
      open(unit=unitNum, file=fileName, status='unknown', &
        form='unformatted', access='stream', convert='big_endian')
    else if (spcFormat == 0 .and. append == 1) then
      open(unit=unitNum, file=fileName, position='append', status='unknown', &
        form='unformatted', access='stream', convert='big_endian')
    else if (spcFormat == 1 .and. append == 0) then
      open(unit=unitNum, file=fileName, status='unknown')
    else if (spcFormat == 1 .and. append == 1) then
      open(unit=unitNum, file=fileName, position='append', status='unknown')
    else
      write(*,*) "WARNING: spcFormat must be 0 or 1"
    end if
  end subroutine

  subroutine write1dble(unitNum, spcFormat, dble1)
    implicit none
    integer, intent(in) :: unitNum  ! Unit number of output spc file.
    integer, intent(in) :: spcFormat  ! Format of output spc file (0:binary, 1:ascii).
    real(8), intent(in) :: dble1

    if (spcFormat == 0) then
      write(unitNum) dble1
    else if (spcFormat == 1) then
      write(unitNum,*) dble1
    else
      write(*,*) "WARNING: spcFormat must be 0 or 1"
    end if
  end subroutine

  subroutine write2dble(unitNum, spcFormat, dble1, dble2)
    implicit none
    integer, intent(in) :: unitNum  ! Unit number of output spc file.
    integer, intent(in) :: spcFormat  ! Format of output spc file (0:binary, 1:ascii).
    real(8), intent(in) :: dble1, dble2

    if (spcFormat == 0) then
      write(unitNum) dble1, dble2
    else if (spcFormat == 1) then
      write(unitNum,*) dble1, dble2
    else
      write(*,*) "WARNING: spcFormat must be 0 or 1"
    end if
  end subroutine

  subroutine write3dble(unitNum, spcFormat, dble1, dble2, dble3)
    implicit none
    integer, intent(in) :: unitNum  ! Unit number of output spc file.
    integer, intent(in) :: spcFormat  ! Format of output spc file (0:binary, 1:ascii).
    real(8), intent(in) :: dble1, dble2, dble3

    if (spcFormat == 0) then
      write(unitNum) dble1, dble2, dble3
    else if (spcFormat == 1) then
      write(unitNum,*) dble1, dble2, dble3
    else
      write(*,*) "WARNING: spcFormat must be 0 or 1"
    end if
  end subroutine

  subroutine write3int(unitNum, spcFormat, int1, int2, int3)
    implicit none
    integer, intent(in) :: unitNum  ! Unit number of output spc file.
    integer, intent(in) :: spcFormat  ! Format of output spc file (0:binary, 1:ascii).
    integer, intent(in) :: int1, int2, int3

    if (spcFormat == 0) then
      write(unitNum) int1, int2, int3
    else if (spcFormat == 1) then
      write(unitNum,*) int1, int2, int3
    else
      write(*,*) "WARNING: spcFormat must be 0 or 1"
    end if
  end subroutine

  subroutine write1int2dble(unitNum, spcFormat, int1, dble1, dble2)
    implicit none
    integer, intent(in) :: unitNum  ! Unit number of output spc file.
    integer, intent(in) :: spcFormat  ! Format of output spc file (0:binary, 1:ascii).
    integer, intent(in) :: int1
    real(8), intent(in) :: dble1, dble2

    if (spcFormat == 0) then
      write(unitNum) int1, dble1, dble2
    else if (spcFormat == 1) then
      write(unitNum,*) int1, dble1, dble2
    else
      write(*,*) "WARNING: spcFormat must be 0 or 1"
    end if
  end subroutine

  subroutine closeSPCFile(unitNum)
    implicit none
    integer, intent(in) :: unitNum  ! Unit number of output spc file.

    close(unitNum)
  end subroutine

end module
