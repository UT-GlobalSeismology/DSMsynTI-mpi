
!------------------------------------------------------------------------
! Reads the input parameters, which are given from standard input.
! The input is temporarily written in a work file, excluding the comments.
! Then, that temporary file is read in.
!------------------------------------------------------------------------
subroutine readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegaI, imin, imax, nZone, rminOfZone, rmaxOfZone, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
  r0, eqlat, eqlon, mt, nReceiver, lat, lon, theta, phi, output)
!------------------------------------------------------------------------
  implicit none
  character(len=80), parameter :: tmpfile = 'workPSV.txt'

  integer, intent(in) :: maxNZone, maxNReceiver  ! Maximum number of zones and receivers.
  real(8), intent(out) :: tlen  ! Time length [s].
  integer, intent(out) :: np  ! Number of points in frequency domain.
  real(8), intent(out) :: re  ! Desired relative error due to vertical gridding.
  real(8), intent(out) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8), intent(out) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8), intent(out) :: omegaI  ! Imaginary part of angular frequency for artificial damping [1/s].
  integer, intent(out) :: imin, imax  ! Index of minimum and maximum frequency.
  integer, intent(out) :: nZone  ! Number of zones.
  real(8), intent(out) :: rminOfZone(*), rmaxOfZone(*)  ! Lower and upper radii of each zone [km].
  real(8), intent(out) :: rhoPolynomials(4,*), vpvPolynomials(4,*), vphPolynomials(4,*)
  real(8), intent(out) :: vsvPolynomials(4,*), vshPolynomials(4,*), etaPolynomials(4,*)
  !:::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta structure.
  real(8), intent(out) :: qmuOfZone(*), qkappaOfZone(*)  ! Qmu and Qkappa of each zone.
  real(8), intent(out) :: r0, eqlat, eqlon, mt(3,3)  ! Depth [km], coordinates [deg], and moment tensor [10^25 dyn cm] of source.
  integer, intent(out) :: nReceiver  ! Number of receivers.
  real(8), intent(out) :: lat(*), lon(*)  ! Coordinates [deg] of receivers.
  real(8), intent(out) :: theta(*), phi(*)  ! Colatitude and longitude of receivers with event at north pole [rad].
  character(len=80), intent(out) :: output(*)  ! Output file names.

  integer :: i
  character(len=80) :: dummy

  ! Open temporary file.
  open(unit=11, file=tmpfile, status='unknown')

  ! Write to the temporary file.
  do
    read(5,'(a80)') dummy
    if (dummy(1:1) == 'c' .or. dummy(1:1) == '!') cycle
    if (dummy(1:3) == 'end') exit
    write(11,'(a80)') dummy
  end do

  ! Close temporary file.
  close(11)

  ! Re-open temporary file.
  open(unit=11, file=tmpfile, status='unknown')

  ! Read parameters.
  read(11,*) tlen, np
  read(11,*) re      ! relative error (vertical grid)
  read(11,*) ratc    ! ampratio (vertical grid cut-off)
  read(11,*) ratl    ! ampratio (for l-cutoff)
  read(11,*) omegaI  ! adamp (for artificial damping)
  omegaI = -log(omegaI) / tlen  ! omegai [1/s]
  read(11,*) imin, imax  ! index of minimum and maximum frequency

  ! earth structure
  read(11,*) nZone
  if (nZone > maxNZone) stop 'nZone is too large. (readInput)'
  do i = 1, nZone
    read(11,*) rminOfZone(i), rmaxOfZone(i), &
      rhoPolynomials(1,i), rhoPolynomials(2,i), rhoPolynomials(3,i), rhoPolynomials(4,i)
    read(11,*) vpvPolynomials(1,i), vpvPolynomials(2,i), vpvPolynomials(3,i), vpvPolynomials(4,i)
    read(11,*) vphPolynomials(1,i), vphPolynomials(2,i), vphPolynomials(3,i), vphPolynomials(4,i)
    read(11,*) vsvPolynomials(1,i), vsvPolynomials(2,i), vsvPolynomials(3,i), vsvPolynomials(4,i)
    read(11,*) vshPolynomials(1,i), vshPolynomials(2,i), vshPolynomials(3,i), vshPolynomials(4,i)
    read(11,*) etaPolynomials(1,i), etaPolynomials(2,i), etaPolynomials(3,i), etaPolynomials(4,i), qmuOfZone(i), qkappaOfZone(i)
  end do

  ! source parameter
  read(11,*) r0, eqlat, eqlon
  read(11,*) mt(1,1), mt(1,2), mt(1,3), mt(2,2), mt(2,3), mt(3,3)

  ! receivers
  read(11,*) nReceiver
  if (nReceiver > maxNReceiver) stop 'nReceiver is too large. (readInput)'
  do i = 1, nReceiver
    read(11,*) lat(i), lon(i)
    call computeThetaPhi(eqlat, eqlon, lat(i), lon(i), theta(i), phi(i))
  end do

  do i = 1, nReceiver
    read(11,'(a80)') output(i)
  end do

  ! Close temporary file.
  close(11)

  return
end subroutine


!------------------------------------------------------------------------!!Common
! Converts geodetic latitude to geocentric latitude.
!------------------------------------------------------------------------
subroutine transformLatitude(geodetic, geocentric)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: flattening = 1.d0 / 298.25d0
  real(8), parameter :: pi = 3.1415926535897932d0

  real(8), intent(in) :: geodetic  ! Input geodetic latitude [deg].
  real(8), intent(out) :: geocentric  ! Output geocentric latitude [deg].
  real(8) :: latitude  ! Latitude variable for use in computation.

  if (geodetic < -90.d0 .or. 90.d0 < geodetic) stop 'Latitude is out of range. (transformLatitude)'

  ! degrees to radians
  latitude = geodetic / 180.d0 * pi
  ! gedetic to geocentric
  latitude = atan((1.d0 - flattening) * (1.d0 - flattening) * tan(latitude))
  ! radians to degrees
  geocentric = latitude * 180.d0 / pi

  return
end subroutine


!------------------------------------------------------------------------!!Common
! Computes the colatitude (theta) and longitude (phi) of a receiver
! when the source is shifted to the north pole.
! Note that the longitude of the original source is set as 0 after the shift,
! so the shifted longitude of receiver [rad] is (pi - azimuth).
!------------------------------------------------------------------------
subroutine computeThetaPhi(iEvLat, iEvLon, iStLat, iStLon, theta, phi)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  real(8), intent(in) :: iEvLat, iEvLon, iStLat, iStLon  ! Input latitudes and longitudes of source and receiver [deg].
  real(8), intent(out) :: theta, phi  ! Colatitude and longitude of receiver with event at north pole [rad].
  real(8) :: evColat, evLon, stColat, stLon  ! Colatitudes and longitudes of source and receiver [rad].
  real(8) :: cosAlpha, sinAlpha
  real(8) :: tmp

  ! Transform geographic latitudes [deg] to geocentric colatitudes [rad].
  call transformLatitude(iEvLat, tmp)
  evColat = (90.d0 - tmp) / 180.d0 * pi
  call transformLatitude(iStLat, tmp)
  stColat = (90.d0 - tmp) / 180.d0 * pi

  ! Transform longitudes from degrees to radians.
  evLon = iEvLon / 180.d0 * pi
  stLon = iStLon / 180.d0 * pi

  ! Compute epicentral distance [rad], which will directly be the colatitude of receiver after shift.
  cosAlpha = cos(evColat) * cos(stColat) + sin(evColat) * sin(stColat) * cos(evLon - stLon)
  if (1.d0 < cosAlpha) cosAlpha = 1.d0
  if (cosAlpha < -1.d0) cosAlpha = -1.d0
  theta = acos(cosAlpha)

  ! Compute shifted longitude of receiver [rad], which is (pi - azimuth).
  if (theta == 0.d0) then
    phi = 0.d0
  else
    cosAlpha = (cos(stColat) * sin(evColat) - sin(stColat) * cos(evColat) * cos(stLon - evLon)) / sin(theta)
    if (1.d0 < cosAlpha) cosAlpha = 1.d0
    if (cosAlpha < -1.d0) cosAlpha = -1.d0
    sinAlpha = sin(stColat) * sin(stLon - evLon) / sin(theta)
    ! pi - azimuth
    if (sinAlpha >= 0.d0) then
      phi = pi - acos(cosAlpha)
    else
      phi = pi + acos(cosAlpha)
    end if
  end if

  return
end subroutine

