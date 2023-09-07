
!------------------------------------------------------------------------
! Reads the input parameters, which are given from standard input.
! The input is temporarily written in a work file, excluding the comments.
! Then, that temporary file is read in.
!------------------------------------------------------------------------
subroutine pinput2(maxnzone, maxnr, &
  re, ratc, ratl, &
  tlen, np, omegai, imin, imax, &
  nzone, vrmin, vrmax, rho, vsv, vsh, qmu, &
  r0, eqlat, eqlon, mt, nr, theta, phi, lat, lon, output)
!------------------------------------------------------------------------
  implicit none
  integer :: maxnzone, maxnr  ! Maximum number of zones and stations.
  integer :: np  ! Number of points in frequency domain.
  integer :: imin, imax  ! Index of minimum and maximum frequency.
  integer :: nzone, nr
  real(8) :: tlen, omegai, re, ratc, ratl
  real(8), dimension(*) :: vrmin, vrmax, qmu
  real(8), dimension(4,*) :: rho, vsv, vsh
  real(8) :: r0, eqlat, eqlon, eqlattmp
  real(8), dimension(3,3) :: mt
  real(8), dimension(*) :: theta, phi, lat, lon
  character(len=80), dimension(*) :: output
  integer :: i
  character(len=80) :: dummy, tmpfile

  data tmpfile / 'workSH.txt' /

! Temporary file open
  open(unit=11, file=tmpfile, status='unknown')
! Writing to the temporary file
  do
    read(5,'(a80)') dummy
    if (dummy(1:1) == 'c') cycle
    if (dummy(1:3) == 'end') exit
    write(11,'(a80)') dummy
  end do

! Temporary file close
  close(11)

! Temporary file open
  open(unit=11, file=tmpfile, status='unknown')
! Reading the parameter
  read(11,*) tlen, np
  read(11,*) re      ! relative error (vertical grid)
  read(11,*) ratc    ! ampratio (vertical grid cut-off)
  read(11,*) ratl    ! ampratio (for l-cutoff)
  read(11,*) omegai  ! omegai
  omegai = -dlog(omegai) / tlen

  read(11,*) imin, imax

  read(11,*) nzone
  if (nzone > maxnzone) stop 'nzone is too large. (pinput)'
  do i = 1, nzone
    read(11,*) vrmin(i), vrmax(i), &
      rho(1,i), rho(2,i), rho(3,i), rho(4,i), &
      vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), &
      vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), qmu(i)
  end do

! Source parameter
  read(11,*) r0, eqlat, eqlon
  eqlattmp = eqlat
  call translat(eqlattmp, eqlattmp)
  read(11,*) mt(1,1), mt(1,2), mt(1,3), mt(2,2), mt(2,3), mt(3,3)

! Stations
  read(11,*) nr
  if (nr > maxnr) stop 'nr is too large. (pinput)'
  do i = 1, nr
    read(11,*) lat(i), lon(i)
    call translat(lat(i), lat(i))
    call calthetaphi(eqlattmp, eqlon, lat(i), lon(i), theta(i), phi(i))
  end do

  do i = 1, nr
    read(11,'(a80)') output(i)
  end do

! Temporary file close
  close(11)

  return
end subroutine


!------------------------------------------------------------------------
! Converts geodetic latitude to geocentric latitude.
!------------------------------------------------------------------------
subroutine translat(geodetic, geocentric)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: flattening = 1.d0 / 298.25d0
  real(8), parameter :: pi = 3.1415926535897932d0
  real(8), intent(in) :: geodetic  ! Input geodetic latitude.
  real(8), intent(out) :: geocentric  ! Output geocentric latitude.
  real(8) :: latitude  ! Latitude variable for use in computation.

  if (geodetic > 90.d0) stop 'Latitude is over 90. (pinput)'

  ! degrees to radians
  latitude = geodetic / 180.d0 * pi
  ! gedetic to geocentric
  latitude = atan((1.d0 - flattening) * (1.d0 - flattening) * tan(latitude))
  ! radians to degrees
  geocentric = latitude * 180.d0 / pi

  return
end subroutine


!------------------------------------------------------------------------
! Computes distance and azimuth from event to station.
!------------------------------------------------------------------------
subroutine calthetaphi(iEvLat, iEvLon, iStLat, iStLon, theta, phi)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0
  real(8), intent(in) :: iEvLat, iEvLon, iStLat, iStLon  ! Input latitudes and longitudes of event and station (in degrees).
  real(8) :: evColat, evLon, stColat, stLon  ! Colatitudes and longitudes of event and station (in radians).
  real(8), intent(out) :: theta, phi  ! Resulting distance and azimuth.
  real(8) :: gcarc, azimuth
  real(8) :: cosAzimuth, sinAzimuth

  ! transformation to spherical coordinates
  evColat = 90.d0 - iEvLat
  stColat = 90.d0 - iStLat

  ! degrees to radians
  evColat = evColat / 180.d0 * pi
  evLon = iEvLon / 180.d0 * pi
  stColat = stColat / 180.d0 * pi
  stLon = iStLon / 180.d0 * pi

  gcarc = acos(cos(evColat) * cos(stColat) + sin(evColat) * sin(stColat) * cos(evLon - stLon))

  cosAzimuth = (cos(stColat) * sin(evColat) - sin(stColat) * cos(evColat) * cos(stLon - evLon)) / sin(gcarc)
  sinAzimuth = sin(stColat) * sin(stLon - evLon) / sin(gcarc)

  azimuth = acos(cosAzimuth)
  if (sinAzimuth < 0.d0) azimuth = -1.d0 * azimuth

  ! radians to degrees
  theta = gcarc * 180.d0 / pi
  phi = azimuth * 180.d0 / pi

  phi = 180.d0 - phi
  return
end subroutine


!------------------------------------------------------------------------
! Function to compute a + bx + cx^2 + dx^3, where x = r/R
!------------------------------------------------------------------------
subroutine valueatradius(coefficients, radius, rmax, result)
  implicit none
  real(8), intent(in) :: coefficients(4)  ! Coefficients of cubic function. [a, b, c, d] in a + bx + cx^2 + dx^3.
  real(8), intent(in) :: radius  ! r
  real(8), intent(in) :: rmax  ! R
  real(8), intent(out) :: result
  integer :: j
  real(8) :: x_n  ! Power of x = r/R.
  real(8) :: accumulatedValue  ! Variable to store the temporary result of accumulation.

  x_n = 1.d0
  accumulatedValue = coefficients(1)
  do j = 2, 4
    x_n = x_n * (radius / rmax)
    accumulatedValue = accumulatedValue + coefficients(j) * x_n
  end do

  result = accumulatedValue
  return
end subroutine


!------------------------------------------------------------------------
! Computes vertical wavenumber k_z at each zone. (See section 3.2 of Kawai et al. 2006.)
!------------------------------------------------------------------------
subroutine computeKz(nzone, vrmin, vrmax, vs, rmin, rmax, imax, lmin, tlen, kz)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nzone  ! Number of zones.
  integer, intent(in) :: imax  ! Index of maximum frequency.
  integer, intent(in) :: lmin  ! Smallest angular order l.
  real(8), intent(in) :: vrmin(:), vrmax(:), vs(:,:), rmin, rmax, tlen
  real(8), intent(out) :: kz(:)  ! Computed value of vertical wavenumber k_z at each zone.
  integer :: izone
  real(8) :: v(4), vs1, vs2, vmin, rh, omega, kx, gtmp

  do izone = 1, nzone
    v(:) = vs(:, izone)
    ! compute Vs at bottom (vs1) and top (vs2) of zone
    call valueatradius(v, vrmin(izone), rmax, vs1)
    call valueatradius(v, vrmax(izone), rmax, vs2)
    ! zone thickness
    rh = vrmax(izone) - vrmin(izone)
    ! compute omega
    omega = 2.d0 * pi * real(imax,8) / tlen
    ! get smaller Vs value (This is to get larger k_z value.)
    vmin = min(vs1, vs2)
    ! compute k_x (eq. 30 of Kawai et al. 2006)
    kx = (real(lmin,8) + 0.5d0) / vrmax(izone)
    ! compute k_z^2 (eq. 32 of Kawai et al. 2006)
    gtmp = ( omega * omega ) / ( vmin * vmin ) - ( kx * kx )
    ! compute k_z (When it is not real, it is set to 0.)
    if ( gtmp > 0.d0 ) then
      kz(izone) = sqrt(gtmp)
    else
      kz(izone) = 0.d0
    endif
  enddo

  return
end subroutine computekz


!------------------------------------------------------------------------
! Deciding the distribution of grid points.
!------------------------------------------------------------------------
subroutine computeGridRadii(maxnlayer, maxnzone, nzone, kz, vrmin, vrmax, rmin, rmax, re, nlayer, nlayerinzone, radii)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: maxnlayer, maxnzone
  integer, intent(in) :: nzone  ! Number of zones.
  real(8), intent(in) :: kz(*)  ! Vertical wavenumber k_z at each zone.
  real(8), intent(in) :: vrmin(*), vrmax(*), rmin, rmax
  real(8), intent(in) :: re  ! Desired relative error.
  integer, intent(out) :: nlayer  ! Total number of grid points.
  integer, intent(out) :: nlayerinzone(maxnzone)  ! Number of grid points in each zone.
  real(8), intent(out) :: radii(maxnlayer + maxnzone + 1)  ! Radius at each grid point.
  integer :: izone, igrid, i, ntmp
  real(8) :: rh

! initializing variables
  radii = 0.d0
  nlayerinzone = 0

! computing the distribution of grid points
  igrid = 1
  radii(1) = rmin
  do izone = 1, nzone
    ! zone thickness
    rh = vrmax(izone) - vrmin(izone)
    ! decide the number of layers in this zone
    if (kz(izone) == 0.d0) then
      ntmp = 1
    else
      ! rh / dz = rh * (lambda_z / dz) / lambda_z = rh * sqrt(3.3 / re) * (k_z / 2 pi)
      !  (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
      !  The "/0.7 +1" is to increase the number of grids a bit.
      ntmp = int(sqrt(3.3d0 / re) * rh * kz(izone) / 2.d0 / pi / 7.d-1 + 1)
    endif
    nlayerinzone(izone) = min(ntmp, 5)
    ! compute radius at each grid point
    do i = 1, nlayerinzone(izone)
      igrid = igrid + 1
      radii(igrid) = vrmin(izone) + dble(i) * rh / dble(nlayerinzone(izone))
    enddo
  enddo

! recounting the total number of grid points
  nlayer = sum(nlayerinzone)

  return
end subroutine computeGridRadii



