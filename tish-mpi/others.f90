
!------------------------------------------------------------------------
! Reads the input parameters, which are given from standard input.
! The input is temporarily written in a work file, excluding the comments.
! Then, that temporary file is read in.
!------------------------------------------------------------------------
subroutine readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegai, imin, imax, &
  nZone, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, qmuOfZone, &
  r0, eqlat, eqlon, mt, nReceiver, theta, phi, lat, lon, output)
!------------------------------------------------------------------------
  implicit none
  character(len=80), parameter :: tmpfile = 'workSH.txt'

  integer, intent(in) :: maxNZone, maxNReceiver  ! Maximum number of zones and receivers.
  real(8), intent(out) :: tlen  ! Time length.
  integer, intent(out) :: np  ! Number of points in frequency domain.
  real(8), intent(out) :: re  ! Desired relative error due to vertical gridding.
  real(8), intent(out) :: ratc  ! Amplitude ratio for vertical grid cut-off.
  real(8), intent(out) :: ratl  ! Amplitude ratio for angular order cut-off.
  real(8), intent(out) :: omegai  ! omegai
  integer, intent(out) :: imin, imax  ! Index of minimum and maximum frequency.
  integer, intent(out) :: nZone  ! Number of zones.
  real(8), intent(out) :: rminOfZone(:), rmaxOfZone(:)  ! Lower and upper radii of each zone.
  real(8), intent(out) :: rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho, vsv, and vsh structure.
  real(8), intent(out) :: qmuOfZone(:)  ! Qmu of each zone.
  integer, intent(out) :: nReceiver  ! Number of receivers.
  real(8), intent(out) :: r0, eqlat, eqlon, mt(3,3)
  real(8), intent(out) :: theta(:), phi(:), lat(:), lon(:)
  character(len=80), intent(out) :: output(:)

  integer :: i
  character(len=80) :: dummy
  real(8) :: eqlattmp

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

  ! Reading parameters
  read(11,*) tlen, np
  read(11,*) re      ! relative error (vertical grid)
  read(11,*) ratc    ! ampratio (vertical grid cut-off)
  read(11,*) ratl    ! ampratio (for l-cutoff)
  read(11,*) omegai  ! omegai
  omegai = -dlog(omegai) / tlen
  read(11,*) imin, imax  ! index of minimum and maximum frequency

  ! Earth structure
  read(11,*) nZone
  if (nZone > maxNZone) stop 'nZone is too large. (pinput)'
  do i = 1, nZone
    read(11,*) rminOfZone(i), rmaxOfZone(i), &
      rhoPolynomials(1,i), rhoPolynomials(2,i), rhoPolynomials(3,i), rhoPolynomials(4,i), &
      vsvPolynomials(1,i), vsvPolynomials(2,i), vsvPolynomials(3,i), vsvPolynomials(4,i), &
      vshPolynomials(1,i), vshPolynomials(2,i), vshPolynomials(3,i), vshPolynomials(4,i), qmuOfZone(i)
  end do

  ! Source parameter
  read(11,*) r0, eqlat, eqlon
  eqlattmp = eqlat
  call translat(eqlattmp, eqlattmp)
  read(11,*) mt(1,1), mt(1,2), mt(1,3), mt(2,2), mt(2,3), mt(3,3)

  ! Receivers
  read(11,*) nReceiver
  if (nReceiver > maxNReceiver) stop 'nReceiver is too large. (pinput)'
  do i = 1, nReceiver
    read(11,*) lat(i), lon(i)
    call translat(lat(i), lat(i))
    call calthetaphi(eqlattmp, eqlon, lat(i), lon(i), theta(i), phi(i))
  end do

  do i = 1, nReceiver
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

  if (geodetic < -90.d0 .or. 90.d0 < geodetic) stop 'Latitude is out of range. (pinput)'

  ! degrees to radians
  latitude = geodetic / 180.d0 * pi
  ! gedetic to geocentric
  latitude = atan((1.d0 - flattening) * (1.d0 - flattening) * tan(latitude))
  ! radians to degrees
  geocentric = latitude * 180.d0 / pi

  return
end subroutine


!------------------------------------------------------------------------
! Computes distance and azimuth from source to receiver.
!------------------------------------------------------------------------
subroutine calthetaphi(iEvLat, iEvLon, iStLat, iStLon, theta, phi)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  real(8), intent(in) :: iEvLat, iEvLon, iStLat, iStLon  ! Input latitudes and longitudes of source and receiver (in degrees).
  real(8) :: evColat, evLon, stColat, stLon  ! Colatitudes and longitudes of source and receiver (in radians).
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
subroutine valueAtRadius(coefficients, radius, rmax, result)
!------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: coefficients(4)  ! Coefficients of cubic function. [a, b, c, d] in a + bx + cx^2 + dx^3.
  real(8), intent(in) :: radius  ! r : The radius to compute the value at.
  real(8), intent(in) :: rmax  ! R: Maximum radius of region considered.
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
subroutine computeKz(nZone, rminOfZone, rmaxOfZone, vsPolynomials, rmax, imax, lmin, tlen, kz)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: imax  ! Index of maximum frequency.
  integer, intent(in) :: lmin  ! Smallest angular order l.
  real(8), intent(in) :: rminOfZone(:), rmaxOfZone(:)  ! Lower and upper radii of each zone.
  real(8), intent(in) :: vsPolynomials(:,:)  ! Polynomial functions of vs structure.
  real(8), intent(in) :: rmax  ! Maximum radius of region considered.
  real(8), intent(in) :: tlen  ! Time length.
  real(8), intent(out) :: kz(:)  ! Computed value of vertical wavenumber k_z at each zone.
  integer :: iZone
  real(8) :: v(4), vs1, vs2, vmin, omega, kx, kz2

  do iZone = 1, nZone
    v(:) = vsPolynomials(:, iZone)
    ! compute Vs at bottom (vs1) and top (vs2) of zone
    call valueAtRadius(v, rminOfZone(iZone), rmax, vs1)
    call valueAtRadius(v, rmaxOfZone(iZone), rmax, vs2)
    ! get smaller Vs value (This is to get larger k_z value.)
    vmin = min(vs1, vs2)
    ! compute largest omega (This is to get larger k_z value.)
    omega = 2.d0 * pi * dble(imax) / tlen
    ! compute smallest k_x (See eq. 30 of Kawai et al. 2006.) (This is to get larger k_z value.)
    kx = (dble(lmin) + 0.5d0) / rmaxOfZone(iZone)
    ! compute k_z^2 (See eq. 32 of Kawai et al. 2006.)
    kz2 = (omega ** 2) / (vmin ** 2) - (kx ** 2)
    ! compute k_z (When it is not real, it is set to 0.)
    if (kz2 > 0.d0) then
      kz(iZone) = sqrt(kz2)
    else
      kz(iZone) = 0.d0
    end if
  end do

  return
end subroutine


!------------------------------------------------------------------------
! Deciding the distribution of grid points.
!------------------------------------------------------------------------
subroutine computeGridRadii(nZone, kz, rminOfZone, rmaxOfZone, rmin, re, nLayer, nLayerInZone, gridRadii)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: kz(:)  ! Vertical wavenumber k_z at each zone.
  real(8), intent(in) :: rminOfZone(:), rmaxOfZone(:)  ! Lower and upper radii of each zone.
  real(8), intent(in) :: rmin  ! Minimum radius of region considered.
  real(8), intent(in) :: re  ! Desired relative error due to vertical gridding.
  integer, intent(out) :: nLayer  ! Total number of layers (= number of grid points + 1).
  integer, intent(out) :: nLayerInZone(:)  ! Number of layers in each zone.
  real(8), intent(out) :: gridRadii(:)  ! Radius at each grid point.
  integer :: iZone, iGrid, i, nTemp
  real(8) :: rh

! initializing variables
  gridRadii = 0.d0
  nLayerInZone = 0

! computing the distribution of grid points
  iGrid = 1
  gridRadii(1) = rmin
  do iZone = 1, nZone
    ! zone thickness
    rh = rmaxOfZone(iZone) - rminOfZone(iZone)
    ! decide the number of layers in this zone
    if (kz(iZone) == 0.d0) then
      nTemp = 1
    else
      ! rh / dz = rh * (lambda_z / dz) / lambda_z = rh * sqrt(3.3 / re) * (k_z / 2 pi)
      !  (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
      !  The "/0.7 +1" is to increase the number of grids a bit.
      nTemp = int(sqrt(3.3d0 / re) * rh * kz(iZone) / 2.d0 / pi / 7.d-1 + 1)
    end if
    nLayerInZone(iZone) = min(nTemp, 5)
    ! compute radius at each grid point
    do i = 1, nLayerInZone(iZone)
      iGrid = iGrid + 1
      gridRadii(iGrid) = rminOfZone(iZone) + dble(i) * rh / dble(nLayerInZone(iZone))
    end do
  end do

! recounting the total number of layers
  nLayer = sum(nLayerInZone)

  return
end subroutine


!------------------------------------------------------------------------
! Computing the indices of the first grid point and the first (iLayer, k', k)-pair in each zone.
!------------------------------------------------------------------------
subroutine computeFirstIndices(nZone, nLayerInZone, iFirstGrid, jFirstComponent)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: nLayerInZone(:)  ! Number of layers in each zone.
  integer, intent(out) :: iFirstGrid(:)  ! Index of the first grid point in each zone.
  integer, intent(out) :: jFirstComponent(:)  ! Index of the first component of (iLayer, k', k)-pairs in each zone.
  integer :: i

  iFirstGrid(1) = 1
  jFirstComponent(1) = 1
  do i = 1, nZone - 1
    iFirstGrid(i+1) = iFirstGrid(i) + nLayerInZone(i)
    jFirstComponent(i+1) = jFirstComponent(i) + 4 * nLayerInZone(i)
  end do

end subroutine


!------------------------------------------------------------------------
! Computing the source position.
!------------------------------------------------------------------------
subroutine computeSourcePosition(nLayer, rmaxOfZone, rmin, rmax, gridRadii, isp, r0, iZoneOfSource, qLayerOfSource)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayer  ! Total number of layers.
  real(8), intent(in) :: rmaxOfZone(:)  ! Upper radius of each zone.
  real(8), intent(in) :: rmin, rmax  ! Minimum and maximum radii of region considered.
  real(8), intent(in) :: gridRadii(:)  ! Radii of grid points.
  integer, intent(in) :: isp(:)  ! Index of the first grid point in each zone.
  real(8), intent(inout) :: r0  ! Source radius. Its value may be fixed in this subroutine.
  integer, intent(out) :: iZoneOfSource  ! Which zone the source is in.
  real(8), intent(out) :: qLayerOfSource  ! A double-value index of source position in its zone.
  !::::::::::::::::::::::::::::::::::::::::: (0 at bottom of zone, nlayer(iZone) at top of zone.)
  integer :: iLayer  ! Index of layer. (1 at rmin, nLayer+1 at rmax.)
  real(8) :: xLayerOfSource  ! A double-value index of source position. (0 at rmin, nLayer at rmax.)

  ! checking the parameter
  if (r0 < rmin .or. rmax < r0) stop 'The source location is improper.(calspo)'

  ! computing a double-value index of source position
  if (r0 == rmax) then
    xLayerOfSource = dble(nLayer) - 0.01d0
    r0 = gridRadii(nLayer) + (xLayerOfSource - dble(nLayer-1)) * (gridRadii(nLayer+1) - gridRadii(nLayer))
    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Note that radii(nLayer+1) = rmax.
    !!!TODO (xLayerOfSource - dble(nLayer-1)) = 0.99d0 ?

  else
    ! find the layer that the source is in (Note that the index of the lowermost layer is 1, not 0.)
    iLayer = 1
    do
      if (r0 < gridRadii(iLayer + 1)) exit
      iLayer = iLayer + 1
    end do

    xLayerOfSource = dble(iLayer - 1) + (r0 - gridRadii(iLayer)) / (gridRadii(iLayer + 1) - gridRadii(iLayer))

    ! fix source position when it is too close to a grid point
    if ((xLayerOfSource - dble(iLayer-1)) < 0.01d0) then
      xLayerOfSource = dble(iLayer-1) + 0.01d0
      r0 = gridRadii(iLayer) + (xLayerOfSource - dble(iLayer-1)) * (gridRadii(iLayer+1) - gridRadii(iLayer))
      !!!TODO (xLayerOfSource - dble(iLayer-1)) = 0.01d0 ?
    elseif ((xLayerOfSource - dble(iLayer-1)) > 0.99d0) then
      xLayerOfSource = dble(iLayer-1) + 0.99d0
      r0 = gridRadii(iLayer) + (xLayerOfSource - dble(iLayer-1)) * (gridRadii(iLayer+1) - gridRadii(iLayer))
      !!!TODO (xLayerOfSource - dble(iLayer-1)) = 0.99d0 ?
    end if
  end if

  ! find the zone that the source is in
  iZoneOfSource = 1
  do
    if (r0 <= rmaxOfZone(iZoneOfSource)) exit
    iZoneOfSource = iZoneOfSource + 1
  end do

  ! compute the double-value index of source position in its zone
  qLayerOfSource = xLayerOfSource - dble(isp(iZoneOfSource) - 1)

end subroutine


!------------------------------------------------------------------------
! Compute grids for the source based on input parameters.
!------------------------------------------------------------------------
subroutine computeSourceGrid(isp, gridRadii, r0, iZoneOfSource, qLayerOfSource, gridRadiiForSource)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: isp(:)  ! Index of the first grid point in each zone.
  real(8), intent(in) :: gridRadii(:)  ! Radii of grid points.
  real(8), intent(in) :: r0  ! Input source radius.
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  real(8), intent(in) :: qLayerOfSource  ! A double-value index of source position in its zone.
  !:::::::::::::::::::::::::::::::::::::::: (0 at bottom of zone, nlayer(iZone) at top of zone.)
  real(8), intent(out) :: gridRadiiForSource(:)  ! Radii to use for source-related computations.
  integer :: iLayer

  ! find the layer that the source is in
  iLayer = isp(iZoneOfSource) + int(qLayerOfSource)  ! Note that int(x) rounds down the value x.

  ! assign grid radii
  gridRadiiForSource(1) = gridRadii(iLayer)  ! radius of bottom of layer
  gridRadiiForSource(2) = r0 ! radius of source
  gridRadiiForSource(3) = gridRadii(iLayer + 1) ! radius of top of layer

end subroutine


!------------------------------------------------------------------------
! Computing variable values at grid points.
!------------------------------------------------------------------------
subroutine computeStructureValues(nZone, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, &
  nLayerInZone, gridRadii, nValue, valuedRadii, rhoValues, ecLValues, ecNValues)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmax  ! Maximum radius of region considered.
  integer, intent(in) :: nLayerInZone(:)  ! Number of layers in each zone.
  real(8), intent(in) :: rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho, vsv, and vsh structure.
  real(8), intent(in) :: gridRadii(:)  ! Radii of grid points.
  integer, intent(out) :: nValue  ! Total number of values for each variable.
  real(8), intent(out) :: valuedRadii(:)  ! Radii corresponding to each variable value.
  real(8), intent(out) :: rhoValues(:), ecLValues(:), ecNValues(:)  ! Values of rho, L, and N at each point
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: (with 2 values at boundaries).
  real(8) :: rhoTemp, vsvTemp, vshTemp
  integer :: iZone, iLayer, iValue, iGrid

  ! Initializing the data
  valuedRadii = 0.d0
  rhoValues = 0.d0
  ecLValues = 0.d0
  ecNValues = 0.d0
  iValue = 0
  iGrid = 0

  ! Computing variable values at grid points
  do iZone = 1, nZone
    do iLayer = 1, nLayerInZone(iZone) + 1
      iValue = iValue + 1
      iGrid = iGrid + 1
      valuedRadii(iValue) = gridRadii(iGrid)

      ! Evaluating the density and elastic constants at this point
      call valueAtRadius(rhoPolynomials(:, iZone), valuedRadii(iValue), rmax, rhoTemp)
      call valueAtRadius(vsvPolynomials(:, iZone), valuedRadii(iValue), rmax, vsvTemp)
      call valueAtRadius(vshPolynomials(:, iZone), valuedRadii(iValue), rmax, vshTemp)
      rhoValues(iValue) = rhoTemp
      ecLValues(iValue) = rhoValues(iValue) * vsvTemp**2
      ecNValues(iValue) = rhoValues(iValue) * vshTemp**2
    end do

    iGrid = iGrid - 1
  end do

  nValue = iValue

end subroutine


!------------------------------------------------------------------------
! Computing variable values near the source.
!------------------------------------------------------------------------
subroutine computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, gridRadiiForSource, &
  rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, mu0)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  real(8), intent(in) :: rmax  ! Maximum radius of region considered.
  real(8), intent(in) :: rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho, vsv, and vsh structure.
  real(8), intent(in) :: gridRadiiForSource(:)  ! Radii to use for source-related computations.
  real(8), intent(out) :: rhoValuesForSource(:), ecLValuesForSource(:), ecNValuesForSource(:), mu0
  real(8) :: rhoTemp, vsvTemp, vshTemp
  integer :: i

  ! Computing the structure grid points
  do i = 1, 3
    ! Evaluating the density and elastic constants at this point
    call valueAtRadius(rhoPolynomials(:, iZoneOfSource), gridRadiiForSource(i), rmax, rhoTemp)
    call valueAtRadius(vsvPolynomials(:, iZoneOfSource), gridRadiiForSource(i), rmax, vsvTemp)
    call valueAtRadius(vshPolynomials(:, iZoneOfSource), gridRadiiForSource(i), rmax, vshTemp)
    rhoValuesForSource(i) = rhoTemp
    ecLValuesForSource(i) = rhoValuesForSource(i) * vsvTemp**2
    ecNValuesForSource(i) = rhoValuesForSource(i) * vshTemp**2
  end do

  mu0 = ecLValuesForSource(2)

end subroutine


!------------------------------------------------------------------------
! Initialize complex vector of size n with zeros.
!------------------------------------------------------------------------
subroutine initComplexVector(n, b)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: n
  complex(8), intent(out) :: b(n)
  integer :: i

  ! Initialize vector 'b' with zeros
  do i = 1, n
    b(i) = dcmplx(0.0d0, 0.0d0)
  end do

end subroutine


!------------------------------------------------------------------------
! Initialize complex matrix of shape (n1, n2) with zeros.
!------------------------------------------------------------------------
subroutine initComplexMatrix(n1, n2, a)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: n1, n2
  complex(8), intent(out) :: a(n1, n2)
  integer :: i, j

  ! Initialize matrix 'a' with zeros
  do j = 1, n2
    do i = 1, n1
      a(i, j) = dcmplx(0.d0, 0.d0)
    end do
  end do

end subroutine


!------------------------------------------------------------------------
! Computes the lsuf parameter based on given input.
!------------------------------------------------------------------------
subroutine computeLsuf(omega, nZone, rmaxOfZone, vsvPolynomials, lsuf)
!------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: omega  ! Angular frequency.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmaxOfZone(nZone)  ! Upper radii of each zone.
  real(8), intent(in) :: vsvPolynomials(4, nZone)  ! Polynomial functions of vsv structure.
  integer, intent(out) :: lsuf  !!TODO probably The angular order that is sufficient to compute the slowest phase velocity.
  real(8) :: vsAtSurface

  ! Compute vs at planet surface
  call valueAtRadius(vsvPolynomials(:, nZone), 1.d0, 1.d0, vsAtSurface)

  ! Calculate lsuf (See eq. 29 of Kawai et al. 2006.)
  ! The slowest velocity (vs at surface) and largest radius (planet radius) is used to gain larger bound of angular order.
  lsuf = int(omega * rmaxOfZone(nZone) / vsAtSurface - 0.5d0) + 1

end subroutine


!------------------------------------------------------------------------
! Computes the coefficient based on the given input.   !!TODO probably coefficient of attenuation
!------------------------------------------------------------------------
subroutine computeCoef(nZone, omega, qmuOfZone, coef)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: omega  ! Angular frequency.
  real(8), intent(in) :: qmuOfZone(nZone)  ! Qmu of each zone.
  complex(8), intent(out) :: coef(nZone)  !!TODO probably Coefficient derived from attenuation for each zone.
  real(8) :: aa, bb
  integer :: iZone

  ! Calculate coefficients
  do iZone = 1, nZone
    if (omega == 0.d0) then
      aa = 1.d0
    else
      aa = 1.d0 + log(omega / (2.d0 * pi)) / (pi * qmuOfZone(iZone))
    end if
    bb = 1.d0 / (2.d0 * qmuOfZone(iZone))
    coef(iZone) = dcmplx(aa, bb) * dcmplx(aa, bb)
  end do

end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
! Accumulates the value of 'u' at a certain receiver for a certain (l, m)-pair (= for a certain trial function).
! (See eq. 1 of Kawai et al. 2006.)
! The trial function is specified by (k=k_max (at surface of planet), l, m, 3 (the T spherical harmonic)).
! (See eqs. 12 & 13 of Kawai et al. 2006.)
!------------------------------------------------------------------------
subroutine computeU(c0, l, trialFunctionValues, u)
!------------------------------------------------------------------------
  implicit none

  complex(8), intent(in) :: c0  ! Expansion coefficent corresponding to this trial function (k=k_max, l, m, 3).
  real(8), intent(in) :: l  ! Angular order.
  complex(8), intent(in) :: trialFunctionValues(3)  ! Trial function term. The coefficient 1/largeL is not multiplied yet.
  complex(8), intent(inout) :: u(3)
  real(8) :: largeL

  ! compute largeL (See the part after eq. 12 of Kawai et al. 2006.)
  largeL = sqrt(dble(l * (l + 1)))

  ! Accumulate value of u. (See eq. 1 of Kawai et al. 2006.)
  ! The coefficient 1/largeL, which is not included in the trial function term, is multiplied here. (See eq. 12 of Kawai et al. 2006.)
  u(1) = dcmplx(0.d0, 0.d0)
  u(2) = u(2) + c0 * trialFunctionValues(2) / dcmplx(largeL, 0.d0)
  u(3) = u(3) + c0 * trialFunctionValues(3) / dcmplx(largeL, 0.d0)

end subroutine
