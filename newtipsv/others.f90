
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
  if (sin(theta) == 0.d0) then
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


!------------------------------------------------------------------------
! Decide if each zone is solid or liquid.
!------------------------------------------------------------------------
subroutine judgeSolidOrLiquid(nZone, vsPolynomials, phaseOfZone, nZoneSolid, nZoneLiquid)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: vsPolynomials(4,nZone)  ! Polynomial functions of vs structure [km/s].
  integer, intent(out) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: liquid).
  integer, intent(out) :: nZoneSolid, nZoneLiquid  ! Number of solid and liquid zones.
  integer :: i

  nZoneSolid = 0
  nZoneLiquid = 0

  do i = 1, nZone
    if (vsPolynomials(1,i) == 0.0d0 .and. vsPolynomials(2,i) == 0.0d0 &
      .and. vsPolynomials(3,i) == 0.0d0 .and. vsPolynomials(4,i) == 0.0d0) then
      ! liquid
      nZoneLiquid = nZoneLiquid + 1
      phaseOfZone(i) = 2
    else
      ! solid
      nZoneSolid = nZoneSolid + 1
      phaseOfZone(i) = 1
    end if
  end do

end subroutine


!------------------------------------------------------------------------!!Common
! Function to compute a + bx + cx^2 + dx^3, where x = r/R.
!------------------------------------------------------------------------
subroutine valueAtRadius(coefficients, radius, rmax, result)
!------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: coefficients(4)  ! Coefficients of cubic function. [a, b, c, d] in a + bx + cx^2 + dx^3.
  real(8), intent(in) :: radius  ! r : The radius to compute the value at [km].
  real(8), intent(in) :: rmax  ! R: Maximum radius of region considered [km].
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
subroutine computeKz(nZone, rminOfZone, rmaxOfZone, phaseOfZone, vpPolynomials, vsPolynomials, rmax, imax, lmin, tlen, kzAtZone)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: liquid).
  real(8), intent(in) :: vpPolynomials(4,nZone), vsPolynomials(4,nZone)  ! Polynomial functions of vp and vs structure [km/s].
  real(8), intent(in) :: rmax  ! Maximum radius of region considered [km].
  integer, intent(in) :: imax  ! Index of maximum frequency.
  integer, intent(in) :: lmin  ! Smallest angular order l.
  real(8), intent(in) :: tlen  ! Time length [s].
  real(8), intent(out) :: kzAtZone(*)  ! Computed value of vertical wavenumber k_z at each zone [1/km].
  integer :: iZone
  real(8) :: v(4), vBottom, vTop, vmin, omega, kx, kz2

  do iZone = 1, nZone
    ! Use Vs in solid, Vp in liquid.
    if (phaseOfZone(iZone) == 1) then
      v(:) = vsPolynomials(:, iZone)
    else
      v(:) = vpPolynomials(:, iZone)
    end if
    ! Compute velocity [km/s] at bottom and top of zone.
    call valueAtRadius(v, rminOfZone(iZone), rmax, vBottom)
    call valueAtRadius(v, rmaxOfZone(iZone), rmax, vTop)
    ! Get smaller velocity value [km/s]. (This is to get larger k_z value.)
    vmin = min(vBottom, vTop)
    ! largest omega [1/s] (This is to get larger k_z value.)
    omega = 2.d0 * pi * dble(imax) / tlen
    ! smallest k_x [1/km] (See eq. 30 of Kawai et al. 2006.) (This is to get larger k_z value.)
    kx = (dble(lmin) + 0.5d0) / rmaxOfZone(iZone)
    ! k_z^2 [1/km^2] (See eq. 32 of Kawai et al. 2006.)
    kz2 = (omega ** 2) / (vmin ** 2) - (kx ** 2)
    ! k_z [1/km] (When it is not real, it is set to 0.)
    if (kz2 > 0.d0) then
      kzAtZone(iZone) = sqrt(kz2)
    else
      kzAtZone(iZone) = 0.d0
    end if
  end do

  return
end subroutine


!------------------------------------------------------------------------
! Deciding the distribution of grid points.
!------------------------------------------------------------------------
subroutine computeGridRadii(nZone, kzAtZone, rminOfZone, rmaxOfZone, phaseOfZone, rmin, re, &
  nGrid, nLayerInZone, nLayerSolid, nLayerLiquid, gridRadii)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: kzAtZone(nZone)  ! Vertical wavenumber k_z at each zone [1/km].
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: liquid).
  real(8), intent(in) :: rmin  ! Minimum radius of region considered [km].
  real(8), intent(in) :: re  ! Desired relative error due to vertical gridding.
  integer, intent(out) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(out) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  integer, intent(out) :: nLayerSolid, nLayerLiquid  ! Number of layers in solid and liquid regions.
  real(8), intent(out) :: gridRadii(*)  ! Radius at each grid point [km].
  integer :: iZone, iGrid, i, nTemp
  real(8) :: rh

  ! Compute the distribution of grid points.
  iGrid = 1
  nLayerSolid = 0
  nLayerLiquid = 0
  gridRadii(1) = rmin
  do iZone = 1, nZone
    ! zone thickness [km]
    rh = rmaxOfZone(iZone) - rminOfZone(iZone)
    ! Decide the number of layers in this zone.
    if (kzAtZone(iZone) == 0.d0) then
      ! We usually do not compute for the evanescent regime
      !  (unless they can be seen on the surface, which is the case of shallow sources).
      nTemp = 1
    else
      ! rh / dz = rh * (lambda_z / dz) / lambda_z = rh * sqrt(3.3 / re) * (k_z / 2 pi)
      !  (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
      !  The "/0.7 +1" is to increase the number of grids a bit.
      nTemp = int(sqrt(3.3d0 / re) * rh * kzAtZone(iZone) / 2.d0 / pi / 7.d-1 + 1)
    end if
    nLayerInZone(iZone) = max(nTemp, 5)
    ! Accumulate number of layers in solid and liquid regions.
    if (phaseOfZone(iZone) == 1) then
      nLayerSolid = nLayerSolid + nLayerInZone(iZone)
    else
      nLayerLiquid = nLayerLiquid + nLayerInZone(iZone)
    end if
    ! Compute radius at each grid point [km].
    do i = 1, nLayerInZone(iZone)
      iGrid = iGrid + 1
      gridRadii(iGrid) = rminOfZone(iZone) + dble(i) * rh / dble(nLayerInZone(iZone))
    end do
  end do

  ! Register the total number of grid points.
  nGrid = iGrid

  return
end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------




!------------------------------------------------------------------------
! Computing the source position.
!------------------------------------------------------------------------
subroutine computeSourcePosition(nGrid, rmaxOfZone, phaseOfZone, gridRadii, r0, iZoneOfSource, iLayerOfSource, oRowOfSource)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nGrid  ! Total number of grid points.
  real(8), intent(in) :: rmaxOfZone(*)  ! Upper radius of each zone [km].
  integer, intent(in) :: phaseOfZone(*)  ! Phase of each zone (1: solid, 2: liquid).
  real(8), intent(in) :: gridRadii(*)  ! Radii of grid points [km].
  real(8), intent(inout) :: r0  ! Source radius [km]. Its value may be fixed in this subroutine.
  integer, intent(out) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(out) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(out) :: oRowOfSource  ! Index of the first row in the vector of (iLayer, k', k)-pairs for the source layer.
  integer :: iLayer  ! Index of layer. (1 at rmin, nGrid-1 just below rmax.)
  real(8) :: xLayerOfSource  ! A double-value index of source position. (1 at rmin, nGrid at rmax.)

  ! Compute a double-value index of source position.
  if (r0 >= gridRadii(nGrid)) then
    ! Fix source position when it is at planet surface.
    xLayerOfSource = dble(nGrid) - 0.01d0
    r0 = gridRadii(nGrid) - 0.01d0 * (gridRadii(nGrid) - gridRadii(nGrid-1))

  else
    ! Find the layer that the source is in. (Note that the index of the lowermost layer is 1, not 0.)
    iLayer = 1
    do
      if (r0 < gridRadii(iLayer + 1)) exit
      iLayer = iLayer + 1
    end do

    ! Compute the double-value index of source position.
    xLayerOfSource = dble(iLayer) + (r0 - gridRadii(iLayer)) / (gridRadii(iLayer + 1) - gridRadii(iLayer))

    ! Fix source position when it is too close to a grid point.
    if ((xLayerOfSource - dble(iLayer)) < 0.01d0) then
      xLayerOfSource = dble(iLayer) + 0.01d0
      r0 = gridRadii(iLayer) + 0.01d0 * (gridRadii(iLayer + 1) - gridRadii(iLayer))
    elseif ((xLayerOfSource - dble(iLayer)) > 0.99d0) then
      xLayerOfSource = dble(iLayer) + 0.99d0
      r0 = gridRadii(iLayer) + 0.99d0 * (gridRadii(iLayer + 1) - gridRadii(iLayer))
    end if
  end if

  ! Find the zone that the source is in.
  iZoneOfSource = 1
  do
    if (r0 < rmaxOfZone(iZoneOfSource)) exit
    iZoneOfSource = iZoneOfSource + 1
  end do
  if (phaseOfZone(iZoneOfSource) /= 1) stop 'The source is in a liquid zone. (computeSourcePosition)'

  ! Find the layer that the source is in.
  iLayerOfSource = int(xLayerOfSource)  ! Note that int(x) rounds down the value x.
  ! Find the first index of (iLayer, k', k)-pair corresponding to the layer that the source is in.
  oRowOfSource = 4 * iLayerOfSource - 3

end subroutine


!------------------------------------------------------------------------
! Computing variable values at grid points.
!------------------------------------------------------------------------
subroutine computeStructureValues(nZone, rmax, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, nLayerInZone, gridRadii, &
  nValue, valuedRadii, rhoValues, kappaValues, ecKxValues, ecKyValues, ecKzValues, ecLValues, ecNValues)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmax  ! Maximum radius of region considered [km].
  real(8), intent(in) :: rhoPolynomials(4,nZone), vpvPolynomials(4,nZone), vphPolynomials(4,nZone)
  real(8), intent(in) :: vsvPolynomials(4,nZone), vshPolynomials(4,nZone), etaPolynomials(4,nZone)
  !:::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta structure.
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  real(8), intent(in) :: gridRadii(*)  ! Radii of grid points [km].
  integer, intent(out) :: nValue  ! Total number of values for each variable.
  real(8), intent(out) :: valuedRadii(*)  ! Radii corresponding to each variable value [km].
  real(8), intent(out) :: rhoValues(*), kappaValues(*), ecKxValues(*), ecKyValues(*), ecKzValues(*), ecLValues(*), ecNValues(*)
  ! Values of rho [g/cm^3], L, and N [10^10 dyn/cm^2 = GPa]
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: at each point (with 2 values at boundaries).  !!!!TODO
  real(8) :: rhoTemp, vpvTemp, vphTemp, vsvTemp, vshTemp, etaTemp, ecA, ecC, ecF
  integer :: iZone, iLayer, iValue, iGrid

  ! Initialize variables.
  iValue = 0
  iGrid = 0

  ! Compute variable values at grid points.
  do iZone = 1, nZone
    do iLayer = 1, nLayerInZone(iZone) + 1
      iValue = iValue + 1
      iGrid = iGrid + 1
      valuedRadii(iValue) = gridRadii(iGrid)

      ! Evaluate the density and elastic constants at this point.
      call valueAtRadius(rhoPolynomials(:, iZone), valuedRadii(iValue), rmax, rhoTemp)
      call valueAtRadius(vpvPolynomials(:, iZone), valuedRadii(iValue), rmax, vpvTemp)
      call valueAtRadius(vphPolynomials(:, iZone), valuedRadii(iValue), rmax, vphTemp)
      call valueAtRadius(vsvPolynomials(:, iZone), valuedRadii(iValue), rmax, vsvTemp)
      call valueAtRadius(vshPolynomials(:, iZone), valuedRadii(iValue), rmax, vshTemp)
      call valueAtRadius(etaPolynomials(:, iZone), valuedRadii(iValue), rmax, etaTemp)
      rhoValues(iValue) = rhoTemp
      ecLValues(iValue) = rhoTemp * vsvTemp * vsvTemp
      ecNValues(iValue) = rhoTemp * vshTemp * vshTemp
      ecA = rhoTemp * vphTemp * vphTemp
      ecC = rhoTemp * vpvTemp * vpvTemp
      ecF = etaTemp * (ecA - 2.d0 * ecLValues(iValue))
      kappaValues(iValue) = (4.d0 * ecA + ecC + 4.d0 * ecF - 4.d0 * ecNValues(iValue) ) / 9.d0
      ecKxValues(iValue) = ecA - 4.d0 / 3.d0 * ecNValues(iValue)
      ecKyValues(iValue) = ecF + 2.d0 / 3.d0 * ecNValues(iValue)
      ecKzValues(iValue) = (ecC + 2.d0 * ecF) / 3.d0
    end do

    iGrid = iGrid - 1
  end do

  nValue = iValue

end subroutine


!------------------------------------------------------------------------
! Computing variable values near the source.
!------------------------------------------------------------------------
subroutine computeSourceStructureValues(iZoneOfSource, r0, rmax, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, ecC0, ecF0, ecL0)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  real(8), intent(in) :: r0  ! Input source radius [km].
  real(8), intent(in) :: rmax  ! Maximum radius of region considered [km].
  real(8), intent(in) :: rhoPolynomials(4,*), vpvPolynomials(4,*), vphPolynomials(4,*)
  real(8), intent(in) :: vsvPolynomials(4,*), vshPolynomials(4,*), etaPolynomials(4,*)
  !:::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta structure.
  real(8), intent(out) :: ecC0, ecF0, ecL0  ! Elastic moduli C, F, L at source position [10^10 dyn/cm^2 = GPa].
  real(8) :: rhoTemp, vpvTemp, vphTemp, vsvTemp, vshTemp, etaTemp, ecA0

  ! Evaluate the density and elastic constants at source.
  call valueAtRadius(rhoPolynomials(:, iZoneOfSource), r0, rmax, rhoTemp)
  call valueAtRadius(vpvPolynomials(:, iZoneOfSource), r0, rmax, vpvTemp)
  call valueAtRadius(vphPolynomials(:, iZoneOfSource), r0, rmax, vphTemp)
  call valueAtRadius(vsvPolynomials(:, iZoneOfSource), r0, rmax, vsvTemp)
  call valueAtRadius(vshPolynomials(:, iZoneOfSource), r0, rmax, vshTemp)
  call valueAtRadius(etaPolynomials(:, iZoneOfSource), r0, rmax, etaTemp)
  ecL0 = rhoTemp * vsvTemp * vsvTemp
  ecA0 = rhoTemp * vphTemp * vphTemp
  ecC0 = rhoTemp * vpvTemp * vpvTemp
  ecF0 = etaTemp * (ecA0 - 2.d0 * ecL0)

end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------





!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------









