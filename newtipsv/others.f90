
!----------------------------------------------------------------------------------------------------------------------------
! Reads the input parameters, which are given from standard input.
! The input is temporarily written in a work file, excluding the comments.
! Then, that temporary file is read in.
!----------------------------------------------------------------------------------------------------------------------------
subroutine readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegaI, imin, imax, nZone, rminOfZone, rmaxOfZone, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
  r0, eqlat, eqlon, mt, nReceiver, lat, lon, theta, phi, output)
!----------------------------------------------------------------------------------------------------------------------------
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


!----------------------------------------------------------------------------------------------------------------------------
! Decide if each zone is solid or fluid.
!----------------------------------------------------------------------------------------------------------------------------
subroutine judgeSolidOrFluid(nZone, vsPolynomials, phaseOfZone, nZoneSolid, nZoneFluid)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: vsPolynomials(4,nZone)  ! Polynomial functions of vs structure [km/s].
  integer, intent(out) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(out) :: nZoneSolid, nZoneFluid  ! Number of solid and fluid zones.
  integer :: iZone

  nZoneSolid = 0
  nZoneFluid = 0

  do iZone = 1, nZone
    if (vsPolynomials(1,iZone) == 0.d0 .and. vsPolynomials(2,iZone) == 0.d0 &
      .and. vsPolynomials(3,iZone) == 0.d0 .and. vsPolynomials(4,iZone) == 0.d0) then
      ! fluid
      nZoneFluid = nZoneFluid + 1
      phaseOfZone(iZone) = 2
    else
      ! solid
      nZoneSolid = nZoneSolid + 1
      phaseOfZone(iZone) = 1
    end if
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computes vertical wavenumber k_z at each zone. (See section 3.2 of Kawai et al. 2006.)
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeKz(nZone, rminOfZone, rmaxOfZone, phaseOfZone, vpPolynomials, vsPolynomials, rmax, imax, lmin, tlen, kzAtZone)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  real(8), intent(in) :: vpPolynomials(4,nZone), vsPolynomials(4,nZone)  ! Polynomial functions of vp and vs structure [km/s].
  real(8), intent(in) :: rmax  ! Maximum radius of region considered [km].
  integer, intent(in) :: imax  ! Index of maximum frequency.
  integer, intent(in) :: lmin  ! Smallest angular order l.
  real(8), intent(in) :: tlen  ! Time length [s].
  real(8), intent(out) :: kzAtZone(*)  ! Computed value of vertical wavenumber k_z at each zone [1/km].
  integer :: iZone
  real(8) :: v(4), vBottom, vTop, vmin, omega, kx, kz2

  do iZone = 1, nZone
    ! Use Vs in solid, Vp in fluid.
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


!----------------------------------------------------------------------------------------------------------------------------
! Deciding the distribution of grid points.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeGridRadii(maxNGrid, maxNGridSolid, maxNGridFluid, nZone, kzAtZone, rminOfZone, rmaxOfZone, phaseOfZone, &
  rmin, re, nGrid, nLayerInZone, gridRadii)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: maxNGrid  ! Maximum number of grid points.
  integer, intent(in) :: maxNGridSolid, maxNGridFluid  ! Maximum number of grid points in solid and fluid regions.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: kzAtZone(nZone)  ! Vertical wavenumber k_z at each zone [1/km].
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  real(8), intent(in) :: rmin  ! Minimum radius of region considered [km].
  real(8), intent(in) :: re  ! Desired relative error due to vertical gridding.
  integer, intent(out) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(out) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  real(8), intent(out) :: gridRadii(*)  ! Radius at each grid point [km].
  integer :: nGridSolid, nGridFluid  ! Number of grid points in solid and fluid regions.
  integer :: iZone, iGrid, i, nTemp
  real(8) :: rh

  ! Compute the distribution of grid points.
  iGrid = 1
  nGridSolid = 0
  nGridFluid = 0
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
    ! Accumulate number of layers in solid and fluid regions.
    if (phaseOfZone(iZone) == 1) then
      nGridSolid = nGridSolid + nLayerInZone(iZone) + 1
    else
      nGridFluid = nGridFluid + nLayerInZone(iZone) + 1
    end if
    ! Compute radius at each grid point [km].
    do i = 1, nLayerInZone(iZone)
      iGrid = iGrid + 1
      gridRadii(iGrid) = rminOfZone(iZone) + dble(i) * rh / dble(nLayerInZone(iZone))
    end do
  end do

  ! Register the total number of grid points.
  nGrid = iGrid
  if (nGrid > maxNGrid) stop 'The number of grid points is too large. (computeGridRadii)'
  if (nGridSolid > maxNGridSolid) stop 'The number of solid grid points is too large. (computeGridRadii)'
  if (nGridFluid > maxNGridFluid) stop 'The number of fluid grid points is too large. (computeGridRadii)'

  return
end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the first indices of each zone for vectors and matrices used later in the program.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeFirstIndices(nZone, nLayerInZone, phaseOfZone, oGridOfZone, oValueOfZone, oValueOfZoneSolid, &
  oPairOfZoneSolid, oPairOfZoneFluid, oElementOfZone, oColumnOfZone, nColumn)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(out) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(out) :: oValueOfZone(nZone)  ! Index of the first value in each zone.
  integer, intent(out) :: oValueOfZoneSolid(nZone)  ! Index of the first value in each zone, when counting only solid zones.
  integer, intent(out) :: oPairOfZoneSolid(nZone), oPairOfZoneFluid(nZone)
  !::::::::::::::::::::: Index of the first (iLayer, k', k)-pair in each zone, counted separately for solid and fluid zones.
  integer, intent(out) :: oElementOfZone(nZone)  ! Index of the first (iLayer, k'-gamma', k-gamma)-pair in each zone.
  integer, intent(out) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  integer, intent(out) :: nColumn  ! Total number of columns in the band matrix.
  integer :: iZone, iSolid, iFluid

  ! Set first index.
  oGridOfZone(1) = 1
  oValueOfZone(1) = 1
  oValueOfZoneSolid(1) = 1
  oPairOfZoneSolid(1) = 1
  oPairOfZoneFluid(1) = 1
  oElementOfZone(1) = 1
  oColumnOfZone(1) = 1

  ! Accumulate indices for each zone.
  iSolid = 0
  iFluid = 0
  do iZone = 1, nZone - 1
    oGridOfZone(iZone + 1) = oGridOfZone(iZone) + nLayerInZone(iZone)
    oValueOfZone(iZone + 1) = oValueOfZone(iZone) + nLayerInZone(iZone) + 1

    if (phaseOfZone(iZone) == 1) then
      ! solid
      iSolid = iSolid + 1
      oPairOfZoneSolid(iSolid + 1) = oPairOfZoneSolid(iSolid) + 4 * nLayerInZone(iZone)
      oValueOfZoneSolid(iSolid + 1) = oValueOfZoneSolid(iSolid) + nLayerInZone(iZone) + 1
      oElementOfZone(iZone + 1) = oElementOfZone(iZone) + 16 * nLayerInZone(iZone)
      if (phaseOfZone(iZone + 1) == 1) then
        oColumnOfZone(iZone + 1) = oColumnOfZone(iZone) + 2 * nLayerInZone(iZone)
      else
        oColumnOfZone(iZone + 1) = oColumnOfZone(iZone) + 2 * nLayerInZone(iZone) + 2
      end if

    else
      ! fluid
      iFluid = iFluid + 1
      oPairOfZoneFluid(iFluid + 1) = oPairOfZoneFluid(iFluid) + 4 * nLayerInZone(iZone)
      oElementOfZone(iZone + 1) = oElementOfZone(iZone) + 4 * nLayerInZone(iZone)
      if (phaseOfZone(iZone + 1) == 1) then
        oColumnOfZone(iZone + 1) = oColumnOfZone(iZone) + nLayerInZone(iZone) + 1
      else
        oColumnOfZone(iZone + 1) = oColumnOfZone(iZone) + nLayerInZone(iZone)
      end if

    end if
  end do

  ! Find end of last zone.
  if (phaseOfZone(nZone) == 1) then
    ! solid
    oColumnOfZone(nZone + 1) = oColumnOfZone(nZone) + 2 * nLayerInZone(nZone) + 2
  else
    ! fluid
    oColumnOfZone(nZone + 1) = oColumnOfZone(nZone) + nLayerInZone(nZone) + 1
  end if
  nColumn = oColumnOfZone(nZone + 1) - 1

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the source position.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeSourcePosition(nGrid, rmaxOfZone, phaseOfZone, gridRadii, r0, iZoneOfSource, iLayerOfSource, oPairOfSource)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nGrid  ! Total number of grid points.
  real(8), intent(in) :: rmaxOfZone(*)  ! Upper radius of each zone [km].
  integer, intent(in) :: phaseOfZone(*)  ! Phase of each zone (1: solid, 2: fluid).
  real(8), intent(in) :: gridRadii(*)  ! Radii of grid points [km].
  real(8), intent(inout) :: r0  ! Source radius [km]. Its value may be fixed in this subroutine.
  integer, intent(out) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(out) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(out) :: oPairOfSource  ! Index of the first (iLayer, k', k)-pair for the layer with the source.
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
  if (phaseOfZone(iZoneOfSource) /= 1) stop 'The source is in a fluid zone. (computeSourcePosition)'

  ! Find the layer that the source is in.
  iLayerOfSource = int(xLayerOfSource)  ! Note that int(x) rounds down the value x.
  ! Find the first index of (iLayer, k', k)-pair corresponding to the layer that the source is in.
  oPairOfSource = 4 * iLayerOfSource - 3

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing variable values at grid points.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeStructureValues(nZone, rmax, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, nLayerInZone, gridRadii, &
  nValue, valuedRadii, rhoValues, kappaValues, ecKxValues, ecKyValues, ecKzValues, ecLValues, ecNValues)
!----------------------------------------------------------------------------------------------------------------------------
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


!----------------------------------------------------------------------------------------------------------------------------
! Computing variable values near the source.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeSourceStructureValues(iZoneOfSource, r0, rmax, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, ecC0, ecF0, ecL0)
!----------------------------------------------------------------------------------------------------------------------------
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


!----------------------------------------------------------------------------------------------------------------------------
! Computing the inverse of density and elastic constant.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeReciprocals(nValue, rhoValues, kappaValues, rhoReciprocals, kappaReciprocals)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nValue
  real(8), intent(in) :: rhoValues(nValue), kappaValues(nValue)
  !::::::::::::::::::::::: Values of rho [g/cm^3] and kappa [10^10 dyn/cm^2 = GPa] at each point (with 2 values at boundaries).
  real(8), intent(out) :: rhoReciprocals(nValue), kappaReciprocals(nValue)
  !:::::::::: Values of 1/rho [cm^3/g] and 1/kappa [(10^10 dyn/cm^2)^(-1) = 1/GPa] at each point (with 2 values at boundaries).

  ! Compute inverses using array operations
  rhoReciprocals(:) = 1.d0 / rhoValues(:)
  kappaReciprocals(:) = 1.d0 / kappaValues(:)

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computes the coefficients to multiply to elastic moduli for anelastic attenuation.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeCoef(nZone, omega, qmuOfZone, qkappaOfZone, coefQmu, coefQkappa, coefQfluid)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: omega  ! Angular frequency [1/s].
  real(8), intent(in) :: qmuOfZone(nZone)  ! Qmu of each zone.
  real(8), intent(in) :: qkappaOfZone(nZone)  ! Qkappa of each zone.
  complex(8), intent(out) :: coefQmu(nZone), coefQkappa(nZone), coefQfluid(nZone)
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  real(8) :: aa, bb
  integer :: iZone

  do iZone = 1, nZone

    ! Compute coefficient for Qmu.
    if (qmuOfZone(iZone) <= 0.d0) then
      coefQmu(iZone) = dcmplx(1.d0)
    else
      if (omega == 0.d0) then
        aa = 1.d0
      else
        aa = 1.d0 + log(omega / (2.d0 * pi)) / (pi * qmuOfZone(iZone))
      end if
      bb = 1.d0 / (2.d0 * qmuOfZone(iZone))
      coefQmu(iZone) = dcmplx(aa, bb) ** 2
    end if

    ! Compute coefficient for Qkappa.
    if (qkappaOfZone(iZone) <= 0.d0) then
      coefQkappa(iZone) = dcmplx(1.d0)
    else
      if (omega == 0.d0) then
        aa = 1.d0
      else
        aa = 1.d0 + dlog(omega / (2.d0 * pi)) / (pi * qkappaOfZone(iZone))
      end if
      bb = 1.d0 / (2.d0 * qkappaOfZone(iZone))
      coefQkappa(iZone) = dcmplx(aa, bb) ** 2
    end if

    ! Compute coefficient for Q in fluid.
    coefQfluid(iZone) = dcmplx(1.d0) / coefQkappa(iZone)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Evaluate the cut-off depth based on the relative amplitude at each depth.
! (See the end of section 3.2 of Kawai et al. 2006.)
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeCutoffColumn(nZone, phaseOfZone, nGrid, oGridOfZone, nColumn, oColumnOfZone, &
  amplitudeAtColumn, ratc, cutoffColumn)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(in) :: nGrid  ! Total number of grid points.
  integer, intent(in) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(in) :: nColumn  ! Total number of columns in the band matrix.
  integer, intent(in) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  real(8), intent(in) :: amplitudeAtColumn(nColumn)  ! Estimate of the amplitude at each c vector component [km].
  real(8), intent(in) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  integer, intent(out) :: cutoffColumn  ! Index of column at cut-off depth.
  real(8) :: amplitudeAtGrid(nGrid)  ! Estimate of the amplitude at each grid point [km].
  real(8) :: amplitudeThreshold  ! Threshold to decide at which grid to cut off.
  integer :: cutoffGrid  ! Index of grid at cut-off depth.
  integer :: iZone, iGrid, iColumn, iZoneOfCutoff

  ! Transplant amplitudeAtColumn(:) into amplitudeAtGrid(:), picking only the vertical component for solid zones.
  iZone = 1
  iGrid = 1
  do iColumn = 1, nColumn

    ! When next zone starts.
    if (iColumn == oColumnOfZone(iZone + 1)) then
      ! At solid-fluid boundary, the same grid is referenced in both zones, so ignore one of them.
      if (phaseOfZone(iZone + 1) /= phaseOfZone(iZone)) iGrid = iGrid - 1
      ! Increment.
      iZone = iZone + 1
    end if

    if (phaseOfZone(iZone) == 1) then
      ! solid
      ! Use only the vertical component.
      if (mod((iColumn - oColumnOfZone(iZone)), 2) == 1) then
        amplitudeAtGrid(iGrid) = amplitudeAtColumn(iColumn)
        iGrid = iGrid + 1
      end if
    else
      ! fluid
      amplitudeAtGrid(iGrid) = amplitudeAtColumn(iColumn)
      iGrid = iGrid + 1
    end if
  end do

  if (iGrid - 1 /= nGrid) stop 'Computation error: nGrid does not match. (computeCutoffDepth)'

  ! Set the threshold amplitude as ratc * the maximum amplitude.
  amplitudeThreshold = maxval(amplitudeAtGrid(1:nGrid)) * ratc

  ! If maxamp is zero, set cutoffColumn to 1 and return.
  if (amplitudeThreshold == 0.d0) then
    cutoffColumn = 1
    return
  end if

  ! Identify the first grid with amplitude greater than threshold value.
  do iGrid = 1, nGrid
    if (amplitudeAtGrid(iGrid) > amplitudeThreshold) then
      cutoffGrid = iGrid
      exit
    end if
  end do

  ! Find the zone that cutoff grid is in.
  do iZone = 1, nZone - 1
    if (oGridOfZone(iZone + 1) > cutoffGrid) then
      iZoneOfCutoff = iZone
      exit
    end if
  end do

  ! Find the column corresponding to the cutoff grid.
  if (phaseOfZone(iZoneOfCutoff) == 1) then
    ! solid
    cutoffColumn = oColumnOfZone(iZoneOfCutoff) + 2 * (cutoffGrid - oGridOfZone(iZoneOfCutoff))
  else if (phaseOfZone(iZoneOfCutoff) == 2) then
    ! fluid
    cutoffColumn = oColumnOfZone(iZoneOfCutoff) + cutoffGrid - oGridOfZone(iZoneOfCutoff)
  end if

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Checks whether the expansion coefficient amplitude has decayed enough to stop computing subsequent l's for this frequency.
! When the amplitude ratio is smaller than the threshold and l has surpassed the accuracy threshold, a counter is incremented.
! The maximum amplitude encountered so far is also recorded using this subroutine.
!----------------------------------------------------------------------------------------------------------------------------
subroutine checkAmplitudeDecay(c0, l, lsuf, ratl, recordAmplitude, decayCounter)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  complex(8), intent(in) :: c0(2)  ! Expansion coefficients at topmost grid [km].
  integer, intent(in) :: l  ! Angular order.
  integer, intent(in) :: lsuf  ! Accuracy threshold of angular order.
  real(8), intent(in) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8), intent(inout) :: recordAmplitude  ! Maximum amplitude encountered, updated if the current amplitude is larger [km].
  integer, intent(inout) :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  real(8) :: amp, ampratio

  ! Calculate the amplitude of expansion coefficients.
  amp = sqrt(abs(c0(1)) ** 2 + abs(c0(2)) ** 2)

  ! Update maximum amplitude if current amplitude is greater.
  if (amp > recordAmplitude) recordAmplitude = amp

  ! Calculate the amplitude ratio.
  ampratio = 0.d0
  if (amp /= 0.d0 .and. recordAmplitude /= 0.d0) then
    ampratio = amp / recordAmplitude
  endif

  ! Increment the counter if amplitude ratio is smaller than its threshold and l has surpassed the accuracy threshold.
  !  Reset the counter otherwise.
  if (l > lsuf .and. ampratio < ratl) then
    decayCounter = decayCounter + 1
  else
    decayCounter = 0
  endif

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Accumulates the value of 'u' at a certain receiver on planet surface for a certain (l, m)-pair (= for a certain trial function).
!  (See eq. 1 of Kawai et al. 2006.)
! The trial function is specified by (k=k_max (at surface of planet), l, m, 1:2 (the S^1 and S^2 spherical harmonics)).
!  (See eqs. 12 & 13 of Kawai et al. 2006.)
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeU(c0, largeL2, harmonicsValues, u)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  complex(8), intent(in) :: c0(2)  ! Expansion coefficent corresponding to this trial function [km] (k=k_max, l, m, 3).
  real(8), intent(in) :: largeL2  ! L^2 = l(l+1).
  complex(8), intent(in) :: harmonicsValues(3)  ! Vector harmonics term. The coefficient 1/largeL is not multiplied yet.
  complex(8), intent(inout) :: u(3)  ! Displacement velocity - the unit is [km] in the frequency domain.
  complex(8) :: largeLc

  ! Compute L.
  largeLc = dcmplx(sqrt(largeL2))

  ! Accumulate value of u. (See eq. 1 of Kawai et al. 2006.)
  ! The coefficient 1/largeL is not yet multiplied to the vector harmonics term, so is multiplied here.
  !  (See eq. 12 of Kawai et al. 2006.)
  u(1) = u(1) + c0(1) * harmonicsValues(1)
  u(2) = u(2) + c0(2) * harmonicsValues(2) / largeLc
  u(3) = u(3) + c0(2) * harmonicsValues(3) / largeLc

end subroutine


