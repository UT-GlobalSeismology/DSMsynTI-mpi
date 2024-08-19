
!----------------------------------------------------------------------------------------------------------------------------
! Reads the input parameters, which are given from standard input.
! The input is temporarily written in a work file, excluding the comments.
! Then, that temporary file is read in.
!----------------------------------------------------------------------------------------------------------------------------
subroutine readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegaI, imin, imax, &
  nZone, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, qmuOfZone, &
  r0, eqlat, eqlon, mt, nReceiver, lat, lon, theta, phi, output)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=80), parameter :: tmpfile = 'workSH.txt'

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
  real(8), intent(out) :: rhoPolynomials(4,*), vsvPolynomials(4,*), vshPolynomials(4,*)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vsv, and vsh [km/s] structure.
  real(8), intent(out) :: qmuOfZone(*)  ! Qmu of each zone.
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
    read(11,*) vsvPolynomials(1,i), vsvPolynomials(2,i), vsvPolynomials(3,i), vsvPolynomials(4,i)
    read(11,*) vshPolynomials(1,i), vshPolynomials(2,i), vshPolynomials(3,i), vshPolynomials(4,i), qmuOfZone(i)
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

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computes vertical wavenumber k_z at each zone. (See section 3.2 of Kawai et al. 2006.)
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeKz(nZone, rminOfZone, rmaxOfZone, vsPolynomials, rmax, imax, lmin, tlen, kzAtZone)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  real(8), intent(in) :: vsPolynomials(4,nZone)  ! Polynomial functions of vs structure [km/s].
  real(8), intent(in) :: rmax  ! Maximum radius of region considered [km].
  integer, intent(in) :: imax  ! Index of maximum frequency.
  integer, intent(in) :: lmin  ! Smallest angular order l.
  real(8), intent(in) :: tlen  ! Time length [s].
  real(8), intent(out) :: kzAtZone(*)  ! Computed value of vertical wavenumber k_z at each zone [1/km].
  integer :: iZone
  real(8) :: v(4), vBottom, vTop, vmin, omega, kx, kz2

  do iZone = 1, nZone
    v(:) = vsPolynomials(:, iZone)
    ! Compute velocity [km/s] at bottom and top of zone.
    call valueAtRadius(v(:), rminOfZone(iZone), rmax, vBottom)
    call valueAtRadius(v(:), rmaxOfZone(iZone), rmax, vTop)
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

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Deciding the distribution of grid points.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeGridRadii(maxNGrid, nZone, kzAtZone, rminOfZone, rmaxOfZone, rmin, re, nGrid, nLayerInZone, gridRadii)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: maxNGrid  ! Maximum number of grid points.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: kzAtZone(nZone)  ! Vertical wavenumber k_z at each zone [1/km].
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  real(8), intent(in) :: rmin  ! Minimum radius of region considered [km].
  real(8), intent(in) :: re  ! Desired relative error due to vertical gridding.
  integer, intent(out) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(out) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  real(8), intent(out) :: gridRadii(*)  ! Radius at each grid point [km].
  integer :: iZone, iGrid, i, nTemp
  real(8) :: rh

  ! Compute the distribution of grid points.
  iGrid = 1
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
    ! Compute radius at each grid point [km].
    do i = 1, nLayerInZone(iZone)
      iGrid = iGrid + 1
      gridRadii(iGrid) = rminOfZone(iZone) + dble(i) * rh / dble(nLayerInZone(iZone))
    end do
  end do

  ! Register the total number of grid points.
  nGrid = iGrid
  if (nGrid > maxNGrid) stop 'The number of grid points is too large. (computeGridRadii)'

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the indices of the first grid point and the first (iLayer, k', k)-pair in each zone.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeFirstIndices(nZone, nLayerInZone, oGridOfZone, oValueOfZone, oPairOfZone)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  integer, intent(out) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(out) :: oValueOfZone(nZone)  ! Index of the first value in each zone.
  integer, intent(out) :: oPairOfZone(nZone)  ! Index of the first (iLayer, k', k)-pair in each zone.
  integer :: iZone

  oGridOfZone(1) = 1
  oValueOfZone(1) = 1
  oPairOfZone(1) = 1
  do iZone = 1, nZone - 1
    oGridOfZone(iZone + 1) = oGridOfZone(iZone) + nLayerInZone(iZone)
    oValueOfZone(iZone + 1) = oValueOfZone(iZone) + nLayerInZone(iZone) + 1
    oPairOfZone(iZone + 1) = oPairOfZone(iZone) + 4 * nLayerInZone(iZone)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the source position.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeSourcePosition(nGrid, rmaxOfZone, gridRadii, r0, iZoneOfSource, iLayerOfSource, oPairOfSource)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nGrid  ! Total number of grid points.
  real(8), intent(in) :: rmaxOfZone(*)  ! Upper radius of each zone [km].
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

  ! Find the layer that the source is in.
  iLayerOfSource = int(xLayerOfSource)  ! Note that int(x) rounds down the value x.
  ! Find the first index of (iLayer, k', k)-pair corresponding to the layer that the source is in.
  oPairOfSource = 4 * iLayerOfSource - 3

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Compute grids for the source based on input parameters.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeSourceGrid(gridRadii, r0, iLayerOfSource, gridRadiiForSource)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: gridRadii(*)  ! Radii of grid points [km].
  real(8), intent(in) :: r0  ! Input source radius [km].
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  real(8), intent(out) :: gridRadiiForSource(3)  ! Radii to use for source-related computations [km].

  ! Assign grid radii.
  gridRadiiForSource(1) = gridRadii(iLayerOfSource)  ! radius of bottom of layer
  gridRadiiForSource(2) = r0 ! radius of source
  gridRadiiForSource(3) = gridRadii(iLayerOfSource + 1) ! radius of top of layer

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing variable values at grid points.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeStructureValues(nZone, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, nLayerInZone, gridRadii, &
  nValue, valuedRadii, rhoValues, ecLValues, ecNValues)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmax  ! Maximum radius of region considered [km].
  real(8), intent(in) :: rhoPolynomials(4,nZone), vsvPolynomials(4,nZone), vshPolynomials(4,nZone)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vsv, and vsh [km/s] structure.
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  real(8), intent(in) :: gridRadii(*)  ! Radii of grid points [km].
  integer, intent(out) :: nValue  ! Total number of values for each variable.
  real(8), intent(out) :: valuedRadii(*)  ! Radii corresponding to each variable value [km].
  real(8), intent(out) :: rhoValues(*), ecLValues(*), ecNValues(*)  ! Values of rho [g/cm^3], L, and N [10^10 dyn/cm^2 = GPa]
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: at each point (with 2 values at boundaries).
  real(8) :: rhoTemp, vsvTemp, vshTemp
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
      call valueAtRadius(vsvPolynomials(:, iZone), valuedRadii(iValue), rmax, vsvTemp)
      call valueAtRadius(vshPolynomials(:, iZone), valuedRadii(iValue), rmax, vshTemp)
      rhoValues(iValue) = rhoTemp
      ecLValues(iValue) = rhoTemp * vsvTemp * vsvTemp
      ecNValues(iValue) = rhoTemp * vshTemp * vshTemp
    end do

    iGrid = iGrid - 1
  end do

  nValue = iValue

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing variable values near the source.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, gridRadiiForSource, &
  rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, ecL0)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  real(8), intent(in) :: rmax  ! Maximum radius of region considered [km].
  real(8), intent(in) :: rhoPolynomials(4,*), vsvPolynomials(4,*), vshPolynomials(4,*)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vsv, and vsh [km/s] structure.
  real(8), intent(in) :: gridRadiiForSource(3)  ! Radii to use for source-related computations [km].
  real(8), intent(out) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::  Values of rho [g/cm^3], L, and N [10^10 dyn/cm^2 = GPa] at each point.
  real(8), intent(out) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  real(8) :: rhoTemp, vsvTemp, vshTemp
  integer :: i

  ! Compute variable values at grid points.
  do i = 1, 3
    ! Evaluate the density and elastic constants at this point.
    call valueAtRadius(rhoPolynomials(:, iZoneOfSource), gridRadiiForSource(i), rmax, rhoTemp)
    call valueAtRadius(vsvPolynomials(:, iZoneOfSource), gridRadiiForSource(i), rmax, vsvTemp)
    call valueAtRadius(vshPolynomials(:, iZoneOfSource), gridRadiiForSource(i), rmax, vshTemp)
    rhoValuesForSource(i) = rhoTemp
    ecLValuesForSource(i) = rhoTemp * vsvTemp * vsvTemp
    ecNValuesForSource(i) = rhoTemp * vshTemp * vshTemp
  end do

  ! Record the value of L at source position.
  ecL0 = ecLValuesForSource(2)

end subroutine


!!TODO understand
!----------------------------------------------------------------------------------------------------------------------------
! Computes the coefficients to multiply to elastic moduli for anelastic attenuation.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeCoef(nZone, omega, qmuOfZone, coefQmu)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: omega  ! Angular frequency [1/s].
  real(8), intent(in) :: qmuOfZone(nZone)  ! Qmu of each zone.
  complex(8), intent(out) :: coefQmu(nZone)  ! Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
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
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Evaluate the cut-off depth based on the relative amplitude at each depth.
! (See the end of section 3.2 of Kawai et al. 2006.)
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeCutoffGrid(nGrid, amplitudeAtGrid, ratc, cutoffGrid)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nGrid  ! Total number of grid points.
  real(8), intent(in) :: amplitudeAtGrid(nGrid)  ! Estimate of the amplitude at each grid point [km].
  real(8), intent(in) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  integer, intent(out) :: cutoffGrid  ! Index of grid at cut-off depth.
  real(8) :: amplitudeThreshold  ! Threshold to decide at which grid to cut off.
  integer :: iGrid

  ! Set the threshold amplitude as ratc * the maximum amplitude.
  amplitudeThreshold = maxval(amplitudeAtGrid(1:nGrid)) * ratc

  ! If maxamp is zero, set cutoffGrid to 1 and return.
  if (amplitudeThreshold == 0.d0) then
    cutoffGrid = 1
    return
  end if

  ! Identify the first grid with amplitude greater than threshold value.
  do iGrid = 1, nGrid
    if (amplitudeAtGrid(iGrid) > amplitudeThreshold) then
      cutoffGrid = iGrid
      exit
    end if
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Checks whether the expansion coefficient amplitude has decayed enough to stop computing subsequent l's for this frequency.
! When the amplitude ratio is smaller than the threshold and l has surpassed the accuracy threshold, a counter is incremented.
! The maximum amplitude encountered so far is also recorded using this subroutine.
!----------------------------------------------------------------------------------------------------------------------------
subroutine checkAmplitudeDecay(c0, l, lsuf, ratl, recordAmplitude, decayCounter)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  complex(8), intent(in) :: c0  ! Expansion coefficient at topmost grid [km].
  integer, intent(in) :: l  ! Angular order.
  integer, intent(in) :: lsuf  ! Accuracy threshold of angular order.
  real(8), intent(in) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8), intent(inout) :: recordAmplitude  ! Maximum amplitude encountered, updated if the current amplitude is larger [km].
  integer, intent(inout) :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  real(8) :: amp, ampratio

  ! Calculate the amplitude of expansion coefficient.
  amp = abs(c0)

  ! Update maximum amplitude if current amplitude is greater.
  if (amp > recordAmplitude) recordAmplitude = amp

  ! Calculate the amplitude ratio.
  ampratio = 0.d0
  if (amp /= 0.d0 .and. recordAmplitude /= 0.d0) then
    ampratio = amp / recordAmplitude
  end if

  ! Increment the counter if amplitude ratio is smaller than its threshold and l has surpassed the accuracy threshold.
  !  Reset the counter otherwise.
  if (l > lsuf .and. ampratio < ratl) then
    decayCounter = decayCounter + 1
  else
    decayCounter = 0
  end if

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Accumulates the value of 'u' at a certain receiver on planet surface for a certain (l, m)-pair (= for a certain trial function).
!  (See eq. 1 of Kawai et al. 2006.)
! The trial function is specified by (k=k_max (at surface of planet), l, m, 3 (the T spherical harmonic)).
!  (See eqs. 12 & 13 of Kawai et al. 2006.)
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeU(c0, largeL2, harmonicsValues, u)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  complex(8), intent(in) :: c0  ! Expansion coefficent corresponding to this trial function [km] (k=k_max, l, m, 3).
  real(8), intent(in) :: largeL2  ! L^2 = l(l+1).
  complex(8), intent(in) :: harmonicsValues(3)  ! Vector harmonics term. The coefficient 1/largeL is not multiplied yet.
  complex(8), intent(inout) :: u(3)  ! Displacement velocity - the unit is [km] in the frequency domain.
  complex(8) :: largeLc

  ! Compute L.
  largeLc = dcmplx(sqrt(largeL2))

  ! Accumulate value of u. (See eq. 1 of Kawai et al. 2006.)
  ! The coefficient 1/largeL is not yet multiplied to the vector harmonics term, so is multiplied here.
  !  (See eq. 12 of Kawai et al. 2006.)
  u(1) = dcmplx(0.d0, 0.d0)
  u(2) = u(2) + c0 * harmonicsValues(3) / largeLc
  u(3) = u(3) - c0 * harmonicsValues(2) / largeLc

end subroutine
