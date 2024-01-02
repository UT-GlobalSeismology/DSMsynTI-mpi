!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ************** tipsv ****************
!  Computation of PSV synthetic seismograms
!  in transversely isotropic media for anisotropic PREM
!  using modified DSM operators & modified source representation.
!  Synthetics for shallow events can be computed.
!
!  Main historical authors: K.Kawai, N.Takeuchi, R.J.Geller
!  (C) 2002 - 2023  University of Tokyo
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program tipsv
  use FileIO
  implicit none

  !----------------------------<<constants>>----------------------------
  real(8), parameter :: pi = 3.1415926535897932d0
  integer, parameter :: maxNGrid = 88300  ! Maximum number of grid points.
  integer, parameter :: maxNGridSolid = 48840  ! Maximum number of grid points in solid region.
  integer, parameter :: maxNGridFluid = 32040  ! Maximum number of grid points in fluid region.
  integer, parameter :: maxNZone = 15  ! Maximum number of zones.
  integer, parameter :: maxNReceiver = 1500  ! Maximum number of receivers.
  integer, parameter :: maxL = 80000  ! Maximum of angular order to loop for.
  real(8), parameter :: lmaxdivf = 2.d4  !    !!TODO where did this value come from?
  real(8), parameter :: shallowdepth = 100.d0  ! Threshold to consider evanescent regime for shallow events [km].
  integer, parameter :: spcFormat = 1  ! Format of output spc file (0:binary, 1:ascii).
  integer, parameter :: ilog = 1

  !----------------------------<<variables>>----------------------------
  ! Variables for the structure
  integer :: nZone  ! Number of zones.
  integer :: nZoneSolid, nZoneFluid  ! Number of solid and fluid zones.
  integer :: phaseOfZone(maxNZone)  ! Phase of each zone (1: solid, 2: fluid).
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled [km].
  real(8) :: rminOfZone(maxNZone), rmaxOfZone(maxNZone)  ! Minimum and maximum radii of each zone [km].
  real(8) :: rhoPolynomials(4, maxNZone), vpvPolynomials(4, maxNZone), vphPolynomials(4, maxNZone)
  real(8) :: vsvPolynomials(4, maxNZone), vshPolynomials(4, maxNZone), etaPolynomials(4, maxNZone)
  !::::: Polynomial functions (coefficients of cubic function) of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta in each zone.
  real(8) :: qmuOfZone(maxNZone), qkappaOfZone(maxNZone)  ! Qmu and Qkappa of each zone.
  integer :: i, iSolid, iFluid

  ! Variables for the source
  real(8) :: r0, eqlat, eqlon, mt(3, 3)  ! Depth [km], coordinates [deg], and moment tensor [10^25 dyn cm] of source.
  real(8) :: ecC0, ecF0, ecL0  ! Elastic moduli C, F, L at source position [10^10 dyn/cm^2 = GPa].

  ! Variables for the receivers
  integer :: nReceiver  ! Number of receivers.
  real(8) :: lat(maxNReceiver), lon(maxNReceiver)  ! Coordinates [deg] of receivers.
  real(8) :: theta(maxNReceiver), phi(maxNReceiver)  ! Colatitude and longitude of receivers with event at north pole [rad].
  integer :: ir

  ! Variables for the periodic range
  real(8) :: tlen  ! Time length [s].
  integer :: np  ! Number of points in frequency domain.
  real(8) :: omega  ! Angular frequency (real part) [1/s].
  real(8) :: omegaI  ! Imaginary part of angular frequency for artificial damping [1/s].
  !:::::::::::::::::::::::::::::::::::::::: (See section 5.1 of Geller & Ohminato 1994.)
  integer :: imin, imax  ! Index of minimum and maximum frequency.
  integer :: iFreq, iCount
  integer :: ltmp(2), imaxFixed

  ! Variables for grid spacing and cut-off
  real(8) :: kzAtZone(maxNZone)  ! Vertical wavenumber k_z at each zone [1/km]. (See section 3.2 of Kawai et al. 2006.)
  real(8) :: re  ! Desired relative error due to vertical gridding. (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
  real(8) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8) :: amplitudeAtGrid(maxNGrid)  ! Estimate of the amplitude at each grid point [km], used for vertical grid cut-off.
  integer :: cutoffGrid  ! Index of grid at cut-off depth.
  integer :: lsuf  ! Accuracy threshold of angular order. (Corresponds to l_d; see eq. 29 of Kawai et al. 2006.)
  real(8) :: recordAmplitude    ! Maximum amplitude encountered [km], used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: llog

  ! Variables for the vertical grid
  integer :: nGrid  ! Total number of grid points.
  real(8) :: gridRadii(maxNGrid)  ! Radii of each grid point [km].
  integer :: nLayerInZone(maxNZone)  ! Number of layers in each zone.
  real(8) :: gridRadiiForSource(3)  ! Radii to use for source-related computations [km].
  integer :: iZoneOfSource  ! Which zone the source is in.
  integer :: iLayerOfSource  ! Index of layer that the source is in.

  ! Variables for the values
  integer :: nValue  ! Total number of values.
  real(8) :: valuedRadii(maxNGrid + maxNZone - 1)  ! Radii corresponding to each variable value [km].
  integer :: oValueOfZone(maxNZone)  ! Index of the first value in each zone.
  integer :: oValueOfZoneSolid(maxNZone)  ! Index of the first value in each zone, when counting only solid zones.
  real(8) :: rhoValues(maxNGrid + maxNZone - 1)  ! Rho at each grid point (with 2 values at boundaries) [g/cm^3].
  real(8) :: kappaValues(maxNGrid + maxNZone - 1)  !
  real(8) :: ecKxValues(maxNGrid + maxNZone - 1)  ! 3*Kx=3A-4N
  real(8) :: ecKyValues(maxNGrid + maxNZone - 1)  ! 3*Ky=3F+2N
  real(8) :: ecKzValues(maxNGrid + maxNZone - 1)  ! 3*Kz=2F+C
  real(8) :: ecLValues(maxNGrid + maxNZone - 1)  ! L at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: ecNValues(maxNGrid + maxNZone - 1)  ! N at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: rhoReciprocals(maxNGrid + maxNZone - 1)  !
  real(8) :: kappaReciprocals(maxNGrid + maxNZone - 1)  !
!  real(8) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)  ! Rho, L, and N at each source-related grid.
  complex(8) :: coefQmu(maxNZone), coefQkappa(maxNZone), coefQfluid(maxNZone)
  !::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  integer :: oV, oVS

  ! Variables for the trial function
  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: largeL2, largeL  ! L^2 = l(l+1).
  real(8) :: plm(3, 0:3, maxNReceiver)  ! Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  !::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8) :: harmonicsValues(3, -2:2, maxNReceiver)  ! Values of vector harmonics terms at each receiver, computed for each l.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Arguments: term (1:3), m (-2:2), iReceiver.

  ! Variables for the matrix elements
  real(8) :: t(4 * maxNGridSolid - 4)
  real(8) :: h1x(4 * maxNGridSolid - 4), h2L(4 * maxNGridSolid - 4), h2N(4 * maxNGridSolid - 4)
  real(8) :: hUn3y(4 * maxNGridSolid - 4), hResid3y(4 * maxNGridSolid - 4), h3Mod2y(-2:1, maxNGridSolid)
  real(8) :: hUn4L(4 * maxNGridSolid - 4), hResid4L(4 * maxNGridSolid - 4), h4Mod1L(-1:2, maxNGridSolid)
  real(8) :: hUn4N(4 * maxNGridSolid - 4), hResid4N(4 * maxNGridSolid - 4), h4Mod2N(-2:1, maxNGridSolid)
  real(8) :: hUn5y(4 * maxNGridSolid - 4), hResid5y(4 * maxNGridSolid - 4), h5Mod1y(-1:2, maxNGridSolid)
  real(8) :: hUn6L(4 * maxNGridSolid - 4), hResid6L(4 * maxNGridSolid - 4), h6Mod2L(-2:1, maxNGridSolid)
  real(8) :: hUn6N(4 * maxNGridSolid - 4), hResid6N(4 * maxNGridSolid - 4), h6Mod1N(-1:2, maxNGridSolid)
  real(8) :: h7y(4 * maxNGridSolid - 4), h7z(4 * maxNGridSolid - 4), h8L(4 * maxNGridSolid - 4), h8N(4 * maxNGridSolid - 4)
  real(8) :: p1(4 * maxNGridFluid - 4), p2(4 * maxNGridFluid - 4), p3(4 * maxNGridFluid - 4)
  real(8) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  integer :: oRowOfZoneSolid(maxNZone), oRowOfZoneFluid(maxNZone)
  !:: Index of the first row in the vector of (iLayer, k', k)-pairs in each zone. Vectors are separate for solid and fluid zones.
  integer :: oRowOfSource  ! Index of the first row in the vector of (iLayer, k', k)-pairs for the layer with the source.
  complex(8) :: a0(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a1(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a2(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: aaParts(4), aSourceParts(8), aSource(2, 3)
  complex(8) :: g_or_c(maxNGrid)  ! This holds either vector g [10^15 N] or c [km], depending on where in the code it is. CAUTION!!
  complex(8) :: u(3, maxNReceiver)  ! Displacement velocity - the unit is [km] in the frequency domain,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: but when converted to the time domain, the unit becomes [km/s].
  integer :: oElementOfZone(maxNZone)  ! Index of the first (iLayer, k'-gamma', k-gamma)-pair in each zone.
  integer :: oColumnOfZone(maxNZone + 1)  ! Index of the first column in the band matrix for each zone.
  integer :: nColumn  ! Total number of columns in the band matrix.
  integer :: oR, oElement, oColumn

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Other variables
  real(8) :: work(4 * maxNGrid - 4)  ! Working array for matrix computations.
  complex(8) :: cwork(16 * maxNGridSolid - 16 + 4 * maxNGridFluid - 4)  ! Working array for matrix computations.


  ! Efficiency improvement variables
  !     When the values to be output use memory over outputMemory MB,
  !     they are written in output files. The interval is outputInterval.
  real(8) :: outputMemory = 10  ! Approximate amount of memory to write at a time. [MB]
  real(8) :: memoryPerOmega  ! The amount of memory used for one omega step. [MB]
  integer :: outputInterval  ! Interval of omega that data should be written out.
  integer :: outputCounter  ! Counter to keep track of how many omegas have been computed after previous output.
  integer :: iOut
  integer, allocatable :: outputi(:)
  complex(8), allocatable :: outputu(:,:,:)


  ! ************************** Inputting parameters **************************
  ! --- read parameters ---
  call readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegaI, imin, imax, nZone, rminOfZone, rmaxOfZone, &
    rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
    r0, eqlat, eqlon, mt, nReceiver, lat, lon, theta, phi, output)

  ! --- computing the required parameters ---
  rmin = rminOfZone(1)
  rmax = rmaxOfZone(nZone)
  if (r0 < rmin .or. r0 > rmax) stop 'The source position is improper.'

  call judgeSolidOrFluid(nZone, vsvPolynomials, phaseOfZone, nZoneSolid, nZoneFluid)

  ! Find the amount of memory that is written in 1 omega step.
  !  For each omega and receiver, 3 complex numbers (16 B each) are written. 1 B = 0.000001 MB.
  memoryPerOmega = 3 * 16 * nReceiver * 0.000001
  ! Find how many omegas can be written within outputMemory.
  outputInterval = int(outputMemory / memoryPerOmega)
  ! Allocate arrays to store output.
  allocate(outputi(outputInterval))
  allocate(outputu(3, nReceiver, outputInterval))

  if (imin == 0) imin = 1
  ! Decide which omega to use when deciding grid spacing. Usually, this is just the upper limit of omega range.
  imaxFixed = imax


  ! ************************** Files handling **************************
  do ir = 1, nReceiver
    call openSPCFile(output(ir), 11, spcFormat, 0)
    call writeSPCFile(11, spcFormat, tlen)
    call writeSPCFile(11, spcFormat, np, 1, 3)
    call writeSPCFile(11, spcFormat, omegaI, lat(ir), lon(ir))
    call writeSPCFile(11, spcFormat, eqlat, eqlon, r0)
    call closeSPCFile(11)
  end do

  if (ilog == 1) then
    open(unit = 11, file = 'llog.log', status = 'unknown')
    write(11, *) 'iFreq, llog, nGrid-1'
    close(11)
  end if


  ! ########################## Option for shallow events ##########################
  ! Here, we find the maximum angular order needed for our frequency range. (See fig. 7 of Kawai et al. 2006.)










  ! ########################## Main computation ##########################

  ! ******************* Computing the matrix elements *******************

  !<<< operation from here >>>>

  ! ------------------- Computing parameters -------------------
  ! Design the number and position of grid points.
  call computeKz(nZone, rminOfZone(:), rmaxOfZone(:), phaseOfZone(:), vpvPolynomials(:,:), vsvPolynomials(:,:), &
    rmax, imaxFixed, 1, tlen, kzAtZone(:))
  call computeGridRadii(maxNGrid, maxNGridSolid, maxNGridFluid, nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), phaseOfZone(:), &
    rmin, re, nGrid, nLayerInZone(:), gridRadii(:))

  ! Compute the first indices in each zone.
  call computeFirstIndices(nZone, nLayerInZone(:), phaseOfZone(:), oValueOfZone(:), oValueOfZoneSolid(:), &
    oRowOfZoneSolid(:), oRowOfZoneFluid(:), oElementOfZone(:), oColumnOfZone(:), nColumn)

  ! Compute the source position.
  call computeSourcePosition(nGrid, rmaxOfZone(:), phaseOfZone(:), gridRadii(:), r0, iZoneOfSource, iLayerOfSource, oRowOfSource)


  ! ------------------- Computing the matrix elements -------------------
  ! Compute variable values at grid points.
  call computeStructureValues(nZone, rmax, rhoPolynomials(:,:), vpvPolynomials(:,:), vphPolynomials(:,:), &
    vsvPolynomials(:,:), vshPolynomials(:,:), etaPolynomials(:,:), nLayerInZone(:), gridRadii(:), &
    nValue, valuedRadii(:), rhoValues(:), kappaValues(:), ecKxValues(:), ecKyValues(:), ecKzValues(:), ecLValues(:), ecNValues(:))
  call computeSourceStructureValues(iZoneOfSource, r0, rmax, rhoPolynomials(:,:), vpvPolynomials(:,:), vphPolynomials(:,:), &
    vsvPolynomials(:,:), vshPolynomials(:,:), etaPolynomials(:,:), ecC0, ecF0, ecL0)

  ! Compute 1/rho and 1/kappa.
  call computeReciprocals(nValue, rhoValues(:), kappaValues(:), rhoReciprocals(:), kappaReciprocals(:))

  ! Compute mass and rigitidy matrices.
  iSolid = 0
  iFluid = 0
  do i = 1, nZone
    oV = oValueOfZone(i)
    if (phaseOfZone(i) == 1) then
      ! solid
      iSolid = iSolid + 1
      oR = oRowOfZoneSolid(iSolid)

      ! Compute unmodified matrices.
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), 2, 0, 0, t(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), 0, 0, 0, h1x(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 0, 0, 0, h2L(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 0, 0, 0, h2N(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), 1, 0, 1, hUn5y(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 1, 0, 1, hUn6L(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 1, 0, 1, hUn6N(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), 2, 1, 1, h7y(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), 2, 1, 1, h7z(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 2, 1, 1, h8L(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 2, 1, 1, h8N(oR:))
      call computeTranspose(nLayerInZone(i), hUn5y(oR:), hUn3y(oR:))
      call computeTranspose(nLayerInZone(i), hUn6L(oR:), hUn4L(oR:))
      call computeTranspose(nLayerInZone(i), hUn6N(oR:), hUn4N(oR:))
      ! Modify matrices for I0 and I2 using lumped matrices.
      call computeLumpedT(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), t(oR:), work(oR:), t(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h1x(oR:), work(oR:), h1x(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h2L(oR:), work(oR:), h2L(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h2N(oR:), work(oR:), h2N(oR:))

    else
      ! fluid
      iFluid = iFluid + 1
      oR = oRowOfZoneFluid(iFluid)

      ! Compute unmodified matrices.
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), 2, 1, 1, p1(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), 0, 0, 0, p2(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), kappaReciprocals(oV:), 2, 0, 0, p3(oR:))
      ! Modify matrices for I0 and I2 using lumped matrices.
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), p2(oR:), work(oR:), p2(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), kappaReciprocals(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), p3(oR:), work(oR:), p3(oR:))

    end if
  end do

  ! Compute the modified operator of the 1st derivative.
  iSolid = 0
  do i = 1, nZone
    oV = oValueOfZone(i)
    if (phaseOfZone(i) == 1) then
      ! solid
      iSolid = iSolid + 1
      oR = oRowOfZoneSolid(iSolid)
      oVS = oValueOfZoneSolid(iSolid)

      ! Compute residual after subtracting step-wise matrix from unmodified matrix.
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn5y(oR:), work(oR:), hResid5y(oR:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn6L(oR:), work(oR:), hResid6L(oR:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn6N(oR:), work(oR:), hResid6N(oR:))
      call computeTranspose(nLayerInZone(i), hResid5y(oR:), hResid3y(oR:))
      call computeTranspose(nLayerInZone(i), hResid6L(oR:), hResid4L(oR:))
      call computeTranspose(nLayerInZone(i), hResid6N(oR:), hResid4N(oR:))

      ! Compute modified matrices.
      call computeModifiedH1(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), h5Mod1y(-1:2, oVS:))
      call computeModifiedH1(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), h6Mod1N(-1:2, oVS:))
      call computeModifiedH2(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), h6Mod2L(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -1, 2, h5Mod1y(-1:2, oVS:), h3Mod2y(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -1, 2, h6Mod1N(-1:2, oVS:), h4Mod2N(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -2, 1, h6Mod2L(-2:1, oVS:), h4Mod1L(-1:2, oVS:))

    end if
  end do

  !<<< operation up to here >>>>

  !******************** Computing the displacement *********************
  outputCounter = 1  !!! difference from shallow-source section

  do iFreq = imin, imax  ! omega-loop
    omega = 2.d0 * pi * dble(iFreq) / tlen

    ! Initialize matrices.
    a0(:, :nColumn) = dcmplx(0.d0, 0.d0)
    a1(:, :nColumn) = dcmplx(0.d0, 0.d0)
    a2(:, :nColumn) = dcmplx(0.d0, 0.d0)
    u(:, :nReceiver) = dcmplx(0.d0, 0.d0)
    ! Plm must be cleared for each omega.  !!! difference from shallow-source section
    plm(:, :, :nReceiver) = 0.d0

    ! Compute the angular order that is sufficient to compute the slowest phase velocity.
    call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

    ! Compute coefficients to multiply to elastic moduli for anelastic attenuation.
    call computeCoef(nZone, omega, qmuOfZone(:), qkappaOfZone(:), coefQmu(:), coefQkappa(:), coefQfluid(:))


    ! call calabnum
    !!TODO



    ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
    iSolid = 0
    iFluid = 0
    do i = 1, nZone
      oElement = oElementOfZone(i)
      oColumn = oColumnOfZone(i)

      if (phaseOfZone(i) == 1) then
        ! solid
        iSolid = iSolid + 1
        oR = oRowOfZoneSolid(iSolid)
        oVS = oValueOfZoneSolid(iSolid)

        ! All parts of A0 are either unmodified or already modified using lumped matrix.
        call computeA0Solid(nLayerInZone(i), omega, omegaI, t(oR:), h1x(oR:), h2L(oR:), h2N(oR:), &
          hUn3y(oR:), hUn4L(oR:), hUn4N(oR:), hUn5y(oR:), hUn6L(oR:), hUn6N(oR:), h7y(oR:), h7z(oR:), h8L(oR:), h8N(oR:), &
          coefQmu(i), coefQkappa(i), cwork(oElement:))
        call overlapASolid(nLayerInZone(i), cwork(oElement:), a0(:, oColumn:))
        ! All parts of A2 are either unmodified or already modified using lumped matrix.
        call computeA2Solid(nLayerInZone(i), h1x(oR:), h2L(oR:), h2N(oR:), coefQmu(i), coefQkappa(i), cwork(oElement:))
        call overlapASolid(nLayerInZone(i), cwork(oElement:), a2(:,oColumn:))
        ! Unmodified residual part of A1.
        call computeA1Solid(nLayerInZone(i), h1x(oR:), h2L(oR:), h2N(oR:), hResid3y(oR:), hResid4L(oR:), hResid4N(oR:), &
          hResid5y(oR:), hResid6L(oR:), hResid6N(oR:), coefQmu(i), coefQkappa(i), cwork(oElement:))
        call overlapASolid(nLayerInZone(i), cwork(oElement:), a1(:, oColumn:))
        ! Modified step part of A1.
        call addModifiedHToA1(nLayerInZone(i), coefQmu(i), coefQkappa(i), &
          h3Mod2y(-2:1, oVS:), h4Mod1L(-1:2, oVS:), h4Mod2N(-2:1, oVS:), &
          h5Mod1y(-1:2, oVS:), h6Mod1N(-1:2, oVS:), h6Mod2L(-2:1, oVS:), a1(:, oColumn:))

      else
        ! fluid
        iFluid = iFluid + 1
        oR = oRowOfZoneFluid(iFluid)

        ! All parts of A0 are either unmodified or already modified using lumped matrix.
        call computeA0Fluid(nLayerInZone(i), omega, omegaI, p1(oR:), p3(oR:), coefQfluid(i), cwork(oElement:))
        call overlapAFluid(nLayerInZone(i), cwork(oElement:), a0(:,oColumn:))
        ! All parts of A2 are either unmodified or already modified using lumped matrix.
        call computeA2Fluid(nLayerInZone(i), omega, omegaI, p2(oR:), cwork(oElement:))
        call overlapAFluid(nLayerInZone(i), cwork(oElement:), a2(:,oColumn:))

      end if
    end do

    ! Initially, no depth cut-off, so set to the index of deepest grid, which is 1.
    cutoffGrid = 1
    ! Clear counter.
    decayCounter = 0
    ! Clear amplitude record.
    recordAmplitude = -1.d0

    do l = 0, maxL  ! l-loop
      ! When the counter detecting the decay of amplitude has reached a threshold, stop l-loop for this frequency.
      if (decayCounter > 20) exit

      ! L^2 and L. (See the part after eq. 12 of Kawai et al. 2006.)
      ! NOTE that integers are casted with dble() before multiplying, because the product can exceed the size of integer(4).
      largeL2 = dble(l) * dble(l + 1)
!      largeL = sqrt(largeL2)

      ! Initialize matrices.
      a(:, :nColumn) = dcmplx(0.d0, 0.d0)
!      aSource(:, :) = dcmplx(0.d0, 0.d0)
      ! Clear the amplitude accumulated for all m's.
      if (mod(l, 100) == 0) amplitudeAtGrid(:nGrid) = 0.d0

      ! Compute trial functions.  !!! difference from shallow-source section
      do ir = 1, nReceiver
        call computeHarmonicsValues(l, theta(ir), phi(ir), plm(:, :, ir), harmonicsValues(:, :, ir))
      end do

      ! Assemble A matrix from parts that have already been computed.
      call assembleAWhole(nZone, phaseOfZone(:), oColumnOfZone(:), largeL2, a0(:,:), a1(:,:), a2(:,:), a(:,:))







      do m = -2, 2  ! m-loop
        if (abs(m) > abs(l)) cycle







      end do  ! m-loop

      ! Decide cut-off depth (at a certain interval of l).
      if (mod(l, 100) == 0) then
        call computeCutoffDepth(nGrid, amplitudeAtGrid(:), ratc, cutoffGrid)
      end if

    end do  ! l-loop

    ! Register the final l (or maxL instead of maxL-1 when all loops are completed).
    llog = min(l, maxL)

    ! Store results.  !!! difference from shallow-source section
    outputi(outputCounter) = iFreq
    do ir = 1, nReceiver
      outputu(:, ir, outputCounter) = u(:, ir)
    end do


    ! ************************** Files Handling **************************
    ! Write to file when the output interval is reached, or when this is the last omega.
    if (outputCounter >= outputInterval .or. iFreq == imax) then
      write(*,*) "kakikomimasu"

      do ir = 1, nReceiver
        call openSPCFile(output(ir), 11, spcFormat, 1)
        do iOut = 1, outputCounter
          call writeSPCFile(11, spcFormat, outputi(iOut), dble(outputu(1, ir, iOut)), imag(outputu(1, ir, iOut)))
          call writeSPCFile(11, spcFormat, dble(outputu(2, ir, iOut)), imag(outputu(2, ir, iOut)))
          call writeSPCFile(11, spcFormat, dble(outputu(3, ir, iOut)), imag(outputu(3, ir, iOut)))
        end do
        call closeSPCFile(11)
      end do

      outputCounter = 0
    end if

    if (ilog == 1) then
      open(unit=11, file='llog.log', position='append', status='old')
      write(11,*) iFreq, llog, nGrid-1
      close(11)
    end if

    outputCounter = outputCounter + 1

  end do  ! omega-loop

  ! Deallocate arrays.
  deallocate(outputi)
  deallocate(outputu)

  write(*,*) "Ivalice looks to the horizon"

  stop
end program tipsv
