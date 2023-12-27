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
  integer, parameter :: maxNLayerSolid = 48840  ! Maximum number of layers in solid region.
  integer, parameter :: maxNLayerLiquid = 32040  ! Maximum number of layers in liquid region.
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
  integer :: nZoneSolid, nZoneLiquid  ! Number of solid and liquid zones.
  integer :: phaseOfZone(maxNZone)  ! Phase of each zone (1: solid, 2: liquid).
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled [km].
  real(8) :: rminOfZone(maxNZone), rmaxOfZone(maxNZone)  ! Minimum and maximum radii of each zone [km].
  real(8) :: rhoPolynomials(4, maxNZone), vpvPolynomials(4, maxNZone), vphPolynomials(4, maxNZone)
  real(8) :: vsvPolynomials(4, maxNZone), vshPolynomials(4, maxNZone), etaPolynomials(4, maxNZone)
  !::::: Polynomial functions (coefficients of cubic function) of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta in each zone.
  real(8) :: qmuOfZone(maxNZone), qkappaOfZone(maxNZone)  ! Qmu and Qkappa of each zone.
  integer :: i, iS, iL

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
  integer :: nLayerSolid, nLayerLiquid  ! Number of layers in solid and liquid regions.
  integer :: oGridOfZone(maxNZone)  ! Index of the first grid point in each zone.
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
  complex(8) :: qCoef(maxNZone)  ! Coefficient to multiply to elastic moduli for attenuation at each zone.
  integer :: oV, oVS

  ! Variables for the trial function
  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: largeL2  ! L^2 = l(l+1).
  real(8) :: plm(3, 0:3, maxNReceiver)  ! Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  !::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8) :: harmonicsValues(3, -2:2, maxNReceiver)  ! Values of vector harmonics terms at each receiver, computed for each l.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Arguments: term (1:3), m (-2:2), iReceiver.

  ! Variables for the matrix elements
  real(8) :: t(4 * maxNGrid - 4)  !!!TODO change to 4*maxNLayerSolid
  real(8) :: h1x(4 * maxNGrid - 4), h1y(4 * maxNGrid - 4), h1z(4 * maxNGrid - 4), h2L(4 * maxNGrid - 4), h2N(4 * maxNGrid - 4)
  real(8) :: hUn3x(4 * maxNGrid - 4), hUn3y(4 * maxNGrid - 4), hUn3z(4 * maxNGrid - 4)
  real(8) :: hUn4L(4 * maxNGrid - 4), hUn4N(4 * maxNGrid - 4)
  real(8) :: hResid3x(4 * maxNGrid - 4), hResid3y(4 * maxNGrid - 4), hResid3z(4 * maxNGrid - 4)
  real(8) :: hResid4L(4 * maxNGrid - 4), hResid4N(4 * maxNGrid - 4)
  real(8) :: h3Mod2x(-2:1, maxNLayerSolid + maxNZone), h3Mod2y(-2:1, maxNLayerSolid + maxNZone), &
    h3Mod2z(-2:1, maxNLayerSolid + maxNZone)
  real(8) :: h4Mod1L(-1:2, maxNLayerSolid + maxNZone), h4Mod1N(-1:2, maxNLayerSolid + maxNZone)
  real(8) :: h4Mod2L(-2:1, maxNLayerSolid + maxNZone), h4Mod2N(-2:1, maxNLayerSolid + maxNZone)
  real(8) :: hUn5x(4 * maxNGrid - 4), hUn5y(4 * maxNGrid - 4), hUn5z(4 * maxNGrid - 4)
  real(8) :: hUn6L(4 * maxNGrid - 4), hUn6N(4 * maxNGrid - 4)
  real(8) :: hResid5x(4 * maxNGrid - 4), hResid5y(4 * maxNGrid - 4), hResid5z(4 * maxNGrid - 4)
  real(8) :: hResid6L(4 * maxNGrid - 4), hResid6N(4 * maxNGrid - 4)
  real(8) :: h5Mod1x(-1:2, maxNLayerSolid + maxNZone), h5Mod1y(-1:2, maxNLayerSolid + maxNZone), &
    h5Mod1z(-1:2, maxNLayerSolid + maxNZone)
  real(8) :: h6Mod1L(-1:2, maxNLayerSolid + maxNZone), h6Mod1N(-1:2, maxNLayerSolid + maxNZone)
  real(8) :: h6Mod2L(-2:1, maxNLayerSolid + maxNZone), h6Mod2N(-2:1, maxNLayerSolid + maxNZone)
  real(8) :: h7x(4 * maxNGrid - 4), h7y(4 * maxNGrid - 4), h7z(4 * maxNGrid - 4), h8L(4 * maxNGrid - 4), h8N(4 * maxNGrid - 4)
  real(8) :: p1(4 * maxNGrid - 4), p2(4 * maxNGrid - 4), p3(4 * maxNGrid - 4)  !!!TODO change to 4*maxNLayerLiquid
  real(8) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  integer :: oRowOfZoneSolid(maxNZone), oRowOfZoneLiquid(maxNZone)
  !:: Index of the first row in the vector of (iLayer, k', k)-pairs in each zone. Vectors are separate for solid and liquid zones.
  integer :: oRowOfSource  ! Index of the first row in the vector of (iLayer, k', k)-pairs for the layer with the source.
  complex(8) :: a0(2, maxNGrid), a2(2, maxNGrid)
  complex(8) :: a(2, maxNGrid)
  complex(8) :: aaParts(4), aSourceParts(8), aSource(2, 3)
  complex(8) :: g_or_c(maxNGrid)  ! This holds either vector g [10^15 N] or c [km], depending on where in the code it is. CAUTION!!
  complex(8) :: u(3, maxNReceiver)  ! Displacement velocity - the unit is [km] in the frequency domain,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: but when converted to the time domain, the unit becomes [km/s].
  integer :: oR

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Other variables
  real(8) :: work(4 * maxNGrid - 4)  ! Working array for matrix computations.


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

  call judgeSolidOrLiquid(nZone, vsvPolynomials, phaseOfZone, nZoneSolid, nZoneLiquid)

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
  call computeGridRadii(nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), phaseOfZone(:), rmin, re, &
    nGrid, nLayerInZone(:), nLayerSolid, nLayerLiquid, gridRadii(:))
  if (nGrid > maxNGrid) stop 'The number of grid points is too large.'
  if (nLayerSolid > maxNLayerSolid) stop 'The number of solid layers is too large.'
  if (nLayerLiquid > maxNLayerLiquid) stop 'The number of liquid layers is too large.'

  ! Compute the first indices in each zone.
!  call computeFirstIndices(nZone, nLayerInZone(:), phaseOfZone(:), oGridOfZone(:), oValueOfZone(:), oValueOfZoneSolid(:), &
!    oRowOfZoneSolid(:), oRowOfZoneLiquid(:))

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
  iS = 0
  iL = 0
  do i = 1, nZone
    oV = oValueOfZone(i)
    if (phaseOfZone(i) == 1) then
      ! solid
      iS = iS + 1
      oR = oRowOfZoneSolid(iS)

      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), 2, 0, 0, t(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), 0, 0, 0, h1x(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), 0, 0, 0, h1y(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), 0, 0, 0, h1z(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 0, 0, 0, h2L(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 0, 0, 0, h2N(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), 1, 0, 1, hUn5x(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), 1, 0, 1, hUn5y(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), 1, 0, 1, hUn5z(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 1, 0, 1, hUn6L(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 1, 0, 1, hUn6N(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), 2, 1, 1, h7x(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), 2, 1, 1, h7y(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), 2, 1, 1, h7z(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 2, 1, 1, h8L(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 2, 1, 1, h8N(oR:))
      call computeLumpedT(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), t(oR:), work(oR:), t(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h1x(oR:), work(oR:), h1x(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h1y(oR:), work(oR:), h1y(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h1z(oR:), work(oR:), h1z(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h2L(oR:), work(oR:), h2L(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), h2N(oR:), work(oR:), h2N(oR:))
      call computeTranspose(nLayerInZone(i), hUn5x(oR:), hUn3x(oR:))
      call computeTranspose(nLayerInZone(i), hUn5y(oR:), hUn3y(oR:))
      call computeTranspose(nLayerInZone(i), hUn5z(oR:), hUn3z(oR:))
      call computeTranspose(nLayerInZone(i), hUn6L(oR:), hUn4L(oR:))
      call computeTranspose(nLayerInZone(i), hUn6N(oR:), hUn4N(oR:))

    else
      !liquid
      iL = iL + 1
      oR = oRowOfZoneLiquid(iL)

      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), 2, 1, 1, p1(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), 0, 0, 0, p2(oR:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), kappaReciprocals(oV:), 2, 0, 0, p3(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), p2(oR:), work(oR:), p2(oR:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), kappaReciprocals(oV:), work(oR:))
      call computeAverage(nLayerInZone(i), p3(oR:), work(oR:), p3(oR:))

    end if
  end do

  ! Compute the modified operator of the 1st derivative.
  iS = 0
  do i = 1, nZone
    oV = oValueOfZone(i)
    if (phaseOfZone(i) == 1) then
      ! solid
      iS = iS + 1
      oR = oRowOfZoneSolid(iS)
      oVS = oValueOfZoneSolid(iS)

      ! Compute residual after subtracting step-wise matrix from unmodified matrix.
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn5x(oR:), work(oR:), hResid5x(oR:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn5y(oR:), work(oR:), hResid5y(oR:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn5z(oR:), work(oR:), hResid5z(oR:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn6L(oR:), work(oR:), hResid6L(oR:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), work(oR:))
      call subtractMatrix(nLayerInZone(i), hUn6N(oR:), work(oR:), hResid6N(oR:))
      call computeTranspose(nLayerInZone(i), hResid5x(oR:), hResid3x(oR:))
      call computeTranspose(nLayerInZone(i), hResid5y(oR:), hResid3y(oR:))
      call computeTranspose(nLayerInZone(i), hResid5z(oR:), hResid3z(oR:))
      call computeTranspose(nLayerInZone(i), hResid6L(oR:), hResid4L(oR:))
      call computeTranspose(nLayerInZone(i), hResid6N(oR:), hResid4N(oR:))

      ! Compute modified matrices.
      call computeModifiedH1(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), h5Mod1x(-1:2, oVS:))
      call computeModifiedH1(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), h5Mod1y(-1:2, oVS:))
      call computeModifiedH1(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), h5Mod1z(-1:2, oVS:))
      call computeModifiedH1(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), h6Mod1L(-1:2, oVS:))
      call computeModifiedH1(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), h6Mod1N(-1:2, oVS:))
      call computeModifiedH2(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), h6Mod2L(-2:1, oVS:))
      call computeModifiedH2(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), h6Mod2N(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -1, 2, h5Mod1x(-1:2, oVS:), h3Mod2x(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -1, 2, h5Mod1y(-1:2, oVS:), h3Mod2y(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -1, 2, h5Mod1z(-1:2, oVS:), h3Mod2z(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -1, 2, h6Mod1L(-1:2, oVS:), h4Mod2L(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -1, 2, h6Mod1N(-1:2, oVS:), h4Mod2N(-2:1, oVS:))
      call computeTransposeMod(nLayerInZone(i), -2, 1, h6Mod2L(-2:1, oVS:), h4Mod1L(-1:2, oVS:))
      call computeTransposeMod(nLayerInZone(i), -2, 1, h6Mod2N(-2:1, oVS:), h4Mod1N(-1:2, oVS:))

    end if
  end do

  !<<< operation up to here >>>>

  !******************** Computing the displacement *********************
  outputCounter = 1  !!! difference from shallow-source section

  do iFreq = imin, imax  ! omega-loop
    omega = 2.d0 * pi * dble(iFreq) / tlen









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
