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
  real(8) :: amplitudeAtColumn(2 * maxNGridSolid + maxNGridFluid)
  !:::::::::::::::::::::::::::: Estimate of the amplitude at each column [km], used for vertical grid cut-off.
  integer :: lsuf  ! Accuracy threshold of angular order. (Corresponds to l_d; see eq. 29 of Kawai et al. 2006.)
  real(8) :: recordAmplitude    ! Maximum amplitude encountered [km], used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: llog

  ! Variables for the vertical grid
  integer :: nGrid  ! Total number of grid points.
  real(8) :: gridRadii(maxNGrid)  ! Radii of each grid point [km].
  integer :: nLayerInZone(maxNZone)  ! Number of layers in each zone.
  integer :: oGridOfZone(maxNZone)  ! Index of the first grid point in each zone.
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
  integer :: oVS

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
  real(8) :: hUn3y(4 * maxNGridSolid - 4), hResid3y(4 * maxNGridSolid - 4), hModL3y(-2:1, maxNGridSolid)
  real(8) :: hUn4L(4 * maxNGridSolid - 4), hResid4L(4 * maxNGridSolid - 4), hModR4L(-1:2, maxNGridSolid)
  real(8) :: hUn4N(4 * maxNGridSolid - 4), hResid4N(4 * maxNGridSolid - 4), hModL4N(-2:1, maxNGridSolid)
  real(8) :: hUn5y(4 * maxNGridSolid - 4), hResid5y(4 * maxNGridSolid - 4), hModR5y(-1:2, maxNGridSolid)
  real(8) :: hUn6L(4 * maxNGridSolid - 4), hResid6L(4 * maxNGridSolid - 4), hModL6L(-2:1, maxNGridSolid)
  real(8) :: hUn6N(4 * maxNGridSolid - 4), hResid6N(4 * maxNGridSolid - 4), hModR6N(-1:2, maxNGridSolid)
  real(8) :: h7y(4 * maxNGridSolid - 4), h7z(4 * maxNGridSolid - 4), h8L(4 * maxNGridSolid - 4), h8N(4 * maxNGridSolid - 4)
  real(8) :: p1(4 * maxNGridFluid - 4), p2(4 * maxNGridFluid - 4), p3(4 * maxNGridFluid - 4)
  integer :: oPairOfZoneSolid(maxNZone), oPairOfZoneFluid(maxNZone)
  !::::::::::::::::::::: Index of the first (iLayer, k', k)-pair in each zone, counted separately for solid and fluid zones.
  complex(8) :: a0(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a1(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a2(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: aSmall(2, maxNGridSolid + maxNGridFluid)
  complex(8) :: g_or_c(2 * maxNGridSolid + maxNGridFluid)
  !::::::::::::::::::::::::::::: This holds either vector g [10^15 N] or c [km], depending on where in the code it is. CAUTION!!
  complex(8) :: g_or_c_Small(maxNGridSolid + maxNGridFluid)
  !::::::::::::::::::::::::::::: This holds either vector g [10^15 N] or c [km], depending on where in the code it is. CAUTION!!
  complex(8) :: u(3, maxNReceiver)  ! Displacement velocity - the unit is [km] in the frequency domain,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: but when converted to the time domain, the unit becomes [km/s].
  integer :: oElementOfZone(maxNZone)  ! Index of the first (iLayer, k'-gamma', k-gamma)-pair in each zone.
  integer :: oColumnOfZone(maxNZone + 1)  ! Index of the first column in the band matrix for each zone.
  integer :: nColumn, nQuasiColumn  ! Total number of columns in the band matrix.
  integer :: cutoffColumn  ! Index of column at cut-off depth.
  integer :: oP, oElement, oColumn

  !!TODO ???
  complex(8) :: anum(4, 4, 10), bnum(4, 4, 10)
  complex(8) :: ya(4), yb(4), yc(4), yd(4)
  integer :: oColumnOfSource  ! Index of the first column in the band matrix for the layer with source.

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Other variables
  real(8) :: work(4 * maxNGrid - 4)  ! Working array for matrix computations.
  complex(8) :: cwork(16 * maxNGridSolid - 16 + 4 * maxNGridFluid - 4)  ! Working array for matrix computations.
  complex(8) :: z(2 * maxNGridSolid + maxNGridFluid), w(2 * maxNGridSolid + maxNGridFluid)
  !::::::::::::::::::::::::::::::::::: Working arrays used when solving linear equations.

  ! Constants
  real(8) :: eps = -1.d0

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
  call computeMatrixElements(maxNGrid, maxNGridSolid, maxNGridFluid, tlen, re, imaxFixed, r0, &
    nZone, rmin, rmax, rminOfZone, rmaxOfZone, phaseOfZone, &
    rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, &
    kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oValueOfZoneSolid, &
    oPairOfZoneSolid, oPairOfZoneFluid, oElementOfZone, oColumnOfZone, nColumn, &
    iZoneOfSource, iLayerOfSource, oColumnOfSource, &
    nValue, valuedRadii, rhoValues, kappaValues, ecKxValues, ecKyValues, ecKzValues, ecLValues, ecNValues, ecC0, ecF0, ecL0, &
    rhoReciprocals, kappaReciprocals, &
    t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
    hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, work)


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

    !!TODO organize
    call calabnum(omega, omegaI, rmax, &
      rhoPolynomials(:, iZoneOfSource), vpvPolynomials(:, iZoneOfSource), vphPolynomials(:, iZoneOfSource), &
      vsvPolynomials(:, iZoneOfSource), vshPolynomials(:, iZoneOfSource), etaPolynomials(:, iZoneOfSource), &
      gridRadii(iLayerOfSource), r0, coefQmu(iZoneOfSource), coefQkappa(iZoneOfSource), anum(:, :, :), bnum(:, :, :) )

    ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
    iSolid = 0
    iFluid = 0
    do i = 1, nZone
      oElement = oElementOfZone(i)
      oColumn = oColumnOfZone(i)

      if (phaseOfZone(i) == 1) then
        ! solid
        iSolid = iSolid + 1
        oP = oPairOfZoneSolid(iSolid)
        oVS = oValueOfZoneSolid(iSolid)

        ! All parts of A0 are either unmodified or already modified using lumped matrix.
        call computeA0Solid(nLayerInZone(i), omega, omegaI, t(oP:), h1x(oP:), h2L(oP:), h2N(oP:), &
          hUn3y(oP:), hUn4L(oP:), hUn4N(oP:), hUn5y(oP:), hUn6L(oP:), hUn6N(oP:), h7y(oP:), h7z(oP:), h8L(oP:), h8N(oP:), &
          coefQmu(i), coefQkappa(i), cwork(oElement:))
        call overlapASolid(nLayerInZone(i), cwork(oElement:), a0(:, oColumn:))
        ! All parts of A2 are either unmodified or already modified using lumped matrix.
        call computeA2Solid(nLayerInZone(i), h1x(oP:), h2L(oP:), h2N(oP:), coefQmu(i), coefQkappa(i), cwork(oElement:))
        call overlapASolid(nLayerInZone(i), cwork(oElement:), a2(:,oColumn:))
        ! Unmodified residual part of A1.
        call computeA1Solid(nLayerInZone(i), h1x(oP:), h2L(oP:), h2N(oP:), hResid3y(oP:), hResid4L(oP:), hResid4N(oP:), &
          hResid5y(oP:), hResid6L(oP:), hResid6N(oP:), coefQmu(i), coefQkappa(i), cwork(oElement:))
        call overlapASolid(nLayerInZone(i), cwork(oElement:), a1(:, oColumn:))
        ! Modified step part of A1.
        call addModifiedHToA1(nLayerInZone(i), coefQmu(i), coefQkappa(i), &
          hModL3y(-2:1, oVS:), hModR4L(-1:2, oVS:), hModL4N(-2:1, oVS:), &
          hModR5y(-1:2, oVS:), hModL6L(-2:1, oVS:), hModR6N(-1:2, oVS:), a1(:, oColumn:))

      else
        ! fluid
        iFluid = iFluid + 1
        oP = oPairOfZoneFluid(iFluid)

        ! All parts of A0 are either unmodified or already modified using lumped matrix.
        call computeA0Fluid(nLayerInZone(i), omega, omegaI, p1(oP:), p3(oP:), coefQfluid(i), cwork(oElement:))
        call overlapAFluid(nLayerInZone(i), cwork(oElement:), a0(:, oColumn:))
        ! All parts of A2 are either unmodified or already modified using lumped matrix.
        call computeA2Fluid(nLayerInZone(i), omega, omegaI, p2(oP:), cwork(oElement:))
        call overlapAFluid(nLayerInZone(i), cwork(oElement:), a2(:, oColumn:))

      end if
    end do

    ! Initially, no depth cut-off, so set to the column of deepest grid, which is 1.
    cutoffColumn = 1
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
      largeL = sqrt(largeL2)

      ! Initialize matrices.
      a(:, :nColumn) = dcmplx(0.d0, 0.d0)
      ! Clear the amplitude accumulated for all m's.
      if (mod(l, 100) == 0) amplitudeAtColumn(:nColumn) = 0.d0

      ! Compute trial functions.  !!! difference from shallow-source section
      do ir = 1, nReceiver
        call computeHarmonicsValues(l, theta(ir), phi(ir), plm(:, :, ir), harmonicsValues(:, :, ir))
      end do

      ! Assemble A matrix from parts that have already been computed.
      call assembleAWhole(nZone, phaseOfZone(:), oColumnOfZone(:), largeL2, a0(:,:), a1(:,:), a2(:,:), a(:,:))
      ! Set boundary condition elements
      call setBoundaryConditions(nZone, rmaxOfZone(:), phaseOfZone(:), oColumnOfZone(:), a(:,:))

      !!TODO organize
      call calya(anum(:,:,:), bnum(:,:,:), largeL2, gridRadii(iLayerOfSource:), r0, ya(:), yb(:), yc(:), yd(:))

      do m = -2, 2  ! m-loop
        if (abs(m) > abs(l)) cycle

        ! Form and solve the linear equation Ac=-g.
        call formAndSolveEquation(l, m, largeL, iZoneOfSource, iLayerOfSource, oColumnOfSource, r0, mt, ecC0, ecF0, ecL0, &
          ya, yb, yc, yd, rmin, nZone, phaseOfZone, oGridOfZone, oColumnOfZone, coefQmu, coefQkappa, &
          nGrid, gridRadii, nColumn, cutoffColumn, a, aSmall, g_or_c, g_or_c_Small, amplitudeAtColumn, nQuasiColumn, eps, z, w)

        if (l == 0) then
          ! Record u.  !!! difference from shallow-source section
          do ir = 1, nReceiver
            u(1, ir) = g_or_c_Small(nQuasiColumn) * harmonicsValues(1, m, ir)
          end do

        else
          ! Check whether the amplitude has decayed enough to stop the l-loops.
          !  This is checked for the topmost-grid expansion coefficent of each m individually.
          call checkAmplitudeDecay(g_or_c(nColumn-1:nColumn), l, lsuf, ratl, recordAmplitude, decayCounter)

          ! Accumulate u.  !!! difference from shallow-source section
          do ir = 1, nReceiver
            call computeU(g_or_c(nColumn-1:nColumn), largeL2, harmonicsValues(:, m, ir), u(:, ir))
          end do

        end if

      end do  ! m-loop

      ! Decide cut-off depth (at a certain interval of l).
      if (mod(l, 100) == 0) then
        call computeCutoffColumn(nZone, phaseOfZone(:), nGrid, oGridOfZone(:), nColumn, oColumnOfZone(:), &
          amplitudeAtColumn(:), ratc, cutoffColumn)
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
