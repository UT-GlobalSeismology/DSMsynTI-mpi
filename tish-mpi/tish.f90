!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ************** tish ****************
!  Computation of SH synthetic seismograms
!  in transversely isotropic media for anisotropic PREM
!  using modified DSM operators & modified source representation.
!  Synthetics for shallow events can be computed.
!
!  Main historical authors: K.Kawai, N.Takeuchi, R.J.Geller
!  (C) 2002.10  University of Tokyo
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program tish
  implicit none

  !----------------------------<<constants>>----------------------------
  real(8), parameter :: pi = 3.1415926535897932d0
  integer, parameter :: maxNLayer = 88300  ! Maximum number of layers.
  integer, parameter :: maxNZone = 15  ! Maximum number of zones.
  integer, parameter :: maxNReceiver = 1500  ! Maximum number of receivers.
  integer, parameter :: maxL = 80000  ! Maximum of angular order to loop for.
  real(8), parameter :: lmaxdivf = 2.d4
  real(8), parameter :: shallowdepth = 100.d0
  integer, parameter :: spcFormat = 0  ! Format of output spc file (0:binary, 1:ascii).
  integer, parameter :: ilog = 0

  !----------------------------<<variables>>----------------------------
  ! Variables for the structure
  integer :: nZone  ! Number of zones.
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled.
  real(8) :: rminOfZone(maxNZone), rmaxOfZone(maxNZone)  ! Minimum and maximum radii of each zone.
  real(8) :: rhoPolynomials(4, maxNZone)  ! Rho of each zone (coefficients of cubic function).
  real(8) :: vsvPolynomials(4, maxNZone)  ! Vsv of each zone (coefficients of cubic function).
  real(8) :: vshPolynomials(4, maxNZone)  ! Vsh of each zone (coefficients of cubic function).
  real(8) :: qmuOfZone(maxNZone)  ! Qmu of each zone.

  ! Variables for the source
  real(8) :: r0, mt(3, 3), eqlat, eqlon, mu0

  ! Variables for the receivers
  integer :: nReceiver  ! Number of receivers.
  real(8) :: theta(maxNReceiver), phi(maxNReceiver)
  real(8) :: lat(maxNReceiver), lon(maxNReceiver)
  integer :: ir

  ! Variables for the periodic range
  real(8) :: tlen  ! Time length.
  integer :: np  ! Number of points in frequency domain.
  real(8) :: omega  ! Angular frequency (real part).
  real(8) :: omegai  ! Imaginary part of angular frequency for artificial damping. (See section 5.1 of Geller & Ohminato 1994.)
  integer :: imin, imax  ! Index of minimum and maximum frequency.

  ! Variables for grid spacing and cut-off
  real(8) :: kzAtZone(maxNZone)  ! Vertical wavenumber k_z at each zone. (See section 3.2 of Kawai et al. 2006.)
  real(8) :: re  ! Desired relative error due to vertical gridding. (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
  real(8) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8) :: amplitudeAtGrid(maxNLayer + 1)  ! Estimate of the amplitude at each grid point, used for vertical grid cut-off.
  integer :: cutoffGrid  ! Index of grid at cut-off depth.
  integer :: lsuf  ! Accuracy threshold of angular order. (Corresponds to l_d; see eq. 29 of Kawai et al. 2006.)
  real(8) :: recordAmplitude    ! Maximum amplitude encountered, used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: llog

  ! Variables for the vertical grid
  integer :: nLayer  ! Total number of layers.   !!!!!TODO switch to use nGrid?
  integer :: nGrid  ! Total number of grid points.
  real(8) :: gridRadii(maxNLayer + 1)  ! Radii of each grid point.
  integer :: nLayerInZone(maxNZone)  ! Number of layers in each zone.
  integer :: oGridOfZone(maxNZone)  ! Index of the first grid point in each zone.
  real(8) :: gridRadiiForSource(3)  ! Radii to use for source-related computations.
  integer :: iZoneOfSource  ! Which zone the source is in.
  integer :: iLayerOfSource  ! Index of layer that the source is in.

  ! Variables for the values
  integer :: nValue  ! Total number of values.
  real(8) :: valuedRadii(maxNLayer + maxNZone)  ! Radii corresponding to each variable value.
  real(8) :: rhoValues(maxNLayer + maxNZone)  ! Rho at each grid point (with 2 values at boundaries).
  real(8) :: ecLValues(maxNLayer + maxNZone)  ! L at each grid point (with 2 values at boundaries).
  real(8) :: ecNValues(maxNLayer + maxNZone)  ! N at each grid point (with 2 values at boundaries).
  real(8) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)  ! Rho, L, and N at each source-related grid.
  complex(8) :: coef(maxNZone)

  ! Variables for the trial function
  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: plm(3, 0:3, maxNReceiver)  ! Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  !::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8) :: trialFunctionValues(3, -2:2, maxNReceiver)  ! Values of trial function at each receiver, computed for each l.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Arguments: component (1:3), m (-2:2), iReceiver.

  ! Variables for the stack points
  integer :: oRowOfZone(maxNZone)  ! Index of the first row in the vector of (iLayer, k', k)-pairs in each zone.
  integer :: oRowOfSource  ! Index of the first row in the vector of (iLayer, k', k)-pairs for the layer with the source.

  ! Variables for the matrix elements
  real(8) :: t(4 * maxNLayer)
  real(8) :: h1(4 * maxNLayer), h2(4 * maxNLayer), h3(4 * maxNLayer), h4(4 * maxNLayer)
  real(8) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  complex(8) :: a0(2, maxNLayer + 1), a2(2, maxNLayer + 1)
  complex(8) :: a(2, maxNLayer + 1)
  complex(8) :: aaParts(4), aSourceParts(8), aSource(2, 3)
  complex(8) :: g_or_c(maxNLayer + 1)  ! This holds either the vector g or c, depending on where in the code it is (CAUTION!).
  complex(8) :: u(3, maxNReceiver)

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Other variables
  integer :: i, j, ii
  integer :: ltmp(2), iimax
  real(8) :: work(4 * maxNLayer)  ! Working array for matrix computations.
  complex(8) :: cwork(4 * maxNLayer)  ! Working array for matrix computations.
  integer :: ier  ! Error code from subroutine solving linear equations.
  complex(8) :: dr(maxNLayer + 1), z(maxNLayer + 1), gdr(3)  ! Working arrays used when solving linear equations.

  ! Constants
  integer :: lda = 2
  real(8) :: eps = -1.d0

  ! Efficiency improvement variables
  integer :: outputmemory = 10  ! MB
  integer :: outputinterval
  real(8) :: memoryperomega ! MB
  integer :: outputindex, mpii
  integer, allocatable :: outputi(:)
  complex(8), allocatable :: outputu(:,:,:)
  !     When the values to be output use memory over outputmemory MB,
  !     they are written in output files. The interval is outputinterval.
  !     memoryperomega is the quantity of memory used for one omega step.


  ! ************************** Inputting parameters **************************
  ! --- read parameters ---
  call readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegai, imin, imax, &
    nZone, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, qmuOfZone, &
    r0, eqlat, eqlon, mt, nReceiver, theta, phi, lat, lon, output)

  memoryperomega = 3 * 16 * nReceiver * 0.000001
  outputinterval = outputmemory / memoryperomega
  allocate(outputi(outputinterval))
  allocate(outputu(3, nReceiver, outputinterval))

  ! --- computing the required parameters ---
  rmin = rminOfZone(1)
  rmax = rmaxOfZone(nZone)

  do ir = 1, nReceiver
    theta(ir) = theta(ir) / 180.0d0 * pi ! Convert theta from degrees to radians
    phi(ir) = phi(ir) / 180.0d0 * pi     ! Convert phi from degrees to radians
  end do

  if (r0 < rmin .or. r0 > rmax) then
    stop 'Location of the source is improper.'
  end if

  if (imin == 0) imin = 1

  ! ************************** Files handling **************************
  if (spcFormat == 0) then
    do ir = 1, nReceiver
      open(unit = 11, file = output(ir), status = 'unknown', &
        form = 'unformatted', access = 'stream', convert = 'big_endian')
      write(11) tlen
      write(11) np, 1, 3
      write(11) omegai, lat(ir), lon(ir)
      write(11) eqlat, eqlon, r0
      close(11)
    end do
  else if (spcFormat == 1) then
    do ir = 1, nReceiver
      open(unit = 11, file = output(ir), status = 'unknown')
      write(11, *) tlen
      write(11, *) np, 1, 3
      write(11, *) omegai, lat(ir), lon(ir)
      write(11, *) eqlat, eqlon, r0
      close(11)
    end do
  else
    write(*, *) 'WARNING:(tish.f)  set spcFormat 0 or 1'
  end if

  if (ilog == 1) then
    open(unit = 11, file = 'llog.log', status = 'unknown')
    close(11)
  end if



  ! Decide which omega to use when deciding grid spacing. Usually, this is just the upper limit of omega range.
  iimax = imax
  ! ************************** Option for shallow events **************************
  ! Here, we find the maximum angular order needed for our frequency range. (See fig. 7 of Kawai et al. 2006.)
  if ((rmax - r0) < shallowdepth) then
    ! ******************* Computing parameters *******************
    ! Set a large value so that we can compute using fine grids for this process.
    iimax = int(tlen * 2.d0)  !!!diff

    ! Design the number and position of grid points.
    call computeKz(nZone, rminOfZone(:), rmaxOfZone(:), vsvPolynomials(:,:), rmax, iimax, 1, tlen, kzAtZone(:))
    call computeGridRadii(nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), rmin, re, nLayer, nLayerInZone(:), gridRadii(:))
    if (nLayer > maxNLayer) stop 'The number of grid points is too large.'

    ! Compute the first indices in each zone.
    call computeFirstIndices(nZone, nLayerInZone(:), oGridOfZone(:), oRowOfZone(:))

    ! Compute the source position.
    call computeSourcePosition(nLayer, rmaxOfZone(:), rmin, rmax, gridRadii(:), r0, iZoneOfSource, iLayerOfSource)

    ! Design grids for source computations.
    call computeSourceGrid(gridRadii(:), r0, iLayerOfSource, gridRadiiForSource(:))

    ! ******************* Computing the matrix elements *******************
    ! Compute variable values at grid points.
    call computeStructureValues(nZone, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), nLayerInZone(:), &
      gridRadii(:), nValue, valuedRadii(:), rhoValues(:), ecLValues(:), ecNValues(:))
    call computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), &
      gridRadiiForSource(:), rhoValuesForSource(:), ecLValuesForSource(:), ecNValuesForSource(:), mu0)

    ! Compute mass and rigitidy matrices.
    do i = 1, nZone
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), rhoValues(:), 2, 0, 0, gridRadii(oGridOfZone(i):), &
        t(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), 2, 1, 1, gridRadii(oGridOfZone(i):), &
        h1(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), 1, 1, 0, gridRadii(oGridOfZone(i):), &
        h2(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), 0, 0, 0, gridRadii(oGridOfZone(i):), &
        h3(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecNValues(:), 0, 0, 0, gridRadii(oGridOfZone(i):), &
        h4(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeLumpedT(nLayerInZone(i), nValue, valuedRadii(:), rhoValues(:), gridRadii(oGridOfZone(i):), work(oRowOfZone(i):))
      call computeAverage(nLayerInZone(i), t(oRowOfZone(i):), work(oRowOfZone(i):), t(oRowOfZone(i):))
      call computeLumpedH(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), gridRadii(oGridOfZone(i):), work(oRowOfZone(i):))
      call computeAverage(nLayerInZone(i), h3(oRowOfZone(i):), work(oRowOfZone(i):), h3(oRowOfZone(i):))
      call computeLumpedH(nLayerInZone(i), nValue, valuedRadii(:), ecNValues(:), gridRadii(oGridOfZone(i):), work(oRowOfZone(i):))
      call computeAverage(nLayerInZone(i), h4(oRowOfZone(i):), work(oRowOfZone(i):), h4(oRowOfZone(i):))
    end do

    ! Compute mass and rigitidy matrices near source.
    call computeIntermediateIntegral(2, 3, gridRadiiForSource, rhoValuesForSource, 2, 0, 0, gridRadiiForSource, gt, work)
    call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecLValuesForSource, 2, 1, 1, gridRadiiForSource, gh1, work)
    call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecLValuesForSource, 1, 1, 0, gridRadiiForSource, gh2, work)
    call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecLValuesForSource, 0, 0, 0, gridRadiiForSource, gh3, work)
    call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecNValuesForSource, 0, 0, 0, gridRadiiForSource, gh4, work)
    call computeLumpedT(2, 3, gridRadiiForSource, rhoValuesForSource, gridRadiiForSource, work)
    call computeAverage(2, gh3, work, gh3)
    call computeLumpedH(2, 3, gridRadiiForSource, ecLValuesForSource, gridRadiiForSource, work)
    call computeAverage(2, gh3, work, gh3)
    call computeLumpedH(2, 3, gridRadiiForSource, ecNValuesForSource, gridRadiiForSource, work)
    call computeAverage(2, gh4, work, gh4)

    nGrid = nLayer + 1
    ! Find the first index of (iLayer, k', k)-pair corresponding to the layer that the source is in.
    oRowOfSource = 4 * iLayerOfSource - 3

    llog = 0

    ! Find the maximum angular order needed for the lowest and highest frequencies. (See fig. 7 of Kawai et al. 2006.)
    do ii = 1, 2  ! omega-loop
      if (ii == 1) then  !!!diff
        i = imin
      else
        i = imax
      end if
      omega = 2.d0 * pi * dble(i) / tlen

      ! Initialize matrices.
      call initComplexMatrix(lda, nGrid, a0(:,:))
      call initComplexMatrix(lda, nGrid, a2(:,:))

      ! Compute the angular order that is sufficient to compute the slowest phase velocity.
      call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

      ! Compute coefficient related to attenuation.
      call computeCoef(nZone, omega, qmuOfZone(:), coef(:))

      ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
      do j = 1, nZone
        call computeA0(nLayerInZone(j), omega, omegai, t(oRowOfZone(j):), &
          h1(oRowOfZone(j):), h2(oRowOfZone(j):), h3(oRowOfZone(j):), h4(oRowOfZone(j):), coef(j), cwork(oRowOfZone(j):))
        call overlapMatrixBlocks(nLayerInZone(j), cwork(oRowOfZone(j):), a0(1, oGridOfZone(j):))

        call computeA2(nLayerInZone(j), h4(oRowOfZone(j):), coef(j), cwork(oRowOfZone(j):))
        call overlapMatrixBlocks(nLayerInZone(j), cwork(oRowOfZone(j):), a2(1, oGridOfZone(j):))
      end do

      ! Initially, no depth cut-off, so set to the index of deepest grid, which is 1.
      cutoffGrid = 1
      ! Clear counter.
      decayCounter = 0
      ! Clear amplitude record.
      recordAmplitude = -1.d0
      !TODO
      ltmp(ii) = maxL  !!!diff

      do l = 0, maxL  ! l-loop
        ! When the counter detecting the decay of amplitude has reached a threshold, stop l-loop for this frequency.
        if (decayCounter > 20) then
          if (ltmp(ii) > l) ltmp(ii) = l
          exit
        end if

        ! Clear the amplitude accumulated for all m's.
        amplitudeAtGrid(:) = 0.d0

        ! Initialize matrices.
        call initComplexMatrix(lda, nGrid, a(:,:))
        call initComplexMatrix(lda, 3, aSource(:,:))

        ! Assemble A matrix from parts that have already been computed.
        call assembleA(nGrid, l, a0(:,:), a2(:,:), a(:,:))

        ! Compute A matrix in layer near source.   TODO: Can't part of a(:,:) be used?
        call computeA(1, omega, omegai, l, t(oRowOfSource:), &
          h1(oRowOfSource:), h2(oRowOfSource:), h3(oRowOfSource:), h4(oRowOfSource:), coef(iZoneOfSource), aaParts(:))

        ! Compute A matrix near source.
        call computeA(2, omega, omegai, l, gt(:), gh1(:), gh2(:), gh3(:), gh4(:), coef(iZoneOfSource), aSourceParts(:))
        call overlapMatrixBlocks(2, aSourceParts(:), aSource(:,:))

        do m = -2, 2  ! m-loop
          if (m == 0 .or. abs(m) > abs(l)) cycle

          call initComplexVector(nGrid, g_or_c)

          ! Computate excitation vector g.
          call computeG(l, m, iLayerOfSource, r0, mt, mu0, coef(iZoneOfSource), aSourceParts(:), aaParts(:), aSource(:,:), &
            gdr, g_or_c(:))

          if (mod(l, 100) == 0) then
            ! Once in a while, compute for all grids to decide the cut-off depth.

            ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
            if ((m == -2) .or. (m == -l)) then
              ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
              call dclisb0(a(:,:), nGrid, 1, lda, g_or_c(:), eps, dr, z, ier)
            else
              ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
              call dcsbsub0(a(:,:), nGrid, 1, lda, g_or_c(:), eps, dr, z, ier)
            end if

            ! Accumulate the absolute values of expansion coefficent c for all m at each grid point.
            !  This is to be used as an estimate of the amplitude at each depth when deciding the cut-off depth.
            amplitudeAtGrid(:) = amplitudeAtGrid(:) + abs(g_or_c(:))

          else
            ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
            if ((m == -2) .or. (m == -l)) then
              ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
              call dclisb(a(:, cutoffGrid:), nGrid - cutoffGrid + 1, 1, lda, iLayerOfSource - cutoffGrid + 1, &
                g_or_c(cutoffGrid:), eps, dr, z, ier)
            else
              ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
              call dcsbsub(a(:, cutoffGrid:), nGrid - cutoffGrid + 1, 1, lda, iLayerOfSource - cutoffGrid + 1, &
                g_or_c(cutoffGrid:), eps, dr, z, ier)
            end if
          end if

          ! Check whether the amplitude has decayed enough to stop the l-loops.
          !  This is checked for the topmost-grid expansion coefficent of each m individually.
          call checkAmplitudeDecay(g_or_c(nGrid), l, lsuf, ratl, recordAmplitude, decayCounter)

        end do  ! m-loop

        ! Decide cut-off depth (at a certain interval of l).
        if (mod(l, 100) == 0) then
          call computeCutoffDepth(nGrid, amplitudeAtGrid(:), ratc, cutoffGrid)
        end if

      end do  ! l-loop
    end do  ! omega-loop

    iimax = int(dble(max(ltmp(1), ltmp(2))) * tlen / lmaxdivf)
  end if  ! option for shallow events


  ! ******************* Computing parameters *******************
  ! Design the number and position of grid points.
  call computeKz(nZone, rminOfZone(:), rmaxOfZone(:), vsvPolynomials(:,:), rmax, iimax, 1, tlen, kzAtZone(:))
  call computeGridRadii(nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), rmin, re, nLayer, nLayerInZone(:), gridRadii(:))
  if (nLayer > maxNLayer) stop 'The number of grid points is too large.'

  ! Compute the first indices in each zone.
  call computeFirstIndices(nZone, nLayerInZone(:), oGridOfZone(:), oRowOfZone(:))

  ! Compute the source position.
  call computeSourcePosition(nLayer, rmaxOfZone(:), rmin, rmax, gridRadii(:), r0, iZoneOfSource, iLayerOfSource)

  ! Design grids for source computations.
  call computeSourceGrid(gridRadii(:), r0, iLayerOfSource, gridRadiiForSource(:))


  ! ******************* Computing the matrix elements *******************
  ! Compute variable values at grid points.
  call computeStructureValues(nZone, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), nLayerInZone(:), &
    gridRadii(:), nValue, valuedRadii(:), rhoValues(:), ecLValues(:), ecNValues(:))
  call computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), &
    gridRadiiForSource(:), rhoValuesForSource(:), ecLValuesForSource(:), ecNValuesForSource(:), mu0)

  ! Compute mass and rigitidy matrices.
  do i = 1, nZone
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), rhoValues(:), 2, 0, 0, gridRadii(oGridOfZone(i):), &
      t(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), 2, 1, 1, gridRadii(oGridOfZone(i):), &
      h1(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), 1, 1, 0, gridRadii(oGridOfZone(i):), &
      h2(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), 0, 0, 0, gridRadii(oGridOfZone(i):), &
      h3(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii(:), ecNValues(:), 0, 0, 0, gridRadii(oGridOfZone(i):), &
      h4(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeLumpedT(nLayerInZone(i), nValue, valuedRadii(:), rhoValues(:), gridRadii(oGridOfZone(i):), work(oRowOfZone(i):))
    call computeAverage(nLayerInZone(i), t(oRowOfZone(i):), work(oRowOfZone(i):), t(oRowOfZone(i):))
    call computeLumpedH(nLayerInZone(i), nValue, valuedRadii(:), ecLValues(:), gridRadii(oGridOfZone(i):), work(oRowOfZone(i):))
    call computeAverage(nLayerInZone(i), h3(oRowOfZone(i):), work(oRowOfZone(i):), h3(oRowOfZone(i):))
    call computeLumpedH(nLayerInZone(i), nValue, valuedRadii(:), ecNValues(:), gridRadii(oGridOfZone(i):), work(oRowOfZone(i):))
    call computeAverage(nLayerInZone(i), h4(oRowOfZone(i):), work(oRowOfZone(i):), h4(oRowOfZone(i):))
  end do

  ! Compute mass and rigitidy matrices near source.
  call computeIntermediateIntegral(2, 3, gridRadiiForSource, rhoValuesForSource, 2, 0, 0, gridRadiiForSource, gt, work)
  call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecLValuesForSource, 2, 1, 1, gridRadiiForSource, gh1, work)
  call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecLValuesForSource, 1, 1, 0, gridRadiiForSource, gh2, work)
  call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecLValuesForSource, 0, 0, 0, gridRadiiForSource, gh3, work)
  call computeIntermediateIntegral(2, 3, gridRadiiForSource, ecNValuesForSource, 0, 0, 0, gridRadiiForSource, gh4, work)
  call computeLumpedT(2, 3, gridRadiiForSource, rhoValuesForSource, gridRadiiForSource, work)
  call computeAverage(2, gt, work, gt)
  call computeLumpedH(2, 3, gridRadiiForSource, ecLValuesForSource, gridRadiiForSource, work)
  call computeAverage(2, gh3, work, gh3)
  call computeLumpedH(2, 3, gridRadiiForSource, ecNValuesForSource, gridRadiiForSource, work)
  call computeAverage(2, gh4, work, gh4)


  !******************** Computing the displacement *********************
  outputindex = 1
  nGrid = nLayer + 1
  ! Find the first index of (iLayer, k', k)-pair corresponding to the layer that the source is in.
  oRowOfSource = 4 * iLayerOfSource - 3

  llog = 0
  do i = imin, imax  ! omega-loop
    omega = 2.d0 * pi * dble(i) / tlen

    ! Initialize matrices.
    call initComplexMatrix(3, nReceiver, u(:,:))
    call initComplexMatrix(lda, nGrid, a0(:,:))
    call initComplexMatrix(lda, nGrid, a2(:,:))
    ! Plm must be cleared for each omega.
    do ir = 1, nReceiver
      call matinit(3, 4, plm(:, :, ir))
    end do

    ! Compute the angular order that is sufficient to compute the slowest phase velocity.
    call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

    ! Compute coefficient related to attenuation.
    call computeCoef(nZone, omega, qmuOfZone(:), coef(:))

    ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
    do j = 1, nZone
      call computeA0(nLayerInZone(j), omega, omegai, t(oRowOfZone(j):), &
        h1(oRowOfZone(j):), h2(oRowOfZone(j):), h3(oRowOfZone(j):), h4(oRowOfZone(j):), coef(j), cwork(oRowOfZone(j):))
      call overlapMatrixBlocks(nLayerInZone(j), cwork(oRowOfZone(j):), a0(:, oGridOfZone(j):))

      call computeA2(nLayerInZone(j), h4(oRowOfZone(j):), coef(j), cwork(oRowOfZone(j):))
      call overlapMatrixBlocks(nLayerInZone(j), cwork(oRowOfZone(j):), a2(:, oGridOfZone(j):))
    end do

    ! Initially, no depth cut-off, so set to the index of deepest grid, which is 1.
    cutoffGrid = 1
    ! Clear counter.
    decayCounter = 0
    ! Clear amplitude record.
    recordAmplitude = -1.d0
    !TODO erase
    llog = maxL

    do l = 0, maxL  ! l-loop
      ! When the counter detecting the decay of amplitude has reached a threshold, stop l-loop for this frequency.
      if (decayCounter > 20) then
        llog = min(llog, l)
        cycle
      end if

      ! Clear the amplitude accumulated for all m's.
      amplitudeAtGrid(:) = 0.d0

      ! Compute trial functions.
      do ir = 1, nReceiver
        call computeTrialFunctionValues(l, theta(ir), phi(ir), plm(:, :, ir), trialFunctionValues(:, :, ir))
      end do

      ! Initialize matrices.
      call initComplexMatrix(lda, nGrid, a(:,:))
      call initComplexMatrix(lda, 3, aSource(:,:))

      ! Assemble A matrix from parts that have already been computed.
      call assembleA(nGrid, l, a0(:,:), a2(:,:), a(:,:))

      ! Compute A matrix in layer near source.   TODO: Can't part of a(:,:) be used?
      call computeA(1, omega, omegai, l, t(oRowOfSource:), &
        h1(oRowOfSource:), h2(oRowOfSource:), h3(oRowOfSource:), h4(oRowOfSource:), coef(iZoneOfSource), aaParts(:))

      ! Compute A matrix near source.
      call computeA(2, omega, omegai, l, gt(:), gh1(:), gh2(:), gh3(:), gh4(:), coef(iZoneOfSource), aSourceParts(:))
      call overlapMatrixBlocks(2, aSourceParts(:), aSource(:,:))

      do m = -2, 2  ! m-loop
        if (m == 0 .or. abs(m) > abs(l)) cycle

        call initComplexVector(nGrid, g_or_c)

        ! Computate excitation vector g.
        call computeG(l, m, iLayerOfSource, r0, mt, mu0, coef(iZoneOfSource), aSourceParts(:), aaParts(:), aSource(:,:), &
          gdr(:), g_or_c(:))

        if (mod(l, 100) == 0) then
          ! Once in a while, compute for all grids to decide the cut-off depth.

          ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
          if (m == -2 .or. m == -l) then
            ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
            call dclisb0(a(:,:), nGrid, 1, lda, g_or_c(:), eps, dr, z, ier)
          else
            ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
            call dcsbsub0(a(:,:), nGrid, 1, lda, g_or_c(:), eps, dr, z, ier)
          end if

          ! Accumulate the absolute values of expansion coefficent c for all m's at each grid point.
          !  This is to be used as an estimate of the amplitude at each depth when deciding the cut-off depth.
          amplitudeAtGrid(:) = amplitudeAtGrid(:) + abs(g_or_c(:))

        else
          ! Otherwise, compute for just the grids above the cut-off depth.

          ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
          if (m == -2 .or. m == -l) then
            ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
            call dclisb(a(:, cutoffGrid:), nGrid - cutoffGrid + 1, 1, lda, iLayerOfSource - cutoffGrid + 1, &
              g_or_c(cutoffGrid:), eps, dr, z, ier)
          else
            ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
            call dcsbsub(a(:, cutoffGrid:), nGrid - cutoffGrid + 1, 1, lda, iLayerOfSource - cutoffGrid + 1, &
              g_or_c(cutoffGrid:), eps, dr, z, ier)
          end if
        end if

        ! Check whether the amplitude has decayed enough to stop the l-loops.
        !  This is checked for the topmost-grid expansion coefficent of each m individually.
        call checkAmplitudeDecay(g_or_c(nGrid), l, lsuf, ratl, recordAmplitude, decayCounter)

        ! Accumulate u.
        do ir = 1, nReceiver
          call computeU(g_or_c(nGrid), l, trialFunctionValues(:, m, ir), u(:, ir))
        end do

      end do  ! m-loop

      ! Decide cut-off depth (at a certain interval of l).
      if (mod(l, 100) == 0) then
        call computeCutoffDepth(nGrid, amplitudeAtGrid(:), ratc, cutoffGrid)
      end if

    end do  ! l-loop

    ! Store results.
    outputi(outputindex) = i
    do ir = 1, nReceiver
      outputu(:, ir, outputindex) = u(:, ir)
    end do

    ! ************************** Files Handling **************************
    ! Write to file when the output interval is reached, or when this is the last omega.
    if (outputindex >= outputinterval .or. i == imax) then
      write(*,*) "kakikomimasu"
      if (spcFormat == 0) then
        do ir = 1, nReceiver
          open(unit=10, file=output(ir), position='append', status='unknown', &
            form='unformatted', access='stream', convert='big_endian')
          do mpii = 1, outputindex
            write(10) outputi(mpii), dble(outputu(1, ir, mpii)), imag(outputu(1, ir, mpii))
            write(10) dble(outputu(2, ir, mpii)), imag(outputu(2, ir, mpii))
            write(10) dble(outputu(3, ir, mpii)), imag(outputu(3, ir, mpii))
          end do
          close(10)
        end do
      else if (spcFormat == 1) then
        do ir = 1, nReceiver
          open(unit=10, file=output(ir), position='append', status='unknown')
          do mpii = 1, outputindex
            write(10,*) outputi(mpii), dble(outputu(1, ir, mpii)), imag(outputu(1, ir, mpii))
            write(10,*) dble(outputu(2, ir, mpii)), imag(outputu(2, ir, mpii))
            write(10,*) dble(outputu(3, ir, mpii)), imag(outputu(3, ir, mpii))
          end do
          close(10)
        end do
      else
        write(*,*) "WARNING: set spcFormat 0 or 1"
      end if
      outputindex = 0
    end if

    if (ilog == 1) then
      open(unit=11, file='llog.log', position='append', status='old')
      write(11,*) i, llog, nLayer
      close(11)
    end if

    outputindex = outputindex + 1

  end do  ! omega-loop

  write(*,*) "Ivalice looks to the horizon"

  stop
end program tish


