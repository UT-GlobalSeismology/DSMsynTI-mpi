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
  integer, parameter :: maxNGrid = 88300  ! Maximum number of grid points.
  integer, parameter :: maxNZone = 15  ! Maximum number of zones.
  integer, parameter :: maxNReceiver = 1500  ! Maximum number of receivers.
  integer, parameter :: maxL = 80000  ! Maximum of angular order to loop for.
  real(8), parameter :: lmaxdivf = 2.d4
  real(8), parameter :: shallowdepth = 100.d0
  integer, parameter :: spcFormat = 1  ! Format of output spc file (0:binary, 1:ascii).
  integer, parameter :: ilog = 1

  !----------------------------<<variables>>----------------------------
  ! Variables for the structure
  integer :: nZone  ! Number of zones.
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled.
  real(8) :: rminOfZone(maxNZone), rmaxOfZone(maxNZone)  ! Minimum and maximum radii of each zone.
  real(8) :: rhoPolynomials(4, maxNZone)  ! Rho of each zone (coefficients of cubic function).
  real(8) :: vsvPolynomials(4, maxNZone)  ! Vsv of each zone (coefficients of cubic function).
  real(8) :: vshPolynomials(4, maxNZone)  ! Vsh of each zone (coefficients of cubic function).
  real(8) :: qmuOfZone(maxNZone)  ! Qmu of each zone.
  integer :: i

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
  real(8) :: omegaI  ! Imaginary part of angular frequency for artificial damping. (See section 5.1 of Geller & Ohminato 1994.)
  integer :: imin, imax  ! Index of minimum and maximum frequency.
  integer :: iFreq, iCount
  integer :: ltmp(2), imaxFixed

  ! Variables for grid spacing and cut-off
  real(8) :: kzAtZone(maxNZone)  ! Vertical wavenumber k_z at each zone. (See section 3.2 of Kawai et al. 2006.)
  real(8) :: re  ! Desired relative error due to vertical gridding. (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
  real(8) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8) :: amplitudeAtGrid(maxNGrid)  ! Estimate of the amplitude at each grid point, used for vertical grid cut-off.
  integer :: cutoffGrid  ! Index of grid at cut-off depth.
  integer :: lsuf  ! Accuracy threshold of angular order. (Corresponds to l_d; see eq. 29 of Kawai et al. 2006.)
  real(8) :: recordAmplitude    ! Maximum amplitude encountered, used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: llog

  ! Variables for the vertical grid
  integer :: nGrid  ! Total number of grid points.
  real(8) :: gridRadii(maxNGrid)  ! Radii of each grid point.
  integer :: nLayerInZone(maxNZone)  ! Number of layers in each zone.
  integer :: oGridOfZone(maxNZone)  ! Index of the first grid point in each zone.
  real(8) :: gridRadiiForSource(3)  ! Radii to use for source-related computations.
  integer :: iZoneOfSource  ! Which zone the source is in.
  integer :: iLayerOfSource  ! Index of layer that the source is in.

  ! Variables for the values
  integer :: nValue  ! Total number of values.
  real(8) :: valuedRadii(maxNGrid + maxNZone - 1)  ! Radii corresponding to each variable value.
  integer :: oValueOfZone(maxNZone)  ! Index of the first value in each zone.
  real(8) :: rhoValues(maxNGrid + maxNZone - 1)  ! Rho at each grid point (with 2 values at boundaries).
  real(8) :: ecLValues(maxNGrid + maxNZone - 1)  ! L at each grid point (with 2 values at boundaries).
  real(8) :: ecNValues(maxNGrid + maxNZone - 1)  ! N at each grid point (with 2 values at boundaries).
  real(8) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)  ! Rho, L, and N at each source-related grid.
  complex(8) :: coef(maxNZone)

  ! Variables for the trial function
  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: largeL2  ! L^2 = l(l+1).
  real(8) :: plm(3, 0:3, maxNReceiver)  ! Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  !::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8) :: trialFunctionValues(3, -2:2, maxNReceiver)  ! Values of trial function at each receiver, computed for each l.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Arguments: component (1:3), m (-2:2), iReceiver.

  ! Variables for the stack points
  integer :: oRowOfZone(maxNZone)  ! Index of the first row in the vector of (iLayer, k', k)-pairs in each zone.
  integer :: oRowOfSource  ! Index of the first row in the vector of (iLayer, k', k)-pairs for the layer with the source.

  ! Variables for the matrix elements
  real(8) :: t(4 * maxNGrid - 4)
  real(8) :: h1(4 * maxNGrid - 4), h2(4 * maxNGrid - 4), h3(4 * maxNGrid - 4), h4(4 * maxNGrid - 4)
  real(8) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  complex(8) :: a0(2, maxNGrid), a2(2, maxNGrid)
  complex(8) :: a(2, maxNGrid)
  complex(8) :: aaParts(4), aSourceParts(8), aSource(2, 3)
  complex(8) :: g_or_c(maxNGrid)  ! This holds either the vector g or c, depending on where in the code it is (CAUTION!).
  complex(8) :: u(3, maxNReceiver)

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Other variables
  real(8) :: work(4 * maxNGrid - 4)  ! Working array for matrix computations.
  complex(8) :: cwork(4 * maxNGrid - 4)  ! Working array for matrix computations.
  integer :: ier  ! Error code from subroutine solving linear equations.
  complex(8) :: dr(maxNGrid), z(maxNGrid), gdr(3)  ! Working arrays used when solving linear equations.

  ! Constants
  integer :: lda = 2
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
  call readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegaI, imin, imax, &
    nZone, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, qmuOfZone, &
    r0, eqlat, eqlon, mt, nReceiver, theta, phi, lat, lon, output)

  ! Find the amount of memory that is written in 1 omega step.
  !  For each omega and receiver, 3 complex numbers (16 B each) are written. 1 B = 0.000001 MB.
  memoryPerOmega = 3 * 16 * nReceiver * 0.000001
  ! Find how many omegas can be written within outputMemory.
  outputInterval = int(outputMemory / memoryPerOmega)
  ! Allocate arrays to store output.
  allocate(outputi(outputInterval))
  allocate(outputu(3, nReceiver, outputInterval))

  ! --- computing the required parameters ---
  rmin = rminOfZone(1)
  rmax = rmaxOfZone(nZone)

  write(*, *) 'theta:', theta(1:nReceiver)  !TODO erase
  write(*, *) 'phi:', phi(1:nReceiver)  !TODO erase

  if (r0 < rmin .or. r0 > rmax) then
    stop 'Location of the source is improper.'
  end if

  if (imin == 0) imin = 1
  ! Decide which omega to use when deciding grid spacing. Usually, this is just the upper limit of omega range.
  imaxFixed = imax


  ! ************************** Files handling **************************
  if (spcFormat == 0) then
    do ir = 1, nReceiver
      open(unit = 11, file = output(ir), status = 'unknown', &
        form = 'unformatted', access = 'stream', convert = 'big_endian')
      write(11) tlen
      write(11) np, 1, 3
      write(11) omegaI, lat(ir), lon(ir)
      write(11) eqlat, eqlon, r0
      close(11)
    end do
  else if (spcFormat == 1) then
    do ir = 1, nReceiver
      open(unit = 11, file = output(ir), status = 'unknown')
      write(11, *) tlen
      write(11, *) np, 1, 3
      write(11, *) omegaI, lat(ir), lon(ir)
      write(11, *) eqlat, eqlon, r0
      close(11)
    end do
  else
    write(*, *) 'WARNING:(tish.f)  set spcFormat 0 or 1'
  end if

  if (ilog == 1) then
    open(unit = 11, file = 'llog.log', status = 'unknown')
    write(11, *) 0
    close(11)
  end if


  ! ************************** Option for shallow events **************************
  ! Here, we find the maximum angular order needed for our frequency range. (See fig. 7 of Kawai et al. 2006.)
  if ((rmax - r0) < shallowdepth) then

    write(*, *) 'Shallow event!'  !TODO erase

    ! ******************* Computing parameters *******************
    ! Set a large value so that we can compute using fine grids for this process.
    imaxFixed = int(tlen * 2.d0)  !!!diff

    ! Design the number and position of grid points.
    call computeKz(nZone, rminOfZone(:), rmaxOfZone(:), vsvPolynomials(:,:), rmax, imaxFixed, 1, tlen, kzAtZone(:))
    write(*, *) 'imaxFixed, tlen:', imaxFixed, tlen  !TODO erase
    call computeGridRadii(nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), rmin, re, nGrid, nLayerInZone(:), gridRadii(:))
    if (nGrid > maxNGrid) stop 'The number of grid points is too large.'
    write(*, *) 'nGrid:', nGrid  !TODO erase
    write(*, *) 'nLayerInZone:', nLayerInZone(1:nZone)  !TODO erase

    ! Compute the first indices in each zone.
    call computeFirstIndices(nZone, nLayerInZone(:), oGridOfZone(:), oValueOfZone(:), oRowOfZone(:))

    ! Compute the source position.
    call computeSourcePosition(nGrid, rmaxOfZone(:), rmin, rmax, gridRadii(:), r0, iZoneOfSource, iLayerOfSource)
    write(*, *) 'iLayerOfSource:', iLayerOfSource  !TODO erase

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
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), rhoValues(oValueOfZone(i):), 2, 0, 0, &
        t(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), 2, 1, 1, &
        h1(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), 1, 1, 0, &
        h2(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), 0, 0, 0, &
        h3(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecNValues(oValueOfZone(i):), 0, 0, 0, &
        h4(oRowOfZone(i):), work(oRowOfZone(i):))
      call computeLumpedT(nLayerInZone(i), valuedRadii(oValueOfZone(i):), rhoValues(oValueOfZone(i):), work(oRowOfZone(i):))
      call computeAverage(nLayerInZone(i), t(oRowOfZone(i):), work(oRowOfZone(i):), t(oRowOfZone(i):))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), work(oRowOfZone(i):))
      call computeAverage(nLayerInZone(i), h3(oRowOfZone(i):), work(oRowOfZone(i):), h3(oRowOfZone(i):))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecNValues(oValueOfZone(i):), work(oRowOfZone(i):))
      call computeAverage(nLayerInZone(i), h4(oRowOfZone(i):), work(oRowOfZone(i):), h4(oRowOfZone(i):))
    end do

    ! Compute mass and rigitidy matrices near source.
    call computeIntermediateIntegral(2, gridRadiiForSource, rhoValuesForSource, 2, 0, 0, gt, work)
    call computeIntermediateIntegral(2, gridRadiiForSource, ecLValuesForSource, 2, 1, 1, gh1, work)
    call computeIntermediateIntegral(2, gridRadiiForSource, ecLValuesForSource, 1, 1, 0, gh2, work)
    call computeIntermediateIntegral(2, gridRadiiForSource, ecLValuesForSource, 0, 0, 0, gh3, work)
    call computeIntermediateIntegral(2, gridRadiiForSource, ecNValuesForSource, 0, 0, 0, gh4, work)
    call computeLumpedT(2, gridRadiiForSource, rhoValuesForSource, work)
    call computeAverage(2, gt, work, gt)
    call computeLumpedH(2, gridRadiiForSource, ecLValuesForSource, work)
    call computeAverage(2, gh3, work, gh3)
    call computeLumpedH(2, gridRadiiForSource, ecNValuesForSource, work)
    call computeAverage(2, gh4, work, gh4)


    !******************** Computing the expansion coefficients *********************
    ! Find the first index of (iLayer, k', k)-pair corresponding to the layer that the source is in.
    oRowOfSource = 4 * iLayerOfSource - 3

    ! Find the maximum angular order needed for the lowest and highest frequencies. (See fig. 7 of Kawai et al. 2006.)
    do iCount = 1, 2  ! omega-loop
      if (iCount == 1) then  !!!diff
        iFreq = imin
      else
        iFreq = imax
      end if
      omega = 2.d0 * pi * dble(iFreq) / tlen

      ! Initialize matrices.
      a0(:lda, :nGrid) = dcmplx(0.0d0, 0.0d0)
      a2(:lda, :nGrid) = dcmplx(0.0d0, 0.0d0)

      ! Compute the angular order that is sufficient to compute the slowest phase velocity.
      call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

      ! Compute coefficient related to attenuation.
      call computeCoef(nZone, omega, qmuOfZone(:), coef(:))

      ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
      do i = 1, nZone
        call computeA0(nLayerInZone(i), omega, omegaI, t(oRowOfZone(i):), &
          h1(oRowOfZone(i):), h2(oRowOfZone(i):), h3(oRowOfZone(i):), h4(oRowOfZone(i):), coef(i), cwork(oRowOfZone(i):))
        call overlapMatrixBlocks(nLayerInZone(i), cwork(oRowOfZone(i):), a0(:, oGridOfZone(i):))

        call computeA2(nLayerInZone(i), h4(oRowOfZone(i):), coef(i), cwork(oRowOfZone(i):))
        call overlapMatrixBlocks(nLayerInZone(i), cwork(oRowOfZone(i):), a2(:, oGridOfZone(i):))
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

        ! L^2. (See the part after eq. 12 of Kawai et al. 2006.)
        ! NOTE that integers are casted with dble() before multiplying, because the product can exceed the size of integer(4).
        largeL2 = dble(l) * dble(l + 1)

        ! Initialize matrices.
        a(:lda, :nGrid) = dcmplx(0.0d0, 0.0d0)
        aSource(:lda, :) = dcmplx(0.0d0, 0.0d0)
        ! Clear the amplitude accumulated for all m's.
        amplitudeAtGrid(:nGrid) = 0.d0

        ! Assemble A matrix from parts that have already been computed.
        call assembleA(nGrid, largeL2, a0(:,:), a2(:,:), a(:,:))

        ! Compute A matrix in layer near source.   TODO: Can't part of a(:,:) be used?
        call computeA(1, omega, omegaI, largeL2, t(oRowOfSource:), &
          h1(oRowOfSource:), h2(oRowOfSource:), h3(oRowOfSource:), h4(oRowOfSource:), coef(iZoneOfSource), aaParts(:))

        ! Compute A matrix near source.
        call computeA(2, omega, omegaI, largeL2, gt(:), gh1(:), gh2(:), gh3(:), gh4(:), coef(iZoneOfSource), aSourceParts(:))
        call overlapMatrixBlocks(2, aSourceParts(:), aSource(:,:))

        do m = -2, 2  ! m-loop
          if (m == 0 .or. abs(m) > abs(l)) cycle

          ! Initialize vector.
          g_or_c(:nGrid) = dcmplx(0.d0, 0.d0)

          ! Computate excitation vector g.
          call computeG(l, m, iLayerOfSource, r0, mt, mu0, coef(iZoneOfSource), aSourceParts(:), aaParts(:), aSource(:,:), &
            gdr(:), g_or_c(:))

          if (mod(l, 100) == 0) then
            ! Once in a while, compute for all grids to decide the cut-off depth.
            ! CAUTION: In this case, all values of g_or_c(:) are computed.

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
            amplitudeAtGrid(1:nGrid) = amplitudeAtGrid(1:nGrid) + abs(g_or_c(1:nGrid))

          else
            ! Otherwise, compute for just the grids above the cut-off depth.
            ! CAUTION: In this case, only g_or_c(nGrid) is computed. Other values of g_or_c(:nGrid-1) still hold values of g!!!

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

        end do  ! m-loop

        ! Decide cut-off depth (at a certain interval of l).
        if (mod(l, 100) == 0) then
          call computeCutoffDepth(nGrid, amplitudeAtGrid(:), ratc, cutoffGrid)
        end if

      end do  ! l-loop

      ! Register the final l (or maxL instead of maxL-1 when all loops are completed).
      ltmp(iCount) = min(l, maxL)
      write(*, *) 'ltmp:', ltmp(iCount), iCount  !TODO erase

    end do  ! omega-loop

    imaxFixed = int(dble(max(ltmp(1), ltmp(2))) * tlen / lmaxdivf)
  end if  ! option for shallow events


  ! ******************* Computing parameters *******************
  ! Design the number and position of grid points.
  call computeKz(nZone, rminOfZone(:), rmaxOfZone(:), vsvPolynomials(:,:), rmax, imaxFixed, 1, tlen, kzAtZone(:))
  write(*, *) 'imaxFixed, tlen:', imaxFixed, tlen  !TODO erase

  call computeGridRadii(nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), rmin, re, nGrid, nLayerInZone(:), gridRadii(:))
  if (nGrid > maxNGrid) stop 'The number of grid points is too large.'
  write(*, *) 'nGrid:', nGrid  !TODO erase
  write(*, *) 'nLayerInZone:', nLayerInZone(1:nZone)  !TODO erase

  ! Compute the first indices in each zone.
  call computeFirstIndices(nZone, nLayerInZone(:), oGridOfZone(:), oValueOfZone(:), oRowOfZone(:))

  ! Compute the source position.
  call computeSourcePosition(nGrid, rmaxOfZone(:), rmin, rmax, gridRadii(:), r0, iZoneOfSource, iLayerOfSource)

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
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), rhoValues(oValueOfZone(i):), 2, 0, 0, &
      t(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), 2, 1, 1, &
      h1(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), 1, 1, 0, &
      h2(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), 0, 0, 0, &
      h3(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecNValues(oValueOfZone(i):), 0, 0, 0, &
      h4(oRowOfZone(i):), work(oRowOfZone(i):))
    call computeLumpedT(nLayerInZone(i), valuedRadii(oValueOfZone(i):), rhoValues(oValueOfZone(i):), work(oRowOfZone(i):))
    call computeAverage(nLayerInZone(i), t(oRowOfZone(i):), work(oRowOfZone(i):), t(oRowOfZone(i):))
    call computeLumpedH(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecLValues(oValueOfZone(i):), work(oRowOfZone(i):))
    call computeAverage(nLayerInZone(i), h3(oRowOfZone(i):), work(oRowOfZone(i):), h3(oRowOfZone(i):))
    call computeLumpedH(nLayerInZone(i), valuedRadii(oValueOfZone(i):), ecNValues(oValueOfZone(i):), work(oRowOfZone(i):))
    call computeAverage(nLayerInZone(i), h4(oRowOfZone(i):), work(oRowOfZone(i):), h4(oRowOfZone(i):))
  end do
  write(*, *) 't:', t(1:4), t((nGrid-1)*4-3:(nGrid-1)*4)  !TODO erase
  write(*, *) 'h1:', h1(1:4), h1((nGrid-1)*4-3:(nGrid-1)*4)  !TODO erase
  write(*, *) 'h2:', h2(1:4), h2((nGrid-1)*4-3:(nGrid-1)*4)  !TODO erase
  write(*, *) 'h3:', h3(1:4), h3((nGrid-1)*4-3:(nGrid-1)*4)  !TODO erase
  write(*, *) 'h4:', h4(1:4), h4((nGrid-1)*4-3:(nGrid-1)*4)  !TODO erase

  ! Compute mass and rigitidy matrices near source.
  call computeIntermediateIntegral(2, gridRadiiForSource, rhoValuesForSource, 2, 0, 0, gt, work)
  call computeIntermediateIntegral(2, gridRadiiForSource, ecLValuesForSource, 2, 1, 1, gh1, work)
  call computeIntermediateIntegral(2, gridRadiiForSource, ecLValuesForSource, 1, 1, 0, gh2, work)
  call computeIntermediateIntegral(2, gridRadiiForSource, ecLValuesForSource, 0, 0, 0, gh3, work)
  call computeIntermediateIntegral(2, gridRadiiForSource, ecNValuesForSource, 0, 0, 0, gh4, work)
  call computeLumpedT(2, gridRadiiForSource, rhoValuesForSource, work)
  call computeAverage(2, gt, work, gt)
  call computeLumpedH(2, gridRadiiForSource, ecLValuesForSource, work)
  call computeAverage(2, gh3, work, gh3)
  call computeLumpedH(2, gridRadiiForSource, ecNValuesForSource, work)
  call computeAverage(2, gh4, work, gh4)
  write(*, *) 'gt:', gt  !TODO erase
  write(*, *) 'gh1:', gh1  !TODO erase
  write(*, *) 'gh2:', gh2  !TODO erase
  write(*, *) 'gh3:', gh3  !TODO erase
  write(*, *) 'gh4:', gh4  !TODO erase
  write(*, *) '--------'  !TODO erase

  !******************** Computing the displacement *********************
  outputCounter = 1
  ! Find the first index of (iLayer, k', k)-pair corresponding to the layer that the source is in.
  oRowOfSource = 4 * iLayerOfSource - 3

  do iFreq = imin, imax  ! omega-loop
    omega = 2.d0 * pi * dble(iFreq) / tlen

    ! Initialize matrices.
    a0(:lda, :nGrid) = dcmplx(0.0d0, 0.0d0)
    a2(:lda, :nGrid) = dcmplx(0.0d0, 0.0d0)
    u(:, :nReceiver) = dcmplx(0.0d0, 0.0d0)
    ! Plm must be cleared for each omega.
    plm(:, :, :nReceiver) = 0.d0

    ! Compute the angular order that is sufficient to compute the slowest phase velocity.
    call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

    ! Compute coefficient related to attenuation.
    call computeCoef(nZone, omega, qmuOfZone(:), coef(:))

    ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
    do i = 1, nZone
      call computeA0(nLayerInZone(i), omega, omegaI, t(oRowOfZone(i):), &
        h1(oRowOfZone(i):), h2(oRowOfZone(i):), h3(oRowOfZone(i):), h4(oRowOfZone(i):), coef(i), cwork(oRowOfZone(i):))
      call overlapMatrixBlocks(nLayerInZone(i), cwork(oRowOfZone(i):), a0(:, oGridOfZone(i):))

      call computeA2(nLayerInZone(i), h4(oRowOfZone(i):), coef(i), cwork(oRowOfZone(i):))
      call overlapMatrixBlocks(nLayerInZone(i), cwork(oRowOfZone(i):), a2(:, oGridOfZone(i):))
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

      ! L^2. (See the part after eq. 12 of Kawai et al. 2006.)
      ! NOTE that integers are casted with dble() before multiplying, because the product can exceed the size of integer(4).
      largeL2 = dble(l) * dble(l + 1)

      ! Initialize matrices.
      a(:lda, :nGrid) = dcmplx(0.0d0, 0.0d0)
      aSource(:lda, :) = dcmplx(0.0d0, 0.0d0)
      ! Clear the amplitude accumulated for all m's.
      amplitudeAtGrid(:nGrid) = 0.d0

      ! Compute trial functions.
      do ir = 1, nReceiver
        call computeTrialFunctionValues(l, theta(ir), phi(ir), plm(:, :, ir), trialFunctionValues(:, :, ir))
      end do

      ! Assemble A matrix from parts that have already been computed.
      call assembleA(nGrid, largeL2, a0(:,:), a2(:,:), a(:,:))

      ! Compute A matrix in layer near source.   TODO: Can't part of a(:,:) be used?
      call computeA(1, omega, omegaI, largeL2, t(oRowOfSource:), &
        h1(oRowOfSource:), h2(oRowOfSource:), h3(oRowOfSource:), h4(oRowOfSource:), coef(iZoneOfSource), aaParts(:))

      ! Compute A matrix near source.
      call computeA(2, omega, omegaI, largeL2, gt(:), gh1(:), gh2(:), gh3(:), gh4(:), coef(iZoneOfSource), aSourceParts(:))
      call overlapMatrixBlocks(2, aSourceParts(:), aSource(:,:))

      do m = -2, 2  ! m-loop
        if (m == 0 .or. abs(m) > abs(l)) cycle

        ! Initialize vector.
        g_or_c(:nGrid) = dcmplx(0.d0, 0.d0)

        ! Computate excitation vector g.
        call computeG(l, m, iLayerOfSource, r0, mt, mu0, coef(iZoneOfSource), aSourceParts(:), aaParts(:), aSource(:,:), &
          gdr(:), g_or_c(:))

        if (mod(l, 100) == 0) then
          ! Once in a while, compute for all grids to decide the cut-off depth.
          ! CAUTION: In this case, all values of g_or_c(:) are computed.

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
          amplitudeAtGrid(1:nGrid) = amplitudeAtGrid(1:nGrid) + abs(g_or_c(1:nGrid))

        else
          ! Otherwise, compute for just the grids above the cut-off depth.
          ! CAUTION: In this case, only g_or_c(nGrid) is computed. Other values of g_or_c(:nGrid-1) still hold values of g!!!

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
          call computeU(g_or_c(nGrid), largeL2, trialFunctionValues(:, m, ir), u(:, ir))
        end do

      end do  ! m-loop

      ! Decide cut-off depth (at a certain interval of l).
      if (mod(l, 100) == 0) then
        call computeCutoffDepth(nGrid, amplitudeAtGrid(:), ratc, cutoffGrid)
      end if

    end do  ! l-loop

    ! Register the final l (or maxL instead of maxL-1 when all loops are completed).
    llog = min(l, maxL)

    ! Store results.
    outputi(outputCounter) = iFreq
    do ir = 1, nReceiver
      outputu(:, ir, outputCounter) = u(:, ir)
    end do

    ! ************************** Files Handling **************************
    ! Write to file when the output interval is reached, or when this is the last omega.
    if (outputCounter >= outputInterval .or. iFreq == imax) then
      write(*,*) "kakikomimasu"
      if (spcFormat == 0) then
        do ir = 1, nReceiver
          open(unit=10, file=output(ir), position='append', status='unknown', &
            form='unformatted', access='stream', convert='big_endian')
          do iOut = 1, outputCounter
            write(10) outputi(iOut), dble(outputu(1, ir, iOut)), imag(outputu(1, ir, iOut))
            write(10) dble(outputu(2, ir, iOut)), imag(outputu(2, ir, iOut))
            write(10) dble(outputu(3, ir, iOut)), imag(outputu(3, ir, iOut))
          end do
          close(10)
        end do
      else if (spcFormat == 1) then
        do ir = 1, nReceiver
          open(unit=10, file=output(ir), position='append', status='unknown')
          do iOut = 1, outputCounter
            write(10,*) outputi(iOut), dble(outputu(1, ir, iOut)), imag(outputu(1, ir, iOut))
            write(10,*) dble(outputu(2, ir, iOut)), imag(outputu(2, ir, iOut))
            write(10,*) dble(outputu(3, ir, iOut)), imag(outputu(3, ir, iOut))
          end do
          close(10)
        end do
      else
        write(*,*) "WARNING: set spcFormat 0 or 1"
      end if
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
end program tish


