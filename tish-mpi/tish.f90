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
  integer, parameter :: maxlmax = 80000
  integer, parameter :: ilog = 0
  real(8), parameter :: lmaxdivf = 2.d4
  real(8), parameter :: shallowdepth = 100.d0
  integer, parameter :: spcform = 0  ! 0:binary, 1:ascii

  !----------------------------<<variables>>----------------------------
  ! Variables for the trial function
  integer :: nLayer, nLayerInZone(maxNZone)  ! Number of layers (total, and in each zone).
  real(8) :: gridRadii(maxNLayer + maxNZone + 1)  ! Radii of each grid point.
  real(8) :: gridRadiiForSource(3)  ! Radii to use for source-related computations.
  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: plm(3, 0:3, maxNReceiver)  ! Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  !::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8) :: trialFunctionValues(3, -2:2, maxNReceiver)  ! Values of trial function at each receiver, computed for each l.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Arguments: component (1:3), m (-2:2), iReceiver.

  ! Variables for the structure
  integer :: nZone  ! Number of zones.
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled.
  real(8) :: rminOfZone(maxNZone), rmaxOfZone(maxNZone)  ! Minimum and maximum radii of each zone.
  real(8) :: rhoPolynomials(4, maxNZone)  ! Rho of each zone (coefficients of cubic function).
  real(8) :: vsvPolynomials(4, maxNZone)  ! Vsv of each zone (coefficients of cubic function).
  real(8) :: vshPolynomials(4, maxNZone)  ! Vsh of each zone (coefficients of cubic function).
  real(8) :: qmuOfZone(maxNZone)  ! Qmu of each zone.
  integer :: nValue  ! Total number of values.
  real(8) :: valuedRadii(maxNLayer + 2 * maxNZone + 1)  ! Radii corresponding to each variable value.
  real(8) :: rhoValues(maxNLayer + 2 * maxNZone + 1)  ! Rho at each grid point (with 2 values at boundaries).
  real(8) :: ecLValues(maxNLayer + 2 * maxNZone + 1)  ! L at each grid point (with 2 values at boundaries).
  real(8) :: ecNValues(maxNLayer + 2 * maxNZone + 1)  ! N at each grid point (with 2 values at boundaries).
  real(8) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)  ! Rho, L, and N at each source-related grid.
  complex(8) :: coef(maxNZone)

  ! Variables for the periodic range
  real(8) :: tlen  ! Time length.
  integer :: np  ! Number of points in frequency domain.
  real(8) :: omega, omegai
  integer :: imin, imax  ! Index of minimum and maximum frequency.
  complex(8) :: u(3,maxNReceiver)

  ! Variables for the source
  integer :: ns
  real(8) :: r0, mt(3, 3), eqlat, eqlon, mu0
  integer :: iZoneOfSource  ! Which zone the source is in.
  real(8) :: qLayerOfSource  ! A double-value index of source position in its zone.
  !:::::::::::::::::::::::::::: (0 at bottom of zone, nLayerOfZone(iZone) at top of zone.)

  ! Variables for the receivers
  integer :: nReceiver  ! Number of receivers.
  real(8) :: theta(maxNReceiver), phi(maxNReceiver)
  real(8) :: lat(maxNReceiver), lon(maxNReceiver)
  integer :: ir

  ! Variables for the matrix elements
  complex(8) :: a0(2, maxNLayer+1), a2(2, maxNLayer+1)
  complex(8) :: a(2, maxNLayer+1)
  real(8) :: t(4 * maxNLayer)
  real(8) :: h1(4 * maxNLayer), h2(4 * maxNLayer)
  real(8) :: h3(4 * maxNLayer), h4(4 * maxNLayer)
  real(8) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  complex(8) :: aa(4), ga(8), ga2(2, 3), gdr(3)
  complex(8) :: g(maxNLayer + 1)

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Variables for grid spacing
  real(8) :: tmpr(maxNLayer + 1)
  real(8) :: kz(maxNZone)  ! Vertical wavenumber k_z. (See section 3.2 of Kawai et al. 2006.)
  real(8) :: re  ! Desired relative error due to vertical gridding. (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
  real(8) :: ratc  ! Amplitude ratio for vertical grid cut-off.
  real(8) :: ratl  ! Amplitude ratio for angular order cut-off.
  real(8) :: maxamp
  integer :: kc, lsuf, ismall, llog

  ! Variables for the stack points
  integer :: isp(maxNZone), jsp(maxNZone), ins

  ! Other variables
  integer :: i, j, ii, jj, nn, ier
  real(8) :: work(4 * maxNLayer)
  complex(8) :: dr(maxNLayer + 1), z(maxNLayer + 1)
  complex(8) :: cwork(4 * maxNLayer)
  integer :: ltmp(2), iimax

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


  ! ************************** Files handling **************************
  if (spcform == 0) then
    do ir = 1, nReceiver
      open(unit = 11, file = output(ir), status = 'unknown', &
        form = 'unformatted', access = 'stream', convert = 'big_endian')
      write(11) tlen
      write(11) np, 1, 3
      write(11) omegai, lat(ir), lon(ir)
      write(11) eqlat, eqlon, r0
      close(11)
    end do
  else if (spcform == 1) then
    do ir = 1, nReceiver
      open(unit = 11, file = output(ir), status = 'unknown')
      write(11, *) tlen
      write(11, *) np, 1, 3
      write(11, *) omegai, lat(ir), lon(ir)
      write(11, *) eqlat, eqlon, r0
      close(11)
    end do
  else
    write(*, *) 'WARNING:(tish.f)  set spcform 0 or 1'
  end if

  if (ilog == 1) then
    open(unit = 11, file = 'llog.log', status = 'unknown')
    close(11)
  end if



  ! ************************** Option for shallow events **************************
  iimax = imax
  if ((rmax - r0) < shallowdepth) then
    ! ******************* Computing parameters *******************
    ! computing of the number and position of grid points
    iimax = int(tlen * 2.d0)
    call computeKz(nZone, rminOfZone, rmaxOfZone, vsvPolynomials, rmax, iimax, 1, tlen, kz)
    call computeGridRadii(nZone, kz, rminOfZone, rmaxOfZone, rmin, re, nLayer, nLayerInZone, gridRadii)
    if (nLayer > maxNLayer) stop 'The number of grid points is too large.'

    ! computing the first indices in each zone
    call computeFirstIndices(nZone, nLayerInZone, isp, jsp)

    ! computing the source position
    call computeSourcePosition(nLayer, rmaxOfZone, rmin, rmax, gridRadii, isp, r0, iZoneOfSource, qLayerOfSource)

    ! computing grids for source computations
    call computeSourceGrid(isp, gridRadii, r0, iZoneOfSource, qLayerOfSource, gridRadiiForSource)

    ! ******************* Computing the matrix elements *******************
    ! computing variable values at grid points
    call computeStructureValues(nZone, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, nLayerInZone, gridRadii, &
      nValue, valuedRadii, rhoValues, ecLValues, ecNValues)
    call computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, gridRadiiForSource, &
      rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, mu0)

    ! computing mass and rigitidy matrices
    do i = 1, nZone
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, rhoValues, 2, 0, 0, gridRadii(isp(i)), &
        t(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecLValues, 2, 1, 1, gridRadii(isp(i)), &
        h1(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecLValues, 1, 1, 0, gridRadii(isp(i)), &
        h2(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecLValues, 0, 0, 0, gridRadii(isp(i)), &
        h3(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecNValues, 0, 0, 0, gridRadii(isp(i)), &
        h4(jsp(i)), work(jsp(i)))
      call computeLumpedT(nLayerInZone(i), nValue, valuedRadii, rhoValues, gridRadii(isp(i)), work(jsp(i)))
      call computeAverage(nLayerInZone(i), t(jsp(i)), work(jsp(i)), t(jsp(i)))
      call computeLumpedH(nLayerInZone(i), nValue, valuedRadii, ecLValues, gridRadii(isp(i)), work(jsp(i)))
      call computeAverage(nLayerInZone(i), h3(jsp(i)), work(jsp(i)), h3(jsp(i)))
      call computeLumpedH(nLayerInZone(i), nValue, valuedRadii, ecNValues, gridRadii(isp(i)), work(jsp(i)))
      call computeAverage(nLayerInZone(i), h4(jsp(i)), work(jsp(i)), h4(jsp(i)))
    end do

    ! computing mass and rigitidy matrices near source
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

    nn = nLayer + 1
    ns = isp(iZoneOfSource) + int(qLayerOfSource)
    ins = 4 * ns - 3


    llog = 0
    do ii = 1,2  ! omega-loop
      if (ii == 1) then
        if (imin == 0) then
          i = 1
        else
          i = imin
        end if
      end if
      if (ii == 2) i = imax
      omega = 2.d0 * pi * dble(i) / tlen
      call computeLsuf(omega, nZone, rmaxOfZone, vsvPolynomials, lsuf)
      call computeCoef(nZone, omega, qmuOfZone, coef)

      call initComplexMatrix(lda, nn, a0)
      call initComplexMatrix(lda, nn, a2)
      do j = 1, nZone
        call computeA0(nLayerInZone(j), omega, omegai, t(jsp(j)), h1(jsp(j)), h2(jsp(j)), h3(jsp(j)), h4(jsp(j)), coef(j), &
          cwork(jsp(j)))
        call overlapMatrixBlocks(nLayerInZone(j), cwork(jsp(j)), a0(1, isp(j)))
        call computeA2(nLayerInZone(j), h4(jsp(j)), coef(j), cwork(jsp(j)))
        call overlapMatrixBlocks(nLayerInZone(j), cwork(jsp(j)), a2(1, isp(j)))
      end do

      kc = 1
      ismall = 0
      maxamp = -1.d0
      ltmp(ii) = maxlmax
      do l = 0, maxlmax  ! l-loop
        if (ismall > 20) then
          if (ltmp(ii) > l) ltmp(ii) = l
          exit
        end if

        do jj = 1, maxnlayer + 1  ! initialize
          tmpr(jj) = 0.d0
        end do
        ! computing the coefficient matrix elements
        ! --- Initializing the matrix elements
        call initComplexMatrix(lda, nn, a)
        call initComplexMatrix(lda, 3, ga2)
        call assembleA(nn, l, a0, a2, a)
        call computeA(1, omega, omegai, l, t(ins), h1(ins), h2(ins), h3(ins), h4(ins), coef(iZoneOfSource), aa)
        call computeA(2, omega, omegai, l, gt, gh1, gh2, gh3, gh4, coef(iZoneOfSource), ga)
        call overlapMatrixBlocks(2, ga, ga2)

        do m = -2, 2  ! m-loop
          if ((m /= 0) .and. (iabs(m) <= iabs(l))) then
            call initComplexVector(nn, g)
            call computeG(l, m, qLayerOfSource, r0, mt, mu0, coef(iZoneOfSource), ga, aa, ga2, gdr, g(isp(iZoneOfSource)))
            if (mod(l, 100) == 0) then
              if ((m == -2) .or. (m == -l)) then
                call dclisb0(a, nn, 1, lda, g, eps, dr, z, ier)
              else
                call dcsbsub0(a, nn, 1, lda, g, eps, dr, z, ier)
              end if
              do jj = 1, nn  ! sum up c of the same l
                tmpr(jj) = tmpr(jj) + abs(g(jj))
              end do
            else
              if ((m == -2) .or. (m == -l)) then
                call dclisb(a(1, kc), nn - kc + 1, 1, lda, ns - kc + 1, g(kc), eps, dr, z, ier)
              else
                call dcsbsub(a(1, kc), nn - kc + 1, 1, lda, ns - kc + 1, g(kc), eps, dr, z, ier)
              end if
            end if

            if (mod(l, 100) == 0) then
              call calcutd(nZone, nLayerInZone, tmpr, ratc, nn, gridRadii, kc)
            end if

            call calamp(g(nn), l, lsuf, maxamp, ismall, ratl)
          end if
        end do  ! m-loop
      end do  ! l-loop
    end do  ! omega-loop

    iimax = dble(max(ltmp(1),ltmp(2))) * tlen / lmaxdivf
  end if  ! option for shallow events


  ! ******************* Computing parameters *******************
  ! computing of the number and position of grid points
  call computeKz(nZone, rminOfZone, rmaxOfZone, vsvPolynomials, rmax, iimax, 1, tlen, kz)
  call computeGridRadii(nZone, kz, rminOfZone, rmaxOfZone, rmin, re, nLayer, nLayerInZone, gridRadii)
  if (nLayer > maxNLayer) stop 'The number of grid points is too large.'

  ! computing the first indices in each zone
  call computeFirstIndices(nZone, nLayerInZone, isp, jsp)

  ! computing the source position
  call computeSourcePosition(nLayer, rmaxOfZone, rmin, rmax, gridRadii, isp, r0, iZoneOfSource, qLayerOfSource)

  ! computing grids for source computations
  call computeSourceGrid(isp, gridRadii, r0, iZoneOfSource, qLayerOfSource, gridRadiiForSource)


  ! ******************* Computing the matrix elements *******************
  ! computing variable values at grid points
  call computeStructureValues(nZone, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, nLayerInZone, gridRadii, &
    nValue, valuedRadii, rhoValues, ecLValues, ecNValues)
  call computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials, vsvPolynomials, vshPolynomials, gridRadiiForSource, &
    rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, mu0)

  ! computing mass and rigitidy matrices
  do i = 1, nZone
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, rhoValues, 2, 0, 0, gridRadii(isp(i)), &
      t(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecLValues, 2, 1, 1, gridRadii(isp(i)), &
      h1(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecLValues, 1, 1, 0, gridRadii(isp(i)), &
      h2(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecLValues, 0, 0, 0, gridRadii(isp(i)), &
      h3(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nLayerInZone(i), nValue, valuedRadii, ecNValues, 0, 0, 0, gridRadii(isp(i)), &
      h4(jsp(i)), work(jsp(i)))
    call computeLumpedT(nLayerInZone(i), nValue, valuedRadii, rhoValues, gridRadii(isp(i)), work(jsp(i)))
    call computeAverage(nLayerInZone(i), t(jsp(i)), work(jsp(i)), t(jsp(i)))
    call computeLumpedH(nLayerInZone(i), nValue, valuedRadii, ecLValues, gridRadii(isp(i)), work(jsp(i)))
    call computeAverage(nLayerInZone(i), h3(jsp(i)), work(jsp(i)), h3(jsp(i)))
    call computeLumpedH(nLayerInZone(i), nValue, valuedRadii, ecNValues, gridRadii(isp(i)), work(jsp(i)))
    call computeAverage(nLayerInZone(i), h4(jsp(i)), work(jsp(i)), h4(jsp(i)))
  end do

  ! computing mass and rigitidy matrices near source
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
  nn = nLayer + 1
  ns = isp(iZoneOfSource) + nint(qLayerOfSource)
  ins = 4 * ns - 3

  llog = 0
  do i = imin, imax  ! omega-loop
    if (i == 0) cycle
    omega = 2.d0 * pi * dble(i) / tlen

    call initComplexMatrix(3, nReceiver, u)

    ! compute the angular order that is sufficient to compute the slowest phase velocity
    call computeLsuf(omega, nZone, rmaxOfZone, vsvPolynomials, lsuf)

    do ir = 1, nReceiver
      call matinit(3, 4, plm(:, :, ir))
    end do

    ! compute coefficient related to attenuation
    call computeCoef(nZone, omega, qmuOfZone, coef)

    call initComplexMatrix(lda, nn, a0)
    call initComplexMatrix(lda, nn, a2)

    ! compute parts of A matrix (omega^2 T - H) (It is split into parts to exclude l-dependence.)
    do j = 1, nZone
      call computeA0(nLayerInZone(j), omega, omegai, t(jsp(j)), h1(jsp(j)), h2(jsp(j)), h3(jsp(j)), h4(jsp(j)), coef(j), &
        cwork(jsp(j)))
      call overlapMatrixBlocks(nLayerInZone(j), cwork(jsp(j)), a0(1, isp(j)))

      call computeA2(nLayerInZone(j), h4(jsp(j)), coef(j), cwork(jsp(j)))
      call overlapMatrixBlocks(nLayerInZone(j), cwork(jsp(j)), a2(1, isp(j)))
    end do

    kc = 1
    ismall = 0
    maxamp = -1.d0
    llog = maxlmax

    do l = 0, maxlmax  ! l-loop
      if (ismall > 20) then
        llog = min(llog, l)
        cycle
      end if

      tmpr(:) = 0.d0

      ! compute trial functions
      do ir = 1, nReceiver
        call computeTrialFunctionValues(l, theta(ir), phi(ir), plm(:, :, ir), trialFunctionValues(:, :, ir))
      end do

      ! initializing the matrix elements
      call initComplexMatrix(lda, nn, a)
      call initComplexMatrix(lda, 3, ga2)

      ! assemble A matrix from parts that have already been computed
      call assembleA(nn, l, a0, a2, a)

      ! TODO ???
      call computeA(1, omega, omegai, l, t(ins), h1(ins), h2(ins), h3(ins), &
        h4(ins), coef(iZoneOfSource), aa)

      ! compute A matrix near source
      call computeA(2, omega, omegai, l, gt, gh1, gh2, gh3, gh4, coef(iZoneOfSource), ga)
      call overlapMatrixBlocks(2, ga, ga2)

      do m = -2, 2  ! m-loop
        if (m == 0 .or. abs(m) > abs(l)) cycle

        call initComplexVector(nn, g)

        call computeG(l, m, qLayerOfSource, r0, mt, mu0, coef(iZoneOfSource), ga, aa, ga2, gdr, &
          g(isp(iZoneOfSource)))

        if (mod(l, 100) == 0) then
          if (m == -2 .or. m == -l) then
            ! in the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed
            call dclisb0(a, nn, 1, lda, g, eps, dr, z, ier)
          else
            ! in consecutive m-loops, start from forward substitution (decomposition is skipped)
            call dcsbsub0(a, nn, 1, lda, g, eps, dr, z, ier)
          end if

          tmpr(:) = tmpr(:) + abs(g(:))

        else
          if (m == -2 .or. m == -l) then
            ! in the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed
            call dclisb(a(1, kc), nn-kc+1, 1, lda, ns-kc+1, g(kc), eps, dr, z, ier)
          else
            ! in consecutive m-loops, start from forward substitution (decomposition is skipped)
            call dcsbsub(a(1, kc), nn-kc+1, 1, lda, ns-kc+1, g(kc), eps, dr, z, ier)
          end if
        end if

        if (mod(l, 100) == 0) then
          call calcutd(nzone, nLayerInZone, tmpr, ratc, nn, gridRadii, kc)
        end if

        call calamp(g(nn), l, lsuf, maxamp, ismall, ratl)

        ! accumulate u
        do ir = 1, nReceiver
          call computeU(g(nn), l, trialFunctionValues(:, m, ir), u(:, ir))
        end do

      end do  ! m-loop
    end do  ! l-loop

    ! store results
    outputi(outputindex) = i
    do ir = 1, nReceiver
      outputu(:, ir, outputindex) = u(:, ir)
    end do

    ! ************************** Files Handling **************************
    ! Write to file when the output interval is reached, or when this is the last omega.
    if (outputindex >= outputinterval .or. i == imax) then
      write(*,*) "kakikomimasu"
      if (spcform == 0) then
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
      else if (spcform == 1) then
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
        write(*,*) "WARNING: set spcform 0 or 1"
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


