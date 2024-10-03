
!------------------------------------------------------------------------
! Computes matrix elements common for all omega, l, and m.
!------------------------------------------------------------------------
subroutine computeMatrixElements(maxNGrid, shallowDepth, tlen, re, imin, imax, lmax, r0, &
  nZone, rmin, rmax, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, &
  kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oPairOfZone, &
  iZoneOfSource, iLayerOfSource, oPairOfSource, gridRadiiForSource, &
  nValue, valuedRadii, rhoValues, ecLValues, ecNValues, rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, ecL0, &
  t, h1, h2sum, h3, h4, gt, gh1, gh2sum, gh3, gh4, work)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: maxNGrid  ! Maximum number of grid points.
  real(8), intent(in) :: shallowDepth  ! Threshold to consider evanescent regime for shallow events [km].
  real(8), intent(in) :: tlen  ! Time length [s].
  real(8), intent(in) :: re  ! Desired relative error due to vertical gridding.
  integer, intent(in) :: imin, imax  ! Index of minimum and maximum frequency.
  integer, intent(in) :: lmax  ! Largest angular order l, if it is computed. Otherwise, 0.
  real(8), intent(inout) :: r0  ! Source radius [km]. Its value may be fixed in this subroutine.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmin, rmax  ! Minimum and maximum radii of region considered [km].
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  real(8), intent(in) :: rhoPolynomials(4,nZone), vsvPolynomials(4,nZone), vshPolynomials(4,nZone)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vsv, and vsh [km/s] structure.
  real(8), intent(out) :: kzAtZone(nZone)  ! Computed value of vertical wavenumber k_z at each zone [1/km].
  integer, intent(out) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(out) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  real(8), intent(out) :: gridRadii(maxNGrid)  ! Radius at each grid point [km].
  integer, intent(out) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(out) :: oValueOfZone(nZone)  ! Index of the first value in each zone.
  integer, intent(out) :: oPairOfZone(nZone)  ! Index of the first (iLayer, k', k)-pair in each zone.
  integer, intent(out) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(out) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(out) :: oPairOfSource  ! Index of the first (iLayer, k', k)-pair for the layer with the source.
  real(8), intent(out) :: gridRadiiForSource(3)  ! Radii to use for source-related computations [km].
  integer, intent(out) :: nValue  ! Total number of values for each variable.
  real(8), intent(out) :: valuedRadii(maxNGrid+nZone-1)  ! Radii corresponding to each variable value [km].
  real(8), intent(out) :: rhoValues(maxNGrid+nZone-1), ecLValues(maxNGrid+nZone-1), ecNValues(maxNGrid+nZone-1)
  !::::::::::::::::::::::::: Values of rho [g/cm^3], L, and N [10^10 dyn/cm^2 = GPa] at each point (with 2 values at boundaries).
  real(8), intent(out) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::  Values of rho [g/cm^3], L, and N [10^10 dyn/cm^2 = GPa] at each point.
  real(8), intent(out) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  real(8), intent(out) :: t(4*maxNGrid-4), h1(4*maxNGrid-4), h2sum(4*maxNGrid-4), h3(4*maxNGrid-4), h4(4*maxNGrid-4)
  real(8), intent(out) :: gt(8), gh1(8), gh2sum(8), gh3(8), gh4(8)
  real(8), intent(out) :: work(4*maxNGrid-4)  ! Working matrix.
  integer :: i, oV, oP
  real(8) :: rShallowThreshold

  ! ------------------- Computing parameters -------------------
  ! Design the number and position of grid points.
  rShallowThreshold = rmax - shallowDepth
  call computeKz(rShallowThreshold, nZone, rminOfZone(:), rmaxOfZone(:), vsvPolynomials(:,:), rmax, imin, imax, 1, lmax, &
    tlen, kzAtZone(:))
  call computeGridRadii(maxNGrid, nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), rmin, re, nGrid, nLayerInZone(:), gridRadii(:))

  ! Compute the first indices in each zone.
  call computeFirstIndices(nZone, nLayerInZone(:), oGridOfZone(:), oValueOfZone(:), oPairOfZone(:))

  ! Compute the source position.
  call computeSourcePosition(nGrid, rmaxOfZone(:), gridRadii(:), r0, iZoneOfSource, iLayerOfSource, oPairOfSource)

  ! Design grids for source computations.
  call computeSourceGrid(gridRadii(:), r0, iLayerOfSource, gridRadiiForSource(:))


  ! ------------------- Computing the matrix elements -------------------
  ! Compute variable values at grid points.
  call computeStructureValues(nZone, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), nLayerInZone(:), &
    gridRadii(:), nValue, valuedRadii(:), rhoValues(:), ecLValues(:), ecNValues(:))
  call computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), &
    gridRadiiForSource(:), rhoValuesForSource(:), ecLValuesForSource(:), ecNValuesForSource(:), ecL0)

  ! Compute mass and rigitidy matrices.
  do i = 1, nZone
    oV = oValueOfZone(i)
    oP = oPairOfZone(i)
    ! Compute unmodified matrices.
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), 2, 0, 0, t(oP:))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 2, 1, 1, h1(oP:))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 1, 1, 0, work(oP:))
    call addTranspose(nLayerInZone(i), work(oP:), h2sum(oP:))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 0, 0, 0, h3(oP:))
    call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 0, 0, 0, h4(oP:))
    ! Modify matrices for I^(0) using lumped matrices.
    call computeLumpedT(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), work(oP:))
    call averageMatrix(nLayerInZone(i), t(oP:), work(oP:), t(oP:))
    call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), work(oP:))
    call averageMatrix(nLayerInZone(i), h3(oP:), work(oP:), h3(oP:))
    call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), work(oP:))
    call averageMatrix(nLayerInZone(i), h4(oP:), work(oP:), h4(oP:))
  end do

  ! Compute mass and rigitidy matrices near source.
  ! Compute unmodified matrices.
  call computeIntermediateIntegral(2, gridRadiiForSource(:), rhoValuesForSource(:), 2, 0, 0, gt(:))
  call computeIntermediateIntegral(2, gridRadiiForSource(:), ecLValuesForSource(:), 2, 1, 1, gh1(:))
  call computeIntermediateIntegral(2, gridRadiiForSource(:), ecLValuesForSource(:), 1, 1, 0, work(:))
  call addTranspose(2, work(:), gh2sum(:))
  call computeIntermediateIntegral(2, gridRadiiForSource(:), ecLValuesForSource(:), 0, 0, 0, gh3(:))
  call computeIntermediateIntegral(2, gridRadiiForSource(:), ecNValuesForSource(:), 0, 0, 0, gh4(:))
  ! Modify matrices for I^(0) using lumped matrices.
  call computeLumpedT(2, gridRadiiForSource(:), rhoValuesForSource(:), work(:))
  call averageMatrix(2, gt(:), work(:), gt(:))
  call computeLumpedH(2, gridRadiiForSource(:), ecLValuesForSource(:), work(:))
  call averageMatrix(2, gh3(:), work(:), gh3(:))
  call computeLumpedH(2, gridRadiiForSource(:), ecNValuesForSource(:), work(:))
  call averageMatrix(2, gh4(:), work(:), gh4(:))

end subroutine


!------------------------------------------------------------------------
! Content of loop for shallow event for 1 value of omega.
!------------------------------------------------------------------------
subroutine omegaLoopForShallowEvent(omega, omegaI, maxL, nZone, rmaxOfZone, vsvPolynomials, qmuOfZone, &
  r0, mt, ecL0, ratc, ratl, amplitudeAtGrid, &
  nGrid, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, coefQmu, &
  t, h1, h2sum, h3, h4, gt, gh1, gh2, gh3, gh4, oPairOfZone, oPairOfSource, aaParts, aSourceParts, aSource, &
  a0, a2, a, g_or_c, cwork, dr, z, gdr, eps, ltmpI)
!------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  integer, intent(in) :: maxL  ! Maximum of angular order to loop for.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmaxOfZone(nZone)  ! Upper radii of each zone [km].
  real(8), intent(in) :: vsvPolynomials(4,nZone)  ! Polynomial function of vsv [km/s] structure.
  real(8), intent(in) :: qmuOfZone(nZone)  ! Qmu of each zone.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  real(8), intent(in) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8), intent(in) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8), intent(out) :: amplitudeAtGrid(nGrid)  ! Estimate of the amplitude at each grid point [km].
  integer, intent(in) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  integer, intent(in) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  complex(8), intent(out) :: coefQmu(nZone)  ! Coefficients to multiply to elastic moduli for attenuation at each zone.
  real(8), intent(in) :: t(4 * nGrid - 4)
  real(8), intent(in) :: h1(4 * nGrid - 4), h2sum(4 * nGrid - 4), h3(4 * nGrid - 4), h4(4 * nGrid - 4)
  real(8), intent(in) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  integer, intent(in) :: oPairOfZone(nZone)  ! Index of the first (iLayer, k', k)-pair in each zone.
  integer, intent(in) :: oPairOfSource  ! Index of the first (iLayer, k', k)-pair for the layer with the source.
  complex(8), intent(out) :: aaParts(4), aSourceParts(8)  ! Unassembled A matrix [10^12 kg/s^2].
  complex(8), intent(out) :: aSource(2,3)  ! Assembled A matrix [10^12 kg/s^2].
  complex(8), intent(out) :: a0(2, nGrid), a2(2, nGrid)
  complex(8), intent(out) :: a(2,nGrid)  ! Assembled A matrix [10^12 kg/s^2].
  complex(8), intent(out) :: g_or_c(nGrid)  ! This holds either vector g [10^15 N] or c [km], depending on where in the code it is.
  complex(8), intent(out) :: cwork(4 * nGrid - 4)  ! Working array for matrix computations.
  complex(8), intent(inout) :: dr(nGrid), z(nGrid), gdr(3)  ! Working arrays used when solving linear equations.
  real(8), intent(inout) :: eps
  integer :: ltmpI

  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: largeL2  ! L^2 = l(l+1).
  integer :: cutoffGrid  ! Index of grid at cut-off depth.
  integer :: lsuf  ! Accuracy threshold of angular order. (Corresponds to l_d; see eq. 29 of Kawai et al. 2006.)
  real(8) :: recordAmplitude  ! Maximum amplitude encountered [km], used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: i
  integer :: oP

  ! Initialize matrices.
  a0(:, :nGrid) = dcmplx(0.d0, 0.d0)
  a2(:, :nGrid) = dcmplx(0.d0, 0.d0)

  ! Compute the angular order that is sufficient to compute the slowest phase velocity.
  call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

  ! Compute coefficients to multiply to elastic moduli for anelastic attenuation.
  call computeCoef(nZone, omega, qmuOfZone(:), coefQmu(:))

  ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
  do i = 1, nZone
    oP = oPairOfZone(i)
    call computeA0(nLayerInZone(i), omega, omegaI, t(oP:), h1(oP:), h2sum(oP:), h3(oP:), h4(oP:), coefQmu(i), cwork(oP:))
    call overlapA(nLayerInZone(i), cwork(oP:), a0(:, oGridOfZone(i):))
    call computeA2(nLayerInZone(i), h4(oP:), coefQmu(i), cwork(oP:))
    call overlapA(nLayerInZone(i), cwork(oP:), a2(:, oGridOfZone(i):))
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
    a(:, :nGrid) = dcmplx(0.d0, 0.d0)
    aSource(:, :) = dcmplx(0.d0, 0.d0)
    ! Clear the amplitude accumulated for all m's.
    if (mod(l, 100) == 0) amplitudeAtGrid(:nGrid) = 0.d0

    ! Assemble A matrix from parts that have already been computed.
    call assembleA(nGrid, largeL2, a0(:,:), a2(:,:), a(:,:))

    ! Compute UNASSEMBLED A matrix in layer with source.
    ! NOTE that a(:,:) cannot be used instead of aaParts(:), because a(:,:) is already assembled.
    call computeA(1, omega, omegaI, largeL2, t(oPairOfSource:), &
      h1(oPairOfSource:), h2sum(oPairOfSource:), h3(oPairOfSource:), h4(oPairOfSource:), coefQmu(iZoneOfSource), aaParts(:))

    ! Compute A matrix near source.
    call computeA(2, omega, omegaI, largeL2, gt(:), gh1(:), gh2(:), gh3(:), gh4(:), coefQmu(iZoneOfSource), aSourceParts(:))
    call overlapA(2, aSourceParts(:), aSource(:,:))

    do m = -2, 2  ! m-loop
      if (m == 0 .or. abs(m) > abs(l)) cycle

      ! Form and solve the linear equation Ac=-g.
      call formAndSolveEquation(l, m, iZoneOfSource, iLayerOfSource, r0, mt, ecL0, coefQmu, aaParts, aSourceParts, aSource, &
        nGrid, cutoffGrid, a, eps, g_or_c, amplitudeAtGrid, dr, z, gdr)

      ! Check whether the amplitude has decayed enough to stop the l-loops.
      !  This is checked for the topmost-grid expansion coefficent of each m individually.
      call checkAmplitudeDecay(g_or_c(nGrid), l, lsuf, ratl, recordAmplitude, decayCounter)

    end do  ! m-loop

    ! Decide cut-off depth (at a certain interval of l).
    if (mod(l, 100) == 0) then
      call computeCutoffGrid(nGrid, amplitudeAtGrid(:), ratc, cutoffGrid)
    end if

  end do  ! l-loop

  ! Register the final l (or maxL instead of maxL-1 when all loops are completed).  !!! difference from main section
  ltmpI = min(l, maxL)

end subroutine


!------------------------------------------------------------------------
! Content of main loop for 1 value of omega.
!------------------------------------------------------------------------
subroutine omegaLoop(omega, omegaI, maxL, nZone, rmaxOfZone, vsvPolynomials, qmuOfZone, &
  r0, mt, ecL0, nReceiver, theta, phi, ratc, ratl, amplitudeAtGrid, &
  nGrid, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, coefQmu, plm, harmonicsValues, &
  t, h1, h2sum, h3, h4, gt, gh1, gh2, gh3, gh4, oPairOfZone, oPairOfSource, aaParts, aSourceParts, aSource, &
  a0, a2, a, g_or_c, u, cwork, dr, z, gdr, eps, llog)
!------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  integer, intent(in) :: maxL  ! Maximum of angular order to loop for.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmaxOfZone(nZone)  ! Upper radii of each zone [km].
  real(8), intent(in) :: vsvPolynomials(4,nZone)  ! Polynomial function of vsv [km/s] structure.
  real(8), intent(in) :: qmuOfZone(nZone)  ! Qmu of each zone.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  integer, intent(in) :: nReceiver  ! Number of receivers.
  real(8), intent(in) :: theta(nReceiver), phi(nReceiver)  ! Colatitude and longitude of receiver with event at north pole [rad].
  real(8), intent(in) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8), intent(in) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8), intent(out) :: amplitudeAtGrid(nGrid)  ! Estimate of the amplitude at each grid point [km].
  integer, intent(in) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  integer, intent(in) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  complex(8), intent(out) :: coefQmu(nZone)  ! Coefficients to multiply to elastic moduli for attenuation at each zone.
  real(8), intent(out) :: plm(3, 0:3, nReceiver)
  !:::::::::::::::::::::::::::::::::::::: Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  complex(8), intent(out) :: harmonicsValues(3, -2:2, nReceiver)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Vector harmonics term. The coefficient 1/largeL is not multiplied yet.
  real(8), intent(in) :: t(4 * nGrid - 4)
  real(8), intent(in) :: h1(4 * nGrid - 4), h2sum(4 * nGrid - 4), h3(4 * nGrid - 4), h4(4 * nGrid - 4)
  real(8), intent(in) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  integer, intent(in) :: oPairOfZone(nZone)  ! Index of the first (iLayer, k', k)-pair in each zone.
  integer, intent(in) :: oPairOfSource  ! Index of the first (iLayer, k', k)-pair for the layer with the source.
  complex(8), intent(out) :: aaParts(4), aSourceParts(8)  ! Unassembled A matrix [10^12 kg/s^2].
  complex(8), intent(out) :: aSource(2,3)  ! Assembled A matrix [10^12 kg/s^2].
  complex(8), intent(out) :: a0(2, nGrid), a2(2, nGrid)
  complex(8), intent(out) :: a(2,nGrid)  ! Assembled A matrix [10^12 kg/s^2].
  complex(8), intent(out) :: g_or_c(nGrid)  ! This holds either vector g [10^15 N] or c [km], depending on where in the code it is.
  complex(8), intent(out) :: u(3, nReceiver)  ! Displacement velocity - the unit is [km] in the frequency domain,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: but when converted to the time domain, the unit becomes [km/s].
  complex(8), intent(out) :: cwork(4 * nGrid - 4)  ! Working array for matrix computations.
  complex(8), intent(inout) :: dr(nGrid), z(nGrid), gdr(3)  ! Working arrays used when solving linear equations.
  real(8), intent(inout) :: eps
  integer, intent(out) :: llog

  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: largeL2  ! L^2 = l(l+1).
  integer :: cutoffGrid  ! Index of grid at cut-off depth.
  integer :: lsuf  ! Accuracy threshold of angular order. (Corresponds to l_d; see eq. 29 of Kawai et al. 2006.)
  real(8) :: recordAmplitude  ! Maximum amplitude encountered [km], used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: i
  integer :: ir
  integer :: oP

  ! Initialize matrices.
  a0(:, :nGrid) = dcmplx(0.d0, 0.d0)
  a2(:, :nGrid) = dcmplx(0.d0, 0.d0)
  u(:, :nReceiver) = dcmplx(0.d0, 0.d0)
  ! Plm must be cleared for each omega.  !!! difference from shallow-source section
  plm(:, :, :nReceiver) = 0.d0

  ! Compute the angular order that is sufficient to compute the slowest phase velocity.
  call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

  ! Compute coefficients to multiply to elastic moduli for anelastic attenuation.
  call computeCoef(nZone, omega, qmuOfZone(:), coefQmu(:))

  ! Compute parts of A matrix (omega^2 T - H). (It is split into parts to exclude l-dependence.)
  do i = 1, nZone
    oP = oPairOfZone(i)
    call computeA0(nLayerInZone(i), omega, omegaI, t(oP:), h1(oP:), h2sum(oP:), h3(oP:), h4(oP:), coefQmu(i), cwork(oP:))
    call overlapA(nLayerInZone(i), cwork(oP:), a0(:, oGridOfZone(i):))
    call computeA2(nLayerInZone(i), h4(oP:), coefQmu(i), cwork(oP:))
    call overlapA(nLayerInZone(i), cwork(oP:), a2(:, oGridOfZone(i):))
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
    a(:, :nGrid) = dcmplx(0.d0, 0.d0)
    aSource(:, :) = dcmplx(0.d0, 0.d0)
    ! Clear the amplitude accumulated for all m's.
    if (mod(l, 100) == 0) amplitudeAtGrid(:nGrid) = 0.d0

    ! Compute trial functions.  !!! difference from shallow-source section
    do ir = 1, nReceiver
      call computeHarmonicsValues(l, theta(ir), phi(ir), plm(:, :, ir), harmonicsValues(:, :, ir))
    end do

    ! Assemble A matrix from parts that have already been computed.
    call assembleA(nGrid, largeL2, a0(:,:), a2(:,:), a(:,:))

    ! Compute UNASSEMBLED A matrix in layer with source.
    ! NOTE that a(:,:) cannot be used instead of aaParts(:), because a(:,:) is already assembled.
    call computeA(1, omega, omegaI, largeL2, t(oPairOfSource:), &
      h1(oPairOfSource:), h2sum(oPairOfSource:), h3(oPairOfSource:), h4(oPairOfSource:), coefQmu(iZoneOfSource), aaParts(:))

    ! Compute A matrix near source.
    call computeA(2, omega, omegaI, largeL2, gt(:), gh1(:), gh2(:), gh3(:), gh4(:), coefQmu(iZoneOfSource), aSourceParts(:))
    call overlapA(2, aSourceParts(:), aSource(:,:))

    do m = -2, 2  ! m-loop
      if (m == 0 .or. abs(m) > abs(l)) cycle

      ! Form and solve the linear equation Ac=-g.
      call formAndSolveEquation(l, m, iZoneOfSource, iLayerOfSource, r0, mt, ecL0, coefQmu, aaParts, aSourceParts, aSource, &
        nGrid, cutoffGrid, a, eps, g_or_c, amplitudeAtGrid, dr, z, gdr)

      ! Check whether the amplitude has decayed enough to stop the l-loops.
      !  This is checked for the topmost-grid expansion coefficent of each m individually.
      call checkAmplitudeDecay(g_or_c(nGrid), l, lsuf, ratl, recordAmplitude, decayCounter)

      ! Accumulate u.  !!! difference from shallow-source section
      do ir = 1, nReceiver
        call computeU(g_or_c(nGrid), largeL2, harmonicsValues(:, m, ir), u(:, ir))
      end do

    end do  ! m-loop

    ! Decide cut-off depth (at a certain interval of l).
    if (mod(l, 100) == 0) then
      call computeCutoffGrid(nGrid, amplitudeAtGrid(:), ratc, cutoffGrid)
    end if

  end do  ! l-loop

  ! Register the final l (or maxL instead of maxL-1 when all loops are completed).
  llog = min(l, maxL)

end subroutine


!------------------------------------------------------------------------
! Form and solve the linear equation Ac=-g.
!------------------------------------------------------------------------
subroutine formAndSolveEquation(l, m, iZoneOfSource, iLayerOfSource, r0, mt, ecL0, coefQmu, aaParts, aSourceParts, aSource, &
  nGrid, cutoffGrid, a, eps, g_or_c, amplitudeAtGrid, dr, z, gdr)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: l  ! Angular order.
  integer, intent(in) :: m  ! Azimuthal order.
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  complex(8), intent(in) :: coefQmu(*)  ! Coefficients to multiply to elastic moduli for attenuation at each zone.
  complex(8), intent(in) :: aaParts(4), aSourceParts(8)  ! Unassembled A matrix [10^12 kg/s^2].
  complex(8), intent(inout) :: aSource(2,3)  ! Assembled A matrix [10^12 kg/s^2].
  integer, intent(in) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(in) :: cutoffGrid  ! Index of grid at cut-off depth.
  complex(8), intent(inout) :: a(2,nGrid)  ! Assembled A matrix [10^12 kg/s^2].
  complex(8), intent(out) :: g_or_c(nGrid)  ! This holds either vector g [10^15 N] or c [km], depending on where in the code it is.
  real(8), intent(inout) :: amplitudeAtGrid(nGrid)  ! Estimate of the amplitude at each grid point [km].
  real(8), intent(inout) :: eps
  complex(8), intent(inout) :: dr(nGrid), z(nGrid), gdr(3)  ! Working arrays used when solving linear equations.
  integer :: ier  ! Error code from subroutine solving linear equations.

  ! Compute excitation vector g.
  g_or_c(:) = dcmplx(0.d0, 0.d0)
  call computeG(l, m, iLayerOfSource, r0, mt, ecL0, coefQmu(iZoneOfSource), aaParts(:), aSourceParts(:), aSource(:,:), &
    gdr(:), g_or_c(:))

  if (mod(l, 100) == 0) then
    ! Once in a while, compute for all grids to decide the cut-off depth.
    ! CAUTION: In this case, all values of g_or_c(:) are computed.

    ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
    ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
    if (m == -2 .or. m == -l) then
      call decomposeAByCholesky(a(:,:), nGrid, 1, 2, eps, dr, ier)
    end if
    ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
    call solveWholeCAfterCholesky(a(:,:), nGrid, 1, 2, g_or_c(:), dr, z)

    ! Accumulate the absolute values of expansion coefficent c for all m's at each grid point.
    !  This is to be used as an estimate of the amplitude at each depth when deciding the cut-off depth.
    amplitudeAtGrid(1:nGrid) = amplitudeAtGrid(1:nGrid) + abs(g_or_c(1:nGrid))

  else
    ! Otherwise, compute for just the grids above the cut-off depth.
    ! CAUTION: In this case, only g_or_c(nGrid) is computed. Other values of g_or_c(:nGrid-1) still hold values of g!!!

    ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
    ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
    if (m == -2 .or. m == -l) then
      call decomposeAByCholesky(a(:, cutoffGrid:), nGrid - cutoffGrid + 1, 1, 2, eps, dr, ier)
    end if
    ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
    call solveSurfaceCAfterCholesky(a(:, cutoffGrid:), nGrid - cutoffGrid + 1, 1, 2, iLayerOfSource - cutoffGrid + 1, &
      g_or_c(cutoffGrid:), dr, z)

  end if

end subroutine
