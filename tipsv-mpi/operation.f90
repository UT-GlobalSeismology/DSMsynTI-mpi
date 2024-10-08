
!------------------------------------------------------------------------
! Computes matrix elements common for all omega, l, and m.
!------------------------------------------------------------------------
subroutine computeMatrixElements(maxNGrid, maxNGridSolid, maxNGridFluid, shallowDepth, tlen, re, imin, imax, lmax, r0, &
  nZone, rmin, rmax, rminOfZone, rmaxOfZone, phaseOfZone, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, &
  kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oValueOfZoneSolid, &
  oPairOfZoneSolid, oPairOfZoneFluid, oElementOfZone, oColumnOfZone, nColumn, &
  iZoneOfSource, iLayerOfSource, oColumnOfSource, &
  nValue, valuedRadii, rhoValues, kappaValues, ecKxValues, ecKyValues, ecKzValues, ecLValues, ecNValues, ecC0, ecF0, ecL0, &
  rhoReciprocals, kappaReciprocals, &
  t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
  hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, work)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: maxNGrid  ! Maximum number of grid points.
  integer, intent(in) :: maxNGridSolid, maxNGridFluid  ! Maximum number of grid points in solid and fluid regions.
  real(8), intent(in) :: shallowDepth  ! Threshold to consider evanescent regime for shallow events [km].
  real(8), intent(in) :: tlen  ! Time length [s].
  real(8), intent(in) :: re  ! Desired relative error due to vertical gridding.
  integer, intent(in) :: imin, imax  ! Index of minimum and maximum frequency.
  integer, intent(in) :: lmax  ! Largest angular order l, if it is computed. Otherwise, 0.
  real(8), intent(inout) :: r0  ! Source radius [km]. Its value may be fixed in this subroutine.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmin, rmax  ! Minimum and maximum radii of region considered [km].
  real(8), intent(in) :: rminOfZone(nZone), rmaxOfZone(nZone)  ! Lower and upper radii of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  real(8), intent(in) :: rhoPolynomials(4,nZone), vpvPolynomials(4,nZone), vphPolynomials(4,nZone)
  real(8), intent(in) :: vsvPolynomials(4,nZone), vshPolynomials(4,nZone), etaPolynomials(4,nZone)
  !:::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta structure.
  real(8), intent(out) :: kzAtZone(nZone)  ! Computed value of vertical wavenumber k_z at each zone [1/km].
  integer, intent(out) :: nGrid  ! Total number of grid points (= number of layers + 1).
  integer, intent(out) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  real(8), intent(out) :: gridRadii(maxNGrid)  ! Radius at each grid point [km].
  integer, intent(out) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(out) :: oValueOfZone(nZone)  ! Index of the first value in each zone.
  integer, intent(out) :: oValueOfZoneSolid(nZone)  ! Index of the first value in each zone, when counting only solid zones.
  integer, intent(out) :: oPairOfZoneSolid(nZone), oPairOfZoneFluid(nZone)
  !::::::::::::::::::::: Index of the first (iLayer, k', k)-pair in each zone, counted separately for solid and fluid zones.
  integer, intent(out) :: oElementOfZone(nZone)  ! Index of the first (iLayer, k'-gamma', k-gamma)-pair in each zone.
  integer, intent(out) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  integer, intent(out) :: nColumn  ! Total number of columns in the band matrix.
  integer, intent(out) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(out) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(out) :: oColumnOfSource  ! Index of the first column in the band matrix for the layer with source.
  integer, intent(out) :: nValue  ! Total number of values for each variable.
  real(8), intent(out) :: valuedRadii(maxNGrid+nZone-1)  ! Radii corresponding to each variable value [km].
  real(8), intent(out) :: rhoValues(maxNGrid+nZone-1), kappaValues(maxNGrid+nZone-1)
  real(8), intent(out) :: ecKxValues(maxNGrid+nZone-1), ecKyValues(maxNGrid+nZone-1), ecKzValues(maxNGrid+nZone-1)
  real(8), intent(out) :: ecLValues(maxNGrid+nZone-1), ecNValues(maxNGrid+nZone-1)
  !::::::::::::::::::::::::::::::::::: Values of rho [g/cm^3], kappa, A-4N/3, F+2N/3, (C+2F)/3, L, and N [10^10 dyn/cm^2 = GPa]
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: at each point (with 2 values at boundaries).
  real(8), intent(out) :: ecC0, ecF0, ecL0  ! Elastic moduli C, F, and L at source position [10^10 dyn/cm^2 = GPa].
  real(8), intent(out) :: rhoReciprocals(maxNGrid+nZone-1), kappaReciprocals(maxNGrid+nZone-1)
  !:::::::::: Values of 1/rho [cm^3/g] and 1/kappa [(10^10 dyn/cm^2)^(-1) = 1/GPa] at each point (with 2 values at boundaries).
  real(8), intent(out) :: t(4*maxNGridSolid-4)
  real(8), intent(out) :: h1x(4*maxNGridSolid-4), h2L(4*maxNGridSolid-4), h2N(4*maxNGridSolid-4)
  real(8), intent(out) :: hUn3y(4*maxNGridSolid-4), hResid3y(4*maxNGridSolid-4), hModL3y(-2:1,maxNGridSolid)
  real(8), intent(out) :: hUn4L(4*maxNGridSolid-4), hResid4L(4*maxNGridSolid-4), hModR4L(-1:2,maxNGridSolid)
  real(8), intent(out) :: hUn4N(4*maxNGridSolid-4), hResid4N(4*maxNGridSolid-4), hModL4N(-2:1,maxNGridSolid)
  real(8), intent(out) :: hUn5y(4*maxNGridSolid-4), hResid5y(4*maxNGridSolid-4), hModR5y(-1:2,maxNGridSolid)
  real(8), intent(out) :: hUn6L(4*maxNGridSolid-4), hResid6L(4*maxNGridSolid-4), hModL6L(-2:1,maxNGridSolid)
  real(8), intent(out) :: hUn6N(4*maxNGridSolid-4), hResid6N(4*maxNGridSolid-4), hModR6N(-1:2,maxNGridSolid)
  real(8), intent(out) :: h7y(4*maxNGridSolid-4), h7z(4*maxNGridSolid-4), h8L(4*maxNGridSolid-4), h8N(4*maxNGridSolid-4)
  real(8), intent(out) :: p1(4*maxNGridFluid-4), p2(4*maxNGridFluid-4), p3(4*maxNGridFluid-4)
  real(8), intent(out) :: work(4*maxNGrid-4)  ! Working matrix.
  integer :: i, iSolid, iFluid, oV, oVS, oP
  real(8) :: rShallowThreshold

  ! ------------------- Computing parameters -------------------
  ! Design the number and position of grid points.
  rShallowThreshold = rmax - shallowDepth
  call computeKz(rShallowThreshold, nZone, rminOfZone(:), rmaxOfZone(:), phaseOfZone(:), &
    vpvPolynomials(:,:), vsvPolynomials(:,:), rmax, imin, imax, 1, lmax, tlen, kzAtZone(:))
  call computeGridRadii(maxNGrid, maxNGridSolid, maxNGridFluid, nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), phaseOfZone(:), &
    rmin, re, nGrid, nLayerInZone(:), gridRadii(:))

  ! Compute the first indices in each zone.
  call computeFirstIndices(nZone, nLayerInZone(:), phaseOfZone(:), oGridOfZone(:), oValueOfZone(:), oValueOfZoneSolid(:), &
    oPairOfZoneSolid(:), oPairOfZoneFluid(:), oElementOfZone(:), oColumnOfZone(:), nColumn)

  ! Compute the source position.
  call computeSourcePosition(nGrid, rmaxOfZone(:), phaseOfZone(:), gridRadii(:), r0, iZoneOfSource, iLayerOfSource)
  oColumnOfSource = oColumnOfZone(iZoneOfSource) + 2 * (iLayerOfSource - oGridOfZone(iZoneOfSource))


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
      oP = oPairOfZoneSolid(iSolid)
      oVS = oValueOfZoneSolid(iSolid)

      ! Compute unmodified matrices.
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), 2, 0, 0, t(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), 0, 0, 0, h1x(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 0, 0, 0, h2L(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 0, 0, 0, h2N(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), 1, 0, 1, hUn5y(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 1, 0, 1, hUn6L(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 1, 0, 1, hUn6N(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), 2, 1, 1, h7y(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecKzValues(oV:), 2, 1, 1, h7z(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), 2, 1, 1, h8L(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), 2, 1, 1, h8N(oP:))
      call transposeMatrix(nLayerInZone(i), hUn5y(oP:), hUn3y(oP:))
      call transposeMatrix(nLayerInZone(i), hUn6L(oP:), hUn4L(oP:))
      call transposeMatrix(nLayerInZone(i), hUn6N(oP:), hUn4N(oP:))
      ! Modify matrices for I^(0) using lumped matrices.
      call computeLumpedT(nLayerInZone(i), valuedRadii(oV:), rhoValues(oV:), work(oP:))
      call averageMatrix(nLayerInZone(i), t(oP:), work(oP:), t(oP:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecKxValues(oV:), work(oP:))
      call averageMatrix(nLayerInZone(i), h1x(oP:), work(oP:), h1x(oP:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), work(oP:))
      call averageMatrix(nLayerInZone(i), h2L(oP:), work(oP:), h2L(oP:))
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), work(oP:))
      call averageMatrix(nLayerInZone(i), h2N(oP:), work(oP:), h2N(oP:))
      ! Compute residual after subtracting step-wise matrix from unmodified I^(1) matrix.
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), work(oP:))
      call subtractMatrix(nLayerInZone(i), hUn5y(oP:), work(oP:), hResid5y(oP:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), work(oP:))
      call subtractMatrix(nLayerInZone(i), hUn6L(oP:), work(oP:), hResid6L(oP:))
      call computeStepH(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), work(oP:))
      call subtractMatrix(nLayerInZone(i), hUn6N(oP:), work(oP:), hResid6N(oP:))
      call transposeMatrix(nLayerInZone(i), hResid5y(oP:), hResid3y(oP:))
      call transposeMatrix(nLayerInZone(i), hResid6L(oP:), hResid4L(oP:))
      call transposeMatrix(nLayerInZone(i), hResid6N(oP:), hResid4N(oP:))
      ! Compute modified matrices for I^(1).
      call computeModifiedHR(nLayerInZone(i), valuedRadii(oV:), ecKyValues(oV:), hModR5y(-1:2, oVS:))
      call computeModifiedHR(nLayerInZone(i), valuedRadii(oV:), ecNValues(oV:), hModR6N(-1:2, oVS:))
      call computeModifiedHL(nLayerInZone(i), valuedRadii(oV:), ecLValues(oV:), hModL6L(-2:1, oVS:))
      call transposeMatrixMod(nLayerInZone(i), -1, 2, hModR5y(-1:2, oVS:), hModL3y(-2:1, oVS:))
      call transposeMatrixMod(nLayerInZone(i), -1, 2, hModR6N(-1:2, oVS:), hModL4N(-2:1, oVS:))
      call transposeMatrixMod(nLayerInZone(i), -2, 1, hModL6L(-2:1, oVS:), hModR4L(-1:2, oVS:))

    else
      ! fluid
      iFluid = iFluid + 1
      oP = oPairOfZoneFluid(iFluid)

      ! Compute unmodified matrices.
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), 2, 1, 1, p1(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), 0, 0, 0, p2(oP:))
      call computeIntermediateIntegral(nLayerInZone(i), valuedRadii(oV:), kappaReciprocals(oV:), 2, 0, 0, p3(oP:))
      ! Modify matrices for I^(0) using lumped matrices.
      call computeLumpedH(nLayerInZone(i), valuedRadii(oV:), rhoReciprocals(oV:), work(oP:))
      call averageMatrix(nLayerInZone(i), p2(oP:), work(oP:), p2(oP:))
      call computeLumpedT(nLayerInZone(i), valuedRadii(oV:), kappaReciprocals(oV:), work(oP:))
      call averageMatrix(nLayerInZone(i), p3(oP:), work(oP:), p3(oP:))

    end if
  end do

end subroutine


!------------------------------------------------------------------------
! Content of loop for shallow event for 1 value of omega.
!------------------------------------------------------------------------
subroutine omegaLoopForShallowEvent(omega, omegaI, maxL, maxNGridSolid, maxNGridFluid, &
  nZone, rmin, rmax, rmaxOfZone, phaseOfZone, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
  r0, mt, ecC0, ecF0, ecL0, ratc, ratl, amplitudeAtColumn, &
  nGrid, gridRadii, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, &
  oValueOfZoneSolid, coefQmu, coefQkappa, coefQfluid, &
  t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
  hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, &
  oPairOfZoneSolid, oPairOfZoneFluid, a0, a1, a2, a, aSmall, g_or_c, g_or_c_Small, &
  oElementOfZone, oColumnOfZone, nColumn, nQuasiColumn, anum, bnum, ya, yb, yc, yd, oColumnOfSource, cwork, z, w, eps, ltmpI)
!------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  integer, intent(in) :: maxL  ! Maximum of angular order to loop for.
  integer, intent(in) :: maxNGridSolid, maxNGridFluid  ! Maximum number of grid points in solid and fluid regions.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmin, rmax  ! Minimum and maximum radii of region considered [km].
  real(8), intent(in) :: rmaxOfZone(nZone)  ! Upper radii of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  real(8), intent(in) :: rhoPolynomials(4,nZone), vpvPolynomials(4,nZone), vphPolynomials(4,nZone)
  real(8), intent(in) :: vsvPolynomials(4,nZone), vshPolynomials(4,nZone), etaPolynomials(4,nZone)
  !:::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta structure.
  real(8), intent(in) :: qmuOfZone(nZone), qkappaOfZone(nZone)  ! Qmu and Qkappa of each zone.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecC0, ecF0, ecL0  ! Elastic moduli C, F, and L at source position [10^10 dyn/cm^2 = GPa].
  real(8), intent(in) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8), intent(in) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8), intent(out) :: amplitudeAtColumn(nColumn)  ! Estimate of the amplitude at each column [km].
  integer, intent(in) :: nGrid  ! Total number of grid points (= number of layers + 1).
  real(8), intent(in) :: gridRadii(nGrid)  ! Radius at each grid point [km].
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  integer, intent(in) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(in) :: oValueOfZoneSolid(nZone)  ! Index of the first value in each zone, when counting only solid zones.
  complex(8), intent(out) :: coefQmu(nZone), coefQkappa(nZone), coefQfluid(nZone)
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  real(8), intent(in) :: t(4*maxNGridSolid-4)
  real(8), intent(in) :: h1x(4*maxNGridSolid-4), h2L(4*maxNGridSolid-4), h2N(4*maxNGridSolid-4)
  real(8), intent(in) :: hUn3y(4*maxNGridSolid-4), hResid3y(4*maxNGridSolid-4), hModL3y(-2:1,maxNGridSolid)
  real(8), intent(in) :: hUn4L(4*maxNGridSolid-4), hResid4L(4*maxNGridSolid-4), hModR4L(-1:2,maxNGridSolid)
  real(8), intent(in) :: hUn4N(4*maxNGridSolid-4), hResid4N(4*maxNGridSolid-4), hModL4N(-2:1,maxNGridSolid)
  real(8), intent(in) :: hUn5y(4*maxNGridSolid-4), hResid5y(4*maxNGridSolid-4), hModR5y(-1:2,maxNGridSolid)
  real(8), intent(in) :: hUn6L(4*maxNGridSolid-4), hResid6L(4*maxNGridSolid-4), hModL6L(-2:1,maxNGridSolid)
  real(8), intent(in) :: hUn6N(4*maxNGridSolid-4), hResid6N(4*maxNGridSolid-4), hModR6N(-1:2,maxNGridSolid)
  real(8), intent(in) :: h7y(4*maxNGridSolid-4), h7z(4*maxNGridSolid-4), h8L(4*maxNGridSolid-4), h8N(4*maxNGridSolid-4)
  real(8), intent(in) :: p1(4*maxNGridFluid-4), p2(4*maxNGridFluid-4), p3(4*maxNGridFluid-4)
  integer, intent(in) :: oPairOfZoneSolid(nZone), oPairOfZoneFluid(nZone)
  !::::::::::::::::::::: Index of the first (iLayer, k', k)-pair in each zone, counted separately for solid and fluid zones.
  complex(8), intent(out) :: a0(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8), intent(out) :: a1(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8), intent(out) :: a2(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8), intent(out) :: a(4, nColumn), aSmall(2, nColumn)  ! Assembled A matrix.
  !::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2] for solid, [m^5/N] for fluid, [10^6 m^2] for boundary elements.
  complex(8), intent(out) :: g_or_c(nColumn), g_or_c_Small(nColumn)
  !::::::::::::::: These hold either vector g ([10^15 N] & [10^9 m^3]) or c ([km] & [GPa]), depending on where in the code it is.
  integer, intent(in) :: oElementOfZone(nZone)  ! Index of the first (iLayer, k'-gamma', k-gamma)-pair in each zone.
  integer, intent(in) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  integer, intent(in) :: nColumn  ! Total number of columns in the band matrix.
  integer, intent(out) :: nQuasiColumn  ! Total number of columns in the rearranged band matrix.
  complex(8), intent(out) :: anum(4, 4, 10), bnum(4, 4, 10)
  complex(8), intent(out) :: ya(4), yb(4), yc(4), yd(4)
  integer, intent(in) :: oColumnOfSource  ! Index of the first column in the band matrix for the layer with source.
  complex(8), intent(out) :: cwork(16 * maxNGridSolid - 16 + 4 * maxNGridFluid - 4)  ! Working array for matrix computations.
  complex(8), intent(inout) :: z(nColumn), w(nColumn)  ! Working arrays used when solving linear equations.
  real(8), intent(inout) :: eps
  integer :: ltmpI

  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: largeL  ! L = sqrt(l(l+1)).
  real(8) :: largeL2  ! L^2 = l(l+1).
  integer :: lsuf  ! Accuracy threshold of angular order.
  real(8) :: recordAmplitude  ! Maximum amplitude encountered [km], used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: cutoffColumn  ! Index of column at cut-off depth.
  integer :: i, iSolid, iFluid
  integer :: oVS
  integer :: oP, oElement, oColumn

  ! Initialize matrices.
  a0(:, :nColumn) = dcmplx(0.d0, 0.d0)
  a1(:, :nColumn) = dcmplx(0.d0, 0.d0)
  a2(:, :nColumn) = dcmplx(0.d0, 0.d0)

  ! Compute the angular order that is sufficient to compute the slowest phase velocity.
  call computeLsuf(omega, nZone, rmaxOfZone(:), vsvPolynomials(:,:), lsuf)

  ! Compute coefficients to multiply to elastic moduli for anelastic attenuation.
  call computeCoef(nZone, omega, qmuOfZone(:), qkappaOfZone(:), coefQmu(:), coefQkappa(:), coefQfluid(:))

  !!TODO organize
  call computeabnum(omega, omegaI, rmax, &
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

    ! Assemble A matrix from parts that have already been computed.
    call assembleAWhole(nZone, phaseOfZone(:), oColumnOfZone(:), largeL2, a0(:,:), a1(:,:), a2(:,:), a(:,:))
    ! Set boundary condition elements
    call setBoundaryConditions(nZone, rmaxOfZone(:), phaseOfZone(:), oColumnOfZone(:), a(:,:))

    !!TODO organize
    call computeya(anum(:,:,:), bnum(:,:,:), largeL2, gridRadii(iLayerOfSource:), r0, ya(:), yb(:), yc(:), yd(:))

    do m = -2, 2  ! m-loop
      if (abs(m) > abs(l)) cycle

      ! Form and solve the linear equation Ac=-g.
      call formAndSolveEquation(l, m, largeL, iZoneOfSource, iLayerOfSource, oColumnOfSource, r0, mt, ecC0, ecF0, ecL0, &
        ya, yb, yc, yd, rmin, nZone, phaseOfZone, oGridOfZone, oColumnOfZone, coefQmu, coefQkappa, &
        nGrid, gridRadii, nColumn, cutoffColumn, a, aSmall, g_or_c, g_or_c_Small, amplitudeAtColumn, nQuasiColumn, eps, z, w)

      if (l > 0) then
        ! Check whether the amplitude has decayed enough to stop the l-loops.
        !  This is checked for the topmost-grid expansion coefficent of each m individually.
        call checkAmplitudeDecay(g_or_c(nColumn-1:nColumn), l, lsuf, ratl, recordAmplitude, decayCounter)

      end if

    end do  ! m-loop

    ! Decide cut-off depth (at a certain interval of l).
    if (mod(l, 100) == 0) then
      call computeCutoffColumn(nZone, phaseOfZone(:), nGrid, oGridOfZone(:), nColumn, oColumnOfZone(:), &
        amplitudeAtColumn(:), ratc, cutoffColumn)
    end if

  end do  ! l-loop

  ! Register the final l (or maxL instead of maxL-1 when all loops are completed).  !!! difference from main section
  ltmpI = min(l, maxL)

end subroutine


!------------------------------------------------------------------------
! Content of main loop for 1 value of omega.
!------------------------------------------------------------------------
subroutine omegaLoop(omega, omegaI, maxL, maxNGridSolid, maxNGridFluid, &
  nZone, rmin, rmax, rmaxOfZone, phaseOfZone, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
  r0, mt, ecC0, ecF0, ecL0, nReceiver, theta, phi, ratc, ratl, amplitudeAtColumn, &
  nGrid, gridRadii, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, &
  oValueOfZoneSolid, coefQmu, coefQkappa, coefQfluid, plm, harmonicsValues, &
  t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
  hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, &
  oPairOfZoneSolid, oPairOfZoneFluid, a0, a1, a2, a, aSmall, g_or_c, g_or_c_Small, u, &
  oElementOfZone, oColumnOfZone, nColumn, nQuasiColumn, anum, bnum, ya, yb, yc, yd, oColumnOfSource, cwork, z, w, eps, llog)
!------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  integer, intent(in) :: maxL  ! Maximum of angular order to loop for.
  integer, intent(in) :: maxNGridSolid, maxNGridFluid  ! Maximum number of grid points in solid and fluid regions.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmin, rmax  ! Minimum and maximum radii of region considered [km].
  real(8), intent(in) :: rmaxOfZone(nZone)  ! Upper radii of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  real(8), intent(in) :: rhoPolynomials(4,nZone), vpvPolynomials(4,nZone), vphPolynomials(4,nZone)
  real(8), intent(in) :: vsvPolynomials(4,nZone), vshPolynomials(4,nZone), etaPolynomials(4,nZone)
  !:::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta structure.
  real(8), intent(in) :: qmuOfZone(nZone), qkappaOfZone(nZone)  ! Qmu and Qkappa of each zone.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecC0, ecF0, ecL0  ! Elastic moduli C, F, and L at source position [10^10 dyn/cm^2 = GPa].
  integer, intent(in) :: nReceiver  ! Number of receivers.
  real(8), intent(in) :: theta(nReceiver), phi(nReceiver)  ! Colatitude and longitude of receiver with event at north pole [rad].
  real(8), intent(in) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8), intent(in) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8), intent(out) :: amplitudeAtColumn(nColumn)  ! Estimate of the amplitude at each column [km].
  integer, intent(in) :: nGrid  ! Total number of grid points (= number of layers + 1).
  real(8), intent(in) :: gridRadii(nGrid)  ! Radius at each grid point [km].
  integer, intent(in) :: nLayerInZone(nZone)  ! Number of layers in each zone.
  integer, intent(in) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(in) :: oValueOfZoneSolid(nZone)  ! Index of the first value in each zone, when counting only solid zones.
  complex(8), intent(out) :: coefQmu(nZone), coefQkappa(nZone), coefQfluid(nZone)
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  real(8), intent(out) :: plm(3, 0:3, nReceiver)
  !:::::::::::::::::::::::::::::::::::::: Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  complex(8), intent(out) :: harmonicsValues(3, -2:2, nReceiver)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Vector harmonics term. The coefficient 1/largeL is not multiplied yet.
  real(8), intent(in) :: t(4*maxNGridSolid-4)
  real(8), intent(in) :: h1x(4*maxNGridSolid-4), h2L(4*maxNGridSolid-4), h2N(4*maxNGridSolid-4)
  real(8), intent(in) :: hUn3y(4*maxNGridSolid-4), hResid3y(4*maxNGridSolid-4), hModL3y(-2:1,maxNGridSolid)
  real(8), intent(in) :: hUn4L(4*maxNGridSolid-4), hResid4L(4*maxNGridSolid-4), hModR4L(-1:2,maxNGridSolid)
  real(8), intent(in) :: hUn4N(4*maxNGridSolid-4), hResid4N(4*maxNGridSolid-4), hModL4N(-2:1,maxNGridSolid)
  real(8), intent(in) :: hUn5y(4*maxNGridSolid-4), hResid5y(4*maxNGridSolid-4), hModR5y(-1:2,maxNGridSolid)
  real(8), intent(in) :: hUn6L(4*maxNGridSolid-4), hResid6L(4*maxNGridSolid-4), hModL6L(-2:1,maxNGridSolid)
  real(8), intent(in) :: hUn6N(4*maxNGridSolid-4), hResid6N(4*maxNGridSolid-4), hModR6N(-1:2,maxNGridSolid)
  real(8), intent(in) :: h7y(4*maxNGridSolid-4), h7z(4*maxNGridSolid-4), h8L(4*maxNGridSolid-4), h8N(4*maxNGridSolid-4)
  real(8), intent(in) :: p1(4*maxNGridFluid-4), p2(4*maxNGridFluid-4), p3(4*maxNGridFluid-4)
  integer, intent(in) :: oPairOfZoneSolid(nZone), oPairOfZoneFluid(nZone)
  !::::::::::::::::::::: Index of the first (iLayer, k', k)-pair in each zone, counted separately for solid and fluid zones.
  complex(8), intent(out) :: a0(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8), intent(out) :: a1(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8), intent(out) :: a2(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8), intent(out) :: a(4, nColumn), aSmall(2, nColumn)  ! Assembled A matrix.
  !::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2] for solid, [m^5/N] for fluid, [10^6 m^2] for boundary elements.
  complex(8), intent(out) :: g_or_c(nColumn), g_or_c_Small(nColumn)
  !::::::::::::::: These hold either vector g ([10^15 N] & [10^9 m^3]) or c ([km] & [GPa]), depending on where in the code it is.
  complex(8), intent(out) :: u(3, nReceiver)  ! Displacement velocity - the unit is [km] in the frequency domain,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: but when converted to the time domain, the unit becomes [km/s].
  integer, intent(in) :: oElementOfZone(nZone)  ! Index of the first (iLayer, k'-gamma', k-gamma)-pair in each zone.
  integer, intent(in) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  integer, intent(in) :: nColumn  ! Total number of columns in the band matrix.
  integer, intent(out) :: nQuasiColumn  ! Total number of columns in the rearranged band matrix.
  complex(8), intent(out) :: anum(4, 4, 10), bnum(4, 4, 10)
  complex(8), intent(out) :: ya(4), yb(4), yc(4), yd(4)
  integer, intent(in) :: oColumnOfSource  ! Index of the first column in the band matrix for the layer with source.
  complex(8), intent(out) :: cwork(16 * maxNGridSolid - 16 + 4 * maxNGridFluid - 4)  ! Working array for matrix computations.
  complex(8), intent(inout) :: z(nColumn), w(nColumn)  ! Working arrays used when solving linear equations.
  real(8), intent(inout) :: eps
  integer, intent(out) :: llog

  integer :: l, m  ! Angular order and azimuthal order of spherical harmonics.
  real(8) :: largeL  ! L = sqrt(l(l+1)).
  real(8) :: largeL2  ! L^2 = l(l+1).
  integer :: lsuf  ! Accuracy threshold of angular order.
  real(8) :: recordAmplitude  ! Maximum amplitude encountered [km], used for angular order cut-off.
  integer :: decayCounter  ! Counter detecting the decay of amplitude, used for angular order cut-off.
  integer :: cutoffColumn  ! Index of column at cut-off depth.
  integer :: i, iSolid, iFluid
  integer :: ir
  integer :: oVS
  integer :: oP, oElement, oColumn

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
  call computeabnum(omega, omegaI, rmax, &
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
    call computeya(anum(:,:,:), bnum(:,:,:), largeL2, gridRadii(iLayerOfSource:), r0, ya(:), yb(:), yc(:), yd(:))

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

end subroutine


!------------------------------------------------------------------------
! Form and solve the linear equation Ac=-g.
!------------------------------------------------------------------------
subroutine formAndSolveEquation(l, m, largeL, iZoneOfSource, iLayerOfSource, oColumnOfSource, r0, mt, ecC0, ecF0, ecL0, &
  ya, yb, yc, yd, rmin, nZone, phaseOfZone, oGridOfZone, oColumnOfZone, coefQmu, coefQkappa, &
  nGrid, gridRadii, nColumn, cutoffColumn, a, aSmall, g_or_c, g_or_c_Small, amplitudeAtColumn, nQuasiColumn, eps, z, w)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: l  ! Angular order.
  integer, intent(in) :: m  ! Azimuthal order.
  real(8), intent(in) :: largeL  ! L = sqrt(l(l+1)).
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(in) :: oColumnOfSource  ! Index of the first column in the band matrix for the layer with source.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecC0, ecF0, ecL0  ! Elastic moduli C, F, and L at source position [10^10 dyn/cm^2 = GPa].
  complex(8), intent(in) :: ya(4), yb(4), yc(4), yd(4)
  real(8), intent(in) :: rmin  ! Minimum radius of region considered [km].
  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(in) :: oGridOfZone(nZone)  ! Index of the first grid point in each zone.
  integer, intent(in) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  complex(8), intent(in) :: coefQmu(nZone), coefQkappa(nZone)
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  integer, intent(in) :: nGrid  ! Total number of grid points (= number of layers + 1).
  real(8), intent(in) :: gridRadii(nGrid)  ! Radius at each grid point [km].
  integer, intent(in) :: nColumn  ! Total number of columns in the band matrix.
  integer, intent(in) :: cutoffColumn  ! Index of column at cut-off depth.
  complex(8), intent(inout) :: a(4, nColumn), aSmall(2, nColumn)  ! Assembled A matrix.
  !::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2] for solid, [m^5/N] for fluid, [10^6 m^2] for boundary elements.
  complex(8), intent(out) :: g_or_c(nColumn), g_or_c_Small(nColumn)
  !::::::::::::::: These hold either vector g ([10^15 N] & [10^9 m^3]) or c ([km] & [GPa]), depending on where in the code it is.
  real(8), intent(inout) :: amplitudeAtColumn(nColumn)  ! Estimate of the amplitude at each column [km].
  integer, intent(out) :: nQuasiColumn  ! Total number of columns in the rearranged band matrix.
  real(8), intent(inout) :: eps
  complex(8), intent(inout) :: z(nColumn), w(nColumn)  ! Working arrays used when solving linear equations.
  integer :: oQuasiColumnOfZoneWithSource  ! Index of the first column in rearraged matrix for the zone with source.
  integer :: iColumnBeforeSource, startColumn, ll(12), lli(12), llj(12)
  integer :: ier  ! Error code from subroutine solving linear equations.

  ! Compute excitation vector g.
  g_or_c(:) = dcmplx(0.d0, 0.d0)
  call calg(l, m, coefQmu(iZoneOfSource), coefQkappa(iZoneOfSource), largeL, ecC0, ecF0, ecL0, &
    ya(:), yb(:), yc(:), yd(:), gridRadii(iLayerOfSource:), r0, mt(:,:), g_or_c(oColumnOfSource:))

  if (l == 0) then
    ! Rearrange matrix elements that are used in l=0 case.
    call rearrangeAForL0(nZone, phaseOfZone(:), oColumnOfZone(:), iZoneOfSource, a(:,:), g_or_c(:), &
      aSmall(:,:), g_or_c_Small(:), oQuasiColumnOfZoneWithSource, nQuasiColumn)

    ! Use all columns, except for those at planet center to impose essential boundary condition u=0.
    if (rmin > 0.d0) then
      startColumn = 1
    else
      startColumn = 2
    end if

    iColumnBeforeSource = oQuasiColumnOfZoneWithSource + (iLayerOfSource - oGridOfZone(iZoneOfSource)) - 1

    ! Decompose matrix A.
    call decomposeAByGauss(aSmall(:, startColumn:), 1, nQuasiColumn - startColumn + 1, 1, eps, z(startColumn:), &
      w(startColumn:), ll, lli, llj, ier)
    ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
    call solveSurfaceCAfterGauss(aSmall(:, startColumn:), g_or_c_Small(startColumn:), 1, nQuasiColumn - startColumn + 1, &
      iColumnBeforeSource - startColumn + 1, z(startColumn:))

  else if (mod(l, 100) == 0) then
    ! Once in a while, compute for all grids to decide the cut-off depth.
    ! CAUTION: In this case, all values of g_or_c(:) are computed.

    ! Use all columns, except for those at planet center to impose essential boundary condition u=0.
    if (rmin > 0.d0) then
      startColumn = 1
    else if (phaseOfZone(1) == 1) then
      ! solid
      startColumn = 3
    else
      ! liquid
      startColumn = 2
    end if

    ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
    ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
    if (m == -2 .or. m == -l) then
      call decomposeAByGauss(a(:, startColumn:), 3, nColumn - startColumn + 1, 6, eps, z(startColumn:), &
        w(startColumn:), ll, lli, llj, ier)
    end if
    ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
    call solveWholeCAfterGauss(a(:, startColumn:), g_or_c(startColumn:), 3, nColumn - startColumn + 1, z(startColumn:))

    ! Accumulate the absolute values of expansion coefficent c for all m's at each grid point.
    !  This is to be used as an estimate of the amplitude at each depth when deciding the cut-off depth.
    amplitudeAtColumn(1:nColumn) = amplitudeAtColumn(1:nColumn) + abs(g_or_c(1:nColumn))

  else
    ! Otherwise, compute for just the grids above the cut-off depth.
    ! CAUTION: In this case, only g_or_c(nColumn-1:nColumn) is computed.
    !   Other values of g_or_c(:nColumn-2) still hold values of g!!!

    ! Consider cutoff + exclude columns at planet center to impose essential boundary condition u=0.
    if (rmin > 0.d0 .or. cutoffColumn > 1) then
      startColumn = cutoffColumn
    else if (phaseOfZone(1) == 1) then
      ! solid
      startColumn = 3
    else
      ! liquid
      startColumn = 2
    end if

    iColumnBeforeSource = oColumnOfZone(iZoneOfSource) + 2 * (iLayerOfSource - oGridOfZone(iZoneOfSource)) - 2

    ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
    ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
    if (m == -2 .or. m == -l) then
      call decomposeAByGauss(a(:, startColumn:), 3, nColumn - startColumn + 1, 6, eps, z(startColumn:), &
        w(startColumn:), ll, lli, llj, ier)
    end if
    ! Solve Ac=g (i.e. (omega^2 T - H) c = -g).
    call solveSurfaceCAfterGauss(a(:, startColumn:), g_or_c(startColumn:), 3, nColumn - startColumn + 1, &
      iColumnBeforeSource - startColumn + 1, z(startColumn:))

  end if

end subroutine
