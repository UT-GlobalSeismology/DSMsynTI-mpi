!------------------------------------------------------------------------
! Computes matrix elements common for all omega, l, and m.
!------------------------------------------------------------------------
subroutine computeMatrixElements(maxNGrid, tlen, re, imaxFixed, r0, &
  nZone, rmin, rmax, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, &
  kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oRowOfZone, &
  iZoneOfSource, iLayerOfSource, oRowOfSource, gridRadiiForSource, &
  nValue, valuedRadii, rhoValues, ecLValues, ecNValues, rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, ecL0, &
  t, h1, h2, h3, h4, gt, gh1, gh2, gh3, gh4, work)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: maxNGrid! Maximum number of grid points.
  real(8), intent(in) :: tlen  ! Time length [s].
  real(8), intent(in) :: re  ! Desired relative error due to vertical gridding.
  integer, intent(in) :: imaxFixed  ! Index of maximum frequency.
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
  integer, intent(out) :: oRowOfZone(nZone)  ! Index of the first row in the vector of (iLayer, k', k)-pairs in each zone.
  integer, intent(out) :: iZoneOfSource  ! Which zone the source is in.
  integer, intent(out) :: iLayerOfSource  ! Which layer the source is in.
  integer, intent(out) :: oRowOfSource  ! Index of the first row in the vector of (iLayer, k', k)-pairs for the source layer.
  real(8), intent(out) :: gridRadiiForSource(3)  ! Radii to use for source-related computations [km].
  integer, intent(out) :: nValue  ! Total number of values for each variable.
  real(8), intent(out) :: valuedRadii(maxNGrid+nZone-1)  ! Radii corresponding to each variable value [km].
  real(8), intent(out) :: rhoValues(maxNGrid+nZone-1), ecLValues(maxNGrid+nZone-1), ecNValues(maxNGrid+nZone-1)
  !::::::::::::::::::::::::: Values of rho [g/cm^3], L, and N [10^10 dyn/cm^2 = GPa] at each point (with 2 values at boundaries).
  real(8), intent(out) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::  Values of rho [g/cm^3], L, and N [10^10 dyn/cm^2 = GPa] at each point.
  real(8), intent(out) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  real(8), intent(out) :: t(4*maxNGrid-4), h1(4*maxNGrid-4), h2(4*maxNGrid-4), h3(4*maxNGrid-4), h4(4*maxNGrid-4)
  real(8), intent(out) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  real(8), intent(out) :: work(4*maxNGrid-4)  ! Working matrix.
  integer :: i


  ! ******************* Computing parameters *******************
  ! Design the number and position of grid points.
  call computeKz(nZone, rminOfZone(:), rmaxOfZone(:), vsvPolynomials(:,:), rmax, imaxFixed, 1, tlen, kzAtZone(:))
  call computeGridRadii(nZone, kzAtZone(:), rminOfZone(:), rmaxOfZone(:), rmin, re, nGrid, nLayerInZone(:), gridRadii(:))
  if (nGrid > maxNGrid) stop 'The number of grid points is too large.'

  ! Compute the first indices in each zone.
  call computeFirstIndices(nZone, nLayerInZone(:), oGridOfZone(:), oValueOfZone(:), oRowOfZone(:))

  ! Compute the source position.
  call computeSourcePosition(nGrid, rmaxOfZone(:), gridRadii(:), r0, iZoneOfSource, iLayerOfSource, oRowOfSource)

  ! Design grids for source computations.
  call computeSourceGrid(gridRadii(:), r0, iLayerOfSource, gridRadiiForSource(:))


  ! ******************* Computing the matrix elements *******************
  ! Compute variable values at grid points.
  call computeStructureValues(nZone, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), nLayerInZone(:), &
    gridRadii(:), nValue, valuedRadii(:), rhoValues(:), ecLValues(:), ecNValues(:))
  call computeSourceStructureValues(iZoneOfSource, rmax, rhoPolynomials(:,:), vsvPolynomials(:,:), vshPolynomials(:,:), &
    gridRadiiForSource(:), rhoValuesForSource(:), ecLValuesForSource(:), ecNValuesForSource(:), ecL0)

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

end subroutine
