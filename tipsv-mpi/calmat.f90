
!----------------------------------------------------------------------------------------------------------------------------
! Computing the transpose of a tridiagonal matrix in a certain zone.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine transposeMatrix(nLayerInZoneI, matIn, matOut)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: matIn(4*nLayerInZoneI)  ! Input tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  real(8), intent(out) :: matOut(4*nLayerInZoneI)  ! Transposed tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  integer :: i

  do i = 4, 4 * nLayerInZoneI, 4
    matOut(i - 3) = matIn(i - 3)
    matOut(i - 2) = matIn(i - 1)
    matOut(i - 1) = matIn(i - 2)
    matOut(i) = matIn(i)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the step-wise part of the unmodified operator for a certain zone.
! (See eqs. 3.44 and 3.45 of Geller & Takeuchi 1995; eqs. 11 and 12 of Takeuchi et al. 1996.)
! Note that in the above papers, the average of ec*r within each layer is taken,
!  but in this program, its value at the grid point is used.
!      / -3D0-D1  3D0+D1           0     ..    \
! 1/8 |  -D0-3D1  D0+3D1-3D1-D2  3D1+D2  0  ..  |
!      \    :            :           :   ..    /
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeStepH(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, h5)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: ecValuesInZoneI(nLayerInZoneI+1)  ! Modulus values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(out) :: h5(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  integer :: iLayer

  do iLayer = 1, nLayerInZoneI
    h5(4 * iLayer - 3) = - 3.d0 / 8.d0 * ecValuesInZoneI(iLayer) * valuedRadiiInZoneI(iLayer) &
      - 1.d0 / 8.d0 * ecValuesInZoneI(iLayer + 1) * valuedRadiiInZoneI(iLayer + 1)
    h5(4 * iLayer - 2) = - h5(4 * iLayer - 3)
    h5(4 * iLayer - 1) = - 1.d0 / 8.d0 * ecValuesInZoneI(iLayer) * valuedRadiiInZoneI(iLayer) &
      - 3.d0 / 8.d0 * ecValuesInZoneI(iLayer + 1) * valuedRadiiInZoneI(iLayer + 1)
    h5(4 * iLayer) = - h5(4 * iLayer - 1)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Subtracting two tridiagonal matrices for a certain zone.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine subtractMatrix(nLayerInZoneI, mat1, mat2, difference)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: mat1(4*nLayerInZoneI), mat2(4*nLayerInZoneI)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::: Input tridiagonal matrices, stored for each (iLayer, k', k)-pair.
  real(8), intent(out) :: difference(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  integer :: i

  do i = 1, 4 * nLayerInZoneI
    difference(i) = mat1(i) - mat2(i)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the step-wise part of the modified matrix for a certain zone.
! (See eq. 16 of Takeuchi et al. 1996.)
! Note that in the paper, the average of ec*r within each layer is taken,
!  but in this program, its value at the grid point is used.
!       / -7D0  8D0  -D0   0     ..    \
!      |  -5D1  -3D1  9D1  -D1  0  ..  |
! 1/12 |    :     :    :    :    ..    |
!      |   ..  0  -5DN-1 -3DN-1 8DN-1  |
!       \           ..  0  -5DN  5DN  /
! The result is a band matrix, stored by (offset from diagonal, row number) for each zone.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeModifiedHR(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, hm1)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: ecValuesInZoneI(nLayerInZoneI+1)  ! Modulus values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(inout) :: hm1(-1:2, nLayerInZoneI+1)  ! Resulting band matrix, stored by (offset from diagonal, row number)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  integer :: iGrid

  iGrid = 1
  hm1(0, iGrid) = -7.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm1(1, iGrid) = 8.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm1(2, iGrid) = -1.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)

  do iGrid = 2, nLayerInZoneI - 1
    hm1(-1, iGrid) = -5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
    hm1(0, iGrid) = -3.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
    hm1(1, iGrid) = 9.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
    hm1(2, iGrid) = -1.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  end do

  iGrid = nLayerInZoneI
  hm1(-1, iGrid) = -5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm1(0, iGrid) = -3.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm1(1, iGrid) = 8.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)

  iGrid = nLayerInZoneI + 1
  hm1(-1, iGrid) = -5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm1(0, iGrid) = 5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the step-wise part of the modified matrix for a certain zone.
! (See eq. 17 of Takeuchi et al. 1996.)
! Note that in the paper, the average of ec*r within each layer is taken,
!  but in this program, its value at the grid point is used.
!       / -5D0  5D0   0   ..           \
!      |  -8D1  3D1  5D1  0   ..        |
! 1/12 |    :    :    :   :   ..        |
!      |  .. 0  DN-1 -9DN-1 3DN-1 5DN-1 |
!       \    ..   0    DN   -8DN   7DN /
! The result is a band matrix, stored by (offset from diagonal, row number) for each zone.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeModifiedHL(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, hm2)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: ecValuesInZoneI(nLayerInZoneI+1)  ! Modulus values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(inout) :: hm2(-2:1, nLayerInZoneI+1)  ! Resulting band matrix, stored by (offset from diagonal, row number)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  integer :: iGrid

  iGrid = 1
  hm2(0, iGrid) = -5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm2(1, iGrid) = 5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)

  iGrid = 2
  hm2(-1, iGrid) = -8.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm2(0, iGrid) = 3.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm2(1, iGrid) = 5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)

  do iGrid = 3, nLayerInZoneI
    hm2(-2, iGrid) = 1.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
    hm2(-1, iGrid) = -9.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
    hm2(0, iGrid) = 3.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
    hm2(1, iGrid) = 5.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  end do

  iGrid = nLayerInZoneI + 1
  hm2(-2, iGrid) = 1.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm2(-1, iGrid) = -8.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)
  hm2(0, iGrid) = 7.d0 / 12.d0 * ecValuesInZoneI(iGrid) * valuedRadiiInZoneI(iGrid)

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the transpose of a band matrix stored by (offset from diagonal, row number) in a certain zone.
!----------------------------------------------------------------------------------------------------------------------------
subroutine transposeMatrixMod(nLayerInZoneI, pMin, pMax, hModIn, hModOut)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: pMin, pMax  ! Minimum and maximum value of offset from diagonal for the input band matrix.
  real(8), intent(in) :: hModIn(pMin:pMax, nLayerInZoneI+1)  ! Input band matrix, stored by (offset from diagonal, row number).
  real(8), intent(out) :: hModOut(-pMax:-pMin, nLayerInZoneI+1)  ! Resulting band matrix,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: stored by (offset from diagonal, row number).
  integer :: pOut, iOut, pIn, iIn

  do pOut = -pMax, -pMin
    do iOut = 1, nLayerInZoneI + 1
      pIn = -pOut
      iIn = pOut + iOut
      if (1 <= iIn .and. iIn <= nLayerInZoneI + 1) then
        hModOut(pOut, iOut) = hModIn(pIn, iIn)
      end if
    end do
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes the part of H_(k'1k1) and H_(k'2k2) without largeL coefficients. (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a block tridiagonal matrix, stored for each (iLayer, k'-gamma', k-gamma)-pair.
! Only elements in the upper triangle of A is computed.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA0Solid(nLayerInZoneI, omega, omegaI, t, h1x, h2L, h2N, h3ay, h4aL, h4aN, h5ay, h6aL, h6aN, h7y, h7z, h8L, h8N, &
  coefQmu, coefQkappa, a0Tmp)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: t(4*nLayerInZoneI)  ! T matrix stored for each (iLayer, k', k)-pair [10^12 kg].
  real(8), intent(in) :: h1x(4*nLayerInZoneI), h2L(4*nLayerInZoneI), h2N(4*nLayerInZoneI)
  real(8), intent(in) :: h3ay(4*nLayerInZoneI), h4aL(4*nLayerInZoneI), h4aN(4*nLayerInZoneI)
  real(8), intent(in) :: h5ay(4*nLayerInZoneI), h6aL(4*nLayerInZoneI), h6aN(4*nLayerInZoneI)
  real(8), intent(in) :: h7y(4*nLayerInZoneI), h7z(4*nLayerInZoneI), h8L(4*nLayerInZoneI), h8N(4*nLayerInZoneI)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(in) :: coefQmu, coefQkappa
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  complex(8), intent(out) :: a0Tmp(16*nLayerInZoneI)  ! Resulting block tridiagonal matrix,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: stored for each (iLayer, k'-gamma', k-gamma)-pair [10^12 kg/s^2].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  complex(8) :: hh0
  integer :: iElement, iBlock

  ! Initialize all elements, since they are referenced when overlapping.
  a0Tmp(:) = dcmplx(0.d0)

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  do iBlock = 1, 4 * nLayerInZoneI
    if (mod(iBlock, 4) /= 3) then

      ! Compute the part of H_(k'1k1) without largeL coefficients.
      iElement = iBlock * 4 - 3
      hh0 = coefQkappa * dcmplx(4.d0 * h1x(iBlock)) + coefQmu * dcmplx(16.d0/3.d0 * h2N(iBlock)) &  ! 4 * [(A-4N/3) + 4N/3] = 4*A
        - coefQmu * dcmplx(4.d0 * h2N(iBlock)) &  ! -4*N
        + coefQkappa * dcmplx(2.d0 * (h3ay(iBlock) + h5ay(iBlock))) - coefQmu * dcmplx(4.d0/3.d0 * (h4aN(iBlock) + h6aN(iBlock))) &
      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 2 * [(F+2N/3) - 2N/3] = 2*F
        + coefQkappa * dcmplx(3.d0 * h7z(iBlock) - 2.d0 * h7y(iBlock)) + coefQmu * dcmplx(4.d0/3.d0 * h8N(iBlock))
      !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 3 * (C+2F)/3 - 2 * (F+2N/3) + 4N/3 = C
      a0Tmp(iElement) = omegaDamped2 * dcmplx(t(iBlock)) - hh0

      ! Compute the part of H_(k'2k2) without largeL coefficients.
      iElement = iBlock * 4
      hh0 = coefQmu * dcmplx(h2L(iBlock)) &  ! L
        - coefQmu * dcmplx(2.d0 * h2N(iBlock)) &  ! -2*N
        - coefQmu * dcmplx(h4aL(iBlock) + h6aL(iBlock)) &  ! -L
        + coefQmu * dcmplx(h8L(iBlock))  ! L
      a0Tmp(iElement) = omegaDamped2 * dcmplx(t(iBlock)) - hh0

    end if
  end do

end subroutine



!----------------------------------------------------------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes H_(k'1k2) and H_(k'2k1) (largeL is not multiplied). (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a block tridiagonal matrix, stored for each (iLayer, k'-gamma', k-gamma)-pair.
! Only elements in the upper triangle of A is computed.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA1Solid(nLayerInZoneI, h1x, h2L, h2N, h3y, h4L, h4N, h5y, h6L, h6N, coefQmu, coefQkappa, a1Tmp)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: h1x(4*nLayerInZoneI), h2L(4*nLayerInZoneI), h2N(4*nLayerInZoneI)
  real(8), intent(in) :: h3y(4*nLayerInZoneI), h4L(4*nLayerInZoneI), h4N(4*nLayerInZoneI)
  real(8), intent(in) :: h5y(4*nLayerInZoneI), h6L(4*nLayerInZoneI), h6N(4*nLayerInZoneI)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(in) :: coefQmu, coefQkappa
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  complex(8), intent(out) :: a1Tmp(16*nLayerInZoneI)  ! Resulting block tridiagonal matrix,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: stored for each (iLayer, k'-gamma', k-gamma)-pair [10^12 kg/s^2].
  complex(8) :: hh1
  integer :: iElement, iBlock

  ! Initialize all elements, since they are referenced when overlapping.
  a1Tmp(:) = dcmplx(0.d0)

  do iBlock = 1, 4 * nLayerInZoneI

    ! Compute H_(k'1k2).
    iElement = iBlock * 4 - 2
    if (mod(iBlock, 4) /= 3) then
      hh1 = -coefQkappa * dcmplx(2.d0 * h1x(iBlock)) - coefQmu * dcmplx(8.d0/3.d0 * h2N(iBlock)) &  ! -2 * [(A-4N/3) + 4N/3] = -2*A
        - coefQmu * dcmplx(h2L(iBlock)) &  ! -L
        + coefQmu * dcmplx(2.d0 * h2N(iBlock)) &  ! 2N
        - coefQkappa * dcmplx(h3y(iBlock)) + coefQmu * dcmplx(2.d0/3.d0 * h4N(iBlock)) &  ! -[(F+2N/3) - 2N/3] = -F
        + coefQmu * dcmplx(h6L(iBlock))  ! L
      a1Tmp(iElement) = -hh1
    end if

    ! Compute H_(k'2k1).
    iElement = iBlock * 4 - 1
    if (mod(iBlock, 4) == 2) then
      hh1 = -coefQkappa * dcmplx(2.d0 * h1x(iBlock)) - coefQmu * dcmplx(8.d0/3.d0 * h2N(iBlock)) &  ! -2 * [(A-4N/3) + 4N/3] = -2*A
        - coefQmu * dcmplx(h2L(iBlock)) &  ! -L
        + coefQmu * dcmplx(2.d0 * h2N(iBlock)) &  ! 2N
        - coefQkappa * dcmplx(h5y(iBlock)) + coefQmu * dcmplx(2.d0/3.d0 * h6N(iBlock)) &  ! -[(F+2N/3) - 2N/3] = -F
        + coefQmu * dcmplx(h4L(iBlock))  ! L
      a1Tmp(iElement) = -hh1
    end if

  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes the part of H_(k'1k1) and H_(k'2k2) with coefficient largeL^2. (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a block tridiagonal matrix, stored for each (iLayer, k'-gamma', k-gamma)-pair.
! Only elements in the upper triangle of A is computed.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA2Solid(nLayerInZoneI, h1x, h2L, h2N, coefQmu, coefQkappa, a2Tmp)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: h1x(4*nLayerInZoneI), h2L(4*nLayerInZoneI), h2N(4*nLayerInZoneI)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(in) :: coefQmu, coefQkappa
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  complex(8), intent(out) :: a2Tmp(16*nLayerInZoneI)  ! Resulting block tridiagonal matrix,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: stored for each (iLayer, k'-gamma', k-gamma)-pair [10^12 kg/s^2].
  complex(8) :: hh2
  integer :: iElement, iBlock

  ! Initialize all elements, since they are referenced when overlapping.
  a2Tmp(:) = dcmplx(0.d0)

  do iBlock = 1, 4 * nLayerInZoneI
    if (mod(iBlock, 4) /= 3) then

      ! Compute the part of H_(k'1k1) with coefficient largeL^2.
      iElement = iBlock * 4 - 3
      hh2 = coefQmu * dcmplx(h2L(iBlock))  ! L
      a2Tmp(iElement) = -hh2

      ! Compute the part of H_(k'2k2) with coefficient largeL^2.
      iElement = iBlock * 4
      hh2 = coefQkappa * dcmplx(h1x(iBlock)) + coefQmu * dcmplx(4.d0/3.d0 * h2N(iBlock))  ! (A-4N/3) + 4N/3 = A
      a2Tmp(iElement) = -hh2

    end if
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Overlapping the coefficient matrix elements for a certain zone in the solid part.
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine overlapASolid(nLayerInZoneI, aTmp, aOut)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  complex(8), intent(in) :: aTmp(16*nLayerInZoneI)  ! Input block tridiagonal matrix,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: stored for each (iLayer, k'-gamma', k-gamma)-pair [10^12 kg/s^2].
  complex(8), intent(inout) :: aOut(4, 2*nLayerInZoneI+2)  ! Upper band of the overlapped matrix [10^12 kg/s^2].
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Should be initialized with 0s beforehand.
  integer :: iGrid

  do iGrid = 1, nLayerInZoneI
    ! (i,i)-block
    if (iGrid == 1) then
      ! This overlaps with previous zone if phase is same.
      aOut(4, 2 * iGrid - 1) = aOut(4, 2 * iGrid - 1) + aTmp(16 * iGrid - 15)
      aOut(3, 2 * iGrid) = aOut(3, 2 * iGrid) + aTmp(16 * iGrid - 14)
      aOut(4, 2 * iGrid) = aOut(4, 2 * iGrid) + aTmp(16 * iGrid - 12)
    else
      aOut(4, 2 * iGrid - 1) = aTmp(16 * iGrid - 19) + aTmp(16 * iGrid - 15)
      aOut(3, 2 * iGrid) = aTmp(16 * iGrid - 18) + aTmp(16 * iGrid - 14)
      aOut(4, 2 * iGrid) = aTmp(16 * iGrid - 16) + aTmp(16 * iGrid - 12)
    end if
    ! (i,i+1)-block
    aOut(2, 2 * iGrid + 1) = aTmp(16 * iGrid - 11)
    aOut(1, 2 * iGrid + 2) = aTmp(16 * iGrid - 10)
    aOut(3, 2 * iGrid + 1) = aTmp(16 * iGrid - 9)
    aOut(2, 2 * iGrid + 2) = aTmp(16 * iGrid - 8)
  end do
  ! (N,N)-block
  iGrid = nLayerInZoneI + 1
  aOut(4, 2 * iGrid - 1) = aTmp(16 * iGrid - 19)
  aOut(3, 2 * iGrid) = aTmp(16 * iGrid - 18)
  aOut(4, 2 * iGrid) = aTmp(16 * iGrid - 16)

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Adding the modified part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes H_(k'1k2) and H_(k'2k1) (largeL is not multiplied). (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine addModifiedHToA1(nLayerInZoneI, coefQmu, coefQkappa, hModL3y, hModR4L, hModL4N, hModR5y, hModL6L, hModR6N, a1)
  !--------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  complex(8), intent(in) :: coefQmu, coefQkappa
  !::::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  real(8), intent(in) :: hModL3y(-2:1, nLayerInZoneI+1), hModR4L(-1:2, nLayerInZoneI+1), hModL4N(-2:1, nLayerInZoneI+1)
  real(8), intent(in) :: hModR5y(-1:2, nLayerInZoneI+1), hModL6L(-2:1, nLayerInZoneI+1), hModR6N(-1:2, nLayerInZoneI+1)
  !::::::::::::::::::::::::::::::: Modified operator band matrices, stored by (offset from diagonal, row number) [10^12 kg/s^2].
  complex(8), intent(inout) :: a1(4, 2*nLayerInZoneI+2)  ! Resulting matrix, containing upper band elements.
  integer :: q, iColumn, p, iGrid

  ! Add modified H_(k'1k2) to A.
  do iColumn = 2, 2 * nLayerInZoneI + 2, 2
    do q = 1, 3, 2
      ! Skip if out of bounds.
      if (iColumn == 2 .and. q == 1) cycle
      ! Find indices of modified matrix corresponding to this element.
      p = (-q + 3) / 2
      iGrid = (q + iColumn - 3) / 2
      ! -H_(k'1k2)mod = -( -L (I[F] - I[L]) ) = L (I[F] - I[L])  (L is not multiplied here.)
      a1(q, iColumn) = a1(q, iColumn) &
        + (coefQkappa * dcmplx(hModL3y(p, iGrid)) - coefQmu * dcmplx(2.d0 / 3.d0 * hModL4N(p, iGrid)) &  ! -[(F+2N/3) - 2N/3] = -F
        - coefQmu * dcmplx(hModL6L(p, iGrid)))  ! L
    end do
  end do

  ! Add modified H_(k'2k1) to A.
  do iColumn = 3, 2 * nLayerInZoneI + 1, 2
    do q = 1, 3, 2
      ! Skip if out of bounds.
      if (iColumn == 3 .and. q == 1) cycle
      ! Find indices of modified matrix corresponding to this element.
      p = (-q + 5) / 2
      iGrid = (q + iColumn - 4) / 2
      ! -H_(k'2k1)mod = -( -L (I[F]mod - I[L]mod) ) = L (I[F]mod - I[L]mod)  (L is not multiplied here.)
      a1(q, iColumn) = a1(q, iColumn) &
        + (coefQkappa * dcmplx(hModR5y(p, iGrid)) - coefQmu * dcmplx(2.d0 / 3.d0 * hModR6N(p, iGrid)) &  ! -[(F+2N/3) - 2N/3] = -F
        - coefQmu * dcmplx(hModR4L(p, iGrid)))  ! L
    end do
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the fluid part.
! This computes the part without largeL coefficients. (See eqs. 2 & 21 of Kawai et al. 2006.)
! Note that a factor of 1/(omega^2) is multiplied.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA0Fluid(nLayerInZoneI, omega, omegaI, p1, p3, coefQfluid, a0Tmp)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: p1(4*nLayerInZoneI), p3(4*nLayerInZoneI)
  !:::::::::::::::::::::::::::::::: Parts of T and H matrix stored for each (iLayer, k', k)-pair. I^F2: [m^5/N], others: [m^4/kg].
  complex(8), intent(in) :: coefQfluid  ! Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  complex(8), intent(out) :: a0Tmp(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [m^5/N].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  ! Compute the part of H_(k'k) without largeL coefficients.
  do i = 1, 4 * nLayerInZoneI
    a0Tmp(i) = -dcmplx(p1(i)) / omegaDamped2 + coefQfluid * dcmplx(p3(i))
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the fluid part.
! This computes the part with coefficient largeL^2. (See eqs. 2 & 21 of Kawai et al. 2006.)
! Note that a factor of 1/(omega^2) is multiplied.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA2Fluid(nLayerInZoneI, omega, omegaI, p2, a2Tmp)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: p2(4*nLayerInZoneI)  ! Parts of H matrix stored for each (iLayer, k', k)-pair [m^4/kg].
  complex(8), intent(out) :: a2Tmp(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [m^5/N].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  ! Compute the part of H_(k'k) with coefficient largeL^2.
  do i = 1, 4 * nLayerInZoneI
    a2Tmp(i) = -dcmplx(p2(i)) / omegaDamped2
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Overlapping the coefficient matrix elements for a certain zone in the fluid part.
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine overlapAFluid(nLayerInZoneI, aTmp, aOut)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  complex(8), intent(in) :: aTmp(4*nLayerInZoneI)  ! Input tridiagonal matrix, stored for each (iLayer, k', k)-pair [m^5/N].
  complex(8), intent(inout) :: aOut(4, 2*nLayerInZoneI+2)  ! Upper band of the overlapped matrix [m^5/N].
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Should be initialized with 0s beforehand.
  integer :: iGrid

  do iGrid = 1, nLayerInZoneI
    ! (i,i)-component
    if (iGrid == 1) then
      ! This overlaps with previous zone if phase is same.
      aOut(4, 1) = aOut(4, 1) + aTmp(1)
    else
      aOut(4, iGrid) = aTmp(4 * iGrid - 4) + aTmp(4 * iGrid - 3)
    end if
    ! (i,i+1)-component
    aOut(3, iGrid + 1) = aTmp(4 * iGrid - 2)
  end do
  ! (N,N)-component
  aOut(4, nLayerInZoneI + 1) = aTmp(4 * nLayerInZoneI)

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Assembling the coefficient matrix 'A' for the whole region from several parts.
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine assembleAWhole(nZone, phaseOfZone, oColumnOfZone, largeL2, a0, a1, a2, a)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(in) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  real(8), intent(in) :: largeL2  ! L^2 = l(l+1).
  complex(8), intent(in) :: a0(4, *), a1(4, *), a2(4, *)  ! Parts of the A matrix, containing upper band elements.
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: For solid, [10^12 kg/s^2]. For fluid, [m^5/N].
  complex(8), intent(out) :: a(4, *)  ! Assembled A matrix, containing upper band elements.
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: For solid, [10^12 kg/s^2]. For fluid, [m^5/N].
  integer :: iZone, q, iColumn, iStart, iEnd
  complex(8) :: largeL2c, largeLc

  ! These are changed to complex beforehand to reduce computation amount inside loop.
  largeL2c = dcmplx(largeL2)
  largeLc = dcmplx(sqrt(largeL2))

  do iZone = 1, nZone
    iStart = oColumnOfZone(iZone)
    iEnd = oColumnOfZone(iZone + 1) - 1

    if (phaseOfZone(iZone) == 1) then
      ! solid
      do iColumn = iStart, iEnd
        do q = 2, 4, 2
          a(q, iColumn) = a0(q, iColumn) + largeL2c * a2(q, iColumn)
        end do
        do q = 1, 3, 2
          a(q, iColumn) = largeLc * a1(q, iColumn)
        end do
      end do

    else
      ! fluid
      do iColumn = iStart, iEnd
        do q = 3, 4
          a(q, iColumn) = a0(q, iColumn) + largeL2c * a2(q, iColumn)
        end do
      end do

    end if
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Setting the boundary condition elements to the coefficient matrix 'A', with elements
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
! Note that although matrix A should be (omega^2 T - H + omega R), a factor of 1/omega is applied to this part of the matrix.
!----------------------------------------------------------------------------------------------------------------------------
subroutine setBoundaryConditions(nZone, rmaxOfZone, phaseOfZone, oColumnOfZone, a)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmaxOfZone(nZone)  ! Upper radius of each zone [km].
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(in) :: oColumnOfZone(nZone + 1)  ! Index of the first column in the band matrix for each zone.
  complex(8), intent(inout) :: a(4,*)  ! A matrix, containing upper band elements. [10^6 m^2] for boundary elements.
  integer :: iZone

  do iZone = 1, nZone - 1
    if (phaseOfZone(iZone) == 1 .and. phaseOfZone(iZone + 1) == 2) then
      a(2, oColumnOfZone(iZone + 1)) = dcmplx(rmaxOfZone(iZone) ** 2)
    endif
    if (phaseOfZone(iZone) == 2 .and. phaseOfZone(iZone + 1) == 1) then
      a(3, oColumnOfZone(iZone + 1)) = -dcmplx(rmaxOfZone(iZone) ** 2)
    endif
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Rearranging the elements of matrix A for the case of l=0.
! When l=0, only the vertical component can exist, so we extract the matrix elements for gamma'=gamma=1.
!----------------------------------------------------------------------------------------------------------------------------
subroutine rearrangeAForL0(nZone, phaseOfZone, oColumnOfZone, iZoneOfSource, aIn, gIn, &
  aSmall, gSmall, oQuasiColumnOfZoneWithSource, nQuasiColumn)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(in) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  integer, intent(in) :: iZoneOfSource  ! Which zone the source is in.
  complex(8), intent(in) :: aIn(4,*), gIn(*)  ! Input A matrix and g vector.
  ! A: [10^12 kg/s^2] for solid, [m^5/N] for fluid, [10^6 m^2] for boundary elements. g: [10^15 N] for solid, [10^9 m^3] for fluid.
  complex(8), intent(out) :: aSmall(2,*), gSmall(*)  ! Rearranged A matrix and g vector.
  ! A: [10^12 kg/s^2] for solid, [m^5/N] for fluid, [10^6 m^2] for boundary elements. g: [10^15 N] for solid, [10^9 m^3] for fluid.
  integer, intent(out) :: oQuasiColumnOfZoneWithSource  ! Index of the first column in rearraged matrix for the zone with source.
  integer, intent(out) :: nQuasiColumn  ! Number of columns in rearranged A matrix.
  integer :: iZone, firstQuasiColumn, lastQuasiColumn, iColumn, iQuasiColumn

  lastQuasiColumn = 0
  do iZone = 1, nZone

    if (phaseOfZone(iZone) == 1) then
      ! solid
      if (iZone == iZoneOfSource) then
        oQuasiColumnOfZoneWithSource = lastQuasiColumn + 1
      endif
      firstQuasiColumn = lastQuasiColumn + 1
      lastQuasiColumn = firstQuasiColumn + (oColumnOfZone(iZone + 1) - oColumnOfZone(iZone)) / 2 - 1

      gSmall(firstQuasiColumn) = gIn(oColumnOfZone(iZone))
      ! When previous zone exists, get the element above the top corner of the block for this zone.
      if (iZone /= 1) then
        if (phaseOfZone(iZone - 1) == 1) then
          ! If previous zone is solid, get gamma'=gamma=1 element as usual.
          aSmall(1, firstQuasiColumn) = aIn(2, oColumnOfZone(iZone))
        else
          ! If previous zone is fluid, get the continuity condition at subdiagonal.
          aSmall(1, firstQuasiColumn) = aIn(3, oColumnOfZone(iZone))
        endif
      endif
      ! First diagonal element.
      aSmall(2, firstQuasiColumn) = aIn(4, oColumnOfZone(iZone))

      ! Loop for the rest of the columns in zone.
      do iQuasiColumn = firstQuasiColumn + 1, lastQuasiColumn
        iColumn = oColumnOfZone(iZone) + 2 * (iQuasiColumn - firstQuasiColumn)
        gSmall(iQuasiColumn) = gIn(iColumn)
        ! Get gamma'=gamma=1 elements.
        aSmall(1, iQuasiColumn) = aIn(2, iColumn)
        aSmall(2, iQuasiColumn) = aIn(4, iColumn)
      end do

    else
      ! fluid
      firstQuasiColumn = lastQuasiColumn + 1
      lastQuasiColumn = firstQuasiColumn + (oColumnOfZone(iZone + 1) - oColumnOfZone(iZone)) - 1

      gSmall(firstQuasiColumn) = gIn(oColumnOfZone(iZone))
      ! When previous zone exists, get the element above the top corner of the block for this zone.
      if (iZone /= 1) then
        if (phaseOfZone(iZone - 1) == 1) then
          ! If previous zone is solid, get the continuity condition at subdiagonal.
          aSmall(1, firstQuasiColumn) = aIn(2, oColumnOfZone(iZone))
        else
          ! If previous zone is fluid, get subdiagonal element as usual.
          aSmall(1, firstQuasiColumn) = aIn(3, oColumnOfZone(iZone))
        endif
      endif
      ! First diagonal element.
      aSmall(2, firstQuasiColumn) = aIn(4, oColumnOfZone(iZone))

      ! Loop for the rest of the columns in zone.
      do iQuasiColumn = firstQuasiColumn + 1, lastQuasiColumn
        iColumn = oColumnOfZone(iZone) + (iQuasiColumn - firstQuasiColumn)
        gSmall(iQuasiColumn) = gIn(iColumn)
        ! Get subdiagonal and diagonal elements.
        aSmall(1, iQuasiColumn) = aIn(3, iColumn)
        aSmall(2, iQuasiColumn) = aIn(4, iColumn)
      end do

    endif
  end do
  nQuasiColumn = lastQuasiColumn

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computes numerical values for matrices anum and bnum.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeabnum(omega, omegaI, rmax, &
  rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, &
  gridRadiiForSource, r0, coefQmu, coefQkappa, anum, bnum)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: rmax
  real(8), intent(in) :: rhoPolynomials(4), vpvPolynomials(4), vphPolynomials(4)
  real(8), intent(in) :: vsvPolynomials(4), vshPolynomials(4), etaPolynomials(4)
  !:::::::::::::::::::::::::::::::::::::::::: Polynomial functions of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta structure.
  real(8), intent(in) :: gridRadiiForSource(2), r0
  complex(8), intent(in) :: coefQmu, coefQkappa  ! Coefficients to multiply to elastic moduli for attenuation.
  complex(8), intent(out) :: anum(4, 4, 10), bnum(4, 4, 10)
  integer :: i
  real(8) :: r, r2
  real(8) :: trho, tvpv, tvph, tvsv, tvsh, teta
  real(8) :: tmu, tecA, tecC, tecL, tecN, tAkappa, tCkappa
  complex(8) :: xrho, xecA, xecC, xecF, xecL, xecN
  complex(8) :: xmu, xFdC, xAmF2dC, xAkappa, xCkappa
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.

  anum(:,:,:) = dcmplx(0.d0)
  bnum(:,:,:) = dcmplx(0.d0)

  do i = 1, 10
    ! compute radius
    if (i <= 5) then
      r = gridRadiiForSource(1) + dble(i-1) / 4.d0 * (r0 - gridRadiiForSource(1))
    else
      r = gridRadiiForSource(2) - dble(10-i) / 4.d0 * (gridRadiiForSource(2) - r0)
    endif

    ! compute variable values at this radius
    call valueAtRadius(rhoPolynomials(:), r, rmax, trho)
    call valueAtRadius(vpvPolynomials(:), r, rmax, tvpv)
    call valueAtRadius(vphPolynomials(:), r, rmax, tvph)
    call valueAtRadius(vsvPolynomials(:), r, rmax, tvsv)
    call valueAtRadius(vshPolynomials(:), r, rmax, tvsh)
    call valueAtRadius(etaPolynomials(:), r, rmax, teta)
    tecA = trho * tvph * tvph
    tecC = trho * tvpv * tvpv
    tecL = trho * tvsv * tvsv
    tecN = trho * tvsh * tvsh
    tmu  = trho * tvsv * tvsv

    tAkappa = tecA - tmu * 4.d0 / 3.d0
    tCkappa = tecC - tmu * 4.d0 / 3.d0
    xAkappa = dcmplx(tAkappa) * coefQkappa
    xCkappa = dcmplx(tCkappa) * coefQkappa
    xrho = dcmplx(trho)
    xmu  = dcmplx(tmu) * coefQmu
    xecL = dcmplx(tecL) * coefQmu
    xecN = dcmplx(tecN) * coefQmu
    xecA = xAkappa + xmu * dcmplx(4.d0 / 3.d0)
    xecC = xCkappa + xmu * dcmplx(4.d0 / 3.d0)
    xecF = dcmplx(teta) * (xecA - dcmplx(2.d0) * xecL)
    xFdC = xecF / xecC
    xAmF2dC = xecA - xecF * xFdC

    r2 = r * r
    omegaDamped2 = dcmplx(omega, -omegaI) ** 2

    ! compute matrix elements
    anum(1,1,i) = - dcmplx(2.d0 / r) * xFdC
    anum(1,2,i) =   dcmplx(1.d0) / xecC
    anum(2,1,i) = - xrho * omegaDamped2 + dcmplx(4.d0 / r2) * (xAmF2dC - xecN)
    anum(2,2,i) =   dcmplx(2.d0 / r) * (xFdC - 1.d0)
    anum(3,1,i) = - dcmplx(1.d0 / r)
    anum(3,3,i) =   dcmplx(1.d0 / r)
    anum(3,4,i) =   dcmplx(1.d0) / xecL
    anum(4,1,i) = - dcmplx(2.d0 / r2) * (xAmF2dC - xecN)
    anum(4,2,i) = - dcmplx(1.d0 / r) * xFdC
    anum(4,3,i) = - xrho * omegaDamped2 - dcmplx(2.d0 / r2) * xecN
    anum(4,4,i) = - dcmplx(3.d0 / r)
    bnum(1,3,i) =   dcmplx(1.d0 / r) * xFdC
    bnum(2,3,i) = - dcmplx(2.d0 / r2) * (xAmF2dC - xecN)
    bnum(2,4,i) =   dcmplx(1.d0 / r)
    bnum(4,3,i) =   xAmF2dC / dcmplx(r2)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
 ! Computing the excitation vector using Geller and Hatori (1995).
 ! Parameters by N.Takeuchi 1995.7
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeya(aa, bb, largeL2, gridRadiiForSource, r0, ya, yb, yc, yd)
!----------------------------------------------------------------------------------------------------------------------------
  common a,b,cl2,itmp
  external eqmotion1
  real(8), parameter :: pi = 3.1415926535897932d0
  complex(8) :: a(64), b(64)
  real(8) :: cl2
  integer :: itmp

  real(8), intent(in) :: r0, largeL2, gridRadiiForSource(2)
  complex(8), intent(in) :: aa(160), bb(160)
  complex(8), intent(out) :: ya(4), yb(4), yc(4), yd(4)
  complex(8) :: work(4, 2), yn(4)
  real(8) :: xs, xe, dr
  integer :: i, k
  integer :: ktmp1, ktmp2, ktmp3, ktmp4

  ! -----------------------<< common variables >>-----------------------
  cl2 = largeL2
  ! ---------------------<< numerical integration >>---------------------
  ! integration from the lower boundary
  ya(1) = dcmplx(0.d0)
  ya(2) = dcmplx(1.d0)
  ya(3) = dcmplx(0.d0)
  ya(4) = dcmplx(0.d0)
  yb(1) = dcmplx(0.d0)
  yb(2) = dcmplx(0.d0)
  yb(3) = dcmplx(0.d0)
  yb(4) = dcmplx(1.d0)
  dr = (r0 - gridRadiiForSource(1)) / 2.d0
  if (dr >= 0.d0) then
    xs = gridRadiiForSource(1)
    do k = 1, 2
      ktmp1 = 32 * (k - 1)
      ktmp2 = ktmp1 - 16
      do i = 1, 64
        if (i <= 32) then
          a(i) = aa(i + ktmp1)
          b(i) = bb(i + ktmp1)
        else
          a(i) = aa(i + ktmp2)
          b(i) = bb(i + ktmp2)
        end if
      end do
      itmp = 0
      xe = xs + dr
      call rk3(4,eqmotion1,xs,xe,1,ya,yn,4,work)
      itmp = 0
      xs = xe - dr
      call rk3(4,eqmotion1,xs,xe,1,yb,yn,4,work)
    end do
  end if

  ! integration from the upper boundary
  yc(1) = dcmplx(0.d0)
  yc(2) = dcmplx(1.d0)
  yc(3) = dcmplx(0.d0)
  yc(4) = dcmplx(0.d0)
  yd(1) = dcmplx(0.d0)
  yd(2) = dcmplx(0.d0)
  yd(3) = dcmplx(0.d0)
  yd(4) = dcmplx(1.d0)
  dr = (gridRadiiForSource(2) - r0) / 2.d0
  if (dr >= 0.d0) then
    xs = gridRadiiForSource(2)
    do k = 1, 2
      ktmp1 = 144 - 32 * (k - 1)
      ktmp2 = ktmp1 - 32
      ktmp3 = ktmp1 - 48
      ktmp4 = ktmp1 - 80
      do i = 1, 64
        if (i <= 32) then
          if (i <= 16) then
            a(i) = aa(i + ktmp1)
            b(i) = bb(i + ktmp1)
          else
            a(i) = aa(i + ktmp2)
            b(i) = bb(i + ktmp2)
          end if
        else
          if (i <= 48) then
            a(i) = aa(i + ktmp3)
            b(i) = bb(i + ktmp3)
          else
            a(i) = aa(i + ktmp4)
            b(i) = bb(i + ktmp4)
          end if
        end if
      end do
      itmp = 0
      xe = xs - dr
      call rk3(4,eqmotion1,xs,xe,1,yc,yn,4,work)
      itmp = 0
      xs = xe + dr
      call rk3(4,eqmotion1,xs,xe,1,yd,yn,4,work)
    end do
  end if

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the excitation vector using Geller and Hatori(1995).
! Parameters by N.Takeuchi 1995.7
!----------------------------------------------------------------------------------------------------------------------------
subroutine calg(l, m, coef1, coef2, lsq, ecC0, ecF0, ecL0, ya, yb, yc, yd, ra, r0, mt, g)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: l, m
  real(8), intent(in) :: r0, ecC0, ecF0, ecL0, ra(2), lsq
  complex(8), intent(in) :: coef1, coef2
  complex(8), intent(in) :: ya(4), yb(4), yc(4), yd(4)
  real(8), intent(in) :: mt(3,3)
  complex(8), intent(out) :: g(4)
  complex(8) :: dd, ee, s1, s2, s(4)
  integer :: ip(4), ier, i
  complex(8) :: a(4,4), b(4), wk(4), xtmp
  real(8) :: eps, sgn, b1, b2, r03, dtmp(4)

  ! Initialization
  eps = -1.d0

  ! ---------------------<< parameter computation >>---------------------
  ! computation of the discontinuity
  if ( m.ge.0 ) then
    sgn = 1.d0
  else
    sgn = -1.d0
  end if
  b1 = sqrt( dble(2*l+1)/(16.d0*pi) )
  if ( l.ne.0 ) then
    b2 = sqrt( dble(2*l+1)*dble(l-1)*dble(l+2)/(64.d0*pi) )
  end if
  if ( iabs(m).eq.2 ) then
    dd = dcmplx( 0.d0 )
    ee = dcmplx( 0.d0 )
  end if
  if ( iabs(m).eq.1 ) then
    dd = dcmplx( 0.d0 )
    ee = dcmplx( - b1 * sgn * mt(1,2), b1 * mt(1,3) ) / ( dcmplx( r0 * r0 * ecL0 ) * coef1 )
  end if
  if ( iabs(m).eq.0 ) then
    dd = dcmplx( 2.d0 * b1 * mt(1,1) / ( r0 * r0 ) ) / (   dcmplx( ecC0 - 4.d0 / 3.d0 * ecL0 ) * coef2 &
      + dcmplx( 4.d0/3.d0*ecL0 ) * coef1 )
    ee = dcmplx( 0.d0 )
  end if

  if ( iabs(m).eq.0 ) then
    r03 = r0 * r0 * r0
    xtmp =   ( ( ecF0 + 2.d0 / 3.d0 * ecL0) * coef2 &
      - 2.d0 / 3.d0 * ecL0 * coef1 ) / ( (ecC0 - 4.d0 / 3.d0 * ecL0)* coef2 &
      + 4.d0 / 3.d0 * ecL0 * coef1 )
    s1 = dcmplx( - 2.d0 * b1 * ( mt(2,2) + mt(3,3) ) / r03 ) + dcmplx( 4.d0 * b1 * mt(1,1) / r03 ) * xtmp
    s2 = dcmplx( b1 * lsq * ( mt(2,2) + mt(3,3) ) / r03 ) - dcmplx( 2.d0 * b1 * lsq * mt(1,1) / r03 ) * xtmp
  end if
  if ( iabs(m).eq.1 ) then
    s1 = dcmplx( 0.d0 )
    s2 = dcmplx( 0.d0 )
  end if
  if ( iabs(m).eq.2 ) then
    r03 = r0 * r0 * r0
    s1 = dcmplx( 0.d0 )
    s2 = dcmplx( - b2 * ( mt(2,2) - mt(3,3) ) / r03 , sgn * 2.d0 * b2 * mt(2,3) / r03 )
  end if

  s(1) = dd
  s(2) = s1
  if ( l.ne.0 ) then
    s(3) = dcmplx( dble(ee)/lsq, dimag(ee)/lsq )
    s(4) = dcmplx( dble(s2)/lsq, dimag(s2)/lsq )
  else
    s(3) = dcmplx( 0.d0 )
    s(4) = dcmplx( 0.d0 )
  end if
  ! consideration of the boundary conditions
  ! determination of the analytical solution
  if ( l.ne.0 ) then
    call sab1( ya,yb,yc,yd,s,a,b )
    call glu(a,4,4,b,eps,wk,ip,ier)
  else
    call sab2( ya,yc,s,a,b )
    call glu(a,2,4,b,eps,wk,ip,ier)
    b(3) = b(2)
    b(2) = dcmplx( 0.d0 )
    b(4) = dcmplx( 0.d0 )
  end if
  ! computation of the excitation vector
  dtmp(1) = - ra(1) * ra(1)
  dtmp(2) = dtmp(1) * lsq
  dtmp(3) = ra(2) * ra(2)
  dtmp(4) = dtmp(3) * lsq
  do i=1,4
    xtmp = b(i)
    g(i) = dcmplx( dble(xtmp)*dtmp(i), dimag(xtmp)*dtmp(i) )
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Equation of the motion of the solid medium.
!----------------------------------------------------------------------------------------------------------------------------
subroutine eqmotion1(r, y, f)
!----------------------------------------------------------------------------------------------------------------------------
  common a, b, l2, itmp
  complex(8) :: a(64), b(64)
  real(8) :: l2
  integer :: itmp

  real(8), intent(in) :: r
  complex(8), intent(in) :: y(4)
  complex(8), intent(out) :: f(4)
  complex(8) :: c(4,4)
  integer :: i, j, mtmp

  itmp = itmp + 1

  ! Computation of the differential coefficients
  mtmp = 16 * (itmp - 1)
  do j = 1, 4
    do i = 1, 4
      mtmp = mtmp + 1
      c(i, j) = a(mtmp) + l2 * b(mtmp)
    end do
  end do

  ! Compute function values
  do i = 1, 4
    f(i) = dcmplx(0.d0)
    do j = 1, 4
      f(i) = f(i) + c(i, j) * y(j)
    end do
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Determination of the matrix imposing the boundary conditions.
!----------------------------------------------------------------------------------------------------------------------------
subroutine sab1(ya, yb, yc, yd, s, a, b)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  complex(8), intent(in) :: ya(4), yb(4), yc(4), yd(4), s(4)
  complex(8), intent(out) :: a(4,4), b(4)
  integer :: i

  do i = 1, 4
    a(i, 1) = -ya(i)
    a(i, 2) = -yb(i)
    a(i, 3) =  yc(i)
    a(i, 4) =  yd(i)
    b(i) = s(i)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Determination of the matrix imposing the boundary conditions.
! Only the first two elements of the arrays are used.
!----------------------------------------------------------------------------------------------------------------------------
subroutine sab2(ya, yc, s, a, b)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  complex(8), intent(in) :: ya(4), yc(4), s(4)
  complex(8), intent(out) :: a(4,4), b(4)
  integer :: i

  do i = 1, 2
    a(i, 1) = -ya(i)
    a(i, 2) =  yc(i)
    b(i) = s(i)
  end do

end subroutine
