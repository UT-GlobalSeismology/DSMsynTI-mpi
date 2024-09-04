
!----------------------------------------------------------------------------------------------------------------------------
! Add transpose of tridiagonal matrix to itself.
! This is to replace "I_(k'k)^4" with "I_(k'k)^4 + I_(kk')^4", because "I_(k'k)^4" appears only in this form.
!  (See eq. 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine addTranspose(nLayerInZoneI, matIn, matSum)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: matIn(4*nLayerInZoneI)  ! Input tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg].
  real(8), intent(out) :: matSum(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg].
  integer :: iLayer

  ! Replace "I_(k'k)^4" with "I_(k'k)^4 + I_(kk')^4". (See eq. 19 of Kawai et al. 2006.)
  do iLayer = 1, 4 * nLayerInZoneI
    select case(mod(iLayer, 4))
     case (0, 1)  ! (k1, k2) = (i, i), (i+1, i+1) -> double the value
      matSum(iLayer) = 2.d0 * matIn(iLayer)
     case (2)  ! (k1, k2) = (i, i+1) -> add "(k1, k2) = (i+1, i)" case
      matSum(iLayer) = matIn(iLayer) + matIn(iLayer + 1)
     case (3)  ! (k1, k2) = (i+1, i) -> add "(k1, k2) = (i, i+1)" case
      matSum(iLayer) = matIn(iLayer - 1) + matIn(iLayer)
    end select
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes (omega^2 T - (I2 - I4 - I4' + I6 - 2*I7)). (See eqs. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA0(nLayerInZoneI, omega, omegaI, t, h1, h2sum, h3, h4, coefQmu, a0Tmp)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: t(4*nLayerInZoneI)  ! T matrix stored for each (iLayer, k', k)-pair [10^12 kg].
  real(8), intent(in) :: h1(4*nLayerInZoneI), h2sum(4*nLayerInZoneI), h3(4*nLayerInZoneI), h4(4*nLayerInZoneI)
  !::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2]. Note that h2sum is (I4 + I4').
  complex(8), intent(in) :: coefQmu  ! Coefficient to multiply to elastic moduli for attenuation.
  complex(8), intent(out) :: a0Tmp(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  real(8) :: h
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  do i = 1, 4 * nLayerInZoneI
    ! I2 - I4 - I4' + I6 - 2*I7. (See eq. 19 of Kawai et al. 2006.)
    h = h1(i) - h2sum(i) + h3(i) - 2.d0 * h4(i)
    ! omega^2 T - (I2 - I4 - I4' + I6 - 2*I7). (See eq. 2 of Kawai et al. 2006.)
    a0Tmp(i) = omegaDamped2 * dcmplx(t(i)) - coefQmu * dcmplx(h)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes (- (I7)). (See eqs. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA2(nLayerInZoneI, h4, coefQmu, a2Tmp)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: h4(4*nLayerInZoneI)  ! Part of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(in) :: coefQmu  ! Coefficient to multiply to elastic moduli for attenuation.
  complex(8), intent(out) :: a2Tmp(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  integer :: i

  do i = 1, 4 * nLayerInZoneI
    ! -(I7). (See eqs. 2 & 19 of Kawai et al. 2006.)
    a2Tmp(i) = - coefQmu * dcmplx(h4(i))
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Overlapping the coefficient matrix elements for a certain zone in the solid part.
! The results are the diagonal and subdiagonal components of the tridiagonal matrix,
!  stored for each (k', k) = (1,1), (1,2),(2,2), (2,3),(3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine overlapA(nLayerInZoneI, aTmp, aOut)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  complex(8), intent(in) :: aTmp(4*nLayerInZoneI)  ! Tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(inout) :: aOut(2, nLayerInZoneI+1)  ! Diagonal and subdiagonal components of the overlapped matrix
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2]. Should be initialized with 0s beforehand.
  integer :: i

  do i = 1, nLayerInZoneI
    ! (i,i)-component
    if (i == 1) then
      ! This overlaps with previous zone (if present).
      aOut(2, i) = aOut(2, i) + aTmp(1)
    else
      aOut(2, i) = aTmp(4 * i - 4) + aTmp(4 * i - 3)
    end if
    ! (i,i+1)-component
    aOut(1, i + 1) = aTmp(4 * i - 2)
  end do
  ! (N,N)-component
  aOut(2, nLayerInZoneI + 1) = aTmp(4 * nLayerInZoneI)

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Assembling the coefficient matrix 'A' in the solid part from several parts.
!----------------------------------------------------------------------------------------------------------------------------
subroutine assembleA(nGrid, largeL2, a0, a2, a)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nGrid  ! Total number of grid points.
  real(8), intent(in) :: largeL2  ! L^2 = l(l+1).
  complex(8), intent(in) :: a0(2,nGrid), a2(2,nGrid)  ! Parts of the A matrix, in diagonal and subdiagonal component format
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  complex(8), intent(out) :: a(2,nGrid)  ! Assembled A matrix [10^12 kg/s^2].
  integer :: iGrid
  complex(8) :: largeL2c

  ! This is changed to complex beforehand to reduce computation amount inside loop.
  largeL2c = dcmplx(largeL2)

  do iGrid = 1, nGrid
    ! A = A0 + L^2 A2.
    a(1, iGrid) = a0(1, iGrid) + largeL2c * a2(1, iGrid)
    a(2, iGrid) = a0(2, iGrid) + largeL2c * a2(2, iGrid)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the coefficient matrix 'A' for a certain zone in the solid part. (See eqs. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeA(nLayerInZoneI, omega, omegaI, largeL2, t, h1, h2sum, h3, h4, coefQmu, a)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: largeL2  ! L^2 = l(l+1).
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: t(4*nLayerInZoneI)  ! T matrix stored for each (iLayer, k', k)-pair [10^12 kg].
  real(8), intent(in) :: h1(4*nLayerInZoneI), h2sum(4*nLayerInZoneI), h3(4*nLayerInZoneI), h4(4*nLayerInZoneI)
  !::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2]. Note that h2sum is (I4 + I4').
  complex(8), intent(in) :: coefQmu  ! Coefficient to multiply to elastic moduli for attenuation.
  complex(8), intent(out) :: a(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  real(8) :: h, largeL2m2
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  ! L^2 - 2
  ! This is computed beforehand to reduce computation amount inside loop.
  largeL2m2 = largeL2 - 2.d0

  do i = 1, 4 * nLayerInZoneI
    ! I2 - I4 - I4' + I6 + (L^2 - 2)*I7. (See eq. 19 of Kawai et al. 2006.)
    h = h1(i) - h2sum(i) + h3(i) + largeL2m2 * h4(i)
    ! omega^2 T - H. (See eq. 2 of Kawai et al. 2006.)
    a(i) = omegaDamped2 * dcmplx(t(i)) - coefQmu * dcmplx(h)
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the excitation vector g.
! Here, the source time function is set as the delta function (in the time domain).
! (Note that the Fourier transform of the delta function is a constant function.
!  The g set in this subroutine does not depend on omega, meaning it is a constant function in the frequency domain.)
! The unit of g is [10^15 N] in the frequency domain, which corresponds to [10^15 N/s] in the time domain.
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeG(l, m, iLayerOfSource, r0, mt, ecL0, coefQmu, aaParts, aSourceParts, aSource, gdr, g)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: l  ! Angular order.
  integer, intent(in) :: m  ! Azimuthal order.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  complex(8), intent(in) :: coefQmu  ! Coefficient to multiply to elastic moduli for attenuation.
  complex(8), intent(in) :: aaParts(4), aSourceParts(8)  ! Unassembled A matrix [10^12 kg/s^2].
  complex(8), intent(inout) :: aSource(2,3)  ! Assembled A matrix [10^12 kg/s^2].
  complex(8), intent(inout) :: gdr(3)  ! Working array.
  complex(8), intent(out) :: g(*)  ! The vector -g [10^15 N].
  real(8) :: b, sgnM, eps
  complex(8) :: dd, gS_or_cS(3)
  integer :: i, ier
  complex(8) :: z(3)

  ! Initialize.
  gS_or_cS(:) = dcmplx(0.d0, 0.d0)
  dd = dcmplx(0.d0, 0.d0)
  eps = -1.d0

  ! Record sign of m.
  if (m >= 0) then
    sgnM = 1.d0
  else
    sgnM = -1.d0
  end if

  if (abs(m) == 1) then
    ! b1 in eq. (26) of Kawai et al. (2006).
    b = sqrt(dble(2 * l + 1) / (16.d0 * pi))
    ! D3 [km] of eq. (26) of Kawai et al. (2006).
    dd = dcmplx(b) * dcmplx(sgnM * mt(1, 3), mt(1, 2)) / (dcmplx(r0 * r0 * ecL0) * coefQmu)

    !TODO ??
    do i = 2, 3
      gS_or_cS(i) = -dd * (aSourceParts(i * 2 + 1) + aSourceParts(i * 2 + 2))  ! i=2 -> 5, 6 ; i=3 -> 7, 8
    end do

  else if (abs(m) == 2) then
    ! b2 in eq. (27) of Kawai et al. (2006).
    ! NOTE that integers are casted with dble() before multiplying, because the product can exceed the size of integer(4).
    b = sqrt(dble(2 * l + 1) * dble(l - 1) * dble(l + 2) / (64.d0 * pi))
    ! -1 * gk3 [10^15 N] of eq. (27) of Kawai et al. (2006).
    gS_or_cS(2) = dcmplx(b / r0) * dcmplx(2.d0 * mt(2, 3), sgnM * (mt(2, 2) - mt(3, 3)))
  end if

  ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
  ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
  if ((m == -2) .or. (m == -l)) then
    call decomposeAByCholesky(aSource(:,:), 3, 1, 2,eps, gdr, ier)
  end if
  ! Solve Ac=g (i.e. (omega^2 T - H) c = -g) for grids near source.
  call solveWholeCAfterCholesky(aSource(:,:), 3, 1, 2, gS_or_cS(:), gdr, z)

  ! Add displacement to c.
  gS_or_cS(3) = gS_or_cS(3) + dd

  ! Compute excitation vector g.
  g(iLayerOfSource) = aaParts(1) * gS_or_cS(1) + aaParts(2) * gS_or_cS(3)
  g(iLayerOfSource + 1) = aaParts(3) * gS_or_cS(1) + aaParts(4) * gS_or_cS(3)

end subroutine
