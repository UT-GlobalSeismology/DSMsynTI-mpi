
!------------------------------------------------------------------------
! Computing \int dr con r^rpow X_k1^(dot1) X_k2^(dot2). (See eq. 16 of Kawai et al. 2006.)
! "I_(k'k)^4" is replaced with "I_(k'k)^4 + I_(kk')^4". (See eq. 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeIntermediateIntegral(nLayerInZoneI, valuedRadiiInZoneI, conInZoneI, rpow, dot1, dot2, mat, work)
!------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxrpow = 2  ! Maximum value of rpow to allow.

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: conInZoneI(nLayerInZoneI+1)  ! Values of a variable at each point (with 2 values at boundaries).
  integer, intent(in) :: rpow  ! The exponent of r.
  integer, intent(in) :: dot1, dot2  ! Whether or not to differentiate X_k1 and X_k2 (1: differentiate, 0: do not differentiate).
  real(8), intent(out) :: mat(4*nLayerInZoneI)  ! Resulting integrals, "I_(k'k)^4" replaced with "I_(k'k)^4 + I_(kk')^4"
  !::::::::::::::::::::::::::::::::::::::::::::::: (See eq. 19 of Kawai et al. 2006.)
  real(8), intent(out) :: work(4*nLayerInZoneI)  ! Resulting integrals. I^0 is in [10^12 kg], others are in [10^12 kg/s^2].
  integer :: iLayer, j1, j2, i, iRow
  real(8) :: a(2,2), b(2,2), c(5), rh

  ! Check input validity.
  if (rpow > maxrpow) stop "Invalid arguments. (computeIntermediateIntegral)"

  ! Compute matrix elements.
  do iLayer = 1, nLayerInZoneI
    ! layer thickness [km]
    rh = valuedRadiiInZoneI(iLayer + 1) - valuedRadiiInZoneI(iLayer)

    ! Set X_k1^(dot1), for both k1=i and k1=i+1.
    select case(dot1)
     case (0)
      a(:, 1) = [valuedRadiiInZoneI(iLayer + 1) / rh, -1.d0 / rh]
      a(:, 2) = [-valuedRadiiInZoneI(iLayer) / rh, 1.d0 / rh]
     case (1)
      a(:, 1) = [-1.d0 / rh, 0.d0]
      a(:, 2) = [1.d0 / rh, 0.d0]
     case default
      stop "Invalid arguments. (computeIntermediateIntegral)"
    end select

    ! Set X_k2^(dot2), for both k2=i and k2=i+1.
    select case(dot2)
     case (0)
      b(:, 1) = [valuedRadiiInZoneI(iLayer + 1) / rh, -1.d0 / rh]
      b(:, 2) = [-valuedRadiiInZoneI(iLayer) / rh, 1.d0 / rh]
     case (1)
      b(:, 1) = [-1.d0 / rh, 0.d0]
      b(:, 2) = [1.d0 / rh, 0.d0]
     case default
      stop "Invalid arguments. (computeIntermediateIntegral)"
    end select

    do j1 = 1, 2  ! k1=i and k1=i+1
      do j2 = 1, 2  ! k2=i and k2=i+1
        c(:) = 0.d0
        ! Multiply X_k1^(dot1) and X_k2^(dot2).
        call multiplyPolynomials(2, a(:, j1), 2, b(:, j2), 3, c(:))
        ! Multiply by r^rpow.
        if (rpow > 0) then
          do i = 3, 1, -1
            c(i + rpow) = c(i)
            c(i) = 0.d0
          end do
        end if
        ! Integrate; the result is saved for each (iLayer, k1, k2)-pair.
        iRow = 4 * (iLayer - 1) + 2 * (j1 - 1) + j2
        call integrateProduct(5, c, valuedRadiiInZoneI(iLayer), valuedRadiiInZoneI(iLayer + 1), &
          valuedRadiiInZoneI(iLayer:), conInZoneI(iLayer:), work(iRow))
      end do
    end do
  end do

  ! Replace "I_(k'k)^4" with "I_(k'k)^4 + I_(kk')^4". (See eq. 19 of Kawai et al. 2006.)
  if (dot1 /= dot2) then
    do iLayer = 1, 4 * nLayerInZoneI
      select case(mod(iLayer, 4))
       case (0, 1)  ! (k1, k2) = (i, i), (i+1, i+1) -> double the value
        mat(iLayer) = 2.d0 * work(iLayer)
       case (2)  ! (k1, k2) = (i, i+1) -> add "(k1, k2) = (i+1, i)" case
        mat(iLayer) = work(iLayer) + work(iLayer + 1)
       case (3)  ! (k1, k2) = (i+1, i) -> add "(k1, k2) = (i, i+1)" case
        mat(iLayer) = work(iLayer - 1) + work(iLayer)
      end select
    end do
  else
    mat(:) = work(:)
  end if

end subroutine


!------------------------------------------------------------------------
! Computing the (nc-1)-degree polynomial c(x) which is the product of
! the (na-1)-degree polynomial a(x) and the (nb-1)-degree polynomial b(x).
!------------------------------------------------------------------------
subroutine multiplyPolynomials(na, a, nb, b, nc, c)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: na, nb, nc  ! Size of the arrays of a, b, and c.
  real(8), intent(in) :: a(na), b(nb)  ! Coefficients of polynimials a and b in ascending order (a(x) = a1 + a2 x + a3 x^2 + ...).
  real(8), intent(out) :: c(nc)  ! Coefficients of polynimial c in ascending order.

  integer :: i, j

  ! Check input validity.
  if (na + nb - 1 /= nc) stop "Invalid arguments. (multiplyPolynomials)"

  ! Initialize the polynomial c.
  c(:) = 0.d0

  ! Compute the product polynomial.
  do i = 1, na
    do j = 1, nb
      c(i + j - 1) = c(i + j - 1) + a(i) * b(j)
    end do
  end do

end subroutine


!------------------------------------------------------------------------
! Evaluating the integrated value of p(r)*con(r) from 'lowerRadius' to 'upperRadius'.
! Here, p(r) is an (n-1)-degree polynomial, and con(r) is the profile of a variable.
! The range [lowerRadius, upperRadius] must be within a certain layer [valuedRadii(1), valuedRadii(2)].
!------------------------------------------------------------------------
subroutine integrateProduct(n, p, lowerRadius, upperRadius, valuedRadii, con, result)
!------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxn = 5  ! Maximum number of polynomial degrees.

  integer, intent(in) :: n  ! Size of the array of p.
  real(8), intent(in) :: p(n)  ! Coefficients of the polynimial in ascending order (p(r) = p1 + p2 r + p3 r^2 + ...).
  real(8), intent(in) :: lowerRadius, upperRadius  ! Radius range to integrate [km].
  real(8), intent(in) :: valuedRadii(2)  ! Radii at both ends of an interval containing integration range [km].
  real(8), intent(in) :: con(2)  ! Values of a variable at both ends of an interval containing integration range.
  real(8), intent(out) :: result
  real(8) :: q(2), pq(maxn+1)

  ! Check input validity.
  if (n > maxn) stop 'Degree of polynomial is too large. (integrateProduct)'

  ! Express con(r) as a polynomial (linear) function q1+q2*r.
  q(2) = (con(2) - con(1)) / (valuedRadii(2) - valuedRadii(1))  ! slope
  q(1) = con(1) - q(2) * valuedRadii(1)  ! intercept
  ! Compute p(r)*con(r).
  call multiplyPolynomials(n, p(:), 2, q(:), n + 1, pq(:))
  ! Evaluate integrated value within subrange [lowerRadius, upperRadius].
  call integratePolynomial(n + 1, pq(:), lowerRadius, upperRadius, result)

end subroutine


!------------------------------------------------------------------------
! Evaluating the integrated value of an (n-1)-degree polynomial 'p(x)' from 'x1' to 'x2'.
!------------------------------------------------------------------------
subroutine integratePolynomial(n, p, x1, x2, result)
!------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxn = 6  ! Maximum number of polynomial degrees

  integer, intent(in) :: n  ! Size of the array of p.
  real(8), intent(in) :: p(n)  ! Coefficients of the polynimial in ascending order (p(r) = p1 + p2 r + p3 r^2 + ...).
  real(8), intent(in) :: x1, x2  ! X range to integrate.
  real(8), intent(out) :: result
  integer :: i, j
  real(8) :: a(maxn), b(maxn), dx, xx

  ! Check input validity.
  if (n > maxn) stop 'Degree of polynomial is too large. (integratePolynomial)'

  ! Initialize: a=(1 x1 x1^2 ...), b=(1 x2 x2^2 ...).
  a(1) = 1.d0
  b(1) = 1.d0
  if (n >= 2) then
    do i = 2, n
      a(i) = a(i - 1) * x1
      b(i) = b(i - 1) * x2
    end do
  end if
  dx = x2 - x1

  ! Evaluate the integrated value.
  result = 0.d0
  do i = 1, n  ! Loop for each term of p(x).
    xx = 0.d0
    ! (x2^(i-1) + x1 x2^(i-2) + x1^2 x2^(i-3) + ... + x1^(i-1)) / i
    do j = 1, i
      xx = xx + a(j) * b(i - j + 1) / dble(i)
    end do
    ! i=1 : p1(x2-x1) ; i=2: p2(x2-x1)(x2+x1)/2 ; i=3: p3(x2-x1)(x2^2+x1x2+x1^2)/3 ; ...
    result = result + p(i) * dx * xx
  end do

end subroutine


!------------------------------------------------------------------------
! Computing the lumped mass matrix for a certain zone in the solid part. (See eqs. 15-17 of Cummins et al. 1994.)
!  T_kk^lumped = m_k r_k^2.
!  T_k'k^lumped = 0 when k' /= k.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeLumpedT(nLayerInZoneI, valuedRadiiInZoneI, rhoValuesInZoneI, tl)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: rhoValuesInZoneI(nLayerInZoneI+1)  ! Rho values at each point (with 2 values at boundaries) [g/cm^3].
  real(8), intent(out) :: tl(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg].
  integer :: i, nn
  real(8) :: c(3), lowerRadius, upperRadius

  ! Initialize. This is c(r) = r^2.
  c = [0.d0, 0.d0, 1.d0]

  do i = 1, nLayerInZoneI
    nn = 4 * (i - 1)

    ! Right side of m_k r_k^2 for k=i. Integrate rho*r^2.
    lowerRadius = valuedRadiiInZoneI(i)
    upperRadius = (valuedRadiiInZoneI(i) + valuedRadiiInZoneI(i + 1)) / 2.d0
    call integrateProduct(3, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), rhoValuesInZoneI(i:), tl(nn + 1))

    tl(nn + 2) = 0.d0
    tl(nn + 3) = 0.d0

    ! Left side of m_k r_k^2 for k=i+1. Integrate rho*r^2.
    lowerRadius = upperRadius
    upperRadius = valuedRadiiInZoneI(i + 1)
    call integrateProduct(3, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), rhoValuesInZoneI(i:), tl(nn + 4))
  end do

end subroutine


!------------------------------------------------------------------------
! Computing the lumped rigidity matrix for a certain zone in the solid part. (See eqs. 15-17 of Cummins et al. 1994.)
!  H_kk^lumped = s_k.
!  H_k'k^lumped = 0 when k' /= k.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeLumpedH(nLayerInZoneI, valuedRadiiInZoneI, muValuesInZoneI, hl)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: muValuesInZoneI(nLayerInZoneI+1)  ! Mu values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(out) :: hl(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  integer :: i, nn
  real(8) :: c(1), lowerRadius, upperRadius

  ! Initialize. This is c(r) = 1 (constant).
  c = [1.d0]

  do i = 1, nLayerInZoneI
    nn = 4 * (i - 1)

    ! Right side of s_k for k=i. Integrate mu (ecL or ecM).
    lowerRadius = valuedRadiiInZoneI(i)
    upperRadius = (valuedRadiiInZoneI(i) + valuedRadiiInZoneI(i + 1)) / 2.d0
    call integrateProduct(1, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), muValuesInZoneI(i:), hl(nn + 1))

    hl(nn + 2) = 0.d0
    hl(nn + 3) = 0.d0

    ! Left side of s_k for k=i+1. Integrate mu (ecL or ecM).
    lowerRadius = upperRadius
    upperRadius = valuedRadiiInZoneI(i + 1)
    call integrateProduct(1, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), muValuesInZoneI(i:), hl(nn + 4))
  end do

end subroutine


!------------------------------------------------------------------------
! Averaging the values of two tridiagonal matrices for a certain zone in the solid part. (See eq. 17 of Cummins et al. 1994.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeAverage(nLayerInZoneI, mat1, mat2, average)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: mat1(4*nLayerInZoneI), mat2(4*nLayerInZoneI)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::: Input tridiagonal matrices, stored for each (iLayer, k', k)-pair.
  real(8), intent(out) :: average(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  integer :: i

  do i = 1, 4 * nLayerInZoneI
    average(i) = (mat1(i) + mat2(i)) / 2.d0
  end do

end subroutine


!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes (omega^2 T - (I2 - I4 - I4' + I6 - 2*I7)). (See eq. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeA0(nLayerInZoneI, omega, omegaI, t, h1, h2, h3, h4, qCoef, a0)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: t(4*nLayerInZoneI)  ! T matrix stored for each (iLayer, k', k)-pair [10^12 kg].
  real(8), intent(in) :: h1(4*nLayerInZoneI), h2(4*nLayerInZoneI), h3(4*nLayerInZoneI), h4(4*nLayerInZoneI)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(in) :: qCoef  ! Coefficient to multiply to elastic moduli for attenuation.
  complex(8), intent(out) :: a0(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  real(8) :: h
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  do i = 1, 4 * nLayerInZoneI
    ! I2 - I4 - I4' + I6 - 2*I7. (See eq. 19 of Kawai et al. 2006.)
    h = h1(i) - h2(i) + h3(i) - 2.d0 * h4(i)
    ! omega^2 T - (I2 - I4 - I4' + I6 - 2*I7). (See eq. 2 of Kawai et al. 2006.)
    a0(i) = omegaDamped2 * dcmplx(t(i)) - qCoef * dcmplx(h)
  end do

end subroutine


!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes (- (I7)). (See eq. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeA2(nLayerInZoneI, h4, qCoef, a2)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: h4(4*nLayerInZoneI)  ! Part of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(in) :: qCoef  ! Coefficient to multiply to elastic moduli for attenuation.
  complex(8), intent(out) :: a2(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  integer :: i

  do i = 1, 4 * nLayerInZoneI
    ! -(I7). (See eq. 2 & 19 of Kawai et al. 2006.)
    a2(i) = - qCoef * dcmplx(h4(i))
  end do

end subroutine


!------------------------------------------------------------------------
! Assembling the coefficient matrix 'A' in the solid part from several parts.
!------------------------------------------------------------------------
subroutine assembleA(nGrid, largeL2, a0, a2, a)
!------------------------------------------------------------------------
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


!------------------------------------------------------------------------
! Computing the coefficient matrix 'A' for a certain zone in the solid part. (See eq. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeA(nLayerInZoneI, omega, omegaI, largeL2, t, h1, h2, h3, h4, qCoef, a)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: largeL2  ! L^2 = l(l+1).
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: t(4*nLayerInZoneI)  ! T matrix stored for each (iLayer, k', k)-pair [10^12 kg].
  real(8), intent(in) :: h1(4*nLayerInZoneI), h2(4*nLayerInZoneI), h3(4*nLayerInZoneI), h4(4*nLayerInZoneI)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(in) :: qCoef  ! Coefficient to multiply to elastic moduli for attenuation.
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
    h = h1(i) - h2(i) + h3(i) + largeL2m2 * h4(i)
    ! omega^2 T - H. (See eq. 2 of Kawai et al. 2006.)
    a(i) = omegaDamped2 * dcmplx(t(i)) - qCoef * dcmplx(h)
  end do

end subroutine


!------------------------------------------------------------------------
! Overlapping the coefficient matrix elements for a certain zone in the solid part.
! The results are the diagonal and subdiagonal components of the tridiagonal matrix,
!  stored for each (k', k) = (1,1), (1,2),(2,2), (2,3),(3,3), ...
!------------------------------------------------------------------------
subroutine overlapMatrixBlocks(nLayerInZoneI, aIn, aOut)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  complex(8), intent(in) :: aIn(4*nLayerInZoneI)  ! Tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  complex(8), intent(out) :: aOut(2, nLayerInZoneI+1)  ! Diagonal and subdiagonal components of the overlapped matrix
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2]. Should be initialized with 0s beforehand.
  integer :: j

  do j = 1, nLayerInZoneI
    ! (j,j)-component
    if (j == 1) then
      aOut(2, j) = aOut(2, j) + aIn(1)
    else
      aOut(2, j) = aOut(2, j) + aIn(4 * j - 4) + aIn(4 * j - 3)
    end if
    ! (j,j+1)-component
    aOut(1, j + 1) = aOut(1, j + 1) + aIn(4 * j - 2)
  end do
  ! (N,N)-component
  aOut(2, nLayerInZoneI + 1) = aOut(2, nLayerInZoneI + 1) + aIn(4 * nLayerInZoneI)

end subroutine


!------------------------------------------------------------------------
! Computing the excitation vector g.
!------------------------------------------------------------------------
subroutine computeG(l, m, iLayerOfSource, r0, mt, ecL0, qCoef, aaParts, aSourceParts, aSource, dr, g)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: l  ! Angular order.
  integer, intent(in) :: m  ! Azimuthal order.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  real(8), intent(in) :: r0, mt(3,3)  ! Depth [km] and moment tensor [10^25 dyn cm] of source.
  real(8), intent(in) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].
  complex(8), intent(in) :: qCoef  ! Coefficient to multiply to elastic moduli for attenuation.
  complex(8), intent(in) :: aaParts(4), aSourceParts(8)  ! Unassembled A matrix [10^12 kg/s^2].
  complex(8), intent(in) :: aSource(2,3)  ! Assembled A matrix [10^12 kg/s^2].
  complex(8), intent(inout) :: dr(3)  ! Working array.
  complex(8), intent(out) :: g(*)  ! The vector -g [10^15 N].
  real(8) :: b, sgnM
  complex(8) :: dd, gS_or_cS(3)
  integer :: i
  complex(8) :: z(3)
  real(8) :: eps, ier

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
    dd = dcmplx(b) * dcmplx(sgnM * mt(1, 3), mt(1, 2)) / (dcmplx(r0 * r0 * ecL0) * qCoef)

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

  ! Solve Ac=g (i.e. (omega^2 T - H) c = -g) for grids near source.
  if ((m == -2) .or. (m == -l)) then
    ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
    call solveWholeCFromStart(aSource(:,:), 3, 1, 2, gS_or_cS(:), eps, dr, z, ier)
  else
    ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
    call solveWholeCFromMiddle(aSource(:,:), 3, 1, 2, gS_or_cS(:), eps, dr, z, ier)
  end if

  ! Add displacement to c.
  gS_or_cS(3) = gS_or_cS(3) + dd

  ! Compute excitation vector g.
  g(iLayerOfSource) = aaParts(1) * gS_or_cS(1) + aaParts(2) * gS_or_cS(3)
  g(iLayerOfSource + 1) = aaParts(3) * gS_or_cS(1) + aaParts(4) * gS_or_cS(3)

end subroutine


