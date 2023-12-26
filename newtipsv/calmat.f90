
!------------------------------------------------------------------------!!Common
! Computing \int dr con r^rpow X_k1^(dot1) X_k2^(dot2). (See eq. 16 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeIntermediateIntegral(nLayerInZoneI, valuedRadiiInZoneI, valuesInZoneI, rpow, dot1, dot2, mat)
!------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxrpow = 2  ! Maximum value of rpow to allow.

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: valuesInZoneI(nLayerInZoneI+1)  ! Values of a variable at each point (with 2 values at boundaries).
  integer, intent(in) :: rpow  ! The exponent of r.
  integer, intent(in) :: dot1, dot2  ! Whether or not to differentiate X_k1 and X_k2 (1: differentiate, 0: do not differentiate).
  real(8), intent(out) :: mat(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: I^0 is in [10^12 kg], others are in [10^12 kg/s^2].
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
          valuedRadiiInZoneI(iLayer:), valuesInZoneI(iLayer:), mat(iRow))
      end do
    end do
  end do

end subroutine


!------------------------------------------------------------------------!!Common
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


!------------------------------------------------------------------------!!Common
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


!------------------------------------------------------------------------!!Common
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


!------------------------------------------------------------------------!!Common
! Computing the lumped mass matrix for a certain zone. (See eqs. 15-17 of Cummins et al. 1994.)
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


!------------------------------------------------------------------------!!Common
! Computing the lumped rigidity matrix for a certain zone. (See eqs. 15-17 of Cummins et al. 1994.)
!  H_kk^lumped = s_k.
!  H_k'k^lumped = 0 when k' /= k.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeLumpedH(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, hl)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: ecValuesInZoneI(nLayerInZoneI+1)  ! Modulus values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(out) :: hl(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  integer :: i, nn
  real(8) :: c(1), lowerRadius, upperRadius

  ! Initialize. This is c(r) = 1 (constant).
  c = [1.d0]

  do i = 1, nLayerInZoneI
    nn = 4 * (i - 1)

    ! Right side of s_k for k=i. Integrate elastic modulus.
    lowerRadius = valuedRadiiInZoneI(i)
    upperRadius = (valuedRadiiInZoneI(i) + valuedRadiiInZoneI(i + 1)) / 2.d0
    call integrateProduct(1, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), ecValuesInZoneI(i:), hl(nn + 1))

    hl(nn + 2) = 0.d0
    hl(nn + 3) = 0.d0

    ! Left side of s_k for k=i+1. Integrate elastic modulus.
    lowerRadius = upperRadius
    upperRadius = valuedRadiiInZoneI(i + 1)
    call integrateProduct(1, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), ecValuesInZoneI(i:), hl(nn + 4))
  end do

end subroutine


!------------------------------------------------------------------------!!Common
! Averaging the values of two tridiagonal matrices for a certain zone. (See eq. 17 of Cummins et al. 1994.)
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
! Computing the transpose of a tridiagonal matrix in a certain zone.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeTranspose(nLayerInZoneI, matIn, matOut)
!------------------------------------------------------------------------
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


!------------------------------------------------------------------------
! Subtracting two tridiagonal matrices for a certain zone.
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine subtractMatrix(nLayerInZoneI, mat1, mat2, difference)
!------------------------------------------------------------------------
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


!------------------------------------------------------------------------
! Computing the step-wise part of the unmodified operator for a certain zone.
! (See eqs. 3.44 and 3.45 of Geller & Takeuchi 1995; eqs. 11 and 12 of Takeuchi et al. 1996.)
! Note that in the above papers, the average of ec*r within each layer is taken,
!  but in this program, its value at the grid point is used.
!      / -3D0-D1  3D0+D1           0     ..    \
! 1/8 |  -D0-3D1  D0+3D1-3D1-D2  3D1+D2  0  ..  |
!      \    :            :           :   ..    /
!------------------------------------------------------------------------
subroutine computeStepH(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, h5)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: ecValuesInZoneI(nLayerInZoneI+1)  ! Modulus values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(out) :: h5(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair [10^12 kg/s^2].
  integer :: iLayer

  do iLayer = 1, nLayerInZoneI
    h5(4 * iLayer - 3) = - 3.0d0 / 8.0d0 * ecValuesInZoneI(iLayer) * valuedRadiiInZoneI(iLayer) &
      - 1.0d0 / 8.0d0 * ecValuesInZoneI(iLayer + 1) * valuedRadiiInZoneI(iLayer + 1)
    h5(4 * iLayer - 2) = - h5(4 * iLayer - 3)
    h5(4 * iLayer - 1) = - 1.0d0 / 8.0d0 * ecValuesInZoneI(iLayer) * valuedRadiiInZoneI(iLayer) &
      - 3.0d0 / 8.0d0 * ecValuesInZoneI(iLayer + 1) * valuedRadiiInZoneI(iLayer + 1)
    h5(4 * iLayer)   = - h5(4 * iLayer - 1)
  end do

end subroutine


!------------------------------------------------------------------------
! Computing the step-wise part of the modified matrix for a certain zone.
! (See eq. 16 of Takeuchi et al. 1996.)
! Note that in the paper, the average of ec*r within each layer is taken,
!  but in this program, its value at the grid point is used.
!       / -7D0  8D0  -D0   0     ..    \
!      |  -5D1  -3D1  9D1  -D1  0  ..  |
! 1/12 |    :     :    :    :    ..    |
!      |   ..  0  -5DN-1 -3DN-1 8DN-1  |
!       \           ..  0  -5DN  5DN  /
!------------------------------------------------------------------------
subroutine computeModifiedH1(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, hm1)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: ecValuesInZoneI(nLayerInZoneI+1)  ! Modulus values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(inout) :: hm1(-1:2, nLayerInZoneI+1)  ! Resulting matrix, stored by (offset from diagonal, row number)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  integer :: i

  i = 1
  hm1(0, i) = hm1(0, i) - 7.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm1(1, i) = hm1(1, i) + 8.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm1(2, i) = hm1(2, i) - 1.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)

  do i = 2, nLayerInZoneI - 1
    hm1(-1, i) = hm1(-1, i) - 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
    hm1(0, i) = hm1(0, i) - 3.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
    hm1(1, i) = hm1(1, i) + 9.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
    hm1(2, i) = hm1(2, i) - 1.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  end do

  i = nLayerInZoneI
  hm1(-1, i) = hm1(-1, i) - 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm1(0, i) = hm1(0, i) - 3.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm1(1, i) = hm1(1, i) + 8.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)

  i = nLayerInZoneI + 1
  hm1(-1, i) = hm1(-1, i) - 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm1(0, i) = hm1(0, i) + 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)

end subroutine


!------------------------------------------------------------------------
! Computing the step-wise part of the modified matrix for a certain zone.
! (See eq. 17 of Takeuchi et al. 1996.)
! Note that in the paper, the average of ec*r within each layer is taken,
!  but in this program, its value at the grid point is used.
!       / -5D0  5D0   0   ..           \
!      |  -8D1  3D1  5D1  0   ..        |
! 1/12 |    :    :    :   :   ..        |
!      |  .. 0  DN-1 -9DN-1 3DN-1 5DN-1 |
!       \    ..   0    DN   -8DN   7DN /
!------------------------------------------------------------------------
subroutine computeModifiedH2(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, hm2)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: ecValuesInZoneI(nLayerInZoneI+1)  ! Modulus values at each point (with 2 values at boundaries) [GPa].
  real(8), intent(inout) :: hm2(-2:1, nLayerInZoneI+1)  ! Resulting matrix, stored by (offset from diagonal, row number)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2].
  integer :: i

  i = 1
  hm2(0, i) = hm2(0, i) - 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm2(1, i) = hm2(1, i) + 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)

  i = 2
  hm2(-1, i) = hm2(-1, i) - 8.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm2(0, i) = hm2(0, i) + 3.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm2(1, i) = hm2(1, i) + 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)

  do i = 3, nLayerInZoneI
    hm2(-2, i) = hm2(-2, i) + 1.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
    hm2(-1, i) = hm2(-1, i) - 9.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
    hm2(0, i) = hm2(0, i) + 3.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
    hm2(1, i) = hm2(1, i) + 5.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  end do

  i = nLayerInZoneI + 1
  hm2(-2, i) = hm2(-2, i) + 1.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm2(-1, i) = hm2(-1, i) - 8.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)
  hm2(0, i) = hm2(0, i) + 7.0d0 / 12.0d0 * ecValuesInZoneI(i) * valuedRadiiInZoneI(i)

end subroutine


!------------------------------------------------------------------------
! Computing the transpose of a band matrix stored by (offset from diagonal, row number) in a certain zone.
!------------------------------------------------------------------------
subroutine computeTransposeMod(nLayerInZoneI, pMin, pMax, hModIn, hModOut)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: pMin, pMax
  real(8), intent(in) :: hModIn(pMin:pMax, nLayerInZoneI+1)  ! Input matrix, stored by (offset from diagonal, row number).
  real(8), intent(out) :: hModOut(-pMax:-pMin, nLayerInZoneI+1)  ! Resulting matrix, stored by (offset from diagonal, row number).
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



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------



