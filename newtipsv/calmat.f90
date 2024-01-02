
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
subroutine averageMatrix(nLayerInZoneI, mat1, mat2, average)
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
subroutine transposeMatrix(nLayerInZoneI, matIn, matOut)
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
! Computing the step-wise part of the unmodified operator for a certain zone.
! (See eqs. 3.44 and 3.45 of Geller & Takeuchi 1995; eqs. 11 and 12 of Takeuchi et al. 1996.)
! Note that in the above papers, the average of ec*r within each layer is taken,
!  but in this program, its value at the grid point is used.
!      / -3D0-D1  3D0+D1           0     ..    \
! 1/8 |  -D0-3D1  D0+3D1-3D1-D2  3D1+D2  0  ..  |
!      \    :            :           :   ..    /
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
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
    h5(4 * iLayer - 3) = - 3.d0 / 8.d0 * ecValuesInZoneI(iLayer) * valuedRadiiInZoneI(iLayer) &
      - 1.d0 / 8.d0 * ecValuesInZoneI(iLayer + 1) * valuedRadiiInZoneI(iLayer + 1)
    h5(4 * iLayer - 2) = - h5(4 * iLayer - 3)
    h5(4 * iLayer - 1) = - 1.d0 / 8.d0 * ecValuesInZoneI(iLayer) * valuedRadiiInZoneI(iLayer) &
      - 3.d0 / 8.d0 * ecValuesInZoneI(iLayer + 1) * valuedRadiiInZoneI(iLayer + 1)
    h5(4 * iLayer) = - h5(4 * iLayer - 1)
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
!------------------------------------------------------------------------
subroutine computeModifiedHR(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, hm1)
!------------------------------------------------------------------------
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
! The result is a band matrix, stored by (offset from diagonal, row number) for each zone.
!------------------------------------------------------------------------
subroutine computeModifiedHL(nLayerInZoneI, valuedRadiiInZoneI, ecValuesInZoneI, hm2)
!------------------------------------------------------------------------
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


!------------------------------------------------------------------------
! Computing the transpose of a band matrix stored by (offset from diagonal, row number) in a certain zone.
!------------------------------------------------------------------------
subroutine transposeMatrixMod(nLayerInZoneI, pMin, pMax, hModIn, hModOut)
!------------------------------------------------------------------------
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


!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes the part of H_(k'1k1) and H_(k'2k2) without largeL coefficients. (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a block tridiagonal matrix, stored for each (iLayer, k'-gamma', k-gamma)-pair.
! Only elements in the upper triangle of A is computed.
!------------------------------------------------------------------------
subroutine computeA0Solid(nLayerInZoneI, omega, omegaI, t, h1x, h2L, h2N, h3ay, h4aL, h4aN, h5ay, h6aL, h6aN, h7y, h7z, h8L, h8N, &
  coefQmu, coefQkappa, a0Tmp)
!------------------------------------------------------------------------
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



!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes H_(k'1k2) and H_(k'2k1) (largeL is not multiplied). (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a block tridiagonal matrix, stored for each (iLayer, k'-gamma', k-gamma)-pair.
! Only elements in the upper triangle of A is computed.
!------------------------------------------------------------------------
subroutine computeA1Solid(nLayerInZoneI, h1x, h2L, h2N, h3y, h4L, h4N, h5y, h6L, h6N, coefQmu, coefQkappa, a1Tmp)
!------------------------------------------------------------------------
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
        + coefQmu * dcmplx(2.d0 * h2N(iBlock)) &  ! 2N
        - coefQmu * dcmplx(h2L(iBlock)) &  ! -L
        + coefQmu * dcmplx(h4L(iBlock)) &  ! L
        - coefQkappa * dcmplx(h5y(iBlock)) + coefQmu * dcmplx(2.d0/3.d0 * h6N(iBlock))  ! -[(F+2N/3) - 2N/3] = -F
      a1Tmp(iElement) = -hh1
    end if

  end do

end subroutine


!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes the part of H_(k'1k1) and H_(k'2k2) with coefficient largeL^2. (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a block tridiagonal matrix, stored for each (iLayer, k'-gamma', k-gamma)-pair.
! Only elements in the upper triangle of A is computed.
!------------------------------------------------------------------------
subroutine computeA2Solid(nLayerInZoneI, h1x, h2L, h2N, coefQmu, coefQkappa, a2Tmp)
!------------------------------------------------------------------------
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


!------------------------------------------------------------------------
! Overlapping the coefficient matrix elements for a certain zone in the solid part.
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!------------------------------------------------------------------------
subroutine overlapASolid(nLayerInZoneI, aTmp, aOut)
!------------------------------------------------------------------------
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
    aOut(3, 2 * iGrid + 1) = aTmp(16 * iGrid - 10)
    aOut(1, 2 * iGrid + 2) = aTmp(16 * iGrid - 9)
    aOut(2, 2 * iGrid + 2) = aTmp(16 * iGrid - 8)
  end do
  ! (N,N)-block
  iGrid = nLayerInZoneI + 1
  aOut(4, 2 * iGrid - 1) = aTmp(16 * iGrid - 19)
  aOut(3, 2 * iGrid) = aTmp(16 * iGrid - 18)
  aOut(4, 2 * iGrid) = aTmp(16 * iGrid - 16)

end subroutine


!------------------------------------------------------------------------
! Adding the modified part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes H_(k'1k2) and H_(k'2k1) (largeL is not multiplied). (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!------------------------------------------------------------------------
subroutine addModifiedHToA1(nLayerInZoneI, coefQmu, coefQkappa, hModL3y, hModR4L, hModL4N, hModR5y, hModL6L, hModR6N, a1)
  !------------------------------------------------------------------------
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
        + coefQkappa * dcmplx(hModL3y(p, iGrid)) - coefQmu * dcmplx(2.d0 / 3.d0 * hModL4N(p, iGrid)) &  ! -[(F+2N/3) - 2N/3] = -F
        - coefQmu * dcmplx(hModL6L(p, iGrid))  ! L
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
        + coefQkappa * dcmplx(hModR5y(p, iGrid)) - coefQmu * dcmplx(2.d0 / 3.d0 * hModR6N(p, iGrid)) &  ! -[(F+2N/3) - 2N/3] = -F
        - coefQmu * dcmplx(hModR4L(p, iGrid))  ! L
    end do
  end do

end subroutine


!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the fluid part.
! This computes the part without largeL coefficients. (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeA0Fluid(nLayerInZoneI, omega, omegaI, p1, p3, coefQfluid, a0Tmp)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: p1(4*nLayerInZoneI), p3(4*nLayerInZoneI)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of T and H matrix stored for each (iLayer, k', k)-pair [TODO].
  complex(8), intent(in) :: coefQfluid  ! Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.
  complex(8), intent(out) :: a0Tmp(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [TODO].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  ! Calculate b0 for each relevant index
  do i = 1, 4 * nLayerInZoneI
    a0Tmp(i) = -dcmplx(p1(i)) / omegaDamped2 + coefQfluid * dcmplx(p3(i))
  end do

end subroutine


!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the fluid part.
! This computes the part with coefficient largeL^2. (See eqs. 2 & 17-18 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeA2Fluid(nLayerInZoneI, omega, omegaI, p2, a2Tmp)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency [1/s] (real and imaginary). Imaginary part is for artificial damping.
  real(8), intent(in) :: p2(4*nLayerInZoneI)  ! Parts of H matrix stored for each (iLayer, k', k)-pair [TODO].
  complex(8), intent(out) :: a2Tmp(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: [TODO].
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping [1/s^2]. (omega - i omega_I)^2.
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  ! Calculate b0 for each relevant index
  do i = 1, 4 * nLayerInZoneI
    a2Tmp(i) = -dcmplx(p2(i)) / omegaDamped2
  end do

end subroutine


!------------------------------------------------------------------------
! Overlapping the coefficient matrix elements for a certain zone in the fluid part.
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!------------------------------------------------------------------------
subroutine overlapAFluid(nLayerInZoneI, aTmp, aOut)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  complex(8), intent(in) :: aTmp(4*nLayerInZoneI)  ! Input tridiagonal matrix, stored for each (iLayer, k', k)-pair [TODO].
  complex(8), intent(inout) :: aOut(4, 2*nLayerInZoneI+2)  ! Upper band of the overlapped matrix [TODO].
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


!------------------------------------------------------------------------
! Assembling the coefficient matrix 'A' for the whole region from several parts.
! The results are the components in the upper band of the A matrix,
!  stored for each (iRow, iColumn) = (1,1), (1,2),(2,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4),(4,4), (2,5),(3,5),(4,5),(5,5), ...
!------------------------------------------------------------------------
subroutine assembleAWhole(nZone, phaseOfZone, oColumnOfZone, largeL2, a0, a1, a2, a)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nZone  ! Number of zones.
  integer, intent(in) :: phaseOfZone(nZone)  ! Phase of each zone (1: solid, 2: fluid).
  integer, intent(in) :: oColumnOfZone(nZone+1)  ! Index of the first column in the band matrix for each zone.
  real(8), intent(in) :: largeL2  ! L^2 = l(l+1).
  complex(8), intent(in) :: a0(4, *), a1(4, *), a2(4, *)  ! Parts of the A matrix, containing upper band elements.
  complex(8), intent(out) :: a(4, *)  ! Assembled A matrix, containing upper band elements.
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



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------



