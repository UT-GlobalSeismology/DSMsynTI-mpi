
!------------------------------------------------------------------------
! Computing \int dr con r^rpow X_k1^(dot1) X_k2^(dot2). (See eq. 16 of Kawai et al. 2006.)
! "I_(k'k)^4" is replaced with "I_(k'k)^4 + I_(kk')^4". (See eq. 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeIntermediateIntegral(nLayerInZoneI, nValue, valuedRadii, con, rpow, dot1, dot2, gridRadiiInZoneI, mat, work)
!------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxrpow = 2  ! Maximum value of rpow to allow.

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: nValue  ! Total number of values for each variable.
  integer, intent(in) :: rpow  ! The exponent of r.
  integer, intent(in) :: dot1, dot2  ! Whether or not to differentiate X_k1 and X_k2 (1: differentiate, 0: do not differentiate).
  real(8), intent(in) :: valuedRadii(nValue)  ! Radii corresponding to each variable value.
  real(8), intent(in) :: con(nValue)  ! Values of a variable at each point (with 2 values at boundaries).
  real(8), intent(in) :: gridRadiiInZoneI(nLayerInZoneI+1)  ! Radii of grid points in zone of interest.
  real(8), intent(out) :: mat(4*nLayerInZoneI)  ! Resulting integrals, "I_(k'k)^4" replaced with "I_(k'k)^4 + I_(kk')^4"
  !::::::::::::::::::::::::::::::::::::::::::::::: (See eq. 19 of Kawai et al. 2006.)
  real(8), intent(out) :: work(4*nLayerInZoneI)  ! Resulting integrals.
  integer :: iLayer, j, k, l, nn
  real(8) :: a(2,2), b(2,2), c(5), rh

  ! parameter check
  if (rpow > maxrpow) stop "Invalid arguments.(computeIntermediateIntegral)"

  ! computing matrix elements
  do iLayer = 1, nLayerInZoneI
    ! layer thickness
    rh = gridRadiiInZoneI(iLayer + 1) - gridRadiiInZoneI(iLayer)

    ! set X_k1^(dot1), for both k1=i and k1=i+1
    select case(dot1)
     case (0)
      a(:, 1) = [gridRadiiInZoneI(iLayer + 1) / rh, -1.d0 / rh]
      a(:, 2) = [-gridRadiiInZoneI(iLayer) / rh, 1.d0 / rh]
     case (1)
      a(:, 1) = [-1.d0 / rh, 0.d0]
      a(:, 2) = [1.d0 / rh, 0.d0]
     case default
      stop "Invalid arguments.(computeIntermediateIntegral)"
    end select

    ! set X_k2^(dot2), for both k2=i and k2=i+1
    select case(dot2)
     case (0)
      b(:, 1) = [gridRadiiInZoneI(iLayer + 1) / rh, -1.d0 / rh]
      b(:, 2) = [-gridRadiiInZoneI(iLayer) / rh, 1.d0 / rh]
     case (1)
      b(:, 1) = [-1.d0 / rh, 0.d0]
      b(:, 2) = [1.d0 / rh, 0.d0]
     case default
      stop "Invalid arguments.(computeIntermediateIntegral)"
    end select

    do j = 1, 2  ! k1=i and k1=i+1
      do k = 1, 2  ! k2=i and k2=i+1
        c = 0.d0
        ! multiply X_k1^(dot1) and X_k2^(dot2)
        call multiplyPolynomials(2, a(:, j), 2, b(:, k), 3, c(:))
        ! multiply by r^rpow
        if (rpow > 0) then
          do l = 3, 1, -1
            c(l + rpow) = c(l)
            c(l) = 0.d0
          end do
        end if
        ! integrate; the result is saved for each (iLayer, k1, k2)-pair
        nn = 4 * (iLayer - 1) + 2 * (j - 1) + k
        call integrateProduct(5, c, gridRadiiInZoneI(iLayer), gridRadiiInZoneI(iLayer + 1), nValue, valuedRadii(:), con(:), &
          work(nn))
      end do
    end do
  end do

  ! replace "I_(k'k)^4" with "I_(k'k)^4 + I_(kk')^4" (See eq. 19 of Kawai et al. 2006.)
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
! Computing the (l-1)-degree polynomial c(x) which is the product of
! the (n-1)-degree polynomial a(x) and the (m-1)-degree polynomial b(x).
!------------------------------------------------------------------------
subroutine multiplyPolynomials(n, a, m, b, l, c)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: n, m, l  ! Size of the arrays of a, b, and c.
  real(8), intent(in) :: a(n), b(m)  ! Coefficients of polynimials a and b in ascending order (a(x) = a1 + a2 x + a3 x^2 + ...).
  real(8), intent(out) :: c(l)  ! Coefficients of polynimial c in ascending order.

  integer :: i, j

  ! Check for invalid arguments
  if (n + m - 1 /= l) stop "Invalid arguments.(multiplyPolynomials)"

  ! Initialize the polynomial c
  c(:) = 0.d0

  ! Compute the product polynomial
  do i = 1, n
    do j = 1, m
      c(i + j - 1) = c(i + j - 1) + a(i) * b(j)
    end do
  end do

end subroutine


!------------------------------------------------------------------------
! Evaluating the integrated value of p(r)*con(r) from 'lowerX' to 'upperX'.
! Here, p(r) is an (n-1)-degree polynomial, and con(r) is the profile of a variable.
!------------------------------------------------------------------------
subroutine integrateProduct(n, p, lowerRadius, upperRadius, nValue, valuedRadii, con, result)
!------------------------------------------------------------------------
  implicit none

  integer, parameter :: maxn = 5  ! Maximum number of polynomial degrees.
  integer, intent(in) :: n  ! Size of the array of p.
  real(8), intent(in) :: p(n)  ! Coefficients of the polynimial in ascending order (p(r) = p1 + p2 r + p3 r^2 + ...).
  real(8), intent(in) :: lowerRadius, upperRadius  ! Radius range to integrate.
  integer, intent(in) :: nValue  ! Total number of values for each variable.
  real(8), intent(in) :: valuedRadii(nValue)  ! Radii corresponding to each variable value.
  real(8), intent(in) :: con(nValue)  ! Values of a variable at each point (with 2 values at boundaries).
  real(8), intent(out) :: result
  real(8) :: r1, r2, q(2), pq(maxn+1), dS
  integer :: iValue

  ! Check the number of polynomial degrees
  if (n > maxn) stop 'Degree of polynomial is too large.(integrateProduct)'

  ! Initialization
  result = 0.d0
  r1 = lowerRadius
  iValue = 1

  do
    ! find the first value within range of integration
    if (valuedRadii(iValue + 1) <= r1) then
      iValue = iValue + 1
      cycle
    end if

    r2 = min(upperRadius, valuedRadii(iValue + 1))

    ! express con(r) as a polynomial (linear) function q1+q2*r
    q(2) = (con(iValue + 1) - con(iValue)) / (valuedRadii(iValue + 1) - valuedRadii(iValue))  ! slope
    q(1) = con(iValue) - q(2) * valuedRadii(iValue)  ! intercept
    ! compute p(r)*con(r)
    call multiplyPolynomials(n, p(:), 2, q(:), n + 1, pq(:))
    ! evaluate integrated value within subrange [r1,r2]
    call integratePolynomial(n + 1, pq(:), r1, r2, dS)
    result = result + dS

    if (r2 == upperRadius) exit
    r1 = r2
  end do

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

  ! Check the number of polynomial degrees.
  if (n > maxn) stop 'Degree of polynomial is too large.(integratePolynomial)'

  ! Initialization: a=(1 x1 x1^2 ...), b=(1 x2 x2^2 ...)
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
  ! loop for each term of p(x)
  do i = 1, n
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
! Computing the lumped mass matrix for a certain zone in the solid part. (See eq. 15 of Cummins et al. 1994.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeLumpedT(nLayerInZoneI, nValue, valuedRadii, rhoValues, gridRadiiInZoneI, tl)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: nValue  ! Total number of values for each variable.
  real(8), intent(in) :: valuedRadii(nValue)  ! Radii corresponding to each variable value.
  real(8), intent(in) :: rhoValues(nValue)  ! Rho values at each point (with 2 values at boundaries).
  real(8), intent(in) :: gridRadiiInZoneI(nLayerInZoneI+1)  ! Radii of grid points in zone of interest.
  real(8), intent(out) :: tl(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  integer :: i, nn
  real(8) :: c(3), lowerRadius, upperRadius

  ! Initialization
  c = [0.d0, 0.d0, 1.d0]

  do i = 1, nLayerInZoneI
    lowerRadius = gridRadiiInZoneI(i)
    upperRadius = (gridRadiiInZoneI(i) + gridRadiiInZoneI(i + 1)) / 2.d0
    nn = 4 * (i - 1)

    call integrateProduct(3, c(:), lowerRadius, upperRadius, nValue, valuedRadii(:), rhoValues(:), tl(nn + 1))

    tl(nn + 2) = 0.d0
    tl(nn + 3) = 0.d0

    lowerRadius = upperRadius
    upperRadius = gridRadiiInZoneI(i + 1)

    call integrateProduct(3, c(:), lowerRadius, upperRadius, nValue, valuedRadii(:), rhoValues(:), tl(nn + 4))
  end do

end subroutine


!------------------------------------------------------------------------
! Computing the lumped rigidity matrix for a certain zone in the solid part. (See eq. 15 of Cummins et al. 1994.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeLumpedH(nLayerInZoneI, nValue, valuedRadii, muValues, gridRadiiInZoneI, hl)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: nValue  ! Total number of values for each variable.
  real(8), intent(in) :: valuedRadii(nValue)  ! Radii corresponding to each variable value.
  real(8), intent(in) :: muValues(nValue)  ! Mu values at each point (with 2 values at boundaries).
  real(8), intent(in) :: gridRadiiInZoneI(nLayerInZoneI+1)  ! Radii of grid points in zone of interest.
  real(8), intent(out) :: hl(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  integer :: i, nn
  real(8) :: c(1), lowerRadius, upperRadius

  ! Initialization
  c = [1.d0]

  do i = 1, nLayerInZoneI
    lowerRadius = gridRadiiInZoneI(i)
    upperRadius = (gridRadiiInZoneI(i) + gridRadiiInZoneI(i + 1)) / 2.d0
    nn = 4 * (i - 1)

    call integrateProduct(1, c(:), lowerRadius, upperRadius, nValue, valuedRadii(:), muValues(:), hl(nn + 1))

    hl(nn + 2) = 0.d0
    hl(nn + 3) = 0.d0

    lowerRadius = upperRadius
    upperRadius = gridRadiiInZoneI(i + 1)

    call integrateProduct(1, c(:), lowerRadius, upperRadius, nValue, valuedRadii(:), muValues(:), hl(nn + 4))
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
subroutine computeA0(nLayerInZoneI, omega, omegaI, t, h1, h2, h3, h4, coef, a0)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency (real and imaginary parts). Imaginary part is for artificial damping.
  real(8), intent(in) :: t(4*nLayerInZoneI)  ! T matrix stored for each (iLayer, k', k)-pair.
  real(8), intent(in) :: h1(4*nLayerInZoneI), h2(4*nLayerInZoneI), h3(4*nLayerInZoneI), h4(4*nLayerInZoneI)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair.
  complex(8), intent(in) :: coef
  complex(8), intent(out) :: a0(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping. (omega - i omega_I)^2.
  real(8) :: h
  integer :: i

  ! Introduce artificial damping into angular frequency. (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  do i = 1, 4 * nLayerInZoneI
    ! I2 - I4 - I4' + I6 - 2*I7. (See eq. 19 of Kawai et al. 2006.)
    h = h1(i) - h2(i) + h3(i) - 2.d0 * h4(i)
    ! omega^2 T - (I2 - I4 - I4' + I6 - 2*I7). (See eq. 2 of Kawai et al. 2006.)
    a0(i) = omegaDamped2 * dcmplx(t(i)) - coef * dcmplx(h)
  end do

end subroutine


!------------------------------------------------------------------------
! Computing part of the coefficient matrix 'A' for a certain zone in the solid part.
! This computes (- (I7)). (See eq. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeA2(nLayerInZoneI, h4, coef, a2)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: h4(4*nLayerInZoneI)  ! Part of H matrix stored for each (iLayer, k', k)-pair.
  complex(8), intent(in) :: coef
  complex(8), intent(out) :: a2(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  integer :: i

  do i = 1, 4 * nLayerInZoneI
    ! -(I7). (See eq. 2 & 19 of Kawai et al. 2006.)
    a2(i) = - coef * dcmplx(h4(i))
  end do

end subroutine


!------------------------------------------------------------------------
! Assembling the coefficient matrix 'A' in the solid part from several parts.
!------------------------------------------------------------------------
subroutine assembleA(nGrid, l, a0, a2, a)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nGrid  ! Total number of grid points.
  integer, intent(in) :: l  ! Angular order.
  complex(8), intent(in) :: a0(2,nGrid), a2(2,nGrid)  ! Parts of the A matrix, in diagonal and subdiagonal component format.
  complex(8), intent(out) :: a(2,nGrid)  ! Assembled A matrix.
  integer :: i, j

  do j = 1, nGrid
    do i = 1, 2
      a(i, j) = a0(i, j) + dcmplx(l * (l + 1)) * a2(i, j)
    end do
  end do

end subroutine


!------------------------------------------------------------------------
! Computing the coefficient matrix 'A' for a certain zone in the solid part. (See eq. 2 & 19 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!------------------------------------------------------------------------
subroutine computeA(nLayerInZoneI, omega, omegaI, l, t, h1, h2, h3, h4, coef, a)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: l  ! Angular order.
  real(8), intent(in) :: omega, omegaI  ! Angular frequency (real and imaginary parts). Imaginary part is for artificial damping.
  real(8), intent(in) :: t(4*nLayerInZoneI)  ! T matrix stored for each (iLayer, k', k)-pair.
  real(8), intent(in) :: h1(4*nLayerInZoneI), h2(4*nLayerInZoneI), h3(4*nLayerInZoneI), h4(4*nLayerInZoneI)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Parts of H matrix stored for each (iLayer, k', k)-pair.
  complex(8), intent(in) :: coef
  complex(8), intent(out) :: a(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  complex(8) :: omegaDamped2  ! Squared angular frequency with artificial damping. (omega - i omega_I)^2.
  real(8) :: h
  integer :: i

  ! introduce artificial damping into angular frequency (See section 5.1 of Geller & Ohminato 1994.)
  omegaDamped2 = dcmplx(omega, -omegaI) ** 2

  do i = 1, 4 * nLayerInZoneI
    ! compute I2 - I4 - I4' + I6 + (L^2 - 2)*I7 (See eq. 19 of Kawai et al. 2006.)
    h = h1(i) - h2(i) + h3(i) + dble(l * (l + 1) - 2) * h4(i)
    ! compute (omega^2 T - H) (See eq. 2 of Kawai et al. 2006.)
    a(i) = omegaDamped2 * dcmplx(t(i)) - coef * dcmplx(h)
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
  complex(8), intent(in) :: aIn(4*nLayerInZoneI)  ! Tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  complex(8), intent(out) :: aOut(2, nLayerInZoneI+1)  ! Diagonal and subdiagonal components of the overlapped matrix.
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::: Should be initialized with 0s beforehand.
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
subroutine computeG(l, m, iLayerOfSource, r0, mt, mu0, coef, aSourceParts, aaParts, aSource, dr, g)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: l  ! Angular order.
  integer, intent(in) :: m  ! Azimuthal order.
  integer, intent(in) :: iLayerOfSource  ! Which layer the source is in.
  real(8), intent(in) :: r0, mu0, mt(3,3)
  complex(8), intent(in) :: coef  !!TODO probably Coefficient derived from attenuation for each zone.
  complex(8), intent(in) :: aSourceParts(8), aaParts(4), aSource(2,3)
  complex(8), intent(inout) :: dr(3)
  complex(8), intent(out) :: g(*)  ! The vector -g
  real(8) :: b, sgnM
  complex(8) :: dd, gS_or_cS(3)
  integer :: i
  complex(8) :: z(3)
  real(8) :: eps
  real(8) :: ier

  ! Initialize.
  call initComplexVector(3, gS_or_cS(:))
  dd = dcmplx(0.d0, 0.d0)
  eps = -1.d0

  ! Record sign of m.
  if (m >= 0) then
    sgnM = 1.d0
  else
    sgnM = -1.d0
  end if

  if (abs(m) == 1) then
    ! b1 in eq. (26) of Kawai et al. (2006)
    b = sqrt((2 * l + 1) / (16.d0 * pi))
    ! D3 of eq. (26) of Kawai et al. (2006)
    dd = dcmplx(b) * dcmplx(sgnM * mt(1, 3), mt(1, 2)) / (dcmplx(r0 * r0 * mu0) * coef)

    !TODO ??
    do i = 2, 3
      gS_or_cS(i) = -dd * (aSourceParts(i * 2 + 1) + aSourceParts(i * 2 + 2))  ! i=2 -> 5, 6 ; i=3 -> 7, 8
    end do

  else if (abs(m) == 2) then
    ! b2 in eq. (27) of Kawai et al. (2006)
    b = sqrt((2 * l + 1) * (l - 1) * (l + 2) / (64.d0 * pi))
    ! -1 * gk3 of eq. (27) of Kawai et al. (2006)
    gS_or_cS(2) = dcmplx(b / r0) * dcmplx(2.d0 * mt(2, 3), sgnM * (mt(2, 2) - mt(3, 3)))
  end if

  ! Solve Ac=g (i.e. (omega^2 T - H) c = -g) for grids near source.
  if ((m == -2) .or. (m == -l)) then
    ! In the first m-loop (m=-1 for l=1; m=-2 otherwise), matrix A must be decomposed.
    call dclisb0(aSource(:,:), 3, 1, 2, gS_or_cS(:), eps, dr, z, ier)
  else
    ! In consecutive m-loops, start from forward substitution (decomposition is skipped).
    call dcsbsub0(aSource(:,:), 3, 1, 2, gS_or_cS(:), eps, dr, z, ier)
  end if

  ! Add displacement to c.
  gS_or_cS(3) = gS_or_cS(3) + dd

  ! Computate excitation vector g.
  g(iLayerOfSource) = aaParts(1) * gS_or_cS(1) + aaParts(2) * gS_or_cS(3)
  g(iLayerOfSource + 1) = aaParts(3) * gS_or_cS(1) + aaParts(4) * gS_or_cS(3)

end subroutine


