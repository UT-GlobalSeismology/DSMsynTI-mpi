
!------------------------------------------------------------------------
! Computing \int dr con r^rpow X_k1^(dot1) X_k2^(dot2). (See eq. 16 of Kawai et al. 2006.)
! "I_(k'k)^4" is replaced with "I_(k'k)^4 + I_(kk')^4". (See eq. 19 of Kawai et al. 2006.)
!------------------------------------------------------------------------
subroutine computeIntermediateIntegral(nLayerInZoneI, nValue, valuedRadii, con, rpow, dot1, dot2, gridRadiiInZoneI, m, work)
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
  real(8), intent(out) :: m(4*nLayerInZoneI), work(4*nLayerInZoneI)
  integer :: iLayer, j, k, l, nn
  real(8) :: a(2,2), b(2,2), c(5), rh

  ! parameter check
  if (rpow > maxrpow) stop "Invalid arguments.(computeIntermediateIntegral)"

  ! computing matrix elements
  do iLayer = 1, nLayerInZoneI
    ! layer thickness
    rh = gridRadiiInZoneI(iLayer+1) - gridRadiiInZoneI(iLayer)

    ! set X_k1^(dot1), for both k1=i and k1=i+1
    select case(dot1)
     case (0)
      a(:, 1) = [gridRadiiInZoneI(iLayer+1)/rh, -1.d0/rh]
      a(:, 2) = [-gridRadiiInZoneI(iLayer)/rh, 1.d0/rh]
     case (1)
      a(:, 1) = [-1.d0/rh, 0.d0]
      a(:, 2) = [1.d0/rh, 0.d0]
     case default
      stop "Invalid arguments.(computeIntermediateIntegral)"
    end select

    ! set X_k2^(dot2), for both k2=i and k2=i+1
    select case(dot2)
     case (0)
      b(:, 1) = [gridRadiiInZoneI(iLayer+1)/rh, -1.d0/rh]
      b(:, 2) = [-gridRadiiInZoneI(iLayer)/rh, 1.d0/rh]
     case (1)
      b(:, 1) = [-1.d0/rh, 0.d0]
      b(:, 2) = [1.d0/rh, 0.d0]
     case default
      stop "Invalid arguments.(computeIntermediateIntegral)"
    end select

    do j = 1, 2  ! k1=i and k1=i+1
      do k = 1, 2  ! k2=i and k2=i+1
        c = 0.d0
        ! multiply X_k1^(dot1) and X_k2^(dot2)
        call multiplyPolynomials(2, a(:, j), 2, b(:, k), 3, c)
        ! multiply by r^rpow
        if (rpow > 0) then
          do l = 3, 1, -1
            c(l+rpow) = c(l)
            c(l) = 0.d0
          enddo
        endif
        ! integrate; the result is saved for each (iLayer, j, k)-pair
        nn = 4 * (iLayer-1) + 2 * (j-1) + k
        call integrateProduct(1, 5, c, gridRadiiInZoneI(iLayer), gridRadiiInZoneI(iLayer+1), nValue, valuedRadii, con, work(nn))
      enddo
    enddo
  enddo

  ! replace "I_(k'k)^4" with "I_(k'k)^4 + I_(kk')^4" (See eq. 19 of Kawai et al. 2006.)
  if (dot1 /= dot2) then
    do iLayer = 1, 4*nLayerInZoneI
      select case(mod(iLayer,4))
       case (0, 1)  ! (k1, k2) = (i, i), (i+1, i+1) -> double the value
        m(iLayer) = 2.d0 * work(iLayer)
       case (2)  ! (k1, k2) = (i, i+1) -> add "(k1, k2) = (i+1, i)" case
        m(iLayer) = work(iLayer) + work(iLayer+1)
       case (3)  ! (k1, k2) = (i+1, i) -> add "(k1, k2) = (i, i+1)" case
        m(iLayer) = work(iLayer-1) + work(iLayer)
      end select
    enddo
  else
    m = work
  endif

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
  c = 0.d0

  ! Compute the product polynomial
  do i = 1, n
    do j = 1, m
      c(i+j-1) = c(i+j-1) + a(i) * b(j)
    enddo
  enddo

end subroutine


!------------------------------------------------------------------------
! Evaluating the integrated value of p(r)*con(r) from 'lowerX' to 'upperX'.
! Here, p(r) is an (n-1)-degree polynomial, and con(r) is the profile of a variable.
!------------------------------------------------------------------------
subroutine integrateProduct(snpIn, n, p, lowerRadius, upperRadius, nValue, valuedRadii, con, result)
!------------------------------------------------------------------------
  implicit none

  integer, parameter :: maxn = 5  ! Maximum number of polynomial degrees.
  integer, intent(in) :: snpIn  !!!TODO ERASE
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
    call multiplyPolynomials(n, p, 2, q, n+1, pq)
    ! evaluate integrated value within subrange [r1,r2]
    call integratePolynomial(n+1, pq, r1, r2, dS)
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

  ! Check the number of polynomial degrees
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

  ! Evaluate the integrated value
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
! Computing the lumped mass matrix. (See eq. 15 of Cummins et al. 1994.)
!------------------------------------------------------------------------
subroutine computeLumpedT(nLayerInZoneI, nValue, valuedRadii, rhoValues, gridRadiiInZoneI, tl)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: nValue  ! Total number of values for each variable.
  real(8), intent(in) :: valuedRadii(nValue)  ! Radii corresponding to each variable value.
  real(8), intent(in) :: rhoValues(nValue)  ! Rho values at each point (with 2 values at boundaries).
  real(8), intent(in) :: gridRadiiInZoneI(nLayerInZoneI+1)  ! Radii of grid points in zone of interest.
  real(8), intent(out) :: tl(4*nLayerInZoneI)
  integer :: i, nn
  real(8) :: c(3), lowerRadius, upperRadius

  ! Initialization
  c = [0.d0, 0.d0, 1.d0]

  do i = 1, nLayerInZoneI
    lowerRadius = gridRadiiInZoneI(i)
    upperRadius = (gridRadiiInZoneI(i) + gridRadiiInZoneI(i+1)) / 2.d0
    nn = 4 * (i-1)

    call integrateProduct(1, 3, c, lowerRadius, upperRadius, nValue, valuedRadii, rhoValues, tl(nn+1))

    tl(nn+2) = 0.d0
    tl(nn+3) = 0.d0

    lowerRadius = upperRadius
    upperRadius = gridRadiiInZoneI(i+1)

    call integrateProduct(1, 3, c, lowerRadius, upperRadius, nValue, valuedRadii, rhoValues, tl(nn+4))
  enddo

end subroutine


!------------------------------------------------------------------------
! Computing the lumped rigidity matrix. (See eq. 15 of Cummins et al. 1994.)
!------------------------------------------------------------------------
subroutine computeLumpedH(nLayerInZoneI, nValue, valuedRadii, muValues, gridRadiiInZoneI, hl)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: nValue  ! Total number of values for each variable.
  real(8), intent(in) :: valuedRadii(nValue)  ! Radii corresponding to each variable value.
  real(8), intent(in) :: muValues(nValue)  ! Mu values at each point (with 2 values at boundaries).
  real(8), intent(in) :: gridRadiiInZoneI(nLayerInZoneI+1)  ! Radii of grid points in zone of interest.
  real(8), intent(out) :: hl(4*nLayerInZoneI)
  integer :: i, nn
  real(8) :: c(1), lowerRadius, upperRadius

  ! Initialization
  c = [1.d0]

  do i = 1, nLayerInZoneI
    lowerRadius = gridRadiiInZoneI(i)
    upperRadius = (gridRadiiInZoneI(i) + gridRadiiInZoneI(i+1)) / 2.d0
    nn = 4 * (i-1)

    call integrateProduct(1, 1, c, lowerRadius, upperRadius, nValue, valuedRadii, muValues, hl(nn+1))

    hl(nn+2) = 0.d0
    hl(nn+3) = 0.d0

    lowerRadius = upperRadius
    upperRadius = gridRadiiInZoneI(i+1)

    call integrateProduct(1, 1, c, lowerRadius, upperRadius, nValue, valuedRadii, muValues, hl(nn+4))
  enddo

end subroutine


!------------------------------------------------------------------------
! Averaging the values of two matrices. (See eq. 17 of Cummins et al. 1994.)
!------------------------------------------------------------------------
subroutine computeAverage(nLayerInZoneI, v1, v2, average)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: v1(4*nLayerInZoneI), v2(4*nLayerInZoneI)
  real(8), intent(out) :: average(4*nLayerInZoneI)

  integer :: i

  do i = 1, 4*nLayerInZoneI
    average(i) = (v1(i) + v2(i)) / 2.d0
  enddo

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
