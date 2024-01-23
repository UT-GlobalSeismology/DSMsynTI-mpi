
!----------------------------------------------------------------------------------------------------------------------------
! Converts geodetic latitude to geocentric latitude.
!----------------------------------------------------------------------------------------------------------------------------
subroutine transformLatitude(geodetic, geocentric)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: flattening = 1.d0 / 298.25d0
  real(8), parameter :: pi = 3.1415926535897932d0

  real(8), intent(in) :: geodetic  ! Input geodetic latitude [deg].
  real(8), intent(out) :: geocentric  ! Output geocentric latitude [deg].
  real(8) :: latitude  ! Latitude variable for use in computation.

  if (geodetic < -90.d0 .or. 90.d0 < geodetic) stop 'Latitude is out of range. (transformLatitude)'

  ! degrees to radians
  latitude = geodetic / 180.d0 * pi
  ! gedetic to geocentric
  latitude = atan((1.d0 - flattening) * (1.d0 - flattening) * tan(latitude))
  ! radians to degrees
  geocentric = latitude * 180.d0 / pi

  return
end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computes the colatitude (theta) and longitude (phi) of a receiver
! when the source is shifted to the north pole.
! Note that the longitude of the original source is set as 0 after the shift,
! so the shifted longitude of receiver [rad] is (pi - azimuth).
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeThetaPhi(iEvLat, iEvLon, iStLat, iStLon, theta, phi)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  real(8), intent(in) :: iEvLat, iEvLon, iStLat, iStLon  ! Input latitudes and longitudes of source and receiver [deg].
  real(8), intent(out) :: theta, phi  ! Colatitude and longitude of receiver with event at north pole [rad].
  real(8) :: evColat, evLon, stColat, stLon  ! Colatitudes and longitudes of source and receiver [rad].
  real(8) :: cosAlpha, sinAlpha
  real(8) :: tmp

  ! Transform geographic latitudes [deg] to geocentric colatitudes [rad].
  call transformLatitude(iEvLat, tmp)
  evColat = (90.d0 - tmp) / 180.d0 * pi
  call transformLatitude(iStLat, tmp)
  stColat = (90.d0 - tmp) / 180.d0 * pi

  ! Transform longitudes from degrees to radians.
  evLon = iEvLon / 180.d0 * pi
  stLon = iStLon / 180.d0 * pi

  ! Compute epicentral distance [rad], which will directly be the colatitude of receiver after shift.
  cosAlpha = cos(evColat) * cos(stColat) + sin(evColat) * sin(stColat) * cos(evLon - stLon)
  if (1.d0 < cosAlpha) cosAlpha = 1.d0
  if (cosAlpha < -1.d0) cosAlpha = -1.d0
  theta = acos(cosAlpha)

  ! Compute shifted longitude of receiver [rad], which is (pi - azimuth).
  if (sin(theta) == 0.d0) then
    phi = 0.d0
  else
    cosAlpha = (cos(stColat) * sin(evColat) - sin(stColat) * cos(evColat) * cos(stLon - evLon)) / sin(theta)
    if (1.d0 < cosAlpha) cosAlpha = 1.d0
    if (cosAlpha < -1.d0) cosAlpha = -1.d0
    sinAlpha = sin(stColat) * sin(stLon - evLon) / sin(theta)
    ! pi - azimuth
    if (sinAlpha >= 0.d0) then
      phi = pi - acos(cosAlpha)
    else
      phi = pi + acos(cosAlpha)
    end if
  end if

  return
end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Function to compute a + bx + cx^2 + dx^3, where x = r/R.
!----------------------------------------------------------------------------------------------------------------------------
subroutine valueAtRadius(coefficients, radius, rmax, result)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: coefficients(4)  ! Coefficients of cubic function. [a, b, c, d] in a + bx + cx^2 + dx^3.
  real(8), intent(in) :: radius  ! r : The radius to compute the value at [km].
  real(8), intent(in) :: rmax  ! R: Maximum radius of region considered [km].
  real(8), intent(out) :: result
  integer :: j
  real(8) :: x_n  ! Power of x = r/R.
  real(8) :: accumulatedValue  ! Variable to store the temporary result of accumulation.

  x_n = 1.d0
  accumulatedValue = coefficients(1)
  do j = 2, 4
    x_n = x_n * (radius / rmax)
    accumulatedValue = accumulatedValue + coefficients(j) * x_n
  end do

  result = accumulatedValue
  return
end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computes the accuracy threshold of angular order that is sufficient to compute the slowest phase velocity.
! (See eq. 29 of Kawai et al. 2006.)
! This corresponds to l_d in Kawai et al. (2006).
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeLsuf(omega, nZone, rmaxOfZone, vsvPolynomials, lsuf)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  real(8), intent(in) :: omega  ! Angular frequency.
  integer, intent(in) :: nZone  ! Number of zones.
  real(8), intent(in) :: rmaxOfZone(nZone)  ! Upper radii of each zone [km].
  real(8), intent(in) :: vsvPolynomials(4,nZone)  ! Polynomial functions of vsv structure [km/s].
  integer, intent(out) :: lsuf  ! Accuracy threshold of angular order.
  real(8) :: vsAtSurface

  ! Compute Vs at planet surface [km/s].
  call valueAtRadius(vsvPolynomials(:, nZone), 1.d0, 1.d0, vsAtSurface)

  ! Compute lsuf. (See eq. 29 of Kawai et al. 2006.)
  !  The slowest velocity (vs at surface) and largest radius (planet radius) is used to gain larger bound of angular order.
  lsuf = int(omega * rmaxOfZone(nZone) / vsAtSurface - 0.5d0) + 1

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing \int var r^rpow X_k1^(dot1) X_k2^(dot2) dr. (See eqs. 16 and 20 of Kawai et al. 2006.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeIntermediateIntegral(nLayerInZoneI, valuedRadiiInZoneI, valuesInZoneI, rpow, dot1, dot2, mat)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxrpow = 2  ! Maximum value of rpow to allow.

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: valuesInZoneI(nLayerInZoneI+1)  ! Values of a variable at each point (with 2 values at boundaries).
  integer, intent(in) :: rpow  ! The exponent of r.
  integer, intent(in) :: dot1, dot2  ! Whether or not to differentiate X_k1 and X_k2 (1: differentiate, 0: do not differentiate).
  real(8), intent(out) :: mat(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  !:::::::::::::::::::::::::::: For solid, I^0: [10^12 kg], others: [10^12 kg/s^2]. For fluid, I^F2: [m^5/N], others: [m^4/kg].
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


!----------------------------------------------------------------------------------------------------------------------------
! Computing the (nc-1)-degree polynomial c(x) which is the product of
! the (na-1)-degree polynomial a(x) and the (nb-1)-degree polynomial b(x).
!----------------------------------------------------------------------------------------------------------------------------
subroutine multiplyPolynomials(na, a, nb, b, nc, c)
!----------------------------------------------------------------------------------------------------------------------------
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


!----------------------------------------------------------------------------------------------------------------------------
! Evaluating the integrated value of p(r)*var(r) from 'lowerRadius' to 'upperRadius'.
! Here, p(r) is an (n-1)-degree polynomial, and var(r) is the profile of a variable.
! The range [lowerRadius, upperRadius] must be within a certain layer [valuedRadii(1), valuedRadii(2)].
!----------------------------------------------------------------------------------------------------------------------------
subroutine integrateProduct(n, p, lowerRadius, upperRadius, valuedRadii, values, result)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxn = 5  ! Maximum number of polynomial degrees.

  integer, intent(in) :: n  ! Size of the array of p.
  real(8), intent(in) :: p(n)  ! Coefficients of the polynimial in ascending order (p(r) = p1 + p2 r + p3 r^2 + ...).
  real(8), intent(in) :: lowerRadius, upperRadius  ! Radius range to integrate [km].
  real(8), intent(in) :: valuedRadii(2)  ! Radii at both ends of an interval containing integration range [km].
  real(8), intent(in) :: values(2)  ! Values of a variable at both ends of an interval containing integration range.
  real(8), intent(out) :: result
  real(8) :: q(2), pq(maxn+1)

  ! Check input validity.
  if (n > maxn) stop 'Degree of polynomial is too large. (integrateProduct)'

  ! Express var(r) as a polynomial (linear) function q1+q2*r.
  q(2) = (values(2) - values(1)) / (valuedRadii(2) - valuedRadii(1))  ! slope
  q(1) = values(1) - q(2) * valuedRadii(1)  ! intercept
  ! Compute p(r)*var(r).
  call multiplyPolynomials(n, p(:), 2, q(:), n + 1, pq(:))
  ! Evaluate integrated value within subrange [lowerRadius, upperRadius].
  call integratePolynomial(n + 1, pq(:), lowerRadius, upperRadius, result)

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Evaluating the integrated value of an (n-1)-degree polynomial 'p(x)' from 'x1' to 'x2'.
!----------------------------------------------------------------------------------------------------------------------------
subroutine integratePolynomial(n, p, x1, x2, result)
!----------------------------------------------------------------------------------------------------------------------------
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


!----------------------------------------------------------------------------------------------------------------------------
! Computing the lumped mass matrix for a certain zone. (See eqs. 15-17 of Cummins et al. 1994.)
!  T_kk^lumped = m_k r_k^2 = \int var r^2 dr,
!  T_k'k^lumped = 0 (when k' /= k).
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeLumpedT(nLayerInZoneI, valuedRadiiInZoneI, valuesInZoneI, tl)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: valuesInZoneI(nLayerInZoneI+1)  ! Values of a variable at each point (with 2 values at boundaries).
  real(8), intent(out) :: tl(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: For solid, [10^12 kg]. For fluid, [m^5/N].
  integer :: i, nn
  real(8) :: c(3), lowerRadius, upperRadius

  ! Initialize. This is c(r) = r^2.
  c = [0.d0, 0.d0, 1.d0]

  do i = 1, nLayerInZoneI
    nn = 4 * (i - 1)

    ! Right side of m_k r_k^2 for k=i. Integrate rho*r^2.
    lowerRadius = valuedRadiiInZoneI(i)
    upperRadius = (valuedRadiiInZoneI(i) + valuedRadiiInZoneI(i + 1)) / 2.d0
    call integrateProduct(3, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), valuesInZoneI(i:), tl(nn + 1))

    tl(nn + 2) = 0.d0
    tl(nn + 3) = 0.d0

    ! Left side of m_k r_k^2 for k=i+1. Integrate rho*r^2.
    lowerRadius = upperRadius
    upperRadius = valuedRadiiInZoneI(i + 1)
    call integrateProduct(3, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), valuesInZoneI(i:), tl(nn + 4))
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Computing the lumped rigidity matrix for a certain zone. (See eqs. 15-17 of Cummins et al. 1994.)
!  H_kk^lumped = s_k = \int var dr,
!  H_k'k^lumped = 0 (when k' /= k).
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine computeLumpedH(nLayerInZoneI, valuedRadiiInZoneI, valuesInZoneI, hl)
!----------------------------------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  real(8), intent(in) :: valuedRadiiInZoneI(nLayerInZoneI+1)  ! Radii corresponding to each variable value [km].
  real(8), intent(in) :: valuesInZoneI(nLayerInZoneI+1)  ! Values of a variable at each point (with 2 values at boundaries).
  real(8), intent(out) :: hl(4*nLayerInZoneI)  ! Resulting tridiagonal matrix, stored for each (iLayer, k', k)-pair.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: For solid, [10^12 kg/s^2]. For fluid, [m^4/kg].
  integer :: i, nn
  real(8) :: c(1), lowerRadius, upperRadius

  ! Initialize. This is c(r) = 1 (constant).
  c = [1.d0]

  do i = 1, nLayerInZoneI
    nn = 4 * (i - 1)

    ! Right side of s_k for k=i. Integrate elastic modulus.
    lowerRadius = valuedRadiiInZoneI(i)
    upperRadius = (valuedRadiiInZoneI(i) + valuedRadiiInZoneI(i + 1)) / 2.d0
    call integrateProduct(1, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), valuesInZoneI(i:), hl(nn + 1))

    hl(nn + 2) = 0.d0
    hl(nn + 3) = 0.d0

    ! Left side of s_k for k=i+1. Integrate elastic modulus.
    lowerRadius = upperRadius
    upperRadius = valuedRadiiInZoneI(i + 1)
    call integrateProduct(1, c(:), lowerRadius, upperRadius, valuedRadiiInZoneI(i:), valuesInZoneI(i:), hl(nn + 4))
  end do

end subroutine


!----------------------------------------------------------------------------------------------------------------------------
! Averaging the values of two tridiagonal matrices for a certain zone. (See eq. 17 of Cummins et al. 1994.)
! The result is a tridiagonal matrix,
!  stored for each (iLayer, k', k) = (1,1,1),(1,1,2),(1,2,1),(1,2,2), (2,2,2),(2,2,3),(2,3,2),(2,3,3), ...
!----------------------------------------------------------------------------------------------------------------------------
subroutine averageMatrix(nLayerInZoneI, mat1, mat2, average)
!----------------------------------------------------------------------------------------------------------------------------
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


