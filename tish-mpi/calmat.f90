
!------------------------------------------------------------------------
! Computing \int r^rpow con W_p^(w1dn) W_q^(w2dn) dr. (See eq. 16 of Kawai et al. 2006.)
!------------------------------------------------------------------------
subroutine computeIntermediateIntegral(nLayerInZoneI, nValue, valuedRadii, con, rpow, w1dn, w2dn, gridRadiiInZoneI, m, work)
!------------------------------------------------------------------------
  implicit none
  integer, parameter :: maxrpow = 2  ! Maximum value of rpow to allow.

  integer, intent(in) :: nLayerInZoneI  ! Number of layers in zone of interest.
  integer, intent(in) :: nValue  ! Total number of values for each variable.
  integer, intent(in) :: rpow, w1dn, w2dn
  real(8), intent(in) :: valuedRadii(nValue)  ! Radii corresponding to each variable value.
  real(8), intent(in) :: con(nValue)  ! Values of a variable at each point (with 2 values at boundaries).
  real(8), intent(in) :: gridRadiiInZoneI(nLayerInZoneI+1)  ! Radii of grid points in zone of interest.
  real(8), intent(out) :: m(4*nLayerInZoneI), work(4*nLayerInZoneI)
  integer :: i, j, k, l, nn, snp
  real(8) :: a(2,2), b(2,2), c(5), rh

  ! Parameter check
  if (rpow > maxrpow) stop "Invalid arguments.(calmatc)"

  ! Computing matrix elements
  snp = 1
  do i = 1, nLayerInZoneI
    ! layer thickness
    rh = gridRadiiInZoneI(i+1) - gridRadiiInZoneI(i)

    select case(w1dn)
     case (0)
      a(:, 1) = [gridRadiiInZoneI(i+1)/rh, -1.d0/rh]
      a(:, 2) = [-gridRadiiInZoneI(i)/rh, 1.d0/rh]
     case (1)
      a(:, 1) = [-1.d0/rh, 0.d0]
      a(:, 2) = [1.d0/rh, 0.d0]
     case default
      stop "Invalid arguments.(calmatc)"
    end select

    select case(w2dn)
     case (0)
      b(:, 1) = [gridRadiiInZoneI(i+1)/rh, -1.d0/rh]
      b(:, 2) = [-gridRadiiInZoneI(i)/rh, 1.d0/rh]
     case (1)
      b(:, 1) = [-1.d0/rh, 0.d0]
      b(:, 2) = [1.d0/rh, 0.d0]
     case default
      stop "Invalid arguments.(calmatc)"
    end select

    do j = 1, 2
      do k = 1, 2
        c = 0.d0
        call multiplyPolynomials(2, a(:, j), 2, b(:, k), 3, c)
        do l = 3, 1, -1
          c(l+rpow) = c(l)
          if (rpow > 0) c(l) = 0.d0
        enddo
        nn = 4 * (i-1) + 2 * (j-1) + k
        call pinteg(snp, 5, c, gridRadiiInZoneI(i), gridRadiiInZoneI(i+1), nValue, valuedRadii, con, work(nn))
      enddo
    enddo
  enddo

  if (w1dn /= w2dn) then
    do i = 1, 4*nLayerInZoneI
      select case(mod(i,4))
       case (0, 1)
        m(i) = 2.d0 * work(i)
       case (2)
        m(i) = work(i) + work(i+1)
       case (3)
        m(i) = work(i-1) + work(i)
      end select
    enddo
  else
    m = work
  endif

end subroutine


!------------------------------------------------------------------------
! Computing the (l-1)-degree polynomial c(n) which is the product of
! the (n-1)-degree polynomial a(n) and the (m-1)-degree polynomial b(n).
!------------------------------------------------------------------------
subroutine multiplyPolynomials(n, a, m, b, l, c)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: n, m, l  ! Size of the arrays of a, b, and c.
  real(8), intent(in) :: a(n), b(m)  ! Coefficients of polynimials a and b.
  real(8), intent(out) :: c(l)  ! Coefficients of polynimial c.

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
  integer :: i, nn, snp
  real(8) :: c(3), from, to

  ! Initialization
  snp = 1
  c = [0.d0, 0.d0, 1.d0]

  do i = 1, nLayerInZoneI
    from = gridRadiiInZoneI(i)
    to = (gridRadiiInZoneI(i) + gridRadiiInZoneI(i+1)) / 2.d0
    nn = 4 * (i-1)

    call pinteg(snp, 3, c, from, to, nValue, valuedRadii, rhoValues, tl(nn+1))

    tl(nn+2) = 0.d0
    tl(nn+3) = 0.d0

    from = to
    to = gridRadiiInZoneI(i+1)

    call pinteg(snp, 3, c, from, to, nValue, valuedRadii, rhoValues, tl(nn+4))
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
  integer :: i, nn, snp
  real(8) :: c(1), from, to

  ! Initialization
  snp = 1
  c = [1.d0]

  do i = 1, nLayerInZoneI
    from = gridRadiiInZoneI(i)
    to = (gridRadiiInZoneI(i) + gridRadiiInZoneI(i+1)) / 2.d0
    nn = 4 * (i-1)

    call pinteg(snp, 1, c, from, to, nValue, valuedRadii, muValues, hl(nn+1))

    hl(nn+2) = 0.d0
    hl(nn+3) = 0.d0

    from = to
    to = gridRadiiInZoneI(i+1)

    call pinteg(snp, 1, c, from, to, nValue, valuedRadii, muValues, hl(nn+4))
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
