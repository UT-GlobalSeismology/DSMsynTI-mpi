
!------------------------------------------------------------------------
! Evaluating the value of each term of the vector harmonics (fully normalized) for a certain receiver and l value.
! Term 1 is the r-component of S_lm^1.
! Term 2 is the theta-component of S_lm^2, or (-1) times the phi-component of T_lm.
! Term 3 is the phi-component of S_lm^2, or the theta-component of T_lm.
! The coefficient 1/largeL is not multiplied here. (See eq. 12 of Kawai et al. 2006.)
! For each receiver, this subroutine must be called for each l in order, since results from previous l's are referenced.
!------------------------------------------------------------------------
subroutine computeHarmonicsValues(l, theta, phi, plm, harmonicsValues)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: pi = 3.1415926535897932d0

  integer, intent(in) :: l  ! Angular order.
  real(8), intent(in) :: theta, phi  ! Colatitude and longitude of receiver with event at north pole [rad].
  real(8), intent(inout) :: plm(3, 0:3)  ! Values of the associated Legendre polynomials for this receiver.
  !:::::::::::::::::::::::::::::::::::::::: Values for previous l's must be provided.
  !:::::::::::::::::::::::::::::::::::::::: The value for this l will be placed in plm(1,m).
  !:::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8), intent(out) :: harmonicsValues(3, -2:2)  ! Values of vector harmonics terms, computed for each l and receiver.
  !::::::::::::::::::::::::::::::::::::::::::: The coefficient 1/largeL is not multiplied yet.
  !::::::::::::::::::::::::::::::::::::::::::: Arguments: term (1:3), m (-2:2).
  integer :: m, i
  real(8) :: x, fact, coef
  complex(8) :: expimp

  x = cos(theta)

  ! Compute values of the associated Legendre polynomials.
  do m = 0, min(l, 3)
    call computePlm(l, m, x, plm(1, m))
  end do

  ! Compute values of vector harmonics. Computations for m and -m is done at the same time.
  do m = 0, min(l, 2)
    ! Compute (l + |m|)! / (l - |m|)! = (l+m)(l+m-1)...(l-m+1).  (This is 1 when m=0.)
    fact = 1.d0
    if (m /= 0) then
      do i = l - m + 1, l + m
        fact = fact * dble(i)
      end do
    endif
    ! Compute the factor with the square root.
    coef = sqrt(dble(2 * l + 1) / (4.d0 * pi) / fact)
    ! Compute e^{i m phi}.
    expimp = exp(dcmplx(0.d0, dble(m) * phi))

    ! Y_lm(theta,phi) = coef P_lm(x) e^{i m phi}.
    harmonicsValues(1, m) = coef * plm(1, m) * expimp
    harmonicsValues(1, -m) = conjg(harmonicsValues(1, m))
    ! [del Y_lm(theta,phi) / del theta] = coef [ mx/sin(theta) P_lm(x) + P_{l}{m+1}(x) ] e^{i m phi}.
    harmonicsValues(2, m) = coef * (dble(m) * x / sin(theta) * plm(1, m) + plm(1, m+1)) * expimp
    harmonicsValues(2, -m) = conjg(harmonicsValues(2, m))
    ! 1/sin(theta) [del Y_lm(theta,phi) / del phi] = im/sin(theta) coef P_lm(x) e^{i m phi}.
    harmonicsValues(3, m) = dcmplx(0.d0, dble(m)) / sin(theta) * coef * plm(1, m) * expimp
    harmonicsValues(3, -m) = conjg(harmonicsValues(3, m))

    ! When m is odd, sign of P_{l}{-m}(x) must be flipped from that of P_lm(x).
    if (mod(m, 2) == 1) then
      harmonicsValues(1, -m) = -harmonicsValues(1, -m)
      harmonicsValues(2, -m) = -harmonicsValues(2, -m)
      harmonicsValues(3, -m) = -harmonicsValues(3, -m)
    endif
  end do

end subroutine


!------------------------------------------------------------------------
! Compute the value of the associated Legendre polynomial for given degree l, order m (> 0), and value of x.
! The result is stored in plm(1).
! For each x and m, this subroutine must be called for each l in order, since results from previous l's are referenced.
! When l > m, the resulting values for previous l's must be provided in plm(1:3).
!------------------------------------------------------------------------
subroutine computePlm(l, m, x, plm)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: l, m  ! Degree and order (m > 0).
  real(8), intent(in) :: x  ! Value of x.
  real(8), intent(inout) :: plm(3)  ! Values of the associated Legendre polynomial for this receiver and m-value.
  !::::::::::::::::::::::::::::::::::: The value for this l will be computed from values for previous l's and stored in plm(1).
  !::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before).
  integer :: i
  real(8) :: pmm, somx2, fact

  ! Check input validity
  if (m < 0 .or. m > l .or. abs(x) > 1.d0) then
    stop "Invalid arguments. (computePlm)"
  endif

  if (m == l) then
    ! When this is the first time to compute for this m-value, compute directly.
    !  This is when (l,m) = (0,0), (1,1), (2,2), (3,3).
    !  P_00(x) = 1, P_11(x) = -(1-x^2)^{1/2}, P_22(x) = 3(1-x^2), P_33(x) = -15(1-x^2)^{3/2}.
    pmm = 1.d0
    if (m > 0) then
      somx2 = sqrt((1.d0 - x) * (1.d0 + x))
      fact = 1.d0
      do i = 1, m
        pmm = -pmm * fact * somx2
        fact = fact + 2.d0
      end do
    endif
    plm = [pmm, 0.d0, 0.d0]

  else
    ! shift previous values
    plm(3) = plm(2)
    plm(2) = plm(1)

    if (l == m + 1) then
      ! When this is the second time to compute for this m-value, compute from the previous Plm value.
      !  This is when (l,m) = (1,0), (2,1), (3,2), (4,3).
      !  P_{m+1}{m}(x) = x(2m+1)P_mm(x).
      plm(1) = x * dble(2 * m + 1) * plm(2)

    else
      ! When this is after the second time to compute for this m-value, compute from the previous two Plm values.
      !  P{l+1}{m}(x) = [ x(2l+1)P_lm(x) - (l+m)P_{l-1}{m}(x) ] / (l-m+1),
      !  thus P{l}{m}(x) = [ x(2l-1)P_{l-1}{m}(x) - (l+m-1)P_{l-2}{m}(x) ] / (l-m).
      plm(1) = (x * dble(2 * l - 1) * plm(2) - dble(l + m - 1) * plm(3)) / dble(l - m)
    endif
  endif

end subroutine
