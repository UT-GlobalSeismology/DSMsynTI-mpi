
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Program to numerically integrate a system of NEQ first order ordinary differential equations of the form
! DY(I)/DX = F(X, Y(1),..., Y(NEQ))
! by the classical Runge-Kutta formula.
! This code is written based on rk3.f created by M. Sugihara in 1989.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------------------------------------------------------------
SUBROUTINE RK3(NEQ, FUNC, X0, XE, N, Y0, YN, N1, WORK)
!------------------------------------------------------------------------
  implicit none
  ! External function declaration
  external :: FUNC  ! Subroutine FUNC(X,Y,F) to evaluate derivatives F(I)=DY(I)/DX.

  integer, intent(in) :: NEQ  ! Number of equations to be integrated.
  real(8), intent(inout) :: X0  ! Initial value of independent variable.
  real(8), intent(in) :: XE  ! Output point at which the solution is desired.
  integer, intent(in) :: N  ! Number of subintervals to divide the interval (X0, XE).
  !::::::::::::::::::::::::::: The classical Runge-Kutta formula is used in each subinterval.
  integer, intent(in) :: N1
  complex(8), intent(inout) :: Y0(NEQ)  ! Initial value at X0.
  complex(8), intent(out) :: YN(NEQ)  ! Approximate solution at XE.
  complex(8), intent(out) :: WORK(N1, 2)  ! Working array.
  real(8) :: H
  integer :: I, J

  ! Calculate the step size.
  H = (XE - X0) / dble(N)

  ! Main loop for the RK3 method.
  do I = 1, N
    call RKSTEP(NEQ, FUNC, X0, H, Y0, YN, WORK(1, 1), WORK(1, 2))
    X0 = X0 + H
    do J = 1, NEQ
      Y0(J) = YN(J)
    end do
  end do

  ! Set the final value of X0 to XE.
  X0 = XE

end subroutine


!------------------------------------------------------------------------
SUBROUTINE RKSTEP(NEQ, FUNC, X, H, Y0, YN, AK, W)
!------------------------------------------------------------------------
  implicit none
  real(8), parameter :: A2 = 0.5d0, A3 = 0.5d0
  real(8), parameter :: B2 = 0.5d0, B3 = 0.5d0
  real(8), parameter :: C1 = 1.d0 / 6.d0
  real(8), parameter :: C2 = 1.d0 / 3.d0
  real(8), parameter :: C3 = 1.d0 / 3.d0
  real(8), parameter :: C4 = 1.d0 / 6.d0
  ! External function declaration
  external :: FUNC  ! Subroutine FUNC(X,Y,F) to evaluate derivatives F(I)=DY(I)/DX.

  integer, intent(in) :: NEQ  ! Number of equations to be integrated.
  real(8), intent(in) :: X, H
  complex(8), intent(in) :: Y0(NEQ)
  complex(8), intent(out) :: YN(NEQ)
  complex(8), intent(inout) :: AK(NEQ), W(NEQ)
  integer :: I

  call FUNC(X, Y0, AK)
  do I = 1, NEQ
    YN(I) = Y0(I) + dcmplx(H * C1) * AK(I)
  end do

  do I = 1, NEQ
    W(I) = Y0(I) + dcmplx(H * B2) * AK(I)
  end do
  call FUNC(X + A2 * H, W, AK)
  do I = 1, NEQ
    YN(I) = YN(I) + dcmplx(H * C2) * AK(I)
  end do

  do I = 1, NEQ
    W(I) = Y0(I) + dcmplx(H * B3) * AK(I)
  end do
  call FUNC(X + A3 * H, W, AK)
  do I = 1, NEQ
    YN(I) = YN(I) + dcmplx(H * C3) * AK(I)
  end do

  do I = 1, NEQ
    W(I) = Y0(I) + dcmplx(H) * AK(I)
  end do
  call FUNC(X + H, W, AK)
  do I = 1, NEQ
    YN(I) = YN(I) + dcmplx(H * C4) * AK(I)
  end do

end subroutine
