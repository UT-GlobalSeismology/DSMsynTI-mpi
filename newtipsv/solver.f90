
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Program to solve the simultaneous linear equation Ac=g using the Gauss method.
! This code is written based on dcsymbdl.f created by Fumiko Nagahori in 1991.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------------------------------------------------------------
! Conducts LU decomposition of matrix A.
! This subroutine must be called when using a matrix A for the first time.
! In subsequent usages, this decomposition can be skipped.
!------------------------------------------------------------------------
subroutine decomposeAByGauss(A, M, N, NN, EPS, Z, W, L, LI, LJ, IER)
!------------------------------------------------------------------------
  implicit none

  complex(8), intent(inout) :: A((M+1)*N)
  integer, intent(in) :: M, N, NN
  real(8), intent(in) :: EPS
  integer, intent(inout) :: L(NN), LI(NN), LJ(NN)
  complex(8), intent(inout) :: Z(M+1), W(M+1)
  integer, intent(out) :: IER  ! Error code.
  integer :: I, J, K, IJ, KK, NK, NKK, NKI, MM
  complex(8) :: PIV

  IER = 0
  IJ = 0
  do I = 1, M
    IJ = IJ + 1
    LI(IJ) = I + 1
    LJ(IJ) = 2
    L(IJ) = M * (I - 1) + 1
    do J = I + 1, M
      IJ = IJ + 1
      L(IJ) = L(IJ - 1) + M + 1
      LI(IJ) = LI(IJ - 1) + 1
      LJ(IJ) = LJ(IJ - 1) + 1
    end do
  end do

  MM = (M + 1) * M / 2

  do K = 1, N - M
    NK = (M + 1) * (K - 1) + 1
    if (cdabs(A(NK + M)) < EPS) then
      write(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
      IER = 1
      return
    endif
    PIV = dcmplx(1.d0) / A(NK + M)
    do J = 2, M + 1
      Z(J) = -A(NK + M * J)
      W(J) = A(NK + M * J) * PIV
      A(NK + M * J) = W(J)
    end do

    KK = NK + M + M
    ! VORTION VEC
    do I = 1, MM
      A(KK + L(I)) = A(KK + L(I)) + W(LJ(I)) * Z(LI(I))
    end do
  end do

  do K = N - M + 1, N - 1
    NK = (M + 1) * (K - 1) + 1
    NKK = (M + 1) * K - 1
    if (cdabs(A(NK + M)) < EPS) then
      write(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
      IER = 1
      return
    endif
    PIV = dcmplx(1.d0) / A(NK + M)
    do J = 2, N - K + 1
      Z(J) = -A(NK + M * J)
      W(J) = A(NK + M * J) * PIV
      A(NK + M * J) = W(J)

      NKI = NKK + M * (J - 1)
      do I = 2, J
        A(NKI + I) = A(NKI + I) + W(J) * Z(I)
      end do
    end do
  end do

end subroutine

!------------------------------------------------------------------------
! Solves the simultaneous linear equation Ac=g.
! All values of c are computed. Cut-off grid is not used.
! Matrix A must be decomposed beforehand.
!------------------------------------------------------------------------
subroutine solveWholeCAfterGauss(A, B, M, N, Z)
  implicit none

  complex(8), intent(in) :: A(M+1,N)
  integer, intent(in) :: M, N
  complex(8), intent(inout) :: B(N), Z(N)
  integer :: MM, J, K, I1, J1
  complex(8) :: SUM

  ! Forward substitution
  MM = M + 1
  if (MM < 3) then
    Z(1) = B(1)
    do J = 2, N
      Z(J) = B(J) - A(1,J) * Z(J-1)
    end do
    do J = 1, N
      Z(J) = Z(J) / A(M + 1, J)
    end do
    B(N) = Z(N)
    do J = 1, N-1
      B(N-J) = Z(N-J) - A(1,N-J+1) * B(N-J+1)
    end do

  else
    Z(1) = B(1)
    Z(2) = B(2) - A(M-1,2) * Z(1)

    do J = 3, N
      if (J > MM) then
        I1 = 1
      else
        I1 = MM - J + 1
      end if
      SUM = dcmplx(0.d0, 0.d0)
      do K = I1, MM-1
        SUM = SUM + A(K,J) * Z(J-MM+K)
      end do
      Z(J) = B(J) - SUM
    end do

    do J = 1, N
      Z(J) = Z(J) / A(M + 1, J)
    end do

    B(N) = Z(N)
    B(N-1) = Z(N-1) - A(MM-1,N) * Z(N)

    do J = 3, N
      J1 = N - J + 1
      I1 = 1
      if (J < MM) I1 = MM - J + 1
      SUM = dcmplx(0.d0, 0.d0)
      do K = I1, MM-1
        SUM = SUM + A(K,MM-K+J1) * B(MM-K+J1)
      end do
      B(J1) = Z(J1) - SUM
    end do
  endif

end subroutine


!------------------------------------------------------------------------
! Solves the simultaneous linear equation Ac=g.
! Only the values of c(n-1) and c(n) are computed!!!
! Computation can be sped-up by only considering parts above a certain cut-off grid.
! Matrix A must be decomposed beforehand.
!------------------------------------------------------------------------
subroutine solveSurfaceCAfterGauss(A, B, M, N, NP, Z)
  implicit none

  complex(8), intent(in) :: A(M+1,N)
  integer, intent(in) :: M, N, NP
  complex(8), intent(inout) :: B(N), Z(N)
  integer :: MM, J, K, I1
  complex(8) :: SUM

  ! Forward substitution
  MM = M + 1

  if (MM < 3) then
    Z(NP) = B(NP)
    do J = NP + 1, N
      Z(J) = B(J) - A(1, J) * Z(J - 1)
    end do
    B(N) = Z(N) / A(M + 1, N)

  else
    Z(NP) = B(NP)
    Z(NP + 1) = B(NP + 1) - A(MM - 1, NP + 1) * Z(NP)
    do J = NP + 2, N
      if (J > NP - 1 + MM) then
        I1 = 1
      else
        I1 = NP - 1 + MM - J + 1
      end if
      SUM = dcmplx(0.d0, 0.d0)
      do K = I1, MM - 1
        SUM = SUM + A(K, J) * Z(J - MM + K)
      end do
      Z(J) = B(J) - SUM
    end do

    do J = N - 1, N
      Z(J) = Z(J) / A(M + 1, J)
    end do

    B(N) = Z(N)
    B(N - 1) = Z(N - 1) - A(MM - 1, N) * Z(N)
  end if

end subroutine
