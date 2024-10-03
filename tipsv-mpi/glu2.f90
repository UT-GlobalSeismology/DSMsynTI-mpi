
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Program to solve the simultaneous linear equation Ac=g using the Gaussian elimination method.
! This code is written based on glu2.f created by H. Hasegawa in 1989.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------------------------------------------------------------
subroutine GLU(A, N, N1, B, EPS, WK, IP, IER)
!------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: N, N1  ! Order of matrix and size of matrix A.
  complex(8), intent(inout) :: A(N1, N)  ! Coefficient matrix. At the end, this will contain result of Gaussian elimination.
  complex(8), intent(inout) :: B(N)  ! Right-hand side vector. At the end, this will contain the solution.
  real(8), intent(inout) :: EPS  ! Parameter to check the singularity of the matrix.
  complex(8), intent(out) :: WK(N)  ! Working array.
  integer, intent(out) :: IP(N)  ! Pivot number.
  integer, intent(out) :: IER  ! 0: normal execution, 1: singular matrix, 2: singular original matrix, 3: invalid argument.
  integer :: I, J, K, IPK
  complex(8) :: AMAX, AIK, W, T

  ! Set the default value for EPS if not provided correctly.
  if (EPS < 0.0d0) EPS = 3.52d-15

  ! Validate the input dimensions.
  if (N1 < N .or. N <= 0) then
    IER = 3
    write(*,*) '  (SUBR. GLU)  INVALID ARGUMENT.  N1, N =', N1, N
    return
  end if

  ! Check the original matrix.
  do I = 1, N
    WK(I) = cdabs(A(I, 1))
  end do
  do J = 2, N
    do I = 1, N
      WK(I) = max(cdabs(WK(I)), cdabs(A(I, J)))
    end do
  end do
  do I = 1, N
    if (cdabs(WK(I)) < EPS) then
      IER = 2
      write(*,*) '  (SUBR. GLU)  ORIGINAL MATRIX IS SINGULAR.'
      return
    end if
  end do

  ! Initialize IER to 0 (indicating no errors initially).
  IER = 0

  do K = 1, N
    ! Find maximum element in the K-th column.
    AMAX = cdabs(A(K, K))
    IPK = K

    do I = K+1, N
      AIK = cdabs(A(I, K))
      if (cdabs(AIK) > cdabs(AMAX)) then
        IPK = I
        AMAX = AIK
      end if
    end do

    IP(K) = IPK

    if (cdabs(AMAX) > EPS) then
      ! If the maximum is not on the diagonal, swap rows
      if (IPK /= K) then
        W = A(IPK, K)
        A(IPK, K) = A(K, K)
        A(K, K) = W
      end if

      ! Compute ALFA
      do I = K+1, N
        A(I, K) = -A(I, K) / A(K, K)
        WK(I) = A(I, K)
      end do

      ! Perform Gaussian Elimination
      do J = K+1, N
        if (IPK /= K) then
          W = A(IPK, J)
          A(IPK, J) = A(K, J)
          A(K, J) = W
        end if

        T = A(K, J)
        do I = K+1, N
          A(I, J) = A(I, J) + WK(I) * T
        end do
      end do

    else
      ! Matrix is singular
      IER = 1
      IP(K) = K

      do I = K+1, N
        A(I, K) = 0.d0
      end do

      write(*,*) '  (SUBR. GLU)  MATRIX IS SINGULAR AT K =', K
      return
    end if
  end do

!------------------------------------------------------------------------
 ENTRY GLUSUB( A, N, N1, B, EPS, WK, IP, IER )
!------------------------------------------------------------------------

  ! FORWARD ELIMINATION PROCESS
  do K = 1, N
    if (IP(K) /= K) then
      W = B(IP(K))
      B(IP(K)) = B(K)
      B(K) = W
    end if

    T = B(K)
    do I = K+1, N
      B(I) = B(I) + A(I, K) * T
    end do
  end do

  ! BACKWARD SUBSTITUTION PROCESS
  B(N) = B(N) / A(N, N)
  do K = N-1, 1, -1
    T = B(K + 1)
    do I = 1, K
      B(I) = B(I) - A(I, K + 1) * T
    end do
    B(K) = B(K) / A(K, K)
  end do

end subroutine GLU
