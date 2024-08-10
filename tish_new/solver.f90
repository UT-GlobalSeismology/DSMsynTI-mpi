
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Program to solve the simultaneous linear equation Ac=g using the modified Cholesky decomposition,
! for a real symmetric positive definite band matrix.
! This code is written based on dclisb.f created by T. Oguni in 1989.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------------------------------------------------------------
! Conducts modified Cholesky decomposition of matrix A.
! This subroutine must be called when using a matrix A for the first time.
! In subsequent usages, this decomposition can be skipped.
!------------------------------------------------------------------------
subroutine decomposeAByCholesky(a, n, nud, n1, eps, dr, ier)
!------------------------------------------------------------------------
  implicit none

  complex(8), intent(inout) :: a(n1,n)  ! 2-D array containing matrix A.
  integer, intent(in) :: n  ! Order of matrix A.
  integer, intent(in) :: nud  ! Size of band's half width.
  integer, intent(in) :: n1  ! Row size of the array a(:,:).
  real(8), intent(inout) :: eps  ! Parameter to check singurarity of the matrix. When a negative value is set, 1.0d-14.
  complex(8), intent(inout) :: dr(n)  ! Working array.
  integer, intent(out) :: ier  ! Error code.
  complex(8) :: xx, s, sum, au, t
  real(8) :: eps1
  integer :: i, m, j, k1, mj, k

  ! Initialize.
  ier = 0
  eps1 = 1.0d-14
  m = nud + 1

  ! Check input validity.
  if ((n <= 0) .or. (nud <= 0 ) .or. (n1 < m)) then
    ier = 2
    write(*,*) ' n, nud, n1: ', n, nud, n1
    stop 'Invalid argument. (solver)'
  endif
  if (eps <= 0.0) eps = eps1

  ! Modified Cholesky decomposition.
  j = 1
  if (abs(a(m,1)) <= eps) then
    ier = 1
    write(*,*) '(solver) Singular at step # ', j
    return
  endif
  dr(1) = dcmplx(1.0d0) / a(m,1)
  xx = a(m-1,2)
  a(m-1,2) = a(m-1,2) * dr(1)
  s = a(m,2) - xx * a(m-1,2)
  j = 2
  if (abs(s) <= eps) then
    ier = 1
    write(*,*) '(solver) Singular at step # ', j
    return
  endif
  dr(2) = dcmplx(1.0d0) / s

  if (m < 3) then
    do j = 3, n
      xx = a(1,j)
      a(1,j) = xx * dr(j-1)
      s = a(2,j) - xx * a(1,j)

      if (abs(s) <= eps) then
        ier = 1
        write(*,*) '(solver) Singular at step # ', j
        return
      end if

      dr(j) = dcmplx(1.0d0) / s
    end do
  else
    do j = 3, n
      k1 = 1
      if (j >= m) k1 = j - m + 1
      mj = m - j

      do i = k1+1, j-1
        sum = dcmplx(0.0d0)

        do k = k1, i-1
          sum = sum + a(m-i+k,i) * a(mj+k,j)
        end do

        a(mj+i,j) = a(mj+i,j) - sum
      end do

      sum = dcmplx(0.0d0)

      do i = k1, j-1
        xx = a(mj+i,j)
        au = xx * dr(i)
        sum = sum + xx * au
        a(mj+i,j) = au
      end do

      t = a(m,j) - sum

      if (abs(t) <= eps) then
        ier = 1
        write(*,*) '(solver) Singular at step # ', j
        return
      end if

      dr(j) = dcmplx(1.0d0) / t
    end do
  endif

end subroutine


!------------------------------------------------------------------------
! Solves the simultaneous linear equation Ac=g.
! All values of c are computed. Cut-off grid is not used.
! Matrix A must be decomposed beforehand.
!------------------------------------------------------------------------
subroutine solveWholeCAfterCholesky(a, n, nud, n1, b, dr, z)
!------------------------------------------------------------------------
  implicit none

  complex(8), intent(inout) :: a(n1,n)  ! 2-D array containing matrix A.
  integer, intent(in) :: n  ! Order of matrix A.
  integer, intent(in) :: nud  ! Size of band's half width.
  integer, intent(in) :: n1  ! Row size of the array a(:,:).
  complex(8), intent(inout) :: b(n)  ! For input, vector g. For output, vector c.
  complex(8), intent(inout) :: dr(n), z(n)  ! Working arrays.
  complex(8) :: sum
  integer :: m, j, i1, k, j1

  ! forward substitution
  m = nud + 1
  if (m < 3) then
    z(1) = b(1)
    do j = 2, n
      z(j) = b(j) - a(1,j) * z(j-1)
    end do

    do j = 1, n
      z(j) = z(j) * dr(j)
    end do

    b(n) = z(n)
    do j = 1, n-1
      b(n-j) = z(n-j) - a(1,n-j+1) * b(n-j+1)
    end do

  else
    z(1) = b(1)
    z(2) = b(2) - a(m-1,2) * z(1)

    do j = 3, n
      if (j > m) then
        i1 = 1
      else
        i1 = m - j + 1
      end if
      sum = dcmplx(0.0d0)
      do k = i1, m-1
        sum = sum + a(k, j) * z(j-m+k)
      end do
      z(j) = b(j) - sum
    end do

    do j = 1, n
      z(j) = z(j) * dr(j)
    end do

    b(n) = z(n)
    b(n-1) = z(n-1) - a(m-1,n) * z(n)
    do j = 3, n
      j1 = n - j + 1
      if (j < m) then
        i1 = m - j + 1
      else
        i1 = 1
      end if
      sum = dcmplx(0.0d0)
      do k = i1, m-1
        sum = sum + a(k, m-k+j1) * b(m-k+j1)
      end do
      b(j1) = z(j1) - sum
    end do

  endif

  return
end subroutine


!------------------------------------------------------------------------
! Solves the simultaneous linear equation Ac=g.
! Only the value of c(n) is computed!!!
! Computation can be sped-up by only considering parts above a certain cut-off grid.
! Matrix A must be decomposed beforehand.
!------------------------------------------------------------------------
subroutine solveSurfaceCAfterCholesky(a, n, nud, n1, np, b, dr, z)
!------------------------------------------------------------------------
  implicit none

  complex(8), intent(inout) :: a(n1,n)  ! 2-D array containing matrix A.
  integer, intent(in) :: n  ! Order of matrix A.
  integer, intent(in) :: nud  ! Size of band's half width.
  integer, intent(in) :: n1  ! Row size of the array a(:,:).
  integer, intent(in) :: np  ! Index of first grid to consider in computation.
  complex(8), intent(inout) :: b(n)  ! For input, vector g. For output, the value c(n) is stored in last element.
  complex(8), intent(inout) :: dr(n), z(n)  ! Working arrays.
  complex(8) :: sum
  integer :: m, j, i1, k

  ! forward substitution
  m = nud + 1
  if (m < 3) then
    z(np) = b(np)
    do j = np+1, n
      z(j) = b(j) - a(1,j) * z(j-1)
    end do
    b(n) = z(n) * dr(n)

  else
    z(np) = b(np)
    z(np+1) = b(np+1) - a(m-1, np+1) * z(np)

    do j = np+2, n
      if (j > np-1+m) then
        i1 = 1
      else
        i1 = np-1+m - j + 1
      end if
      sum = dcmplx(0.0d0)
      do k = i1, m-1
        sum = sum + a(k, j) * z(j-m+k)
      end do
      z(j) = b(j) - sum
    end do

    do j = n-1, n
      z(j) = z(j) * dr(j)
    end do

    b(n) = z(n)
    b(n-1) = z(n-1) - a(m-1, n) * z(n)

  end if

  return
end subroutine
