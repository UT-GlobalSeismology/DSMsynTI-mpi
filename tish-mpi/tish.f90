!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ************** tish ****************
!  Computation of SH synthetic seismograms
!  in transversely isotropic media for anisotropic PREM
!  using modified DSM operators & modified source representation.
!  Synthetics for shallow events can be computed.
!
!  Main historical authors: K.Kawai, N.Takeuchi, R.J.Geller
!  (C) 2002.10  University of Tokyo
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program tish
  implicit none

  !----------------------------<<constants>>----------------------------
  real(8) :: pi, lmaxdivf, shallowdepth
  integer :: maxnlayer, maxnzone, maxnr, maxlmax, ilog, spcform
  parameter (pi = 3.1415926535897932d0)
  parameter (maxnlayer = 88300)  ! Maximum number of layers.
  parameter (maxnzone = 15)  ! Maximum number of zones.
  parameter (maxnr = 1500)
  parameter (maxlmax = 80000)
  parameter (ilog = 0)
  parameter (lmaxdivf = 2.d4)
  parameter (shallowdepth = 100.d0)
  parameter (spcform = 0)  ! 0:binary, 1:ascii

  !----------------------------<<variables>>----------------------------
  ! Variables for the trial function
  integer :: nlayer, nlayerinzone(maxnzone)  ! Number of layers (total, and in each zone).
  integer :: l, m
  real(8) :: radii(maxnlayer + maxnzone + 1)  ! Radii of each grid point.
  real(8) :: gra(3)

  ! Variables for the structure
  integer :: nzone  ! Number of zones.
  integer :: ndc, vnp
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled.
  real(8) :: vrmin(maxnzone), vrmax(maxnzone)  ! Minimum and maximum radii of each zone.
  real(8) :: rrho(4, maxnzone), vsv(4, maxnzone), vsh(4, maxnzone)  ! Rho, vsv, and vsh of each zone (coefficients of cubic function).
  real(8) :: qmu(maxnzone)
  real(8) :: vra(maxnlayer + 2 * maxnzone + 1)
  real(8) :: rho(maxnlayer + 2 * maxnzone + 1)
  real(8) :: ecL(maxnlayer + 2 * maxnzone + 1)
  real(8) :: ecN(maxnlayer + 2 * maxnzone + 1)
  real(8) :: gvra(3), grho(3), gecL(3), gecN(3)
  complex(16) :: coef(maxnzone)

  ! Variables for the periodic range
  integer :: np  ! Number of points in frequency domain.
  integer :: imin, imax  ! Index of minimum and maximum frequency.
  real(8) :: tlen, omega, omegai

  ! Variables for the source
  integer :: spn, ns
  real(8) :: r0, mt(3, 3), spo, mu0, eqlat, eqlon

  ! Variables for the stations
  integer :: nr, ir
  real(8) :: theta(maxnr), phi(maxnr)
  real(8) :: lat(maxnr), lon(maxnr)

  ! Variables for the matrix elements
  complex(16) :: a0(2, maxnlayer+1), a2(2, maxnlayer+1)
  complex(16) :: a(2, maxnlayer+1)
  real(8) :: t(4 * maxnlayer)
  real(8) :: h1(4 * maxnlayer), h2(4 * maxnlayer)
  real(8) :: h3(4 * maxnlayer), h4(4 * maxnlayer)
  real(8) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  complex(16) :: aa(4), ga(8), ga2(2, 3), gdr(3)
  complex(16) :: g(maxnlayer + 1)

  ! Variables for the output file
  character(len=80) :: output(maxnr)

  ! Variables for grid spacing
  real(8) :: tmpr(maxnlayer + 1)
  real(8) :: kz(maxnzone)  ! Vertical wavenumber k_z. (See section 3.2 of Kawai et al. 2006.)
  real(8) :: re, ratc, ratl, maxamp
  integer :: kc, lsuf, ismall, llog

  ! Variables for the stack points
  integer :: isp(maxnzone), jsp(maxnzone), ins

  ! Other variables
  integer :: i, j, ii, jj, nn, ier
  real(8) :: work(4 * maxnlayer), lsq
  complex(16) :: dr(maxnlayer + 1), z(maxnlayer + 1)
  complex(16) :: cwork(4 * maxnlayer)
  integer :: ltmp(2), iimax

  ! Constants with data statements
  integer :: lda = 2
  real(8) :: eps = -1.d0

  ! Efficiency improvement variables
  integer :: mpios
  integer :: outputmemory = 10  ! MB
  integer :: outputinterval
  real(8) :: memoryperomega ! MB
  integer :: outputindex, mpii
  integer, allocatable :: outputi(:)
  complex(16), allocatable :: outputu(:,:,:)
!     When the values to be output use memory over outputmemory MB,
!     they are written in output files. The interval is outputinterval.
!     memoryperomega is the quantity of memory used for one omega step.
  character(len=2) :: char_rank
  real(8) :: ark, angel


  ! ************************** Inputting parameters **************************
  ! --- read parameters ---
  call pinput2(maxnzone, maxnr, re, ratc, ratl, tlen, np, omegai, &
    imin, imax, nzone, vrmin, vrmax, rrho, vsv, vsh, qmu, &
    r0, eqlat, eqlon, mt, nr, theta, phi, lat, lon, output)

  memoryperomega = 3 * 16 * nr * 0.000001
  outputinterval = outputmemory / memoryperomega
  allocate(outputi(outputinterval))
  allocate(outputu(3, nr, outputinterval))

  ! --- computing the required parameters ---
  rmin = vrmin(1)
  rmax = vrmax(nzone)
  ndc = nzone - 1

  do ir = 1, nr
    theta(ir) = theta(ir) / 180.0d0 * pi ! Convert theta from degrees to radians
    phi(ir) = phi(ir) / 180.0d0 * pi     ! Convert phi from degrees to radians
  end do

  if (r0 < rmin .or. r0 > rmax) then
    stop 'Location of the source is improper.'
  end if

  ! ************************** Files handling **************************
  if (spcform == 0) then
    do ir = 1, nr
      open(unit = 11, file = output(ir), status = 'unknown', &
        form = 'unformatted', access = 'stream', convert = 'big_endian')
      write(11) tlen
      write(11) np, 1, 3
      write(11) omegai, lat(ir), lon(ir)
      write(11) eqlat, eqlon, r0
      close(11)
    end do
  else if (spcform == 1) then
    do ir = 1, nr
      open(unit = 11, file = output(ir), status = 'unknown')
      write(11, *) tlen
      write(11, *) np, 1, 3
      write(11, *) omegai, lat(ir), lon(ir)
      write(11, *) eqlat, eqlon, r0
      close(11)
    end do
  else
    write(*, *) 'WARNING:(tish.f)  set spcform 0 or 1'
  end if

  if (ilog == 1) then
    open(unit = 11, file = 'llog.log', status = 'unknown')
    close(11)
  end if



  ! ************************** Option for shallow events **************************
  iimax = imax
  if ((rmax-r0) < shallowdepth) then
    ! computing of the number and the location of grid points
    iimax = int(tlen * 2.d0)
    call computeKz(nzone, vrmin, vrmax, vsv, rmin, rmax, iimax, 1, tlen, kz)
    call computeGridRadii(maxnlayer, maxnzone, nzone, kz, vrmin, vrmax, rmin, rmax, re, nlayer, nlayerinzone, radii)

    ! checking the parameter
    if (nlayer > maxnlayer) stop 'The number of grid points is too large.'

    ! computing the stack points
    call calsp(ndc, nlayerinzone, isp, jsp)

    ! computing the source location
    call calspo(ndc, vrmax, nlayer, r0, rmin, rmax, radii, isp, spo, spn)

    ! computing grids for source computations
    call calgra(isp, radii, r0, spn, spo, gra)

    ! ******************* Computing the matrix elements *******************
    ! computing the structure grid points
    call calstg(nzone, rrho, vsv, vsh, nlayer, nlayerinzone, radii, rmax, vnp, vra, rho, ecL, ecN)
    call calgstg(spn, rrho, vsv, vsh, gra, gvra, rmax, grho, gecL, gecN, r0, mu0)

    do i = 1, ndc + 1
      call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, rho, 2, 0, 0, radii(isp(i)), t(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecL, 2, 1, 1, radii(isp(i)), h1(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecL, 1, 1, 0, radii(isp(i)), h2(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecL, 0, 0, 0, radii(isp(i)), h3(jsp(i)), work(jsp(i)))
      call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecN, 0, 0, 0, radii(isp(i)), h4(jsp(i)), work(jsp(i)))
      call computeLumpedT(nlayerinzone(i), vnp, vra, rho, radii(isp(i)), work(jsp(i)))
      call computeAverage(nlayerinzone(i), t(jsp(i)), work(jsp(i)), t(jsp(i)))
      call computeLumpedH(nlayerinzone(i), vnp, vra, ecL, radii(isp(i)), work(jsp(i)))
      call computeAverage(nlayerinzone(i), h3(jsp(i)), work(jsp(i)), h3(jsp(i)))
      call computeLumpedH(nlayerinzone(i), vnp, vra, ecN, radii(isp(i)), work(jsp(i)))
      call computeAverage(nlayerinzone(i), h4(jsp(i)), work(jsp(i)), h4(jsp(i)))
    end do

    call computeIntermediateIntegral(2, 3, gvra, grho, 2, 0, 0, gra, gt, work)
    call computeIntermediateIntegral(2, 3, gvra, gecL, 2, 1, 1, gra, gh1, work)
    call computeIntermediateIntegral(2, 3, gvra, gecL, 1, 1, 0, gra, gh2, work)
    call computeIntermediateIntegral(2, 3, gvra, gecL, 0, 0, 0, gra, gh3, work)
    call computeIntermediateIntegral(2, 3, gvra, gecN, 0, 0, 0, gra, gh4, work)
    call computeLumpedT(2, 3, gvra, grho, gra, work)
    call computeAverage(2, gh3, work, gh3)
    call computeLumpedH(2, 3, gvra, gecL, gra, work)
    call computeAverage(2, gh3, work, gh3)
    call computeLumpedH(2, 3, gvra, gecN, gra, work)
    call computeAverage(2, gh4, work, gh4)

    nn = nlayer + 1
    ns = isp(spn) + dint(spo)
    ins = 4 * ns - 3


    llog = 0
    do ii = 1,2  ! omega-loop
      if (ii == 1) then
        if (imin == 0) then
          i = 1
        else
          i = imin
        endif
      endif
      if (ii == 2) i = imax
      omega = 2.d0 * pi * dble(i) / tlen
      call callsuf(omega, nzone, vrmax, vsv, lsuf)
      call calcoef(nzone, omega, qmu, coef)

      call cmatinit(lda, nn, a0)
      call cmatinit(lda, nn, a2)
      do j = 1, ndc + 1
        call cala0(nlayerinzone(j), omega, omegai, t(jsp(j)), h1(jsp(j)), h2(jsp(j)), h3(jsp(j)), h4(jsp(j)), coef(j), cwork(jsp(j)))
        call overlap(nlayerinzone(j), cwork(jsp(j)), a0(1, isp(j)))
        call cala2(nlayerinzone(j), h4(jsp(j)), coef(j), cwork(jsp(j)))
        call overlap(nlayerinzone(j), cwork(jsp(j)), a2(1, isp(j)))
      enddo

      kc = 1
      ismall = 0
      maxamp = -1.d0
      ltmp(ii) = maxlmax
      do l = 0, maxlmax  ! l-loop
        if (ismall > 20) then
          if (ltmp(ii) > l) ltmp(ii) = l
          exit
        endif

        do jj = 1, maxnlayer + 1  ! initialize
          tmpr(jj) = 0.d0
        enddo
        lsq = dsqrt(dble(l) * dble(l + 1))
        ! computing the coefficient matrix elements
        ! --- Initializing the matrix elements
        call cmatinit(lda, nn, a)
        call cmatinit(lda, 3, ga2)
        call cala(nn, l, lda, a0, a2, a)
        call calga(1, omega, omegai, l, t(ins), h1(ins), h2(ins), h3(ins), h4(ins), coef(spn), aa)
        call calga(2, omega, omegai, l, gt, gh1, gh2, gh3, gh4, coef(spn), ga)
        call overlap(2, ga, ga2)

        do m = -2, 2  ! m-loop
          if ((m /= 0) .and. (iabs(m) <= iabs(l))) then
            call cvecinit(nn, g)
            call calg2(l, m, spo, r0, mt, mu0, coef(spn), ga, aa, ga2, gdr, g(isp(spn)))
            if (mod(l, 100) == 0) then
              if ((m == -2) .or. (m == -l)) then
                call dclisb0(a, nn, 1, lda, g, eps, dr, z, ier)
              else
                call dcsbsub0(a, nn, 1, lda, g, eps, dr, z, ier)
              endif
              do jj = 1, nn  ! sum up c of the same l
                tmpr(jj) = tmpr(jj) + abs(g(jj))
              enddo
            else
              if ((m == -2) .or. (m == -l)) then
                call dclisb(a(1, kc), nn - kc + 1, 1, lda, ns - kc + 1, g(kc), eps, dr, z, ier)
              else
                call dcsbsub(a(1, kc), nn - kc + 1, 1, lda, ns - kc + 1, g(kc), eps, dr, z, ier)
              endif
            endif

            if (mod(l, 100) == 0) then
              call calcutd(nzone, nlayerinzone, tmpr, ratc, nn, radii, kc)
            endif

            call calamp(g(nn), l, lsuf, maxamp, ismall, ratl)
          endif
        enddo  ! m-loop
      enddo  ! l-loop
    enddo  ! omega-loop

    iimax = dble(max(ltmp(1),ltmp(2))) * tlen / lmaxdivf
  endif  ! option for shallow events


! ******************* Computing parameters *******************
! computing of the number and the location of grid points
  call computeKz(nzone, vrmin, vrmax, vsv, rmin, rmax, iimax, 1, tlen, kz)
  call computeGridRadii(maxnlayer, maxnzone, nzone, kz, vrmin, vrmax, rmin, rmax, re, nlayer, nlayerinzone, radii)

! checking the parameter
  if (nlayer > maxnlayer) stop 'The number of grid points is too large.'

! computing the stack points
  call calsp(ndc, nlayerinzone, isp, jsp)

! computing the source location
  call calspo(ndc, vrmax, nlayer, r0, rmin, rmax, radii, isp, spo, spn)

! computing grids for source computations
  call calgra(isp, radii, r0, spn, spo, gra)

! ******************* Computing the matrix elements *******************
! computing the structure grid points
  call calstg(nzone, rrho, vsv, vsh, nlayer, nlayerinzone, radii, rmax, vnp, vra, rho, ecL, ecN)
  call calgstg(spn, rrho, vsv, vsh, gra, gvra, rmax, grho, gecL, gecN, r0, mu0)

  do i = 1, ndc + 1
    call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, rho, 2, 0, 0, radii(isp(i)), t(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecL, 2, 1, 1, radii(isp(i)), h1(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecL, 1, 1, 0, radii(isp(i)), h2(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecL, 0, 0, 0, radii(isp(i)), h3(jsp(i)), work(jsp(i)))
    call computeIntermediateIntegral(nlayerinzone(i), vnp, vra, ecN, 0, 0, 0, radii(isp(i)), h4(jsp(i)), work(jsp(i)))
    call computeLumpedT(nlayerinzone(i), vnp, vra, rho, radii(isp(i)), work(jsp(i)))
    call computeAverage(nlayerinzone(i), t(jsp(i)), work(jsp(i)), t(jsp(i)))
    call computeLumpedH(nlayerinzone(i), vnp, vra, ecL, radii(isp(i)), work(jsp(i)))
    call computeAverage(nlayerinzone(i), h3(jsp(i)), work(jsp(i)), h3(jsp(i)))
    call computeLumpedH(nlayerinzone(i), vnp, vra, ecN, radii(isp(i)), work(jsp(i)))
    call computeAverage(nlayerinzone(i), h4(jsp(i)), work(jsp(i)), h4(jsp(i)))
  enddo

  call computeIntermediateIntegral(2, 3, gvra, grho, 2, 0, 0, gra, gt, work)
  call computeIntermediateIntegral(2, 3, gvra, gecL, 2, 1, 1, gra, gh1, work)
  call computeIntermediateIntegral(2, 3, gvra, gecL, 1, 1, 0, gra, gh2, work)
  call computeIntermediateIntegral(2, 3, gvra, gecL, 0, 0, 0, gra, gh3, work)
  call computeIntermediateIntegral(2, 3, gvra, gecN, 0, 0, 0, gra, gh4, work)
  call computeLumpedT(2, 3, gvra, grho, gra, work)
  call computeAverage(2, gt, work, gt)
  call computeLumpedH(2, 3, gvra, gecL, gra, work)
  call computeAverage(2, gh3, work, gh3)
  call computeLumpedH(2, 3, gvra, gecN, gra, work)
  call computeAverage(2, gh4, work, gh4)






end program tish


