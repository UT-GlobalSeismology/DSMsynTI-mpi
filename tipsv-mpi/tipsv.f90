!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ************** tipsv ****************
!  Computation of PSV synthetic seismograms
!  in transversely isotropic media for anisotropic PREM
!  using modified DSM operators & modified source representation.
!  Synthetics for shallow events can be computed.
!
!  Main historical authors: K.Kawai, N.Takeuchi, R.J.Geller
!  (C) 2002 - 2024  University of Tokyo
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program tipsv
  use FileIO
  implicit none

  !----------------------------<<constants>>----------------------------
  real(8), parameter :: pi = 3.1415926535897932d0
  integer, parameter :: maxNGrid = 88300  ! Maximum number of grid points.
  integer, parameter :: maxNGridSolid = 48840  ! Maximum number of grid points in solid region.
  integer, parameter :: maxNGridFluid = 32040  ! Maximum number of grid points in fluid region.
  integer, parameter :: maxNZone = 15  ! Maximum number of zones.
  integer, parameter :: maxNReceiver = 1500  ! Maximum number of receivers.
  integer, parameter :: maxL = 80000  ! Maximum of angular order to loop for.
  real(8), parameter :: lmaxdivf = 2.d4  !    !!TODO where did this value come from?
  real(8), parameter :: shallowdepth = 100.d0  ! Threshold to consider evanescent regime for shallow events [km].
  integer, parameter :: spcFormat = 1  ! Format of output spc file (0:binary, 1:ascii).
  integer, parameter :: ilog = 1

  !----------------------------<<variables>>----------------------------
  ! Variables for the structure
  integer :: nZone  ! Number of zones.
  integer :: nZoneSolid, nZoneFluid  ! Number of solid and fluid zones.
  integer :: phaseOfZone(maxNZone)  ! Phase of each zone (1: solid, 2: fluid).
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled [km].
  real(8) :: rminOfZone(maxNZone), rmaxOfZone(maxNZone)  ! Minimum and maximum radii of each zone [km].
  real(8) :: rhoPolynomials(4, maxNZone), vpvPolynomials(4, maxNZone), vphPolynomials(4, maxNZone)
  real(8) :: vsvPolynomials(4, maxNZone), vshPolynomials(4, maxNZone), etaPolynomials(4, maxNZone)
  !::::: Polynomial functions (coefficients of cubic function) of rho [g/cm^3], vpv, vph, vsv, vsh [km/s], and eta in each zone.
  real(8) :: qmuOfZone(maxNZone), qkappaOfZone(maxNZone)  ! Qmu and Qkappa of each zone.

  ! Variables for the source
  real(8) :: r0, eqlat, eqlon, mt(3, 3)  ! Depth [km], coordinates [deg], and moment tensor [10^25 dyn cm] of source.
  real(8) :: ecC0, ecF0, ecL0  ! Elastic moduli C, F, and L at source position [10^10 dyn/cm^2 = GPa].

  ! Variables for the receivers
  integer :: nReceiver  ! Number of receivers.
  real(8) :: lat(maxNReceiver), lon(maxNReceiver)  ! Coordinates [deg] of receivers.
  real(8) :: theta(maxNReceiver), phi(maxNReceiver)  ! Colatitude and longitude of receivers with event at north pole [rad].
  integer :: ir

  ! Variables for the periodic range
  real(8) :: tlen  ! Time length [s].
  integer :: np  ! Number of points in frequency domain.
  real(8) :: omega  ! Angular frequency (real part) [1/s].
  real(8) :: omegaI  ! Imaginary part of angular frequency for artificial damping [1/s].
  !:::::::::::::::::::::::::::::::::::::::: (See section 5.1 of Geller & Ohminato 1994.)
  integer :: imin, imax  ! Index of minimum and maximum frequency.
  integer :: iFreq, iCount
  integer :: ltmp(2), imaxFixed

  ! Variables for grid spacing and cut-off
  real(8) :: kzAtZone(maxNZone)  ! Vertical wavenumber k_z at each zone [1/km]. (See section 3.2 of Kawai et al. 2006.)
  real(8) :: re  ! Desired relative error due to vertical gridding. (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
  real(8) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8) :: amplitudeAtColumn(2 * maxNGridSolid + maxNGridFluid)
  !:::::::::::::::::::::::::::: Estimate of the amplitude at each column [km], used for vertical grid cut-off.
  integer :: llog

  ! Variables for the vertical grid
  integer :: nGrid  ! Total number of grid points.
  real(8) :: gridRadii(maxNGrid)  ! Radii of each grid point [km].
  integer :: nLayerInZone(maxNZone)  ! Number of layers in each zone.
  integer :: oGridOfZone(maxNZone)  ! Index of the first grid point in each zone.
  integer :: iZoneOfSource  ! Which zone the source is in.
  integer :: iLayerOfSource  ! Index of layer that the source is in.

  ! Variables for the values
  integer :: nValue  ! Total number of values.
  real(8) :: valuedRadii(maxNGrid + maxNZone - 1)  ! Radii corresponding to each variable value [km].
  integer :: oValueOfZone(maxNZone)  ! Index of the first value in each zone.
  integer :: oValueOfZoneSolid(maxNZone)  ! Index of the first value in each zone, when counting only solid zones.
  real(8) :: rhoValues(maxNGrid + maxNZone - 1)  ! Rho at each grid point (with 2 values at boundaries) [g/cm^3].
  real(8) :: kappaValues(maxNGrid + maxNZone - 1)  ! Kappa at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: ecKxValues(maxNGrid + maxNZone - 1)  ! Kx=A-4N/3 at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: ecKyValues(maxNGrid + maxNZone - 1)  ! Ky=F+2N/3 at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: ecKzValues(maxNGrid + maxNZone - 1)  ! Kz=(C+2F)/3 at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: ecLValues(maxNGrid + maxNZone - 1)  ! L at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: ecNValues(maxNGrid + maxNZone - 1)  ! N at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: rhoReciprocals(maxNGrid + maxNZone - 1)  ! 1/rho at each grid point (with 2 values at boundaries) [cm^3/g].
  real(8) :: kappaReciprocals(maxNGrid + maxNZone - 1)  ! 1/kappa at each grid point (with 2 values at boundaries) [1/GPa].
  complex(8) :: coefQmu(maxNZone), coefQkappa(maxNZone), coefQfluid(maxNZone)
  !::::::::::::::::::::::::::::::::::::::::: Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.

  ! Variables for the trial function
  real(8) :: plm(3, 0:3, maxNReceiver)  ! Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  !::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8) :: harmonicsValues(3, -2:2, maxNReceiver)  ! Values of vector harmonics terms at each receiver, computed for each l.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Arguments: term (1:3), m (-2:2), iReceiver.

  ! Variables for the matrix elements
  real(8) :: t(4 * maxNGridSolid - 4)
  real(8) :: h1x(4 * maxNGridSolid - 4), h2L(4 * maxNGridSolid - 4), h2N(4 * maxNGridSolid - 4)
  real(8) :: hUn3y(4 * maxNGridSolid - 4), hResid3y(4 * maxNGridSolid - 4), hModL3y(-2:1, maxNGridSolid)
  real(8) :: hUn4L(4 * maxNGridSolid - 4), hResid4L(4 * maxNGridSolid - 4), hModR4L(-1:2, maxNGridSolid)
  real(8) :: hUn4N(4 * maxNGridSolid - 4), hResid4N(4 * maxNGridSolid - 4), hModL4N(-2:1, maxNGridSolid)
  real(8) :: hUn5y(4 * maxNGridSolid - 4), hResid5y(4 * maxNGridSolid - 4), hModR5y(-1:2, maxNGridSolid)
  real(8) :: hUn6L(4 * maxNGridSolid - 4), hResid6L(4 * maxNGridSolid - 4), hModL6L(-2:1, maxNGridSolid)
  real(8) :: hUn6N(4 * maxNGridSolid - 4), hResid6N(4 * maxNGridSolid - 4), hModR6N(-1:2, maxNGridSolid)
  real(8) :: h7y(4 * maxNGridSolid - 4), h7z(4 * maxNGridSolid - 4), h8L(4 * maxNGridSolid - 4), h8N(4 * maxNGridSolid - 4)
  real(8) :: p1(4 * maxNGridFluid - 4), p2(4 * maxNGridFluid - 4), p3(4 * maxNGridFluid - 4)
  integer :: oPairOfZoneSolid(maxNZone), oPairOfZoneFluid(maxNZone)
  !::::::::::::::::::::: Index of the first (iLayer, k', k)-pair in each zone, counted separately for solid and fluid zones.
  complex(8) :: a0(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a1(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a2(4, 2 * maxNGridSolid + maxNGridFluid)
  complex(8) :: a(4, 2 * maxNGridSolid + maxNGridFluid), aSmall(2, maxNGridSolid + maxNGridFluid)
  !::::::::::::::::::::::::::::::::::::::::::::::: [10^12 kg/s^2] for solid, [m^5/N] for fluid, [10^6 m^2] for boundary elements.
  complex(8) :: g_or_c(2 * maxNGridSolid + maxNGridFluid), g_or_c_Small(maxNGridSolid + maxNGridFluid)
  !::::: These hold either vector g ([10^15 N] & [10^9 m^3]) or c ([km] & [GPa]), depending on where in the code it is. CAUTION!!
  complex(8) :: u(3, maxNReceiver)  ! Displacement velocity - the unit is [km] in the frequency domain,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: but when converted to the time domain, the unit becomes [km/s].
  integer :: oElementOfZone(maxNZone)  ! Index of the first (iLayer, k'-gamma', k-gamma)-pair in each zone.
  integer :: oColumnOfZone(maxNZone + 1)  ! Index of the first column in the band matrix for each zone.
  integer :: nColumn, nQuasiColumn  ! Total number of columns in the band matrix.

  !!TODO ???
  complex(8) :: anum(4, 4, 10), bnum(4, 4, 10)
  complex(8) :: ya(4), yb(4), yc(4), yd(4)
  integer :: oColumnOfSource  ! Index of the first column in the band matrix for the layer with source.

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Other variables
  real(8) :: work(4 * maxNGrid - 4)  ! Working array for matrix computations.
  complex(8) :: cwork(16 * maxNGridSolid - 16 + 4 * maxNGridFluid - 4)  ! Working array for matrix computations.
  complex(8) :: z(2 * maxNGridSolid + maxNGridFluid), w(2 * maxNGridSolid + maxNGridFluid)
  !::::::::::::::::::::::::::::::::::: Working arrays used when solving linear equations.

  ! Constants
  real(8) :: eps = -1.d0

  ! Efficiency improvement variables
  !     When the values to be output use memory over outputMemory MB,
  !     they are written in output files. The interval is outputInterval.
  real(8) :: outputMemory = 10  ! Approximate amount of memory to write at a time. [MB]
  real(8) :: memoryPerOmega  ! The amount of memory used for one omega step. [MB]
  integer :: outputInterval  ! Interval of omega that data should be written out.
  integer :: outputCounter  ! Counter to keep track of how many omegas have been computed after previous output.
  integer :: iOut
  integer, allocatable :: outputi(:)
  complex(8), allocatable :: outputu(:,:,:)


  ! ########################## Setup ##########################

  ! ************************** Inputting parameters **************************
  ! --- read parameters ---
  call readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegaI, imin, imax, nZone, rminOfZone, rmaxOfZone, &
    rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
    r0, eqlat, eqlon, mt, nReceiver, lat, lon, theta, phi, output)

  ! --- computing the required parameters ---
  rmin = rminOfZone(1)
  rmax = rmaxOfZone(nZone)
  if (r0 < rmin .or. r0 > rmax) stop 'The source position is improper.'

  call judgeSolidOrFluid(nZone, vsvPolynomials, phaseOfZone, nZoneSolid, nZoneFluid)

  ! Find the amount of memory that is written in 1 omega step.
  !  For each omega and receiver, 3 complex numbers (16 B each) are written. 1 B = 0.000001 MB.
  memoryPerOmega = 3 * 16 * nReceiver * 0.000001
  ! Find how many omegas can be written within outputMemory.
  outputInterval = int(outputMemory / memoryPerOmega)
  ! Allocate arrays to store output.
  allocate(outputi(outputInterval))
  allocate(outputu(3, nReceiver, outputInterval))

  if (imin == 0) imin = 1
  ! Decide which omega to use when deciding grid spacing. Usually, this is just the upper limit of omega range.
  imaxFixed = imax


  ! ************************** Files handling **************************
  do ir = 1, nReceiver
    call openSPCFile(output(ir), 11, spcFormat, 0)
    call writeSPCFile(11, spcFormat, tlen)
    call writeSPCFile(11, spcFormat, np, 1, 3)
    call writeSPCFile(11, spcFormat, omegaI, lat(ir), lon(ir))
    call writeSPCFile(11, spcFormat, eqlat, eqlon, r0)
    call writeSPCFile(11, spcFormat, 0, 0.d0, 0.d0)  !TODO erase when other software become compatible
    call writeSPCFile(11, spcFormat, 0.d0, 0.d0)  !TODO erase when other software become compatible
    call writeSPCFile(11, spcFormat, 0.d0, 0.d0)  !TODO erase when other software become compatible
    call closeSPCFile(11)
  end do

  if (ilog == 1) then
    open(unit = 11, file = 'llog_PSV.log', status = 'unknown')
    write(11, *) 'iFreq, llog, nGrid-1'
    close(11)
  end if


  ! ########################## Option for shallow events ##########################
  ! Here, we find the maximum angular order needed for our frequency range. (See fig. 7 of Kawai et al. 2006.)
  if ((rmax - r0) < shallowdepth) then

    write(*, *) 'Shallow event!'  !TODO erase

    ! Set a large value so that we can compute using fine grids for this process.
    imaxFixed = int(tlen * 2.d0)  !!! difference from main section

    ! ******************* Computing the matrix elements *******************
    call computeMatrixElements(maxNGrid, maxNGridSolid, maxNGridFluid, tlen, re, imaxFixed, r0, &
      nZone, rmin, rmax, rminOfZone, rmaxOfZone, phaseOfZone, &
      rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, &
      kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oValueOfZoneSolid, &
      oPairOfZoneSolid, oPairOfZoneFluid, oElementOfZone, oColumnOfZone, nColumn, &
      iZoneOfSource, iLayerOfSource, oColumnOfSource, &
      nValue, valuedRadii, rhoValues, kappaValues, ecKxValues, ecKyValues, ecKzValues, ecLValues, ecNValues, ecC0, ecF0, ecL0, &
      rhoReciprocals, kappaReciprocals, &
      t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
      hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, work)

    ! ******************** Computing the expansion coefficients *********************
    ! Find the maximum angular order needed for the lowest and highest frequencies. (See fig. 7 of Kawai et al. 2006.)
    do iCount = 1, 2  ! omega-loop
      if (iCount == 1) then  !!! difference from main section
        iFreq = imin
      else
        iFreq = imax
      end if
      omega = 2.d0 * pi * dble(iFreq) / tlen

      call omegaLoopForShallowEvent(omega, omegaI, maxL, maxNGridSolid, maxNGridFluid, &
        nZone, rmin, rmax, rmaxOfZone, phaseOfZone, &
        rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
        r0, mt, ecC0, ecF0, ecL0, ratc, ratl, amplitudeAtColumn, &
        nGrid, gridRadii, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, &
        oValueOfZoneSolid, coefQmu, coefQkappa, coefQfluid, &
        t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
        hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, &
        oPairOfZoneSolid, oPairOfZoneFluid, a0, a1, a2, a, aSmall, g_or_c, g_or_c_Small, &
        oElementOfZone, oColumnOfZone, nColumn, nQuasiColumn, anum, bnum, ya, yb, yc, yd, oColumnOfSource, cwork, z, w, eps, &
        ltmp(iCount))

    end do  ! omega-loop

    imaxFixed = max(imax, int(dble(max(ltmp(1), ltmp(2))) * tlen / lmaxdivf))  !!! difference from main section

    write(*, *) ltmp(1), ltmp(2), imax, imaxFixed   !TODO erase
    write(*, *) 'Ending shallow-event section.'  !TODO erase
  end if  ! option for shallow events


  ! ########################## Main computation ##########################

  ! ******************* Computing the matrix elements *******************
  call computeMatrixElements(maxNGrid, maxNGridSolid, maxNGridFluid, tlen, re, imaxFixed, r0, &
    nZone, rmin, rmax, rminOfZone, rmaxOfZone, phaseOfZone, &
    rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, &
    kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oValueOfZoneSolid, &
    oPairOfZoneSolid, oPairOfZoneFluid, oElementOfZone, oColumnOfZone, nColumn, &
    iZoneOfSource, iLayerOfSource, oColumnOfSource, &
    nValue, valuedRadii, rhoValues, kappaValues, ecKxValues, ecKyValues, ecKzValues, ecLValues, ecNValues, ecC0, ecF0, ecL0, &
    rhoReciprocals, kappaReciprocals, &
    t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
    hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, work)

  !******************** Computing the displacement *********************
  outputCounter = 1  !!! difference from shallow-source section

  do iFreq = imin, imax  ! omega-loop
    omega = 2.d0 * pi * dble(iFreq) / tlen

    call omegaLoop(omega, omegaI, maxL, maxNGridSolid, maxNGridFluid, &
      nZone, rmin, rmax, rmaxOfZone, phaseOfZone, &
      rhoPolynomials, vpvPolynomials, vphPolynomials, vsvPolynomials, vshPolynomials, etaPolynomials, qmuOfZone, qkappaOfZone, &
      r0, mt, ecC0, ecF0, ecL0, nReceiver, theta, phi, ratc, ratl, amplitudeAtColumn, &
      nGrid, gridRadii, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, &
      oValueOfZoneSolid, coefQmu, coefQkappa, coefQfluid, plm, harmonicsValues, &
      t, h1x, h2L, h2N, hUn3y, hResid3y, hModL3y, hUn4L, hResid4L, hModR4L, hUn4N, hResid4N, hModL4N, &
      hUn5y, hResid5y, hModR5y, hUn6L, hResid6L, hModL6L, hUn6N, hResid6N, hModR6N, h7y, h7z, h8L, h8N, p1, p2, p3, &
      oPairOfZoneSolid, oPairOfZoneFluid, a0, a1, a2, a, aSmall, g_or_c, g_or_c_Small, u, &
      oElementOfZone, oColumnOfZone, nColumn, nQuasiColumn, anum, bnum, ya, yb, yc, yd, oColumnOfSource, cwork, z, w, eps, llog)

    ! Store results.  !!! difference from shallow-source section
    outputi(outputCounter) = iFreq
    do ir = 1, nReceiver
      outputu(:, ir, outputCounter) = u(:, ir)
    end do

    ! ************************** Files Handling **************************
    ! Write to file when the output interval is reached, or when this is the last omega.
    if (outputCounter >= outputInterval .or. iFreq == imax) then
      write(*,*) "kakikomimasu"

      do ir = 1, nReceiver
        call openSPCFile(output(ir), 11, spcFormat, 1)
        do iOut = 1, outputCounter
          call writeSPCFile(11, spcFormat, outputi(iOut), dble(outputu(1, ir, iOut)), imag(outputu(1, ir, iOut)))
          call writeSPCFile(11, spcFormat, dble(outputu(2, ir, iOut)), imag(outputu(2, ir, iOut)))
          call writeSPCFile(11, spcFormat, dble(outputu(3, ir, iOut)), imag(outputu(3, ir, iOut)))
        end do
        call closeSPCFile(11)
      end do

      outputCounter = 0
    end if

    if (ilog == 1) then
      open(unit=11, file='llog_PSV.log', position='append', status='old')
      write(11,*) iFreq, llog, nGrid-1
      close(11)
    end if

    outputCounter = outputCounter + 1

  end do  ! omega-loop


  ! ########################## Finishing ##########################

  ! Deallocate arrays.
  deallocate(outputi)
  deallocate(outputu)

  write(*,*) "Ivalice looks to the horizon"

  stop
end program tipsv
