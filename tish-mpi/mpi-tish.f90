!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ************** tish ****************
!  Computation of SH synthetic seismograms
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

program tish
  use FileIO
  implicit none

  !----------------------------<<constants>>----------------------------
  real(8), parameter :: pi = 3.1415926535897932d0
  integer, parameter :: maxNGrid = 88300  ! Maximum number of grid points.
  integer, parameter :: maxNZone = 15  ! Maximum number of zones.
  integer, parameter :: maxNReceiver = 1500  ! Maximum number of receivers.
  integer, parameter :: maxL = 80000  ! Maximum of angular order to loop for.
  real(8), parameter :: shallowDepth = 100.d0  ! Threshold to consider evanescent regime for shallow events [km].
  integer, parameter :: spcFormat = 1  ! Format of output spc file (0:binary, 1:ascii).
  integer, parameter :: ilog = 0

  !----------------------------<<variables>>----------------------------
  ! Variables for the structure
  integer :: nZone  ! Number of zones.
  real(8) :: rmin, rmax  ! Minimum and maximum radii of region that will be handled [km].
  real(8) :: rminOfZone(maxNZone), rmaxOfZone(maxNZone)  ! Minimum and maximum radii of each zone [km].
  real(8) :: rhoPolynomials(4, maxNZone), vsvPolynomials(4, maxNZone), vshPolynomials(4, maxNZone)
  !:::::::::::::::::::: Polynomial functions (coefficients of cubic function) of rho [g/cm^3], vsv, and vsh [km/s] in each zone.
  real(8) :: qmuOfZone(maxNZone)  ! Qmu of each zone.

  ! Variables for the source
  real(8) :: r0, eqlat, eqlon, mt(3, 3)  ! Depth [km], coordinates [deg], and moment tensor [10^25 dyn cm] of source.
  real(8) :: ecL0  ! Elastic modulus L at source position [10^10 dyn/cm^2 = GPa].

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
  integer :: ltmp(2), lmax, imaxFixed

  ! Variables for grid spacing and cut-off
  real(8) :: kzAtZone(maxNZone)  ! Vertical wavenumber k_z at each zone [1/km]. (See section 3.2 of Kawai et al. 2006.)
  real(8) :: re  ! Desired relative error due to vertical gridding. (See eqs. 6.1-6.3 of Geller & Takeuchi 1995.)
  real(8) :: ratc  ! Threshold amplitude ratio for vertical grid cut-off.
  real(8) :: ratl  ! Threshold amplitude ratio for angular order cut-off.
  real(8) :: amplitudeAtGrid(maxNGrid)  ! Estimate of the amplitude at each grid point [km], used for vertical grid cut-off.
  integer :: llog

  ! Variables for the vertical grid
  integer :: nGrid  ! Total number of grid points.
  real(8) :: gridRadii(maxNGrid)  ! Radii of each grid point [km].
  integer :: nLayerInZone(maxNZone)  ! Number of layers in each zone.
  integer :: oGridOfZone(maxNZone)  ! Index of the first grid point in each zone.
  real(8) :: gridRadiiForSource(3)  ! Radii to use for source-related computations [km].
  integer :: iZoneOfSource  ! Which zone the source is in.
  integer :: iLayerOfSource  ! Index of layer that the source is in.

  ! Variables for the values
  integer :: nValue  ! Total number of values.
  real(8) :: valuedRadii(maxNGrid + maxNZone - 1)  ! Radii corresponding to each variable value [km].
  integer :: oValueOfZone(maxNZone)  ! Index of the first value in each zone.
  real(8) :: rhoValues(maxNGrid + maxNZone - 1)  ! Rho at each grid point (with 2 values at boundaries) [g/cm^3].
  real(8) :: ecLValues(maxNGrid + maxNZone - 1)  ! L at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: ecNValues(maxNGrid + maxNZone - 1)  ! N at each grid point (with 2 values at boundaries) [GPa].
  real(8) :: rhoValuesForSource(3), ecLValuesForSource(3), ecNValuesForSource(3)  ! Rho, L, and N at each source-related grid.
  complex(8) :: coefQmu(maxNZone)  ! Coefficients to multiply to elastic moduli for anelastic attenuation at each zone.

  ! Variables for the trial function
  real(8) :: plm(3, 0:3, maxNReceiver)  ! Values of the associated Legendre polynomials at each receiver and m, stored for 3 l's.
  !::::::::::::::::::::::::::::::::::::::: Arguments: previous l's (1 before : 3 before), m (0:3).
  complex(8) :: harmonicsValues(3, -2:2, maxNReceiver)  ! Values of vector harmonics terms at each receiver, computed for each l.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Arguments: term (1:3), m (-2:2), iReceiver.

  ! Variables for the matrix elements
  real(8) :: t(4 * maxNGrid - 4)
  real(8) :: h1(4 * maxNGrid - 4), h2sum(4 * maxNGrid - 4), h3(4 * maxNGrid - 4), h4(4 * maxNGrid - 4)
  real(8) :: gt(8), gh1(8), gh2(8), gh3(8), gh4(8)
  integer :: oPairOfZone(maxNZone)  ! Index of the first (iLayer, k', k)-pair in each zone.
  integer :: oPairOfSource  ! Index of the first (iLayer, k', k)-pair for the layer with the source.
  complex(8) :: a0(2, maxNGrid), a2(2, maxNGrid)
  complex(8) :: a(2, maxNGrid)
  complex(8) :: aaParts(4), aSourceParts(8), aSource(2, 3)
  complex(8) :: g_or_c(maxNGrid)  ! This holds either vector g [10^15 N] or c [km], depending on where in the code it is. CAUTION!!
  complex(8) :: u(3, maxNReceiver)  ! Displacement velocity - the unit is [km] in the frequency domain,
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: but when converted to the time domain, the unit becomes [km/s].

  ! Variables for the output file
  character(len=80) :: output(maxNReceiver)

  ! Other variables
  real(8) :: work(4 * maxNGrid - 4)  ! Working array for matrix computations.
  complex(8) :: cwork(4 * maxNGrid - 4)  ! Working array for matrix computations.
  complex(8) :: dr(maxNGrid), z(maxNGrid), gdr(3)  ! Working arrays used when solving linear equations.

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

  ! Variables for MPI   !!!diff from non-mpi
  include 'mpif.h'
  integer :: petot, my_rank, ierr
  integer :: unitNum
  integer, allocatable, dimension (:) :: mpimin, mpimax


  ! ########################## Setup ##########################

  ! ************************** MPI **************************   !!!diff from non-mpi
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, petot, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  allocate(mpimin(petot), mpimax(petot))
  unitNum = 10 + my_rank


  ! ************************** Inputting parameters **************************
  ! --- read parameters ---
  if (my_rank == 0) then   !!!diff from non-mpi
    call readInput(maxNZone, maxNReceiver, tlen, np, re, ratc, ratl, omegaI, imin, imax, &
      nZone, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, qmuOfZone, &
      r0, eqlat, eqlon, mt, nReceiver, lat, lon, theta, phi, output)
  end if
  ! share input values
  call MPI_BCAST(re, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)   !!!diff from non-mpi
  call MPI_BCAST(ratc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ratl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tlen, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(np, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(omegaI, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(imin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(imax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nZone, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(rminOfZone, maxNZone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(rmaxOfZone, maxNZone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(rhoPolynomials, 4 * maxNZone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(vsvPolynomials, 4 * maxNZone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(vshPolynomials, 4 * maxNZone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(qmuOfZone, maxNZone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(r0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eqlat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eqlon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(mt, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nReceiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(lat, maxNReceiver, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(lon, maxNReceiver, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(theta, maxNReceiver, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(phi, maxNReceiver, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(output, 80 * maxNReceiver, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  ! --- computing the required parameters ---
  rmin = rminOfZone(1)
  rmax = rmaxOfZone(nZone)
  if (r0 < rmin .or. r0 > rmax) stop 'The source position is improper.'

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
  lmax = 0


  ! ************************** Files handling **************************
  if (my_rank == 0) then   !!!diff from non-mpi
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
      open(unit = 11, file = 'llog_SH.log', status = 'unknown')
      write(11, *) 'iFreq, llog, nGrid-1'
      close(11)
    end if
  end if


  ! ########################## Option for shallow events ##########################
  ! Here, we find the maximum angular order needed for our frequency range. (See fig. 7 of Kawai et al. 2006.)
  if ((rmax - r0) < shallowDepth) then

    ! Set a large value so that we can compute using fine grids for this process.
    imaxFixed = int(tlen * 2.d0)  !!! difference from main section

    ! ******************* Computing the matrix elements *******************  !!! difference from main section
    call computeMatrixElements(maxNGrid, shallowDepth, tlen, re, imin, imaxFixed, 0, r0, &
      nZone, rmin, rmax, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, &
      kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oPairOfZone, &
      iZoneOfSource, iLayerOfSource, oPairOfSource, gridRadiiForSource, &
      nValue, valuedRadii, rhoValues, ecLValues, ecNValues, rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, ecL0, &
      t, h1, h2sum, h3, h4, gt, gh1, gh2, gh3, gh4, work)

    ! ******************** Computing the expansion coefficients *********************
    ! Find the maximum angular order needed for the lowest and highest frequencies. (See fig. 7 of Kawai et al. 2006.)
    do iCount = 1, 2  ! omega-loop
      if (iCount == 1) then  !!! difference from main section
        iFreq = imin
      else
        iFreq = imax
      end if
      omega = 2.d0 * pi * dble(iFreq) / tlen

      call omegaLoopForShallowEvent(omega, omegaI, maxL, nZone, rmaxOfZone, vsvPolynomials, qmuOfZone, &
        r0, mt, ecL0, ratc, ratl, amplitudeAtGrid, &
        nGrid, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, coefQmu, &
        t, h1, h2sum, h3, h4, gt, gh1, gh2, gh3, gh4, oPairOfZone, oPairOfSource, aaParts, aSourceParts, aSource, &
        a0, a2, a, g_or_c, cwork, dr, z, gdr, eps, ltmp(iCount))

    end do  ! omega-loop

    lmax = max(ltmp(1), ltmp(2))  !!! difference from main section
  end if  ! option for shallow events


  ! ########################## Main computation ##########################

  ! ******************* Computing the matrix elements *******************
  call computeMatrixElements(maxNGrid, shallowDepth, tlen, re, imin, imax, lmax, r0, &
    nZone, rmin, rmax, rminOfZone, rmaxOfZone, rhoPolynomials, vsvPolynomials, vshPolynomials, &
    kzAtZone, nGrid, nLayerInZone, gridRadii, oGridOfZone, oValueOfZone, oPairOfZone, &
    iZoneOfSource, iLayerOfSource, oPairOfSource, gridRadiiForSource, &
    nValue, valuedRadii, rhoValues, ecLValues, ecNValues, rhoValuesForSource, ecLValuesForSource, ecNValuesForSource, ecL0, &
    t, h1, h2sum, h3, h4, gt, gh1, gh2, gh3, gh4, work)

  ! ******************** Computing the displacement *********************
  outputCounter = 1  !!! difference from shallow-source section

  !call trapezoidSplit(imin, imax, petot, mpimin, mpimax)   !!!diff from non-mpi
  !do iFreq = mpimin(my_rank + 1), mpimax(my_rank + 1)  ! omega-loop   !!!diff from non-mpi

  do iFreq = imin + my_rank, imax, petot  ! omega-loop   !!!diff from non-mpi
    omega = 2.d0 * pi * dble(iFreq) / tlen

    call omegaLoop(omega, omegaI, maxL, nZone, rmaxOfZone, vsvPolynomials, qmuOfZone, &
      r0, mt, ecL0, nReceiver, theta, phi, ratc, ratl, amplitudeAtGrid, &
      nGrid, nLayerInZone, oGridOfZone, iZoneOfSource, iLayerOfSource, coefQmu, plm, harmonicsValues, &
      t, h1, h2sum, h3, h4, gt, gh1, gh2, gh3, gh4, oPairOfZone, oPairOfSource, aaParts, aSourceParts, aSource, &
      a0, a2, a, g_or_c, u, cwork, dr, z, gdr, eps, llog)

    ! Store results.  !!! difference from shallow-source section
    outputi(outputCounter) = iFreq
    do ir = 1, nReceiver
      outputu(:, ir, outputCounter) = u(:, ir)
    end do

    ! ************************** Files Handling **************************
    ! Write to file when the output interval is reached, or when this is the last omega.   !!!diff from non-mpi
    if (outputCounter >= outputInterval .or. iFreq + petot > imax) then

      do ir = 1, nReceiver
        call openSPCFileMPI(output(ir), unitNum, spcFormat)
        do iOut = 1, outputCounter
          call writeSPCFile(unitNum, spcFormat, outputi(iOut), dble(outputu(1, ir, iOut)), imag(outputu(1, ir, iOut)))
          call writeSPCFile(unitNum, spcFormat, dble(outputu(2, ir, iOut)), imag(outputu(2, ir, iOut)))
          call writeSPCFile(unitNum, spcFormat, dble(outputu(3, ir, iOut)), imag(outputu(3, ir, iOut)))
        end do
        call closeSPCFile(unitNum)
      end do

      outputCounter = 0
    end if

    if (ilog == 1) then
      open(unit=11, file='llog_SH.log', position='append', status='old')
      write(11,*) iFreq, llog, nGrid-1
      close(11)
    end if

    outputCounter = outputCounter + 1

  end do  ! omega-loop


  ! ########################## Finishing ##########################

  ! Deallocate arrays.
  deallocate(outputi)
  deallocate(outputu)
  deallocate(mpimin)   !!!diff from non-mpi
  deallocate(mpimax)

  write(*,*) my_rank, "Ivalice looks to the horizon"   !!!diff from non-mpi
  call MPI_FINALIZE(ierr)

  stop
end program tish
