program glowbasicssj

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.98 beta, 1/2017

! Adapted from glowbasic by Liam Kilcommons 3/2017

! Basic single-processor driver for the GLOW model.
! Uses MSIS/IRI for input.
! Runs GLOW for designated inputs once, or multiple times.
! MPI and netCDF libraries not required.

! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.

! Other definitions:
! f107p   Solar 10.7 cm flux for previous day
! ap      Ap index of geomagnetic activity
! z       altitude array, km

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! lmax    number of wavelength intervals for solar flux
! nmaj    number of major species
! nst     number of states produced by photoionization/dissociation
! nei     number of states produced by electron impact
! nex     number of ionized/excited species
! nw      number of airglow emission wavelengths
! nc      number of component production terms for each emission
  
  use cglow,only: pedcond,hallcond,zz,jmax

  implicit none

  integer,parameter :: kmax=1 ! Maximum number of datapoints to be read 
  integer,parameter :: nssjch=19 ! Number of SSJ channels (should not change) 
  
  integer :: kread

  real :: idate(kmax),ut(kmax),glat(kmax),glong(kmax),f107a(kmax),f107(kmax),f107p(kmax),ap(kmax),ef(kmax),ec(kmax)
  real :: phissj(kmax,nssjch)
  integer :: i,iostatus,j,k
  !real,allocatable :: pedcond(:),hallcond(:)

!  do instance=1,10000
!
! Read input data
  do k=1,kmax
!    write(6,"('Enter date, UT, lat, lon, F107a, F107, F107p, Ap, Ef, Ec, PhiSSJ( 19 values)')")
    read(5,*,iostat=iostatus) idate(k),ut(k),glat(k),glong(k),f107a(k),f107(k),f107p(k),ap(k),ef(k),ec(k),(phissj(k,i),i=1,19)
    if (iostatus /= 0) stop
    kread=k
  enddo
!
! Call subroutine
!
  do k=1,kread
	call glowssjcond(idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,kmax,phissj,1)

	write(6,"(1x,10f8.1)") idate(k),ut(k),glat(k),glong(k),f107a(k),f107(k),f107p(k),ap(k),ef(k),ec(k)
	do j=1,jmax
		write(6,"(f8.1,2e10.2)") zz(j)/1.e5,pedcond(k,j),hallcond(k,j)
	enddo
  enddo

!  enddo

stop

end program glowbasicssj
