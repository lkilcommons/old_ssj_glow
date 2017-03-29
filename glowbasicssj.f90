program glowbasicssj

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.98 beta, 1/2017

! Adapted from glowdriver by Stan Solomon, 2/2016

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
  
  use cglow,only: pedcond,hallcond

  implicit none

  real :: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
  real :: phissj(1,19)
  integer :: i,iostatus,j
  !real,allocatable :: pedcond(:),hallcond(:)

!  do instance=1,10000
!
! Get input values:
!
!    write(6,"('Enter date, UT, lat, lon, F107a, F107, F107p, Ap, Ef, Ec, PhiSSJ( 19 values)')")
    read(5,*,iostat=iostatus) idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,(phissj(1,i),i=1,19)
!    if (iostatus /= 0) stop
!
! Call subroutine
!
    call glowssjcond(idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,phissj,1)

!  enddo

stop

end program glowbasicssj
