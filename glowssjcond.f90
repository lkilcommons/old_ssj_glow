subroutine glowssjcond(idates,uts,glats,glongs,f107as,f107s,f107ps,aps,efs,ecs,xncond,phij,usephij)

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.98 beta, 1/2017

! Adapted from glowbasic by Liam Kilcommons, 3/2017

! Basic single-processor driver for the GLOW model.
! Uses MSIS/IRI for input.
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

!f2py depend(xncond) idates
!f2py depend(xncond) uts
!f2py depend(xncond) glats
!f2py depend(xncond) glongs
!f2py depend(xncond) f107as
!f2py depend(xncond) f107s
!f2py depend(xncond) f107ps
!f2py depend(xncond) aps,efs,ecs
!f2py depend(xncond) phij

  use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst,ncond
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
  use cglow,only: iscale,jlocal,kchem,xuvfac
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,tei,tpi,tir,ecalc,zxden,zeta,zceta,zlbh
  use cglow,only: pedcond,hallcond
  use cglow,only: cglow_init, cglow_init_backup, cglow_init_cond, cglow_free_cond
  use cglow,only: cglow_clear, cglow_clear_exsect
  use cglow,only: data_dir

  implicit none

  character(len=1024) :: iri90_dir
  
  !Input arguments are arrays to fill module data from cglow for each ssj spectra
  integer, intent(in) :: xncond
  real, intent(in) :: idates(xncond),uts(xncond),glats(xncond),glongs(xncond)
  real, intent(in) :: f107as(xncond),f107s(xncond),f107ps(xncond),aps(xncond)
  real, intent(in) :: efs(xncond),ecs(xncond)
  real, intent(in) :: phij(xncond,19)
  integer, intent(in):: usephij
  real,allocatable :: z(:)                    ! glow height coordinate in km (jmax)
  real,allocatable :: zun(:), zvn(:)          ! neutral wind components (not in use)
  real,allocatable :: outf(:,:)               ! iri output (11,jmax)

  real :: rz12,stl,fmono,emono
  real :: x0,x1,y0,y1,x,y
  real :: d(8), t(2), sw(25), oarr(30), enerj(19), delj(19)
  integer :: l,j,jj,ijf,jmag,iday,mmdd,i,ii,n,k,ix,m,itail
  integer :: instance,iostatus
  logical :: jf(12)
  data sw/25*1./
!
! DMSP SSJ Channel Center Energies in eV     
!
  data enerj/  30000.,20400.,13900.,9450.,&
               6460.,4400.,3000.,2040.,&
               1392.,949.,646.,440.,&
               300.,204.0,139.,95.,&
               65.,44.,30./
!
! DMSP SSJ Channel Effective Energy Widths in eV (as in Hardy 1984)
!
  data delj/ 9600.0,8050.0,5475.0,3720.0,&
               2525.0,1730.0,1180.0,804.0,&
               545.5,373.0,254.5,173.0,&
               118.0,80.5,54.5,37.0,&
               25.5,17.5,14.0/

!
! Initialize standard switches:
!
  iscale=1
  xuvfac=3.
  kchem=4
  jlocal=0
  itail=0
  fmono=0.
  emono=0.

!
! Initialize DMSP SSJ specific switches
!
!  usephij=0
  
!
! Set data directories:
!
  data_dir    = 'data/'
  iri90_dir   = 'data/iri90/'
!
! Set number of altitude levels:
!
  jmax = 102

!
! Allocate local arrays:
!
  allocate(z(jmax))
  allocate(zun(jmax))
  allocate(zvn(jmax))
  allocate(outf(11,jmax))
!
! Call CGLOW_INIT (module CGLOW) to set array dimensions and allocate use-associated variables:
! (This was formerly done using common blocks, including common block /cglow/.)
!
  ncond = xncond !Set module variable to passed dimension
  if (.not. allocated(zz)) call cglow_init
  if (allocated(pedcond)) then
    call cglow_free_cond
  end if 
  call cglow_init_cond
    
  !call cglow_init_backup
  
!
! Loop over inputs
!
  
  pedcond(:,:)=0.
  hallcond(:,:)=0.

  write(6,*) "Computing GLOW conductivities for ",ncond," points."      

  do instance=1,ncond

    z(:)=0.;zun(:)=0.;zvn(:)=0.;outf(:,:)=0.;
    call cglow_clear
    !call cglow_clear_exsect
    !if (instance .gt. 1) call cglow_load_exsect
    
! 
! Set correspond module data parameters to parameters passed in
!
    idate=idates(instance)
    ut=uts(instance)
    glat=glats(instance)
    glong=glongs(instance)
    f107a=f107as(instance)
    f107=f107s(instance)
    f107p=f107ps(instance)
    ap=aps(instance)
    ef=efs(instance)
    ec=ecs(instance)
    !write(6,*) idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec

!
! Call EGRID to set up electron energy grid:
!
    call egrid (ener, del, nbins)

!
! Calculate local solar time:
!
    stl = ut/3600. + glong/15.
    if (stl < 0.) stl = stl + 24.
    if (stl >= 24.) stl = stl - 24.
!
! Call MZGRID to use MSIS/NOEM/IRI inputs on default altitude grid:
!
    call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                 z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
!
! Call MAXT to put auroral electron flux specified by namelist input into phitop array:
!
    phitop(:) = 0.
    if (ef>.001 .and. ec>1.) call maxt (ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
    !Fill in appropirate bins with SSJ flux
    n=19
    if (usephij.eq.1) then
      !write(6,*) "Inserting DMSP fluxes into phitop"
      do m=1,nbins
        if ((ener(m).ge. (enerj(19)-delj(19))).and.(ener(m).le.enerj(1)+delj(1))) then
          !If within bounds of SSJ energy grid
          if (ener(m).gt.enerj(1)) then
            !Grid energy higher than highest center energy but w/in highest channel width
            phitop(m) = phij(instance,1)
          else if (ener(m).lt.enerj(19)) then
            !Grid energy lower than lowest center energy but w/in lowest channel width
            phitop(m) = phij(instance,19)
          else
            !Check if need to advance index
            if (ener(m).gt.enerj(n-1)) then
              n=n-1
            endif
            !Grid energy m is within SSJ coverage, linearly interpolate
            x = ener(m); x0 = enerj(n); x1 = enerj(n-1);
            y0 = phij(instance,n); y1 = phij(instance,n-1);
            phitop(m) = y0 + (x-x0)*(y1-y0)/(x1-x0)

          endif
          ! if (mod(instance,300).eq.0) then
          !   write(6,*) ener(m),'jchannel: ',n,'jflux: ',phij(instance,n),'gflux: ',phitop(m)
          ! endif
        else if (ener(m) .gt. enerj(1)+delj(1)) then 
          exit
        endif
      enddo
    endif


!
! Fill altitude array, converting to cm:
!
    zz(:) = z(:) * 1.e5     ! km to cm at all jmax levels

    !write(6,*) 'before: ', instance, z(j), zo(j), zo2(j), zn2(j), &
    !           zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j)

!
! Call GLOW to calculate ionized and excited species, airglow emission rates,
! and vertical column brightnesses:
!
    call glow
!
! Call CONDUCT to calculate Pederson and Hall conductivities:
!
    do j=1,jmax
      call conduct (glat, glong, z(j), zo(j), zo2(j), zn2(j), &
                    zxden(3,j), zxden(6,j), zxden(7,j), ecalc(j), &
                    ztn(j), zti(j), zte(j), &
                    pedcond(instance,j), hallcond(instance,j))
    enddo
    !j=30
    !write(6,*) 'after: ',instance, z(j), zo(j), zo2(j), zn2(j), &
    !           zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
    !           pedcond(instance,j), hallcond(instance,j)

!
! Output section:
!

    !write(6,*) ener
    !write(6,"(1x,i7,9f8.1)") idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
    !write(6,"('   Z     Tn       O        N2        NO      Ne(in)    Ne(out)  Ionrate      O+       O2+      NO+       N(2D)    Pederson   Hall')")
    !write(6,"(1x,0p,f5.1,f6.0,1p,12e10.2)") (z(j),ztn(j),zo(j),zn2(j),zno(j),ze(j), &
    !   ecalc(j),tir(j),zxden(3,j),zxden(6,j),zxden(7,j),zxden(10,j),pedcond(instance,j),hallcond(instance,j),j=1,jmax)
    !write(6,"('   Z      3371    4278    5200    5577    6300    7320   10400    3644    7774    8446    3726    LBH     1356    1493    1304')")
    !write(6,"(1x,f5.1,15f8.2)")(z(j),(zeta(ii,j),ii=1,15),j=1,jmax)

  enddo

  deallocate(z)
  deallocate(zun)
  deallocate(zvn)
  deallocate(outf)
  
return

end subroutine glowssjcond
