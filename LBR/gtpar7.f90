subroutine gtpar7(ius,nl,nptc,ins,iver,cnp,ispcvp,rfracp,asphr,rop, &
     dryap,nawp,awcrp,nv,nwlv,wlv,rfi)
! get parameters of particle polydispersions
!--- history
! 95. 1.25  created from getph3
! 95. 1.24  dintpm -> aerd (now includs also aerosol water uptake data)
! 95  9.13  with modified aerd file format
! 96. 3.11  added non-spherical parameters
!     5. 5  awcrp(*,*,2)
! 08.10.31  modified fortran 90 version.
! 12.02.28  include nonspherical flag
!--- input
! ius       i         read unit number of particle parameter file
! nl        i         number of layers (come from mlatm)
!--- output
! nptc      i         number of particle polydispersions
! ins    i(kpcl)      nonspherical flag
! iver   i(kpcl)      type of vertical profile
! cnp   r(knl,kver)   dry volume concentration profile (relative unit)
! ispcvp  i(3,kptc)   fundamental materials for 3 comp. internal mixture (1-8)
! rfracp  r(3,kptc)   dry volume fraction of the dry mixture
! asphr  r(3,kptc)    non-spherical parameters (x0, g, r) (in case of ins=3)
! rop     r(kptc)     paricle density relative to water
! dryap r(6,4,kptc)   dv/dlnr parameters for the dry mixture
!                     see vlspc2
!                     c-values (coefficients of volume spectrum) are ralative
! nawp    i(kptc)      number of aw
! awcrp r(kaw,kptc,2)  1: aw      water activity (see shettle and fenn.)
!                     2: rmmd    rmmd
! nv        i         number of fundamental species (1-8)
! nwlv      i         number of wavelengths for refractive index tables
! wlv    r(kwl)       wavelengths (micron)
! rfi   r(kwl,knv,2)  refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                     with log-regular wavelength interval
!---
  use paras,only: kver,kptc,knl,kaw,kwlv,knv
  integer,intent(in):: ius
  integer,intent(in):: nl
  integer,intent(out):: nptc
  integer,intent(out):: ins(kptc)
  integer,intent(out):: iver(kptc)
  real(8),intent(out):: cnp(knl,kver)
  real(8),intent(out):: rfracp(3,kptc)
  integer,intent(out):: ispcvp(3,kptc)
  real(8),intent(out):: asphr(3,kptc)
  real(8),intent(out):: rop(kptc)
  real(8),intent(out):: dryap(6,4,kptc)
  integer,intent(out):: nawp(kptc)
  real(8),intent(out):: awcrp(kaw,kptc,2)
  integer,intent(out):: nv
  integer,intent(out):: nwlv
  real(8),intent(out):: wlv(kwlv)
  real(8),intent(out):: rfi(kwlv,knv,2)

! work
  character(1)::ch
  integer::i,j,l,iv,iptc
  integer:: nmode, nver
  real(8)::rmin,rmax
!---
  rewind ius
!  number concentration profiles (relative unit)
!  number is conservative in growing aerosols
  read(ius,'(a1)') ch
  read(ius,*) rmin,rmax
  read(ius,*) nver, nptc
  do iv=1,nver
     read(ius,'(a1)') ch
     read(ius,*) (cnp(l,iv),l=1,nl)
  enddo
! loop for aerosol mptcs
  do iptc=1,nptc
     read(ius,'(a1)') ch
     read(ius,*) ins(iptc)
     read(ius,*) iver(iptc)
     if(ins(iptc)==0.or.ins(iptc)==1.or.ins(iptc)==3) then
        read(ius,*) (ispcvp (i,iptc),i=1,3)
        read(ius,*) (rfracp (i,iptc),i=1,3)
     endif

     if(ins(iptc)==3) &
!! Cuzzi and Pollack
     read(ius,*) (asphr(i,iptc),i=1,3)

     read(ius,*) rop(iptc)
     read(ius,*) nmode
     dryap(1,2,iptc)=nmode
     dryap(1,3,iptc)=rmin
     dryap(1,4,iptc)=rmax
     do j=1,nmode
        read(ius,*) dryap(2:6,j,iptc)
     enddo
     read(ius,*) nawp(iptc)
     if(nawp(iptc) > 0) then
        read(ius,'(a1)') ch
        do i=1,nawp(iptc)
           read(ius,*) awcrp(i,iptc,1:2)
        enddo
     endif
  enddo

! get refractive index table
  read(ius,'(a1)') ch
  read(ius,*) nv,nwlv
  read(ius,*) wlv(1:nwlv)
  do iv=1,nv
     read(ius,'(a1)') ch
     read(ius,*) rfi(1:nwlv,iv,1)
     read(ius,*) rfi(1:nwlv,iv,2)
  enddo
  return
end subroutine gtpar7
