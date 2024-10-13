!mm220221
subroutine abgasx(init,iug,nln,nch,nbnd,dp,np,amtpb,t,wgt,tau)
! calculate absorption coefficients.
!--- history
! 06.03.17 created (miho)
! 07.11.30 fix dimension (nln => knln, nch => kch)
! 08.08.14 add dintvl
! 08. 9.26 modified fortran 90 free-style form
! 08.12.22 bug fixed
! 10. 1.15 modified for parag file (mstrnx)
!--- input
! iug      i     read unit number of gas line absorption para file
! nln      i     number of layers
! wla      r     the lower boundary of wavelength
! wlb      r     the upper boundary of wavelength
! dp     r(nln)  residual of the pressure
! np     i(nln)  grid number of the pressure
! t      r(nln)  temperature of the layer
!--- output
! tau  r(nln,kch,kww) optical thickness
!---
  use paras, only: kmol,kt,kmol1,kww,knln,kch,kbnd,kp,kt

  implicit none
! input
  integer,intent(inout):: init   !! initialization flag
  integer,intent(in):: iug       !! device number
  integer,intent(in):: nln       !! device number
  integer,intent(in):: nbnd       !! device number
!x  integer :: nbnd       !! device number
  integer,intent(inout):: nch     !! device number
  real(8),intent(in):: dp(knln)  !! device number
  integer,intent(in):: np(knln)  !! device number
  real(8),intent(in):: amtpb(kmol1,knln) !! gas concentration
  real(8),intent(in):: t(knln)   !! temperature in layer

! output
  real(8),intent(out):: wgt(kch)          !! weight
  real(8),intent(out):: tau(knln,kch,kww) !! optical thickness
!  character,intent(inout):: err*64

! work
!  integer:: irec
  integer::iww,l,ich,imol,it
  integer:: nww
  real(8):: akt(kt,knln),knu(knln)

! gttbl
  integer:: i,iw,ia,ip
  integer,parameter:: kflg=7
  integer,save:: nch0(kbnd)
  integer,save:: iflgb(kflg,kbnd)
  real(8),save:: wgt0(kch,kbnd)
  integer,save:: nabs(kbnd),iabs(kmol,kbnd)
  real(8),save:: akd(kch,kp,kt,kmol,kbnd)
  real(8),save:: skd(kch,kp,kt,kbnd)

!--- exec
  nww=1
  tau(1:nln,1:kch,1:nww)=0.d0
  if(init==1) then
     init=0
!!!                  !!!
!!! from gttbl2.f90  !!!
!!!                  !!!
! read gas parameters
  do i=1,9
     read(iug,*)
  enddo

! quantities for each band
  do iw=1,kbnd
!! optical properties flag
     read(iug,*)
     read(iug,*) iflgb(1:kflg,iw)
!! number of subintervals
     read(iug,*)
     read(iug,*) nch0(iw)

!! weights for channels
     read(iug,*)
     read(iug,*) wgt0(1:nch0(iw),iw)

!! molecules
     read(iug,*)
     read(iug,*) nabs(iw)

!!! major gas absorption
     if(nabs(iw)>0) then
        do ia=1,nabs(iw)
           read(iug,*) iabs(ia,iw)
           do it=1,kt; do ip=1,kp
              read(iug,*) akd(1:nch0(iw),ip,it,ia,iw)
           enddo; enddo
        enddo
     endif
!!! h2o continuum
     if(iflgb(5,iw)>0) then
        read(iug,*)
        do it=1,kt; do ip=1,kp
           read(iug,*) skd(1:nch0(iw),ip,it,iw)
        enddo; enddo
     endif

!!! cfc absorption
     if(iflgb(7,iw)>0) then
        do l=1,7
           read(iug,*)
        enddo
     endif
  enddo
  close(iug)
!!!                !!!
!!! gttbl2.f90 end !!!
!!!                !!!

  endif
  nww=1

!x lamda=0.556 micron nbnd=24
!x  nbnd=24

  nch=nch0(nbnd)
  wgt(1:nch)=wgt0(1:nch,nbnd)

! find the position of the beginning of the data.
! wlb = 10+10**(irec/dintvl)

  do iww=1,nww
     do ich=1,nch
        i=1
        do imol=1,kmol
           if(iabs(i,nbnd)==imol) then
              do l=1,nln
                 do it=1,kt
                    akt(it,l)=akd(ich,np(l)-1,it,i,nbnd) &
                         +(akd(ich,np(l),it,i,nbnd)-akd(ich,np(l)-1,it,i,nbnd))*dp(l)
                 enddo
              enddo
              call tdok2(nln,akt,t,knu)
              tau(1:nln,ich,iww)=tau(1:nln,ich,iww)+knu(1:nln)*amtpb(imol,1:nln)
              i=i+1
           endif
        enddo
        if(iabs(1,nbnd)==1) then
           do l=1,nln
              do it=1,kt
                 akt(it,l)=skd(ich,np(l)-1,it,nbnd) &
                      +(skd(ich,np(l),it,nbnd)-skd(ich,np(l)-1,it,nbnd))*dp(l)
              enddo
           enddo
           call tdok2(nln,akt,t,knu)
          tau(1:nln,ich,iww)=tau(1:nln,ich,iww)+knu(1:nln)*amtpb(kmol1,1:nln)
        endif
     enddo
  enddo
  return
end subroutine abgasx
