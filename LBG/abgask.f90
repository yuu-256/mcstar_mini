subroutine abgask(iug,nln,nch,dintvl,nww,wlb,dp,np,amtpb,t,tau,err)
! calculate absorption coefficients.
!--- history
! 06.03.17 created (miho)
! 07.11.30 fix dimension (nln => knln, nch => kch) 
! 08.08.14 add dintvl
! 08. 9.26 modified Fortran 90 free-style form
! 08.12.22 Bug fixed
! 21. 1. 6 parameter knln
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
  use paras, only: kmol,kp,kt,kmol1,kww,knln,kch
  implicit none
! input
  integer,intent(in):: iug       !! device number
  integer,intent(in):: nln       !! number of layer
  integer,intent(in):: nch       !! number of channel
  real(8),intent(in):: dintvl    !! interval of wavenumber
  integer,intent(in):: nww       !! number of intervals
  real(8),intent(in):: wlb       !! wavelength at larger boundary 
  real(8),intent(in):: dp(knln)  !! ratio of pressure grid
  integer,intent(in):: np(knln)  !! pressure grid point
  real(8),intent(in):: amtpb(kmol1,knln) !! gas concentration
  real(8),intent(in):: t(knln)   !! temperature in layer

! output
  real(8),intent(out):: tau(knln,kch,kww) !! optical thickness
  character,intent(inout):: err*64

! work
  integer:: irec,iww,l,ich,imol,it
  real(4):: ak(kch,kp,kt,kmol1)
  real(8):: akt(kt,knln),knu(knln)
!--- exec
  tau(1:nln,1:nch,1:nww)=0.d0

! find the position of the beginning of the data.
! wlb = 10+10**(irec/dintvl)
  irec=int((log10(1.d4/wlb)-1.d0)*dintvl)
  do iww=1,nww
     irec=irec+1
     read(iug,rec=irec) ak(1:nch,1:kp,1:kt,1:kmol1)
     if(irec.lt.1) then
        err='wavenumber out of the range [0.2-1000 micron]'
        return
     endif

     do ich=1,nch
        do imol=1,kmol1
           do l=1,nln
              do it=1,kt
                 akt(it,l)=ak(ich,np(l)-1,it,imol) &
                      +(ak(ich,np(l),it,imol)-ak(ich,np(l)-1,it,imol))*dp(l)
              enddo
           enddo
           call tdok2(nln,akt,t,knu)
           tau(1:nln,ich,iww)=tau(1:nln,ich,iww)+knu(1:nln)*amtpb(imol,1:nln)
        enddo
     enddo
  enddo
  return
end subroutine abgask
