function rawb(rh1,r,rho,naw,aw,rmmd)
! growth of wet aerosols
! raw is r(aw) in eq.(2) in shettle and fenn.
!--- history
! 94.12.26 created by i. lensky
! 95. 9.13 modified by t. nakajima
! 96. 5. 5 drop caw, number of iterations 10 -> 6
! 08. 9.19 modified Fortran 90 free-style form
!--- input
! rh1    r           relative humidity 0 <= rh1 < 1
! r      r           dry aerosol radius (cm)
! rho    r           paricle density relative to water
! naw    i           number of aw (if 0 then no growth)
! aw     r(kaw)      water activity
! rmmd   r(kaw)      Hanel's water uptake data
!--- output
! raw    r           r(aw) mode radius of wet aerosols (cm)
!--
  use paras, only: kaw
  implicit none

! input & output
  real(8),intent(in):: rh1       !! relative humidity
  real(8),intent(in):: r         !! dry aerosol radius [cm]
  real(8),intent(in):: rho       !! particle density relative to water
  integer,intent(in):: naw       !! number of aw
  real(8),intent(in):: aw(kaw)   !! water activity
  real(8),intent(in):: rmmd(kaw) !! Hanel's water uptake data
  real(8):: rawb                 !! mode radius of wet aerosols [cm]

! work
  integer:: it,iaw
  real(8):: aw1,caw,r1
!--exec
  if (naw < 0 .or. rh1 < aw(1)) then
     rawb = r
     return
  endif
! iteration
  aw1 = rh1
  do it=1,6
! eq.(5)  of shettle & fenn
     if(aw1 <= aw(1)) then
        iaw=1
     else if(aw1 >= aw(naw)) then
        iaw=naw
     else
        do iaw=1,naw-1
           if(aw1 >= aw(iaw) .and. aw1 < aw(iaw+1)) exit
        enddo
     endif
     caw=log(rmmd(iaw+1)/rmmd(iaw))/log((1-aw(iaw+1))/(1-aw(iaw)))
     r1=rmmd(iaw)*((1-aw1)/(1-aw(iaw)))**caw
     rawb=r*(1.d0+rho*r1)**(1.d0/3.d0)
! aw1 is aw, eq.(4) in shettle & fenn. rawb in cm.
     aw1=rh1*exp(-0.001056d0/(1.d4*rawb))
  enddo
  return
end function rawb
