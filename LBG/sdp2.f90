subroutine sdp2(nln,p,dp,np)
! derive absorption coefficient from pt table and fitting.
!--history
! 02.02.27  Modified knupt in ckdng.f (Zhang Hua)
! 02.11.25  Modified from knpt1.      (Miho Sekiguchi)
! 03.09.09  Modified from sdp for free source form.(Miho Sekiguchi)
! 07.03.28  Modified about dimension nln => knln. (Miho Sekiguchi)
! 08. 9.25 modified Fortran 90 free-style form
!--
  use paras, only: kp,knln
  implicit none

! input & output
  integer,intent(in):: nln       !! number of layer
  real(8),intent(in):: p(knln)   !! pressure in layer
  real(8),intent(out):: dp(knln) !! ratio of pressure grid
  integer,intent(out):: np(knln) !! number of pressure grid

! works
  integer:: l,ip
  real(8):: prs(1:26)=(/ &
       1.0130d+3,  0.63096d+3, 0.39811d+3,  0.25119d+3,  0.15849d+3, &
       1.0000d+2,  0.63096d+2, 0.39811d+2,  0.25119d+2,  0.15849d+2, &
       1.0000d+1,  0.63096d+1, 0.39811d+1,  0.25119d+1,  0.15849d+1, &
       1.0000d+0,  0.63096d+0, 0.39811d+0,  0.25119d+0,  0.15849d+0, &
       1.0000d-1,  0.63096d-1, 0.39811d-1,  0.25119d-1,  0.15849d-1, &
       1.0000d-2/)               !! pressure grid
!--exec
  do l=1,nln
     if(p(l) < prs(kp)) then
        np(l)=kp
        dp(l)=1.d0
        cycle
     endif
     do ip = 1, kp
        if (p(l) >= prs(ip)) exit
     enddo
     np(l)=ip
     if(np(l) > kp) np(l)=kp
     if(np(l) <= 1) np(l)=2
     dp(l)=log10(p(l)/prs(np(l)-1))/log10(prs(np(l))/prs(np(l)-1))
  enddo
end subroutine sdp2
