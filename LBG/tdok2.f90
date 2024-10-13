subroutine tdok2(nln,akt,t,knu)
!----------------------------------------------------------------------c
!                using Shi's formula: (t/to)**(a+bt)
!----------------------------------------------------------------------c
  use paras, only: kt,knln
  implicit none
  integer,intent(in):: nln       !! number of layer

! 09.05.18 Bug fixed (Miho, reported by Hiro Masunaga)
  real(8),intent(inout):: akt(kt,knln)   !! absorption coefficients at T grid
  real(8),intent(in):: t(knln)   !! temperature in layer
  real(8),intent(out):: knu(knln) !! absorption coefficients
! works
  integer:: l
  real(8):: a,b,at1,at2,bt1,bt2
!--exec
  at1 = log10(200.d0/260.d0); at2 = log10(320.d0/260.d0)
  
  do l=1,nln
     if(akt(2,l) <= -14.d0.and.akt(1,l) > -14.d0.and. &
          akt(3,l) > -14.d0) akt(2,l)=(akt(1,l)+akt(3,l))*0.5d0
     bt1 = akt(1,l)-akt(2,l)
     bt2 = akt(3,l)-akt(2,l)
     if(bt1 == 0.d0 .and. bt2 == 0.d0) then
        knu(l)=10**akt(2,l); cycle
     endif
     b   = (bt2/at2-bt1/at1)/120.d0
     a   = bt2/at2 - b*320.d0
     
     knu(l)=10.d0**(akt(2,l))*(t(l)/260.d0)**(a+b*t(l))
  enddo
  return
end subroutine tdok2
