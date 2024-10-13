!mm220221
subroutine fbox(ps,bx,nb,ib,erc)
!c find the box containing the point ps
!c--- history
!c 88. 9.16  created
!c--- input
!c ps      r(3)        (x,y,z)
!c bx     r(kb1,3)     coordinates of box interfaces.
!c nb      i(3)        nbr of interfaces for (x,y,z)-axes.
!c--- output
!c ib      i(3)        box position.
!c erc     c*64        ' ' no error
!c--- parameter
!c kb1        i       max(kx,ky,kz)+1.
!c$endi
!      parameter (kb1   =11)
  use paras
  implicit none
  character(len=64) ::  erc
  real(8) :: ps(3),bx(kb1,3)
  integer :: nb(3),ib(3)
  integer :: l,i

  out : do l=1,3
     do i=1,nb(l)
        if(ps(l).ge.bx(i,l) .and. ps(l).le.bx(i+1,l)) then
           ib(l)=i
           cycle out
        endif
     enddo
     erc='fbox'
     return
  enddo out
  erc=' '
  return
end subroutine fbox