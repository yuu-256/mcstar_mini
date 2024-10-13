!mm220221
subroutine boxss(dbx,ds,ib,ps,is,pe,eps)
! ray trace in a box ( the box system case)
!--- history
! 88. 9.21 created
!tr211120     if(ak .le. 0) -> ak<0
!--- input
! dbx       r(kb1,3,2)  x,y,z-coordinates of interfaces.
! ds        r(3)     incident direction vector.
! ib        i(3)     current box location.
! ps        r(3)     current location.
! eps        r       convergence criterion for hitting wall
!--- output
! is         i       surface of emergence (+- 1,2,3 for x,y,z)
!                    >0 then positive x,y,z side relative to the box center
!                    <0 then negative x,y,z side relative to the box center
!                    =0 error
! pe        r(3)     point of interception.
!removed erc        i       ' '    no error.
!removed                    ' ='   p1=p2
!removed                    'box no'   no interception
!--- parameter
! kb1        i       max(kx,ky,kz)+1.
!$endi
!  parameter (kb1   =11)
  use paras
  implicit none

  real(8), intent(in) :: dbx(kb1,3,2), ds(3), ps(3), eps
  integer, intent(in) :: ib(3)
  real(8), intent(out) :: pe(3)
  integer, intent(out) :: is
!  character(len=64), intent(out) ::  erc
! work area
  real(8) :: bxc(3),p1(3),p2(3),w2(3)
  integer :: i,j1,j2,i1,i2,i3
  real(8) :: ak

! sizes of the box
  w2(1)  = dbx(ib(1),1,1)
  w2(2)  = dbx(ib(2),2,1)
  w2(3)  = dbx(ib(3),3,1)
! absolute box center coordinates
  bxc(1)  = dbx(ib(1),1,2)
  bxc(2)  = dbx(ib(2),2,2)
  bxc(3)  = dbx(ib(3),3,2)
! position of incidence on a wall relative to the center
  p1(1:3)=ps(1:3)-bxc(1:3)
  
!  erc = ' '
  
  if(abs(ds(1)) .gt. eps) then
     is = sign(1.,ds(1))
     p2(1)=w2(1)*is
     ak=(p2(1)-p1(1))/ds(1)
     p2(2)=p1(2)+ak*ds(2)
     p2(3)=p1(3)+ak*ds(3)
     if(abs(p2(2)) .le. w2(2) .and. abs(p2(3)) .le. w2(3)) then
! p2=p1
!        if(ak .le. eps*w2(1)) erc=' ='
! retrieve the absolute coordinates
        pe(1:3)=p2(1:3)+bxc(1:3)
        return
     end if
  end if
  
  if(abs(ds(2)) .gt. eps) then
     is = sign(1.,ds(2))
     p2(2)=w2(2)*is
     is = is*2
     ak=(p2(2)-p1(2))/ds(2)
     p2(3)=p1(3)+ak*ds(3)
     p2(1)=p1(1)+ak*ds(1)
     if(abs(p2(3)) .le. w2(3) .and. abs(p2(1)) .le. w2(1)) then
! p2=p1
!        if(ak .le. eps*w2(2)) erc=' ='
! retrieve the absolute coordinates
        pe(1:3)=p2(1:3)+bxc(1:3)
        return
     end if
  end if
  
  if(abs(ds(3)) .gt. eps) then
     is = sign(1.,ds(3))
     p2(3)=w2(3)*is
     is = is*3
     ak=(p2(3)-p1(3))/ds(3)
     p2(1)=p1(1)+ak*ds(1)
     p2(2)=p1(2)+ak*ds(2)
     if(abs(p2(1)) .le. w2(1) .and. abs(p2(2)) .le. w2(2)) then
! p2=p1
!        if(ak .le. eps*w2(3)) erc=' ='
! retrieve the absolute coordinates
        pe(1:3)=p2(1:3)+bxc(1:3)
        return
     end if
  end if
  
  is = 0
!  erc='box no'
  return
end subroutine boxss
