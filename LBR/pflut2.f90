!mm220221
subroutine pflut2(ipol,init,nang,npol,ang,phs,na,angr,cs,ph,angdp,err)
  !mm220118 created from pflut, add subgrid kn
  use paras, only : pi, rad
  implicit none
  integer :: init,nang,ipol,npol
  real(8) :: ang(nang), phs(nang,npol)
  integer :: na
  real(8) :: angr(na), cs(na)
  real(8) :: ph(na,npol), angdp(na)
  character (len=64) :: err
  !work
  !mm220118
  integer, parameter :: kn = 1
  integer :: nb
  real(8) :: angh((na-1)*kn+1), cs1(2)
  !end 220118
  integer :: m, initd,jpol
  real(8) :: a1, pd1, dp2((na-1)*kn+1), da, dx, x1, p1
  real(8) :: pint4c, aintp
  
  da = 180./(na-1)/kn !mm220118
  !mm220118
  nb = (na-1)*kn+1
  do m = 1, nb
    angh(m) = da*(m-1)
  end do
  
! phase function table
  do m=1,na
    a1=angr(m)
    initd=1
    do jpol = 1, npol
      ph(m,jpol)=pint4c(initd,a1,nang,ang,phs(:,jpol))
    end do
  enddo
! distribution table
  pd1=0
  dp2(1)=0
  do m=1,nb-1
    a1=angh(m)+da/2 !mm220118
    cs1(1) = cos(angh(m)*rad) !mm220118
    cs1(2) = cos(angh(m+1)*rad) !mm220118
    initd=1
    p1=pint4c(initd,a1,nang,ang,phs(:,1))
    !mm220118
    !pd1=pd1+p1*abs(cs(m+1)-cs(m))
    pd1=pd1+p1*abs(cs1(2)-cs1(1))
    dp2(m+1)=pd1
  enddo
  if(pd1.le.0) then
    err='pd1=0'
    return
  endif
  dp2(1:nb)=dp2(1:nb)/pd1 !mm220118
  ph(1:na,1:npol)=ph(1:na,1:npol)/pd1/2/pi
  if(init>0) return
! inverse distribution table
  angdp(1)=0
  angdp(na)=180
  dx=1.0/(na-1)
  do m=2,na-1
    x1=dx*(m-1)
    angdp(m)=aintp(x1,nb,dp2,angh) !mm220118
  enddo
end subroutine pflut2