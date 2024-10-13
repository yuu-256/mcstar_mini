subroutine phasa2mc(ipol,mmax1,nang,ang,qrl,bx,bt,trf,b0)
!------ history
! 95. 7.11     created
!     8.10     qlr -> qrl
!    10.31     g <- g.old=g*(2l-1)
!mm210826 kmmax1 => mmax1
!mm210826 changed gsfc
!mm231113 phasa for mc from rpstar.pack
!------ input
! nang      i            number of scattering angles
! ang     r(knang)       scattering angles in degree
! qrl     r(knang,4)     elements of the phase matrix
! mmax1     i            maximum fourier order + 1
!------ output
! g      r(6,kmmax1)    gsf moments of phase matrix
! qfit    r(knang,6)
! f         r            truncation factor
!------
  use paras, only: knang,kpol2,klgn1,rad,pi
  implicit none
  integer,intent(in):: ipol
  integer,intent(in):: mmax1
  integer,intent(in):: nang
  real(8),intent(in):: ang(knang)
  real(8),intent(inout):: qrl(knang,kpol2)
  real(8),intent(out):: bx(klgn1,kpol2)
  real(8),intent(out):: bt(klgn1,kpol2)
  real(8),intent(out):: trf
  real(8),intent(out):: b0

  integer,parameter:: kg=5
  integer,parameter:: knang2=knang*kg
  integer,save::init=1
  integer,save::ng
  real(8),save::gw(kg),gx(kg)

  integer:: i,j,ii,m,n,ip,ip1
  real(8):: qcp(knang,kpol2),xx,ww,yy,qq(klgn1)
  real(8):: x(knang)!,y(knang),gsfc(knang2,kpol2,kmmax1)
!mm210826
  real(8):: gsfc(klgn1)
  
  real(8):: cnrm,p3,p2
  integer:: mmaxt,k1,npol
  real(8):: cintp2

  integer:: mm(1:kpol2)=(/0,2, 2,0,0,0/)
  integer:: nn(1:kpol2)=(/0,2,-2,0,2,2/)
  
  real(8) :: pl(klgn1)
  real(8) :: p(kpol2), bw(klgn1,kpol2)
!--exec
! shifted gaussian quadrature.
  if(init>=1) then
     init=0
     ng=kg
     call qgausn(gw,gx,ng)
     gw(1:ng)=gw(1:ng)*0.5d0
  endif
  
  npol=ipol
  if(npol>1) npol=6

! l,r -> cp representation
  x(1:nang)=cos(ang(1:nang)*rad)
  qcp(1:nang,1)=qrl(1:nang,1)
  if(ipol>1) then
    qcp(1:nang,2)=qrl(1:nang,2)+qrl(1:nang,3)
    qcp(1:nang,3)=qrl(1:nang,2)-qrl(1:nang,3)
    qcp(1:nang,4:npol)=qrl(1:nang,4:npol)
  end if

  bx(1:mmax1,1:kpol2)=0.d0  !mm210826
  bt(1:mmax1,1:kpol2)=0.d0  !mm210826
  ii=0
  do i=1,nang-1; do j=1,ng
     ii=ii+1
     xx=x(i)+gx(j)*(x(i+1)-x(i))
     ww=gw(j)*abs(x(i+1)-x(i))

     do ip=1,ipol
        m=mm(ip)
        n=nn(ip)
        call gsf(m,n,mmax1,xx,qq)
! p=q*(i)**(n-m)
!  =q*(-1)**((n-m)/2)

!mm210826
!        gsfc(ii,ip,1:mmax1)=qq(1:mmax1)*(-1)**((n-m)/2)        
        gsfc(1:mmax1)=qq(1:mmax1)*(-1)**((n-m)/2)
           
        yy=cintp2(xx,nang,x,qcp(1:nang,ip))
        do k1=1,mmax1
!mm210826           
!           bx(k1,ip)=bx(k1,ip)+yy*gsfc(ii,ip,k1)*ww
           bx(k1,ip)=bx(k1,ip)+yy*gsfc(k1)*ww
        enddo
     enddo


! off-diagonal elements
     do ip1=1,ipol/2
        ip=4+ip1
        m=mm(ip)
        n=nn(ip)
        call gsf(m,n,mmax1,xx,qq)
!mm210826        
!        gsfc(ii,ip,1:mmax1)=qq(1:mmax1)*(-1)**((n-m)/2)
        gsfc(1:mmax1)=qq(1:mmax1)*(-1)**((n-m)/2)
        yy=cintp2(xx,nang,x,qcp(1:nang,ip))
        do k1=1,mmax1
!           g(ip,l)=g(ip,l)+dble(2*l-1)*yy*gsfc(ii,ip,l)*ww
!mm210826
!           bx(k1,ip)=bx(k1,ip)+yy*gsfc(ii,ip,k1)*ww
           bx(k1,ip)=bx(k1,ip)+yy*gsfc(k1)*ww           
        enddo
     enddo
  enddo; enddo

  do k1=1,mmax1
     p2=(bx(k1,2)+bx(k1,3))*0.5d0
     p3=(bx(k1,2)-bx(k1,3))*0.5d0
     bx(k1,2)=p2
     bx(k1,3)=p3
  enddo

! normalization
  b0=bx(1,1)
  bx(1:mmax1,1:npol)=bx(1:mmax1,1:npol)/b0
  qrl(1:nang,1:npol)=qrl(1:nang,1:npol)/(b0*4.d0*pi)
   
!---
! Truncation factor
  mmaxt=mmax1-1
  trf=bx(mmax1,1)

! Truncation
  bt(1:mmaxt,1:4)=bx(1:mmaxt,1:4)-trf
  bt(1:mmaxt,5:6)=bx(1:mmaxt,5:6)

  bt(1:2,2:3)=0.d0

  cnrm=1.d0/(1.d0-trf)
  bt(1:mmaxt,1:kpol2)=bt(1:mmaxt,1:kpol2)*cnrm

! for Monte Carlo
  bw(1:mmaxt,1) = bt(1:mmaxt,1)
  if( ipol > 1 ) then
    bw(1:mmaxt,2) = bt(1:mmaxt,2) + bt(1:mmaxt,3)
    bw(1:mmaxt,3) = bt(1:mmaxt,2) - bt(1:mmaxt,3)
    bw(1:mmaxt,4:6) = bt(1:mmaxt,4:6)
  end if
  
  do i = 1, nang
    p = 0.d0
    do ip = 1, npol
      m = mm(ip)
      n = nn(ip)
      call gsf(m,n,mmaxt,x(i),pl(1:mmaxt))
      do k1 = 1, mmaxt
        p(ip) = p(ip) + (2*k1-1)*bw(k1,ip)*pl(k1)*(-1)**((n-m)/2.d0)
      end do
    end do
    qrl(i,1) = p(1)
    if( ipol > 1 ) then
      qrl(i,2) = ( p(2) + p(3) ) * 0.5d0
      qrl(i,3) = ( p(2) - p(3) ) * 0.5d0
      qrl(i,4:6) = p(4:6)
    end if
  end do 
  if(abs(bw(1,1))>0.d0) qrl(1:nang,1:npol) = qrl(1:nang,1:npol) / bw(1,1) / 4.d0 / pi

  return
end subroutine phasa2mc

