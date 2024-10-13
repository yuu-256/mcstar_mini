subroutine lgndf3(lmax1,n,x,y,g)
! Legendre expansion (same as lgndf2 but generating g-moments).
!--- history
! 90. 1.20 created from lgndf2, use expfn.
!          direction of integration from x(1) to x(n).
! 91. 2.16 strip ng from save statement.
! 92. 4. 3 kill the statement of gw=gw/2 after qgausn
!     6.22 add gw/2 again because the above change is mistake
! 08.10.20 modified fortran 90 version.
!--- input
! lmax1      i      maximum order + 1.
! n          i      number of data on (-1, 1). .ge. 4.
! x        r(na)    independent variables on (-1, 1)
! y        r(na)    y(x).
!--- output
! g    r(lmax1)     y = sum(l1=1,lmax1) (2*l1-1)*g(l1)*pl(l1)
!                      where pl(l1) is (l1-1)th order legendre
!                      polynomial.
!--
  implicit none
  integer,intent(in)::lmax1,n
  real(8),intent(in)::x(n),y(n)
  real(8),intent(out)::g(lmax1)
! working areas.
  integer,save::init=1
  integer,parameter::ng=5
  real(8),save::gw(ng),gx(ng)
  integer::i,j,i1,i2,i3,i4,isign,l1
  real(8)::x1,x2,x3,x4,alp1,alp2,alp3,alp4,xx,ww,pp,pl,pl1,pl2
  real(8)::expfn
!--exec
! shifted gaussian quadrature.
  if(init >= 1) then
     init=0
     call qgausn(gw,gx,ng)
     gw(1:ng)=gw(1:ng)/2.0d0
  endif
! clear
  g(1:lmax1)=0.d0
! loop for angle
  do i=1,n-1
! cubic interpolation
     if(i <= 2) then
        i1=1; i4=4
     else
        if(i <= n-2) then
           i1=i-1; i4=i+2
        else
           i1=n-3; i4=n
        endif
     endif

     i2=i1+1; i3=i2+1
     x1=x(i1); x2=x(i2); x3=x(i3); x4=x(i4)
     if( (y(i1) <= 0.) .or. (y(i2) <= 0.) .or. (y(i3) <= 0.) .or.  &
          (y(i4) <= 0.) ) then
        isign=-1; alp1=y(i1); alp2=y(i2); alp3=y(i3); alp4=y(i4)
     else
        isign=1;  alp1=log(y(i1)); alp2=log(y(i2))
        alp3=log(y(i3)); alp4=log(y(i4))
     endif

! loop for gaussian integration
     do j=1,ng
!! interpolated value of y
        xx=x(i)+gx(j)*(x(i+1)-x(i))
        ww=gw(j)*(x(i+1)-x(i))
        pp =  (xx-x2)*(xx-x3)*(xx-x4)/(x1-x2)/(x1-x3)/(x1-x4)*alp1 &
             +(xx-x1)*(xx-x3)*(xx-x4)/(x2-x1)/(x2-x3)/(x2-x4)*alp2 &
             +(xx-x1)*(xx-x2)*(xx-x4)/(x3-x1)/(x3-x2)/(x3-x4)*alp3 &
             +(xx-x1)*(xx-x2)*(xx-x3)/(x4-x1)/(x4-x2)/(x4-x3)*alp4
        if(isign == 1) pp=expfn(pp)

! legendre sum
        pl1=0.d0; pl=1.d0
        g(1)=g(1)+pp*ww
        if(lmax1 >= 2) then
           do l1=2,lmax1
              pl2=pl1; pl1=pl
              pl=(dble(2*l1-3)*xx*pl1-dble(l1-2)*pl2)/dble(l1-1)
              g(l1)=g(l1)+pp*pl*ww
           enddo
        endif
     enddo
  enddo
  return
end subroutine lgndf3
