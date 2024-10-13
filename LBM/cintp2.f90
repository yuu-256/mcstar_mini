function cintp2(xx,n,x,y)
! cubic interpolation (log-interpolation)
!--- history
! 87.12.21 created
! 88. 1. 9 change calculation algorithm of k1
! 96. 2. 7 log-interpolation by a. higurashi
!--- input
! xx    r     interpolation point
! n     r     nbr of data
! x   r(nn)   independent variable
! y   r(nn)   dependent   variable
!--- output
! cintp r     output
!--exec
  implicit none
! input & output
  integer,intent(in):: n
  real(8),intent(in):: xx,x(n),y(n)
  real(8):: cintp2

! for work
  integer:: i,i1,i2,i3,k1,k2,i4,i5,isign
  real(8):: xi,yy,dinc
  real(8):: expfn
!--exec
! 1 point
  if(n<=1) then
     cintp2=y(1)
     return
  endif
! 2 points
  if(n==2) then
     cintp2=y(1)+(y(2)-y(1))*(xx-x(1))/(x(2)-x(1))
     return
  endif
! 3 points
  if(n==3) then
     cintp2=0.d0
     i2=1
     i3=2
     do i=1,3
        i1=i2
        i2=i3
        i3=mod(i+1,3)+1
        cintp2=cintp2+(xx-x(i2))*(xx-x(i3))/(x(i1)-x(i2))/(x(i1)-x(i3))*y(i1)
     enddo
     return
  endif

! cubic interpolation
  dinc=x(2)-x(1)
  do i=1,n
     k1=i-2
     if(dinc >= 0.d0 .and. xx<=x(i)) exit
     if(dinc < 0.d0 .and. xx>=x(i)) exit
  enddo
  k1=max(1,i-2)
  k2=k1+3
  if(k2>n) then
     k2=n
     k1=k2-3
  endif

  cintp2=0.d0
  i2=k1
  i3=k1+1
  i4=k1+2
  i5=k1+3
  if((y(i2)<=0.d0) .or. (y(i3)<=0.d0) .or. (y(i4)<=0.d0) .or. (y(i5)<=0.d0)) then
     isign=-1
  else
     isign=1
  endif
  do i=1,4
     i1=i2
     i2=i3
     i3=i4
     i4=mod(i+2,4)+k1
     xi=x(i1)
     if(isign>0)then
        yy=log(y(i1))
     else
        yy=y(i1)
     endif
     cintp2=cintp2 &
          +(xx-x(i2))*(xx-x(i3))*(xx-x(i4))/(xi-x(i2))/(xi-x(i3))/(xi-x(i4))*yy
  enddo
  if(isign == 1) cintp2=expfn(cintp2)
  return
end function cintp2

