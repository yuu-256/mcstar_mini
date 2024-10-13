function bintp(xx,n,x,y)
! 2nd order polynomial interpolation
!--- history
! 88. 1.14 created
!     3. 7 linear extrapolation for out-ou-boundary data
! 08. 9.19 modified Fortran 90 free-style form
!--- input
! xx    r     interpolation point
! n     r     nbr of data
! x   r(nn)   independent variable,  x(i) .lt. x(i+1)
! y   r(nn)   dependent   variable
!--- output
! bintp r     output
!--
  implicit none

! for input & output
  integer,intent(in)::n
  real(8),intent(in)::xx,x(n),y(n)
  real(8):: bintp
! for work
  integer::i,i1,i2,i3,k1,k2
!--exec
! 1 point
  if(n <= 1) then
     bintp=y(1)
     return
  endif

! 2 points or xx.le.x(1)
  if(xx <= x(1) .or. n==2) then
     bintp=y(1)+(y(2)-y(1))*(xx-x(1))/(x(2)-x(1))
     return
  endif
! xx.ge.xx(n)
  if(xx > x(n)) then
     bintp=y(n-1)+(y(n)-y(n-1))*(xx-x(n-1))/(x(n)-x(n-1))
     return
  endif

! 3 points
  do i=1,n
     if(xx < x(i)) goto 10
  enddo
  i=n
10 k1=max(1,i-1); k2=k1+2
  if(k2 > n) k2=n; k1=k2-2

  bintp=0.0d0; i2=k1; i3=k1+1
  do i=1,3
     i1=i2; i2=i3; i3=mod(i+1,3)+k1
     bintp=bintp+(xx-x(i2))*(xx-x(i3))/(x(i1)-x(i2))/(x(i1)-x(i3))*y(i1)
  enddo
  return
end function bintp
