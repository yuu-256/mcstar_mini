function plgd(init,x)
! Legendre polynomials
!--- history
! 87.11.12 created
! 90. 1.16 save statement
!--- input
! init   i    if 1 then l=0
!             if 0 then l=l+1
! x      r    (-1,1)
!--- out
! init=0
  implicit none
  integer,intent(inout)::init
  real(8),intent(in)::x
  integer,save::l
  real(8),save::pl,pl1,pl2
  real(8)::plgd
!--
  if(init > 0) then
     init=0; l=-1
  endif
  l=l+1
  if(l == 0) then
     pl=1
  elseif(l == 1) then
     pl1=pl; pl=x
  else
     pl2=pl1; pl1=pl
     pl=(dble(2*l-1)*x*pl1-dble(l-1)*pl2)/dble(l)
  endif
  plgd=pl
  return
end function plgd
