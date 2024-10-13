function pint4c(init,ang1,na,ang,p)
! interpolation of the phase function.
!--- history
! 89.11. 8 created from pint3. change interpolation x-> ang
! 90. 1.23 debug
! 08.10.20 modified fortran 90 version.
! 21. 1. 7 p as one dimensional variable
!--- input
! init       i      1 then search ang1-interval else not search.
! ang1       r      scattering angle in degree for interpolation.
! na         i      no. of scattering angles.
! ang      r(na)    scattering angles in degree.
! p     r(kna)    phase function
!--- output variables
! init       i      0
! pint4      r      interpolated value of p at ang1.
!--- variables for the routine.
!  use paras, only: knang, knln, pi, rad
!  use paras, only: knang, kpol
  implicit none

! input & output  
  integer,intent(inout):: init
  real(8),intent(in):: ang1
  integer,intent(in):: na
!  integer,intent(in):: l
!  real(8),intent(in):: ang(knang)
!  real(8),intent(in):: p(knang,kpol)
  real(8),intent(in):: ang(na)
  real(8),intent(in):: p(na)
  real(8):: pint4c

! work
  integer,save:: i1,i2,i3
  real(8):: alp1,alp2,alp3
  real(8):: x1,x2,x3,xx
  integer:: i,isign
  real(8):: pp,expfn
!--exec
! search near angles  
  if(init >= 1) then
     init=0
     do i=1,na-1
        if((ang1-ang(i))*(ang1-ang(i+1)) <= 0.d0) then
           if(i<=1) then
              i1=1; i3=3
           elseif(i<na-1) then
              i1=i-1; i3=i+1
           else
              i1=na-2; i3=na
           endif
           i2=i1+1
           exit
        endif
     enddo
     if(i==na) then
        i1=na-2; i2=na-1; i3=na
     endif
  endif
  xx=ang1
  x1=ang(i1)
  x2=ang(i2)
  x3=ang(i3)
!  alp1=p(i1,l)
!  alp2=p(i2,l)
!  alp3=p(i3,l)
  alp1=p(i1)
  alp2=p(i2)
  alp3=p(i3)
  isign=-1
  if(alp1 > 0.d0 .and. alp2 > 0.d0 .and. alp3 > 0.d0) then
     isign=1
     alp1=log(alp1)
     alp2=log(alp2)
     alp3=log(alp3)
  endif
  pp=(xx-x2)*(xx-x3)/(x1-x2)/(x1-x3)*alp1 &
       +(xx-x1)*(xx-x3)/(x2-x1)/(x2-x3)*alp2 &
       +(xx-x1)*(xx-x2)/(x3-x1)/(x3-x2)*alp3
  if(isign >= 1) pp=expfn(pp)
  pint4c=pp
  return
end function pint4c
