function bplnk(w,t)
! plank function with respect to wavelength
!--- history
! 89. 8. 1  created
! 08. 9.25 modified Fortran 90 free-style form
!--- input 
! w      r      wavelength (micron)
! t      r      absolute temperature (K)
!--- output
! bplnk  rf     plank function (W/m2/str/micron)
!--
  implicit none

! input 
  real(8),intent(in)::w,t

  real(8)::x,bplnk
!--exec
  if(w*t <= 0.d0) then
     bplnk=0.d0
  else
     x=1.438786d4/w/t
     bplnk=1.1911d8/w**5/(exp(x)-1.d0)
  endif
  return
end function bplnk
