subroutine smllop(x,cr,ci,qext,qsca)
! optical cross sections for small particles (x<0.1)
! in this approximation
!     q(ang,1)=3.0/8.0/pi*qsca
!     q(ang,2)=q(ang,1)*cos(ang)**2
!     q(ang,3)=q(ang,1)*cos(ang)
!     q(ang,4)=0
!--- Reference
! H.C. van de Hulst, 1981: Light Scattering by Small Particles, Dover. 
! Section 14.21. Series Expansion (P.270) 
!--- history
! 95. 9.19 created
! 08. 9.25 modified Fortran 90 free-style form
!--- input
! x       r       size parameter (<1.5)
! cr      r       m= cr - i ci
! ci      r       ci>0
!--- output
! qext    r       extinction efficiency factor
! qsca    r       scattering efficiency factor
!---
  implicit none

! input & output
  real(8),intent(in):: x         !! size parameter
  real(8),intent(in):: cr        !! refractive index (real part)
  real(8),intent(in):: ci        !! refractive index (imaginary part)
  real(8),intent(out):: qext     !! extinction efficiency factor
  real(8),intent(out):: qsca     !! scattering efficiency factor

! work
  real(8):: a,b,g,z1,z2,e1,e3,e4,s4
!--
  a=(cr**2+ci**2)**2
  b=cr**2-ci**2
  g=cr*ci
  z1=a+4.d0*b+4.d0
  z2=4.d0*a+12.d0*b+9.d0
  e1=24.d0*g/z1
  e3=g*(4.d0/15.d0+20.d0/3.d0/z2+24.d0/5.d0/z1**2*(7.d0*a+4.d0*b-20.d0))
  e4=8.d0/3.d0/z1**2*((a+b-2.d0)**2-36.d0*g**2)
  s4=8.d0/3.d0/z1**2*((a+b-2.d0)**2+36.d0*g**2)
  qext=e1*x+e3*x**3+e4*x**4
  qsca=s4*x**4
  return
end subroutine smllop
