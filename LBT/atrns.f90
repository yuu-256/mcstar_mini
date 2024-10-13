!mm220221
subroutine atrns(sc,fi,ds,eps)
!c get scattering direction vector
!c--- history
!c 88. 9.16 created
!c--- input
!c sc     r      scattering angle (degree)
!c fi     r      azimuthal angle (degree)
!c ds    r(3)    incident direction vector
!c eps    r      criterion for abs(ds(3))=1
!c--- output
!c ds            scattering direction vector
!c$endi
  use paras
  implicit none

  real(8), intent(inout) :: ds(3)
  real(8), intent(in) :: sc,fi,eps
  real(8) :: ct1,st1,cf1,sf1,ct2,st2,cf2,sf2
  real(8) :: vn

  ct1=cos(sc*rad)
  st1=sin(sc*rad)
  cf1=cos(fi*rad)
  sf1=sin(fi*rad)
  ct2=ds(3)
  st2=sqrt(sum(ds(1:2)**2))
  if(st2.le.eps) then
     cf2=1.d0
     sf2=0.d0
  else
     cf2=ds(1)/st2
     sf2=ds(2)/st2
     vn=sqrt(cf2**2+sf2**2)
     cf2=cf2/vn
     sf2=sf2/vn
  endif
  ds(1)=st1*(cf1*ct2*cf2-sf1*sf2)+ct1*st2*cf2
  ds(2)=st1*(cf1*ct2*sf2+sf1*cf2)+ct1*st2*sf2
  ds(3)=-st1*cf1*st2+ct1*ct2
 return
end subroutine atrns