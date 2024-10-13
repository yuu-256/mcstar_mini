subroutine scatmat(ipol,phsf,ds,drs,out,eps)
  use paras, only : kpol2, rad, kpol
  implicit none
  integer, intent(in) :: ipol
  real(8), intent(in) :: phsf(kpol2)
  real(8), intent(in) :: ds(3), drs(3)
  real(8), intent(out) :: out(kpol,kpol)
  real(8), intent(in) :: eps
  ! workspace
  real(8) :: fis, firs, rotc1, rotc2, rots1, rots2
  real(8) :: c(kpol,kpol), phmx(kpol,kpol)
  
  out = 0.d0
  if(ipol==1) then
    out(1,1) = phsf(1)
    return
  end if
  phmx(1:ipol,1:ipol) = 0.d0
  phmx(1,1) = phsf(1)
  phmx(1,2) = phsf(5)
  phmx(2,1) = phsf(5)
  phmx(2,2) = phsf(2)
  phmx(3,3) = phsf(3)
  if(ipol==4) then
    phmx(4,4) = phsf(4)
    phmx(3,4) = phsf(6)
    phmx(4,3) = -phsf(6)
  end if
  
  fis = sqrt(sum(ds(1:2)**2.))
  if(fis.lt.eps) then
    fis = 0.d0
  else
    fis = sign(acos(ds(1)/fis)/rad,ds(2))
  end if
  
  firs = sqrt(sum(drs(1:2)**2.))
  if(firs.lt.eps) then
    fis = 0.d0
  else
    fis = fis-sign(acos(drs(1)/firs)/rad,drs(2))
  end if
  
  call rotsp(ds(3),drs(3),fis,rotc1,rots1,rotc2,rots2,eps)
  
  c(1:ipol,1:ipol) = 0.d0
  c(1,1) = 1.d0
  if(ipol==4) c(4,4) = 1.d0
  c(2,2) = rotc2
  c(3,3) = rotc2
  c(2,3) = rots2
  c(3,2) = -rots2
  phmx(1:ipol,1:ipol) = matmul(c(1:ipol,1:ipol),phmx(1:ipol,1:ipol))
  c(2,2) = rotc1
  c(3,3) = rotc1
  c(2,3) = rots1
  c(3,2) = -rots1
  out(1:ipol,1:ipol) = matmul(phmx(1:ipol,1:ipol),c(1:ipol,1:ipol))
  
end subroutine scatmat