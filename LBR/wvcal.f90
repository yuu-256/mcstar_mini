subroutine wvcal(ind,p,t,rh,ppmv,gm3,e,gm3s,es)
! conversion of units for water vapor
! ppmv =n/na =e/(p-e)*1.0e6, rh =e/es =gm3/gm3s
! ev = nrt, (p-e)v = na rt
!--- history
! 95. 9.30  created
! 08. 9.20  modified Fortran 90 free-style form
!tr210105 ind out of range action 
!--- input
! ind   i    1: give rh  (relative humidity, 0-1)
!            2: give ppmv (volume mixing ratio in ppmv)
!            3: give gm3 (mass mixing ratio in g/m3)
!            4: give e    (water vapor pressure in hpa)
! p     r    atmospheric pressure (hPa)
! t     r    temperature (k)
! one of rh, ppmv, gm3, or e
!--- output
! other than input variable among rh, ppmv and gm3.
! gm3s  r    saturation water vapor content (g/m3)
! es    r    saturation water vapor pressure (hPa)
!--- gas constant (cgs) and molecular weights of water vapor
!           and air (g/mol)
!--
  implicit none
! for input & output
  integer,intent(in)::ind        !! indicator
  real(8),intent(in)::p,t        !! pressure and temperature
  real(8),intent(inout)::rh,ppmv,gm3,e,gm3s,es !! quantities of watere vaour

! work
  real(8):: gcm3,gcm3s
  real(8):: wvsat                !! function for saturated
  real(8),parameter::r=8.314d7, w=18.02d0, air=28.964d0
!--exec
! satuaration water vapor (g/cm3)
  gcm3s=wvsat(t)*1.0d-6
  es=gcm3s/w*r*t/1.0d3
  if(ind==1) then
     gcm3=rh*gcm3s
     e=rh*es
     ppmv=e/(p-e)*1.0d6
  elseif(ind==2) then
     e=p*ppmv/(1.0d6+ppmv)
     rh=e/es
     gcm3=rh*gcm3s
  elseif(ind==3) then
     gcm3=gm3*1.0d-6
     rh=gcm3/gcm3s
     e=rh*es
     ppmv=e/(p-e)*1.0d6
  elseif(ind==4) then
     rh=e/es
     gcm3=rh*gcm3s
     ppmv=e/(p-e)*1.0d6
!tr210105
  else
    rh=9.9e9
    gcm3=9.9e9
    gcm3s=9.9e9
    return
  endif
  gm3=gcm3*1.0d6
  gm3s=gcm3s*1.0d6
  if(rh>1.d0) rh=1.d0
  return
end subroutine wvcal

function wvsat(t)
! saturation water vapor amount (LOWTRAN)
!  ws=a*exp(18.9766-14.9595a-2.4388a**2) g/m3
!  a=t0/t, t0=273.15, t in k (t0-50,t0+50)
!--- history
! 93. 3. 3 created
! 08. 9.20  modified Fortran 90 free-style form
!--- input
! t       r      absolute temperature (K)
!                t0-50 to t0+50
!--- output
! wvsat  rf      water vapor amount (g/m3)
!---
  implicit none
  real(8),intent(in)::t
  real(8)::a,wvsat
!--exec
  a=273.15d0/t
  wvsat=a*dexp(18.9766d0-14.9595d0*a-2.4388d0*a**2)
  return
end function wvsat
