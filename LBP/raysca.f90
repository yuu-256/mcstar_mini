subroutine raysca(nln,wl,co2,dpf,taur)
! get Rayleigh scattering coefficient
  use paras, only: knln,knln1,kpol2,knang,pi
  implicit none
  integer,intent(in):: nln
  real(8),intent(in):: wl
  real(8),intent(in):: co2(knln)
  real(8),intent(out):: dpf(knln)
  real(8),intent(out):: taur(knln)
  
  integer:: l
  real(8):: fn2,fo2,ri300,wlc,co2vm,fair,air,fk,riair,c1,c2

  real(8),parameter:: vn2=78.084d0, vo2=20.946d0, var=0.934d0
  real(8),parameter:: dm = 2.546899d+19  ! [1/cm3]
!--exec
!! Bodhaine et al. 1999
  fn2=1.034d0+3.17d0*1.d-4/wl**2
  fo2=1.096d0+1.385d0*1.d-3/wl**2+1.448d0*1.d-4/wl**4

  if(wl > 0.23d0) then
     ri300 = 5791817.0d0/(238.0185d0 - 1.0d0/wl**2)  &
            + 167909.0d0/(57.362d0 - 1.0d0/wl**2)
  else
     ri300 = 8060.51d0 + 2480990.0d0/(132.274d0 - 1.0d0/WL**2)  &
            + 17455.7d0/(39.32957d0 - 1.0d0/WL**2)
  endif
  ri300 = ri300 * 1.d-8
  wlc = wl * 1.d-4

  do l=1,nln
     co2vm=co2(l)*1.d-4
     fair = vn2*fn2 + vo2*fo2 + var + co2vm*1.15d0
     air = vn2 + vo2 + var + co2vm
     fk=fair/air

     dpf(l)=6.d0*(fk - 1.d0)/(7.d0*fk + 3.d0)
     
     riair = 1.d0 + 0.54d0 * ( co2(l)*1.d-6 - 0.0003d0 )
     riair = riair * ri300 + 1.d0

     c1 = 24.d0 * pi**3 / (wlc**4 * dm**2)
     c2 = (riair**2 - 1.d0)**2 / (riair**2 + 2.d0)**2
     taur(l) = c1 * c2 * fk *dm
!     xsec = c1 * c2 * fk
!     taur(l)= xsec * dm * (pl(l)/pstd) * (tstd/tl(l))
  enddo
  end subroutine raysca
