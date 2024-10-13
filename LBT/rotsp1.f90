subroutine rotsp(amui,amuj,fi,rotc1,rots1,rotc2,rots2,eps)
!------ history
! 94. 6.25 created
!     6.29 modified
! 95. 5.30 kndi1->kndi,etc
!     8.14 rotc=1-2g**2, rots=sqrt(1-rotc**2)
!          if(abs(amuj(j)).eq.1 .and. abs(amui(i)).eq.1) then....
! 21. 7. 5  introduction of SF for negative fi as ps4
!mm231003 from rot1 to single way
!------ memo
! p=c(-kai)*p(thi)*c(kai')
! thi - emerget direction, thi' - incident direction
!  f(thi,thi')=(cos(thi)-cos(thi')*cos(thi))/sin(thi')/sin(thi)
!  g(thi,thi')=sin(thi)*sin(fi-fi')/sin(thi)
! rotc1(thi,thi')=2f(thi,thi')**2-1=1-2g(thi,thi')**2
! rots1(thi,thi')=2f(thi,thi')*g(thi,thi')=sqrt(1-rotc1(thi,thi')**2)
! rotc2(thi,thi')=rotc1(thi',thi)
! rots2(thi,thi')=rots1(thi',thi)
!
! c(kai')=|1                                   |
!         |   rotc1(thi,thi') rots1(thi,thi')  |
!         |  -rots1(thi,thi') rotc1(thi,thi')  |
!         |                                   1|
!
! c(-kai)=|1                                   |
!         |   rotc2(thi,thi') rots2(thi,thi')  |
!         |  -rots2(thi,thi') rotc2(thi,thi')  |
!         |                                   1|
!
!------ input
! amui  r     cosine of the zenith angle
! amuj  r     cosine of the zenith angle
! nfi       i       number of azimuthal angles
! fi    r(knfi)     azimuthal angle in degree.
!------ output
! cs       r(kndi,kndj,knfi)     scattring angle
! rotc1    r(kndi,kndj,knfi)
! rots1    r(kndi,kndj,knfi)
! rotc2    r(kndi,kndj,knfi)
! rots2    r(kndi,kndj,knfi)
!--
  use paras, only: rad
  implicit none
  real(8),intent(in):: amui
  real(8),intent(in):: amuj
  real(8),intent(in):: fi

  real(8),intent(out):: rotc1
  real(8),intent(out):: rots1
  real(8),intent(out):: rotc2
  real(8),intent(out):: rots2

! work
  integer:: ir,i,j,k
  real(8):: cs1,cs2,cf,g1,g2,ccf,ccs,ss,f1,f2
  real(8):: smui,smuj
!tr210705
  real(8):: sf
  real(8):: eps
!--exec

  smui=sqrt(1.d0-amui**2)
  smuj=sqrt(1.d0-amuj**2)
  cs1=amui*amuj
  cs2=smui*smuj
        
  cf=cos(fi*rad)
!tr210705
  sf=sin(fi*rad)
  if(abs(cf) < eps) cf=0.d0
  ccf=cf

  if(abs(cs1)<eps) then
  !if(abs(cs1)<eps) then
  !if(abs(cs1)==1.d0) then
     rotc1=1.d0; rots1=0.d0
     rotc2=1.d0; rots2=0.d0
  else
     if(abs(amui-amuj)<eps*100. .and. 1.d0-cf<eps) then
       ccf=cos((180.01)*rad)
       ! cf==1.d0; sin(thi-thj)
       !ccs = cs1+cs2 = amui*amuj + smui*smuj = cos(thi-thj)
       !ss = sin(thi-thj)
       !amui*smuj -smui*amuj = -sin(thi-thj)
       !f1 = -sin(thi-thj)/sin(thi-thj)
       !g1 = 0.0
       !f2 = sin(thi-thj)/sin(thi-thj)
       !g2 = 0.0
       rotc1=1.d0; rots1=0.d0
       rotc2=1.d0; rots2=0.d0
     else if(abs(amui+amuj)<eps*100. .and. 1.d0+cf<eps) then
       ccf=cos((179.99)*rad)
       ! cf==-1.d0; sin(thi-thj)
       !ccs = cs1-cs2 = amui*amuj - smui*smuj = cos(thi+thj)
       !ss = sin(thi+thj)
       !amui*smuj +smui*amuj = sin(thi+thj)
       !f1 = 1.d0
       !g1 = 0.0
       !f2 = 1.d0
       !g2 = 0.0
       rotc1=1.d0; rots1=0.d0
       rotc2=1.d0; rots2=0.d0
     else
       ccs=cs1+cs2*ccf
       ss=sqrt(1.d0-ccs**2)
       if((1.0-amui)<eps .or. (1.0+amui)<eps .or. ss<eps*100.) then
         rotc1=1.d0; rots1=0.d0
       else
         f1=(amui*smuj -smui*amuj*cf)/ss
         g1=smui*sf/ss
         rotc1=2.d0*f1**2-1.d0
         rots1=2.d0*f1*g1
       end if
       if((1.0-amuj)<eps .or. (1.0+amuj)<eps .or. ss<eps*100.) then
         rotc2=1.d0; rots2=0.d0
       else
         f2=(smui*amuj-amui*smuj*cf)/ss
         g2=smuj*sf/ss
         rotc2=2.d0*f2**2-1.d0
         rots2=2.d0*f2*g2
       end if
     end if
  endif
  return
end subroutine rotsp

