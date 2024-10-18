subroutine Preprocess()

    !!! DECLARATION:

    !!! EXECUTION:
    ! Angles settings 
    call AngleSettings()

    ! Viewing positions settings
    call ViewSettings()

    !!! Vertical property settings
    call VerticalPropertySettings()
    
    !!! 3D Structure settings
    call FieldSettings()

end subroutine Preprocess


subroutine AngleSettings(na1,nfi,th1,fi,nfi12)

    !!! DECLARATION:
    use params
    implicit none
    integer,intent(inout) :: na1,nfi
    real(8),intent(inout) :: th1(kna1u),fi(knfi)
    integer,intent(out) :: nfi12

    integer :: i

    !!! EXECUTION:    
    ! emerging zenith angles  i=1:na1 th1<0 (upgoing), na1+1-2*na1 th1>0 (downgoing)
    do i=1,na1
      th1(i+na1)=th1(i)
      th1(i)=180-th1(i)
    enddo

    na1=2*na1

    do i=1,nfi
      fi(i+nfi)=fi(i)+180
    enddo

    nfi=2*nfi
    nfi12=nfi/2

end subroutine AngleSettings


subroutine ViewSettings(npv,eps,dseed,rpvw)

    !!! DECLARATION:
    use params
    implicit none
    integer,intent(in) :: npv
    real(8),intent(in) :: eps,dseed
    real(8),intent(inout) :: rpvw(knpv,2)

    integer :: lpv,nrnd
    real(8) :: rndv(2)

    !!! EXECUTION:
    do lpv=1,npv
      nrnd=2
      call rdmu2(dseed,nrnd,rndv)
      rpvw(lpv,1:2)=rpvw(lpv,1:2)*(1-100*eps*rndv(1:2))
    enddo

end subroutine ViewSettings


subroutine VerticalPropertySettings()

    !!! DECLARATION:
    use params
    implicit none
    integer,intent(in) :: knm1,nptc
    real(8),intent(in) :: galbl(nx,ny),galbs(nx,ny),wbnd(nbnd+1)
    real(8),intent(inout) :: wl(nbnd),dw(nbnd),cnp1(nl,nptc),pmatm(knl,natm),tmatm(knl,natm),amol(knl,knm1,natm),trac(nl,21)
    real(8),intent(inout) :: cnp(nl,nptc)
    real(8),intent(inout) :: gtmp(natm)

    integer :: iw

    !!! EXECUTION:
    !-- set wl and dw from wbnd
    do iw=1,nbnd
       wla=1.d4/wbnd(iw+1)
       wlb=1.d4/wbnd(iw)
       wl(iw)=0.5d0*(wla+wlb)
       dw(iw)=wlb-wla

    !-- setting the ground same albedo for all the band
       if(wbnd(iw)<2500.d0) then
          galb(1:nx,1:ny,iw)=galbl(1:nx,1:ny)
         else
          galb(1:nx,1:ny,iw)=galbs(1:nx,1:ny)
       endif
    enddo

    !-- initialization
    co2ppm = 330.

    prs(1:knl,1:natm)=pmatm(1:knl,1:natm)
    tmp(1:knl,1:natm)=tmatm(1:knl,1:natm)
    cng(1:knl,1:knm1,1:natm)=amol(1:knl,1:knm1,1:natm)
    do iatm=1,natm
      cng(:,8:28,iatm)=trac(:,1:21)
    end do

    nmol=knm1

    ! set vertical profiles of polydispersion
    do iptc=1,nptc 
       cnp(1:nl,iptc)=cnp1(1:nl,iver(iptc))
    enddo
    
  !-- ground temperature given by user
    gtmp(1:natm)=tmp(1,1:natm)

end subroutine VerticalPropertySettings


subroutine FieldSettings()

    !!! DECLARATION:
    use params
    implicit none
    

    !!! EXECUTION:
    if(ifrh==1) then
      do iatm=1,natm
        if( count(jatm(1:nx,1:ny)==iatm)==0 ) cycle
        m=1
        do l=1,nl
          p1=prs(l,iatm)
    ! don't apply trh for pressure smaller than 50.0
          if(p1>50.0) then
             ppmv=cng(l,m,iatm)
             t1=tmp(l,iatm)
             rh=trh
             ind=1
             call wvcal(ind,p1,t1,rh,ppmv,gm3,e,gm3s,es)
             cng(l,m,iatm)=ppmv
          endif
        enddo
      enddo
    endif

    init0=1


    nz=nln
    if(nx>kx .or. ny>ky .or. nz>kz) then
      err='nx>kx .or. ny>ky .or. nz>kz'
      goto 1
    endif

    nb(1)=nx; nb(2)=ny; nb(3)=nz
    dbx(1)=dx; dbx(2)=dy

    do j=1,2
      do i=1,nb(j)+1
      bx(i,j)=(i-1)*dbx(j)   ! km unit
    enddo; enddo

    do lz=1,nz+1
      ltau=nz+2-lz
      bx(lz,3)=alt(iabs(ipbf(ltau)))         ! km unit
    enddo

    nx=nb(1); ny=nb(2); nz=nb(3)
    na12=na1/2

end subroutine FieldSettings
