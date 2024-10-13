!mm220221
subroutine imsbxb1(init,ipol,icycl,dseed,pbnd,drs,dvw,pvw,nb,bx,dbx, &
    ce,omg,phs,ph2,angdp,uh,utrn0,lim1,mlim,ind)
! 3d-ims correction in box system.
! ray tracing til mlim1
! inverse monte-carlo. one photon case.
!--- reference
! M.Momoi et al., in prep.
!--- history
!mm220109 created pims from imsbxb1
! proc is flow in MCims_doc (M. Momoi)
!mm220131 invalid flag
!mm220213 pims => ims M.Momoi
!--- parameter
! kx,ky,kz            declared numbers of boxes in x/y/z directions
! kb1                 kb1=max(kx, ky, kz) + 1
! knang               declared number of scattering angles

  use paras
  implicit none
  character(len=64) :: erc
  real(8), save ::  eps,u1min,da
  integer, save :: mlim1,mlim2
!  double precision dseed
  real(8) :: dseed, pbnd
  real(8) :: drs(3),dvw(3),pvw(3),bx(kb1,3),dbx(kb1,3,2)
  integer :: nb(3),mlim(*)
!c working area
  real(8) :: ps(3),pe(3),ds(3),pe2(3), &
       cpc(3)
  integer :: ibs(3),ibe2(3), ibe(3)
  integer :: init,ind,ibx,iby,ibz,i, &
       icycl,lim1,lim2,is
  real(8) :: ut(kpol,kpol),uh(kpol),u1(kpol),u2(kpol),u3(kpol),trne,tau, &
       cst,cet,dst,taue,phs1,dtau,w0t, taus, dtaus
  real(8) :: angdp(na4,kx,ky,kz),phs(na4,kx,ky,kz,kpol2),ph2(na4,kx,ky,kz)
  integer :: nrnd,nx,ny,nz
  real(8) :: taus1,trns1,ce(kx,ky,kz),omg(kx,ky,kz),scat,sfi
  real(8) :: rndv(2)
  integer:: ind1, jpol
  real(8) :: tr11, tr12, tr13, tr15, trh
  !mm231005
  real(8) :: utrn(kpol),utrn0(kpol),phmx(kpol,kpol), ds2(3)
  integer :: ipol

! initialization
  if(init.ne.0) then
     init=0
     call cpcon(cpc)
     eps=cpc(1)*100.0
     u1min=10.0**(cpc(2)*0.8)
     mlim1=mlim(1)
     mlim2=mlim(2)
     da=180.0/(na4-1)
  endif

! start
  ind=0
  nx=nb(1); ny=nb(2); nz=nb(3)
  ut(1:ipol,1:ipol) = 0.0
  do jpol = 1, ipol
    ut(jpol,jpol) = 1.
  end do
  uh(1:ipol)=0
  utrn(1:ipol)=0.
  utrn0(1:ipol)=0.
! backward tracing from the receiver
  ps(1:3)=pvw(1:3)
  ds(1:3)=dvw(1:3)
! find box, ibs: found position
  call fbox(ps,bx,nb,ibs,erc)
  if (erc.ne.' ') then
    ind=9
    return
  endif
  
! extinction of the ray
  lim1=0
  do10 : do
    lim1=lim1+1
    if(lim1.gt.mlim1) then
      !ind=8
      return
    endif

    ! transmittance
    ! proc-1
    call gttau(dseed,trne,taue)
    if(lim1.eq.1) then
      trne = 1. - (1.-trne)*(1.-pbnd)
      taue = - log(trne)
    end if
    if(trne.le.0) return
    ! ray tracing along the incident direction
    ! ray trace in a box (number of intercepting walls)
    call raytrace2(icycl,nb,bx,dbx,ds,ps,ibs,pe2,ibe2,&
        ce,omg,taue,tau,dtau,taus,dtaus,eps,mlim2,ind)
    if( taue > tau + dtau .or. iabs(ind).gt.3 ) return
    ! proc-1A
    ps(1:3) = pe2(1:3)
    ibs(1:3) = ibe2(1:3)
    ibx = ibs(1); iby = ibs(2); ibz = ibs(3)
    if( omg(ibx,iby,ibz) .eq. 0. ) return
    !! scattering point
    dst=(taue-tau)/ce(ibx,iby,ibz)
    ps(1:3)=dst*ds(1:3)+ps(1:3)
    trh=taus+omg(ibx,iby,ibz)*(taue-tau)
    trh=exp(-trh)
    ! proc-1B
    !! scattering angle
    i=acos(drs(1)*ds(1)+drs(2)*ds(2)+drs(3)*ds(3))/rad/da+1.5
    if(i<1 .or. i>na4) then
      ind=9
      return
    endif
    if(ipol==1) then
      u1(1)=omg(ibx,iby,ibz)*ut(1,1)*phs(i,ibx,iby,ibz,1)
    else
      call scatmat(ipol,phs(i,ibx,iby,ibz,1:kpol2),ds,drs,phmx,eps)
      phmx(1:ipol,1:ipol) = matmul(ut(1:ipol,1:ipol),phmx(1:ipol,1:ipol))
      u1(1:ipol)=omg(ibx,iby,ibz)*phmx(1:ipol,1)
    end if
    ! proc-1C
    ut(1:ipol,1:ipol)=ut(1:ipol,1:ipol)*omg(ibx,iby,ibz)*trh
    ! proc-1D
    call gttau(dseed,tr11,taue)
    ! ray tracing along the solar direction
    ! ray trace in a box (number of intercepting walls)
    call raytrace2(icycl,nb,bx,dbx,drs,ps,ibs,pe2,ibe2,&
        ce,omg,taue,tau,dtau,taus,dtaus,eps,mlim2,ind1)
    if( taue <= tau + dtau .and. iabs(ind1).le.3 ) then
      pe(1:3) = pe2(1:3)
      ibe(1:3) = ibe2(1:3)
      ibx = ibe(1); iby = ibe(2); ibz = ibe(3)
      ! proc-1E
      u2(1:ipol)=-u1(1:ipol)*trh*omg(ibx,iby,ibz) !/2/pi
      ! proc-1D: trh11
      !! scattering point
      dst=(taue-tau)/ce(ibx,iby,ibz)
      pe(1:3)=dst*drs(1:3)+pe(1:3)
      trh=taus+omg(ibx,iby,ibz)*(taue-tau)
      trh=exp(-trh)
      ! proc-1F
      call gttau(dseed,tr12,taue)
      ! ray tracing along the solar direction
      ! ray trace in a box (number of intercepting walls)
      call raytrace2(icycl,nb,bx,dbx,drs,pe,ibe,pe2,ibe2,&
          ce,omg,taue,tau,dtau,taus,dtaus,eps,mlim2,ind1)
      if( taue <= tau + dtau .and. iabs(ind1).le.3 ) then
        pe(1:3) = pe2(1:3)
        ibe(1:3) = ibe2(1:3)
        ibx = ibe(1); iby = ibe(2); ibz = ibe(3)
        !! scattering point
        dst=(taue-tau)/ce(ibx,iby,ibz)
        pe(1:3)=dst*drs(1:3)+pe(1:3)
        call strnm(icycl,drs,ibe,pe,ibe2,pe2,nb,bx,dbx, &
            ce,taue,tr13,eps,mlim2,ind1)
        ! proc-1G
        u3(1:ipol)=-u2(1:ipol)*trh*omg(ibx,iby,ibz) !/2/pi
        ! proc-1H
        if(knni-lim1.ge.2) uh(1:ipol)=uh(1:ipol)+u3(1:ipol)*tr13
        if(knni-lim1.ge.3) utrn(1:ipol)=utrn(1:ipol)+u3(1:ipol)*tr13
      else if(iabs(ind1).le.3) then
        taue = tau+dtau
        tr12=exp(-taue)
        tr13=1.
      else
        tr13=0.
      end if
      ! proc-1I
      if(knni-lim1.ge.1) uh(1:ipol)=uh(1:ipol)+u2(1:ipol)*tr12*tr13
      if(knni-lim1.ge.2) utrn(1:ipol)=utrn(1:ipol)+u2(1:ipol)*tr12*tr13
    else if(iabs(ind1).le.3) then
      taue = tau+dtau
      tr11=exp(-taue)
      tr12=1.; tr13=1.
    else
      tr12=0.; tr13=0.
    end if
    ! proc-1J
    uh(1:ipol)=uh(1:ipol)+u1(1:ipol)*tr11*tr12*tr13
    if(knni-lim1.ge.1) utrn(1:ipol)=utrn(1:ipol)+u1(1:ipol)*tr11*tr12*tr13
    
    !call strnm(icycl,drs,ibs,ps,ibe2,pe2,nb,bx, &
    !        ce,taue,trne,eps,mlim2,ind1)
    !uh=u1*tr13
    !if( tr11*tr12*tr13.ne.trne ) then
    !  write(*,*) tr11*tr12*tr13, tr11,tr12,tr13,trne
    !  stop
    !end if
    
    !if(mlim1-lim1.lt.1 .or. ut.le.u1min) return
    if(knni-lim1.lt.1 .or. abs(ut(1,1)).le.u1min) then
      utrn0(1:ipol)=uh(1:ipol)-utrn(1:ipol)
      return
    end if
    
    ! proc-2
    call gttau(dseed,trne,taue)
    if(trne.gt.0) then
      ! ray tracing along the incident direction
      ! ray trace in a box (number of intercepting walls)
      call raytrace2(icycl,nb,bx,dbx,ds,ps,ibs,pe2,ibe2,&
          ce,omg,taue,tau,dtau,taus,dtaus,eps,mlim2,ind1)
    end if
    if( taue <= tau + dtau .and. trne.gt.0 .and. &
      omg(ibe2(1),ibe2(2),ibe2(3)).ne.0. .and. iabs(ind1).le.3) then
      ! proc-2A
      pe(1:3) = pe2(1:3)
      ibe(1:3) = ibe2(1:3)
      ibx = ibe(1); iby = ibe(2); ibz = ibe(3)
      !! scattering point
      dst=(taue-tau)/ce(ibx,iby,ibz)
      pe(1:3)=dst*ds(1:3)+pe(1:3)
      call strnm(icycl,drs,ibe,pe,ibe2,pe2,nb,bx,dbx, &
          ce,taue,tr15,eps,mlim2,ind1)
      ! proc-2B
      !! scattering angle
      i=acos(drs(1)*ds(1)+drs(2)*ds(2)+drs(3)*ds(3))/rad/da+1.5
      if(i<1 .or. i>na4) then
        ind=9
        return
      endif
      if(ipol==1) then
        u1(1:ipol)=-omg(ibx,iby,ibz)*ut(1,1)*phs(i,ibx,iby,ibz,1)
      else
        call scatmat(ipol,phs(i,ibx,iby,ibz,1:kpol2),ds,drs,phmx,eps)
        phmx(1:ipol,1:ipol) = matmul(ut(1:ipol,1:ipol),phmx(1:ipol,1:ipol))
        u1(1:ipol)=-omg(ibx,iby,ibz)*phmx(1:ipol,1)
      end if
      ! proc-2C
      uh(1:ipol)=uh(1:ipol)+u1(1:ipol)*tr15 !/2/pi
      if(knni-lim1.eq.1) utrn(1:ipol)=utrn(1:ipol)+u1(1:ipol)*tr15
    end if
    
    ! proc-3
    ibx = ibs(1); iby = ibs(2); ibz = ibs(3) !mm220126
    nrnd=2
    call rdmu2(dseed,nrnd,rndv)
    !mm220118 interpolartion
    !i=rndv(1)*(na4-1)+1.5
    !if(i<1 .or. i>na4) then
    !  ind=9
    !  return
    !endif
    !scat=angdp(i,ibx,iby,ibz)
    call angintp(rndv(1),na4,angdp(1:na4,ibx,iby,ibz),ind,scat)
    if(ind.eq.9) return
    !end interpolation
    sfi=rndv(2)*360.0
    ds2(1:3) = ds(1:3)
    call atrns(scat,sfi,ds,eps)
    
    i=scat/da+1.5
    if( i < 0 .or. i > na4 ) then
      ind=9
      return
    end if
    if(ipol==1) then
      ut(1,1) = ut(1,1) * phs(i,ibx,iby,ibz,1) / ph2(i,ibx,iby,ibz)
    else
      call scatmat(ipol,phs(i,ibx,iby,ibz,1:kpol2)/ph2(i,ibx,iby,ibz),ds2,ds,phmx,eps)
      ut(1:ipol,1:ipol) = matmul(ut(1:ipol,1:ipol),phmx(1:ipol,1:ipol))
    end if
    if( abs(phs(i,ibx,iby,ibz,1)).le.eps .or. abs(ut(1,1)).le.u1min ) return
  end do do10
  
  return
end subroutine imsbxb1
