subroutine MCST()
    !!! DECLARATION:
    implicit none

    !!! EXECUTION:
  ! CKD band integration
    do iw=1,nbnd

      if(wl(iw)<wlmin .or. wl(iw)>wlmax) cycle
  ! ground albedo horizontal distributions
      ag(1:nx,1:ny)=galb(1:nx,1:ny,iw)
      if(icont(2)==1) ngas(1)=iw

      !mm220221 added iup,iukd,ipol, ins/ deleted dw, galb
      call mcst2(init0,icont,iuk,iup,iukd,iug,indg,inda,isol,ipol,imthd,nda,&
         fis,na0,th0,na1,th1,nfi,fi,npv,llu,llr,rpvw,ngas,wl(iw), &
         icn,wlcn,npoly,ncomp,mptc,vptc,cnpt, &
         nl,alt,prs,tmp,gtmp,cng,ins,cnp,ispcvp,rfracp,asphr, &
         rop,dryap,nawp,awcrp,nwlv,wlv,rfi,nln,ipbf, &
         icycl,dseed,np0,nb,bx,cconc,ag,jatm, &
         fsol,tauav,tausav,taurz,tckdav,gav,flxs1,flxs01,flxe1,flxuxy1, &
         flxd001,ais1,aie1,flxsp1,flxep1,nsos,err)

       if(err/=' ') then
         write(*,*) 'error in mcst1'
         cycle ! iw
       endif

      call IncrementBackwards()

      call IncrementForwards()

    end do

end subroutine MCST
