subroutine MCST()

    !!! DECLARATION:
    implicit none

    !!! EXECUTION:
  ! CKD band integration
    do iw=1,nbnd
 
        if(wl(iw)<wlmin .or. wl(iw)>wlmax) cycle

        ag(1:nx,1:ny)=galb(1:nx,1:ny,iw)
        
        if(icont(2)==1) ngas(1)=iw

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

        call IncrementBackwards(na12,na0,nfi,npv,ipol,inda,fsol,dw(iw),sol,ais1,aie1,flxsp1,flxep1,ais2,aie2,flxsp2,flxep2)

        call IncrementForwards(na0, nx, ny, dw(iw), flxs01, flxs1, flxe1, flxuxy1, flxd001, flxs02, flxs2, flxe2, flxuxy2, flxd002)

    end do

end subroutine MCST
