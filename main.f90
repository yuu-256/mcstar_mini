module params
    implicit none

    ! global I/O units
    integer :: iuo=77, iuow=79
end module params


program main
    use params
    implicit none

    !!! CLOCK START
    call cpu_time(time_beg)

    !!! COLD START
    call ColdStart() ! done

    !!! CALCULATE PRECISION SET
    call SetPrecision(eps) ! done

    !!! OPEN CONFIGURATION FILE
    call OpenConfig(iui) ! done

    !!! READ CONFIGURATION FILE
    call  ReadConfig(iui,
    fi,wlmin,wlmax,galbs,galbl,matm,nln,ipbf,ifrh,trh,      &
    npoly,icn,wlcn,ncomp,cnpt,mptc,vptc,icycl,np0,nx,ny,    &
    dx,dy,fis,indg,jatm,cconc,npv,llr,llu,rpvw)

    !!! OPEN DATA FILES
    call OpenData(icont,                                    &                       
    iua,ius,iuk,iup,iukd,iug,iuo,iuow)

    !!! READ DATA FILES
    call ReadData(iua,ius,iug,icont,wlmin,wlmax,            &
    natm,nl,alt,prs,tmp,nmol,cng,nptc,ins,cnp,ispcvp,       &
    rfracp,asphr,rop,dryap,nawp,awcrp,nv,nwlv,wlv,rfi,      &
    nbnd,wbnd)

    !!! CLOSE DATA FILES
    !!! ON DEVELOPMENT !!!
    ! call CloseData(iui, iua, ius, iuk, iup, iukd, iug)

    !!! PREPROCESSING
    call Preprocess

    !!! DISPLAY SETTINGS AND DATA
    call DisplaySettings(icont,nfi,fi,kback) ! done

    !!! INITIALIZE VARIABLES
    call Initialize(na12, na0, nfi, npv, ipol, nx, ny, sol, &
    ais2, aie2, flxsp2, flxep2, flxs2, flxs02,              &   
    flxuxy2, flxd002) ! done

    !!! MONTE CARLO SIMULATION INCLUDING CKD BANDS INTEGRATION
    call MCST
    
    !!! WRITE OUTPUT FILES
    call WriteOutput !ON DEVELOPMENT!!!

    !!! CLOSE OUTPUT FILES
    call CloseOutput

    !!! DISPLAY OUTPUT
    call DisplayOutput() ! これいる？

    !!! CLOCK STOP
    call cpu_time(time_end)
    print *, 'Elapsed time: ', time_end-time_beg, ' seconds'

end program main
