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

    !!! OPEN DATA FILES
    call OpenData(iui, iua, ius, iuk, iup, iukd, iuo, iuow) ! done

    !!! READ DATA FILES
    call ReadData(iui, iua, iug, iuk, ius, iup, iukd)

    !!! READ CONFIGURATION FILE
    call  ReadConfig(iui,fi,wlmin,wlmax,galbs,galbl,matm,   &
    nln,ipbf,ifrh,trh,npoly,icn,wlcn,ncomp,cnpt,mptc,vptc,  &
    icycl,np0,nx,ny,dx,dy,fis,indg,jatm,cconc,npv,llr,llu,  &
    rpvw)
    
    do lpv=1,npv
      nrnd=2
      call rdmu2(dseed,nrnd,rndv)
      rpvw(lpv,1:2)=rpvw(lpv,1:2)*(1-100*eps*rndv(1:2))
    enddo

    !!! CLOSE DATA FILES
    call CloseData

    !!! PREPROCESSING
    call Preprocess

    !!! DISPLAY SETTINGS AND DATA
    call DisplaySettings

    !!! INITIALIZE VARIABLES
    call Initialize

    !!! GET 3D STRUCTURE OF THE SYSTEM
    call GetStructure

    !!! MONTE CARLO SIMULATION INCLUDING CKD BANDS INTEGRATION
    call MCST
    
    !!! WRITE OUTPUT FILES
    call WriteOutput

    !!! CLOSE OUTPUT FILES
    call CloseOutput

    !!! DISPLAY OUTPUT
    call DisplayOutput

    !!! CLOCK STOP
    call cpu_time(time_end)
    print *, 'Elapsed time: ', time_end-time_beg, ' seconds'

end program main
