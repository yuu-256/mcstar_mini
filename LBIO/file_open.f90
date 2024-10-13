subroutine file_open(iui,iua,ius,iuk,iup,iukd,iuo,iuow)
    
    implicit none
    integer, intent(in) :: iui,iua,ius,iuk,iup,iukd,iuo,iuow

    ! file open
    fnm='mcdata'
    open (iui,file=fnm)
    fnm='01MCD1/MLATMD'              ! atmospheric library
    open(iua,file=fnm,status='old')
    fnm='01MCD1/AERDB7_mc3'           ! aerosol library
    open(ius,file=fnm,status='old')
    fnm='01MCD1/PKRNL.OUT'            ! Mie kernel
    open(iuk,file=fnm,status='old')
    fnm='01MCD1/py13phsf_solid_column.dat'
    open(iup,file=fnm,status='old',access='direct',form='unformatted', &
          recl=4*knang2*kpol2)
    fnm='01MCD1/dkrnl/PKRNL.OUT_asp207'
    open(iukd,file=fnm,status='old',form='unformatted')
 
    open(iuo,file='out')
    open(iuow,file='outw')

end subroutine file_open
