subroutine OpenData(iui,iua,ius,iuk,iup,iukd,iuo,iuow)
    !!! OPEN DATA FILES

    !!! DESCRIPTION:
    implicit none
    integer, intent(out) :: iui,iua,ius,iuk,iup,iukd,iuo,iuow
    character(len=80) :: fnm 

    ! data file i/o units
    iui= 1    ! user data file
    iua= 3    ! atmospheric model library read
    iug= 7    ! CKD gas absorption table
    iuk=11    ! Mie kernel file
    ius=12    ! Aerosol library
    
    !!! EXECUTION:
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

    open (iuo,file='out' )           ! standard output file
    open(iuow,file='outw')           ! work file for photon info

end subroutine OpenData
