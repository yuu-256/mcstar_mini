subroutine OpenData(icont,iua,ius,iuk,iup,iukd,iug,iuo,iuow)
    !!! OPEN DATA FILES

    !!! DECLARATION:
    use params
    implicit none
    integer, intent(in) :: icont(10)
    integer, intent(out) :: iua,ius,iuk,iup,iukd,iug,iuo,iuow
    character(len=80) :: fnm
    integer :: ngas(2)

    ! data file i/o units
    iua= 3    ! atmospheric model library read
    iug= 7 ! CKD gas absorption table for narrow band
    iuk=11    ! Mie kernel file
    ius=12    ! Aerosol library
    iup=13    ! photon scattering phase function
    iukd=14

    ! work file i/o units
    iuo=77
    iuow=79
    

    !!! EXECUTION:
    write(*,*) 'Open data files'

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

    if(icont(2)<=1) then
       fnm='01MCD1/PARAG.ch111'         ! CKD gas absorption table
       open(iug,file=fnm,status='old')
    else if(icont(2)==2) then
       ngas(1:2)=(/2,2000/)
       fnm='01MCD1/ckd.g.ch_2_2e3_ltl'     ! CKD gas absorption table
       open(iug,file=fnm,status='old',access='direct',form='unformatted', &
            recl=4*kp*kt*ngas(1)*kmol1)
    end if

    open (iuo,file='out' )           ! standard output file
    open(iuow,file='outw')           ! work file for photon info

    write(*,*) 'Open data files done'

end subroutine OpenData
