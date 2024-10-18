subroutine OpenData(iua,ius,iuk,iup,iukd,iug_n,iug_b,iuo,iuow)
    !!! OPEN DATA FILES

    !!! DESCRIPTION:
    use params
    implicit none
    integer, intent(out) :: iua,ius,iuk,iup,iukd,iug_n,iug_b,iuo,iuow
    character(len=80) :: fnm 

    ! data file i/o units
    iua= 3    ! atmospheric model library read
    iug_n= 8 ! CKD gas absorption table for narrow band
    iug_b= 9 ! CKD gas absorption table for broad band
    iuk=11    ! Mie kernel file
    ius=12    ! Aerosol library

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
    fnm='01MCD1/PARAG.ch111'         ! CKD gas absorption table
    open(iug_b,file=fnm,status='old')
    fnm='01MCD1/ckd.g.ch_2_2e3_ltl'         ! CKD gas absorption table
    open(iug_n,file=fnm,status='old',access='direct',form='unformatted', &
         recl=4*kp*kt*ngas(1)*kmol1)
    open (iuo,file='out' )           ! standard output file
    open(iuow,file='outw')           ! work file for photon info

    write(*,*) 'Open data files done'

end subroutine OpenData
