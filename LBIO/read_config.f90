subroutine read_config(iui,icont,isol,inda,indg,imthd,ipol,nda,nds,na0,th0,na1,th1)
    ! Read configuration file
    ! IN : iui : mcdata file number

    implicit none
    integer, intent(in) :: iui
    integer, intent(out) :: icont,isol,inda,imthd,ipol,nda,nds,na0,na1
    real(8), intent(out) :: th0(0:),th1(0:),indg(0:0)

    ! Read experimenal parameters
    read(iui,*,err=1,end=1) icont(1:4)
    read(iui,*) isol,inda,indg(1,1),imthd,ipol,nda,nds
    read(iui,*) na0, (th0(i),i=1,na0)
    read(iui,*) na1,(th1(i),i=1,na1)

end subroutine read_config
