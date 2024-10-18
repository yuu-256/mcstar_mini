subroutine OpenConfig(iui)

    !!! DECLARATION:
    implicit none
    integer, intent(out) :: iui
    character(len=80) :: fnm
    iui=1

    !!! EXECUTION:
    fnm='mcdata'
    open (iui,file=fnm)

end subroutine OpenConfig
