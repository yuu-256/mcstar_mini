subroutine CloseData(iui, iua, ius, iuk, iup, iukd)

    !!! DECLARATION:
    implicit none
    integer, intent(in) :: iui, iua, ius, iuk, iup, iukd

    !!! EXECUTION:
    close(iui)
    close(iua)
    close(ius)
    close(iuk)
    close(iup)
    close(iukd)

end subroutine CloseData
