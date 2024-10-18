subroutine CloseData(iui, iua, ius, iuk, iup, iukd, iug)

    !!! DECLARATION:
    implicit none
    integer, intent(in) :: iui, iua, ius, iuk, iup, iukd, iug

    !!! EXECUTION:
    close(iui)
    close(iua)
    close(ius)
    close(iuk)
    close(iup)
    close(iukd)
    close(iug)

end subroutine CloseData
