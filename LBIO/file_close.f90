subroutine file_close(iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow)

    implicit none
    integer, intent(in) :: iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow

    ! Close all open files
    close(iui)
    close(iua)
    close(iug)
    close(iuk)
    close(ius)
    close(iup)
    close(iukd)
    close(iuo)
    close(iuow)

end subroutine file_close
