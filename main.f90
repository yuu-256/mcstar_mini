module params
    implicit none

    ! global I/O units
    integer :: iuo=77, iuow=79
end module params


program main
    use params
    implicit none

    ! 
    ! file name
    integer :: iui, iua, iug, iuk, ius
    integer :: iup, iukd
    character(len=100) :: fnm
    character(len=64) :: err

    ! data file I/O units
    iui = 1      ! user data file
    iua = 3      ! atmospheric model library read
    iug = 7      ! CKD gas absorption table
    iuk = 11     ! Mie kernel file
    ius = 12     ! Aerosol library

    ! Open data files
    call file_open(iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow)

    ! Read configuration file
    call read_config(iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow)

    ! Close data files
    call file_close(iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow)

    

end program main
