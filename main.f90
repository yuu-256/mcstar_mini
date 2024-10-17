module params
    implicit none

    ! global I/O units
    integer :: iuo=77, iuow=79
end module params


program main
    use params
    implicit none

    !!! CLOCK START
    call cpu_time(time_beg)

    !!! OPEN DATA FILES
    call OpenData()

    !!! READ DATA FILES
    call ReadData(iui, iua, iug, iuk, ius, iup, iukd)

    !!! READ CONFIGURATION FILE
    call ReadConfig(iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow)

    !!! CLOSE DATA FILES
    call CloseData

    !!! PREPROCESSING
    call Preprocess

    !!! DISPLAY SETTINGS AND DATA
    call DisplaySettings

    !!! INITIALIZE VARIABLES
    call Initialize

    !!! GET 3D STRUCTURE OF THE SYSTEM
    call GetStructure

    !!! MONTE CARLO SIMULATION INCLUDING CKD BANDS INTEGRATION
    call MCST
    
    !!! WRITE OUTPUT FILES
    call WriteOutput

    !!! CLOSE OUTPUT FILES
    call CloseOutput

    !!! DISPLAY OUTPUT
    call DisplayOutput

    !!! CLOCK STOP
    call cpu_time(time_end)
    print *, 'Elapsed time: ', time_end-time_beg, ' seconds'

end program main
