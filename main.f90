module params
    implicit none

    ! global I/O units
    integer :: iuo=77, iuow=79
end module params


program main
    use params
    implicit none

    !!! OPEN DATA FILES
    call OpenData

    !!! READ DATA FILES
    call ReadData(iui, iua, iug, iuk, ius, iup, iukd)

    !!! READ CONFIGURATION FILE
    call ReadConfig(iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow)

    !!! CLOSE DATA FILES
    call CloseData(iui, iua, iug, iuk, ius, iup, iukd, iuo, iuow)

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
    
    !!! DISPLAY OUTPUT
    call DisplayOutput

    !!! CLOSE OUTPUT FILES
    call CloseOutput

    !!! DISPLAY OUTPUT
    call DisplayOutput

end program main
