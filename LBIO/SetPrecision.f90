subroutine SetPrecision(eps)

    !!! DECLARATION:
    real(8) :: eps
    real(8) :: cpc(3)

    !!! EXECUTION:    
    call cpcon(cpc)
    eps=cpc(1)*10

end subroutine SetPrecision
