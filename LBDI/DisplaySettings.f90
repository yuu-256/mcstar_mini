subroutine DisplaySettings(icont,nfi,fi,kback)

    !!! DECLARATION:
    implicit none
    integer, intent(in) :: icont(4)
    integer, intent(in) :: nfi
    real(8), intent(in) :: fi(nfi)
    integer, intent(in) :: kback

    !!! EXECUTION:
    write(*,*) 'Simulation settings:'
    write(*,*) icont(1), ' : on/off Rayleigh (1/0)'
    write(*,*) icont(2), ' : on/off CKD gas absorption (2: narrow/1: broad/0: no)'
    write(*,*) icont(3), ' : on/off 0th Fourier L(0) output (0/1)'
    write(*,*) icont(4), ' : on/off photon info output to outw-file (0/1)'
    write(*,*) 'this routine set paired fi and fi+180 as'
    write(*,'(10f8.1)') fi(1:nfi)
    write(*,'(a,i3,a)') ' reduced np by factor of ',kback,' for backward MC'

end subroutine DisplaySettings
