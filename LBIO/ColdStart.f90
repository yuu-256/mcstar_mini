subroutine ColdStart()
  
    !!! DECLARATION:
    implicit none
    integer :: i
    real(8) :: dseed
    integer :: nrnd
    real(8) :: rndv(10)

    !!! EXECUTION:
    print *, 'Cold start of random generator'

    !-- cold start of random generator
    dseed=1.0d0
    nrnd=1
    do i=1,10
        call rdmu2(dseed,nrnd,rndv)
    enddo

    print *, 'Cold start of random generator done'

end subroutine ColdStart
