subroutine gtgas(iug,nbnd,wbnd)

    !!! DECLARATION:
    use params
    implicit none
    integer, intent(in) :: iug
    integer, intent(out) :: nbnd
    real(8), intent(out) :: wbnd(knw)
    integer :: iw, ibnd
    integer :: ngas(2)
    real(8) :: wlmin, wlmax


    
    !!! EXECUTION:
    if(icont(2)<=1) then ! broad band
        read(iug,*) nbnd
        read(iug,*)
        read(iug,*) wbnd(1:nbnd+1)
    else if(icont(2)==2) then ! narrow band
        nbnd=int(log10(1.d4/wlmin)*ngas(2)+1.d0)-int(log10(1.d4/wlmax)*ngas(2)+1.d0)+1
        ibnd=int(log10(1.d4/wlmax)*dble(ngas(2)))
        do iw = 0, nbnd
          if(iw==0) then
            wbnd(iw+1)=1.d4/wlmax
          else if(iw==nbnd) then
            wbnd(iw+1)=1.d4/wlmin
          else
           wbnd(iw+1)=10**((ibnd+iw)/dble(ngas(2)))
          end if
        end do
    end if

end subroutine gtgas
