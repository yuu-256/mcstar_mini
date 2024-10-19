subroutine Initialize(na12, na0, nfi, npv, ipol, nx, ny, sol, ais2, aie2, flxsp2, flxep2, flxs2, flxs02, flxuxy2, flxd002)
  
    !!! DECLARATION:
    implicit none
    ! input
    integer, intent(in) :: na12, na0, nfi, npv, ipol, nx, ny
    ! output
    real(8), intent(out) :: sol
    real(8), dimension(na12,na0,nfi,npv,2,ipol,3) :: ais2
    real(8), dimension(na12,nfi,npv,2,ipol,3) :: aie2
    real(8), dimension(na0,npv,3,3) :: flxsp2
    real(8), dimension(npv,3,3) :: flxep2
    real(8), dimension(na0,6,3) :: flxs2
    real(8), dimension(na0,6) :: flxs02
    real(8), dimension(na0,nx,ny) :: flxuxy2
    real(8), dimension(na0) :: flxd002

    !!!EXECUTION: 
    sol=0
    ais2(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)=0
    aie2(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)=0
    flxsp2(1:na0,1:npv,1:3,1:3)=0
    flxep2(1:npv,1:3,1:3)=0
    flxs2(1:na0,1:6,1:3)=0
    flxs02(1:na0,1:6)=0
    flxuxy2(1:na0,1:nx,1:ny)=0
    flxd002(1:na0)=0

end subroutine Initialize
