subroutine IncrementForwards(na0, nx, ny, dw, flxs01, flxs1, flxe1, flxuxy1, flxd001, flxs02, flxs2, flxe2, flxuxy2, flxd002)

    !!! DECLARATION:
    implicit none
    integer, intent(in) :: na0, nx, ny
    real(8), intent(in) :: dw
    real(8), intent(in) :: flxs01(na0,6), flxs1(na0,6,3), flxe1(6,3), flxuxy1(na0,nx,ny), flxd001(na0)
    real(8), intent(inout) :: flxs02(na0,6), flxs2(na0,6,3), flxe2(6,3), flxuxy2(na0,nx,ny), flxd002(na0)

    !!! EXECUTION:
    flxs02(1:na0,1:6)=flxs02(1:na0,1:6)+flxs01(1:na0,1:6)*dw
    flxs2(1:na0,1:6,1:3)=flxs2(1:na0,1:6,1:3)+flxs1(1:na0,1:6,1:3)*dw
    flxe2(1:6,1:3)=flxe2(1:6,1:3)+flxe1(1:6,1:3)*dw
    flxuxy2(1:na0,1:nx,1:ny)=flxuxy2(1:na0,1:nx,1:ny)  &
                            +flxuxy1(1:na0,1:nx,1:ny)*dw
    flxd002(1:na0)=flxd002(1:na0)+flxd001(1:na0)*dw

end subroutine IncrementForwards
