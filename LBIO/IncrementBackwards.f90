subroutine IncrementBackwards(na12,na0,nfi,npv,ipol,inda,fsol,dw,sol,ais1,aie1,flxsp1,flxep1,ais2,aie2,flxsp2,flxep2)

    !!! DECLARATION:
    implicit none
    integer,intent(in) :: na12,na0,nfi,npv,ipol,inda
    real(8),intent(in) :: fsol,dw
    real(8),intent(inout) :: sol
    real(8),intent(in) :: ais1(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3),aie1(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3),flxsp1(1:na0,1:npv,1:3,1:3),flxep1(1:npv,1:3,1:3)
    real(8),intent(inout) :: ais2(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3),aie2(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3),flxsp2(1:na0,1:npv,1:3,1:3),flxep2(1:npv,1:3,1:3)

    !!! EXECUTION:
    sol=sol+fsol*dw
    if(inda==1) then
    ais2(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)=ais2(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)  &
        +ais1(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)*dw
    aie2(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)=aie2(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)  &
        +aie1(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)*dw
    end if
    flxsp2(1:na0,1:npv,1:3,1:3)=flxsp2(1:na0,1:npv,1:3,1:3)  &
        +flxsp1(1:na0,1:npv,1:3,1:3)*dw
    flxep2(1:npv,1:3,1:3)=flxep2(1:npv,1:3,1:3)  &
        +flxep1(1:npv,1:3,1:3)*dw

end subroutine IncrementBackwards
