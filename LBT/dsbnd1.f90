!mm220221
subroutine dsbnd1(icycl,dvw,pvw,nb,bx,dbx,pe2,ce,pbnd,mlim)

! transmittance in incident direction in box system.
!--- history
!mm220109 created
!--- parameter
! kx,ky,kz            declared numbers of boxes in x/y/z directions
! kb1                 kb1=max(kx, ky, kz) + 1
! knang               declared number of scattering angles

  use paras
  implicit none
  character(len=64) :: erc
  real(8), save ::  eps
  integer, save :: mlim1,mlim2
  real(8) :: dvw(3),pvw(3),bx(kb1,3),dbx(kb1,3,2),pbnd
  integer :: nb(3),mlim(*)
  real(8) :: ps(3),ds(3),pe2(3), &
       cpc(3), ce(kx,ky,kz), taue
  integer :: ibs(3),ibe2(3)
  integer :: ind1,icycl

! initialization
  !if(init.ne.0) then
  !   init=0
     call cpcon(cpc)
     eps=cpc(1)*100.0
     mlim1=mlim(1)
     mlim2=mlim(2)
  !endif
  
! backward tracing from the receiver
  ps(1:3)=pvw(1:3)
  ds(1:3)=dvw(1:3)
! find box, ibs: found position
  call fbox(ps,bx,nb,ibs,erc)
  if (erc.ne.' ') then
    pbnd=0.0
    ind1=9
    return
  endif
  
  call strnm(icycl,ds,ibs,ps,ibe2,pe2,nb,bx,dbx, &
            ce,taue,pbnd,eps,mlim2,ind1)
  
  return
end subroutine dsbnd1