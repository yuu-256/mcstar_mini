subroutine ReadData(iui,fi,wlmin,wlmax,galbs,galbl,matm,nln,ipbf,ifrh,trh,npoly,icn,wlcn,ncomp,cnpt,mptc,vptc,icycl,np0,nx,ny,dx,dy,fis,indg,jatm,cconc,npv,llr,llu,rpvw)
!!! 宣言部未完
  !!! DECLARATIONS:
  implicit none
  integer, intent(in) :: iui
  real(8), intent(out) :: fi(:)
  real(8), intent(out) :: wlmin,wlmax
  real(8), intent(out) :: galbs(:,:),galbl(:,:)
  integer, intent(out) :: matm,nln
  integer, intent(out) :: ipbf(:)
  integer, intent(out) :: ifrh
  real(8), intent(out) :: trh
  integer, intent(out) :: npoly,icn
  real(8), intent(out) :: wlcn
  integer, intent(out) :: ncomp(:),cnpt(:)
  integer, intent(out) :: mptc(:,:),vptc(:,:)
  integer, intent(out) :: icycl,np0
  integer, intent(out) :: nx,ny
  real(8), intent(out) :: dx,dy
  real(8), intent(out) :: fis
  integer, intent(out) :: indg(:,:)
  integer, intent(out) :: jatm(:,:)
  real(8), intent(out) :: cconc(:,:,:,npoly)
  integer, intent(out) :: npv,llr
  integer, intent(out) :: llu(:),rpvw(:,:)

  integer :: i,j,k,lz,k1,l1

  !!! EXECUTION:
  read(iui,*) nfi,(fi(i),i=1,nfi)
  read(iui,*) wlmin,wlmax
  read(iui,*) galbs(1,1),galbl(1,1)
  read(iui,*) matm,nln
  read(iui,*) (ipbf(i),i=1,nln1)
  read(iui,*) ifrh,trh
  read(iui,*) npoly,icn,wlcn

  do  i=1,npoly
    read(iui,*) ncomp(i),cnpt(i)
    read(iui,*) (mptc(i,j),j=1,ncomp(i))
    read(iui,*) (vptc(i,j),j=1,ncomp(i))
  enddo

  read (iui, *, end=1, err=1) icycl,np0
  read(iui,*) nx,ny,dx,dy  !tr for mc7d
  read(iui,*) fis

  if( indg(1,1) < 0 ) then
    do i = 1, nx
      read(iui,*) (indg(i,j),j=1,ny)
    end do
  else
    indg(1:nx,1:ny)=indg(1,1)
  end if

  if( galbs(1,1) < 0 ) then
    do i = 1, nx
      read(iui,*) (galbs(i,j),j=1,ny)
    end do
  else
    galbs(1:nx,1:ny)=galbs(1,1)
  end if

  if( galbl(1,1) < 0 ) then
    do i = 1, nx
      read(iui,*) (galbl(i,j),j=1,ny)
    end do
  else
      galbl(1:nx,1:ny)=galbl(1,1)
  end if

  if( matm <= 0 ) then
    do i = 1, nx
      read(iui,*) (jatm(i,j),j=1,ny)
    end do
  else
    jatm(1:nx,1:ny)=matm
  end if

  do k=1,npoly
    do lz=1,nln  ! bottom to top
      read(iui,*) k1,l1
      do i=1,nx
        read(iui,*) (cconc(i,j,lz,k),j=1,ny)
      enddo
    enddo
  enddo

! viewing (receiver) positions on a box boundary plane
! llu: boundary plane number of viewing (receiver) positions (1-6)
!  1 xy+z (toa), 2 xy-z (boa),  3 yz+x,  4 yz-x,  5 zx+y,  6 zx-y
! llr: 0 for position prescribed by rpvw, 1 by random position
! rpvw: relative location of the viewing position on the boundary plane

  read(iui,*) npv,llr
  do lpv=1,npv
    read(iui,*) llu(lpv),rpvw(lpv,1:2)
    nrnd=2
    call rdmu2(dseed,nrnd,rndv)  !　これの宣言を書いておく
    rpvw(lpv,1:2)=rpvw(lpv,1:2)*(1-100*eps*rndv(1:2))
  enddo

end subroutine ReadData
