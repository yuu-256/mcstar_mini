subroutine WriteOutput(OutFile, OutData, OutDataSize)
    !!! DECLARATION:
    implicit none

    !!! EXECUTION: ON DEVELOPMENT!
    !-- Write simulation parameters
    write(iuo,*) icont(1:3),' icont (itray itkd immax)'
    write(iuo,'(7i5,a)') isol,inda,indg(1,1),imthd,ipol,nda,nds,  &
                       ' isol inda indg imthd ipol nda nds'
    write(iuo,'(i5,a)') na0,' na0/th0'
    write(iuo,'(8f10.3)') th0(1:na0)
    write(iuo,'(i5,a)') na1,' na1/th1'
    write(iuo,'(8f10.3)') th1(1:na1)
    write(iuo,'(i5,a)') nfi,' nfi/fi'
    write(iuo,'(8f10.3)') fi(1:nfi)
    write(iuo,'(2f10.3,a)') wlmin,wlmax,' wl-min -max'
    write(iuo,'(2f10.3,a)') galbs(1,1),galbl(1,1),' ground albedos for SW and LW bands'
    write(iuo,*) matm,nln,' matm nln'
    write(iuo,*) 'ipbf'
    write(iuo,'(15i5)') ipbf(1:nln1)
    write(iuo,'(i5,f10.3," ifrh trh")') ifrh,trh
    write(iuo,'(2i5,1pe11.3," npoly icn wlcn")') npoly,icn,wlcn
    do  i=1,npoly
      write(iuo,'(i5,f10.3," ncomp cnpt")') ncomp(i),cnpt(i)
      write(iuo,*) 'mptc vptc'
      write(iuo,'(10(i5,f10.3))') (mptc(i,j),vptc(i,j),j=1,ncomp(i))
    enddo
    write(iuo,*)

    !-- Write MC parameters
    write(iuo,*) icycl,np0,' icycl np0'
    write(iuo,'(2i4,2f10.3,a)') nx,ny,dx,dy,' nx ny dx dy'
    write(iuo,'(f8.2,a)') fis,' fis solar azimuthal angle relative to x-axis'

    if( indg(1,1) < 0 ) then
      write(iuo,'(a)') ' surface type inhomogeneous'
      do i = 1, nx
        write(iuo,'(10i3)') (indg(i,j),j=1,ny)
      end do
    else

    if( galbs(1,1) < 0 ) then
      write(iuo,'(a)') ' surface albedo inhomogeneous (SW)'
      do i = 1, nx
        write(iuo,'(10f10.3)') (galbs(i,j),j=1,ny)
      end do
    end if

    if( galbl(1,1) < 0 ) then
      write(iuo,'(a)') ' surface albedo inhomogeneous (LW)'
      do i = 1, nx
        write(iuo,'(10f10.3)') (galbl(i,j),j=1,ny)
      end do
    end if

    if( matm <= 0 ) then
      write(iuo,'(a)') ' atmospheric inhomogeneous'
      do i = 1, nx
        write(iuo,'(10i3)') (jatm(i,j),j=1,ny)
      end do
    end if
    write(iuo,*) 'cconc(ix,iy,iz,ipoly)'
    do k=1,npoly
      do lz=1,nln  ! bottom to top
        write(iuo,'(2i4,a)') k,lz,' ipoly lz'
        do i=1,nx
          write(iuo,'(10f8.4)') (cconc(i,j,lz,k),j=1,ny)
        enddo
      enddo
    enddo

    !-- Write viewing positions
    write(iuo,'(/2i5,a)') npv,llr,' number of viewing positions and prescribed/random'
    write(iuo,*) 'llu : boundary plane 1 xy+z (toa), 2 xy-z (boa), 3 yz+x, 4 yz-x, 5 zx+y, 6 zx-y'
    write(iuo,*) llr,' llr: 0 for user defined viewing position on the boundary plane, 1 randome allocation'
    write(iuo,*) 'rpvw values are slightly randomized to avoid falling between adjacent boxes.'
    write(iuo,'(a)') '  lpv  llu   rpvw1   rpvw2'
    do lpv=1,npv
      write(iuo,'(2i5,2f8.3)') lpv,llu(lpv),rpvw(lpv,1:2)
    enddo

    !-- Write output data
    
end subroutine WriteOutput
