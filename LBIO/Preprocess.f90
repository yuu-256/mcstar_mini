subroutine Preprocess()

    !!! DECLARATION:

    !!! EXECUTION:

    ! emerging zenith angles  i=1:na1 th1<0 (upgoing), na1+1-2*na1 th1>0 (downgoing)
    do i=1,na1
      th1(i+na1)=th1(i)
      th1(i)=180-th1(i)
    enddo

    na1=2*na1
    do i=1,nfi
      fi(i+nfi)=fi(i)+180
    enddo
    nfi=2*nfi
    nfi12=nfi/2

    nln1=nln+1

    do lpv=1,npv
      nrnd=2
      call rdmu2(dseed,nrnd,rndv)
      rpvw(lpv,1:2)=rpvw(lpv,1:2)*(1-100*eps*rndv(1:2))
    enddo


    
end subroutine Preprocess
