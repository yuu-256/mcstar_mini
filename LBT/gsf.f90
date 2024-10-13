subroutine gsf(m,n,lmax1,amu,q)
! generalized spherical functions
!--- history
! 95. 7. 1     created
!     8. 5     add  if(...) return
!mm210825 kmmax1 => lmax1  
!--- input
! m      i
! n      i
! lmax1  i     maximux order of l + 1
! amu    r     -1 <,= amu <,= 1
!--- output
! q    r(klmax1)   gsf function
!--- memo ---
! p(l,m,n)=q*(i)**(n-m)
!------------
!  use paras, only: kmmax1 !mm210825
  implicit none
  integer,intent(in):: m,n,lmax1
  real(8),intent(in):: amu
  real(8),intent(out):: q(lmax1) !mm210825
  
  integer:: nn,mm,lmin,lmin1,lmax,k,l,l1
  real(8):: fct,fct1
  real(8):: c1,c2,c3
!--
  q(1:lmax1)=0.d0  !mm210825

  nn=abs(n)
  mm=abs(m)
  lmin=max(mm,nn)
  lmin1=lmin+1
  lmax=lmax1-1

  if(lmin > mm)then
! fct=(2*nn)!/(nn-m)!/(nn+m)!
!    = exp{ sum;k=0,nn-m-1 (log(2*nn-k)-log(nn-m-k)) }
     fct1=0
     do k=1,nn-m
        fct1=fct1+log(dble(2*nn-k+1))-log(dble(nn-m-k+1))
     enddo
     fct=exp(fct1)

     if(n > 0)then
        q(lmin1)=(-1)**(nn-m)*sqrt(fct) &
             *(1.d0-amu)**((nn-m)/2.)*(1.d0+amu)**((nn+m)/2.)/2.**nn
     else
        q(lmin1)=sqrt(fct) &
             *(1.d0-amu)**((nn+m)/2.)*(1.d0+amu)**((nn-m)/2.)/2.**nn
     endif
  else
     fct1=0.d0
     do k=1,mm-n
        fct1=fct1+log(dble(2*mm-k+1))-log(dble(mm-n-k+1))
     enddo
     fct=exp(fct1)
     if(m > 0)then
        q(lmin1)=sqrt(fct) &
             *(1-amu)**((mm-n)/2.)*(1+amu)**((mm+n)/2.)/2.**mm
     else
        q(lmin1)=(-1)**(mm-n)*sqrt(fct) &
             *(1-amu)**((mm+n)/2.)*(1+amu)**((mm-n)/2.)/2.**mm
     endif
  endif

  if(lmin1==lmax1) return

  l=lmin
  c1=sqrt(dble((l+m+1)*(l-m+1)*(l+n+1)*(l-n+1)))/dble((l+1)*(2*l+1))
  if(l/=0)then
     c3=m*n/dble(l*(l+1))
  else
     c3=0
  endif
  q(lmin1+1)=(amu-c3)*q(lmin1)/c1

  if(lmin1+1==lmax1) return

  do l=lmin1,lmax-1
     l1=l+1
     c1=sqrt(dble((l+m+1)*(l-m+1)*(l+n+1)*(l-n+1)))/dble((l+1)*(2*l+1))
     c2=sqrt(dble((l+m)*(l-m)*(l+n)*(l-n)))/dble(l*(2*l+1))
     c3=m*n/dble(l*(l+1))
     q(l1+1)=((amu-c3)*q(l1)-c2*q(l1-1))/c1
  enddo

  return
end subroutine gsf

