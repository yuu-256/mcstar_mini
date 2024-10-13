subroutine cpcon(c)
! machine constants of computer
!--- output
! c     r(3)      (1)  minimum positive x for  1+x      .ne. 1
!                      average as a result of complex arithmetic
!                      operations.
!                 (2)  minimum exponent y for  10.0**y  .ne. 0
!                 (3)  maximum exponent z for  10.0**z  is max. value
!                  if init=1 (data statement) then  set as z=y
!                  if init=2 then this routine gets actual value of z.
!                  - see note -
!--- history
! 90. 1.20  created
!     6.27  change the algorithm to get x and y taking into account the
!           high accuracy co-processor and graceful underflow.
! 92. 3.21  n=1000 from n=200
!     7. 9  bug in x-definition
! 03. 9. 9   changed fpr free source term
!--- note
! this program will generate -underflow error- and -overflow error-
!  messages.  on some computer -overflow error- message may be
!  fatal error.  in that case, please set init = 1 in the data
!  satement for suppressing the procedure of getting c(3).
!--
  use paras, only: pi
  implicit none
  real(8),intent(out)::c(3)
  real(8),save::x,y,z
  integer,save::init=1
  real(8)::y1,y2,y3,z1,z2
  integer::n,m0,m,k,i
!  character ch*80
!--exec
! resolution of computation
  if(init <= 0) then
     c(1)=x; c(2)=y; c(3)=z
     return
  endif

!! test sum(k=1,m) cos((2k-1)*pi/(2m+1)) = 0.5
!! simple check, x = x + e, is not valid when the computer
!!  use a high accurate co-processor.
  n=500; m0=10; x=0.d0
  do m=1,m0
     y=0.d0
     do k=1,m
        y=y+cos((2.d0*k-1)*pi/(2.d0*m+1))
     enddo
     y=abs(2.d0*y-1)
     x=x+y
  enddo
  x=x/m0
  c(1)=x

! exponent for minimum positive value
! this procedure will generate -underflow error message-
  y2=1.d0; n=1000
  do i=1,n
     y1=y2; y3=y1/10.d0
!! for graceful underflow
!! even y2 becomes 0 as output, y2 is not 0 inside
!! computer when graceful underflow is applied.
!! so we replace the value of y2 by output.
!         ch='0'
!         write(ch,107) y3
!  107    format(1p,e12.5)
!         y2=0
!         read(ch,*) y2
     y2=y3
     if(abs(10.d0*y2/y1-1.d0) > 5.d0*x) goto 104
  enddo
  i=n+1
104 y=1-i
  c(2)=y

! exponent for maximum positive value
! this procedure will generate -overflow message-
  if(init <= 1) then
     z=-y
  else
     z2=1
     do i=1,n
        z1=z2; z2=z1*10.d0
        if(abs(z2/z1/10.d0-1.d0) > 5.d0*x) goto 106
     enddo
     i=n+1
106  z=i-1
  endif
  c(3)=z

  init=0
  return
end subroutine cpcon
