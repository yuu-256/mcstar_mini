!mm220221
subroutine ggdir(dseed,ds)
!c random generator of direction after hitting ground.
!c--- history
!c 88. 9.21  created
!c--- input
!c dseed       d      initial condition for radom generator
!c--- output
!c ds        r(3)     emergent direction vector.
!  parameter (pi=3.141592654)
  use paras
  implicit none
  real(8) :: dseed
  real(8) :: ds(3)
  real(8) :: rndv(2)
  integer :: nrnd
  real(8) :: cth,sth,fi,cfi,sfi

  nrnd=2
  call rdmu2(dseed,nrnd,rndv)
  cth=sqrt(rndv(1))
  sth=sqrt(abs(1.0-cth**2))
  fi=2.0*pi*rndv(2)
  cfi=cos(fi)
  sfi=sin(fi)
  ds(1)=sth*cfi
  ds(2)=sth*sfi
  ds(3)=cth
  return
end subroutine ggdir

subroutine ggdir1(dseed,ds)
!c random generator of direction after hitting ground.
!c--- history
!c 88. 9.21  created
!c--- input
!c dseed       d      initial condition for radom generator
!c--- output
!c ds        r(3)     emergent direction vector.
!  parameter (pi=3.141592654)
  use paras
  implicit none
  real(8) :: dseed
  real(8) :: ds(3)
  real(8) :: rndv(2)
  integer :: nrnd
  real(8) :: cth,sth,fi,cfi,sfi

  nrnd=2
  call rdmu2(dseed,nrnd,rndv)
  cth=rndv(1)
  sth=sqrt(abs(1.0-cth**2))
  fi=2.0*pi*rndv(2)
  cfi=cos(fi)
  sfi=sin(fi)
  ds(1)=sth*cfi
  ds(2)=sth*sfi
  ds(3)=cth
  return
end subroutine ggdir1

subroutine gttau(dseed,trn,tau)
!c random generator of trial transmissivity
!c--- history
!c 88. 9.21 created
!c--- input
!c dseed       d      initial condition for random generator
!c--- output
!c trn         r      transmissivity
!c tau         r      -log(trns)
!c$endif
  implicit none
  real(8) :: dseed
  real(8) :: rndv(1)
  integer :: nrnd
  real(8) :: trn,tau
  nrnd=1
  call rdmu2(dseed,nrnd,rndv)
  trn=rndv(1)
  if(trn.le.0.0) then
     trn=0
     tau=1.0e10
  else
     tau=-log(trn)
  endif
  return
end subroutine gttau