!mm220221 not use
subroutine phdns(s,nang,sca,phs,phs0,erc)
!c get phase function from scattering angle
!c--- history
!c 88. 9.16 created
!c--- input
!c s         r         scattering angle
!c nang      i         number of scattering angles
!c sca     r(nang)     scattering angles
!c phs     r(nang)     phase functions
!c--- output
!c phs0      r         phase function values
!c erc     c*64        ' ' no error
!c$endi
  implicit none
  character(len=64) :: erc
  real(8) ::  sca(nang),phs(nang)
  integer :: nang
  real(8) :: s,phs0
  real(8) :: bintp
  if(s.lt.0.0 .or. s.gt.180.0) then
     erc='phase'
     return
  endif
  phs0=bintp(s,nang,sca,phs)
  erc=' '
  return
end subroutine phdns