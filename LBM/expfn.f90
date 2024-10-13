function expfn(x)
! exponential function with over/under flow setting.
!--- history
! 89. 5. 4   created by T. Nakajima.
! 90. 1.17   updated with cpcon
! 03. 9. 9   changed fpr free source term
!--- input
! x        R         independent variable.
!--- output
! expfn    F         exp(x).
!                     if x.le. vmn then exp(x) is reset as   0.
!                     if x.ge. vmx then exp(x) is reset as exp(vmx).
!--- parameters
! system set the -vmn- and -vmx- by the function r1mach.
!--
  use paras, only: ccp
  implicit none
  integer,save::init
  real(8),save::vmn,vmx,expmn,expmx
  real(8),intent(in)::x
  real(8):: expfn
!--exec
  init=1

! set vmn      r         enderflow limit.
!     vmx      r         overflow limt.
  if(init > 0) then
     init=0
     vmn=ccp(2)*0.8d0*2.3d0; vmx=ccp(3)*0.8d0*2.3d0
     expmn=0.0d0; expmx=exp(vmx)
  endif

  if(x <= vmn) then
     expfn=expmn
  else
     if(x >= vmx) then
        expfn=expmx
     else
        expfn=exp(x)
     endif
  endif
  return
end function expfn
