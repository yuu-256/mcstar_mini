!mm220221
subroutine angintp(rd,na,angdp,ind,scat)
  !mm220118 created
  implicit none
  real(8) :: rd
  integer :: na, ind
  real(8) :: angdp(na), scat
  real(8) :: dsca
  integer :: i
  
  dsca = rd*(na-1)
  i = dsca
  if(i<0.or.i>na-1) then
    ind=9
    return
  end if
  scat = angdp(i+1)+(angdp(i+2)-angdp(i+1))*(dsca-i)
  !dsca = 180./(na-1)
  !scat = int(scat/dsca+.5)*dsca
  
  return
end subroutine angintp
