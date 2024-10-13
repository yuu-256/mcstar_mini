function vlspc2(pr)
! volume spectrum of partile polydisperison: v(x) = dv / d ln r
!--- history
! 95. 9.14 created with parameter packet only
!--- input
! pr    r(6,4)      parameter packet
!
!    pr(1,1)=r       particle radius in cm
!    pr(1,2)=nmode   number of mode radius
!    pr(1,3)=rmin    minimum particle radius in cm
!    pr(1,4)=rmax    maximum particle radius in cm
!
!    for each j-th mode (<= 4)
!
!  pr(2,j): type of function (itp) for the mode.
!   itp=1: power law
!     pr(3,j)=c, pr(4,j)=r0,  pr(5,j)=p
!     vj = c * (r/r0)**(4-p) if r>r0; = c * (r/r0)**4 if r<r0
!   itp=2: log-normal
!     pr(3,j)=c, pr(4,j)=s,   pr(5,j)=rm
!     vj = c * exp((ln(r/rm)/ln(s))**2 / 2)
!   itp=3: modified gamma
!     pr(3,j)=c, pr(4,j)=alfa, pr(5,j)=beta, pr(6,j)=gamma
!     vj = c * (r1)**(alfa+4) exp (-beta*r1**gamma) where r1=r*1.0e4
!--- output
! vlspc2   rf       dv/d ln r
!--
  implicit none
  real(8),intent(in)::pr(6,4)    !! parameter packet
  integer::nmode,itp,m
  real(8)::r,rmin,rmax,e1,r1
  real(8)::c,rc,pdndr,pn,s,rm,alf,bet,gam
  real(8)::vlspc2
!--exec
  r=pr(1,1); nmode=int(pr(1,2))
  rmin =pr(1,3); rmax =pr(1,4)
  vlspc2=0.0d0
  if(r < rmin .or. r > rmax) return

  do m=1,nmode
     itp=int(pr(2,m))

     if(itp==1) then
! power law
        c=pr(3,m); rc=pr(4,m); pdndr=pr(5,m)
        if(r <= rc) then
           pn=4.d0
        else
           pn=4.d0-pdndr
        endif
        e1=pn*log(r/rc)

     elseif(itp==2) then
! log-normal
        c =pr(3,m); s =pr(4,m); rm=pr(5,m)
        e1=-0.5d0*(log(r/rm)/log(s))**2

     else
! modified gamma
        r1=r*1.0d4; c=pr(3,m); alf=pr(4,m); bet=pr(5,m); gam=pr(6,m)
        e1=(alf+4.d0)*log(r1)-bet*r1**gam
     endif
     if(e1 > -100.d0) vlspc2=vlspc2+c*exp(e1)
  enddo
  return
end function vlspc2
