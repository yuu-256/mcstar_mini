subroutine ReadData(iua,ius,iug,icont,wlmin,wlmax,natm,nl,alt,prs,tmp,nmol,cng,nptc,ins,cnp,ispcvp,rfracp,asphr,rop,dryap,nawp,awcrp,nv,nwlv,wlv,rfi,nbnd,wbnd)

    !!! DECLARATION:
    implicit none

    ! input
    integer,intent(in):: iua,ius,iug   !! device number
    integer,intent(in):: icont(10)     !! control parameters
    real(8),intent(in):: wlmin,wlmax  !! min and max wavelength [micron]
    
    ! output
    integer,intent(inout):: natm,nl       !! number of atmospheric layer
    real(8),intent(out):: alt(knl) !! altitude [km]
    real(8),intent(out):: prs(knl,katm) !! pressure [mb]
    real(8),intent(out):: tmp(knl,katm) !! temperature [K]
    integer,intent(out):: nmol          !! number of gases
    real(8),intent(inout):: cng(knl,knm0,katm) !! gas concentration [ppmv]

    integer,intent(inout):: nptc          !! number of particle polydispersions
    integer,intent(inout):: ins(kptc)     !! nonspherical flag
    real(8),intent(inout):: cnp(knl,kptc) !! dry volume concentration
    integer,intent(inout):: ispcvp(3,kptc) !! internal mixture of fundamental
    real(8),intent(inout):: rfracp(3,kptc) !! dry component volume fractions
    real(8),intent(inout):: asphr(3,kptc)  !! non-spherical parameters
    real(8),intent(inout):: rop(kptc)     !! particle density of dry mixture [g/cm3]
    real(8),intent(inout):: dryap(6,4,kptc) !! param to define volume size dist.
    integer,intent(inout):: nawp(kptc)    !! number of RH to define particle growth
    real(8),intent(inout):: awcrp(kaw,kptc,2) !! growth parameters

    integer,intent(inout):: nv            !! number of fundamental species
    integer,intent(inout):: nwlv          !! number of wavelenggths for rfi tables
    real(8),intent(inout):: wlv(kwlv)     !! wavelengths [micron] for rfi tables
    real(8),intent(inout):: rfi(kwlv,knv,2)  !! refractive index of fundamental
    
    ! CKD parameters
    integer,intent(out):: nbnd
    real(8),intent(out):: wbnd(knw)

    !!! EXECUTION:
    ! model atmospheres
    call mlatm(iua,nl,natm,nm1,nm2,idm,idms,wmol,rams,airm, &
          alt,pmatm,tmatm,dnsty,amol,trac)

    ! get aerosol parameters
    call gtpar7(ius,nl,nptc,ins,iver,cnp1,ispcvp,rfracp,asphr,rop, &
          dryap,nawp, awcrp,nv,nwlv,wlv,rfi)

    ! get CKD parameters
    call gtgas(iug,icont,wlmin,wlmax,nbnd,wbnd)
    
end subroutine ReadData
