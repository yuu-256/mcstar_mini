subroutine  qgausn( gwt, gmu, m )
!  compute weights and abscissae for ordinary gaussian quadrature
!   (no weight function inside integral) on the interval (0,1)
!--- history
! 90. 1.17  registered
! 08. 9.25 modified Fortran 90 free-style form
!--- input
! m        i       order of quadrature rule
!--- output
! gmu    r(m)      array of abscissae (0, 1)
! gwt    r(m)      array of weights   sum=1
!--- notes
! reference:  Davis,p.j. and p. Rabinowitz, methods of numerical
!             integration,academic press, new york, 1975, pp. 87
! method:     compute the abscissae as roots of the legendre
!             polynomial p-sub-n using a cubically convergent
!             refinement of newton's method.  compute the
!             weights from eq. 2.7.3.8 of davis/rabinowitz.
!             accuracy:  at least 13 significant digits.
!--- internal variables
! pm2,pm1,p : 3 successive legendre polynomials
! ppr       : derivative of legendre polynomial
! p2pri     : 2nd derivative of legendre polynomial
! tol       : convergence criterion
! x,xi      : successive iterates in cubically-
!             convergent version of newton's method
!            ( seeking roots of legendre polynomial )
!--
  use paras, only:pi
  implicit none

! input & output
  integer,intent(inout)::m
  real(8),intent(out)::gmu(m),gwt(m)

! work
  integer::lim,np1,k,nn
  real(8)::cona, t, en, nnp1, p, pm1, pm2, ppr, p2pri, prod, tmp, x, xi
  real(8)::tol=1.0d-13
!--exec
  if ( m <= 1 )  then
     m = 1
     gmu( 1 ) = 0.5d0
     gwt( 1 ) = 1.0d0
     return
  end if

  en   = m
  np1  = m + 1
  nnp1 = m * np1
  cona = dble( m-1 ) / ( 8.d0 * m**3 )
!+---------------------------------------------------------------------+
!!         initial guess for k-th root of legendre polynomial,         |
!!         from davis/rabinowitz  eq. (2.7.3.3a)                       |
!+---------------------------------------------------------------------+
  lim  = m / 2
  do k = 1, lim
     t = dble( 4*k - 1 ) * pi / dble( 4*m + 2 )
     x = cos ( t + cona / tan( t ) )

!+---------------------------------------------------------------------+
!!             recursion relation for legendre polynomials             |
!!       initialize legendre polynomials: (pm2) p-sub-0, (pm1) p-sub-1 |
!+---------------------------------------------------------------------+
10   pm2 = 1.d0
     pm1 = x
     do nn = 2, m
!            p   = ( ( 2*nn - 1 ) * x * pm1 - ( nn-1 ) * pm2 ) / nn
        p   =  2* x * pm1 - pm2 - ( x * pm1 - pm2 ) / nn
        pm2 = pm1
        pm1 = p
     enddo

     tmp   = 1.d0 / ( 1.d0 - x**2 )
     ppr   = en * ( pm2 - x * p ) * tmp
     p2pri = ( 2.d0 * x * ppr - nnp1 * p ) * tmp
     xi    = x - ( p / ppr ) * ( 1.d0 + ( p / ppr ) &
          * p2pri / ( 2.d0 * ppr ) )

     if ( dabs(xi-x) > tol ) then
!!          check for convergence
        x = xi
        go to 10
     end if

!       ** iteration finished--calc. weights, abscissae for (-1,1)
     gmu( k ) = - x
     gwt( k ) = 2.d0 / ( tmp * ( en * pm2 )**2 )
     gmu( np1 - k ) = - gmu( k )
     gwt( np1 - k ) =   gwt( k )
  enddo

  if ( mod( m,2 ) /= 0 )  then
!!       set middle abscissa and weight for rules of odd order
     gmu( lim + 1 ) = 0.d0
     prod = 1.d0
     do k = 3, m, 2
        prod = prod * k / ( k-1 )
     enddo
     gwt( lim + 1 ) = 2.d0 / prod**2
  end if

  do k = 1, m
!!       convert from (-1,1) to (0,1)
     gmu( k ) = 0.5d0 * gmu( k ) + 0.5d0
     gwt( k ) = 0.5d0 * gwt( k )
  enddo
  return
end subroutine qgausn
