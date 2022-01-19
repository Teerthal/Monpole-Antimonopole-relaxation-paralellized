      subroutine fdflux(f,i,j,k,r,hatn,dxhatn,dyhatn,dzhatn,lxb,lyl,lzd)
      implicit none
      include 'parameters.inc'
      integer i,j,k,n,nn,lxb,lyl,lzd,ilocal,jlocal,klocal
      real*8 f,r
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 cd,g
      real*8 x,y,z,rr
      real*8 rmserror
      real*8 dfddx,dfddy,dfddz
      real*8 res
      real*8 laplacian

      real*8 hatn,dxhatn,dyhatn,dzhatn

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension hatn(3,1:localsize,1:localsize,1:localsize)
      dimension dxhatn(3,1:localsize,1:localsize,1:localsize)
      dimension dyhatn(3,1:localsize,1:localsize,1:localsize)
      dimension dzhatn(3,1:localsize,1:localsize,1:localsize)

      dimension r(nf)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
      dimension laplacian(nf)

!=======================================================================
! fs=field strength; first index is the number of gauge fields.
!=======================================================================

      dimension fs(3,0:3,0:3),dfddx(nf),dfddy(nf),dfddz(nf)

!=======================================================================
! cd=covariant derivative; second index is the gauge field index.
!=======================================================================

      dimension cd(0:3,4),g(nf)
      dimension res(nf)

!=======================================================================
! fields for SO(3):
!     f(1)=phi^1,f(2)=phi^2,f(3)=phi^3,
!     f(4)=W^1_0, f(5)=W^1_1, f(6)=W^1_2, f(7)=W^1_3,
!     f(8)=W^2_0, f(9)=W^2_1, f(10)=W^2_2, f(11)=W^2_3,
!     f(12)=W^3_0, f(13)=W^3_1, f(14)=W^3_2, f(15)=W^3_3,
!     gauge fixing functions (Num Rel): f(16)=\partial_i W^1_i, 
!     f(17)=\partial_i W^2_i, f(18)=\partial_i W^3_i, 
!     f(19)=|{\vec \phi}|
!     time derivatives:
!     fd(i)=time derivative of f(i) for i=1-3 (scalar fields),
!     fd(i)=electric fields(=W^a_{0j}) for i=4-15 (gauge fields).
!     fd(19)=time derivative of |phi|
!=======================================================================

!local here seem to imply global coordinates.
            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

! spatial derivatives to 4th order:

      call derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)
!     call derivativesPBC(f,i,j,k,dfdx,dfdy,dfdz,
!    1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

! covariant derivatives:

      call covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)

! find gauge field strengths (uses derivatives):

      call fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)

! for convenience of writing:

      do n=1,nf
      g(n)=f(n,i,j,k)
      enddo

!=======================================================================
! define laplacian() such that the expression is accurate near the
! boundaries as well. For this relaxation code, we will use the boundary 
! conditions: $D_i\phi^a =0$ and $W^a_{ij}=0$. These relations then let
! us write second order derivatives in terms of first order derivatives.
! if on boundaries:

      if(abs(ilocal).eq.latx.or.abs(jlocal)
     1                          .eq.laty.or.abs(klocal).eq.latz) then

          do n=1,nf

! laplacian on the boundary:
! scalar fields;

          laplacian(1)=gw*(
     1                     dfdx(2)*g(13)+dfdy(2)*g(14)+dfdz(2)*g(15)
     2                    -dfdx(3)*g(9)-dfdy(3)*g(10)-dfdz(3)*g(11)
     3                    +g(2)*g(18)-g(3)*g(17))
          laplacian(2)=gw*(
     1                     dfdx(3)*g(5)+dfdy(3)*g(6)+dfdz(3)*g(7)
     2                    -dfdx(1)*g(13)-dfdy(1)*g(14)-dfdz(1)*g(15)
     3                    +g(3)*g(16)-g(1)*g(18))
          laplacian(3)=gw*(
     1                     dfdx(1)*g(9)+dfdy(1)*g(10)+dfdz(1)*g(11)
     2                    -dfdx(2)*g(5)-dfdy(2)*g(6)-dfdz(2)*g(7)
     3                    +g(1)*g(17)-g(2)*g(16))

! gauge fields;

          laplacian(4)=0.
          laplacian(5)=dfdx(16)-gw*(g(17)*g(13)-g(18)*g(9))
     1   +gw*(g(9)*dfdx(13)+g(10)*dfdy(13)+g(11)*dfdz(13)
     2       -g(13)*dfdx(9)-g(14)*dfdy(9)-g(15)*dfdz(9))
          laplacian(6)=dfdy(16)-gw*(g(17)*g(14)-g(18)*g(10))
     1   +gw*(g(9)*dfdx(14)+g(10)*dfdy(14)+g(11)*dfdz(14)
     2       -g(13)*dfdx(10)-g(14)*dfdy(10)-g(15)*dfdz(10))
          laplacian(7)=dfdz(16)-gw*(g(17)*g(15)-g(18)*g(11))
     1   +gw*(g(9)*dfdx(15)+g(10)*dfdy(15)+g(11)*dfdz(15)
     2       -g(13)*dfdx(11)-g(14)*dfdy(11)-g(15)*dfdz(11))

          laplacian(8)=0.
          laplacian(9)=dfdx(17)-gw*(g(18)*g(5)-g(16)*g(13))
     1   +gw*(g(13)*dfdx(5)+g(14)*dfdy(5)+g(15)*dfdz(5)
     2       -g(5)*dfdx(13)-g(6)*dfdy(13)-g(7)*dfdz(13))
          laplacian(10)=dfdy(17)-gw*(g(18)*g(6)-g(16)*g(14))
     1   +gw*(g(13)*dfdx(6)+g(14)*dfdy(6)+g(15)*dfdz(6)
     2       -g(5)*dfdx(14)-g(6)*dfdy(14)-g(7)*dfdz(14))
          laplacian(11)=dfdz(17)-gw*(g(18)*g(7)-g(16)*g(15))
     1   +gw*(g(13)*dfdx(7)+g(14)*dfdy(7)+g(15)*dfdz(7)
     2       -g(5)*dfdx(15)-g(6)*dfdy(15)-g(7)*dfdz(15))

          laplacian(12)=0.
          laplacian(13)=dfdx(18)-gw*(g(16)*g(9)-g(17)*g(5))
     1   +gw*(g(5)*dfdx(9)+g(6)*dfdy(9)+g(7)*dfdz(9)
     2       -g(9)*dfdx(5)-g(10)*dfdy(5)-g(11)*dfdz(5))
          laplacian(14)=dfdy(18)-gw*(g(16)*g(10)-g(17)*g(6))
     1   +gw*(g(5)*dfdx(10)+g(6)*dfdy(10)+g(7)*dfdz(10)
     2       -g(9)*dfdx(6)-g(10)*dfdy(6)-g(11)*dfdz(6))
          laplacian(15)=dfdz(18)-gw*(g(16)*g(11)-g(17)*g(7))
     1   +gw*(g(5)*dfdx(11)+g(6)*dfdy(11)+g(7)*dfdz(11)
     2       -g(9)*dfdx(7)-g(10)*dfdy(7)-g(11)*dfdz(7))

! gauge constraints:

          laplacian(16)=0.
          laplacian(17)=0.
          laplacian(18)=0. 

!=======================================================================
! Higgs field magnitude:
! Here using the kink BPS equation: $df/dx = -\sqrt{2V}$, and
! generalizing it to our case: 
! $\partial_i f = \pm n_i\sqrt{2V} = n_i \sqrt{\lambda/2}(vev^2-f^2)$
! (where the sign of the square root has been chosen to that f=|\phi| 
! grows to its vev at infinity). This then gives the Laplacian as
! $-\sqrt{2*lambda}*f*n_i\partial_i f 
!        + \sqrt{lambda/2}(vev^2-f^2)*\partial_i n_i$.
!=======================================================================

          laplacian(19)=-sqrt(2.*lambda)*g(19)*(hatn(1,i,j,k)*dfdx(19)
     1                 +hatn(2,i,j,k)*dfdy(19)+hatn(3,i,j,k)*dfdz(19))
     2          +sqrt(lambda/2.)*(vev**2-g(19)**2)
     3              *(dxhatn(1,i,j,k)+dyhatn(2,i,j,k)+dzhatn(3,i,j,k))

          enddo

! if not on boundaries:
      else
          do n=1,nf
          laplacian(n)=(1.0-vel**2)*d2fdx2(n)+d2fdy2(n)+d2fdz2(n)
          enddo
      endif

!=======================================================================
! Scalar fluxes:
!=======================================================================

! phi^1:
      r(1)= laplacian(1)
     1 +gw*((g(9)*dfdx(3)+g(10)*dfdy(3)+g(11)*dfdz(3))
     2     -(g(13)*dfdx(2)+g(14)*dfdy(2)+g(15)*dfdz(2)))
     3 +gw*((g(9)*cd(1,3)+g(10)*cd(2,3)+g(11)*cd(3,3))
     4     -(g(13)*cd(1,2)+g(14)*cd(2,2)+g(15)*cd(3,2)))
     5 -lambda*(g(1)**2+g(2)**2+g(3)**2-vev**2)*g(1)
     6 +gw*(g(17)*g(3)-g(18)*g(2))


! phi^2:
      r(2)= laplacian(2)
     1 +gw*((g(13)*dfdx(1)+g(14)*dfdy(1)+g(15)*dfdz(1))
     2     -(g(5)*dfdx(3)+g(6)*dfdy(3)+g(7)*dfdz(3)))
     3 +gw*((g(13)*cd(1,1)+g(14)*cd(2,1)+g(15)*cd(3,1))
     4     -(g(5)*cd(1,3)+g(6)*cd(2,3)+g(7)*cd(3,3)))
     5 -lambda*(g(1)**2+g(2)**2+g(3)**2-vev**2)*g(2)
     6 +gw*(g(18)*g(1)-g(16)*g(3))
! phi^3:
      r(3)= laplacian(3)
     1 +gw*((g(5)*dfdx(2)+g(6)*dfdy(2)+g(7)*dfdz(2))
     2     -(g(9)*dfdx(1)+g(10)*dfdy(1)+g(11)*dfdz(1)))
     3 +gw*((g(5)*cd(1,2)+g(6)*cd(2,2)+g(7)*cd(3,2))
     4     -(g(9)*cd(1,1)+g(10)*cd(2,1)+g(11)*cd(3,1)))
     5 -lambda*(g(1)**2+g(2)**2+g(3)**2-vev**2)*g(3)
     6 +gw*(g(16)*g(2)-g(17)*g(1))

!
! W-fluxes (assumes $W^a_0=0$ gauge):
!
! W^1_{00}:
      r(4)=0.
! W^1_{01}:
      r(5)= laplacian(5)
     1 +gw*(
     2   -(dfdx(9)*f(13,i,j,k)-dfdx(13)*f(9,i,j,k))
     2   -(dfdy(9)*f(14,i,j,k)-dfdy(13)*f(10,i,j,k))
     2   -(dfdz(9)*f(15,i,j,k)-dfdz(13)*f(11,i,j,k))
     3   -(f(9,i,j,k)*fs(3,1,1)-f(13,i,j,k)*fs(2,1,1))
     3   -(f(10,i,j,k)*fs(3,1,2)-f(14,i,j,k)*fs(2,1,2))
     3   -(f(11,i,j,k)*fs(3,1,3)-f(15,i,j,k)*fs(2,1,3)))
     4   -gw*(g(2)*cd(1,3)-g(3)*cd(1,2))
     5   -dfdx(16)+gw*(g(17)*g(13)-g(18)*g(9))

! W^1_{02}:
      r(6)= laplacian(6)
     1 +gw*( 
     2   -(dfdx(10)*f(13,i,j,k)-dfdx(14)*f(9,i,j,k))
     2   -(dfdy(10)*f(14,i,j,k)-dfdy(14)*f(10,i,j,k))
     2   -(dfdz(10)*f(15,i,j,k)-dfdz(14)*f(11,i,j,k))
     3   -(f(9,i,j,k)*fs(3,2,1)-f(13,i,j,k)*fs(2,2,1))
     3   -(f(10,i,j,k)*fs(3,2,2)-f(14,i,j,k)*fs(2,2,2))
     3   -(f(11,i,j,k)*fs(3,2,3)-f(15,i,j,k)*fs(2,2,3)))
     4   -gw*(g(2)*cd(2,3)-g(3)*cd(2,2))
     5   -dfdy(16)+gw*(g(17)*g(14)-g(18)*g(10))
! W^1_{03}:
      r(7)= laplacian(7)
     1 +gw*( 
     2   -(dfdx(11)*f(13,i,j,k)-dfdx(15)*f(9,i,j,k))
     2   -(dfdy(11)*f(14,i,j,k)-dfdy(15)*f(10,i,j,k))
     2   -(dfdz(11)*f(15,i,j,k)-dfdz(15)*f(11,i,j,k))
     3   -(f(9,i,j,k)*fs(3,3,1)-f(13,i,j,k)*fs(2,3,1))
     3   -(f(10,i,j,k)*fs(3,3,2)-f(14,i,j,k)*fs(2,3,2))
     3   -(f(11,i,j,k)*fs(3,3,3)-f(15,i,j,k)*fs(2,3,3)))
     4   -gw*(g(2)*cd(3,3)-g(3)*cd(3,2))
     5   -dfdz(16)+gw*(g(17)*g(15)-g(18)*g(11))
!
! W^2_{00}:
      r(8)=0.
! W^2_{01}:
      r(9)= laplacian(9)
     1 +gw*( 
     2   -(dfdx(13)*f(5,i,j,k)-dfdx(5)*f(13,i,j,k))
     2   -(dfdy(13)*f(6,i,j,k)-dfdy(5)*f(14,i,j,k))
     2   -(dfdz(13)*f(7,i,j,k)-dfdz(5)*f(15,i,j,k))
     3   -(f(13,i,j,k)*fs(1,1,1)-f(5,i,j,k)*fs(3,1,1))
     3   -(f(14,i,j,k)*fs(1,1,2)-f(6,i,j,k)*fs(3,1,2))
     3   -(f(15,i,j,k)*fs(1,1,3)-f(7,i,j,k)*fs(3,1,3)))
     4   -gw*(g(3)*cd(1,1)-g(1)*cd(1,3))
     5   -dfdx(17)+gw*(g(18)*g(5)-g(16)*g(13))

! W^2_{02}:
      r(10)= laplacian(10)
     1 +gw*( 
     2   -(dfdx(14)*f(5,i,j,k)-dfdx(6)*f(13,i,j,k))
     2   -(dfdy(14)*f(6,i,j,k)-dfdy(6)*f(14,i,j,k))
     2   -(dfdz(14)*f(7,i,j,k)-dfdz(6)*f(15,i,j,k))
     3   -(f(13,i,j,k)*fs(1,2,1)-f(5,i,j,k)*fs(3,2,1))
     3   -(f(14,i,j,k)*fs(1,2,2)-f(6,i,j,k)*fs(3,2,2))
     3   -(f(15,i,j,k)*fs(1,2,3)-f(7,i,j,k)*fs(3,2,3)))
     4   -gw*(g(3)*cd(2,1)-g(1)*cd(2,3))
     5   -dfdy(17)+gw*(g(18)*g(6)-g(16)*g(14))
! W^2_{03}:
      r(11)= laplacian(11)
     1 +gw*( 
     2   -(dfdx(15)*f(5,i,j,k)-dfdx(7)*f(13,i,j,k))
     2   -(dfdy(15)*f(6,i,j,k)-dfdy(7)*f(14,i,j,k))
     2   -(dfdz(15)*f(7,i,j,k)-dfdz(7)*f(15,i,j,k))
     3   -(f(13,i,j,k)*fs(1,3,1)-f(5,i,j,k)*fs(3,3,1))
     3   -(f(14,i,j,k)*fs(1,3,2)-f(6,i,j,k)*fs(3,3,2))
     3   -(f(15,i,j,k)*fs(1,3,3)-f(7,i,j,k)*fs(3,3,3)))
     4   -gw*(g(3)*cd(3,1)-g(1)*cd(3,3))
     5   -dfdz(17)+gw*(g(18)*g(7)-g(16)*g(15))
! W^3_{00}:
      r(12)=0.
! W^3_{01}:
      r(13)= laplacian(13)
     1 +gw*( 
     2   -(dfdx(5)*f(9,i,j,k)-dfdx(9)*f(5,i,j,k))
     2   -(dfdy(5)*f(10,i,j,k)-dfdy(9)*f(6,i,j,k))
     2   -(dfdz(5)*f(11,i,j,k)-dfdz(9)*f(7,i,j,k))
     3   -(f(5,i,j,k)*fs(2,1,1)-f(9,i,j,k)*fs(1,1,1))
     3   -(f(6,i,j,k)*fs(2,1,2)-f(10,i,j,k)*fs(1,1,2))
     3   -(f(7,i,j,k)*fs(2,1,3)-f(11,i,j,k)*fs(1,1,3)))
     4   -gw*(g(1)*cd(1,2)-g(2)*cd(1,1))
     5   -dfdx(18)+gw*(g(16)*g(9)-g(17)*g(5))
! W^3_{02}:
      r(14)= laplacian(14)
     1 +gw*( 
     2   -(dfdx(6)*f(9,i,j,k)-dfdx(10)*f(5,i,j,k))
     2   -(dfdy(6)*f(10,i,j,k)-dfdy(10)*f(6,i,j,k))
     2   -(dfdz(6)*f(11,i,j,k)-dfdz(10)*f(7,i,j,k))
     3   -(f(5,i,j,k)*fs(2,2,1)-f(9,i,j,k)*fs(1,2,1))
     3   -(f(6,i,j,k)*fs(2,2,2)-f(10,i,j,k)*fs(1,2,2))
     3   -(f(7,i,j,k)*fs(2,2,3)-f(11,i,j,k)*fs(1,2,3)))
     4   -gw*(g(1)*cd(2,2)-g(2)*cd(2,1))
     5   -dfdy(18)+gw*(g(16)*g(10)-g(17)*g(6))
! W^3_{03}:
      r(15)= laplacian(15)
     1 +gw*( 
     2   -(dfdx(7)*f(9,i,j,k)-dfdx(11)*f(5,i,j,k))
     2   -(dfdy(7)*f(10,i,j,k)-dfdy(11)*f(6,i,j,k))
     2   -(dfdz(7)*f(11,i,j,k)-dfdz(11)*f(7,i,j,k))
     3   -(f(5,i,j,k)*fs(2,3,1)-f(9,i,j,k)*fs(1,3,1))
     3   -(f(6,i,j,k)*fs(2,3,2)-f(10,i,j,k)*fs(1,3,2))
     3   -(f(7,i,j,k)*fs(2,3,3)-f(11,i,j,k)*fs(1,3,3)))
     4   -gw*(g(1)*cd(3,2)-g(2)*cd(3,1))
     5   -dfdz(18)+gw*(g(16)*g(11)-g(17)*g(7))

! Gamma^a:
! gauge functions obey first order equations (so there is no
! need of fd(21),fd(22),fd(23))::

        r(16)=0.
        r(17)=0.
        r(18)=0.

! |phi| equation (on 3 October 2016):

      r(19)= laplacian(19)
     1-(dxhatn(1,i,j,k)**2+dxhatn(2,i,j,k)**2+dxhatn(3,i,j,k)**2
     2 +dyhatn(1,i,j,k)**2+dyhatn(2,i,j,k)**2+dyhatn(3,i,j,k)**2
     3 +dzhatn(1,i,j,k)**2+dzhatn(2,i,j,k)**2+dzhatn(3,i,j,k)**2)*g(19)
     4  -2.*gw*g(19)*
     5  (g(5)*(hatn(2,i,j,k)*dxhatn(3,i,j,k)
     6                            -hatn(3,i,j,k)*dxhatn(2,i,j,k))
     7  +g(6)*(hatn(2,i,j,k)*dyhatn(3,i,j,k)
     8                            -hatn(3,i,j,k)*dyhatn(2,i,j,k))
     9  +g(7)*(hatn(2,i,j,k)*dzhatn(3,i,j,k)
     1                            -hatn(3,i,j,k)*dzhatn(2,i,j,k))
     2  +g(9)*(hatn(3,i,j,k)*dxhatn(1,i,j,k)
     3                            -hatn(1,i,j,k)*dxhatn(3,i,j,k))
     4  +g(10)*(hatn(3,i,j,k)*dyhatn(1,i,j,k)
     5                            -hatn(1,i,j,k)*dyhatn(3,i,j,k))
     6  +g(11)*(hatn(3,i,j,k)*dzhatn(1,i,j,k)
     7                            -hatn(1,i,j,k)*dzhatn(3,i,j,k))
     8  +g(13)*(hatn(1,i,j,k)*dxhatn(2,i,j,k)
     9                            -hatn(2,i,j,k)*dxhatn(1,i,j,k))
     1  +g(14)*(hatn(1,i,j,k)*dyhatn(2,i,j,k)
     2                            -hatn(2,i,j,k)*dyhatn(1,i,j,k))
     3  +g(15)*(hatn(1,i,j,k)*dzhatn(2,i,j,k)
     4                            -hatn(2,i,j,k)*dzhatn(1,i,j,k)))
     5  -gw**2*g(19)*(g(5)**2+g(6)**2+g(7)**2+g(9)**2+g(10)**2+g(11)**2
     6                +g(13)**2+g(14)**2+g(15)**2
     7 -(hatn(1,i,j,k)*g(5)+hatn(2,i,j,k)*g(9)+hatn(3,i,j,k)*g(13))**2
     8 -(hatn(1,i,j,k)*g(6)+hatn(2,i,j,k)*g(10)+hatn(3,i,j,k)*g(14))**2
     9 -(hatn(1,i,j,k)*g(7)+hatn(2,i,j,k)*g(11)+hatn(3,i,j,k)*g(15))**2)
     1 -lambda*(g(19)**2-vev**2)*g(19)

      return
      end
