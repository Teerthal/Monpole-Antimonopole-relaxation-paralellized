      subroutine fdflux(f,i,j,k,r,lxb,lyl,lzd)
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

      dimension f(nf,1:localsize,1:localsize,1:localsize)
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
! fields for SU(1)xU(1):
!     f(1)=phi^1,f(2)=phi^2,f(3)=phi^3,f(4)=phi^4
!     f(5)=W^1_0, f(6)=W^1_1, f(7)=W^1_2, f(8)=W^1_3,
!     f(9)=W^2_0, f(10)=W^2_1, f(11)=W^2_2, f(12)=W^2_3,
!     f(13)=W^3_0, f(14)=W^3_1, f(15)=W^3_2, f(16)=W^3_3,
!     f(17)=B_0, f(18)=B_1, f(19)=B_2, f(20)=B_3 
!     gauge fixing functions (Num Rel): f(21)=\partial_i W^1_i, 
!     f(22)=\partial_i W^2_i, f(23)=\partial_i W^3_i, 
!     f(24)=\partial_iB_i
!     time derivatives:
!     fd(i)=time derivative of f(i) for i=1-4 (scalar fields),
!     fd(i)=electric fields(=W^a_{0j}) for i=5-20 (gauge fields).
!     fd(24)=time derivative of |phi|
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

          laplacian(1)=(gw/2.)*((-g(6)*dfdx(4)-g(7)*dfdy(4)-g(8)*dfdz(4))
     1                     +(g(10)*dfdx(3)+g(11)*dfdy(3)+g(12)*dfdz(3))
     2                     -(g(14)*dfdx(2)+g(15)*dfdy(2)+g(16)*dfdz(2))
     3                     -(g(21)*g(4)+g(22)*g(3)+g(23)*g(2)))
     4                     +(gy/2.)*(
     5                     -(g(18)*dfdx(2)+g(19)*dfdy(2)+g(20)*dfdz(2))
     6                     -(g(24)*g(2)))

          laplacian(2)=(gw/2.)*((g(6)*dfdx(1)+g(7)*dfdy(1)+g(8)*dfdz(1))
     1                     +(g(10)*dfdx(4)+g(11)*dfdy(4)+g(12)*dfdz(4))
     2                     +(g(14)*dfdx(1)+g(15)*dfdy(1)+g(16)*dfdz(1))
     3                     +(g(21)*g(3)+g(22)*g(3)+g(23)*g(1)))
     4                     +(gy/2.)*(
     5                     +(g(18)*dfdx(1)+g(19)*dfdy(1)+g(20)*dfdz(1))
     6                     -(g(24)*g(1)))

          laplacian(3)=(gw/2.)*(-(g(6)*dfdx(2)+g(7)*dfdy(2)+g(8)*dfdz(2))
     1                     -(g(10)*dfdx(1)+g(11)*dfdy(1)+g(12)*dfdz(1))
     2                     +(g(14)*dfdx(4)+g(15)*dfdy(4)+g(16)*dfdz(4))
     3                     +(-g(21)*g(2)-g(22)*g(1)+g(23)*g(4)))
     4                     +(gy/2.)*(
     5                     -(g(18)*dfdx(4)+g(19)*dfdy(4)+g(20)*dfdz(4))
     6                     -(g(24)*g(4)))

          laplacian(4)=(gw/2.)*((g(6)*dfdx(1)+g(7)*dfdy(1)+g(8)*dfdz(1))
     1                     -(g(10)*dfdx(2)+g(11)*dfdy(2)+g(12)*dfdz(2))
     2                     -(g(14)*dfdx(3)+g(15)*dfdy(3)+g(16)*dfdz(3))
     3                     +(g(21)*g(1)-g(22)*g(2)-g(23)*g(3)))
     4                     +(gy/2.)*(
     5                     -(g(18)*dfdx(4)+g(19)*dfdy(4)+g(20)*dfdz(4))
     6                     -(g(24)*g(4)))

! gauge fields;
! W^1_i
          laplacian(5)=0.

          laplacian(6)=-dfdx(21)+gw*((dfdx(10)*g(14)+dfdy(10)*g(15)
     1                 +dfdz(10)*g(16))-(dfdx(14)*g(10)+dfdy(14)*g(11)
     2                 +dfdz(14)*g(12))+gy*((g(10)*g(23))-(g(14)*g(22))))

          laplacian(7)=-dfdy(21)+gw*((dfdx(11)*g(14)+dfdy(11)*g(15)
     1                 +dfdz(11)*g(16))-(dfdx(15)*g(10)+dfdy(15)*g(11)
     2                 +dfdz(15)*g(12))+gy*((g(11)*g(23))-(g(15)*g(22))))

          laplacian(8)=-dfdz(21)+gw*((dfdx(12)*g(14)+dfdy(12)*g(15)
     1                 +dfdz(12)*g(16))-(dfdx(16)*g(10)+dfdy(16)*g(11)
     2                 +dfdz(16)*g(12))+gy*((g(12)*g(23))-(g(16)*g(22))))

! W^2_i
          laplacian(9)=0.

          laplacian(10)=-dfdx(22)+gw*((dfdx(14)*g(6)+dfdy(14)*g(7)
     1                 +dfdz(14)*g(8))-(dfdx(6)*g(14)+dfdy(6)*g(15)
     2                 +dfdz(6)*g(16))+gy*((g(14)*g(21))-(g(6)*g(23))))

          laplacian(11)=-dfdy(22)+gw*((dfdx(15)*g(6)+dfdy(15)*g(7)
     1                 +dfdz(15)*g(8))-(dfdx(7)*g(14)+dfdy(7)*g(15)
     2                 +dfdz(7)*g(16))+gy*((g(15)*g(21))-(g(7)*g(23))))

          laplacian(12)=-dfdz(22)+gw*((dfdx(16)*g(6)+dfdy(16)*g(7)
     1                 +dfdz(16)*g(8))-(dfdx(8)*g(14)+dfdy(8)*g(15)
     2                 +dfdz(8)*g(16))+gy*((g(16)*g(21))-(g(8)*g(24))))

! W^3_i
          laplacian(13)=0.

          laplacian(14)=-dfdx(23)+gw*((dfdx(6)*g(10)+dfdy(6)*g(11)
     1                 +dfdz(6)*g(12))-(dfdx(10)*g(6)+dfdy(10)*g(7)
     2                 +dfdz(10)*g(8))+gy*((g(6)*g(22))-(g(10)*g(21))))

          laplacian(15)=-dfdy(23)+gw*((dfdx(7)*g(10)+dfdy(7)*g(11)
     1                 +dfdz(7)*g(12))-(dfdx(11)*g(6)+dfdy(11)*g(7)
     2                 +dfdz(11)*g(8))+gy*((g(7)*g(22))-(g(11)*g(21))))

          laplacian(16)=-dfdz(23)+gw*((dfdx(8)*g(10)+dfdy(8)*g(11)
     1                 +dfdz(8)*g(12))-(dfdx(12)*g(6)+dfdy(12)*g(7)
     2                 +dfdz(12)*g(8))+gy*((g(8)*g(22))-(g(12)*g(21))))

!B_i
          laplacian(17)=0.

          laplacian(18)=dfdx(24)

          laplacian(19)=dfdy(24)

          laplacian(20)=dfdz(24)

! gauge constraints:

          laplacian(21)=0.
          laplacian(22)=0.
          laplacian(23)=0. 
          laplacian(24)=0.


!------------Omitting phi magnitude for now----------------------------!          
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

!          laplacian(19)=-sqrt(2.*lambda)*g(19)*(hatn(1,i,j,k)*dfdx(19)
!     1                 +hatn(2,i,j,k)*dfdy(19)+hatn(3,i,j,k)*dfdz(19))
!     2          +sqrt(lambda/2.)*(vev**2-g(19)**2)
!     3              *(dxhatn(1,i,j,k)+dyhatn(2,i,j,k)+dzhatn(3,i,j,k))

          enddo
!------------'End' Omitting phi magnitude for now----------------------!


! if not on boundaries:
      else
          do n=1,nf
          laplacian(n)=d2fdx2(n)+d2fdy2(n)+d2fdz2(n)
          enddo
      endif

!=======================================================================
! Scalar fluxes:
!=======================================================================

! phi^1:
      r(1)= laplacian(1)-(gw/2.)*((-g(6)*dfdx(4)-g(7)*dfdy(4)-g(8)*dfdz(4))
     1      +(g(10)*dfdx(3)+g(11)*dfdy(3)+g(12)*dfdz(3))-(g(14)*dfdx(2)
     2      +g(15)*dfdy(2)+g(16)*dfdz(2))-(g(21)*g(4)+g(22)*g(3)+g(23)*g(2)))
     3      -(gy/2.)*(-(g(18)*dfdx(2)+g(19)*dfdy(2)+g(20)*dfdz(2))
     4      -(g(24)*g(2)))-(gw/2.)*(-(g(6)*cd(1,4)+g(7)*cd(2,4)+g(8)*cd(3,4))
     5      +(g(10)*cd(1,4)+g(11)*cd(2,4)+g(12)*cd(3,4))-(g(14)*cd(1,2)               
     6      +g(15)*cd(2,2)+g(16)*cd(3,2)))+(gy/2.)*(f(18)*cd(1,2)
     7      +f(19)*cd(2,2)+f(20)*cd(3,2))-2*lambda*(g(1)**2+g(2)**2+g(3)**2
     8      +g(4)**2-vev**2)*g(1)

! phi^2:
      r(2)= laplacian(2)-(gw/2.)*((g(6)*dfdx(1)+g(7)*dfdy(1)+g(8)*dfdz(1))
     1      +(g(10)*dfdx(4)+g(11)*dfdy(4)+g(12)*dfdz(4))+(g(14)*dfdx(1)
     2      +g(15)*dfdy(1)+g(16)*dfdz(1))+(g(21)*g(3)+g(22)*g(3)+g(23)*g(1)))
     3      -(gy/2.)*(+(g(18)*dfdx(1)+g(19)*dfdy(1)+g(20)*dfdz(1))
     4      -(g(24)*g(1)))-(gw/2.)*( (g(6)*cd(1,1)+g(7)*cd(2,1)+g(8)*cd(3,1))
     5      +(g(10)*cd(1,4)+g(11)*cd(2,4)+g(12)*cd(3,4))+(g(14)*cd(1,1)               
     6      +g(15)*cd(2,1)+g(16)*cd(3,1)))-(gy/2.)*(f(18)*cd(1,1)
     7      +f(19)*cd(2,1)+f(20)*cd(3,1))-2*lambda*(g(1)**2+g(2)**2+g(3)**2
     8      +g(4)**2-vev**2)*g(2)
                    

! phi^3:
      r(3)= laplacian(3)-(gw/2.)*(-(g(6)*dfdx(2)+g(7)*dfdy(2)+g(8)*dfdz(2))
     1      -(g(10)*dfdx(1)+g(11)*dfdy(1)+g(12)*dfdz(1))+(g(14)*dfdx(4)
     2      +g(15)*dfdy(4)+g(16)*dfdz(4))+(-g(21)*g(2)-g(22)*g(1)
     3      +g(23)*g(4)))-(gy/2.)*(-(g(18)*dfdx(4)+g(19)*dfdy(4)
     4      +g(20)*dfdz(4))-(g(24)*g(4)))-(gw/2.)*(-(g(6)*cd(1,2)
     5      +g(7)*cd(2,2)+g(8)*cd(3,2))-(g(10)*cd(1,1)+g(11)*cd(2,1)               
     6      +g(12)*cd(3,1))+(g(14)*cd(1,4)+g(15)*cd(2,4)+g(16)*cd(3,4)))
     7      +(gy/2.)*(f(18)*cd(1,4)+f(19)*cd(2,4)+f(20)*cd(3,4))-2*lambda
     8      *(g(1)**2+g(2)**2+g(3)**2+g(4)**2-vev**2)*g(3)

! phi^4
      r(4) = laplacian(4)-(gw/2.)*((g(6)*dfdx(1)+g(7)*dfdy(1)+g(8)*dfdz(1))
     1       -(g(10)*dfdx(2)+g(11)*dfdy(2)+g(12)*dfdz(2))-(g(14)*dfdx(3)
     2       +g(15)*dfdy(3)+g(16)*dfdz(3))+(g(21)*g(1)-g(22)*g(2)
     3       -g(23)*g(3)))+(gy/2.)*(-(g(18)*dfdx(4)+g(19)*dfdy(4)            
     4       +g(20)*dfdz(4))-(g(24)*g(4)))-(gw/2.)*((g(6)*cd(1,1)
     5      +g(7)*cd(2,1)+g(8)*cd(3,1))-(g(10)*cd(1,2)+g(11)*cd(2,2)               
     6      +g(12)*cd(3,2))-(g(14)*cd(1,3)+g(15)*cd(2,3)+g(16)*cd(3,3)))
     7      -(gy/2.)*(f(18)*cd(1,3)+f(19)*cd(2,3)+f(20)*cd(3,3))-2*lambda
     8      *(g(1)**2+g(2)**2+g(3)**2+g(4)**2-vev**2)*g(4)


! W-fluxes (assumes $W^a_0=0$ gauge):
!
! W^1_{0}:
      r(5)=0.

! W^1_{1}:
      r(6)= laplacian(6)-(-dfdx(21)+gw*((dfdx(10)*g(14)+dfdy(10)*g(15)
     1      +dfdz(10)*g(16))-(dfdx(14)*g(10)+dfdy(14)*g(11)+dfdz(14)*g(12))
     2      +gy*((g(10)*g(23))-(g(14)*g(22)))))+gw*(g(10)*fs(3,1,1)
     3      +g(11)*fs(3,1,2)+g(12)*fs(3,1,3)-(g(14)*fs(2,1,1)
     4      +g(15)*fs(2,1,2)+g(16)*fs(2,1,3)))+gw*(g(1)*cd(4,1)-g(2)*cd(3,1)
     5      +g(3)*cd(3,1)-g(4)*cd(1,1))

! W^1_{2}:
      r(7)= laplacian(7)-(-dfdy(21)+gw*((dfdx(11)*g(14)+dfdy(11)*g(15)
     1      +dfdz(11)*g(16))-(dfdx(15)*g(10)+dfdy(15)*g(11)+dfdz(15)*g(12))
     2      +gy*((g(11)*g(23))-(g(15)*g(22)))))+gw*(g(10)*fs(3,2,1)
     3      +g(11)*fs(3,2,2)+g(12)*fs(3,2,3)-(g(14)*fs(2,2,1)
     4      +g(15)*fs(2,2,2)+g(16)*fs(2,2,3)))+gw*(g(1)*cd(4,2)-g(2)*cd(3,2)
     5      +g(3)*cd(3,2)-g(4)*cd(1,2))

! W^1_{3}:
      r(8)= laplacian(8)-(-dfdz(21)+gw*((dfdx(12)*g(14)+dfdy(12)*g(15)
     1      +dfdz(12)*g(16))-(dfdx(16)*g(10)+dfdy(16)*g(11)+dfdz(16)*g(12))
     2      +gy*((g(12)*g(23))-(g(16)*g(22)))))+gw*(g(10)*fs(3,3,1)
     3      +g(11)*fs(3,3,2)+g(12)*fs(3,3,3)-(g(14)*fs(2,3,1)
     4      +g(15)*fs(2,3,2)+g(16)*fs(2,3,3)))+gw*(g(1)*cd(4,3)-g(2)*cd(3,3)
     5      +g(3)*cd(3,3)-g(4)*cd(1,3))
!
! W^2_{0}:
      r(9)=0.

! W^2_{1}:
      r(10)= laplacian(10)-(-dfdx(22)+gw*((dfdx(14)*g(6)+dfdy(14)*g(7)
     1       +dfdz(14)*g(8))-(dfdx(6)*g(14)+dfdy(6)*g(15)+dfdz(6)*g(16))
     2       +gy*((g(14)*g(21))-(g(6)*g(23)))))+gw*(g(14)*fs(1,1,1)
     3      +g(15)*fs(1,1,2)+g(16)*fs(1,1,3)-(g(6)*fs(3,1,1)
     4      +g(7)*fs(3,1,2)+g(8)*fs(3,1,3)))+gw*(-g(1)*cd(3,1)
     5      -g(2)*cd(4,1)+g(3)*cd(1,1)+g(4)*cd(2,1))

! W^2_{2}:
      r(11)= laplacian(11)-(-dfdy(22)+gw*((dfdx(15)*g(6)+dfdy(15)*g(7)
     1       +dfdz(15)*g(8))-(dfdx(7)*g(14)+dfdy(7)*g(15)+dfdz(7)*g(16))
     2       +gy*((g(15)*g(21))-(g(7)*g(23)))))+gw*(g(14)*fs(1,2,1)
     3      +g(15)*fs(1,2,2)+g(16)*fs(1,2,3)-(g(6)*fs(3,2,1)
     4      +g(7)*fs(3,2,2)+g(8)*fs(3,2,3)))+gw*(-g(1)*cd(3,2)
     5      -g(2)*cd(4,2)+g(3)*cd(1,2)+g(4)*cd(2,2))

! W^2_{3}:
      r(12)= laplacian(12)-(-dfdz(22)+gw*((dfdx(16)*g(6)+dfdy(16)*g(7)
     1       +dfdz(16)*g(8))-(dfdx(8)*g(14)+dfdy(8)*g(15)+dfdz(8)*g(16))
     2       +gy*((g(16)*g(21))-(g(8)*g(24)))))+gw*(g(14)*fs(1,3,1)
     3      +g(15)*fs(1,3,2)+g(16)*fs(1,3,3)-(g(6)*fs(3,3,1)
     4      +g(7)*fs(3,3,2)+g(8)*fs(3,3,3)))+gw*(-g(1)*cd(3,3)
     5      -g(2)*cd(4,3)+g(3)*cd(1,3)+g(4)*cd(2,3))

! W^3_{0}:
      r(13)=0.

! W^3_{1}:
      r(14)= laplacian(13)-(-dfdx(23)+gw*((dfdx(6)*g(10)+dfdy(6)*g(11)
     1       +dfdz(6)*g(12))-(dfdx(10)*g(6)+dfdy(10)*g(7)+dfdz(10)*g(8))
     2       +gy*((g(6)*g(22))-(g(10)*g(21)))))+gw*(g(6)*fs(2,1,1)
     3      +g(7)*fs(2,1,2)+g(8)*fs(2,1,3)-(g(10)*fs(1,1,1)
     4      +g(11)*fs(1,1,2)+g(12)*fs(1,1,3)))+gw*(g(1)*cd(2,1)
     5      -g(2)*cd(1,1)-g(3)*cd(4,1)+g(4)*cd(3,1))

! W^3_{2}:
      r(15)= laplacian(15)-(-dfdy(23)+gw*((dfdx(7)*g(10)+dfdy(7)*g(11)
     1       +dfdz(7)*g(12))-(dfdx(11)*g(6)+dfdy(11)*g(7)+dfdz(11)*g(8))
     2       +gy*((g(7)*g(22))-(g(11)*g(21)))))+gw*(g(6)*fs(2,2,1)
     3      +g(7)*fs(2,2,2)+g(8)*fs(2,2,3)-(g(10)*fs(1,2,1)
     4      +g(11)*fs(1,2,2)+g(12)*fs(1,2,3)))+gw*(g(1)*cd(2,2)
     5      -g(2)*cd(1,2)-g(3)*cd(4,2)+g(4)*cd(3,2))

! W^3_{3}:
      r(16)= laplacian(16)-(-dfdz(23)+gw*((dfdx(8)*g(10)+dfdy(8)*g(11)
     1       +dfdz(8)*g(12))-(dfdx(12)*g(6)+dfdy(12)*g(7)+dfdz(12)*g(8))
     2       +gy*((g(8)*g(22))-(g(12)*g(21)))))+gw*(g(6)*fs(2,3,1)
     3      +g(7)*fs(2,3,2)+g(8)*fs(2,3,3)-(g(10)*fs(1,3,1)
     4      +g(11)*fs(1,3,2)+g(12)*fs(1,3,3)))+gw*(g(1)*cd(2,3)
     5      -g(2)*cd(1,3)-g(3)*cd(4,3)+g(4)*cd(3,3))

! B_0:
      r(17)-0.

! B_1:
      r(18)=laplacian(18)-(dfdx(24))+gy*((cd(1,2)*f(1))-(cd(1,4)*f(2))
     1      +(cd(1,4)*f(3))-(cd(1,3)*f(4)))

! B_2:
      r(19)=laplacian(19)-(dfdx(24))+gy*((cd(1,2)*f(1))-(cd(1,4)*f(2))
     1      +(cd(1,4)*f(3))-(cd(1,3)*f(4)))

! B_3:
      r(20)=laplacian(20)-(dfdx(24))+gy*((cd(1,2)*f(1))-(cd(1,4)*f(2))
     1      +(cd(1,4)*f(3))-(cd(1,3)*f(4)))

