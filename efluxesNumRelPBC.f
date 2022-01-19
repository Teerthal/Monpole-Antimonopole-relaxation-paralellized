      subroutine efdfluxPBC(f,fd,i,j,k,r,lxb,lyl,lzd)
      implicit none
      include 'parameters.inc'
      integer i,j,k,n,nn,lxb,lyl,lzd,ilocal,jlocal,klocal
      real*8 f,fd,r
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 cd,g,gd
      real*8 x,y,z,rr
      real*8 rmserror
      real*8 dfddx,dfddy,dfddz
      real*8 res
      real*8 px,py,pz,dist
      real*8 laplacian
c
        dimension f(nf,1:localsize,1:localsize,1:localsize)
        dimension fd(nf,1:localsize,1:localsize,1:localsize)

      dimension r(nf)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
c fs=field strength; first index is the number of gauge fields.
      dimension fs(ngauge/4,0:3,0:3),d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
c cd=covariant derivative; second index is the gauge field index.
      dimension cd(0:3,nscalar),g(nf),gd(nf)
      dimension res(nf)
      dimension laplacian(nf)
c
c fields:
c     f(1)=Real(phi(up)),f(2)=Imag(phi(down)),f(3)=Real(phi(down)),
c     f(4)=Imag(phi(down)).
c     f(5)=W^1_0, f(6)=W^1_1, f(7)=W^1_2, f(8)=W^1_3,
c     f(9)=W^2_0, f(10)=W^2_1, f(11)=W^2_2, f(12)=W^2_3,
c     f(13)=W^3_0, f(14)=W^3_1, f(15)=W^3_2, f(16)=W^3_3,
c     f(17)=Y_0, f(18)=Y_1, f(19)=Y_2, f(20)=Y_3;
c     gauge fixing functions (Num Rel): 
c     f(21)=\partial_i W^1_i, f(22)=\partial_i W^2_i, 
c     f(23)=\partial_i W^3_i, f(24)=\partial_i Y_i. 
c     time derivatives:
c     fd(i)=time derivative of f(i) for i=1-4 (scalar fields),
c     fd(i)=electric fields(=W^a_{0j}) for i=4-20 (gauge fields).
c

            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

c spatial derivatives to 4th order:

      call ederivativesPBC(f,fd,i,j,k,dfdx,dfdy,dfdz,dfddx,dfddy,dfddz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

c covariant derivatives:

      call ecovderiv(f,fd,i,j,k,dfdx,dfdy,dfdz,cd)

c find gauge field strengths (uses derivatives):

      call efieldstrength(f,fd,i,j,k,dfdx,dfdy,dfdz,fs)

c for convenience of writing:
      do n=1,nf
      g(n)=f(n,i,j,k)
      gd(n)=fd(n,i,j,k)
      enddo
c
c if on boundaries:
!        if(abs(ilocal).eq.latx.or.abs(jlocal)
!     1                          .eq.laty.or.abs(klocal).eq.latz) then
c radial unit vector:
!        px=float(ilocal)/sqrt(float(ilocal**2+jlocal**2+klocal**2))
!        py=float(jlocal)/sqrt(float(ilocal**2+jlocal**2+klocal**2))
!        pz=float(klocal)/sqrt(float(ilocal**2+jlocal**2+klocal**2))
!        dist=sqrt(float(ilocal**2+jlocal**2+klocal**2))*dx

!         do n=1,nf
c see README file for explanation of this formula for the
c laplacian on the boundary:

!          laplacian(n)=-(px*dfddx(n)+py*dfddy(n)+pz*dfddz(n))

c
!         enddo

c if not on boundaries:

!      else
          do n=1,nf
          laplacian(n)=d2fdx2(n)+d2fdy2(n)+d2fdz2(n)
          enddo
!      endif
c
c Higgs fluxes:
c
c equation for scalar is:
c $$
c \partial_t^2\phi^a = \nabla^2\phi^a - g \epsilon^{abc}\partial_i\phi^b W_i^c
c       - g \epsilon^{abc} (D_i \phi)^b W_i^c - \lambda (\phi^b\phi^b-vev^2)\phi^a
c       - g \epsilon^{abc} \phi^b \Gamma^c
c $$
c where $\Gamma^a=\partial_i W_i^a$ and we have chosen $W_0^a=0$.
c
      r(1)=laplacian(1)
c $- g \epsilon^{abc}\partial_i\phi^b W_i^c$ terms:
     1  -gw*((dfdx(2)*g(13)-dfdx(3)*g(9))+(dfdy(2)*g(14)-dfdy(3)*g(10))
     2       +(dfdz(2)*g(15)-dfdz(3)*g(11)))
c $- g \epsilon^{abc} (D_i \phi)^b W_i^c$ terms:
     3  -gw*((cd(1,2)*g(13)-cd(1,3)*g(9))+(cd(2,2)*g(14)-cd(2,3)*g(10))
     4       +(cd(3,2)*g(15)-cd(3,3)*g(11)))
c $- \lambda (\phi^b\phi^b-vev^2)\phi^a$ terms:
     5  -lambda*(g(1)**2+g(2)**2+g(3)**2-vev**2)*g(1)
c $- g \epsilon^{abc} \phi^b \Gamma^c$ terms:
     6  -gw*(g(2)*g(18)-g(3)*g(17))
c
      r(2)=laplacian(2)
     1  -gw*((dfdx(3)*g(5)-dfdx(1)*g(13))+(dfdy(3)*g(6)-dfdy(1)*g(14))
     2     +(dfdz(3)*g(7)-dfdz(1)*g(15)))
     3  -gw*((cd(1,3)*g(5)-cd(1,1)*g(13))+(cd(2,3)*g(6)-cd(2,1)*g(14))
     4     +(cd(3,3)*g(7)-cd(3,1)*g(15)))
     5  -lambda*(g(1)**2+g(2)**2+g(3)**2-vev**2)*g(2)
     6  -gw*(g(3)*g(16)-g(1)*g(18))
c
      r(3)=laplacian(3)
     1  -gw*((dfdx(1)*g(9)-dfdx(2)*g(5))+(dfdy(1)*g(10)-dfdy(2)*g(6))
     2     +(dfdz(1)*g(11)-dfdz(2)*g(7)))
     3  -gw*((cd(1,1)*g(9)-cd(1,2)*g(5))+(cd(2,1)*g(10)-cd(2,2)*g(6))
     4     +(cd(3,1)*g(11)-cd(3,2)*g(7)))
     5  -lambda*(g(1)**2+g(2)**2+g(3)**2-vev**2)*g(3)
     6  -gw*(g(1)*g(17)-g(2)*g(16))
c
c gauge field fluxes:
c
      r(4)=0.
c
      r(5)=laplacian(5)
c following is: $+g \epsilon^{abc} W^b_j \partial_j W^c_i$:
     1 +gw*( 
     2   -(dfdx(9)*g(13)-dfdx(13)*g(9))
     2   -(dfdy(9)*g(14)-dfdy(13)*g(10))
     2   -(dfdz(9)*g(15)-dfdz(13)*g(11))
c following is: $-g \epsilon^{abc} W^b_j W^c_{ij}$:
     3   -(g(9)*fs(3,1,1)-g(13)*fs(2,1,1))
     3   -(g(10)*fs(3,1,2)-g(14)*fs(2,1,2))
     3   -(g(11)*fs(3,1,3)-g(15)*fs(2,1,3)))
c gauge function $\Gamma^a$ terms:
     4   -dfdx(16)-gw*(g(9)*g(18)-g(13)*g(17))
c current is: $j^a_\mu = -g \epsilon^{abc} \phi^b (D_\mu \phi)^c$
     5   -gw*(g(2)*cd(1,3)-g(3)*cd(1,2))
c
      r(6)=laplacian(6)+gw*( 
     2   -(dfdx(10)*g(13)-dfdx(14)*g(9))
     2   -(dfdy(10)*g(14)-dfdy(14)*g(10))
     2   -(dfdz(10)*g(15)-dfdz(14)*g(11))
     3   -(g(9)*fs(3,2,1)-g(13)*fs(2,2,1))
     3   -(g(10)*fs(3,2,2)-g(14)*fs(2,2,2))
     3   -(g(11)*fs(3,2,3)-g(15)*fs(2,2,3)))
     4   -dfdy(16)-gw*(g(10)*g(18)-g(14)*g(17))
     5   -gw*(g(2)*cd(2,3)-g(3)*cd(2,2))
c
      r(7)=laplacian(7)+gw*( 
     2   -(dfdx(11)*g(13)-dfdx(15)*g(9))
     2   -(dfdy(11)*g(14)-dfdy(15)*g(10))
     2   -(dfdz(11)*g(15)-dfdz(15)*g(11))
     3   -(g(9)*fs(3,3,1)-g(13)*fs(2,3,1))
     3   -(g(10)*fs(3,3,2)-g(14)*fs(2,3,2))
     3   -(g(11)*fs(3,3,3)-g(15)*fs(2,3,3)))
     4   -dfdz(16)-gw*(g(11)*g(18)-g(15)*g(17))
     5   -gw*(g(2)*cd(3,3)-g(3)*cd(3,2))
c
      r(8)=0.
c
      r(9)=laplacian(9)+gw*( 
     2   -(dfdx(13)*g(5)-dfdx(5)*g(13))
     2   -(dfdy(13)*g(6)-dfdy(5)*g(14))
     2   -(dfdz(13)*g(7)-dfdz(5)*g(15))
     3   -(g(13)*fs(1,1,1)-g(5)*fs(3,1,1))
     3   -(g(14)*fs(1,1,2)-g(6)*fs(3,1,2))
     3   -(g(15)*fs(1,1,3)-g(7)*fs(3,1,3)))
     4   -dfdx(17)-gw*(g(13)*g(16)-g(5)*g(18))
     5   -gw*(g(3)*cd(1,1)-g(1)*cd(1,3))
c
      r(10)=laplacian(10)+gw*( 
     2   -(dfdx(14)*g(5)-dfdx(6)*g(13))
     2   -(dfdy(14)*g(6)-dfdy(6)*g(14))
     2   -(dfdz(14)*g(7)-dfdz(6)*g(15))
     3   -(g(13)*fs(1,2,1)-g(5)*fs(3,2,1))
     3   -(g(14)*fs(1,2,2)-g(6)*fs(3,2,2))
     3   -(g(15)*fs(1,2,3)-g(7)*fs(3,2,3)))
     4   -dfdy(17)-gw*(g(14)*g(16)-g(6)*g(18))
     5   -gw*(g(3)*cd(2,1)-g(1)*cd(2,3))
c
      r(11)=laplacian(11)+gw*( 
     2   -(dfdx(15)*g(5)-dfdx(7)*g(13))
     2   -(dfdy(15)*g(6)-dfdy(7)*g(14))
     2   -(dfdz(15)*g(7)-dfdz(7)*g(15))
     3   -(g(13)*fs(1,3,1)-g(5)*fs(3,3,1))
     3   -(g(14)*fs(1,3,2)-g(6)*fs(3,3,2))
     3   -(g(15)*fs(1,3,3)-g(7)*fs(3,3,3)))
     4   -dfdz(17)-gw*(g(15)*g(16)-g(7)*g(18))
     5   -gw*(g(3)*cd(3,1)-g(1)*cd(3,3))
c
      r(12)=0.
c
      r(13)=laplacian(13)+gw*( 
     2   -(dfdx(5)*g(9)-dfdx(9)*g(5))
     2   -(dfdy(5)*g(10)-dfdy(9)*g(6))
     2   -(dfdz(5)*g(11)-dfdz(9)*g(7))
     3   -(g(5)*fs(2,1,1)-g(9)*fs(1,1,1))
     3   -(g(6)*fs(2,1,2)-g(10)*fs(1,1,2))
     3   -(g(7)*fs(2,1,3)-g(11)*fs(1,1,3)))
     4   -dfdx(18)-gw*(g(5)*g(17)-g(9)*g(16))
     5   -gw*(g(1)*cd(1,2)-g(2)*cd(1,1))
c
      r(14)=laplacian(14)+gw*( 
     2   -(dfdx(6)*g(9)-dfdx(10)*g(5))
     2   -(dfdy(6)*g(10)-dfdy(10)*g(6))
     2   -(dfdz(6)*g(11)-dfdz(10)*g(7))
     3   -(g(5)*fs(2,2,1)-g(9)*fs(1,2,1))
     3   -(g(6)*fs(2,2,2)-g(10)*fs(1,2,2))
     3   -(g(7)*fs(2,2,3)-g(11)*fs(1,2,3)))
     4   -dfdy(18)-gw*(g(6)*g(17)-g(10)*g(16))
     5   -gw*(g(1)*cd(2,2)-g(2)*cd(2,1))
c
      r(15)=laplacian(15)+gw*( 
     2   -(dfdx(7)*g(9)-dfdx(11)*g(5))
     2   -(dfdy(7)*g(10)-dfdy(11)*g(6))
     2   -(dfdz(7)*g(11)-dfdz(11)*g(7))
     3   -(g(5)*fs(2,3,1)-g(9)*fs(1,3,1))
     3   -(g(6)*fs(2,3,2)-g(10)*fs(1,3,2))
     3   -(g(7)*fs(2,3,3)-g(11)*fs(1,3,3)))
     4   -dfdz(18)-gw*(g(7)*g(17)-g(11)*g(16))
     5   -gw*(g(1)*cd(3,2)-g(2)*cd(3,1))
c
c gauge constraint functions obey first order equations (so there is no
c need of fd(16)-fd(18)):
        r(16)=0.
        r(17)=0.
        r(18)=0.
c
      return
      end
