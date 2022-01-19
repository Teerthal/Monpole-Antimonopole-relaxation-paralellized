c 21 June 2016:
c This is for SO(3) model.
c
c includes absorbing boundary conditions right now --
c may be better to move these to a separate subroutine.
c
      subroutine effluxes(f,fd,i,j,k,s,lxb,lyl,lzd)
      implicit none
      include 'parameters.inc'
      integer i,j,k,n,nn,lxb,lyl,lzd,ilocal,jlocal,klocal
      real*8 f,fd,s
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 cd,g,gd
      real*8 x,y,z,rr
      real*8 dfddx,dfddy,dfddz
      real*8 px,py,pz
c
        dimension f(nf,1:localsize,1:localsize,1:localsize)
        dimension fd(nf,1:localsize,1:localsize,1:localsize)
      dimension s(nf)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
c fs=field strength; first index is the number of gauge fields.
      dimension fs(ngauge/4,0:3,0:3),d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
c cd=covariant derivative; second index is the gauge field index.
      dimension cd(0:3,nscalar),g(nf),gd(nf)
c
c fields:
c     f(1)=Real(phi(up)),f(2)=Imag(phi(down)),f(3)=Real(phi(down)),
c     f(4)=Imag(phi(down)).
c     f(5)=W^1_0, f(6)=W^1_1, f(7)=W^1_2, f(8)=W^1_3,
c     f(9)=W^2_0, f(10)=W^2_1, f(11)=W^2_2, f(12)=W^2_3,
c     f(13)=W^3_0, f(14)=W^3_1, f(15)=W^3_2, f(16)=W^3_3,
c     f(17)=Y_0, f(18)=Y_1, f(19)=Y_2, f(20)=Y_3;
c     f(21)=\partial_\mu W^{1\mu}, f(22)=\partial_\mu W^{2\mu},
c     f(23)=\partial_\mu W^{3\mu}, f(24)=\partial_\mu Y^\mu;
c     fd(i)=time derivative of f(i)
c

            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

c spatial derivatives to 4th order:

      call ederivatives(f,fd,i,j,k,dfdx,dfdy,dfdz,dfddx,dfddy,dfddz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

c covariant derivatives:

      call ecovderiv(f,fd,i,j,k,dfdx,dfdy,dfdz,cd)

c find gauge field strengths (uses derivatives):

      call efieldstrength(f,fd,i,j,k,dfdx,dfdy,dfdz,fs)
c
c The equation we are solving here is $\partial_t (f)=fd$
c and so there are no spatial derivatives on the rhs. Therefore
c we don't need any absorbing boundary conditions for this equation.
c
c for the scalar field, $fd()=D_0\Phi$ which is the same as
c $fd()=\partial_0\Phi$ in the $W^a_0=0=Y_0$ gauge (which we
c are adopting). 
       s(1)=fd(1,i,j,k)
       s(2)=fd(2,i,j,k)
       s(3)=fd(3,i,j,k)
c
c Eq. (2.11) of Baumgarte&Shapiro is $\partial_t A_i = -E_i -...$ so 
c we are taking fd(...)=+\partial_t A_i = -E_i (note the sign).
       s(4)=0.
       s(5)=fd(5,i,j,k)
       s(6)=fd(6,i,j,k)
       s(7)=fd(7,i,j,k)
c
       s(8)=0.
       s(9)=fd(9,i,j,k)
       s(10)=fd(10,i,j,k)
       s(11)=fd(11,i,j,k)
c
       s(12)=0.
       s(13)=fd(13,i,j,k)
       s(14)=fd(14,i,j,k)
       s(15)=fd(15,i,j,k)
c
c fluxes for gauge functions:
c The equation for the gauge functions are:
c \partial_t g(16) = \partial_i (W^a_{0i}) 
c                    - gp2*(\partial_i (W^a_{0i})+g\epsilon^{abc}W^b_i W^c_{0i} -j^a_0)
       s(16)=
c This line is $(1-gp2)* \partial_i (W^a_{0i})$ Note: we are working in the $W^a_0=0=Y_0$ gauge.
     1  (1.-gp2)*(dfddx(5)+dfddy(6)+dfddz(7))
c This line is: - gp2*(g\epsilon^{abc}W^b_i W^c_{0i})
     2    +gp2*gw*(
     3   -(f(9,i,j,k)*fs(3,0,1)-f(13,i,j,k)*fs(2,0,1))
     3   -(f(10,i,j,k)*fs(3,0,2)-f(14,i,j,k)*fs(2,0,2))
     3   -(f(11,i,j,k)*fs(3,0,3)-f(15,i,j,k)*fs(2,0,3)))
c This line is: -gp2*(-j^a_0) with j^a_0 = -g\epsilon^{abc}\phi^b (D_t \phi)^c
     4   -gp2*gw*(f(2,i,j,k)*cd(0,3)-f(3,i,j,k)*cd(0,2))
c
       s(17)=
     1  (1.-gp2)*(dfddx(9)+dfddy(10)+dfddz(11))
     2    +gp2*gw*(
     3   -(f(13,i,j,k)*fs(1,0,1)-f(5,i,j,k)*fs(3,0,1))
     3   -(f(14,i,j,k)*fs(1,0,2)-f(6,i,j,k)*fs(3,0,2))
     3   -(f(15,i,j,k)*fs(1,0,3)-f(7,i,j,k)*fs(3,0,3)))
     4   -gp2*gw*(f(3,i,j,k)*cd(0,1)-f(1,i,j,k)*cd(0,3))
c
       s(18)=
     1  (1.-gp2)*(dfddx(13)+dfddy(14)+dfddz(15))
     2    +gp2*gw*(
     3   -(f(5,i,j,k)*fs(2,0,1)-f(9,i,j,k)*fs(1,0,1))
     3   -(f(6,i,j,k)*fs(2,0,2)-f(10,i,j,k)*fs(1,0,2))
     3   -(f(7,i,j,k)*fs(2,0,3)-f(11,i,j,k)*fs(1,0,3)))
     4   -gp2*gw*(f(1,i,j,k)*cd(0,2)-f(2,i,j,k)*cd(0,1))
c
      return
      end
