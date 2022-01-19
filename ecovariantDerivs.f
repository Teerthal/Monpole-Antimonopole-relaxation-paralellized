c 21 June 2016:
c These are covariant derivatives of the SO(3) (tHooft-Polyakov) model.
c I used Mathematica notebook "so3CovariantDerivs.nb" to find
c the expressions here. 
c
c 19 December 2015:
c
      subroutine ecovderiv(f,fd,i,j,k,dfdx,dfdy,dfdz,cd)
c
      implicit none
      include 'parameters.inc'
      integer i,j,k
      real*8 f,fd
      real*8 dfdx,dfdy,dfdz,cd
c
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension fd(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension cd(0:3,nscalar)
c
c covariant derivative: (D_\mu \Phi)^a = cd(mu,a).
c
        cd(0,1)=fd(1,i,j,k)
     1          +gw*f(3,i,j,k)*f(8,i,j,k)-gw*f(2,i,j,k)*f(12,i,j,k)
        cd(1,1)=dfdx(1)
     1          +gw*f(3,i,j,k)*f(9,i,j,k)-gw*f(2,i,j,k)*f(13,i,j,k)
        cd(2,1)=dfdy(1)
     1          +gw*f(3,i,j,k)*f(10,i,j,k)-gw*f(2,i,j,k)*f(14,i,j,k)
        cd(3,1)=dfdz(1)
     1          +gw*f(3,i,j,k)*f(11,i,j,k)-gw*f(2,i,j,k)*f(15,i,j,k)
        cd(0,2)=fd(2,i,j,k)
     1          -gw*f(3,i,j,k)*f(4,i,j,k)+gw*f(1,i,j,k)*f(12,i,j,k)
        cd(1,2)=dfdx(2)
     1          -gw*f(3,i,j,k)*f(5,i,j,k)+gw*f(1,i,j,k)*f(13,i,j,k)
        cd(2,2)=dfdy(2)
     1          -gw*f(3,i,j,k)*f(6,i,j,k)+gw*f(1,i,j,k)*f(14,i,j,k)
        cd(3,2)=dfdz(2)
     1          -gw*f(3,i,j,k)*f(7,i,j,k)+gw*f(1,i,j,k)*f(15,i,j,k)
        cd(0,3)=fd(3,i,j,k)
     1          +gw*f(2,i,j,k)*f(4,i,j,k)-gw*f(1,i,j,k)*f(8,i,j,k) 
        cd(1,3)=dfdx(3)
     1          +gw*f(2,i,j,k)*f(5,i,j,k)-gw*f(1,i,j,k)*f(9,i,j,k)
        cd(2,3)=dfdy(3)
     1          +gw*f(2,i,j,k)*f(6,i,j,k)-gw*f(1,i,j,k)*f(10,i,j,k)
        cd(3,3)=dfdz(3)
     1          +gw*f(2,i,j,k)*f(7,i,j,k)-gw*f(1,i,j,k)*f(11,i,j,k)
      return
      end
