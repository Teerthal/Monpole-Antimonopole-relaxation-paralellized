!=======================================================================
! 18 January 2016 -- for SO(3)
!=======================================================================
      
      subroutine covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)

      implicit none
      include 'parameters.inc'
      integer i,j,k
      real*8 f,fd
      real*8 dfdx,dfdy,dfdz,cd

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension cd(0:3,3)

!=======================================================================
! covariant derivative: (D_\mu \Phi)^a = cd(mu,a) where
! a=1 means real part of upper component of doublet, a=2 is imaginary
! part of upper component, a=3 is real part of lower component,
! and a=4 is imaginary part of lower component..
!
!
! SO(3) case: 
! Dictionary: phi^1=f(1),phi^2=f(2),phi^3=f(3); W^1_0=f(4)...;
! W^2_0=f(8)...; W^3_0=f(12)....
!=======================================================================


       cd(1,1)=dfdx(1)
     1         +gw*(f(9,i,j,k)*f(3,i,j,k)-f(13,i,j,k)*f(2,i,j,k))
       cd(1,2)=dfdx(2)
     1         +gw*(f(13,i,j,k)*f(1,i,j,k)-f(5,i,j,k)*f(3,i,j,k))
       cd(1,3)=dfdx(3)
     1         +gw*(f(5,i,j,k)*f(2,i,j,k)-f(9,i,j,k)*f(1,i,j,k))

       cd(2,1)=dfdy(1)
     1         +gw*(f(10,i,j,k)*f(3,i,j,k)-f(14,i,j,k)*f(2,i,j,k))
       cd(2,2)=dfdy(2)
     1         +gw*(f(14,i,j,k)*f(1,i,j,k)-f(6,i,j,k)*f(3,i,j,k))
       cd(2,3)=dfdy(3)
     1         +gw*(f(6,i,j,k)*f(2,i,j,k)-f(10,i,j,k)*f(1,i,j,k))

       cd(3,1)=dfdz(1)
     1         +gw*(f(11,i,j,k)*f(3,i,j,k)-f(15,i,j,k)*f(2,i,j,k))
       cd(3,2)=dfdz(2)
     1         +gw*(f(15,i,j,k)*f(1,i,j,k)-f(7,i,j,k)*f(3,i,j,k))
       cd(3,3)=dfdz(3)
     1         +gw*(f(7,i,j,k)*f(2,i,j,k)-f(11,i,j,k)*f(1,i,j,k))

       cd(0,1)=-vel*cd(1,1)
       cd(0,2)=-vel*cd(1,2)
       cd(0,3)=-vel*cd(1,3)

      return
      end
