c 18 January 2016 -- for SO(3)
c
      subroutine covd(f,i,j,k,hatn,dxhatn,dyhatn,dzhatn,cd)
c
      implicit none
      include 'parameters.inc'
      integer i,j,k
      real*8 f,cd
      real*8 hatn,dxhatn,dyhatn,dzhatn


      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
      dimension cd(0:3,3)
      dimension hatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dxhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dyhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dzhatn(3,-latx:latx,-laty:laty,-latz:latz)
      
c
c covariant derivative: (D_\mu \Phi)^a = cd(mu,a) where
c a=1 means real part of upper component of doublet, a=2 is imaginary
c part of upper component, a=3 is real part of lower component,
c and a=4 is imaginary part of lower component..
c
c
c SO(3) case: 
c Dictionary: phi^1=f(1),phi^2=f(2),phi^3=f(3); W^1_0=f(4)...;
c W^2_0=f(8)...; W^3_0=f(12)....
c
       cd(0,1)=0.
       cd(0,2)=0.
       cd(0,3)=0.
c
       cd(1,1)=dxhatn(1,i,j,k)
     1         +gw*(f(9,i,j,k)*hatn(3,i,j,k)-f(13,i,j,k)*hatn(2,i,j,k))
       cd(1,2)=dxhatn(2,i,j,k)
     1         +gw*(f(13,i,j,k)*hatn(1,i,j,k)-f(5,i,j,k)*hatn(3,i,j,k))
       cd(1,3)=dxhatn(3,i,j,k)
     1         +gw*(f(5,i,j,k)*hatn(2,i,j,k)-f(9,i,j,k)*hatn(1,i,j,k))
c
       cd(2,1)=dyhatn(1,i,j,k)
     1         +gw*(f(10,i,j,k)*hatn(3,i,j,k)-f(14,i,j,k)*hatn(2,i,j,k))
       cd(2,2)=dyhatn(2,i,j,k)
     1         +gw*(f(14,i,j,k)*hatn(1,i,j,k)-f(6,i,j,k)*hatn(3,i,j,k))
       cd(2,3)=dyhatn(3,i,j,k)
     1         +gw*(f(6,i,j,k)*hatn(2,i,j,k)-f(10,i,j,k)*hatn(1,i,j,k))
c
       cd(3,1)=dzhatn(1,i,j,k)
     1         +gw*(f(11,i,j,k)*hatn(3,i,j,k)-f(15,i,j,k)*hatn(2,i,j,k))
       cd(3,2)=dzhatn(2,i,j,k)
     1         +gw*(f(15,i,j,k)*hatn(1,i,j,k)-f(7,i,j,k)*hatn(3,i,j,k))
       cd(3,3)=dzhatn(3,i,j,k)
     1         +gw*(f(7,i,j,k)*hatn(2,i,j,k)-f(11,i,j,k)*hatn(1,i,j,k))
c
      return
      end
