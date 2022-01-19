c 19 December 2015:
c This subroutine calculates the gauge field strengths for all fields
c at a given spatial point. At the moment, all the equations are
c explicitly listed. It would be better to organize in a do loop over
c fields. However, that will need some thinking, e.g. parameters
c that specify which fields are scalar and gauge fields (so that
c the do loop can run over only the gauge fields), e.g. instead
c writing x-derivative as dfdx(...) write it as df(1,...) so that
c a loop over derivatives can be performed..
c
      subroutine efieldstrength(f,fd,i,j,k,dfdx,dfdy,dfdz,fs)
c
      implicit none
      include 'parameters.inc'
      integer i,j,k
      real*8 f,fd,r
      real*8 dfdx,dfdy,dfdz,fs
c
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension fd(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension fs(ngauge/4,0:3,0:3)
c
c field strength: W^a_{\mu\nu}=fs(a,mu,nu):
c
c a=1:
      fs(1,0,0)=0.
      fs(1,0,1)=fd(5,i,j,k)
      fs(1,0,2)=fd(6,i,j,k)
      fs(1,0,3)=fd(7,i,j,k)
      fs(1,1,0)=-fs(1,0,1)
      fs(1,1,1)=0.
      fs(1,1,2)=dfdx(6)-dfdy(5)
     1         +gw*(f(9,i,j,k)*f(14,i,j,k)-f(10,i,j,k)*f(13,i,j,k))
      fs(1,1,3)=dfdx(7)-dfdz(5)
     1         +gw*(f(9,i,j,k)*f(15,i,j,k)-f(11,i,j,k)*f(13,i,j,k))
      fs(1,2,0)=-fs(1,0,2)
      fs(1,2,1)=-fs(1,1,2)
      fs(1,2,2)=0.
      fs(1,2,3)=dfdy(7)-dfdz(6)
     1         +gw*(f(10,i,j,k)*f(15,i,j,k)-f(11,i,j,k)*f(14,i,j,k))
      fs(1,3,0)=-fs(1,0,3)
      fs(1,3,1)=-fs(1,1,3)
      fs(1,3,2)=-fs(1,2,3)
      fs(1,3,3)=0.
c a=2:
      fs(2,0,0)=0.
      fs(2,0,1)=fd(9,i,j,k)
      fs(2,0,2)=fd(10,i,j,k)
      fs(2,0,3)=fd(11,i,j,k)
      fs(2,1,0)=-fs(2,0,1)
      fs(2,1,1)=0.
      fs(2,1,2)=dfdx(10)-dfdy(9)
     1         +gw*(f(13,i,j,k)*f(6,i,j,k)-f(14,i,j,k)*f(5,i,j,k))
      fs(2,1,3)=dfdx(11)-dfdz(9)
     1         +gw*(f(13,i,j,k)*f(7,i,j,k)-f(15,i,j,k)*f(5,i,j,k))
      fs(2,2,0)=-fs(2,0,2)
      fs(2,2,1)=-fs(2,1,2)
      fs(2,2,2)=0.
      fs(2,2,3)=dfdy(11)-dfdz(10)
     1         +gw*(f(14,i,j,k)*f(7,i,j,k)-f(15,i,j,k)*f(6,i,j,k))
      fs(2,3,0)=-fs(2,0,3)
      fs(2,3,1)=-fs(2,1,3)
      fs(2,3,2)=-fs(2,2,3)
      fs(2,3,3)=0.
c a=3:
      fs(3,0,0)=0.
      fs(3,0,1)=fd(13,i,j,k)
      fs(3,0,2)=fd(14,i,j,k)
      fs(3,0,3)=fd(15,i,j,k)
      fs(3,1,0)=-fs(3,0,1)
      fs(3,1,1)=0.
      fs(3,1,2)=dfdx(14)-dfdy(13)
     1         +gw*(f(5,i,j,k)*f(10,i,j,k)-f(6,i,j,k)*f(9,i,j,k))
      fs(3,1,3)=dfdx(15)-dfdz(13)
     1         +gw*(f(5,i,j,k)*f(11,i,j,k)-f(7,i,j,k)*f(9,i,j,k))
      fs(3,2,0)=-fs(3,0,2)
      fs(3,2,1)=-fs(3,1,2)
      fs(3,2,2)=0.
      fs(3,2,3)=dfdy(15)-dfdz(14)
     1         +gw*(f(6,i,j,k)*f(11,i,j,k)-f(7,i,j,k)*f(10,i,j,k))
      fs(3,3,0)=-fs(3,0,3)
      fs(3,3,1)=-fs(3,1,3)
      fs(3,3,2)=-fs(3,2,3)
      fs(3,3,3)=0.
c
      return
      end
