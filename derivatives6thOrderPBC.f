! 2 December 2018: derivatives using periodic boundary conditions.
! Modified version of su2u1 6th order PBC subroutine for 
      subroutine derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1            d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)
      implicit none
      include 'parameters.inc'
      integer i,j,k,n,is
      integer ip3,ip2,ip1,im1,im2,im3
      integer jp3,jp2,jp1,jm1,jm2,jm3
      integer kp3,kp2,kp1,km1,km2,km3
      real*8 f,fd
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      integer ilocal,jlocal,klocal,lxb,lyl,lzd
!
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
!
!     Global coordinates
            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1
!
      do 35 n=1,nf
!
      ip3=i+3
      ip2=i+2
      ip1=i+1
      im1=i-1
      im2=i-2
      im3=i-3
! Coordinate shift from end of lattice to beginning for PBC
      if(ip3.gt.latx) ip3=i+3-2*latx
      if(ip2.gt.latx) ip2=i+2-2*latx
      if(ip1.gt.latx) ip1=i+1-2*latx
      if(im1.lt.-latx) im1=i-1+2*latx
      if(im2.lt.-latx) im2=i-2+2*latx
      if(im3.lt.-latx) im3=i-3+2*latx
c x-derivatives:
c using periodic boundary conditions:
        dfdx(n)=(f(n,ip3,j,k)-9.*f(n,ip2,j,k)+45.*f(n,ip1,j,k)
     1   -45.*f(n,im1,j,k)+9.*f(n,im2,j,k)-f(n,im3,j,k))/(60.*dx)
        d2fdx2(n)=(f(n,ip3,j,k)/90.-3.*f(n,ip2,j,k)/20.
     1   +3.*f(n,ip1,j,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,im1,j,k)/2.
     2   -3.*f(n,im2,j,k)/20.+f(n,im3,j,k)/90.)/(dx**2)
c
      jp3=j+3
      jp2=j+2
      jp1=j+1
      jm1=j-1
      jm2=j-2
      jm3=j-3
      if(jp3.gt.laty) jp3=j+3-2*laty
      if(jp2.gt.laty) jp2=j+2-2*laty
      if(jp1.gt.laty) jp1=j+1-2*laty
      if(jm1.lt.-laty) jm1=j-1+2*laty
      if(jm2.lt.-laty) jm2=j-2+2*laty
      if(jm3.lt.-laty) jm3=j-3+2*laty
c y-derivatives:
c using periodic boundary conditions:
        dfdy(n)=(f(n,i,jp3,k)-9.*f(n,i,jp2,k)+45.*f(n,i,jp1,k)
     1   -45.*f(n,i,jm1,k)+9.*f(n,i,jm2,k)-f(n,i,jm3,k))/(60.*dx)
        d2fdy2(n)=(f(n,i,jp3,k)/90.-3.*f(n,i,jp2,k)/20.
     1   +3.*f(n,i,jp1,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,jm1,k)/2.
     2   -3.*f(n,i,jm2,k)/20.+f(n,i,jm3,k)/90.)/(dx**2)
c
      kp3=k+3
      kp2=k+2
      kp1=k+1
      km1=k-1
      km2=k-2
      km3=k-3
      if(kp3.gt.latz) kp3=k+3-2*latz
      if(kp2.gt.latz) kp2=k+2-2*latz
      if(kp1.gt.latz) kp1=k+1-2*latz
      if(km1.lt.-latz) km1=k-1+2*latz
      if(km2.lt.-latz) km2=k-2+2*latz
      if(km3.lt.-latz) km3=k-3+2*latz
c z-derivatives:
c using periodic boundary conditions:
        dfdz(n)=(f(n,i,j,kp3)-9.*f(n,i,j,kp2)+45.*f(n,i,j,kp1)
     1   -45.*f(n,i,j,km1)+9.*f(n,i,j,km2)-f(n,i,j,km3))/(60.*dx)
        d2fdz2(n)=(f(n,i,j,kp3)/90.-3.*f(n,i,j,kp2)/20.
     1   +3.*f(n,i,j,kp1)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,j,km1)/2.
     2   -3.*f(n,i,j,km2)/20.+f(n,i,j,km3)/90.)/(dx**2)
c
35    continue
      return
      end
