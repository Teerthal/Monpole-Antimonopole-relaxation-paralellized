       subroutine coefficient(coeff)
       implicit none
       include 'parameters.inc'
       integer i,j,k
       real*8 coeffx,coeffy,coeffz,coeff
       dimension coeff(-latx:latx,-laty:laty,-latz:latz)
c
       do 10 i=-latx,latx
       do 20 j=-laty,laty
       do 30 k=-latz,latz
c
c 2nd derivatives (4th order):       
       if(abs(i).le.latx-2) then
c coeff is the coefficient of the f(n,i,j,k) term in the formula
c for the second derivative.
        coeffx=-30.
c        d2fdx2(n)=(-f(n,i+2,j,k)+16.*f(n,i+1,j,k)-30.*f(n,i,j,k)
c     1         +16.*f(n,i-1,j,k)-f(n,i-2,j,k))/(12.*dx**2)
       endif
c
c one sided derivatives at 4th order near the boundary:
       if(i.eq.-(latx-1)) then
        coeffx=-20.
c        d2fdx2(n)=(-f(n,i+3,j,k)+4.*f(n,i+2,j,k)+6.*f(n,i+1,j,k)
c     1       -20.*f(n,i,j,k)+11.*f(n,i-1,j,k))/(12.*dx**2)
       endif
c
       if(i.eq.+(latx-1)) then
        coeffx=-20.
c        d2fdx2(n)=(11.*f(n,i+1,j,k)-20.*f(n,i,j,k)+6.*f(n,i-1,j,k)
c     1       +4.*f(n,i-2,j,k)-f(n,i-3,j,k))/(12.*dx**2)
       endif
c
       if(i.eq.latx) then
        coeffx=45.
c       d2fdx2(n)=(45.*f(n,i,j,k)-154.*f(n,i-1,j,k)
c     1         +214.*f(n,i-2,j,k)-156.*f(n,i-3,j,k)
c     2         +61.*f(n,i-4,j,k)-10.*f(n,i-5,j,k))/(12.*dx**2)
       endif
c
       if(i.eq.-latx) then
        coeffx=45.
c       d2fdx2(n)=(45.*f(n,i,j,k)-154.*f(n,i+1,j,k)
c     1         +214.*f(n,i+2,j,k)-156.*f(n,i+3,j,k)
c     2         +61.*f(n,i+4,j,k)-10.*f(n,i+5,j,k))/(12.*dx**2)
       endif
c
c y-derivatives:
c
c 2nd derivatives (4th order):       
       if(abs(j).le.laty-2) then 
        coeffy=-30.
c        d2fdy2(n)=(-f(n,i,j+2,k)+16.*f(n,i,j+1,k)-30.*f(n,i,j,k)
c     1         +16.*f(n,i,j-1,k)-f(n,i,j-2,k))/(12.*dx**2)
       endif
c       
       if(j.eq.-(laty-1)) then
        coeffy=-20.
c        d2fdy2(n)=(-f(n,i,j+3,k)+4.*f(n,i,j+2,k)+6.*f(n,i,j+1,k)
c     1       -20.*f(n,i,j,k)+11.*f(n,i,j-1,k))/(12.*dx**2)
       endif
c
       if(j.eq.+(laty-1)) then
        coeffy=-20.
c        d2fdy2(n)=(11.*f(n,i,j+1,k)-20.*f(n,i,j,k)+6.*f(n,i,j-1,k)
c     1       +4.*f(n,i,j-2,k)-f(n,i,j-3,k))/(12.*dx**2)
       endif
c
       if(j.eq.laty) then
        coeffy=45.
c        d2fdy2(n)=(45.*f(n,i,j,k)-154.*f(n,i,j-1,k)
c     1         +214.*f(n,i,j-2,k)-156.*f(n,i,j-3,k)
c     2         +61.*f(n,i,j-4,k)-10.*f(n,i,j-5,k))/(12.*dx**2)
       endif
c
       if(j.eq.-laty) then
        coeffy=45.
c        d2fdy2(n)=(45.*f(n,i,j,k)-154.*f(n,i,j+1,k)
c     1         +214.*f(n,i,j+2,k)-156.*f(n,i,j+3,k)
c     2         +61.*f(n,i,j+4,k)-10.*f(n,i,j+5,k))/(12.*dx**2)
       endif
c
c z-derivatives:
c
       if(abs(k).le.latz-2) then 
        coeffz=-30.
c      d2fdz2(n)=(-f(n,i,j,k+2)+16.*f(n,i,j,k+1)-30.*f(n,i,j,k)
c     1         +16.*f(n,i,j,k-1)-f(n,i,j,k-2))/(12.*dx**2)
       endif
c
       if(k.eq.-(latz-1)) then
        coeffz=-20.
c        d2fdz2(n)=(-f(n,i,j,k+3)+4.*f(n,i,j,k+2)+6.*f(n,i,j,k+1)
c     1       -20.*f(n,i,j,k)+11.*f(n,i,j,k-1))/(12.*dx**2)
       endif
c
       if(k.eq.+(latz-1)) then
        coeffz=-20.
c        d2fdz2(n)=(11.*f(n,i,j,k+1)-20.*f(n,i,j,k)+6.*f(n,i,j,k-1)
c     1       +4.*f(n,i,j,k-2)-f(n,i,j,k-3))/(12.*dx**2)
       endif
c
       if(k.eq.latz) then
        coeffz=45.
c       d2fdz2(n)=(45.*f(n,i,j,k)-154.*f(n,i,j,k-1)
c     1         +214.*f(n,i,j,k-2)-156.*f(n,i,j,k-3)
c     2         +61.*f(n,i,j,k-4)-10.*f(n,i,j,k-5))/(12.*dx**2)
       endif
c
       if(k.eq.-latz) then
        coeffz=45.
c        d2fdz2(n)=(45.*f(n,i,j,k)-154.*f(n,i,j,k+1)
c     1         +214.*f(n,i,j,k+2)-156.*f(n,i,j,k+3)
c     2         +61.*f(n,i,j,k+4)-10.*f(n,i,j,k+5))/(12.*dx**2)
       endif
c
35    continue
c
      coeff(i,j,k)=-(coeffx+coeffy+coeffz)
c
30    continue
20    continue
10    continue
c
      return
      end
