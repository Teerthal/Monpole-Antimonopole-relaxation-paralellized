c here we define our differencing scheme, in this case 4th order.
c
      subroutine derivativesPBC(f,i,j,k,dfdx,dfdy,dfdz,
     1            d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)
      implicit none
      include 'parameters.inc'
      integer i,j,k,n,is,three1,two1,one1,three2,two2,one2
      integer ilocal,jlocal,klocal,lxb,lyl,lzd
      real*8 f
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
c
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
c
            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

      do n=1,nf
c
c x-derivatives:
!      if(abs(ilocal).le.latx-3) then
         three1=3
         three2=3
         two1=2
         two2=2
         one1=1
         one2=1
!      endif

!      if(ilocal.eq.latx-2) then
!        three1=-2*latx+3
!        three2=3
!        two1=2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(ilocal.eq.latx-1) then
!        three1=-2*latx+3
!        three2=3
!        two1=-2*latx+2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(ilocal.eq.latx) then
!        three1=-2*latx+3
!        three2=3
!        two1=-2*latx+2
!        two2=2
!        one1=-2*latx+1
!        one2=1
!      endif

!      if(ilocal.eq.-latx+2) then
!        three1=3
!        three2=-2*latx+3
!        two1=2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(ilocal.eq.-latx+1) then
!        three1=3
!        three2=-2*latx+3
!        two1=2
!        two2=-2*latx+2
!        one1=1
!        one2=1
!      endif

!      if(ilocal.eq.-latx) then
!        three1=3
!        three2=-2*latx+3
!        two1=2
!        two2=-2*latx+2
!        one1=1
!        one2=-2*latx+1
!      endif
c
c 6th order 1st- and 2nd- derivatives:
c central derivatives in lattice interior:
        dfdx(n)=(f(n,i+three1,j,k)-9.*f(n,i+two1,j,k)+45.
     1   *f(n,i+one1,j,k)-45.*f(n,i-one2,j,k)+9.*f(n,i-two2,j,k)
     2   -f(n,i-three2,j,k))/(60.*dx)
        d2fdx2(n)=(f(n,i+three1,j,k)/90.-3.*f(n,i+two1,j,k)/20.
     1   +3.*f(n,i+one1,j,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i-one2,j,k)/2.
     2   -3.*f(n,i-two2,j,k)/20.+f(n,i-three2,j,k)/90.)/(dx**2)
c
c y-derivatives:
!      if(abs(jlocal).le.laty-3) then
!        three1=3
!        three2=3
!        two1=2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(jlocal.eq.laty-2) then
!        three1=-2*laty+3
!        three2=3
!        two1=2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(jlocal.eq.laty-1) then
!        three1=-2*laty+3
!        three2=3
!        two1=-2*laty+2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(jlocal.eq.laty) then
!        three1=-2*laty+3
!        three2=3
!        two1=-2*laty+2
!        two2=2
!        one1=-2*laty+1
!        one2=1
!      endif

!      if(jlocal.eq.-laty+2) then
!        three1=3
!        three2=-2*laty+3
!        two1=2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(jlocal.eq.-laty+1) then
!        three1=3
!        three2=-2*laty+3
!        two1=2
!        two2=-2*laty+2
!        one1=1
!        one2=1
!      endif

!      if(jlocal.eq.-laty) then
!        three1=3
!        three2=-2*laty+3
!        two1=2
!        two2=-2*laty+2
!        one1=1
!        one2=-2*laty+1
!      endif
c
c 6th order 1st- and 2nd- derivatives:
c central derivatives in lattice interior:
        dfdy(n)=(f(n,i,j+three1,k)-9.*f(n,i,j+two1,k)+45.
     1   *f(n,i,j+one1,k)-45.*f(n,i,j-one2,k)+9.*f(n,i,j-two2,k)
     2   -f(n,i,j-three2,k))/(60.*dx)
        d2fdy2(n)=(f(n,i,j+three1,k)/90.-3.*f(n,i,j+two1,k)/20.
     1   +3.*f(n,i,j+one1,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,j-one2,k)/2.
     2   -3.*f(n,i,j-two2,k)/20.+f(n,i,j-three2,k)/90.)/(dx**2)
c
c
c z-derivatives:
      if(latz.le.3) then
        dfdz(n)=0
        d2fdz2(n)=0
      endif

      if(latz.gt.3) then
!      if(abs(klocal).le.latz-3) then
         three1=3
         three2=3
         two1=2
         two2=2
         one1=1
         one2=1
!      endif

!      if(klocal.eq.latz-2) then
!        three1=-2*latz+3
!        three2=3
!        two1=2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(klocal.eq.latz-1) then
!        three1=-2*latz+3
!        three2=3
!        two1=-2*latz+2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(klocal.eq.latz) then
!        three1=-2*latz+3
!        three2=3
!        two1=-2*latz+2
!        two2=2
!        one1=-2*latz+1
!        one2=1
!      endif

!      if(klocal.eq.-latz+2) then
!        three1=3
!        three2=-2*latz+3
!        two1=2
!        two2=2
!        one1=1
!        one2=1
!      endif

!      if(klocal.eq.-latz+1) then
!        three1=3
!        three2=-2*latz+3
!        two1=2
!        two2=-2*latz+2
!        one1=1
!        one2=1
!      endif

!      if(klocal.eq.-latz) then
!        three1=3
!        three2=-2*latz+3
!        two1=2
!        two2=-2*latz+2
!        one1=1
!        one2=-2*latz+1
!      endif

c 6th order 1st- and 2nd- derivatives:
c central derivatives in lattice interior:

        dfdz(n)=(f(n,i,j,k+three1)-9.*f(n,i,j,k+two1)+45.
     1   *f(n,i,j,k+one1)-45.*f(n,i,j,k-one2)+9.*f(n,i,j,k-two2)
     2   -f(n,i,j,k-three2))/(60.*dx)
        d2fdz2(n)=(f(n,i,j,k+three1)/90.-3.*f(n,i,j,k+two1)/20.
     1   +3.*f(n,i,j,k+one1)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,j,k-one2)/2.
     2   -3.*f(n,i,j,k-two2)/20.+f(n,i,j,k-three2)/90.)/(dx**2)
      endif
 
        enddo
      return
      end
