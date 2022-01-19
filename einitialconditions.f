c 8 July 2016: I first tried initial conditions with only f(3)=vev
c and f(1)=0=f(2) and the gauge waves only in W^3. This led to trivial 
c classical scattering because the gauge because the W^3 in this case
c is the photon and it can easily be checked that photon-photon collisions 
c do not excite any of the other fields in the classical equations of motion.
c I then tried f(1)=vev, f(2)=0=f(3), keeping the gauge waves still in
c W^3. Here I got non-trivial scattering but f(3), W^1 and W^2 remained
c zero throughout the evolution. I then checked that this evolution
c is consistent with the equations of motion and we are effectively
c evolving in a U(1) sub-space of the SO(3) model. However, I did find
c that the evolution led to zeros of the Higgs field if I chose large
c values of amp and omega (e.g. amp=2,omega=10). The zeros were 
c propagating along the z-direction but slightly off the z-axis.
c (The monopole winding around these zeros trivially vanished because
c f(3)=0.) Thinking back, the only way to reconcile these zeros of the 
c Higgs is that the collision produced a loop of string moving along
c the z-direction. This is important to revisit because the production
c of string loops is also relevant in the electroweak model.
c To study monopole production, I had to generalize to the case when
c both f(1) and f(3) are initially non-zero. So I chose f(1)=f(2)=vev/sqrt(2.).
c With omega=10, amp=3, I am finding zeros of the Higgs at late times
c that carry non-trivial monopole winding.
c
c 28 June 2016: Make sure that initial conditions satisfy Gauss
c law constraints.
c
c 26 June 2016: Now evolving using Numerical Relativity methods.
c One issue is that the gauge used is W^a_0=0. However, if we
c simply boost a static monopole solution we get $W^a_0 \ne 0$.
c So then we need to gauge transform back to $W^a_0=0$. I haven't
c yet figured out how to do this. This complication doesn't apply
c to the case when the mmbar start at rest as then no boost is
c needed. For now I am focussing on getting the code to work 
c when there is no boost.
c
c 28 September 2015: Fixed sign of W^a_0 and now the energy of
c boosted monopole-antimonopole agrees approximatedly with that 
c calculated from the BPS formula multiplied by the Lorentz boost factor.
c
c 16 September 2015: Trying to implement mmbar initial conditions
c in which monopole and antimonopole are boosted in opposite
c directions. More specifically, monopole and antimonopole are
c located along the z-axis and the monopole is boosted (in a general
c direction) and the antimonopole in the opposite direction, so
c that we are always in the center of momentum frame. As of now,
c the v=0 initial conditions are working but the v!=0 initial
c conditions are not.
c
c 10 August 2015: 
c Gauged SO(3) model. Initial conditions.
c
c 15 April 2015:

      subroutine einitialconditions(f,fd,fnew,fdnew,pid,
     1                            lxb,lyl,lzd,np,bb,fb,lb,rb,db,ub)
      implicit none
      include 'parameters.inc'
      include 'einitialParameters.inc'
      integer i,j,k,n
c xm,ym,zm,vxm,vym,vzm are now declared in initialParameters.inc
c      real*8 xm,ym,zm,vxm,vym,vzm
      real*8 plusminus
      real*8 xa,ya,za,vxa,vya,vza
      real*8 uvxm,uvym,uvzm,uvxa,uvya,uvza
      real*8 x,y,z,xgm,ygm,zgm,xga,yga,zga,rm,ra,rhom,rhoa
      real*8 sm,sa,ms
      real*8 uvdotxmrest,uvdotxarest
      real*8 pprofilem,pprofilea,wprofilem,wprofilea
      real*8 uvdotxm,uvdotxa,vdotxm,vdotxa
      real*8 uvcrossxm,uvcrossxa,vcrossxm,vcrossxa 
      real*8 pdashm,wdashm,pdasha,wdasha
      real*8 dtrhom, dtrhoa
      real*8 f,fd,fnew,fdnew
      real*8 dfdx,dfdy,dfdz,fs,cd,d2fdx2,d2fdy2,d2fdz2
      real*8 dfddx,dfddy,dfddz
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub
      integer ilocal,jlocal,klocal,np,pid
c
      dimension vcrossxm(3),uvcrossxm(3),vcrossxa(3),uvcrossxa(3) 
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension fd(nf,1:localsize,1:localsize,1:localsize)
      dimension fnew(nf,1:localsize,1:localsize,1:localsize)
      dimension fdnew(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
      dimension fs(3,0:3,0:3)
      dimension cd(0:3,3)
c =====================
c begin{dictionary}
c Gauged SO(3) model contains 3 real scalar fields:
c f(1)= phi(1), f(2)=phi(2), f(3)=phi(3), 
c and gauge fields: 
c f(4)=W^1_0, f(5)=W^1_1, f(6)=W^1_2, f(7)=W^1_3,
c f(8)=W^2_0, f(9)=W^2_1, f(10)=W^2_2, f(11)=W^2_3,
c f(12)=W^3_0, f(13)=W^3_1, f(14)=W^3_2, f(15)=W^3_3.
c Also: fd(n)=dot(f(n)). 
c end{dictionary}
c =====================
c scalar mass in model:


c start: loop over lattice:

        fd=0.0
        fnew=0.0
        fdnew=0.0

        do k=db,ub
            do j=lb,rb
                do i=bb,fb
c (x,y,z) for lattice site:
                ilocal = i+lxb-1
                jlocal = j+lyl-1
                klocal = k+lzd-1
                if(klocal.lt.0) plusminus = 1.0
                if(klocal.gt.0) plusminus = -1.0
                if(klocal.eq.0) plusminus = 0.0
!               plusminus=-1.0

                x=float(ilocal)*dx
                y=float(jlocal)*dx
                z=float(klocal)*dx

      call ederivatives(f,fd,i,j,k,dfdx,dfdy,dfdz,dfddx,dfddy,dfddz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

      call covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)

      call fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)

c give initial fields first:
c scalar field is in vacuum, pointing in the 1 & 3 direction.
c Note: if scalar points in pure 1 or 2 or 3 direction, dynamics
c is in a subspace and monopoles cannot be produced. To see this
c for phi in 3 direction is easy. To see this when phi is in
c 1 direction (for example), is a bit harder but one can check
c from the equations of motion that then phi^3=0=W^1_i=W^2_i
c  at all times. So need to take 2 components of phi to be
c non-zero.

c {\dot \phi}^a:
       fd(1,i,j,k)=plusminus*vel*cd(1,1)
       fd(2,i,j,k)=plusminus*vel*cd(1,2)
       fd(3,i,j,k)=plusminus*vel*cd(1,3)
c {\dot W}^1_\mu:
       fd(4,i,j,k)=0.0
       fd(5,i,j,k)=plusminus*vel*fs(1,1,1)
       fd(6,i,j,k)=plusminus*vel*fs(1,1,2)
       fd(7,i,j,k)=plusminus*vel*fs(1,1,3)
c {\dot W}^2_\mu:
       fd(8,i,j,k)=0.0
       fd(9,i,j,k)=plusminus*vel*fs(2,1,1)
       fd(10,i,j,k)=plusminus*vel*fs(2,1,2)
       fd(11,i,j,k)=plusminus*vel*fs(2,1,3)
c {\dot W}^3_\mu:
       fd(12,i,j,k)=0.0
       fd(13,i,j,k)=plusminus*vel*fs(3,1,1)
       fd(14,i,j,k)=plusminus*vel*fs(3,1,2)
       fd(15,i,j,k)=plusminus*vel*fs(3,1,3)

       fd(16,i,j,k)=0.0
       fd(17,i,j,k)=plusminus*vel*dfdx(17)
       fd(18,i,j,k)=plusminus*vel*dfdx(18)
       fd(19,i,j,k)=plusminus*vel*dfdx(19)

!         if(ilocal.eq.-latx.and.jlocal.eq.0.and.klocal.eq.0) then
!!!       print *, "latx"
!!        print *, "cd   ",cd(1,1),cd(2,1),cd(3,1)
!!        print *, "cd   ",cd(1,2),cd(2,2),cd(3,2)
!!        print *, "cd   ",cd(1,3),cd(2,3),cd(3,3)
!!        print *, f(1,i,j,k),f(2,i,j,k),f(3,i,j,k) 
!         print*,ilocal,jlocal,klocal,pid
!        print *, "cd   ",d(1,1),fd(1,i,j,k)
!!!       print*, dfdx(1),dfdx(2),dfdx(3)
!!!       print*,f(1,i,j,k)
!!!       print*,f(1,i+1,j,k),f(1,i+2,j,k),f(1,i+3,j,k)
!!!       print*,f(1,i-1,j,k),f(1,i-2,j,k),f(1,i-3,j,k)
!!!       print*,cd(1,1),dfdx(1),f(9,i,j,k),f(3,i,j,k),f(13,i,j,k),
!!!    1                      f(2,i,j,k),dfdx(2),dfdx(3)
!!!       print *, gradientEnergyPhi
!         endif


c initialize fnew and fdnew arrays:
!              do n=1,nf
!                 fnew(n,i,j,k)=0.
!                 fdnew(n,i,j,k)=0.
!              enddo

                  enddo
              enddo
          enddo

      return
      end
