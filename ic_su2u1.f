!=======================================================================
! Gauged SO(3) model. Initial conditions.
!
! 15 April 2015:
! The initial BPS monopole profile is taken from Vilenkin & Shellard, 
! Sec. 14.1 (1994):
! phi=vev*h(r)*{\hat x}
! A_i^a = (1-K(r)) \epsilon^{aij} x^j/(er^2) (note: overall sign in V&S is incorrect)
! h = 1./tanh(\xi)- 1/\xi
! K = \xi /sinh(\xi)
! \xi = m_v r = gw*vev*r
! For the non-BPS case, this is modified to:
! h = 1./tanh(\xi)- (1+m \xi)e^(-m \xi)/\xi
! K = \xi /sinh(\xi)
! where $m^2 = 2\lambda \eta^2$ ($\lambda$ occurs in the Higgs potential).
! Note: $m$ is really the scalar mass $m_s$, here written as "ms".
!=======================================================================

!----------Variable Definitions---------------------------
!--f = field array 
!--localsize = local array size
!--indices run from 1:localsize
!--h:profile functions
!--

      subroutine initialconditions(f,pid,
     1                            lxb,lyl,lzd,bb,fb,lb,rb,db,ub,np)
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'
      integer i,j,k,n
      integer itw,izm

!      real*8 zm (now declared in initialParameters.inc)
      real*8 xa,ya,za
      real*8 x,y,z,xgm,ygm,zgm,xga,yga,zga,rm,ra,rhom,rhoa
      real*8 sm,sa,ms,mv
      real*8 pprofilem,pprofilea,wprofilem,wprofilea
      real*8 yprofilem, yprofilea
      real*8 f
!      real*8 hatn
!      real*8 dxhatn,dyhatn,dzhatn
      real*8 phisq
!      real*8 dhdx,dhdy,dhdz,d2hdx2,d2hdy2,d2hdz2
      real*8 twist,ctw,stw,ctw2,stw2
      real*8 correctionm,correctiona
      real*8 theta_m,theta_a, sph_phi
      real*8 jac
      real*8 xua,yua,zua,xum,yum,zum,rum,rua
      integer bx,fx,ly,ry,dz,uz,lxb,lyl,lzd
      integer iglobal,jglobal,kglobal,np,pid,bb,fb,lb,rb,db,ub
    
      dimension jac(2,3,3)
      dimension f(nf,1:localsize,1:localsize,1:localsize)
!      dimension hatn(4,1:localsize,1:localsize,1:localsize)
!      dimension dxhatn(4,1:localsize,1:localsize,1:localsize)
!      dimension dyhatn(4,1:localsize,1:localsize,1:localsize)
!      dimension dzhatn(4,1:localsize,1:localsize,1:localsize)
!      dimension dhdx(nscalar),dhdy(nscalar),dhdz(nscalar)
!      dimension d2hdx2(nscalar),d2hdy2(nscalar),d2hdz2(nscalar)

!=======================================================================
! This section is for running the code on several nodes at the
! same time, each code with a different value of the parameters.
! read scan parameters and create files to write out successful runs 
!=======================================================================

        character*500 filename2
        character (len=200) winfile
        character (len=200) itwist
        character (len=200) izmono
      
        call getarg(1,winfile)
        read(winfile,'(a500)')filename2

        call getarg(2,itwist)
        read(itwist,*) itw
        twist=float(itw)*3.1416/6.

        call getarg(3,izmono)
        read(izmono,*) izm
        zm=(float(izm)+0.5)*dx
     
        if(pid==0) then
        print*,' ================= '
        print*,' lambda, twist, zm ', lambda,twist,zm
        print*,' ================= '
        endif

!=======================================================================
!       initialize:
!=======================================================================

        do k=1,localsize
            do j=1,localsize
                do i=1,localsize

                do n=1,nf
                  f(n,i,j,k)=0.
                enddo
                enddo
            enddo
        enddo



!=======================================================================
! fields for SU(1)xU(1):
!     f(1)=phi^1,f(2)=phi^2,f(3)=phi^3,f(4)=phi^4
!     f(5)=W^1_0, f(6)=W^1_1, f(7)=W^1_2, f(8)=W^1_3,
!     f(9)=W^2_0, f(10)=W^2_1, f(11)=W^2_2, f(12)=W^2_3,
!     f(13)=W^3_0, f(14)=W^3_1, f(15)=W^3_2, f(16)=W^3_3,
!     f(17)=B_0, f(18)=B_1, f(19)=B_2, f(20)=B_3 
!     gauge fixing functions (Num Rel): f(21)=\partial_i W^1_i, 
!     f(22)=\partial_i W^2_i, f(23)=\partial_i W^3_i, 
!     f(24)=\partial_iB_i
!     time derivatives:
!     fd(i)=time derivative of f(i) for i=1-4 (scalar fields),
!     fd(i)=electric fields(=W^a_{0j}) for i=5-20 (gauge fields).
!     fd(24)=time derivative of |phi|
!=======================================================================

!=======================================================================
! scalar mass in model:
!=======================================================================

        ms=sqrt(2.*lambda)*vev
        mv=gw*vev

        if(pid==0) then
        print*, ' ms, mv, & ms/mv: ', ms,mv,ms/mv
        print*,' ======================= '
        print*, ' '
        endif

!=======================================================================
! Initial location of monopole is specified in initialParameters.inc
! (or passed on as an argument when executing) and is chosen so that 
! it is not on a lattice site.
! Antimonopole will be at (-xm,-ym,-zm):
!=======================================================================

        xa=-xm
        ya=-ym
        za=-zm

! cosine and sine of twist:

        ctw=cos(twist)
        stw=sin(twist)
        ctw2=cos(twist/2.)
        stw2=sin(twist/2.)


!=======================================================================
! start: loop over larger lattice (+1) for scalar fields:
!=======================================================================

        do k=1,localsize
            do j=1,localsize
                do i=1,localsize

! ------------------ begin: setup variables -----------------
!               (x,y,z) for lattice site:
!--local labels are actually global coordinates
!--x,y,z are global coordinates
!--xgm, ... xga : global coordinates for monopole(m) and antimonopole(a) 
!---Boosts are in the x-direction
!--h are the higgs isovector direction as in eq(12) of 1705.0309

                iglobal = i+lxb-1
                jglobal = j+lyl-1
                kglobal = k+lzd-1

                x=float(iglobal)*dx
                y=float(jglobal)*dx
                z=float(kglobal)*dx

! Unboosted coordinates

                xum = (x-xm)
                yum = (y-ym)
                zum = (z-zm)
                xua = (x-xa)
                yua = (y-ya)
                zua = (z-za)
                rum = xum**2+yum**2+zum**2
                rua = xua**2+yua**2+zua**2

! Spherical coordinates in terms of non-boosted
! cartesian coordinates                

                theta_m = acos(zum/rum)
                theta_a = acos(zua/rua)
! Since y and x coordinates are the same for the monpole and
! antimonopole, we use monopole-centered coordinates here to
! define the azimuthal angle. 4.D0*DATAN(1.D0) = pi

                if(xgm.gt.0.) then
                sph_phi = atan(yum/xum)
                else if(xgm.eq.0.) then
                sph_phi = atan(yum/xum) + 4.D0*DATAN(1.D0)
                else
                sph_phi = 2.D0*DATAN(1.D0)
                endif

! dimensionless distances (in units of vector mass):
!Need to correct for su2u1!

                sm=gw*vev*rum
                sa=gw*vev*rua

! -- end: setup variables ------
!
!--begin scalar profile -----
! monopole and antimonopole scalar profile functions:

                if(sm.ne.0.) then
                pprofilem=vev*(1./tanh(sm)-(1.+ms*sm)*exp(-ms*sm)/sm)
                else
                pprofilem=0.
                endif

                if(sa.ne.0.) then
                pprofilea=vev*(1./tanh(sa)-(1.+ms*sa)*exp(-ms*sa)/sa)
                else
                pprofilea=0.
                endif

!---------------------end scalar profiles and derivatives --------------

! Profiles for the scalar

        if(rum*rua.ne.0.) then
    
        f(1,i,j,k)=pprofilem*pprofilea*(sin(theta_m/2)*sin(theta_a/2)*
    1   ctw+cos(theta_m/2)*cos(theta_a/2))/vev
    
        f(2,i,j,k)=pprofilem*pprofilea*(sin(theta_m/2)*sin(theta_a/2)*
    1   stw)/vev

        f(3,i,j,k)=pprofilem*pprofilea*(sin(theta_m/2)*cos(theta_a/2)*
    1   cos(sph_phi)-cos(theta_m/2)*sin(theta_a/2)*cos(sph_phi-twist))/vev

        f(4,i,j,k)=pprofilem*pprofilea*(sin(theta_m/2)*cos(theta_a/2)*
    1   sin(sph_phi)-cos(theta_m/2)*sin(theta_a/2)*sin(sph_phi-twist))/vev

        f(1,i,j,k)=0.
        f(2,i,j,k)=0.
        f(3,i,j,k)=0.
        f(4,i,j,k)=0.
        
        endif

                 enddo
             enddo
         enddo



! Removed the piece of code for computing dxhatn

!=======================================================================
! now for the initial gauge fields and the electric fields:
!=======================================================================

        do k=1,localsize
            do j=1,localsize
                do i=1,localsize
 
! -- begin: setup variables ------
! (x,y,z) for lattice site:
  
                iglobal = i+lxb-1
                jglobal = j+lyl-1
                kglobal = k+lzd-1

                x=float(iglobal)*dx
                y=float(jglobal)*dx
                z=float(kglobal)*dx

! Unboosted coordinates

                xum = (x-xm)
                yum = (y-ym)
                zum = (z-zm)
                xua = (x-xa)
                yua = (y-ya)
                zua = (z-za)
                rum = xum**2+yum**2+zum**2
                rua = xua**2+yua**2+zua**2

! Spherical coordinates in terms of non-boosted
! cartesian coordinates                

                theta_m = acos(zum/rum)
                theta_a = acos(zua/rua)
! Since y and x coordinates are the same for the monpole and
! antimonopole, we use monopole-centered coordinates here to
! define the azimuthal angle. 4.D0*DATAN(1.D0) = pi

                if(xgm.gt.0.) then
                sph_phi = atan(yum/xum)
                else if(xgm.eq.0.) then
                sph_phi = atan(yum/xum) + 4.D0*DATAN(1.D0)
                else
                sph_phi = 2.D0*DATAN(1.D0)
                endif

! dimensionless distances (in units of vector mass):

                sm=gw*vev*rum
                sa=gw*vev*rua

! Jacobian elements
! First index:1/2 -> monopole/antimonopole centered coordinates
! These are the elements corresponding to the partial derivatives of
! the spherical coordinates w.r.t the un-boosted cartesian coordinates.

                jac(1,1,1) = xum/rum
                jac(1,1,2) = yum/rum
                jac(1,1,3) = zum/rum
                jac(1,2,1) = xum*zum/(rum**2.*sqrt(xum**2+yum**2.))
                jac(1,2,2) = yum*zum/(rum**2.*sqrt(xum**2.+yum**2.))
                jac(1,2,3) = -sqrt(xum**2.+yum**2.)/(rum**2.)
                jac(1,3,1) = -yum/(xum**2.+yum**2.)
                jac(1,3,2) = xum/(xum**2.+yum**2.)
                jac(1,3,3) = 0.

                jac(2,1,1) = xua/rua
                jac(2,1,2) = yua/rua
                jac(2,1,3) = zua/rua
                jac(2,2,1) = xua*zua/(rua**2.*sqrt(xua**2+yua**2.))
                jac(2,2,2) = yua*zua/(rua**2.*sqrt(xua**2.+yua**2.))
                jac(2,2,3) = -sqrt(xua**2.+yua**2.)/(rua**2.)
                jac(2,3,1) = -yua/(xua**2.+yua**2.)
                jac(2,3,2) = xua/(xua**2.+yua**2.)
                jac(2,3,3) = 0.


!-----------------------------------------------------------------------
! monopole and antimonopole gauge profile function (1-K(r)):
! Note -- no factor of 1/r as in icMMbar.f and icMMbarTwisted.f
! Correction factors for profiles as evaluated by Ayush Saurabh:
!-----------------------------------------------------------------------

        correctionm=(1.+(1.-sqrt(lambda))*sm**2/4.+sm**4/16.)/
     1                (1.+sm**2/4.+sm**4/16.)
        correctiona=(1.+(1.-sqrt(lambda))*sa**2/4.+sa**4/16.)/
     1                (1.+sa**2/4.+sa**4/16.)
        if(sm.ne.0.) then
        wprofilem=(1.-(sm/sinh(sm))*correctionm)/gw
        else
        wprofilem=0.
        endif
        if(sa.ne.0.) then
        wprofilea=(1.-(sa/sinh(sa))*correctiona)/gw
        else
        wprofilea=0.
        endif

! square of scalar field su2u1:

       phisq=f(1,i,j,k)**2+f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2

       if(phisq.ne.0.) then
! W_\mu^1:
       f(5,i,j,k)=0.

       f(6,i,j,k)=(gw*wprofilem*wprofilea
     1  *(jac(1,2,1)*(cos(theta_a)*cos(twist-sph_phi)*sin(twist)-cos(twist)*
     2 sin(twist-sph_phi))+(1./4.)*jac(2,2,1)*(2.*(gw**2-1.)*cos(twist-
     3 sph_phi)*(cos(theta_a)*sin(2.*twist)*(sin(theta_m))**2-sin(twist)
     4 *sin(2.*theta_m)*sin(theta_a))+(sin(twist-sph_phi)*(3.+gw**2-
     5 (gw**2-1.)*cos(2.*twist)-2.*(gw**2-1)*cos(2.*theta_m)*stw**2)))
     6 +(-gw**2+(gw**2-1.)*cos(theta_m)*cos(theta_a)+(gw**2-1.)*ctw
     7 *sin(theta_m)*sin(theta_a))*(cos(twist-sph_phi)*(cos(theta_m)
     8 *sin(theta_a)-ctw*cos(theta_a)*sin(theta_m))-stw*sin(theta_m)
     9 *sin(twist-sph_phi))*jac(1,3,1)))

       f(7,i,j,k)=(gw*wprofilem*wprofilea
     1  *(jac(1,2,2)*(cos(theta_a)*cos(twist-sph_phi)*sin(twist)-cos(twist)*
     2 sin(twist-sph_phi))+(1./4.)*jac(2,2,2)*(2.*(gw**2-1.)*cos(twist-
     3 sph_phi)*(cos(theta_a)*sin(2.*twist)*(sin(theta_m))**2-sin(twist)
     4 *sin(2.*theta_m)*sin(theta_a))+(sin(twist-sph_phi)*(3.+gw**2-
     5 (gw**2-1.)*cos(2.*twist)-2.*(gw**2-1)*cos(2.*theta_m)*stw**2)))
     6 +(-gw**2+(gw**2-1.)*cos(theta_m)*cos(theta_a)+(gw**2-1.)*ctw
     7 *sin(theta_m)*sin(theta_a))*(cos(twist-sph_phi)*(cos(theta_m)
     8 *sin(theta_a)-ctw*cos(theta_a)*sin(theta_m))-stw*sin(theta_m)
     9 *sin(twist-sph_phi))*jac(1,3,2)))

       f(8,i,j,k)=(gw*wprofilem*wprofilea
     1  *(jac(1,2,3)*(cos(theta_a)*cos(twist-sph_p)*sin(twist)-cos(twist)*
     2 sin(twist-sph_phi))+(1./4.)*jac(2,2,3)*(2.*(gw**2-1.)*cos(twist-
     3 sph_phi)*(cos(theta_a)*sin(2.*twist)*(sin(theta_m))**2-sin(twist)
     4 *sin(2.*theta_m)*sin(theta_a))+(sin(twist-sph_phi)*(3.+gw**2-
     5 (gw**2-1.)*cos(2.*twist)-2.*(gw**2-1)*cos(2.*theta_m)*stw**2)))))

! W_\mu^2:
       f(9,i,j,k)=0.

       f(10,i,j,k)=(gw*wprofilem*wprofilea
     1  *((1./8.)*jac(2,2,1)*(2.*cos(twist-sph_phi)*(3.-gw**2-(gw**2
     2 -1.)*cos(2.*twist) - 2.*(gw**2-1.)*cos(2.*theta_m)*stw**2)-4.*
     3 (gw**2-1.)*sin(twist-sph_phi)*(cos(theta_m)*sin(2.*twist)*
     4 (sin(theta_m)**2)-stw*sin(2.*theta_m)*sin(theta_a)))-jac(1,2,1)
     5 *(ctw*cos(twist-sph_phi)+cos(theta_a)*sin(twist-sph_phi))
     6 -(-gw**2+(gw**2-1.)*cos(theta_m)*cos(theta_a)+(gw**2-1.)*ctw
     7 *sin(theta_m)*sin(theta_a))*(cos(twist-sph_phi)*(cos(theta_m)
     8 *sin(theta_a)-ctw*cos(theta_a)*sin(theta_m))-stw*sin(theta_m)
     9 *sin(twist-sph_phi))*jac(1,3,1)))
     
       f(11,i,j,k)=(gw*wprofilem*wprofilea
     1  *((1./8.)*jac(2,2,2)*(2.*cos(twist-sph_phi)*(3.-gw**2-(gw**2
     2 -1.)*cos(2.*twist) - 2.*(gw**2-1.)*cos(2.*theta_m)*stw**2)-4.*
     3 (gw**2-1.)*sin(twist-sph_phi)*(cos(theta_m)*sin(2.*twist)*
     4 (sin(theta_m)**2)-stw*sin(2.*theta_m)*sin(theta_a)))-jac(1,2,2)
     5 *(ctw*cos(twist-sph_phi)+cos(theta_a)*sin(twist-sph_phi))
     6 -(-gw**2+(gw**2-1.)*cos(theta_m)*cos(theta_a)+(gw**2-1.)*ctw
     7 *sin(theta_m)*sin(theta_a))*(cos(twist-sph_phi)*(cos(theta_m)
     8 *sin(theta_a)-ctw*cos(theta_a)*sin(theta_m))-stw*sin(theta_m)
     9 *sin(twist-sph_phi))*jac(1,3,2)))

       f(12,i,j,k)=(gw*wprofilem*wprofilea
     1  *((1./8.)*jac(2,2,3)*(2.*cos(twist-sph_phi)*(3.-gw**2-(gw**2
     2 -1.)*cos(2.*twist) - 2.*(gw**2-1.)*cos(2.*theta_m)*stw**2)-4.*
     3 (gw**2-1.)*sin(twist-sph_phi)*(cos(theta_m)*sin(2.*twist)*
     4 (sin(theta_m)**2)-stw*sin(2.*theta_m)*sin(theta_a)))-jac(1,2,3)
     5 *(ctw*cos(twist-sph_phi)+cos(theta_a)*sin(twist-sph_phi))))

! W_\mu^3:
       f(13,i,j,k)=0.
       f(14,i,j,k)=(gw*wprofilem*wprofilea
     1  *(stw*sin(theta_a)*jac(1,2,1) + (cos(theta_m)*cos(theta_a)+ctw
     2 *sin(theta_m)*sin(theta_a))*((gw**2-1.)*sin(theta_m)*jac(2,2,3)
     3 +(gw**2-(gw**2-1.)*cos(theta_m)*cos(theta_a)-(gw**2-1.)*ctw
     4 *sin(theta_m)*sin(theta_a))*jac(1,3,1))-jac(1,3,1)))

       f(15,i,j,k)=(gw*wprofilem*wprofilea
     1  *(stw*sin(theta_a)*jac(1,2,2) + (cos(theta_m)*cos(theta_a)+ctw
     2 *sin(theta_m)*sin(theta_a))*((gw**2-1.)*sin(theta_m)*jac(2,2,3)
     3 +(gw**2-(gw**2-1.)*cos(theta_m)*cos(theta_a)-(gw**2-1.)*ctw
     4 *sin(theta_m)*sin(theta_a))*jac(1,3,2))-jac(1,3,2)))

       f(16,i,j,k)=(gw*wprofilem*wprofilea
     1  *(stw*sin(theta_a)*jac(1,2,3) + (cos(theta_m)*cos(theta_a)+ctw
     2 *sin(theta_m)*sin(theta_a))*((gw**2-1.)*sin(theta_m)*jac(2,2,3)
     3 +(gw**2-(gw**2-1.)*cos(theta_m)*cos(theta_a)-(gw**2-1.)*ctw
     4 *sin(theta_m)*sin(theta_a))*jac(1,3,3))-jac(1,3,3)))

! Guessing profile from yiyand and T.
      pprofilem = 0.
      pprofilea = 0.

! B_0 = 0:
       f(17,i,j,k)=0.

! B_1:
       f(18,i,j,k)=(gy*yprofilem*yprofilea*gy**2*(stw*sin(theta_m)*
     1             jac(2,2,1)-(-1.+cos(theta_m)*cos(theta_a)+ctw*
     2             sin(theta_m)*sin(theta_a))*jac(1,3,1)))

! B_2
       f(19,i,j,k)=(gy*yprofilem*yprofilea*gy**2*(stw*sin(theta_m)*
     1             jac(2,2,2)-(-1.+cos(theta_m)*cos(theta_a)+ctw*
     2             sin(theta_m)*sin(theta_a))*jac(1,3,2)))

! B_3
       f(20,i,j,k)=(gy*yprofilem*yprofilea*gy**2*(stw*sin(theta_m)*
     1             jac(2,2,3)-(-1.+cos(theta_m)*cos(theta_a)+ctw*
     2             sin(theta_m)*sin(theta_a))*jac(1,3,3)))

       else
       f(5,i,j,k)=0.
       f(6,i,j,k)=0.
       f(7,i,j,k)=0.
       f(8,i,j,k)=0.
       f(9,i,j,k)=0.
       f(10,i,j,k)=0.
       f(11,i,j,k)=0.
       f(12,i,j,k)=0.
       f(13,i,j,k)=0.
       f(14,i,j,k)=0.
       f(15,i,j,k)=0.
       f(16,i,j,k)=0.
       f(17,i,j,k)=0.
       f(18,i,j,k)=0.
       f(19,i,j,k)=0.
       f(20,i,j,k)=0.

       endif

                enddo
            enddo
        enddo
        
!=======================================================================
! gauge functions $\Gamma^a = \partial_i W^a_i$,$\Xi_i = \partial_iB_i$:
! This uses a call to derivatives. That is why it is necessary to
! evaluate the gauge functions in a separate do loop (after all the
! gauge fields have been specified).
!=======================================================================

!-------------Set local indices for gamma -x----------------------------

                bx=4
                fx=localsize-3

!-------------Set local indices for gamma -y----------------------------

                ly=4
                ry=localsize-3

!-------------Set local indices for gamma -z----------------------------

                dz=4
                uz=localsize-3

        call icgamma(f,bx,fx,ly,ry,dz,uz,lxb,lyl,lzd)

        return
        end
