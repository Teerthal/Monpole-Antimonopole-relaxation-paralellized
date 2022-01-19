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

      subroutine initialconditions(f,hatn,dxhatn,dyhatn,dzhatn,pid,
     1                            lxb,lyl,lzd,np)
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
      real*8 f
      real*8 hatn
      real*8 dxhatn,dyhatn,dzhatn
      real*8 phisq
      real*8 dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2
      real*8 twist,ctw,stw,ctw2,stw2
      real*8 correctionm,correctiona
      integer bx,fx,ly,ry,dz,uz,lxb,lyl,lzd
      integer ilocal,jlocal,klocal,np,pid
 

      dimension f(nf,1:localsize,1:localsize,1:localsize)

      dimension hatn(3,1:localsize,1:localsize,1:localsize)
      dimension dxhatn(3,1:localsize,1:localsize,1:localsize)
      dimension dyhatn(3,1:localsize,1:localsize,1:localsize)
      dimension dzhatn(3,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)

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
                do n=1,3
                hatn(n,i,j,k)=0.
                dxhatn(n,i,j,k)=0.
                dyhatn(n,i,j,k)=0.
                dzhatn(n,i,j,k)=0.
                enddo
                enddo
            enddo
        enddo



!=======================================================================
! begin{dictionary}
! Gauged SO(3) model contains 3 real scalar fields:
! f(1)= phi(1), f(2)=phi(2), f(3)=phi(3), 
! and gauge fields: 
! f(4)=W^1_0, f(5)=W^1_1, f(6)=W^1_2, f(7)=W^1_3,
! f(8)=W^2_0, f(9)=W^2_1, f(10)=W^2_2, f(11)=W^2_3,
! f(12)=W^3_0, f(13)=W^3_1, f(14)=W^3_2, f(15)=W^3_3.
! Also: fd(n)=dot(f(n)). 
! end{dictionary}
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

                ilocal = i+lxb-1
                jlocal = j+lyl-1
                klocal = k+lzd-1

                x=float(ilocal)*dx
                y=float(jlocal)*dx
                z=float(klocal)*dx

                xgm=boost*(x-xm)
                ygm=y-ym
                zgm=(z-zm)

                xga=boost*(x-xa)
                yga=(y-ya)
                zga=(z-za)


! radial coordinate (with boosted x) with origin on monopole:

                rm=sqrt(xgm**2+ygm**2+zgm**2)
! radial coordinate (with boosted x) with origin on antimonopole:

                ra=sqrt(xga**2+yga**2+zga**2)

! xy-radial boosted coordinate with origin on monopole:
       
                rhom=sqrt(rm**2-zgm**2)

! xy-radial boosted coordinate with origin on antimonopole:

                rhoa=sqrt(ra**2-zga**2)

! dimensionless distances (in units of vector mass):

                sm=gw*vev*rm
                sa=gw*vev*ra

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
!
! \phi^a:
! (note: monopole-antimonopole separation = 2.*z1)
! the direction of the Higgs at a point whose angular coordinates
! are (thetam,phi) and (thetaa,phia=phi) w.r.t monopole and antimonopole
! is: (sin(thetam-thetaa)*(cos(phi),sin(phi)),cos(thetam-thetaa)).
! Some algebra shows this to be:(xg*2*z1,yg*2*z1,zg*zga+rho*rhoa)/(r*ra).
! This implies that the antimonpole \phi^1 and \phi^2 are the same
! as those for the monopole. However, \phi^3 for the antimonopole
! is -\phi^3 for the monopole, i.e. only the third component is
! reversed. Caution: this is not a simple "product ansatz".
!-----------------------------------------------------------------------

        if(rm*ra.ne.0.) then
        hatn(1,i,j,k)=(zga*(x*ctw2+y*stw2)*ctw-zgm*(x*ctw2+y*stw2)
     1     -ra*(y*ctw2-x*stw2)*stw)/(rm*ra)
        f(1,i,j,k)=pprofilem*pprofilea*hatn(1,i,j,k)/vev

        hatn(2,i,j,k)=(zga*(y*ctw2-x*stw2)*ctw-zgm*(y*ctw2-x*stw2)
     1     +ra*(x*ctw2+y*stw2)*stw)/(rm*ra)
        f(2,i,j,k)=pprofilem*pprofilea*hatn(2,i,j,k)/vev

        hatn(3,i,j,k)=(zgm*zga+ctw*(x**2+y**2))/(rm*ra)
        f(3,i,j,k)=pprofilem*pprofilea*hatn(3,i,j,k)/vev

        f(19,i,j,k)=pprofilem*pprofilea/vev

        else
        hatn(1,i,j,k)=0.
        hatn(2,i,j,k)=0.
        hatn(3,i,j,k)=1.
        
        f(1,i,j,k)=0.
        f(2,i,j,k)=0.
        f(3,i,j,k)=0.
        f(19,i,j,k)=0.
        endif
        
!-----------------------------------------------------------------------
! Following formulae assume there monopoles are not boosted (initial velocity=0).
! Also assumes monopole and antimonopole are initially on z-axis (xm=0=xa=ym=ya).
!-----------------------------------------------------------------------
        if(rm*ra.ne.0.) then
         dxhatn(1,i,j,k)=(ctw2*(zga*ctw-zgm)+stw2*ra*stw
     1                 -(y*ctw2-x*stw2)*(x/ra)*stw)/(rm*ra)
     2                              -(x/rm**2+x/ra**2)*hatn(1,i,j,k)
         dyhatn(1,i,j,k)=(stw2*(zga*ctw-zgm)-ctw2*ra*stw
     1                 -(y*ctw2-x*stw2)*(y/ra)*stw)/(rm*ra)
     2                              -(y/rm**2+y/ra**2)*hatn(1,i,j,k)
         dzhatn(1,i,j,k)=((x*ctw2+y*stw2)*(ctw-1.)
     1                 -(y*ctw2-x*stw2)*(zga/ra)*stw)/(rm*ra)
     2                          -(zgm/rm**2+zga/ra**2)*hatn(1,i,j,k)
         dxhatn(2,i,j,k)=(-stw2*(zga*ctw-zgm)+ctw2*ra*stw
     1                 +(x*ctw2+y*stw2)*(x/ra)*stw)/(rm*ra)
     2                              -(x/rm**2+x/ra**2)*hatn(2,i,j,k)
         dyhatn(2,i,j,k)=(ctw2*(zga*ctw-zgm)+stw2*ra*stw
     1                 +(x*ctw2+y*stw2)*(y/ra)*stw)/(rm*ra)
     2                              -(y/rm**2+y/ra**2)*hatn(2,i,j,k)
         dzhatn(2,i,j,k)=((y*ctw2-x*stw2)*(ctw-1.)
     1                  +(x*ctw2+y*stw2)*(zga/ra)*stw)/(rm*ra)
     2                              -(zgm/rm**2+zga/ra**2)*hatn(2,i,j,k)
         dxhatn(3,i,j,k)=2.*x*ctw/(rm*ra)
     1                    -(x/rm**2+x/ra**2)*hatn(3,i,j,k)
         dyhatn(3,i,j,k)=2.*y*ctw/(rm*ra)
     1                    -(y/rm**2+y/ra**2)*hatn(3,i,j,k)
         dzhatn(3,i,j,k)=(zgm+zga)/(rm*ra)
     1                    -(zgm/rm**2+zga/ra**2)*hatn(3,i,j,k)
      else
         dxhatn(1,i,j,k)=0.
         dyhatn(1,i,j,k)=0.
         dzhatn(1,i,j,k)=0.
         dxhatn(2,i,j,k)=0.
         dyhatn(2,i,j,k)=0.
         dzhatn(2,i,j,k)=0.
         dxhatn(3,i,j,k)=0.
         dyhatn(3,i,j,k)=0.
         dzhatn(3,i,j,k)=0.
      endif

                enddo
            enddo
       enddo

!=======================================================================
! now for the initial gauge fields and the electric fields:
!=======================================================================

        do k=1,localsize
            do j=1,localsize
                do i=1,localsize
 
! -- begin: setup variables ------
! (x,y,z) for lattice site:
  
                ilocal = i+lxb-1
                jlocal = j+lyl-1
                klocal = k+lzd-1

                x=float(ilocal)*dx
                y=float(jlocal)*dx
                z=float(klocal)*dx

                xgm=boost*(x-xm)
                ygm=(y-ym)
                zgm=(z-zm)

                xga=boost*(x-xa)
                yga=(y-ya)
                zga=(z-za)

 
! radial coordinate (with boosted x) with origin on monopole:

                rm=sqrt(xgm**2+ygm**2+zgm**2)

! radial coordinate (with boosted x) with origin on antimonopole:

                ra=sqrt(xga**2+yga**2+zga**2)

! xy-radial boosted coordinate with origin on monopole:

                rhom=sqrt(rm**2-zgm**2)

! xy-radial boosted coordinate with origin on antimonopole:

                rhoa=sqrt(ra**2-zga**2)

! dimensionless distances (in units of vector mass):

                sm=gw*vev*rm
                sa=gw*vev*ra

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

! square of scalar field:

       phisq=f(1,i,j,k)**2+f(2,i,j,k)**2+f(3,i,j,k)**2

       if(phisq.ne.0.) then
! W_\mu^1:
       f(4,i,j,k)=0.
       f(5,i,j,k)=2*boost*(-gw*wprofilem*wprofilea
     1  *(hatn(2,i,j,k)*dxhatn(3,i,j,k)-hatn(3,i,j,k)*dxhatn(2,i,j,k)))
       f(6,i,j,k)=(-gw*wprofilem*wprofilea
     1  *(hatn(2,i,j,k)*dyhatn(3,i,j,k)-hatn(3,i,j,k)*dyhatn(2,i,j,k)))
       f(7,i,j,k)=(-gw*wprofilem*wprofilea
     1  *(hatn(2,i,j,k)*dzhatn(3,i,j,k)-hatn(3,i,j,k)*dzhatn(2,i,j,k)))

! W_\mu^2:
       f(8,i,j,k)=0.
       f(9,i,j,k)=2*boost*(-gw*wprofilem*wprofilea
     1  *(hatn(3,i,j,k)*dxhatn(1,i,j,k)-hatn(1,i,j,k)*dxhatn(3,i,j,k)))
       f(10,i,j,k)=(-gw*wprofilem*wprofilea
     1  *(hatn(3,i,j,k)*dyhatn(1,i,j,k)-hatn(1,i,j,k)*dyhatn(3,i,j,k)))
       f(11,i,j,k)=(-gw*wprofilem*wprofilea
     1  *(hatn(3,i,j,k)*dzhatn(1,i,j,k)-hatn(1,i,j,k)*dzhatn(3,i,j,k)))

! W_\mu^3:
       f(12,i,j,k)=0.
       f(13,i,j,k)=2*boost*(-gw*wprofilem*wprofilea
     1  *(hatn(1,i,j,k)*dxhatn(2,i,j,k)-hatn(2,i,j,k)*dxhatn(1,i,j,k)))
       f(14,i,j,k)=(-gw*wprofilem*wprofilea
     1  *(hatn(1,i,j,k)*dyhatn(2,i,j,k)-hatn(2,i,j,k)*dyhatn(1,i,j,k)))
       f(15,i,j,k)=(-gw*wprofilem*wprofilea
     1  *(hatn(1,i,j,k)*dzhatn(2,i,j,k)-hatn(2,i,j,k)*dzhatn(1,i,j,k)))

       else
       f(4,i,j,k)=0.
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
       endif

!        if(ilocal.eq.latx+1.and.jlocal.eq.-1.and.klocal.eq.1) then
!        print *, "latx+1"
!        print *,f(1,i,j,k),f(2,i,j,k),f(3,i,j,k),f(4,i,j,k)
!        print *,f(5,i,j,k),f(6,i,j,k),f(7,i,j,k),f(8,i,j,k)
!        print *,f(9,i,j,k),f(10,i,j,k),f(11,i,j,k),f(12,i,j,k)
!        print *,f(13,i,j,k),f(14,i,j,k),f(15,i,j,k)
!        endif
!        if(ilocal.eq.-latx+1.and.jlocal.eq.-1.and.klocal.eq.1) then
!        print *, "-latx+1"
!        print *,f(1,i,j,k),f(2,i,j,k),f(3,i,j,k),f(4,i,j,k)
!        print *,f(5,i,j,k),f(6,i,j,k),f(7,i,j,k),f(8,i,j,k)
!        print *,f(9,i,j,k),f(10,i,j,k),f(11,i,j,k),f(12,i,j,k)
!        print *,f(13,i,j,k),f(14,i,j,k),f(15,i,j,k)
!        endif



                enddo
            enddo
        enddo
        
!=======================================================================
! gauge functions $\Gamma^a = \partial_i W^a_i$:
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
