!=======================================================================
! Just a monopole.
!
! Gauged SO(3) model. Initial conditions.
!
! 15 April 2015:
! 
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
     1                            lxl,lxr)
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'
      integer i,j,k,n
      integer itw,izm

!      real*8 zm (now declared in initialParameters.inc)
      real*8 x,y,z,xgm,ygm,zgm,rm,rhom
      real*8 sm,ms,mv
      real*8 pprofilem,wprofilem
      real*8 f,fgamma
      real*8 hatn
      real*8 dxhatn,dyhatn,dzhatn
      real*8 phisq
      real*8 dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2
      real*8 twist
      real*8 correctionm
      integer pid,lxl,lxr,ilocal
 

      dimension f(nf,1:localsize,-laty:laty,-latz:latz)
      dimension fgamma(nf,6,-laty:laty,-latz:latz)

      dimension hatn(3,1:localsize,-laty:laty,-latz:latz)
      dimension dxhatn(3,1:localsize,-laty:laty,-latz:latz)
      dimension dyhatn(3,1:localsize,-laty:laty,-latz:latz)
      dimension dzhatn(3,1:localsize,-laty:laty,-latz:latz)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)

!=======================================================================
! This section is for running the code on several nodes at the
! same time, each code with a different value of the parameters.
! read scan parameters and create files to write out successful runs 
!=======================================================================

        character*500 filename2
        character (len=200) winfile
        character (len=200) izmono
      
        call getarg(1,winfile)
        read(winfile,'(a500)')filename2

        call getarg(2,izmono)
        read(izmono,*) izm
        zm=(float(izm)+0.5)*dx
     
        if(pid==0) then
        print*,' ======================= '
        print*,' lambda, zm, dx, latx*dx ', lambda,zm,dx,float(latx)*dx
        print *,' '
        endif

!=======================================================================
!       initialize:
!=======================================================================

        do k=-latz,latz
            do j=-laty,laty
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
! start: loop over larger lattice (+1) for scalar fields:
!=======================================================================

        do k=-latz,latz
            do j=-laty,laty
                do i=1,localsize

! ------------------ begin: setup variables -----------------
!               (x,y,z) for lattice site:

                ilocal = i+lxl-1
                x=float(ilocal)*dx
                y=float(j)*dx
                z=float(k)*dx

                xgm=x-xm
                ygm=y-ym
                zgm=z-zm

! radial coordinate (with boosted x) with origin on monopole:

                rm=sqrt(xgm**2+ygm**2+zgm**2)

! xy-radial boosted coordinate with origin on monopole:
       
                rhom=sqrt(rm**2-zgm**2)

! dimensionless distances (in units of vector mass):

                sm=gw*vev*rm

! -- end: setup variables ------
!
!--begin scalar profile -----
! monopole and antimonopole scalar profile functions:

                pprofilem=vev*(1./tanh(sm)-(1.+ms*sm)*exp(-ms*sm)/sm)

!--end scalar profile -----
!
! \phi^a:

        if(rm.gt.0.001) then
                hatn(1,i,j,k)=xgm/rm
                hatn(2,i,j,k)=ygm/rm
                hatn(3,i,j,k)=zgm/rm

                dxhatn(1,i,j,k)=1./rm-xgm**2/rm**3
                dyhatn(1,i,j,k)=-xgm*ygm/rm**3
                dzhatn(1,i,j,k)=-xgm*zgm/rm**3
                dxhatn(2,i,j,k)=-ygm*xgm/rm**3
                dyhatn(2,i,j,k)=1./rm-ygm**2/rm**3
                dzhatn(2,i,j,k)=-ygm*zgm/rm**3
                dxhatn(3,i,j,k)=-zgm*xgm/rm**3
                dyhatn(3,i,j,k)=-zgm*ygm/rm**3
                dzhatn(3,i,j,k)=1./rm-zgm**2/rm**3

                f(1,i,j,k)=pprofilem*hatn(1,i,j,k)
                f(2,i,j,k)=pprofilem*hatn(2,i,j,k)
                f(3,i,j,k)=pprofilem*hatn(3,i,j,k)

                f(19,i,j,k)=pprofilem/vev
       else
                hatn(1,i,j,k)=0.
                hatn(2,i,j,k)=0.
                hatn(3,i,j,k)=1.

                f(1,i,j,k)=0.
                f(2,i,j,k)=0.
                f(3,i,j,k)=0.
                f(19,i,j,k)=0.
       endif
                enddo
            enddo
       enddo

!=======================================================================
! now for the initial gauge fields and the electric fields:
!=======================================================================

        do k=-latz,latz
            do j=-laty,laty
                do i=1,localsize
 
! -- begin: setup variables ------
! (x,y,z) for lattice site:

                ilocal = i+lxl-1
                x=float(ilocal)*dx
                y=float(j)*dx
                z=float(k)*dx

                xgm=(x-xm)
                ygm=(y-ym)
                zgm=(z-zm)
 
! radial coordinate (with boosted x) with origin on monopole:

                rm=sqrt(xgm**2+ygm**2+zgm**2)

! xy-radial boosted coordinate with origin on monopole:

                rhom=sqrt(rm**2-zgm**2)

! dimensionless distances (in units of vector mass):

                sm=gw*vev*rm

! monopole gauge profile function (1-K(r)):

        correctionm=(1.+(1.-sqrt(lambda))*sm**2/4.+sm**4/16.)/
     1                (1.+sm**2/4.+sm**4/16.)
        wprofilem=(1.-(sm/sinh(sm))*correctionm)/gw

! W_\mu^1:

       f(4,i,j,k)=0.
       f(5,i,j,k)=-wprofilem
     1  *(hatn(2,i,j,k)*dxhatn(3,i,j,k)-hatn(3,i,j,k)*dxhatn(2,i,j,k))
       f(6,i,j,k)=-wprofilem
     1  *(hatn(2,i,j,k)*dyhatn(3,i,j,k)-hatn(3,i,j,k)*dyhatn(2,i,j,k))
       f(7,i,j,k)=-wprofilem
     1  *(hatn(2,i,j,k)*dzhatn(3,i,j,k)-hatn(3,i,j,k)*dzhatn(2,i,j,k))

!       f(5,i,j,k)=0.
!       f(6,i,j,k)=+gw*wprofilem*zgm/rm**2
!       f(7,i,j,k)=-gw*wprofilem*ygm/rm**2
! W_\mu^2:

       f(8,i,j,k)=0.
       f(9,i,j,k)=-wprofilem
     1  *(hatn(3,i,j,k)*dxhatn(1,i,j,k)-hatn(1,i,j,k)*dxhatn(3,i,j,k))
       f(10,i,j,k)=-wprofilem
     1  *(hatn(3,i,j,k)*dyhatn(1,i,j,k)-hatn(1,i,j,k)*dyhatn(3,i,j,k))
       f(11,i,j,k)=-wprofilem
     1  *(hatn(3,i,j,k)*dzhatn(1,i,j,k)-hatn(1,i,j,k)*dzhatn(3,i,j,k))
!       f(9,i,j,k)=-gw*wprofilem*zgm/rm**2
!       f(10,i,j,k)=0.
!       f(11,i,j,k)=+gw*wprofilem*xgm/rm**2
! W_\mu^3:

       f(12,i,j,k)=0.
       f(13,i,j,k)=-wprofilem
     1  *(hatn(1,i,j,k)*dxhatn(2,i,j,k)-hatn(2,i,j,k)*dxhatn(1,i,j,k))
       f(14,i,j,k)=-wprofilem
     1  *(hatn(1,i,j,k)*dyhatn(2,i,j,k)-hatn(2,i,j,k)*dyhatn(1,i,j,k))
       f(15,i,j,k)=-wprofilem
     1  *(hatn(1,i,j,k)*dzhatn(2,i,j,k)-hatn(2,i,j,k)*dzhatn(1,i,j,k))
!       f(13,i,j,k)=+gw*wprofilem*ygm/rm**2
!       f(14,i,j,k)=-gw*wprofilem*xgm/rm**2
!       f(15,i,j,k)=0.

                enddo
            enddo
        enddo
        
!=======================================================================
! gauge functions $\Gamma^a = \partial_i W^a_i$:
! This uses a call to derivatives. That is why it is necessary to
! evaluate the gauge functions in a separate do loop (after all the
! gauge fields have been specified).
!=======================================================================

        
        if(pid==0) then
        do k=-latz,latz
            do j=-laty,laty
                do i=1,localsize-6

       call derivatives(f,i,j,k,dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2,lxl)

                f(16,i,j,k)=dfdx(5)+dfdy(6)+dfdz(7)
                f(17,i,j,k)=dfdx(9)+dfdy(10)+dfdz(11)
                f(18,i,j,k)=dfdx(13)+dfdy(14)+dfdz(15)
               
                enddo
            enddo
        enddo

        elseif(pid==np-1) then
        do k=-latz,latz
            do j=-laty,laty
                do i=7,localsize

       call derivatives(f,i,j,k,dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2,lxl)

                f(16,i,j,k)=dfdx(5)+dfdy(6)+dfdz(7)
                f(17,i,j,k)=dfdx(9)+dfdy(10)+dfdz(11)
                f(18,i,j,k)=dfdx(13)+dfdy(14)+dfdz(15)
               
                enddo
            enddo
        enddo

        else
        do k=-latz,latz
            do j=-laty,laty
                do i=4,localsize-3

       call derivatives(f,i,j,k,dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2,lxl)

                f(16,i,j,k)=dfdx(5)+dfdy(6)+dfdz(7)
                f(17,i,j,k)=dfdx(9)+dfdy(10)+dfdz(11)
                f(18,i,j,k)=dfdx(13)+dfdy(14)+dfdz(15)
               
                enddo
            enddo
        enddo

        endif



!       need ghostfields here if(pid.ne.0.and.pid.ne.np-1) then
!      if(pid==0.and.i==18.and.j==-29.and.k==-29) then
!      print *,pid,f(17,i-1,j,k),f(17,i+3,j,k),dfdx(9)
!      endif
!      if(pid==1.and.i==5.and.j==-29.and.k==-29) then
!      print *,pid,f(17,i,j,k),f(17,i+3,j,k),dfdx(9)
!      endif
 

        return
        end
