!======================================================================
! For electroweak model:
! (following conventions in Achucarro and Vachaspati,
!  Phys. Rept. 327 (2000) 347-426.)
!
! nf=number of fields; nt=number of time steps;
! lat=semi-lattice size; dx=spatial spacing;
! alpha, beta etc. are parameters of the model
!
! Allows for different lattice sizes along different axes. 
!
! Fiducial lattice only goes up to lat-1. The last layer on the 
! lattice is a boundary layer where only some variables are calculated.
!======================================================================

      integer nf,latx,laty,latz,nt,nte,nsnaps,outerlayer,localsize
      integer nscalar,ngauge,nconstraints,nprocessx
      real*8 dx,dt,gw,gy,gp2,lambda,vev,relaxparam
      real*8 pi
      real*8 vel,boost

!======================================================================
! For abelian-Higgs model: gy is the gauge coupling.
! Note: conventions like those in electroweak but truncated to hypercharge.
! model parameters:
!
! The number of fields (nf) includes real scalar fields, gauge fields,
! and the gauge constraints. E.g. for electroweak there are 4 real
! scalars, 4x4=16 gauge fields, and 4 gauge constraints (Lorenz gauge
! conditions. (gauge constraint equations are first order.)
!======================================================================

      parameter(pi=3.14159265359)
      parameter(nf=24,nscalar=4,ngauge=4*4,nconstraints=ngauge/4)
      parameter(gw=1.0,gy=1.0,lambda=1.0,vev=1.0)

!======================================================================
! Simulation parameters next:
! gp2=gauge function parameter ($gp2=a^2$ of Num. Rel., Baumgarte & Shapiro,
! see Chapter 11). $0 < gp2 < 1$. 
! relaxparam = relaxation parameter in "successive over relaxation". (Similar
! to dt in the evolution code.)
! e.g. eq.(107) in http://www.damtp.cam.ac.uk/lab/people/sd/lectures/nummeth98/pdes.htm
! And, np = number of processors.
!======================================================================

      parameter(latx=32,laty=latx,latz=latx,nprocessx=2,
     1 nt=10,nte=1200,dx=0.1,dt=dx/8.,localsize=(2*latx/nprocessx)+13,
     2  	gp2=0.75,relaxparam=1.0,vel=0.0,boost=1.0/sqrt(1.0-vel*vel))

!======================================================================
! Output parameter -- how many snapshots of the energy distribution?
!======================================================================

      parameter(nsnaps=nt/2)

!======================================================================
! Calculate energy in a sub-box with dimension -lat+outerlayer,+lat-outerlayer.
! outerlayer must be greater than equal to 1 (see energy subroutine).
!======================================================================

      parameter(outerlayer=2)

