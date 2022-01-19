c For SO(3) ('t Hooft-Polyakov) model: conventions in 
c MMbar scattering paper arXiv:1511.05095
c
c nf=number of fields; nt=number of time steps;
c lat=semi-lattice size; dx=spatial spacing;
c dt=temporal spacing; alpha, beta etc. are parameters of the model
c
c Allows for different lattice sizes along different axes but will need
c to check code where latx=laty=latz may have been assumed. 
c
      integer nf,latx,laty,latz,nt,nsnaps,outerlayer
      integer nscalar,ngauge,nconstraints,nprocessx,localsize
      real*8 dx,dt,gw,gy,gp2,lambda,vev
      real*8 pi
c For abelian-Higgs model: gy is the gauge coupling.
c Note: conventions like those in electroweak but truncated to hypercharge.
c model parameters:
c the number of fields (nf) includes real scalar fields, gauge fields,
c and the gauge constraints. E.g. for electroweak there are 4 real
c scalars, 4x4=16 gauge fields, and 4 gauge constraints (Lorenz gauge
c conditions. (gauge constraint equations are first order.)
c
      parameter(pi=3.14159265359)
c also specify how many real scalar fields there are, how many gauge
c fields (e.g. electroweak has 16), and how many gauge constraints:
c (this is used when implementing absorbing boundary conditions)
      parameter(nf=18,nscalar=3,ngauge=3*4,nconstraints=ngauge/4)
      parameter(gw=1.0,gy=0.,lambda=1.0,vev=1.0)
c
c simulation parameters next:
c gp2=gauge function parameter ($gp2=a^2$ of Num. Rel., Baumgarte & Shapiro,
c see Chapter 11). $0 < gp2 < 1$.
c
      parameter(latx=32,laty=latx,latz=latx,nprocessx=2,nt=0,
     1      dx=0.1,dt=dx/4.,localsize=(2*latx/nprocessx)+13,gp2=0.75)
c
c Output parameter -- how many snapshots of the energy distribution?
      parameter(nsnaps=20)
c
c Calculate helicity in a sub-box with dimension -lat+outerlayer,+lat-outerlayer.
c outerlayer must be greater than equal to 1.
      parameter(outerlayer=2)
c
