c Here we list the problem-specific parameters. These will be used
c in the initial conditions, in detecting features, and in main.f.
c latx, dx have been defined in parameters.inc and we use those here.
      real*8 xm,ym,zm,vxm,vym,vzm,twist
      real*8 x0,y0,z0,vx0,vy0,vz0
      real*8 width,omega,amp
c      parameter(xm=0.0*dx,ym=0.0*dx,zm=(0.7*float(latx)+0.5)*dx)
c      parameter(xm=0.0*dx,ym=0.0*dx,zm=(12.0+0.5)*dx)
c      parameter(vxm=0.,vym=0.,vzm=0.0,twist=0.0*3.1416)
c
c parameters for gauge collisions (to investigate mmbar production):
c location of initial gauge packets:
      parameter(x0=0.0*dx,y0=0.0*dx,z0=16.0*dx)
c width of gauge packet:
      parameter(width=4.0*dx)
c velocity of gauge packets (not needed for so3.icGaugeWaves.f):
      parameter(vx0=0.,vy0=0.,vz0=0.6)
c circular polarization frequency:
      parameter(omega=4.)
c amplitude of incoming waves:
      parameter(amp=10.0)
c
