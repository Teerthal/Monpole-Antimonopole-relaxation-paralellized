c Here we list the problem-specific parameters. These will be used
c in the initial conditions, in detecting features, and in main.f.
c latx, dx have been defined in parameters.inc and we use those here.
      integer im,jm
      real*8 xm,ym,zm,vxm,vym
      parameter(im=0,jm=0)
c location of monopole (zm is specified as an argument when executing):
      parameter(xm=float(im)*dx,ym=float(jm)*dx)
      parameter(vxm=0.,vym=0.)
