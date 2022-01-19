c SO(3) energy:
c
c 15 December 2015:
c Gives same values as energyElectroweak.f for electricEnergyW
c and magneticEnergyW. So staying with the Mathematica generated
c energyElectroweak.f. However, advantage here is that we can 
c also print out the field strengths.
c
      subroutine magneticfield(f,it,hatn,dxhatn,dyhatn,dzhatn,gauge,mf)
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'
      integer i,j,k,it
      integer lx, ly, lz
      real*8 x, y, z, r, dV, totalHelicity, helicityDensity
      real*8 magneticEnergy, magneticEnergyDensity
      real*8 f, gauge, mf, cd
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 dxgauge,dygauge,dzgauge
      real*8 hatn,dxhatn,dyhatn,dzhatn

c
      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
c cd(0:3,number of scalar fields):
      dimension gauge(3,-latx:latx,-laty:laty,-latz:latz)
      dimension mf(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dxgauge(3), dygauge(3), dzgauge(3)
      dimension hatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dxhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dyhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dzhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension fs(3,0:3,0:3)
      dimension cd(0:3,3)
c
c total scalar KE, gradient energy; total gauge electric, magnetic; 
c total scalar potential; total energy.
        totalHelicity=0.
        helicityDensity=0.
        magneticEnergy=0.
        magneticEnergyDensity=0.
c
        lx=latx-1
        ly=laty-1
        lz=latz-1
c
c     volume element:
               dV=dx**3
c
c        print*, ' energy computation sub-box semi-size = ', lx 
c     

       open(unit=300, file='helicityDensity.dat',
     1                         access='append')       
       open(unit=301, file='magneticEnergyDensity.dat',
     1                         access='append')       
 
        do 30 k=-lz,lz
         do 20 j=-ly,ly
            do 10 i=-lx,lx

               x=float(i)*dx
               y=float(j)*dx
               z=float(k)*dx
c     
      r=sqrt(x**2+y**2+z**2)

c     if one wants to only count the energy within a small region:
c     
      call derivatives(f,i,j,k,dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2)
      
      gauge(1,i,j,k)=f(5,i,j,k)*hatn(1,i,j,k)
     1              +f(9,i,j,k)*hatn(2,i,j,k)
     2              +f(13,i,j,k)*hatn(3,i,j,k)

      gauge(2,i,j,k)=f(6,i,j,k)*hatn(1,i,j,k)
     1              +f(10,i,j,k)*hatn(2,i,j,k)
     2              +f(14,i,j,k)*hatn(3,i,j,k)

      gauge(3,i,j,k)=f(7,i,j,k)*hatn(1,i,j,k)
     1              +f(11,i,j,k)*hatn(2,i,j,k)
     2              +f(15,i,j,k)*hatn(3,i,j,k)

      dxgauge(1)=dfdx(5)*hatn(1,i,j,k)+f(5,i,j,k)*dxhatn(1,i,j,k)
     1          +dfdx(9)*hatn(2,i,j,k)+f(9,i,j,k)*dxhatn(2,i,j,k)
     2          +dfdx(13)*hatn(3,i,j,k)+f(13,i,j,k)*dxhatn(3,i,j,k)

      dygauge(1)=dfdy(5)*hatn(1,i,j,k)+f(5,i,j,k)*dyhatn(1,i,j,k)
     1          +dfdy(9)*hatn(2,i,j,k)+f(9,i,j,k)*dyhatn(2,i,j,k)
     2          +dfdy(13)*hatn(3,i,j,k)+f(13,i,j,k)*dyhatn(3,i,j,k)

      dzgauge(1)=dfdz(5)*hatn(1,i,j,k)+f(5,i,j,k)*dzhatn(1,i,j,k)
     1          +dfdz(9)*hatn(2,i,j,k)+f(9,i,j,k)*dzhatn(2,i,j,k)
     2          +dfdz(13)*hatn(3,i,j,k)+f(13,i,j,k)*dzhatn(3,i,j,k)

      dxgauge(2)=dfdx(6)*hatn(1,i,j,k)+f(6,i,j,k)*dxhatn(1,i,j,k)
     1          +dfdx(10)*hatn(2,i,j,k)+f(10,i,j,k)*dxhatn(2,i,j,k)
     2          +dfdx(14)*hatn(3,i,j,k)+f(14,i,j,k)*dxhatn(3,i,j,k)

      dygauge(2)=dfdy(6)*hatn(1,i,j,k)+f(6,i,j,k)*dyhatn(1,i,j,k)
     1          +dfdy(10)*hatn(2,i,j,k)+f(10,i,j,k)*dyhatn(2,i,j,k)
     2          +dfdy(14)*hatn(3,i,j,k)+f(14,i,j,k)*dyhatn(3,i,j,k)

      dzgauge(2)=dfdz(6)*hatn(1,i,j,k)+f(6,i,j,k)*dzhatn(1,i,j,k)
     1          +dfdz(10)*hatn(2,i,j,k)+f(10,i,j,k)*dzhatn(2,i,j,k)
     2          +dfdz(14)*hatn(3,i,j,k)+f(14,i,j,k)*dzhatn(3,i,j,k)

      dxgauge(3)=dfdx(7)*hatn(1,i,j,k)+f(7,i,j,k)*dxhatn(1,i,j,k)
     1          +dfdx(11)*hatn(2,i,j,k)+f(11,i,j,k)*dxhatn(2,i,j,k)
     2          +dfdx(15)*hatn(3,i,j,k)+f(15,i,j,k)*dxhatn(3,i,j,k)

      dygauge(3)=dfdy(7)*hatn(1,i,j,k)+f(7,i,j,k)*dyhatn(1,i,j,k)
     1          +dfdy(11)*hatn(2,i,j,k)+f(11,i,j,k)*dyhatn(2,i,j,k)
     2          +dfdy(15)*hatn(3,i,j,k)+f(15,i,j,k)*dyhatn(3,i,j,k)

      dzgauge(3)=dfdz(7)*hatn(1,i,j,k)+f(7,i,j,k)*dzhatn(1,i,j,k)
     1          +dfdz(11)*hatn(2,i,j,k)+f(11,i,j,k)*dzhatn(2,i,j,k)
     2          +dfdz(15)*hatn(3,i,j,k)+f(15,i,j,k)*dzhatn(3,i,j,k)


c     call fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)
c     call covd(f,i,j,k,hatn,dxhatn,dyhatn,dzhatn,cd)
c     call covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)

c     mf(1,i,j,k)=-(fs(1,2,3)*hatn(1,i,j,k)
c    1            +fs(2,2,3)*hatn(2,i,j,k)
c    2            +fs(3,2,3)*hatn(3,i,j,k)
c    3           - abs(f(19,i,j,k)**3)*(hatn(1,i,j,k)*cd(2,2)*cd(3,3)
c    4           - hatn(1,i,j,k)*cd(2,3)*cd(3,2)
c    5           + hatn(2,i,j,k)*cd(2,3)*cd(3,1)
c    6           - hatn(2,i,j,k)*cd(2,1)*cd(3,3)
c    7           + hatn(3,i,j,k)*cd(2,1)*cd(3,2)
c    8           - hatn(3,i,j,k)*cd(2,2)*cd(3,1)))
 
c     mf(2,i,j,k)=-(fs(1,3,1)*hatn(1,i,j,k)
c    1            +fs(2,3,1)*hatn(2,i,j,k)
c    2            +fs(3,3,1)*hatn(3,i,j,k)
c    3           - abs(f(19,i,j,k)**3)*(hatn(1,i,j,k)*cd(3,2)*cd(1,3)
c    4           - hatn(1,i,j,k)*cd(3,3)*cd(1,2)
c    5           + hatn(2,i,j,k)*cd(3,3)*cd(1,1)
c    6           - hatn(2,i,j,k)*cd(3,1)*cd(1,3)
c    7           + hatn(3,i,j,k)*cd(3,1)*cd(1,2)
c    8           - hatn(3,i,j,k)*cd(3,2)*cd(1,1)))

c     mf(3,i,j,k)=-(fs(1,1,2)*hatn(1,i,j,k)
c    1            +fs(2,1,2)*hatn(2,i,j,k)
c    2            +fs(3,1,2)*hatn(3,i,j,k) 
c    3           - abs(f(19,i,j,k)**3)*(hatn(1,i,j,k)*cd(1,2)*cd(2,3)
c    4           - hatn(1,i,j,k)*cd(1,3)*cd(2,2)
c    5           + hatn(2,i,j,k)*cd(1,3)*cd(2,1)
c    6           - hatn(2,i,j,k)*cd(1,1)*cd(2,3)
c    7           + hatn(3,i,j,k)*cd(1,1)*cd(2,2)
c    8           - hatn(3,i,j,k)*cd(1,2)*cd(2,1)))

c     mf(1,i,j,k)=-(fs(1,2,3)*hatn(1,i,j,k)
c    1            +fs(2,2,3)*hatn(2,i,j,k)
c    2            +fs(3,2,3)*hatn(3,i,j,k)
c    3           - (hatn(1,i,j,k)*cd(2,2)*cd(3,3)
c    4           - hatn(1,i,j,k)*cd(2,3)*cd(3,2)
c    5           + hatn(2,i,j,k)*cd(2,3)*cd(3,1)
c    6           - hatn(2,i,j,k)*cd(2,1)*cd(3,3)
c    7           + hatn(3,i,j,k)*cd(2,1)*cd(3,2)
c    8           - hatn(3,i,j,k)*cd(2,2)*cd(3,1)))
 
c     mf(2,i,j,k)=-(fs(1,3,1)*hatn(1,i,j,k)
c    1            +fs(2,3,1)*hatn(2,i,j,k)
c    2            +fs(3,3,1)*hatn(3,i,j,k)
c    3           - (hatn(1,i,j,k)*cd(3,2)*cd(1,3)
c    4           - hatn(1,i,j,k)*cd(3,3)*cd(1,2)
c    5           + hatn(2,i,j,k)*cd(3,3)*cd(1,1)
c    6           - hatn(2,i,j,k)*cd(3,1)*cd(1,3)
c    7           + hatn(3,i,j,k)*cd(3,1)*cd(1,2)
c    8           - hatn(3,i,j,k)*cd(3,2)*cd(1,1)))

c     mf(3,i,j,k)=-(fs(1,1,2)*hatn(1,i,j,k)
c    1            +fs(2,1,2)*hatn(2,i,j,k)
c    2            +fs(3,1,2)*hatn(3,i,j,k) 
c    3           - (hatn(1,i,j,k)*cd(1,2)*cd(2,3)
c    4           - hatn(1,i,j,k)*cd(1,3)*cd(2,2)
c    5           + hatn(2,i,j,k)*cd(1,3)*cd(2,1)
c    6           - hatn(2,i,j,k)*cd(1,1)*cd(2,3)
c    7           + hatn(3,i,j,k)*cd(1,1)*cd(2,2)
c    8           - hatn(3,i,j,k)*cd(1,2)*cd(2,1)))


c     mf(1,i,j,k)=-(fs(1,2,3)*f(1,i,j,k)
c    1            +fs(2,2,3)*f(2,i,j,k)
c    2            +fs(3,2,3)*f(3,i,j,k)
c    3           - (f(1,i,j,k)*cd(2,2)*cd(3,3)
c    4           - f(1,i,j,k)*cd(2,3)*cd(3,2)
c    5           + f(2,i,j,k)*cd(2,3)*cd(3,1)
c    6           - f(2,i,j,k)*cd(2,1)*cd(3,3)
c    7           + f(3,i,j,k)*cd(2,1)*cd(3,2)
c    8           - f(3,i,j,k)*cd(2,2)*cd(3,1)))
 
c     mf(2,i,j,k)=-(fs(1,3,1)*f(1,i,j,k)
c    1            +fs(2,3,1)*f(2,i,j,k)
c    2            +fs(3,3,1)*f(3,i,j,k)
c    3           - (f(1,i,j,k)*cd(3,2)*cd(1,3)
c    4           - f(1,i,j,k)*cd(3,3)*cd(1,2)
c    5           + f(2,i,j,k)*cd(3,3)*cd(1,1)
c    6           - f(2,i,j,k)*cd(3,1)*cd(1,3)
c    7           + f(3,i,j,k)*cd(3,1)*cd(1,2)
c    8           - f(3,i,j,k)*cd(3,2)*cd(1,1)))

c     mf(3,i,j,k)=-(fs(1,1,2)*f(1,i,j,k)
c    1            +fs(2,1,2)*f(2,i,j,k)
c    2            +fs(3,1,2)*f(3,i,j,k) 
c    3           - (f(1,i,j,k)*cd(1,2)*cd(2,3)
c    4           - f(1,i,j,k)*cd(1,3)*cd(2,2)
c    5           + f(2,i,j,k)*cd(1,3)*cd(2,1)
c    6           - f(2,i,j,k)*cd(1,1)*cd(2,3)
c    7           + f(3,i,j,k)*cd(1,1)*cd(2,2)
c    8           - f(3,i,j,k)*cd(1,2)*cd(2,1)))


      mf(1,i,j,k)=  -(dygauge(3)-dzgauge(2)
     1           - (hatn(1,i,j,k)*dyhatn(2,i,j,k)*dzhatn(3,i,j,k)
     2           - hatn(1,i,j,k)*dyhatn(3,i,j,k)*dzhatn(2,i,j,k)
     3           + hatn(2,i,j,k)*dyhatn(3,i,j,k)*dzhatn(1,i,j,k)
     4           - hatn(2,i,j,k)*dyhatn(1,i,j,k)*dzhatn(3,i,j,k)
     5           + hatn(3,i,j,k)*dyhatn(1,i,j,k)*dzhatn(2,i,j,k)
     6           - hatn(3,i,j,k)*dyhatn(2,i,j,k)*dzhatn(1,i,j,k)))

      mf(2,i,j,k)=-(dzgauge(1)-dxgauge(3)
     1           - (hatn(1,i,j,k)*dzhatn(2,i,j,k)*dxhatn(3,i,j,k)
     2           - hatn(1,i,j,k)*dzhatn(3,i,j,k)*dxhatn(2,i,j,k)
     3           + hatn(2,i,j,k)*dzhatn(3,i,j,k)*dxhatn(1,i,j,k)
     4           - hatn(2,i,j,k)*dzhatn(1,i,j,k)*dxhatn(3,i,j,k)
     5           + hatn(3,i,j,k)*dzhatn(1,i,j,k)*dxhatn(2,i,j,k)
     6           - hatn(3,i,j,k)*dzhatn(2,i,j,k)*dxhatn(1,i,j,k)))

      mf(3,i,j,k)=(dygauge(1)-dxgauge(2)
     1           - (hatn(1,i,j,k)*dyhatn(2,i,j,k)*dxhatn(3,i,j,k)
     2           - hatn(1,i,j,k)*dyhatn(3,i,j,k)*dxhatn(2,i,j,k)
     3           + hatn(2,i,j,k)*dyhatn(3,i,j,k)*dxhatn(1,i,j,k)
     4           - hatn(2,i,j,k)*dyhatn(1,i,j,k)*dxhatn(3,i,j,k)
     5           + hatn(3,i,j,k)*dyhatn(1,i,j,k)*dxhatn(2,i,j,k)
     6           - hatn(3,i,j,k)*dyhatn(2,i,j,k)*dxhatn(1,i,j,k)))

c     mf(1,i,j,k)=  -(dygauge(3)-dzgauge(2)
c    1           - (f(1,i,j,k)*dfdy(2)*dfdz(3)
c    2           - f(1,i,j,k)*dfdy(3)*dfdz(2)
c    3           + f(2,i,j,k)*dfdy(3)*dfdz(1)
c    4           - f(2,i,j,k)*dfdy(1)*dfdz(3)
c    5           + f(3,i,j,k)*dfdy(1)*dfdz(2)
c    6           - f(3,i,j,k)*dfdy(2)*dfdz(1)))

c     mf(2,i,j,k)=-(dzgauge(1)-dxgauge(3)
c    1           - (f(1,i,j,k)*dfdz(2)*dfdx(3)
c    2           - f(1,i,j,k)*dfdz(3)*dfdx(2)
c    3           + f(2,i,j,k)*dfdz(3)*dfdx(1)
c    4           - f(2,i,j,k)*dfdz(1)*dfdx(3)
c    5           + f(3,i,j,k)*dfdz(1)*dfdx(2)
c    6           - f(3,i,j,k)*dfdz(2)*dfdx(1)))

c     mf(3,i,j,k)=(dygauge(1)-dxgauge(2)
c    1           - (f(1,i,j,k)*dfdy(2)*dfdx(3)
c    2           - f(1,i,j,k)*dfdy(3)*dfdx(2)
c    3           + f(2,i,j,k)*dfdy(3)*dfdx(1)
c    4           - f(2,i,j,k)*dfdy(1)*dfdx(3)
c    5           + f(3,i,j,k)*dfdy(1)*dfdx(2)
c    6           - f(3,i,j,k)*dfdy(2)*dfdx(1)))


        helicityDensity=gauge(1,i,j,k)*mf(1,i,j,k)
     1                 +gauge(2,i,j,k)*mf(2,i,j,k)
     2                 +gauge(3,i,j,k)*mf(3,i,j,k)
        magneticEnergyDensity= 0.5*(mf(1,i,j,k)**2+mf(2,i,j,k)**2
     1                          +mf(3,i,j,k)**2)  
       
        totalHelicity = totalHelicity+helicityDensity*dV
        magneticEnergy= magneticEnergy+magneticEnergyDensity*dV


        if(it.eq.nt) then
        write(300,*) r, helicityDensity
        endif
         if(it.eq.nt) then
        write(301,*) r, magneticEnergyDensity
        endif
 
10         continue
20       continue
30     continue
c     
c    
       
        close(300)
        close(301)

       open(unit=200, file='totalHelicity.dat',
     1                         access='append')       
       open(unit=201, file='magneticEnergy.dat',
     1                         access='append')       
 
        write(200,*) it, totalHelicity
        write(201,*) it, magneticEnergy

        close(200)
        close(201)


      return
      end
