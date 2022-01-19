!=======================================================================
! SO(3) energy:
!
! 15 December 2015:
! Gives same values as energyElectroweak.f for electricEnergyW
! and magneticEnergyW. So staying with the Mathematica generated
! energyElectroweak.f. However, advantage here is that we can 
! also print out the field strengths.
!=======================================================================

      subroutine uenergy(f,energy,itime,pid,lxb,lyl,lzd,
     1                          bb,fb,lb,rb,db,ub)
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'

      integer n,i,j,k,itime,pid,proci,procj,prock
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub
      integer ilocal,jlocal,klocal,enbb,enfb,enlb,enrb,endb,enub
      integer lx, ly, lz
!     integer energyfile,potentialfile,rmagneticfile,fieldfile
      real*8 x, y, z, dV
      real*8 kineticEnergyPhi, gradientEnergyPhi, electricEnergyW
      real*8 magneticEnergyW, electricEnergyY, magneticEnergyY
      real*8 potentialEnergy, energyDensity
      real*8 totalKEPhi, totalGEPhi, totalEEW, totalMEW, totalEEY
      real*8 totalMEY, totalPE, totalEnergy
      real*8 f
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 dfddx,dfddy,dfddz
      real*8 cd
      real*8 energyInSphere,r, energy(12), energyBoundary
      character(len=20) :: fileTE,fileME,fileGE,fileEE,filePE,fileKE

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
      dimension fs(3,0:3,0:3),d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)

!=======================================================================
! Covariant Derivatives:  cd(0:3,number of scalar fields):
!=======================================================================

      dimension cd(0:3,3)

!---Identifying the coordinates of the sublattice
!---Prefix l : sublattice coordinates
!---Prefix en : global lattice coordinates to constructed to omit
!---last lattice site 

      call processcoord(pid,proci,procj,prock)

       enbb=bb
       if(proci==nprocessx-1) then
              enfb=fb
       else
              enfb=fb-1
       endif

       enlb=lb
       if(procj==nprocessx-1) then
              enrb=rb
       else
              enrb=rb-1
       endif

       endb=db
       if(prock==nprocessx-1) then
              enub=ub
       else
              enub=ub-1
       endif
      
        
!=======================================================================
! total scalar KE, gradient energy; total gauge electric, magnetic; 
! total scalar potential; total energy.
!=======================================================================

        totalKEPhi=0.
        totalGEPhi=0.
        totalEEW=0.
        totalMEW=0.
        totalEEY=0.
        totalMEY=0.
        totalPE=0.
        totalEnergy=0.
        energyInSphere=0.
        energyBoundary=0.

!---lx,ly,lz seem to be phased out. Omitting
!        lx=latx-1
!        ly=laty-1
!        lz=latz-1

!=======================================================================
!     volume element:
!=======================================================================

        dV=dx**3

!        print*, ' energy computation sub-box semi-size = ', lx 

!       if(pid==1) then

!       open(unit=100, file='energyDensity.dat',
!    1                         access='append')       
!       open(unit=101, file='kineticEnergyPhi.dat',
!    1                         access='append')       
!       open(unit=102, file='gradientEnergyPhi.dat',
!    1                         access='append')       
!       open(unit=103, file='electricEnergyW.dat',
!    1                         access='append')       
!       open(unit=104, file='magneticEnergyW.dat',
!    1                         access='append')       
!       open(unit=105, file='electricEnergyY.dat',
!    1                         access='append')       
!       open(unit=106, file='magneticEnergyY.dat',
!    1                         access='append')       
!       open(unit=107, file='potentialEnergy.dat',
!    1                         access='append')       
!
!     endif
!        if(itime.eq.nt) then
!       write(fileTE,"(A3,I2.2,A1,I4.4)") "TE.",pid,".",itime
!       write(filePE,"(A3,I2.2,A1,I4.4)") "PE.",pid,".",itime
!       write(fileME,"(A3,I2.2,A1,I4.4)") "ME.",pid,".",itime
!       write(fileGE,"(A3,I2.2,A1,I4.4)") "GE.",pid,".",itime
!       write(fileEE,"(A3,I2.2,A1,I4.4)") "EE.",pid,".",itime
!       write(fileKE,"(A3,I2.2,A1,I4.4)") "KE.",pid,".",itime

!       open(unit=10,file="data/"//trim(adjustl(fileTE)),
!    1          form="unformatted",access="stream")
!       open(unit=11,file="data/"//trim(adjustl(filePE)),
!    1          form="unformatted",access="stream")
!       open(unit=12,file="data/"//trim(adjustl(fileME)),
!    1          form="unformatted",access="stream")
!       open(unit=13,file="data/"//trim(adjustl(fileGE)),
!    1          form="unformatted",access="stream")
!       open(unit=14,file="data/"//trim(adjustl(fileEE)),
!    1          form="unformatted",access="stream")
!       open(unit=15,file="data/"//trim(adjustl(fileKE)),
!    1          form="unformatted",access="stream")

!       endif
        do k=endb,enub
          do j=enlb,enrb
            do i=enbb,enfb
     
!     if one wants to only count the energy within a small region:
               ilocal = i+lxb-1
               jlocal = j+lyl-1
               klocal = k+lzd-1

               x=float(ilocal)*dx
               y=float(jlocal)*dx
               z=float(klocal)*dx

      call derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)
!     call derivativesPBC(f,i,j,k,dfdx,dfdy,dfdz,
!    1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)


!=======================================================================
! covariant derivative: cd(mu,a)=(D_\mu\Phi)^a (a=1 is the real 
! part of the upper component of the doublet; a=2 is the imag
! part; a=3 is real part of the lower component; a=4 is the
! imag part of the lower component.
!=======================================================================
         
        call covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)

!=======================================================================     
!find gauge field strengths (uses derivatives):
!=======================================================================

        call fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)



        kineticEnergyPhi= 0.5*(cd(0,1)**2+cd(0,2)**2+cd(0,3)**2)

        gradientEnergyPhi= 0.5*(cd(1,1)**2+cd(2,1)**2+cd(3,1)**2  
     1                       +cd(1,2)**2+cd(2,2)**2+cd(3,2)**2  
     2                       +cd(1,3)**2+cd(2,3)**2+cd(3,3)**2)

        electricEnergyW=0.5*(fs(1,0,1)**2+fs(1,0,2)**2+fs(1,0,3)**2
     1                    +fs(2,0,1)**2+fs(2,0,2)**2+fs(2,0,3)**2
     1                    +fs(3,0,1)**2+fs(3,0,2)**2+fs(3,0,3)**2)

        magneticEnergyW=0.5*(fs(1,1,2)**2+fs(1,1,3)**2+fs(1,2,3)**2
     1                    +fs(2,1,2)**2+fs(2,1,3)**2+fs(2,2,3)**2
     2                    +fs(3,1,2)**2+fs(3,1,3)**2+fs(3,2,3)**2)

        electricEnergyY=0.
        magneticEnergyY=0.

        potentialEnergy=0.25*lambda*(f(1,i,j,k)**2+f(2,i,j,k)**2
     1                       +f(3,i,j,k)**2-vev**2)**2

! ======================================================================
     
       energyDensity=kineticEnergyPhi+gradientEnergyPhi+electricEnergyW+
     1   magneticEnergyW+electricEnergyY+magneticEnergyY+potentialEnergy

!       if(gradientEnergyPhi.ge.10.0) then
!       if(ilocal.eq.latx.and.jlocal.eq.0.and.klocal.eq.0) then
!       print *, "latx"
!       print *, "cd   ",cd(1,1),cd(2,1),cd(3,1)
!       print *, "cd   ",cd(1,2),cd(2,2),cd(3,2)
!       print *, "cd   ",cd(1,3),cd(2,3),cd(3,3)
!       print *, f(1,i,j,k),f(2,i,j,k),f(3,i,j,k)
!       print*,ilocal,jlocal,klocal,pid
!       print*, dfdx(1),dfdx(2),dfdx(3)
!       print*,f(1,i,j,k)
!       print*,f(1,i+1,j,k),f(1,i+2,j,k),f(1,i+3,j,k)
!       print*,f(1,i-1,j,k),f(1,i-2,j,k),f(1,i-3,j,k)
!        print*,cd(1,1),dfdx(1),f(9,i,j,k),f(3,i,j,k),f(13,i,j,k),
!    1                      f(2,i,j,k),dfdx(2),dfdx(3)
!       print *, gradientEnergyPhi
!       endif
!        if(ilocal.eq.-latx.and.jlocal.eq.0.and.klocal.eq.0) then
!!       print *, "latx"
!        print *, "cd   ",cd(1,1),cd(2,1),cd(3,1)
!        print *, "cd   ",cd(1,2),cd(2,2),cd(3,2)
!        print *, "cd   ",cd(1,3),cd(2,3),cd(3,3)
!        print *, f(1,i,j,k),f(2,i,j,k),f(3,i,j,k) 
!        print*,ilocal,jlocal,klocal,pid
!       print *, "cd   ",cd(0,1),cd(0,2),cd(0,3),cd(1,1)
!!       print*, dfdx(1),dfdx(2),dfdx(3)
!!       print*,f(1,i,j,k)
!!       print*,f(1,i+1,j,k),f(1,i+2,j,k),f(1,i+3,j,k)
!!       print*,f(1,i-1,j,k),f(1,i-2,j,k),f(1,i-3,j,k)
!!       print*,cd(1,1),dfdx(1),f(9,i,j,k),f(3,i,j,k),f(13,i,j,k),
!!    1                      f(2,i,j,k),dfdx(2),dfdx(3)
!!       print *, gradientEnergyPhi
!        print *, kineticEnergyPhi
!        endif
!       if(ilocal.eq.-latx.and.jlocal.eq.-1.and.klocal.eq.1) then
!       print *, "-latx"
!       print *, "cd   ",cd(1,1),cd(2,1),cd(3,1)
!       print *, "cd   ",cd(1,2),cd(2,2),cd(3,2)
!       print *, "cd   ",cd(1,3),cd(2,3),cd(3,3)
!       print *, f(1,i,j,k),f(2,i,j,k),f(3,i,j,k)
!       print *, gradientEnergyPhi
!       endif

     
        totalKEPhi=totalKEPhi+kineticEnergyPhi*dV
        totalGEPhi=totalGEPhi+gradientEnergyPhi*dV
        totalEEW=totalEEW+electricEnergyW*dV
        totalMEW=totalMEW+magneticEnergyW*dV
        totalEEY=totalEEY+electricEnergyY*dV
        totalMEY=totalMEY+magneticEnergyY*dV
        totalPE=totalPE+potentialEnergy*dV
        totalEnergy=totalEnergy+energyDensity*dV

!=======================================================================
! energy within sphere -- useful for single monopole:
!        r=sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)
! energy within sphere -- for mmbar the center is at origin:
!=======================================================================

        r=sqrt(x**2+y**2+z**2)

        if(r.lt.float(latx-1)*dx) then
        energyInSphere=energyInSphere+energyDensity*dV
        endif

        if((abs(ilocal).ge.(latx)).or.(abs(jlocal).ge.(latx))
     1              .and.(abs(klocal).ge.(latx))) then
        energyBoundary=energyBoundary+energyDensity*dV
        endif

!       if(pid==1) then

!       if(itime.eq.nt) then
!       write(100,*) r, energyDensity
!       endif
!       
!       if(itime.eq.nt) then
!       write(101,*) r, kineticEnergyPhi
!       endif

!       if(itime.eq.nt) then
!       write(102,*) r, gradientEnergyPhi
!       endif

!       if(itime.eq.nt) then
!       write(103,*) r, electricEnergyW
!       endif
!       
!       if(itime.eq.nt) then
!       write(104,*) r, magneticEnergyW
!       endif
!       
!       if(itime.eq.nt) then
!       write(105,*) r, electricEnergyY
!       endif
!      
!       if(itime.eq.nt) then
!       write(106,*) r, magneticEnergyY
!       endif
!       
!       if(itime.eq.nt) then
!       write(107,*) r, potentialEnergy
!       endif
!       
!      endif 
!
!        if(itime.eq.nt) then
!         write(10) energyDensity
!         write(11) potentialEnergy
!         write(12) magneticEnergyW
!         write(13) gradientEnergyPhi
!         write(14) electricEnergyW
!         write(15) kineticEnergyPhi
!        endif
                enddo
            enddo
        enddo

        energy(1)=itime
        energy(2)=totalKEPhi
        energy(3)=totalGEPhi
        energy(4)=totalEEW
        energy(5)=totalMEW
        energy(6)=totalEEY
        energy(7)=totalMEY
        energy(8)=totalPE
        energy(9)=totalEnergy
        energy(10)=float(latx-1)*dx
        energy(11)=energyInSphere
        energy(12)=energyBoundary

!        if(itime.eq.nt) then
!       close(10)
!       close(11)
!       close(12)
!       close(13)
!       close(14)
!       close(15)
!        endif

          return
      end
