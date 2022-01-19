c 24 June 2016: Energies for SO(3) model.
c
      subroutine eenergy(f,fd,energy,itime,pid,lxb,lyl,lzd,
     1                          bb,fb,lb,rb,db,ub)
      implicit none
      include 'parameters.inc'
      integer n,i,j,k,itime,pid,proci,procj,prock
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub
      integer ilocal,jlocal,klocal,enbb,enfb,enlb,enrb,endb,enub
      integer lx, ly, lz
      real*8 x, y, z, dV
      real*8 kineticEnergyPhi, gradientEnergyPhi, electricEnergyW
      real*8 magneticEnergyW, electricEnergyY, magneticEnergyY
      real*8 potentialEnergy, energyDensity
      real*8 totalKEPhi, totalGEPhi, totalEEW, totalMEW
c      real*8 totalMEY, totalEEY
      real*8 totalPE, totalEnergy
      real*8 f, fd
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 dfddx,dfddy,dfddz
      real*8 cd
      real*8 magneticCharge,totalCharge,chargeInOctant
      real*8 energy(7)
      character(len=20) :: fileTE,filePE,fileTEE

c       ,fileME,fileGE,fileEE,fileKE
c
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension fd(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
      dimension fs(ngauge/4,0:3,0:3),d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
      dimension cd(0:3,nscalar)

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

c     
c total scalar KE, gradient energy; total gauge electric, magnetic; 
c total scalar potential; total energy.
        totalKEPhi=0.
        totalGEPhi=0.
        totalEEW=0.
        totalMEW=0.
        totalPE=0.
        totalEnergy=0.
c
        totalCharge=0.
        chargeInOctant=0.
c
        lx=latx
        ly=laty
        lz=latz
c
c     volume element:
               dV=dx**3
c
c        print*, ' energy computation sub-box semi-size = ', lx
c
!       write(fileTE,"(A3,I3.3)") "TE.",itime
!       write(filePE,"(A3,I3.3)") "PE.",itime
c       write(fileME,"(A3,I3.3)") "ME.",itime
c       write(fileGE,"(A3,I3.3)") "GE.",itime
c       write(fileEE,"(A3,I3.3)") "EE.",itime
c       write(fileKE,"(A3,I3.3)") "KE.",itime
!       open(unit=10,file="data/"//trim(adjustl(fileTE)),
!    1          form="unformatted",access="stream")
!       open(unit=11,file="data/"//trim(adjustl(filePE)),
!    1          form="unformatted",access="stream")
c       open(unit=12,file="data/"//trim(adjustl(fileME)),
c    1          form="unformatted",access="stream")
c       open(unit=13,file="data/"//trim(adjustl(fileGE)),
c    1          form="unformatted",access="stream")
c       open(unit=14,file="data/"//trim(adjustl(fileEE)),
c    1          form="unformatted",access="stream")
c       open(unit=15,file="data/"//trim(adjustl(fileKE)),
c    1          form="unformatted",access="stream")
!        if(itime.le.nte) then 

!       write(fileTEE,"(A4,I2.2,A1,I4.4)") "TEE.",pid,".",itime
!       write(filePE,"(A3,I2.2,A1,I3.3)") "PE.",pid,".",itime
c       write(fileME,"(A3,I3.3)") "ME.",itime
c       write(fileGE,"(A3,I3.3)") "GE.",itime
c       write(fileEE,"(A3,I3.3)") "EE.",itime
c       write(fileKE,"(A3,I3.3)") "KE.",itime
!       open(unit=10,file="data/"//trim(adjustl(fileTEE)),
!    1          form="unformatted",access="stream")
!       open(unit=11,file="data/"//trim(adjustl(filePE)),
!    1          form="unformatted",access="stream")
c       open(unit=12,file="data/"//trim(adjustl(fileME)),
c    1          form="unformatted",access="stream")
c       open(unit=13,file="data/"//trim(adjustl(fileGE)),
c    1          form="unformatted",access="stream")
c       open(unit=14,file="data/"//trim(adjustl(fileEE)),
c    1          form="unformatted",access="stream")
c       open(unit=15,file="data/"//trim(adjustl(fileKE)),
c    1          form="unformatted",access="stream")
!       endif

        do k=endb,enub
          do j=enlb,enrb
            do i=enbb,enfb
c     c     if one wants to only count the energy within a small region:
               ilocal = i+lxb-1
               jlocal = j+lyl-1
               klocal = k+lzd-1

               x=float(ilocal)*dx
               y=float(jlocal)*dx
               z=float(klocal)*dx

      call ederivatives(f,fd,i,j,k,dfdx,dfdy,dfdz,dfddx,dfddy,dfddz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

c covariant derivative: cd(mu,a)=(D_\mu\Phi)^a (a=1 is the real 
c part of the upper component of the doublet; a=2 is the imag
c part; a=3 is real part of the lower component; a=4 is the
c imag part of the lower component.

      call ecovderiv(f,fd,i,j,k,dfdx,dfdy,dfdz,cd)

c find gauge field strengths (uses derivatives):

      call efieldstrength(f,fd,i,j,k,dfdx,dfdy,dfdz,fs)

c
      kineticEnergyPhi= 0.5*(cd(0,1)**2+cd(0,2)**2+cd(0,3)**2)
c
      gradientEnergyPhi= 0.5*(cd(1,1)**2+cd(2,1)**2+cd(3,1)**2
     1                       +cd(1,2)**2+cd(2,2)**2+cd(3,2)**2
     2                       +cd(1,3)**2+cd(2,3)**2+cd(3,3)**2)

      electricEnergyW=0.5*(fs(1,0,1)**2+fs(1,0,2)**2+fs(1,0,3)**2
     1                    +fs(2,0,1)**2+fs(2,0,2)**2+fs(2,0,3)**2
     1                    +fs(3,0,1)**2+fs(3,0,2)**2+fs(3,0,3)**2)
c
      magneticEnergyW=0.5*(fs(1,1,2)**2+fs(1,1,3)**2+fs(1,2,3)**2
     1                    +fs(2,1,2)**2+fs(2,1,3)**2+fs(2,2,3)**2
     2                    +fs(3,1,2)**2+fs(3,1,3)**2+fs(3,2,3)**2)
c
      potentialEnergy=0.25*lambda*(f(1,i,j,k)**2+f(2,i,j,k)**2
     1                                    +f(3,i,j,k)**2-vev**2)**2
c
c ==============
c     
c       energyDensity=kineticEnergyPhi+gradientEnergyPhi+electricEnergyW+
c     1   magneticEnergyW+electricEnergyY+magneticEnergyY+potentialEnergy
       energyDensity=kineticEnergyPhi+gradientEnergyPhi+electricEnergyW+
     1   magneticEnergyW+potentialEnergy
c     
        totalKEPhi=totalKEPhi+kineticEnergyPhi*dV
        totalGEPhi=totalGEPhi+gradientEnergyPhi*dV
        totalEEW=totalEEW+electricEnergyW*dV
        totalMEW=totalMEW+magneticEnergyW*dV
c        totalEEY=totalEEY+electricEnergyY*dV
c        totalMEY=totalMEY+magneticEnergyY*dV
        totalPE=totalPE+potentialEnergy*dV
        totalEnergy=totalEnergy+energyDensity*dV

!         if(ilocal.eq.-latx.and.jlocal.eq.0.and.klocal.eq.0) then
!!!       print *, "latx"
!!        print *, "cd   ",cd(1,1),cd(2,1),cd(3,1)
!!        print *, "cd   ",cd(1,2),cd(2,2),cd(3,2)
!!        print *, "cd   ",cd(1,3),cd(2,3),cd(3,3)
!!        print *, f(1,i,j,k),f(2,i,j,k),f(3,i,j,k) 
!         print*,ilocal,jlocal,klocal,pid
!        print *, "cd   ",cd(0,1),cd(0,2),cd(0,3), cd(1,1),fd(1,i,j,k)
!!!       print*, dfdx(1),dfdx(2),dfdx(3)
!!!       print*,f(1,i,j,k)
!!!       print*,f(1,i+1,j,k),f(1,i+2,j,k),f(1,i+3,j,k)
!!!       print*,f(1,i-1,j,k),f(1,i-2,j,k),f(1,i-3,j,k)
!!!       print*,cd(1,1),dfdx(1),f(9,i,j,k),f(3,i,j,k),f(13,i,j,k),
!!!    1                      f(2,i,j,k),dfdx(2),dfdx(3)
!!!       print *, gradientEnergyPhi
!         print *, kineticEnergyPhi
!         endif

c
!        if(itime.le.nt) then
!         write(10) energyDensity
!         write(11) potentialEnergy
c         write(12) magneticEnergyW
c         write(13) gradientEnergyPhi
c         write(14) electricEnergyW
c         write(15) kineticEnergyPhi
!        endif
c
c
c
!        if(itime.le.nte) then
!         write(10) energyDensity
!         write(11) potentialEnergy
c         write(12) magneticEnergyW
c         write(13) gradientEnergyPhi
c         write(14) electricEnergyW
c         write(15) kineticEnergyPhi
!        endif

        enddo
        enddo
        enddo

        energy(1)=itime
        energy(2)=totalKEPhi
        energy(3)=totalGEPhi
        energy(4)=totalEEW
        energy(5)=totalMEW
        energy(6)=totalPE
        energy(7)=totalEnergy
         if(itime.le.nte) then

!       close(10)
!       close(11)
!       close(10)
!       close(11)
!       close(12)
!       close(13)
!       close(14)
!       close(15)
        endif

      end
