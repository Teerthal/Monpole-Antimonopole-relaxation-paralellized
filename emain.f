c 1 August 2015
c Gauged models require subroutine that checks if Lorenz gauge is satisfied.
c
c     November 16, 2008
c     Evolve classical evolution equations for a field theory. 
c     (See also the directory  ~/programs/ab-vortices)
c
c ====================
c begin{declarations}
      subroutine emain(f,pid,np,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

        use mpi
      implicit none

c ===================
c begin{input}
c model-dependent parameters:
      include 'parameters.inc'
c problem-dependent parameters:
      include 'einitialParameters.inc'
c end{input}
c ===================
c
c     declare all fields and time derivatives -- nf is the number of
c     fields; f denotes field; fd denotes time derivative of f. We should
c     consider if it is better to have just one array whose elements
c     contain both the field and its time derivative.
c
c     A note about boundaries: the box goes from -latx+1 to latx-1 in x, y
c and z. Only the field value is calculated on the boundaries of the box 
c(x=+lat etc.), never the value of the time-derivative. This should be
c kept in mind when calculating energy etc..
c
      integer interval
      real*8 f, fd, fnew, fdnew
      real*8 energyinitial, totalEnergy
      real*8 energy(7), energy0(7)
c     real*8 totalHelicity
c     real*8 csNumber
c
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension fd(nf,1:localsize,1:localsize,1:localsize)
c
c     staggered leapfrog requires storing fields at two time steps:
      dimension fnew(nf,1:localsize,1:localsize,1:localsize)
      dimension fdnew(nf,1:localsize,1:localsize,1:localsize)
c
      integer it, itime,i,j,k,l,ierr,np,pid
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub

c
c end{declarations}
c ==================

      call mpi_barrier(mpi_comm_world,ierr)

      if(pid==0) then

      print*,'latx=',latx,' nt=',nt
      print*, ' dx=', dx,' dt= ', dt,' gp2=',gp2 
      print*, ' nsnaps= ',nsnaps
c theory parameters:
      print*, ' nf,gw,gy,lambda,vev = ',nf,gw, gy, lambda,vev
c for SO(3) gauge wave collisions (mmbar creation problem):
      print*, ' z0,width,vz0,omega,amp=',z0,width,vz0,omega,amp
      endif

c
c ======================
c begin{initial conditions}

        call einitialconditions(f,fd,fnew,fdnew,pid,
     1                          lxb,lyl,lzd,np,bb,fb,lb,rb,db,ub)
!----------periodic boundary conditions--------------------------

!       call eghostcellsx(f,fd,
!    1          pid,np,bb,fb,lb,rb,db,ub)

!       call eghostcellsy(f,fd,
!    1          pid,np,bb,fb,lb,rb,db,ub)

!       call eghostcellsz(f,fd,
!    1          pid,np,bb,fb,lb,rb,db,ub)

c end{initial conditions}
c ======================
c begin{write out initial energies}
c 'it' in subroutine energy is the time step (zero right now):

      it=0
        call eenergy(f,fd,energy,it,
     1          pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

c remember initial energy in box:

        call eprintenergy(energy,pid,np)
        energyinitial=energy(7)

        if(pid==0) then 
        print *,'totalEnergy Initial',energyinitial
        endif

c end{write out initial energies}
cc =======================
c begin{evolution}
c
      do itime=1,nte
c
c Euler method to start evolution (includes calls to fdflux and absorbingbc):
c (Note: the evolution is only done on the lattice minus the boundary layer.
c The value of the field on the boundary layer is determined by the boundary
c condition. The time derivative of the field on the boundary layer is not
c needed. The boundary value of the field is used to determine the flux on
c the inner layer (e.g. x=latx-1), where the field and time derivatives are
c evolved by the euler, etc. methods.)

        call eevolveeuler(f,fd,fnew,fdnew,pid,np,
     1                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub)
!       call eevolveeulerPBC(f,fd,fnew,fdnew,pid,np,
!    1                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

c At this point fnew, fdnew are f,fd euler evolved by dt.

        call eaverageforhalfstep(f,fd,fnew,fdnew)

! At this point fnew, fdnew are at dt/2.
! Use leapfrog -- (includes call to boundary conditions)

        call eleapforward(f,fd,fnew,fdnew,pid,np,
     1                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

!       call eleapforwardPBC(f,fd,fnew,fdnew,pid,np,
!    1                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub)


c At this point fnew, fdnew are at dt.
c Now do second iteration (average then leap):

        call eaverageforhalfstep(f,fd,fnew,fdnew)

        call eleapforward(f,fd,fnew,fdnew,pid,np,
     1                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub)
!       call eleapforwardPBC(f,fd,fnew,fdnew,pid,np,
!    1                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub)


c At this point f,fd are at initial time, fnew,fdnew at final time.
c Rename f to be fnew and fd to be fdnew to go on to the next time step:

        call echangename(f,fd,fnew,fdnew)


cc-----
cccc use if interested in soliton creation:
cccc find min and max of |phi| (small values would be signal for vortices):
c        call phiwinding(f,itime,pid)
cc-----
c take nsnaps (defined in parameters.inc) snapshots of energy distribution:
c (nsnaps shouldn't be too large otherwise memory will fill up.)
c A snapshot at the final time is included.
c        interval=int(nt/nsnaps)

         interval=1
c        if(int(nt/nsnaps).eq.0) interval=1

         if(mod(itime,interval).eq.0.or.itime.eq.nte) then

       call eenergy(f,fd,energy,itime,
     1          pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)
       call eprintenergy(energy,pid,np)

        totalEnergy=energy(7)

c so3.charge.f -- to evaluate magnetic charge -- still needs some work.
c           call charge(f,fd,itime)
c           call chernSimons(f,fd,csNumber,itime) 
c           call helicity(f,fd,totalHelicity,itime)
c           print*, ' timestep, total energy ', itime,totalEnergy,
c
c stop if most of the energy has left the box:

        if(totalEnergy.lt.energyinitial/10.) then

        if(pid==0) then
                print*, ' 90 percent energy lost; stopping run '
        endif
                exit
        endif

        endif

        enddo

      end
