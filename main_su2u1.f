!=======================================================================
! Relax monopole-antimonopole field configuration using Gauss-Seidel 
! method. And then evolve using time Crank-Nicolson
!=======================================================================

! begin{declarations}

      use mpi
      implicit none

!=======================================================================
! begin{input}
! model-dependent parameters:

      include 'parameters.inc'

! problem-dependent parameters:

      include 'initialParameters.inc'

! end{input}
! ======================================================================

!     Declare all fields and time derivatives -- nf is the number of
!     fields; f denotes field. 

      real*8 f, gauge, mf
      real*8 energyinitial, totalEnergy
      real*8 energy(12), energy0(12)
      real*8 energymin
!      real*8 hatn,dxhatn,dyhatn,dzhatn
!      dimension hatn(3,1:localsize,1:localsize,1:localsize)
!      dimension dxhatn(3,1:localsize,1:localsize,1:localsize)
!      dimension dyhatn(3,1:localsize,1:localsize,1:localsize)
!      dimension dzhatn(3,1:localsize,1:localsize,1:localsize)
      dimension f(nf,1:localsize,1:localsize,1:localsize)

      dimension gauge(3,1:localsize,1:localsize,1:localsize)
      dimension mf(3,1:localsize,1:localsize,1:localsize)

      integer it, itime,i,j,k,l,ierr,np,pid,interval
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub
      integer proci,procj,prock

! end{declarations}

        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world, pid, ierr)
        call mpi_comm_size(mpi_comm_world, np, ierr)

!        call processcoord(pid,proci,procj,prock)
!        print *, pid,proci,procj,prock

        call boundindices(pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)
!        print*,pid,lxb,lyl,lzd,bb,fb,lb,rb

        if(pid==0) then
        print *, ' no. of fields =                       ', nf
        print *, ' lattice size =                        ', latx
        print *, ' lattice spacing =                     ', dx
        print *, ' Gauss-Seidel relaxparam =             ', relaxparam
        print *, ' number of snap shots =                ', nsnaps
        print *, ' gauge function parameter =            ', gp2
        print *, ' theory parameters: gw,gy,lambda,vev = ', gw, gy,
     1                  lambda,vev
        print *, 'number of processors = ', np
        print *, ' '
        endif

!=======================================================================
! begin{initial conditions}
!=======================================================================

        call initialconditions(f,hatn,dxhatn,dyhatn,dzhatn,pid,
     1                          lxb,lyl,lzd,bb,fb,lb,rb,db,ub,np)

!         call ghostcellsx(f,
!     1          pid,np,bb,fb,lb,rb,db,ub)
! 
!         call ghostcellsy(f,
!     1          pid,np,bb,fb,lb,rb,db,ub)
! 
!         call ghostcellsz(f,
!     1          pid,np,bb,fb,lb,rb,db,ub)

!=======================================================================
! end{initial conditions}
!=======================================================================

!=======================================================================
! begin{write out initial energies}
! 'it' in subroutine energy is the time step (zero right now):
!=======================================================================

        it=0
        call uenergy(f,energy,it,pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

!=======================================================================
! remember initial energy in box:
!=======================================================================

        call printenergy(energy,pid,np)
        energyinitial=energy(9)

!=======================================================================
! initalize energymin (used to stop after the energy has become minimum)
!=======================================================================

         energymin=energyinitial

!=======================================================================
! begin{evolution}
!=======================================================================

        do itime=1,nt

!=======================================================================
!Subrountines to call relaxation subroutine as well as as mpi send and
!receive relevant ghost cell data
!=======================================================================

        call evolveeuler(f,hatn,dxhatn,dyhatn,dzhatn,pid,np,
     1                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub,itime)

!         call ghostcellsx(f,
!     1          pid,np,bb,fb,lb,rb,db,ub)
! 
!         call ghostcellsy(f,
!     1          pid,np,bb,fb,lb,rb,db,ub)
! 
!         call ghostcellsz(f,
!     1          pid,np,bb,fb,lb,rb,db,ub)

!Compute and print energies at specified intervals
        interval=int(nt/nsnaps)

        if(int(nt/nsnaps).eq.0) interval=1
        if(mod(itime,interval).eq.0.or.itime.eq.nt) then
        call uenergy(f,energy,itime,pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

        call printenergy(energy,pid,np)

        totalEnergy=energy(9)

!=======================================================================
! Stop if energy minimum has been reached. Once the energy
! starts growing, it is seen to grow rapidly, signaling an
! instability. So best to stop if energy grows:
!=======================================================================

        if(totalEnergy.le.energymin) then

                energymin=totalEnergy
        else

        if(pid==0) then
                print*, ' totalEnergy.gt.energymin =  STOP '
        endif
                exit
        endif

        endif

        enddo

        call mpi_barrier(mpi_comm_world,ierr)

!        call emain(f,pid,np,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

        call mpi_finalize(ierr)

      end
