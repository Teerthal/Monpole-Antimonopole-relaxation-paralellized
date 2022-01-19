      subroutine evolveeuler(f,hatn,dxhatn,dyhatn,dzhatn,pid,np,
     1                      lxb,lyl,lzd,bb,fb,lb,rb,db,ub,itime)
              use mpi
      implicit none
      include 'parameters.inc'

!=======================================================================
! need initialParameters to implement the *constrained* relaxation.
! we need to hold the position of the mmbar fixed and these are given
! in the initialParameters. Better would be to impose the constraints
! when calling the evolution subroutines and not in the evolution
! routine itself.
!=======================================================================
      
      include 'initialParameters.inc'

      integer n,i,j,k,itime
      integer bb,fb,lb,rb,db,ub,lxb,lyl,lzd
      real*8 f
      real*8 error
      real*8 hatn,dxhatn,dyhatn,dzhatn
 
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension hatn(3,1:localsize,1:localsize,1:localsize)
      dimension dxhatn(3,1:localsize,1:localsize,1:localsize)
      dimension dyhatn(3,1:localsize,1:localsize,1:localsize)
      dimension dzhatn(3,1:localsize,1:localsize,1:localsize)

      dimension error(nf)
      integer np,pid,proci,procj,prock

      call processcoord(pid,proci,procj,prock)

              
!=======================================================================
! e.g. eq.(107) in http://www.damtp.cam.ac.uk/lab/people/sd/lectures/nummeth98/pdes.htm
!
!
! evolve on lattice (one-sided differencing on boundaries as
! defined in subroutine derivatives):: 
!=======================================================================

           do n=1,nf
           error(n)=0.
           enddo

!=======================================================================
! The Gauss-Seidel method is sensitive to how one loops over the
! lattice. I've tried looping from one end of the lattice to the
! other. But if I start from the center of the lattice, the code
! seems to do better (converge faster). Should explore this further.
! Some Dictionary:
! bb = back boundary...xm
! fb = front boundary...xp
! lb = left boundary...ym
! rb = right boundary...yp
! Above are local coordinates that the evolution subroutines loop over 
!=======================================================================

! In the comments below, only ghosts around the boundaries are sent and received.
! Word ghost is ommitted i.e. any reference to boundary indices is to be understood as
! the relevant ghost values. Additionally, only fb and bb corresponding to x-coordinates are 
! used for brevity. The same comments apply to y and z coordinates.

!------------------------Relaxation-------------------------------------

! For the 0th process, firt evolve, then send the local front boundary (fb) to front process (pidf).
! Then finally receive back the updated front boundary(fb) values from the pidf for next iteration.

        if(proci==0.and.procj==0.and.prock==0) then

        call evolveloop(f,hatn,dxhatn,dyhatn,dzhatn,
     1                  bb,fb,lb,rb,db,ub,lxb,lyl,lzd)

        call ghostsend(f,pid,np,bb,fb,lb,rb,db,ub,itime)

        call ghostrecvl(f,pid,np,bb,fb,lb,rb,db,ub,itime)

! For the end boundary processes, first receive back boundary (bb) from back process (pidb)
! Then evolve local sublattice, followed by sending the updated bb values back to pidb
        elseif(proci==nprocessx-1.and.procj==nprocessx-1
     1  .and.prock==nprocessx-1) then

        call ghostrecvu(f,pid,np,bb,fb,lb,rb,db,ub,itime)

        call evolveloop(f,hatn,dxhatn,dyhatn,dzhatn,
     1                  bb,fb,lb,rb,db,ub,lxb,lyl,lzd)
        call ghostsend(f,pid,np,bb,fb,lb,rb,db,ub,itime)

        else

! For non boundary processes, first receive bb from pidb.
! Then evolve local sublattice. Then send the bb and fb ghosts 
! to pidb and pidf. Finally receive updated fb from pidf for next iteration.
        call ghostrecvu(f,pid,np,bb,fb,lb,rb,db,ub,itime)

        call evolveloop(f,hatn,dxhatn,dyhatn,dzhatn,
     1                  bb,fb,lb,rb,db,ub,lxb,lyl,lzd)

        call ghostsend(f,pid,np,bb,fb,lb,rb,db,ub,itime)

        call ghostrecvl(f,pid,np,bb,fb,lb,rb,db,ub,itime)

        endif

        return
        end
