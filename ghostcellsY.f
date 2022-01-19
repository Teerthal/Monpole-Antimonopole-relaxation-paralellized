        subroutine ghostcellsy(f,
     1          pid,np,bb,fb,lb,rb,db,ub)

!=======================================================================
!              Subroutine to exchange ghostpoints
!=======================================================================

                use mpi
                implicit none
                include 'parameters.inc'
                
        integer ghostsize,ierr,np,status(mpi_status_size)
        integer bb,fb,lb,rb,db,ub
        integer pid,pidl,pidr,proci,procj,prock
        integer msgl,msgr
        integer bx,fx,dz,uz
        integer lysl,lysu,lyrl,lyru
        integer rysl,rysu,ryrl,ryru
        integer pidrbound, pidlbound

        real*8 f
        dimension f(nf,1:localsize,1:localsize,1:localsize)

        call processcoord(pid,proci,procj,prock)

        pidl=pid-nprocessx
        pidr=pid+nprocessx
        pidlbound=pid-nprocessx*(nprocessx-1)
        pidrbound=pid+nprocessx*(nprocessx-1)

        ghostsize=nf*6*localsize*localsize
       
        msgl=1
        msgr=2
!--------------------------Exchange = 1--------------------------------- 
!l & r msg correspond to left & right of message destination in both 
!mpi_send and mpi_recv. A message should be sent and recieved in such a 
!way that no two messages are sent or recieved at the same time .  
!-----------------------------------------------------------------------
               bx=1
               fx=localsize
               dz=1
               uz=localsize
 
               rysl=rb-6
               rysu=rb-1

               lyrl=lb-6
               lyru=lb-1

               lysl=lb+1
               lysu=lb+6

               ryrl=rb+1
               ryru=rb+6


!        if(procj/=nprocessx-1) then
!                call mpi_send(fnew(1:nf-1,bx:fx,rysl:rysu,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidr,
!     2                  msgl,mpi_comm_world,ierr)
!        endif
!        if(procj/=0) then
!                call mpi_recv(fnew(1:nf-1,bx:fx,lyrl:lyru,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidl,
!     2                  msgl,mpi_comm_world,status,ierr)
!                call mpi_send(fnew(1:nf-1,bx:fx,lysl:lysu,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidl,
!     2                  msgr,mpi_comm_world,ierr)
!        endif
!        if(procj/=nprocessx-1) then
!                call mpi_recv(fnew(1:nf-1,bx:fx,ryrl:ryru,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidr,
!     2                  msgr,mpi_comm_world,status,ierr)
!        endif
!
!        call mpi_barrier(mpi_comm_world,ierr)

!---------------------periodic boundary conditions----------------------

        if(procj==0) then
                call mpi_send(f(1:nf,bx:fx,lysl:lysu,dz:uz),
     1                  ghostsize,mpi_double_precision,pidrbound,
     2                  msgr,mpi_comm_world,ierr)
                call mpi_recv(f(1:nf,bx:fx,lyrl:lyru,dz:uz),
     1                  ghostsize,mpi_double_precision,pidrbound,
     2                  msgl,mpi_comm_world,status,ierr)
        endif
        if(procj==nprocessx-1) then
                call mpi_recv(f(1:nf,bx:fx,ryrl:ryru,dz:uz),
     1                  ghostsize,mpi_double_precision,pidlbound,
     2                  msgr,mpi_comm_world,status,ierr)
                call mpi_send(f(1:nf,bx:fx,rysl:rysu,dz:uz),
     1                  ghostsize,mpi_double_precision,pidlbound,
     2                  msgl,mpi_comm_world,ierr)
        endif

!----------------------------------------------------------------

        call mpi_barrier(mpi_comm_world,ierr)


!        if(procj/=nprocessx-1) then
!                call mpi_send(fdnew(1:nf-1,bx:fx,rysl:rysu,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidr,
!     2                  msgl,mpi_comm_world,ierr)
!        endif
!        if(procj/=0) then
!                call mpi_recv(fdnew(1:nf-1,bx:fx,lyrl:lyru,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidl,
!     2                  msgl,mpi_comm_world,status,ierr)
!                call mpi_send(fdnew(1:nf-1,bx:fx,lysl:lysu,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidl,
!     2                  msgr,mpi_comm_world,ierr)
!        endif
!        if(procj/=nprocessx-1) then
!                call mpi_recv(fdnew(1:nf-1,bx:fx,ryrl:ryru,dz:uz),
!     1                  ghostsize,mpi_double_precision,pidr,
!     2                  msgr,mpi_comm_world,status,ierr)
!        endif
!
!        call mpi_barrier(mpi_comm_world,ierr)

!------------------periodic boundary conditions---------------------------



        return
        end
