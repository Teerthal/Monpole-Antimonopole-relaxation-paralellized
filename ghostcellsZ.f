        subroutine ghostcellsz(f,
     1          pid,np,bb,fb,lb,rb,db,ub)

!=======================================================================
!              Subroutine to exchange ghostpoints along z-axis
!=======================================================================

                use mpi
                implicit none
                include 'parameters.inc'
                
        integer ghostsize,ierr,np,status(mpi_status_size)
        integer bb,fb,lb,rb,db,ub
        integer pid,pidd,pidu,proci,procj,prock
        integer msgd,msgu
        integer bx,fx,ly,ry
        integer dzsl,dzsu,dzrl,dzru
        integer uzsl,uzsu,uzrl,uzru
        integer piddbound, pidubound

        real*8 f
        dimension f(nf,1:localsize,1:localsize,1:localsize)

        call processcoord(pid,proci,procj,prock)

        pidd=pid-(nprocessx*nprocessx)
        pidu=pid+(nprocessx*nprocessx)
        piddbound=pid-nprocessx*nprocessx*(nprocessx-1)
        pidubound=pid+nprocessx*nprocessx*(nprocessx-1)


        ghostsize=nf*6*localsize*localsize
       
        msgd=1
        msgu=2
!--------------------------Exchange = 1--------------------------------- 
!l & r msg correspond to left & right of message destination in both 
!mpi_send and mpi_recv. A message should be sent and recieved in such a 
!way that no two messages are sent or recieved at the same time .  
!-----------------------------------------------------------------------
               ly=1
               ry=localsize
               bx=1
               fx=localsize

 
               uzsl=ub-3
               uzsu=ub-1

               dzrl=db-3
               dzru=db-1

               dzsl=db+1
               dzsu=db+3

               uzrl=ub+1
               uzru=ub+3


!        if(prock/=nprocessx-1) then
!                call mpi_send(fnew(1:nf-1,bx:fx,ly:ry,uzsl:uzsu),
!     1                  ghostsize,mpi_double_precision,pidu,
!     2                  msgd,mpi_comm_world,ierr)
!        endif
!        if(prock/=0) then
!                call mpi_recv(fnew(1:nf-1,bx:fx,ly:ry,dzrl:dzru),
!     1                  ghostsize,mpi_double_precision,pidd,
!     2                  msgd,mpi_comm_world,status,ierr)
!                call mpi_send(fnew(1:nf-1,bx:fx,ly:ry,dzsl:dzsu),
!     1                  ghostsize,mpi_double_precision,pidd,
!     2                  msgu,mpi_comm_world,ierr)
!        endif
!        if(prock/=nprocessx-1) then
!                call mpi_recv(fnew(1:nf-1,bx:fx,ly:ry,uzrl:uzru),
!     1                  ghostsize,mpi_double_precision,pidu,
!     2                  msgu,mpi_comm_world,status,ierr)
!        endif
!
!        call mpi_barrier(mpi_comm_world,ierr)

!-----------------periodic boundary conditions-----------------

        if(prock==0) then
                call mpi_send(f(1:nf,bx:fx,ly:ry,dzsl:dzsu),
     1                  ghostsize,mpi_double_precision,pidubound,
     2                  msgu,mpi_comm_world,ierr)
                call mpi_recv(f(1:nf,bx:fx,ly:ry,dzrl:dzru),
     1                  ghostsize,mpi_double_precision,pidubound,
     2                  msgd,mpi_comm_world,status,ierr)
        endif
        if(prock==nprocessx-1) then
                call mpi_recv(f(1:nf,bx:fx,ly:ry,uzrl:uzru),
     1                  ghostsize,mpi_double_precision,piddbound,
     2                  msgu,mpi_comm_world,status,ierr)
                call mpi_send(f(1:nf,bx:fx,ly:ry,uzsl:uzsu),
     1                  ghostsize,mpi_double_precision,piddbound,
     2                  msgd,mpi_comm_world,ierr)
        endif
!---------------------------------------------------------------

        call mpi_barrier(mpi_comm_world,ierr)

!        if(prock/=nprocessx-1) then
!                call mpi_send(fdnew(1:nf-1,bx:fx,ly:ry,uzsl:uzsu),
!     1                  ghostsize,mpi_double_precision,pidu,
!     2                  msgd,mpi_comm_world,ierr)
!        endif
!        if(prock/=0) then
!                call mpi_recv(fdnew(1:nf-1,bx:fx,ly:ry,dzrl:dzru),
!     1                  ghostsize,mpi_double_precision,pidd,
!     2                  msgd,mpi_comm_world,status,ierr)
!                call mpi_send(fdnew(1:nf-1,bx:fx,ly:ry,dzsl:dzsu),
!     1                  ghostsize,mpi_double_precision,pidd,
!     2                  msgu,mpi_comm_world,ierr)
!        endif
!        if(prock/=nprocessx-1) then
!                call mpi_recv(fdnew(1:nf-1,bx:fx,ly:ry,uzrl:uzru),
!     1                  ghostsize,mpi_double_precision,pidu,
!     2                  msgu,mpi_comm_world,status,ierr)
!        endif
!
!        call mpi_barrier(mpi_comm_world,ierr)

!--------------------periodic boundary conditions----------------



        return
        end
