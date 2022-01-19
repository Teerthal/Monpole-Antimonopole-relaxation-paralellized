        subroutine eghostcellsx(fnew,fdnew,
     1          pid,np,bb,fb,lb,rb,db,ub)

!=======================================================================
!              Subroutine to exchange ghostpoints
!=======================================================================

                use mpi
                implicit none
                include 'parameters.inc'
                
        integer ghostsize,ierr,np,status(mpi_status_size)
        integer bb,fb,lb,rb,db,ub
        integer pid,pidb,pidf,proci,procj,prock
        integer msgb,msgf
        integer ly,ry,dz,uz
        integer bxsl,bxsu,bxrl,bxru
        integer fxsl,fxsu,fxrl,fxru
        integer pidfbound, pidbbound


        real*8 fnew,fdnew
        dimension fnew(nf,1:localsize,1:localsize,1:localsize)
        dimension fdnew(nf,1:localsize,1:localsize,1:localsize)

        call processcoord(pid,proci,procj,prock)

        pidb=pid-1
        pidf=pid+1
        pidbbound=pid-(nprocessx-1)
        pidfbound=pid+(nprocessx-1)


        ghostsize=(nf-1)*localsize*localsize*3
       
        msgb=1
        msgf=2
!--------------------------Exchange = 1--------------------------------- 
!l & r msg correspond to left & right of the processor executing
!mpi_send and mpi_recv. A message should be sent and recieved in such a 
!way that no two messages are sent or recieved at the same time.  
!-----------------------------------------------------------------------
               ly=1
               ry=localsize
               dz=1
               uz=localsize
 
               fxsl=fb-3
               fxsu=fb-1

               bxrl=bb-3
               bxru=bb-1

               bxsl=bb+1
               bxsu=bb+3

               fxrl=fb+1
               fxru=fb+3


        if(proci/=nprocessx-1) then
                call mpi_send(fnew(1:nf-1,fxsl:fxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidf,
     2                  msgb,mpi_comm_world,ierr)
        endif
        if(proci/=0) then
                call mpi_recv(fnew(1:nf-1,bxrl:bxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidb,
     2                  msgb,mpi_comm_world,status,ierr)
                call mpi_send(fnew(1:nf-1,bxsl:bxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidb,
     2                  msgf,mpi_comm_world,ierr)
        endif
        if(proci/=nprocessx-1) then
                call mpi_recv(fnew(1:nf-1,fxrl:fxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidf,
     2                  msgf,mpi_comm_world,status,ierr)
        endif

        call mpi_barrier(mpi_comm_world,ierr)

!-------------periodic boundary conditions----------------------

        if(proci==0) then
                call mpi_send(fnew(1:nf,bxsl:bxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidfbound,
     2                  msgf,mpi_comm_world,ierr)
                call mpi_recv(fnew(1:nf,bxrl:bxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidfbound,
     2                  msgb,mpi_comm_world,status,ierr)
        endif
        if(proci==nprocessx-1) then
                call mpi_recv(fnew(1:nf,fxrl:fxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidbbound,
     2                  msgf,mpi_comm_world,status,ierr)
                call mpi_send(fnew(1:nf,fxsl:fxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidbbound,
     2                  msgb,mpi_comm_world,ierr)
        endif
!---------------------------------------------------------------

        call mpi_barrier(mpi_comm_world,ierr)


        if(proci/=nprocessx-1) then
                call mpi_send(fdnew(1:nf-1,fxsl:fxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidf,
     2                  msgb,mpi_comm_world,ierr)
        endif
        if(proci/=0) then
                call mpi_recv(fdnew(1:nf-1,bxrl:bxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidb,
     2                  msgb,mpi_comm_world,status,ierr)
                call mpi_send(fdnew(1:nf-1,bxsl:bxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidb,
     2                  msgf,mpi_comm_world,ierr)
        endif
        if(proci/=nprocessx-1) then
                call mpi_recv(fdnew(1:nf-1,fxrl:fxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidf,
     2                  msgf,mpi_comm_world,status,ierr)
        endif

        call mpi_barrier(mpi_comm_world,ierr)

!------------------periodic conditions--------------------------

        if(proci==0) then
                call mpi_send(fdnew(1:nf,bxsl:bxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidfbound,
     2                  msgf,mpi_comm_world,ierr)
                call mpi_recv(fdnew(1:nf,bxrl:bxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidfbound,
     2                  msgb,mpi_comm_world,status,ierr)
        endif
        if(proci==nprocessx-1) then
                call mpi_recv(fdnew(1:nf,fxrl:fxru,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidbbound,
     2                  msgf,mpi_comm_world,status,ierr)
                call mpi_send(fdnew(1:nf,fxsl:fxsu,ly:ry,dz:uz),
     1                  ghostsize,mpi_double_precision,pidbbound,
     2                  msgb,mpi_comm_world,ierr)
        endif

        call mpi_barrier(mpi_comm_world,ierr)

!--------------------------------------------------------------

        return
        end
