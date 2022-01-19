        subroutine eprintenergy(energy,pid,np)
                use mpi

        implicit none

        real*8 energy(7), energy0(7)
        integer i,msg,ierr,np,pid,status(mpi_status_size)

        msg=1

        if(pid/=0) then
              call mpi_send(energy(1),7,mpi_double_precision,0,msg,
     1          mpi_comm_world, ierr)
        endif

        if(pid==0) then
              do i=1,np-1

              call mpi_recv(energy0(1),7,mpi_double_precision,i,msg,
     1             mpi_comm_world,status, ierr)
              
              energy=energy+energy0
              enddo
        endif

              call mpi_barrier(mpi_comm_world,ierr)
              call mpi_bcast(energy(1),7,mpi_double_precision,0,
     1             mpi_comm_world,ierr)
              call mpi_barrier(mpi_comm_world,ierr)


                
        if(pid==0) then
        print*,' =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* '
        print*,' time:           ', energy(1)/np
        print*,' =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* '
        print*,' KE:             ', energy(2)
        print*,' GE:             ', energy(3)
        print*,' EEW:            ', energy(4)
        print*,' MEW:            ', energy(5)
        print*,' PEtotalE:       ', energy(6)
        print*,' totalE:         ', energy(7)
        print*,' '
        endif

        return
        end
