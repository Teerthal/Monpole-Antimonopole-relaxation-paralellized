        subroutine evolveloop(f,hatn,dxhatn,dyhatn,dzhatn,
     1                          bb,fb,lb,rb,db,ub,lxb,lyl,lzd)

        implicit none
        include 'parameters.inc'
        
        integer i,j,k,bb,fb,lb,rb,db,ub,lxb,lyl,lzd
        integer ilocal,jlocal,klocal
        integer n

        real*8 f,coeff,r
        real*8 hatn,dxhatn,dyhatn,dzhatn

        dimension f(nf,1:localsize,1:localsize,1:localsize)
        dimension hatn(3,1:localsize,1:localsize,1:localsize)
        dimension dxhatn(3,1:localsize,1:localsize,1:localsize)
        dimension dyhatn(3,1:localsize,1:localsize,1:localsize)
        dimension dzhatn(3,1:localsize,1:localsize,1:localsize)

        dimension r(nf)

        coeff=49./6.

! Looping over local coordinates
        do k=db,ub
            do j=lb,rb
                do i=bb,fb
! Global coordinates (suffix local)
            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

        call fdflux(f,i,j,k,r,hatn,dxhatn,dyhatn,dzhatn,lxb,lyl,lzd)

! Omitting end boundaries
        if(abs(ilocal).ne.latx.and.abs(jlocal)
     1                          .ne.laty.and.abs(klocal).ne.latz) then

! Scalar field evolutions
        f(19,i,j,k)= f(19,i,j,k)+relaxparam*dx**2*r(19)/coeff
        f(1,i,j,k)=f(19,i,j,k)*hatn(1,i,j,k)
        f(2,i,j,k)=f(19,i,j,k)*hatn(2,i,j,k)
        f(3,i,j,k)=f(19,i,j,k)*hatn(3,i,j,k)

! Gauge field evolutions
                do n=4,nf-1
                f(n,i,j,k)= f(n,i,j,k)+relaxparam*dx**2*r(n)/coeff
                enddo
        endif

                enddo
            enddo
        enddo


        return
        end
