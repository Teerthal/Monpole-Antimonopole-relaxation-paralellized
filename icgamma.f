        subroutine icgamma(f,bx,fx,ly,ry,dz,uz,lxb,lyl,lzd)

        implicit none
        include 'parameters.inc'
        
        integer i,j,k,bx,fx,ly,ry,dz,uz,lxb,lyl,lzd

        real*8 f
        real*8 dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2

        dimension f(nf,1:localsize,1:localsize,1:localsize)
        dimension dfdx(nf),dfdy(nf),dfdz(nf)
        dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)

        do k=dz,uz
            do j=ly,ry
                do i=bx,fx

        call derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

                f(16,i,j,k)=dfdx(5)+dfdy(6)+dfdz(7)
                f(17,i,j,k)=dfdx(9)+dfdy(10)+dfdz(11)
                f(18,i,j,k)=dfdx(13)+dfdy(14)+dfdz(15)

                enddo
            enddo
        enddo


        return
        end
