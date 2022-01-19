      subroutine echangename(f,fd,fnew,fdnew)
      implicit none
      include 'parameters.inc'
      integer n,i,j,k
      real*8 f, fd, fnew, fdnew

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension fd(nf,1:localsize,1:localsize,1:localsize)
      dimension fnew(nf,1:localsize,1:localsize,1:localsize)
      dimension fdnew(nf,1:localsize,1:localsize,1:localsize)
c     
      do k=1,localsize
         do j=1,localsize
            do i=1,localsize

               do n=1,nf-1
                  f(n,i,j,k)=fnew(n,i,j,k)
                  fd(n,i,j,k)=fdnew(n,i,j,k)
               enddo

            enddo
         enddo
      enddo

      return
      end
