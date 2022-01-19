      subroutine eleapforward(f,fd,fnew,fdnew,pid,np,
     1                      lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

      implicit none
      include 'parameters.inc'
      integer n,i,j,k,pid,np
      integer bb,fb,lb,rb,db,ub,lxb,lyl,lzd

      real*8 f,fd,fnew,fdnew,ftemp,fdtemp,r,s

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension fd(nf,1:localsize,1:localsize,1:localsize)
      dimension fnew(nf,1:localsize,1:localsize,1:localsize)
      dimension fdnew(nf,1:localsize,1:localsize,1:localsize)
      dimension ftemp(nf,1:localsize,1:localsize,1:localsize)
      dimension fdtemp(nf,1:localsize,1:localsize,1:localsize)
      dimension r(nf),s(nf)

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
!=======================================================================

        do k=db,ub
            do j=lb,rb
                do i=bb,fb

        call efdflux(fnew,fdnew,i,j,k,r,lxb,lyl,lzd)
        call effluxes(fnew,fdnew,i,j,k,s,lxb,lyl,lzd)

                do n=1,nf-1

              ftemp(n,i,j,k)=f(n,i,j,k)+dt*s(n)
              fdtemp(n,i,j,k)=fd(n,i,j,k)+dt*r(n)

                enddo


                enddo
            enddo
        enddo

        do k=db,ub
            do j=lb,rb
                do i=bb,fb

                do n=1,nf-1
                  fnew(n,i,j,k)=ftemp(n,i,j,k)
                  fdnew(n,i,j,k)=fdtemp(n,i,j,k)
                enddo

                enddo
            enddo
        enddo

        call eghostcellsx(fnew,fdnew,
     1          pid,np,bb,fb,lb,rb,db,ub)

        call eghostcellsy(fnew,fdnew,
     1          pid,np,bb,fb,lb,rb,db,ub)

        call eghostcellsz(fnew,fdnew,
     1          pid,np,bb,fb,lb,rb,db,ub)

      return
      end
