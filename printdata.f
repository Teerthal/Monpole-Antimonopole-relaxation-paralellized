      subroutine printdata(f, mf)


c begin{declarations}
      implicit none
c ===================
c model-dependent parameters:
      include 'parameters.inc'
c problem-dependent parameters:
      include 'initialParameters.inc'
c
        integer i, j, k, l

        real*8 f, mf 
      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
      dimension mf(3,-latx:latx,-laty:laty,-latz:latz)


       open(unit=1, file='mfield.dat',
     1                         access='append')       
       open(unit=2, file='Phi1.dat',
     1                         access='append')       
       open(unit=3, file='Phi2.dat',
     1                         access='append') 
       open(unit=4, file='Phi3.dat',
     1                         access='append') 
       open(unit=5, file='W_1^1.dat',
     1                         access='append') 
       open(unit=6, file='W_2^1.dat',
     1                         access='append') 
       open(unit=7, file='W_3^1.dat',
     1                         access='append') 
       open(unit=8, file='W_1^2.dat',
     1                         access='append') 
       open(unit=9, file='W_2^2.dat',
     1                         access='append') 
       open(unit=10, file='W_3^2.dat',
     1                         access='append') 
       open(unit=11, file='W_1^3.dat',
     1                         access='append') 
       open(unit=12, file='W_2^3.dat',
     1                         access='append') 
       open(unit=13, file='W_3^3.dat',
     1                         access='append') 
 


       do 23 i = -latx+1,latx-1
          do 22 j = -laty+1,laty-1
           do 21 k = -latz+1,latz-1

        write(1,*) mf(1,i,j,k),mf(2,i,j,k),mf(3,i,j,k)
        write(2,*) f(1,i,j,k)
        write(3,*) f(2,i,j,k)
        write(4,*) f(3,i,j,k)
        write(5,*) f(5,i,j,k)
        write(6,*) f(6,i,j,k)
        write(7,*) f(7,i,j,k)
        write(8,*) f(9,i,j,k)
        write(9,*) f(10,i,j,k)
        write(10,*) f(11,i,j,k)
        write(11,*) f(13,i,j,k)
        write(12,*) f(14,i,j,k)
        write(13,*) f(15,i,j,k)

21      continue
22      continue
23      continue 

         do 200 l=1,13
             close(l)  
200     continue
        
        return
        end
