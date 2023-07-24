      program elabora
c
      implicit none
      integer i,num5
      double precision num1,num2,num3,num4,num6,num7,num8
c
      open(unit=10,file='riferimento1.dat',status='unknown')
      open(unit=11,file='i1.dat',status='unknown')
      open(unit=12,file='fileout.dat',status='unknown')
      do i=1,1000
         read(10,*) num1,num2,num3,num4
         read(11,*) num5,num6,num7,num8
         write(12,*) num3,num6,num4
      end do
      close(10)
      close(11)
      close(12)
c
      end
