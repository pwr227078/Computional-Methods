      Program FTABLE  ! tabulating a function (P.Scharoch 15.09.2010)

      implicit none
      real*8 xl,xr,h,x
      real*8 myfunction
      integer Np,ip
c
      open(10,file='ftable.inp',status='unknown') ! input file
      open(11,file='ftable.out',status='unknown') ! output file
c
      read(10,*)xl,xr,Np
c
      h=(xr-xl)/(Np-1)    ! setting the step
      do ip=1,Np
	x=xl+(ip-1)*h
	write(11,'(2F15.9)') x, myfunction(x)
      enddo
c
      stop
      end
c
!=======================================================================
c
      real*8 FUNCTION myfunction(x)   !
c
      implicit none
      real*8 x
c
      myfunction=x*sin(x)		! define here your own function
c
      return
      end
c
!=======================================================================
