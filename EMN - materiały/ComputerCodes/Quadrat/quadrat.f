      Program quadrature   !P.Scharoch; 15.09.2010r
c
      implicit none
c
      real*8 xl,xr,ifact
      real*8 func, simpson
c
      external func
      integer Nh1,Nh,ih,Nstep
c
      open(10,file='quadrat.inp',status='unknown')
      open(11,file='quadrat.out',status='unknown')
c
      read(10,*)xl,xr,Nh1,ifact,Nstep
c
      Nh=Nh1
      do ih=1,Nstep
      	 write(11,'(i4,F15.6)')Nh,simpson(xl,xr,Nh,func)
      Nh=Nh*ifact
      enddo
c
      stop
      end
c
c====================================================================
c
      real*8 FUNCTION simpson(xl,xr,Nh,fun)
c
      implicit none
c
      real*8 x,h,fun,xl,xr
      integer ih,Nh
c
      h=(xr-xl)/Nh/2.d0
c
      simpson=0.d0
      x=xl+h
c
      do ih=1,Nh
        simpson=simpson + h*(fun(x+h)+4.d0*fun(x)+fun(x-h))/3.d0
	x=x+2.d0*h
      enddo
c
      return
      end
c
c====================================================================
c
      real*8 FUNCTION func(x)
c
      implicit none
c
      real*8 x
c
      func=sin(x)*sin(x)		! tu wstaw swoj¹ funkcjê
c
      return
	end
c
c====================================================================		 
