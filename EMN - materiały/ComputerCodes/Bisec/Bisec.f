      Program BISEC  ! finding root    (P.Scharoch 15.09.2010)
c
      implicit none
      real*8 xl,xr,eps
      real*8 func,zero
c
      external func
c
      open(10,file='bisec.inp',status='unknown')
      open(11,file='bisec.out',status='unknown')
c
      read(10,*)xl,xr,eps
c
      write(11,'(2F20.15)')ZERO(xl,xr,eps,func), eps
c
      stop
      end
c
!=======================================================================
c
      real*8 function ZERO(xl,xr,eps,fun)  ! zero must be within (xl,xr)
c
      implicit none
      real*8 xl,xr,xm,eps
      real*8 fun
c
100   xm=(xl+xr)/2.d0
      if((xr-xl).lt.eps) goto 200
c
      if(fun(xl)*fun(xm).lt.0.0d0)then
	xr=xm
        goto 100
      else
        xl=xm
        goto 100
      endif
c
200   ZERO=xm
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION FUNC(x)
c
      implicit none
      real*8 x
c
      func=sin(x)	!
c
      return
      end
c
!=======================================================================
