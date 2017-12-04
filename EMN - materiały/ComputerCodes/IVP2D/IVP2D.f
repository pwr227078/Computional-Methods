      Program InitialValueProplem2D    !(P.Scharoch 15.09.2010)
c
      implicit none
c
      real*8 h,y10,y20,x0,xE,x(0:1000)
      real*8 y1E(0:1000),y1AB3(0:1000),y1RK5(0:1000)
      real*8 y2E(0:1000),y2AB3(0:1000),y2RK5(0:1000)
c
      integer Np, n
c
      external drfun1, drfun2
c
      open(10,file='IVP2D.inp',status='unknown')
      open(11,file='IVP2D.out',status='unknown')
c
      read(10,*)y10,y20,x0,xE,Np
c
      h=(xE-x0)/Np
c
      x(0)=x0
      do n=1,Np
	x(n)=x0+n*h
      enddo
c
!------------------------Euler------------------------------------------
c
      y1E(0)=y10
      y2E(0)=y20
	do n=0,Np-1
      call EulerD2(x(n),h,y1E(n),y2E(n),drfun1,drfun2,y1E(n+1),y2E(n+1))
	enddo
c
!------------------------Runge-Kutta 5----------------------------------
c
      y1RK5(0)=y10
      y2RK5(0)=y20
      do n=0,Np-1
        call RK5D2
     &(x(n),h,y1RK5(n),y2RK5(n),drfun1,drfun2,y1RK5(n+1),y2RK5(n+1))
      enddo
c
!--------------------Adams-Bashforth 3-punktowy-------------------------
c
      y1AB3(0)=y10
      y1AB3(1)=y1E(1)
      y1AB3(2)=y1E(2)
      y2AB3(0)=y20
      y2AB3(1)=y2E(1)
      y2AB3(2)=y2E(2)
c
      do n=2,Np-1
	 call AB3D2
     &  (x(n),h,y1AB3(n-2),y2AB3(n-2),y1AB3(n-1),y2AB3(n-1),
     &  y1AB3(n),y2AB3(n),drfun1,drfun2,y1AB3(n+1),y2AB3(n+1))
      enddo
!-----------------------------------------------------------------------

      do n=0,Np
        write(11,'(5F12.6)')
     &  x(n),y1E(n),y1AB3(n),y1RK5(n),y10*dexp(-x(n))
      enddo
c
      stop
      end
c
!=======================================================================
c
      SUBROUTINE EulerD2(x,h,y1,y2,fun1,fun2,y1p1,y2p1)
c
      implicit none
c
      real*8 x,y1,y2,h,fun1,fun2,y1p1,y2p1
c
      y1p1=y1+h*fun1(x,y1,y2)
      y2p1=y2+h*fun2(x,y1,y2)
c
      return
      end
c
!=======================================================================
c
      SUBROUTINE 
     & AB3D2(x,h,y1m2,y2m2,y1m1,y2m1,y1,y2,fun1,fun2,y1p1,y2p1)
c
      implicit none
c
      real*8 x,h
      real*8 y1m2,y1m1,y1,fun1,y1p1
      real*8 y2m2,y2m1,y2,fun2,y2p1
c
      y1p1=y1+ h/12.*(5.*fun1(x-2.*h,y1m2,y2m2)
     &              -16.*fun1(x-h,y1m1,y2m1)+23.*fun1(x,y1,y2))
      y2p1=y2+ h/12.*(5.*fun2(x-2.*h,y1m2,y2m2)
     &              -16.*fun2(x-h,y1m1,y2m1)+23.*fun2(x,y1,y2))

c
      return
      end
c
!=======================================================================
c
      SUBROUTINE RK5D2(x,h,y1,y2,fun1,fun2,y1p1,y2p1)
c
      implicit none
c
      real*8 x,h,y1,y2,fun1,fun2,y1p1,y2p1
      real*8 k11,k12,k13,k14
      real*8 k21,k22,k23,k24
c
      k11=h*fun1(x,y1,y2)	 ! k1
      k21=h*fun2(x,y1,y2)
c
      k12=h*fun1(x+h/2.,y1+k11/2.,y2+k21/2.)	 ! k2
      k22=h*fun2(x+h/2.,y1+k11/2.,y2+k21/2.)
c
      k13=h*fun1(x+h/2.,y1+k12/2.,y2+k22/2.)
      k23=h*fun2(x+h/2.,y1+k12/2.,y2+k22/2.)
c
      k14=h*fun1(x+h,y1+k13,y2+k23)
      k24=h*fun2(x+h,y1+k13,y2+k23)
c
      y1p1=y1+(k11+2.*k12+2.*k13+k14)/6.
      y2p1=y2+(k21+2.*k22+2.*k23+k24)/6.
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfun1(x,y1,y2)
c
      implicit none
c
      real*8 x,y1,y2
c
      drfun1=y2
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION drfun2(x,y1,y2)
c
      implicit none
c
      real*8 x,y1,y2
c
      drfun2=y1
c
      return
      end
c
!=======================================================================
c
