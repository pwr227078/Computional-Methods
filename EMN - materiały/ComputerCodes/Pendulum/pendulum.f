      Program Pendulum      ! (P.Scharoch 15.09.2010)
c
      implicit none
c
      real*8 l,g, m
      real*8 h,alfa0,omega0,t0,tE, E0,frEtop
      real*8 t(0:100000),alfa(0:100000),omega(0:100000),Etot(0:100000)
c
      integer Np, n
c
      external dfalfa, dfomega
c
      parameter(g=10.00d0)
c
      common /l/l
c
      open(10,file='pendulum.inp',status='unknown')
      open(11,file='pendulum.out',status='unknown')
c
      read(10,*)frEtop,tE,Np
c
      l=1.0d0
      m=1.0d0
      E0=frEtop*m*(2.0d0*l)*g
      t0=0.0d0
      alfa0=0.0d0
      omega0=dsqrt(2.0d0*E0/m)/l
c
      h=(tE-t0)/Np
c
      t(0)=t0
      do n=1,Np
	  t(n)=t0+n*h
      enddo
c
!------------------------Runge-Kutta 5----------------------------------
c
      alfa(0)=alfa0
      omega(0)=omega0
      do n=0,Np-1
        call RK5D2
     & (t(n),h,alfa(n),omega(n),dfalfa,dfomega,alfa(n+1),omega(n+1))
      enddo
c
      do n=0,Np
        Etot(n)=0.5d0*m*l*l*omega(n)**2+m*g*l*(1.d0-dcos(alfa(n)))
      enddo
c
!-----------------------------------------------------------------------
c
      do n=0,Np
       write(11,'(4F12.6)')
     & t(n),alfa(n),omega(n),Etot(n)
      enddo
c
      stop
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
      k12=h*fun1(x+h/2.d0,y1+k11/2.d0,y2+k21/2.d0)	 ! k2
      k22=h*fun2(x+h/2.d0,y1+k11/2.d0,y2+k21/2.d0)
c
      k13=h*fun1(x+h/2.d0,y1+k12/2.d0,y2+k22/2.d0)
      k23=h*fun2(x+h/2.d0,y1+k12/2.d0,y2+k22/2.d0)
c
      k14=h*fun1(x+h,y1+k13,y2+k23)
      k24=h*fun2(x+h,y1+k13,y2+k23)
c
      y1p1=y1+(k11+2.d0*k12+2.d0*k13+k14)/6.d0
      y2p1=y2+(k21+2.d0*k22+2.d0*k23+k24)/6.d0
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION dfomega(t,alfa,omega)
c
      implicit none
c
      real*8 t,alfa,omega,l,g
c
      parameter(g=10.0d0)
c
      common /l/l
c
      dfomega=-(g/l)*dsin(alfa)
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION dfalfa(t,alfa,omega)
c
      implicit none
c
      real*8 t,alfa,omega
c
      dfalfa=omega
c
      return
      end
c
!=======================================================================
c
