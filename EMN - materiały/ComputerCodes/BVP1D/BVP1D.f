      Program BoundaryValueProblem1D   !(P.Scharoch 15.09.2010)
!-----------gravitational field inside a star---------------------------
      implicit none
c
      real*8 h,xE
      real*8 x(0:100000), y(0:100000)
c
      integer Np,n
c
      external S, k2
c
      open(10,file='BVP1D.inp',status='unknown')
      open(11,file='BVP1D.out',status='unknown')
c
      read(10,*)y(0),x(0),xE,Np
c
      h=(xE-x(0))/Np
c
      do n=1,Np
       x(n)=x(0)+n*h
       y(n)=0.0
      enddo
c
c
!------------------------Numerow-Cawling--------------------------------
c
      y(1)=0.5*h
c
	do n=1,Np-1
        call Numerow(x(n),h,y(n),y(n-1),y(n+1),k2,S)
 	enddo
c
!-----------------------------------------------------------------------

      do n=0,Np
      write(11,'(3F12.6)')x(n),y(n),
     & (1.-.5*(x(n)+2.)*dexp(-x(n)))
	enddo
c
      close(10)
      close(11)
c
      stop
      end
c
!=======================================================================
c
      SUBROUTINE Numerow(x,h,y,ym1,yp1,k2,S)
c
      implicit none
c
      real*8 x,h,y,ym1,k2,S,yp1
      real*8 hh,Ap1,Am1,A,SS
c
      hh=h**2/12.
	Ap1=1.0+hh*k2(x+h)
	Am1=1.0+hh*k2(x-h)
	A= 2.0*(1.0-5.0*hh*k2(x))
	SS=hh*(S(x+h)+10.0*S(x)+S(x-h))
c
      if(yp1.eq.0.0)then
	yp1=A/Ap1*y-Am1/Ap1*ym1+SS/Ap1
	endif
c
      if(ym1.eq.0.0)then
	ym1=A/Am1*y-Ap1/Am1*yp1+SS/Am1
	endif
c
      return
	end
c
!=======================================================================
c
      real*8 FUNCTION S(x)
c
      implicit none
c
      real*8 x
c
      S=-x*dexp(-x)/2.0
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION k2(x)
c
      implicit none
c
      real*8 x
c
      k2=0.0*x
c
      return
      end
c
!=======================================================================
