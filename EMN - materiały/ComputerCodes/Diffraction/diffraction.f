      Program diffraction   !(P.Scharoch, 15.09.2010r.)
c
      implicit none
c
      real*8 yl,yr,a,d,lambda,k,Ampl
      real*8 ReCompA,ImCompA,simpson
      real*8 hy,y
      real*8 I
      real*8 Pi
      integer Ny,Nint,iy
c
      external ReCompA,ImCompA
c
      parameter (Nint=100, Pi=3.141592654d0)
c
      common /dyfrpar/Ampl,k,d,y
c
      open(10,file='diffraction.inp',status='unknown')
      open(11,file='diffraction.out',status='unknown')
c
      read(10,*)Ampl,a,d,yl,yr,Ny ! a,d,yl,yr are in units of wavelength
      lambda=1.0d0
c
      hy=(yr-yl)/(Ny-1)
      k=2.d0*Pi/lambda
c
      do iy=1,Ny
      y=yl+(iy-1)*hy
      I=simpson(-a/2.0d0,a/2.0d0,Nint,ReCompA)**2 +
     &	simpson(-a/2.0d0,a/2.0d0,Nint,ImCompA)**2
c
      write(11,'(2F15.6)')y,I
      enddo
c
      close(10)
      close(11)
      stop
      end
c
!=======================================================================
c
      real*8 FUNCTION simpson(xl,xr,Nh,fun)
c
      implicit none
c
      real*8 x,h,fun,xl,xr
      integer ih, Nh
c
      h=(xr-xl)/Nh/2.d0
c
      simpson=0.0d0
	x=xl+h
c
      do ih=1,Nh
        simpson=simpson + h*(fun(x+h)+4.0d0*fun(x)+fun(x-h))/3.d0
	x=x+2.0d0*h
      enddo
c
      return
	end
c
!=======================================================================
c
      real*8 FUNCTION ReCompA(x) ! real part ofthe complex amplitude
c
      implicit none
c
      real*8 x,k,r,d,y,Ampl
c
      common /dyfrpar/Ampl,k,d,y
c
      r=dsqrt(d*d+(y-x)**2)
      ReCompA=Ampl*dcos(k*r)/dsqrt(r)		
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION ImCompA(x) !imaginary part ofthe complex amplitude
c
      implicit none
c
      real*8 x,k,r,d,y,Ampl
c
      common /dyfrpar/Ampl,k,d,y
c
      r=dsqrt(d*d+(y-x)**2)
      ImCompA=Ampl*dsin(k*r)/dsqrt(r)		
c
      return
      end
c
!=======================================================================
