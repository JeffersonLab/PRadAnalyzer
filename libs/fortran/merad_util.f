      subroutine grid_init
      implicit none
      integer i
      include 'merad_grid.inc'
      do i=1,30
        grv(i)=dble(i-1)/dble(29)/4d0
      enddo
      do i=31,45
        grv(i)=0.25d0+dble(i-30)/dble(45-30)/4d0
      enddo
      do i=46,nv
        grv(i)=0.5d0+dble(i-45)/dble(nv-45)/2d0
      enddo
      do i=1,7
        grt1(i)=0.1d0*dble(i)**2/49d0/2d0
        grt1(31-i)=1d0+grt1(1)-grt1(i)
      enddo
      do i=8,15
        grt1(i)=(0.1d0+0.9d0*(dble(i-7))/8d0)/2d0
        grt1(31-i)=1d0+grt1(1)-grt1(i)
      enddo
      do i=1,30
        grz(i)=0.5d0*dble(i**2)/30d0**2
      enddo
      do i=1,30
        grz(61-i)=1d0-0.49d0*dble(i**2)/30d0**2
      enddo
      end

      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
! simps
! a1,b1 -the limits of integration
! h1 -an initial step of integration
! reps1,aeps1 - relative and absolute precision of integration
! funct -a name of function subprogram for calculation of integrand +
! x - an argument of the integrand
! ai - the value of integral
! aih- the value of integral with the step of integration
! aiabs- the value of integral for module of the integrand
! this subrogram calculates the definite integral with the relative or
! absolute precision by simpson+s method with the automatical choice
! of the step of integration
! if aeps1    is very small(like 1.e-17),then calculation of integral
! with reps1,and if reps1 is very small (like 1.e-10),then calculation
! of integral with aeps1
! when aeps1=reps1=0. then calculation with the constant step h1
!
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpsx(a,b,np,ep,func,res)
      implicit real*8 (a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

