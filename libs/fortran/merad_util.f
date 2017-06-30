! Subroutines copied from MERADGEN 1.0
! Reference: A. Afanasev, E. Chudakov, A. Ilyichev, V. Zykunov,
!            Comput. Phys. Commun. 176, 218 (2007)
! Adjusted grid_init, using nv, nt1, nz instead of numbers, Chao Peng, 05/13/2017

!===============================================================================
      subroutine merad_init(si)
! Initialize MERADGEN by Mandelstam variable s in MeV^2
     &bind(C, name = "merad_init")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: si
      include 'merad_const.inc'
      data pi/3.14159265359d0/, alfa/.7297352568d-2/,
     .     m/0.510998918d0/

      m2=m*m
      s=si
      als=s*(s-4d0*m2)
      coeb=4d0*pi*alfa**2/als
      coer=alfa**3/als/pi/4d0
      call grid_init

      end subroutine merad_init


!===============================================================================
      subroutine grid_init
! grid initialization
      implicit none
      integer i, n1, n2
      include 'merad_grid.inc'

      ! for grid nt1
      n1 = nt1/4
      n2 = nt1/2
      do i=1,n1
        grt1(i)=0.1d0*dble(i**2)/dble(n1**2)/2d0
        grt1(nt1+1-i)=1d0+grt1(1)-grt1(i)
      enddo
      do i=n1+1,n2
        grt1(i)=(0.1d0+0.9d0*(dble(i-n1))/dble(n2-n1))/2d0
        grt1(nt1+1-i)=1d0+grt1(1)-grt1(i)
      enddo
      ! for grid nz
      n1 = nz/2
      do i=1,n1
        grz(i)=0.5d0*dble(i**2)/dble(n1**2)
      enddo
      do i=1,n1
        grz(nz+1-i)=1d0-0.49d0*dble(i**2)/dble(n1**2)
      enddo

      end subroutine grid_init

!===============================================================================
      real*8 function fspen(x)
! spence function
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/

      if(x) 8,1,1
    1 if(x-.5d0) 2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0) 4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0) 6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0) 10,9,9
    9 fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
   10 fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return

      end function fspen

!===============================================================================
      real*8 function fspens(x)
! core calculation for fspen
      implicit real*8(a-h,o-z)

      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
    1 an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch) 2,2,1
    2 fspens=f
      return

      end function fspens

!===============================================================================
      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
! simpson integration
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
    3 f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
    8 f(k)=10.d16
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
   19 f(k)=10.d16
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

      end subroutine simps

!===============================================================================
      subroutine simpsx(a,b,np,ep,func,res)
! interface for simpson integration
      implicit real*8 (a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)

      end subroutine simpsx

