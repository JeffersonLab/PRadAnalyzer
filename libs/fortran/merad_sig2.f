! Copy from Mehdi Meziane's coding
! Moller radiative cross section beyond ultra-relativistic approximation
! Reference: Eur. Phys. J. A51, 1 (2015)
!            I. Akushevich, H. Gao, A. Ilyichev, and M. Meziane
! Added C interfaces, changed unit from GeV to MeV, Chao Peng, 05/13/2017

!===============================================================================
      real(C_DOUBLE) function merad_sig2(t, xs_type)
! return the Born level cross section or self-energy and vertex contribution
! tin: the Mandelstam variable t (MeV^2)
! xs_type:
!           0 = Born,
!           1 = self-energy and vertex contribution,
!           others = box diagram contribution
     &bind(C, name = "merad_sig2")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: t
      integer(C_INT), intent(IN), VALUE :: xs_type
      real*8 u, sig2, verself2, sigF0, sigFB1, sigFB2
      include 'merad_const.inc'

      u=4d0*m2-t-s
      if(xs_type.eq.0) then
        merad_sig2=sig2(s,t,u)+sig2(s,u,t)
      else if(xs_type.eq.1) then
        merad_sig2=verself2(s,t,u)+verself2(s,u,t)
      else
        merad_sig2=sigF0(s,t,u)+sigFB1(s,t,u)+sigFB2(s,t,u)
     .             +sigF0(s,u,t)+sigFB1(s,u,t)+sigFB2(s,u,t)
      endif
      return

      end function merad_sig2

!===============================================================================
      real(C_DOUBLE) function merad_sigir2(vmin, t)
! return the infrared divergence part of Bremsstrahlung
! vmin: upper limit for "soft" Bremsstrahlung
! t: the Mandelstam variable t (MeV^2)
! pl: degree of polarization
     &bind(C, name = "merad_sigir2")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: vmin, t
      real*8 u, sig2, dcanc2
      include 'merad_const.inc'

      u=4d0*m2-s-t
      merad_sigir2=(sig2(s,t,u)+sig2(s,u,t))*dcanc2(vmin,s,t,u)
      return

      end function merad_sigir2

!*************************BEYOND URA ROUTINES*******************
!ccccccccccccccccccccccBORN ccccccccccccccccccccccccccccccccccc
      real*8 function sig2(xs,xt,xu)
      implicit none
      real*8 xt,xu,xs,u1,u2,u3,pl1,pl2,pl3
      real*8 dsvt,dsvu,du1,du2,du3,dp1,dp2,dp3,vacpol,l1f
      real*8 dbornu,dbornt,um,sm,tm,xsis,xsiu
      integer i
      include 'merad_const.inc'
!        xu=4*m2-xs-xt
      um=(xu-4d0*m2)
      tm=(xt-4d0*m2)
      sm=(xs-4d0*m2)
      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
!        print *,xu,m2,xsiu
      dbornu=xu*xu/4d0/xs/xsis/xsis*(4d0*xsiu**4-(2d0+xt/xu)*
     .(1-xsiu**2)*(1-xsiu**2))-xs*xs*xsis*xsis*xsis*xsis/xu
      sig2=2d0*pi*alfa*alfa*dbornu/(xs*xt**2)
      return

      end function sig2

!ccccccccccccccccc   BOX ccccccccccccccccccccccccccccccccccc

      real*8 function sigF0(xs,xt,xu)
      implicit none
      real*8 xs,xt,xu,fspen,dbornu
      real*8 xsis,xsiu,xsit,au,as
      real*8 um,sm,tm

      include 'merad_const.inc'

      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
      xsit=sqrt(-xt+4d0*m2)/sqrt(-xt)

      um=(xu-4d0*m2)
      tm=(xt-4d0*m2)
      sm=(xs-4d0*m2)

      dbornu= (um**2 - 4d0*m2*m2*(2d0+xt/xu))/sm -sm**2/xu

      as= (1d0/xsis)*(xsis**2+1d0)*(-2d0*pi**2
     . +log((xsis+1d0)/(1d0-xsis))*log((xsis+1d0)/(1d0-xsis))
     . + 4d0*fspen(2d0*xsis/(xsis+1d0)))

      au= (1d0/xsiu)*(xsiu**2+1d0)*(-log((xsiu+1d0)/(-1d0+xsiu))
     . *log((xsiu+1d0)/(-1d0+xsiu))+2d0*log((xsiu+1d0)/(2d0*xsiu))
     . *log((xsiu+1d0)/(2d0*xsiu))+4d0*fspen((xsiu+1d0)/(2d0*xsiu))
     . -pi**2/3d0)

      sigF0=alfa**3/(xs*xt**2)*dbornu*(as+au)
      return

      end function sigF0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function sigFB1(xs,xt,xu)
      implicit none
      real*8 xs,xt,xu,fspen
      real*8 xsis,xsiu,xsit
      real*8 um,sm,tm
      real*8 bf1t1,bf1t2,bf1t3,bf1t4

      include 'merad_const.inc'

      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
      xsit=sqrt(-xt+4d0*m2)/sqrt(-xt)

      um=(xu-4d0*m2)
      tm=(xt-4d0*m2)
      sm=(xs-4d0*m2)

      bf1t1=1d0/(12*xsis*xt)*((xsis**2+1d0)*(xsis**4-6d0*xsis**2-3d0)*xs
     . *xs*xt-2d0*xsis**2*(xsis**2+1d0)**3*xs**3
     .-12d0*xsis**2*xs*xt**2-4d0*
     . xt**3)*(4d0*pi**2+3d0*(log((xsis+1d0)/(-xsis+1d0)))**2
     . - 6d0*(log((xsis+1d0)/(2d0*xsis)))**2
     .- 12d0*fspen((xsis-1d0)/(2d0*xsis))
     .- 6d0*log((xsis+1d0)/(-xsis+1d0))*log(-xsis*xsis*xs/xt))

       bf1t2=1d0/(12d0*xt*xsit**3)*(xsit**2*(xsit**2-3d0)*(3d0*xsit**2
     . +1d0)*xt**3-2d0*xt**2*xu*(3d0*xsit**6+2d0*xsit**4+10d0*xsit**2
     . -1d0)-4d0*xu**2*xt*(5d0*xsit**4+4d0*xsit**2-1d0)
     . -16d0*xsit**2*xu**3)
     . *(4d0*pi**2-6d0*(log((xsit+1d0)/2d0))**2
     . +3d0*(log((xsit+1d0)/(xsit-1d0)))**2
     . -12d0*fspen((-xsit+1d0)/2d0))

       bf1t3=(1d0/xsis)*log((xsis+1d0)/(-xsis+1d0))*(xs*xsis**2+xt)*
     . (xsis**2*(xsis**2+1d0)*xs-2d0*(xsis**2-2d0)*xt)

       bf1t4=log((xsit**2-1d0)/4d0)*(2*xt**2-xs**2*(xsis**4+xsis**2)
     . +(3*xsis**2-1d0)*xs*xt-2d0*xs*(xt+2d0*xu)/xsit**2)

       sigFB1=alfa**3/(xu*xt*(xs**2)*xsis**2)*(bf1t1+bf1t2+bf1t3+bf1t4)
       return

       end function sigFB1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function sigFB2(xs,xt,xu)
      implicit none
      real*8 xs,xt,xu,fspen
      real*8 xsis,xsiu,xsit
      real*8 um,sm,tm
      real*8 bf2t1,bf2t2,bf2t3,bf2t4

      include 'merad_const.inc'

      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
      xsit=sqrt(-xt+4d0*m2)/sqrt(-xt)

      um=(xu-4d0*m2)
      tm=(xt-4d0*m2)
      sm=(xs-4d0*m2)

      bf2t1=1d0/(12*xsiu*xt)*(4d0*xt**3-2d0*(xsiu**4-6d0*xsiu**2-1d0)
     .      *xt**2*xu+(-xsiu**6+xsiu**4+9d0*xsiu**2+7d0)*xt*xu**2
     .      +2d0*xu**3*(xsiu**2+1d0)**3)
     .      *(-6d0*(log((xsiu-1d0)/(2d0*xsiu)))**2
     .      +3d0*(log((xsiu+1d0)/(xsiu-1d0)))**2
      ! CHANGED
      ! the original part (comment below) missed a pair of bracket in the
      ! denominator, but the error is small
!     .     -12d0*fspen((xsiu+1d0)/2d0*xsiu)
     .      -12d0*fspen((xsiu+1d0)/(2d0*xsiu))
     .      +6d0*log((xsiu+1d0)/(xsiu-1d0))*log((xsiu**2)*xu/xt)+pi**2)

      bf2t2=1d0/(12d0*xt*xsit**3)
     .      *((-xsit**4+2d0*xsit**2+3d0)*(xsit**2)*xt**3
     .      +2d0*(xsit**6-4d0*xsit**4+8d0*xsit**2+1d0)*(xt**2)*xu
     .      +4d0*(3d0*xsit**4+1d0)*xt*xu**2+16d0*xsit**2*xu**3)
     .      *(-6d0*(log((xsit+1d0)/2d0))**2
     .      +3d0*(log((xsit+1d0)/(xsit-1d0)))**2
     .      -12d0*fspen((1d0-xsit)/2d0)+4d0*pi**2)

      bf2t3=log((xsit**2-1d0)/4d0)
     .      *(2*xu*(xsit**2*xt+xt+2d0*xu)/xsit**2
     .      +(xt-xu)*(2d0*xt+xsiu**2*xu+xu))

      bf2t4=-(1/xsiu)*log((xsiu+1d0)/(xsiu-1d0))
     .      *(xsiu**2*(xt-xu)-2d0*xt)*(2d0*xt+xsiu**2*xu+xu)

      sigFB2= alfa**3/(xu*xt*xs**2*xsis**2)*(bf2t1+bf2t2+bf2t3+bf2t4)
      return

      end function sigFB2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ************* Cancellation of Infrared Divergence **********

      real*8 function dcanc2(xvmin,xs,xt,xu)
      implicit none
      real*8 xvmin,xs,xt,xu,lm,lr,del1s,eps,fspen
      real*8 um,tm,sm,lambm,sqlm,lmm
      real*8 j01,j02,j03,j0,xsis,xsit,xsiu
      real*8 delt1s, delt1h
      include 'merad_const.inc'

      um=(xu-4d0*m2)
      tm=(xt-4d0*m2)
      sm=(xs-4d0*m2)

      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
      xsit=sqrt(-xt+4d0*m2)/sqrt(-xt)

!      j01=(xs-2d0*m2)/(sqrt(xs*sm))*log((sqrt(sm)+sqrt(xs))**2/4d0/m2)
!      j02=(xt-2d0*m2)/(sqrt(xt*tm))*log((sqrt(-tm)+sqrt(-xt))**2/4d0/m2)
!      j03=(xu-2d0*m2)/(sqrt(xu*um))*log((sqrt(-um)+sqrt(-xu))**2/4d0/m2)

      j01=((xsis*xsis+1d0)/xsis)*log((1d0+xsis)/(1d0-xsis))
      j02=-((xsit*xsit+1d0)/xsit)*log((1d0+xsit)/(-1d0+xsit))
      j03=-((xsiu*xsiu+1d0)/xsiu)*log((1d0+xsiu)/(-1d0+xsiu))
      ! CHANGED
      ! the original formula commented out the j0 expression beyond URA,
      ! and used URA form, here I commented out the URA form and used beyond
      ! URA form
      j0=-2d0*(j01+j02+j03+2d0)
!      j0=-4d0*(dlog(m2*xs/xt/xu)+1d0)
      dcanc2=alfa/pi*(j0*log(xvmin/m2)+delt1h(xvmin,xs,xt,xu)+delt1s(xs,xt,xu))
      return

      end function dcanc2

      real*8 function delt1h(xvmin,xs,xt,xu)
      implicit none
      real*8 xv,xs,xt,xu,fspen,xvmin
      real*8 j01,j02,j03,j0,xsis,xsit,xsiu
      real*8 z1u,z2u,hs,ht,zs,zt,z1,z2,z3,z4,li2
!      external delt1hv
      include 'merad_const.inc'
!      include 'merad_xx.inc'
      xv=xvmin
      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
      xsit=sqrt(-xt+4d0*m2)/sqrt(-xt)

      z1u=sqrt(xv-xsiu*xsiu*(xv+xu))/sqrt(-xu)/xsiu
      z2u=sqrt((xv+xsiu*xsiu*xu)/(xv+xu))/xsiu
      zs=(xsis/xv)*(sqrt(xsis*xsis*xs*xs-2d0*xs*xv+xv*xv)-xsis*xs)+1d0
      zt=(xsit/xv)*(sqrt(xsit*xsit*xt*xt-2d0*xt*xv+xv*xv)-xsit*xt)+1d0
      z1=1d0+xsis
      z2=(1d0+xsis)*(1d0+xsis)/(1d0-xsis)
      z3=1d0-xsis
      z4=(1d0-xsis)*(1d0-xsis)/(1d0+xsis)

      hs=((xsis*xsis+1d0)/2d0/xsis)*(fspen(zs/z1)+fspen(zs/z2)
     .-fspen(zs/z3)-fspen(zs/z4)
     .-dlog(abs((1d0+xsis)*(1d0+xsis)/(xsis-1d0)/(xsis-1d0)))
     .*dlog(abs(((zs-1d0)*(zs-1d0)-xsis*xsis)/(1d0-xsis*xsis))))

      z1=1d0+xsit
      z2=(1d0+xsit)*(1d0+xsit)/(1d0-xsit)
      z3=1d0-xsit
      z4=(1d0-xsit)*(1d0-xsit)/(1d0+xsit)
      ht=((xsit*xsit+1d0)/2d0/xsit)*(fspen(zt/z1)+fspen(zt/z2)
     .-fspen(zt/z3)-fspen(zt/z4)
     .-dlog(abs((1d0+xsit)*(1d0+xsit)/(xsit-1d0)/(xsit-1d0)))
     .*dlog(abs(((zt-1d0)*(zt-1d0)-xsit*xsit)/(1d0-xsit*xsit))))

      li2=fspen(4d0*xsiu/(1d0+xsiu)/(1d0+xsiu))
     . -fspen(-4d0*xsiu/(-1d0+xsiu)/(-1d0+xsiu))
     .-2d0*fspen(2d0*xsiu/(-1d0+xsiu))+2d0*fspen(2d0*xsiu/(1d0+xsiu))
     .+fspen(2d0*xsiu*(z1u-1d0)/(-1d0+xsiu)/(-1d0+xsiu))
     .+fspen(-2d0*xsiu*(z1u+1d0)/(-1d0+xsiu)/(-1d0+xsiu))
     .-fspen(-2d0*xsiu*(z1u-1d0)/(1d0+xsiu)/(1d0+xsiu))
     .-fspen(2d0*xsiu*(z1u+1d0)/(1d0+xsiu)/(1d0+xsiu))
     .+2d0*fspen(-xsiu*(z2u-1d0)/(-1d0+xsiu))
     .+2d0*fspen(xsiu*(z2u+1d0)/(-1d0+xsiu))
     .-2d0*fspen(xsiu*(z2u+1d0)/(1d0+xsiu))
     .-2d0*fspen(xsiu*(1d0-z2u)/(1d0+xsiu))
     .+2d0*dlog(abs((xsiu+1d0)/(xsiu-1d0)))
     .*dlog(abs((xsiu*xsiu*z2u*z2u-1d0)/(xsiu*xsiu-1d0)))


      delt1h=dlog(1d0+xv/m2)
     .+hs-ht+(xsiu*xsiu+1d0)*li2/2d0/xsiu
!      print *,li2
      return

      end function delt1h

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function H1t(xvmin,xs,xt,xu)
      implicit none
      real*8 xv,xs,xt,xu,fspen,xvmin
      real*8 j01,j02,j03,j0,xsis,xsit,xsiu
      real*8 z1u,z2u,hs,ht,zs,zt,z1,z2,z3,z4,li2
!      external delt1hv
      include 'merad_const.inc'
!      include 'merad_xx.inc'
      xv=xvmin
      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
      xsit=sqrt(-xt+4d0*m2)/sqrt(-xt)

      z1u=sqrt(xv-xsiu*xsiu*(xv+xu))/sqrt(-xu)/xsiu
      z2u=sqrt((xv+xsiu*xsiu*xu)/(xv+xu))

      zs=(xsis/xv)*(sqrt(xsis*xsis*xs*xs-2d0*xs*xv+xv*xv)-xsis*xs)+1d0
      zt=(xsit/xv)*(sqrt(xsit*xsit*xt*xt-2d0*xt*xv+xv*xv)-xsit*xs)+1d0

      z1=1d0+xsit
      z2=(1d0+xsit)*(1d0+xsit)/(1d0-xsit)
      z3=1d0-xsit
      z4=(1d0-xsit)*(1d0-xsit)/(1d0+xsit)


      H1t=((xsit*xsit+1d0)/2d0/xsit)*(fspen(zt/z1)+fspen(zt/z2)
     .-fspen(zt/z3)-fspen(zt/z4)
     .-dlog(abs((1d0+xsit)*(1d0+xsit)/(xsit-1d0)/(xsit-1d0)))
     .*dlog(abs(((zt-1d0)*(zt-1d0)-xsit*xsit)/(1d0-xsit*xsit))))

!      print *,li2
      return

      end function H1t


      real*8 function delt1s(xs,xt,xu)
      implicit none
      real*8 xs,xt,xu,fspen,del1st
      real*8 del1ss,del1su,xsis,xsit,xsiu,sphi
      include 'merad_const.inc'

      xsis=sqrt(xs-4d0*m2)/sqrt(xs)
      xsiu=sqrt(-xu+4d0*m2)/sqrt(-xu)
      xsit=sqrt(-xt+4d0*m2)/sqrt(-xt)

      del1ss=dlog((xsis+1d0)/(1d0-xsis))
     .*(1d0+dlog((xsis+1d0)/(1d0-xsis)))
     .+fspen(4d0*xsis/(xsis+1d0)/(xsis+1d0))

      del1st=dlog((xsit+1d0)/(xsit-1d0))
     .*(dlog((xsit+1d0)/(xsit-1d0))-1d0)
     .+fspen(4d0*xsit/(xsit+1d0)/(xsit+1d0))

      del1su=dlog((xsiu+1d0)/(xsiu-1d0))
     .*(dlog((xsiu+1d0)/(xsiu-1d0))-1d0)
     .+fspen(4d0*xsiu/(xsiu+1d0)/(xsiu+1d0))

      delt1s=(xsis*xsis+1d0)*del1ss/2d0/xsis
     .     -(xsit*xsit+1d0)*del1st/2d0/xsit
     .     -(xsiu*xsiu+1d0)*del1su/2d0/xsiu+1d0
     .-sphi(-xu*(1d0+xsiu**2),xs*(1d0+xsis**2),-xt*(1d0+xsit**2))
     .+sphi(-xu*(1d0+xsiu**2),-xt*(1d0+xsit**2),xs*(1d0+xsis**2))
     .-sphi(-xt*(1d0+xsit**2),xs*(1d0+xsis**2),-xu*(1d0+xsiu**2))
      return

      end function delt1s

      real*8 function sphi(xs,xt,xu)
      implicit none
      real*8 sj(4),delij(4,4),z(4),sl(4),z2(6,6)
      real*8 sphiu,sphid,sphiu4,sphid4,zu,zd
      real*8 xs,xt,xu,fspen
      include 'merad_const.inc'
      integer*4 i,j
       do i=1,4
        sj(i)=0d0
        z(i)=0d0
        sl(i)=0d0
      enddo

      do i=1,4
       do j=1,4
       delij(i,j)=0d0
       enddo
      enddo

      sj(1)=1d0
      sj(2)=1d0
      sj(3)=-1d0
      sj(4)=-1d0

      do i=1,4
       do j=1,4
        if (j .eq. i) then
          delij(i,j)=1d0
          else
          delij(i,j)=0d0
        endif
        enddo
      enddo

      sl(1)=sqrt(xs*xs-16d0*m2*m2)
      sl(2)=sqrt(xt*xt-16d0*m2*m2)
      sl(3)=sqrt(xu*xu-16d0*m2*m2)


      z(1)=(4d0*m2*(xt+sl(2))/(xu+sl(3))-xs)/sl(2)-1d0
      z(2)=((xt+sl(2))*(xu+sl(3))/4d0/m2-xs)/sl(2)-1d0
      z(3)=(-4d0*m2*(xu+sl(3))/(xt+sl(2))+xs)/sl(2)-1d0
      z(4)=(-64d0*m2**3/(xu+sl(3))/(xt+sl(2))+xs)/sl(2)-1d0


      zu=sl(1)/sl(2)-1d0

      zd=-1d0+(xs*xt-4d0*m2*xu)/sl(2)**2
      do i=1,4
       do j=1,4
       z2(i,j)=z(i)-z(j)
       enddo
      enddo
      z2(5,1)=zu-z(1)
      z2(5,2)=zu-z(2)
      z2(5,3)=zu-z(3)
      z2(5,4)=16d0*m2**2*(4d0*m2/(xt+sl(2))/(xu+sl(3))
     .-1d0/(xs+sl(1)))/sl(2)
      z2(6,1)=zd-z(1)
      z2(6,2)=zd-z(2)
      z2(6,3)=4d0*m2*(4d0*m2*xs-xu*xt+sl(2)*sl(3))/(xt+sl(2))/sl(2)**2
      z2(6,4)=zd-z(4)


      sphiu=dlog(16d0*m2**2/(xt+sl(2))**2)
     .*dlog(abs(z2(5,1)*z2(5,3)/z2(5,2)/z2(5,4)))


      sphid=dlog(16d0*m2**2/(xt+sl(2))**2)
     .*dlog(abs(z2(6,1)*z2(6,3)/z2(6,2)/z2(6,4)))


      sphiu4=0.0
      sphid4=0.0

      do i=1,4
       do j=1,4
        if (i .eq. j) then
        sphiu4=sphiu4+sj(j)*(-1d0)**(i+1)*0.5*dlog(abs(z2(5,i)))**2
        sphid4=sphid4+sj(j)*(-1d0)**(i+1)*0.5*dlog(abs(z2(6,i)))**2


         else

         sphiu4=sphiu4+sj(j)*((-1d0)**(i+1))
     .   *(dlog(abs(z2(5,i)))*dlog(abs(z2(i,j)))
     .   -fspen(z2(5,i)/z2(j,i)))


         sphid4=sphid4+sj(j)*((-1d0)**(i+1))
     .   *(dlog(abs(z2(6,i)))*dlog(abs(z2(i,j)))
     .   -fspen(z2(6,i)/z2(j,i)))
        endif

        enddo
      enddo

      sphi= xu*(sphiu+sphiu4-sphid-sphid4)/2d0/sl(3)
      return

      end function sphi

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ****************** anomalous magnetic moment **************************
      real*8 function damag(xs,xt,xu)
      implicit none
      include 'merad_const.inc'
      real*8 xs,xt,xu,amag1,xsit

      ! CHANGED
      ! the original formula (comments below) is different with the paper's
      ! formula, and the dimension is incorrect
      ! its difference may be huge due to the missing t, but the final
      ! contribution to the total cross section is still very small
!      amag1=log((sqrt(4d0*m2-xt)+sqrt(-xt))**2/(4d0*m2))*
!     .(3d0*(xs-2d0*m2)/xu+(10d0*m2-3d0*xu)/(xs-4d0*m2))/xt
!      damag=-8d0*m2*alfa**3*(amag1)/xs

      xsit=sqrt(1-4d0*m2/xt)
      amag1=log((xsit+1)/(xsit-1))
     .      *(3d0*(xs-2d0*m2)/xu+(10d0*m2-3d0*xu)/(xs-4d0*m2))/xsit

      damag=-4d0*m2*alfa**3*amag1/(xs*xt**2)
      return

      end function damag

!ccccccccccccccccccccccvertex + self energy ccccccccccccccccccccccccccccc


      real*8 function verself2(xs,xt,xu)
      implicit none
      real*8 xt,xu,xs
      real*8 vacpol2,deltvert,fspen,sig2,damag
      real*8 lambm,sqlm,lmm,tmm
      integer i
      include 'merad_const.inc'

      tmm=-xt+2d0*m2
      lambm=xt**2-4d0*m2*xt
      sqlm=dsqrt(lambm)
      lmm=(dlog((sqlm-xt)/(sqlm+xt)))/sqlm
      deltvert=-2d0+(-1.5d0*xt+4d0*m2)*lmm-(tmm/sqlm)*(.5d0*lambm*lmm**2
     . +2d0*fspen((2d0*sqlm)/(-xt+sqlm))-pi**2/2d0)

      ! CHANGED
      ! using sig2 instead of calculating dbornt
      ! add AMM contribution to it
      ! using vacpol2 instead of vacpol (hadron part discarded)
      verself2= sig2(xs,xt,xu)*alfa/pi*(vacpol2(-xt)+2d0*deltvert)
     .          +damag(xs,xt,xu)

      return

      end function verself2

!===============================================================================
      real*8 function vacpol2(t)
! vacuum polarization for leptons
! hadrons contribution is discarded since the parameterisation is not working
! well for low Q^2
! input Q^2 in MeV^2
! NOTICE: parameterisation for hadron contribution is not accurate at low Q^2
      implicit none
      real*8 t,am2,a2,sqlmi,allmi,suml
      integer i
      include 'merad_const.inc'
      dimension am2(3)
      !am2 : squared masses of charge leptons
      data am2/0.2611198942d0,1.11636921d4,3.157089d6/

      suml=0.
      do i=1,3
        a2=2.*am2(i)
        sqlmi=dsqrt(t*t+2.*a2*t)
        allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
        suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      enddo
      vacpol2=suml
      return

      end function vacpol2

!cccccccccccccccccccccccEND OF ROUTINE BEYOND URAcccccccccccccccccccccccc

