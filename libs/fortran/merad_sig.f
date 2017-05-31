! Subroutines copied from MERADGEN 1.0
! Reference: A. Afanasev, E. Chudakov, A. Ilyichev, V. Zykunov,
!            Comput. Phys. Commun. 176, 218 (2007)
! Added C interfaces, changed unit from GeV to MeV
! Chao Peng, 05/13/2017

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
      real(C_DOUBLE) function merad_sig(t, pl, xs_type)
! return the Born level cross section or self-energy and vertex contribution
! tin: the Mandelstam variable t (MeV^2)
! plin: degree of polarization
! xs_type:
!           0 = Born,
!           1 = self-energy and vertex contribution,
!           others = box diagram contribution
     &bind(C, name = "merad_sig")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: t, pl
      integer(C_INT), intent(IN), VALUE :: xs_type
      real*8 u,ss,u1,u2,u3,pl1,pl2,pl3
      real*8 dsvt,dsvu,du1,du2,du3,dp1,dp2,dp3,vacpol,l1f,xsBt
      include 'merad_const.inc'

      u=4*m2-s-t
      ss=s-2d0*m2
      u1=(ss**2+u**2)/2d0+2d0*m2*(s+2d0*t-3d0*m2)
      u2=(ss**2+t**2)/2d0+2d0*m2*(s+2d0*u-3d0*m2)
      u3=(s-2d0*m2)*(s-6*m2)
      pl1=-t*(-t*s**2/2d0/als-ss)
      pl2=-t*(-t*s**2/2d0/als-2d0*m2)+als/2d0-ss*s
      pl3=-(ss**4+4d0*m2*(ss**2*t+m2*(-ss**2+2d0*t**2-4d0*m2*t)))/als
      ! born
      if(xs_type .eq. 0) then
        merad_sig=coeb*(u1/t**2+u2/u**2+u3/u/t
     .            +pl*(pl1/t**2+pl2/u**2+pl3/u/t))
      ! self energy and vertex contribution
      else if(xs_type .eq. 1) then
        dsvt=vacpol(-t)+L1f(t,s,m2)
        dsvu=vacpol(-u)+L1f(u,s,m2)
        du1=dsvt
        du2=dsvu
        du3=dsvt+dsvu
        dp1=dsvt
        dp2=dsvu
        dp3=dsvt+dsvu
        merad_sig=alfa/pi*coeb*(du1*u1/t**2+du2*u2/u**2+du3*u3/u/t/2d0
     .            +pl*(dp1*pl1/t**2+dp2*pl2/u**2+dp3*pl3/u/t/2d0))
      ! box diagram contribution
      else
        merad_sig=(xsBt(pl,s,t,u)+xsBt(pl,s,u,t))/ss
      endif
      return

      end function merad_sig

!===============================================================================
      real(C_DOUBLE) function merad_sigir(vmin, t, pl)
! return the infrared divergence part of Bremsstrahlung
! vmin: upper limit for "soft" Bremsstrahlung
! t: the Mandelstam variable t (MeV^2)
! pl: degree of polarization
     &bind(C, name = "merad_sigir")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: vmin, t, pl
      real*8 u0, sigborn, dcanc
      include 'merad_const.inc'

      u0=4d0*m2 - s - t
      merad_sigir=sigborn(t,pl)*dcanc(vmin,s,t,u0)
      return

      end function merad_sigir

!===============================================================================
      real*8 function sigborn(t,pl)
! retrun the Born cross section
! t: the Mandelstam variable t (MeV^2)
! pl: degree of polarization
      implicit none
      real*8 t,pl,u,ss,u1,u2,u3,pl1,pl2,pl3
      include 'merad_const.inc'

      u=4d0*m2-s-t
      ss=s-2d0*m2
      u1=(ss**2+u**2)/2d0+2d0*m2*(s+2d0*t-3d0*m2)
      u2=(ss**2+t**2)/2d0+2d0*m2*(s+2d0*u-3d0*m2)
      u3=(s-2d0*m2)*(s-6d0*m2)
      pl1=-t*(-t*s**2/2d0/als-ss)
      pl2=-t*(-t*s**2/2d0/als-2d0*m2)+als/2d0-ss*s
      pl3=-(ss**4+4d0*m2*(ss**2*t+m2*(-ss**2+2d0*t**2-4d0*m2*t)))/als
      sigborn=coeb*(u1/t**2+u2/u**2+u3/u/t+pl*(pl1/t**2+pl2/u**2+pl3/u/t))
      return

      end function sigborn

!===============================================================================
      real*8 function vacpol(t)
! vacuum polarization for leptons and hadrons
! input Q^2 in MeV^2
! NOTICE: parameterisation for hadron contribution is not accurate at low Q^2
      implicit none
      real*8 t,am2,a2,sqlmi,allmi,suml
      real*8 aaa,bbb,ccc,sumh
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
      if(t.lt.1.d6)then
        aaa = -1.345d-9
        bbb = -2.302d-3
        ccc = 4.091*1d-6
      elseif(t.lt.64d6)then
        aaa = -1.512d-3
        bbb =  -2.822d-3
        ccc = 1.218*1d-6
      else
        aaa = -1.1344d-3
        bbb = -3.0680d-3
        ccc = 9.9992d-1*1d-6
      endif
      sumh = -(aaa+bbb*dlog(1d0+ccc*t)) *2d0*pi/alfa
      vacpol=suml + sumh
      return

      end function vacpol

!===============================================================================
      real*8 function L1f(xt,xs,xm2)
      implicit none
      real*8 xt,xs,xm2
      include 'merad_const.inc'

      L1f= - 2d0*dlog(abs(xt)/xs)*( dlog(abs(xt)/xm2)-1d0 )
     .     + dlog(abs(xt)/xm2) + (dlog(abs(xt)/xm2))**2
     .     + 4d0*( pi**2/12d0-1d0 )
      return

      end function L1f

!===============================================================================
      real*8 function xsBt(pl,xs,xt,xu)
! box diagram contribution
      implicit none
      real*8 xs,xt,xu,l1f,dgg1,dgg2,pl
      include 'merad_const.inc'

      xsBt=2d0*alfa**3/xt**2
     .     *((1d0+pl)*xu**2/xs*dgg1(xs,xt,xu)
     .       -(1d0-pl)*xs**2/xu*dgg2(xs,xt,xu))
      return

      end function xsBt

!===============================================================================
      real*8 function dgg1(xs,xt,xu)
! box diagram contribution dgg1
      implicit none
      real*8 xs,xt,xu,ls,lx,dgg
      include 'merad_const.inc'

      LS=dlog(xs/abs(xt))
      LX=dlog(xu/xt)
      dgg = LS**2*(xS**2+xU**2)/2./xt-LS*xU-(LX**2+PI**2)*xU**2/xt
      dgg1=2d0*dlog(xs/abs(xu))*dlog(dsqrt(abs(xu/xs)))-xt/xu**2*dgg
      return

      end function dgg1

!===============================================================================
      real*8 function dgg2(xs,xt,xu)
! box diagram contribution dgg2
      implicit none
      real*8 xs,xt,xu,ls,lx,dgg
      include 'merad_const.inc'

      LS=dlog(xs/abs(xt))
      LX=dlog(xu/xt)
      dgg = LS**2*xS**2/xt+LX*xS-(LX**2+PI**2)*(xS**2+xU**2)/2./xt
      dgg2=2d0*dlog(xs/abs(xu))*dlog(dsqrt(abs(xu/xs)))-xt/xs**2*dgg
      return

      end function dgg2

!===============================================================================
      real*8 function dcanc(xvmin,xs,xt,xu)
! Bremsstrahlung's infrared divergence contribution
      implicit none
      real*8 xvmin,xs,xt,xu,lm,lr,del1s,del1h,del1inf,fspen
      include 'merad_const.inc'
      lm=dlog(-xt/m2)
      lr=dlog(-xu/xs)
      del1s=-2.5d0*lm**2+(3d0-2d0*lr)*lm-lr**2/2d0
     .      -(lm-1)*dlog(xs*(xs+xt)/xt**2)-pi**2/3d0+1d0

      del1h=-lm**2/2d0+(dlog(xt**2*(xs+xt)**2*(xs-xvmin)
     .      /xs/(xvmin-xt)/xvmin/(xs+xt-xvmin)**2)+1d0)*lm
     .      -dlog(-xvmin/xt)**2/2d0-dlog(1d0-xvmin/xt)**2
     .      +dlog((xs+xt)/(xs+xt-xvmin))*dlog((xs+xt)*(xs+xt-xvmin)/xt**2)
     .      +dlog((xvmin-xs)/xt)*dlog(1d0-xvmin/xs)+dlog(-xvmin/xt)
     .      +fspen((xs-xvmin)/xs)-fspen((xt-xvmin)/xt)
     .      +2d0*(fspen(xvmin/xs)-fspen(xvmin/xt)-fspen(xvmin/(xt+xs)))
     .      -pi**2/6d0
      del1inf=4d0*dlog(xvmin/dsqrt(m2*xs))*(dlog(xt*xu/m2/xs)-1d0)
      dcanc=alfa/pi*(del1inf + del1s + del1h)
      return

      end function dcanc

