!===============================================================================
      real(C_DOUBLE) function merad_sig(t, pl, xs_type)
! return the Born level cross section or self-energy and vertex contribution
! tin: the Mandelstam variable t (MeV^2)
! plin: degree of polarization
! xs_type: 0 = Born, others = self-energy and vertex contribution
     &bind(C, name = "merad_sig")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: t, pl
      integer(C_INT), intent(IN), VALUE :: xs_type
      real*8 u,ss,u1,u2,u3,pl1,pl2,pl3
      real*8 dsvt,dsvu,du1,du2,du3,dp1,dp2,dp3,vacpol,l1f
      include 'merad_const.inc'

      u=4*m2-s-t
      ss=s-2d0*m2
      u1=(ss**2+u**2)/2d0+2d0*m2*(s+2d0*t-3d0*m2)
      u2=(ss**2+t**2)/2d0+2d0*m2*(s+2d0*u-3d0*m2)
      u3=(s-2d0*m2)*(s-6*m2)
      pl1=-t*(-t*s**2/2d0/als-ss)
      pl2=-t*(-t*s**2/2d0/als-2d0*m2)+als/2d0-ss*s
      pl3=-(ss**4+4d0*m2*(ss**2*t+m2*(-ss**2+2d0*t**2-4d0*m2*t)))
     .    /als
      if(xs_type.eq.0)then
        merad_sig=coeb*(u1/t**2+u2/u**2+u3/u/t
     .            +pl*(pl1/t**2+pl2/u**2+pl3/u/t))
        return
      else
        dsvt=vacpol(-t)+L1f(t,s,m2)
        dsvu=vacpol(-u)+L1f(u,s,m2)
        PRINT *, vacpol(-t), vacpol(-u)
        du1=dsvt
        du2=dsvu
        du3=dsvt+dsvu
        dp1=dsvt
        dp2=dsvu
        dp3=dsvt+dsvu
        merad_sig=alfa/pi*coeb*(du1*u1/t**2+du2*u2/u**2+du3*u3/u/t/2d0
     .            +pl*(dp1*pl1/t**2+dp2*pl2/u**2+dp3*pl3/u/t/2d0))
      endif

      end function merad_sig

      ! vacuum polarization for leptons and hadrons
      ! input Q^2 in MeV^2
      ! NOTICE: parameterisation for hadron contribution is not accurate at low
      ! Q^2
      real*8 function vacpol(t)
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
      end function vacpol


      real*8 function L1f(xt,xs,xm2)
      implicit none
      real*8 xt,xs,xm2
      include 'merad_const.inc'
      L1f= - 2d0*dlog(abs(xt)/xs)*( dlog(abs(xt)/xm2)-1d0 )
     .     + dlog(abs(xt)/xm2) + (dlog(abs(xt)/xm2))**2
     .     + 4d0*( pi**2/12d0-1d0 )
      return
      end function L1f
