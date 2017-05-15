! Subroutines copied from MERADGEN 1.0
! Reference: A. Afanasev, E. Chudakov, A. Ilyichev, V. Zykunov,
!            Comput. Phys. Commun. 176, 218 (2007)
! Added C interfaces, changed unit from GeV to MeV
! Chao Peng, 05/13/2017

!===============================================================================
      subroutine merad_init(elab)
! Initialize MERADGEN by beam energy, input is in MeV
     &bind(C, name = "merad_init")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: elab
      include 'merad_const.inc'
      data pi/3.14159265359d0/, alfa/.7297352568d-2/,
     .     m/0.510998918d0/

      m2 = m*m
      En=elab
      s=2d0*(En*m+m2)
      als=s*(s-4d0*m2)
      coeb=4d0*pi*alfa**2*(s-2d0*m2)/als
      coer=alfa**3/als/pi*(s - 2d0*m2)/4d0
      call grid_init

      end

!===============================================================================
      real(C_DOUBLE) function merad_sigfs(vmin, tin, plin)
! return the "soft" Bremsstrahlung cross section
! vmin: the minimum v for photons to be detected (MeV^2)
! tin: the Mandelstam variable t (MeV^2)
! plin: degree of polarization
     &bind(C, name = "merad_sigfs")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: vmin, tin, plin
      real*8 fsirsoft, xs, sigborn
      integer nn
      include 'merad_tv.inc'

      t = tin
      pl = plin
      sig0 = sigborn(t, pl)
      ! redundant, but it prevents crash
      xs = fsirsoft(vmin)
      call simpsx(1d-22,vmin,10000,1d-3,fsirsoft,xs)
      merad_sigfs = xs

      end

!===============================================================================
      real*8 function fsirsoft(v)
! a helper function to use the simpson integration subroutine
! kinematics t, v are shared through a common block
      implicit none
      real*8 v,fsir
      integer nn
      include 'merad_tv.inc'

      fsirsoft=fsir(t,0d0,v,0d0,pl,nn,-1,sig0)

      end

!===============================================================================
      real(C_DOUBLE) function merad_sigfh(vmin, vmax, tin, plin)
! retrun the "hard" Bremsstrahlung cross section
! vmin: the minimum v for photons to be detected (MeV^2)
! vmax: v cuts on the hard photons to be calculated, cannot exceed the
!       kinematics limit (MeV^2)
! tin: the Mandelstam variable t (MeV^2)
! plin: degree of polarization
     &bind(C, name = "merad_sigfh")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: vmin, vmax, tin, plin
      real*8 va, vb, sia, sib, fsir
      integer nn, iv
      include 'merad_tv.inc'
      include 'merad_grid.inc'

      ! this hard photon part will be revisited again
      ! so use grid to calculate it and store the values
      ! Trapezoid rule
      vb = vmin
      sib = 0d0
      do iv = 1, nv
        va = vb
        sia = sib
        vb = vmin + (vmax - vmin)*grv(iv)
        sib = fsir(tin, 0d0, vb, 0d0, plin, nn, 2, 0d0)
        distsiv(iv) = distsiv(iv-1) + (sib + sia)*(vb-va)/2d0
        distarv(iv) = vb
      enddo
      merad_sigfh = distsiv(nv)

      end

!===============================================================================
      double precision function sigborn(t,pl)
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

      end

!===============================================================================
      subroutine zd(t,t1,v)
! helper subroutine to calculate some frequently used variables
! should be called before fsir when the ikey to fsir is not 2 or -1
      implicit none
      real*8 t,t1,v,u
      include 'merad_const.inc'
      include 'merad_gr.inc'

      u=v-s-t+4d0*m2

      az=(v-t)**2-4d0*m2*t
      bz=-(v*(2d0*m2*(t+t1)+t1*(v-t)))
     .   +s*(-t**2+t1*v+t*(t1+v))
      cz=(s*(t-t1)+t1*v)**2
     .   -4d0*m2*(s*(t-t1)**2+t1*v**2)

      az1=az
      bz1=-(t*(s+t-4d0*m2)*(t-t1))
     .    +(t*(2*t-t1)-2d0*m2*(t+t1)+s*(t+t1))*v - t*v**2
      cz1=((s+t)*(t-t1)-t*v)**2+4d0*m2
     .    *(-((s+t)*(t-t1)**2)+(t-t1)*(t+t1)*v-t1*v**2)

      az2=az
      bz2=(4d0*m2-s-t)*t*(4d0*m2-t1)
     .    +(6d0*m2*t-s*t-2d0*m2*t1+s*t1-t*t1)*v+(-4d0*m2 + s)*v**2
      cz2=u*(-4d0*m2*(s+t1-4d0*m2)**2
     .    +(16d0*m2**2+t1**2-4d0*m2*(s+2d0*t1))*u)
     .    -2d0*(2d0*m2-t1)*(4d0*m2-s-t1)*u*v+(s+t1-4d0*m2)**2*v**2

      end

!===============================================================================
      double precision function fsir(t,t1,v,z,pl,nn,ikey,sig0)
! The cross section of real photon emission
! t,t1,v,z are kinematic invariant
! pl is degree of polarization
! nn number of point where cross section is negative
! ikey=0 radiative cross section ds/dt/dv/dt1/dz
! ikey=1 radiative cross section integrated over z ds/dt/dv/dt1
! ikey=2 radiative cross section integrated over z and t1 ds/dt/dv
! ikey=-1 finite  part of
! radiative cross section integrated over z and t1 ds/dt/dv
      implicit none
      real*8 t,t1,v,z,pl,u,uu,vv,v1,vv1,z1,z2,zz,zz1,zz2,ds
      real*8 dz,dz1,vt,pls,tt,tt1,dt,sd,fir,sig0
      integer nn,ikey
      real*8 sr1,sr2,sr3,sr4,sr5,sr6,sr7,sr8,sr9,sr10
      real*8 aj1,aj2,aj3,aj4,aj5,aj6,aj7,aj8,aj9,aj10
      real*8 aj11,aj12,aj13,aj14,aj15,aj16,aj17,aj18,aj19,aj20
      real*8 aj21,aj22,aj23,aj24,aj25,aj26,aj27,aj28,aj29,aj30
      real*8 aj31,aj32,aj33,aj34,aj35,aj36,aj37,aj38,aj39,aj40
      real*8 aj41,aj42,aj43,aj44,aj45,aj46,aj47,aj48,aj49,aj50
      real*8 aj51
      include 'merad_const.inc'
      include 'merad_gr.inc'
      u=v-s-t+4d0*m2
      uu=1d0/u
      v1=s+u+t1-4d0*m2
      z1=z-t1+t
      z2=z-s-t1+4.*m2
      zz=1d0/z
      zz1=1d0/z1
      zz2=1d0/z2
      vv1=1d0/v1
      vv=1d0/v
      pls=-pl/als
      tt=1d0/t
      tt1=1d0/t1
      dt=1d0/(t-t1)
      dz=1d0/(s+t1-4*m2)
      dz1=1d0/(s+t-4*m2)
      ds=1d0/(s-4*m2)
      vt=1d0/(v-t)
      if(ikey.eq.0)then
       sd=1d0/dsqrt(-az*z**2-2d0*bz*z-cz)
!        sd=1d0/dsqrt(az*(zmax-z)*(z-zmin))
      aj1=sd
      aj2=z*sd
      aj3=z**2*sd
      aj4=zz*sd
      aj5=m2*zz**2*sd
      aj6=zz1*sd
      aj7=m2*zz1**2*sd
      aj8=zz2*sd
      aj9=m2*zz2**2*sd
      aj10=t1*aj1
      aj11=tt1*aj1
      aj12=m2*tt1**2*aj1
      aj13=vv1*aj1
      aj14=m2*vv1**2*aj1
      aj15=tt1*aj2
      aj16=vv1*aj2
      aj17=m2*vv1**2*aj2
      aj18=vv1*aj3
      aj19=m2*vv1**2*aj3
      aj20=t1*aj4
      aj21=tt1*aj4
      aj22=m2*tt1**2*aj4
      aj23=dt*(aj4-aj6)
      aj24=vv1*aj4
      aj25=dz*aj4
      aj26=m2*dz**2*aj4
      aj27=m2*dz**3*aj4
      aj28=tt1*aj5
      aj29=tt1**2*aj5
      aj30=dz*aj5
      aj31=dz**2*aj5
      aj32=t1*aj6
      aj33=t1**2*aj6
      aj34=tt1*aj6
      aj35=m2*tt1**2*aj6
      aj36=vv1*aj6
      aj37=t1*aj7
      aj38=t1**2*aj7
      aj39=tt1*aj7
      aj40=tt1**2*aj7
      aj41=t1*aj8
      aj42=tt1*aj8
      aj43=vv1*aj8
      aj44=m2*vv1**2*aj8
      aj45=dz*aj8
      aj46=m2*dz**2*aj8
      aj47=m2*dz**3*aj8
      aj48=vv1*aj9
      aj49=vv1**2*aj9
      aj50=dz*aj9
      aj51=dz**2*aj9
      elseif(ikey.eq.1)then
      aj1=pi/sqrt(az)
      aj2=-pi*bz/az**(3d0/2d0)
      aj3=pi*(3d0*bz**2-az*cz)/2d0/az**(5d0/2d0)
      aj4=pi/sqrt(cz)
      aj5=-m2*pi*bz/cz**(3d0/2d0)
      aj6=pi/sqrt(cz1)
      aj7=-m2*pi*bz1/cz1**(3d0/2d0)
      aj8=-pi/sqrt(cz2)
      aj9=m2*pi*bz2/cz2**(3d0/2d0)
      aj10=t1*aj1
      aj11=tt1*aj1
      aj12=m2*tt1**2*aj1
      aj13=vv1*aj1
      aj14=m2*vv1**2*aj1
      aj15=tt1*aj2
      aj16=vv1*aj2
      aj17=m2*vv1**2*aj2
      aj18=vv1*aj3
      aj19=m2*vv1**2*aj3
      aj20=t1*aj4
      aj21=tt1*aj4
      aj22=m2*tt1**2*aj4
      aj23=dt*(aj4-aj6)
      aj24=vv1*aj4
      aj25=dz*aj4
      aj26=m2*dz**2*aj4
      aj27=m2*dz**3*aj4
      aj28=tt1*aj5
      aj29=tt1**2*aj5
      aj30=dz*aj5
      aj31=dz**2*aj5
      aj32=t1*aj6
      aj33=t1**2*aj6
      aj34=tt1*aj6
      aj35=m2*tt1**2*aj6
      aj36=vv1*aj6
      aj37=t1*aj7
      aj38=t1**2*aj7
      aj39=tt1*aj7
      aj40=tt1**2*aj7
      aj41=t1*aj8
      aj42=tt1*aj8
      aj43=vv1*aj8
      aj44=m2*vv1**2*aj8
      aj45=dz*aj8
      aj46=m2*dz**2*aj8
      aj47=m2*dz**3*aj8
      aj48=vv1*aj9
      aj49=vv1**2*aj9
      aj50=dz*aj9
      aj51=dz**2*aj9
      elseif(ikey.eq.2.or.ikey.eq.-1)then
      aj1=pi*v/(m2+v)
      aj2=-pi*v**2*(v-s+2d0*m2)/2d0/(m2+v)**2
      aj3=0d0!pi*(3d0*bz**2-az*cz)/2d0/az**(5d0/2d0)
      aj4=pi/sqrt((s-v)**2-4d0*m2*s)*
     .    dlog((s-v-2d0*m2+sqrt((s-v)**2-4d0*m2*s))**2/4d0/m2/(m2+v))
      aj5=pi/v
      aj6=pi/sqrt((v-u)**2-4d0*m2*u)*
     .    dlog((v-u+2d0*m2+sqrt((v-u)**2-4d0*m2*u))**2/4d0/m2/(m2+v))
      aj7=pi/v
      aj8=pi*dlog(4d0*m2*u**2*(m2+v)/(v*(v-u+sqrt((v-u)**2-4d0*m2*u))
     .    -2d0*m2*u)**2)/dsqrt((v-u)**2-4d0*m2*u)
      aj9=pi*v/u**2
      aj10=pi*v*(2d0*m2*t+(t-v)*v)/2d0/(m2+v)**2
      aj11=pi*dlog(4d0*m2*t**2*(m2+v)/(v*(v-t+sqrt((v-t)**2-4d0*m2*t))
     .     -2d0*m2*t)**2)/dsqrt((v-t)**2-4d0*m2*t)
      aj12=pi*v/t**2
      aj13=pi*dlog((v-t+4d0*m2+sqrt((v-t)**2-4d0*m2*t))**2
     .     /4d0/m2/(m2+v))/sqrt((v-t)**2-4d0*m2*t)
      aj14=pi/v
      aj15=pi*v*(v*(2d0*m2-t+v)-s*(t+v))/((v-t)**2-4d0*m2*t)/(m2+v)
     .     +pi*t*(s*(t-v)+2d0*m2*v)/((v-t)**2-4d0*m2*t)**(1.5d0)
     .     *dlog(4d0*m2*t**2*(m2+v)/(v*(v-t+sqrt((v-t)**2-4d0*m2*t))
     .     -2d0*m2*t)**2)
      aj16=pi*v*(v*(2d0*m2-t+v)-s*(t+v))/((v-t)**2-4d0*m2*t)/(m2+v)
     .     -pi*v*(u*(t-v)+2d0*m2*v)/((v-t)**2-4d0*m2*t)**(1.5d0)
     .     *dlog(4d0*m2*(m2+v)/(v-t+sqrt((v-t)**2-4d0*m2*t)-2d0*m2)**2)
      aj17=-pi*(1d0+(s*(t-v)+2d0*m2*v)/((v-t)**2-4d0*m2*t))
     .     +m2*pi*(s*(t+v)+v*(t-v-2d0*m2))/((v-t)**2-4d0*m2*t)**(1.5d0)
     .     *dlog(4d0*m2*(m2+v)/(v-t+sqrt((v-t)**2-4d0*m2*t)-2d0*m2)**2)
      aj18=pi*((v**2*(4d0*m2**3*(4d0*s*t+(8d0*t-3d0*v)*v)
     .     -2d0*m2**2*(2d0*s*t*(s+7d0*t)-2d0*(s - 6d0*t)*t*v
     .     -(4d0*s+31d0*t)*v**2+13*v**3)+(t-v)*(3d0*(t-v)**2*v**2
     .     +2d0*s*(t-v)*v*(2d0*t+3d0*v)+s**2*(-t**2+4d0*t*v+3d0*v**2))
     .     +2d0*m2*(s**2*(4d0*t**2-t*v-v**2)+s*(t+v)
     .     *(3d0*t**2-15d0*t*v+8d0*v**2)+(t-v)*v*(2d0*t**2-13d0*t*v
     .     +8d0*v**2))))/(2d0*(-4d0*m2*t+(t-v)**2)**2*(m2+v)**2)
     .     -(v**2*(u**2*(t-v)**2+6d0*m2**2*v**2
     .     -2d0*m2*u*(s*t+2d0*v*(v-t)))
     .     *dlog((4*m2*(m2+v))/(2d0*m2-t
     .     +dsqrt(-4d0*m2*t+(t-v)**2)+v)**2))
     .     /(-4d0*m2*t+(t-v)**2)**2.5)
      aj19=pi*((2d0*v*((u-4d0*m2)**2*(t-v)**2*v
     .     -4d0*m2**3*(4d0*(s-t)*t
     .     +4d0*t*v-3d0*v**2)+4d0*m2**2*((s-2d0*t)*t*(s+t)
     .     -3d0*(s-3d0*t)*t*v-(2d0*s+9d0*t)*v**2+4d0*v**3)
     .     +m2*(2d0*s**2*(t+v)**2+2d0*s*(t-v)*(t**2-3d0*t*v+4d0*v**2)
     .     +(t-v)**2*(t**2-10*t*v+6*v**2))))
     .     /(((t-v)**2-4d0*m2*t)**2*(m2+v))+m2
     .     *(4d0*v*(-2d0*m2**2*(4d0*s*t+(4d0*t-3d0*v)*v)+m2*(2d0*s*t*(s+t)
     .     +2d0*t*(s+t)*v-3d0*t*v**2+v**3)+u*(t-v)*((t-v)*v
     .     +s*(2d0*t+v)))*dlog((4d0*m2*(m2+v))/(2d0*m2-t
     .     +dsqrt(-4*m2*t+(t-v)**2)+v)**2))/((t-v)**2-4d0*m2*t)**2.5)/2d0
      aj20=-((s*(t+v)-v*(v-t+2d0*m2))*aj1
     .     +(s*t*(v-s)+m2*(4d0*s*t-2d0*v**2))*aj4)/((s-v)**2-4d0*m2*s)
      aj21=pi/t*dlog(((s-2d0*m2)*(s-2d0*m2+dsqrt(als)))/2d0/m2**2
     .     -1d0)/dsqrt(als)
      aj22=-pi/t**3/als*((s*(v-t)-2d0*m2*v)*v
     .     +dlog(((s-2d0*m2)*((s-2d0*m2)+dsqrt(als)))/2d0/m2**2
     .     -1d0)*(t*(s*v-als)-2d0*m2*v**2)*m2/dsqrt(als))
      aj23=2d0*pi*dlog((dsqrt(t**2-4d0*m2*t)-t+2d0*m2)/2d0/m2)
     .     /v/dsqrt(t**2-4d0*m2*t)
      aj24=2d0*pi*dlog((dsqrt(u**2-4d0*m2*u)-u+2d0*m2)/2d0/m2)
     .     /v/dsqrt(u**2-4d0*m2*u)
      aj25=-pi*dlog((m2**2*(v**2-2d0*u**2)+(m2*v-(s-2d0*m2)*u)**2
     .     +dsqrt(als)*u*(s*u-2d0*m2*(u+v)))/2d0/m2/(m2*(u+v)**2-s*u*v))
     .     /u/dsqrt(als)
      aj26=(-pi/als/u**2*(v*(s*u*(s+u-2d0*v)+2d0*m2*(v*(u+v)-2d0*s*u))
     .     /(s*u*v-m2*(u+v)**2))
     .     +(s*u*(v-s)+m2*(4d0*s*u-2d0*v**2))*aj25/u**2/als)*m2
      aj27=(-pi*v*(s**3*u**3*v*(u**2-3d0*s**2-2d0*s*u+2d0*(5d0*s+u)*v
     .     -6d0*v**2)+2d0*m2*s**2*u**2*(2d0*s*u**2*(s+u)-u*(s*u-15d0*s**2
     .     +4d0*u**2)*v+(s**2-34d0*s*u+2d0*u**2)*v**2
     .     -(11d0*s - 5d0*u)*v**3 + 12d0*v**4) + 2d0*m2**2*s*u
     .     *(2d0*u**2*(13d0*s*u-24d0*s**2+u**2)*v-8d0*s*u**3*(2d0*s + u)
     .     -2d0*s*(4d0*s-37d0*u)*u*v**2+(62d0*s-11d0*u)*u*v**3
     .     +6d0*(s-4d0*u)*v**4-15d0*v**5)+4d0*m2**3*(u+v)
     .     *(3d0*v**3*(u+v)**2-12d0*s*u*v*(u+v)**2+8d0*s**2*u**2*(2d0*u
     .     +v))))*m2/(2d0*u**4*(s*u*v-m2*(u+v)**2)**2*als**2)
     .     +(aj25*(2d0*m2*s*u*(3d0*s+u-3d0*v)*v**2+6d0*m2**2*v**2
     .     *(v**2-4d0*s*u)+u**2*(als-s*v)**2))*m2/(u**4*als**2)
      aj28=m2*pi/t**2/als*((s*(v-t)-2d0*m2*v)
     .      *dlog((s-2d0*m2)*(s-2d0*m2+dsqrt(als))/2d0/m2
     .      -1d0)/dsqrt(als)-(t*(s*v-als)-2d0*m2*v**2)/m2/v)
      aj29=-pi/t**4/als**(2.5d0)*(dsqrt(als)*(als*t**2*(2d0*s*v-als)
     .     +v**2*(2d0*s*t*v*(s+2d0*m2)-v**2*((s-2d0*m2)**2+8d0*m2)
     .     -2d0*s*t*(s*t+2d0*(2d0*s+t)*m2-16d0*m2**2)))/v
     .     +2d0*m2*(2d0*m2**2*v*(8d0*s*t-3d0*v**2)+s*v*m2*(3d0*v**2
     .     +2d0*t*(v-2d0*s-t))-s*t*(2d0*s*u*v+t*als))
     .     *dlog((s-2d0*m2)*(s-2d0*m2+dsqrt(als))/2d0/m2-1d0))
      aj30=-(pi*(als*u-s*u*v+2d0*m2*v**2)/v/als
     .     -m2*(m2*(4d0*s*u-2d0*v*(u+v))-s*u*(s+u-2d0*v))*aj25)/u**2
      aj31=pi*(u*(u**2/v+(u**2*(v-4d0*m2-2d0*s))/als
     .     +12d0*v*(-s*t*u+m2*v**2)*m2/als**2+(u**2*((s-u)**2*v
     .     +4d0*m2*(u**2+(u-s)*v)))*m2/(als*(-s*u*v+m2*(u+v)**2)))
     .     +(2d0*(3d0*s*u*(s*u-8d0*m2**2+2d0*m2*s)*v**2-9d0*m2*s*u*v**3
     .     +s*als*u**2*(v-t)+6d0*m2**2*v**3*(u+v)
     .     -2d0*s*(s-m2)*u**2*v*(2d0*(v-t)-u))*dlog((dsqrt(als)*u+s*u
     .     -2d0*m2*(u+v))**2/(4d0*m2*(-s*u*v+m2*(u+v)**2))))*m2
     .     /als**2.5)/u**5
      aj32=((u*(t+v)-v*(v-t+2d0*m2))*aj1
     .     -(u*t*(v-u)+m2*(4d0*u*t-2d0*v**2))*aj6)/((u-v)**2-4d0*m2*u)
      aj33=(-(pi*v*(4d0*m2*(4d0*m2-u)*u**2*(t-v)**2+(t-v)*u
     .     *(32d0*m2**3-4d0*m2**2*(5d0*s-24d0*u)+3d0*u**2*(s+u)
     .     -2d0*m2*u*(6d0*s+19d0*u))*v+u*(128d0*m2**3
     .     -4d0*m2**2*(33d0*s-40d0*u)-(s-9d0*u)*u*(s+u)+2d0*m2*(11d0*s**2
     .     -6d0*s*u-42d0*u**2))*v**2-(20d0*m2**3-2d0*m2**2*(8d0*s+5d0*u)
     .     +2d0*m2*(s**2+19d0*s*u-21d0*u**2)+u*(9d0*u**2-5d0*s**2
     .     -2d0*s*u))*v**3-(14d0*m2**2+s**2+4d0*s*u-3d0*u**2
     .     +6d0*m2*(u-2d0*s))*v**4+2d0*m2*v**5))/(2d0*(m2+v)**2)
     .     +aj6*(v**2*(6d0*m2**2*v**2-(s*t*u*(2d0*m2+u)))
     .     +t*u*(u-4d0*m2)*(v**3-2d0*t*u*v+(u-4d0*m2)*(t*u - v**2))))
     .     /((u - v)**2-4d0*m2*u)**2
      aj34=-pi*2d0*dlog(2d0*m2/(2d0*m2-u+dsqrt(u*(u-4d0*m2))))/t
     .     /dsqrt(u*(u-4d0*m2))
      aj35=-pi*(v*(u*(t-v)+2d0*m2*v)+2d0*m2*(t*u*(s+t)-2d0*m2*v**2)
     .     *dlog((dsqrt(u*(u-4d0*m2))-u)/2d0/m2-1d0)/dsqrt(u*(u-4d0*m2)))
     .     /t**3/u/(u-4d0*m2)
      aj36=pi*2d0*dlog((s-2d0*m2+dsqrt(als))/2d0/m2)/v/dsqrt(als)
      aj37=pi*((t/v+(s*(u-v)+2d0*m2*v)/(4d0*m2*u-(u-v)**2)-1d0)
     .     +m2*(t*u-s*v+2d0*m2*v)/((u-v)**2-4d0*m2*u)**(1.5)
     .     *dlog((v-u+2d0*m2+dsqrt((v-u)**2-4d0*m2*u))**2/4d0/m2/(m2+v)))
      aj38=pi*((2d0*(4d0*m2*u-(u-v)**2)*(m2*t**2*u**2*(u-4d0*m2)**2
     .     +t**2*u**2*(24d0*m2**2-10d0*m2*u+u**2)*v+2*t*u
     .     *(m2*(22d0*m2-5d0*s)*u-2d0*m2**2*(4d0*m2+s)
     .     +(s-9d0*m2)*u**2+u**3)*v**2-t*u*(16d0*m2**2+6d0*m2*(s-3d0*u)
     .     +u*(s+3d0*u))*v**3+(12d0*m2**3-u**2*(s+u)-4d0*m2**2*(s+4d0*u)
     .     +m2*(s**2+4d0*s*u+8*u**2))*v**4+(8d0*m2**2-4d0*m2*u+u**2)*v**5))
     .     /(v*(m2+v))+4d0*m2*dsqrt((u-v)**2-4d0*m2*u)*(t**2*(4d0*m2-u)*u**2
     .     +t*u*(8d0*m2**2+(s-u)*u+m2*(2d0*u-6d0*s))*v
     .     +t*u*(u-2d0*s)*v**2-m2*(2d0*(m2-s)+u)*v**3+m2*v**4)
     .     *dlog((2d0*m2-u+dsqrt((u-v)**2-4d0*m2*u)+v)**2/(4d0*m2*(m2+v))))
     .     /(2d0*(4d0*m2*u-(u-v)**2)**3)
      aj39=pi/t**2/(u*(u-4d0*m2))*((u*(v-t)-2d0*m2*v)
     .     *2d0*dlog(2d0*m2/(2d0*m2-u+dsqrt(u*(u-4d0*m2))))*m2
     .     /dsqrt(u*(u-4d0*m2))-(u*t*(v-u)+m2*(4d0*u*t-2d0*v**2))/v)
      aj40=pi*((u*(t**2*u**2*(u-4d0*m2)*(u-4d0*m2-2d0*v)
     .     -2d0*t*u*(8d0*m2**2+2d0*m2*(s-3d0*u)+u*(s+u))*v**2
     .     +(12d0*m2**2-4d0*m2*u+u**2)*v**4))/((u-4d0*m2)**2*v)
     .     +m2*(2d0*dsqrt(-u)*(-(t**2*(4d0*m2-u)*u**2)
     .     +2d0*t*u*(-m2*(4d0*m2+s)+(m2+s)*u)*v+3d0*m2*(2d0*m2-u)*v**3)
     .     *dlog(((2d0*m2-u)*(2d0*m2+dsqrt(4d0*m2-u)*dsqrt(-u)-u))
     .     /(2d0*m2**2)-1d0))/(4d0*m2-u)**2.5)/(t**4*u**3)
      aj41=-pi*(v*(s*v-t*u-2d0*m2*v)/(m2+v)
     .     -((v**2-4d0*u*m2)*(s-4d0*m2)-u*(s*v-2d0*m2*(3d0*v-2d0*u)))
     .     *dlog((v*(v-u+dsqrt((v-u)**2-4d0*m2*u))-2d0*m2*u)**2
     .     /4d0/m2/u**2/(m2+v))/dsqrt((v-u)**2-4d0*m2*u))/((v-u)**2-4d0*m2*u)
      aj42=2d0*pi*dlog((v*dsqrt(s-4d0*m2)+dsqrt(4d0*u*t*m2
     .     +v**2*(s-4d0*m2)))**2/(4d0*m2*t*u))
     .     /dsqrt((s-4d0*m2)*(4d0*u*t*m2+v**2*(s-4d0*m2)))
      aj43=2d0*pi*dlog((dsqrt(t*(t-4d0*m2))+2d0*m2-t)/(2d0*m2))
     .     /u/dsqrt(t*(t-4d0*m2))
      aj44=pi*(dsqrt(t*(t-4d0*m2))*(2d0*m2*v**2-t*u*(s+u))/v
     .     +2d0*m2*(t*(v-u)-2d0*m2*v)*dlog((dsqrt(t*(t-4d0*m2))
     .     +2d0*m2-t)/(2d0*m2)))/(t*(t-4d0*m2))**1.5/u**2
      aj45=pi*dlog((v*(dsqrt(als)+s)-2d0*m2*(u+v))**2
     .     /4d0/m2/(m2*(u+v)**2-s*u*v))/dsqrt(als)/u
      aj46=m2*(-pi*v*(s*u*(t+v)-2d0*m2*v*(u+v))/u/(m2*(u+v)**2-s*u*v)
     .     +(s*(u-v)+2d0*m2*v)*aj45)/u/als
      aj47=m2*(2d0*aj45*(s**2*(u-v)**2+6d0*m2**2*v**2-2d0*m2*s*(t*u
     .     +2d0*v*(v-v)))-pi*(v*(-(s**3*u**2*v*(s**2-3d0*u**2+10d0*u*v
     .     -6d0*v**2+2d0*s*(v-u)))+4d0*m2**3*(u+v)*(8d0*s**2*u**2
     .     +(u+v)**2*(4d0*s*u-3d0*v**2))+2d0*m2*s**2*u*(2d0*u**3*v-u**4
     .     +13d0*u**2*v**2-7d0*u*v**3-6d0*v**4+s**2*u*(u+5d0*v)
     .     +s*v*(8d0*u*v-5d0*u**2+v**2))-2d0*m2**2*s*(8d0*s*u**2*v*(s+2d0*v)
     .     +(u+v)*(8d0*s**2*u**2+2d0*u**3*(s+u)+6d0*u**3*v + 6d0*s*u*v**2
     .     +2d0*u**2*v**2-15d0*u*v**3-3d0*v**4))))
     .     /(u*(s*u*v-m2*(u+v)**2)**2))/(2d0*als**2*u**2)
      aj48=pi*(t*v*(4d0*m2-t)*(t*(u-v)+2d0*m2*v)
     .     -2d0*m2*dsqrt(t*(t-4d0*m2))*(t*u*(s+u)-2d0*m2*v**2)
     .     *dlog((2d0*m2-t+dsqrt(t*(t-4d0*m2)))/2d0/m2))
     .     /(t*(t-4d0*m2))**2/u**3
      aj49=pi*((t*(t-4d0*m2)*(-2d0*m2*(t-6d0*m2)*v**4+t*((-64d0*m2**3
     .     +t*((s+t)**2-(s+3d0*t)*u)+2d0*m2*(2d0*s**2+5d0*t*(2d0*u-v)))
     .     *v**2+(t*u-8d0*m2**2-2d0*m2*(3d0*s-2d0*t+u))*v**3)
     .     -t**2*(t-4d0*m2)*(2d0*(2d0*m2-u)*v**2+u**2*(s+u+v))))/v
     .     -4d0*m2*dsqrt(t*(t-4d0*m2))*(-3d0*m2*(t-2d0*m2)*v**3
     .     +t*u*(2d0*s*(t-m2)*v+(t-4d0*m2)*(t*u+2d0*m2*v)))
     .     *dlog((2d0*m2-t+dsqrt(t*(t-4d0*m2)))
     .     /(2d0*m2)))/(t**3*(t-4d0*m2)**3*u**4)
      aj50=-pi/u**3/als*(v*(s*(v-u)-2d0*m2*v)
     .     +m2*(s*u*(t+v)-2d0*m2*v*(v+u))/dsqrt(als)
     .     *dlog((v*dsqrt(als)+(s*v-2d0*m2*(u+v)))**2
     .     /(4d0*m2*(m2*(u+v)**2-s*u*v))))
      aj51=pi*((v*(-(s**3*u*(u-v)**2*v)-4d0*m2**2*s
     .     *((2d0*s-u)*u**2*(s+u)-u**2*(11d0*s+3d0*u)*v+u*(u-2d0*s)*v**2
     .     +6d0*u*v**3+v**4)+m2*s**2*(u**2*(u**2+(s+u)**2-8d0*(s+u)*v)
     .     +v**2*(2d0*u**2+4d0*u*v+v**2))+4d0*m2**3*((u+v)**2*(3d0*v**2
     .     -4d0*s*u)+4d0*s*u*(s*u-v*(u+v)))))/(als**2*(m2*(u+v)**2
     .     -s*u*v))+(2d0*(-3d0*(m2**2+(s-m2)**2)*u*v**2+3d0*m2*(s-2d0*m2)
     .     *v**3-s*(2d0*m2+s)*u**2*(v-t)+2d0*s*u*v*(s*u+(s-2d0*m2)*(v-t)))
     .     *dlog((dsqrt(als)*v+s*v-2d0*m2*(u+v))**2/(4d0*m2
     .     *(m2*(u+v)**2-s*u*v)))*m2)/als**2.5)/u**4
      endif

      sr1=4.0*(2.0*(2.0*(2.0*t-3.0*v+5.0*s-2.0*m2)*m2*pls-((3.0*(t-v)
     . +2.0*s)*pls*s+1.0))*aj7-(6.0*(2.0*m2-s)*m2*pls+pls*s**2+1.0)*
     . aj6+2.0*(pls*s**2-1.0)*aj5+4.0*((4.0*(t-v+s)+pls*s*t**2)*m2-(
     . 2.0*(pls*t**2+2.0)*m2**2+(s+t-v)**2))*aj40-(pls*s**2-1.0)*aj4-
     . 2.0*(2.0*(((t-v)*(t-2.0*v)+6.0*s**2+2.0*(3.0*t-2.0*v)*s)*pls+
     . 4.0-4.0*(2.0*s+t)*m2*pls)*m2-((2.0*((t-v)**2+s**2)+(3.0*t-2.0*
     . v)*s)*pls*s+2.0*(t-v+s)))*aj39-2.0*(2.0*m2-s)*aj37*pls-2.0*(
     . 2.0*(4.0*(tt*v-1.0-2.0*s*tt)-(2.0*t-v)*pls*t+8.0*m2*tt)*m2+4.0
     . *s**2*tt-t-4.0*(tt*v-1.0)*s+pls*s*t**2)*aj35-(2.0*(((2.0*tt*v-
     . 9.0-4.0*s*tt)*s+(2.0*tt*v-3.0)*v)*pls*s-(4.0*(tt*v-3.0-s*tt)*s
     . *tt+6.0*tt*v-5.0)-16.0*(2.0*pls*s-tt)*m2**2*tt-2.0*(((2.0*tt*v
     . -5.0)*v+t+4.0*(tt*v-2.0-3.0*s*tt)*s)*pls-4.0*(tt*v-1.0-2.0*s*
     . tt)*tt)*m2)*m2+2.0*(2.0*tt*v-3.0-2.0*s*tt)*s-(2.0*(tt*v-2.0)*v
     . +3.0*t)+(2.0*s+t)*pls*s**2)*aj34-4.0*(2.0*m2-s)**2*aj29+4.0*(
     . 2.0*(4.0*m2-3.0*s)*m2*pls+pls*s**2-1.0)*aj28*s+2.0*((4.0*(((
     . 6.0*s**2-v**2)*tt-2.0*(tt*v-2.0)*s)*pls+2.0*(tt*v-1.0-2.0*s*tt
     . )*tt-4.0*(2.0*pls*s-tt)*m2*tt)*m2-(2.0*(2.0*(tt*v-3.0-s*tt)*s*
     . tt+3.0*tt*v-1.0)-(2.0*(tt*v-7.0-2.0*s*tt)*s+(2.0*tt*v+3.0)*v)*
     . pls*s))*m2-((tt*v-1.0)*v+t-2.0*(tt*v-1.0-s*tt)*s-(t-v+2.0*s)*
     . pls*s**2))*aj23+2.0*(8.0*(tt*v-1.0-2.0*s*tt+2.0*m2*tt)*m2+4.0*
     . s**2*tt-t-4.0*(tt*v-1.0)*s)*aj22+(2.0*(4.0*(((6.0*s**2-v**2)*
     . tt-2.0*(tt*v-2.0)*s)*pls+2.0*(tt*v-1.0-2.0*s*tt)*tt-4.0*(2.0*
     . pls*s-tt)*m2*tt)*m2+((2.0*tt*v-9.0-4.0*s*tt)*s+(2.0*tt*v+1.0)*
     . v)*pls*s-(4.0*(tt*v-3.0-s*tt)*s*tt+6.0*tt*v-1.0))*m2+2.0*(2.0*
     . tt*v-1.0-2.0*s*tt)*s-(2.0*(tt*v-1.0)*v+t)+(t-2.0*v+2.0*s)*pls*
     . s**2)*aj21+4.0*(pls*s*t-1.0-2.0*m2*pls*t)*aj12-2.0*(4.0*m2**2*
     . pls-2.0*m2*pls*s+1.0)*aj11)
      sr2=-4.0*((2.0*(tt*v-1.0+t*vt-2.0*(tt-vt-6.0*s*tt*vt)*s-2.0*(t*
     . vt+tt*v-4.0*(tt-vt)*s)*pls*s**2-4.0*(5.0*(tt-vt)*pls*s+6.0*tt*
     . vt)*m2*s)*m2-(tt*v-1.0+t*vt-2.0*(tt-vt-2.0*s*tt*vt)*s-(tt*v-
     . 1.0+t*vt-2.0*(tt-vt)*s)*pls*s**2)*s+16.0*((tt*v+1.0+t*vt+2.0*(
     . tt-vt)*s)*pls+2.0*tt*vt)*m2**3)*aj36+(2.0*(((2.0*tt-11.0*vv-
     . 8.0*s*tt*vv)*s-(3.0*t*vv-4.0*tt*v-2.0))*pls*s-(tt+6.0*vv+6.0*s
     . *tt*vv))*m2-(((tt-4.0*vv-2.0*s*tt*vv)*s-2.0*(t*vv-1.0))*(pls*s
     . **2+1.0)+32.0*(s*tt+1.0)*m2**3*pls*vv)-8.0*((2.0*tt*v+1.0-2.0*
     . t*vv-(5.0*s*tt+4.0)*s*vv)*pls-2.0*tt*vv)*m2**2)*aj6+(8.0*(2.0*
     . (tt+vv)-5.0*s*tt*vv+4.0*m2*tt*vv)*m2**2*pls*s-((2.0*s*vv-3.0)*
     . s*tt+tt*v-1.0)*(pls*s**2-1.0)+2.0*((8.0*s**2*vv-7.0*s+2.0*v)*
     . pls*s-(2.0*s*vv+1.0))*m2*tt)*aj4+(2.0*(4.0*(((2.0*(t*vv+1.0)-
     . 5.0*(vt-vv)*s)*s-(t**2*vv-2.0*t+v))*pls-2.0*(tt+2.0*vv-3.0*(vt
     . -vv)*s*tt)-2.0*(((vt-vv)*t+1.0-2.0*(vt-vv)*s)*pls+2.0*(vt-vv)*
     . tt)*m2)*m2-(t*vt-7.0*t*vv+4.0*tt*v+2.0+12.0*(vt-vv)*s**2*tt+
     . 2.0*(vt-9.0*vv-2.0*tt)*s+(3.0*t-2.0*v-8.0*(vt-vv)*s**2-((2.0*
     . vt-7.0*vv)*t-2.0)*s)*pls*s))*m2+2.0*(2.0*(vt-vv)*s*tt+vt-5.0*
     . vv)*s**2-(3.0*t**2*vv-4.0*t+2.0*v)+((vt-9.0*vv)*t+6.0)*s-(s*vt
     . -s*vv-t*vv)*(2.0*s+t)*pls*s**2)*aj34-(((pls*s**2+1.0)*(s*tt+
     . 1.0)-16.0*m2**3*pls*tt)*vv+8.0*(tt+vv+2.0*s*tt*vv)*m2**2*pls-
     . 2.0*((2.0*tt+3.0*vv+4.0*s*tt*vv)*pls*s+tt*vv)*m2)*aj32-(2.0*(
     . 4.0*(2.0*((tt*v+1.0+t*vt+2.0*(tt-vt)*s)*pls+2.0*tt*vt)*m2+(2.0
     . *(tt*v-3.0-t*vt)*s-(5.0*(tt-vt)*s**2-tt*v**2))*pls+2.0*(2.0*tt
     . +vt-3.0*s*tt*vt))*m2+7.0*tt*v-1.0+t*vt-6.0*(3.0*tt+vt-2.0*s*tt
     . *vt)*s-((7.0*tt*v-18.0-4.0*t*vt)*s-(8.0*(tt-vt)*s**2-3.0*v))*
     . pls*s)*m2+(3.0*tt*v-4.0)*v+2.0*t+2.0*(5.0*tt+vt-2.0*s*tt*vt)*s
     . **2-(9.0*tt*v-5.0+t*vt)*s-((tt*v-4.0)*v+2.0*t+2.0*(tt-vt)*s**2
     . -(3.0*tt*v-7.0-t*vt)*s)*pls*s**2)*aj24+((2.0*(vt+vv-2.0*(vt-vv
     . )*s*tt)*s-((vt-vv)*t+2.0)-((vt+vv)*t-2.0*(vt-vv)*s)*pls*s**2)*
     . s+16.0*((t*vt+1.0-2.0*(vt-vv)*s)*pls+2.0*(vt-vv)*tt)*m2**3-8.0
     . *(((2.0*vt+vv)*t+2.0-5.0*(vt-vv)*s)*pls*s-2.0*(tt+vt-3.0*(vt-
     . vv)*s*tt))*m2**2+2.0*((vt-vv)*t+2.0+12.0*(vt-vv)*s**2*tt-2.0*(
     . 3.0*vt+vv+2.0*tt)*s+((4.0*vt+vv)*t+2.0-8.0*(vt-vv)*s)*pls*s**2
     . )*m2)*aj21-((s*vv-1.0-2.0*m2*vv)*(pls*s**2-1.0)+8.0*m2**2*pls*
     . s*vv)*aj20*tt-2.0*(4.0*(tt-2.0*vt)*m2**2*pls+tt-2.0*((tt-2.0*
     . vt)*pls*s+2.0*tt*vt)*m2)*aj16-2.0*(4.0*(2.0*vt+vv)*m2**2*pls+
     . vv-2.0*((2.0*vt+vv)*pls*s-2.0*tt*vt)*m2)*aj15-2.0*(2.0*tt*v-
     . 1.0-(3.0*tt+vt)*s+pls*s**2+8.0*(tt+vt)*m2**3*pls+4.0*((2.0*tt*
     . v-1.0-t*vt-2.0*s*tt)*pls-2.0*tt*vt)*m2**2-((3.0*tt*v-1.0-2.0*t
     . *vt-2.0*s*tt)*pls*s-2.0*(2.0*s*vt+1.0)*tt)*m2)*aj13-2.0*(2.0*(
     . 2.0*(((vt-vv)*t+3.0)*pls+2.0*(vt-vv)*tt-2.0*(vt-vv)*m2*pls)*m2
     . -(2.0*((vt-vv)*s*tt-2.0*vv)+(t*vt+2.0)*pls*s))*m2-(t*vv-1.0-(
     . vt-3.0*vv)*s-pls*s**2))*aj11-2.0*(4.0*(tt+vv+2.0*s*tt*vv-2.0*
     . m2*tt*vv)*m2**2*pls+(s*vv-1.0)*tt+(tt+vv)*pls*s**2-((tt+3.0*vv
     . +2.0*s*tt*vv)*pls*s+2.0*tt*vv)*m2)*aj1)
      sr3=-4.0*(2.0*((((2.0*m2*tt-s*vv)*m2*pls+(pls*s**2-1.0)*vv)*aj2
     . -2.0*(pls*s-tt-2.0*m2*pls)*aj19+(4.0*m2**2*pls-2.0*m2*pls*s+
     . 1.0)*aj18*vv)*tt-((4.0*tt*v-1.0-2.0*s*tt)*pls*s-2.0*(2.0*tt*v-
     . 1.0-2.0*s*tt)*tt+2.0*((4.0*s-3.0*v)*pls-4.0*tt)*m2*tt)*aj17+(
     . 2.0*tt-vv-2.0*s*tt*vv+pls*s**2*vv-(2.0*(2.0*tt-3.0*vv-2.0*s*tt
     . *vv)*tt-(tt-3.0*vv-2.0*s*tt*vv)*pls*s)*m2-2.0*((tt*v-3.0-4.0*s
     . *vv)*pls+4.0*tt*vv)*m2**2*tt)*aj16+(2.0*((6.0*s**2+v**2-(tt*v+
     . 4.0)*s*v)*pls-4.0*(s-v)*tt-2.0*((4.0*s-tt*v**2)*pls-2.0*tt)*m2
     . )*m2*tt-(((2.0*tt*v-1.0)*v+2.0*s**2*tt-(3.0*tt*v-1.0)*s)*pls*s
     . +2.0*(2.0*tt*v-1.0-s*tt)*s*tt-(2.0*tt**2*v**2-2.0*tt*v+1.0)))*
     . aj14)-(2.0*((2.0*(6.0*tt-7.0*vv-2.0*s*tt*vv)*s-(tt*v+3.0)*(tt*
     . v-1.0))*pls*s-(4.0*(2.0*tt-3.0*vv-s*tt*vv)*s*tt-(3.0*tt**2*v-
     . 6.0*tt+2.0*vv))-8.0*((4.0*s*vv-tt*v)*pls-2.0*tt*vv)*m2**2*tt+
     . 2.0*((((2.0*tt*v-3.0)*v+12.0*s**2*vv)*tt-2.0*(tt**2*v+6.0*tt-
     . 4.0*vv)*s)*pls+4.0*(2.0*tt-vv-2.0*s*tt*vv)*tt)*m2)*m2-(3.0*tt*
     . v-4.0+2.0*t*vv-2.0*(3.0*tt-2.0*vv-2.0*s*tt*vv)*s-(tt*v-4.0+2.0
     . *t*vv-2.0*(tt-2.0*vv)*s)*pls*s**2))*aj13-((pls*s**2-1.0)*vv+
     . 2.0*(pls*s-vv)*m2*tt)*aj10*tt+(4.0*((2.0*(3.0*s**2*vv**2+tt*v)
     . -(tt+8.0*vv)*s)*pls-4.0*(s*vv-1.0)*tt*vv)*m2**2*tt-(8.0*((4.0*
     . s*vv**2-tt)*pls-2.0*tt*vv**2)*m2**3*tt-(tt-vv-2.0*s*tt*vv)*(
     . pls*s**2-1.0))+2.0*(((9.0*tt-vv-2.0*s*tt*vv)*s*vv-(tt*v+3.0)*
     . tt)*pls*s-(2.0*(2.0*tt-vv-s*tt*vv)*s*tt*vv-(3.0*tt**2-tt*vv+vv
     . **2)))*m2)*aj1)
      sr4=-4.0*(((2.0*(2.0*t*vv-3.0+s*vv)*s+2.0*t**2*vv-5.0*t+4.0*v)*
     . (pls*s**2+1.0)-16.0*(2.0*t*vv-3.0+6.0*s*vv)*m2**3*pls-2.0*(((
     . 11.0*s*vv+18.0*t*vv-21.0)*s+9.0*t**2*vv-20.0*t+12.0*v)*pls*s+
     . 9.0*t*vv-10.0+7.0*s*vv)*m2+4.0*(((22.0*t*vv-25.0+20.0*s*vv)*s+
     . 4.0*t**2*vv-18.0*t+13.0*v)*pls+8.0*vv)*m2**2)*aj6-2.0*(((11.0*
     . t-8.0*v+5.0*s)*s+2.0*(4.0*t**2-6.0*t*v+3.0*v**2))*pls*s+4.0*(t
     . -v+s)+16.0*(2.0*t-v+4.0*s)*m2**2*pls-2.0*(((26.0*t-17.0*v+18.0
     . *s)*s+2.0*(4.0*t**2-6.0*t*v+3.0*v**2))*pls+8.0)*m2)*aj7-2.0*((
     . 2.0*(t-v)-3.0*s**2*vv+2.0*(t*vv+1.0)*s)*pls*s-(2.0*(t*vv-1.0)-
     . s*vv)+4.0*(5.0*(s*vv-1.0)*pls*s+vv-8.0*m2*pls*s*vv)*m2)*aj4*m2
     . -2.0*(2.0*((2.0*(7.0*t-6.0*v+4.0*s)*s**2+(3.0*t**2-3.0*t*v+2.0
     . *v**2)*(t-v)+(13.0*t**2-17.0*t*v+8.0*v**2)*s)*pls+8.0*(t-v+s))
     . *m2-(8.0*(((4.0*(t-v)+5.0*s)*s+2.0*t**2-2.0*t*v+v**2)*pls+3.0)
     . *m2**2+((2.0*s**2+3.0*s*t-2.0*s*v+4.0*t**2-4.0*t*v+2.0*v**2)*
     . pls*s+2.0*(s+t-v))*(s+t-v)-32.0*m2**3*pls*s))*aj39+4.0*(2.0*m2
     . -s)*aj38*pls-4.0*((7.0*t-6.0*v+9.0*s-8.0*m2)*m2*pls-((3.0*(t-v
     . )+2.0*s)*pls*s+1.0))*aj37-(((2.0*(2.0*(t*vv-1.0)+s*vv)*s+3.0*(
     . t**2*vv-2.0*t+v))*s+t**3*vv-3.0*t**2+3.0*t*v-v**2)*(pls*s**2+
     . 1.0)+64.0*(tt+vv)*m2**4*pls*s+16.0*(((2.0*(tt*v+1.0-2.0*t*vv)-
     . (5.0*tt+7.0*vv)*s)*s-((t*vv-3.0)*t+(tt*v+1.0)*v))*pls-3.0*(tt+
     . vv))*m2**3-4.0*((6.0*(tt*v+2.0-3.0*t*vv)*s**2-(8.0*s**3*tt+
     . 18.0*s**3*vv+2.0*t**3*vv-10.0*t**2+13.0*t*v-5.0*v**2)-(3.0*(
     . 2.0*t*vv-7.0)*t+2.0*(2.0*tt*v+3.0)*v)*s)*pls+2.0*(tt*v+3.0-7.0
     . *t*vv-(4.0*tt+7.0*vv)*s))*m2**2-2.0*((2.0*t**3*vv-10.0*t**2+
     . 11.0*t*v-3.0*v**2+2.0*(tt+5.0*vv)*s**3-(2.0*tt*v+13.0-15.0*t*
     . vv)*s**2+((7.0*t*vv-18.0)*t+2.0*(tt*v+3.0)*v)*s)*pls*s+7.0*t**
     . 2*vv-10.0*t+3.0*v+2.0*(tt+5.0*vv)*s**2-(2.0*(tt*v+5.0)-15.0*t*
     . vv)*s)*m2)*aj34-(8.0*(t*vv-4.0+2.0*s*vv-2.0*m2*vv)*m2**2*pls+(
     . t*vv-2.0+s*vv)*(pls*s**2+1.0)-2.0*(3.0*(t*vv-2.0+s*vv)*pls*s+
     . 2.0*vv)*m2)*aj32+(8.0*((4.0*s**3*tt+v**2-3.0*(tt*v-1.0)*s**2+(
     . (2.0*tt*v-1.0)*v-6.0*t)*s)*pls-(tt*v-3.0-4.0*s*tt))*m2**2-(
     . 16.0*(((5.0*s**2+v**2)*tt-2.0*(tt*v+1.0)*s)*pls+3.0*tt)*m2**3-
     . ((pls*s**2+1.0)*(2.0*s**2-2.0*s*v+v**2)+64.0*m2**4*pls*s*tt))+
     . 2.0*((2.0*(tt*v-4.0-s*tt)*s**2+4.0*t**2-3.0*t*v-3.0*v**2-2.0*(
     . (tt*v-3.0)*v-2.0*t)*s)*pls*s-(2.0*s**2*tt-5.0*v-2.0*(tt*v-4.0)
     . *s))*m2)*aj23-(16.0*((((5.0*tt+7.0*vv)*s-2.0*tt*v)*s+(tt*v+1.0
     . )*v)*pls+3.0*(tt+vv))*m2**3-((pls*s**2+1.0)*(2.0*s**2*vv-2.0*s
     . +v)+64.0*(tt+vv)*m2**4*pls)*s-2.0*((3.0*t**2*vv-5.0*t-2.0*tt*v
     . **2-2.0*(tt+5.0*vv)*s**2+(2.0*tt*v+5.0+t*vv)*s)*pls*s**2-(t**2
     . *vv+v+2.0*(tt+5.0*vv)*s**2-(2.0*(tt*v+1.0)+t*vv)*s))*m2+8.0*((
     . (3.0*tt*v+1.0+t*vv-(4.0*tt+9.0*vv)*s)*s+(t*vv-1.0)*t-(2.0*tt*v
     . +1.0)*v)*pls*s-(t*vv-tt*v+(4.0*tt+7.0*vv)*s))*m2**2)*aj21+2.0*
     . ((3.0*s*vv-2.0)*pls*s-vv-4.0*m2*pls*s*vv)*aj20*m2-(4.0*(2.0*m2
     . -s)*(t*vv-2.0)*m2*pls+(pls*s**2+1.0)*(t*vv-1.0))*aj15+(4.0*((
     . 2.0*(4.0*(t*vv-1.0)+3.0*s*vv)*s+2.0*t**2*vv-7.0*t+7.0*v)*pls+
     . 6.0*vv)*m2**2+(t**2*vv-3.0*t+2.0*v+(3.0*t*vv-2.0)*s)*(pls*s**2
     . +1.0)-16.0*(2.0*s+t)*m2**3*pls*vv-2.0*(((9.0*t*vv-8.0+2.0*s*vv
     . )*s+4.0*t**2*vv-11.0*t+9.0*v)*pls*s+7.0*t*vv-8.0+2.0*s*vv)*m2)
     . *aj11-(8.0*(t*vv-3.0+3.0*s*vv-2.0*m2*vv)*m2**2*pls+(t*vv-2.0+s
     . *vv)*(pls*s**2+1.0)-2.0*((5.0*t*vv-7.0+4.0*s*vv)*pls*s+4.0*vv)
     . *m2)*aj1)*uu
      sr5=-4.0*(((16.0*(((2.0*(t*uu+3.0)-7.0*s*uu)*s-2.0*(t**2*uu+t+v
     . ))*pls-3.0*uu)*m2**3-((2.0*(t*uu+1.0-s*uu)*s-(t**2*uu+t+v))*(
     . pls*s**2+1.0)-64.0*m2**4*pls*uu)*s-8.0*(((5.0*t*uu+11.0-9.0*s*
     . uu)*s-(3.0*t**2*uu+5.0*t+4.0*v))*pls*s-7.0*(s*uu-1.0))*m2**2+
     . 2.0*((5.0*t*uu+12.0-10.0*s*uu)*s-(2.0*t**2*uu+2.0*t+v)-(2.0*(
     . 5.0*s*uu-4.0*t*uu-6.0)*s+4.0*t**2*uu+6.0*t+5.0*v)*pls*s**2)*m2
     . -4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls*s**3+s)*(2.0*
     . m2-s+v)*dz1*m2)*ds+2.0*(2.0*((t*uu+3.0)*t+3.0*s**2*uu-(3.0*t*
     . uu+2.0)*s)*pls*s-3.0*(t*uu+1.0-2.0*s*uu))*m2+(2.0*(t*uu+1.0-s*
     . uu)*s-(t**2*uu+t+v))*(pls*s**2+1.0)-8.0*(2.0*(s-t)*s*uu+t**2*
     . uu+t+v)*m2**2*pls+((pls*s**2+1.0)*(2.0*s**2-2.0*s*v+v**2)-32.0
     . *m2**3*pls*s-2.0*(2.0*(4.0*s**2-3.0*s*v+v**2)*pls*s+8.0*s-3.0*
     . v)*m2+8.0*((5.0*s**2-2.0*s*v+v**2)*pls+3.0)*m2**2)*dz1+4.0*(
     . 8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls*s**3+s)*(2.0*m2-s)*
     . ds**2*m2)*aj45+(((pls*s**2+1.0)*(s+t)-48.0*m2**3*pls+4.0*(7.0*
     . t-4.0*v+21.0*s)*m2**2*pls-2.0*((11.0*t-4.0*v+11.0*s)*pls*s+4.0
     . )*m2)*dz1-2.0*(2.0*(8.0*m2-5.0*s)*m2*pls+pls*s**2+1.0))*aj8-((
     . pls*s**2+1.0)*(s+t)-48.0*m2**3*pls+4.0*(7.0*t-4.0*v+13.0*s)*m2
     . **2*pls-2.0*((7.0*t-4.0*v+7.0*s)*pls*s+4.0)*m2)*aj6*dz1-16.0*
     . aj5*m2*pls*s+4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls*s
     . **3+s)*(2.0*m2-s)*aj46*ds-(2.0*(4.0*((3.0*(3.0*s*vt+2.0)*s+t**
     . 2*vt+3.0*t+v)*pls*s+t*vt+1.0+7.0*s*vt)*m2**2+((pls*s**2+1.0)*s
     . **2+32.0*m2**4*pls)*s*vt-8.0*(((7.0*s*vt+4.0)*s+t**2*vt+t+v)*
     . pls+3.0*vt)*m2**3+((t*vt+1.0-10.0*s*vt)*s-(t**2*vt+t+v)-2.0*(
     . 5.0*s**2*vt+2.0*s+2.0*t)*pls*s**2)*m2)*dz1+16.0*(6.0*t*vt-1.0+
     . 6.0*s*vt-8.0*m2*vt)*m2**3*pls-((2.0*s**2+t**2)*uu-(2.0*uu-vt)*
     . s*t)*(pls*s**2+1.0)+4.0*((t-2.0*v-2.0*(uu+vt)*t**2+((4.0*uu-
     . 15.0*vt)*t-4.0*(uu+2.0*vt)*s)*s)*pls-2.0*vt)*m2**2-2.0*((2.0*(
     . 3.0*(uu-vt)*t-1.0-(3.0*uu+vt)*s)*s+2.0*t-v-2.0*(uu+vt)*t**2)*
     . pls*s+(3.0*uu+vt)*t-1.0-6.0*s*uu)*m2+(((pls*s**2+1.0)*(2.0*s**
     . 2-2.0*s*t+t**2)+64.0*m2**4*pls)*s-16.0*((7.0*s**2-2.0*s*t+2.0*
     . t**2)*pls+3.0)*m2**3+8.0*((9.0*s**2-5.0*s*t+3.0*t**2)*pls+7.0)
     . *m2**2*s-2.0*(10.0*s**2-5.0*s*t+2.0*t**2+2.0*(5.0*s**2-4.0*s*t
     . +2.0*t**2)*pls*s**2)*m2)*(uu-vt)*ds)*aj43-((((t-v+3.0*s)*(t-v)
     . -2.0*(s*vt-2.0)*s**2)*(pls*s**2+1.0)-64.0*m2**4*pls*s*vt+16.0*
     . (((7.0*s*vt-8.0)*s+t**2*vt-2.0*t+2.0*v)*pls+3.0*vt)*m2**3+2.0*
     . (((10.0*(s*vt-2.0)*s-(14.0*t-13.0*v))*s-(6.0*t**2-9.0*t*v+4.0*
     . v**2))*pls*s-((t*vt+20.0-10.0*s*vt)*s-(t**2*vt-8.0*t+9.0*v)))*
     . m2-4.0*((2.0*(9.0*s*vt-16.0)*s**2-(5.0*t**2-7.0*t*v+4.0*v**2)+
     . (2.0*t**2*vt-17.0*t+16.0*v)*s)*pls+2.0*(t*vt-10.0+7.0*s*vt))*
     . m2**2)*dz1-(2.0*(((6.0*t*vt-1.0+2.0*s*vt)*s+2.0*t**2*vt-6.0*t+
     . 9.0*v)*pls*s-t*vt)*m2-(((t*vt-1.0)*s-(t-2.0*v))*(pls*s**2+1.0)
     . +128.0*m2**4*pls*vt-16.0*(6.0*t*vt+5.0+6.0*s*vt)*m2**3*pls)-
     . 4.0*(((15.0*t*vt+7.0+8.0*s*vt)*s+2.0*t**2*vt-3.0*t+7.0*v)*pls+
     . 2.0*vt)*m2**2)+4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls
     . *s**3+s)*(2.0*m2-s)*ds**2*m2+(16.0*(((2.0*(t*vt+3.0)-7.0*s*vt)
     . *s-2.0*(t**2*vt+t+v))*pls-3.0*vt)*m2**3-((2.0*(t*vt+1.0-s*vt)*
     . s-(t**2*vt+t+v))*(pls*s**2+1.0)-64.0*m2**4*pls*vt)*s-8.0*(((
     . 5.0*t*vt+11.0-9.0*s*vt)*s-(3.0*t**2*vt+5.0*t+4.0*v))*pls*s-7.0
     . *(s*vt-1.0))*m2**2+2.0*((5.0*t*vt+12.0-10.0*s*vt)*s-(2.0*t**2*
     . vt+2.0*t+v)-(2.0*(5.0*s*vt-4.0*t*vt-6.0)*s+4.0*t**2*vt+6.0*t+
     . 5.0*v)*pls*s**2)*m2-4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*
     . m2+pls*s**3+s)*(2.0*m2-s+v)*dz1*m2)*ds)*aj42-4.0*(2.0*m2-s)*
     . aj41*dz1*m2*pls+8.0*(s+t-4.0*m2)*aj4*dz1*m2*pls*s+2.0*(4.0*((
     . 3.0*(3.0*s*vt+2.0)*s+t**2*vt+3.0*t+v)*pls*s+t*vt+1.0+7.0*s*vt)
     . *m2**2+((pls*s**2+1.0)*s**2+32.0*m2**4*pls)*s*vt-8.0*(((7.0*s*
     . vt+4.0)*s+t**2*vt+t+v)*pls+3.0*vt)*m2**3+((t*vt+1.0-10.0*s*vt)
     . *s-(t**2*vt+t+v)-2.0*(5.0*s**2*vt+2.0*s+2.0*t)*pls*s**2)*m2)*
     . aj36*dz1+(((t-v+3.0*s)*(t-v)-2.0*(s*vt-2.0)*s**2)*(pls*s**2+
     . 1.0)+64.0*(tt-vt)*m2**4*pls*s+16.0*(((2.0*(tt*v-4.0)-(5.0*tt-
     . 7.0*vt)*s)*s+t**2*vt-2.0*t+2.0*v)*pls-3.0*(tt-vt))*m2**3+2.0*(
     . ((2.0*(tt*v-10.0-(tt-5.0*vt)*s)*s-(14.0*t-13.0*v))*s-(6.0*t**2
     . -9.0*t*v+4.0*v**2))*pls*s+t**2*vt-8.0*t+9.0*v-2.0*(tt-5.0*vt)*
     . s**2+(2.0*(tt*v-10.0)-t*vt)*s)*m2-4.0*((2.0*(3.0*tt*v-16.0-(
     . 4.0*tt-9.0*vt)*s)*s**2-(5.0*t**2-7.0*t*v+4.0*v**2)+(2.0*t**2*
     . vt-17.0*t+16.0*v)*s)*pls+2.0*(3.0*tt*v-10.0+t*vt-(4.0*tt-7.0*
     . vt)*s))*m2**2)*aj34*dz1+4.0*(2.0*m2-s)*aj32*dz1*m2*pls-4.0*(
     . 8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls*s**3+s)*(2.0*m2-s)*
     . aj30*ds+4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls*s**3+s
     . )*(2.0*m2-s)*aj28*ds-4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*
     . m2+pls*s**3+s)*(2.0*m2-s)*aj26*ds-((16.0*(((2.0*(t*uu+3.0)-7.0
     . *s*uu)*s-2.0*(t**2*uu+t+v))*pls-3.0*uu)*m2**3-((2.0*(t*uu+1.0-
     . s*uu)*s-(t**2*uu+t+v))*(pls*s**2+1.0)-64.0*m2**4*pls*uu)*s-8.0
     . *(((5.0*t*uu+11.0-9.0*s*uu)*s-(3.0*t**2*uu+5.0*t+4.0*v))*pls*s
     . -7.0*(s*uu-1.0))*m2**2+2.0*((5.0*t*uu+12.0-10.0*s*uu)*s-(2.0*t
     . **2*uu+2.0*t+v)-(2.0*(5.0*s*uu-4.0*t*uu-6.0)*s+4.0*t**2*uu+6.0
     . *t+5.0*v)*pls*s**2)*m2-4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-
     . 6.0*m2+pls*s**3+s)*(2.0*m2-s+v)*dz1*m2)*ds+2.0*(2.0*((t*uu+3.0
     . )*t+3.0*s**2*uu-(3.0*t*uu+2.0)*s)*pls*s-3.0*(t*uu+1.0-2.0*s*uu
     . ))*m2+(2.0*(t*uu+1.0-s*uu)*s-(t**2*uu+t+v))*(pls*s**2+1.0)-8.0
     . *(2.0*(s-t)*s*uu+t**2*uu+t+v)*m2**2*pls+((pls*s**2+1.0)*(2.0*s
     . **2-2.0*s*v+v**2)-32.0*m2**3*pls*s-2.0*(2.0*(4.0*s**2-3.0*s*v+
     . v**2)*pls*s+8.0*s-3.0*v)*m2+8.0*((5.0*s**2-2.0*s*v+v**2)*pls+
     . 3.0)*m2**2)*dz1+4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+
     . pls*s**3+s)*(2.0*m2-s)*ds**2*m2)*aj25+((((pls*s**2+1.0)*(2.0*s
     . **2-2.0*s*t+t**2)+64.0*m2**4*pls)*s-16.0*((7.0*s**2-2.0*s*t+
     . 2.0*t**2)*pls+3.0)*m2**3+8.0*((9.0*s**2-5.0*s*t+3.0*t**2)*pls+
     . 7.0)*m2**2*s-2.0*(10.0*s**2-5.0*s*t+2.0*t**2+2.0*(5.0*s**2-4.0
     . *s*t+2.0*t**2)*pls*s**2)*m2)*(uu-vt)*ds-(2.0*((2.0*(3.0*t*uu+
     . 1.0-3.0*s*uu)*s-(2.0*t**2*uu+6.0*t-v))*pls*s+3.0*(t*uu+1.0-2.0
     . *s*uu))*m2-((2.0*(t*uu+1.0-s*uu)*s-(t**2*uu+t+v))*(pls*s**2+
     . 1.0)-8.0*(2.0*(s-t)*s*uu+t**2*uu+t+v)*m2**2*pls)))*aj24-(8.0*(
     . (4.0*s**3*tt+v**2-3.0*(tt*v-1.0)*s**2-(6.0*t+v)*s)*pls-(3.0*(
     . tt*v-1.0)-4.0*s*tt))*m2**2+16.0*((2.0*(tt*v+1.0)-5.0*s*tt)*pls
     . *s-3.0*tt)*m2**3+(pls*s**2+1.0)*(2.0*s**2-2.0*s*v+v**2)+64.0*
     . m2**4*pls*s*tt+2.0*(((2.0*(tt*v-4.0-s*tt)*s+4.0*t+5.0*v)*s+4.0
     . *t**2-t*v-2.0*v**2)*pls*s-(2.0*s**2*tt-3.0*v-2.0*(tt*v-4.0)*s)
     . )*m2)*aj23*dz1+((16.0*(((2.0*(t*vt+3.0)-7.0*s*vt)*s-2.0*(t**2*
     . vt+t+v))*pls-3.0*vt)*m2**3-((2.0*(t*vt+1.0-s*vt)*s-(t**2*vt+t+
     . v))*(pls*s**2+1.0)-64.0*m2**4*pls*vt)*s-8.0*(((5.0*t*vt+11.0-
     . 9.0*s*vt)*s-(3.0*t**2*vt+5.0*t+4.0*v))*pls*s-7.0*(s*vt-1.0))*
     . m2**2+2.0*((5.0*t*vt+12.0-10.0*s*vt)*s-(2.0*t**2*vt+2.0*t+v)-(
     . 2.0*(5.0*s*vt-4.0*t*vt-6.0)*s+4.0*t**2*vt+6.0*t+5.0*v)*pls*s**
     . 2)*m2-4.0*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls*s**3+s)*
     . (2.0*m2-s+v)*dz1*m2)*ds-4.0*((2.0*m2-s+v)*dz1*tt-(2.0*m2-s)*ds
     . **2)*(8.0*m2**2*pls*s-6.0*m2*pls*s**2-6.0*m2+pls*s**3+s)*m2)*
     . aj21+4.0*(2.0*m2-s)*aj16*m2*pls*vt-4.0*(2.0*m2-s)*aj15*m2*pls*
     . vt-(32.0*m2**3*pls*vt-8.0*m2**2*pls*s*vt-16.0*m2**2*pls*t*vt+
     . 4.0*m2**2*pls+8.0*m2*pls*s*t*vt-2.0*m2*pls*s-4.0*m2*vt+pls*s**
     . 2+1.0)*aj13+(32.0*m2**3*pls*vt-8.0*m2**2*pls*s*vt-16.0*m2**2*
     . pls*t*vt-20.0*m2**2*pls+8.0*m2*pls*s*t*vt+6.0*m2*pls*s-4.0*m2*
     . vt-pls*s**2-1.0)*aj11)
      sr6=4.0*((((2.0*(s*tt+2.0)*s*vv+3.0*(t*vv-1.0))*s+t**2*vv-2.0*t
     . +v)*(pls*s**2+1.0)+64.0*m2**4*pls*s*tt*vv-16.0*((tt*v-1.0+2.0*
     . t*vv+(7.0*s*vv+2.0)*s*tt)*pls+3.0*tt*vv)*m2**3+4.0*(((2.0*(3.0
     . *tt+4.0*vv+9.0*s*tt*vv)*s-(2.0*tt*v-1.0))*s+(2.0*t*vv-1.0)*t+(
     . 2.0*tt*v-1.0)*v)*pls+2.0*(2.0*tt+3.0*vv+7.0*s*tt*vv))*m2**2+
     . 2.0*(((2.0*tt*v+3.0-3.0*t*vv-2.0*(tt+6.0*vv+5.0*s*tt*vv)*s)*s+
     . (t*vv-1.0)*t-(2.0*tt*v-1.0)*v)*pls*s+tt*v+3.0-5.0*t*vv-(3.0*tt
     . +13.0*vv+10.0*s*tt*vv)*s)*m2)*aj6-2.0*(((pls*s**2+1.0)*s**2+
     . 32.0*m2**4*pls)*s*tt-8.0*((7.0*s**2*tt-4.0*s+tt*v**2)*pls+3.0*
     . tt)*m2**3-(2.0*((5.0*s*tt-2.0)*s-2.0*(t-v))*pls*s**2+(10.0*s**
     . 2+s*v+v**2)*tt)*m2+4.0*(((tt*v+2.0)*v-2.0*t+3.0*(3.0*s*tt-2.0)
     . *s)*pls*s+(7.0*s-v)*tt)*m2**2)*aj36+(((pls*s**2+1.0)*(s*tt+1.0
     . )-32.0*m2**3*pls*tt)*vv+8.0*(tt+vv+3.0*s*tt*vv)*m2**2*pls-4.0*
     . ((tt+vv+2.0*s*tt*vv)*pls*s+tt*vv)*m2)*aj33+(((tt-4.0*vv-2.0*s*
     . tt*vv)*s-2.0*(t*vv-1.0))*(pls*s**2+1.0)-16.0*(tt-4.0*vv-2.0*s*
     . tt*vv)*m2**3*pls+4.0*(4.0*tt*v-1.0-4.0*t*vv+(tt-18.0*vv-10.0*s
     . *tt*vv)*s)*m2**2*pls-2.0*(((2.0*tt-15.0*vv-8.0*s*tt*vv)*s-(5.0
     . *t*vv-4.0*tt*v-2.0))*pls*s-(2.0*tt+7.0*vv+5.0*s*tt*vv))*m2)*
     . aj32+(8.0*(4.0*s*tt+1.0-4.0*m2*tt)*m2**2*pls+(pls*s**2+1.0)*(s
     . *tt+1.0)-2.0*(5.0*pls*s**2+1.0)*m2*tt)*aj2*vv+(4.0*(2.0*m2-s)*
     . (2.0*tt-vv)*m2*pls+(pls*s**2+1.0)*(tt-vv))*aj18-(4.0*((2.0*(
     . 8.0*tt-vv-3.0*s*tt*vv)*s-(7.0*tt*v-3.0))*pls-6.0*tt*vv)*m2**2-
     . ((2.0*(tt*v-1.0)-(3.0*tt-2.0*vv)*s)*(pls*s**2+1.0)-32.0*(s*vv-
     . 1.0)*m2**3*pls*tt)+2.0*((9.0*tt*v-4.0-t*vv-(13.0*tt-5.0*vv-2.0
     . *s*tt*vv)*s)*pls*s-(3.0*(3.0*tt-vv)-2.0*s*tt*vv))*m2)*aj16+(((
     . tt*v-1.0)*v+2.0*(2.0*tt-vv)*s**2-(3.0*tt*v-2.0)*s)*(pls*s**2+
     . 1.0)-64.0*m2**4*pls*s*tt*vv+16.0*(((5.0*s**2*vv+2.0*v)*tt-2.0*
     . (4.0*tt+vv)*s)*pls+3.0*tt*vv)*m2**3+4.0*((2.0*(16.0*tt-3.0*vv-
     . 4.0*s*tt*vv)*s**2+(5.0*tt*v-3.0)*v-(17.0*tt*v+6.0-12.0*t*vv)*s
     . )*pls+2.0*(10.0*tt-3.0*vv-4.0*s*tt*vv))*m2**2+2.0*(((14.0*tt*v
     . -3.0-4.0*t*vv+2.0*(s*tt*vv-10.0*tt+4.0*vv)*s)*s-((4.0*t*vv-7.0
     . )*t+(6.0*tt*v-1.0)*v))*pls*s+2.0*(s*tt*vv-10.0*tt+4.0*vv)*s+
     . 8.0*tt*v-5.0)*m2)*aj13+(((pls*s**2+1.0)*(s*tt+1.0)-32.0*m2**3*
     . pls*tt)*vv+8.0*(tt+vv+5.0*s*tt*vv)*m2**2*pls-4.0*(3.0*(s*tt+
     . 1.0)*pls*s+tt)*m2*vv)*aj10-(16.0*((2.0*(tt-vv)+5.0*s**2*tt*vv
     . **2-2.0*(3.0*tt+2.0*vv)*s*vv)*pls+3.0*tt*vv**2)*m2**3-(((tt-
     . 3.0*vv-2.0*s*tt*vv)*s-(t*vv-2.0))*(pls*s**2+1.0)+64.0*m2**4*
     . pls*s*tt*vv**2)+2.0*(((2.0*(s*vv-6.0)*s*tt*vv+8.0*tt-5.0*vv)*s
     . +t*vv-3.0*tt*v)*pls*s+2.0*(s*vv-6.0)*s*tt*vv+6.0*tt-5.0*vv)*m2
     . +4.0*((2.0*(11.0*tt+2.0*vv-4.0*s*tt*vv)*s**2*vv+2.0*t*vv-1.0-(
     . 13.0*tt+4.0*vv-4.0*t*vv**2)*s)*pls-2.0*(4.0*s*vv-7.0)*tt*vv)*
     . m2**2)*aj1)*uu
      sr7=-4.0*((2.0*((2.0*uu-3.0*vv+7.0*tt-12.0*s*tt*vv)*s-(2.0*(tt*
     . v-1.0)+t*vv)-((2.0*uu-5.0*vv-5.0*tt+2.0*(uu+4.0*vv)*s*tt)*s+
     . 4.0*tt*v+3.0+(2.0*uu-3.0*vv)*t)*pls*s**2)*m2+((pls*s**2+1.0)*(
     . 2.0*s**2*vv-2.0*s+v)*s-256.0*m2**5*pls*uu)*tt+64.0*(tt+2.0*uu+
     . (5.0*uu-vv)*s*tt)*m2**4*pls-16.0*(((4.0*(2.0*uu-vv)+5.0*tt+(
     . 10.0*uu-vv)*s*tt)*s+2.0*tt*v+1.0+t*uu)*pls+(uu+3.0*vv)*tt)*m2
     . **3+8.0*(((5.0*uu-7.0*vv+tt+(5.0*uu+4.0*vv)*s*tt)*s+3.0*(tt*v+
     . 1.0)+(2.0*uu-vv)*t)*pls*s-(uu-vv+2.0*tt-(uu+9.0*vv)*s*tt))*m2
     . **2)*aj45-(16.0*(6.0*tt+5.0*vv+13.0*s*tt*vv-12.0*m2*tt*vv)*m2
     . **3*pls+((s*vv-2.0)*s*tt-(t*vv+1.0))*(pls*s**2+1.0)-4.0*(2.0*
     . tt*v+3.0+2.0*t*vv+(16.0*s*tt*vv+15.0*tt+18.0*vv)*s)*m2**2*pls+
     . 2.0*((2.0*tt*v+3.0+4.0*t*vv+(8.0*tt+7.0*vv+s*tt*vv)*s)*pls*s+
     . tt-vv-3.0*s*tt*vv)*m2)*aj8-2.0*(2.0*(((2.0*tt*v-17.0-2.0*s*tt)
     . *s**2-(3.0*t**2-3.0*t*v+2.0*v**2)-((2.0*tt*v-13.0)*v+20.0*t)*s
     . +16.0*(tt*v-5.0-3.0*s*tt+4.0*m2*tt)*m2**2)*pls-4.0*(((3.0*tt*v
     . -16.0-4.0*s*tt)*s-((tt*v-4.0)*v+7.0*t))*pls-tt)*m2)*m2+((5.0*t
     . -3.0*v+3.0*s)*s+2.0*(2.0*t**2-2.0*t*v+v**2))*pls*s-2.0*t)*aj44
     . -(2.0*((2.0*t**2*vv+4.0*t-3.0*v-2.0*(uu-vv)*s**3*tt-(2.0*uu-
     . 7.0*vv+2.0*tt)*s**2-((2.0*uu-5.0*vv)*t-2.0)*s)*pls*s+(2.0*uu+
     . vv+tt)*s-(t*vv-4.0))*m2-(16.0*((tt*v-6.0+(uu-5.0*vv)*t-(tt-8.0
     . *uu+14.0*vv-10.0*(uu-vv)*s*tt)*s)*pls+(uu-vv)*tt)*m2**3-(64.0*
     . (2.0*(uu-2.0*vv-tt)+5.0*(uu-vv)*s*tt)*m2**4*pls-((pls*s**2+1.0
     . )*(s**2+t**2)*vv+256.0*(uu-vv)*m2**5*pls*tt)))-4.0*((2.0*t**2*
     . vv+4.0*t-v-10.0*(uu-vv)*s**3*tt-5.0*(2.0*(uu-2.0*vv)+tt)*s**2-
     . (2.0*tt*v-11.0+4.0*(uu-3.0*vv)*t)*s)*pls+2.0*(uu-vv+tt-(uu-vv)
     . *s*tt))*m2**2)*aj43+(((pls*s**2+1.0)*(s*tt+1.0)-32.0*m2**3*pls
     . *tt)*vv+8.0*(tt+vv+3.0*s*tt*vv)*m2**2*pls-4.0*((tt+vv+2.0*s*tt
     . *vv)*pls*s+tt*vv)*m2)*aj41+4.0*((s**2*tt*vv+s*vv+1.0)*pls*s-(s
     . *tt+1.0)*vv-2.0*(3.0*s*vv-1.0-4.0*m2*vv)*m2*pls*s*tt)*aj4*m2-(
     . 2.0*((2.0*uu-3.0*vv+7.0*tt-12.0*s*tt*vv)*s-(2.0*(tt*v-1.0)+t*
     . vv)-((2.0*uu-5.0*vv-5.0*tt+2.0*(uu+4.0*vv)*s*tt)*s+4.0*tt*v+
     . 3.0+(2.0*uu-3.0*vv)*t)*pls*s**2)*m2+((pls*s**2+1.0)*(2.0*s**2*
     . vv-2.0*s+v)*s-256.0*m2**5*pls*uu)*tt+64.0*(tt+2.0*uu+(5.0*uu-
     . vv)*s*tt)*m2**4*pls-16.0*(((4.0*(2.0*uu-vv)+5.0*tt+(10.0*uu-vv
     . )*s*tt)*s+2.0*tt*v+1.0+t*uu)*pls+(uu+3.0*vv)*tt)*m2**3+8.0*(((
     . 5.0*uu-7.0*vv+tt+(5.0*uu+4.0*vv)*s*tt)*s+3.0*(tt*v+1.0)+(2.0*
     . uu-vv)*t)*pls*s-(uu-vv+2.0*tt-(uu+9.0*vv)*s*tt))*m2**2)*aj25+(
     . ((pls*s**2+1.0)*(2.0*s**2-2.0*s*v+v**2)-256.0*m2**5*pls*uu)*tt
     . +64.0*(tt+2.0*uu+5.0*s*tt*uu)*m2**4*pls+8.0*((5.0*s**3*tt*uu+
     . 2.0*s*t*uu+tt*v**2+(6.0*tt+5.0*uu)*s**2)*pls+2.0*tt-uu+s*tt*uu
     . )*m2**2-16.0*((tt*v+1.0+t*uu+(5.0*tt+8.0*uu+10.0*s*tt*uu)*s)*
     . pls+tt*uu)*m2**3-2.0*(((2.0*tt*v+5.0)*v-4.0*t+2.0*(4.0*tt+uu+s
     . *tt*uu)*s**2-(3.0*tt*v+2.0-2.0*t*uu)*s)*pls*s+2.0*(3.0*tt-uu)*
     . s-5.0*tt*v)*m2)*aj24-2.0*((3.0*s*vv-2.0)*pls*s-vv-4.0*m2*pls*s
     . *vv)*aj20*m2*tt+4.0*(2.0*m2-s)*aj19*pls*tt-4.0*((2.0*tt*v+1.0-
     . s*tt)*pls*s-tt-(3.0*(tt*v+1.0)-2.0*s*tt-4.0*m2*tt)*m2*pls)*
     . aj17+((pls*s**2+1.0)*(tt+vv)-16.0*m2**3*pls*tt*vv+2.0*(tt-3.0*
     . vv-2.0*s*tt*vv)*m2*pls*s+8.0*(tt+vv+2.0*s*tt*vv)*m2**2*pls)*
     . aj16+2.0*(((3.0*tt*v-1.0-2.0*s*tt)*s-2.0*((2.0*tt*v-1.0)*v+2.0
     . *t))*pls*s+2.0*(tt*v+1.0-s*tt)-8.0*(3.0*s*tt+4.0-4.0*m2*tt)*m2
     . **2*pls-2.0*(((4.0*tt*v-5.0-6.0*s*tt)*s-3.0*(t+tt*v**2))*pls-
     . 4.0*tt)*m2)*aj14-(2.0*(8.0*(3.0*(tt+vv+s*tt*vv)-4.0*m2*tt*vv)*
     . m2**2*pls+tt-3.0*vv-2.0*s*tt*vv+(s+6.0*t)*pls*s*vv)*m2-(2.0*tt
     . *v+1.0+t*vv-(2.0*tt+vv)*s)*(pls*s**2+1.0)-4.0*((2.0*tt*v+3.0+
     . 2.0*t*vv+(5.0*(tt+2.0*vv)+2.0*s*tt*vv)*s)*pls-2.0*tt*vv)*m2**2
     . )*aj13+(((pls*s**2+1.0)*(s*tt+1.0)-32.0*m2**3*pls*tt)*vv+8.0*(
     . tt+vv+3.0*s*tt*vv)*m2**2*pls-2.0*((5.0*tt+4.0*vv+2.0*s*tt*vv)*
     . pls*s+tt*vv)*m2)*aj1)
      sr8=-4.0*((2.0*(2.0*(4.0*((2.0*(7.0*t*vv-4.0+10.0*s*vv)*s+2.0*t
     . **2*vv-7.0*t+5.0*v)*pls+6.0*vv-16.0*m2*pls*s*vv)*m2-(((4.0*(
     . 9.0*s*vv+12.0*t*vv-8.0)*s+22.0*t**2*vv-37.0*t+16.0*v)*s+2.0*t
     . **3*vv-9.0*t**2+12.0*t*v-4.0*v**2)*pls+4.0*(5.0*t*vv-2.0+5.0*s
     . *vv)))*m2+(((14.0*s*vv+25.0*t*vv-20.0)*s+16.0*t**2*vv-29.0*t+
     . 12.0*v)*s+5.0*t**3*vv-13.0*t**2+12.0*t*v-3.0*v**2)*pls*s+2.0*(
     . 5.0*s*vv+9.0*t*vv-6.0)*s+12.0*t**2*vv-14.0*t+5.0*v)*m2-((2.0*(
     . 2.0*(t*vv-1.0)+s*vv)*s+4.0*t**2*vv-6.0*t+3.0*v)*s+2.0*t**3*vv-
     . 4.0*t**2+3.0*t*v-v**2-((2.0*t-v)*(t-v)-2.0*s**3*vv-4.0*(t*vv-
     . 1.0)*s**2-(2.0*t**2*vv-6.0*t+3.0*v)*s)*pls*s**2))*aj6-2.0*(2.0
     . *(((2.0*t-v)*(t-v)**2+6.0*s**3+2.0*(3.0*t-4.0*v)*s**2+(6.0*t**
     . 2-11.0*t*v+6.0*v**2)*s)*pls+4.0*(t-v+s))*m2-(4.0*((4.0*(2.0*(t
     . -v)+3.0*s)*s+4.0*t**2-6.0*t*v+3.0*v**2)*pls+6.0)*m2**2+((s**2-
     . s*v+t**2-2.0*t*v+v**2)*pls*s+s+t-v)*(s+t-v)-64.0*m2**3*pls*s))
     . *aj7-4.0*(2.0*(t-v+3.0*s-4.0*m2)*m2*pls-((t-v+s)*pls*s+1.0))*
     . aj38+2.0*(2.0*((4.0*t-3.0*v)*(t-v)+4.0*s**2+(12.0*t-11.0*v)*s)
     . *m2*pls-(8.0*(4.0*t-3.0*v+2.0*s)*m2**2*pls+(pls*s**2+3.0*pls*s
     . *t-3.0*pls*s*v+2.0)*(s+t-v)))*aj37+2.0*(4.0*m2**2*pls-2.0*m2*
     . pls*s+1.0)*(4.0*m2*vv-s*vv-t*vv+1.0)*aj33+2.0*(((2.0*(t*vv-1.0
     . )+s*vv)*s+t**2*vv-2.0*t+v)*(pls*s**2+1.0)-8.0*(4.0*t*vv-7.0+
     . 2.0*s*vv)*m2**3*pls+2.0*((20.0*t*vv-21.0+12.0*s*vv)*s+4.0*t**2
     . *vv-11.0*t+6.0*v)*m2**2*pls-(((16.0*t*vv-15.0+9.0*s*vv)*s+7.0*
     . t**2*vv-13.0*t+6.0*v)*pls*s+2.0*(3.0*(t*vv-1.0)+s*vv))*m2)*
     . aj32-((t*vv-1.0+s*vv)*(pls*s**2+1.0)-32.0*m2**3*pls*vv+4.0*(
     . 2.0*t*vv-3.0+8.0*s*vv)*m2**2*pls-2.0*((3.0*t*vv-4.0+5.0*s*vv)*
     . pls*s+3.0*vv)*m2)*aj2-((t*vv-1.0+s*vv)*(pls*s**2+1.0)-32.0*m2
     . **3*pls*vv+8.0*(t*vv-1.0+3.0*s*vv)*m2**2*pls-2.0*((2.0*t*vv-
     . 1.0+4.0*s*vv)*pls*s+vv)*m2)*aj10+(((3.0*t*vv-4.0+2.0*s*vv)*s+t
     . **2*vv-3.0*t+2.0*v)*(pls*s**2+1.0)-128.0*m2**4*pls*s*vv**2+8.0
     . *((4.0*(t*vv-6.0+3.0*s*vv)*s*vv-(4.0*t*vv-11.0))*pls+6.0*vv**2
     . )*m2**3+2.0*(((s*vv-14.0)*s**2*vv-(5.0*t**2*vv-12.0*t+8.0*v)-(
     . t**2*vv**2+15.0*t*vv-21.0)*s)*pls*s+(s*vv-10.0)*s*vv+t**2*vv**
     . 2-9.0*t*vv+11.0)*m2-4.0*((2.0*((t*vv-16.0+3.0*s*vv)*s*vv-(11.0
     . *t*vv-16.0))*s-(2.0*t**2*vv-7.0*t+6.0*v))*pls+4.0*(t*vv-4.0+s*
     . vv)*vv)*m2**2)*aj1)*uu**2
      sr9=4.0*(2.0*(((s*vv-2.0)*s-(t**2*vv-v)-(t-2.0*v+3.0*s)*pls*s**
     . 2+8.0*(5.0*t*vv+4.0+11.0*s*vv-12.0*m2*vv)*m2**3*pls+((9.0*t-
     . 13.0*v+2.0*s**2*vv+(2.0*t*vv+21.0)*s)*pls*s+4.0*(t*vv+2.0))*m2
     . -4.0*(((5.0*t*vv+12.0+6.0*s*vv)*s+t**2*vv+3.0*t-5.0*v)*pls+2.0
     . *vv)*m2**2+(4.0*(((17.0*t-12.0*v+20.0*s)*s+3.0*t**2-5.0*t*v+v
     . **2)*pls+6.0)*m2**2+(t-v+s)*t+(2.0*s-v)*(s+t)*pls*s**2-8.0*(
     . 6.0*t-5.0*v+14.0*s)*m2**3*pls-(((20.0*t-13.0*v+21.0*s)*s+7.0*t
     . **2-13.0*t*v+2.0*v**2)*pls*s+4.0*(3.0*t-v+2.0*s))*m2)*dz1)*aj8
     . -(s**2+s*t-t*v+(s+t-v)*(s+t)*pls*s**2-8.0*(6.0*t-5.0*v+6.0*s)*
     . m2**3*pls+4.0*((13.0*t-9.0*v+10.0*s)*s+3.0*t**2-5.0*t*v+v**2)*
     . m2**2*pls-(((11.0*s+18.0*t-10.0*v)*s+7.0*t**2-10.0*t*v+2.0*v**
     . 2)*pls*s+2.0*(t-2.0*v+3.0*s))*m2)*aj6*dz1)-(((3.0*s*vv-1.0)*s+
     . t**2*vv-t+v+((s*vv-3.0)*s-(t**2*vv-t-v))*pls*s**2)*s-16.0*((
     . 6.0*s**2*vv-2.0*s*t*vv+v)*pls+4.0*vv)*m2**3+8.0*(((t*vv-8.0+
     . 8.0*s*vv)*s-(t**2*vv-2.0*t-v))*pls*s+2.0*t*vv-1.0+8.0*s*vv)*m2
     . **2-2.0*((4.0*t*vv-3.0+11.0*s*vv)*s+t**2*vv-t+v+((2.0*(t*vv-
     . 7.0)+7.0*s*vv)*s-(t**2*vv-5.0*v))*pls*s**2)*m2+(8.0*((5.0*s-
     . 2.0*v)*pls*s*v-2.0*(3.0*s-v)-2.0*((2.0*s-v)*pls*v-2.0)*m2)*m2
     . **2-(4.0*s**2-2.0*s*v+v**2-(2.0*s-v)*pls*s**2*v)*s+2.0*(12.0*s
     . **2-6.0*s*v+v**2-4.0*(2.0*s-v)*pls*s**2*v)*m2)*dz1)*aj45+(16.0
     . *(((11.0*s*vv+16.0*t*vv-8.0)*s+5.0*t**2*vv-2.0*t-2.0*v)*pls+
     . 4.0*vv-8.0*(3.0*s+2.0*t-2.0*m2)*m2*pls*vv)*m2**3-(((t*vv-2.0)*
     . s+t)*s+3.0*t**3*vv-5.0*t**2+3.0*t*v-v**2+(((t*vv-2.0)*s+t)*s-(
     . t**3*vv-3.0*t**2+3.0*t*v-v**2))*pls*s**2)+2.0*(((3.0*t-v)*(t-v
     . )+s**3*vv+(6.0*t*vv-11.0)*s**2+(t**2*vv+2.0*t+v)*s)*pls*s+(4.0
     . *t*vv-7.0+s*vv)*s+11.0*t**2*vv-9.0*t+5.0*v)*m2-8.0*(((t*vv-1.0
     . )*t**2+4.0*s**3*vv+(5.0*t*vv-2.0)*s*t+(10.0*t*vv-11.0)*s**2)*
     . pls+8.0*t*vv-5.0+2.0*s*vv)*m2**2-(8.0*(3.0*(2.0*t-v+6.0*s)*pls
     . *s**2+2.0*(2.0*t-v+5.0*s)+32.0*m2**2*pls*s)*m2**2+(2.0*t**2-
     . 2.0*t*v+v**2+2.0*s**2-(2.0*t**2-2.0*t*v+v**2-2.0*s**2)*pls*s**
     . 2)*s-16.0*((20.0*s**2+v**2+2.0*(2.0*t-v)*s)*pls+6.0)*m2**3-2.0
     . *(2.0*(2.0*(2.0*t-v)+5.0*s)*s+2.0*t**2-2.0*t*v+v**2+2.0*((2.0*
     . t-v+7.0*s)*s-(t**2-t*v+v**2))*pls*s**2)*m2)*dz1)*aj43-2.0*(2.0
     . *(2.0*(t*vv-2.0+3.0*s*vv-4.0*m2*vv)*m2*pls-((t*vv-2.0+s*vv)*
     . pls*s+2.0*vv))*m2+t*vv-1.0+s*vv-pls*s**2+(((s+t)*s**2-48.0*m2
     . **3+4.0*(3.0*t-2.0*v+9.0*s)*m2**2)*pls-2.0*((3.0*t-2.0*v+5.0*s
     . )*pls*s+2.0)*m2)*dz1)*aj41+((pls*s**2-1.0)*(s**2*vv+s-t)-64.0*
     . m2**3*pls*s*vv+8.0*((5.0*s+2.0*t)*pls*s+1.0)*m2**2*vv-2.0*((
     . 5.0*s**2*vv+3.0*s-v)*pls*s+s*vv+1.0)*m2+2.0*(16.0*m2**2*pls*s-
     . 6.0*m2*pls*s**2+3.0*m2*pls*s*v-6.0*m2+pls*s**3-pls*s**2*t-s+t)
     . *(4.0*m2-s-t)*dz1)*aj4+(8.0*(3.0*(2.0*t-v+6.0*s)*pls*s**2+2.0*
     . (2.0*t-v+5.0*s)+32.0*m2**2*pls*s)*m2**2+(2.0*t**2-2.0*t*v+v**2
     . +2.0*s**2-(2.0*t**2-2.0*t*v+v**2-2.0*s**2)*pls*s**2)*s-16.0*((
     . 20.0*s**2+v**2+2.0*(2.0*t-v)*s)*pls+6.0)*m2**3-2.0*(2.0*(2.0*(
     . 2.0*t-v)+5.0*s)*s+2.0*t**2-2.0*t*v+v**2+2.0*((2.0*t-v+7.0*s)*s
     . -(t**2-t*v+v**2))*pls*s**2)*m2)*aj36*dz1-2.0*(2.0*((3.0*t-2.0*
     . v+3.0*s)*pls*s+4.0)*m2+48.0*m2**3*pls-36.0*m2**2*pls*s-12.0*m2
     . **2*pls*t+8.0*m2**2*pls*v-s-t)*aj32*dz1+(((3.0*s*vv-1.0)*s+t**
     . 2*vv-t+v+((s*vv-3.0)*s-(t**2*vv-t-v))*pls*s**2)*s-16.0*((6.0*s
     . **2*vv-2.0*s*t*vv+v)*pls+4.0*vv)*m2**3+8.0*(((t*vv-8.0+8.0*s*
     . vv)*s-(t**2*vv-2.0*t-v))*pls*s+2.0*t*vv-1.0+8.0*s*vv)*m2**2-
     . 2.0*((4.0*t*vv-3.0+11.0*s*vv)*s+t**2*vv-t+v+((2.0*(t*vv-7.0)+
     . 7.0*s*vv)*s-(t**2*vv-5.0*v))*pls*s**2)*m2+(8.0*((5.0*s-2.0*v)*
     . pls*s*v-2.0*(3.0*s-v)-2.0*((2.0*s-v)*pls*v-2.0)*m2)*m2**2-(4.0
     . *s**2-2.0*s*v+v**2-(2.0*s-v)*pls*s**2*v)*s+2.0*(12.0*s**2-6.0*
     . s*v+v**2-4.0*(2.0*s-v)*pls*s**2*v)*m2)*dz1)*aj25+((2.0*(s-v)*s
     . +2.0*t**2+v**2+(2.0*(s-v)*s-(2.0*t**2-v**2))*pls*s**2)*t+256.0
     . *m2**4*pls*s-16.0*((12.0*s**2+v**2+4.0*(3.0*t-v)*s)*pls+6.0)*
     . m2**3+8.0*((6.0*s**3+t*v**2+2.0*(2.0*t-v)*(t-v)*s+2.0*(7.0*t-
     . 2.0*v)*s**2)*pls+2.0*(5.0*t-v+2.0*s))*m2**2+2.0*((3.0*(t-v)*t*
     . v-2.0*s**3-2.0*(6.0*t-v)*s**2-(2.0*t**2-9.0*t*v+v**2)*s)*pls*s
     . -(2.0*(4.0*t-v+s)*s+10.0*t**2-4.0*t*v+v**2))*m2)*aj23*dz1-((s*
     . vv+1.0-2.0*m2*vv)*(pls*s**2-1.0)+8.0*m2**2*pls*s*vv-2.0*(s+t-
     . 4.0*m2)*(pls*s**2-1.0)*dz1)*aj20-((pls*s**2+1.0)*(t*vv-1.0)-
     . 16.0*m2**3*pls*vv+8.0*(t*vv-2.0+2.0*s*vv)*m2**2*pls-2.0*((3.0*
     . t*vv-5.0+s*vv)*pls*s+vv)*m2)*aj16+(16.0*(3.0*t*vv+1.0+s*vv-4.0
     . *m2*vv)*m2**3*pls-(t**2*vv-v-(t*vv-2.0)*s)*(pls*s**2+1.0)+8.0*
     . (((s*vv-6.0)*s-(t**2*vv-3.0*v))*pls+vv)*m2**2+2.0*((4.0*(t-2.0
     . *v)-s**2*vv-3.0*(t*vv-3.0)*s)*pls*s-(t*vv-6.0+s*vv))*m2)*aj13-
     . ((t*vv-2.0+s*vv)*(pls*s**2+1.0)-32.0*m2**3*pls*vv+8.0*(t*vv-
     . 2.0+2.0*s*vv)*m2**2*pls-2.0*((2.0*t*vv-3.0+5.0*s*vv)*pls*s+2.0
     . *vv)*m2)*aj1)*uu
      sr10=4.0*(2.0*(2.0*(2.0*(t-v+3.0*s-4.0*m2)*m2*pls-((t-v+s)*pls*
     . s+1.0))*aj9-(4.0*m2**2*pls-2.0*m2*pls*s+1.0)*aj8-2.0*(2.0*m2-s
     . )**2*aj51+(4.0*(2.0*t*uu-1.0+2.0*s*uu-4.0*m2*uu)*m2-((4.0*t*uu
     . -1.0)*s-(t-v)))*aj50+(pls*s**2-1.0)*aj5+2.0*(16.0*(t-v+2.0*s)*
     . m2**3*pls-(32.0*m2**4*pls+t**2)+((s+t-v)**2*pls*s+4.0*t)*m2-
     . 2.0*((5.0*s+t-v)*(s+t-v)*pls+2.0)*m2**2)*aj49-(2.0*(4.0*((4.0*
     . t-3.0*v+6.0*s)*pls-2.0*uu-8.0*m2*pls)*m2+2.0*(2.0*t*uu-1.0+2.0
     . *s*uu)-(6.0*s+2.0*t-v)*(s+t-v)*pls)*m2-((4.0*t*uu-1.0)*s-(t-v)
     . -(s+t-v)**2*pls*s))*aj48+4.0*(2.0*m2-s)**2*aj47-(4.0*(2.0*t*uu
     . -1.0+2.0*s*uu+3.0*pls*s**2-4.0*(pls*s+uu)*m2)*m2-((4.0*t*uu-
     . 3.0)*s-(t-v)+2.0*pls*s**3))*aj46)-(2.0*(2.0*(2.0*t*uu-1.0)*s*
     . uu-(2.0*t*uu-3.0)-(6.0*s*t*uu-3.0*s+2.0*t**2*uu+2.0*t+3.0*v)*
     . pls*s+64.0*m2**3*pls*uu-16.0*((2.0*t*uu+1.0+4.0*s*uu)*pls-uu**
     . 2)*m2**2)*m2+(2.0*s*uu-1.0)*s+2.0*t**2*uu+t+v+(t+v-s)*pls*s**2
     . +8.0*(((8.0*t*uu+1.0+3.0*s*uu)*s+t**2*uu+t+v)*pls-2.0*(t*uu+
     . 1.0+s*uu)*uu)*m2**2)*aj45+2.0*((s**2+2.0*t**2+(t+v)*s)*pls*s-
     . 2.0*t+2.0*(4.0*(t+2.0*v+s)*m2-(3.0*s+t)*(s+t+v))*m2*pls)*aj44+
     . (2.0*(2.0*(2.0*t*uu-1.0)*s*uu-(2.0*t*uu-3.0)-(6.0*s*t*uu-3.0*s
     . +2.0*t**2*uu+2.0*t-v)*pls*s+64.0*m2**3*pls*uu-8.0*((4.0*t*uu+
     . 1.0+8.0*s*uu)*pls-2.0*uu**2)*m2**2)*m2+(2.0*s*uu+1.0)*s+2.0*t
     . **2*uu-t+v+(t-v-s)*pls*s**2+4.0*(((16.0*t*uu+1.0+6.0*s*uu)*s+
     . 2.0*t**2*uu+t-2.0*v)*pls-4.0*(t*uu+1.0+s*uu)*uu)*m2**2)*aj43-(
     . pls*s**2-1.0)*aj4-4.0*(2.0*m2-s)**2*aj31-4.0*(2.0*(4.0*m2-3.0*
     . s)*m2*pls+pls*s**2-1.0)*aj30*s-8.0*(2.0*m2-s)**2*aj27+2.0*(4.0
     . *(2.0*t*uu-1.0+2.0*s*uu+3.0*pls*s**2-4.0*(pls*s+uu)*m2)*m2-((
     . 4.0*t*uu-3.0)*s-(t-v)+2.0*pls*s**3))*aj26+(2.0*(2.0*(2.0*t*uu-
     . 1.0)*s*uu-(2.0*t*uu-3.0)-(6.0*s*t*uu-3.0*s+2.0*t**2*uu+2.0*t+
     . 3.0*v)*pls*s+64.0*m2**3*pls*uu-16.0*((2.0*t*uu+1.0+4.0*s*uu)*
     . pls-uu**2)*m2**2)*m2+(2.0*s*uu-1.0)*s+2.0*t**2*uu+t+v+(t+v-s)*
     . pls*s**2+8.0*(((8.0*t*uu+1.0+3.0*s*uu)*s+t**2*uu+t+v)*pls-2.0*
     . (t*uu+1.0+s*uu)*uu)*m2**2)*aj25+2.0*(16.0*((2.0*t*uu+1.0+4.0*s
     . *uu)*pls-uu**2-4.0*m2*pls*uu)*m2**3-(t**2*uu+v+s**2*uu-(s-t)*
     . pls*s**2)+((2.0*t**2*uu+2.0*t+5.0*v+6.0*(t*uu-1.0)*s)*pls*s-
     . 2.0*((2.0*t*uu-1.0)*s*uu-(t*uu-2.0)))*m2-4.0*(((8.0*t*uu+1.0+
     . 3.0*s*uu)*s+t**2*uu+t+v)*pls-2.0*(t*uu+1.0+s*uu)*uu)*m2**2)*
     . aj24-2.0*(2.0*m2-s)*aj17*pls+2.0*((2.0*t+v)*pls*s-1.0-2.0*(t+
     . 2.0*v-2.0*m2)*m2*pls)*aj14-(4.0*(2.0*m2-s)*m2*pls+pls*s**2+1.0
     . )*aj13)
       fsir=coer*(sr1+sr2+sr3+sr4+sr5+sr6+sr7+sr8+sr9+sr10)
!     fsir=aj6
!     print *,'fsir',fsir
      if(fsir.lt.0)nn=nn+1
      if(ikey.eq.-1)then
       fir=-4d0*((aj5+aj7+m2*aj1*vv**2+aj14)+(t-2d0*m2)*
     .      (aj23+aj13*vv)+(s-2d0*m2)*
     .      (aj36+vv*aj4)-(s+t-2d0*m2)*(vv*aj6+aj24))
       fsir=fsir-alfa/4d0/pi**2*fir*sig0
      endif
      end
