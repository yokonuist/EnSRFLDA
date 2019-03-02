subroutine vtm_cal(ref_x,vtm,hr)

! From Zeigler, 1978

real,parameter :: zfrez=3000
real,parameter :: zice=8000
real,parameter :: h0=7000
real,parameter :: denom=1./(zice-zfrez)

real  :: ref_x
real  :: vtm
real  :: hr
real  :: rhofact
real  :: refz
real  :: s1,s2

refz=10**(0.1*ref_x)
rhofact=exp(0.4*hr/h0)

if( hr < zfrez ) then
  vtm=2.6*refz**0.107*rhofact
else if( hr < zice ) then
  s1=(zice-hr)*denom
  s2=2*(hr-zfrez)*denom
  vtm=s1*2.6*refz**0.107*rhofact+s2
else
  vtm=2.0
endif

end subroutine vtm_cal


SUBROUTINE terv1D(irsg,rhobar,q,vtr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the terminal velocities (vtr) of rain water qr, snow qs and
!  hail qg.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Goddard Cumulus Ensemble Modeling Group, NASA
!
!  MODIFICATION HISTORY:
!  12/14/2004
!  Mingjing Tong
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  irsg    Flag identifying hydrometeor field as rain, snow or hail
!          =  0 for rain;  qpt is rain
!          =  1 for snow;  qpt is snow
!          =  2 for hail;  qpt is hail
!  q       Hydrometeor field defined at the scalar point
!  rhobar  Air density defined at scalar point (g/cm**3)
!
!  OUTPUT:
!
!  vtr     Vertical velocity defined at the scalar point (cm/s)
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: irsg
  REAL :: q
  REAL :: rhobar
  REAL :: vtr
  REAL :: tema, temb, temc, temd, rho0cgs
  REAL :: temp, interp, f1, f2, rstep
  INTEGER :: INDEX
  real :: rho0
  REAL    :: pwr2
  REAL    :: pwr0625
  real    :: vrc,vsc,vgc

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  
  rho0=1.2250
  rho0cgs = rho0*0.001
  rstep = 1.0/50.0E-10

  temc = 0.0
  
  call stcstice1(vrc,vsc,vgc)

  IF (irsg == 0) THEN    ! irsg = 0    for rain water (qr)

    tema=SQRT(rho0cgs/rhobar)
    temb=rhobar*q

    IF (temb > 1.e-16) THEN

     !  original Lin scheme

        temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
        index = int(temp)
        
        call initlktb1(index,pwr2,pwr0625,1)
        f1 = pwr2
        call initlktb1(index+1,pwr2,pwr0625,1)
        f2 = pwr2

        vtr = vrc * tema * (f1 + (f2-f1)*(temp-index))

    ELSE
      vtr = 0.0
    END IF
    temc = max(temc,vtr)

!   print *,'max rain fall speed is=',temc,irsg

  ELSE IF( irsg == 1) THEN  !irsg = 1 for snow (qs)

    tema=SQRT(rho0cgs/rhobar)
    temb=rhobar*q

    IF (temb > 1.e-16) THEN

    ! original Lin scheme.......
        temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
        index = int(temp)
        
        call initlktb1(index,pwr2,pwr0625,2)        
        f1 = pwr0625
        call initlktb1(index+1,pwr2,pwr0625,2)        
        f2 = pwr0625
        vtr = vsc * tema * ( f1 + (f2-f1)*(temp-index) )

    ELSE
      vtr = 0.0
    END IF

    temc = max(temc,vtr)

!   print *,'max snow fall speed is=',temc,irsg

  ELSE IF( irsg == 2) THEN  ! irsg = 2   for graupel (qg)

    tema=sqrt(rho0cgs/rhobar)
    temb=rhobar*q
    IF (temb > 1.e-16) THEN

     ! Original method
        temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
        index = int(temp)
        
        call initlktb1(index,pwr2,pwr0625,2)           
        f1 = pwr0625
        call initlktb1(index+1,pwr2,pwr0625,2) 
        f2 = pwr0625
        interp = f1 + (f2 - f1) * (temp - index)
        vtr = vgc / sqrt(rhobar) * interp * interp
   
        temc = max(temc,vtr)

   ELSE
     vtr = 0.0
   END IF
   temc = max(temc,vtr)

!   print *,'max hail fall speed is=',temc,irsg

  END IF

  RETURN
END SUBROUTINE terv1D

SUBROUTINE initlktb1(index,pwr2,pwr0625,ireg)

  INTEGER :: index
  REAL    :: rhoqx
  REAL    :: pwr2
  REAL    :: pwr0625
  integer :: ireg

  rhoqx = 0.05 * index * 1.0E-7 ! EMK cgs units
  if(ireg == 1) then
    pwr2    =  ( rhoqx ) ** 0.2
  else if(ireg == 2) then
    pwr0625 =  ( rhoqx ) ** 0.0625
  endif
  
  RETURN
END SUBROUTINE initlktb1

SUBROUTINE stcstice1(vrc,vsc,vgc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set constants used by the ice microphyscs parameterization routine
!  ICECVT
!
!  Lin et.al.  J. Clim. Appl. Meteor.  22, 1065-1092
!  Modified and coded by tao and simpson (JAS, 1989; Tao, 1993)
!
!-----------------------------------------------------------------------

 real   :: roqr,tnw,roqs,tns,roqg,tng,cpi,cpi2,grvt
 real   :: ag,bg,as,bs,aww,bww,ga4b,ga4g,ga4d
 real   :: ac1,bc1,cc1,dc1,cd1,cd2,zrc,zsc,zgc
 real   :: vrc,vsc,vgc

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Define the density and size distribution of precipitation
!
!-----------------------------------------------------------------------
!
  roqr = 1.
  tnw = .08
  roqs = .1
  tns = .03
  roqg = .913
  tng = .0004
!
  cpi = 4.*ATAN(1.)
  cpi2 = cpi*cpi
  grvt = 980.

!
!-----------------------------------------------------------------------
!
!  Define the coefficients used in terminal velocity
!
!-----------------------------------------------------------------------
!

  ag = 1400.
  bg = .5
  as = 152.93
  bs = .25
  aww= 2115.
  bww= .8

! original Lin configuration
  ga4b = 17.83779
  ga4g = 11.63177
  ga4d = 8.285063

!
!-----------------------------------------------------------------------
!
!  Lin et al., 1983
!
!-----------------------------------------------------------------------
!
  ac1 = aww
  bc1 = bww
  cc1 = as
  dc1 = bs
  cd1 = 6.e-1
  cd2 = 4.*grvt/(3.*cd1)
  zrc = (cpi*roqr*tnw)**0.25
  zsc = (cpi*roqs*tns)**0.25
  zgc = (cpi*roqg*tng)**0.25
  vrc = ac1*ga4b/(6.*zrc**bww)
  vsc = cc1*ga4d/(6.*zsc**bs)
  vgc = ga4g*SQRT(cd2*roqg/zgc)/6.  
  
  RETURN
  
END SUBROUTINE stcstice1  

