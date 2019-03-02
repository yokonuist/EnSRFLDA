subroutine cal_mapfactor_lambert(lat,mapfact)

real                 ::   lat
real                 ::   mapfact

real,parameter       ::   pi=3.1415926
real,parameter       ::   a=6371
real,parameter       ::   le=11423.37
real,parameter       ::   hd2jd=pi/180.0

real                 ::   l
real                 ::   k

k = (log(sin(60*hd2jd))-log(sin(30*hd2jd)))/(log(tan(60*hd2jd)/2)-log(tan(30*hd2jd/2)))

l =  le*(cos(lat*hd2jd)/(1+sin(lat*hd2jd)))**(k)

mapfact = k*l/(a*(1-sin(lat*hd2jd)**2)**(0.5))

end subroutine cal_mapfactor_lambert

subroutine cal_rho(mub,mu,qv,phl,phu,rdnw,p_top,znu,t,       &
                   nx,ny,nz,rho,p                            &                  
                  )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine cal_rho
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer         ::  nx,ny,nz                  
real            ::  mub,mu,qv,phl,phu,rdnw,znu,p
real            ::  p_top
real            ::  rho  

REAL, PARAMETER    :: pi = 3.1415926535897932346
REAL, PARAMETER    :: gas_constant = 287.     ! Value used in WRF.
REAL, PARAMETER    :: gas_constant_v = 461.6  ! Value used in WRF.
REAL, PARAMETER    :: cp = 7.*gas_constant/2. ! Value used in WRF.
REAL, PARAMETER    :: t_kelvin = 273.15
REAL, PARAMETER    :: kappa = gas_constant / cp
REAL, PARAMETER    :: rd_over_rv = gas_constant / gas_constant_v  
REAL, PARAMETER    :: gravity = 9.81        ! m/s - value used in MM5.


real            :: albn,aln,ppb,ttb,qvf1,t
real            :: cvpm,cpovcv,ps0,ts0,tis0,tlp



ps0  = 100000.0    ! Base sea level pressure
ts0  = 300.0       ! Base potential temperature for all levels.
tis0 = 290.0       ! Base sea level temperature
tlp  = 50.0        ! temperature difference from 1000mb to 300mb


cvpm =  - (1. - gas_constant/cp)
cpovcv = cp / (cp - gas_constant)

   ppb  = znu * mub + p_top
   ttb  = (tis0 + tlp*log(ppb/ps0)) * (ps0/ppb)**kappa
   albn = (gas_constant/ps0) * ttb * (ppb/ps0)**cvpm

   qvf1 = 1. + qv / rd_over_rv
   aln  = -1. / (mub+mu) * ( albn*mu + rdnw *(phu - phl) )
               
   p = ps0 * ( (gas_constant*(ts0+t)*qvf1) / &
                       (ps0*(aln+albn)) )**cpovcv            

   rho= 1.0 / (albn+aln)
                                                    
end subroutine cal_rho        

subroutine cal_tc(t,p,tc,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine cal_tc
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer         ::  nx,ny,nz 
real            ::  t,p,tc

REAL, PARAMETER    :: gas_constant = 287.     ! Value used in WRF.
REAL, PARAMETER    :: cp = 7.*gas_constant/2. ! Value used in WRF. 
REAL, PARAMETER    :: t_kelvin = 273.15
real            :: ps0
REAL, PARAMETER    :: kappa = gas_constant / cp

ps0  = 100000.0    ! Base sea level pressure


  tc=(t+300)*(p/ps0)**kappa-t_kelvin 


end subroutine cal_tc

subroutine rv_operator(u,v,w,rv_x,azimuth,evl,rf_x,ref_r,ref_s,ref_h,   &
                       qr,qs,qgr,rho,uctb,vctb,wctb                     &
                       )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine rv_operator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
include 'namelist.inc'

real                      :: qr,qs,qgr,rho
real                      :: rf_x,ref_r,ref_s,ref_h
real                      :: u,v,w
real                      :: rv_x
real                      :: azimuth,evl
real                      :: vtm,vtr,vts,vth
real                      :: uctb,vctb,wctb


i=0
call terv1D(i,rho,qr,vtr)
i=1
call terv1D(i,rho,qs,vts)
i=2
call terv1D(i,rho,qgr,vth)

if(rf_x >= 10 ) then
vtm=(vtr*ref_r+vts*ref_s+vth*ref_h)/rf_x
vtm=vtm*0.01
else
vtm=0
endif	

uctb=u*cos(evl)*sin(azimuth)
vctb=v*cos(evl)*cos(azimuth)
wctb=(w-vtm)*sin(evl)


rv_x =  uctb + vctb + wctb
     

end subroutine rv_operator


subroutine ref_operator(qr,qs,qgr,ref_x,ref_r,ref_s,ref_h,rho,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine ref_operator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real                      :: qr,qs,qgr
real                      :: t
real                      :: ref_x,ref_r,ref_s,ref_h
real,parameter            :: pi=3.1415926
real,parameter            :: rhor=1000
real,parameter            :: rhos=100
real,parameter            :: rhoi=917
real,parameter            :: rhoh=913
real,parameter            :: nr=8*1e6
real,parameter            :: ns=3*1e6
real,parameter            :: nh=4*1e4
real,parameter            :: Ki=0.176
real,parameter            :: Kr=0.93
real,parameter            :: para1=1e18*720
real,parameter            :: para2=0.9

ref_r=0
ref_s=0
ref_h=0

ref_r=para1*(rho*qr)**1.75/(pi**1.75*nr**0.75*rhor**1.75)

if(t<0) then
ref_s=para1*Ki*rhos**0.25*(rho*qs)**1.75/(pi**1.75*Kr*ns**0.75*rhoi**2)
else
ref_s=para1*(rho*qs)**1.75/(pi**1.75*ns**0.75*rhos**0.175)
endif

ref_h=(para1/(pi**1.75*nh**0.75*rhoh**1.75))**0.95*(rho*qgr)**1.6625


if(ref_r+ref_s+ref_h > 1.0) then
ref_x=10*para2*log10(ref_r+ref_s+ref_h)
else
ref_x=0
endif

if(ref_r >= 1.0) then
ref_r=10*para2*log10(ref_r)*para2
else
ref_r=1e-8
endif

if(ref_s >= 1.0) then
ref_s=10*para2*log10(ref_s)*para2
else
ref_r=1e-8
endif

if(ref_h >= 1.0) then
ref_h=10*para2*log10(ref_h)*para2
else
ref_r=1e-8
endif

end subroutine

subroutine innovation_cal(innovationmean,obs,hxmean)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine innovation_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real        :: obs
real        :: innovationmean,hxmean
    
  innovationmean=obs-hxmean
  
end subroutine innovation_cal

subroutine zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine zzh_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer     :: ens_num,numstat
real        :: xstat(numstat,ens_num),hx(ens_num),zzh(numstat)
real        :: hxmean
integer     :: temp1

zzh=0


do i=1,numstat
 do j=1,ens_num
  zzh(i)=zzh(i)+xstat(i,j)*(hx(j)-hxmean)
 enddo
  zzh(i)=zzh(i)/(ens_num-1)
enddo



end subroutine zzh_cal

subroutine hph_cal(hx,hxmean,hph,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine hph_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer     :: ens_num
real        :: hx(ens_num)
real        :: hxmean
real        :: hph
hph=0
do i=1,ens_num
  hph=hph+(hx(i)-hxmean)*(hx(i)-hxmean)
enddo
  
hph=hph/(ens_num-1)

end subroutine hph_cal

subroutine kalman_cal(zzh,hph,obs_err,numstat,kalman)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine kalman_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer     :: numstat
real        :: zzh(numstat)
real        :: hph,obs_err
real        :: kalman(numstat)


do i=1,numstat
  
  kalman(i)=zzh(i)/(hph+obs_err**2)
 !print*,'kalman(i),zzh(i),hph',kalman(i),zzh(i),hph
enddo


end subroutine kalman_cal

subroutine kalman_alpha_cal(kalman_alpha,hph,obs_err)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine kalman_alpha_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real        :: hph,obs_err
real        :: kalman_alpha

kalman_alpha=1./(1+sqrt(obs_err**2/(hph+obs_err**2)))
!print*,'kalman_alpha',kalman_alpha

end subroutine kalman_alpha_cal

subroutine innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,    &
                          xstat,xstat_mean,kalman,hx_pert                 &
                          )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine innovate_xstat
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                          
integer       :: numstat,ens_num
real          :: xstat(numstat,ens_num),xstat_mean(numstat),xstatf(numstat,ens_num)
real          :: innovationmean,kalman_alpha
real          :: kalman(numstat),hx_pert(ens_num)
integer       :: temp1,temp2,temp3



  do i=1,numstat  
   xstat_mean(i)=xstat_mean(i)+kalman(i)*innovationmean        !  innovate ens mean    
   do j=1,ens_num 
    xstat(i,j)=xstat(i,j)-kalman_alpha*kalman(i)*hx_pert(j)    !  innovate ens perturbation
  !   print*,'xstat(i,j),xstat_mean(i)',xstat(i,j),xstat_mean(i)
   enddo   
  enddo



end subroutine innovate_xstat  


subroutine xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine xstat_mean_seperate_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
integer       ::  numstat,ens_num
real          ::  xstat(numstat,ens_num)
real          ::  xstat_mean(numstat)

xstat_mean=0


do i=1,numstat
 do j=1,ens_num
   xstat_mean(i)=xstat_mean(i)+xstat(i,j)
 enddo
   xstat_mean(i)=xstat_mean(i)/ens_num
enddo



do i=1,numstat
 do j=1,ens_num
   xstat(i,j)=xstat(i,j)-xstat_mean(i)

!*******************test code********************* 
! if(i > 4032000 .and. xstat(i,j) .ne. 0)then
! print*,'ok,find good things? i,xstat(i,j)',i,xstat(i,j)
! endif
!**********************xstat(i,j)*****************
 
enddo
enddo

  

end subroutine xstat_mean_seperate_cal

subroutine xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine xstat_mean_merge_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
integer       ::  numstat,ens_num
real          ::  xstat(numstat,ens_num)
real          ::  xstat_mean(numstat)


!print*,'11xstat(3920001,1)',xstat(3920001,1)
do i=1,numstat
 do j=1,ens_num
   xstat(i,j)=xstat(i,j)+xstat_mean(i)
 enddo
enddo

!print*,'11xstat(3920001,1)',xstat(3920001,1)

end subroutine xstat_mean_merge_cal

subroutine check_xstat_mean__cal(xstat,xstat_mean,numstat,ens_num,xstatf,xstatff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine xstat_mean_seperate_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
integer       ::  numstat,ens_num
real          ::  xstat(numstat,ens_num)
real          ::  xstatf(numstat,ens_num)
real          ::  xstatff(numstat,ens_num)
real          ::  xstat_mean(numstat)
real          ::  temp

xstat_mean=0


do i=1,numstat
 do j=1,ens_num
   xstat_mean(i)=xstat_mean(i)+xstat(i,j)

 enddo
   xstat_mean(i)=xstat_mean(i)/ens_num
enddo


open(8888,file='./check.log')
write(8888,*)'numstat,ens_num,pert_ana,ana,inv_pert'
do i=1,numstat
 do j=1,ens_num
   xstatf(i,j)=xstat(i,j)-xstat_mean(i)
   temp=xstat(i,j)-xstatff(i,j)
!*******************test code********************* 
 if(abs(xstatf(i,j)) > 10 .or. abs(temp) > 8)then
! print*,'ok,find good things? i,xstat(i,j)',i,xstatf(i,j)
write(8888,*)i,j,xstatf(i,j),xstat(i,j),temp
 endif
!**********************xstat(i,j)*****************
 
enddo
enddo

  

end subroutine check_xstat_mean__cal


subroutine dis_correlation(hor_distance,vert_distance,local_dis,verti_dis,dis_cor)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine dis_correlation
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
real     :: hor_distance,vert_distance
real     :: dis_cor,dis_cor_h,dis_cor_v
real     :: a,b,z
real     :: local_dis,verti_dis

a=local_dis
b=hor_distance
z=b/a

if(b >= 0 .and. b < 0.5*a ) then
 dis_cor_h=-28./33.*z**5+8./11.*z**4+20./11.*z**3-80./33.*z**2+1
endif
if(b >= 0.5*a .and. b < a ) then
 dis_cor_h=20./33.*z**5-16./11.*z**4+100./33.*z**2-45./11.*z+51./22.-7./44./z
endif
if(b >= a .and. b < 1.5*a ) then
 dis_cor_h=-4./11.*z**5+16./11.*z**4-10./11.*z**3-100./33.*z**2+5*z-61./22.+115./132./z
endif
if(b >= 1.5*a .and. b< 2*a) then	
 dis_cor_h=4./33.*z**5-8./11.*z**4+10./11.*z**3+80./33.*z**2-80./11.*z+64./11.-32./33./z
endif
if(b>2*a) then
 dis_cor_h=0
endif

a=verti_dis
b=vert_distance
z=b/a

if(b >= 0 .and. b < 0.5*a ) then
 dis_cor_v=-28./33.*z**5+8./11.*z**4+20./11.*z**3-80./33.*z**2+1
endif
if(b >= 0.5*a .and. b < a ) then
 dis_cor_v=20./33.*z**5-16./11.*z**4+100./33.*z**2-45./11.*z+51./22.-7./44./z
endif
if(b >= a .and. b < 1.5*a ) then
 dis_cor_v=-4./11.*z**5+16./11.*z**4-10./11.*z**3-100./33.*z**2+5*z-61./22.+115./132./z
endif
if(b >= 1.5*a .and. b< 2*a) then	
 dis_cor_v=4./33.*z**5-8./11.*z**4+10./11.*z**3+80./33.*z**2-80./11.*z+64./11.-32./33./z
endif
if(b>2*a) then
 dis_cor_v=0
endif

dis_cor=dis_cor_h*dis_cor_v

end subroutine dis_correlation


subroutine ijk_to_istat(ivar,ii,jj,kk,istat,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine ijk_to_istat
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
implicit none

integer :: nx,ny,nz
integer :: ivar,ii,jj,kk,istat

istat=(ivar-1)*nx*ny*nz+(kk-1)*nx*ny+(jj-1)*nx+ii

end subroutine ijk_to_istat

subroutine ladijk_to_istat(ivar,ii,jj,kk,istat,nx,ny,nz,nlev)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine ladijk_to_istat
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
implicit none

integer :: nx,ny,nz,nlev
integer :: ivar,ii,jj,kk,istat

if(ivar < 11)then
istat=(ivar-1)*nx*ny*nz+(kk-1)*nx*ny+(jj-1)*nx+ii
!print*,'ii,jj,kk,istat',ii,jj,kk,istat
endif

if(ivar == 11)then
istat=(ivar-1)*nx*ny*nz+(kk-1)*nx*ny+(jj-1)*nx+ii
endif

if(ivar == 12)then
istat=(ivar-2)*nx*ny*nz+nx*ny*4+(kk-1)*nx*ny+(jj-1)*nx+ii
endif

if(ivar == 13)then
istat=(ivar-3)*nx*ny*nz+nx*ny*4*2+(jj-1)*nx+ii  !TSK
endif

!print*,'ivar,istat(istatmp from obs_real_sub)from subroutine ladijk_to_istat',ivar,istat 

end subroutine ladijk_to_istat

subroutine istat_to_ijk(ivar,ii,jj,kk,istat,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine istat_to_ijk
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
implicit none

integer :: nx,ny,nz
integer :: ivar,ii,jj,kk,istat
integer :: tem1,tem2,tem3,tem4,tem5

ivar=int(istat/(nx*ny*nz))+1

tem1=istat-(ivar-1)*nx*ny*nz
kk  =int(tem1/(nx*ny))+1

tem2=istat-(ivar-1)*nx*ny*nz-(kk-1)*nx*ny
jj  =int(tem2/nx)+1

ii  =istat-(ivar-1)*nx*ny*nz-(kk-1)*nx*ny-(jj-1)*nx


end subroutine istat_to_ijk

subroutine ladistat_to_ijk(ivar,ii,jj,kk,istat,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine istat_to_ijk
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!

!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
implicit none

integer :: nx,ny,nz
integer :: ivar,ii,jj,kk,istat
integer :: tem1,tem2,tem3,tem4,tem5

ivar=int(istat/(nx*ny*4))+1

tem1=istat-(ivar-1)*nx*ny*4
kk  =int(tem1/(nx*ny))+1

tem2=istat-(ivar-1)*nx*ny*4-(kk-1)*nx*ny
jj  =int(tem2/nx)+1

ii  =istat-(ivar-1)*nx*ny*4-(kk-1)*nx*ny-(jj-1)*nx


end subroutine ladistat_to_ijk

subroutine check_hydr(nx,ny,nz,ens_num,ensqv,ensqr,ensqi,ensqs,ensqgr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine check_hydr
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer          :: nx,ny,nz
integer          :: ens_num
integer          :: nn
real             :: ensqv(nx,ny,nz,ens_num)
real             :: ensqr(nx,ny,nz,ens_num)
real             :: ensqi(nx,ny,nz,ens_num)
real             :: ensqs(nx,ny,nz,ens_num)
real             :: ensqgr(nx,ny,nz,ens_num)


do k=1,nz
do j=1,ny
do i=1,nx
do nn=1,ens_num 
 if(ensqv (i,j,k,nn)<0) ensqv (i,j,k,nn)=0
 if(ensqr (i,j,k,nn)<0) ensqr (i,j,k,nn)=0
 if(ensqi (i,j,k,nn)<0) ensqi (i,j,k,nn)=0
 if(ensqs (i,j,k,nn)<0) ensqs( i,j,k,nn)=0
 if(ensqgr(i,j,k,nn)<0) ensqgr(i,j,k,nn)=0
enddo
enddo
enddo
enddo

end subroutine check_hydr

subroutine check_soillayer(nx,ny,nz,ens_num,enssmois,enstslb,enstsk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine check_hydr
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer          :: nx,ny,nz
integer          :: ens_num
integer          :: nn
real             :: enssmois(nx,ny,4,ens_num)
real             :: enstslb(nx,ny,4,ens_num)
real             :: enstsk(nx,ny,ens_num)

do k=1,4
do j=1,ny
do i=1,nx
do nn=1,ens_num 
 if(enssmois (i,j,k,nn)<0) then
   enssmois (i,j,k,nn)=0
 !print*,'wrong enssmois (i,j,k,nn)',nn,enssmois (i,j,k,nn)
 endif 
 if(enssmois (i,j,k,nn)> 1) then
   enssmois(i,j,k,nn)=0.99
 endif
if(enstslb (i,j,k,nn)< 233.15)then
 enstslb (i,j,k,nn)=233.15
 endif
 if(enstslb (i,j,k,nn)> 343.15)then
 enstslb (i,j,k,nn)=343.15
 endif
enddo
enddo
enddo
enddo

!do j=1,ny
!do i=1,nx
!do nn=1,ens_num 
!if(enstsk (i,j,nn)< 243.15)then
! print*,'ckeck enstsk',enstsk(i,j,nn)
! enstsk (i,j,nn)=243.15
! endif
! if(enstsk (i,j,nn)> 313.15)then
! print*,'ckeck enstsk',enstsk(i,j,nn)
! enstsk (i,j,nn)=313.15
! endif
!enddo
!enddo
!enddo
do j=1,ny
do i=1,nx
do nn=1,ens_num 

!************************
!*** fill data like this:
!  ---------
!  1 2 3 4 |
!  2 3 4 5 |
!****************************

if( i == nx  .and. enstsk(i,j,nn) == 0)then
 enstsk(i,j,nn)=enstsk(i-1,j,nn)
elseif(i == 1)then
 enstsk(i,j,nn)=enstsk(i+1,j,nn)
elseif( j == ny .and. enstsk(i,j,nn) == 0)then
 enstsk(i,j,nn)=enstsk(i,j-1,nn)
elseif( j == 1 )then
 enstsk(i,j,nn) = enstsk(i,j+1,nn)  
elseif(enstsk (i,j,nn) < 233.15)then
 print*,'ckeck enstsk',enstsk(i,j,nn),i,j,nn
 enstsk (i,j,nn) = 233.15
elseif(enstsk (i,j,nn)> 343.15)then
 print*,'ckeck enstsk',enstsk(i,j,nn)
 enstsk (i,j,nn)=343.15
 endif
enddo
enddo
enddo

return

end subroutine check_soillayer

subroutine rm_correlate(zzh,numstat,ivar,nx,ny,nz,ireg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine rm_correlate
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer         :: nx,ny,nz
integer         :: numstat
real            :: zzh(numstat)
integer         :: ivar
integer         :: ipoint
integer         :: ireg
integer         :: istat

if (ivar < 11)then
ipoint=(ivar-1)*nx*ny*nz+1

do k=1,1
do j=1,ny
do i=1,nx
  zzh(ipoint)=0
  ipoint=ipoint+1
enddo
enddo
enddo

endif

!!for soil layer!!
if (ivar >= 11)then
ipoint=(ivar-1)*nx*ny*4+1

do k=1,1
do j=1,ny
do i=1,nx
  zzh(ipoint)=0
  ipoint=ipoint+1
enddo
enddo
enddo

endif


if(ireg == 2 .and. (ivar < 7 ) ) then
ipoint=(ivar-1)*nx*ny*nz+1


do istat=ipoint,ipoint+nx*ny*nz-1
  zzh(istat)=0
enddo
  

endif  


if(ireg == 1 .and. ivar > 6 ) then
ipoint=(ivar-1)*nx*ny*nz+1



do istat=ipoint,ipoint+nx*ny*nz-1
  zzh(istat)=0
enddo

  

endif  

end subroutine rm_correlate


subroutine relax_inflation_scheme(xstat,xstatf,numstat)

include 'namelist.inc'

real    :: xstat(numstat,ens_num)
real    :: xstatf(numstat,ens_num)

do i=1,numstat
 do j=1,ens_num
   xstat(i,j)=(1-ratio_af)*xstatf(i,j)+ratio_af*xstat(i,j)
 enddo
enddo

end subroutine relax_inflation_scheme

!***************************************************************
!***************************************************************
!***************************************************************
!       The following codes are used for new local tactic
!***************************************************************  
!***************************************************************
!***************************************************************

subroutine rm_correlate2(zzh,numstat,gridsgn,nx,ny,nz,ireg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine rm_correlate
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer         :: nx,ny,nz
integer         :: numstat
real            :: zzh(numstat)
integer         :: ivar
integer         :: ipoint
integer         :: ireg
integer         :: istat
integer         :: gridsgn(numstat)



if( ireg == 2  ) then
ivar=6
ipoint=(ivar-1)*nx*ny*nz+1 

do istat=1,numstat
 if(gridsgn(istat) < ipoint) then
  zzh(istat)=0
 endif
enddo
 
endif  

if( ireg == 1  ) then
ivar=7
ipoint=(ivar-1)*nx*ny*nz+1 

do istat=1,numstat
 if(gridsgn(istat) >= ipoint) then
  zzh(istat)=0
 endif
enddo
 
endif  


end subroutine rm_correlate2

subroutine merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)

integer       ::  gridnum,ens_num,numstat
integer       ::  gridsgn(gridnum)
real          ::  xstat_sub(gridnum,ens_num)
real          ::  xstat(numstat,ens_num)


do i=1,gridnum
 do j=1,ens_num
   xstat(gridsgn(i),j)=xstat_sub(i,j)
 enddo
enddo

end subroutine merge_local
!***************************************************************
!***************************************************************

!***************************************************************
!***************************************************************

subroutine cal_minvalue_1d(x,min,dim1)

integer     ::  dim1
real        ::  x(dim1)
real        ::  min

real        ::  temr01,temr02,temr03


temr01=9999999

do i=1,dim1
  if( x(i) < temr01 ) then
      temr01=x(i)
  endif
enddo

min=temr01

end subroutine cal_minvalue_1d

subroutine cal_maxvalue_1d(x,max,dim1)

integer     ::  dim1
real        ::  x(dim1)
real        ::  max

real        ::  temr01,temr02,temr03


temr01=-9999999

do i=1,dim1
  if( x(i) > temr01 ) then
      temr01=x(i)
  endif
enddo

max=temr01

end subroutine cal_maxvalue_1d
