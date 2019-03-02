
subroutine obs_sim_grid_u(uwest,ueast,hx,ens_num)
!##########################################################
!####                      			    #######
!####         hx for sim obs at grid point          #######
!####      all analysis vars will be assimilated    #######
!####         grid is define at mass grid           #######
!####        u,v,w,ph are averaged at mass grid     #######
!####                                               #######   
!##########################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine obs_sim_grid_u
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
integer    :: ens_num
real       :: uwest(ens_num),ueast(ens_num)
real       :: hx(ens_num)

do nn=1,ens_num

 hx(nn)=0.5*(uwest(nn)+ueast(nn))

enddo


end subroutine obs_sim_grid_u

subroutine obs_sim_grid_v(vsouth,vnorth,hx,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine obs_sim_grid_v
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
integer    :: ens_num
real       :: vsouth(ens_num),vnorth(ens_num)
real       :: hx(ens_num)

do nn=1,ens_num

 hx(nn)=0.5*(vsouth(nn)+vnorth(nn))

enddo


end subroutine obs_sim_grid_v

subroutine obs_sim_grid_w(whigh,wlow,hx,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine obs_sim_grid_w
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
integer    :: ens_num
real       :: whigh(ens_num),wlow(ens_num)
real       :: hx(ens_num)

do nn=1,ens_num

 hx(nn)=0.5*(whigh(nn)+wlow(nn))

enddo


end subroutine obs_sim_grid_w

subroutine obs_sim_grid_ph(phhigh,phlow,hx,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine obs_sim_grid_ph
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
integer    :: ens_num
real       :: phhigh(ens_num),phlow(ens_num)
real       :: hx(ens_num)

do nn=1,ens_num

 hx(nn)=0.5*(phhigh(nn)+phlow(nn))

enddo


end subroutine obs_sim_grid_ph

subroutine obs_sim_grid_mass(mass_grid_var,hx,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine obs_sim_grid_mass
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
integer    :: ens_num
real       :: mass_grid_var(ens_num)
real       :: hx(ens_num)

do nn=1,ens_num

 hx(nn)=mass_grid_var(nn)

enddo

end subroutine obs_sim_grid_mass

subroutine point_distance_sim(iobs,jobs,kobs,ivar,igrid,jgrid,kgrid,analysis_var_num,  &
                              hor_distance,vert_distance,phb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine point_distance_sim
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
real,parameter   :: gravity = 9.81        ! m/s - value used in MM5.
integer          :: analysis_var_num
integer          :: iobs,jobs,kobs
integer          :: igrid,jgrid,kgrid
integer          :: ivar                          
real             :: hor_distance,vert_distance
real             :: xdis,ydis,zdis
integer          :: temp1
real             :: phb(nx,ny,nz)


! step 1 cal the horizontal distance
  
  if( ivar == 1 ) then  !for u point
   if(iobs >= igrid ) then
    xdis=abs(iobs-igrid)*dx+0.5*dx
   else
    xdis=abs(iobs-igrid)*dx-0.5*dx
   endif
    ydis=abs(jobs-jgrid)*dy   
  endif
  
  if( ivar == 2 ) then  !for v point
    xdis=abs(iobs-igrid)*dx
   if(jobs >= jgrid ) then   
    ydis=abs(jobs-jgrid)*dy+0.5*dy
   else
    ydis=abs(jobs-jgrid)*dy-0.5*dy
   endif 
  endif    

  if( ivar > 2) then ! for non u v point
    xdis=abs(iobs-igrid)*dx
    ydis=abs(jobs-jgrid)*dy
  endif  
     
  hor_distance=sqrt(xdis**2+ydis**2)  
  
  
! step 2 cal vertial distance
  
  if (ivar == 3 .or. ivar == 4 )  then ! for w and ph point

   zdis= abs(  (phb(iobs,jobs,kobs)+phb(iobs,jobs,kobs+1))/2     &
             -  phb(igrid,jgrid,kgrid)                           &
             )/gravity                                                        
  endif
  
  if (ivar .ne. 3 .and. ivar .ne. 4 ) then  ! for non w ph point

   zdis= abs(  (phb(iobs,jobs,kobs)+phb(iobs,jobs,kobs+1))/2          &
             - (phb(igrid,jgrid,kgrid)+phb(igrid,jgrid,kgrid+1))/2    &
             )/gravity 
             
  endif  
  
  vert_distance=zdis 
  
   
! finished cal distance
end subroutine point_distance_sim  

subroutine localization_cal_sim(zzh,numstat,iobs,jobs,kobs,phb, &
                                analysis_var_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine localization_cal_sim
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
integer     :: numstat
integer     :: analysis_var_num
integer     :: iobs,jobs,kobs
real        :: zzh(numstat)
real        :: pbh(nx,ny,nz)
real        :: hor_distance,vert_distance
real        :: dis_cor

integer     :: ivar
integer     :: ipoint
integer     :: itemp

ipoint=1

do k=1,nz
do j=1,ny
do i=1,nx

   do ivar=1,analysis_var_num
   
 ! step 1 cal the distance
   call point_distance_sim(iobs,jobs,kobs,ivar,i,j,k,analysis_var_num,  &
                           hor_distance,vert_distance,phb)
                      
 ! step 2 cal the correlation
   if(optlocal  == 1 ) then
     if( hor_distance <= local_dis .and. vert_distance <= verti_dis)  then
       dis_cor=1
     else 
       dis_cor=0
     endif
   else if(optlocal  == 2 ) then
     call dis_correlation(hor_distance,vert_distance,local_dis,verti_dis,dis_cor) 
   endif
                             
 ! step 3 add the correlation to zzh
  
   itemp=(ivar-1)*nx*ny*nz+ipoint
   zzh(itemp)=zzh(itemp)*dis_cor
     
   enddo
   
   ipoint=ipoint+1                           
enddo
enddo
enddo
                              

end subroutine localization_cal_sim  


subroutine radardata_kr_cal(nx,ny,nz,ir,jr,hr,phb,krw,krm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine radardata_kr_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real,parameter :: gravity=9.81
integer        :: nx,ny,nz
integer        :: ir,jr
real           :: hr
real           :: phb(nx,ny,nz) 
integer        :: krw,krm
real           :: temp1,temp2,temp3

do k=1,nz
 
 temp1=phb(ir,jr,k)/gravity
 if( temp1 > hr ) then
    krw = k -1                            !  the w level just lower than hr
  exit
 endif                
  
enddo

do k=1,nz-1
  
 temp1=0.5*(phb(ir,jr,k)+phb(ir,jr,k+1))/gravity
 if( temp1 > hr ) then
    krm = k-1                             !  the mass level just lower than hr
  exit
 endif   
   
enddo

end subroutine radardata_kr_cal


subroutine interpolation_sim_radar_wind(nx,ny,nz,ir,jr,krw,krm,hr,phb,     &
                                        uuw,uue,ulw,ule,                   &
                                        vus,vun,vls,vln,                   &
                                        wu,wl,                             &
                                        ur,vr,wr                           &
                                        )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine interpolation_sim_radar_wind
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real,parameter :: gravity=9.81
integer        :: nx,ny,nz
integer        :: ir,jr
real           :: hr
real           :: phb(nx,ny,nz) 

real           :: uuw,uue,ulw,ule,                   &
                  vus,vun,vls,vln,                   &
                  wu,wl
                  
integer        :: krw,krm
real           :: ur,vr,wr
real           :: uu,ul,vu,vl
real           :: huuv,hluv,huw,hlw
real           :: ww,wwu,wwl                 ! weight for up and low point w
real           :: wm,wmu,wml                 ! weight for up and low point uv
real           :: temp1,temp2,temp3


huuv=0.5*(phb(ir,jr,krm+1)+phb(ir,jr,krm+2))/gravity
hluv=0.5*(phb(ir,jr,krm)+phb(ir,jr,krm+1))/gravity
huw=phb(ir,jr,krw+1)/gravity
hlw=phb(ir,jr,krw  )/gravity

ww=(hr-hlw)/(huw-hlw)
wwl=1-ww
wwu=ww

wm=(hr-hluv)/(huuv-hluv)
wml=1-wm
wmu=wm

uu=0.5*(uuw+uue)
ul=0.5*(ulw+ule)

vu=0.5*(vus+vun)
vl=0.5*(vls+vln)

ur=wml*ul+wmu*uu
vr=wml*vl+wmu*vu
wr=wwl*wl+wwu*wu


!if(ir==nx/2.and.jr==ny/2) then
!print*,'hr',hr
!print*,'huuv',huuv,hluv,huw,hlw
!print*,'(huuvw-hluvw)',(huuv-hluv)
!print*,'(huw-hlw)',(huw-hlw)
!print*,'ww',ww,wwl,wwu
!print*,'wm',wm,wml,wmu
!print*,'uu',uu,ul,vu,vl,wu,wl
!endif

end subroutine interpolation_sim_radar_wind

subroutine interpolation_sim_radar_mass(nx,ny,nz,ir,jr,krm,hr,phb,         &
                                        qru,qrl,qsu,qsl,qgru,qgrl,         &
                                        rhou,rhol,tcu,tcl,                 &
                                        qrr,qsr,qgrr,rhor,tcr              &
                                        )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine interpolation_sim_radar_mass
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real,parameter :: gravity=9.81
integer        :: nx,ny,nz
integer        :: ir,jr
real           :: hr
real           :: qru,qrl,qsu,qsl,qgru,qgrl,tcu,tcl
real           :: phb(nx,ny,nz) 
integer        :: krm
real           :: qrr,qsr,qgrr,tcr
real           :: rhou,rhol,rhor
real           :: hum,hlm
real           :: wm,wmu,wml                 ! weight for up and low point m
real           :: temp1,temp2,temp3


hum=0.5*(phb(ir,jr,krm+1)+phb(ir,jr,krm+2))/gravity
hlm=0.5*(phb(ir,jr,krm)+phb(ir,jr,krm+1))/gravity

wm=(hr-hlm)/(hum-hlm)
wml=1-wm
wmu=wm

qrr=wml*qrl+wmu*qru
qsr=wml*qsl+wmu*qsu
qgrr=wml*qgrl+wmu*qgru
rhor=wml*rhol+wmu*rhou
tcr=wml*tcl+wmu*tcu

end subroutine interpolation_sim_radar_mass

subroutine radar_distance_sim(iobs,jobs,hr,ivar,igrid,jgrid,kgrid,analysis_var_num,  &
                              hor_distance,vert_distance,phb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine radar_distance_sim
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
real,parameter   :: gravity = 9.81        ! m/s - value used in MM5.
integer          :: analysis_var_num
integer          :: iobs,jobs,kobs
integer          :: igrid,jgrid
real             :: hr
integer          :: ivar                          
real             :: hor_distance,vert_distance
real             :: xdis,ydis,zdis
integer          :: temp1
real             :: phb(nx,ny,nz)


! step 1 cal the horizontal distance
  
  if( ivar == 1 ) then  !for u point
   if(iobs >= igrid ) then
    xdis=abs(iobs-igrid)*dx+0.5*dx
   else
    xdis=abs(iobs-igrid)*dx-0.5*dx
   endif
    ydis=abs(jobs-jgrid)*dy   
  endif
  
  if( ivar == 2 ) then  !for v point
    xdis=abs(iobs-igrid)*dx
   if(jobs >= jgrid ) then   
    ydis=abs(jobs-jgrid)*dy+0.5*dy
   else
    ydis=abs(jobs-jgrid)*dy-0.5*dy
   endif 
  endif    

  if( ivar > 2) then ! for non u v point
    xdis=abs(iobs-igrid)*dx
    ydis=abs(jobs-jgrid)*dy
  endif  
     
  hor_distance=sqrt(xdis**2+ydis**2)  
  
! step 2 cal vertial distance
  
  if (ivar == 3 .or. ivar == 4 )  then ! for w and ph point

   zdis= abs(  hr -  phb(igrid,jgrid,kgrid)/gravity )                                                       
  endif
  
  if (ivar .ne. 3 .and. ivar .ne. 4 ) then  ! for non w ph point

   zdis= abs(  hr - (phb(igrid,jgrid,kgrid)+phb(igrid,jgrid,kgrid+1))/2/gravity)
             
  endif  
  
  vert_distance=zdis 
   
! finished cal distance
end subroutine radar_distance_sim  

subroutine localization_cal_radar_sim(zzh,numstat,iobs,jobs,hr,phb, &
                                      analysis_var_num,ipoint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine localization_cal_radar_sim
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
integer     :: numstat
integer     :: analysis_var_num
integer     :: iobs,jobs
real        :: hr
real        :: zzh(numstat)
real        :: pbh(nx,ny,nz)
real        :: hor_distance,vert_distance
real        :: dis_cor
integer     :: ii,jj,kk
integer     :: ivar
integer     :: ipoint
integer     :: itemp

real        :: local_h_dis,local_v_dis

real        :: kkr,nzr

   call istat_to_ijk(ivar,ii,jj,kk,ipoint,nx,ny,nz)


 ! step 1 cal the distance
   call radar_distance_sim(iobs,jobs,hr,ivar,ii,jj,kk,analysis_var_num,  &
                           hor_distance,vert_distance,phb)   

   if(ivar <= 2) then
     local_h_dis=local_dis
     local_v_dis=verti_dis  
   endif
   if(ivar == 3) then
     local_h_dis=local_dis_w
     local_v_dis=verti_dis_w   
   endif   
   if(ivar == 4) then
     local_h_dis=local_dis_ph
     local_v_dis=verti_dis_ph   
   endif 
   if(ivar == 5) then
     local_h_dis=local_dis_t
     local_v_dis=verti_dis_t   
   endif    
   if(ivar >5 ) then
     local_h_dis=local_dis_hydro
     local_v_dis=verti_dis_hydro   
   endif      
                  
 ! step 2 cal the correlation
   
   if(optlocal  == 1 ) then
     if( hor_distance <= local_h_dis .and. vert_distance <= local_v_dis)  then
       dis_cor=1
     else 
       dis_cor=0
     endif
   else if(optlocal  == 2 ) then
     call dis_correlation(hor_distance,vert_distance,local_h_dis,local_v_dis,dis_cor) 
   endif
                        
 ! step 3 add the correlation to zzh
    
   zzh(ipoint)=zzh(ipoint)*dis_cor                              



end subroutine localization_cal_radar_sim 

subroutine rv_maker(u,v,w,rv_x,azimuth,evl,rf_x,hr,ref_r,ref_s,ref_h,   &
                    qr,qs,qgr,rho                                       &
                    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine rv_maker
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
real                      :: u,v,w
real                      :: qr,qs,qgr,rho
real                      :: rf_x,ref_r,ref_s,ref_h
real                      :: rv_x
real                      :: azimuth,evl
real                      :: hr
real                      :: vtm,vtr,vts,vth

if(rf_x >= 5) then

i=0
call terv1D(i,rho,qr,vtr)
i=1
call terv1D(i,rho,qs,vts)
i=2
call terv1D(i,rho,qgr,vth)

vtm=(vtr*ref_r+vts*ref_s+vth*ref_h)/rf_x
vtm=vtm*0.01

rv_x =  u*cos(evl)*sin(azimuth)  &
      + v*cos(evl)*cos(azimuth)  &
      +(w-vtm)*sin(evl)   
else
rv_x = -99999
endif	  


end subroutine rv_maker


subroutine ref_maker(qr,qs,qgr,ref_x,ref_r,ref_s,ref_h,rho,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine ref_maker
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
real,parameter            :: para2=1

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

if(ref_r+ref_s+ref_h >= 1.0) then
ref_x=10*para2*log10(ref_r+ref_s+ref_h)
else
ref_x=-99999
endif

if(ref_r >= 1.0) then
ref_r=10*para2*log10(ref_r)
else
ref_r=1e-8
endif

if(ref_s >= 1.0) then
ref_s=10*para2*log10(ref_s)
else
ref_r=1e-8
endif

if(ref_h >= 1.0) then
ref_h=10*para2*log10(ref_h)
else
ref_r=1e-8
endif

end subroutine ref_maker

subroutine inflation_sim_cal(xstat,xstat_mean,numstat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine inflation_sim_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
integer       :: numstat
real          :: xstat(numstat,ens_num),xstat_mean(numstat)

call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)

do i=1,numstat
do j=1,ens_num
  xstat(i,j)=xstat(i,j)*(1+coef_inflat)
enddo
enddo

call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)

end subroutine inflation_sim_cal


!***************************************************************
!***************************************************************
!***************************************************************
!       The following codes are used for new local tactic
!***************************************************************  
!***************************************************************
!***************************************************************

subroutine localization_cal_radar_sim2(dis_rho,numstat,iobs,jobs,hr,phb,       &
                                       analysis_var_num,ii,jj,kk,ivar,         &
                                       local_h_dis,local_v_dis,ipoint)
include 'namelist.inc'
integer     :: numstat
integer     :: analysis_var_num
integer     :: iobs,jobs
real        :: hr
real        :: dis_rho(numstat)
real        :: pbh(nx,ny,nz)
real        :: hor_distance,vert_distance
real        :: dis_cor
integer     :: ii,jj,kk
integer     :: ivar
integer     :: ipoint
integer     :: itemp

real        :: local_h_dis,local_v_dis


 ! step 1 cal the distance
   call radar_distance_sim(iobs,jobs,hr,ivar,ii,jj,kk,analysis_var_num,  &
                           hor_distance,vert_distance,phb)   

         
 ! step 2 cal the correlation
   
   if(optlocal  == 1 ) then
     if( hor_distance <= local_h_dis .and. vert_distance <= local_v_dis)  then
       dis_rho(ipoint)=1
     else 
       dis_rho(ipoint)=0
     endif
   else if(optlocal  == 2 ) then
     call dis_correlation(hor_distance,vert_distance,local_h_dis,local_v_dis,dis_rho(ipoint)) 
   endif                              



end subroutine localization_cal_radar_sim2 

subroutine prepare_info4local(analysis_var_num,ra_data_num,numstat,                          &                               
                              ir,jr,hr,xradar,yradar,hradar,phb,       			     &
                              dis_rho,gridsgn,gridnum                                        &
                             )
include 'namelist.inc'
                             
  integer               :: analysis_var_num,ra_data_num,numstat
  integer               :: iobs,ivar
  integer               :: nn  

  real                  :: phb(nx,ny,nz)                 
  real                  :: hr
  integer               :: ir,jr

  real                  :: dis_rho(numstat)
  integer               :: istatmp
  integer               :: gridnum,gridsgn(numstat)
  integer               :: ii,jj,kk  
  integer               :: istart,iend,jstart,jend
  character*8           :: timestart,timeend
  real                  :: local_h_dis,local_v_dis
  
  real                  :: local_v_dis_temp,local_h_dis_temp

do ivar=1,analysis_var_num

!***************************************************************  
!    mark the area for calculate distance correlation coeffient
!***************************************************************  

   if(ivar <= 2) then
     local_h_dis=local_dis
     local_v_dis=verti_dis
   endif
   if(ivar == 3) then
     local_h_dis=local_dis_w
     local_v_dis=verti_dis_w
   endif   
   if(ivar == 4) then
     local_h_dis=local_dis_ph
     local_v_dis=verti_dis_ph
   endif 
   if(ivar == 5) then
     local_h_dis=local_dis_t
     local_v_dis=verti_dis_t
   endif    
   if(ivar >5 ) then
     local_h_dis=local_dis_hydro
     local_v_dis=verti_dis_hydro
   endif       
   
istart=ir-(int(2*local_h_dis/dx)+1)
if(istart < 1) istart=1
iend=ir+(int(2*local_h_dis/dx)+1)
if(iend > nx) iend=nx
jstart=jr-(int(2*local_h_dis/dy)+1)
if(jstart < 1) jstart=1
jend=jr+(int(2*local_h_dis/dy)+1)
if(jend > ny) jend=ny
   
!***************************************************************  
!       calculate distance correlation coeffient
!*************************************************************** 

 
do k=1,nz-8

do j=jstart,jend
do i=istart,iend
  
  call ijk_to_istat(ivar,i,j,k,istatmp,nx,ny,nz)

  call localization_cal_radar_sim2(dis_rho,numstat,ir,jr,hr,phb,       &
                                   analysis_var_num,i,j,k,ivar,        &
                                   local_h_dis,local_v_dis,istatmp)   
                                                                                                        
enddo    
enddo
enddo                      

enddo  

!***************************************************************  
!      count an record the grid number and grid subscript
!*************************************************************** 
gridnum=0
do istatmp=1,numstat 
  if( dis_rho(istatmp) > 0.0001 ) then
    gridsgn(istatmp)=istatmp
    gridnum=gridnum+1
  endif
enddo

end subroutine prepare_info4local


subroutine simple_local(zzh,gridnum,dis_rho)

integer    :: gridnum
real       :: zzh(gridnum)
real       :: dis_rho(gridnum)

do i=1,gridnum
 zzh(i)=zzh(i)*dis_rho(i)
enddo

end subroutine simple_local