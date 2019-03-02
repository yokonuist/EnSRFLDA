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


subroutine interpolation_radar_wind(nx,ny,nz,ir,jr,krw,krm,hr,phb,     &
                                    uuw,uue,ulw,ule,                   &
                                    vus,vun,vls,vln,                   &
                                    wu,wl,                             &
                                    ur,vr,wr                           &
                                    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine interpolation_radar_wind
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.03.27  
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

end subroutine interpolation_radar_wind

subroutine interpolation_radar_mass(nx,ny,nz,ir,jr,krm,hr,phb,         &
                                    qru,qrl,qsu,qsl,qgru,qgrl,         &
                                    rhou,rhol,tcu,tcl,                 &
                                    qrr,qsr,qgrr,rhor,tcr              &
                                   )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine interpolation_radar_mass
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


huuv=0.5*(phb(ir,jr,krm+1)+phb(ir,jr,krm+2))/gravity
hluv=0.5*(phb(ir,jr,krm)+phb(ir,jr,krm+1))/gravity

wm=(hr-hluv)/(huuv-hluv)
wml=1-wm
wmu=wm

qrr=wml*qrl+wmu*qru
qsr=wml*qsl+wmu*qsu
qgrr=wml*qgrl+wmu*qgru
rhor=wml*rhol+wmu*rhou
tcr=wml*tcl+wmu*tcu

end subroutine interpolation_radar_mass

subroutine radar_distance(iobs,jobs,hr,ivar,igrid,jgrid,kgrid,analysis_var_num,  &
                          hor_distance,vert_distance,phb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine radar_distance
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
end subroutine radar_distance  

subroutine localization_cal_radar(zzh,numstat,iobs,jobs,hr,phb, &
                                  analysis_var_num,ipoint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine localization_cal_radar
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.03.27   
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
   call radar_distance(iobs,jobs,hr,ivar,ii,jj,kk,analysis_var_num,  &
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



end subroutine localization_cal_radar


subroutine inflation_cal(xstat,xstat_mean,numstat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine inflation_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.03.27  
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

!if(i > 4032001 .and. xstat(i,j) > 0) then
!if(i == 4045586 .and. xstat(i,j) > 0) then
!print*,i,xstat(i,j),j
!endif
 
  xstat(i,j)=xstat(i,j)*(1+coef_inflat)

!if(i > 4032001 .and. xstat(i,j) > 0) then
!if(i == 4045586 .and. xstat(i,j) > 0) then
!print*,i,xstat(i,j),j
!endif
!***** ok , tested on 2013/5/12

enddo
enddo

call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)

end subroutine inflation_cal


!***************************************************************
!***************************************************************
!***************************************************************
!       The following codes are used for new local tactic
!***************************************************************  
!***************************************************************
!***************************************************************

subroutine localization_cal_radar2(dis_rho,numstat,iobs,jobs,hr,phb,       &
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
   call radar_distance(iobs,jobs,hr,ivar,ii,jj,kk,analysis_var_num,  &
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



end subroutine localization_cal_radar2

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

  call localization_cal_radar2(dis_rho,numstat,ir,jr,hr,phb,       &
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

!**********************************************************************************************************
!**********************************************************************************************************
!**********************************************************************************************************
!**********************************************************************************************************
!**********************************************************************************************************
!**********************************************************************************************************
!**********************************************************************************************************
!**********************************************************************************************************


subroutine snddata_ks_cal(nx,ny,nz,ism,jsm,hgtsnd,phb,ksm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine snddata_ks_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real,parameter :: gravity=9.81
integer        :: nx,ny,nz
integer        :: ism,jsm
real           :: hgtsnd
real           :: phb(nx,ny,nz) 
integer        :: ksm
real           :: temp1,temp2,temp3

ksm=-1

do k=1,nz-2
  
 temp1=0.5*(phb(ism,jsm,k)+phb(ism,jsm,k+1))/gravity
 if( temp1 > hgtsnd ) then
    ksm = k-1                             !  the mass level just lower than hgtsnd
	if(ksm<1) ksm=-1
  exit
 endif   
   
enddo

end subroutine snddata_ks_cal

subroutine snddata_ks_lad_cal(nx,ny,nz,ism,jsm,hgtsnd,ksm,dzs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine snddata_ks_lad_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real,parameter :: gravity=9.81
integer        :: nx,ny,nz
integer        :: ism,jsm
real           :: hgtsnd
real           :: phb(nx,ny,nz)
real           :: dzs(4) 
integer        :: ksm
real           :: temp1,temp2,temp3

ksm=-1
print*,dzs(1)

if(dzs(1) .ne. 0.1)then
print*,'dzs(1)=',dzs(1),'there is something wrong in your code '
endif

do k=1,4
 print*,'k',k,'dzs=',dzs(k) 
 !temp1=0.5*(phb(ism,jsm,k)+phb(ism,jsm,k+1))/gravity
 temp1=dzs(k)
print*,'k',k,'temp1=',temp1
!stop
 if( temp1 >= hgtsnd ) then
    ksm = k-1                             !  the mass level just lower than hgtsnd
	if(ksm<1) ksm=-1
  exit
 endif   
enddo

end subroutine snddata_ks_lad_cal

subroutine snddata_isjs_cal(nx,ny,nz,ism,jsm,xlon,xlat,latsnd,lonsnd,minlon,minlat,maxlon,maxlat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine radardata_kr_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real,parameter :: gravity=9.81
integer        :: nx,ny,nz
integer        :: ism,jsm
real           :: latsnd,lonsnd
real           :: xlon(nx,ny),xlat(nx,ny)
real           :: temp1,temp2,temp3
real           :: minlon,minlat,maxlon,maxlat


ism=-1
jsm=-1

if(latsnd < minlat .or. latsnd > maxlat  .and.&
   lonsnd < minlon .or. lonsnd > maxlon) then
  
  ism=-1
  jms=-1

else

  do j=1,ny-1
  do i=1,nx-1
        
    if(lonsnd >= xlon(i,j) .and. xlon(i+1,j) >= lonsnd  .and.&
       latsnd >= xlat(i,j) .and. xlat(i,j+1) >= latsnd ) then
        ism=i
        jsm=j
print*,'ism,jsm,i,j,lonsnd,latsnd',ism,jsm,i,j,lonsnd,latsnd
         exit
    endif

  enddo
  enddo
  

endif


end subroutine snddata_isjs_cal

subroutine obs_latlonhgt_to_ijk(latsnd,lonsnd,hgtsnd,isnd,jsnd,ksnd,    &
                                phb,nx,ny,nz,xlat,xlon,                 &
				                         minlon,minlat,maxlon,maxlat,            &
                                ism,jsm,ksm,wwsn,wsew,wul               &
				)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine obs_latlonhgt_to_ijk
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

integer          ::  nx,ny,nz
real             ::  phb(nx,ny,nz)
real             ::  latsnd,lonsnd,hgtsnd
real             ::  xlat(nx,ny),xlon(nx,ny)
real             ::  isnd,jsnd,ksnd
integer          ::  ism,jsm,ksm
real             ::  minlon,minlat,maxlon,maxlat

real             ::  wwsn,wsew,wul
real             ::  temr01,temr02,temr03
real             ::  gravity=9.81

  
  call snddata_isjs_cal(nx,ny,nz,ism,jsm,xlon,xlat,latsnd,lonsnd,minlon,minlat,maxlon,maxlat)
  call snddata_ks_cal(nx,ny,nz,ism,jsm,hgtsnd,phb,ksm)

  if( ism == -1 .or. jsm== -1 .or. ksm == -1 ) return
       
  wwsn=(xlat(ism,jsm+1)-latsnd)/(xlat(ism,jsm+1)-xlat(ism,jsm))
  
  jsnd=jsm+1-wwsn
   
  wsew=(xlon(ism+1,jsm)-lonsnd)/(xlon(ism+1,jsm)-xlon(ism,jsm)) 

  isnd=ism+1-wsew

  temr01=0.5*(phb(ism,jsm,ksm+2)+phb(ism,jsm,ksm+1))/gravity
  temr02=0.5*(phb(ism,jsm,ksm  )+phb(ism,jsm,ksm+1))/gravity

  wul=(temr01-hgtsnd)/(temr01-temr02)

  ksnd=ksm+1-wul
  
  print*,ksm,ksnd,wul,temr01,temr02,phb(ism,jsm,ksm+2)
  pause
  
end subroutine obs_latlonhgt_to_ijk

subroutine obs_latlonhgtlad_to_ijk(latsnd,lonsnd,hgtsnd,isnd,jsnd,ksnd,    &
                                phb,dzs,nx,ny,nz,xlat,xlon,                 &
				                         minlon,minlat,maxlon,maxlat,            &
                                ism,jsm,ksm,wwsn,wsew,wul               &
				)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine obs_latlonhgtlad_to_ijk
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

integer          ::  nx,ny,nz
real             ::  phb(nx,ny,nz)
real             ::  dzs(4)
!!!
real             ::  latsnd,lonsnd,hgtsnd
real*4           ::  xlat(nx,ny),xlon(nx,ny)
real             ::  isnd,jsnd,ksnd
integer          ::  ism,jsm,ksm
real             ::  minlon,minlat,maxlon,maxlat

real             ::  wwsn,wsew,wul
real             ::  temr01,temr02,temr03,temp1,temp2
real             ::  gravity=9.81

!print*,dzs(1) ok! 
!  call snddata_isjs_cal(nx,ny,nz,ism,jsm,xlon,xlat,latsnd,lonsnd,minlon,minlat,maxlon,maxlat)
!  call snddata_ks_lad_cal(nx,ny,nz,ism,jsm,hgtsnd,phb,ksm,dzs)

ism=-1
jsm=-1

if(latsnd < minlat .or. latsnd > maxlat  .and.&
   lonsnd < minlon .or. lonsnd > maxlon) then
  
  ism=-1
  jms=-1

else

  do j=1,ny-1
  do i=1,nx-1
        
    if(lonsnd >= real(xlon(i,j)) .and. real(xlon(i+1,j)) >= lonsnd  .and.&
       latsnd >= real(xlat(i,j)) .and. real(xlat(i,j+1)) >= latsnd ) then
        ism=i
        jsm=j
         exit
    endif

  enddo
  enddo
  

endif

!-----------------------------------------------------------------------------!
!-------for test  as for snddata_ks_lad_cal or snddata_ks_cal ------!
!-----------------------------------------------------------------------------!
ksm=-1
!print*,'dzs(1),ism,jsm,',dzs(1),ism,jsm

if(dzs(1) .ne. 0.1)then
print*,'dzs(1)=',dzs(1),'there is something wrong in your code '
endif

dzs(0)=0.0 !!! first level
dzs(5)=1.5 !!! lost level
dzs(6)=2.0 !!! proper?

!do k=1,4
! if( hgtsnd >= dzs(k) .and. dzs(k+1) > hgtsnd) then
!    ksm = k 
! !   temp2 = (hgtsnd-dzs(k))/(dzs(k+1)-dzs(k)) 
! !   temp1=ksm+(hgtsnd-dzs(k))/(dzs(k+1)-dzs(k))                           !  the mass level just lower than hgtsnd
!	if(ksm<1) ksm=-1
!  exit
! endif   
!enddo
!temp1=0.5*dzs(1)
do k=1,4
 temp2=dzs(k-1)+0.5*dzs(k-1)+0.5*dzs(k)
 !temp1=0.5*(dzs(k-1)+dzs(k))
 !temp2=0.5*(dzs(k)+dzs(k+1))
 if( temp2 > hgtsnd) then
 !if(dzs(k+1) > hgtsnd .and. hgtsnd >= dzs(k) ) then
    ksm = k-1
 !   temp2 = (hgtsnd-dzs(k))/(dzs(k+1)-dzs(k)) 
                              !  the mass level just lower than hgtsnd
	if(ksm<1) ksm=-1
  exit
 endif   
enddo

  if( ism == -1 .or. jsm== -1 .or. ksm == -1 ) return 

!print*,'ism,jsm,ksm,hgtsnd,latsnd,lonsnd,',ism,jsm,ksm,hgtsnd,latsnd,lonsnd
    
  wwsn=(xlat(ism,jsm+1)-latsnd)/(xlat(ism,jsm+1)-xlat(ism,jsm))
  
  jsnd=jsm+1-wwsn
   
  wsew=(xlon(ism+1,jsm)-lonsnd)/(xlon(ism+1,jsm)-xlon(ism,jsm)) 

  isnd=ism+1-wsew

  temr01=dzs(ksm)+0.5*dzs(ksm)+0.5*dzs(ksm+1)
  temr02=dzs(ksm-1)+0.5*dzs(ksm-1)+0.5*dzs(ksm)

 ! temr01=0.5*(dzs(ksm)+dzs(ksm+1))
 ! temr02=0.5*(dzs(ksm-1)+dzs(ksm))

 ! wul=(hgtsnd-dzs(ksm))/(dzs(ksm+1)-dzs(ksm))

!   wul=temp2

!  ksnd=ksm+wul
  
!   ksnd=temp1

wul=(temr01-hgtsnd)/(temr01-temr02)

  ksnd=ksm+1-wul
  
end subroutine obs_latlonhgtlad_to_ijk



subroutine snd_distance(isnd,jsnd,ksnd,hgtsnd,ivar,igrid,jgrid,kgrid,analysis_var_num,  &
                        hor_distance,vert_distance,phb,ism,jsm,ksm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine snd_distance
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
real,parameter   :: gravity = 9.81        ! m/s - value used in MM5.
integer          :: analysis_var_num
real             :: isnd,jsnd,ksnd
integer          :: ism,jsm,ksm
integer          :: igrid,jgrid
real             :: hr
integer          :: ivar                          
real             :: hor_distance,vert_distance
real             :: xdis,ydis,zdis
integer          :: temp1
real             :: phb(nx,ny,nz)


! step 1 cal the horizontal distance
  
  if( ivar == 1 ) then  !for u point

   if(ism >= igrid ) then
    xdis=abs(ism-igrid)*dx+0.5*dx+(isnd-real(ism))*dx
   else
    xdis=abs(ism-igrid)*dx-0.5*dx-(isnd-real(ism))*dx
   endif

   if(jsm >= jgrid ) then
    ydis=abs(jsm-jgrid)*dy+(jsnd-real(jsm))*dy 
   else	
	ydis=abs(jsm-jgrid)*dy-(jsnd-real(jsm))*dy  
   endif	
	 
  endif
  
  if( ivar == 2 ) then  !for v point

   if(ism >= igrid ) then
    xdis=abs(ism-igrid)*dx+(isnd-real(ism))*dx
   else	
    xdis=abs(ism-igrid)*dx-(isnd-real(ism))*dx
   endif

   if(jsm >= jgrid ) then   
    ydis=abs(jsm-jgrid)*dy+0.5*dy+(jsnd-real(jsm))*dy
   else
    ydis=abs(jsm-jgrid)*dy-0.5*dy-(jsnd-real(jsm))*dy
   endif 

  endif    

  if( ivar > 2) then ! for non u v point

   if(ism >= igrid ) then
    xdis=abs(ism-igrid)*dx+(isnd-real(ism))*dx
   else	
    xdis=abs(ism-igrid)*dx-(isnd-real(ism))*dx
   endif

   if(jsm >= jgrid ) then
    ydis=abs(jsm-jgrid)*dy+(jsnd-real(jsm))*dy 
   else	
	ydis=abs(jsm-jgrid)*dy-(jsnd-real(jsm))*dy  
   endif

  endif  
     
  hor_distance=sqrt(xdis**2+ydis**2)  
  
! step 2 cal vertial distance
  
  if (ivar == 3 .or. ivar == 4 )  then ! for w and ph point

   zdis= abs(  hgtsnd -  phb(igrid,jgrid,kgrid)/gravity )                                                       
  endif
!!!!for SMOIS & TSLB!!!!
!  if (ivar == 11 .or. ivar == 12 )  then ! for soilm and soilt point

!   zdis= abs(hgtsnd-1.0)   !0.1,0.3,0.6,1.0 VERY simple                                                    
!  endif
!!!!
  if (ivar .ne. 3 .and. ivar .ne. 4 .and. ivar .ne. 11 .and. ivar .ne. 12) then  ! for non w ph point

   zdis= abs(  hgtsnd - (phb(igrid,jgrid,kgrid)+phb(igrid,jgrid,kgrid+1))/2/gravity)
             
  endif  
  
  vert_distance=zdis 
   
! finished cal distance
end subroutine snd_distance  


subroutine localization_cal_snd2(dis_rho,numstat,isnd,jsnd,ksnd,hgtsnd,ism,jsm,ksm,      &
                                 phb,analysis_var_num,ii,jj,kk,ivar,                     &
                                 local_h_dis,local_v_dis,ipoint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine localization_cal_snd2
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'

integer     :: numstat
integer     :: analysis_var_num

integer     :: ism,jsm,ksm
real        :: isnd,jsnd,ksnd,hgtsnd

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
   call snd_distance(isnd,jsnd,ksnd,hgtsnd,ivar,ii,jj,kk,analysis_var_num,  &
                     hor_distance,vert_distance,phb,ism,jsm,ksm)
!!!!for SMOIS & TSLB!!!!
!  if (ivar == 11 .or. ivar == 12 )  then ! for soilm and soilt point


!   zdis= abs(hgtsnd-1.0)   !0.1,0.3,0.6,1.0 VERY simple                                                    
!  endif
!!!!

 !  print*, 'distance, isnd,jsnd,ksnd,ii,jj,kk',hor_distance  ,vert_distance, isnd,jsnd,ksnd,ii,jj,kk
      
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

 !  print*,dis_rho(ipoint)

end subroutine localization_cal_snd2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine localization_cal_lad2(dis_rho,numstat,isnd,jsnd,ksnd,hgtsnd,ism,jsm,ksm,      &
                                 phb,dzs,analysis_var_num,ii,jj,kk,ivar,                 &
                                 local_h_dis,local_v_dis,ipoint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine localization_cal_snd2
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'

integer     :: numstat
integer     :: analysis_var_num

integer     :: ism,jsm,ksm
real        :: isnd,jsnd,ksnd,hgtsnd

real        :: dis_rho(numstat)
real        :: pbh(nx,ny,nz)
real        :: dzs(4)
real        :: hor_distance,vert_distance
real        :: dis_cor
integer     :: ii,jj,kk,mm
integer     :: ivar
integer     :: ipoint
integer     :: itemp

real        :: local_h_dis,local_v_dis


 ! step 1 cal the distance
 !  call snd_distance(isnd,jsnd,ksnd,hgtsnd,ivar,ii,jj,kk,analysis_var_num,  &
 !  hor_distance,vert_distance,phb,ism,jsm,ksm)
 !!!!for SMOIS & TSLB!!!!
 !vert_distance=0.0

  if (ivar == 11 .or. ivar == 12 )  then ! for soilm and soilt point

  if(ism >= ii ) then
    xdis=abs(ism-ii)*dx+(isnd-real(ism))*dx
   else	
    xdis=abs(ism-ii)*dx-(isnd-real(ism))*dx
   endif

   if(jsm >= jj ) then
    ydis=abs(jsm-jj)*dy+(jsnd-real(jsm))*dy 
   else	
	ydis=abs(jsm-jj)*dy-(jsnd-real(jsm))*dy  
   endif

  hor_distance=sqrt(xdis**2+ydis**2)

!   print*,'kk,dzs(kk),hgtsnd',kk,dzs(kk),hgtsnd

  vert_distance = abs(hgtsnd-dzs(kk))   !0.1,0.3,0.6,1.0 VERY simple 
                                                   
!   print*,'vert_distance= abs(hgtsnd-dzs(kk))', vert_distance
  endif !if (ivar == 11 .or. ivar == 12 )

  if (ivar == 13 )  then ! for soilm and soilt point
  if(ism >= ii ) then
    xdis=abs(ism-ii)*dx+(isnd-real(ism))*dx
   else	
    xdis=abs(ism-ii)*dx-(isnd-real(ism))*dx
   endif
   if(jsm >= jj ) then
    ydis=abs(jsm-jj)*dy+(jsnd-real(jsm))*dy 
   else	
	ydis=abs(jsm-jj)*dy-(jsnd-real(jsm))*dy  
   endif

  hor_distance=sqrt(xdis**2+ydis**2)
  vert_distance = 0   ! TSK 2d only 
!   print*,'vert_distance= abs(hgtsnd-dzs(kk))', vert_distance
  endif !if (ivar == 13 ) : 2013/5/12

!   print*, 'distance, isnd,jsnd,ksnd,ii,jj,kk',hor_distance,vert_distance, isnd,jsnd,ksnd,ii,jj,kk
      
! step 2 cal the correlation
   
   if(optlocal  == 1 ) then
     if( hor_distance <= local_h_dis .and. vert_distance <= local_v_dis)  then
       dis_rho(ipoint)=1
!print*,'okokokok'
     else 
       dis_rho(ipoint)=0
     endif
   else if(optlocal  == 2 ) then
     call dis_correlation(hor_distance,vert_distance,local_h_dis,local_v_dis,dis_rho(ipoint)) 
   endif 

!mm=0                             
!if(dis_rho(ipoint).ne. 0)then
!   mm=mm+1
!   print*,'dis_rho(ipoint),ipoint',dis_rho(ipoint),ipoint
!endif

end subroutine localization_cal_lad2


subroutine prepare_info4local_snd(analysis_var_num,ra_data_num,numstat,                          &                               
                                  isnd,jsnd,ksnd,hgtsnd,ism,jsm,ksm,                		 &
                                  phb,dis_rho,gridsgn,gridnum                                    &
                                 )

include 'namelist.inc'
                             
  integer               :: analysis_var_num,ra_data_num,numstat
  integer               :: iobs,ivar
  integer               :: nn  

  real                  :: phb(nx,ny,nz)                 
  real                  :: isnd,jsnd,ksnd,hgtsnd
  integer               :: ism,jsm,ksm

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
!print*,' !!!!!TF FOR SMOIS OR TSLB ? !!!!'  
istart=ism-(int(2*local_h_dis/dx)+1)
if(istart < 1) istart=1
iend=ism+(int(2*local_h_dis/dx)+1)
if(iend > nx) iend=nx
jstart=jsm-(int(2*local_h_dis/dy)+1)
if(jstart < 1) jstart=1
jend=jsm+(int(2*local_h_dis/dy)+1)
if(jend > ny) jend=ny
   
  print*, istart,iend,jstart,jend
   
!***************************************************************  
!       calculate distance correlation coeffient
!*************************************************************** 

 
do k=1,nz-8

do j=jstart,jend
do i=istart,iend
  
  call ijk_to_istat(ivar,i,j,k,istatmp,nx,ny,nz)
! call ijk_to_istat(ivar,i,j,k,istatmp,nx,ny,nz) !this for snd  ::bugs

  call localization_cal_snd2(dis_rho,numstat,isnd,jsnd,ksnd,hgtsnd,ism,jsm,ksm,      &
                             phb,analysis_var_num,i,j,k,ivar,                        &
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

end subroutine prepare_info4local_snd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! for soil layer!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine prepare_info4local_lad(analysis_var_num,ra_data_num,numstat,ivar,                        &                               
                                  isnd,jsnd,ksnd,hgtsnd,ism,jsm,ksm,                		             &
                                  phb,dzs,dis_rho,gridsgn,gridnum                                   &
                                 )

include 'namelist.inc'
                             
  integer               :: analysis_var_num,ra_data_num,numstat
  integer               :: iobs,ivar,nlev
  integer               :: nn  

  real                  :: phb(nx,ny,nz)
  real                  :: dzs(4)              
  real                  :: isnd,jsnd,ksnd,hgtsnd
  integer               :: ism,jsm,ksm

  real                  :: dis_rho(numstat)
  integer               :: istatmp
  integer               :: gridnum,gridsgn(numstat)
  integer               :: ii,jj,kk,itest  
  integer               :: istart,iend,jstart,jend
  character*8           :: timestart,timeend
  real                  :: local_h_dis,local_v_dis
  
  real                  :: local_v_dis_temp,local_h_dis_temp


!!! asume that soil does not have any relationship with above soil layer
!!!
!!! 
!do ivar=11,12

!***************************************************************  
!    mark the area for calculate distance correlation coeffient
!***************************************************************  
   
!   if(ivar >10 ) then

!   endif       
!print*,' !!!!!TF FOR SMOIS OR TSLB ? !!!!'  

!print*, 'ism ,jsm, soil layer:',ism ,jsm,istart,iend,jstart,jend
   
!***************************************************************  
!       calculate distance correlation coeffient
!*************************************************************** 

if( ivar > 10 .and. ivar < 13 )then

istart=ism-(int(2*local_dis_soil/dx)+1)
if(istart < 1) istart=1
iend=ism+(int(2*local_dis_soil/dx)+1)
if(iend > nx) iend=nx
jstart=jsm-(int(2*local_dis_soil/dy)+1)
if(jstart < 1) jstart=1
jend=jsm+(int(2*local_dis_soil/dy)+1)
if(jend > ny) jend=ny

local_h_dis=local_dis_soil
local_v_dis=verti_dis_soil   

nlev=4
!else
!nlev=nz-8
!endif

 
do k=1,nlev

do j=jstart,jend
do i=istart,iend
  
  call ladijk_to_istat(ivar,i,j,k,istatmp,nx,ny,nz,nlev)
! print*,i,j,k,istatmp
! call ijk_to_istat(ivar,i,j,k,istatmp,nx,ny,nz) !this for snd  ::bugs

  call localization_cal_lad2(dis_rho,numstat,isnd,jsnd,ksnd,hgtsnd,ism,jsm,ksm,      &
                             phb,dzs,analysis_var_num,i,j,k,ivar,                        &
                             local_dis_soil,verti_dis_soil,istatmp)							     

! print*,'dis_rho,dzs(k),ism,jsm,ksm',dis_rho,dzs(kk),ism,jsm,ksm,nlev
! stop                                                                                                   

enddo    
enddo
enddo                      

endif  ! if ivar > 10

! print*,'dzs(k),ism,jsm,ksm,nlev,ivar',dzs(k),ism,jsm,ksm,nlev,ivar,numstat,istatmp
!stop    
!print*,'ok!!!'

!enddo  


if( ivar == 13 )then   !!! for TSK
 local_h_dis=local_dis_soil
 local_v_dis=verti_dis_soil
!!
istart=ism-(int(2*local_dis_soil/dx)+1)
if(istart < 1) istart=1
iend=ism+(int(2*local_dis_soil/dx)+1)
if(iend > nx) iend=nx
jstart=jsm-(int(2*local_dis_soil/dy)+1)
if(jstart < 1) jstart=1
jend=jsm+(int(2*local_dis_soil/dy)+1)
if(jend > ny) jend=ny
!!
 nlev=4

do j=jstart,jend
do i=istart,iend
  
  call ladijk_to_istat(ivar,i,j,k,istatmp,nx,ny,nz,nlev)


  call localization_cal_lad2(dis_rho,numstat,isnd,jsnd,ksnd,hgtsnd,ism,jsm,ksm,      &
                             phb,dzs,analysis_var_num,i,j,k,ivar,                    &
                             local_dis_soil,verti_dis_soil,istatmp)							     
                                                                                                   
enddo    
enddo                      

endif  ! if ivar == 13
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

!print*,'dzs(kk),ism,jsm,ksm,nlev,ivar',dzs(kk),ism,jsm,ksm,nlev,ivar,numstat,istatmp
!print*,'gridnum',gridnum
!stop

end subroutine prepare_info4local_lad



