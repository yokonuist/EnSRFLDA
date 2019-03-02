subroutine innoation_rmse_rv(analysis_var_num,ra_data_num,numstat,                          &
                             xstat,                                                         &
                             rv,ir,jr,hr,evl,azimuth,                                       &
                             xradar,yradar,hradar,lev_num,                                  &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                             rdnw,rdn,mub,mu,iradar,ireg,                                   &
                             rmse_rv,spd_rv                                                 &
                             )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine innoation_rmse_rv
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
implicit none                            
include 'namelist.inc'

  integer               :: analysis_var_num,ra_data_num,numstat
  integer               :: iobs,ivar
  integer               :: nn
  
  integer               :: istatnum1,istatnum2,istatnum3,istatnum4,      &
                           istatnum6,istatnum5,istatnum7,istatnum8,      &
                           istatnum9,istatnum10,istatnum11,istatnum12,   &
                           istatnum13,istatnum14,istatnum15,istatnum16,  &
                           istatnum17,istatnum18,istatnum19,istatnum20   
                           
  real                  :: uuw,uue,ulw,ule,                   &
                           vus,vun,vls,vln,                   &
                           wu,wl
  real                  :: qru,qrl,qsu,qsl,qgru,qgrl,         &
                           rhou,rhol,tcu,tcl   
  real                  :: pu,pl
  
  real                  :: ur,vr,wr                         
  real                  :: qrr,qsr,qgrr,tcr,rhor                      
  real                  :: ref_x,ref_r,ref_s,ref_h
  real                  :: uctb,vctb,wctb
  
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat,hgt
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real                  :: xstat(numstat,ens_num)
  
  real,parameter        :: gravity=9.81                         
  real                  :: rv(ra_data_num),hr(ra_data_num)
  integer               :: ir(ra_data_num),jr(ra_data_num)
  real                  :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num  
  real                  :: azimuth(lev_num)
  real                  :: azimuth_temp
  real                  :: evl(ra_data_num)
  real                  :: evl_temp
  real                  :: hor_dis
  integer               :: krm,krw
                           
  real,allocatable      :: hx(:),innovation(:)  
  real                  :: hxmean  
  
  real,allocatable      :: spd(:) 
  real                  :: rmse_rv,spd_rv
  real                  :: innovationmean
  real                  :: kalman_alpha
  
  integer               :: temp1,temp2,temp3
                           
  real                  :: tem1,tem2,tem3,tem4,tem5
  real                  :: ratio_r2s_rv
  
  real                  :: rv_grid(nx,ny,lev_num)
  integer               :: ilv
  character(len=2)      :: lab_radar_num    
  character(len=3)      :: time_lab
  real                  :: evl_record(lev_num)  
  integer               :: i,j,k,ireg
  integer               :: iradar
  
  allocate(spd(ens_num))
  allocate(hx(ens_num)) 
  allocate(innovation(ens_num)) 
  
  rmse_rv=0
  spd_rv =0
  temp1=0
  spd=0

  temp1=int((time_num)/100     )-int((time_num)/1000    )*10  
  temp2=int((time_num)/10      )-int((time_num)/100     )*10 
  temp3=int((time_num)/1       )-int((time_num)/10      )*10
  time_lab=char(48+temp1)//char(48+temp2)//char(48+temp3)

  rv_grid = miss_data
  open(1101,file='../../obs_real/angleinfo2_'//time_lab//'.dat')
  do i=1,lev_num
   read(1101,*) evl_record(i)
  enddo
  close(1101)


! counting the number of valid radial velocity data
  do iobs=1,ra_data_num
  if(rv(iobs) .ne. miss_data) then ! for valid radial velocity 
  temp1=temp1+1  
  endif          ! for valid radial velocity    
  enddo
    
  do iobs=1,ra_data_num
  
  azimuth_temp=azimuth(iobs)*3.1415926/180
  evl_temp=evl(iobs)*3.1415926/180
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr(iobs),phb,krw,krm)   
  
  if(rv(iobs) .ne. miss_data .and. hr(iobs) > 100) then ! for valid radial velocity  
  hxmean=0 
  
  do i=1,ens_num   
  
  ivar=1   ! for u 
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum1 ,nx,ny,nz)  ! for uuw
  call ijk_to_istat(ivar,ir(iobs)+1,jr(iobs)  ,krm+1,istatnum2 ,nx,ny,nz)  ! for uue
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum3 ,nx,ny,nz)  ! for ulw
  call ijk_to_istat(ivar,ir(iobs)+1,jr(iobs)  ,krm  ,istatnum4 ,nx,ny,nz)  ! for ule
  ivar=2   ! for v 
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum5 ,nx,ny,nz)  ! for vus
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)+1,krm+1,istatnum6 ,nx,ny,nz)  ! for vun
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum7 ,nx,ny,nz)  ! for vls
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)+1,krm  ,istatnum8 ,nx,ny,nz)  ! for vln  
  ivar=3   ! for w
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krw+1,istatnum9 ,nx,ny,nz)  ! for wu
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krw  ,istatnum10,nx,ny,nz)  ! for wl
  
! prepare for cal vtm using method copy from arps(based on NASA )
  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw  ,istatnum11,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum12,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum13,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum14,nx,ny,nz)  
  
  
  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                                            &
               xstat(istatnum14,i),                                                                         &
               xstat(istatnum11,i),xstat(istatnum12,i),                                                     &
               rdnw(krw),p_top,znu(krm),                                                                    &
               xstat(istatnum13,i),                                                                         &
               nx,ny,nz,rhol,pl                                                                             &                  
              )  
  call cal_tc(xstat(istatnum13,i),pl,tcl,nx,ny,nz) 


  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum11,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+2,istatnum12,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum13,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum14,nx,ny,nz)  

  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                                            &
               xstat(istatnum14,i),                                                                         &
               xstat(istatnum11,i),xstat(istatnum12,i),                                                     &
               rdnw(krw),p_top,znu(krm),                                                                    &
               xstat(istatnum13,i),                                                                         &
               nx,ny,nz,rhou,pu                                                                             &                  
              )  
  call cal_tc(xstat(istatnum13,i),pu,tcu,nx,ny,nz)

                                             
  ivar=7   ! for qr 
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum11 ,nx,ny,nz)  ! for qru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum12 ,nx,ny,nz)  ! for qrl
  ivar=9   ! for qs
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum13 ,nx,ny,nz)  ! for qsu
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum14 ,nx,ny,nz)  ! for qsl  
  ivar=10  ! for qgr
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum15 ,nx,ny,nz)  ! for qgru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum16 ,nx,ny,nz)  ! for qgrl


  uuw=xstat(istatnum1,i)
  uue=xstat(istatnum2,i)
  ulw=xstat(istatnum3,i)
  ule=xstat(istatnum4,i)
  vus=xstat(istatnum5,i)
  vun=xstat(istatnum6,i)
  vls=xstat(istatnum7,i)
  vln=xstat(istatnum8,i)
  wu =xstat(istatnum9,i)
  wl =xstat(istatnum10,i)
  
  
      call interpolation_radar_wind(nx,ny,nz,ir(iobs),jr(iobs),krw,krm,hr(iobs),phb,     &
                                    uuw,uue,ulw,ule,                                     &
                                    vus,vun,vls,vln,                                     &
                                    wu,wl,                                               &
                                    ur,vr,wr                                             &
                                   ) 

  qru =xstat(istatnum11 ,i)
  qrl =xstat(istatnum12 ,i)
  qsu =xstat(istatnum13 ,i)
  qsl =xstat(istatnum14 ,i)
  qgru=xstat(istatnum15 ,i)
  qgrl=xstat(istatnum16 ,i)
 
                                       
      call interpolation_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr(iobs),phb,     &
                                    qru,qrl,qsu,qsl,qgru,qgrl,                       &
                                    rhou,rhol,tcu,tcl,                               &
                                    qrr,qsr,qgrr,rhor,tcr                            &
                                    )
  
  call ref_operator(qrr,qsr,qgrr,ref_x,ref_r,ref_s,ref_h,rhor,tcr)                                   
                                               
  call rv_operator(ur,vr,wr,hx(i),azimuth_temp,evl_temp,ref_x,ref_r,ref_s,ref_h,   &
                   qrr,qsr,qgrr,rhor,uctb,vctb,wctb                                &
                   ) 
  
  hxmean=hxmean+hx(i)/ens_num

  enddo        !  for ensemble loop
  
  !-----------------------------------
  ! record the hxmean to a new array
  !-----------------------------------
  call cal_lv_num_real(evl(iobs),evl_record,ilv,lev_num)
  rv_grid(ir(iobs),jr(iobs),ilv)=hxmean
  !-----------------------------------  
!  print*,ir(iobs),jr(iobs),ilv,'position'
!  print*,rv_grid(ir(iobs),jr(iobs),ilv),hxmean
 
  
  call innovation_cal(innovationmean,rv(iobs),hxmean) 
  
  do i=1,ens_num
  spd(i)=spd(i) + (hx(i)-hxmean)**2/temp1
  enddo
  
  rmse_rv = rmse_rv + innovationmean**2/temp1
  


!  temp1=temp1+1  
  endif          ! for valid radial velocity  
  
  enddo
!  rmse_rv=sqrt(rmse_rv/temp1)
  rmse_rv=sqrt(rmse_rv)
  
  do i=1,ens_num
!  spd_rv=spd_rv+sqrt(spd(i)/temp1)/ens_num
  spd_rv=spd_rv+sqrt(spd(i))/ens_num
  enddo

  call ratio_rmse2spd(rmse_rv,spd_rv,ratio_r2s_rv)
  
  print*,'rmse_rv is : ',rmse_rv
  print*,'spd_rv is : ',spd_rv
  print*,'the ratio of rmse_rv to spd_rv is:',ratio_r2s_rv
  
  
  !-----------------------------------
  ! record the hxmean to a new array
  !-----------------------------------  
  temp1=int((time_num)/100     )-int((time_num)/1000    )*10  
  temp2=int((time_num)/10      )-int((time_num)/100     )*10 
  temp3=int((time_num)/1       )-int((time_num)/10      )*10
  time_lab=char(48+temp1)//char(48+temp2)//char(48+temp3) 

  temp1=int(iradar/10      )-int(iradar/100     )*10  
  temp2=int(iradar/1       )-int(iradar/10      )*10
  lab_radar_num = char(48+temp1)//char(48+temp2) 
  
  
  if( ireg == 0 ) then
      open(1100,file='./rv_grid_real_'//lab_radar_num//'_'//time_lab//'_before.dat',form='binary')  !
  else if( ireg == 1 ) then
      open(1100,file='./rv_grid_real_'//lab_radar_num//'_'//time_lab//'_after.dat',form='binary')    !
  endif
  do k=1,lev_num	
  	do j=1,ny
  	  do i=1,nx
       write(1100) rv_grid(i,j,k)            !,*,i,j,k
!       print*,rv_grid(i,j,k)
!       if(rv_grid(i,j,k) /= miss_data)  then
!        pause
!       endif 
      enddo
    enddo
  enddo       
  close(1100)    
  
!stop
end subroutine innoation_rmse_rv






subroutine innoation_rmse_rf(analysis_var_num,ra_data_num,numstat,                          &
                             xstat,                                                         &
                             rf,ir,jr,hr,evl,azimuth,                                       &
                             xradar,yradar,hradar,lev_num,                                  &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                             rdnw,rdn,mub,mu,iradar,ireg,                                   &
                             rmse_rf,spd_rf                                                 &
                             )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine innoation_rmse_rf
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
implicit none                              
include 'namelist.inc'

  integer               :: analysis_var_num,ra_data_num,numstat
  integer               :: iobs,ivar
  integer               :: nn
  
  integer               :: istatnum1,istatnum2,istatnum3,istatnum4,      &
                           istatnum6,istatnum5,istatnum7,istatnum8,      &
                           istatnum9,istatnum10,istatnum11,istatnum12,   &
                           istatnum13,istatnum14,istatnum15,istatnum16,  &
                           istatnum17,istatnum18,istatnum19,istatnum20   
                           
  real                  :: uuw,uue,ulw,ule,                   &
                           vus,vun,vls,vln,                   &
                           wu,wl
  real                  :: qru,qrl,qsu,qsl,qgru,qgrl,         &
                           rhou,rhol,tcu,tcl   
  real                  :: pu,pl
  
  real                  :: ur,vr,wr                         
  real                  :: qrr,qsr,qgrr,tcr,rhor                      
  real                  :: ref_x,ref_r,ref_s,ref_h
  
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat,hgt
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real                  :: xstat(numstat,ens_num)
  

  real,parameter        :: gravity=9.81                         
  real                  :: rf(ra_data_num),hr(ra_data_num)
  integer               :: ir(ra_data_num),jr(ra_data_num)
  integer               :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num  
  real                  :: azimuth(lev_num)
  real                  :: azimuth_temp
  real                  :: evl(ra_data_num)
  real                  :: evl_temp
  real                  :: hor_dis
  integer               :: krm,krw
                           
  real,allocatable      :: hx(:),innovation(:)  
  real                  :: hxmean  
  
  real,allocatable      :: spd(:) 
  real                  :: rmse_rf,spd_rf
  real                  :: innovationmean
  real                  :: kalman_alpha
  
  integer               :: temp1,temp2,temp3
                           
  real                  :: tem1,tem2,tem3,tem4,tem5
  real                  :: ratio_r2s_rf

  real                  :: rf_grid(nx,ny,lev_num)
  integer               :: ilv
  character(len=2)      :: lab_radar_num    
  character(len=3)      :: time_lab
  real                  :: evl_record(lev_num) 
  integer               :: i,j,k,ireg
  integer               :: iradar
    
  allocate(spd(ens_num))
  allocate(hx(ens_num)) 
  allocate(innovation(ens_num)) 
  
  rmse_rf=0
  spd_rf =0
  temp1=0
  spd=0
  
  temp1=int((time_num)/100     )-int((time_num)/1000    )*10  
  temp2=int((time_num)/10      )-int((time_num)/100     )*10 
  temp3=int((time_num)/1       )-int((time_num)/10      )*10
  time_lab=char(48+temp1)//char(48+temp2)//char(48+temp3)

  rf_grid = miss_data
  open(1101,file='../../obs_real/angleinfo2_'//time_lab//'.dat')
  do i=1,lev_num
   read(1101,*) evl_record(i)
  enddo
  close(1101)
  
  
! counting the number of valid reflectivity data
  do iobs=1,ra_data_num
  if(rf(iobs) .ne. miss_data)  then ! for valid reflectivity
  temp1=temp1+1  
  endif          ! for valid reflectivity   
  enddo
  
  do iobs=1,ra_data_num

  azimuth_temp=azimuth(iobs)*3.1415926/180
  evl_temp=evl(iobs)*3.1415926/180
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr(iobs),phb,krw,krm)   
    
  if(rf(iobs) .ne. miss_data .and. hr(iobs) > 100)  then ! for valid reflectivity
     
  hxmean=0
  
  do i=1,ens_num
  
  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw  ,istatnum1,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum2,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum3,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum4,nx,ny,nz)  
    
  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                                                     &
               xstat(istatnum4,i),                                                                                   &
               xstat(istatnum1,i),xstat(istatnum2,i),                                                                &
               rdnw(krw),p_top,znu(krm),                                                                             &
               xstat(istatnum3,i),                                                                                   &
               nx,ny,nz,rhol,pl                                                                                      &                  
              )  
  call cal_tc(xstat(istatnum3,i),pl,tcl,nx,ny,nz) 
  
  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum1,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+2,istatnum2,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum3,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum4,nx,ny,nz)  

  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                                                     &
               xstat(istatnum4,i),                                                                                   &
               xstat(istatnum1,i),xstat(istatnum2,i),                                                                &
               rdnw(krw),p_top,znu(krm),                                                                             &
               xstat(istatnum3,i),                                                                                   &
               nx,ny,nz,rhou,pu                                                                                      &                  
              )  
  call cal_tc(xstat(istatnum3,i),pu,tcu,nx,ny,nz) 

  
!*********************for pert reflectivity  hx cal*********************** 
                     
  ivar=7   ! for qr 
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum1 ,nx,ny,nz)  ! for qru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum2 ,nx,ny,nz)  ! for qrl
  ivar=9   ! for qs
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum3 ,nx,ny,nz)  ! for qsu
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum4 ,nx,ny,nz)  ! for qsl  
  ivar=10  ! for qgr
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum5 ,nx,ny,nz)  ! for qgru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum6 ,nx,ny,nz)  ! for qgrl
  
  qru =xstat(istatnum1 ,i)  
  qrl =xstat(istatnum2 ,i)
  qsu =xstat(istatnum3 ,i)
  qsl =xstat(istatnum4 ,i)
  qgru=xstat(istatnum5 ,i)
  qgrl=xstat(istatnum6 ,i)

                                  
  call interpolation_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr(iobs),phb,     &
                                qru,qrl,qsu,qsl,qgru,qgrl,                       &
                                rhou,rhol,tcu,tcl,                               &
                                qrr,qsr,qgrr,rhor,tcr                            &
                                )

 
  call ref_operator(qrr,qsr,qgrr,hx(i),ref_r,ref_s,ref_h,rhor,tcr)  
  
  hxmean=hxmean+hx(i)/ens_num  

  enddo

  !-----------------------------------
  ! record the hxmean to a new array
  !-----------------------------------
  call cal_lv_num_real(evl(iobs),evl_record,ilv,lev_num)
  rf_grid(ir(iobs),jr(iobs),ilv)=hxmean
  !-----------------------------------  

  call innovation_cal(innovationmean,rf(iobs),hxmean) 
  
  rmse_rf = rmse_rf + innovationmean**2/temp1
  
  do i=1,ens_num
  spd(i)=spd(i) + (hx(i)-hxmean)**2/temp1
  enddo
  
!  temp1=temp1+1  
  endif          ! for valid reflectivity  

  enddo
!  rmse_rf=sqrt(rmse_rf/temp1)
  rmse_rf=sqrt(rmse_rf)
  
  do i=1,ens_num
!  spd_rf=spd_rf+sqrt(spd(i)/temp1)/ens_num
  spd_rf=spd_rf+sqrt(spd(i))/ens_num
  enddo

  call ratio_rmse2spd(rmse_rf,spd_rf,ratio_r2s_rf)

  print*,'rmse_rf is : ',rmse_rf
  print*,'spd_rf is : ',spd_rf 
  print*,'the ratio of rmse_rf to spd_rf is:',ratio_r2s_rf
  
  !-----------------------------------
  ! record the hxmean to a new array
  !-----------------------------------  
  temp1=int((time_num)/100     )-int((time_num)/1000    )*10  
  temp2=int((time_num)/10      )-int((time_num)/100     )*10 
  temp3=int((time_num)/1       )-int((time_num)/10      )*10
  time_lab=char(48+temp1)//char(48+temp2)//char(48+temp3)

  temp1=int(iradar/10      )-int(iradar/100     )*10  
  temp2=int(iradar/1       )-int(iradar/10      )*10
  lab_radar_num = char(48+temp1)//char(48+temp2) 
  
  if( ireg == 0 ) then
      open(1100,file='./rf_grid_real_'//lab_radar_num//'_'//time_lab//'_before.dat',form='binary')  !
  else if( ireg == 1 ) then
      open(1100,file='./rf_grid_real_'//lab_radar_num//'_'//time_lab//'_after.dat',form='binary')   !
  endif
  do k=1,lev_num	
  	do j=1,ny
  	  do i=1,nx
       write(1100) rf_grid(i,j,k)                       !,*,i,j,k
!       print*,rv_grid(i,j,k)
!       if(rf_grid(i,j,k) /= miss_data)  then
!        pause
!       endif 
      enddo
    enddo
  enddo    
  close(1100)      
  
end subroutine innoation_rmse_rf


subroutine ratio_rmse2spd(rmse,spd,ratio_r2s)

real      :: rmse
real      :: spd
real      :: ratio_r2s

ratio_r2s=rmse/spd


end subroutine ratio_rmse2spd


subroutine output_radio_r2s_rvrf(rmse_rv,spd_rv,rmse_rf,spd_rf,iradar,lab_time_num,ireg)

integer               :: iradar
real                  :: rmse_rv,spd_rv,rmse_rf,spd_rf
character(len=2)      :: lab_radar_num
character(len=2)      :: lab_time_num
integer               :: temp1,temp2,temp3
integer               :: ireg

temp1=int(iradar/10      )-int(iradar/100     )*10  
temp2=int(iradar/1       )-int(iradar/10      )*10
lab_radar_num = char(48+temp1)//char(48+temp2) 
  
  if(ireg == 0 ) then 
  open(8553,file='radio_r2s_rvrf_'//lab_radar_num//'_'//lab_time_num//'_before_assim.dat')
  else if(ireg == 1 ) then 
  open(8553,file='radio_r2s_rvrf_'//lab_radar_num//'_'//lab_time_num//'_after_assim.dat')
  endif
  
  write(8553,'(2A8)') 'RV','RF'
  write(8553,'(2F8.3)') rmse_rv/spd_rv,rmse_rf/spd_rf
  write(8553,'(2A8)') 'rmse_rv','spd_rv'
  write(8553,'(2F8.3)') rmse_rv,spd_rv
  write(8553,'(2A8)') 'rmse_rf','spd_rf'
  write(8553,'(2F8.3)') rmse_rf,spd_rf  
  close(8553) 
  
  
end subroutine output_radio_r2s_rvrf

subroutine cal_lv_num_real(evlr,evl,ilv,lev_num)

implicit none
integer       :: lev_num
real          :: evlr
real          :: evl(lev_num)
integer       :: ilv

real          :: temr01
integer       :: i,j,k

!print*,'******'
do i=1,lev_num
!  print*,evlr,evl(i)
  if(abs(evlr-evl(i))<0.01) then
		ilv=i
  endif

enddo

end subroutine cal_lv_num_real