subroutine assimilate_radar_cal(analysis_var_num,ra_data_num,numstat,                          &
                                xstat,xstat_mean,                                              &
                                rf,rv,ir,jr,hr_rf,hr_rv,evl,azimuth,                           &
                                xradar,yradar,hradar,lev_num,                                  &
                                xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                                rdnw,rdn,mub,mu                                                &
                               ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine assimilate_radar_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.03.27   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                                 
include 'namelist.inc'
  
  integer               :: analysis_var_num,ra_data_num,numstat
  integer               :: iobs,ivar
  integer               :: nn
  
  integer               :: istatnum1,istatnum2,istatnum3,istatnum4,      &
                           istatnum6,istatnum5,istatnum7,istatnum8,      &
                           istatnum9,istatnum10,istatnum11,istatnum12,   &
                           istatnum13,istatnum14,istatnum15,istatnum16,  &
                           istatnum17,istatnum18,istatnum19,istatnum20   
                           
  real*8                  :: uuw,uue,ulw,ule,                   &
                           vus,vun,vls,vln,                   &
                           wu,wl
  real*8                  :: qru,qrl,qsu,qsl,qgru,qgrl,         &
                           rhou,rhol,tcu,tcl   
  real*8                 :: pu,pl
  
  real*8                  :: ur,vr,wr                         
  real*8                  :: qrr,qsr,qgrr,tcr,rhor                      
  real                  :: ref_x,ref_r,ref_s,ref_h
  real                  :: uctb,vctb,wctb  
  real*8,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real*8                  :: xstat(numstat,ens_num),xstat_mean(numstat)
  
  real,parameter        :: gravity=9.81                         
  real                  :: rf(ra_data_num),rv(ra_data_num),hr_rv(ra_data_num),hr_rf(ra_data_num)
  integer               :: ir(ra_data_num),jr(ra_data_num)
  real                  :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num  
  real                  :: azimuth(ra_data_num)
  real                  :: azimuth_temp
  real                  :: evl(ra_data_num)
  real                  :: evl_temp
  real                  :: hor_dis
  integer               :: krm,krw
                           
  real,allocatable      :: hx_pert(:),hx(:),innovation(:)  
  real                  :: hxmean  

  real,allocatable      :: zzh(:),kalman(:)
  real                  :: hph,innovationmean
  real                  :: kalman_alpha
  
  real                  :: hor_distance,vert_distance,dis_cor
  
  real                  :: temp1(ens_num),temp2(ens_num),temp3(ens_num),temp4(ens_num), &
                           temp5(ens_num),temp6(ens_num),temp7(ens_num),temp8(ens_num), &
                           temp9(ens_num),temp10(ens_num)
                           
  real                  :: tem1,tem2,tem3,tem4,tem5
  real                  :: time1,time2,use_time
  integer               :: istatmp,imp,jmp,kmp
  integer               :: itemp1,itemp2,itemp3
  
!*****************************************************************
  real,allocatable      :: xstat_sub(:,:),xstat_mean_sub(:)
  real,allocatable      :: dis_rho(:)
  integer,allocatable   :: gridsgn(:)      
  integer               :: gridnum 
  real,allocatable      :: temp_dis_rho(:)
  integer,allocatable   :: temp_gridsgn(:)
!*****************************************************************


!  allocate(zzh(numstat))
!  allocate(kalman(numstat))


print*,'start the assimilation of sim radar obs'
print*,ra_data_num


  do iobs=1,ra_data_num
   
!*********************prepare for HX cal*********************************
  azimuth_temp=azimuth(iobs)*3.1415926/180
  evl_temp=evl(iobs)*3.1415926/180
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr_rv(iobs),phb,krw,krm)  
!************************************************************************

  
if(assim_rv == 1 .and. rv(iobs) .ne. miss_data .and. hr_rv(iobs) > 100) then ! for radial velocity

!***************************************************************
!           prepare information for local
!***************************************************************  
  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))

do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 

!  print*,'check whether the allocate is normal'
!  print*,temp_dis_rho(1),temp_dis_rho(1)
!  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)  
  
call prepare_info4local(analysis_var_num,ra_data_num,numstat,                          &                               
                        ir(iobs),jr(iobs),hr_rv(iobs),xradar,yradar,hradar,phb,        &
                        temp_dis_rho,temp_gridsgn,gridnum                              &
                        )
                                           

  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
    
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))

  print*,'allocate finished'
  
  
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo

  
!  print*,'check whether the allocate is normal'
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum)
  
  nn=1
  do i=1,numstat
  if(temp_dis_rho(i) > 0.05 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
     do j=1,ens_num
        xstat_sub(nn,j)=xstat(gridsgn(nn),j) 
      enddo
   nn=nn+1
  endif
  enddo 

  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)
  print*,'there are ',gridnum,'grid being influented by obs'
!************************************************************************
!***********************for radial velocity HX cal***********************
!************************************************************************    
  hxmean=0
  hx_chi_mean=0
  
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
                                                                   
  ivar=7   ! for qr 
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum11 ,nx,ny,nz)  ! for qru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum12 ,nx,ny,nz)  ! for qrl
  ivar=9   ! for qs
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum13 ,nx,ny,nz)  ! for qsu
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum14 ,nx,ny,nz)  ! for qsl  
  ivar=10  ! for qgr
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum15 ,nx,ny,nz)  ! for qgru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum16 ,nx,ny,nz)  ! for qgrl

!*****************************************************************      
  do i=1,ens_num   

! prepare for cal vtm using method copy from arps(based on NASA )
  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw  ,istatnum17,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum18,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum19,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum20,nx,ny,nz)  
    
  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                    &
               xstat(istatnum20,i),                                                 &
               xstat(istatnum17,i),xstat(istatnum18,i),                             &
               rdnw(krw),p_top,znu(krm),                                            &
               xstat(istatnum19,i),                                                 &
               nx,ny,nz,rhol,pl                                                     &                  
              )  
  call cal_tc(xstat(istatnum13,i),pl,tcl,nx,ny,nz) 
  
  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum17,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+2,istatnum18,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum19,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum20,nx,ny,nz)  

  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                    &
               xstat(istatnum20,i),                                                 &
               xstat(istatnum17,i),xstat(istatnum18,i),                             &
               rdnw(krw),p_top,znu(krm),                                            &
               xstat(istatnum19,i),                                                 &
               nx,ny,nz,rhou,pu                                                     &                  
              )  
  call cal_tc(xstat(istatnum13,i),pu,tcu,nx,ny,nz)

  
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
  
      call interpolation_radar_wind(nx,ny,nz,ir(iobs),jr(iobs),krw,krm,hr_rv(iobs),phb,  &
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
 
                                       
      call interpolation_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr_rv(iobs),phb,  &
                                    qru,qrl,qsu,qsl,qgru,qgrl,                       &
                                    rhou,rhol,tcu,tcl,                               &
                                    qrr,qsr,qgrr,rhor,tcr                            &
                                    )


  call ref_operator(qrr,qsr,qgrr,ref_x,ref_r,ref_s,ref_h,rhor,tcr)                                   
                                                                        
  call rv_operator(ur,vr,wr,hx(i),azimuth_temp,evl_temp,ref_x,ref_r,ref_s,ref_h,   &
                   qrr,qsr,qgrr,rhor,uctb,vctb,wctb                                &
                  )
  
  hx(i)=hx(i)
  
  hxmean=hxmean+hx(i)/ens_num

  enddo
  
  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean  
  enddo

     
!************************************************************************   
!************************************************************************  

! cal the innovation vecter    
  call innovation_cal(innovationmean,rv(iobs),hxmean)
  
!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 
  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
    
  
! quality control if innovationmean is great than threhold , assimilation will not be done
  if( abs(innovationmean) > ob_rv_threshold ) goto 101    
  print*,'innovationmean,rv(iobs),hxmean',innovationmean,rv(iobs),hxmean    

  
! cal the ZZH 
!  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)
  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)

  
if(enable_rm_correlation_rv == 1 ) then   
! If the correlation between RV and state variation is not good or even lead to a worse result
! it should be remove.currently,it is found that rv will worsen the 1st level for all var


!  do ivar=1,analysis_var_num   
!  call rm_correlate(zzh,numstat,ivar,nx,ny,nz,1)
!  enddo  


  call rm_correlate2(zzh,gridnum,gridsgn,nx,ny,nz,1) 
 
endif

! cal localization
  if(enable_local ==1) then
  

!do istatmp=1,numstat 
!       call localization_cal_radar(zzh,numstat,ir(iobs),jr(iobs),hr_rv(iobs),phb, &
!                                   analysis_var_num,istatmp)                                      
!enddo                          


  call simple_local(zzh,gridnum,dis_rho)
      
  endif  
  
! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
   
! cal kalman maxtria                                                                          
!  call kalman_cal(zzh,hph,obs_err_rv,numstat,kalman)
  call kalman_cal(zzh,hph,obs_err_rv,gridnum,kalman)                                               
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_rv) 

  
!  print*,'ir(iobs),jr(iobs),krm,krw',ir(iobs),jr(iobs),krm,krw
!  print*,'zzh',zzh(itemp1),zzh(itemp2),zzh(itemp3)
!  print*, 'kalman_alpha,hph',kalman_alpha,hph

! cal analysis field 
!  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,   &
!                      xstat,xstat_mean,kalman,hx_pert                &
!                     )  
                     
  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 
                     
! finished innovation for a rv obs
  print*,'finished innovation for a rv obs'

101 continue 
!  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)  
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman)  
  deallocate(hx) 
  deallocate(hx_pert)  

print*,'assimilate the obs at the',iobs,'th point'     
  endif       ! for radial velocity
                               
  enddo  !  end loop for obs      


  do iobs=1,ra_data_num
   
!*********************prepare for HX cal*********************************
  azimuth_temp=azimuth(iobs)*3.1415926/180
  evl_temp=evl(iobs)*3.1415926/180
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr_rf(iobs),phb,krw,krm)  
!************************************************************************
!print*,'azimuth,evlr',azimuth*180/3.1415926,evlr*180/3.1415926  
  
if(assim_rf == 1 .and. rf(iobs) .ne. miss_data .and. rf(iobs) >= 10 .and. hr_rf(iobs) > 100)  then ! for reflectivity

!***************************************************************
!           prepare information for local
!***************************************************************  
  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))


do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 


!  print*,'check whether the allocate is normal'
!  print*,temp_dis_rho(1),temp_dis_rho(1)
!  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)
      
call prepare_info4local(analysis_var_num,ra_data_num,numstat,                          &                               
                        ir(iobs),jr(iobs),hr_rf(iobs),xradar,yradar,hradar,phb,        &
                        temp_dis_rho,temp_gridsgn,gridnum                              &
                        )
                          

  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))  
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))  
  print*,'allocate finished'
  
   
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo
 
 
!  print*,'check whether the allocate is normal'
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum) 
  
  nn=1
  do i=1,numstat
  if(temp_gridsgn(i) > 0 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
     do j=1,ens_num
        xstat_sub(nn,j)=xstat(gridsgn(nn),j) 
      enddo
   nn=nn+1
  endif
  enddo 
   
  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)  
!************************************************************************
!***********************for reflectivity HX cal**************************
!************************************************************************     
  hxmean=0
                     
  ivar=7   ! for qr 
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum1 ,nx,ny,nz)  ! for qru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum2 ,nx,ny,nz)  ! for qrl
  ivar=9   ! for qs
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum3 ,nx,ny,nz)  ! for qsu
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum4 ,nx,ny,nz)  ! for qsl  
  ivar=10  ! for qgr
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm+1,istatnum5 ,nx,ny,nz)  ! for qgru
  call ijk_to_istat(ivar,ir(iobs)  ,jr(iobs)  ,krm  ,istatnum6 ,nx,ny,nz)  ! for qgrl
  
!************************************************************************    

  do i=1,ens_num
  
  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw  ,istatnum17,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum18,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum19,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm,istatnum20,nx,ny,nz)  
    
  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                              &
               xstat(istatnum20,i),                                                           &
               xstat(istatnum17,i),xstat(istatnum18,i),                                       &
               rdnw(krw),p_top,znu(krm),                                                      &
               xstat(istatnum19,i),                                                           &
               nx,ny,nz,rhol,pl                                                               &                  
              )  
  call cal_tc(xstat(istatnum3,i),pl,tcl,nx,ny,nz) 
  
  ivar=4   ! for ph 
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+1,istatnum17,nx,ny,nz)
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krw+2,istatnum18,nx,ny,nz)
  ivar=5   ! for t
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum19,nx,ny,nz)
  ivar=6   ! for qv
  call ijk_to_istat(ivar,ir(iobs),jr(iobs),krm+1,istatnum20,nx,ny,nz)  

  call cal_rho(mub(ir(iobs),jr(iobs),i),mu(ir(iobs),jr(iobs),i),                              &
               xstat(istatnum20,i),                                                           &
               xstat(istatnum17,i),xstat(istatnum18,i),                                       &
               rdnw(krw),p_top,znu(krm),                                                      &
               xstat(istatnum19,i),                                                           &
               nx,ny,nz,rhou,pu                                                               &                  
              )  
  call cal_tc(xstat(istatnum3,i),pu,tcu,nx,ny,nz) 


  qru =xstat(istatnum1 ,i)  
  qrl =xstat(istatnum2 ,i)
  qsu =xstat(istatnum3 ,i)
  qsl =xstat(istatnum4 ,i)
  qgru=xstat(istatnum5 ,i)
  qgrl=xstat(istatnum6 ,i)
 
                                       
      call interpolation_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr_rf(iobs),phb,  &
                                    qru,qrl,qsu,qsl,qgru,qgrl,                       &
                                    rhou,rhol,tcu,tcl,                               &
                                    qrr,qsr,qgrr,rhor,tcr                            &
                                    )

  call ref_operator(qrr,qsr,qgrr,hx(i),ref_r,ref_s,ref_h,rhor,tcr)
  
  
  hxmean=hxmean+hx(i)/ens_num
 
  enddo
  
  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean  
  enddo    
!************************************************************************   
!************************************************************************  
! cal the innovation vecter    
  call innovation_cal(innovationmean,rf(iobs),hxmean)  

!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 
  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
  
            
! quality control if innovationmean is much great than threhold , assimilation will not be done
  if( abs(innovationmean) > ob_rf_threshold ) goto 102
  print*,'innovationmean,rf(iobs),hxmean',innovationmean,rf(iobs),hxmean  

  
! cal the ZZH 
!  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)
  
if(enable_rm_correlation_rf == 1 ) then  
! If the correlation between RF and state variation is not good or even lead to a worse result
! it should be remove.currently,it is found that rf will worsen the 1st level for all var


!  do ivar=1,analysis_var_num    
!  call rm_correlate(zzh,numstat,ivar,nx,ny,nz,2)
!  enddo
 

  call rm_correlate2(zzh,gridnum,gridsgn,nx,ny,nz,2) 
 
endif
    
! cal localization
  if(enable_local ==1) then
  
  
!do istatmp=1,numstat 
!       call localization_cal_radar(zzh,numstat,ir(iobs),jr(iobs),hr_rf(iobs),phb, &
!                                   analysis_var_num,istatmp)                                      
!enddo                                 


  call simple_local(zzh,gridnum,dis_rho)
  
  endif  

! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 
! cal kalman maxtria                                                                          
!  call kalman_cal(zzh,hph,obs_err_rf,numstat,kalman) 
  call kalman_cal(zzh,hph,obs_err_rf,gridnum,kalman)
                                            
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_rf) 
   
! cal analysis field 
!  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,   &
!                      xstat,xstat_mean,kalman,hx_pert                &
!                     )  
                     
  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 
                                   
    
! finished innovation for a rf obs
  print*,'finished innovation for a rf obs'


                       
102 continue                      
!  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)

  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman) 
  deallocate(hx) 
  deallocate(hx_pert)  

print*,'assimilate the obs at the',iobs,'th point'      
  
  endif       ! for reflectivity 
                          
  enddo  !  end loop for obs                         

                             
end subroutine assimilate_radar_cal         




subroutine assimilate_sounding_cal(analysis_var_num,obs_snd_num,numstat,                          &
                                   xstat,xstat_mean,                                              &
                                   usnd,vsnd,temperaturesnd,latsnd,lonsnd,hgtsnd,                 &
                                   xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                                   rdnw,rdn,mub,mu                                                &
                                   ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!            subroutine assimilate_snd_cal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                                 
include 'namelist.inc'
  
  integer               :: analysis_var_num,obs_snd_num,numstat
  integer               :: iobs,ivar
  integer               :: nn
  
  integer               :: istatnum1,istatnum2,istatnum3,istatnum4,      &
                           istatnum6,istatnum5,istatnum7,istatnum8,      &
                           istatnum9,istatnum10,istatnum11,istatnum12,   &
                           istatnum13,istatnum14,istatnum15,istatnum16,  &
                           istatnum17,istatnum18,istatnum19,istatnum20   
                           
  real                  :: uusw,uuse,uunw,uune,                   &
                           ulsw,ulse,ulnw,ulne,                   &
                           vusw,vuse,vunw,vune,                   &
                           vlsw,vlse,vlnw,vlne,                   &
                           tusw,tuse,tunw,tune,                   &
                           tlsw,tlse,tlnw,tlne                                              
  
  real                  :: ubsnd,vbsnd,tbsnd
  
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real                  :: xstat(numstat,ens_num),xstat_mean(numstat)
  
  real,parameter        :: gravity=9.81                         
  real                  :: usnd(obs_snd_num),vsnd(obs_snd_num),temperaturesnd(obs_snd_num),hgtsnd(obs_snd_num)
  real                  :: latsnd(obs_snd_num),lonsnd(obs_snd_num)
  integer               :: ism,jsm,ksm
  real                  :: isnd,jsnd,ksnd

                           
  real,allocatable      :: hx_pert(:),hx(:),innovation(:)  
  real                  :: hxmean  

  real,allocatable      :: zzh(:),kalman(:)
  real                  :: hph,innovationmean
  real                  :: kalman_alpha
  
  real                  :: hor_distance,vert_distance,dis_cor
  
  real                  :: temp1(ens_num),temp2(ens_num),temp3(ens_num),temp4(ens_num), &
                           temp5(ens_num),temp6(ens_num),temp7(ens_num),temp8(ens_num), &
                           temp9(ens_num),temp10(ens_num)
                           
  real                  :: tem1,tem2,tem3,tem4,tem5
  real                  :: time1,time2,use_time
  integer               :: istatmp,imp,jmp,kmp
  integer               :: itemp1,itemp2,itemp3

  real                  :: minlon,minlat,maxlon,maxlat
  real                  :: wwsn,wsew,wul

  
!*****************************************************************
  real,allocatable      :: xstat_sub(:,:),xstat_mean_sub(:)
  real,allocatable      :: dis_rho(:)
  integer,allocatable   :: gridsgn(:)      
  integer               :: gridnum 
  real,allocatable      :: temp_dis_rho(:)
  integer,allocatable   :: temp_gridsgn(:)
!*****************************************************************


print*,'find max lat & lon and min lat & lon to exclude the obs out of fcst area'
!********************************************
!  find max lat & lon and min lat & lon
!********************************************
do j=1,ny
  call cal_minvalue_1d(xlon(:,j),tem1,nx)
  if(j==1) then
    minlon=tem1
  else
    if(tem1 < minlon)  minlon= tem1  
  endif
enddo

do j=1,ny
  call cal_minvalue_1d(xlat(:,j),tem1,nx)
  if(j==1) then
    minlat=tem1
  else
    if(tem1 < minlat)  minlat= tem1  
  endif
enddo

do j=1,ny
  call cal_maxvalue_1d(xlon(:,j),tem1,nx)
  if(j==1) then
    maxlon=tem1
  else
    if(tem1 > maxlon)  maxlon= tem1  
  endif
enddo

do j=1,ny
  call cal_maxvalue_1d(xlat(:,j),tem1,nx)
  if(j==1) then
    maxlat=tem1
  else
    if(tem1 > maxlat)  maxlat= tem1  
  endif
enddo

print*,'minlon,minlat,maxlon,maxlat',minlon,minlat,maxlon,maxlat

!********************************************
!  begin the obs loop
!********************************************
print*,'begin the obs loop'

do iobs=1,obs_snd_num

!*********************prepare for HX cal*********************************
 call obs_latlonhgt_to_ijk(latsnd(iobs),lonsnd(iobs),hgtsnd(iobs),isnd,jsnd,ksnd,    &
                          phb,nx,ny,nz,xlat,xlon,                 &
			  minlon,minlat,maxlon,maxlat,            &
                          ism,jsm,ksm,wwsn,wsew,wul               &
						  )
!print*,'xlat,xlon',latsnd(iobs),lonsnd(iobs)						  
!print*,'ism,jsm,ksm,isnd,jsnd,ksnd',ism,jsm,ksm	,isnd,jsnd,ksnd					  
if( ism == -1 .or. jsm== -1 .or. ksm == -1 ) cycle
!************************************************************************

if(assim_u == 1 .and. usnd(iobs) .ne. miss_data) then    

!***************************************************************
!           prepare information for local
!***************************************************************  
  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))


do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 

hx_chi=0
hx_chi_pert=0

  print*,'check whether the allocate is normal'
  print*,temp_dis_rho(1),temp_dis_rho(1)
  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)  
  
call prepare_info4local_snd(analysis_var_num,obs_snd_num,numstat,                        &                               
                            isnd,jsnd,ksnd,hgtsnd(iobs),ism,jsm,ksm,                 	&
                            phb,temp_dis_rho,temp_gridsgn,gridnum                                  &
                           )

                                           

  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
    
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))

  print*,'allocate finished'
  
  
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo

  
!  print*,'check whether the allocate is normal'
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum)
  
  nn=1
  do i=1,numstat
  if(temp_dis_rho(i) > 0.05 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
     do j=1,ens_num
        xstat_sub(nn,j)=xstat(gridsgn(nn),j) 
      enddo
   nn=nn+1
  endif
  enddo 

  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)
  print*,'there are ',gridnum,'grid being influented by obs'

   
!************************************************************************
!***********************for U HX cal*************************************
!************************************************************************ 
hxmean=0

do i=1,ens_num 


  ivar=1   ! for u 

  call ijk_to_istat(ivar,ism  ,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm  ,ksm+1,istatnum2 ,nx,ny,nz)  
     uusw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! uusw

  call ijk_to_istat(ivar,ism+1,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+2,jsm  ,ksm+1,istatnum2 ,nx,ny,nz)  
     uuse=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! uuse

  call ijk_to_istat(ivar,ism  ,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum2 ,nx,ny,nz)  
     uunw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! uunw

  call ijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+2,jsm+1,ksm+1,istatnum2 ,nx,ny,nz)  
     uune=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! uune
  !-------------------------------------------------------------------------------- 
  call ijk_to_istat(ivar,ism  ,jsm  ,ksm ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm  ,ksm ,istatnum2 ,nx,ny,nz)  
     ulsw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! ulsw

  call ijk_to_istat(ivar,ism+1,jsm  ,ksm ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+2,jsm  ,ksm ,istatnum2 ,nx,ny,nz)  
     ulse=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! ulse

  call ijk_to_istat(ivar,ism  ,jsm+1,ksm ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm+1,ksm ,istatnum2 ,nx,ny,nz)  
     ulnw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! ulnw

  call ijk_to_istat(ivar,ism+1,jsm+1,ksm ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+2,jsm+1,ksm ,istatnum2 ,nx,ny,nz)  
     ulne=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        ! ulne
  !--------------------------------------------------------------------------------

  tem1=uusw*(1-wwsn)+uunw*wwsn
  tem2=uuse*(1-wwsn)+uune*wwsn
  tem3=tem1*(1-wsew)+tem2*wsew           !  for lower u

  tem1=ulsw*(1-wwsn)+ulnw*wwsn
  tem2=ulse*(1-wwsn)+ulne*wwsn
  tem4=tem1*(1-wsew)+tem2*wsew           !  for upper u

  hx(i)=tem3*(1-wul)+tem4*wul   
           
  hxmean=hxmean+hx(i)/ens_num

enddo

  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean  
  enddo

!************************************************************************
!************************************************************************

! cal the innovation vecter    
  call innovation_cal(innovationmean,usnd(iobs),hxmean)
  
!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 
  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
    
  
! quality control if innovationmean is great than threhold , assimilation will not be done
  if( abs(innovationmean) > 20 ) goto 103  
  print*,'innovationmean,usnd(iobs),hxmean',innovationmean,usnd(iobs),hxmean    

  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)

! cal localization
  if(enable_local ==1) then

  call simple_local(zzh,gridnum,dis_rho)
      
  endif  
  
! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
  
  call kalman_cal(zzh,hph,obs_err_u,gridnum,kalman)                                               
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_u) 

  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 

! finished innovation for a usnd obs
  print*,'finished innovation for a u obs'

103 continue 
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman)  
  deallocate(hx) 
  deallocate(hx_pert)   

endif

if(assim_v == 1 .and. vsnd(iobs) .ne. miss_data) then    

!***************************************************************
!           prepare information for local
!***************************************************************  
  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))


do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 

hx_chi=0
hx_chi_pert=0

!  print*,'check whether the allocate is normal'
!  print*,temp_dis_rho(1),temp_dis_rho(1)
!  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)  
  
call prepare_info4local_snd(analysis_var_num,obs_snd_num,numstat,                        &                               
                            isnd,jsnd,ksnd,hgtsnd(iobs),ism,jsm,ksm,                 	&
                            phb,temp_dis_rho,temp_gridsgn,gridnum                                  &
                           )

                                           

  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
    
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))

  print*,'allocate finished'
  
  
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo

  
!  print*,'check whether the allocate is normal'
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum)
  
  nn=1
  do i=1,numstat
  if(temp_dis_rho(i) > 0.05 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
     do j=1,ens_num
        xstat_sub(nn,j)=xstat(gridsgn(nn),j) 
      enddo
   nn=nn+1
  endif
  enddo 

  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)
  print*,'there are ',gridnum,'grid being influented by obs'

   
!************************************************************************
!***********************for V HX cal*************************************
!************************************************************************ 
hxmean=0

do i=1,ens_num 

ivar=2   ! for v
  
  call ijk_to_istat(ivar,ism  ,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism  ,jsm+1,ksm+1,istatnum2 ,nx,ny,nz)  
     vusw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vusw

  call ijk_to_istat(ivar,ism+1,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum2 ,nx,ny,nz)  
     vuse=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vuse

  call ijk_to_istat(ivar,ism  ,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism  ,jsm+2,ksm+1,istatnum2 ,nx,ny,nz)  
     vunw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vunw

  call ijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm+2,ksm+1,istatnum2 ,nx,ny,nz)  
     vune=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vune
  !-------------------------------------------------------------------------------- 
  call ijk_to_istat(ivar,ism  ,jsm  ,ksm  ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism  ,jsm+1,ksm  ,istatnum2 ,nx,ny,nz)  
     vlsw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vlsw

  call ijk_to_istat(ivar,ism+1,jsm  ,ksm  ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm+1,ksm  ,istatnum2 ,nx,ny,nz)  
     vlse=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vlse

  call ijk_to_istat(ivar,ism  ,jsm+1,ksm  ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism  ,jsm+2,ksm  ,istatnum2 ,nx,ny,nz)  
     vlnw=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vlnw

  call ijk_to_istat(ivar,ism+1,jsm+1,ksm  ,istatnum1 ,nx,ny,nz)  
  call ijk_to_istat(ivar,ism+1,jsm+2,ksm  ,istatnum2 ,nx,ny,nz)  
     vlne=(xstat(istatnum1,i)+xstat(istatnum2,i))*0.5                        !  vlne
  !--------------------------------------------------------------------------------

  tem1=vusw*(1-wwsn)+vunw*wwsn
  tem2=vuse*(1-wwsn)+vune*wwsn
  tem3=tem1*(1-wsew)+tem2*wsew           !  for lower v

  tem1=vlsw*(1-wwsn)+vlnw*wwsn
  tem2=vlse*(1-wwsn)+vlne*wwsn
  tem4=tem1*(1-wsew)+tem2*wsew           !  for upper v

  hx(i)=tem3*(1-wul)+tem4*wul   
           
  hxmean=hxmean+hx(i)/ens_num

enddo

  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean  
  enddo

!************************************************************************
!************************************************************************

! cal the innovation vecter    
  call innovation_cal(innovationmean,vsnd(iobs),hxmean)
  
!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 
  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
    
  
! quality control if innovationmean is great than threhold , assimilation will not be done
  if( abs(innovationmean) > 20 ) goto 104  
  print*,'innovationmean,vsnd(iobs),hxmean',innovationmean,vsnd(iobs),hxmean    


  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)

! cal localization
  if(enable_local ==1) then

  call simple_local(zzh,gridnum,dis_rho)
      
  endif  
  
! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
  
  call kalman_cal(zzh,hph,obs_err_v,gridnum,kalman)                                               
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_v) 

  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 

! finished innovation for a vsnd obs
  print*,'finished innovation for a v obs'

104 continue 
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman)  
  deallocate(hx) 
  deallocate(hx_pert)  



endif


if(assim_t == 1 .and. temperaturesnd(iobs) .ne. miss_data) then    

!***************************************************************
!           prepare information for local
!***************************************************************  
  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))


do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 

hx_chi=0
hx_chi_pert=0

!  print*,'check whether the allocate is normal'
!  print*,temp_dis_rho(1),temp_dis_rho(1)
!  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)  
  
call prepare_info4local_snd(analysis_var_num,obs_snd_num,numstat,                        &                               
                            isnd,jsnd,ksnd,hgtsnd(iobs),ism,jsm,ksm,                 	&
                            phb,temp_dis_rho,temp_gridsgn,gridnum                                  &
                           )

                                           

  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
    
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))

  print*,'allocate finished'
  
  
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo

  
!  print*,'check whether the allocate is normal'
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum)
  
  nn=1
  do i=1,numstat
  if(temp_dis_rho(i) > 0.05 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
     do j=1,ens_num
        xstat_sub(nn,j)=xstat(gridsgn(nn),j) 
      enddo
   nn=nn+1
  endif
  enddo 

  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)
  print*,'there are ',gridnum,'grid being influented by obs'

   
!************************************************************************
!***********************for T HX cal*************************************
!************************************************************************ 
 
do i=1,ens_num 
 
 ivar=3   ! for t 

  call ijk_to_istat(ivar,ism  ,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
     tusw=xstat(istatnum1,i)                                                 ! tusw

  call ijk_to_istat(ivar,ism+1,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
     tuse=xstat(istatnum1,i)                                                 ! tuse

  call ijk_to_istat(ivar,ism  ,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
     tunw=xstat(istatnum1,i)                                                 ! tunw

  call ijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
     tune=xstat(istatnum1,i)                                                 ! tune
  !-------------------------------------------------------------------------------- 

  call ijk_to_istat(ivar,ism  ,jsm  ,ksm,istatnum1 ,nx,ny,nz)  
     tusw=xstat(istatnum1,i)                                                 ! tlsw

  call ijk_to_istat(ivar,ism+1,jsm  ,ksm,istatnum1 ,nx,ny,nz)  
     tuse=xstat(istatnum1,i)                                                 ! tlse

  call ijk_to_istat(ivar,ism  ,jsm+1,ksm,istatnum1 ,nx,ny,nz)  
     tunw=xstat(istatnum1,i)                                                 ! tlnw

  call ijk_to_istat(ivar,ism+1,jsm+1,ksm,istatnum1 ,nx,ny,nz)  
     tune=xstat(istatnum1,i)                                                 ! tlne


  tem1=uusw*(1-wwsn)+uunw*wwsn
  tem2=uuse*(1-wwsn)+uune*wwsn
  tem3=tem1*(1-wsew)+tem2*wsew           !  for lower t

  tem1=ulsw*(1-wwsn)+ulnw*wwsn
  tem2=ulse*(1-wwsn)+ulne*wwsn
  tem4=tem1*(1-wsew)+tem2*wsew           !  for upper t

  hx(i)=tem3*(1-wul)+tem4*wul   
           
  hxmean=hxmean+hx(i)/ens_num

enddo

  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean  
  enddo

!************************************************************************
!************************************************************************

! cal the innovation vecter    
  call innovation_cal(innovationmean,temperaturesnd(iobs),hxmean)
  
!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 
  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
    
  
! quality control if innovationmean is great than threhold , assimilation will not be done
  if( abs(innovationmean) > 10 ) goto 105  
  print*,'innovationmean,temperaturesnd(iobs),hxmean',innovationmean,temperaturesnd(iobs),hxmean    


  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)

! cal localization
  if(enable_local ==1) then

  call simple_local(zzh,gridnum,dis_rho)
      
  endif  
  
! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 
  call kalman_cal(zzh,hph,obs_err_t,gridnum,kalman)                                               
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_t) 

  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 

! finished innovation for a temperaturesnd obs
  print*,'finished innovation for a temperaturesnd obs'

105 continue 
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman)  
  deallocate(hx) 
  deallocate(hx_pert)   


endif




print*,'assimilate the obs at the',iobs,'th point'  

enddo  ! do iobs=1,obs_snd_num
                               

end subroutine assimilate_sounding_cal 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! land info assimilate with EnSRF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assimilate_land_cal(analysis_var_num,obs_lad_num,numstat,                                     &
                                   xstat,xstat_mean,                                                     &
                                   soil_m,soil_t,soil_st,latsnd,lonsnd,hgtsnd,                           &
                                   xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,hgt,znw,znu,p_top,              &
                                   rdnw,rdn,mub,mu                                                       &
                                   ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!            subroutine assimilate_land_cal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.05.26   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                                 
include 'namelist.inc'
  
  integer               :: analysis_var_num,obs_lad_num,numstat
  integer               :: iobs,ivar
  integer               :: nn,inum
  
  integer               :: istatnum1,istatnum2,istatnum3,istatnum4,      &
                           istatnum6,istatnum5,istatnum7,istatnum8,      &
                           istatnum9,istatnum10,istatnum11,istatnum12,   &
                           istatnum13,istatnum14,istatnum15,istatnum16,  &
                           istatnum17,istatnum18,istatnum19,istatnum20   
                           
  real                 :: uusw,uuse,uunw,uune,                   &
                           ulsw,ulse,ulnw,ulne,                   &
                           vusw,vuse,vunw,vune,                   &
                           vlsw,vlse,vlnw,vlne,                   &
                           tusw,tuse,tunw,tune,                   &
                           tlsw,tlse,tlnw,tlne 
  
  real                  :: stusw,stuse,stunw,stune,                   &
                           stlsw,stlse,stlnw,stlne,                   &
                           smusw,smuse,smunw,smune,                   &
                           smlsw,smlse,smlnw,smlne,                   &                        
                           strsw,strse,strnw,strne                   !! str* -> soil_st 
                           
  real                  :: ubsnd,vbsnd,tbsnd
  
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
!!!for soil depth!!!
  real                  :: dzs(4)
  real                  :: xstat(numstat,ens_num),xstat_mean(numstat)
  
  real,parameter        :: gravity=9.81                         
  real                  :: usnd(obs_lad_num),vsnd(obs_lad_num),soil_t(obs_lad_num),soil_m(obs_lad_num),hgtsnd(obs_lad_num)
  real                  :: soil_st(obs_lad_num),temp_soil_t(obs_lad_num),temp_soil_m(obs_lad_num),temp_hgtsnd(obs_lad_num),&
                           temp_soil_st(obs_lad_num)
  real                  :: latsnd(obs_lad_num),lonsnd(obs_lad_num),temp_latsnd(obs_lad_num),temp_lonsnd(obs_lad_num)
  integer               :: ism,jsm,ksm
  real                  :: isnd,jsnd,ksnd

                           
  real,allocatable      :: hx_pert(:),hx(:),innovation(:)  
  real                  :: hxmean  

  real,allocatable      :: zzh(:),kalman(:)
  real                  :: hph,innovationmean
  real                  :: kalman_alpha
  
  real                  :: hor_distance,vert_distance,dis_cor
  
  real                  :: temp1(ens_num),temp2(ens_num),temp3(ens_num),temp4(ens_num), &
                           temp5(ens_num),temp6(ens_num),temp7(ens_num),temp8(ens_num), &
                           temp9(ens_num),temp10(ens_num)
                           
  real                  :: tem1,tem2,tem3,tem4,tem5
  real                  :: time1,time2,use_time
  integer               :: istatmp,imp,jmp,kmp
  integer               :: itemp1,itemp2,itemp3
!  integer               :: th_num

  real                  :: minlon,minlat,maxlon,maxlat
  real                  :: wwsn,wsew,wul

  
!*****************************************************************
  real,allocatable      :: xstat_sub(:,:),xstat_mean_sub(:)
  real,allocatable      :: dis_rho(:)
  integer,allocatable   :: gridsgn(:)      
  integer               :: gridnum 
  real,allocatable      :: temp_dis_rho(:)
  integer,allocatable   :: temp_gridsgn(:)
!*****************************************************************

!open(1001,file='/home/guoyk/test.txt')

print*,'find max lat & lon and min lat & lon to exclude the obs out of fcst area'
!print*,xlon(1,1),xlat(1,1)
!**********************************************************
!  find max lat & lon and min lat & lon of landinfo
!**********************************************************
do j=1,ny
  call cal_minvalue_1d(xlon(:,j),tem1,nx)
  if(j==1) then
    minlon=tem1
  else
    if(tem1 < minlon)  minlon= tem1  
  endif
enddo

do j=1,ny
  call cal_minvalue_1d(xlat(:,j),tem1,nx)
  if(j==1) then
    minlat=tem1
  else
    if(tem1 < minlat)  minlat= tem1  
  endif
enddo

do j=1,ny
  call cal_maxvalue_1d(xlon(:,j),tem1,nx)

  if(j==1) then
    maxlon=tem1
  else
    if(tem1 > maxlon)  maxlon= tem1  
  endif
enddo

do j=1,ny
  call cal_maxvalue_1d(xlat(:,j),tem1,nx)
  if(j==1) then
    maxlat=tem1
  else
    if(tem1 > maxlat)  maxlat= tem1  
  endif
enddo

!write(1001,*)'minlon,minlat,maxlon,maxlat',minlon,minlat,maxlon,maxlat

print*,'minlon,minlat,maxlon,maxlat',minlon,minlat,maxlon,maxlat

!********************************************
!  begin the obs loop of u
!********************************************


!write(1001,*)'iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul,soil_t(iobs) which came from obs_latlonhgtlad_to_ijk'

if(enable_thinning == 1) then

 print*,'BEGIN THINNING OBS,the thinning Num is',th_num !,obs_lad_num,obs_lad_num/4

 open(100,file='./Thinning_data.txt')
  inum=1  !new var
    do i=1,obs_lad_num/th_num
!      do j=i+4,
     temp_latsnd(inum)=latsnd((i-1)*th_num+1)
     temp_lonsnd(inum)=lonsnd((i-1)*th_num+1)
     temp_hgtsnd(inum)=hgtsnd((i-1)*th_num+1)
     temp_soil_m(inum)=soil_m((i-1)*th_num+1)
     temp_soil_t(inum)=soil_t((i-1)*th_num+1)
     temp_soil_st(inum)=soil_st((i-1)*th_num+1)
     inum=inum+1
! print*,i,inum,i+4,soil_st(i+4)
! print*,temp_soil_m(inum-1),temp_soil_t(inum-1),temp_soil_st(inum-1),temp_latsnd(inum-1),temp_lonsnd(inum-1),temp_hgtsnd(inum-1)
 write(100,'(6f11.6)')temp_soil_m(inum-1),temp_soil_t(inum-1),temp_soil_st(inum-1),temp_latsnd(inum-1),temp_lonsnd(inum-1),temp_hgtsnd(inum-1)
    enddo

 print*,'After thinning, there is',inum,'obs left'

else

temp_latsnd = latsnd
temp_lonsnd = lonsnd
temp_hgtsnd = hgtsnd
temp_soil_m  = soil_m
temp_soil_t  = soil_t
temp_soil_st = soil_st
inum = obs_lad_num+1  ! +1 is for next iobs=1,inum-1;because count one more number when do thinning.

endif !if(enable_thinning == 1)

print*,'*********************************'
print*,'begin the obs landinfomation loop'
print*,'*********************************'

do iobs=1,inum-1


!print*, iobs,dzs(1)
!print*,'latsnd(1),lonsnd(1),xlat(1),xlon(1)',latsnd(1),lonsnd(1),xlat(1,1),xlon(1,1)

!*********************prepare for HX cal*********************************
 call obs_latlonhgtlad_to_ijk(temp_latsnd(iobs),temp_lonsnd(iobs),temp_hgtsnd(iobs),isnd,jsnd,ksnd,    &
                          phb,dzs,nx,ny,nz,xlat,xlon,                                   &
			                     minlon,minlat,maxlon,maxlat,                                  &
                          ism,jsm,ksm,wwsn,wsew,wul                                     &
						  )

			  
if( ism == -1 .or. jsm== -1 .or. ksm == -1 ) cycle

!************************************************************************

!print*,'iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul',iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul

!stop  2013/5/12 ok
!write(1001,*)iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul,soil_t(iobs)
!************************************************************************
!print*,'assim_tslb,soil_t(iobs),iobs ',assim_tslb,soil_t(iobs),iobs

if(assim_smois == 1 .and. temp_soil_m(iobs) .ne. miss_data) then    

!print*,'!***************************************************************'
!print*,'!           prepare information for local           !!!soil_m   '
!print*,'!***************************************************************'

ivar=11
  
  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))


do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 

hx_chi=0
hx_chi_pert=0

!  print*,'check whether the allocate is normal'
!  print*,temp_dis_rho(1),temp_dis_rho(1)
!  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)  
  
 call prepare_info4local_lad(analysis_var_num,inum-1,numstat,ivar,                        &                               
                            isnd,jsnd,ksnd,temp_hgtsnd(iobs),ism,jsm,ksm,                 &
                            phb,dzs,temp_dis_rho,temp_gridsgn,gridnum                     &
                           )

!  print*,gridnum                                         
!stop
  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
    
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))

  print*,'allocate finished'
  
  
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo

  
!  print*,'check whether the allocate is normal'
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum)
  
  nn=1
  do i=1,numstat
  if(temp_dis_rho(i) > 0.05 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
     do j=1,ens_num
        xstat_sub(nn,j)=xstat(gridsgn(nn),j) 
      enddo
!print*,'xstat_sub(nn,j),gridsgn(nn)',xstat_sub(nn,j),nn,gridsgn(nn),xstat(3925536,j),xstat_mean(3925536)
   nn=nn+1
  endif
  enddo 

  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)
  print*,'there are ',gridnum,'grid being influented by soilmobs'

   
!************************************************************************
!***********************for SOIL_M HX cal*************************************
!************************************************************************ 
!print*,'iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul',iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul
 
hxmean=0

do i=1,ens_num 
 
 ivar=11   ! for soil_m 

  call ladijk_to_istat(ivar,ism  ,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
     smusw=xstat(istatnum1,i)                                                 ! tusw
!print*,'istatnum1,tusw',istatnum1,tusw
  call ladijk_to_istat(ivar,ism+1,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
     smuse=xstat(istatnum1,i)                                                 ! tuse
!print*,'istatnum1,tuse',istatnum1,tuse
  call ladijk_to_istat(ivar,ism  ,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
     smunw=xstat(istatnum1,i)                                                 ! tunw
!print*,'istatnum1,tunw',istatnum1,tunw
  call ladijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
     smune=xstat(istatnum1,i)                                                 ! tune
  !-------------------------------------------------------------------------------- 
!print*,'istatnum1,tune',istatnum1,tune
  call ladijk_to_istat(ivar,ism  ,jsm  ,ksm,istatnum1 ,nx,ny,nz)  
     smlsw=xstat(istatnum1,i)                                                 ! tlsw
!print*,'istatnum1,tlsw,xstat(istatnum1,i)',istatnum1,tlsw,xstat(istatnum1,i)
  call ladijk_to_istat(ivar,ism+1,jsm  ,ksm,istatnum1 ,nx,ny,nz)  
     smlse=xstat(istatnum1,i)                                                 ! tlse
!print*,'istatnum1,tlse,xstat(istatnum1,i)',istatnum1,tlse,xstat(istatnum1,i)
  call ladijk_to_istat(ivar,ism  ,jsm+1,ksm,istatnum1 ,nx,ny,nz)  
     smlnw=xstat(istatnum1,i)                                                 ! tlnw
!print*,'istatnum1,tlnw,xstat(istatnum1,i)',istatnum1,tlnw,xstat(istatnum1,i)
  call ladijk_to_istat(ivar,ism+1,jsm+1,ksm,istatnum1 ,nx,ny,nz)  
     smlne=xstat(istatnum1,i)                                                 ! tlne
!print*,'istatnum1,tlne,xstat(istatnum1,i)',istatnum1,tlne,xstat(istatnum1,i)

!  tem1=(tusw*wwsn+tunw*(1-wwsn))*1000
!  tem2=(tuse*wwsn+tune*(1-wwsn))*1000
!  tem3=tem1*wsew+tem2*(1-wsew)           !  for upper soilm

!  tem1=(tlsw*wwsn+tlnw*(1-wwsn))*1000
!  tem2=(tlse*wwsn+tlne*(1-wwsn))*1000
!  tem4=tem1*wsew+tem2*(1-wsew)           !  for lower soilm

!  hx(i)=(tem3*wul+tem4*(1-wul))/1000   
           
 
! hxmean=hxmean+hx(i)/ens_num


  tem1=smwsw*(1-wwsn)+smunw*wwsn
  tem2=smwse*(1-wwsn)+smune*wwsn
  tem3=tem1*(1-wsew)+tem2*wsew           !  for lower soilm

  tem1=smlsw*(1-wwsn)+smlnw*wwsn
  tem2=smlse*(1-wwsn)+smlne*wwsn
  tem4=tem1*(1-wsew)+tem2*wsew           !  for upper soilm

  hx(i)=tem3*(1-wul)+tem4*wul   
           
  hxmean=hxmean+hx(i)/ens_num

!print*,'i,hx(i),hxmean',i,hx(i),hxmean
!stop


enddo

  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean 
!print*, hx_pert(i)
  enddo
!stop
!************************************************************************
!print*,'tem1,tem2,tem3,tem4,hxmean,tusw,tuse,tunw,tune',tem1,tem2,tem3,tem4,hxmean,tusw,tuse,tunw,tune!xstat(3920001,1),&
!xstat(3920001,2),xstat(3934001,1),xstat(3934001,2)
!************************************************************************
!************************************************************************

! cal the innovation vecter   no change 

  call innovation_cal(innovationmean,temp_soil_m(iobs),hxmean)
  
!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) no change

  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
    
!  print*,'innovationmean,soil_m(iobs),hxmean',innovationmean,soil_m(iobs),hxmean

! quality control if innovationmean is great than threhold , assimilation will not be done
!  stop

! if( abs(innovationmean) > 0.8 ) goto 200   
  if( abs(innovationmean) > 3*obs_err_smois ) goto 200 !3*obs_error
write(*,'(a30,3f11.6)')'innovmean,soil_m(iobs),hxmean',innovationmean,temp_soil_m(iobs),hxmean    


  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)

! cal localization
  if(enable_local ==1) then

  call simple_local(zzh,gridnum,dis_rho)
      
  endif  
  
! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 
  call kalman_cal(zzh,hph,obs_err_smois,gridnum,kalman)                                               
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_smois) 

  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 

! finished innovation for a temperaturesnd obs
!  print*,'finished innovation for a soil_m obs'

200 continue 
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman)  
  deallocate(hx) 
  deallocate(hx_pert)   


endif

!**************************************************
!**************************************************
!**************************************************

if(assim_tslb == 1 .and. temp_soil_t(iobs) .ne. miss_data) then    

! print*,'xstat(3976001,2),enstslb(1,1,1,1),soil_t(1)',xstat(3976001,2),soil_t(1)
!print*,'!***************************************************************'
!print*,'!           prepare information for local           !!!soil_t   '
!print*,'!***************************************************************' 
 
ivar=12

  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))


do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 

hx_chi=0
hx_chi_pert=0

!write(1001,*)'first temp_dis_rho == 0 & temp_gridsgn == 0 '
!  print*,'check whether the soil allocate is normal'
!  print*,temp_dis_rho(1),temp_dis_rho(1)
!  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)  
!print*,'xlat,xlon',latsnd(iobs),lonsnd(iobs),hgtsnd(iobs)						  
!print*,'ism,jsm,ksm,isnd,jsnd,ksnd',ism,jsm,ksm	,isnd,jsnd,ksnd	 

 call prepare_info4local_lad(analysis_var_num,inum-1,numstat,ivar,                        &                               
                            isnd,jsnd,ksnd,temp_hgtsnd(iobs),ism,jsm,ksm,                 &
                            phb,dzs,temp_dis_rho,temp_gridsgn,gridnum                     &
                           )
!open(100,file='/home/guoyk/dis_rho.txt')
!write(100,*)temp_dis_rho
!stop
!write(1001,*)'after call prepare_info4local_lad, write gridnum,temp_dis_rho,temp_gridsgn',gridnum,temp_dis_rho,temp_gridsgn

  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
    
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))

!  print*,'allocate finished'
  
  
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo

  
!  print*,'check whether the allocate is normal,gridnum',gridnum
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum)

!write(1001,*)'first dis_rho(nn) & gridsgn(nn) are 0 ,then dis = temp_dis, grid = temp_grid '  
  nn=1
  do i=1,numstat
  if(temp_dis_rho(i) > 0.05 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
   
     do j=1,ens_num
!!!!!!!!!! -30 & 30 degree of TLSB is impossible !!!!!!!!!!!!!
    !   if(gridsgn(nn) > 3976001 .and. xstat(gridsgn(nn),j) > 243.15 .and. 303.15 > xstat(gridsgn(nn),j)) then 
        xstat_sub(nn,j)=xstat(gridsgn(nn),j)
    !   endif 
      
     enddo
!write(1001,*)'j,dis_rho(nn),gridsgn(nn),nn,xstat_sub(nn,10),xstat(gridsgn(nn),10)',j-1,dis_rho(nn),gridsgn(nn),nn,xstat_sub(nn,10),xstat(gridsgn(nn),10)
! print*,'gridsgn(nn),xstat(gridsgn(nn),j),xstat(3920001,j)',gridsgn(nn),xstat(gridsgn(nn),j),xstat(3920001,j)
! stop
  nn=nn+1
  endif
  enddo
 



  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)
!  print*,'there are ',gridnum,'grid being influented by soiltobs'

   
!************************************************************************
!***********************for SOIL_T HX cal*************************************
!************************************************************************ 
!write(1001,*)'HX:i,hxmean,hx(i),xstat(istatnum1,i),tlsw,tlse,wul,iobs,wwsn,wsew ' 

hxmean=0

do i=1,ens_num 
 
 ivar=12   ! for soil_t 

  call ladijk_to_istat(ivar,ism  ,jsm  ,ksm+1,istatnum1 ,nx,ny,nz,4)  
     stusw=xstat(istatnum1,i)                                                 ! tusw
  !print*,istatnum1
  call ladijk_to_istat(ivar,ism+1,jsm  ,ksm+1,istatnum1 ,nx,ny,nz,4)  
     stuse=xstat(istatnum1,i)                                                 ! tuse
  !print*,istatnum1
  call ladijk_to_istat(ivar,ism  ,jsm+1,ksm+1,istatnum1 ,nx,ny,nz,4)  
     stunw=xstat(istatnum1,i)                                                 ! tunw
  !print*,istatnum1
  call ladijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum1 ,nx,ny,nz,4)  
     stune=xstat(istatnum1,i)                                                 ! tune
  !-------------------------------------------------------------------------------- 
  !   print*,istatnum1
  call ladijk_to_istat(ivar,ism  ,jsm  ,ksm,istatnum1 ,nx,ny,nz,4)  
     stlsw=xstat(istatnum1,i)                                                 ! tlsw
!   print*,'istatnum1,tlsw,xstat(istatnum1,i)',istatnum1,tlsw,xstat(istatnum1,i)
  call ladijk_to_istat(ivar,ism+1,jsm  ,ksm,istatnum1 ,nx,ny,nz,4)  
     stlse=xstat(istatnum1,i)                                                 ! tlse
!   print*,'istatnum1,tlse,xstat(istatnum1,i)',istatnum1,tlse,xstat(istatnum1,i)

  call ladijk_to_istat(ivar,ism  ,jsm+1,ksm,istatnum1 ,nx,ny,nz,4)  
     stlnw=xstat(istatnum1,i)                                                 ! tlnw
!    print*,'istatnum1,tlnw,xstat(istatnum1,i)',istatnum1,tlnw,xstat(istatnum1,i)
  call ladijk_to_istat(ivar,ism+1,jsm+1,ksm,istatnum1 ,nx,ny,nz,4)  
     stlne=xstat(istatnum1,i)                                                 ! tlne
!     print*,'istatnum1,tlne,xstat(istatnum1,i)',istatnum1,tlne,xstat(istatnum1,i)

  tem1=stusw*(1-wwsn)+stunw*wwsn
  tem2=stuse*(1-wwsn)+stune*wwsn
  tem3=tem1*(1-wsew)+tem2*wsew           !  for lower soilt

  tem1=stlsw*(1-wwsn)+stlnw*wwsn
  tem2=stlse*(1-wwsn)+stlne*wwsn
  tem4=tem1*(1-wsew)+tem2*wsew           !  for upper soilt

  hx(i)=tem3*(1-wul)+tem4*wul
 
!  tem1=tusw*wwsn+tunw*(1-wwsn)
!  tem2=tuse*wwsn+tune*(1-wwsn)
!  tem3=tem1*wsew+tem2*(1-wsew)           !  for upper soilm

!  tem1=tlsw*wwsn+tlnw*(1-wwsn)
!  tem2=tlse*wwsn+tlne*(1-wwsn)
!  tem4=tem1*wsew+tem2*(1-wsew)           !  for lower soilm

!  hx(i)=tem3*wul+tem4*(1-wul)
           
  hxmean=hxmean+hx(i)/ens_num

!print*,'tem1,tem2,tem3,tem4',tem1,tem2,tem3,tem4
!print*, 'hx(i),xstat(istatnum1,i),tlsw,tlse,wul,iobs,wwsn,wsew',hx(i),xstat(istatnum1,i),tlsw,tlse,wul,iobs,wwsn,wsew 
!stop
!write(1001,*)i,hxmean,hx(i),xstat(istatnum1,i),tlsw,tlse,wul,iobs,wwsn,wsew 

enddo


 

  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean  
  enddo

!************************************************************************
!************************************************************************
!print*,'soil_t(iobs),hxmean',soil_t(iobs),hxmean
!stop
! cal the innovation vecter    
 
 call innovation_cal(innovationmean,temp_soil_t(iobs),hxmean)

! write(1001,*)'innovationmean'
! write(1001,*)innovationmean  

!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 
  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
    
  
! quality control if innovationmean is great than threhold , assimilation will not be done
!  if( abs(innovationmean) > 10 ) goto 205 
 if( abs(innovationmean) > 3*obs_err_tslb ) goto 205 !3*obs_error 
!  print*,'innovationmean,soil_t(iobs),hxmean',innovationmean,temp_soil_t(iobs),hxmean 
 write(*,'(a30,3f11.6)')'innovmean,soil_t(iobs),hxmean',innovationmean,temp_soil_t(iobs),hxmean 
!  write(1001,*)'abs(innovationmean) > 10 then write,soil_t(iobs),hxmean'
!  write(1001,*)innovationmean,soil_t(iobs),hxmean

  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)
 
!  write(1001,*)'zzh',zzh

! cal localization
  if(enable_local ==1) then

  call simple_local(zzh,gridnum,dis_rho)
      
  endif  
  
! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)
  
! write(1001,*)'hph',hph

  call kalman_cal(zzh,hph,obs_err_tslb,gridnum,kalman)                                               
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_tslb) 

! write(1001,*)'kalman,kalman_alpha',kalman,kalman_alpha

  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 
! write(1001,*)'final xa == xstat_sub xstat_mean_sub',xstat_sub,xstat_mean_sub

! finished innovation for a temperaturesnd obs
! print*,'finished innovation for a soil_t obs'

205 continue 
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman)  
  deallocate(hx) 
  deallocate(hx_pert)   


endif

!********************************************************************************
!*****************************soil_st**********************************
!********************************************************************************

if(assim_tsk == 1 .and. temp_soil_st(iobs) .ne. miss_data) then    

!print*,'!***************************************************************'
print*,'!           prepare information for local           !!!soil_st   '
!print*,'!***************************************************************'
!stop

 ivar=13
 
  allocate(temp_dis_rho(numstat))
  allocate(temp_gridsgn(numstat))


do i=1,numstat
  temp_dis_rho(i)=0
  temp_gridsgn(i)=0
enddo 

hx_chi=0
hx_chi_pert=0

!  print*,'check whether the allocate is normal'
!  print*,temp_dis_rho(1),temp_dis_rho(1)
!  print*,temp_gridsgn(numstat),temp_gridsgn(numstat)  
  
 call prepare_info4local_lad(analysis_var_num,inum-1,numstat,ivar,                      &                               
                            isnd,jsnd,ksnd,temp_hgtsnd(iobs),ism,jsm,ksm,               &
                            phb,dzs,temp_dis_rho,temp_gridsgn,gridnum                   &
                           )

!  print*,gridnum                                         
!  stop

if(gridnum == 0) then
print*,'Maybe your local distance is too small, or not propoer!'
endif

  allocate(dis_rho(gridnum))
  allocate(gridsgn(gridnum)) 
  allocate(zzh(gridnum))
  allocate(kalman(gridnum)) 
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
    
  allocate(xstat_sub(gridnum,ens_num))
  allocate(xstat_mean_sub(gridnum))

!  print*,'allocate finished'
  
  
  do i=1,gridnum 
   xstat_mean_sub(i)=0
   dis_rho(i)=0
   gridsgn(i)=0
   zzh(i)=0
   kalman(i)=0
   do j=1,ens_num
    xstat_sub(i,j)=0
    hx(j)=0
    hx_pert(j)=0
   enddo
 enddo

  
!  print*,'check whether the allocate is normal'
!  print*,xstat_sub(1,1),xstat_sub(1,ens_num),gridsgn(1),dis_rho(1)
!  print*,xstat_sub(gridnum,1),xstat_sub(gridnum,ens_num),gridsgn(gridnum),dis_rho(gridnum)
!  stop

  nn=1
  do i=1,numstat
  if(temp_dis_rho(i) > 0.05 ) then
   gridsgn(nn)=temp_gridsgn(i) 
   dis_rho(nn)=temp_dis_rho(i)
     do j=1,ens_num
        xstat_sub(nn,j)=xstat(gridsgn(nn),j) 
!print*,'xstat_sub(nn,j),gridsgn(nn)',j,xstat_sub(nn,j),nn,gridsgn(nn),xstat(gridsgn(nn),j),xstat_mean(gridsgn(nn))
      enddo
!print*,'xstat_sub(nn,j),gridsgn(nn)',j,xstat_sub(nn,j),nn,gridsgn(nn),xstat(gridsgn(nn),j),xstat_mean(gridsgn(nn))
   nn=nn+1
  endif
  enddo 

!stop
  deallocate(temp_dis_rho)
  deallocate(temp_gridsgn)

  print*,'there are ',gridnum,'grid being influented by soilstobs'

   
!************************************************************************
!***********************for SOIL_ST HX cal*************************************
!************************************************************************ 
!print*,'iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul',iobs,isnd,jsnd,ksnd,ism,jsm,ksm,wwsn,wsew,wul
 
hxmean=0

do i=1,ens_num 
 
 ivar=13   ! for soil_st 

  call ladijk_to_istat(ivar,ism  ,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
     strsw=xstat(istatnum1,i)                                                 ! strsw
!print*,'istatnum1,tusw',istatnum1,tusw
  call ladijk_to_istat(ivar,ism+1,jsm  ,ksm+1,istatnum1 ,nx,ny,nz)  
     strse=xstat(istatnum1,i)                                                 ! strse
!print*,'istatnum1,tuse',istatnum1,tuse
  call ladijk_to_istat(ivar,ism  ,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
     strnw=xstat(istatnum1,i)                                                 ! strnw
!print*,'istatnum1,tunw',istatnum1,tunw
  call ladijk_to_istat(ivar,ism+1,jsm+1,ksm+1,istatnum1 ,nx,ny,nz)  
     strne=xstat(istatnum1,i)                                                 ! strne
  !-------------------------------------------------------------------------------- 


  tem1=strsw*(1-wwsn)+strnw*wwsn
  tem2=strse*(1-wwsn)+strne*wwsn
  hx(i)=tem1*(1-wsew)+tem2*wsew           !  for surface temp

 !! hx(i)=tem3*(1-wul)+tem4*wul            !?? 
           
  hxmean=hxmean+hx(i)/ens_num

!print*,'i,hx(i),hxmean',i,hx(i),hxmean
!stop


enddo

  do i=1,ens_num  
    hx_pert(i)=hx(i)-hxmean 
!print*, hx_pert(i)
  enddo
!stop
!************************************************************************
!print*,'tem1,tem2,tem3,tem4,hxmean,tusw,tuse,tunw,tune',tem1,tem2,tem3,tem4,hxmean,tusw,tuse,tunw,tune!xstat(3920001,1),&
!xstat(3920001,2),xstat(3934001,1),xstat(3934001,2)
!************************************************************************
!************************************************************************

! cal the innovation vecter   no change 

  call innovation_cal(innovationmean,temp_soil_st(iobs),hxmean)
  
!  call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) no change

  call xstat_mean_seperate_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num) 
    
!  print*,'innovationmean,soil_st(iobs),hxmean',innovationmean,soil_st(iobs),hxmean
!  quality control if innovationmean is great than threhold , assimilation will not be done
!  stop
  
  if( abs(innovationmean) > 5*obs_err_tsk ) goto 300  !5*obs_error
  write(*,'(a30,3f11.6)')'innovmean,soil_st(iobs),hxmean',innovationmean,temp_soil_st(iobs),hxmean    
 
 write(*,*),'*******************************************'
 write(*,'(a30,I4,a8)'),'Assimilate The Soil Obs At The',iobs,'Th Point'
 write(*,*),'*******************************************'

  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)

! cal localization

  if(enable_local ==1) then

  call simple_local(zzh,gridnum,dis_rho)
      
  endif  
  
! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 
  call kalman_cal(zzh,hph,obs_err_smois,gridnum,kalman)                                               
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_smois) 

  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 

! finished innovation for a temperaturesnd obs
!  print*,'finished innovation for a soil_st obs'

300 continue 
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman)  
  deallocate(hx) 
  deallocate(hx_pert)   


endif


!print*,'assimilate the soil obs at the',iobs,'th point'
! write(*,*),'*******************************************'
! write(*,'(a30,I4,a8)'),'Assimilate The Soil Obs At The',iobs,'Th Point'
! write(*,*),'*******************************************'  

enddo  ! do iobs=1,inum
                               

end subroutine assimilate_land_cal    
