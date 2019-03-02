subroutine assimilate_sim_cal(analysis_var_num,obs_sim_num,numstat,                          &
                              xstat,xstat_mean,                                              &
                              xyzu,xyzv,xyzw,xyzph,xyzt,xyzqv,xyzqr,xyzqs,xyzqi,xyzqgr,      &
                              obsu,obsv,obsw,obsph,obst,obsqv,obsqr,obsqs,obsqi,obsqgr,      &
                              xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                              rdnw,rdn,mub,mu                                                &
                             ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine assimilate_sim_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!#############################################################
!#         sim sounding obs assimilation subroutine
!#              all analysis vars have a obs 
!#      this is a very simple test code with lots of bugs
!#############################################################                       
include 'namelist.inc'
  
  integer               :: analysis_var_num,obs_sim_num,numstat
  integer               :: iobs,ivar
  integer               :: istatnum1,istatnum2  
  
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real                  :: xstat(numstat,ens_num),xstat_mean(numstat)
                           
  real                  :: obsu(obs_sim_num),obsv(obs_sim_num),obsw(obs_sim_num),obsph(obs_sim_num),   &
                           obst(obs_sim_num),obsqv(obs_sim_num),obsqr(obs_sim_num),obsqi(obs_sim_num), &
                           obsqs(obs_sim_num),obsqgr(obs_sim_num)

  integer               :: xyzu(3,obs_sim_num),xyzv(3,obs_sim_num),xyzw(3,obs_sim_num),xyzph(3,obs_sim_num),   &
                           xyzt(3,obs_sim_num),xyzqv(3,obs_sim_num),xyzqr(3,obs_sim_num),xyzqi(3,obs_sim_num), &
                           xyzqs(3,obs_sim_num),xyzqgr(3,obs_sim_num) 
                           
  real,allocatable      :: hx_pert(:),hx(:),innovation(:)  
  real                  :: hxmean  

  real,allocatable      :: zzh(:),kalman(:)
  real                  :: hph,innovationmean
  real                  :: kalman_alpha
  
  real                  :: hor_distance,vert_distance,dis_cor
  
  real                  :: temp1(ens_num),temp2(ens_num),temp3(ens_num)
  

  allocate(zzh(numstat))
  allocate(kalman(numstat))
  allocate(hx(ens_num)) 
  allocate(hx_pert(ens_num))
  allocate(innovation(ens_num))  

print*,'start the assimilation of sim obs'
print*,miss_data


 do iobs=1, obs_sim_num    !  loop for obs 
 do ivar=1,analysis_var_num  ! loop for var

 ! assimilate u
 if(ivar == 1 .and. assim_u == 1 .and. obsu(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)

 ! cal the position 
  call ijk_to_istat(ivar,xyzu(1,iobs)  ,xyzu(2,iobs)  ,xyzu(3,iobs)  ,istatnum1,nx,ny,nz)  
  call ijk_to_istat(ivar,xyzu(1,iobs)+1,xyzu(2,iobs)  ,xyzu(3,iobs)  ,istatnum2,nx,ny,nz)           
 ! cal the hx and hxmean for u   
  call obs_sim_grid_u(xstat(istatnum1,:),xstat(istatnum2,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
   temp2(i)=xstat(istatnum2,i)+xstat_mean(istatnum2)
  enddo 
   call obs_sim_grid_u(temp1,temp2,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 

 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsu(iobs),hxmean) 
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)   
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzu(1,iobs),xyzu(2,iobs),xyzu(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)
 ! cal kalman maxtria          
  call kalman_cal(zzh,hph,obs_err_u,numstat,kalman)
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_u) 
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,     &
                      xstat,xstat_mean,kalman,hx_pert                  &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                     
                 
 ! finished innovation for a u obs 
  print*,'finished innovation for a u obs'
  endif


! assimilate v 
 if(ivar == 2 .and. assim_v == 1 .and. obsv(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position   
  call ijk_to_istat(ivar,xyzv(1,iobs)  ,xyzv(2,iobs)  ,xyzv(3,iobs)  ,istatnum1,nx,ny,nz)  
  call ijk_to_istat(ivar,xyzv(1,iobs)  ,xyzv(2,iobs)+1,xyzv(3,iobs)  ,istatnum2,nx,ny,nz)             
 ! cal the hx and hxmean and hx_pert for v   
  call obs_sim_grid_v(xstat(istatnum1,:),xstat(istatnum2,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
   temp2(i)=xstat(istatnum2,i)+xstat_mean(istatnum2)
  enddo 
   call obs_sim_grid_v(temp1,temp2,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsv(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzv(1,iobs),xyzv(2,iobs),xyzv(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 ! cal kalman maxtria          
  call kalman_cal(zzh,hph,obs_err_v,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_v)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,     &
                      xstat,xstat_mean,kalman,hx_pert                  &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                     
 ! finished innovation for a v obs      
  print*,'finished innovation for a v obs'               
  endif
  
! assimilate w 
 if(ivar == 3  .and. assim_w == 1 .and. obsw(iobs) .ne. miss_data) then
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position             
  call ijk_to_istat(ivar,xyzw(1,iobs)  ,xyzw(2,iobs)  ,xyzw(3,iobs)  ,istatnum1,nx,ny,nz)  
  call ijk_to_istat(ivar,xyzw(1,iobs)  ,xyzw(2,iobs)  ,xyzw(3,iobs)+1,istatnum2,nx,ny,nz)    
 ! cal the hx and hxmean for w   
  call obs_sim_grid_w(xstat(istatnum1,:),xstat(istatnum2,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
   temp2(i)=xstat(istatnum2,i)+xstat_mean(istatnum2)
  enddo 
   call obs_sim_grid_w(temp1,temp2,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsw(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzw(1,iobs),xyzw(2,iobs),xyzw(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)   
 ! cal kalman maxtria          
  call kalman_cal(zzh,hph,obs_err_w,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_w)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,     &
                      xstat,xstat_mean,kalman,hx_pert                  &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a w obs    
  print*,'finished innovation for a w obs'                     
  endif
  
! assimilate ph 
 if(ivar == 4 .and. assim_ph == 1 .and. obsph(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position  
  call ijk_to_istat(ivar,xyzph(1,iobs)  ,xyzph(2,iobs)  ,xyzph(3,iobs)  ,istatnum1,nx,ny,nz)  
  call ijk_to_istat(ivar,xyzph(1,iobs)  ,xyzph(2,iobs)  ,xyzph(3,iobs)+1,istatnum2,nx,ny,nz)               
 ! cal the hx and hxmean for ph   
  call obs_sim_grid_ph(xstat(istatnum1,:),xstat(istatnum2,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
   temp2(i)=xstat(istatnum2,i)+xstat_mean(istatnum2)
  enddo 
   call obs_sim_grid_ph(temp1,temp2,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsph(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzph(1,iobs),xyzph(2,iobs),xyzph(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)   
 ! cal kalman maxtria           
  call kalman_cal(zzh,hph,obs_err_ph,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_ph)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,     &
                      xstat,xstat_mean,kalman,hx_pert                  &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a ph obs 
  print*,'finished innovation for a ph obs'                       
  endif
  
! assimilate t 
 if(ivar == 5 .and. assim_t == 1 .and. obst(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position    
  call ijk_to_istat(ivar,xyzt(1,iobs)  ,xyzt(2,iobs)  ,xyzt(3,iobs)  ,istatnum1,nx,ny,nz)           
 ! cal the hx and hxmean for t    
  call obs_sim_grid_mass(xstat(istatnum1,:),hx_pert,ens_num) 
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
  enddo 
   call obs_sim_grid_mass(temp1,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obst(iobs),hxmean)  
  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzt(1,iobs),xyzt(2,iobs),xyzt(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)    
 ! cal kalman maxtria          
  call kalman_cal(zzh,hph,obs_err_t,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_t)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,     &
                      xstat,xstat_mean,kalman,hx_pert                  &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a t obs  
  print*,'finished innovation for a t obs'                        
  endif
  
! assimilate qv 
 if(ivar == 6 .and. assim_qv == 1 .and. obsqv(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position    
  call ijk_to_istat(ivar,xyzqv(1,iobs)  ,xyzqv(2,iobs)  ,xyzqv(3,iobs)  ,istatnum1,nx,ny,nz)             
 ! cal the hx and hxmean for qv   
  call obs_sim_grid_mass(xstat(istatnum1,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
  enddo 
   call obs_sim_grid_mass(temp1,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsqv(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzqv(1,iobs),xyzqv(2,iobs),xyzqv(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)   
 ! cal kalman maxtria          
  call kalman_cal(zzh,hph,obs_err_qv,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_qv)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,     &
                      xstat,xstat_mean,kalman,hx_pert                  &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a qv obs  
  print*,'finished innovation for a qv obs'                       
  endif
  
! assimilate qr 
 if(ivar == 7 .and. assim_qr == 1 .and. obsqr(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position 
  call ijk_to_istat(ivar,xyzqr(1,iobs)  ,xyzqr(2,iobs)  ,xyzqr(3,iobs)  ,istatnum1,nx,ny,nz)           
 ! cal the hx and hxmean for qr   
  call obs_sim_grid_mass(xstat(istatnum1,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
  enddo 
   call obs_sim_grid_mass(temp1,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsqr(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzqr(1,iobs),xyzqr(2,iobs),xyzqr(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 ! cal kalman maxtria              
  call kalman_cal(zzh,hph,obs_err_qr,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_qr)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,    &
                      xstat,xstat_mean,kalman,hx_pert                 &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a qr obs 
  print*,'finished innovation for a qr obs'                       
  endif
  
! assimilate qi 
 if(ivar == 8 .and. assim_qi == 1 .and. obsqi(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position      
  call ijk_to_istat(ivar,xyzqi(1,iobs)  ,xyzqi(2,iobs)  ,xyzqi(3,iobs)  ,istatnum1,nx,ny,nz)                                                                          
 ! cal the hx and hxmean for qi   
  call obs_sim_grid_mass(xstat(istatnum1,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
  enddo 
   call obs_sim_grid_mass(temp1,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsqi(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzqi(1,iobs),xyzqi(2,iobs),xyzqi(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 ! cal kalman maxtria                                                                          
  call kalman_cal(zzh,hph,obs_err_qi,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_qi)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,   &
                      xstat,xstat_mean,kalman,hx_pert                &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a qi obs   
  print*,'finished innovation for a qi obs'                  
  endif
  
! assimilate qs
 if(ivar == 9 .and. assim_qs == 1 .and. obsqs(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position           
  call ijk_to_istat(ivar,xyzqs(1,iobs)  ,xyzqs(2,iobs)  ,xyzqs(3,iobs)  ,istatnum1,nx,ny,nz)                                                                     
 ! cal the hx and hxmean for qs   
  call obs_sim_grid_mass(xstat(istatnum1,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
  enddo 
   call obs_sim_grid_mass(temp1,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo  
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsqs(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzqs(1,iobs),xyzqs(2,iobs),xyzqs(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 ! cal kalman maxtria                                                                          
  call kalman_cal(zzh,hph,obs_err_qs,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_qs)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,    &
                      xstat,xstat_mean,kalman,hx_pert                 &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a qs obs  
  print*,'finished innovation for a qs obs'                        
  endif
  
! assimilate qgr
 if(ivar == 10 .and. assim_qgr == 1 .and. obsqgr(iobs) .ne. miss_data) then  
    
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) 

 ! cal the position  
  call ijk_to_istat(ivar,xyzqgr(1,iobs)  ,xyzqgr(2,iobs)  ,xyzqgr(3,iobs)  ,istatnum1,nx,ny,nz)                                                                              
 ! cal the hx and hxmean for qgr   
  call obs_sim_grid_mass(xstat(istatnum1,:),hx_pert,ens_num)
  do i=1,ens_num                   
   temp1(i)=xstat(istatnum1,i)+xstat_mean(istatnum1)
  enddo 
   call obs_sim_grid_mass(temp1,hx,ens_num)
  hxmean=0
  do i=1,ens_num    
   hxmean=hxmean+hx(i)/ens_num
  enddo 
 ! cal the innovation vecter    
  call innovation_cal(innovationmean,obsqgr(iobs),hxmean)  
 ! cal the ZZH 
  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
 ! cal localization
  if(enable_local ==1) then
   call localization_cal_sim(zzh,numstat,xyzqgr(1,iobs),xyzqgr(2,iobs),xyzqgr(3,iobs), &
                             phb,analysis_var_num) 
  endif  
 ! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
 ! cal kalman maxtria                                                                          
  call kalman_cal(zzh,hph,obs_err_qgr,numstat,kalman)                                             
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_qgr)  
 ! cal analysis field 
  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,   &
                      xstat,xstat_mean,kalman,hx_pert                &
                     )  
  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
 ! finished innovation for a qgr obs
  print*,'finished innovation for a qgr obs'                           
  endif
  
 enddo   !  end loop for var    
print*,'assimilate the obs at the',iobs,'th point' 
 enddo   !  end loop for obs
 
print*,'finish the assimilation' 
 
 deallocate(zzh)
 deallocate(kalman)
 deallocate(hx) 
 deallocate(innovation)
   
end subroutine assimilate_sim_cal


subroutine assimilate_radar_cal(analysis_var_num,ra_data_num,numstat,                          &
                                xstat,xstat_mean,                                              &
                                rf,rv,ir,jr,hr,                                                &
                                xradar,yradar,hradar,lev_num,evl,                              &
                                xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                                rdnw,rdn,mub,mu                                                &
                               ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine assimilate_radar_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!#############################################################
!#          sim radar obs assimilation subroutine
!#             both RV and RF can be processed  
!#      this is a very simple test code with lots of bugs
!#############################################################
                                 
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
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real                  :: xstat(numstat,ens_num),xstat_mean(numstat)
  
  real,parameter        :: gravity=9.81                         
  real                  :: rf(ra_data_num),rv(ra_data_num),hr(ra_data_num)
  integer               :: ir(ra_data_num),jr(ra_data_num)
  integer               :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num  
  real                  :: azimuth
  real                  :: evlr  
  real                  :: evl(lev_num)
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
  tem1=real(jr(iobs)-yradar)/real(ir(iobs)-xradar)
  azimuth=atan(tem1)
  hor_dis= sqrt(((jr(iobs)-yradar)*dy)**2+((ir(iobs)-xradar)*dx)**2)
  evlr=atan( (hr(iobs)-hradar)/hor_dis  )
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr(iobs),phb,krw,krm)  
!************************************************************************
!print*,'azimuth,evlr',azimuth*180/3.1415926,evlr*180/3.1415926
  
if(assim_rv == 1 .and. rv(iobs) .ne. miss_data) then ! for radial velocity

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
                        ir(iobs),jr(iobs),hr(iobs),xradar,yradar,hradar,phb,           &
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
  
  call interpolation_sim_radar_wind(nx,ny,nz,ir(iobs),jr(iobs),krw,krm,hr(iobs),phb,     &
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
 
                                       
  call interpolation_sim_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr(iobs),phb,     &
                                    qru,qrl,qsu,qsl,qgru,qgrl,                       &
                                    rhou,rhol,tcu,tcl,                               &
                                    qrr,qsr,qgrr,rhor,tcr                            &
                                    )


  call ref_operator(qrr,qsr,qgrr,ref_x,ref_r,ref_s,ref_h,rhor,tcr)                                   
                                                                        
  call rv_operator(ur,vr,wr,hx(i),azimuth,evlr,ref_x,ref_r,ref_s,ref_h,   &
                   qrr,qsr,qgrr,rhor,uctb,vctb,wctb                       &
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
!   call localization_cal_radar_sim(zzh,numstat,ir(iobs),jr(iobs),hr(iobs),phb, &
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
  tem1=real(jr(iobs)-yradar)/real(ir(iobs)-xradar)
  azimuth=atan(tem1)
  hor_dis= sqrt(((jr(iobs)-yradar)*dy)**2+((ir(iobs)-xradar)*dx)**2)
  evlr=atan( (hr(iobs)-hradar)/hor_dis  )
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr(iobs),phb,krw,krm)  
!************************************************************************
!print*,'azimuth,evlr',azimuth*180/3.1415926,evlr*180/3.1415926  
  
if(assim_rf == 1 .and. rf(iobs) .ne. miss_data .and. rf(iobs) >= 10)  then ! for reflectivity

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
                        ir(iobs),jr(iobs),hr(iobs),xradar,yradar,hradar,phb,           &
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
 
                                       
  call interpolation_sim_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr(iobs),phb,     &
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
!   call localization_cal_radar_sim(zzh,numstat,ir(iobs),jr(iobs),hr(iobs),phb, &
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





subroutine assimilate_radar_cal2(analysis_var_num,ra_data_num,numstat,                          &
                                 xstat,xstat_mean,                                              &
                                 rf,rv,ir,jr,hr,                                                &
                                 xradar,yradar,hradar,lev_num,evl,                              &
                                 xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                                 rdnw,rdn,mub,mu,                                               &
                                 xstat2,xstat_mean2,mub2,mu2                                    &
                                 ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine assimilate_radar_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!#############################################################
!#          sim radar obs assimilation subroutine
!#             both RV and RF can be processed  
!#      this is a very simple test code with lots of bugs
!#############################################################
                                 
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
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real                  :: xstat(numstat,ens_num),xstat_mean(numstat)
  
  real,parameter        :: gravity=9.81                         
  real                  :: rf(ra_data_num),rv(ra_data_num),hr(ra_data_num)
  integer               :: ir(ra_data_num),jr(ra_data_num)
  integer               :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num  
  real                  :: azimuth
  real                  :: evlr  
  real                  :: evl(lev_num)
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

  real                  :: xstat2(numstat,ens_num),xstat_mean2(numstat)
  real                  :: mub2(nx,ny,ens_num),mu2(nx,ny,ens_num)
  real,allocatable      :: xstat_sub2(:,:),xstat_mean_sub2(:)
  real,allocatable      :: zzh2(:),kalman2(:)


!  allocate(zzh(numstat))
!  allocate(kalman(numstat))


print*,'start the assimilation of sim radar obs'
print*,ra_data_num


  do iobs=1,ra_data_num
   
!*********************prepare for HX cal*********************************
  tem1=real(jr(iobs)-yradar)/real(ir(iobs)-xradar)
  azimuth=atan(tem1)
  hor_dis= sqrt(((jr(iobs)-yradar)*dy)**2+((ir(iobs)-xradar)*dx)**2)
  evlr=atan( (hr(iobs)-hradar)/hor_dis  )
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr(iobs),phb,krw,krm)  
!************************************************************************
!print*,'azimuth,evlr',azimuth*180/3.1415926,evlr*180/3.1415926
  
if(assim_rv == 1 .and. rv(iobs) .ne. miss_data) then ! for radial velocity

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
                        ir(iobs),jr(iobs),hr(iobs),xradar,yradar,hradar,phb,           &
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

  allocate(zzh2(gridnum))
  allocate(kalman2(gridnum)) 
  allocate(xstat_sub2(gridnum,ens_num))
  allocate(xstat_mean_sub2(gridnum))
  
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
   
   
   xstat_mean_sub2(i)=0
   zzh2(i)=0
   kalman2(i)=0
   do j=1,ens_num
    xstat_sub2(i,j)=0
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
        xstat_sub2(nn,j)=xstat2(gridsgn(nn),j) 
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
  
  call interpolation_sim_radar_wind(nx,ny,nz,ir(iobs),jr(iobs),krw,krm,hr(iobs),phb,     &
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
 
                                       
  call interpolation_sim_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr(iobs),phb,     &
                                    qru,qrl,qsu,qsl,qgru,qgrl,                       &
                                    rhou,rhol,tcu,tcl,                               &
                                    qrr,qsr,qgrr,rhor,tcr                            &
                                    )


  call ref_operator(qrr,qsr,qgrr,ref_x,ref_r,ref_s,ref_h,rhor,tcr)                                   
                                                                        
  call rv_operator(ur,vr,wr,hx(i),azimuth,evlr,ref_x,ref_r,ref_s,ref_h,   &
                   qrr,qsr,qgrr,rhor,uctb,vctb,wctb                       &
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
  call xstat_mean_seperate_cal(xstat_sub2,xstat_mean_sub2,gridnum,ens_num)   
  
! quality control if innovationmean is great than threhold , assimilation will not be done
  if( abs(innovationmean) > ob_rv_threshold ) goto 101  
  print*,'innovationmean,rv(iobs),hxmean',innovationmean,rv(iobs),hxmean  
  
! cal the ZZH 
!  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)
  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)
  call zzh_cal(xstat_sub2,hx,hxmean,ens_num,gridnum,zzh2)
    
if(enable_rm_correlation_rv == 1 ) then   
! If the correlation between RV and state variation is not good or even lead to a worse result
! it should be remove.currently,it is found that rv will worsen the 1st level for all var


!  do ivar=1,analysis_var_num   
!  call rm_correlate(zzh,numstat,ivar,nx,ny,nz,1)
!  enddo  


  call rm_correlate2(zzh,gridnum,gridsgn,nx,ny,nz,1) 
  call rm_correlate2(zzh2,gridnum,gridsgn,nx,ny,nz,1)
   
endif

! cal localization
  if(enable_local ==1) then
  

!do istatmp=1,numstat 
!   call localization_cal_radar_sim(zzh,numstat,ir(iobs),jr(iobs),hr(iobs),phb, &
!                                   analysis_var_num,istatmp)                                      
!enddo                          


  call simple_local(zzh,gridnum,dis_rho)
  call simple_local(zzh2,gridnum,dis_rho)
      
  endif  

! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  

  
! cal kalman maxtria                                                                          
!  call kalman_cal(zzh,hph,obs_err_rv,numstat,kalman)
  call kalman_cal(zzh,hph,obs_err_rv,gridnum,kalman) 
  call kalman_cal(zzh2,hph,obs_err_rv,gridnum,kalman2) 
                                                
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
                     
  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,              &
                      xstat_sub2,xstat_mean_sub2,kalman2,hx_pert                &
                     )  
                                                              
! finished innovation for a rv obs
  print*,'finished innovation for a rv obs'

101 continue 
!  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)  
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)
 
  call xstat_mean_merge_cal(xstat_sub2,xstat_mean_sub2,gridnum,ens_num)                        
  call merge_local(xstat_sub2,xstat2,gridsgn,gridnum,ens_num,numstat) 
  
  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman) 
  deallocate(hx) 
  deallocate(hx_pert)  
  
  deallocate(zzh2)
  deallocate(kalman2) 
  deallocate(xstat_sub2)
  deallocate(xstat_mean_sub2)

print*,'assimilate the obs at the',iobs,'th point'     
  endif       ! for radial velocity
                               
  enddo  !  end loop for obs      


  do iobs=1,ra_data_num
   
!*********************prepare for HX cal*********************************
  tem1=real(jr(iobs)-yradar)/real(ir(iobs)-xradar)
  azimuth=atan(tem1)
  hor_dis= sqrt(((jr(iobs)-yradar)*dy)**2+((ir(iobs)-xradar)*dx)**2)
  evlr=atan( (hr(iobs)-hradar)/hor_dis  )
  
! krw and krw is the level just lower than hr      
  call radardata_kr_cal(nx,ny,nz,ir(iobs),jr(iobs),hr(iobs),phb,krw,krm)  
!************************************************************************
!print*,'azimuth,evlr',azimuth*180/3.1415926,evlr*180/3.1415926  
  
if(assim_rf == 1 .and. rf(iobs) .ne. miss_data .and. rf(iobs) >= 10)  then ! for reflectivity

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
                        ir(iobs),jr(iobs),hr(iobs),xradar,yradar,hradar,phb,           &
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
  
  allocate(zzh2(gridnum))
  allocate(kalman2(gridnum)) 
  allocate(xstat_sub2(gridnum,ens_num))
  allocate(xstat_mean_sub2(gridnum))  
  
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
   
   
   xstat_mean_sub2(i)=0
   zzh2(i)=0
   kalman2(i)=0
   do j=1,ens_num
    xstat_sub2(i,j)=0
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
        xstat_sub2(nn,j)=xstat2(gridsgn(nn),j) 
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
 
                                       
  call interpolation_sim_radar_mass(nx,ny,nz,ir(iobs),jr(iobs),krm,hr(iobs),phb,     &
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
  call xstat_mean_seperate_cal(xstat_sub2,xstat_mean_sub2,gridnum,ens_num) 
  
            
! quality control if innovationmean is much great than threhold , assimilation will not be done
  if( abs(innovationmean) > ob_rf_threshold ) goto 102
  print*,'innovationmean,rf(iobs),hxmean',innovationmean,rf(iobs),hxmean

! cal the ZZH 
!  call zzh_cal(xstat,hx,hxmean,ens_num,numstat,zzh)  
  call zzh_cal(xstat_sub,hx,hxmean,ens_num,gridnum,zzh)
  call zzh_cal(xstat_sub2,hx,hxmean,ens_num,gridnum,zzh2)
  
if(enable_rm_correlation_rf == 1 ) then  
! If the correlation between RF and state variation is not good or even lead to a worse result
! it should be remove.currently,it is found that rf will worsen the 1st level for all var


!  do ivar=1,analysis_var_num    
!  call rm_correlate(zzh,numstat,ivar,nx,ny,nz,2)
!  enddo
 

  call rm_correlate2(zzh,gridnum,gridsgn,nx,ny,nz,2) 
  call rm_correlate2(zzh2,gridnum,gridsgn,nx,ny,nz,2) 
  
endif
    
! cal localization
  if(enable_local ==1) then
  
  
!do istatmp=1,numstat 
!   call localization_cal_radar_sim(zzh,numstat,ir(iobs),jr(iobs),hr(iobs),phb, &
!                                   analysis_var_num,istatmp)                                      
!enddo                                 


  call simple_local(zzh,gridnum,dis_rho)
  call simple_local(zzh2,gridnum,dis_rho)
  
  endif  

! cal hph 
  call hph_cal(hx,hxmean,hph,ens_num)  
  
! cal kalman maxtria                                                                          
!  call kalman_cal(zzh,hph,obs_err_rf,numstat,kalman) 
  call kalman_cal(zzh,hph,obs_err_rf,gridnum,kalman)
  call kalman_cal(zzh2,hph,obs_err_rv,gridnum,kalman2) 
                                              
  call kalman_alpha_cal(kalman_alpha,hph,obs_err_rf) 
   
! cal analysis field 
!  call innovate_xstat(ens_num,numstat,innovationmean,kalman_alpha,   &
!                      xstat,xstat_mean,kalman,hx_pert                &
!                     )  
                     
  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,           &
                      xstat_sub,xstat_mean_sub,kalman,hx_pert                &
                     ) 
                     
  call innovate_xstat(ens_num,gridnum,innovationmean,kalman_alpha,              &
                      xstat_sub2,xstat_mean_sub2,kalman2,hx_pert                &
                     )  
                                                                
                     
! finished innovation for a rf obs
  print*,'finished innovation for a rf obs'
                       
102 continue                      
!  call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                      
                    
  call xstat_mean_merge_cal(xstat_sub,xstat_mean_sub,gridnum,ens_num)                        
  call merge_local(xstat_sub,xstat,gridsgn,gridnum,ens_num,numstat)

  call xstat_mean_merge_cal(xstat_sub2,xstat_mean_sub2,gridnum,ens_num)                        
  call merge_local(xstat_sub2,xstat2,gridsgn,gridnum,ens_num,numstat)

  deallocate(xstat_sub)
  deallocate(xstat_mean_sub)
  deallocate(dis_rho)
  deallocate(gridsgn)
  deallocate(zzh)
  deallocate(kalman) 
  deallocate(hx) 
  deallocate(hx_pert)  
  
  deallocate(zzh2)
  deallocate(kalman2) 
  deallocate(xstat_sub2)
  deallocate(xstat_mean_sub2)

print*,'assimilate the obs at the',iobs,'th point'      
  
  endif       ! for reflectivity 
                          
  enddo  !  end loop for obs                         

                             
end subroutine assimilate_radar_cal2