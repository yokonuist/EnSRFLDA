  program ensrf_sim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            program ensrf_sim !!! for simulate storms ,maybe not fit for SMOIS & TSLB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
 !----------------------------------------------------------------
 !Purpose: 
 !        Main program of ENSRF
 !
 !---------------------------------------------------------------- 
  

  implicit none
  include 'netcdf.inc'
  include 'namelist.inc'
   
  character (len=500)   :: input_file_name
  integer               :: length_input,length_file_dir,length_ens_head  
  integer               :: length_obs_dir      
  integer,parameter     :: analysis_var_num=10
  integer,parameter     :: basic_var_num_1=13
  integer,parameter     :: basic_var_num_2=2
  integer               :: i,j,k,ii,jj,kk
  integer               :: istart(4), iend(4)
  character (len=4)     :: lab_ens
  character (len=2)     :: lab_domain
  character (len=2)     :: lab_time_num
  character (len=2)     :: lab_radar_num
  
  integer               :: rec_id_var(analysis_var_num)
  integer,allocatable   :: rec_cdfid(:)
  integer               :: rec_dims(analysis_var_num,3)
  integer,allocatable   :: rec_rcode(:)
  character(len=80)     :: analysis_var_name(analysis_var_num)
  integer               :: file_cdfid,file_rcode
  character(len=80)     :: basic_var_name_1(basic_var_num_1)
  character(len=80)     :: basic_var_name_2(basic_var_num_2)
    
  real,allocatable      :: u(:,:,:),v(:,:,:),w(:,:,:),ph(:,:,:),t(:,:,:),qv(:,:,:),  &
                           qr(:,:,:),qi(:,:,:),qs(:,:,:),qgr(:,:,:) 
                                           
  real,allocatable      :: ensu(:,:,:,:),ensv(:,:,:,:),ensw(:,:,:,:),ensph(:,:,:,:),enst(:,:,:,:),ensqv(:,:,:,:), &
                           ensqr(:,:,:,:),ensqi(:,:,:,:),ensqs(:,:,:,:),ensqgr(:,:,:,:)    

                           
  real,allocatable      :: umean(:,:,:),vmean(:,:,:),wmean(:,:,:),phmean(:,:,:),tmean(:,:,:),qvmean(:,:,:),  &
                           qrmean(:,:,:),qimean(:,:,:),qsmean(:,:,:),qgrmean(:,:,:)       
  
  real,allocatable      :: utrue(:,:,:),vtrue(:,:,:),wtrue(:,:,:),phtrue(:,:,:),ttrue(:,:,:),qvtrue(:,:,:),  &
                           qrtrue(:,:,:),qitrue(:,:,:),qstrue(:,:,:),qgrtrue(:,:,:)
  real,allocatable      :: mubtrue(:,:),mutrue(:,:)                            
  
  real,allocatable      :: xlon(:,:),xlat(:,:),ulon(:,:),ulat(:,:),vlon(:,:),vlat(:,:),phb(:,:,:),hgt(:,:),mub(:,:,:),mu(:,:,:), &
                           znw(:),znu(:),rdnw(:),rdn(:)
  real                  :: p_top
  
  real,allocatable      :: obsu(:),obsv(:),obsw(:),obsph(:),obst(:),obsqv(:),obsqr(:),obsqi(:),obsqs(:),obsqgr(:)
  integer,allocatable   :: xyzu(:,:),xyzv(:,:),xyzw(:,:),xyzph(:,:),xyzt(:,:),xyzqv(:,:),xyzqr(:,:),xyzqi(:,:),xyzqs(:,:),xyzqgr(:,:)  

  character(len=5)      :: radar_name
  integer               :: xradar,yradar !
  real                  :: hradar
  integer               :: lev_num  
  real,allocatable      :: rv(:),rf(:),hr(:)
  integer,allocatable   :: ir(:),jr(:)
  real,allocatable      :: evl(:)
  real                  :: ref_r,ref_s,ref_h
  integer               :: iradar
  real,allocatable      :: rho(:,:,:),ref_x(:,:,:),p(:,:,:),tc(:,:,:)
  
  integer               :: obs_sim_num !
  integer               :: ra_data_num !
  
  real,allocatable      :: xstat(:,:),xstat_mean(:),xstatf(:,:),xstat_temp_mean(:) 

  
  integer               :: numstat 
  integer               :: istat
  real                  :: rmse_rv,spd_rv,rmse_rf,spd_rf
  
  real,allocatable      :: stdu2d(:),stdv2d(:),stdw2d(:),stdph2d(:),stdt2d(:),stdqv2d(:), &
                           stdqr2d(:),stdqi2d(:),stdqs2d(:),stdqgr2d(:)
  real                  :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr
  real                  :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr  
  integer               :: ireg  
  integer               :: temp1,temp2,temp3

!*******************************
!      TEST CODE
!*******************************
  integer               :: rec_id_var2(analysis_var_num)
  integer,allocatable   :: rec_cdfid2(:)
  integer               :: rec_dims2(analysis_var_num,3)
  integer,allocatable   :: rec_rcode2(:)
  character (len=2)     :: lab_time_num2
 
  real,allocatable      :: ensu2(:,:,:,:),ensv2(:,:,:,:),ensw2(:,:,:,:),ensph2(:,:,:,:),enst2(:,:,:,:),ensqv2(:,:,:,:), &
                           ensqr2(:,:,:,:),ensqi2(:,:,:,:),ensqs2(:,:,:,:),ensqgr2(:,:,:,:),mub2(:,:,:),mu2(:,:,:) 
  real,allocatable      :: xstat2(:,:),xstat_mean2(:),xstatf2(:,:),xstat_temp_mean2(:)
  
  real,allocatable      :: utrue2(:,:,:),vtrue2(:,:,:),wtrue2(:,:,:),phtrue2(:,:,:),ttrue2(:,:,:),qvtrue2(:,:,:),  &
                           qrtrue2(:,:,:),qitrue2(:,:,:),qstrue2(:,:,:),qgrtrue2(:,:,:) 
  real,allocatable      :: mubtrue2(:,:),mutrue2(:,:)           
!*******************************
!*******************************  
  
              
  data analysis_var_name/'U','V','W','PH','T','QVAPOR','QRAIN','QICE','QSNOW','QGRAUP'/
  
  data basic_var_name_1/'XLONG','XLAT','XLONG_U','XLAT_U','XLONG_V','XLAT_V','PHB','HGT','P_TOP','ZNW','ZNU',  &
                        'RDNW','RDN'/
                        
  data basic_var_name_2/'MUB','MU'/
  
  open(unit=5, file="ensrf.input", form="formatted", status="old")
  
  read (5,nml = input_file_info)
  write(6,nml = input_file_info)
  read (5,nml =  dimension_info)
  write(6,nml =  dimension_info)
  read (5,nml =   ensemble_info)
  write(6,nml =   ensemble_info)
  read (5,nml =   analysis_info)
  write(6,nml =   analysis_info)  
  read (5,nml =        obs_info)
  write(6,nml =        obs_info)  
  close(5)

  allocate(rec_cdfid(ens_num))  
  allocate(rec_rcode(ens_num))  

  allocate(u  (nx,ny,nz))
  allocate(v  (nx,ny,nz))
  allocate(w  (nx,ny,nz))
  allocate(ph (nx,ny,nz))
  allocate(t  (nx,ny,nz))
  allocate(qv (nx,ny,nz))  
  allocate(qr (nx,ny,nz))
  allocate(qs (nx,ny,nz))
  allocate(qi (nx,ny,nz))
  allocate(qgr(nx,ny,nz))
  
  allocate(ensu  (nx,ny,nz,ens_num))
  allocate(ensv  (nx,ny,nz,ens_num))
  allocate(ensw  (nx,ny,nz,ens_num))
  allocate(ensph (nx,ny,nz,ens_num))
  allocate(enst  (nx,ny,nz,ens_num))
  allocate(ensqv (nx,ny,nz,ens_num))   
  allocate(ensqr (nx,ny,nz,ens_num)) 
  allocate(ensqi (nx,ny,nz,ens_num)) 
  allocate(ensqs (nx,ny,nz,ens_num)) 
  allocate(ensqgr(nx,ny,nz,ens_num))
    
  allocate(umean  (nx,ny,nz))
  allocate(vmean  (nx,ny,nz))
  allocate(wmean  (nx,ny,nz))
  allocate(phmean (nx,ny,nz))
  allocate(tmean  (nx,ny,nz))
  allocate(qvmean (nx,ny,nz))     
  allocate(qrmean (nx,ny,nz))
  allocate(qimean (nx,ny,nz))
  allocate(qsmean (nx,ny,nz))
  allocate(qgrmean(nx,ny,nz))

  allocate(utrue  (nx,ny,nz))
  allocate(vtrue  (nx,ny,nz))
  allocate(wtrue  (nx,ny,nz))
  allocate(phtrue (nx,ny,nz))
  allocate(ttrue  (nx,ny,nz))
  allocate(qvtrue (nx,ny,nz))     
  allocate(qrtrue (nx,ny,nz))
  allocate(qitrue (nx,ny,nz))
  allocate(qstrue (nx,ny,nz))
  allocate(qgrtrue(nx,ny,nz)) 

  allocate(xlon(nx,ny))
  allocate(xlat(nx,ny))
  allocate(ulon(nx,ny))
  allocate(ulat(nx,ny))
  allocate(vlon(nx,ny))
  allocate(vlat(nx,ny)) 
  allocate(phb(nx,ny,nz))
  allocate(hgt(nx,ny)) 
  allocate(mub(nx,ny,ens_num)) 
  allocate(mu (nx,ny,ens_num)) 
  allocate(znw   (nz)) 
  allocate(znu   (nz)) 
  allocate(rdnw  (nz)) 
  allocate(rdn   (nz))   
  
  allocate(stdu2d  (nz))
  allocate(stdv2d  (nz))
  allocate(stdw2d  (nz))
  allocate(stdph2d (nz))
  allocate(stdt2d  (nz))
  allocate(stdqv2d (nz))  
  allocate(stdqr2d (nz))
  allocate(stdqi2d (nz))
  allocate(stdqs2d (nz))  
  allocate(stdqgr2d(nz))
 
  allocate(rho(nx,ny,nz))
  allocate(ref_x(nx,ny,nz))
  allocate(p(nx,ny,nz))
  allocate(tc(nx,ny,nz))  
  allocate(mubtrue(nx,ny)) 
  allocate(mutrue (nx,ny))   

if(enable_rckf_gi==1) then

  allocate(ensu2  (nx,ny,nz,ens_num))
  allocate(ensv2  (nx,ny,nz,ens_num))
  allocate(ensw2  (nx,ny,nz,ens_num))
  allocate(ensph2 (nx,ny,nz,ens_num))
  allocate(enst2  (nx,ny,nz,ens_num))
  allocate(ensqv2 (nx,ny,nz,ens_num))   
  allocate(ensqr2 (nx,ny,nz,ens_num)) 
  allocate(ensqi2 (nx,ny,nz,ens_num)) 
  allocate(ensqs2 (nx,ny,nz,ens_num)) 
  allocate(ensqgr2(nx,ny,nz,ens_num)) 
  allocate(mub2(nx,ny,ens_num)) 
  allocate(mu2 (nx,ny,ens_num))   
  allocate(utrue2  (nx,ny,nz))
  allocate(vtrue2  (nx,ny,nz))
  allocate(wtrue2  (nx,ny,nz))
  allocate(phtrue2 (nx,ny,nz))
  allocate(ttrue2  (nx,ny,nz))
  allocate(qvtrue2 (nx,ny,nz))     
  allocate(qrtrue2 (nx,ny,nz))
  allocate(qitrue2 (nx,ny,nz))
  allocate(qstrue2 (nx,ny,nz))
  allocate(qgrtrue2(nx,ny,nz)) 
  allocate(mubtrue2(nx,ny)) 
  allocate(mutrue2 (nx,ny)) 
  allocate(rec_cdfid2(ens_num))  
  allocate(rec_rcode2(ens_num))    
endif 
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading true files
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------   

if(enable_rckf_gi==1) then

call read_true_file2(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                              &
                     rec_dims,file_rcode,                                                                   &
                     utrue2,vtrue2,wtrue2,phtrue2,ttrue2,qvtrue2,qrtrue2,qitrue2,qstrue2,qgrtrue2,          &  
                     mubtrue2,mutrue2,basic_var_num_2,basic_var_name_2,                                     &
                     xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,rdnw,rdn,                          &
                     basic_var_num_1,basic_var_name_1,                                                      &
                     lab_domain,lab_time_num2                                                               &             
                    )    
  print*,'finished reading true fields2' 
endif  


call read_true_file(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                    &
                    rec_dims,file_rcode,                                                         &
                    utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,          &  
                    mubtrue,mutrue,basic_var_num_2,basic_var_name_2,                             &
                    xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,rdnw,rdn,                &
                    basic_var_num_1,basic_var_name_1,                                            &
                    lab_domain,lab_time_num                                                      &             
                   )      
                                          
                     
  print*,'finished reading true fields'                                 
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading ens files
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------   
                  
if(enable_rckf_gi==1) then

call read_ens_file2(analysis_var_num,analysis_var_name,rec_id_var2,file_cdfid,                   &
                    rec_dims2,file_rcode,rec_cdfid2,rec_rcode2,                                  &
                    ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2,          & 
                    mub2,mu2,basic_var_num_2,basic_var_name_2,                                   &
                    lab_domain,lab_time_num2                                                     &   
                    )
  print*,'finished reading ens fields2'                       
endif

  
call read_ens_file(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                    &
                   rec_dims,file_rcode,rec_cdfid,rec_rcode,                                     &
                   ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr                     & 
                   mub,mu,basic_var_num_2,basic_var_name_2,                                     &
                   lab_domain,lab_time_num                                                      &   
                  )                  
                          
  print*,'finished reading ens fields'                        


!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Processing the analysis  326-379
!
!---------------------------------------------------------------------------------------  
!--------------------------------------------------------------------------------------- 
numstat=nx*ny*nz*analysis_var_num

allocate(xstat(numstat,ens_num))
allocate(xstat_mean(numstat))
allocate(xstatf(numstat,ens_num))
allocate(xstat_temp_mean(numstat) )

if(enable_rckf_gi==1) then
  allocate(xstat2(numstat,ens_num))
  allocate(xstat_mean2(numstat))
  allocate(xstatf2(numstat,ens_num))
  allocate(xstat_temp_mean2(numstat) )
endif
              
if( enable_rms_anal == 1 ) then

!=======================================================================
!          THIS STEP IS USED FOR RMS_CONVECTIVE_REGION_CAL
!=======================================================================
  do k=1,nz-1
  do j=1,ny
  do i=1,nx  
    call cal_rho(mubtrue(i,j),mutrue(i,j),qvtrue(i,j,k),phtrue(i,j,k),phtrue(i,j,k+1),rdnw(k),       &
                 p_top,znu(k),ttrue(i,j,k),nx,ny,nz,rho(i,j,k),p(i,j,k)                              &                  
                )
  enddo
  enddo
  enddo

  do i=1,nx
  do j=1,ny
    rho(i,j,nz)=rho(i,j,nz-1)
    p  (i,j,nz)=p_top
  enddo
  enddo

  do k=1,nz-1
  do j=1,ny
  do i=1,nx                  
     call cal_tc(ttrue(i,j,k),p(i,j,k),tc(i,j,k),nx,ny,nz)  
  enddo
  enddo
  enddo
  
  ref_x=0

  do k=1,nz
  do j=1,ny
  do i=1,nx     
     call ref_operator(qrtrue(i,j,k),qstrue(i,j,k),qgrtrue(i,j,k),ref_x(i,j,k),ref_r,ref_s,ref_h,rho(i,j,k),tc(i,j,k))    
         
  enddo
  enddo
  enddo 


  
!=======================================================================
!                      CALCULATE THE RMS
!=======================================================================
  ireg=0
  print*,'cal the mean '
  call mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,               &
                nx,ny,nz,ens_num,                                                       &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean      &
               )
  
  print*,'cal the rms'   
  if(opt_rms == 1) then
  call rms_cal(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &               
               nx,ny,nz,ens_num,                                                                  &
               umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
               stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
               stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d       &
              )
  else if(opt_rms == 2) then
  call rms_convetive_region(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &
                            nx,ny,nz,ens_num,                                                                  &
                            umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
                            stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                            stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                            ref_x                                                                              &
                            )
  endif         
           
  call output_rmse(stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                   stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                   lab_time_num,lab_domain,ireg,nx,ny,nz)
  
  print*,'cal the spread'
  if(opt_rms == 1) then   
  call ens_spread(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,               &
                  nx,ny,nz,ens_num,                                                       &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,     &
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr     &
                 )
  else if(opt_rms == 2) then 
  call ens_spread_convetive_region(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,               &
                                   nx,ny,nz,ens_num,                                                       &
                                   umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,     &
                                   spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                                   ref_x                                                                   & 
                                   )
  endif
                 
  call output_spd(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                  lab_time_num,lab_domain,ireg,nx,ny,nz)     
  
  call output_ratio_r2s(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                        stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,               &
                        lab_time_num,lab_domain,ireg,nx,ny,nz)                      





  if(enable_rckf_gi==100) then
  
  print*,'cal the mean2 '
  call mean_cal(ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2,               &
                nx,ny,nz,ens_num,                                                                 &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean                &
               )
  
  print*,'cal the rms2'   
  if(opt_rms == 1) then
  call rms_cal(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &               
               nx,ny,nz,ens_num,                                                                  &
               umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
               stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
               stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d       &
              )
  else if(opt_rms == 2) then
  call rms_convetive_region(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &
                            nx,ny,nz,ens_num,                                                                  &
                            umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
                            stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                            stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                            ref_x                                                                              &
                            )
  endif         
           
  call output_rmse(stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                   stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                   lab_time_num2,lab_domain,ireg,nx,ny,nz)
  
  print*,'cal the spread2'
  if(opt_rms == 1) then   
  call ens_spread(ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2,               &
                  nx,ny,nz,ens_num,                                                                 &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,               &
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr                &
                 )
  else if(opt_rms == 2) then 
  call ens_spread_convetive_region(ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2,     &
                                   nx,ny,nz,ens_num,                                                       &
                                   umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,     &
                                   spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                                   ref_x                                                                   & 
                                   )
  endif
                 
  call output_spd(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                  lab_time_num2,lab_domain,ireg,nx,ny,nz)     
  
  call output_ratio_r2s(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                        stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,               &
                        lab_time_num2,lab_domain,ireg,nx,ny,nz) 
                        
  endif                        
     
endif                                      
!---------------------------------------------------------------------------------------  
! analysis for simulated observations
!---------------------------------------------------------------------------------------
if(enable_sim == 1) then  

call fill_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,           &
                ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr   &
               )
               
deallocate(ensu)
deallocate(ensv)  
deallocate(ensw)  
deallocate(ensph)  
deallocate(enst)  
deallocate(ensqv)  
deallocate(ensqr)  
deallocate(ensqi)  
deallocate(ensqs)  
deallocate(ensqgr)              
deallocate(u  )
deallocate(v  )
deallocate(w  )
deallocate(ph )
deallocate(t  )
deallocate(qv )  
deallocate(qr )
deallocate(qs )
deallocate(qi )
deallocate(qgr)   
deallocate(umean  )
deallocate(vmean  )
deallocate(wmean  )
deallocate(phmean )
deallocate(tmean  )
deallocate(qvmean )  
deallocate(qrmean )
deallocate(qsmean )
deallocate(qimean )
deallocate(qgrmean) 


if(enable_rckf_gi==1) then

  call fill_xstat(xstat2,ens_num,numstat,nx,ny,nz,analysis_var_num,                    &
     	          ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2   &
     	          )
     	          
  deallocate(ensu2)
  deallocate(ensv2)  
  deallocate(ensw2)  
  deallocate(ensph2)  
  deallocate(enst2)  
  deallocate(ensqv2)  
  deallocate(ensqr2)  
  deallocate(ensqi2)  
  deallocate(ensqs2)  
  deallocate(ensqgr2)  
	
endif               
               
            
if(obs_type == 0 .or. obs_type == 1) then ! sim sim obs   

!---------------------------------------------------------------------------------------  
!                       Reading observations
!---------------------------------------------------------------------------------------  

call read_sim_sound_dimension(lab_time_num,obs_sim_num)
  
  allocate  (obsu(obs_sim_num))
  allocate  (obsv(obs_sim_num))
  allocate  (obsw(obs_sim_num))
  allocate (obsph(obs_sim_num))
  allocate  (obst(obs_sim_num))
  allocate (obsqv(obs_sim_num))
  allocate (obsqr(obs_sim_num))
  allocate (obsqi(obs_sim_num))
  allocate (obsqs(obs_sim_num))
  allocate(obsqgr(obs_sim_num))
  
  allocate  (xyzu(3,obs_sim_num))
  allocate  (xyzv(3,obs_sim_num))
  allocate  (xyzw(3,obs_sim_num))
  allocate (xyzph(3,obs_sim_num))
  allocate  (xyzt(3,obs_sim_num))
  allocate (xyzqv(3,obs_sim_num))
  allocate (xyzqr(3,obs_sim_num))
  allocate (xyzqi(3,obs_sim_num))
  allocate (xyzqs(3,obs_sim_num))
  allocate(xyzqgr(3,obs_sim_num))
  
  
call read_sim_sound(xyzu,xyzv,xyzw,xyzph,xyzt,xyzqv,xyzqr,xyzqs,xyzqi,xyzqgr,      &
                    obsu,obsv,obsw,obsph,obst,obsqv,obsqr,obsqs,obsqi,obsqgr,      &
                    obs_sim_num)

  close(1001)
  close(1002)
  close(1003)
  close(1004)
  close(1005)
  close(1006)
  close(1007)
  close(1008)
  close(1009)
  close(1010)
  close(1011)





if(enable_inflat==1 ) then
 print*,'calling inflation_sim_cal'
 call inflation_sim_cal(xstat,xstat_mean,numstat)
 print*,'finished inflation_sim_cal'
endif  


call assimilate_sim_cal(analysis_var_num,obs_sim_num,numstat,                          &
                        xstat,xstat_mean,                                              &
                        xyzu,xyzv,xyzw,xyzph,xyzt,xyzqv,xyzqr,xyzqs,xyzqi,xyzqgr,      &
                        obsu,obsv,obsw,obsph,obst,obsqv,obsqr,obsqs,obsqi,obsqgr,      &
                        xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                        rdnw,rdn,mub,mu                                                &                       
                       )    
                       
if( enable_inflat == 1 .and. inflat_opt == 2 ) then

call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
xstat=(1-ratio_af)*xstatf+ratio_af*xstat
call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num)                    

endif     
                       
endif    ! sim sim obs 



if(obs_type == 0 .or. obs_type == 2) then ! sim radar obs
!---------------------------------------------------------------------------------------  
!                       Reading observations
!---------------------------------------------------------------------------------------  

call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
xstatf=xstat
call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 

if(enable_rckf_gi==1) then
call xstat_mean_seperate_cal(xstat2,xstat_mean2,numstat,ens_num)
xstatf2=xstat2
call xstat_mean_merge_cal(xstat2,xstat_mean2,numstat,ens_num) 
endif

do iradar=1,radar_num    !   for radar number loop  

call read_sim_radar_dimension(lab_time_num,iradar,radar_name,xradar,yradar,     &
                              hradar,lev_num,ra_data_num)

if(allocated(rv)) deallocate(rv)
allocate(rv(ra_data_num))
if(allocated(rf)) deallocate(rf)
allocate(rf(ra_data_num))
if(allocated(ir)) deallocate(ir)
allocate(ir(ra_data_num))
if(allocated(jr)) deallocate(jr)
allocate(jr(ra_data_num))
if(allocated(hr)) deallocate(hr)
allocate(hr(ra_data_num))

call read_sim_radar(ra_data_num,rv,rf,ir,jr,hr)

close(1081)
close(1082)
close(1083)
close(1084)
close(1085)

if(opt_rms_rad   == 1) then
call innoation_rmse_rv(analysis_var_num,ra_data_num,numstat,                          &
                       xstat,                                                         &
                       rv,ir,jr,hr,                                                   &
                       xradar,yradar,hradar,lev_num,evl,                              &
                       xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                       rdnw,rdn,mub,mu,                                               &
                       rmse_rv,spd_rv                                                 &
                      )                      
call innoation_rmse_rf(analysis_var_num,ra_data_num,numstat,                          &
                       xstat,                                                         &
                       rf,ir,jr,hr,                                                   &
                       xradar,yradar,hradar,lev_num,evl,                              &
                       xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                       rdnw,rdn,mub,mu,                                               &
                       rmse_rf,spd_rf                                                 &
                       )                      
call output_radio_r2s_rvrf(rmse_rv,spd_rv,rmse_rf,spd_rf,iradar,lab_time_num,0)

endif


if(enable_inflat==1 .and. inflat_opt == 1 ) then
 print*,'calling inflation_sim_cal'
 call inflation_sim_cal(xstat,xstat_mean,numstat)
 print*,'finished inflation_sim_cal'
endif  

if(enable_rckf_gi==0) then
call assimilate_radar_cal(analysis_var_num,ra_data_num,numstat,                          &
                          xstat,xstat_mean,                                              &
                          rf,rv,ir,jr,hr,                                                &
                          xradar,yradar,hradar,lev_num,evl,                              &
                          xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                          rdnw,rdn,mub,mu                                                &
                          ) 
else
call assimilate_radar_cal2(analysis_var_num,ra_data_num,numstat,                          &
                           xstat,xstat_mean,                                              &
                           rf,rv,ir,jr,hr,                                                &
                           xradar,yradar,hradar,lev_num,evl,                              &
                           xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                           rdnw,rdn,mub,mu,                                               &
                           xstat2,xstat_mean2,mub2,mu2                                    &
                           ) 
endif

if(opt_rms_rad   == 1) then
call innoation_rmse_rv(analysis_var_num,ra_data_num,numstat,                          &
                       xstat,                                                         &
                       rv,ir,jr,hr,                                                   &
                       xradar,yradar,hradar,lev_num,evl,                              &
                       xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                       rdnw,rdn,mub,mu,                                               &
                       rmse_rv,spd_rv                                                 &
                      )                      
call innoation_rmse_rf(analysis_var_num,ra_data_num,numstat,                          &
                       xstat,                                                         &
                       rf,ir,jr,hr,                                                   &
                       xradar,yradar,hradar,lev_num,evl,                              &
                       xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                       rdnw,rdn,mub,mu,                                               &
                       rmse_rf,spd_rf                                                 &
                       )                      
call output_radio_r2s_rvrf(rmse_rv,spd_rv,rmse_rf,spd_rf,iradar,lab_time_num,1) 
           
endif   ! opt_rms_rad
         
enddo      !   for radar number loop  

if( enable_inflat == 1 .and. inflat_opt == 2 ) then

call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
!xstat=(1-ratio_af)*xstatf+ratio_af*xstat
call relax_inflation_scheme(xstat,xstatf,numstat)
call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 

if(enable_rckf_gi==1) then 
  call xstat_mean_seperate_cal(xstat2,xstat_mean2,numstat,ens_num)
  xstat2=(1-ratio_af2)*xstatf2+ratio_af2*xstat2
  call xstat_mean_merge_cal(xstat2,xstat_mean2,numstat,ens_num) 
endif

endif

endif   ! sim radar obs                            
                               

  allocate(ensu  (nx,ny,nz,ens_num))
  allocate(ensv  (nx,ny,nz,ens_num))
  allocate(ensw  (nx,ny,nz,ens_num))
  allocate(ensph (nx,ny,nz,ens_num))
  allocate(enst  (nx,ny,nz,ens_num))
  allocate(ensqv (nx,ny,nz,ens_num))   
  allocate(ensqr (nx,ny,nz,ens_num)) 
  allocate(ensqi (nx,ny,nz,ens_num)) 
  allocate(ensqs (nx,ny,nz,ens_num)) 
  allocate(ensqgr(nx,ny,nz,ens_num)) 
 
call dis_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,          &
               ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr  &
              )
                 
 deallocate(xstat) 
 deallocate(xstat_mean)  
 
 
if(enable_rckf_gi==1) then 

  allocate(ensu2  (nx,ny,nz,ens_num))
  allocate(ensv2  (nx,ny,nz,ens_num))
  allocate(ensw2  (nx,ny,nz,ens_num))
  allocate(ensph2 (nx,ny,nz,ens_num))
  allocate(enst2  (nx,ny,nz,ens_num))
  allocate(ensqv2 (nx,ny,nz,ens_num))   
  allocate(ensqr2 (nx,ny,nz,ens_num)) 
  allocate(ensqi2 (nx,ny,nz,ens_num)) 
  allocate(ensqs2 (nx,ny,nz,ens_num)) 
  allocate(ensqgr2(nx,ny,nz,ens_num)) 
 
  call dis_xstat(xstat2,ens_num,numstat,nx,ny,nz,analysis_var_num,                   &
                 ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2  &
                 )
                 
 deallocate(xstat2) 
 deallocate(xstat_mean2)  
  
endif

 
 
  allocate(u  (nx,ny,nz))
  allocate(v  (nx,ny,nz))
  allocate(w  (nx,ny,nz))
  allocate(ph (nx,ny,nz))
  allocate(t  (nx,ny,nz))
  allocate(qv (nx,ny,nz))  
  allocate(qr (nx,ny,nz))
  allocate(qs (nx,ny,nz))
  allocate(qi (nx,ny,nz))
  allocate(qgr(nx,ny,nz))

  allocate(umean  (nx,ny,nz))
  allocate(vmean  (nx,ny,nz))
  allocate(wmean  (nx,ny,nz))
  allocate(phmean (nx,ny,nz))
  allocate(tmean  (nx,ny,nz))
  allocate(qvmean (nx,ny,nz))     
  allocate(qrmean (nx,ny,nz))
  allocate(qimean (nx,ny,nz))
  allocate(qsmean (nx,ny,nz))
  allocate(qgrmean(nx,ny,nz))     
     
endif    
!---------------------------------------------------------------------------------------  
!  end analysis for simulated observations
!---------------------------------------------------------------------------------------
call check_hydr(nx,ny,nz,ens_num,ensqv,ensqr,ensqi,ensqs,ensqgr)


if( enable_rms_anal == 1 ) then
  ireg=1
  print*,'cal the mean '
  call mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,               &
                nx,ny,nz,ens_num,                                                       &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean      &
               )
  
  print*,'cal the rms'   
  if(opt_rms == 1) then
  call rms_cal(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &               
               nx,ny,nz,ens_num,                                                                  &
               umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
               stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
               stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d       &
              )
  else if(opt_rms == 2) then
  call rms_convetive_region(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &
                            nx,ny,nz,ens_num,                                                                  &
                            umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
                            stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                            stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                            ref_x                                                                              &
                            )
  endif         
           
  call output_rmse(stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                   stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                   lab_time_num,lab_domain,ireg,nx,ny,nz)
  
  print*,'cal the spread'
  if(opt_rms == 1) then   
  call ens_spread(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,               &
                  nx,ny,nz,ens_num,                                                       &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,     &
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr     &
                 )
  else if(opt_rms == 2) then 
  call ens_spread_convetive_region(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,               &
                                   nx,ny,nz,ens_num,                                                       &
                                   umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,     &
                                   spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                                   ref_x                                                                   & 
                                   )
  endif
                 
  call output_spd(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                  lab_time_num,lab_domain,ireg,nx,ny,nz)    
                  
  call output_ratio_r2s(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                        stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,               &
                        lab_time_num,lab_domain,ireg,nx,ny,nz)  
                        
                        
                        
                        
                        
  if(enable_rckf_gi==100) then
  
  print*,'cal the mean2 '
  call mean_cal(ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2,               &
                nx,ny,nz,ens_num,                                                                 &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean                &
               )
  
  print*,'cal the rms2'   
  if(opt_rms == 1) then
  call rms_cal(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &               
               nx,ny,nz,ens_num,                                                                  &
               umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
               stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
               stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d       &
              )
  else if(opt_rms == 2) then
  call rms_convetive_region(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,                &
                            nx,ny,nz,ens_num,                                                                  &
                            umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,                &
                            stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                            stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                            ref_x                                                                              &
                            )
  endif         
           
  call output_rmse(stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,                          &
                   stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,      &
                   lab_time_num2,lab_domain,ireg,nx,ny,nz)
  
  print*,'cal the spread2'
  if(opt_rms == 1) then   
  call ens_spread(ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2,               &
                  nx,ny,nz,ens_num,                                                                 &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,               &
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr                &
                 )
  else if(opt_rms == 2) then 
  call ens_spread_convetive_region(ensu2,ensv2,ensw2,ensph2,enst2,ensqv2,ensqr2,ensqi2,ensqs2,ensqgr2,     &
                                   nx,ny,nz,ens_num,                                                       &
                                   umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,     &
                                   spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                                   ref_x                                                                   & 
                                   )
  endif
                 
  call output_spd(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                  lab_time_num2,lab_domain,ireg,nx,ny,nz)     
  
  call output_ratio_r2s(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,     &
                        stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,               &
                        lab_time_num2,lab_domain,ireg,nx,ny,nz) 
                        
  endif                                                   
             
endif          
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                        Writing ens files
! 
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------
print*,'writing analysis field to files'   

  do i=1, ens_num     ! for writing ens file loop
  
  u  = ensu (:,:,:,i)
  v  = ensv (:,:,:,i)
  w  = ensw (:,:,:,i)
  ph = ensph(:,:,:,i)
  t  = enst (:,:,:,i)
  qv = ensqv(:,:,:,i)
  qr = ensqr(:,:,:,i)
  qi = ensqi(:,:,:,i)
  qs = ensqs(:,:,:,i)
  qgr = ensqgr(:,:,:,i)  

  file_cdfid=rec_cdfid(i)
  file_rcode=rec_rcode(i)
  call output_data(u,v,w,ph,t,qv,qr,qi,qs,qgr,nx,ny,nz,                             &
                   analysis_var_num,rec_id_var,file_cdfid,rec_dims,file_rcode       &
                   )
  
  enddo   ! for writing ens file loop
  
  do i=1 , ens_num
  call ncclos(rec_cdfid(i),rec_rcode(i))
  enddo
  
if(enable_rckf_gi==1) then
print*,'writing analysis field to files2' 
  do i=1, ens_num     ! for writing ens file loop
  
  u  = ensu2 (:,:,:,i)
  v  = ensv2 (:,:,:,i)
  w  = ensw2 (:,:,:,i)
  ph = ensph2(:,:,:,i)
  t  = enst2 (:,:,:,i)
  qv = ensqv2(:,:,:,i)
  qr = ensqr2(:,:,:,i)
  qi = ensqi2(:,:,:,i)
  qs = ensqs2(:,:,:,i)
  qgr = ensqgr2(:,:,:,i)  

  file_cdfid=rec_cdfid2(i)
  file_rcode=rec_rcode2(i)
  call output_data(u,v,w,ph,t,qv,qr,qi,qs,qgr,nx,ny,nz,                               &
                   analysis_var_num,rec_id_var2,file_cdfid,rec_dims2,file_rcode       &
                   )
  
  enddo   ! for writing ens file loop
  
  do i=1 , ens_num
  call ncclos(rec_cdfid2(i),rec_rcode2(i))
  enddo

endif 
    
  end program ensrf_sim
