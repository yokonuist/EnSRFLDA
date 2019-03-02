  program ensrf_real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            program ensrf_real
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
 ! integer,parameter     :: analysis_var_num=12
  integer,parameter     :: analysis_var_num=13
  integer,parameter     :: basic_var_num_1=15
  integer,parameter     :: basic_var_num_2=2
  integer               :: i,j,k,ii,jj,kk
  integer               :: istart(4), iend(4)
  integer               :: istartt(3), iendd(3)
  character (len=4)     :: lab_ens
  character (len=2)     :: lab_domain
  character (len=2)     :: lab_time_num
  character (len=2)     :: lab_radar_num
  
  integer               :: rec_id_var(analysis_var_num)
  integer,allocatable   :: rec_cdfid(:)
 ! integer               :: rec_dims(analysis_var_num,3)
  integer               :: rec_dims(12,3)
  integer               :: rec_dims2(2)
  integer,allocatable   :: rec_rcode(:)
  character(len=80)     :: analysis_var_name(analysis_var_num)
  integer               :: file_cdfid,file_rcode
  character(len=80)     :: basic_var_name_1(basic_var_num_1)
  character(len=80)     :: basic_var_name_2(basic_var_num_2)
    
  real,allocatable      :: u(:,:,:),v(:,:,:),w(:,:,:),ph(:,:,:),t(:,:,:),qv(:,:,:),  &
                           qr(:,:,:),qi(:,:,:),qs(:,:,:),qgr(:,:,:),smois(:,:,:),tslb(:,:,:),tsk(:,:) 
                                           
  real,allocatable      :: ensu(:,:,:,:),ensv(:,:,:,:),ensw(:,:,:,:),ensph(:,:,:,:),enst(:,:,:,:),ensqv(:,:,:,:), &
                           ensqr(:,:,:,:),ensqi(:,:,:,:),ensqs(:,:,:,:),ensqgr(:,:,:,:),enssmois(:,:,:,:),enstslb(:,:,:,:),&
                           enstsk(:,:,:)    
                           
  real,allocatable      :: umean(:,:,:),vmean(:,:,:),wmean(:,:,:),phmean(:,:,:),tmean(:,:,:),qvmean(:,:,:),  &
                           qrmean(:,:,:),qimean(:,:,:),qsmean(:,:,:),qgrmean(:,:,:),smoismean(:,:,:),tslbmean(:,:,:),&
                           tskmean(:,:)                          
  
  real,allocatable      :: xlon(:,:),xlat(:,:),ulon(:,:),ulat(:,:),vlon(:,:),vlat(:,:),phb(:,:,:),hgt(:,:),mub(:,:,:),mu(:,:,:), &
                           znw(:),znu(:),rdnw(:),rdn(:),dzs(:),zs(:)
!!!for soil depth !!!
  real                  :: p_top
  

  character(len=5)      :: radar_name
  integer               :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num  
  real,allocatable      :: rv(:),rf(:),hr_rf(:),hr_rv(:)
  integer,allocatable   :: ir(:),jr(:)
  real,allocatable      :: evl(:),azimuth(:)
  integer               :: iradar
  !!!!!!!OBS INFORMATiON!!!!!!
  real*4,allocatable      :: soil_m(:),soil_t(:),soil_st(:)
  integer               :: obs_lad_num
  !!!!!!!!!!!
  integer               :: ra_data_num

  real,allocatable     ::  usnd(:),vsnd(:),temperaturesnd(:),latsnd(:),lonsnd(:),hgtsnd(:)  
  integer              ::  obs_snd_num
  
  real,allocatable      :: xstat(:,:),xstat_mean(:),xstatf(:,:),xstat_meanf(:),xstatff(:,:) 
  integer               :: numstat 
  integer               :: istat
  real                  :: rmse_rv,spd_rv,rmse_rf,spd_rf
  
  real,allocatable      :: stdu2d(:),stdv2d(:),stdw2d(:),stdph2d(:),stdt2d(:),stdqv2d(:), &
                           stdqr2d(:),stdqi2d(:),stdqs2d(:),stdqgr2d(:),stdsmois2d(:),stdtslb2d(:),&
                           stdtsk2d(:)!!!? only one number
  real                  :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,stdtsk
  real                  :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,&
                           spd_tsk  
  integer               :: ireg,ireg_1  
  integer               :: temp1,temp2,temp3

            
  data analysis_var_name/'U','V','W','PH','T','QVAPOR','QRAIN','QICE','QSNOW','QGRAUP','SMOIS','TSLB','TSK'/
  
  data basic_var_name_1/'XLONG','XLAT','XLONG_U','XLAT_U','XLONG_V','XLAT_V','PHB','HGT','P_TOP','ZNW','ZNU',  &
                        'RDNW','RDN','DZS','ZS'/
   !!!!DZS,ZS                     
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
  write(6,nml =        obs_info) !!!! 
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
  allocate(smois (nx,ny,4))
  allocate(tslb (nx,ny,4))
  allocate(tsk (nx,ny))
  
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
  allocate(enssmois (nx,ny,4,ens_num))
  allocate(enstslb (nx,ny,4,ens_num)) 
  allocate(enstsk (nx,ny,ens_num)) 
   
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
  allocate(smoismean (nx,ny,4))
  allocate(tslbmean (nx,ny,4))      
  allocate(tskmean (nx,ny)) !

  allocate(xlon(nx,ny))
  allocate(xlat(nx,ny))
  allocate(ulon(nx,ny))
  allocate(ulat(nx,ny))
  allocate(vlon(nx,ny))
  allocate(vlat(nx,ny)) 
  allocate(phb(nx,ny,nz))
!!!for siol depth!!!
  allocate(dzs(4))
  allocate(zs(4))
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
  allocate(stdsmois2d (4))
  allocate(stdtslb2d (4))  
  allocate(stdtsk2d(1))
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading true files
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------
!stop
!!!!
!print*,'begin the test of ensrf 1'

 call read_mean_file(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                                     &
                    rec_dims,rec_dims2,file_rcode,                                                                 &
                    umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean, & 
                    xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,znw,znu,p_top,rdnw,rdn,                           &
                    basic_var_num_1,basic_var_name_1,                                                              & 
                    lab_domain,lab_time_num                                                                        &            
                   ) 

!write
!  print*,'finished reading mean fields'  
  print*,'******READING SMOIS & TSLB MEAN FIELD!*****'!,tskmean(1,1)                       
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading ens files
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------   
  
 call read_ens_file(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                    &
                   rec_dims,rec_dims2,file_rcode,rec_cdfid,rec_rcode,                            &
                   ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk, & 
                   mub,mu,basic_var_num_2,basic_var_name_2,                                      &
                   lab_domain,lab_time_num                                                       &   
                  )

!open(1000,file='test.txt')
!write(1000,*)enstsk(:,:,1)
!close(1000)

!stop
!  call check_soillayer(nx,ny,nz,ens_num,enssmois,enstslb,enstsk)               
!  print*,'finished reading ens fields'                        
  print*,'******FINISHED READING SMOIS & TSLB ENS FIELD!*****'!!,enstsk(1,1,1)

!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Processing the analysis
!
!---------------------------------------------------------------------------------------  
!--------------------------------------------------------------------------------------- 


numstat=nx*ny*nz*(analysis_var_num-3)+nx*ny*4*2+nx*ny

!write(*,*)'Allocated',numstat*4/1024,'KB space for EnSRF!'

allocate(xstat(numstat,ens_num))
allocate(xstat_mean(numstat))
allocate(xstatf(numstat,ens_num))
allocate(xstat_meanf(numstat) )
allocate(xstatff(numstat,ens_num)) ! For test

!stop
              
if( enable_rms_anal == 1 ) then

!=======================================================================
!                      CALCULATE THE RMS
!=======================================================================
  ireg=0
!  print*,'cal the mean '
!  call mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,                  &
!                nx,ny,nz,ens_num,                                                                                  &
!                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,tslbmean,smoismean,tskmean      &
!               )
  
  
  print*,'cal the spread'!,tskmean(1,1),enstsk(1,1,2)
 
  call ens_spread(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,                 &
                  nx,ny,nz,ens_num,                                                                                 &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,    &
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk     &
                 )
!print*,'begin the test of ensrf 2 ',lab_time_num,'fuck',lab_domain    
  call output_spd(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk,    &
                  lab_time_num,lab_domain,ireg,nx,ny,nz)     
                    
!print*,'begin the test of ensrf 3 '
!print*,' finish RMS analysis '      
endif 


!!! good start ,everything seems ok  2013 /5 /5 !!!
 
!stop                                     
!---------------------------------------------------------------------------------------  
! analysis for real observations
!---------------------------------------------------------------------------------------

if(enable_real == 1 ) then             !  start assimilating real observations

 !print*,'enstsk(1,1,1)',enstsk(1,1,1)
 call fill_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,                                 &
                ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk  &
               )
xstatff=xstat

if(obs_type == 0 .or. obs_type == 1) then ! real snd obs


call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
xstatf=xstat
call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 



call read_real_snd_info(obs_snd_num)

allocate(usnd(obs_snd_num))
allocate(vsnd(obs_snd_num))
allocate(temperaturesnd(obs_snd_num))
allocate(latsnd(obs_snd_num))
allocate(lonsnd(obs_snd_num))
allocate(hgtsnd(obs_snd_num))

call read_real_snd_data(usnd,vsnd,temperaturesnd,latsnd,lonsnd,hgtsnd,obs_snd_num) 

if(enable_inflat==1 .and. inflat_opt == 1 ) then
 print*,'calling inflation_cal'
 call inflation_cal(xstat,xstat_mean,numstat)
 print*,'finished inflation_cal'
endif  


print*,'the step of assimilate'

call assimilate_sounding_cal(analysis_var_num,obs_snd_num,numstat,                          &
                             xstat,xstat_mean,                                              &
                             usnd,vsnd,temperaturesnd,latsnd,lonsnd,hgtsnd,                 &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                             rdnw,rdn,mub,mu                                                &
                            ) 

print*,'finish the step of assimilate'


if( enable_inflat == 1 .and. inflat_opt == 2 ) then
call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
call relax_inflation_scheme(xstat,xstatf,numstat)
call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 
endif


 
 call dis_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,          &
               ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr  &
              )
                 
 !call dis_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,                                &
 !              ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk &
 !             )

 deallocate(xstat) 
 deallocate(xstat_mean)  
 

endif   ! real snd obs


if(obs_type == 0 .or. obs_type == 2) then ! real radar obs
!---------------------------------------------------------------------------------------  
!                       Reading observations
!---------------------------------------------------------------------------------------  

call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
xstatf=xstat
call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 

do iradar=1,radar_num    !   for radar number loop  



call read_real_radar_info(iradar,xradar,yradar,hradar,lev_num,ra_data_num)

if(allocated(rv)) deallocate(rv)
allocate(rv(ra_data_num))
if(allocated(rf)) deallocate(rf)
allocate(rf(ra_data_num))
if(allocated(ir)) deallocate(ir)
allocate(ir(ra_data_num))
if(allocated(jr)) deallocate(jr)
allocate(jr(ra_data_num))
if(allocated(hr_rv)) deallocate(hr_rv)
allocate(hr_rv(ra_data_num))
if(allocated(hr_rf)) deallocate(hr_rf)
allocate(hr_rf(ra_data_num))
if(allocated(evl)) deallocate(evl)
allocate(evl(ra_data_num))
if(allocated(azimuth)) deallocate(azimuth)
allocate(azimuth(ra_data_num))

call read_real_radar_data(ra_data_num,rv,rf,ir,jr,hr_rv,hr_rf,evl,azimuth)




if(opt_rms_rad   == 1) then
      
      ireg_1=0
      
      call innoation_rmse_rv(analysis_var_num,ra_data_num,numstat,                          &
                             xstat,                                                         &
                             rv,ir,jr,hr_rv,evl,azimuth,                                    &
                             xradar,yradar,hradar,lev_num,                                  &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                             rdnw,rdn,mub,mu,iradar,ireg_1,                                   &
                             rmse_rv,spd_rv                                                 &
                             )
                       
      call innoation_rmse_rf(analysis_var_num,ra_data_num,numstat,                          &
                             xstat,                                                         &
                             rf,ir,jr,hr_rf,evl,azimuth,                                    &
                             xradar,yradar,hradar,lev_num,                                  &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                             rdnw,rdn,mub,mu,iradar,ireg_1,                                   &
                             rmse_rf,spd_rf                                                 &
                             )
      call output_radio_r2s_rvrf(rmse_rv,spd_rv,rmse_rf,spd_rf,iradar,lab_time_num,0) 
                                   
endif


if(enable_inflat==1 .and. inflat_opt == 1 ) then
 print*,'calling inflation_cal'
 call inflation_cal(xstat,xstat_mean,numstat)
 print*,'finished inflation_cal'
endif  

call assimilate_radar_cal(analysis_var_num,ra_data_num,numstat,                          &
                          xstat,xstat_mean,                                              &
                          rf,rv,ir,jr,hr_rf,hr_rv,evl,azimuth,                           &
                          xradar,yradar,hradar,lev_num,                                  &
                          xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                          rdnw,rdn,mub,mu                                                &
                         ) 
                         

if(opt_rms_rad   == 1) then

      ireg_1=1

      call innoation_rmse_rv(analysis_var_num,ra_data_num,numstat,                          &
                             xstat,                                                         &
                             rv,ir,jr,hr_rv,evl,azimuth,                                    &
                             xradar,yradar,hradar,lev_num,                                  &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                             rdnw,rdn,mub,mu,iradar,ireg_1,                                   &
                             rmse_rv,spd_rv                                                 &
                             )
                             
      call innoation_rmse_rf(analysis_var_num,ra_data_num,numstat,                          &
                             xstat,                                                         &
                             rf,ir,jr,hr_rf,evl,azimuth,                                    &
                             xradar,yradar,hradar,lev_num,                                  &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,           &
                             rdnw,rdn,mub,mu,iradar,ireg_1,                                   &
                             rmse_rf,spd_rf                                                 &
                             )
      call output_radio_r2s_rvrf(rmse_rv,spd_rv,rmse_rf,spd_rf,iradar,lab_time_num,1)  
                                   
endif   ! opt_rms_rad
         
enddo   !   for radar number loop  

if( enable_inflat == 1 .and. inflat_opt == 2 ) then
call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
call relax_inflation_scheme(xstat,xstatf,numstat)
call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 
endif

                          
                               

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
     
endif ! real radar obs   
!---------------------------------------------------------------------------------------  
!  end analysis for simulated observations
!---------------------------------------------------------------------------------------
call check_hydr(nx,ny,nz,ens_num,ensqv,ensqr,ensqi,ensqs,ensqgr)



if(obs_type == 0 .or. obs_type == 3) then ! real land obs

 write(*,*)'Process land data assimilation,obs_type is',obs_type
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num) ! XSTAT --> pert, XSTAT_MEAN --> mean field
 xstatf=xstat  ! xstaff ---> pert field

 call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 

!***************** test code is ok! *************************************
! print*,'numstat,xstat(3920001,1),xstat_mean(3920001)',numstat,xstat(3920001,1),xstat_mean(3920001)
! stop
! everything is ok ! 2013/5/11 ok
!******************* fine ! **************************************************** 
 call read_real_land_info(obs_lad_num)

!obs_lad_num=10
!if(allocated(usnd)) deallocate(usnd)
!allocate(usnd(obs_snd_num))

if(allocated(soil_t)) deallocate(soil_t)
allocate(soil_t(obs_lad_num))
if(allocated(soil_st)) deallocate(soil_st)
allocate(soil_st(obs_lad_num))
if(allocated(soil_m)) deallocate(soil_m)
allocate(soil_m(obs_lad_num))
if(allocated(latsnd)) deallocate(latsnd)
allocate(latsnd(obs_lad_num))
if(allocated(lonsnd)) deallocate(lonsnd)
allocate(lonsnd(obs_lad_num))
if(allocated(hgtsnd)) deallocate(hgtsnd)
allocate(hgtsnd(obs_lad_num))

 !print*,numstat,obs_lad_num 
 !print*,''

  call read_real_land_data(soil_m,soil_t,soil_st,latsnd,lonsnd,hgtsnd,obs_lad_num)
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print*,obs_lad_num 
! something wrong in reading obs, tslb?!
! print*,'usnd,soil_m,soil_t,latsnd,lonsnd,hgtsnd,obs_lad_num',usnd(1),soil_m(1),soil_t(1),latsnd(1),lonsnd(1),hgtsnd(1),obs_lad_num
! open(1000,file='/home/guoyk/test.dat')
! write(1000,'(6f11.2)')usnd,soil_m,soil_t,latsnd,lonsnd,hgtsnd
! print*,'xstat(3976001,1),enstslb(1,1,1,1),soil_t(1)',xstat(3976001,1),enstslb(1,1,1,1),soil_t(1)
! print*,'xstat(4032001,1),enstsk(1,1,1),soil_st(1)',xstat(4032001,1),enstsk(1,1,1),soil_st(1)
! stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!******************* inflation method 1 *********************************!
! before asimilate: use a number '(1+coef_inflat)' to inflate the pertub !
!**********************************************************************************!!

   if(enable_inflat==1 .and. inflat_opt == 1 ) then

  print*,'calling inflation_cal'

 call inflation_cal(xstat,xstat_mean,numstat)

  print*,'finished inflation_cal'

   endif   !!!if(enable_inflat==1 .and. inflat_opt == 1 )  

  print*,'begin the step of assimilate land obs'!,xlon(1,1)

!stop

 call assimilate_land_cal(analysis_var_num,obs_lad_num,numstat,                          &
                             xstat,xstat_mean,                                           &
                             soil_m,soil_t,soil_st,latsnd,lonsnd,hgtsnd,                 &
                             xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,hgt,znw,znu,p_top,    &
                             rdnw,rdn,mub,mu                                             &
                            ) 

  print*,'finish the step of assimilating soil data! '

!stop

   if( enable_inflat == 1 .and. inflat_opt == 2 ) then
 call xstat_mean_seperate_cal(xstat,xstat_mean,numstat,ens_num)
 call relax_inflation_scheme(xstat,xstatf,numstat)
 call xstat_mean_merge_cal(xstat,xstat_mean,numstat,ens_num) 
  endif  !! if( enable_inflat == 1 .and. inflat_opt == 2 )


! ******check every point of soil data ********************
!********************guoyk 2013/7/2
!******************************************************************

if(1 == 1) then
 call check_xstat_mean__cal(xstat,xstat_mean,numstat,ens_num,xstatf,xstatff)
! print*,'******finish check*******'
endif
deallocate(xstatff)

if(allocated(enstsk))deallocate(enstsk)
allocate(enstsk (nx,ny,ens_num))

 call dis_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,                                &
               ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk &
              )
!stop                 
 deallocate(xstat) 
 deallocate(xstat_mean) 


endif  ! land obs
!-------------------------------------------------------------------------------
! END ASSIMILATE OBS OF SMOIS/ TSLB /TSK
!-------------------------------------------------------------------------------

 call check_soillayer(nx,ny,nz,ens_num,enssmois,enstslb,enstsk)

endif ! if(enable_real == 1)
 
!stop

if( enable_rms_anal == 1 ) then
  ireg=1
  print*,'cal the mean '
  call mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,                &
                nx,ny,nz,ens_num,                                                                                &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean    &
               )
  
  
  print*,'cal the spread'

  call ens_spread(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,                 &
                  nx,ny,nz,ens_num,                                                                                 &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,    &
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk     &
                 )

                 
  call output_spd(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk, &
                  lab_time_num,lab_domain,ireg,nx,ny,nz)    
                     
             
endif  

!stop        
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
  smois = enssmois(:,:,:,i)
  tslb = enstslb(:,:,:,i)  
  tsk = enstsk(:,:,i)

  file_cdfid=rec_cdfid(i)
  file_rcode=rec_rcode(i)
  call output_data(u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz,                     &
                   analysis_var_num,rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode    &
                   )
  
  enddo   ! for writing ens file loop
  
  do i=1 , ens_num
  call ncclos(rec_cdfid(i),rec_rcode(i))
  enddo

 !stop
   
  end program ensrf_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!        obs info
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_real_land_info(obs_lad_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_real_land_info
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.03.27   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

include 'namelist.inc'

character(len=3)      :: lab_time_num
integer               :: obs_lad_num
integer               :: temp1,temp2,temp3

temp1=int(time_num/100     )-int(time_num/1000    )*10  
temp2=int(time_num/10      )-int(time_num/100     )*10  
temp3=int(time_num/1       )-int(time_num/10      )*10
lab_time_num = char(48+temp1)//char(48+temp2)//char(48+temp3) 
   
length_obs_dir=len_trim(obs_file_dir)
 
open(1083,file=''//obs_file_dir(1:length_obs_dir)//'/lad/obs_'//lab_time_num//'.dat')
read(1083,'(i11)')obs_lad_num                             
                                                                         
end subroutine read_real_land_info       

subroutine read_real_land_data(soil_m,soil_t,soil_st,latsnd,lonsnd,hgtsnd,obs_lad_num) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_real_land_data
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
include 'namelist.inc'
integer  :: obs_lad_num
real     :: usnd(obs_lad_num),soil_m(obs_lad_num),soil_t(obs_lad_num),latsnd(obs_lad_num), &
            lonsnd(obs_lad_num),hgtsnd(obs_lad_num),soil_st(obs_lad_num)
integer  :: tem,tem0,i   

tem=0

do i=1,obs_lad_num

  read(1083,'(6f11.6)')soil_m(i),soil_t(i),soil_st(i),latsnd(i),lonsnd(i),hgtsnd(i)

enddo
! print*,soil_m(1),soil_t(1),soil_st(1),latsnd(1),lonsnd(1),hgtsnd(1)
 print*,'finished reading all real lda obs'

close(1083)

end subroutine read_real_land_data
