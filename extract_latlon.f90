  program extract_latlon


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
  
  real,allocatable      :: xlon(:,:),xlat(:,:),ulon(:,:),ulat(:,:),vlon(:,:),vlat(:,:),phb(:,:,:),hgt(:,:),mub(:,:,:),mu(:,:,:), &
                           znw(:),znu(:),rdnw(:),rdn(:)
  real                  :: p_top
  

  character(len=5)      :: radar_name
  integer               :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num  
  real,allocatable      :: rv(:),rf(:),hr_rf(:),hr_rv(:)
  integer,allocatable   :: ir(:),jr(:)
  real,allocatable      :: evl(:),azimuth(:)
  integer               :: iradar
  
  integer               :: ra_data_num
  
  real,allocatable      :: xstat(:,:),xstat_mean(:),xstatf(:,:),xstat_meanf(:) 
  integer               :: numstat 
  integer               :: istat
  real                  :: rmse_rv,spd_rv,rmse_rf,spd_rf
  
  real,allocatable      :: stdu2d(:),stdv2d(:),stdw2d(:),stdph2d(:),stdt2d(:),stdqv2d(:), &
                           stdqr2d(:),stdqi2d(:),stdqs2d(:),stdqgr2d(:)
  real                  :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr
  real                  :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr  
  integer               :: ireg  
  integer               :: temp1,temp2,temp3

            
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
  
  allocate(xlon(nx,ny))
  allocate(xlat(nx,ny))
  allocate(ulon(nx,ny))
  allocate(ulat(nx,ny))
  allocate(vlon(nx,ny))
  allocate(vlat(nx,ny)) 
  allocate(phb(nx,ny,nz))
  allocate(hgt(nx,ny)) 
  allocate(znw   (nz)) 
  allocate(znu   (nz)) 
  allocate(rdnw  (nz)) 
  allocate(rdn   (nz))  


  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)
  
  temp1=int(time_num/10      )-int(time_num/100     )*10  
  temp2=int(time_num/1       )-int(time_num/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)  
  

  length_file_dir =len_trim(input_file_dir)
  length_ens_head=len_trim(ens_file_head)
  
  input_file_name=''//input_file_dir(1:length_file_dir)//'/'//ens_file_head(1:length_ens_head)//'_d'//lab_domain//''
  
  print*,"INPUT FILE IS: ",trim(input_file_name)
  length_input = len_trim(input_file_name)
  
  
! Now read the file
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                &
                          rec_id_var,file_cdfid,rec_dims,file_rcode,                                      &
                          umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,nx,ny,nz     &
                          ) 
                          
  call get_basic_info_1(basic_var_num_1,basic_var_name_1,                    &
                        file_cdfid,file_rcode,nx,ny,nz,                      &
                        xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,               &
                        znw,znu,p_top,rdnw,rdn                               &
                       )  
                                              
                                                 
  call ncclos(file_cdfid,file_rcode)  
  
  print*,'writing lat and lon'                 
  open(1001,file='latlon.dat')
  do i=1,nx
  do j=1,ny
  write(1001,*) xlat(i,j),xlon(i,j)
  enddo
  enddo                 
  close(1001)            
                   
                   
                   
end program extract_latlon    