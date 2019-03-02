program mean_io

!**********************************************
!                NUIST WRF-ENSRF
!  This program is used for reading and generating 
!          a mean file for a certain time
!       only analysis field will be averaged
!**********************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            program mean_io
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  include 'netcdf.inc'
  include 'namelist.inc'
   
  character (len=500)   :: input_file_name
  integer               :: length_input,length_file_dir,length_ens_head      
  !integer,parameter     :: analysis_var_num=12
  integer,parameter     :: analysis_var_num=13
  integer               :: i,j,k,ii,jj,kk
  integer               :: istart(4), iend(4)
  integer               :: istartt(3), iendd(3) !!
  character (len=4)     :: lab_ens
  character (len=2)     :: lab_domain
  character (len=2)     :: lab_time_num
  
  integer               :: rec_id_var(analysis_var_num)
  integer,allocatable   :: rec_cdfid(:)
  integer               :: rec_dims(12,3) ! 3_d
  integer               :: rec_dims2(2)   ! 2_d TSK
  integer,allocatable   :: rec_rcode(:)
  character(len=80)     :: analysis_var_name(analysis_var_num)
  integer               :: file_cdfid,file_rcode

  real,allocatable      :: u(:,:,:),v(:,:,:),w(:,:,:),ph(:,:,:),t(:,:,:),qv(:,:,:), &
                           qr(:,:,:),qi(:,:,:),qs(:,:,:),qgr(:,:,:),smois(:,:,:),tslb(:,:,:),tsk(:,:)   
  real,allocatable      :: ensu(:,:,:,:),ensv(:,:,:,:),ensw(:,:,:,:),ensph(:,:,:,:),enst(:,:,:,:),ensqv(:,:,:,:), &
                           ensqr(:,:,:,:),ensqi(:,:,:,:),ensqs(:,:,:,:),ensqgr(:,:,:,:),enssmois(:,:,:,:),enstslb(:,:,:,:),&
                           enstsk(:,:,:)  
  real,allocatable      :: umean(:,:,:),vmean(:,:,:),wmean(:,:,:),phmean(:,:,:),tmean(:,:,:),qvmean(:,:,:), &
                           qrmean(:,:,:),qimean(:,:,:),qsmean(:,:,:),qgrmean(:,:,:),smoismean(:,:,:),tslbmean(:,:,:),&
                           tskmean(:,:)        

  integer               :: temp1,temp2,temp3

            
  data analysis_var_name/'U','V','W','PH','T','QVAPOR','QRAIN','QICE','QSNOW','QGRAUP','SMOIS','TSLB','TSK'/
  
  open(unit=5, file="ensrf.input", form="formatted", status="old")
  
  read  (5,nml = input_file_info)
  write (6,nml = input_file_info)
  read  (5,nml =  dimension_info)
  write (6,nml =  dimension_info)
  read  (5,nml =   ensemble_info)
  write (6,nml =   ensemble_info)
  read  (5,nml =   analysis_info)
  write (6,nml =   analysis_info)  
  
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
  allocate(tsk (nx,ny))  !

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
  allocate(enstsk (nx,ny,ens_num)) !
    
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
  allocate(tskmean(nx,ny)) !
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading ens files
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------         
  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)
  
  temp1=int(time_num/10      )-int(time_num/100     )*10  
  temp2=int(time_num/1       )-int(time_num/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)  
    
  do i=1, ens_num
  
  temp1=int(i/100     )-int(i/1000    )*10
  temp2=int(i/10      )-int(i/100     )*10
  temp3=int(i/1       )-int(i/10      )*10
  
  lab_ens=char(48)//char(48+temp1)//char(48+temp2)//char(48+temp3)
  

  length_file_dir =len_trim(input_file_dir)
  length_ens_head=len_trim(ens_file_head)
  
  input_file_name=''//input_file_dir(1:length_file_dir)//'/'//ens_file_head(1:length_ens_head)//'_d'//lab_domain//'_'//lab_time_num//'_'//lab_ens//''
  
  print*,"INPUT FILE IS: ",trim(input_file_name)
  length_input = len_trim(input_file_name)
  
  
  
! Now read the file
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                &
                          rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode,                            &
                          u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz                              &
                          )
  rec_cdfid(i)=file_cdfid
  rec_rcode(i)=file_rcode

  ensu  (:,:,:,i) = u
  ensv  (:,:,:,i) = v
  ensw  (:,:,:,i) = w
  ensph (:,:,:,i) = ph
  enst  (:,:,:,i) = t
  ensqv (:,:,:,i) = qv
  ensqr (:,:,:,i) = qr
  ensqi (:,:,:,i) = qi
  ensqs (:,:,:,i) = qs
  ensqgr(:,:,:,i) = qgr
  enssmois (:,:,:,i) = smois
  enstslb (:,:,:,i) = tslb
  enstsk(:,:,i)   = tsk
  
  enddo   ! for ens file loop
    
!! call check_soillayer(nx,ny,nz,ens_num,enssmois,enstslb,enstsk)

  print*,'cal the mean '
  call mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,             &
                nx,ny,nz,ens_num,                                                                             &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean &
               )
               

  do i=1 , ens_num
  call ncclos(rec_cdfid(i),rec_rcode(i))
  enddo
 
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Writing mean file
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------   
  
  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)
  
  temp1=int(time_num/10      )-int(time_num/100     )*10  
  temp2=int(time_num/1       )-int(time_num/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)  
  

  length_file_dir =len_trim(input_file_dir)
  length_ens_head=len_trim(ens_file_head)
  
  input_file_name=''//input_file_dir(1:length_file_dir)//'/'//ens_file_head(1:length_ens_head)//'_d'//lab_domain//'_'//lab_time_num//'_mean'
  
  print*,"INPUT FILE IS: ",trim(input_file_name)
  length_input = len_trim(input_file_name)
    
! Now write mean file

  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                           &
                          rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode,                                       &
                          u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz                                         &
                          ) 
                          
  call output_data(umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,nx,ny,nz,   &
                   analysis_var_num,rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode                                        &
                   )
                   
  call ncclos(file_cdfid,file_rcode)  
 ! print*,'*****ADD SMOIS & TSLB IN THE MEAN FILE*****'
 ! print*,'********* TSK *****************************'
  end program mean_io     
