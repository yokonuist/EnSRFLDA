program anal_result

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          Program anal_result
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: analysis the result files
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
  include 'netcdf.inc'
  include 'namelist.inc'
   
  character (len=500)   :: input_file_name
  integer               :: length_input,length_file_dir,length_ens_head      
  integer,parameter     :: analysis_var_num=13
  integer               :: i,j,k,ii,jj,kk
  integer               :: istart(4), iend(4),istartt(3),iendd(3)
  character (len=4)     :: lab_ens
  character (len=2)     :: lab_domain
  character (len=2)     :: lab_time_num
  
  integer               :: rec_id_var(analysis_var_num)
  integer,allocatable   :: rec_cdfid(:)
  integer               :: rec_dims(12,3)
  integer               :: rec_dims2(2)
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
  real,allocatable      :: utrue(:,:,:),vtrue(:,:,:),wtrue(:,:,:),phtrue(:,:,:),ttrue(:,:,:),qvtrue(:,:,:),  &
                           qrtrue(:,:,:),qitrue(:,:,:),qstrue(:,:,:),qgrtrue(:,:,:),smoistrue(:,:,:),tslbtrue(:,:,:),&
                           tsktrue(:,:)  
                           
  real,allocatable      :: upert(:,:,:,:),vpert(:,:,:,:),wpert(:,:,:,:),phpert(:,:,:,:),tpert(:,:,:,:),qvpert(:,:,:,:), &
                           qrpert(:,:,:,:),qipert(:,:,:,:),qspert(:,:,:,:),qgrpert(:,:,:,:),smoispert(:,:,:,:),tslbpert(:,:,:,:),&
                           tskpert(:,:,:)                                                     
                           
  real,allocatable      :: stdu2d(:),stdv2d(:),stdw2d(:),stdph2d(:),stdt2d(:),stdqv2d(:), &
                           stdqr2d(:),stdqi2d(:),stdqs2d(:),stdqgr2d(:),stdsmois2d(:),stdtslb2d(:) !(nz)
                           
  real                  :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,stdtsk,stdtsk2d ! Alarm
  real                  :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk
  integer               :: temp1,temp2,temp3

!---------------------------------------------------------------------------------------------------------    
!------------------for eignvalue decomposition variables-------------------------------------------------  
!---------------------------------------------------------------------------------------------------------  
  integer               :: ks                !  for eignvalue decomposition option,ks=-1: self; ks=0: depature; ks=1: normalized depature
  integer               :: mnl               !  for eignvalue decomposition mnl=min(m,n),m: number of state variables ,n: number of ensemble
  integer               :: numstat           !  number of state variables
  real,allocatable      :: xstat(:,:),xtp(:,:)
  real,allocatable      :: er(:,:),egvt_stat(:,:),egvt_ens(:,:)
!---------------------------------------------------------------------------------------------------------
  
  real,allocatable      :: normu(:,:,:,:),normv(:,:,:,:),normw(:,:,:,:),normph(:,:,:,:),normt(:,:,:,:),normqv(:,:,:,:), &
                           normqr(:,:,:,:),normqi(:,:,:,:),normqs(:,:,:,:),normqgr(:,:,:,:),normsmois(:,:,:,:),normtslb(:,:,:,:),&
                           normtsk(:,:,:)

            
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
  read  (5,nml =        obs_info)
  write (6,nml =        obs_info)   
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
  allocate(tsk(nx,ny))

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
  allocate(tskmean(nx,ny))

  allocate(upert  (nx,ny,nz,ens_num))
  allocate(vpert  (nx,ny,nz,ens_num))
  allocate(wpert  (nx,ny,nz,ens_num))
  allocate(phpert (nx,ny,nz,ens_num))
  allocate(tpert  (nx,ny,nz,ens_num))
  allocate(qvpert (nx,ny,nz,ens_num))     
  allocate(qrpert (nx,ny,nz,ens_num))
  allocate(qipert (nx,ny,nz,ens_num))
  allocate(qspert (nx,ny,nz,ens_num)) 
  allocate(qgrpert(nx,ny,nz,ens_num))
  allocate(smoispert (nx,ny,4,ens_num))
  allocate(tslbpert (nx,ny,4,ens_num))
  allocate(tskpert(nx,ny,ens_num)) 
  
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
  allocate(smoistrue (nx,ny,4))
  allocate(tslbtrue (nx,ny,4))       
  allocate(tsktrue  (nx,ny))  
  
  allocate(normu  (nx,ny,nz,ens_num))
  allocate(normv  (nx,ny,nz,ens_num))
  allocate(normw  (nx,ny,nz,ens_num))
  allocate(normph (nx,ny,nz,ens_num))
  allocate(normt  (nx,ny,nz,ens_num))
  allocate(normqv (nx,ny,nz,ens_num)) 
  allocate(normqr (nx,ny,nz,ens_num)) 
  allocate(normqi (nx,ny,nz,ens_num)) 
  allocate(normqs (nx,ny,nz,ens_num))  
  allocate(normqgr(nx,ny,nz,ens_num))
  allocate(normsmois (nx,ny,4,ens_num))
  allocate(normtslb (nx,ny,4,ens_num))
  allocate(normtsk  (nx,ny,ens_num))      
    
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

!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading true files
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------   
if(1==1) then  
  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)
  
  temp1=int(time_num/10      )-int(time_num/100     )*10  
  temp2=int(time_num/1       )-int(time_num/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)  
  

  length_file_dir =len_trim(obs_file_dir)
  length_ens_head=len_trim(ens_file_head)
  
  input_file_name=''//obs_file_dir(1:length_file_dir)//'/'//ens_file_head(1:length_ens_head)//'_d'//lab_domain//'_'//lab_time_num//''
  
  print*,"INPUT TRUE FILE IS: ",trim(input_file_name)
  length_input = len_trim(input_file_name)
  
  
! Now read the file
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                &
                          rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode,                             &
                          utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,smoistrue,tslbtrue,tsktrue,nx,ny,nz     &
                          ) 
  call ncclos(file_cdfid,file_rcode) 
endif


!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading mean files
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
  
  print*,"INPUT MEAN FILE IS: ",trim(input_file_name)
  length_input = len_trim(input_file_name)
  
  
! Now read the file
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                &
                          rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode,                            &
                          umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,nx,ny,nz     &
                          )                                              
                                                 
  call ncclos(file_cdfid,file_rcode)                              
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
  enstsk(:,:,i)= tsk
  
  enddo   ! for ens file loop
 
 do i=1 , ens_num
 call ncclos(rec_cdfid(i),rec_rcode(i))
 enddo
   
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Processing the analysis
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  

if(1==1) then   
  print*,'cal the mean '
  call mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,                  &
                nx,ny,nz,ens_num,                                                                                  &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean      &
               )
               
  call pert_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,                   &
                upert,vpert,wpert,phpert,tpert,qvpert,qrpert,qipert,qspert,qgrpert,smoispert,tslbpert,tskpert,      &
                nx,ny,nz,ens_num,                                                                                   &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean       &
                )
                  
 ! print*,'cal the max and min'
 ! call  maxmin(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,                &
 !              nx,ny,nz,ens_num,                                                        &
 !              umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean       &
 !             )  
  print*,'cal the rms'                            

  call rms_cal(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,smoistrue,tslbtrue,tsktrue,              &
               nx,ny,nz,ens_num,                                                                                           &
               umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,              &
               stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,stdtsk,                           &
               stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,stdsmois2d,stdtslb2d,stdtsk2d  &
              )              
  print*,'writing the mean and rms'
  open(8553,file='rms_d'//lab_domain//'_'//lab_time_num//'.dat')
  write(8553,*) 'rms of the ensemble'
  write(8553,*) 'stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,stdtsk'

  write(8553,'(10F20.8)') stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,stdtsk

  write(8553,*) 'vertical rms pattern of the ensemble'
  write(8553,*) 'stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr'
!!!!
    do k=nz,1,-1
      write(8553,'(10F20.8)') stdu2d(k),stdv2d(k),stdw2d(k),stdph2d(k),stdt2d(k),stdqv2d(k),stdqr2d(k),stdqi2d(k),stdqs2d(k),stdqgr2d(k)
    enddo

  write(8553,*) 'vertical rms(lda) pattern of the ensemble '
  write(8553,*) 'stdsmois,stdtslb'
     do k=4,1,-1
     write(8553,'(2F20.8)')stdsmois2d(k),stdtslb2d(k)
    enddo
   write(8553,*) 'stdtsk'
   write(8553,'(F20.8)')stdtsk2d  
  close(8553)
!!!!  
print*,'cal the spread'
  call ens_spread(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,                  &
                  nx,ny,nz,ens_num,                                                                                  &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,     &
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk      &
                 )
                 
  open(8553,file='./ens_spread_d'//lab_domain//'_'//lab_time_num//'.dat')        
  write(8553,*) 'ensemble spread:'
  write(8553,*) 'spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk'    
  write(8553,'(10F20.8)') spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk
  close(8553)     
endif                 

!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Output perturbation field
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
if(0==1)then

  print*,'writing ensemble files for post analysis such as Grads etc.'  
  
  open(8553,file='init_pert_ens_d'//lab_domain//'_'//lab_time_num//'.dat',form='binary')
  
  do ii=1,ens_num
  
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) ensu(i,j,k,ii)-umean(i,j,k)
          enddo
       enddo
    enddo
    
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) ensv(i,j,k,ii)-vmean(i,j,k)
          enddo
       enddo
    enddo    
  
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) ensw(i,j,k,ii)-wmean(i,j,k)
          enddo
       enddo
    enddo

    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) ensph(i,j,k,ii)-phmean(i,j,k)
          enddo
       enddo
    enddo      

    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) enst(i,j,k,ii)-tmean(i,j,k)
          enddo
       enddo
    enddo
    
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) ensqv(i,j,k,ii)-qvmean(i,j,k)
          enddo
       enddo
    enddo  

   do k=1,4
       do j=1,ny
          do i=1,nx
              write(8553) enssmois(i,j,k,ii)-smoismean(i,j,k)
          enddo
       enddo
    enddo 

  do k=1,4
       do j=1,ny
          do i=1,nx
              write(8553) enstslb(i,j,k,ii)-tslbmean(i,j,k)
          enddo
       enddo
    enddo   
  print*,'finished writing ens',ii
  enddo
      
  close(8553)

endif
         
end program anal_result

