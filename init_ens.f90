  program init_ens
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            program init_ens
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
 !        reading the initial data file and generate the ensemble
 !
 !---------------------------------------------------------------- 
  

  implicit none
  include 'netcdf.inc'
  include 'namelist.inc'
   
  character (len=500)   :: input_file_name
  integer               :: length_input,length_file_dir,length_ens_head      
!  integer,parameter     :: analysis_var_num=10
!  integer,parameter     :: analysis_var_num=12
  integer,parameter     :: analysis_var_num=13
!  print*,analysis_var_num,'var num'
  integer,parameter     :: basic_var_num_1=15 
  integer,parameter     :: basic_var_num_2=2   
  integer               :: i,j,k,ii,jj,kk,nt,kt
  integer               :: istart(4), iend(4)
  character (len=4)     :: lab_ens
  character (len=2)     :: lab_domain
  character (len=2)     :: lab_time_num
  
  integer               :: rec_id_var(analysis_var_num)
  integer,allocatable   :: rec_cdfid(:)
  integer               :: rec_dims(analysis_var_num,3),rec_dims2(2)
  integer,allocatable   :: rec_rcode(:),rec_kk(:)
  character(len=80)     :: analysis_var_name(analysis_var_num)
  character(len=80)     :: basic_var_name_1(basic_var_num_1)  
  character(len=80)     :: basic_var_name_2(basic_var_num_2)  
  integer               :: file_cdfid,file_rcode  
  real,allocatable      :: u(:,:,:),v(:,:,:),w(:,:,:),ph(:,:,:),t(:,:,:),qv(:,:,:),qr(:,:,:),qi(:,:,:),qs(:,:,:),qgr(:,:,:)
 ! print*,' var of smois & tslb' 
  real,allocatable      :: smois(:,:,:),tslb(:,:,:),tsk(:,:) ! soil moisture, soil temperature, skin temperature
  real,allocatable      :: ensu(:,:,:,:),ensv(:,:,:,:),ensw(:,:,:,:),ensph(:,:,:,:),enst(:,:,:,:),ensqv(:,:,:,:), &
                           ensqr(:,:,:,:),ensqi(:,:,:,:),ensqs(:,:,:,:),ensqgr(:,:,:,:),enssmois(:,:,:,:),enstslb(:,:,:,:)
  real,allocatable      :: enstsk(:,:,:) 
! print*, 'add ens varname'   
  real,allocatable      :: umean(:,:,:),vmean(:,:,:),wmean(:,:,:),phmean(:,:,:),tmean(:,:,:),qvmean(:,:,:),  &
                           qrmean(:,:,:),qimean(:,:,:),qsmean(:,:,:),qgrmean(:,:,:),smoismean(:,:,:),tslbmean(:,:,:) 
  real,allocatable      :: tskmean(:,:)  !
 !print*, 'add mean varname'
  
  real,allocatable      :: upert(:,:,:,:),vpert(:,:,:,:),wpert(:,:,:,:),phpert(:,:,:,:),tpert(:,:,:,:),qvpert(:,:,:,:),&
                           smoispert(:,:,:,:),tslbpert(:,:,:,:),tskpert(:,:,:) !
 ! print*, 'add pertvar name'
                           
  real,allocatable      :: xlon(:,:),xlat(:,:),ulon(:,:),ulat(:,:),vlon(:,:),vlat(:,:),phb(:,:,:),hgt(:,:),mub(:,:,:),mu(:,:,:), &
                           znw(:),znu(:),rdnw(:),rdn(:),dzs(:),zs(:)  
                                                                             
  real                  :: p_top
                                                         
  integer               :: seed_ens
  real,allocatable      :: ranpert(:,:,:,:),ranpert2(:,:,:,:),ranpertsoil(:,:,:,:),ranpertsoil2(:,:,:,:)
  integer               :: temp1,temp2,temp3,tempkk
  real                  :: radio_to_1,radio2_to_1,spd_ens_ranpert,spd_ens_ranpert2
  real,allocatable      :: spd_ranpert(:),spd_ranpert2(:)
  real                  :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk  
  
  real                  :: ee,alpha,alpha_z 
  integer               :: scalar
  real                  :: hradius_temp
  
  real,allocatable      :: ranpert_mean(:,:,:),ranpert2_mean(:,:,:),ranpertsoil_mean(:,:,:),ranpertsoil2_mean(:,:,:)
  real,allocatable      :: ranpertsoil_2d(:,:,:),ranpertsoil_mean2d(:,:),ranpertsoil2_mean2d(:,:),ranpertsoil2_2d(:,:,:)
            
  data analysis_var_name/'U','V','W','PH','T','QVAPOR','QRAIN','QICE','QSNOW','QGRAUP','SMOIS','TSLB','TSK'/
 ! print*,'***ADD  var******similar as var of PHB or PH!***',analysis_var_name(11),analysis_var_name(12),analysis_var_name(13)
  data basic_var_name_1/'XLONG','XLAT','XLONG_U','XLAT_U','XLONG_V','XLAT_V','PHB','HGT','P_TOP','ZNW','ZNU',  &
                        'RDNW','RDN','DZS','ZS'/
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

  allocate(u(nx,ny,nz))
  allocate(v(nx,ny,nz))
  allocate(w(nx,ny,nz))
  allocate(ph(nx,ny,nz))
  allocate(t(nx,ny,nz))
  allocate(qv(nx,ny,nz))  
  allocate(qr(nx,ny,nz))
  allocate(qi(nx,ny,nz))
  allocate(qs(nx,ny,nz))
  allocate(qgr(nx,ny,nz))
  allocate(smois(nx,ny,4))
  allocate(tslb(nx,ny,4))
  allocate(tsk(nx,ny)) !
!print*,'****OK!****'

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

  allocate(upert  (nx,ny,nz,ens_num))
  allocate(vpert  (nx,ny,nz,ens_num))
  allocate(wpert  (nx,ny,nz,ens_num))
  allocate(phpert (nx,ny,nz,ens_num))
  allocate(tpert  (nx,ny,nz,ens_num))
  allocate(qvpert (nx,ny,nz,ens_num))
  allocate(smoispert (nx,ny,4,ens_num))
  allocate(tslbpert (nx,ny,4,ens_num))
  allocate(tskpert (nx,ny,ens_num))  !

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
  allocate(tskmean (nx,ny))  !!   

  allocate(xlon(nx,ny))
  allocate(xlat(nx,ny))
  allocate(ulon(nx,ny))
  allocate(ulat(nx,ny))
  allocate(vlon(nx,ny))
  allocate(vlat(nx,ny)) 
  allocate(phb(nx,ny,nz))
 ! allocate(smois(nx,ny,nz))
 ! allocate(tslb(nx,ny,nz))
  allocate(hgt(nx,ny)) 
  allocate(mub(nx,ny,ens_num)) 
  allocate(mu (nx,ny,ens_num)) 
  allocate(znw   (nz)) 
  allocate(znu   (nz)) 
  allocate(rdnw  (nz)) 
  allocate(rdn   (nz))
  allocate(dzs   (4))
  allocate(zs    (4)) 
  
  allocate(ranpert(nx,ny,nz,ens_num))
  allocate(ranpert_mean(nx,ny,nz))
  allocate(spd_ranpert(ens_num))
  allocate(ranpert2(nx,ny,nz,ens_num))
  allocate(ranpert2_mean(nx,ny,nz))
  allocate(spd_ranpert2(ens_num))  
  allocate(rec_cdfid(ens_num))  
  allocate(rec_rcode(ens_num))
 
 allocate(ranpertsoil(nx,ny,4,ens_num))
 allocate(ranpertsoil_mean(nx,ny,4))
 allocate(ranpertsoil2(nx,ny,4,ens_num))
 allocate(ranpertsoil2_mean(nx,ny,4))

 allocate(ranpertsoil_2d(nx,ny,ens_num))
 allocate(ranpertsoil_mean2d(nx,ny))
 allocate(ranpertsoil2_2d(nx,ny,ens_num))
 allocate(ranpertsoil2_mean2d(nx,ny))
!stop ok

!print*,'*****ADD SMOIS & TSLB ens, mean name and pert var*****'
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************  
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
!print*,tsk (1,1)
!stop
!  print*,'****ADD VAR SMOIS AND TSLB IN SUBROUTINE get_info_from_cdf!!!****'
  call get_basic_info_2(basic_var_num_2,basic_var_name_2,                    &
                        file_cdfid,file_rcode,nx,ny,                         &
                        mub(:,:,i),mu(:,:,i)                                 &
                        )  
                        
  rec_cdfid(i)=file_cdfid
  rec_rcode(i)=file_rcode
 
  
  ensu (:,:,:,i)   =  u
  ensv (:,:,:,i)   =  v 
  ensw (:,:,:,i)   =  w
  ensph(:,:,:,i)   =  ph
  enst (:,:,:,i)   =  t
  ensqv(:,:,:,i)   =  qv
  ensqr(:,:,:,i)   =  qr
  ensqi(:,:,:,i)   =  qi
  ensqs(:,:,:,i)   =  qs
  ensqgr(:,:,:,i)  =  qgr 
  enssmois(:,:,:,i)   =  smois
  enstslb(:,:,:,i)   =  tslb 
  enstsk(:,:,i) = tsk
  
  
  enddo     ! for ens file loop

!stop
  call get_basic_info_1(basic_var_num_1,basic_var_name_1,                    &
                        rec_cdfid(1),rec_rcode(1),nx,ny,nz,                  &
                        xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,        &
                        znw,znu,p_top,rdnw,rdn                               &
                       )  
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
  ranpert=0
  ranpert2=0
      
  do i=1 , ens_num
  
  seed_ens=i
  
   if(i <= ens_num/2) then
   
    if(opt_pert == 1) then
    
    !print*,rec_id_var

     call addpert(ranpert(:,:,:,i),nx,ny,nz,seed_ens,rec_id_var)
     call addpertsoil(ranpertsoil(:,:,:,i),nx,ny,4,seed_ens)  !for 4 levels soil layers
     call addpertsoil_2d(ranpertsoil_2d(:,:,i),nx,ny,seed_ens)   ! for 2d field, tsk likely? 
!!!!!!!    need rec_idvar for vertical level of soil layer 
    endif ! if(opt_pert == 1)
! stop 
    if(opt_pert == 2) then
     call addpert_local(ranpert(:,:,:,i),nx,ny,nz,seed_ens,pert_cx,pert_cy,pert_radius,dx,dy)  
    endif  
    
    if(opt_pert == 3) then
     call addpert2(ranpert(:,:,:,i),nx,ny,nz,seed_ens,scalar)
     call addpert2(ranpert2(:,:,:,i),nx,ny,nz,seed_ens,scalar)
    endif   
         
   else
     temp1 = i-ens_num/2
     ranpert(:,:,:,i)=-ranpert(:,:,:,temp1)
     ranpertsoil(:,:,:,i)=-ranpertsoil(:,:,:,temp1)
     ranpertsoil_2d(:,:,i)=-ranpertsoil_2d(:,:,temp1) ! for 2d field, tsk likely?

    if(opt_pert == 3) then     
     temp1 = i-ens_num/2
     ranpert2(:,:,:,i)=-ranpert2(:,:,:,temp1)
    endif
     
   endif
   
  enddo 

!!!ok?

  if(enable_recurfilt==1)  then
  print*,'enable_recurfilt'
  
  hradius_temp=hradius
  
  do i=1,ens_num
   call recurfilt_3d(nx,ny,nz,ranpert(:,:,:,i),ipass_filt,hradius_temp,nradius_z)

   
    if(opt_pert == 3) then     
     call recurfilt_3d(nx,ny,nz,ranpert2(:,:,:,i),ipass_filt,hradius_temp*0.5,nradius_z)
    endif
    
  enddo  
  
  ranpert_mean=0
  ranpert2_mean=0
  do k=1,nz
    do j=1,ny
      do i=1,nx
        do ii=1,ens_num
            ranpert_mean(i,j,k)=ranpert_mean(i,j,k)+ranpert(i,j,k,ii)
            if(opt_pert == 3) then
             ranpert2_mean(i,j,k)=ranpert2_mean(i,j,k)+ranpert2(i,j,k,ii)
            endif 
        enddo
            ranpert_mean(i,j,k)=ranpert_mean(i,j,k)/ens_num
            if(opt_pert == 3) then
             ranpert2_mean(i,j,k)=ranpert2_mean(i,j,k)/ens_num
            endif 
      enddo
    enddo
  enddo
  
  spd_ranpert=0
  spd_ranpert2=0
  do ii=1,ens_num
    do k=1,nz
      do j=1,ny
        do i=1,nx  
           spd_ranpert(ii)=spd_ranpert(ii)+(ranpert(i,j,k,ii)-ranpert_mean(i,j,k))**2
           if(opt_pert == 3) then 
            spd_ranpert2(ii)=spd_ranpert2(ii)+(ranpert2(i,j,k,ii)-ranpert2_mean(i,j,k))**2
           endif 
        enddo
      enddo
    enddo  
  enddo
  
  spd_ens_ranpert=0
  spd_ens_ranpert2=0
  do ii=1,ens_num
     spd_ens_ranpert=spd_ens_ranpert+sqrt(spd_ranpert(ii)/(nx*ny*nz))/ens_num
     if(opt_pert == 3) then 
      spd_ens_ranpert2=spd_ens_ranpert2+sqrt(spd_ranpert2(ii)/(nx*ny*nz))/ens_num
     endif 
  enddo

  radio_to_1=1/spd_ens_ranpert
  radio2_to_1=1/spd_ens_ranpert2
  
  print*,'radio_to_1',radio_to_1
  if(opt_pert == 3) then
   print*,'radio2_to_1',radio2_to_1
  endif 
  
  ranpert=ranpert*radio_to_1
  if(opt_pert == 3) then
   ranpert2=ranpert2*radio2_to_1
  endif 
  
  do ii=1,ens_num
    do k=1,nz
      do j=1,ny
        do i=1,nx  
           ranpert(i,j,k,ii)=ranpert(i,j,k,ii)-ranpert_mean(i,j,k)
           if(opt_pert == 3) then
            ranpert2(i,j,k,ii)=ranpert2(i,j,k,ii)-ranpert2_mean(i,j,k)
           endif 
        enddo
      enddo
    enddo  
  enddo

  ranpert_mean=0
  ranpert2_mean=0
  
  do k=1,nz
    do j=1,ny
      do i=1,nx
        do ii=1,ens_num
            ranpert_mean(i,j,k)=ranpert_mean(i,j,k)+ranpert(i,j,k,ii)
            if(opt_pert == 3) then 
             ranpert2_mean(i,j,k)=ranpert2_mean(i,j,k)+ranpert2(i,j,k,ii)
            endif 
        enddo
            ranpert_mean(i,j,k)=ranpert_mean(i,j,k)/ens_num
            if(opt_pert == 3) then
             ranpert2_mean(i,j,k)=ranpert2_mean(i,j,k)/ens_num
            endif
      enddo
    enddo
  enddo
  
  spd_ranpert=0
  spd_ranpert2=0
  do ii=1,ens_num
    do k=1,nz
      do j=1,ny
        do i=1,nx  
           spd_ranpert(ii)=spd_ranpert(ii)+(ranpert(i,j,k,ii)-ranpert_mean(i,j,k))**2
           if(opt_pert == 3) then
            spd_ranpert2(ii)=spd_ranpert2(ii)+(ranpert2(i,j,k,ii)-ranpert2_mean(i,j,k))**2
           endif 
        enddo
      enddo
    enddo  
  enddo
  
  spd_ens_ranpert=0
  spd_ens_ranpert2=0
  do ii=1,ens_num
     spd_ens_ranpert=spd_ens_ranpert+sqrt(spd_ranpert(ii)/(nx*ny*nz))/ens_num
     if(opt_pert == 3) then 
      spd_ens_ranpert2=spd_ens_ranpert2+sqrt(spd_ranpert2(ii)/(nx*ny*nz))/ens_num
     endif 
  enddo
  
  print*,spd_ens_ranpert
  if(opt_pert == 3) then
   print*,spd_ens_ranpert2
  endif
  
  endif    ! if(enable_recurfilt==1)  then
  
  
write(*,'(a40,3f11.6)'),'U center pgrid point before pert:',ensu (nx/2,ny/2,10,1),ensu (nx/2,ny/2,10,ens_num/2),ensu (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'V center pgrid point before pert:',ensv (nx/2,ny/2,10,1),ensv (nx/2,ny/2,10,ens_num/2),ensv (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'W center pgrid point before pert:',ensw (nx/2,ny/2,10,1),ensw (nx/2,ny/2,10,ens_num/2),ensw (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'PH center pgrid point before pert:',ensph (nx/2,ny/2,10,1),ensph (nx/2,ny/2,10,ens_num/2),ensph (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'T center pgrid point before pert:',enst (nx/2,ny/2,10,1),enst (nx/2,ny/2,10,ens_num/2),enst (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'QV center pgrid point before pert:',ensqv (nx/2,ny/2,10,1),ensqv (nx/2,ny/2,10,ens_num/2),ensqv (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'SMOIS center pgrid point before pert:',enssmois (nx/2,ny/2,2,1),enssmois (nx/2,ny/2,2,ens_num/2),enssmois (nx/2,ny/2,2,ens_num)
write(*,'(a40,3f11.6)'),'TSLB center pgrid point before pert:',enstslb (nx/2,ny/2,2,1),enstslb (nx/2,ny/2,2,ens_num/2),enstslb (nx/2,ny/2,2,ens_num)
write(*,'(a40,3f11.6)'),'TSK center pgrid point before pert:',enstsk (nx/2,ny/2,1),enstsk (nx/2,ny/2,ens_num/2),enstsk (nx/2,ny/2,ens_num)

write(*,'(a20,f11.6)'),'  stdu_init=',stdu_init 
write(*,'(a20,f11.6)'),'  stdv_init=',stdv_init
write(*,'(a20,f11.6)'),'  stdw_init=',stdw_init
write(*,'(a20,f11.6)'),'  stdph_init=',stdph_init
write(*,'(a20,f11.6)'),'  stdt_init=',stdt_init
write(*,'(a20,f11.6)'),'  stdqv_init=',stdqv_init
write(*,'(a20,f11.6)'),'  stdqr_init=',stdqr_init
write(*,'(a20,f11.6)'),'  stdqs_init=',stdqs_init
write(*,'(a20,f11.6)'),'  stdqg_init=',stdqg_init
write(*,'(a20,f11.6)'),'  stdphi_init=',stdphi_init
write(*,'(a20,f11.6)'),'  stdpsi_init=',stdpsi_init
write(*,'(a20,f11.6)'),'stdsmois_init=',stdsmois_init
write(*,'(a20,f11.6)'),'stdtslb_init=',stdtslb_init
write(*,'(a20,f11.6)'),' stdtsk_init=',stdtsk_init

!!stop 

!print*,'opt_pert',opt_pert

 if(opt_pert < 3 ) then
  do i=1 , ens_num  
    
    upert(:,:,:,i)  = ranpert(:,:,:,i)* stdu_init
    ensu (:,:,:,i)  = ensu (:,:,:,i) + upert(:,:,:,i)
    
    vpert(:,:,:,i)  = ranpert(:,:,:,i)* stdv_init
    ensv (:,:,:,i)  = ensv (:,:,:,i) + vpert(:,:,:,i)
    
    wpert(:,:,:,i)  = ranpert(:,:,:,i)* stdw_init
    ensw (:,:,:,i)  = ensw (:,:,:,i) + wpert(:,:,:,i)
    
    phpert(:,:,:,i) = ranpert(:,:,:,i)* stdph_init
    ensph (:,:,:,i) = ensph(:,:,:,i) + phpert(:,:,:,i)
    
    tpert(:,:,:,i)  = ranpert(:,:,:,i)* stdt_init
    enst (:,:,:,i)  = enst (:,:,:,i) + tpert(:,:,:,i)
    
    qvpert(:,:,:,i) = ranpert(:,:,:,i)* stdqv_init
    ensqv (:,:,:,i) = ensqv(:,:,:,i) + qvpert(:,:,:,i)
    
    ensqr (:,:,:,i) = ensqr(:,:,:,i) + ranpert(:,:,:,i)* stdqr_init
    ensqs (:,:,:,i) = ensqs(:,:,:,i) + ranpert(:,:,:,i)* stdqs_init
    ensqgr(:,:,:,i) = ensqgr(:,:,:,i) + ranpert(:,:,:,i)* stdqg_init

    smoispert(:,:,:,i) = ranpertsoil(:,:,:,i)*stdsmois_init
    enssmois (:,:,:,i) = enssmois(:,:,:,i) + smoispert(:,:,:,i)
    
    tslbpert(:,:,:,i) = ranpertsoil(:,:,:,i)*stdtslb_init
    enstslb (:,:,:,i) = enstslb(:,:,:,i) + tslbpert(:,:,:,i)

    tskpert(:,:,i) = ranpertsoil_2d(:,:,i)*stdtsk_init
    enstsk (:,:,i) = enstsk(:,:,i) + tskpert(:,:,i)    
   !!!!!!!!!!!!!! 
!    open(10000,file='/home/guoyk/test.txt')
!    write(10000,*)ranpert(:,:,:,i) 
  !!!!!!!!!!!!!!   
 !print*,ranpert(2,5,1,i)
  enddo   

!close(10000)

  do ii=1,ens_num
    
!   if(rec_id_var < 61)then
!      tempkk=nz
!     endif    
!    if(rec_id_var == 62 .or. rec_id_var == 63)then
!       tempkk=4
!    endif 
    
    do k=1,nz
      do j=1,ny
        do i=1,nx 
          if(ensqv (i,j,k,ii) <0) ensqv (i,j,k,ii)=0
          if(ensqr (i,j,k,ii) <0) ensqr (i,j,k,ii)=0
          if(ensqs (i,j,k,ii) <0) ensqs (i,j,k,ii)=0
          if(ensqgr (i,j,k,ii) <0) ensqgr (i,j,k,ii)=0
!!!!!!!!!!! no negtive value of smois
!          if(enssmois (i,j,k,ii) <0) enssmois (i,j,k,ii)=0
!
        enddo
      enddo
    enddo  
  enddo

do ii=1,ens_num
    do k=1,4
      do j=1,ny
        do i=1,nx 
!!!!!!!!!!!! some self-check !!!!!!!!!!!!!!!!!!!!         
          if(abs(tslbpert(i,j,k,ii))>stdtslb_init.or.abs(smoispert(i,j,k,ii))>stdsmois_init)then
            print*,'tslbpert(i,j,k,ii))>stdtslb_init',i,j,k,ii,ranpertsoil(i,j,k,ii),tslbpert(i,j,k,ii)
          endif

!!!!!!!!!!! no negtive value of smois or grater than 1
          if(enssmois (i,j,k,ii) < 0 ) enssmois (i,j,k,ii)=0
          if(enssmois (i,j,k,ii) > 1 ) enssmois (i,j,k,ii)=0.99
!!!!!!!!!!!!! -30 and 30 degree of TSLB is impossible
          if(enstslb (i,j,k,ii) < 233.15 ) enstslb (i,j,k,ii)= 233.15  !!
          if(enstslb (i,j,k,ii) > 313.15 ) enstslb (i,j,k,ii)= 313.15
        enddo
      enddo
    enddo  
  enddo

do ii=1,ens_num
!    do k=1,4
      do j=1,ny
        do i=1,nx 
!!!!!!!!!!!! some self-check !!!!!!!!!!!!!!!!!!!!         
!          if(abs(tskpert(i,j,ii))>stdtsk_init)then
!            print*,'tskpert(i,j,ii))>stdtsk_init',i,j,ii,ranpertsoil_2d(i,j,ii),tskpert(i,j,ii)
!          endif
!!!!!!!!!!! no negtive value of smois or grater than 1
!          if(enssmois (i,j,k,ii) < 0 ) enssmois (i,j,k,ii)=0
!          if(enssmois (i,j,k,ii) > 1 ) enssmois (i,j,k,ii)=0.99
!!!!!!!!!!!!! -30 and 30 degree of TSLB is impossible
          if(enstsk (i,j,ii) < 233.15 ) enstsk (i,j,ii)= 233.15  !!! TSK !!!
          if(enstsk (i,j,ii) > 313.15 ) enstsk (i,j,ii)= 313.15
        enddo
      enddo
!    enddo  
enddo


 
 else
   
   call addpert_with_physical_constrain(xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,                                     &
                                        znw,znu,p_top,rdnw,rdn,                                                    &
                                        mub,mu,								        &
                                        ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb, &
                                        ranpert,ranpert2,                                                          &
                                        nx,ny,nz,ens_num,dx,dy,                                                    &
                                        stdphi_init,stdpsi_init,stdqv_init,stdt_init,                              &
                                        upert,vpert,wpert,phpert,tpert,qvpert,smoispert,tslbpert                   &
                                        )
                                        
    do ii=1,ens_num
      do k=1,nz
        do j=1,ny
          do i=1,nx 
            if(ensqv (i,j,k,ii) <0) ensqv (i,j,k,ii)=0
            if(ensqr (i,j,k,ii) <0) ensqr (i,j,k,ii)=0
            if(ensqs (i,j,k,ii) <0) ensqs (i,j,k,ii)=0
            if(ensqgr (i,j,k,ii) <0) ensqgr (i,j,k,ii)=0
         enddo
       enddo
      enddo  
    enddo   
                                         
 endif  ! if(opt_pert < 3 ) then
 
write(*,'(a40,3f11.6)'),'U center pgrid point after pert:',ensu (nx/2,ny/2,10,1),ensu (nx/2,ny/2,10,ens_num/2),ensu (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'V center pgrid point after pert:',ensv (nx/2,ny/2,10,1),ensv (nx/2,ny/2,10,ens_num/2),ensv (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'W center pgrid point after pert:',ensw (nx/2,ny/2,10,1),ensw (nx/2,ny/2,10,ens_num/2),ensw (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'PH center pgrid point after pert:',ensph (nx/2,ny/2,10,1),ensph (nx/2,ny/2,10,ens_num/2),ensph (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'T center pgrid point after pert:',enst (nx/2,ny/2,10,1),enst (nx/2,ny/2,10,ens_num/2),enst (nx/2,ny/2,10,ens_num)
write(*,'(a40,3f11.6)'),'QV center pgrid point after pert:',ensqv (nx/2,ny/2,10,1),ensqv (nx/2,ny/2,10,ens_num/2),ensqv (nx/2,ny/2,10,ens_num)
!  print*,'SMOIS center pgrid point after pert:',enssmois (nx/2,ny/2,2,1),enssmois (nx/2,ny/2,2,ens_num/2),enssmois (nx/2,ny/2,2,ens_num)
 ! print*,'TSLB center pgrid point after pert:',enstslb (nx/2,ny/2,2,1),enstslb (nx/2,ny/2,2,ens_num/2),enstslb (nx/2,ny/2,2,ens_num)
 ! print*,'TSK center pgrid point after pert:',enstsk (nx/2,ny/2,1),enstsk (nx/2,ny/2,ens_num/2),enstsk (nx/2,ny/2,ens_num)
write(*,'(a40,3f11.6)'),'SMOIS center pgrid point after pert:',enssmois (nx/2,ny/2,2,1),enssmois (nx/2,ny/2,2,ens_num/2),enssmois (nx/2,ny/2,2,ens_num)
write(*,'(a40,3f11.6)'),'TSLB center pgrid point after pert:',enstslb (nx/2,ny/2,2,1),enstslb (nx/2,ny/2,2,ens_num/2),enstslb (nx/2,ny/2,2,ens_num)
write(*,'(a40,3f11.6)'),'TSK center pgrid point after pert:',enstsk (nx/2,ny/2,1),enstsk (nx/2,ny/2,ens_num/2),enstsk (nx/2,ny/2,ens_num)
!stop
   
  call mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,        &
                nx,ny,nz,ens_num,                                                                        &
                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean &
               )  
  call ens_spread(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,      &
                  nx,ny,nz,ens_num,                                                                      &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,&
                  spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk  &
                 )  
                 
 ! print*, spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk                
  
 ! print*,'writing ensemble files'

if(0==1) then
  
  print*,'writing ensemble files for post analysis such as Grads etc.'  
  
  open(8553,file='init_pert_ens.dat',form='binary')
  
  do ii=1,ens_num
  
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) upert(i,j,k,ii)
          enddo
       enddo
    enddo
    
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) vpert(i,j,k,ii)
          enddo
       enddo
    enddo    
  
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) wpert(i,j,k,ii)
          enddo
       enddo
    enddo

    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) phpert(i,j,k,ii)
          enddo
       enddo
    enddo      

    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) tpert(i,j,k,ii)
          enddo
       enddo
    enddo
    
    do k=1,nz
       do j=1,ny
          do i=1,nx
              write(8553) qvpert(i,j,k,ii)
          enddo
       enddo
    enddo    
  enddo

!    do k=1,4
!       do j=1,ny
!          do i=1,nx
!              write(8553) smoispert(i,j,k,ii)
!          enddo
!       enddo
!    enddo 

!    do k=1,4
!       do j=1,ny
!          do i=1,nx
!              write(8553) tslbpert(i,j,k,ii)
!          enddo
!       enddo
!    enddo 
!      print*,'finished writing ens',ii
  close(8553)

endif
  
  print*,'writing ensemble files for EnSRF assimilation cycle'  
  
  do i=1 , ens_num
  
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
!print*,rec_dims(13,0),rec_dims(13,1),rec_dims(13,2),rec_dims(13,3)
  call output_data(u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz,                            &
                   analysis_var_num,rec_id_var,rec_cdfid(i),rec_dims,rec_dims2,rec_rcode(i)       &
                   )
  
  enddo   
!stop  
  do i=1 , ens_num
  call ncclos(rec_cdfid(i),rec_rcode(i))
  enddo
  
  print*,'***********************************'
  print*,'         Finished                  '
  print*,'***********************************'
    
  end program init_ens
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
