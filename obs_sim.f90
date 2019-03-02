  program obs_sim
 !----------------------------------------------------------------
 !Purpose: 
 !        Main program of obs_sim
 !
 !---------------------------------------------------------------- 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            program obs_sim
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
  integer               :: length_obs_dir      
  integer,parameter     :: analysis_var_num=10
  integer,parameter     :: basic_var_num_1=13
  integer,parameter     :: basic_var_num_2=2
  integer               :: i,j,k,ii,jj,kk,kr
  integer               :: istart(4), iend(4)
  character (len=4)     :: lab_ens
  character (len=2)     :: lab_domain
  character (len=2)     :: lab_time_num
  character (len=2)     :: lab_radar_num
  
  integer               :: rec_id_var(analysis_var_num)
  integer               :: rec_cdfid
  integer               :: rec_dims(analysis_var_num,3)
  integer               :: rec_rcode
  character(len=80)     :: analysis_var_name(analysis_var_num)
  integer               :: file_cdfid,file_rcode
  character(len=80)     :: basic_var_name_1(basic_var_num_1)
  character(len=80)     :: basic_var_name_2(basic_var_num_2)
    
  real,allocatable      :: u(:,:,:),v(:,:,:),w(:,:,:),ph(:,:,:),t(:,:,:),qv(:,:,:),qr(:,:,:),qi(:,:,:),qs(:,:,:),qgr(:,:,:)                 
  integer               :: temp1,temp2,temp3,temp4,temp5
  real,allocatable      :: rho(:,:,:),ref_x(:,:,:),p(:,:,:),tc(:,:,:),refcom(:,:)
  real                  :: maxref_vert,iref,jref
  
  
  real,allocatable      :: obsu(:),obsv(:),obsw(:),obsph(:),obst(:),obsqv(:),obsqr(:),obsqi(:),obsqs(:),obsqgr(:)
  integer,allocatable   :: xyzu(:,:),xyzv(:,:),xyzw(:,:),xyzph(:,:),xyzt(:,:),xyzqv(:,:),xyzqr(:,:),xyzqi(:,:),xyzqs(:,:),xyzqgr(:,:)

  
  real,allocatable      :: xlon(:,:),xlat(:,:),ulon(:,:),ulat(:,:),vlon(:,:),vlat(:,:),phb(:,:,:),hgt(:,:),mub(:,:),mu(:,:), &
                           znw(:),znu(:),rdnw(:),rdn(:)
  real                  :: p_top
  real                  :: h_top
  
  character(len=5)      :: radar_name
  integer               :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num
  real,allocatable      :: rv(:),rf(:)
  integer,allocatable   :: ir(:),jr(:)
  real,allocatable      :: hr(:)
  integer               :: radar_data_num
  real                  :: hr_make
  real                  :: rv_make,rf_make
  real                  :: azimuth
  real,allocatable      :: evl(:)
  real                  :: hor_dis
  real                  :: vert_dis
  integer               :: krm,krw
  real                  :: ur,vr,wr,qrr,qsr,qgrr,tcr,rhor
  real                  :: ref_r,ref_s,ref_h
  integer               :: iradar
  real                  :: tem1,tem2,tem3,tem4,tem5
  
  real GASDEV
  integer seed
            
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
  read (5,nml =        obs_info)
  write(6,nml =        obs_info)  
  close(5)


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
  
  allocate(rho(nx,ny,nz))
  allocate(ref_x(nx,ny,nz))
  allocate(p(nx,ny,nz))
  allocate(tc(nx,ny,nz))
  allocate(refcom(nx,ny))
  
  allocate(xlon(nx,ny))
  allocate(xlat(nx,ny))
  allocate(ulon(nx,ny))
  allocate(ulat(nx,ny))
  allocate(vlon(nx,ny))
  allocate(vlat(nx,ny)) 
  allocate(phb(nx,ny,nz))
  allocate(hgt(nx,ny)) 
  allocate(mub(nx,ny)) 
  allocate(mu (nx,ny)) 
  allocate(znw   (nz)) 
  allocate(znu   (nz)) 
  allocate(rdnw  (nz)) 
  allocate(rdn   (nz))   
   
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------  
!                       Reading true files
!
!---------------------------------------------------------------------------------------  
!---------------------------------------------------------------------------------------   
  
  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)
  
  temp1=int(time_num/10      )-int(time_num/100     )*10  
  temp2=int(time_num/1       )-int(time_num/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)  
  

  length_file_dir =len_trim(obs_file_dir)
  length_ens_head=len_trim(ens_file_head)
  
  input_file_name=''//obs_file_dir(1:length_file_dir)//'/'//ens_file_head(1:length_ens_head)//'_d'//lab_domain//'_'//lab_time_num//''
  
  print*,"INPUT FILE IS: ",trim(input_file_name)
  length_input = len_trim(input_file_name)
  
  
! Now read the file
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                &
                          rec_id_var,file_cdfid,rec_dims,file_rcode,                                      &
                          u,v,w,ph,t,qv,qr,qi,qs,qgr,nx,ny,nz                                             &
                          ) 
                          
  print*,'finished reading fields'                        


                          
  call get_basic_info_1(basic_var_num_1,basic_var_name_1,                    &
                        file_cdfid,file_rcode,nx,ny,nz,                      &
                        xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,               &
                        znw,znu,p_top,rdnw,rdn                               &
                       )   
                       
  call get_basic_info_2(basic_var_num_2,basic_var_name_2,                    &
                        file_cdfid,file_rcode,nx,ny,                         &
                        mub,mu                                               &
                        )                                          
  
  print*,'finished reading coordinate info'                         


  call ncclos(file_cdfid,file_rcode)
  
  
  do k=1,nz-1
  do j=1,ny
  do i=1,nx
  
  call cal_rho(mub(i,j),mu(i,j),qv(i,j,k),ph(i,j,k),ph(i,j,k+1),rdnw(k),p_top,znu(k),t(i,j,k),       &
               nx,ny,nz,rho(i,j,k),p(i,j,k)                                                          &                  
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
                  
  call cal_tc(t(i,j,k),p(i,j,k),tc(i,j,k),nx,ny,nz)  

  enddo
  enddo
  enddo
  
  ref_x=0
  refcom=0
  
  do k=1,nz
  do j=1,ny
  do i=1,nx
     
     call ref_operator(qr(i,j,k),qs(i,j,k),qgr(i,j,k),ref_x(i,j,k),ref_r,ref_s,ref_h,rho(i,j,k),tc(i,j,k)) 


     if(i==nx/2 .and. j==ny/2) then
      print*,'ref_X at center is:',ref_x(i,j,k),tc(i,j,k)
     endif
               
  enddo
  enddo
  enddo 
  

  do j=1,ny
  do i=1,nx    
    maxref_vert=0
    do k=1,nz
      if(maxref_vert < ref_x(i,j,k) ) then
        maxref_vert=ref_x(i,j,k)
      endif
    enddo
    refcom(i,j)=maxref_vert 
  enddo
  enddo 
 
if(obs_type == 0 .or. obs_type == 1 ) then  
  open(1001,file =  './obs_u_'//lab_time_num//'.dat')
  open(1002,file =  './obs_v_'//lab_time_num//'.dat')
  open(1003,file =  './obs_w_'//lab_time_num//'.dat')
  open(1004,file = './obs_ph_'//lab_time_num//'.dat')
  open(1005,file =  './obs_t_'//lab_time_num//'.dat')
  open(1006,file = './obs_qv_'//lab_time_num//'.dat')
  open(1007,file = './obs_qr_'//lab_time_num//'.dat')
  open(1008,file = './obs_qi_'//lab_time_num//'.dat')
  open(1009,file = './obs_qs_'//lab_time_num//'.dat')
  open(1010,file ='./obs_num_'//lab_time_num//'.dat')
  open(1011,file ='./obs_qgr_'//lab_time_num//'.dat')
  
  temp1=0
  temp2=0
  temp3=0

  do j=1,ny,sim_obs_den
  do i=1,nx,sim_obs_den
  
    if(refcom(i,j)>=20) then
      do k=2,nz/2        
       if( ref_x(i,j,k) > 10.0 ) then
       
        write(1001,500) (u(i,j,k)+u(i+1,j,k))*0.5  ,i,j,k  
        write(1002,500) (v(i,j,k)+v(i,j+1,k))*0.5  ,i,j,k  
        write(1003,500) (w(i,j,k)+w(i,j,k+1))*0.5  ,i,j,k  
        write(1004,500) (ph(i,j,k)+ph(i,j,k+1))*0.5,i,j,k  
        write(1005,500) t(i,j,k) ,i,j,k  
        write(1006,600) qv(i,j,k),i,j,k  
        write(1007,600) qr(i,j,k),i,j,k  
        write(1008,600) qi(i,j,k),i,j,k  
        write(1009,600) qs(i,j,k),i,j,k
        write(1011,600) qgr(i,j,k),i,j,k 
        
       else
        
        write(1001,500) miss_data,i,j,k  
        write(1002,500) miss_data,i,j,k  
        write(1003,500) miss_data,i,j,k  
        write(1004,500) miss_data,i,j,k  
        write(1005,500) miss_data,i,j,k  
        write(1006,600) miss_data,i,j,k  
        write(1007,600) miss_data,i,j,k  
        write(1008,600) miss_data,i,j,k  
        write(1009,600) miss_data,i,j,k
        write(1011,600) miss_data,i,j,k 
       temp3=temp3+1  
       endif   
     temp1=temp1+1         
     enddo
   endif
   
  enddo
  enddo 

   
  write(1010,*) 'the number of obs grid is:'
  write(1010,*) temp1

  
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
  
length_obs_dir=len_trim(obs_file_dir)  

  open(1001,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_u_'//lab_time_num//'.dat')
  open(1002,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_v_'//lab_time_num//'.dat')
  open(1003,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_w_'//lab_time_num//'.dat')
  open(1004,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_ph_'//lab_time_num//'.dat')
  open(1005,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_t_'//lab_time_num//'.dat')
  open(1006,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qv_'//lab_time_num//'.dat')
  open(1007,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qr_'//lab_time_num//'.dat')
  open(1008,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qi_'//lab_time_num//'.dat')
  open(1009,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qs_'//lab_time_num//'.dat')
  open(1010,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_num_'//lab_time_num//'.dat')
  open(1011,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qgr_'//lab_time_num//'.dat')

  read(1010,'(A)')
  read(1010,*) temp2
  print*,'the number of obs grid is:',temp2
  print*,'the number of missing data obs grid is:'
  print*, temp3  
    
  allocate  (obsu(temp2))
  allocate  (obsv(temp2))
  allocate  (obsw(temp2))
  allocate (obsph(temp2))
  allocate  (obst(temp2))
  allocate (obsqv(temp2))
  allocate (obsqr(temp2))
  allocate (obsqi(temp2))
  allocate (obsqs(temp2))
  allocate(obsqgr(temp2))
  
  allocate  (xyzu(3,temp2))
  allocate  (xyzv(3,temp2))
  allocate  (xyzw(3,temp2))
  allocate (xyzph(3,temp2))
  allocate  (xyzt(3,temp2))
  allocate (xyzqv(3,temp2))
  allocate (xyzqr(3,temp2))
  allocate (xyzqi(3,temp2))
  allocate (xyzqs(3,temp2))
  allocate(xyzqgr(3,temp2))
  
  print*,'allocated normally'
  
  do i=1,temp2
    
    read(1001,500),obsu(i),xyzu(1,i),xyzu(2,i),xyzu(3,i)
    read(1002,500),obsv(i),xyzv(1,i),xyzv(2,i),xyzv(3,i)
    read(1003,500),obsw(i),xyzw(1,i),xyzw(2,i),xyzw(3,i)
    read(1004,500),obsph(i),xyzph(1,i),xyzph(2,i),xyzph(3,i)
    read(1005,500),obst(i),xyzt(1,i),xyzt(2,i),xyzt(3,i)
    read(1006,600),obsqv(i),xyzqv(1,i),xyzqv(2,i),xyzqv(3,i)
    read(1007,600),obsqr(i),xyzqr(1,i),xyzqr(2,i),xyzqr(3,i)
    read(1008,600),obsqi(i),xyzqi(1,i),xyzqi(2,i),xyzu(3,i)
    read(1009,600),obsqs(i),xyzqs(1,i),xyzqs(2,i),xyzqs(3,i)
    read(1011,600),obsqgr(i),xyzqgr(1,i),xyzqgr(2,i),xyzqgr(3,i)
    
 enddo


 print*,'reading the obs files,following are the samples'
 print*,obsu(temp2/2),xyzu(1,temp2/2),xyzu(2,temp2/2),xyzu(3,temp2/2)
 print*,obsv(temp2/2),xyzv(1,temp2/2),xyzv(2,temp2/2),xyzv(3,temp2/2)
 print*,obsw(temp2/2),xyzw(1,temp2/2),xyzw(2,temp2/2),xyzw(3,temp2/2)
 print*,obsph(temp2/2),xyzph(1,temp2/2),xyzph(2,temp2/2),xyzph(3,temp2/2)
 print*,obst(temp2/2),xyzt(1,temp2/2),xyzt(2,temp2/2),xyzt(3,temp2/2)
 print*,obsqv(temp2/2),xyzqv(1,temp2/2),xyzqv(2,temp2/2),xyzqv(3,temp2/2)
 print*,obsqr(temp2/2),xyzqr(1,temp2/2),xyzqr(2,temp2/2),xyzqr(3,temp2/2)
 print*,obsqi(temp2/2),xyzqi(1,temp2/2),xyzqi(2,temp2/2),xyzu(3,temp2/2)
 print*,obsqs(temp2/2),xyzqs(1,temp2/2),xyzqs(2,temp2/2),xyzqs(3,temp2/2)
 print*,obsqgr(temp2/2),xyzqgr(1,temp2/2),xyzqgr(2,temp2/2),xyzqgr(3,temp2/2) 
 
 print*,'the error for each var' 
 print*,obs_err_u  
 print*,obs_err_v 
 print*,obs_err_w  
 print*,obs_err_ph 
 print*,obs_err_t 
 print*,obs_err_qv
 print*,obs_err_qr
 print*,obs_err_qi
 print*,obs_err_qs
 print*,obs_err_qgr 
  
500 FORMAT(F16.4,3I4) 
600 FORMAT(F16.8,3I4) 

endif             !!!!!!!!!!!!!!!!!!!!!!!!!!!!! case for grid to grid sim obs
 

if(obs_type == 0 .or. obs_type == 2 ) then 

print*,obs_err_rf,obs_err_rv
do iradar=1,radar_num

temp1=int(iradar/10      )-int(iradar/100     )*10  
temp2=int(iradar/1       )-int(iradar/10      )*10
lab_radar_num = char(48+temp1)//char(48+temp2) 

print*,'the ',iradar,'th radar'

open(1081,file='radar_sim_info_'//lab_radar_num//'.dat')
open(1082,file='obs_rad_num_'//lab_radar_num//'_'//lab_time_num//'.dat')
open(1083,file='obs_rv_'//lab_radar_num//'_'//lab_time_num//'.dat',form='formatted')
open(1084,file='obs_rf_'//lab_radar_num//'_'//lab_time_num//'.dat',form='formatted')
open(1085,file='obs_radardata_position_'//lab_radar_num//'_'//lab_time_num//'.dat')

read(1081,'(A)')
read(1081,'(5A)') radar_name  
print*,'radar_name:',radar_name
read(1081,'(A)')
read(1081,'(I3,1X,I3)')  xradar,yradar
print*,'radar site:',xradar,yradar
read(1081,'(A)')
read(1081,'(F4.1)') hradar
print*,'radar height:',hradar
read(1081,'(A)')
read(1081,'(I2)') lev_num
print*,'level number:',lev_num
if(allocated(evl)) deallocate(evl)
allocate(evl(lev_num))
print*,'level angle:'
read(1081,'(A)')
do i=1,lev_num
 read(1081,'(F5.2)') evl(i)
 print*,evl(i)
enddo


temp3=0
temp4=0
temp5=0
seed=-777



  do j=1,ny,sim_obs_den
  do i=1,nx,sim_obs_den
  
    if(refcom(i,j)>=20) then 
      
      temp1=j-yradar
      temp2=i-xradar
      tem1=real(temp1)/real(temp2)
      azimuth=atan(tem1)
      hor_dis= sqrt(((temp1)*dy)**2+((temp2)*dx)**2)
      
      do k=1,lev_num

       vert_dis= hor_dis*tan(evl(k)*3.1415926/180)+hradar
             
       call radardata_kr_cal(nx,ny,nz,i,j,vert_dis,phb,krw,krm)    ! krw and krw is the level just lower than vert_dis        

       if(krw >= 1 .and. krm >= 1 .and.  vert_dis < phb(i,j,nz-2)/9.81) then
       
       call interpolation_sim_radar_wind(nx,ny,nz,i ,j ,krw,krm,vert_dis,phb,                            &
                                         u(i,j,krm+1),u(i+1,j,krm+1),u(i,j,krm  ),u(i+1,j,krm  ),        &
                                         v(i,j,krm+1),v(i,j+1,krm+1),v(i,j,krm  ),v(i,j+1,krm),          &
                                         w(i,j,krw+1),w(i,j,krw  ),                                      &
                                         ur,vr,wr                                                        &
                                         )
                                        
       call interpolation_sim_radar_mass(nx,ny,nz,i ,j ,krm,vert_dis,phb,                                     &
                                         qr(i,j,krm+1),qr(i,j,krm  ),qs(i,j,krm+1),qs(i,j,krm  ),             &
                                         qgr(i,j,krm+1),qgr(i,j,krm  ),rho(i,j,krm+1),rho(i,j,krm  ),         & 
                                         tc(i,j,krm+1),tc(i,j,krm  ),                                         &                      
                                         qrr,qsr,qgrr,rhor,tcr                                                &
                                        )


       call ref_maker(qrr,qsr,qgrr,rf_make,ref_r,ref_s,ref_h,rhor,tcr)
                                        
       call rv_maker(ur,vr,wr,rv_make,azimuth,evl(k)*3.1415926/180,rf_make,vert_dis,ref_r,ref_s,ref_h,   &
                    qrr,qsr,qgrr,rhor                                                                    &
                    )
       
      ! rf_make =rf_make+GASDEV(seed)*obs_err_rf
      ! rv_make =rv_make+GASDEV(seed)*obs_err_rv
              
       else
       	
       rv_make=miss_data
       rf_make=miss_data	
       
       endif
       
       if(abs(rv_make) .gt. 100 ) rv_make=miss_data
       if(abs(rf_make) .gt. 100 ) rf_make=miss_data
       
       write(1083,'(F16.3)') rv_make
       write(1084,'(F16.3)') rf_make
       write(1085,'(2I4,F16.3)') i,j,vert_dis
       temp3=temp3+1
       if(rv_make == miss_data) temp4=temp4+1
       if(rf_make == miss_data) temp5=temp5+1
       
       if(i==nx/2.and.j==ny/2) then
       PRINT*,'rv,rf:',rv_make,rf_make,i,j,vert_dis
       endif
       
      enddo
      
    endif
    
  enddo
  enddo
  
  write(1082,*) 'THE NUMBER OF RADAR DATA'
  write(1082,'(I10)'),temp3
  write(*,*) 'THE NUMBER OF MISS RV'
  write(*,'(I10)'),temp4
  write(*,*) 'THE NUMBER OF MISS RF'
  write(*,'(I10)'),temp5
      
close(1081)
close(1082)
close(1083)
close(1084)
close(1085)

enddo                    !     for radar number loop

temp1=0
temp2=0
temp3=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,'Now will test whether the radar data io is working normally'

do iradar=1,radar_num

temp1=int(iradar/10      )-int(iradar/100     )*10  
temp2=int(iradar/1       )-int(iradar/10      )*10
lab_radar_num = char(48+temp1)//char(48+temp2) 

print*,'the ',iradar,'th radar'

open(1081,file='radar_sim_info_'//lab_radar_num//'.dat')
open(1082,file='obs_rad_num_'//lab_radar_num//'_'//lab_time_num//'.dat')
open(1083,file='obs_rv_'//lab_radar_num//'_'//lab_time_num//'.dat')
open(1084,file='obs_rf_'//lab_radar_num//'_'//lab_time_num//'.dat')
open(1085,file='obs_radardata_position_'//lab_radar_num//'_'//lab_time_num//'.dat')

read(1081,'(A)')
read(1081,'(5A)') radar_name  
print*,'radar_name:',radar_name
read(1081,'(A)')
read(1081,'(I3,1X,I3)')  xradar,yradar
print*,'radar site:',xradar,yradar
read(1081,'(A)')
read(1081,'(F4.1)') hradar
print*,'radar height:',hradar
read(1081,'(A)')
read(1081,'(I2)') lev_num
print*,'level number:',lev_num
if(allocated(evl)) deallocate(evl)
allocate(evl(lev_num))
print*,'level angle:'
read(1081,'(A)')
do i=1,lev_num
 read(1081,'(F5.2)') evl(i)
 print*,evl(i)
enddo

read(1082,'(A)')
read(1082,'(I10)'),temp3
write(*,*) 'THE NUMBER OF RADAR DATA:',temp3
if(allocated(rv)) deallocate(rv)
allocate(rv(temp3))
if(allocated(rf)) deallocate(rf)
allocate(rf(temp3))
if(allocated(ir)) deallocate(ir)
allocate(ir(temp3))
if(allocated(jr)) deallocate(jr)
allocate(jr(temp3))
if(allocated(hr)) deallocate(hr)
allocate(hr(temp3))


do i=1,temp3
  read(1083,*) rv(i)
  read(1084,*) rf(i)
  read(1085,'(2I4,F16.3)') ir(i),jr(i),hr(i)
enddo

print*,'reading the obs files,following are the samples'

do i=1,temp3 
  if(ir(i)== nx/2 .and. jr(i) == ny/2) then
   write(*,*) rv(i),ir(i),jr(i),hr(i)
   write(*,*) rf(i),ir(i),jr(i),hr(i)
   print*,'=============================='
  endif
enddo

close(1081)
close(1082)
close(1083)
close(1084)
close(1085)


enddo                    !     for radar number loop

endif    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! case for radar sim obs
       
  end program obs_sim