 program plot_sim_radar
 !----------------------------------------------------------------
 !Purpose: 
 !        Main program of obs_sim
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
 
  real,allocatable      :: xlon(:,:),xlat(:,:),ulon(:,:),ulat(:,:),vlon(:,:),vlat(:,:),phb(:,:,:),hgt(:,:),mub(:,:),mu(:,:), &
                           znw(:),znu(:),rdnw(:),rdn(:)
  real                  :: p_top
  real                  :: h_top
  
  character(len=5)      :: radar_name
  integer               :: xradar,yradar
  real                  :: hradar
  integer               :: lev_num
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
  real,allocatable      :: rv(:,:,:),rf(:,:,:)                 
  
  integer               :: irec,iradar
  real                  :: tem1,tem2,tem3,tem4,tem5
            
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


do iradar=1,radar_num

temp1=int(iradar/10      )-int(iradar/100     )*10  
temp2=int(iradar/1       )-int(iradar/10      )*10
lab_radar_num = char(48+temp1)//char(48+temp2) 

print*,'the ',iradar,'th radar'

temp1=nx*ny
open(1081,file='radar_sim_info_'//lab_radar_num//'.dat')


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

if(allocated(rv)) deallocate(rv)
allocate(rv(nx,ny,lev_num))
if(allocated(rf)) deallocate(rf)
allocate(rf(nx,ny,lev_num))

rv=0
rf=0

temp3=0
temp4=0
temp5=0

  do j=1,ny
  do i=1,nx
  
    if(refcom(i,j)>=20) then 
      
      temp1=j-yradar
      temp2=i-xradar
      tem1=real(temp1)/real(temp2)
      azimuth=atan(tem1)
      hor_dis= sqrt(((temp1)*dy)**2+((temp2)*dx)**2)
      
      do k=1,lev_num

       vert_dis= hor_dis*tan(evl(k)*3.1415926/180)+hradar
             
       call radardata_kr_cal(nx,ny,nz,i,j,vert_dis,phb,krw,krm)    ! krw and krw is the level just lower than vert_dis        

       if(krw >= 1 .and. krm >= 1 .and.  vert_dis < phb(i,j,nz)/9.81) then
       
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
              
       else
       	
       rv_make=miss_data
       rf_make=miss_data	
       
       endif
       
       temp3=temp3+1
       if(rv_make == miss_data) temp4=temp4+1
       if(rf_make == miss_data) temp5=temp5+1
       
       if(i==nx/2.and.j==ny/2) then
       PRINT*,'rv,rf:',rv_make,rf_make,i,j,vert_dis
       endif
       
       if(rv_make .ne. miss_data) then
         rv(i,j,k)=rv_make
       else
       	 rv(i,j,k)=0
       endif
       
       if(rf_make .ne. miss_data) then
         rf(i,j,k)=rf_make
       else
       	 rf(i,j,k)=0
       endif       	
       	
      enddo
      
    endif
    
  enddo
  enddo
  
  write(*,*) 'THE NUMBER OF RADAR DATA'
  write(*,'(I10)'),temp3
  write(*,*) 'THE NUMBER OF MISS RV'
  write(*,'(I10)'),temp4
  write(*,*) 'THE NUMBER OF MISS RF'
  write(*,'(I10)'),temp5

                                                       
     
open(1083,file='plot_radar_'//lab_radar_num//'_'//lab_time_num//'.dat',form='binary')

write(1083) rv
write(1083) rf
print*,'===================RV=========================' 
do j=ny/2-5,ny/2+5 
write(*,'(11F10.3)') (rv(i,j,lev_num/2),i=nx/2-5,nx/2+5) 
enddo
print*,'=============================================='
print*,'===================RF========================='
do j=ny/2-5,ny/2+5 
write(*,'(11F10.3)') (rf(i,j,lev_num/2),i=nx/2-5,nx/2+5) 
enddo 
  
close(1081)
close(1083)  

enddo         !   for radar number loop

end program plot_sim_radar