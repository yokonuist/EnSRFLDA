subroutine read_sim_sound_dimension(lab_time_num,obs_sim_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_sim_sound_dimension
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'

integer            :: length_obs_dir
character(len=2)   :: lab_time_num
integer            :: obs_sim_num

  length_obs_dir=len_trim(obs_file_dir)
  
  open(1001,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_u_'//lab_time_num//'.dat'  )
  open(1002,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_v_'//lab_time_num//'.dat'  )
  open(1003,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_w_'//lab_time_num//'.dat'  )
  open(1004,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_ph_'//lab_time_num//'.dat' )
  open(1005,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_t_'//lab_time_num//'.dat'  )
  open(1006,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qv_'//lab_time_num//'.dat' )
  open(1007,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qr_'//lab_time_num//'.dat' )
  open(1008,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qi_'//lab_time_num//'.dat' )
  open(1009,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qs_'//lab_time_num//'.dat' )
  open(1010,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_num_'//lab_time_num//'.dat')
  open(1011,file =  ''//obs_file_dir(1:length_obs_dir)//'/obs_qgr_'//lab_time_num//'.dat')

  read(1010,'(A)')
  read(1010,*) obs_sim_num
  print*,'the number of obs grid is:',obs_sim_num
    
end subroutine read_sim_sound_dimension

subroutine read_sim_sound(xyzu,xyzv,xyzw,xyzph,xyzt,xyzqv,xyzqr,xyzqs,xyzqi,xyzqgr,      &
                          obsu,obsv,obsw,obsph,obst,obsqv,obsqr,obsqs,obsqi,obsqgr,      &
                          obs_sim_num)

integer      :: obs_sim_num
real         :: obsu(obs_sim_num),obsv(obs_sim_num),obsw(obs_sim_num),obsph(obs_sim_num),   &
                obst(obs_sim_num),obsqv(obs_sim_num),obsqr(obs_sim_num),obsqi(obs_sim_num), &
                obsqs(obs_sim_num),obsqgr(obs_sim_num)

integer      :: xyzu(3,obs_sim_num),xyzv(3,obs_sim_num),xyzw(3,obs_sim_num),xyzph(3,obs_sim_num),   &
                xyzt(3,obs_sim_num),xyzqv(3,obs_sim_num),xyzqr(3,obs_sim_num),xyzqi(3,obs_sim_num), &
                xyzqs(3,obs_sim_num),xyzqgr(3,obs_sim_num) 

  
  do i=1,obs_sim_num
    
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
 
 print*,'finished reading all sim obs'
 
500 FORMAT(F16.4,3I4) 
600 FORMAT(F16.8,3I4) 


end subroutine read_sim_sound



subroutine read_sim_radar_dimension(lab_time_num,iradar,radar_name,xradar,yradar,     &
                                    hradar,lev_num,ra_data_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_sim_radar_dimension
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'

character(len=2)      :: lab_time_num
character(len=2)      :: lab_radar_num
character(len=5)      :: radar_name
integer               :: xradar,yradar
real                  :: hradar
integer               :: lev_num  
integer               :: iradar
real,allocatable      :: evl(:)
integer               :: ra_data_num
integer               :: temp1,temp2,temp3

temp1=int(iradar/10      )-int(iradar/100     )*10  
temp2=int(iradar/1       )-int(iradar/10      )*10
lab_radar_num = char(48+temp1)//char(48+temp2) 

print*,'the ',iradar,'th radar'

length_obs_dir=len_trim(obs_file_dir)
open(1081,file=''//obs_file_dir(1:length_obs_dir)//'/radar_sim_info_'//lab_radar_num//'.dat')
open(1082,file=''//obs_file_dir(1:length_obs_dir)//'/obs_rad_num_'//lab_radar_num//'_'//lab_time_num//'.dat')
open(1083,file=''//obs_file_dir(1:length_obs_dir)//'/obs_rv_'//lab_radar_num//'_'//lab_time_num//'.dat')
open(1084,file=''//obs_file_dir(1:length_obs_dir)//'/obs_rf_'//lab_radar_num//'_'//lab_time_num//'.dat')
open(1085,file=''//obs_file_dir(1:length_obs_dir)//'/obs_radardata_position_'//lab_radar_num//'_'//lab_time_num//'.dat')

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
read(1082,'(I10)'),ra_data_num
write(*,*) 'THE NUMBER OF RADAR DATA:',ra_data_num

end subroutine read_sim_radar_dimension


subroutine read_sim_radar(ra_data_num,rv,rf,ir,jr,hr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_sim_radar
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
integer  :: ra_data_num
real     :: rv(ra_data_num),rf(ra_data_num),hr(ra_data_num)
integer  :: ir(ra_data_num),jr(ra_data_num)
integer  :: tem

tem=0

do i=1,ra_data_num
  read(1083,*) rv(i)
  read(1084,*) rf(i)
  read(1085,'(2I4,F16.3)') ir(i),jr(i),hr(i)
  
  if(rv(i) .ne. miss_data .and. rv(i) .ne. miss_data .and. tem == 0 ) then
   print*,'sample:rv,rf,ir,jr,hr'
   print*,rv(i),rf(i),ir(i),jr(i),hr(i)
   tem=1
  endif
  
enddo
 print*,'finished reading all sim obs'
end subroutine read_sim_radar