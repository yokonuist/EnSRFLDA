subroutine read_real_radar_info(iradar,xradar,yradar,hradar,lev_num,ra_data_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_real_radar_dimension
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
character(len=2)      :: lab_radar_num
character(len=5)      :: radar_name
real                  :: xradar,yradar
real                  :: hradar
integer               :: lev_num  
integer               :: iradar
integer               :: ra_data_num
integer               :: temp1,temp2,temp3

temp1=int(time_num/100     )-int(time_num/1000    )*10  
temp2=int(time_num/10      )-int(time_num/100     )*10  
temp3=int(time_num/1       )-int(time_num/10      )*10
lab_time_num = char(48+temp1)//char(48+temp2)//char(48+temp3) 

temp1=int(iradar/10      )-int(iradar/100     )*10  
temp2=int(iradar/1       )-int(iradar/10      )*10
lab_radar_num = char(48+temp1)//char(48+temp2)
                                     
print*,'the ',iradar,'th radar'   
   
length_obs_dir=len_trim(obs_file_dir)
open(1083,file=''//obs_file_dir(1:length_obs_dir)//'/obs_radar_'//lab_time_num//'_'//lab_radar_num//'.dat') 

read(1083,'(F8.3,F8.3)')  xradar,yradar
print*,'radar site:',xradar,yradar
read(1083,'(F8.3)') hradar
print*,'radar height:',hradar
read(1083,'(I4)')   lev_num 
print*,'level number:',lev_num  
read(1083,'(I10)')  ra_data_num   
write(*,*) 'THE NUMBER OF RADAR DATA:',ra_data_num                                     
                                     

                                     
end subroutine read_real_radar_info       



subroutine read_real_radar_data(ra_data_num,rv,rf,ir,jr,hr_rv,hr_rf,evl,azimuth) 
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
real     :: rv(ra_data_num),rf(ra_data_num),hr_rv(ra_data_num),hr_rf(ra_data_num)
integer  :: ir(ra_data_num),jr(ra_data_num)
real     :: evl(ra_data_num),azimuth(ra_data_num)
integer  :: tem   

tem=0

do i=1,ra_data_num
  read(1083,'(F12.3,I4,I4,F12.3,F12.3,F12.3)') rf(i),ir(i),jr(i),hr_rf(i),evl(i),azimuth(i)
  read(1083,'(F12.3,I4,I4,F12.3,F12.3,F12.3)') rv(i),ir(i),jr(i),hr_rv(i),evl(i),azimuth(i)
  
  if(rv(i) .ne. miss_data .and. rf(i) .ne. miss_data .and. tem == 0 ) then
   print*,'sample:rv,ir,jr,hr_rf,evl,azimuth'
   write(6,'(F12.3,I4,I4,F12.3,F12.3,F12.3)') rv(i),ir(i),jr(i),hr_rv(i),evl(i),azimuth(i)
   print*,'sample:rf,ir,jr,hr_rf,evl,azimuth'   
   write(6,'(F12.3,I4,I4,F12.3,F12.3,F12.3)') rf(i),ir(i),jr(i),hr_rf(i),evl(i),azimuth(i)
   tem=1
  endif
  
enddo
 print*,'finished reading all real radar obs'

close(1083)

end subroutine read_real_radar_data



subroutine read_real_snd_info(obs_snd_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_real_snd_info
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
integer               :: obs_snd_num
integer               :: temp1,temp2,temp3

temp1=int(time_num/100     )-int(time_num/1000    )*10  
temp2=int(time_num/10      )-int(time_num/100     )*10  
temp3=int(time_num/1       )-int(time_num/10      )*10
lab_time_num = char(48+temp1)//char(48+temp2)//char(48+temp3) 
   
length_obs_dir=len_trim(obs_file_dir)
open(1083,file=''//obs_file_dir(1:length_obs_dir)//'/snd/obs_'//lab_time_num//'.dat') 

read(1083,*)    obs_snd_num                             
                                     
                                     
end subroutine read_real_snd_info       

subroutine read_real_snd_data(usnd,vsnd,temperaturesnd,latsnd,lonsnd,hgtsnd,obs_snd_num) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine read_real_snd_data
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
include 'namelist.inc'
integer  :: obs_snd_num
real     :: usnd(obs_snd_num),vsnd(obs_snd_num),temperaturesnd(obs_snd_num),latsnd(obs_snd_num),  &
            lonsnd(obs_snd_num),hgtsnd(obs_snd_num)
integer  :: tem   

tem=0

do i=1,obs_snd_num
  read(1083,'(6f11.6)')  usnd(i),vsnd(i),temperaturesnd(i),latsnd(i),lonsnd(i),hgtsnd(i)
  
  if(tem == 0 ) then
   print*,'sample:usnd,vsnd,temperaturesnd,latsnd,lonsnd,hgtsnd'
   write(*,'(6f11.6)')  usnd(i),vsnd(i),temperaturesnd(i),latsnd(i),lonsnd(i),hgtsnd(i)
   tem=1
  endif
  
enddo
 print*,'finished reading all real snd obs'

close(1083)

end subroutine read_real_snd_data
