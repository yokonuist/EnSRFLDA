  program pert_sounding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!           subroutine pert_sounding
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  implicit none
  include 'pert_sounding.inc'
  
  
  integer,parameter :: n=46
  real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
  integer k
  real GASDEV
  integer seed
  
  open (5,file='pert_sounding.input')
  read (5,nml=pert_info)
  write(6,nml=pert_info)

  open(1001,file='./sounding/input_sounding',form='formatted',status='old')
  read(1001,*) ps, ts, qvs

  do k=1,n
    read(1001,*) h(k), th(k), qv(k), u(k), v(k)
  enddo
  close(1001)
  
  seed=-777
  
  do k=1,n
   if(enable_pert_u==1)  u(k) =u(k)+GASDEV(seed)*stdu
   if(enable_pert_v==1)  v(k) =v(k)+GASDEV(seed)*stdv
   if(enable_pert_t==1)  th(k)=th(k)+GASDEV(seed)*stdt
   if(enable_pert_q==1)  qv(k)=qv(k)+GASDEV(seed)*stdq
  enddo
  

  open(1001,file='./input_sounding',form='formatted')
  write(1001,'(F8.2,5X,F6.2,5X,F5.2)') ps, ts, qvs
  write(6,'(F8.2,5X,F6.2,5X,F5.2)') ps, ts, qvs
  do k=1,n
    write(1001,'(F8.2,5X,F6.2,5X,F5.2,5X,F6.2,5X,F6.2)') h(k), th(k), qv(k), u(k), v(k)
    write(6,'(F8.2,5X,F6.2,5X,F5.2,5X,F6.2,5X,F6.2)') h(k), th(k), qv(k), u(k), v(k)
  enddo
  close(1001)  


  end   program pert_sounding