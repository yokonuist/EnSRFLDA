  subroutine get_info_from_cdf(input_file,length_input,analysis_var_num,analysis_var_name,       &
                               rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode,              &
                               u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz                &
                               ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine get_info_from_cdf
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

  character (len=500)    :: input_file
  integer                :: i,j,k,ivtype, length,length_input
  
  character (len=80)     :: varnam
 ! integer                :: dimids(12)
  integer                :: dimids(13)
  integer                :: cdfid, rcode, id_var, idvar
!  integer                :: nDims, nAtts, dims(4)
  integer                :: nDims(13), nAtts, dims(4) 
 integer                :: dims3

  integer                :: istart(4), iend(4)
 
  integer                :: analysis_var_num
  integer                :: nn
  character(len=80)      :: analysis_var_name(analysis_var_num)  
  integer                :: rec_id_var(analysis_var_num)
  integer                :: file_cdfid
  !integer                :: rec_dims(analysis_var_num,3)
  integer                :: rec_dims(12,3)
 !integer,allocatable    :: rec_dims(:,:)
  integer                :: rec_dims2(2)
  integer                :: file_rcode  
  
  integer                :: nx,ny,nz
  real                   :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ph(nx,ny,nz),t(nx,ny,nz),qv(nx,ny,nz),   &
                            qr(nx,ny,nz),qi(nx,ny,nz),qs(nx,ny,nz),qgr(nx,ny,nz)
  real                   :: smois(nx,ny,4),tslb(nx,ny,4),tsk(nx,ny)
!ADD SMOIS & TSLB
  real,allocatable       :: data_r(:,:,:)
!  real,allocatable       :: data_tskr(:,:)

!--------------------finished var define-----------------------------------
!--------------------------------------------------------------------------  


! Open netCDF file 

   rcode = nf_open(input_file(1:length_input), NF_WRITE, cdfid )
   print*,"Attempting to open netCDF file with write access"

   length = max(1,index(input_file,' ')-1)
   if( rcode == 0) then
     write(6,*) ' '
   else
     write(6,*) ' error opening netcdf file ',input_file(1:length)
     stop
   end if

! start the data reading section
! each var will be read sequentially

  do nn = 1, analysis_var_num
 
  varnam=analysis_var_name(nn)
! Detect the var id according to varnam        
  rcode = nf_inq_varid ( cdfid, varnam, id_var )
       
  dims = 1
  rcode = nf_inq_var( cdfid, id_var, varnam, ivtype, nDims(nn), dimids, nAtts )
!print*,nn,cdfid, id_var, varnam, ivtype, nDims!, dimids             
! Get the dimensions of this field
    do i=1,nDims(nn)
      rcode = nf_inq_dimlen( cdfid, dimids(i), dims(i) )
    enddo
      istart        = 1
      istart(nDims(nn)) = 1
      iend          = 1
    do i = 1,nDims(nn)-1
      iend(i)     = dims(i)
 ! print*,nn,iend(i)
    enddo

    
! record the info of var in nc

  rec_id_var(nn)=id_var
  file_cdfid=cdfid
  if(nn < 13) then
   do i=1,nDims(nn)-1
  rec_dims(nn,i)=dims(i)
!  print*,nn,rec_dims(nn,i)
   enddo
  else if( nn == 13)then
   do i=1,nDims(nn)-1
   rec_dims2(i)=dims(i)
   enddo
  endif
  file_rcode=rcode    
 ! print*,nn,rcode,rec_id_var(nn),file_cdfid,rec_dims2(1),rec_dims2(2)
 ! print*,nn,rcode,rec_id_var(nn),file_cdfid,varnam
  enddo     ! for the var loop  
!stop  
  call fill_data(u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz,                              &
                 analysis_var_num,rec_id_var,file_cdfid,rec_dims,rec_dims2,nDims,file_rcode       &
                ) 
!print*,'get_info_cdf',tsk(1,1) 
return
end subroutine get_info_from_cdf
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
  subroutine get_basic_info_1(basic_var_num_1,basic_var_name_1,                    &
                              file_cdfid,file_rcode,nx,ny,nz,                      &
                              xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,        &
                              znw,znu,p_top,rdnw,rdn                               &
                             )  

  implicit none
  
  include 'netcdf.inc'
  
  integer                :: i,j,k,ivtype 
  character (len=80)     :: varnam
  integer                :: dimids(10)
  
  integer                :: cdfid, rcode, id_var, idvar
  integer                :: nDims, nAtts, dims(4)
  integer                :: dims3

  integer                :: istart(4), iend(4)
 
  integer                :: basic_var_num_1
  integer                :: nn
  character(len=80)      :: basic_var_name_1(basic_var_num_1)  
  integer                :: file_cdfid
  integer                :: file_rcode  
  
  integer                :: nx,ny,nz
  real                   :: xlon(nx,ny),xlat(nx,ny),ulon(nx,ny),ulat(nx,ny),vlon(nx,ny),vlat(nx,ny),phb(nx,ny,nz),hgt(nx,ny),  &
                            znw(nz),znu(nz),rdnw(nz),rdn(nz),dzs(4),zs(4)
  real                   :: p_top
  real,allocatable       :: data_r_2d(:,:)
  real,allocatable       :: data_r_3d(:,:,:)
  real,allocatable       :: data_r_1d(:)
  real                   :: data_r

!--------------------finished var define-----------------------------------
!--------------------------------------------------------------------------  

  do nn = 1, basic_var_num_1
    
  varnam=basic_var_name_1(nn)
! Detect the var id according to varnam        
  rcode = nf_inq_varid ( file_cdfid, varnam, id_var )
       
  dims = 1
  rcode = nf_inq_var( file_cdfid, id_var, varnam, ivtype, nDims, dimids, nAtts )             
! Get the dimensions of this field
    do i=1,ndims
      rcode = nf_inq_dimlen( file_cdfid, dimids(i), dims(i) )
    enddo
      istart        = 1
      istart(nDims) = 1
      iend          = 1
    do i = 1,nDims-1
      iend(i)     = dims(i)
    enddo
    
!  print*,'reading var:',varnam
!  print*,'id_var is',id_var
!  print*,'dimension is',dims
!  print*,'nDims',nDims

    
if(nn .ne. 7 .and. nn .ne. 9 .and. nn .lt. 10) then   
 if(allocated(data_r_2d)) deallocate(data_r_2d)   
 allocate(data_r_2d(dims(1),dims(2)))    
endif
if(nn == 7) then	
 if(allocated(data_r_3d)) deallocate(data_r_3d)   
 allocate(data_r_3d(dims(1),dims(2),dims(3)))
endif
if(nn == 10 .or. nn == 11 .or. nn == 12 .or. nn == 13 ) then	
 if(allocated(data_r_1d)) deallocate(data_r_1d)   
 allocate(data_r_1d(dims(1)))
endif

   if( nn == 1 ) then    !   for XLON

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      xlon(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(xlon,dims(1),dims(2),nx,ny) 

    endif      

   if( nn == 2 ) then    !   for XLAT

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      xlat(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(xlat,dims(1),dims(2),nx,ny) 

    endif   

   if( nn == 3 ) then    !   for ULON

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      ulon(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(ulon,dims(1),dims(2),nx,ny) 

    endif      

   if( nn == 4 ) then    !   for ULAT

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      ulat(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(ulat,dims(1),dims(2),nx,ny) 

    endif   
    
   if( nn == 5 ) then    !   for VLON

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      vlon(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(vlon,dims(1),dims(2),nx,ny) 

    endif      

   if( nn == 6 ) then    !   for VLAT

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      vlat(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(vlat,dims(1),dims(2),nx,ny) 

    endif   


   if( nn == 7 ) then    !   for PHB

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_3d,file_rcode)
    
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      phb(i,j,k)=data_r_3d(i,j,k)    
    enddo
    enddo
    enddo

    call edge_fill(phb,dims(1),dims(2),dims(3),nx,ny,nz) 

    endif
  
   if( nn == 8 ) then    !   for HGT

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      hgt(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(hgt,dims(1),dims(2),nx,ny) 

    endif    
    
   if( nn == 9 ) then    !   for P_TOP

    call ncvgt( file_cdfid,id_var,istart,iend,data_r,file_rcode)
    
      p_top=data_r   

    endif    
    
   if( nn == 10 ) then    !   for ZNW

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_1d,file_rcode)
    
    do i=1,dims(1)
      znw(i)=data_r_1d(i)    
    enddo


    endif                    
    
   if( nn == 11 ) then    !   for ZNU

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_1d,file_rcode)
    
    do i=1,dims(1)
      znu(i)=data_r_1d(i)    
    enddo
    
    znu(nz)=znu(nz-1)
    
    endif        
    
   if( nn == 12 ) then    !   for RDNW

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_1d,file_rcode)
    
    do i=1,dims(1)
      rdnw(i)=data_r_1d(i)    
    enddo
    
    rdnw(nz)=rdnw(nz-1)
    
    endif  
    
   if( nn == 13 ) then    !   for RDN

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_1d,file_rcode)
    
    do i=1,dims(1)
      rdn(i)=data_r_1d(i)    
    enddo
    
    rdn(nz)=rdn(nz-1)
    
    endif
 
if( nn == 14 ) then    !   for DZS

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_1d,file_rcode)
    
    do i=1,dims(1)
      dzs(i)=data_r_1d(i)    
    enddo
    
    !rdn(nz)=rdn(nz-1)
    
    endif
 if( nn == 15) then    !   for ZS

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_1d,file_rcode)
    
    do i=1,dims(1)
      zs(i)=data_r_1d(i)    
    enddo
    
    !rdn(nz)=rdn(nz-1)
    
    endif      

  enddo     ! for the var loop  
  

return
end subroutine get_basic_info_1
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
  subroutine get_basic_info_2(basic_var_num_2,basic_var_name_2,                    &
                              file_cdfid,file_rcode,nx,ny,                         &
                              mub,mu                                               &
                             )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine get_basic_info_2
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
  
  integer                :: i,j,k,ivtype 
  character (len=80)     :: varnam
  integer                :: dimids(10)
  
  integer                :: cdfid, rcode, id_var, idvar
  integer                :: nDims, nAtts, dims(4)
  integer                :: dims3

  integer                :: istart(4), iend(4)
 
  integer                :: basic_var_num_2
  integer                :: nn
  character(len=80)      :: basic_var_name_2(basic_var_num_2)  
  integer                :: file_cdfid
  integer                :: file_rcode  
  
  integer                :: nx,ny
  real                   :: mub(nx,ny),mu(nx,ny)
  real,allocatable       :: data_r_2d(:,:)


!--------------------finished var define-----------------------------------
!--------------------------------------------------------------------------  

  do nn = 1, basic_var_num_2
    
  varnam=basic_var_name_2(nn)
! Detect the var id according to varnam        
  rcode = nf_inq_varid ( file_cdfid, varnam, id_var )
       
  dims = 1
  rcode = nf_inq_var( file_cdfid, id_var, varnam, ivtype, nDims, dimids, nAtts )             
! Get the dimensions of this field
    do i=1,ndims
      rcode = nf_inq_dimlen( file_cdfid, dimids(i), dims(i) )
    enddo
      istart        = 1
      istart(nDims) = 1
      iend          = 1
    do i = 1,nDims-1
      iend(i)     = dims(i)
    enddo

 if(allocated(data_r_2d)) deallocate(data_r_2d)   
 allocate(data_r_2d(dims(1),dims(2)))   
     
   if( nn == 1 ) then    !   for MUB

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      mub(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(mub,dims(1),dims(2),nx,ny) 

    endif    
    
   if( nn == 2 ) then    !   for MU

    call ncvgt( file_cdfid,id_var,istart,iend,data_r_2d,file_rcode)
    
    do j=1,dims(2)
    do i=1,dims(1)
      mu(i,j)=data_r_2d(i,j)    
    enddo
    enddo

    call edge_fill_2d(mu,dims(1),dims(2),nx,ny) 

   endif

  enddo     ! for the var loop  
  
return
end subroutine get_basic_info_2
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine fill_data(u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz,                              &
                     var_num,rec_id_var,file_cdfid,rec_dims,rec_dims2,nDims,file_rcode                &
                    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine fill_data
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

integer                   :: nx,ny,nz
real                      :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ph(nx,ny,nz),t(nx,ny,nz),qv(nx,ny,nz), &
                             qr(nx,ny,nz),qi(nx,ny,nz),qs(nx,ny,nz),qgr(nx,ny,nz),smois(nx,ny,nz),tslb(nx,ny,nz),&
                             tsk(nx,ny)   
real,allocatable          :: data_r(:,:,:)
real,allocatable          :: data_rr(:,:)
integer                   :: i,j,k
integer                   :: nn,ierr
integer                   :: var_num
integer                   :: rec_id_var(var_num)
integer                   :: file_cdfid
integer                   :: rec_dims(var_num-1,3)
integer                   :: rec_dims2(2) 
integer                   :: file_rcode  
integer                   :: istart(4), iend(4),istartt(3),iendd(3),nDims(13)

! print*,rec_dims2(1),rec_dims2(2)
! print*,nDims(13)
do nn=1, var_num

 istart        = 1
 istart(4)     = 1
 iend          = 1

 istartt        = 1
 istartt(3)     = 1
 iendd          = 1

if(nn < 13) then
do j = 1,nDims(nn)-1
   iend(j)     = rec_dims(nn,j)
!print*,nn,iend(j),nDims(nn)
enddo
else if(nn == 13 ) then
do j = 1,nDims(nn)-1
   iendd(j)     = rec_dims2(j)
!print*,nn,iendd(j),nDims(nn)
 enddo  
endif 

!print*,'guoyk',nn,rec_dims(nn,1)     
 if(allocated(data_r)) deallocate(data_r)   
 allocate(data_r(rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3)))
 if(allocated(data_rr)) deallocate(data_rr)   
 allocate(data_rr(rec_dims2(1),rec_dims2(2)))

   if( nn == 1 ) then    !   for u
!print*,file_cdfid,rec_id_var(nn),istart,iend,file_rcode
    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
!print*,'guoyk',nn 
!print*,rec_dims(nn,3)
!print*,rec_dims(nn,2)
!print*,rec_dims(nn,1)   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      u(i,j,k)=data_r(i,j,k) 
!    print*,i,j,k,data_r(i,j,k) 
    enddo
!    print*,j
    enddo
!    print*,k
    enddo
!print*,'guoyk',nn
    call edge_fill(u,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz) 

    endif
!print*,'guoyk',nn
   if( nn == 2 ) then    !   for v

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      v(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(v,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif
    
   if( nn == 3 ) then    !   for w

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      w(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(w,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz) 
            
    endif 
    
   if( nn == 4 ) then    !   for ph

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      ph(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(ph,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif
    




   if( nn == 5 ) then    !   for t

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      t(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(t,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif 
    
   if( nn == 6 ) then    !   for qv

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      qv(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(qv,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif  
    
   if( nn == 7 ) then    !   for qr

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      qr(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(qr,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif  
    
   if( nn == 8 ) then    !   for qi

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      qi(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(qi,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif   
    
   if( nn == 9 ) then    !   for qs

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      qs(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(qs,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif             

   if( nn == 10 ) then    !   for qgr

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      qgr(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
    call edge_fill(qgr,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif    
!print*,'ADD IF FOR SMOIS & TSLB'
    if( nn == 11 ) then    !   for smois

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      SMOIS(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
!    call edge_fill(smois,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif

if( nn == 12 ) then    !   for tslb

    call ncvgt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      tslb(i,j,k)=data_r(i,j,k)    
    enddo
    enddo
    enddo
        
!    call edge_fill(tslb,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3),nx,ny,nz)
            
    endif

if( nn == 13 ) then    !   for tsk

    call ncvgt( file_cdfid,rec_id_var(nn),istartt,iendd,data_rr,file_rcode)
!    print*,rec_id_var(nn),file_cdfid  
!    do k=1,rec_dims(nn,3)
!  ierr=nf_get_var_real(file_cdfid,rec_id_var(nn),data_rr)
    do j=1,rec_dims2(2)
    do i=1,rec_dims2(1)
      tsk(i,j)=data_rr(i,j)    

    enddo
    enddo
  ! enddo
!print*,iendd
!print*,tsk (1,1),data_rr(1,1),ierr       
!    call edge_fill_2d(tsk,rec_dims2(1),rec_dims2(2),nx,ny)
            
    endif
!print*,'guoyk',nn
enddo   
return
end subroutine fill_data   
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine edge_fill(var,dim1,dim2,dim3,nx,ny,nz) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine edge_fill
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        
implicit none

integer  :: nx,ny,nz,i,j,k
integer  :: dim1,dim2,dim3
real     :: var(nx,ny,nz)

if( dim1 == nx .and. dim2 < ny .and. dim3 < nz ) then
  
  do k=1,dim3
  do i=1,nx
   var(i,ny,k)=var(i,ny-1,k)
  enddo
  enddo

  do j=1,ny
  do i=1,nx
   var(i,j,nz)=var(i,j,nz-1)
  enddo
  enddo
  
endif

if( dim1 < nx .and. dim2 == ny .and. dim3 < nz ) then

  do k=1,dim3
  do j=1,ny
   var(nx,j,k)=var(nx-1,j,k)
  enddo
  enddo

  do j=1,ny
  do i=1,nx
   var(i,j,nz)=var(i,j,nz-1)
  enddo
  enddo
  
endif

if( dim1 < nx .and. dim2 < ny .and. dim3 == nz ) then

  do k=1,nz
  do i=1,dim2
   var(i,ny,k)=var(i,ny-1,k)
  enddo
  enddo
  
  do k=1,nz
  do j=1,ny
   var(nx,j,k)=var(nx-1,j,k)
  enddo
  enddo

endif

if( dim1 < nx .and. dim2 < ny .and. dim3 < nz ) then
  do k=1,dim3
  do j=1,dim2
   var(nx,j,k)=var(nx-1,j,k)
  enddo
  enddo
  
  do k=1,dim3
  do i=1,nx
   var(i,ny,k)=var(i,ny-1,k)
  enddo
  enddo
  
  do i=1,nx
  do j=1,ny
   var(i,j,nz)=var(i,j,nz-1)
  enddo
  enddo
    
endif
return
end subroutine edge_fill
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine edge_fill_2d(var,dim1,dim2,nx,ny) 

        
implicit none

integer  :: nx,ny,i,j
integer  :: dim1,dim2
real     :: var(nx,ny)

if( dim1 == nx .and. dim2 < ny) then

  do i=1,nx
   var(i,ny)=var(i,ny-1)
  enddo

endif

if( dim1 < nx .and. dim2 == ny) then

  do j=1,ny
   var(nx,j)=var(nx-1,j)
  enddo

endif

if( dim1 < nx .and. dim2 < ny) then

  do j=1,ny-1
   var(nx,j)=var(nx-1,j)
  enddo
  
  do i=1,nx
   var(i,ny)=var(i,ny-1)
  enddo

endif

return

end subroutine edge_fill_2d
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine output_data(u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk,nx,ny,nz,                            &
                       analysis_var_num,rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode           &
                      )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine output_data
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
                      
integer                   :: nx,ny,nz
integer                   :: xtype,ndims,dimids,natts
character*10              :: vname
real                      :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ph(nx,ny,nz),t(nx,ny,nz),qv(nx,ny,nz), &
                             qr(nx,ny,nz),qi(nx,ny,nz),qs(nx,ny,nz),qgr(nx,ny,nz),smois(nx,ny,4),tslb(nx,ny,4)
real                      :: tsk(nx,ny)
real,allocatable          :: data_r(:,:,:)
real,allocatable          :: data_rr(:,:)
integer                   :: i,j,k,ierr
integer                   :: nn
integer                   :: analysis_var_num
integer                   :: rec_id_var(analysis_var_num)
integer                   :: file_cdfid
integer                   :: rec_dims(12,3)
integer                   :: rec_dims2(2)
integer                   :: file_rcode
integer                   :: istart(4), iend(4)
integer                   :: istartt(3),iendd(3)


 
!stop
do nn=1, analysis_var_num

if(nn < 13) then
 istart        = 1
 istart(4)     = 1
 iend          = 1
 do j = 1,3
   iend(j)     = rec_dims(nn,j)
!print*,iend(j)
 enddo  
!else if(nn == 11 .or. nn == 12) then
! istart        = 1
! istart(4)     = 1
! iend          = 1
! do j = 1,3
!   iend(j)     = rec_dims(nn,j)
! print*,iend(j)
 !enddo 
!   iend(1)     = nx-1
!   iend(2)     = ny-1
!   iend(3)     = 4  !soil layer
else if(nn == 13)then
 istartt        = 1
 istartt(3)     = 1
 iendd          = 1
do j = 1,2
   iendd(j)     = rec_dims2(j)
!print*,iendd(j),j
 enddo
endif

 if(allocated(data_r)) deallocate(data_r) 
 allocate(data_r(rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3)))
 if(allocated(data_rr)) deallocate(data_rr) 
 allocate(data_rr(rec_dims2(1),rec_dims2(2)))

!print*,'ok ?'
   if( nn == 1 ) then    !   for u
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = u(i,j,k)  
    enddo
    enddo
    enddo 
!  print*,rec_dims(nn,1),rec_dims(nn,2),rec_dims(nn,3)
!  print*,file_cdfid,rec_id_var(nn),istart,iend,file_rcode
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)  
 
   endif
!stop
!  print*,'ok ?'  
   if( nn == 2 ) then    !   for v
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = v(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode) 
    
   endif
   
   if( nn == 3 ) then    !   for w
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = w(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
    
   endif
   if( nn == 4 ) then    !   for ph
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = ph(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
    
   endif
   
   if( nn == 5 ) then    !   for t
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = t(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)   
    
   endif   
   
   if( nn == 6 ) then    !   for qv
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = qv(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)  
    
   endif
   
   if( nn == 7 ) then    !   for qr
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = qr(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)  
    
   endif
   
   if( nn == 8 ) then    !   for qi
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = qi(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)  
    
   endif
   
   if( nn == 9 ) then    !   for qs
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = qs(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)  
    
   endif   
    
   if( nn == 10 ) then    !   for qgr
   
    do k=1,rec_dims(nn,3)
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = qgr(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)  
    
   endif   

    if( nn == 11 ) then    !   for smois
   
    do k=1,4
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = smois(i,j,k)  
    enddo
    enddo
    enddo
!print*, file_cdfid,rec_id_var(nn),istart,iend,file_rcode
! stop
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
    
   endif


    if( nn == 12 ) then    !   for tslb
   
    do k=1,4
    do j=1,rec_dims(nn,2)
    do i=1,rec_dims(nn,1)
      data_r(i,j,k) = tslb(i,j,k)  
    enddo
    enddo
    enddo 
 
   call ncvpt( file_cdfid,rec_id_var(nn),istart,iend,data_r,file_rcode)
    
   endif

    if( nn == 13 ) then    !   for tsk
!    ierr = nf_inq_var(file_cdfid,rec_id_var(nn),vname,xtype,ndims,dimids,natts)
!    print*,ierr,file_cdfid,rec_id_var(nn),xtype,ndims,dimids,vname  
!    do k=1,rec_dims(nn,3)
      do j=1,rec_dims2(2)
      do i=1,rec_dims2(1)
!    do j=1,ny-1
!    do i=1,nx-1
     data_rr(i,j) = tsk(i,j)  
    enddo
    enddo
!    enddo 
! print*,data_rr(1,1)
   call ncvpt( file_cdfid,rec_id_var(nn),istartt,iendd,data_rr,file_rcode)
    
   endif
!print*,nn       
enddo
!stop
return
end subroutine output_data

subroutine read_ens_file(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                    &
                         rec_dims,rec_dims2,file_rcode,rec_cdfid,rec_rcode,                           &
                         ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,   & 
                         mub,mu,basic_var_num_2,basic_var_name_2,                                     & 
                         lab_domain,lab_time_num                                                      &                    
                         )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine read_ens_file
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none  
include 'namelist.inc'

character (len=500)   :: input_file_name
integer               :: length_input,length_file_dir,length_ens_head     
integer               :: analysis_var_num
integer               :: basic_var_num_2
integer               :: i,j,k,ii,jj,kk
integer               :: istart(4), iend(4)  !
integer               :: istartt(3),iendd(3) !
character (len=4)     :: lab_ens
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
  
integer               :: rec_id_var(analysis_var_num)
integer               :: rec_cdfid(ens_num)
!integer               :: rec_dims(analysis_var_num,3)
integer               :: rec_dims(12,3)
integer               :: rec_dims2(2)
integer               :: rec_rcode(ens_num)
character(len=80)     :: analysis_var_name(analysis_var_num)
character(len=80)     :: basic_var_name_2(basic_var_num_2)
integer               :: file_cdfid,file_rcode


real                  :: ensu(nx,ny,nz,ens_num)
real                  :: ensv(nx,ny,nz,ens_num)
real                  :: ensw(nx,ny,nz,ens_num)
real                  :: ensph(nx,ny,nz,ens_num)
real                  :: enst(nx,ny,nz,ens_num)
real                  :: ensqv(nx,ny,nz,ens_num)
real                  :: ensqr(nx,ny,nz,ens_num)
real                  :: ensqi(nx,ny,nz,ens_num)
real                  :: ensqs(nx,ny,nz,ens_num)
real                  :: ensqgr(nx,ny,nz,ens_num)
 !!
real                  :: mub(nx,ny,ens_num)
real                  :: mu (nx,ny,ens_num)
!!!
real                  :: enssmois(nx,ny,4,ens_num)
real                  :: enstslb(nx,ny,4,ens_num)

real                  :: enstsk(nx,ny,ens_num)
real                  :: tsk(nx,ny)

integer               :: temp1,temp2,temp3



  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)
  
  temp1=int(time_num/10      )-int(time_num/100     )*10  
  temp2=int(time_num/1       )-int(time_num/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)  

! enstsk = 0
 
  do i=1, ens_num    ! for reading ens file loop
  
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
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                                  &
                          rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode,                                              &
                          ensu(:,:,:,i),ensv(:,:,:,i),ensw(:,:,:,i),ensph(:,:,:,i),enst(:,:,:,i),ensqv(:,:,:,i),            &
                          ensqr(:,:,:,i),ensqi(:,:,:,i),ensqs(:,:,:,i),ensqgr(:,:,:,i),enssmois(:,:,:,i),enstslb(:,:,:,i),  &
                          enstsk(:,:,i),nx,ny,nz                             &
                          ) 


  call get_basic_info_2(basic_var_num_2,basic_var_name_2,                    &
                        file_cdfid,file_rcode,nx,ny,                         &
                        mub(:,:,i),mu(:,:,i)                                 &
                        )       
                                              
  rec_cdfid(i)=file_cdfid
  rec_rcode(i)=file_rcode



  enddo   ! for reading ens file loop

return                
end subroutine read_ens_file 

subroutine read_true_file(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                   &
                         rec_dims,file_rcode,                                                         &
                         utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,          & 
                         mubtrue,mutrue,basic_var_num_2,basic_var_name_2,                             &
                         xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,znw,znu,p_top,rdnw,rdn,                &
                         basic_var_num_1,basic_var_name_1,                                            & 
                         lab_domain,lab_time_num                                                      &            
                         )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine read_true_file
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

include 'namelist.inc'

character (len=500)   :: input_file_name
integer               :: length_input,length_file_dir,length_ens_head     
integer               :: analysis_var_num
integer               :: basic_var_num_1
integer               :: basic_var_num_2
integer               :: i,j,k,ii,jj,kk
integer               :: istart(4), iend(4)
character (len=4)     :: lab_ens
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
  

integer               :: rec_id_var(analysis_var_num)
integer               :: rec_cdfid(ens_num)
integer               :: rec_dims(analysis_var_num,3)
integer               :: rec_rcode(ens_num)
character(len=80)     :: analysis_var_name(analysis_var_num)
character(len=80)     :: basic_var_name_1(basic_var_num_1)
character(len=80)     :: basic_var_name_2(basic_var_num_2)
integer               :: file_cdfid,file_rcode

real                  :: utrue(nx,ny,nz)
real                  :: vtrue(nx,ny,nz)
real                  :: wtrue(nx,ny,nz)
real                  :: phtrue(nx,ny,nz)
real                  :: ttrue(nx,ny,nz)
real                  :: qvtrue(nx,ny,nz)
real                  :: qrtrue(nx,ny,nz)
real                  :: qitrue(nx,ny,nz)
real                  :: qstrue(nx,ny,nz)
real                  :: qgrtrue(nx,ny,nz)
real                  :: mubtrue(nx,ny) 
real                  :: mutrue (nx,ny)  

real                  :: xlon(nx,ny)
real                  :: xlat(nx,ny)
real                  :: ulon(nx,ny)
real                  :: ulat(nx,ny)
real                  :: vlon(nx,ny)
real                  :: vlat(nx,ny) 
real                  :: phb(nx,ny,nz)
real                  :: hgt(nx,ny) 
real                  :: mub(nx,ny,ens_num) 
real                  :: mu (nx,ny,ens_num) 
real                  :: znw   (nz) 
real                  :: znu   (nz) 
real                  :: rdnw  (nz) 
real                  :: rdn   (nz)   
  
integer               :: temp1,temp2,temp3

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
                          utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,nx,ny,nz     &
                          ) 
                          
  call get_basic_info_1(basic_var_num_1,basic_var_name_1,                    &
                        file_cdfid,file_rcode,nx,ny,nz,                      &
                        xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,        &
                        znw,znu,p_top,rdnw,rdn                               &
                       )  
                       
  call get_basic_info_2(basic_var_num_2,basic_var_name_2,                    &
                        file_cdfid,file_rcode,nx,ny,                         &
                        mubtrue,mutrue                                       &
                        )                           
                                                 
  call ncclos(file_cdfid,file_rcode)      
  
end subroutine read_true_file       

subroutine read_mean_file(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                                    &
                         rec_dims,rec_dims2,file_rcode,                                                                &
                         umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,& 
                         xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,znw,znu,p_top,rdnw,rdn,                          &
                         basic_var_num_1,basic_var_name_1,                                                             & 
                         lab_domain,lab_time_num                                                                       &            
                         )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine read_mean_file
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

include 'namelist.inc'

character (len=500)   :: input_file_name
integer               :: length_input,length_file_dir,length_ens_head     
integer               :: analysis_var_num
integer               :: basic_var_num_1

integer               :: i,j,k,ii,jj,kk
integer               :: istart(4), iend(4)
character (len=4)     :: lab_ens
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
  
integer               :: rec_id_var(analysis_var_num)
integer               :: rec_cdfid(ens_num)
!integer               :: rec_dims(analysis_var_num,3)
integer               :: rec_dims(12,3)
integer               :: rec_dims2(2)
integer               :: rec_rcode(ens_num)
character(len=80)     :: analysis_var_name(analysis_var_num)
character(len=80)     :: basic_var_name_1(basic_var_num_1)

integer               :: file_cdfid,file_rcode

real                  :: umean(nx,ny,nz)
real                  :: vmean(nx,ny,nz)
real                  :: wmean(nx,ny,nz)
real                  :: phmean(nx,ny,nz)
real                  :: tmean(nx,ny,nz)
real                  :: qvmean(nx,ny,nz)
real                  :: qrmean(nx,ny,nz)
real                  :: qimean(nx,ny,nz)
real                  :: qsmean(nx,ny,nz)
real                  :: qgrmean(nx,ny,nz)
real                  :: smoismean(nx,ny,nz)
real                  :: tslbmean(nx,ny,nz)  !!
real                  :: tskmean(nx,ny) !! TSK

real                  :: xlon(nx,ny)
real                  :: xlat(nx,ny)
real                  :: ulon(nx,ny)
real                  :: ulat(nx,ny)
real                  :: vlon(nx,ny)
real                  :: vlat(nx,ny) 
real                  :: phb(nx,ny,nz)
!!!!for soil depth!!!!
real                  :: dzs(4)
real                  :: zs(4)
real                  :: hgt(nx,ny) 
real                  :: mub(nx,ny,ens_num) 
real                  :: mu (nx,ny,ens_num) 
real                  :: znw   (nz) 
real                  :: znu   (nz) 
real                  :: rdnw  (nz) 
real                  :: rdn   (nz)   
  
integer               :: temp1,temp2,temp3

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
  
  
! Now read the file
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                              &
                          rec_id_var,file_cdfid,rec_dims,rec_dims2,file_rcode,                                          &
                          umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,&
                          nx,ny,nz     &
                          ) 
                          
  call get_basic_info_1(basic_var_num_1,basic_var_name_1,                    &
                        file_cdfid,file_rcode,nx,ny,nz,                      &
                        xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,               &
                        znw,znu,p_top,rdnw,rdn                               &
                       )  
                                              
                                                 
  call ncclos(file_cdfid,file_rcode)      
  
end subroutine read_mean_file           


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_ens_file2(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                    &
                          rec_dims,file_rcode,rec_cdfid,rec_rcode,                                     &
                          ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,                    & 
                          mub,mu,basic_var_num_2,basic_var_name_2,                                     & 
                          lab_domain,lab_time_num                                                      &                    
                          )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine read_ens_file
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
include 'namelist.inc'

character (len=500)   :: input_file_name
integer               :: length_input,length_file_dir,length_ens_head     
integer               :: analysis_var_num
integer               :: basic_var_num_2
integer               :: i,j,k,ii,jj,kk
integer               :: istart(4), iend(4)
character (len=4)     :: lab_ens
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
  
integer               :: rec_id_var(analysis_var_num)
integer               :: rec_cdfid(ens_num)
integer               :: rec_dims(analysis_var_num,3)
integer               :: rec_rcode(ens_num)
character(len=80)     :: analysis_var_name(analysis_var_num)
character(len=80)     :: basic_var_name_2(basic_var_num_2)
integer               :: file_cdfid,file_rcode


real                  :: ensu(nx,ny,nz,ens_num)
real                  :: ensv(nx,ny,nz,ens_num)
real                  :: ensw(nx,ny,nz,ens_num)
real                  :: ensph(nx,ny,nz,ens_num)
real                  :: enst(nx,ny,nz,ens_num)
real                  :: ensqv(nx,ny,nz,ens_num)
real                  :: ensqr(nx,ny,nz,ens_num)
real                  :: ensqi(nx,ny,nz,ens_num)
real                  :: ensqs(nx,ny,nz,ens_num)
real                  :: ensqgr(nx,ny,nz,ens_num)
real                  :: mub(nx,ny,ens_num)
real                  :: mu (nx,ny,ens_num)

integer               :: temp1,temp2,temp3


  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)  

  temp1=int((time_num-1)/10      )-int((time_num-1)/100     )*10  
  temp2=int((time_num-1)/1       )-int((time_num-1)/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)   
  
  do i=1, ens_num    ! for reading ens file loop
  
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
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                                  &
                          rec_id_var,file_cdfid,rec_dims,file_rcode,                                                        &
                          ensu(:,:,:,i),ensv(:,:,:,i),ensw(:,:,:,i),ensph(:,:,:,i),enst(:,:,:,i),ensqv(:,:,:,i),            &
                          ensqr(:,:,:,i),ensqi(:,:,:,i),ensqs(:,:,:,i),ensqgr(:,:,:,i),nx,ny,nz                             &
                          ) 
            
  call get_basic_info_2(basic_var_num_2,basic_var_name_2,                    &
                        file_cdfid,file_rcode,nx,ny,                         &
                        mub(:,:,i),mu(:,:,i)                                 &
                        )       
                                            
  rec_cdfid(i)=file_cdfid
  rec_rcode(i)=file_rcode

  enddo   ! for reading ens file loop
                          
end subroutine read_ens_file2 

subroutine read_true_file2(analysis_var_num,analysis_var_name,rec_id_var,file_cdfid,                   &
                           rec_dims,file_rcode,                                                         &
                           utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,          & 
                           mubtrue,mutrue,basic_var_num_2,basic_var_name_2,                             &
                           xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,znw,znu,p_top,rdnw,rdn,                &
                           basic_var_num_1,basic_var_name_1,                                            & 
                           lab_domain,lab_time_num                                                      &            
                          )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine read_true_file
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

include 'namelist.inc'

character (len=500)   :: input_file_name
integer               :: length_input,length_file_dir,length_ens_head     
integer               :: analysis_var_num
integer               :: basic_var_num_1
integer               :: basic_var_num_2
integer               :: i,j,k,ii,jj,kk
integer               :: istart(4), iend(4)
character (len=4)     :: lab_ens
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
  
integer               :: rec_id_var(analysis_var_num)
integer               :: rec_cdfid(ens_num)
integer               :: rec_dims(analysis_var_num,3)
integer               :: rec_rcode(ens_num)
character(len=80)     :: analysis_var_name(analysis_var_num)
character(len=80)     :: basic_var_name_1(basic_var_num_1)
character(len=80)     :: basic_var_name_2(basic_var_num_2)
integer               :: file_cdfid,file_rcode

real                  :: utrue(nx,ny,nz)
real                  :: vtrue(nx,ny,nz)
real                  :: wtrue(nx,ny,nz)
real                  :: phtrue(nx,ny,nz)
real                  :: ttrue(nx,ny,nz)
real                  :: qvtrue(nx,ny,nz)
real                  :: qrtrue(nx,ny,nz)
real                  :: qitrue(nx,ny,nz)
real                  :: qstrue(nx,ny,nz)
real                  :: qgrtrue(nx,ny,nz)
real                  :: mubtrue(nx,ny) 
real                  :: mutrue (nx,ny)  

real                  :: xlon(nx,ny)
real                  :: xlat(nx,ny)
real                  :: ulon(nx,ny)
real                  :: ulat(nx,ny)
real                  :: vlon(nx,ny)
real                  :: vlat(nx,ny) 
real                  :: phb(nx,ny,nz)
real                  :: dzs(4)
real                  :: zs(4)
real                  :: hgt(nx,ny) 
real                  :: mub(nx,ny,ens_num) 
real                  :: mu (nx,ny,ens_num) 
real                  :: znw   (nz) 
real                  :: znu   (nz) 
real                  :: rdnw  (nz) 
real                  :: rdn   (nz)   
  
integer               :: temp1,temp2,temp3

  temp1=int(domain_num/10      )-int(domain_num/100     )*10  
  temp2=int(domain_num/1       )-int(domain_num/10      )*10
  lab_domain = char(48+temp1)//char(48+temp2)
  
  temp1=int((time_num-1)/10      )-int((time_num-1)/100     )*10  
  temp2=int((time_num-1)/1       )-int((time_num-1)/10      )*10
  lab_time_num = char(48+temp1)//char(48+temp2)   
 
  length_file_dir =len_trim(obs_file_dir)
  length_ens_head=len_trim(ens_file_head)
  
  input_file_name=''//obs_file_dir(1:length_file_dir)//'/'//ens_file_head(1:length_ens_head)//'_d'//lab_domain//'_'//lab_time_num//''
  
  print*,"INPUT FILE IS: ",trim(input_file_name)
  length_input = len_trim(input_file_name)
  
  
! Now read the file
  call get_info_from_cdf (input_file_name,length_input,analysis_var_num,analysis_var_name,                &
                          rec_id_var,file_cdfid,rec_dims,file_rcode,                                      &
                          utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,nx,ny,nz     &
                          ) 
                          
  call get_basic_info_1(basic_var_num_1,basic_var_name_1,                    &
                        file_cdfid,file_rcode,nx,ny,nz,                      &
                        xlon,xlat,ulon,ulat,vlon,vlat,phb,dzs,zs,hgt,               &
                        znw,znu,p_top,rdnw,rdn                               &
                       )  
                       
  call get_basic_info_2(basic_var_num_2,basic_var_name_2,                    &
                        file_cdfid,file_rcode,nx,ny,                         &
                        mubtrue,mutrue                                       &
                        )                           
                                                 
  call ncclos(file_cdfid,file_rcode)      
  
end subroutine read_true_file2             
 
