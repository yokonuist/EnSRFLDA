
subroutine addpert(ranpert,dim1,dim2,dim3,seed_ens,rec_id_var)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine addpert
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer   :: dim1,dim2,dim3
real      :: ranpert(dim1,dim2,dim3)
real      :: tem(dim1,dim2,dim3)
real      :: std
integer   :: seed_ens,i,j,k
integer   :: tem1
integer   :: ens_num

ranpert=0
tem1=-777-seed_ens 
!print*,rec_id_var
! if (rec_id_var < 60)then
   do k=2,dim3-1
   do j=5,dim2-5 
   do i=5,dim1-5 
    if(seed_ens < ens_num/2)  then
     ranpert(i,j,k)=GASDEV(tem1)
    else
     ranpert(i,j,k)=GASDEV(tem1)
    endif
   enddo
   enddo
   enddo
   
tem=ranpert
   
   do k=2,dim3-1
   do j=2,dim2-1 
   do i=2,dim1-1

    ranpert(i,j,k)= 0.2*(tem(i-1,j,k)+tem(i+1,j,k)+tem(i,j-1,k)+tem(i,j+1,k)) &
                   +0.8*tem(i,j,k)
    
   enddo
   enddo
   enddo
!endif


!!!!!!for smois & tslb pert
! if(rec_id_var == 63 .or. rec_id_var == 62 )then
!print*,'ok',rec_id_var
!   do k=1,4
!   do j=5,dim2-5 
!   do i=5,dim1-5 
!    if(seed_ens < ens_num/2)  then
!     ranpert(i,j,k)=GASDEV(tem1)
!    else
!     ranpert(i,j,k)=GASDEV(tem1)
!    endif
!   enddo
!   enddo
!   enddo
   
!tem=ranpert
   
!   do k=1,4
!   do j=2,dim2-1 
!   do i=2,dim1-1

!    ranpert(i,j,k)= 0.2*(tem(i-1,j,k)+tem(i+1,j,k)+tem(i,j-1,k)+tem(i,j+1,k)) &
!                   +0.8*tem(i,j,k)
    
!   enddo
!   enddo
!   enddo

!endif      
   
end subroutine addpert

subroutine addpertsoil(ranpertsoil,dim1,dim2,dim3,seed_ens)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine addpert
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer   :: dim1,dim2,dim3
real      :: ranpertsoil(dim1,dim2,dim3)
real      :: tem(dim1,dim2,dim3)
real      :: std
integer   :: seed_ens,i,j,k
integer   :: tem1
integer   :: ens_num

ranpert=0
tem1=-777-seed_ens 
!print*,rec_id_var
! if (rec_id_var < 60)then
!   do k=1,dim3-1
   do k=1,4
   do j=2,dim2-5 
   do i=5,dim1-5 
    if(seed_ens < ens_num/2)  then
     ranpertsoil(i,j,k)=GASDEV(tem1)   !0~1 distribution
    else
     ranpertsoil(i,j,k)=GASDEV(tem1)
    endif
   enddo
   enddo
   enddo
   
tem=ranpertsoil
   
  ! do k=1,dim3-1
   do k=1,4
   do j=2,dim2-1 
   do i=2,dim1-1

    ranpertsoil(i,j,k)= 0.2*(tem(i-1,j,k)+tem(i+1,j,k)+tem(i,j-1,k)+tem(i,j+1,k)) &
                   +0.8*tem(i,j,k)
    if(ranpertsoil(i,j,k) > 1)then   !more strict constraint
     ranpertsoil(i,j,k)=1
    endif
    if(ranpertsoil(i,j,k) <-1)then
      ranpertsoil(i,j,k)=-1
    endif
   enddo
   enddo
   enddo
!endif
!endif      
   
end subroutine addpertsoil

subroutine addpertsoil_2d(ranpertsoil_2d,dim1,dim2,seed_ens)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine addpert
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer   :: dim1,dim2
real      :: ranpertsoil_2d(dim1,dim2)
real      :: tem(dim1,dim2)
real      :: std
integer   :: seed_ens,i,j,k
integer   :: tem1
integer   :: ens_num

ranpert=0
tem1=-777-seed_ens 
!print*,rec_id_var
! if (rec_id_var < 60)then
!   do k=1,dim3-1
!   do k=1,4
   do j=2,dim2-5 
   do i=5,dim1-5 
    if(seed_ens .lt. ens_num/2 .or. seed_ens == ens_num/2 )  then
     ranpertsoil_2d(i,j)=GASDEV(tem1)   !0~1 distribution
    else
     ranpertsoil_2d(i,j)=GASDEV(tem1)
    endif
   enddo
   enddo
!   enddo
   
tem=ranpertsoil_2d
   
  ! do k=1,dim3-1
!   do k=1,4
   do j=2,dim2-1 
   do i=2,dim1-1

    ranpertsoil_2d(i,j)= 0.01*(tem(i-1,j)+tem(i+1,j)+tem(i,j-1)+tem(i,j+1)) &
                   +0.99*tem(i,j)
!    if(ranpertsoil(i,j,k) > 1)then   !more strict constraint
!     ranpertsoil(i,j,k)=1
!    endif
!    if(ranpertsoil(i,j,k) <-1)then
!      ranpertsoil(i,j,k)=-1
!    endif
   enddo
   enddo
!  enddo
!endif
!endif      
   
end subroutine addpertsoil_2d


subroutine addpert_local(ranpert,dim1,dim2,dim3,seed_ens,pert_cx,pert_cy,pert_radius,dx,dy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine addpert_local
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
integer   :: dim1,dim2,dim3
real      :: var(dim1,dim2,dim3)
real      :: ranpert(dim1,dim2,dim3)
real      :: tem(dim1,dim2,dim3)
real      :: std
integer   :: seed_ens,i,j,k
integer   :: tem1
integer   :: pert_cx,pert_cy 
real      :: pert_radius
real      :: dx,dy
real      :: temp1,temp2,temp3  
integer   :: ens_num

ranpert=0
tem1=-777-seed_ens   
 
   do k=2,dim3-1   
   do j=1,dim2 
   do i=1,dim1
   
    temp1=abs(i-pert_cx)*dx
    temp2=abs(j-pert_cy)*dy
    temp3=sqrt(temp1**2+temp2**2)
    
    if(temp3 <= pert_radius) then
   
      if(seed_ens < ens_num/2)  then
       ranpert(i,j,k)=GASDEV(tem1)
      else
       ranpert(i,j,k)=GASDEV(tem1)
      endif
    endif
    
   enddo
   enddo   
   enddo
     
tem=ranpert
   
   do k=2,dim3-1
   do j=2,dim2-1 
   do i=2,dim1-1

    ranpert(i,j,k)= 0.2*(tem(i-1,j,k)+tem(i+1,j,k)+tem(i,j-1,k)+tem(i,j+1,k)) &
                   +0.8*tem(i,j,k)
    
   enddo
   enddo
   enddo     
   
end subroutine addpert_local



subroutine addpert2(ranpert,dim1,dim2,dim3,seed_ens,scalar)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!            subroutine addpert
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.04.02   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer   :: dim1,dim2,dim3
real      :: ranpert(dim1,dim2,dim3)
real      :: tem(dim1,dim2,dim3)
real      :: std
integer   :: seed_ens,i,j,k
integer   :: tem1
integer   :: ens_num
integer   :: scalar


ranpert=0
tem1=-777-seed_ens 
 
   do k=2,dim3-1
   do j=5,dim2-5,scalar 
   do i=5,dim1-5,scalar 
     ranpert(i,j,k)=GASDEV(tem1)

   enddo
   enddo
   enddo
   
   
end subroutine addpert2


subroutine addpert_with_physical_constrain(xlon,xlat,ulon,ulat,vlon,vlat,phb,hgt,                                  &
                                           znw,znu,p_top,rdnw,rdn,                                                 &
                                           mub,mu,								   &
                                           ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,&
                                           ranpert,ranpert2,                                                       &
                                           nx,ny,nz,ens_num,dx,dy,                                                 &
                                           stdphi_init,stdpsi_init,stdqv_init,stdt_init,                           &
                                           upert,vpert,wpert,phpert,tpert,qvpert,smoispert,tslbpert                &
                                           )

  integer               :: nx,ny,nz,ens_num
  real                  :: dx,dy
  real                  :: ranpert(nx,ny,nz,ens_num)
  real                  :: ranpert2(nx,ny,nz,ens_num)                                         
  real,dimension(nx,ny) :: xlon,xlat,ulon,ulat,vlon,vlat
  real                  :: mub(nx,ny,ens_num),mu(nx,ny,ens_num),rdnw(nz),rdn(nz),znu(nz),znw(nz)
  real                  :: p_top
  real                  :: phb(nx,ny,nz)
  real                  :: ensu  (nx,ny,nz,ens_num)
  real                  :: ensv  (nx,ny,nz,ens_num)
  real                  :: ensw  (nx,ny,nz,ens_num)
  real                  :: ensph (nx,ny,nz,ens_num)
  real                  :: enst  (nx,ny,nz,ens_num)
  real                  :: ensqv (nx,ny,nz,ens_num)  
  real                  :: ensqr (nx,ny,nz,ens_num) 
  real                  :: ensqi (nx,ny,nz,ens_num)
  real                  :: ensqs (nx,ny,nz,ens_num) 
  real                  :: ensqgr(nx,ny,nz,ens_num)
  real                  :: enssmois (nx,ny,4,ens_num)
  real                  :: enstslb (nx,ny,4,ens_num) 
  
                                               
  real                  :: stdqv_init,stdt_init
                                             
  real                  :: phi(nx,ny,nz,ens_num)
  real                  :: psi(nx,ny,nz,ens_num)
  real                  :: stdphi_init,stdpsi_init
    
  real                  :: temr01,temr02,temr03,temr04,temr05,temr06,temr07,temr08,temr09,temr10
  integer               :: temi01,temi02,temi03,temi04,temi05,temi06,temi07,temi08,temi09,temi10
  
  real                  :: upert(nx,ny,nz,ens_num)
  real                  :: vpert(nx,ny,nz,ens_num)
  real                  :: wpert(nx,ny,nz,ens_num)
  real                  :: phpert(nx,ny,nz,ens_num)
  real                  :: tpert(nx,ny,nz,ens_num)
  real                  :: qvpert(nx,ny,nz,ens_num)
  real                  :: smoispert(nx,ny,4,ens_num)
  real                  :: tslbpert(nx,ny,4,ens_num)
    
  real                  :: upert_mean(nx,ny,nz)
  real                  :: vpert_mean(nx,ny,nz)
  real                  :: wpert_mean(nx,ny,nz)
  real                  :: phpert_mean(nx,ny,nz)
  real                  :: tpert_mean(nx,ny,nz)
  real                  :: smoispert_mean(nx,ny,nz)
  real                  :: tslbpert_mean(nx,ny,nz)  
  
  real                  :: rho(nx,ny,nz),rhobar(nx,ny,nz)
  real                  :: rhou(nx,ny,nz),rhov(nx,ny,nz),rhow(nx,ny,nz)
  real                  :: p(nx,ny,nz),ppb(nx,ny,nz)
  
  real                  :: phpert_laplace(nx,ny,nz)  
  real                  :: phpert_temp(nx,ny,nz)
  
  real                  :: udu_vdu(nx,ny,nz)
  real                  :: udv_vdv(nx,ny,nz)
  real                  :: fu(nx,ny,nz)
  real                  :: fv(nx,ny,nz)
  
  real                  :: curl(nx,ny,nz)
  
  real                  :: pi(nx,ny,nz)
  
  integer               :: ii,jj,kk
  integer               :: iter
  
  real                  :: mapfact_m(nx,ny)
  real                  :: mapfact_u(nx,ny)
  real                  :: mapfact_v(nx,ny)

phi=0
psi=0
upert=0
vpert=0
wpert=0
phpert=0
tpert=0
qvpert=0
smoispert=0
tslbpert=0

upert_mean=0
vpert_mean=0
wpert_mean=0
phpert_mean=0
tpert_mean=0
smoispert_mean=0
tslbpert_mean=0

rho=0
rhobar=0
rhou=0
rhov=0
rhow=0
p=0
ppb=0
phpert_laplace=0
phpert_temp=0
udu_vdu=0
udv_vdv=0
fu=0
fv=0
curl=0


do j=1,ny
 do i=1,nx
   call cal_mapfactor_lambert(xlat(i,j),mapfact_m(i,j))   
 enddo
enddo

do j=1,ny
 do i=1,nx
   call cal_mapfactor_lambert(ulat(i,j),mapfact_u(i,j))   
 enddo
enddo

do j=1,ny
 do i=1,nx
   call cal_mapfactor_lambert(vlat(i,j),mapfact_v(i,j))   
 enddo
enddo
!****************************************************  
!  
!              calculate rho 
!
!**************************************************** 
print*,'cal rho init'
    do k=2,nz-1
      do j=1,ny
        do i=1,nx  
       call cal_rho_init(mub(i,j,1),mu(i,j,1),ensqv(i,j,k,1),ensph(i,j,k,1),ensph(i,j,k+1,1),    &
                         rdnw(k),p_top,znu(k),enst(i,j,k,1),                             &
                         nx,ny,nz,rho(i,j,k),p(i,j,k),ppb(i,j,k)                         &                  
                        )                            
        enddo
      enddo
    enddo 

    do k=2,nz-1
      do j=2,ny-1
        do i=2,nx-1  
           rhou(i,j,k)=(rho(i+1,j,k)+rho(i,j,k))/2
           rhov(i,j,k)=(rho(i,j+1,k)+rho(i,j,k))/2
           rhow(i,j,k)=(rho(i,j,k+1)+rho(i,j,k))/2
        enddo
      enddo
    enddo 
!****************************************************  
!  
!              give the random phi and psi 
!
!**************************************************** 
print*,'set phi and psi pert'
  do i=1,ens_num
     phi(:,:,:,i)= ranpert(:,:,:,i)*stdphi_init
     psi(:,:,:,i)=-ranpert(:,:,:,i)*stdpsi_init
!     upert(:,:,:,i)= ranpert(:,:,:,i)*3
!     vpert(:,:,:,i)=-ranpert(:,:,:,i)*3
  enddo

  do i=1,ens_num
    do k=1,nz
      if(k>=3*nz/4) then
       phi(:,:,k,i)=0
       psi(:,:,k,i)=0
      endif
    enddo 
  enddo
  
!****************************************************  
!  
!          perturb U V will get from  
!          the given random phi and psi
!
!****************************************************   
print*,'cal perturb U V'
  do ii=1,ens_num
    do k=1,nz
      do j=3,ny-2
        do i=3,nx-2
        
           temr01=(phi(i-1,j+1,k,ii)+phi(i,j+1,k,ii))*0.5          ! phi j+1 ave
           temr02=(phi(i-1,j-1,k,ii)+phi(i,j-1,k,ii))*0.5          ! phi j ave
           temr03=-0.5*(temr01-temr02)/dy                          ! -d phi / d y
           
           temr04=0.25*(psi(i+1,j,k,ii)+psi(i,j,k,ii)-psi(i-2,j,k,ii)-psi(i-1,j,k,ii))/dx           !  d psi / d x
           
           temr05=temr03+temr04                                ! ustr= -d phi / d y + d psi / d x          
           upert(i,j,k,ii)=temr05*mapfact_u(i,j)
                   
           temr01=(phi(i-1,j-1,k,ii)+phi(i-1,j,k,ii))*0.5          ! phi i+1 ave
           temr02=(phi(i+1,j-1,k,ii)+phi(i+1,j,k,ii))*0.5          ! phi i-1 ave
           temr03=-0.5*(temr01-temr02)/dx*mapfact_u(i,j)                          !  d phi / d x
           
           temr04=0.25*(psi(i,j+1,k,ii)+psi(i,j,k,ii)-psi(i,j-1,k,ii)-psi(i,j-2,k,ii))/dy           !  d psi / d y
           
           temr05=temr03+temr04                                ! vstr= d phi / d x + d psi / d y          
           vpert(i,j,k,ii)=temr05*mapfact_v(i,j)        
                       
        enddo
      enddo
    enddo
  enddo
  
  print*,'remove bias mean'
  do k=1,nz
    do j=1,ny
      do i=1,nx
        do ii=1,ens_num
            upert_mean(i,j,k)=upert_mean(i,j,k)+upert(i,j,k,ii)
        enddo
            upert_mean(i,j,k)=upert_mean(i,j,k)/ens_num
      enddo
    enddo
  enddo
  
   do ii=1,ens_num
    do k=1,nz
      do j=1,ny
        do i=1,nx  
           upert(i,j,k,ii)=upert(i,j,k,ii)-upert_mean(i,j,k)
        enddo
      enddo
    enddo  
  enddo 
  
  do k=1,nz
    do j=1,ny
      do i=1,nx
        do ii=1,ens_num
            vpert_mean(i,j,k)=vpert_mean(i,j,k)+vpert(i,j,k,ii)
        enddo
            vpert_mean(i,j,k)=vpert_mean(i,j,k)/ens_num
      enddo
    enddo
  enddo
  
   do ii=1,ens_num
    do k=1,nz
      do j=1,ny
        do i=1,nx  
           vpert(i,j,k,ii)=vpert(i,j,k,ii)-vpert_mean(i,j,k)
        enddo
      enddo
    enddo  
  enddo   
!****************************************************  
!  
!          perturb w will get from ustr & vstr
!          
!
!****************************************************   
print*,'cal perturb W' 

  do ii=1,ens_num
  
    do k=nz-2,3,-1
      do j=3,ny-2
        do i=3,nx-2
             
             temr01= ( (upert(i+2,j,k,ii)*rhou(i+2,j,k)+upert(i+1,j,k,ii)*rhou(i+1,j,k))/mapfact_m(i+1,j)      &
                      -(upert(i,j,k,ii)*rhou(i,j,k)+upert(i-1,j,k,ii)*rhou(i-1,j,k))/mapfact_m(i-1,j)   )/dx/4*mapfact_m(i,j)**2                ! d rho*u / dx
             temr02= ( (vpert(i,j+2,k,ii)*rhov(i,j+2,k)+vpert(i,j+1,k,ii)*rhov(i,j+1,k))/mapfact_m(i,j+1)      &
                      -(vpert(i,j,k,ii)*rhov(i,j,k)+vpert(i,j-1,k,ii)*rhov(i,j-1,k))/mapfact_m(i,j-1)   )/dy/4*mapfact_m(i,j)**2                ! d rho*v / dy
             temr03=wpert(i,j,k+1,ii)*rhow(i,j,k+1)-(temr01+temr02)*(phb(i,j,k+1)-phb(i,j,k))/9.8
             wpert(i,j,k,ii)=temr03/rhow(i,j,k)
             
        enddo
      enddo
    enddo
    
  enddo  
  print*,'remove bias mean'
  do k=1,nz
    do j=1,ny
      do i=1,nx
        do ii=1,ens_num
            wpert_mean(i,j,k)=wpert_mean(i,j,k)+upert(i,j,k,ii)
        enddo
            wpert_mean(i,j,k)=wpert_mean(i,j,k)/ens_num
      enddo
    enddo
  enddo
  
   do ii=1,ens_num
    do k=1,nz
      do j=1,ny
        do i=1,nx  
           wpert(i,j,k,ii)=wpert(i,j,k,ii)-wpert_mean(i,j,k)
        enddo
      enddo
    enddo  
  enddo 
  
!****************************************************  
!  
!          perturb PH will get from  
!              the ustr & vstr
!
!****************************************************     

print*,'call perturb PH'
  do ii=1,ens_num
  
  
    do k=2,nz-1
      do j=3,ny-2
        do i=3,nx-2
           

!   udu_vdu(i,j,k)=  										   &
!		     0.5*(ensu(i,j,k,ii)+ensu(i+1,j,k,ii))*(upert(i+1,j,k,ii)-upert(i,j,k,ii))/dx  &
!                   +0.5*(upert(i,j,k,ii)+upert(i+1,j,k,ii))*(ensu(i+1,j,k,ii)-ensu(i,j,k,ii))/dx  &
!                   +0.25*ensv(i,j+1,k,ii)*(upert(i,j+1,k,ii)+upert(i+1,j+1,k,ii)-upert(i,j,k,ii)-upert(i+1,j,k,ii))/dy  &
!                   +0.25*ensv(i,j,k,ii)*(upert(i,j,k,ii)+upert(i+1,j,k,ii)-upert(i,j-1,k,ii)-upert(i+1,j-1,k,ii))/dy   &
!                   +0.25*vpert(i,j+1,k,ii)*(ensu(i,j+1,k,ii)+ensu(i+1,j+1,k,ii)-ensu(i,j,k,ii)-ensu(i+1,j,k,ii))/dy  &
!                   +0.25*vpert(i,j,k,ii)*(ensu(i,j,k,ii)+ensu(i+1,j,k,ii)-ensu(i,j-1,k,ii)-ensu(i+1,j-1,k,ii))/dy   
!&   
!                   +0.5*(upert(i,j,k,ii)+upert(i+1,j,k,ii))*(upert(i+1,j,k,ii)-upert(i,j,k,ii))/dx  &
!                   +0.25*vpert(i,j+1,k,ii)*(upert(i,j+1,k,ii)+upert(i+1,j+1,k,ii)-upert(i,j,k,ii)-upert(i+1,j,k,ii))/dy  &
!                   +0.25*vpert(i,j,k,ii)*(upert(i,j,k,ii)+upert(i+1,j,k,ii)-upert(i,j-1,k,ii)-upert(i+1,j-1,k,ii))/dy    
                                    
!   udu_vdu(i,j,k)= udu_vdu(i,j,k)*rho(i,j,k)
        
!   udv_vdv(i,j,k)=                                                                                                      &
!                    0.25*ensu(i+1,j,k,ii)*(vpert(i+1,j+1,k,ii)+vpert(i+1,j,k,ii)-vpert(i,j+1,k,ii)-vpert(i,j,k,ii))/dx  &
!                   +0.25*ensu(i,j,k,ii)*(vpert(i,j+1,k,ii)+vpert(i,j,k,ii)-vpert(i-1,j+1,k,ii)-vpert(i-1,j,k,ii))/dx  &
!                   +0.25*upert(i+1,j,k,ii)*(ensv(i+1,j+1,k,ii)+ensv(i+1,j,k,ii)-ensv(i,j+1,k,ii)-ensv(i,j,k,ii))/dx  &
!                   +0.25*upert(i,j,k,ii)*(ensv(i,j+1,k,ii)+ensv(i,j,k,ii)-ensv(i-1,j+1,k,ii)-ensv(i-1,j,k,ii))/dx  &                   
!                   +0.5*(vpert(i,j,k,ii)+vpert(i,j+1,k,ii))*(ensv(i,j+1,k,ii)-ensv(i,j,k,ii))/dy                  &
!                   +0.5*(ensv(i,j,k,ii)+ensv(i,j+1,k,ii))*(vpert(i,j+1,k,ii)-vpert(i,j,k,ii))/dy                
! &
!                   +0.5*(vpert(i,j,k,ii)+vpert(i,j+1,k,ii))*(vpert(i,j+1,k,ii)-vpert(i,j,k,ii))/dy                  & 
!                   +0.25*upert(i+1,j,k,ii)*(vpert(i+1,j+1,k,ii)+vpert(i+1,j,k,ii)-vpert(i,j+1,k,ii)-vpert(i,j,k,ii))/dx  &
!                   +0.25*upert(i,j,k,ii)*(vpert(i,j+1,k,ii)+vpert(i,j,k,ii)-vpert(i-1,j+1,k,ii)-vpert(i-1,j,k,ii))/dx                                    
   
!   udv_vdv(i,j,k)= udv_vdv(i,j,k)*rho(i,j,k)


!      curl(i,j,k)=   (vpert(i,j+1,k,ii)+vpert(i,j,k,ii)-vpert(i-1,j+1,k,ii)-vpert(i-1,j,k,ii))/dx  &
!                    -(upert(i,j,k,ii)+upert(i+1,j,k,ii)-upert(i,j-1,k,ii)-upert(i+1,j-1,k,ii))/dy
                

                   
        fu(i,j,k)= (2*sin(ulat(i,j)*3.1415926/180.0)*7.292*(1e-5))*upert(i,j,k,ii)*0.5  &    !+curl(i,j,k)
                  +(2*sin(ulat(i+1,j)*3.1415926/180.0)*7.292*(1e-5))*upert(i+1,j,k,ii)*0.5   !+curl(i+1,j,k)
         
!        fu(i,j,k)= fu(i,j,k)*rho(i,j,k) 
         
        fv(i,j,k)= (2*sin(vlat(i,j)*3.1415926/180.0)*7.292*(1e-5))*vpert(i,j,k,ii)*0.5   &     !+curl(i,j,k)
                  +(2*sin(vlat(i,j+1)*3.1415926/180.0)*7.292*(1e-5))*vpert(i,j+1,k,ii)*0.5 !+curl(i,j+1,k)                
        
!        fv(i,j,k)= fv(i,j,k)*rho(i,j,k)

    udu_vdu(i,j,k)=  0.25*((upert(i+2,j,k,ii)+upert(i+1,j,k,ii))/mapfact_m(i+1,j)-(upert(i,j,k,ii)+upert(i-1,j,k,ii))/mapfact_m(i-1,j))/dx*mapfact_m(i,j)**2  &
                    *0.25*((vpert(i,j+2,k,ii)+vpert(i,j+1,k,ii))/mapfact_m(i,j+1)-(vpert(i,j,k,ii)+vpert(i,j-1,k,ii))/mapfact_m(i,j-1))/dy*mapfact_m(i,j)**2 
    udv_vdv(i,j,k)=  0.25*((upert(i,j+1,k,ii)+upert(i+1,j+1,k,ii))/mapfact_m(i,j+1)-(upert(i,j-1,k,ii)+upert(i+1,j-1,k,ii))/mapfact_m(i,j-1))/dy*mapfact_m(i,j)**2   &
                    *0.25*((vpert(i+1,j+1,k,ii)+vpert(i+1,j,k,ii))/mapfact_m(i+1,j)-(vpert(i-1,j+1,k,ii)+vpert(i-1,j,k,ii))/mapfact_m(i-1,j))/dx*mapfact_m(i,j)**2

            
        enddo        
      enddo
    enddo
    
     do k=2,nz-1
      do j=3,ny-2
        do i=3,nx-2
              
!            phpert_laplace(i,j,k)   = -0.5*(udu_vdu(i+1,j,k)-udu_vdu(i-1,j,k))/dx  &
!                                      -0.5*(udv_vdv(i,j+1,k)-udv_vdv(i,j-1,k))/dy  &
!                                      +0.5*(fv(i+1,j,k)-fv(i-1,j,k))/dx            &
!                                      -0.5*(fu(i,j+1,k)-fu(i,j-1,k))/dy


            
            phpert_laplace(i,j,k)   =  +0.5*(fv(i+1,j,k)-fv(i-1,j,k))/dx            &
                                       -0.5*(fu(i,j+1,k)-fu(i,j-1,k))/dy            &
                                       +2*(udu_vdu(i,j,k)-udv_vdv(i,j,k))            
                                      
        enddo         
      enddo
    enddo
    
    
    do iter=1,50
    
     do k=nz-1,3,-1
      do j=3,ny-2
        do i=3,nx-2
           phpert_temp(i,j,k)=0.25*(phpert(i+1,j,k,ii)+phpert(i-1,j,k,ii)+phpert(i,j+1,k,ii)+phpert(i,j-1,k,ii)-phpert_laplace(i,j,k)*dx*dy)    
        enddo         
      enddo
     enddo  
     
     do k=nz-1,3,-1
      do j=3,ny-2
        do i=3,nx-2     
           phpert(i,j,k,ii)=phpert_temp(i,j,k) 
        enddo         
      enddo
     enddo  
             
    enddo
    
    
   print*,'finished perturb PH menber ',ii    
  enddo  
  print*,'remove bias mean'
  do k=1,nz
    do j=1,ny
      do i=1,nx
        do ii=1,ens_num
            phpert_mean(i,j,k)=phpert_mean(i,j,k)-phpert(i,j,k,ii)
        enddo
            phpert_mean(i,j,k)=phpert_mean(i,j,k)/ens_num
      enddo
    enddo
  enddo
  
   do ii=1,ens_num
    do k=1,nz
      do j=1,ny
        do i=1,nx  
           phpert(i,j,k,ii)=(phpert(i,j,k,ii)-phpert_mean(i,j,k))*10
        enddo
      enddo
    enddo  
  enddo 


  print*,' '
  write(*,*) ((phpert (i,j,10,1),i=nx/2-3,nx/2+3),j=ny/2-3,ny/2+3)
  print*,' '

  
!****************************************************  
!  
!          perturb Theta will get from  
!              the perturb PH
!
!****************************************************  
print*,'call perturb theta'
!  do ii=1,ens_num

!    do k=3,nz-2
!      do j=2,ny-1
!        do i=2,nx-1  
!       call cal_rho_init(mub(i,j,ii),mu(i,j,ii),ensqv(i,j,k,ii),phpert(i,j,k,ii),phpert(i,j,k+1,ii),    &
!                         rdnw(k),p_top,znu(k),enst(i,j,k,ii),                                           &
!                         nx,ny,nz,rho(i,j,k),p(i,j,k),ppb(i,j,k)                                        &                  
!                        )    
!        enddo
!      enddo
!    enddo     
    
!    do k=3,nz-2
!      do j=2,ny-1
!        do i=2,nx-1 
!           rhobar(i,j,k)=-(ppb(i,j,k+1)-ppb(i,j,k-1))/(phb(i,j,k+1)-phb(i,j,k-1))
!        enddo
!      enddo
!    enddo  
    
    
!    do k=3,nz-3
!      do j=2,ny-1
!        do i=2,nx-1 
!        
!           call cal_tc_init(enst(i,j,k,ii),p(i,j,k),tc,nx,ny,nz)            
!           tpert(i,j,k,ii)=(tc+273.15)*((p(i,j,k)-ppb(i,j,k))/ppb(i,j,k)-(rho(i,j,k)-rhobar(i,j,k))/rhobar(i,j,k)  )
!           call cal_theta_init(tpert(i,j,k,ii),p(i,j,k),tpert(i,j,k,ii),nx,ny,nz)
!
!            tpert(i,j,k,ii)=(p(i,j,k)-ppb(i,j,k))/(287.05*rhobar(i,j,k))
!            call cal_theta_init(tpert(i,j,k,ii),p(i,j,k),tpert(i,j,k,ii),nx,ny,nz)
!        enddo
!      enddo
!    enddo  
        
!  enddo

  do i=1,ens_num  
   do k=1,nz
     if(k<3*nz/4) then
       tpert(:,:,k,i) = ranpert2(:,:,k,i)*stdt_init
       enst (:,:,k,i) = enst(:,:,k,i) - tpert(:,:,k,i)
     else
       enst (:,:,k,i) = enst(:,:,k,i)
     endif    
   enddo
  enddo

!  print*,'remove bias mean'  
!  do k=1,nz
!    do j=1,ny
!      do i=1,nx
!        do ii=1,ens_num
!            tpert_mean(i,j,k)=tpert_mean(i,j,k)+tpert(i,j,k,ii)
!        enddo
!            tpert_mean(i,j,k)=tpert_mean(i,j,k)/ens_num
!      enddo
!    enddo
!  enddo
  
!   do ii=1,ens_num
!    do k=1,nz
!      do j=1,ny
!        do i=1,nx  
!           tpert(i,j,k,ii)=tpert(i,j,k,ii)-tpert_mean(i,j,k)
!        enddo
!      enddo
!    enddo  
!  enddo   
!****************************************************  
!  
!          perturb qv will get by  random  number 
!
!****************************************************  
print*,'call perturb qv'
  do i=1,ens_num  
   do k=1,nz
     if(k<3*nz/4) then
       qvpert(:,:,k,i) = ranpert2(:,:,k,i)* stdqv_init
       ensqv (:,:,k,i) = ensqv(:,:,k,i) - qvpert(:,:,k,i)
     else
       ensqv (:,:,k,i) = ensqv(:,:,k,i)
     endif    
   enddo
  enddo
  
!****************************************************  
!  
!         adding perturb to background
!
!**************************************************** 
print*,'adding perturb to background' 
  ensu=ensu+upert
  ensv=ensv+vpert
  ensw=ensw+wpert
!  ensph=ensph+phpert
!  enst=enst+tpert
  
   do ii=1,ens_num
    do k=2,nz-1
      do j=1,ny
        do i=1,nx  
           ensph(i,j,k,ii)=ensph(i,j,k,ii)+(phpert(i,j,k,ii)+phpert(i,j,k-1,ii))*0.5
        enddo
      enddo
    enddo  
  enddo     
  
  
print*,'done' 

end  subroutine addpert_with_physical_constrain           



subroutine cal_rho_init(mub,mu,qv,phl,phu,rdnw,p_top,znu,t,       &
                   nx,ny,nz,rho,p,ppb                             &                  
                  )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine cal_rho
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer         ::  nx,ny,nz                  
real            ::  mub,mu,qv,phl,phu,rdnw,znu,p
real            ::  p_top
real            ::  rho  

REAL, PARAMETER    :: pi = 3.1415926535897932346
REAL, PARAMETER    :: gas_constant = 287.     ! Value used in WRF.
REAL, PARAMETER    :: gas_constant_v = 461.6  ! Value used in WRF.
REAL, PARAMETER    :: cp = 7.*gas_constant/2. ! Value used in WRF.
REAL, PARAMETER    :: t_kelvin = 273.15
REAL, PARAMETER    :: kappa = gas_constant / cp
REAL, PARAMETER    :: rd_over_rv = gas_constant / gas_constant_v  
REAL, PARAMETER    :: gravity = 9.81        ! m/s - value used in MM5.


real            :: albn,aln,ppb,ttb,qvf1,t
real            :: cvpm,cpovcv,ps0,ts0,tis0,tlp



ps0  = 100000.0    ! Base sea level pressure
ts0  = 300.0       ! Base potential temperature for all levels.
tis0 = 290.0       ! Base sea level temperature
tlp  = 50.0        ! temperature difference from 1000mb to 300mb


cvpm =  - (1. - gas_constant/cp)
cpovcv = cp / (cp - gas_constant)

   ppb  = znu * mub + p_top
   ttb  = (tis0 + tlp*log(ppb/ps0)) * (ps0/ppb)**kappa
   albn = (gas_constant/ps0) * ttb * (ppb/ps0)**cvpm

   qvf1 = 1. + qv / rd_over_rv
   aln  = -1. / (mub+mu) * ( albn*mu + rdnw *(phu - phl) )
               
   p = ps0 * ( (gas_constant*(ts0+t)*qvf1) / &
                       (ps0*(aln+albn)) )**cpovcv            

   rho= 1.0 / (albn+aln)
                                                    
end subroutine cal_rho_init        

subroutine cal_tc_init(t,p,tc,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine cal_tc
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer         ::  nx,ny,nz 
real            ::  t,p,tc

REAL, PARAMETER    :: gas_constant = 287.     ! Value used in WRF.
REAL, PARAMETER    :: cp = 7.*gas_constant/2. ! Value used in WRF. 
REAL, PARAMETER    :: t_kelvin = 273.15
real            :: ps0
REAL, PARAMETER    :: kappa = gas_constant / cp

ps0  = 100000.0    ! Base sea level pressure


  tc=(t+300)*(p/ps0)**kappa-t_kelvin 


end subroutine cal_tc_init     

subroutine cal_theta_init(t,p,tc,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine cal_tc
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer         ::  nx,ny,nz 
real            ::  t,p,tc

REAL, PARAMETER    :: gas_constant = 287.     ! Value used in WRF.
REAL, PARAMETER    :: cp = 7.*gas_constant/2. ! Value used in WRF. 
REAL, PARAMETER    :: t_kelvin = 273.15
real            :: ps0
REAL, PARAMETER    :: kappa = gas_constant / cp

ps0  = 100000.0    ! Base sea level pressure


  t=(tc+t_kelvin)*(ps0/p/100)**kappa-300


end subroutine cal_theta_init                     

subroutine cal_mapfactor_lambert(lat,mapfact)

real                 ::   lat
real                 ::   mapfact

real,parameter       ::   pi=3.1415926
real,parameter       ::   a=6371
real,parameter       ::   le=11423.37
real,parameter       ::   hd2jd=pi/180.0

real                 ::   l
real                 ::   k

k = (log(sin(60*hd2jd))-log(sin(30*hd2jd)))/(log(tan(60*hd2jd)/2)-log(tan(30*hd2jd/2)))

l =  le*(cos(lat*hd2jd)/(1+sin(lat*hd2jd)))**(k)

mapfact = k*l/(a*(1-sin(lat*hd2jd)**2)**(0.5))

end subroutine cal_mapfactor_lambert     
