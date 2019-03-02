subroutine mean_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,     &
                    nx,ny,nz,ens_num,                                                                     &
                    umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean &
                    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine mean_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: calculate the ensemble mean
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 implicit none
integer          :: nx,ny,nz,ii,j,k,nn,i 
integer          :: ens_num
real             :: ensu(nx,ny,nz,ens_num)
real             :: ensv(nx,ny,nz,ens_num)
real             :: ensw(nx,ny,nz,ens_num)
real             :: ensph(nx,ny,nz,ens_num)
real             :: enst(nx,ny,nz,ens_num)
real             :: ensqv(nx,ny,nz,ens_num)
real             :: ensqr(nx,ny,nz,ens_num)
real             :: ensqi(nx,ny,nz,ens_num)
real             :: ensqs(nx,ny,nz,ens_num)
real             :: ensqgr(nx,ny,nz,ens_num)
real             :: enssmois(nx,ny,4,ens_num)
real             :: enstslb(nx,ny,4,ens_num)
real             :: enstsk(nx,ny,ens_num)

real             :: umean(nx,ny,nz)
real             :: vmean(nx,ny,nz)
real             :: wmean(nx,ny,nz)
real             :: phmean(nx,ny,nz)
real             :: tmean(nx,ny,nz)
real             :: qvmean(nx,ny,nz)
real             :: qrmean(nx,ny,nz)
real             :: qimean(nx,ny,nz)
real             :: qsmean(nx,ny,nz)
real             :: qgrmean(nx,ny,nz)
real             :: smoismean(nx,ny,4)
real             :: tslbmean(nx,ny,4)
real             :: tskmean(nx,ny)

real             :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13

temp1=0
temp2=0
temp3=0
temp4=0
temp5=0
temp6=0
temp7=0
temp8=0
temp9=0
temp10=0
temp11=0
temp12=0
temp13=0

do k=1,nz
do j=1,ny-1
do i=1,nx-1

temp1=0
temp2=0
temp3=0
temp4=0
temp5=0
temp6=0 
temp7=0 
temp8=0 
temp9=0 
temp10=0
!temp11=0
!temp12=0

 do ii=1,ens_num
  temp1=temp1 +  ensu(i,j,k,ii)
  temp2=temp2 +  ensv(i,j,k,ii)
  temp3=temp3 +  ensw(i,j,k,ii)
  temp4=temp4 + ensph(i,j,k,ii)
  temp5=temp5 +  enst(i,j,k,ii)
  temp6=temp6 + ensqv(i,j,k,ii)
  temp7=temp7 + ensqr(i,j,k,ii)
  temp8=temp8 + ensqi(i,j,k,ii)
  temp9=temp9 + ensqs(i,j,k,ii)
  temp10=temp10 + ensqgr(i,j,k,ii)
!  temp11=temp11 + enssmois(i,j,k,ii)
!  temp12=temp12 + enstslb(i,j,k,ii)
 enddo 
  umean(i,j,k)  =temp1/ens_num  
  vmean(i,j,k)  =temp2/ens_num 
  wmean(i,j,k)  =temp3/ens_num 
  phmean(i,j,k) =temp4/ens_num 
  tmean(i,j,k)  =temp5/ens_num 
  qvmean(i,j,k) =temp6/ens_num 
  qrmean(i,j,k) =temp7/ens_num 
  qimean(i,j,k) =temp8/ens_num 
  qsmean(i,j,k) =temp9/ens_num 
  qgrmean(i,j,k) =temp10/ens_num
!  smoismean(i,j,k) =temp11/ens_num
!  tslbmean(i,j,k) =temp12/ens_num
enddo
enddo
enddo

do k=1,4
do j=1,ny-1
do i=1,nx-1
 temp11=0
 temp12=0

 do ii=1,ens_num

  temp11=temp11 + enssmois(i,j,k,ii)
  temp12=temp12 + enstslb(i,j,k,ii)
 enddo 
  
  smoismean(i,j,k) =temp11/ens_num
  tslbmean(i,j,k) =temp12/ens_num
  !if(smoismean(i,j,k) < 0)then
  ! print*,'wrong data of smois',smoismean(i,j,k)
  ! smoismean(i,j,k)=0
  !endif
  !if(smoismean(i,j,k) > 1)then
  ! print*,'wrong data of smois',smoismean(i,j,k)
  ! smoismean(i,j,k)=0.999
  !endif  
 ! if(tslbmean(i,j,k) < 233.15)then
 !  print*,'wrong small data of tslb',tslbmean(i,j,k)
 !  tslbmean(i,j,k)=233.15
 !endif
 ! if(tslbmean(i,j,k) > 313.15)then
 ! print*,'wrong big data of tslb',tslbmean(i,j,k)
 ! tslbmean(i,j,k)=313.15
!endif
enddo
enddo
enddo


do j=1,ny
do i=1,nx
do nn=1,ens_num 

!************************
!*** fill data like this:
!  ---------
!  1 2 3 4 |
!  2 3 4 5 |
!****************************

if( i == nx .and. enstsk(i,j,nn) == 0)then
 enstsk(i,j,nn)=enstsk(i-1,j,nn)
endif
if(j == ny.and. enstsk(i,j,nn) == 0)then
 enstsk(i,j,nn)=enstsk(i,j-1,nn)
endif  

if(enstsk (i,j,nn) < 233.15)then
 print*,'ckeck enstsk',enstsk(i,j,nn),i,j,nn
 enstsk (i,j,nn) = 233.15
 endif
 if(enstsk (i,j,nn)> 333.15)then
 print*,'ckeck enstsk',enstsk(i,j,nn)
 enstsk (i,j,nn)=333.15
 endif
enddo
enddo
enddo

do j=1,ny-1
do i=1,nx-1
 temp13=0
! temp12=0

 do ii=1,ens_num

!  temp11=temp11 + enssmois(i,j,k,ii)
  temp13=temp13 + enstsk(i,j,ii)
 enddo 
  
!  smoismean(i,j,k) =temp11/ens_num
  tskmean(i,j) =temp13/ens_num
!  if(tskmean(i,j) < 233.15)then
!   print*,'wrong small data of tslb',tskmean(i,j)
!   tskmean(i,j)=233.15
! endif
!  if(tskmean(i,j) > 313.15)then
!  print*,'wrong big data of tslb',tskmean(i,j)
!  tskmean(i,j)=313.15
! endif

enddo
enddo
!print*,tskmean(1,1)
end subroutine mean_cal


subroutine pert_cal(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,              &
                    upert,vpert,wpert,phpert,tpert,qvpert,qrpert,qipert,qspert,qgrpert,smoispert,tslbpert,tskpert, &
                    nx,ny,nz,ens_num,                                                                              &
                    umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean  &
                    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine pert_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: calculate the ensemble perturbation
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer          :: nx,ny,nz 
integer          :: ens_num,nn
real             :: ensu(nx,ny,nz,ens_num)
real             :: ensv(nx,ny,nz,ens_num)
real             :: ensw(nx,ny,nz,ens_num)
real             :: ensph(nx,ny,nz,ens_num)
real             :: enst(nx,ny,nz,ens_num)
real             :: ensqv(nx,ny,nz,ens_num)
real             :: ensqr(nx,ny,nz,ens_num)
real             :: ensqi(nx,ny,nz,ens_num)
real             :: ensqs(nx,ny,nz,ens_num)
real             :: ensqgr(nx,ny,nz,ens_num)
real             :: enssmois(nx,ny,4,ens_num)
real             :: enstslb(nx,ny,4,ens_num)
real             :: enstsk(nx,ny,ens_num)

real             :: upert(nx,ny,nz,ens_num)
real             :: vpert(nx,ny,nz,ens_num)
real             :: wpert(nx,ny,nz,ens_num)
real             :: phpert(nx,ny,nz,ens_num)
real             :: tpert(nx,ny,nz,ens_num)
real             :: qvpert(nx,ny,nz,ens_num)
real             :: qrpert(nx,ny,nz,ens_num)
real             :: qipert(nx,ny,nz,ens_num)
real             :: qspert(nx,ny,nz,ens_num)
real             :: qgrpert(nx,ny,nz,ens_num)
real             :: smoispert(nx,ny,4,ens_num)
real             :: tslbpert(nx,ny,4,ens_num)
real             :: tskpert(nx,ny,ens_num)

real             :: umean(nx,ny,nz)
real             :: vmean(nx,ny,nz)
real             :: wmean(nx,ny,nz)
real             :: phmean(nx,ny,nz)
real             :: tmean(nx,ny,nz)
real             :: qvmean(nx,ny,nz)
real             :: qrmean(nx,ny,nz)
real             :: qimean(nx,ny,nz)
real             :: qsmean(nx,ny,nz)
real             :: qgrmean(nx,ny,nz)
real             :: smoismean(nx,ny,4)
real             :: tslbmean(nx,ny,4)
real             :: tskmean(nx,ny)


do nn=1,ens_num

do k=1,nz
do j=1,ny-1
do i=1,nx-1

  upert(i,j,k,nn)  = ensu(i,j,k,nn)  - umean(i,j,k)
  vpert(i,j,k,nn)  = ensv(i,j,k,nn)  - vmean(i,j,k)
  wpert(i,j,k,nn)  = ensw(i,j,k,nn)  - wmean(i,j,k)
  phpert(i,j,k,nn) = ensph(i,j,k,nn) - phmean(i,j,k)
  tpert(i,j,k,nn)  = enst(i,j,k,nn)  - tmean(i,j,k)
  qvpert(i,j,k,nn) = ensqv(i,j,k,nn) - qvmean(i,j,k)
  qrpert(i,j,k,nn) = ensqr(i,j,k,nn) - qrmean(i,j,k)
  qipert(i,j,k,nn) = ensqi(i,j,k,nn) - qimean(i,j,k)
  qspert(i,j,k,nn) = ensqs(i,j,k,nn) - qsmean(i,j,k)
  qgrpert(i,j,k,nn)= ensqgr(i,j,k,nn)- qgrmean(i,j,k)
 ! smoispert(i,j,k,nn) = enssmois(i,j,k,nn) - smoismean(i,j,k)
 ! tslbpert(i,j,k,nn) = enstslb(i,j,k,nn) - tslbmean(i,j,k)

enddo
enddo
enddo

enddo

do nn=1,ens_num

do k=1,4
do j=1,ny-1
do i=1,nx-1

  smoispert(i,j,k,nn) = enssmois(i,j,k,nn) - smoismean(i,j,k)
  tslbpert(i,j,k,nn) = enstslb(i,j,k,nn) - tslbmean(i,j,k)

enddo
enddo
enddo

enddo
!!!
do nn=1,ens_num

!do k=1,4
do j=1,ny-1
do i=1,nx-1

  tskpert(i,j,nn) = enstsk(i,j,nn) - tskmean(i,j)

enddo
enddo
!enddo

enddo



end subroutine pert_cal

subroutine rms_cal(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,smoistrue,tslbtrue,tsktrue,                &
                   nx,ny,nz,ens_num,                                                                                             &
                   umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,                &
                   stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,stdtsk,                             &
                   stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,stdsmois2d,stdtslb2d,stdtsk2d    &
                  )

!Guoyk 2013/6/14 11:00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine rms_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: calculate the rmse
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer          :: nx,ny,nz 
integer          :: ens_num
real             :: utrue(nx,ny,nz)
real             :: vtrue(nx,ny,nz)
real             :: wtrue(nx,ny,nz)
real             :: phtrue(nx,ny,nz)
real             :: ttrue(nx,ny,nz)
real             :: qvtrue(nx,ny,nz)
real             :: qrtrue(nx,ny,nz)
real             :: qitrue(nx,ny,nz)
real             :: qstrue(nx,ny,nz)
real             :: qgrtrue(nx,ny,nz)
real             :: smoistrue(nx,ny,4)
real             :: tslbtrue(nx,ny,4)
real             :: tsktrue(nx,ny)

real             :: umean(nx,ny,nz)
real             :: vmean(nx,ny,nz)
real             :: wmean(nx,ny,nz)
real             :: phmean(nx,ny,nz)
real             :: tmean(nx,ny,nz)
real             :: qvmean(nx,ny,nz)
real             :: qrmean(nx,ny,nz)
real             :: qimean(nx,ny,nz)
real             :: qsmean(nx,ny,nz)
real             :: qgrmean(nx,ny,nz)
real             :: smoismean(nx,ny,4)
real             :: tslbmean(nx,ny,4)
real             :: tskmean(nx,ny)

real             :: ustd(nx,ny,nz)
real             :: vstd(nx,ny,nz)
real             :: wstd(nx,ny,nz)
real             :: phstd(nx,ny,nz)
real             :: tstd(nx,ny,nz)
real             :: qvstd(nx,ny,nz)
real             :: qrstd(nx,ny,nz)
real             :: qistd(nx,ny,nz)
real             :: qsstd(nx,ny,nz)
real             :: qgrstd(nx,ny,nz)
real             :: smoisstd(nx,ny,4)
real             :: tslbstd(nx,ny,4)
real             :: tskstd(nx,ny)

real             :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,stdtsk
real             :: stdu2d(nz),stdv2d(nz),stdw2d(nz),stdph2d(nz),stdt2d(nz),stdqv2d(nz),stdqr2d(nz),stdqi2d(nz),stdqs2d(nz),stdqgr2d(nz),stdsmois2d(4),stdtslb2d(4),&
                    stdtsk2d 

real             :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13

temp1=0
temp2=0
temp3=0
temp4=0
temp5=0
temp6=0
temp7=0
temp8=0
temp9=0
temp10=0
!temp11=0
!temp12=0

ustd=0
vstd=0
wstd=0
phstd=0
tstd=0
qvstd=0
qrstd=0 
qistd=0 
qsstd=0
qgrstd=0
!smoisstd=0
!tslbstd=0

stdu=0
stdv=0
stdw=0
stdph=0
stdt=0
stdqv=0
stdqr=0
stdqi=0
stdqs=0
stdqgr=0
!stdsmois=0
!stdtslb=0

  do k=1,nz 
  do j=1,ny
  do i=1,nx
 
     ustd(i,j,k)  = ustd(i,j,k)  + ( utrue(i,j,k)  - umean(i,j,k)  )**2
     vstd(i,j,k)  = vstd(i,j,k)  + ( vtrue(i,j,k)  - vmean(i,j,k)  )**2
     wstd(i,j,k)  = wstd(i,j,k)  + ( wtrue(i,j,k)  - wmean(i,j,k)  )**2
     phstd(i,j,k) = phstd(i,j,k) + ( phtrue(i,j,k) - phmean(i,j,k) )**2
     tstd(i,j,k)  = tstd(i,j,k)  + ( ttrue(i,j,k)  - tmean(i,j,k)  )**2
     qvstd(i,j,k) = qvstd(i,j,k) + ( qvtrue(i,j,k) - qvmean(i,j,k) )**2
     qrstd(i,j,k) = qrstd(i,j,k) + ( qstrue(i,j,k) - qrmean(i,j,k) )**2
     qistd(i,j,k) = qistd(i,j,k) + ( qitrue(i,j,k) - qimean(i,j,k) )**2
     qsstd(i,j,k) = qsstd(i,j,k) + ( qstrue(i,j,k) - qsmean(i,j,k) )**2
     qgrstd(i,j,k)= qgrstd(i,j,k)+ ( qgrtrue(i,j,k)- qgrmean(i,j,k))**2
 !    smoisstd(i,j,k)= smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
 !    tslbstd(i,j,k) =  tslbstd(i,j,k) + ( tslbtrue(i,j,k) - tslbmean(i,j,k) )**2

     stdu=stdu+ustd(i,j,k)
     stdv=stdv+vstd(i,j,k)
     stdw=stdw+wstd(i,j,k)
     stdph=stdph+phstd(i,j,k)
     stdt=stdt+tstd(i,j,k)
     stdqv=stdqv+qvstd(i,j,k)   
     stdqr=stdqr+qrstd(i,j,k) 
     stdqi=stdqi+qistd(i,j,k) 
     stdqs=stdqs+qsstd(i,j,k) 
     stdqgr=stdqgr+qgrstd(i,j,k) 
!     stdsmois=stdsmois+smoisstd(i,j,k)
!     stdtslb=stdtslb+tslbstd(i,j,k)

     temp1=temp1+1 
        
  enddo
  enddo
  enddo
    
  stdu=sqrt(stdu/(temp1))
  stdv=sqrt(stdv/(temp1))
  stdw=sqrt(stdw/(temp1))
  stdph=sqrt(stdph/(temp1))
  stdt=sqrt(stdt/(temp1))
  stdqv=sqrt(stdqv/(temp1))  
  stdqr=sqrt(stdqr/(temp1))  
  stdqi=sqrt(stdqi/(temp1))  
  stdqs=sqrt(stdqs/(temp1))    
  stdqgr=sqrt(stdqgr/(temp1)) 
!  stdsmois=sqrt(stdsmois/(temp1))
!  stdtslb=sqrt(stdtslb/(temp1))
  
print*,'rms done'

ustd=0
vstd=0
wstd=0
phstd=0
tstd=0
qvstd=0    
qrstd=0 
qistd=0 
qsstd=0 
qgrstd=0
!smoisstd=0
!tslbstd=0 

stdu2d=0
stdv2d=0
stdw2d=0
stdph2d=0
stdt2d=0
stdqv2d=0
stdqr2d=0
stdqi2d=0
stdqs2d=0
stdqgr2d=0
!stdsmois2d=0
!stdtslb2d=0

  do k=1,nz 
  
    temp1=0
    
  do j=1,ny
  do i=1,nx

     ustd(i,j,k)  = ustd(i,j,k)  + ( utrue(i,j,k)  - umean(i,j,k)  )**2
     vstd(i,j,k)  = vstd(i,j,k)  + ( vtrue(i,j,k)  - vmean(i,j,k)  )**2
     wstd(i,j,k)  = wstd(i,j,k)  + ( wtrue(i,j,k)  - wmean(i,j,k)  )**2
     phstd(i,j,k) = phstd(i,j,k) + ( phtrue(i,j,k) - phmean(i,j,k) )**2
     tstd(i,j,k)  = tstd(i,j,k)  + ( ttrue(i,j,k)  - tmean(i,j,k)  )**2
     qvstd(i,j,k) = qvstd(i,j,k) + ( qvtrue(i,j,k) - qvmean(i,j,k) )**2
     qrstd(i,j,k) = qrstd(i,j,k) + ( qstrue(i,j,k) - qrmean(i,j,k) )**2
     qistd(i,j,k) = qistd(i,j,k) + ( qitrue(i,j,k) - qimean(i,j,k) )**2
     qsstd(i,j,k) = qsstd(i,j,k) + ( qstrue(i,j,k) - qsmean(i,j,k) )**2
     qgrstd(i,j,k)= qgrstd(i,j,k)+ ( qgrtrue(i,j,k)- qgrmean(i,j,k))**2
  !   smoisstd(i,j,k) = smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
  !   tslbstd(i,j,k) = tslbstd(i,j,k) + ( tslbtrue(i,j,k) - tslbmean(i,j,k) )**2
     
     stdu2d(k)=stdu2d(k)+ustd(i,j,k)
     stdv2d(k)=stdv2d(k)+vstd(i,j,k)
     stdw2d(k)=stdw2d(k)+wstd(i,j,k)
     stdph2d(k)=stdph2d(k)+phstd(i,j,k)
     stdt2d(k)=stdt2d(k)+tstd(i,j,k)
     stdqv2d(k)=stdqv2d(k)+qvstd(i,j,k)   
     stdqr2d(k)=stdqr2d(k)+qrstd(i,j,k) 
     stdqi2d(k)=stdqi2d(k)+qistd(i,j,k) 
     stdqs2d(k)=stdqs2d(k)+qsstd(i,j,k) 
     stdqgr2d(k)=stdqgr2d(k)+qgrstd(i,j,k)
 !    stdsmois2d(k)=stdsmois2d(k)+smoisstd(i,j,k)
 !    stdtslb2d(k)=stdtslb2d(k)+tslbstd(i,j,k)               
     
     temp1=temp1+1 

  enddo
  enddo
    
  stdu2d(k)=sqrt(stdu2d(k)/(temp1))
  stdv2d(k)=sqrt(stdv2d(k)/(temp1))
  stdw2d(k)=sqrt(stdw2d(k)/(temp1))
  stdph2d(k)=sqrt(stdph2d(k)/(temp1))
  stdt2d(k)=sqrt(stdt2d(k)/(temp1))
  stdqv2d(k)=sqrt(stdqv2d(k)/(temp1)) 
  stdqr2d(k)=sqrt(stdqr2d(k)/(temp1)) 
  stdqi2d(k)=sqrt(stdqi2d(k)/(temp1)) 
  stdqs2d(k)=sqrt(stdqs2d(k)/(temp1)) 
  stdqgr2d(k)=sqrt(stdqgr2d(k)/(temp1))
!  stdsmois2d(k)=sqrt(stdsmois2d(k)/(temp1))
!  stdtslb2d(k)=sqrt(stdtslb2d(k)/(temp1)) 
  
  enddo
 
print*,'horizontal rms done' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!Soil layer ,guoyk!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp1=0
temp11=0
temp12=0
temp13=0

smoisstd=0
tslbstd=0


stdsmois=0
stdtslb=0

  do k=1,4 
  do j=1,ny-1
  do i=1,nx-1
 
     smoisstd(i,j,k)= smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
     tslbstd(i,j,k) =  tslbstd(i,j,k) + ( tslbtrue(i,j,k) - tslbmean(i,j,k) )**2
 

     stdsmois=stdsmois+smoisstd(i,j,k)
     stdtslb=stdtslb+tslbstd(i,j,k)

     temp1=temp1+1 
        
  enddo
  enddo
  enddo
    
  stdsmois=sqrt(stdsmois/(temp1))
  stdtslb=sqrt(stdtslb/(temp1))

!!!
tskstd=0
stdtsk=0

!  do k=1,4 
  do j=1,ny-1
  do i=1,nx-1
 
!     smoisstd(i,j,k)= smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
     tskstd(i,j) =  tskstd(i,j) + ( tsktrue(i,j) - tskmean(i,j) )**2
 
!     stdsmois=stdsmois+smoisstd(i,j,k)
     stdtsk=stdtsk+tskstd(i,j)

     temp13=temp13+1 
        
  enddo
  enddo
!  enddo
    
!  stdsmois=sqrt(stdsmois/(temp1))
  stdtsk=sqrt(stdtsk/(temp13))
  
print*,'rmslda done'


smoisstd=0
tslbstd=0 


stdsmois2d=0
stdtslb2d=0

  do k=1,4 
  
    temp1=0
    
  do j=1,ny-1
  do i=1,nx-1

     smoisstd(i,j,k) = smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
     tslbstd(i,j,k) = tslbstd(i,j,k) + ( tslbtrue(i,j,k) - tslbmean(i,j,k) )**2
     
     stdsmois2d(k)=stdsmois2d(k)+smoisstd(i,j,k)
     stdtslb2d(k)=stdtslb2d(k)+tslbstd(i,j,k)               
     
     temp1=temp1+1 

  enddo
  enddo
  ! print*,stdsmois2d(k),temp1,smoistrue(100,140,k),smoismean(100,140,k)
  stdsmois2d(k)=sqrt(stdsmois2d(k)/(temp1))
  stdtslb2d(k)=sqrt(stdtslb2d(k)/(temp1)) 
  !print*,stdsmois2d(k),temp1
  enddo
!!!!!!

tskstd=0 
stdtsk2d=0
    
temp13=0
    
  do j=1,ny-1
  do i=1,nx-1

!     smoisstd(i,j,k) = smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
     tskstd(i,j) = tskstd(i,j) + ( tsktrue(i,j) - tskmean(i,j) )**2
     
!     stdsmois2d(k)=stdsmois2d(k)+smoisstd(i,j,k)
     stdtsk2d=stdtsk2d+tskstd(i,j)               
     
     temp13=temp13+1 

  enddo
  enddo
  ! print*,stdsmois2d(k),temp1,smoistrue(100,140,k),smoismean(100,140,k)
!  stdsmois2d(k)=sqrt(stdsmois2d(k)/(temp1))
  stdtsk=sqrt(stdtsk2d/(temp13)) 
  !print*,stdsmois2d(k),temp1
!  enddo
 
print*,'horizontal rmslda done'

end subroutine rms_cal


subroutine ens_spread(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,enstsk,       &
                      nx,ny,nz,ens_num,                                                       &
                      umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,tskmean,  &
                      spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk   &
                      )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine ens_spread
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: calculate the ensemble spread
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer          :: nx,ny,nz 
integer          :: ii,jj,kk,nn
integer          :: ens_num
real             :: ensu(nx,ny,nz,ens_num)
real             :: ensv(nx,ny,nz,ens_num)
real             :: ensw(nx,ny,nz,ens_num)
real             :: ensph(nx,ny,nz,ens_num)
real             :: enst(nx,ny,nz,ens_num)
real             :: ensqv(nx,ny,nz,ens_num)
real             :: ensqr(nx,ny,nz,ens_num)
real             :: ensqi(nx,ny,nz,ens_num)
real             :: ensqs(nx,ny,nz,ens_num)
real             :: ensqgr(nx,ny,nz,ens_num)
real             :: enssmois(nx,ny,4,ens_num)
real             :: enstslb(nx,ny,4,ens_num)
real             :: enstsk(nx,ny,ens_num) !!

real             :: umean(nx,ny,nz)
real             :: vmean(nx,ny,nz)
real             :: wmean(nx,ny,nz)
real             :: phmean(nx,ny,nz)
real             :: tmean(nx,ny,nz)
real             :: qvmean(nx,ny,nz)
real             :: qrmean(nx,ny,nz)
real             :: qimean(nx,ny,nz)
real             :: qsmean(nx,ny,nz)
real             :: qgrmean(nx,ny,nz)
real             :: smoismean(nx,ny,4)
real             :: tslbmean(nx,ny,4)
real             :: tskmean(nx,ny)   !!

real             :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk

real             :: gspd_u(ens_num),gspd_v(ens_num),gspd_w(ens_num),gspd_ph(ens_num),gspd_t(ens_num),gspd_qv(ens_num), &
                    gspd_qr(ens_num),gspd_qi(ens_num),gspd_qs(ens_num),gspd_qgr(ens_num),gspd_smois(ens_num),gspd_tslb(ens_num),&
                    gspd_tsk(ens_num)
integer          :: tem1

do nn=1,ens_num
gspd_u(nn)=0
gspd_v(nn)=0
gspd_w(nn)=0
gspd_ph(nn)=0
gspd_t(nn)=0
gspd_qv(nn)=0
gspd_qr(nn)=0
gspd_qi(nn)=0
gspd_qs(nn)=0
gspd_qgr(nn)=0
gspd_smois(nn)=0
gspd_tslb(nn)=0
gspd_tsk(nn)=0 !!
enddo

spd_u=0
spd_v=0
spd_w=0
spd_ph=0
spd_t=0
spd_qv=0
spd_qr=0
spd_qi=0
spd_qs=0
spd_qgr=0
spd_smois=0
spd_tslb=0
spd_tsk=0

do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_u(nn)=gspd_u(nn)+(ensu(i,j,k,nn)-umean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_u=spd_u+sqrt(gspd_u(nn)/tem1)/ens_num
enddo



do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_v(nn)=gspd_v(nn)+(ensv(i,j,k,nn)-vmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_v=spd_v+sqrt(gspd_v(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_w(nn)=gspd_w(nn)+(ensw(i,j,k,nn)-wmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_w=spd_w+sqrt(gspd_w(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_ph(nn)=gspd_ph(nn)+(ensph(i,j,k,nn)-phmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_ph=spd_ph+sqrt(gspd_ph(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_t(nn)=gspd_t(nn)+(enst(i,j,k,nn)-tmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_t=spd_t+sqrt(gspd_t(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_qv(nn)=gspd_qv(nn)+(ensqv(i,j,k,nn)-qvmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_qv=spd_qv+sqrt(gspd_qv(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_qr(nn)=gspd_qr(nn)+(ensqr(i,j,k,nn)-qrmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_qr=spd_qr+sqrt(gspd_qr(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_qi(nn)=gspd_qi(nn)+(ensqi(i,j,k,nn)-qimean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_qi=spd_qi+sqrt(gspd_qi(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1 
    gspd_qs(nn)=gspd_qs(nn)+(ensqs(i,j,k,nn)-qsmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_qs=spd_qs+sqrt(gspd_qs(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz-1
do j=1,ny-1
do i=1,nx-1
    gspd_qgr(nn)=gspd_qgr(nn)+(ensqgr(i,j,k,nn)-qgrmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_qgr=spd_qgr+sqrt(gspd_qgr(nn)/tem1)/ens_num
enddo

do nn=1,ens_num
tem1=0
do k=1,4
do j=1,ny-1
do i=1,nx-1
    gspd_smois(nn)=gspd_smois(nn)+(enssmois(i,j,k,nn)-smoismean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_smois=spd_smois+sqrt(gspd_smois(nn)/tem1)/ens_num
enddo

do nn=1,ens_num
tem1=0
do k=1,4
do j=1,ny-1
do i=1,nx-1
    gspd_tslb(nn)=gspd_tslb(nn)+(enstslb(i,j,k,nn)-tslbmean(i,j,k))**2
    tem1=tem1+1
enddo
enddo
enddo
  spd_tslb=spd_tslb+sqrt(gspd_tslb(nn)/tem1)/ens_num
enddo

do nn=1,ens_num
! print*,nn,gspd_tsk(nn),enstsk(1,1,nn),tskmean(1,1)
tem1=0
!do k=1,4
do j=1,ny-1
do i=1,nx-1
    gspd_tsk(nn)=gspd_tsk(nn)+(enstsk(i,j,nn)-tskmean(i,j))**2
    tem1=tem1+1

enddo
enddo
!enddo
 
  spd_tsk=spd_tsk+sqrt(gspd_tsk(nn)/tem1)/ens_num
!print*,spd_tsk
enddo


end subroutine ens_spread



subroutine fill_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine fill_xstat
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: write ensemble into vetor X
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
integer          :: nx,ny,nz,analysis_var_num,ens_num,numstat,numstatsoil 
real             :: xstat(numstat,ens_num)
!real             :: xstatsoil(numstatsoil,ens_num)
real             :: u(nx,ny,nz,ens_num)
real             :: v(nx,ny,nz,ens_num)
real             :: w(nx,ny,nz,ens_num)
real             :: ph(nx,ny,nz,ens_num)
real             :: t(nx,ny,nz,ens_num)
real             :: qv(nx,ny,nz,ens_num)
real             :: qr(nx,ny,nz,ens_num)
real             :: qi(nx,ny,nz,ens_num)
real             :: qs(nx,ny,nz,ens_num)
real             :: qgr(nx,ny,nz,ens_num)
real             :: smois(nx,ny,4,ens_num)
real             :: tslb(nx,ny,4,ens_num)
real             :: tsk(nx,ny,ens_num) !!tsk
integer          :: nn,ii,iis
real             :: temp1,temp2,temp3,temp4,temp5

temp1=0
temp2=0
temp3=0

do nn=1,ens_num

ii=0
iis=0

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=u(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=v(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=w(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=ph(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=t(i,j,k,nn)
enddo
enddo
enddo


do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=qv(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=qr(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=qi(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=qs(i,j,k,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=qgr(i,j,k,nn)
enddo
enddo
enddo
!print*,'ismois',ii+1
do k=1,4
!do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=smois(i,j,k,nn)
enddo
enddo
enddo
!print*,'itslb',ii+1

do k=1,4
!do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=tslb(i,j,k,nn)
enddo
enddo
enddo

!do k=1,4
!do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii,nn)=tsk(i,j,nn)
enddo
enddo
!enddo
!print*,ii,'tsk(1,1,nn)',tsk(1,1,1)
enddo
!print*,'ii,xstat(4046001,1),tsk(1,1,1),enssmois(1,1,1,1),enstslb(1,1,1,1)',xstat(3920001,1),smois(1,1,1,1),tslb(1,1,1,1)
!print*,'ii,xstat(4032001,1),tsk(1,1,1)',ii,xstat(4032001,1),tsk(1,1,1)
if(ii==numstat) print*,'fill xstat & xstatsoil finished normally'

end subroutine fill_xstat

subroutine fill_xstat_mean(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine fill_xstat
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2009.01.04  
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: write ensemble mean into vetor X
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
integer          :: nx,ny,nz,analysis_var_num,ens_num,numstat 
real             :: xstat(numstat)
real             :: u(nx,ny,nz)
real             :: v(nx,ny,nz)
real             :: w(nx,ny,nz)
real             :: ph(nx,ny,nz)
real             :: t(nx,ny,nz)
real             :: qv(nx,ny,nz)
real             :: qr(nx,ny,nz)
real             :: qi(nx,ny,nz)
real             :: qs(nx,ny,nz)
real             :: qgr(nx,ny,nz)
real             :: smois(nx,ny,nz)
real             :: tslb(nx,ny,nz)
integer          :: nn,ii
real             :: temp1,temp2,temp3,temp4,temp5

temp1=0
temp2=0
temp3=0

do nn=1,ens_num

ii=0

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=u(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=v(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=w(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=ph(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=t(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=qv(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=qr(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=qi(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=qs(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=qgr(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=smois(i,j,k)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   xstat(ii)=tslb(i,j,k)
enddo
enddo
enddo

enddo

if(ii==numstat) print*,'fill xstat finished normally'

end subroutine fill_xstat_mean

subroutine dis_xstat(xstat,ens_num,numstat,nx,ny,nz,analysis_var_num,u,v,w,ph,t,qv,qr,qi,qs,qgr,smois,tslb,tsk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine fill_xstat
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: write vetor X into ensemble
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
integer          :: nx,ny,nz,analysis_var_num,ens_num,numstat 
real             :: xstat(numstat,ens_num)
real             :: u(nx,ny,nz,ens_num)
real             :: v(nx,ny,nz,ens_num)
real             :: w(nx,ny,nz,ens_num)
real             :: ph(nx,ny,nz,ens_num)
real             :: t(nx,ny,nz,ens_num)
real             :: qv(nx,ny,nz,ens_num)
real             :: qr(nx,ny,nz,ens_num)
real             :: qi(nx,ny,nz,ens_num)
real             :: qs(nx,ny,nz,ens_num)
real             :: qgr(nx,ny,nz,ens_num)
real             :: smois(nx,ny,4,ens_num)
real             :: tslb(nx,ny,4,ens_num)
real             :: tsk(nx,ny,ens_num)
integer          :: nn,ii
real             :: temp1,temp2,temp3,temp4,temp5

temp1=0
temp2=0
temp3=0

do nn=1,ens_num

ii=0

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   u(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   v(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   w(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   ph(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   t(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   qv(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   qr(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   qi(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx

   ii=ii+1
   qs(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,nz
do j=1,ny
do i=1,nx
   ii=ii+1
   qgr(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,4
do j=1,ny
do i=1,nx
   ii=ii+1
   smois(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do k=1,4
do j=1,ny
do i=1,nx
   ii=ii+1
   tslb(i,j,k,nn)=xstat(ii,nn)
enddo
enddo
enddo

do j=1,ny
do i=1,nx
   ii=ii+1
   tsk(i,j,nn)=xstat(ii,nn)  !!tsk
enddo
enddo

enddo

if(ii==numstat) print*,'distributed xstat finished normally'

end subroutine dis_xstat

subroutine transpose_cal(xstat,xtp,ens_num,numstat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine transpose_cal
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: transpose X
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
integer          :: ens_num,numstat 
real             :: xstat(numstat,ens_num)
real             :: xtp(ens_num,numstat)

do i=1,numstat
do j=1,ens_num
  xtp(j,i)=xstat(i,j)
enddo
enddo  

end subroutine transpose_cal

subroutine ens2stat(egvt_stat,egvt_ens,numstat,mnl,er,xstat,xtp,ens_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine ens2stat
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: transform vector from ensemble space to state space
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

integer          :: numstat,mnl,ens_num 
real             :: xstat(numstat,ens_num)
real             :: xtp(ens_num,numstat)
real             :: er(mnl,4)
real             :: egvt_ens(ens_num,mnl)
real             :: egvt_stat(numstat,mnl)

egvt_stat=0
do j=1,mnl
do i=1,numstat
  do ii=1,mnl
     egvt_stat(i,j)=egvt_stat(i,j)+xstat(i,ii)*egvt_ens(ii,j)
  enddo
     egvt_stat(i,j)=egvt_stat(i,j)/sqrt(abs(er(j,1)))
enddo
enddo  

end subroutine ens2stat



subroutine maxmin(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,                  &
                  nx,ny,nz,ens_num,                                                                           &
                  umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean       &
                  )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine maxmin
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: check whether the max and min is normal 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                   
integer          :: nx,ny,nz 
integer          :: ens_num
real             :: ensu(nx,ny,nz,ens_num)
real             :: ensv(nx,ny,nz,ens_num)
real             :: ensw(nx,ny,nz,ens_num)
real             :: ensph(nx,ny,nz,ens_num)
real             :: enst(nx,ny,nz,ens_num)
real             :: ensqv(nx,ny,nz,ens_num)
real             :: ensqr(nx,ny,nz,ens_num)
real             :: ensqi(nx,ny,nz,ens_num)
real             :: ensqs(nx,ny,nz,ens_num)
real             :: ensqgr(nx,ny,nz,ens_num)
real             :: enssmois(nx,ny,4,ens_num)
real             :: enstslb(nx,ny,4,ens_num)

real             :: umean(nx,ny,nz)
real             :: vmean(nx,ny,nz)
real             :: wmean(nx,ny,nz)
real             :: phmean(nx,ny,nz)
real             :: tmean(nx,ny,nz)
real             :: qvmean(nx,ny,nz)
real             :: qrmean(nx,ny,nz)
real             :: qimean(nx,ny,nz)
real             :: qsmean(nx,ny,nz)
real             :: qgrmean(nx,ny,nz)
real             :: smoismean(nx,ny,4)
real             :: tslbmean(nx,ny,4)

real             :: var_max,var_min


 call maxmin_cal(umean,nx,ny,nz,var_max,var_min)
print*,'U max min',var_max,var_min
 call maxmin_cal(vmean,nx,ny,nz,var_max,var_min)
print*,'V max min',var_max,var_min
 call maxmin_cal(wmean,nx,ny,nz,var_max,var_min)
print*,'W max min',var_max,var_min
 call maxmin_cal(phmean,nx,ny,nz,var_max,var_min)
print*,'PH max min',var_max,var_min
 call maxmin_cal(tmean,nx,ny,nz,var_max,var_min)
print*,'T max min',var_max,var_min
 call maxmin_cal(qvmean,nx,ny,nz,var_max,var_min)
print*,'QV max min',var_max,var_min
 call maxmin_cal(qrmean,nx,ny,nz,var_max,var_min)
print*,'QRAIB max min',var_max,var_min
 call maxmin_cal(qimean,nx,ny,nz,var_max,var_min)
print*,'QICE max min',var_max,var_min
 call maxmin_cal(qsmean,nx,ny,nz,var_max,var_min)
print*,'QSNOW max min',var_max,var_min
 call maxmin_cal(qgrmean,nx,ny,nz,var_max,var_min)
print*,'QGRAUP max min',var_max,var_min
 call maxmin_callda(smoismean,nx,ny,nz,var_max,var_min)
print*,'smois max min',var_max,var_min
 call maxmin_callda(tslbmean,nx,ny,nz,var_max,var_min)
print*,'tslb max min',var_max,var_min

end subroutine maxmin

subroutine maxmin_cal(var,nx,ny,nz,var_max,var_min)

implicit none

integer          :: nx,ny,nz 
integer          :: i,j,k,ii,jj,kk,nn
real             :: var(nx,ny,nz)
real             :: var_max,var_min

var_max=0
var_min=0

do k=1,nz
do j=1,ny
do i=1,nx
  if(var_max <= var(i,j,k)) var_max = var(i,j,k)
  if(var_min >= var(i,j,k)) var_min = var(i,j,k)
enddo
enddo
enddo


end subroutine maxmin_cal

subroutine maxmin_callda(var,nx,ny,nz,var_max,var_min)

implicit none

integer          :: nx,ny,nz 
integer          :: i,j,k,ii,jj,kk,nn
real             :: var(nx,ny,4)
real             :: var_max,var_min

var_max=0
var_min=0

do k=1,4
do j=1,ny
do i=1,nx
  if(var_max <= var(i,j,k)) var_max = var(i,j,k)
  if(var_min >= var(i,j,k)) var_min = var(i,j,k)
enddo
enddo
enddo

print*,var(i,j,k)

end subroutine maxmin_callda

subroutine rms_convetive_region(utrue,vtrue,wtrue,phtrue,ttrue,qvtrue,qrtrue,qitrue,qstrue,qgrtrue,smoistrue,tslbtrue,                &
                                nx,ny,nz,ens_num,                                                                                     &
                                umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,                &
                                stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,                            &
                                stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,stdsmois2d,stdtslb2d,    &
                                ref_x                                                                                                 &
                                )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine rms_convetive_region
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: rmse in convetive region
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
integer          :: nx,ny,nz 
integer          :: ens_num
real             :: utrue(nx,ny,nz)
real             :: vtrue(nx,ny,nz)
real             :: wtrue(nx,ny,nz)
real             :: phtrue(nx,ny,nz)
real             :: ttrue(nx,ny,nz)
real             :: qvtrue(nx,ny,nz)
real             :: qrtrue(nx,ny,nz)
real             :: qitrue(nx,ny,nz)
real             :: qstrue(nx,ny,nz)
real             :: qgrtrue(nx,ny,nz)
real             :: smoistrue(nx,ny,nz)
real             :: tslbtrue(nx,ny,nz)

real             :: umean(nx,ny,nz)
real             :: vmean(nx,ny,nz)
real             :: wmean(nx,ny,nz)
real             :: phmean(nx,ny,nz)
real             :: tmean(nx,ny,nz)
real             :: qvmean(nx,ny,nz)
real             :: qrmean(nx,ny,nz)
real             :: qimean(nx,ny,nz)
real             :: qsmean(nx,ny,nz)
real             :: qgrmean(nx,ny,nz)
real             :: smoismean(nx,ny,nz)
real             :: tslbmean(nx,ny,nz)

real             :: ustd(nx,ny,nz)
real             :: vstd(nx,ny,nz)
real             :: wstd(nx,ny,nz)
real             :: phstd(nx,ny,nz)
real             :: tstd(nx,ny,nz)
real             :: qvstd(nx,ny,nz)
real             :: qrstd(nx,ny,nz)
real             :: qistd(nx,ny,nz)
real             :: qsstd(nx,ny,nz)
real             :: qgrstd(nx,ny,nz)
real             :: smoisstd(nx,ny,nz)
real             :: tslbstd(nx,ny,nz)

real             :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb
real             :: stdu2d(nz),stdv2d(nz),stdw2d(nz),stdph2d(nz),stdt2d(nz),stdqv2d(nz),stdqr2d(nz),stdqi2d(nz),stdqs2d(nz),stdqgr2d(nz),stdsmois2d(nz),stdtslb2d(nz)

real             :: ref_x(nx,ny,nz)


ustd=0
vstd=0
wstd=0
phstd=0
tstd=0
qvstd=0
qrstd=0 
qistd=0 
qsstd=0
qgrstd=0
smoisstd=0
tslbstd=0

stdu=0
stdv=0
stdw=0
stdph=0
stdt=0
stdqv=0
stdqr=0
stdqi=0
stdqs=0
stdqgr=0
stdsmois=0
stdtslb=0

  do k=1,nz 
  do j=1,ny
  do i=1,nx

     ustd(i,j,k)  = ustd(i,j,k)  + ( utrue(i,j,k)  - umean(i,j,k)  )**2
     vstd(i,j,k)  = vstd(i,j,k)  + ( vtrue(i,j,k)  - vmean(i,j,k)  )**2
     wstd(i,j,k)  = wstd(i,j,k)  + ( wtrue(i,j,k)  - wmean(i,j,k)  )**2
     phstd(i,j,k) = phstd(i,j,k) + ( phtrue(i,j,k) - phmean(i,j,k) )**2
     tstd(i,j,k)  = tstd(i,j,k)  + ( ttrue(i,j,k)  - tmean(i,j,k)  )**2
     qvstd(i,j,k) = qvstd(i,j,k) + ( qvtrue(i,j,k) - qvmean(i,j,k) )**2
     qrstd(i,j,k) = qrstd(i,j,k) + ( qstrue(i,j,k) - qrmean(i,j,k) )**2
     qistd(i,j,k) = qistd(i,j,k) + ( qitrue(i,j,k) - qimean(i,j,k) )**2
     qsstd(i,j,k) = qsstd(i,j,k) + ( qstrue(i,j,k) - qsmean(i,j,k) )**2
     qgrstd(i,j,k)= qgrstd(i,j,k)+ ( qgrtrue(i,j,k)- qgrmean(i,j,k))**2
     smoisstd(i,j,k) = smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
     tslbstd(i,j,k) = tslbstd(i,j,k) + ( tslbtrue(i,j,k) - tslbmean(i,j,k) )**2

   
   if(ref_x(i,j,k) >= 10 ) then 
     stdu=stdu+ustd(i,j,k)
     stdv=stdv+vstd(i,j,k)
     stdw=stdw+wstd(i,j,k)
     stdph=stdph+phstd(i,j,k)
     stdt=stdt+tstd(i,j,k)
     stdqv=stdqv+qvstd(i,j,k)   
     stdqr=stdqr+qrstd(i,j,k) 
     stdqi=stdqi+qistd(i,j,k) 
     stdqs=stdqs+qsstd(i,j,k) 
     stdqgr=stdqgr+qgrstd(i,j,k)
     stdsmois=stdsmois+smoisstd(i,j,k) 
     stdtslb=stdtslb+tslbstd(i,j,k)
     
     temp1=temp1+1 
   endif
        
  enddo
  enddo
  enddo
    
  stdu=sqrt(stdu/(temp1))
  stdv=sqrt(stdv/(temp1))
  stdw=sqrt(stdw/(temp1))
  stdph=sqrt(stdph/(temp1))
  stdt=sqrt(stdt/(temp1))
  stdqv=sqrt(stdqv/(temp1))  
  stdqr=sqrt(stdqr/(temp1))  
  stdqi=sqrt(stdqi/(temp1))  
  stdqs=sqrt(stdqs/(temp1))    
  stdqgr=sqrt(stdqgr/(temp1))
  stdsmois=sqrt(stdsmois/(temp1))
  stdtslb=sqrt(stdtslb/(temp1)) 
  
  print*,'the number of grids of REF > 10 is :',temp1
  
ustd=0
vstd=0
wstd=0
phstd=0
tstd=0
qvstd=0    
qrstd=0 
qistd=0 
qsstd=0 
qgrstd=0
smoisstd=0
tslbstd=0 

stdu2d=0
stdv2d=0
stdw2d=0
stdph2d=0
stdt2d=0
stdqv2d=0
stdqr2d=0
stdqi2d=0
stdqs2d=0
stdqgr2d=0
stdsmois2d=0
stdtslb2d=0

  do k=1,nz 
  
    temp1=0
    
  do j=1,ny
  do i=1,nx
  
     ustd(i,j,k)  = ustd(i,j,k)  + ( utrue(i,j,k)  - umean(i,j,k)  )**2
     vstd(i,j,k)  = vstd(i,j,k)  + ( vtrue(i,j,k)  - vmean(i,j,k)  )**2
     wstd(i,j,k)  = wstd(i,j,k)  + ( wtrue(i,j,k)  - wmean(i,j,k)  )**2
     phstd(i,j,k) = phstd(i,j,k) + ( phtrue(i,j,k) - phmean(i,j,k) )**2
     tstd(i,j,k)  = tstd(i,j,k)  + ( ttrue(i,j,k)  - tmean(i,j,k)  )**2
     qvstd(i,j,k) = qvstd(i,j,k) + ( qvtrue(i,j,k) - qvmean(i,j,k) )**2
     qrstd(i,j,k) = qrstd(i,j,k) + ( qstrue(i,j,k) - qrmean(i,j,k) )**2
     qistd(i,j,k) = qistd(i,j,k) + ( qitrue(i,j,k) - qimean(i,j,k) )**2
     qsstd(i,j,k) = qsstd(i,j,k) + ( qstrue(i,j,k) - qsmean(i,j,k) )**2
     qgrstd(i,j,k)= qgrstd(i,j,k)+ ( qgrtrue(i,j,k)- qgrmean(i,j,k))**2
     smoisstd(i,j,k) = smoisstd(i,j,k) + ( smoistrue(i,j,k) - smoismean(i,j,k) )**2
     tslbstd(i,j,k) = tslbstd(i,j,k) + ( tslbtrue(i,j,k) - tslbmean(i,j,k) )**2
    
    if(ref_x(i,j,k) >= 10 ) then 
     stdu2d(k)=stdu2d(k)+ustd(i,j,k)
     stdv2d(k)=stdv2d(k)+vstd(i,j,k)
     stdw2d(k)=stdw2d(k)+wstd(i,j,k)
     stdph2d(k)=stdph2d(k)+phstd(i,j,k)
     stdt2d(k)=stdt2d(k)+tstd(i,j,k)
     stdqv2d(k)=stdqv2d(k)+qvstd(i,j,k)   
     stdqr2d(k)=stdqr2d(k)+qrstd(i,j,k) 
     stdqi2d(k)=stdqi2d(k)+qistd(i,j,k) 
     stdqs2d(k)=stdqs2d(k)+qsstd(i,j,k) 
     stdqgr2d(k)=stdqgr2d(k)+qgrstd(i,j,k)
     stdsmois2d(k)=stdsmois2d(k)+smoisstd(i,j,k)
     stdtslb2d(k)=stdtslb2d(k)+tslbstd(i,j,k)               
     
     temp1=temp1+1 
    endif
    
  enddo
  enddo
 
 if( temp1 .ne. 0 ) then  
  stdu2d(k)=sqrt(stdu2d(k)/(temp1))
  stdv2d(k)=sqrt(stdv2d(k)/(temp1))
  stdw2d(k)=sqrt(stdw2d(k)/(temp1))
  stdph2d(k)=sqrt(stdph2d(k)/(temp1))
  stdt2d(k)=sqrt(stdt2d(k)/(temp1))
  stdqv2d(k)=sqrt(stdqv2d(k)/(temp1)) 
  stdqr2d(k)=sqrt(stdqr2d(k)/(temp1)) 
  stdqi2d(k)=sqrt(stdqi2d(k)/(temp1)) 
  stdqs2d(k)=sqrt(stdqs2d(k)/(temp1)) 
  stdqgr2d(k)=sqrt(stdqgr2d(k)/(temp1))
  stdsmois2d(k)=sqrt(stdsmois2d(k)/(temp1))
  stdtslb2d(k)=sqrt(stdtslb2d(k)/(temp1)) 
 else
  stdu2d(k)=-99999
  stdv2d(k)=-99999
  stdw2d(k)=-99999
  stdph2d(k)=-99999
  stdt2d(k)=-99999
  stdqv2d(k)=-99999
  stdqr2d(k)=-99999
  stdqi2d(k)=-99999
  stdqs2d(k)=-99999
  stdqgr2d(k)=-99999 
  stdsmois2d(k)=-99999
  stdtslb2d(k)=-99999
  endif  
  
  enddo
  
end subroutine rms_convetive_region  

subroutine ens_spread_convetive_region(ensu,ensv,ensw,ensph,enst,ensqv,ensqr,ensqi,ensqs,ensqgr,enssmois,enstslb,                 &
                                       nx,ny,nz,ens_num,                                                                          &
                                       umean,vmean,wmean,phmean,tmean,qvmean,qrmean,qimean,qsmean,qgrmean,smoismean,tslbmean,     &
                                       spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,     &
                                       ref_x                                                                   & 
                                       )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine ens_spread_convetive_region
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: ensemble spread in convetive region
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer          :: nx,ny,nz 
integer          :: ii,jj,kk,nn
integer          :: ens_num
real             :: ensu(nx,ny,nz,ens_num)
real             :: ensv(nx,ny,nz,ens_num)
real             :: ensw(nx,ny,nz,ens_num)
real             :: ensph(nx,ny,nz,ens_num)
real             :: enst(nx,ny,nz,ens_num)
real             :: ensqv(nx,ny,nz,ens_num)
real             :: ensqr(nx,ny,nz,ens_num)
real             :: ensqi(nx,ny,nz,ens_num)
real             :: ensqs(nx,ny,nz,ens_num)
real             :: ensqgr(nx,ny,nz,ens_num)
real             :: enssmois(nx,ny,nz,ens_num)
real             :: enstslb(nx,ny,nz,ens_num)
real             :: umean(nx,ny,nz)
real             :: vmean(nx,ny,nz)
real             :: wmean(nx,ny,nz)
real             :: phmean(nx,ny,nz)
real             :: tmean(nx,ny,nz)
real             :: qvmean(nx,ny,nz)
real             :: qrmean(nx,ny,nz)
real             :: qimean(nx,ny,nz)
real             :: qsmean(nx,ny,nz)
real             :: qgrmean(nx,ny,nz)
real             :: smoismean(nx,ny,nz)
real             :: tslbmean(nx,ny,nz)

real             :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb
real             :: gspd_u(ens_num),gspd_v(ens_num),gspd_w(ens_num),gspd_ph(ens_num),gspd_t(ens_num),gspd_qv(ens_num), &
                    gspd_qr(ens_num),gspd_qi(ens_num),gspd_qs(ens_num),gspd_qgr(ens_num),gspd_smois(ens_num),gspd_tslb(ens_num)
                    
real             :: ref_x(nx,ny,nz)   
integer          :: tem1                 

do nn=1,ens_num
gspd_u(nn)=0
gspd_v(nn)=0
gspd_w(nn)=0
gspd_ph(nn)=0
gspd_t(nn)=0
gspd_qv(nn)=0

gspd_qr(nn)=0
gspd_qi(nn)=0
gspd_qs(nn)=0
gspd_qgr(nn)=0
gspd_smois(nn)=0
gspd_tslb(nn)=0
enddo

spd_u=0
spd_v=0
spd_w=0
spd_ph=0
spd_t=0
spd_qv=0
spd_qr=0
spd_qi=0
spd_qs=0
spd_qgr=0
spd_smois=0
spd_tslb=0

do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_u(nn)=gspd_u(nn)+(ensu(i,j,k,nn)-umean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_u=spd_u+sqrt(gspd_u(nn)/tem1)/ens_num
enddo



do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_v(nn)=gspd_v(nn)+(ensv(i,j,k,nn)-vmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_v=spd_v+sqrt(gspd_v(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_w(nn)=gspd_w(nn)+(ensw(i,j,k,nn)-wmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_w=spd_w+sqrt(gspd_w(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_ph(nn)=gspd_ph(nn)+(ensph(i,j,k,nn)-phmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_ph=spd_ph+sqrt(gspd_ph(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_t(nn)=gspd_t(nn)+(enst(i,j,k,nn)-tmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_t=spd_t+sqrt(gspd_t(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_qv(nn)=gspd_qv(nn)+(ensqv(i,j,k,nn)-qvmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_qv=spd_qv+sqrt(gspd_qv(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_qr(nn)=gspd_qr(nn)+(ensqr(i,j,k,nn)-qrmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_qr=spd_qr+sqrt(gspd_qr(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_qi(nn)=gspd_qi(nn)+(ensqi(i,j,k,nn)-qimean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_qi=spd_qi+sqrt(gspd_qi(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_qs(nn)=gspd_qs(nn)+(ensqs(i,j,k,nn)-qsmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_qs=spd_qs+sqrt(gspd_qs(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_qgr(nn)=gspd_qgr(nn)+(ensqgr(i,j,k,nn)-qgrmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_qgr=spd_qgr+sqrt(gspd_qgr(nn)/tem1)/ens_num
enddo


do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_smois(nn)=gspd_smois(nn)+(enssmois(i,j,k,nn)-smoismean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_smois=spd_smois+sqrt(gspd_smois(nn)/tem1)/ens_num
enddo

do nn=1,ens_num
tem1=0
do k=1,nz
do j=1,ny
do i=1,nx
  if(ref_x(i,j,k)>=10) then  
    gspd_tslb(nn)=gspd_tslb(nn)+(enstslb(i,j,k,nn)-tslbmean(i,j,k))**2
    tem1=tem1+1
  endif
enddo
enddo
enddo
  spd_tslb=spd_tslb+sqrt(gspd_tslb(nn)/tem1)/ens_num
enddo

end subroutine ens_spread_convetive_region

subroutine output_rmse(stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb,                              &
                       stdu2d,stdv2d,stdw2d,stdph2d,stdt2d,stdqv2d,stdqr2d,stdqi2d,stdqs2d,stdqgr2d,stdsmois2d,stdtslb2d,      &
                       lab_time_num,lab_domain,ireg,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine output_rmse
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: output rmse
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
integer               :: nx,ny,nz
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
integer               :: ireg
real                  :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb
real                  :: stdu2d(nz),stdv2d(nz),stdw2d(nz),stdph2d(nz),stdt2d(nz),stdqv2d(nz),  &
                         stdqr2d(nz),stdqi2d(nz),stdqs2d(nz),stdqgr2d(nz),stdsmois2d(nz),stdtslb2d(nz)

  print*,'writing the mean and rms'
  
  if(ireg == 0 ) then
  open(8553,file='rms_d'//lab_domain//'_'//lab_time_num//'_forecast.dat')
  else
  open(8553,file='rms_d'//lab_domain//'_'//lab_time_num//'_analysis.dat')
  endif
  
  write(8553,*) 'rms of the ensemble'
  write(8553,*) 'stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb'

     write(8553,'(10F20.8)') stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb

  write(8553,*) 'vertical rms pattern of the ensemble'
  write(8553,*) 'stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtslb'

    write(8553,*) 'vertical rms pattern of the ensemble '
    do k=nz,1,-1
      write(8553,'(10F20.8)') stdu2d(k),stdv2d(k),stdw2d(k),stdph2d(k),stdt2d(k),stdqv2d(k),stdqr2d(k),stdqi2d(k),stdqs2d(k),stdqgr2d(k),stdsmois2d(k),stdtslb2d(k)
    enddo
  close(8553)
  
end subroutine output_rmse

subroutine output_spd(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk,&
                      lab_time_num,lab_domain,ireg,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine output_spd
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: output spread
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        
integer               :: nx,ny,nz                      
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
integer               :: ireg
real                  :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk
  
! print*,lab_domain,'fuck'
! print*,lab_time_num,'ok?'

  if(ireg == 0 ) then
  open(8553,file='./ens_spread_d'//lab_domain//'_'//lab_time_num//'_forecast.dat')   
  else
  open(8553,file='./ens_spread_d'//lab_domain//'_'//lab_time_num//'_analysis.dat') 
  endif       
  
!  print*,'spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk'
!  print*, spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk

write(8553,'(a17)') 'ensemble spread:'
write(8553,'(a96)') 'spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk'    
write(8553,'(13F11.6)')spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,spd_tsk
  close(8553)       
  
end subroutine output_spd

subroutine output_ratio_r2s(spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb,    &
                            stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdsmois,stdtlsb,               &
                            lab_time_num,lab_domain,ireg,nx,ny,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!          subroutine output_ratio_r2s
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   AUTHOR:Wang Shizhang  2008.12.18   
!!   email : brandenburgii@hotmail.com
!!   Guided by Min Jinzhong (NUIST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   purpose: output  ratio of rmse to spread
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                         
integer               :: nx,ny,nz                      
character (len=2)     :: lab_domain
character (len=2)     :: lab_time_num
integer               :: ireg
real                  :: stdu,stdv,stdw,stdph,stdt,stdqv,stdqr,stdqi,stdqs,stdqgr,stdmois,stdtslb
real                  :: spd_u,spd_v,spd_w,spd_ph,spd_t,spd_qv,spd_qr,spd_qi,spd_qs,spd_qgr,spd_smois,spd_tslb
  
  if(ireg == 0 ) then
  open(8553,file='./radio_r2s_d'//lab_domain//'_'//lab_time_num//'_forecast.dat')   
  else
  open(8553,file='./radio_r2s_d'//lab_domain//'_'//lab_time_num//'_analysis.dat') 
  endif       
  
  write(8553,'(10A8)') 'U','V','W','PH','T','QVAPOR','QRAIN','QICE','QSNOW','QGRAUP'
  write(8553,'(10F8.5)') stdu/spd_u,stdv/spd_v,stdw/spd_w,stdph/spd_ph,stdt/spd_t,            &
                         stdqv/spd_qv,stdqr/spd_qr,stdqi/spd_qi,stdqs/spd_qs,stdqgr/spd_qgr,stdsmois/spd_smois,stdtslb/spd_tslb
  close(8553)     
  
end subroutine output_ratio_r2s       
