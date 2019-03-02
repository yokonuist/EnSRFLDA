!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
  FUNCTION GASDEV(IDUM)

!  PURPOSE:
!
!  Generates a random number from normal distribution by feeding
!  a negative integer iseed.
!
!  INPUT :
!
!    IDUM      an arbitrary negative integer as a seed for a
!              sequence of random numbers
!
!  OUTPUT:
!
!    GASDEV    A random number from Gaussian distribution with mean of 0
!              and standard deviation of 1.
!


   IMPLICIT NONE        ! Force explicit declarations

   INTEGER :: IDUM        ! The seed for random number generation
   REAL :: GASDEV         ! The function to generate random number.

   INTEGER,SAVE::ISET
   REAL,SAVE::GSET
   REAL :: V1, V2, R
   REAL :: RAN1
   REAL :: FAC
   DATA ISET/0/

  IF (ISET.EQ.0) THEN
  1 V1=2.*RAN1(IDUM)-1.
    V2=2.*RAN1(IDUM)-1.
    R=V1**2+V2**2
    IF(R.GE.1.)GO TO 1
    FAC=SQRT(-2.*LOG(R)/R)
    GSET=V1*FAC
    GASDEV=V2*FAC
    ISET=1
  ELSE
    GASDEV=GSET
    ISET=0
  ENDIF

  RETURN
  END FUNCTION GASDEV
!-------------------------------------------------------------------------------------------  
  FUNCTION RAN1(IDUM)

!  PURPOSE:
!
!  Generates a random number between 0 and 1 by feeding
!  a negative integer iseed.
!
!  Added by M.Tong
!  Reference: "Seminumerical Algorithms" by Donald Knuth
!
!
!  INPUT :
!
!    IDUM      an arbitrary negative integer as a seed for a
!              sequence of random numbers
!
!  OUTPUT:
!
!    RAN1      A random number between 0 and 1.
!
!
  IMPLICIT NONE        ! Force explicit declarations
  INTEGER :: IDUM        ! The seed for random number generation
  REAL :: RAN1           ! The function to generate random number.

  REAL,SAVE :: R(97)
  INTEGER :: IX1,IX2,IX3,J,IFF
  INTEGER :: M1,M2,M3,IA1,IA2,IA3,IC1,IC2,IC3
  REAL :: RM1,RM2
  SAVE IX1,IX2,IX3

  PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
  PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
  PARAMETER (M3=243000,IA3=4561,IC3=51349)
  DATA IFF /0/
!----------------------------------------------------------------------
!
!  Initialize the sequence of random numbers between 0 and 1,
!  using iseed.
!
!----------------------------------------------------------------------
!
  IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
    IFF=1
    IX1=MOD(IC1-IDUM,M1)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX2=MOD(IX1,M2)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX3=MOD(IX1,M3)
    DO J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
    ENDDO
    IDUM=1
  ENDIF
  IX1=MOD(IA1*IX1+IC1,M1)
  IX2=MOD(IA2*IX2+IC2,M2)
  IX3=MOD(IA3*IX3+IC3,M3)
  J=1+(97*IX3)/M3
  IF(J.GT.97.OR.J.LT.1)THEN
    WRITE(*,*)'J is greater than 97 or less than 1','IDUM=',IDUM
    STOP
  ENDIF
  RAN1=R(J)
  R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1

  RETURN
  END FUNCTION RAN1
!-------------------------------------------------------------------------------------------  
!------------------------------------------------------------------------------------------- 
!------------------------------------------------------------------------------------------- 
!------------------------------------------------------------------------------------------- 
!====================================================================
      SUBROUTINE ranpert(iseed,n,rndn)
!====================================================================

      IMPLICIT NONE
!====================================================================

      INTEGER :: iseed            ! random seed
      INTEGER :: n                ! dimension size
      INTEGER :: i

      REAL :: rndn   (n)          ! return standard normal distri. 
      REAL, ALLOCATABLE :: work1(:),work2(:)
      REAL, PARAMETER   :: pi=3.141592653589
      REAL :: ave,var

      ALLOCATE (work1(n),work2(n))

      do i=1,n
        iseed =  mod(iseed*7141+54773,259200)
        work1(i)=float(iseed)/259199.
      enddo
      do i=1,n
        iseed =  mod(iseed*7141+54773,259200)
        work2(i)=float(iseed)/259199.
      enddo

      rndn = sqrt(-2.0*log(work1+tiny(1.0)))*cos(2.0*pi*work2)
                       ! N(0,1) standard normal distribution

! statistics
      ave = sum(rndn)/float(n)
      work1 = rndn-ave
      var = dot_product(work1,work1)/float(n-1)

      print *,'  AVE,VAR,SQRT(AVR):',ave,var,sqrt(var)

      DEALLOCATE (work1,work2)
      RETURN
      END SUBROUTINE ranpert
!-------------------------------------------------------------------------------------------  
!------------------------------------------------------------------------------------------- 
!-----------------the below are subroutines for EOF-----------------------------------------
!-------------------------------------------------------------------------------------------       
      subroutine eof(m,n,mnl,f,ks,er,egvt)
!-------------------------------------------------------------------------------------------  
!-------------------------------------------------------------------------------------------  
!     EMPIRICAL ORTHOGONAL FUNCTIONS (EOF's)
!     This subroutine applies the EOF approach to analysis time series 
!       of meteorological field f(m,n).
!     input: m,n,mnl,f(m,n),ks
!       m: number of grid-points
!       n: lenth of time series
!       mnl=min(m,n)
!       f(m,n): the raw spatial-temporal seires
!       ks: contral parameter
!           ks=-1: self; ks=0: depature; ks=1: normalized depature
!     output: egvt,ecof,er
!       egvt(m,mnl): array of eigenvactors
!       ecof(mnl,n): array of time coefficients for the respective eigenvectors
!       er(mnl,1): lamda (eigenvalues), its sequence is from big to small.
!       er(mnl,2): accumulated eigenvalues from big to small
!       er(mnl,3): explained variances (lamda/total explain) from big to small
!       er(mnl,4): accumulated explaned variances from big to small
!-------------------------------------------------------------------------------------------  
!-------------------------------------------------------------------------------------------  
      dimension f(m,n)
      dimension er(mnl,4),egvt(m,mnl)
      dimension cov(mnl,mnl),s(mnl,mnl),d(mnl),v(mnl) !work array
!     Preprocessing data
      call transf(m,n,f,ks)
!     Crossed product matrix of the data f(m,n)
      call crossproduct(m,n,mnl,f,cov)
!     Eigenvalues and eigenvectors by the Jacobi method 
      call jacobi(mnl,cov,s,d,0.00001)
!     Specified eigenvalues 
      call arrang(mnl,d,s,er)
!     Fill data in ouyput eignvectors 
      call tcoeff(m,n,mnl,f,s,er,egvt)
      return
      end subroutine eof
!-------------------------------------------------------------------------------------------  
!----Preprocessing data to provide a field by ks.-------------------------------------------  
!------------------------------------------------------------------------------------------- 
      subroutine transf(m,n,f,ks)
      dimension f(m,n)
      dimension fw(n),wn(m)           !work array
	i0=0
	do i=1,m
	  do j=1,n
          fw(j)=f(i,j)
        enddo
	  call meanvar(n,fw,af,sf,vf)
	  if(sf.eq.0.)then
	    i0=i0+1
	    wn(i0)=i
	  endif
	enddo
	if(i0.ne.0)then
	  write(*,*)'****  FAULT  ****'
	  write(*,*)' The program cannot go on because '
	  write(*,*)' The original field has invalid data.'
	  write(*,*)' There are totally ',i0,'  gridpionts with invalid data.'
     	  write(*,*)' The array WN stores the positions of those invalid'
          write(*,*)' grid-points. You must pick off those invalid data'  
          write(*,*)' from the orignal field and then reinput a new'
	  write(*,*)' field to calculate its EOFs.'   
	  write(*,*)'****  FAULT  ****'
	  stop
	endif	    
      if(ks.eq.-1)return
      if(ks.eq.0)then                !anomaly of f
        do i=1,m
          do j=1,n
            fw(j)=f(i,j)
          enddo
          call meanvar(n,fw,af,sf,vf)
          do j=1,n
            f(i,j)=f(i,j)-af
          enddo
        enddo
        return
      endif
      if(ks.eq.1)then                 !normalizing f
        do i=1,m
          do j=1,n
            fw(j)=f(i,j)
          enddo
          call meanvar(n,fw,af,sf,vf)
          do j=1,n
            f(i,j)=(f(i,j)-af)/sf
          enddo
        enddo
      endif
      return
      end subroutine transf
!-------------------------------------------------------------------------------------------------------------------------------------    
!Crossed product martix of the data. It is n times of covariance matrix of the data if ks=0 (i.e. for anomaly).
!------------------------------------------------------------------------------------------------------------------------------------- 
      subroutine crossproduct(m,n,mnl,f,cov)
      dimension f(m,n),cov(mnl,mnl)
      if(n-m) 10,50,50
  10  do 20 i=1,mnl
      do 20 j=i,mnl
        cov(i,j)=0.0
        do is=1,m
          cov(i,j)=cov(i,j)+f(is,i)*f(is,j)
        enddo
        cov(j,i)=cov(i,j)
  20  continue
      return
  50  do 60 i=1,mnl
      do 60 j=i,mnl
        cov(i,j)=0.0
        do js=1,n
          cov(i,j)=cov(i,j)+f(i,js)*f(j,js)
        enddo
        cov(j,i)=cov(i,j)
  60  continue
      return
      end      
!----------------------------------------------------------------------------------------------------
!Computing eigenvalues and eigenvectors of a real symmetric matrix a(m,m) by the Jacobi method.  
!----------------------------------------------------------------------------------------------------
      subroutine jacobi(m,a,s,d,eps)
      dimension a(m,m)
      dimension s(m,m),d(m)
      do 30 i=1,m
      do 30 j=1,i
        if(i-j) 20,10,20
  10    s(i,j)=1.
        go to 30
  20    s(i,j)=0.
        s(j,i)=0.
  30  continue
      g=0.
      do 40 i=2,m
        i1=i-1
        do 40 j=1,i1
  40      g=g+2.*a(i,j)*a(i,j)
      s1=sqrt(g)
      s2=eps/float(m)*s1
      s3=s1
      l=0
  50  s3=s3/float(m)
  60  do 130 iq=2,m
        iq1=iq-1
        do 130 ip=1,iq1
        if(abs(a(ip,iq)).lt.s3) goto 130
        l=1
        v1=a(ip,ip)
        v2=a(ip,iq)
        v3=a(iq,iq)
        u=0.5*(v1-v3)
        if(u.eq.0.0) g=1.
        if(abs(u).ge.1e-10) g=-sign(1.,u)*v2/sqrt(v2*v2+u*u)
        st=g/sqrt(2.*(1.+sqrt(1.-g*g)))
        ct=sqrt(1.-st*st)
        do 110 i=1,m
          g=a(i,ip)*ct-a(i,iq)*st
          a(i,iq)=a(i,ip)*st+a(i,iq)*ct
          a(i,ip)=g
          g=s(i,ip)*ct-s(i,iq)*st
          s(i,iq)=s(i,ip)*st+s(i,iq)*ct
  110     s(i,ip)=g
        do 120 i=1,m
          a(ip,i)=a(i,ip)
  120     a(iq,i)=a(i,iq)
        g=2.*v2*st*ct
        a(ip,ip)=v1*ct*ct+v3*st*st-g
        a(iq,iq)=v1*st*st+v3*ct*ct+g
        a(ip,iq)=(v1-v3)*st*ct+v2*(ct*ct-st*st)
        a(iq,ip)=a(ip,iq)
  130 continue
      if(l-1) 150,140,150
  140 l=0
      go to 60
  150 if(s3.gt.s2) goto 50
      do 160 i=1,m
        d(i)=a(i,i)
  160 continue
      return
      end
!-------------------------------------------------------------------------------------------  
!-----Provides a series of eigenvalues from maximuim to minimuim.-------------------------  
!-------------------------------------------------------------------------------------------  
      subroutine arrang(mnl,d,s,er)
      dimension d(mnl),s(mnl,mnl),er(mnl,4)
      tr=0.0
      do 10 i=1,mnl
        tr=tr+d(i)
        er(i,1)=d(i)
  10  continue
      mnl1=mnl-1
      do 20 k1=mnl1,1,-1
      do 20 k2=k1,mnl1
        if(er(k2,1).lt.er(k2+1,1)) then
          c=er(k2+1,1)
          er(k2+1,1)=er(k2,1)
          er(k2,1)=c
          do 15 i=1,mnl
            c=s(i,k2+1)
            s(i,k2+1)=s(i,k2)
            s(i,k2)=c
  15      continue
        endif
  20  continue
      er(1,2)=er(1,1)
      do 30 i=2,mnl
        er(i,2)=er(i-1,2)+er(i,1)
  30  continue
      do 40 i=1,mnl
        er(i,3)=er(i,1)/tr
        er(i,4)=er(i,2)/tr
  40  continue
      return
      end
!-------------------------------------------------------------------------------------------  
!-----Provides standard eigenvectors and their time coefficients----------------------------  
!-------------------------------------------------------------------------------------------  
      subroutine tcoeff(m,n,mnl,f,s,er,egvt)
      dimension f(m,n)
      dimension s(mnl,mnl),er(mnl,4),egvt(m,mnl)
      real v(mnl)  !work array
      do j=1,mnl
        do i=1,m
          egvt(i,j)=0.
        enddo
      enddo
!     Normalizing the input eignvectors s
      do 10 j=1,mnl
        c=0.
        do i=1,mnl
          c=c+s(i,j)*s(i,j)
        enddo
        c=sqrt(c)
        do i=1,mnl
          s(i,j)=s(i,j)/c
        enddo
  10  continue
!-------------------------------------------
      if(m.le.n) then
        do js=1,mnl
        do i=1,m
          egvt(i,js)=s(i,js)
        enddo
        enddo
      else
        do 40 i=1,m
          do j=1,n
            v(j)=f(i,j)
          enddo
          do js=1,mnl
          do j=1,n
            egvt(i,js)=egvt(i,js)+v(j)*s(j,js)
          enddo
          enddo
  40    continue
!--------------------------------------------------------  
        do 50 js=1,mnl
          do i=1,m
            egvt(i,js)=egvt(i,js)/sqrt(abs(er(js,1)))
          enddo
  50    continue
      endif
      return
      end
!------------------------------------------------------------------------------------------- 
!------Computing the mean ax, standard deviation sx-----------------------------------------     
!-------------------------------------------------------------------------------------------      
      subroutine meanvar(n,x,ax,sx,vx)
      dimension x(n)
      ax=0.
      vx=0.
      sx=0.
      do 10 i=1,n
        ax=ax+x(i)
  10  continue
      ax=ax/float(n)
      do 20 i=1,n
        vx=vx+(x(i)-ax)**2
  20  continue
      vx=vx/float(n)
      sx=sqrt(vx)
      return
      end            
!-------------------------------------------------------------------------------------------
!------Computing the mean ax, standard deviation sx-----------------------------------------     
!-------------------------------------------------------------------------------------------      
!-------------------------------------------------------------------------------------------  
!-----------------the above are subroutines for EOF----------------------------------------- 
!------------------------------------------------------------------------------------------- 
!-------------------------------------------------------------------------------------------     
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SCALE_BCKGRD              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
SUBROUTINE recurfilt_3d(nx,ny,nz,pgrd,ipass_filt,hradius,nradius_z)
!cc
!cc
  IMPLICIT none
  INTEGER :: i,j,k,n
  INTEGER :: nx,ny,nz,ipass_filt,nradius_z
     REAL :: hradius
  REAL :: alpha, alpha_z, ee
  REAL :: temp,temp2

  REAL :: pgrd(nx,ny,nz)
  REAL, DIMENSION(:), allocatable :: temx
  REAL, DIMENSION(:), allocatable :: temy
  REAL, DIMENSION(:), allocatable :: temz
!
  REAL,   DIMENSION (:), allocatable :: varb
!
  allocate ( temx(nx) )
  allocate ( temy(ny) )
  allocate ( temz(nz) )
!
!
  IF( hradius == 0 .and. nradius_z == 0 ) return 
!
  ee    = REAL(ipass_filt) / (hradius*hradius)
  alpha = 1+ee-SQRT( ee*(ee+2.) )
!
  IF( nradius_z /= 0 ) THEN
    ee     = REAL (ipass_filt) / REAL (nradius_z*nradius_z)
    alpha_z = 1 + ee - SQRT (ee*(ee+2.))
  ENDIF
!
!
  DO n = 1, ipass_filt/2
!
!
    allocate ( varb(nx) )
    DO k = 1, nz
      DO j = 1, ny
!
        if(n==1) varb(1) = (1-alpha) * pgrd(1,j,k)
        if(n==2) varb(1) = (1-alpha)/(1-alpha*alpha) * pgrd(1,j,k)
        if(n==3) then
          temp = (1-alpha)/((1-alpha*alpha)*(1-alpha*alpha))
          temp2 =alpha*alpha*alpha
          varb(1) = temp * (pgrd(1,j,k)-temp2*pgrd(1,j,k))
        ENDIF
        if(n>=4) then
          temp2 =alpha*alpha*alpha
          temp = (1-alpha)/(1-3*alpha*alpha+3*temp2*alpha-temp2*temp2)
          varb(1) = temp * (pgrd(1,j,k)-3*temp2*pgrd(2,j,k)+              &
                 temp2*alpha*alpha*pgrd(2,j,k)+temp2*alpha*pgrd(3,j,k))
        ENDIF
!
        DO i = 2, nx, 1
          varb(i) = alpha*varb(i-1) + (1.-alpha)*pgrd(i,j,k)
        END DO
!
        if(n==0) pgrd(nx,j,k) = (1-alpha) * varb(nx)
        if(n==1) pgrd(nx,j,k) = (1-alpha)/(1-alpha*alpha) * varb(nx)
        if(n==2) then
          temp = (1-alpha)/((1-alpha*alpha)*(1-alpha*alpha))
          temp2 =alpha*alpha*alpha
          pgrd(nx,j,k) = temp * (varb(nx)-temp2*varb(nx-1))
        ENDIF
        if(n>=3) then
          temp2 =alpha*alpha*alpha
          temp = (1-alpha)/(1-3*alpha*alpha+3*temp2*alpha-temp2*temp2)
          pgrd(nx,j,k) = temp * (varb(nx)-3*temp2*varb(nx-1)+            &
               temp2*alpha*alpha*varb(nx-1)+temp2*alpha*varb(nx-2))
        ENDIF
!
!
        DO i = nx-1, 1, -1
          pgrd(i,j,k) = alpha*pgrd(i+1,j,k) + (1.-alpha)*varb(i)
        END DO
!
      END DO
    END DO
    deallocate (varb)
!
!
    allocate ( varb(ny) )
    DO k = 1, nz
      DO i = 1, nx
!
        if(n==1) varb(1) = (1-alpha) * pgrd(i,1,k)
        if(n==2) varb(1) = (1-alpha)/(1-alpha*alpha) * pgrd(i,1,k)
        if(n==3) then
          temp = (1-alpha)/((1-alpha*alpha)*(1-alpha*alpha))
          temp2 =alpha*alpha*alpha
          varb(1) = temp * (pgrd(i,1,k)-temp2*pgrd(i,1,k))
        ENDIF
        if(n>=4) then
          temp2 =alpha*alpha*alpha
          temp = (1-alpha)/(1-3*alpha*alpha+3*temp2*alpha-temp2*temp2)
          varb(1) = temp * (pgrd(i,1,k)-3*temp2*pgrd(i,2,k)+              &
                 temp2*alpha*alpha*pgrd(i,2,k)+temp2*alpha*pgrd(i,3,k))
        ENDIF
!
        DO j = 2, ny, 1
          varb(j) = alpha*varb(j-1) + (1.-alpha)*pgrd(i,j,k)
        END DO
!
        if(n==0) pgrd(i,ny,k) = (1-alpha) * varb(ny)
        if(n==1) pgrd(i,ny,k) = (1-alpha)/(1-alpha*alpha) * varb(ny)
        if(n==2) then
          temp = (1-alpha)/((1-alpha*alpha)*(1-alpha*alpha))
          temp2 =alpha*alpha*alpha
          pgrd(i,ny,k) = temp * (varb(ny)-temp2*varb(ny-1))
        ENDIF
        if(n>=3) then
          temp2 =alpha*alpha*alpha
          temp = (1-alpha)/(1-3*alpha*alpha+3*temp2*alpha-temp2*temp2)
          pgrd(i,ny,k) = temp * (varb(ny)-3*temp2*varb(ny-1)+            &
               temp2*alpha*alpha*varb(ny-1)+temp2*alpha*varb(ny-2))
        ENDIF
!
        DO j = ny-1, 1, -1
          pgrd(i,j,k) = alpha*pgrd(i,j+1,k) + (1.-alpha)*varb(j)
        END DO
!
      END DO
    END DO
    deallocate (varb)
!
!
    IF( nradius_z /= 0 ) THEN
      allocate ( varb(nz) )
      DO i = 1, nx
        DO j = 1, ny

        if(n==1) varb(1) = (1-alpha_z) * pgrd(i,j,1)
        if(n==2) varb(1) = (1-alpha_z)/(1-alpha_z*alpha_z) * pgrd(i,j,1)
        if(n==3) then
          temp = (1-alpha_z)/((1-alpha_z*alpha_z)*(1-alpha_z*alpha_z))
          temp2 =alpha_z*alpha_z*alpha_z
          varb(1) = temp * (pgrd(i,j,1)-temp2*pgrd(i,j,1))
        ENDIF
        if(n>=4) then
          temp2 =alpha_z*alpha_z*alpha_z
          temp = (1-alpha_z)/(1-3*alpha_z*alpha_z+3*temp2*alpha_z-temp2*temp2)
          varb(1) = temp * (pgrd(i,j,1)-3*temp2*pgrd(i,j,2)+              &
                 temp2*alpha_z*alpha_z*pgrd(i,j,2)+temp2*alpha_z*pgrd(i,j,3))
        ENDIF
!
        DO k = 2, nz, 1
          varb(k) = alpha_z*varb(k-1) + (1.-alpha_z)*pgrd(i,j,k)
        END DO
!
        if(n==0) pgrd(i,j,nz) = (1-alpha_z) * varb(nz)
        if(n==1) pgrd(i,j,nz) = (1-alpha_z)/(1-alpha_z*alpha_z) * varb(nz)
        if(n==2) then
          temp = (1-alpha_z)/((1-alpha_z*alpha_z)*(1-alpha_z*alpha_z))
          temp2 =alpha_z*alpha_z*alpha_z
          pgrd(i,j,nz) = temp * (varb(nz)-temp2*varb(nz-1))
        ENDIF
        if(n>=3) then
          temp2 =alpha_z*alpha_z*alpha_z
          temp = (1-alpha_z)/(1-3*alpha_z*alpha_z+3*temp2*alpha_z-temp2*temp2)
          pgrd(i,j,nz) = temp * (varb(nz)-3*temp2*varb(nz-1)+            &
               temp2*alpha_z*alpha_z*varb(nz-1)+temp2*alpha_z*varb(nz-2))
        ENDIF
!
        DO k = nz-1, 1, -1
          pgrd(i,j,k) = alpha_z*pgrd(i,j,k+1) + (1.-alpha_z)*varb(k)
        END DO

        END DO
      END DO
      deallocate (varb)
    ENDIF
!
  END DO
!
!
  deallocate (temx)
  deallocate (temy)
  deallocate (temz)
!
!
  RETURN
END SUBROUTINE recurfilt_3d       


