!---------------------------------------------------------------------------
! Description:
!   FFT library with subroutines for 1D FFT, 2D FFT (master), N-dimensional slave FFT.
!   Auxillary subroutines for writing files of fourier transform and zero padding
!
! Method:
!   Followed Daniel Lanczos algorithm as given in Numerical Recipies
!
! Developer: Ashwin Vishnu M
!   Creative Commons (c) 2013
!
! Info:
!   Only +ve frequency range calculated by master subroutines.
!   Sizes of Input matrix (MxN), Output matrix (M/2 x N/2)
!   FFT calculated till Nyquist limit freq = pi/h
!   Least count of frequency spectrum      = pi/(M/2*h)
!---------------------------------------------------------------------------

      subroutine FFT1d(func, M,  transform)    ! Master subroutine
        integer :: Ndim, i, M, zeropadreq
        integer, dimension (1) :: NN
        double precision:: func(M), transform(M/2)
        double precision, dimension(:),allocatable :: data
        double precision :: theta, pi
        double complex ::  Fr, Fi, Ftrans, W, Wp, Hpos, Hneg

        allocate( data(1:2*M) ) 
        Ndim = 1

         do i= 1,M
            data(2*i-1) = func(i)
            data(2*i  ) = 0d0
         enddo

        NN(1) = M 
        print*, 'Starting FFT... '
        call FOURN(data,NN,Ndim,+1)                ! Call slave
        print*, 'FFT Completed... '
        
        pi=4d0*datan(1d0)
        theta = 2d0*pi/M
        Wp   = dcmplx(dcos(theta),dsin(theta))
        W    = Wp

        Fr = data(1)+data(2)
        Fi = data(1)-data(2)
        Ftrans = Fr - Fi*dcmplx(0.0,1.0) 
        transform(1) =Fr !dsqrt(dreal(Ftrans)**2+ dimag(Ftrans)**2)
        do i= 2, M/2
            Hpos =  dcmplx(data(2*i-1), data(2*i))
            Hneg =  dcmplx(data(2*M-2*i+3), data(2*M-2*i+4))
            Fr   =  (Hpos +dconjg(Hneg))  /2d0
            Fi   =  (Hpos -dconjg(Hneg))  /2d0 * W
            Ftrans = Fr - Fi*dcmplx(0.0,1.0) 
            transform(i) =dsqrt(dreal(Ftrans)**2+ dimag(Ftrans)**2)  
            W = W*Wp
        enddo

        transform = transform*2d0 / M       !Normalization

        deallocate(data)
      endsubroutine FFT1d

      subroutine InvFFT1d(transform, M, func)    ! Master subroutine
        integer :: Ndim, i, M, zeropadreq
        integer, dimension (1) :: NN
        double precision:: func(M), transform(M/2)
        double precision, dimension(:),allocatable :: data
        double precision :: theta, pi
        double complex ::  Fr, Fi, Ftrans, W, Wp, Hpos, Hneg

        allocate( data(1:2*M) ) 
        Ndim = 1

        data(1)=0.5d0*transform(1)
        data(2)=0d0
        do i= 2,M/2
           data(2*i-1) = 0.5d0*transform(i)
           data(2*i  ) = 0d0
           data(2*M-2*i+3) = 0.5d0*transform(i)
           data(2*M-2*i+4) = 0d0
        enddo

        NN(1) =M
 
        print*, 'Starting Inverse FFT... '
        call FOURN(data,NN,Ndim,-1)                ! Call slave
        print*, 'Inverse FFT Completed... '
        
        do i= 1,M
            func(i) = dsqrt(data(2*i-1)**2 + data(2*i)**2)
        enddo
        
        deallocate(data)
      endsubroutine InvFFT1d

      subroutine FFT2d(func, M, N, transform)    ! Master subroutine
        implicit none
        integer :: Ndim, i, j, M, N, zeropadreq
        integer, dimension (2) :: NN
        double precision:: func(M,N) , transform(M/2,N/2) 
        double precision, dimension(:,:),allocatable :: data
        double precision :: theta, pi
        integer :: ipos, jpos, ineg, jneg
        double complex ::  Fr, Fi, Ftrans, W, Wp, Hpos, Hneg

        allocate( data(1:2*M,1:N) ) 
        Ndim = 2

        do j= 1,N
         do i= 1,M
            data(2*i-1,j) = func(i,j)
            data(2*i, j ) = 0d0
         enddo
        enddo

        NN(1) = M; NN(2) = N
        print*, 'Starting FFT... '
        call FOURN(data,NN,Ndim,+1)                ! Call slave
        print*, 'FFT Completed... '

        pi    = 4d0*datan(1d0)
        theta = 2d0*pi/M
        Wp   = dcmplx(dcos(theta),dsin(theta))
        transform = 0d0

        W = dcmplx(1.0,0.0)
        do i= 1, M/2
          ipos  =  2*i
          ineg  =  2*M-2*i+4
          if(i==1) then
                  ipos = 2
                  ineg = 2
          endif
          do j=1,N/2 
            jpos  =  j
            jneg  =  N-j+2
            if(j==1) then
                    jpos = 1
                    jneg = 1
            endif
            Hpos =  dcmplx(data(ipos-1,jpos), data(ipos,jpos))
            Hneg =  dcmplx(data(ineg-1,jneg), data(ineg,jneg))
            Fr   =  (Hpos +conjg(Hneg))  /2d0
            Fi   =  (Hpos -conjg(Hneg))  /2d0 * W
            Ftrans = Fr - Fi*dcmplx(0.0,1.0) 
            transform(i,j) =dsqrt(dreal(Ftrans)**2 + dimag(Ftrans)**2)
          enddo
          W = W*Wp
        enddo
        transform = transform*2d0/(M*N)

        deallocate(data)
      endsubroutine FFT2d


!*************************************************************************************!
!            File writing and zeropadding                                             !
!*************************************************************************************!
subroutine WriteFFT1d(transform,M,h,filename)
      implicit none
      integer :: M ,i
      double precision :: transform(M/2), h, k, pi
      character(50) :: filename
      
      pi = 4d0*datan(1d0)
      open(34, file=trim(filename))
      write(34,*) 'Variables= kh,u_fft'
      write(34,*) 'Zone i=',M/2

      do i=1,M/2
         k = (i-1)*2d0*pi/(M*h)*h
         write(34,*) k, transform(i)
      enddo
      close(34)
endsubroutine WriteFFT1d

subroutine WriteFFT2d(transform,M,N,hx,hy,filename)
      implicit none
      integer :: M,N ,i,j
      double precision :: transform(M/2, N/2), hx,hy, kx, ky, pi
      character(50) :: filename
      
      pi = 4d0*datan(1d0)
      open(34, file=trim(filename))
      write(34,*) 'Variables= kx,ky,u_fft'
      write(34,*) 'Zone i=',M/2,' j=', N/2

      do j=1,N/2
        ky = (j-1)*2d0*pi/(N*hy)
        do i=1,M/2
         kx = (i-1)*2d0*pi/(M*hx)
         write(34,*) kx,ky, transform(i,j)
        enddo
      enddo
      close(34)
endsubroutine WriteFFT2d

subroutine ZeroPad1d(data,M, datanew,Mpad)
     implicit none
     integer :: M,Mpad, i, ii
     double precision :: data(1:M), pow, datanew(1:Mpad)
     
     pow = log(real(M))/log(2.0)
     if( mod(pow,1.0) > 0 ) then
        print*, 'Zero padding... ',M,'-->',Mpad
        print*, 'Normalization factor: ', real(Mpad)/M 
        ii = -(Mpad-M)/2  
        do i = 1,Mpad
           if(ii .ge. 1 .and. ii .le. M ) then
               datanew(i) = data(ii)
           else
               datanew(i) = 0d0
           endif
           ii = ii+1
        enddo

     else
        print*, 'Zero padding not required '
        datanew(1:M) = data(1:M)
     endif
endsubroutine ZeroPad1d

subroutine ZeroPad2d(data,M,N, datanew,Mpad,Npad)
     implicit none
     integer :: M,Mpad, i, ii, N,Npad, j, jj
     double precision :: data(1:M,1:N), pow1, pow2, datanew(1:Mpad,1:Npad)
     
     pow1 = log(real(M))/log(2.0)
     pow2 = log(real(N))/log(2.0)
     datanew = 0d0

     if( mod(pow1,1.0) > 0 .or. mod(pow2,1.0)>0 ) then
        print*, 'Zero padding...',M,'x',N,'-->',Mpad,'x',Npad
        print*, 'Normalization factor x: ', real(Mpad)/M 
        print*, 'Normalization factor y: ', real(Npad)/N 

        jj = -(Npad-N)/2
        do  j=1,Npad
         if(jj .ge. 1 .and. jj .le. N ) then
           ii = -(Mpad-M)/2
           do i = 1,Mpad
           if(ii .ge. 1 .and. ii .le. M )then
             datanew(i,j) = data(ii,jj)
           endif
           ii = ii+1
           enddo
         endif
         jj = jj+1
        enddo

     else
        print*, 'Zero padding not required '
        datanew(1:M,1:N) = data(1:M,1:N)
     endif
endsubroutine ZeroPad2d


!*************************************************************************************!
!             Dealiasing for product terms                                            !
!*************************************************************************************!
subroutine Dealias1d(func,func_sq, M)
        integer :: M, i,j,k, ii, zpadindex
        double precision, dimension( 1:M ) :: func
        double precision, dimension(1:M/2) :: transform
        double precision, dimension(1:M) :: convolution, func_sq
        double precision, dimension(1:2*M) :: func_temp
        double precision :: temp, pi, c
        real :: zpad

        pi = 4d0*datan(1d0)
        call FFT1d(func, M, transform)

     ! Convolution Four(f.g) = F(k)*G(k) = Sum(F(k').G(k-k'))
        convolution = 0d0
        do i=1,M
           do j=1,M/2
              k =  i-j+1
              c=0.5d0           ! c to compensate for negative half of the transform 
              if(k<1) then
                    c=1d0
                    k = -k+2
              endif
              convolution(i) = convolution(i)+ transform(j)*transform(k)*c
           enddo
        enddo

      ! Zeropadding
        zpad = 2  ! 3/2 k_m or 2 k_m
        convolution(int(zpad*M/4)+1:M) =0d0

!      call writefft1d(convolution,2*M,6d0/(2*M-1),'afterzeropad.dat')

       call InvFFT1d(convolution,2*M,func_temp)
       do i=1,M
         func_sq(i) = func_temp(2*i-1)
       enddo

endsubroutine Dealias1d                  

subroutine ShiftAlign(func,func_sq,M)
        implicit none
        integer :: M
        integer :: i,ip,im, peak
        double precision, dimension(M):: func, func_sq, temp

        peak = maxloc(func,1)

        temp(peak) = func_sq(1)

        ip = peak+1
        im = peak-1
        do i = 2, M/2
            if(ip == M) then          !When Periodic B.C. is applied
                   temp(ip)=func_sq(i)
                   ip=1
            endif
            if(im == 1) then
                   temp(im)=func_sq(i)
                   im=M
            endif
            temp(ip) = func_sq(i)
            temp(im) = func_sq(i)
            ip = ip+1
            im = im-1
        enddo
        func_sq = temp
endsubroutine ShiftAlign

!*************************************************************************************!
!             Slave subroutine by N.Brenner. Ref: Numerical recipies                  !
!*************************************************************************************!
      subroutine FOURN(data, NN, Ndim, ISIGN)  
         implicit none
         integer :: Ndim, ISIGN                  ! No of dimensions,ISIGN= +1 FFT, -1 Inv FFT
         double precision, dimension(*) :: data  ! I/P and transformed O/P data array
         integer, dimension(Ndim) :: NN          ! No of sampled points in each dimension
         double complex :: W, Wp, Wtemp, temp
         double precision :: theta, pi
         integer :: N, Ntot, Nprev, Nrem
         integer :: i1, i2, i2rev, i3, i3rev, ip1, ip2, ip3,idim, ibit
         integer :: k1, k2, ifp1, ifp2
         
         pi= 4d0*datan(1d0)
         Ntot = 1     ! Count total no. of data
         do idim = 1, Ndim
           Ntot = Ntot * NN(idim)
         enddo
         Nprev=1

         do idim = 1, Ndim   ! Loop over dimensions
           N = NN(idim)
           print*, 'Sweeping along dimension #',idim, N, 'points'
           Nrem =  Ntot /(N*Nprev)

           ip1 = Nprev*2
           ip2 = ip1*N
           ip3 = ip2*Nrem
           i2rev= 1

           !******* BIT REVERSAL ********!
           do i2 = 1, ip2, ip1
                if(i2 .lt. i2rev) then
                  do i1 = i2, i2+ip1-2, 2
                   do i3 = i1,ip3,ip2
                       i3rev = i2rev + i3 - i2
                       call swap(data(i3), data(i3rev))
                       call swap(data(i3+1), data(i3rev+1))
                   enddo
                  enddo
                 endif

                 ibit=ip2/2
                 do while( (ibit .ge. ip1) .and. (i2rev .gt. ibit))
                       i2rev=i2rev-ibit
                       ibit=ibit/2
                 enddo
                 i2rev = i2rev + ibit
           enddo

           !****** DANIELSON LANCZOS ALGO FOR FFT****!
           ifp1=ip1
           do while(ifp1 .lt. ip2)
                ifp2=2*ifp1
                theta= ISIGN*2d0*pi/(ifp2/ip1)
                Wp = dcmplx( -2d0*dsin(0.5d0*theta)**2, dsin(theta))
                W  = dcmplx( 1d0 ,0 )
                do i3=1, ifp1, ip1
                 do i1=i3,i3+ip1-2 ,2
                  do i2=i1,ip3,ifp2
                      k1 = i2
                      k2 = k1+ifp1
                      temp = W*dcmplx(data(k2), data(k2+1))
                      data(k2)   = data(k1)   -dreal(temp)
                      data(k2+1) = data(k1+1) -dimag(temp)
                      data(k1)   = data(k1)   +dreal(temp)
                      data(k1+1) = data(k1+1) +dimag(temp)
                   enddo
                  enddo
                  W = W*Wp + W
                enddo
               ifp1=ifp2
           enddo
           Nprev = N*Nprev
         enddo  
      endsubroutine FOURN

      subroutine swap(a,b)
        implicit none
        double precision :: a,b
        double precision :: temp
             temp = a
             a = b
             b = temp
      endsubroutine swap

      subroutine nextpower(N1,N2)
       implicit none
       real :: pow
       integer:: N1, N2
      
       pow = log(real(N1))/log(2.0)

       if( mod(pow,1.0) > 0.0) then
         pow = pow+1.0
         N2 = 2**int(pow)
       else
         N2 = N1
       endif
      endsubroutine nextpower


