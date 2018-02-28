include 'fftlib.for'
program errori
        implicit none
        integer :: N, i
        double precision,allocatable ,dimension(:):: un ,ue, error, ufft
        double precision :: time, minerr, maxerr, norm, dummy
        character(50) :: fname
        
        print*, 'N?'
        read*, N

        allocate(un(0:N),ue(0:N),error(0:N), ufft(1:N/2))

        call system('ls sol* > fileread')
        call system('mkdir error')
        open(1,file='error/fileread')
        open(2,file='error/error_history.dat')
        !write(2,*) 'Variables = time, min_error, max_error, L2norm_error'
        write(2,*) 'Variables = time, disp_error, phase_error, max_error'
        do
          read(1,*) fname
          open(3, file=trim(fname), status='old')
          print*, fname
          read(3,*)
          read(3,*)
          read(3,*)
          do i=0,N
              read(3,*) dummy, un(i), ue(i)
          enddo
          close(3)
          error = un - ue

         call FFT1d(error(1:N), N, ufft)
         call writefft1d(ufft,N,100d0/N,'error/fft_'//trim(fname))

          print*,  fname(5:13) 
        ! minerr = minval(error)
        ! maxerr = maxval(error)
        ! norm = dsqrt(dot_product(error,error))/N

          minerr = maxval(n/4:n/2)
          maxerr = maxval(1:n/4)
          norm   = maxval(error)
          write(2,*)  fname(5:13),minerr,maxerr,norm
        enddo                              
endprogram errori
