include 'OUCS.f90'
include 'fftlib.for'
module Wave
       integer, parameter :: N = 2048
       double precision,parameter ::  x0=50d0, xmax = 100d0
       double precision :: t, told,x, dx = xmax / (N), dx3,  k, c,amp, Cconv ,Cdisp, L2norm
       double precision,dimension(0:N) :: u, usq, uold, unew, UUx, Uxxx, Lu , uexact, error
endmodule Wave

module RK
      double precision :: dt=1d-5, tmax
      integer :: stage, Niter, Ncycle=5 ,writeafter
      integer :: spatial_conv = 5 ,spatial_disp = 2 ! 1-CD2, 2-UD3,  3-CD4, 4-CD4 non conserv, 5-OUCS3; 1-CD2, 2-CD4, 3-O3CS
      double precision,dimension(1:4) :: fac_dt, fac_rk
endmodule RK

subroutine Initialize
        use Wave
        use RK
        implicit none
        integer :: i
       ! k   = 1d0        
        amp = 2d0*k**2*0.999d0
        c   = 4d0*k**2

        tmax    = xmax/c*Ncycle
        Niter   = tmax/dt
        writeafter = 1/dt/10

        Cconv=6d0
        Cdisp=1d0
        x=0d0
        dx3 = dx**3

!       u=0        ! Delta function
!       u(N/2)=-1d0        

        do i = 0,N ! Soliton wave solution
             u(i) = amp* (1d0- dtanh(k*(x-x0))**2d0)
             x=x+dx
        enddo                

        t = 0d0
        fac_dt(1) = 0.5d0 
        fac_dt(2) = 0.5d0
        fac_dt(3) = 1d0
        fac_dt(4) = 1d0
        
        fac_rk(1) = 1d0/6d0 
        fac_rk(2) = 1d0/3d0
        fac_rk(3) = 1d0/3d0
        fac_rk(4) = 1d0/6d0
endsubroutine Initialize

subroutine Backup
        use Wave
        implicit none
        uold = u
        unew = u
endsubroutine Backup

subroutine Restore
        use Wave
        implicit none
        u = unew
endsubroutine Restore

subroutine SpatialDerivative
        use Wave
        use RK
        implicit none
        integer :: i

        !*** Convection term****!
        select case(spatial_conv)
           case(1)
           !CD2
               do i = 1,N-1
                 UUx(i) = (u(i+1)**2 - u(i-1)**2)/2d0/dx/2d0
               enddo
               UUx(N) = (u(1)**2 - u(N-1)**2)/2d0/dx/2d0
               UUx(0) = UUx(N)
                
           case(2)
           !UD3
              UUx(1) = (3d0*u(1)**2 - 4d0*u(0)**2 + u(N-1)**2)/2d0/dx/2d0
              do i = 2,N
                UUx(i) = (3d0*u(i)**2 - 4d0*u(i-1)**2 + u(i-2)**2)/2d0/dx/2d0
              enddo
              UUx(0) = UUx(N)
         
           case(3)
           !CD4
               do i = 2,N-2
                 UUx(i) = -(u(i+2)**2 - 8d0*u(i+1)**2 + 8d0*u(i-1)**2 - u(i-2)**2)/12d0/dx/2d0
               enddo
               UUx(1)   = -(u(3)**2 - 8d0*u(2)**2 + 8d0*u(N)**2   - u(N-1)**2)/12d0/dx/2d0
               UUx(N)   = -(u(2)**2 - 8d0*u(1)**2 + 8d0*u(N-1)**2 - u(N-2)**2)/12d0/dx/2d0
               UUx(N-1) = -(u(1)**2 - 8d0*u(0)**2 + 8d0*u(N-2)**2 - u(N-3)**2)/12d0/dx/2d0
               UUx(0) = UUx(N)
           case(4)
           !CD4 - Non_conservation form
               do i = 2,N-2
                 UUx(i) = -u(i)*(u(i+2) - 8d0*u(i+1) + 8d0*u(i-1) - u(i-2))/12d0/dx
               enddo
               UUx(1)   = -u(1)*(u(3) - 8d0*u(2)+ 8d0*u(N)   - u(N-1))/12d0/dx
               UUx(N)   = -u(N)*(u(2) - 8d0*u(1)+ 8d0*u(N-1) - u(N-2))/12d0/dx
               UUx(N-1) = -u(N-1)*(u(1) - 8d0*u(0)+ 8d0*u(N-2) - u(N-3))/12d0/dx
               UUx(0) = UUx(N)
           case(5)
              ! call Dealias1d(u, usq,N)
              ! call ShiftAlign(u,usq,N)
               call Periodic_OUCS3(u*u,UUx,N,dx)
               UUx= UUx/2d0
        end select

        ! Dispersion term
        select case(spatial_disp)
        !CD2
            case(1)
               do i = 2,N-2
                 Uxxx(i) = (u(i+2) - 2d0*u(i+1) + 2d0*u(i-1) - u(i-2))/2d0/dx3
               enddo
               Uxxx(1) = (u(3) - 2d0*u(2) + 2d0*u(N) - u(N-1))/2d0/dx3
               Uxxx(N) = (u(2) - 2d0*u(1) + 2d0*u(N-1) - u(N-2))/2d0/dx3
               Uxxx(N-1) = (u(1) - 2d0*u(N) + 2d0*u(N-2) - u(N-3))/2d0/dx3
               Uxxx(0) = Uxxx(N)
        !CD4
            case(2)
               do i = 3,N-3
                 Uxxx(i) = (-u(i+3) +8d0*u(i+2) - 13d0*u(i+1) + 13d0*u(i-1) - 8d0*u(i-2) + u(i-3))/8d0/dx3
               enddo
               Uxxx(2)   = (-u(5) +8d0*u(4) - 13d0*u(3) + 13d0*u(1) - 8d0*u(N)   + u(N-1))/8d0/dx3
               Uxxx(1)   = (-u(4) +8d0*u(3) - 13d0*u(2) + 13d0*u(N) - 8d0*u(N-1) + u(N-2))/8d0/dx3
               Uxxx(N)   = (-u(3) +8d0*u(2) - 13d0*u(1) + 13d0*u(N-1) - 8d0*u(N-2) + u(N-3))/8d0/dx3
               Uxxx(N-1) = (-u(2) +8d0*u(1) - 13d0*u(0) + 13d0*u(N-2) - 8d0*u(N-3) + u(N-4))/8d0/dx3
               Uxxx(N-2) = (-u(1) +8d0*u(0) - 13d0*u(N-1) + 13d0*u(N-3) - 8d0*u(N-4) + u(N-5))/8d0/dx3
               Uxxx(0)   = Uxxx(N)
        !O3CS
             case(3)
               call Periodic_O3CS(u,Uxxx,N,dx3)
       end select
endsubroutine SpatialDerivative

subroutine RKstage
        use Wave
        use RK
        implicit none
        
        call SpatialDerivative
        Lu = -dt *(Cconv*UUx + Cdisp*Uxxx)

        unew = unew + fac_rk(stage)*Lu
        
        u = uold + fac_dt(stage)*Lu
endsubroutine RKstage

subroutine CalcError
        use Wave
        use RK, only: Ncycle
        implicit none
        integer :: i, ii

        x  = 0d0
        ii = 0
        uexact = 0d0

        do i=0,N*(Ncycle+1)
         uexact(ii) = uexact(ii)+amp*(1d0 - dtanh(k*(x-x0 - c*t))**2 )
         if(i/=0 .and. mod(i,N)==0) then
            ii = 0
            uexact(ii) = uexact(N)
         endif
         x  = x+dx
         ii = ii+1
        enddo
        
        error = dabs(u-uexact)
        L2norm = dsqrt(dot_product(error,error)/N)
endsubroutine CalcError

subroutine March
       use Wave
       use RK 
       implicit none 
       integer :: i

       call Initialize
        
       call CalcError
       call SolnWrite
       open(1,file='L2norm.dat',status='replace')
       write(1,*) 'Variables= t, L2norm'

       do i=0,Niter!while(t < tmax)
               call Backup
               do stage=1,4
                 call RKstage
               enddo
               call Restore

               t = t + dt

               if(mod(i,writeafter)==0) then
                call CalcError
                call SolnWrite
                write(*,*) 't=',t, L2norm !'umax=', maxval(u), 'logerrormax=', log10(maxval(error)), 'L2norm=',L2norm
                write(1,*) t, L2norm
               endif
       enddo
       close(3)

endsubroutine March

subroutine SolnWrite
        use Wave
        use RK
       implicit none
       integer :: i
       character(40) :: fname

       if(int(t)<10) then      
        write(fname,11) 'sol_000',t,'.dat'
        11 format(A7,F6.4,A4)
       elseif(int(t)<100) then 
        write(fname,12) 'sol_00',t,'.dat'
        12 format(A6,F7.4,A4)
       elseif(int(t)<1000) then 
        write(fname,13) 'sol_0',t,'.dat'
        13 format(A5,F8.4,A4)
       else
        write(fname,14) 'sol_',t,'.dat'
        14 format(A4,F9.4,A4)
       endif

       open(101, file=trim(fname),status='replace')
       
       write(101,*) 'Variables = x, u'!, u_exact, error'
       write(101,*) 'Zone i=',N+1
       write(101,*) 'Auxdata time="',real(t),'"'
       
       111 format(4ES15.3)
       x=0d0
       do i=0,N
        write(101,111) x,u(i)!, uexact(i), error(i)!, (i), error(i)
        x=x+dx
       enddo

       close(101)
endsubroutine SolnWrite

program KdV
        use Wave
        use RK
        implicit none
        print*, "k, (2/3) CD4/O3CS"
        read*, k, spatial_disp
        call Initialize
        call March
        
endprogram KdV
