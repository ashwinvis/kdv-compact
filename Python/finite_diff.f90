SUBROUTINE periodic_tdma_solver(a,b,c,d,X,N)
    IMPLICIT NONE
    integer, intent(in)::N
    double precision, intent(in), dimension(1:N)::a,b,c,d
    double precision, intent(out),dimension(1:N)::X
    double precision, allocatable,dimension(:) :: p , q , r , w , u , x1     
    INTEGER :: i
    double precision::sum1=0.0d0

    allocate(p(2:N),q(1:N),r(1:N-2),w(1:N-2),u(1:N-1),x1(1:N))

    q(1) = b(1) 
    w(1) = a(1) / b(1) 
    u(1) = c(1) / b(1) 
    r(1) = c(N) 
    sum1 = r(1) * w(1) 

    DO i = 2 , N-2
      p(i) = a(i)            
      q(i) = b(i) - p(i) * u(i-1)    
      u(i) = c(i) / q(i)         
      w(i) = - p(i) * w(i-1) / q(i)  
      r(i) = - r(i-1) * u(i-1)       
      sum1  = sum1 + r(i) * w(i)       
    ENDDO

    i = N-1 
    p(i) = a(i) 
    q(i) = b(i) - p(i) * u(i-1) 
    u(i) = (c(i) - p(i) * w(i-1)) / q(i) 
    p(N) = a(N) - r(N-2) * u(N-2) 
    q(N) = b(N) - sum1 - p(N) * u(N-1) 
    x1(1) = d(1) / q(1)  
!  
    DO i = 2 , N-1
      x1(i) = (d(i) - p(i) * x1(i-1)) / q(i)    
    ENDDO
    sum1 = 0.0d0 
    DO i = 1 , N-2 
      sum1 = sum1 + r(i) * x1(i)    
    ENDDO 
    x1(N) = (d(N) - sum1 - p(N) * x1(N-1)) / q(N)
    X(N) = x1(N)
    
    X(N-1) = x1(N-1) - u(N-1) * X(N)

    DO i = N-2 , 1 , -1
      X(i) = x1(i) - u(i) * X(i+1) - w(i) * X(N)    
    ENDDO  

    deallocate(x1,p,q,u,w,r)

END SUBROUTINE periodic_tdma_solver


subroutine Periodic_OUCS3(temp,tempzi,N,dzi)
    implicit none
    
    integer::j,N
    double precision::alpha1=0.0d0
    double precision::b2=1.0d0
    double precision,dimension(-2:2)::Const
    double precision::b1,b3, dzi
    double precision,dimension(1:N)::an,bn,cn,Rn
    double precision,intent(in),dimension(0:N)::temp
    double precision,intent(out),dimension(0:N)::tempzi

    Const(-2)=-0.183205192/4.0d0+alpha1/300.0d0
    Const( 2)= 0.183205192/4.0d0+alpha1/300.0d0
    Const(-1)=-1.57557379/2.0d0+alpha1/30.0d0
    Const( 1)= 1.57557379/2.0d0+alpha1/30.0d0
    Const( 0)=-11.0d0*alpha1/150.0d0
    
    b2= 1.0d0
    b1= 0.3793894912d0-alpha1/60.0d0
    b3= 0.3793894912d0+alpha1/60.0d0
    
    an=b1
    bn=b2
    cn=b3

    do j=2,N-2
      Rn(j)=(Const(-2)*temp(j-2)+Const(-1)*temp(j-1)+Const(0)*temp(j)&
             &+Const(1)*temp(j+1)+Const(2)*temp(j+2))/dzi        
    enddo
    Rn(1)  =(Const(-2)*temp(N-1)+Const(-1)*temp(N)+Const(0)*temp(1)&
             &+Const(1)*temp(2)+Const(2)*temp(3))/dzi
    Rn(N)  =(Const(-2)*temp(N-2)+Const(-1)*temp(N-1)+Const(0)*temp(N)&
             &+Const(1)*temp(1)+Const(2)*temp(2))/dzi
    Rn(N-1)=(Const(-2)*temp(N-3)+Const(-1)*temp(N-2)+Const(0)*temp(N-1)&
             &+Const(1)*temp(0)+Const(2)*temp(1))/dzi
    call periodic_tdma_solver(an,bn,cn,Rn,tempzi(1:N),N)
    tempzi(0)=tempzi(N)
end subroutine Periodic_OUCS3


subroutine Periodic_O3CS(temp,tempzi3,N,dzi3)     ! for 3rd derivative
    implicit none
    
    integer::j,N
    double precision::b2=1.0d0
    double precision,dimension(-2:2)::Const
    double precision::b1, b3, dzi3
    double precision,dimension(1:N)::an,bn,cn,Rn
    double precision,intent(in),dimension(0:N)::temp
    double precision,intent(out),dimension(0:N)::tempzi3

    Const(-2)= -0.957036d0      !+alpha1/300.0d0
    Const( 2)=  0.957036d0      !+alpha1/300.0d0
    Const(-1)=  1.914072d0      !+alpha1/30.0d0
    Const( 1)= -1.914072d0      !+alpha1/30.0d0
    Const( 0)=  0d0         !-11.0d0*alpha1/150.0d0
    
    b2= 1.0d0
    b1= 0.457036d0     !-alpha1/60.0d0
    b3= 0.457036d0     !+alpha1/60.0d0
    
    an=b1
    bn=b2
    cn=b3

    do j=2,N-2
      Rn(j)=(Const(-2)*temp(j-2)+Const(-1)*temp(j-1)+Const(0)*temp(j)&
             &+Const(1)*temp(j+1)+Const(2)*temp(j+2))/dzi3        
    enddo
    Rn(1)  =(Const(-2)*temp(N-1)+Const(-1)*temp(N)+Const(0)*temp(1)&
             &+Const(1)*temp(2)+Const(2)*temp(3))/dzi3
    Rn(N)  =(Const(-2)*temp(N-2)+Const(-1)*temp(N-1)+Const(0)*temp(N)&
             &+Const(1)*temp(1)+Const(2)*temp(2))/dzi3
    Rn(N-1)=(Const(-2)*temp(N-3)+Const(-1)*temp(N-2)+Const(0)*temp(N-1)&
             &+Const(1)*temp(0)+Const(2)*temp(1))/dzi3
    call periodic_tdma_solver(an,bn,cn,Rn,tempzi3(1:N),N)
    tempzi3(0)=tempzi3(N)
end subroutine Periodic_O3CS


subroutine convection(spatial_conv, UUx, u, N, dx)
    IMPLICIT NONE
    integer::i
    integer, intent(in)::spatial_conv, N
    double precision, intent(in)::dx
    double precision, intent(out), dimension(0:N)::UUx
    double precision, intent(in), dimension(0:N)::u
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
end subroutine convection


subroutine dispersion(spatial_disp, Uxxx, u, N, dx3)
    IMPLICIT NONE
    integer::i
    integer, intent(in)::spatial_disp, N
    double precision, intent(in)::dx3
    double precision, intent(out),dimension(0:N)::Uxxx
    double precision, intent(in), dimension(0:N)::u
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
endsubroutine dispersion
