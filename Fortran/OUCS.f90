SUBROUTINE periodic_tdma_solver(a,b,c,d,X,N)
	IMPLICIT NONE
	integer, intent(in)::N
	double precision, intent(in), dimension(1:N)::a,b,c,d
	double precision, intent(out),dimension(1:N)::X
	double precision, allocatable,dimension(:) :: p , q , r , w , u , x1 	
	INTEGER :: i,j
	double precision::sum1=0.0d0,maxerror
	character(80)::str
	double precision, allocatable, dimension(:)::error

	allocate(p(2:N),q(1:N),r(1:N-2),w(1:N-2),u(1:N-1),x1(1:N),error(1:N))

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

	deallocate(x1,p,q,u,w,r,error)

END SUBROUTINE periodic_tdma_solver

subroutine Periodic_OUCS3(temp,tempzi,N,dzi)
	implicit none
	
	integer::i,j,k,N
	double precision::alpha1=0.0d0
	double precision::b2=1.0d0
	double precision,dimension(-2:2)::Const
	double precision::term1,term2,beta,b1,b3, dzi
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
	
	integer::i,j,k,N
	double precision::alpha1=0.0d0
	double precision::b2=1.0d0
	double precision,dimension(-2:2)::Const
	double precision::term1,term2,beta,b1,b3, dzi3
	double precision,dimension(1:N)::an,bn,cn,Rn
	double precision,intent(in),dimension(0:N)::temp
	double precision,intent(out),dimension(0:N)::tempzi3

	Const(-2)= -0.957036d0      !+alpha1/300.0d0
	Const( 2)=  0.957036d0      !+alpha1/300.0d0
	Const(-1)=  1.914072d0      !+alpha1/30.0d0
	Const( 1)= -1.914072d0      !+alpha1/30.0d0
	Const( 0)=  0d0             !-11.0d0*alpha1/150.0d0
	
	b2= 1.0d0
	b1= 0.457036d0         !-alpha1/60.0d0
	b3= 0.457036d0         !+alpha1/60.0d0
	
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
