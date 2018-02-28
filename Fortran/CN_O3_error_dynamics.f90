MODULE Var

INTEGER:: i,j,k,N=51,p,q,r,m1
DOUBLE PRECISION:: D,E,F,pm1,pp1,qm2,qm1,qp1,qp2,q0
DOUBLE PRECISION::alp,beta2,betaN,Nc,Ne,kh,pi,eta,eps2,eps4,k2,k4,Gr,Gi,a1,b1
DOUBLE PRECISION,DIMENSION(:,:):: A(51,51),B(51,51),C(51,51),G(314,50,50)
DOUBLE PRECISION,DIMENSION(:,:):: A2(51,51),B2(51,51),C2(51,51)
DOUBLE PRECISION,DIMENSION(:,:):: cnbyc1(314,50,50),cnbyc(314,50,50),dcnbyc1(314,50,50)
DOUBLE PRECISION,DIMENSION(:,:):: enbye(314,50,50),dcnbyc(314,50,50),denbye(314,50,50)
DOUBLE PRECISION,DIMENSION(:,:):: AI(51,51),X(51),AA(51,102),Dumm(51,51)
DOUBLE PRECISION,DIMENSION(:,:):: betaj(314,50,50),dbetaj(314,50,50)
DOUBLE PRECISION,DIMENSION(:,:):: term(314,50,50),dterm(314,50,50),vgnbyvg(314,50,50)
DOUBLE PRECISION,DIMENSION(:,:):: phaser(314,50,50),dphaser(314,50,50)
DOUBLE COMPLEX,DIMENSION(:,:):: Aj(314,50,50),Gj(314,50,50)
DOUBLE COMPLEX:: iota,sumC,sumC2
DOUBLE COMPLEX,DIMENSION(:):: C1(51)
character*30:: filename,filename1
END MODULE Var
!-----------------------------------------------------------------------------------------------------------------------------
Subroutine Define_AB
use Var
eta =  0.0d0
D = 0.3793894912d0
E = 1.57557379d0
F = 0.183205192d0
beta2 = -0.025d0
betaN = 0.09d0     ! for N-2
pm1 = D - eta/60.0d0
pp1 = D + eta/60.0d0
qm2 = -(F/4.0d0) + (eta/300.0d0)
qm1 = -(E/2.0d0) + (eta/30.0d0)
q0 = -11.0d0*eta/150.0d0
qp1 = (E/2.0d0) + (eta/30.0d0)
qp2 = (F/4.0d0) + (eta/300.0d0)
DO i = 1,N
DO j = 1,N
   A(i,j) = 0.0d0
   B(i,j) = 0.0d0
END DO
END DO
A(1,1) = 1.0d0
A(2,2) = 1.0d0
A(N-1,N-1) = 1.0d0
A(N,N) = 1.0d0
DO i = 3,N-2
 A(i,i) = 1.0d0
 A(i,i-1) = pm1    ! D - eta/60.0d0
 A(i,i+1) = pp1    ! D + eta/60.d0
END DO

B(1,1) = -1.5d0
B(1,2) = 2.0d0
B(1,3) =  -0.5d0
B(N,N) = 1.5d0
B(N,N-1) = -2.0d0
B(N,N-2) = 0.5d0
B(2,1) = (2.0d0*beta2/3.0d0) - (1.0d0/3.0d0)
B(2,2) = -((8.0d0*beta2/3.0d0) + (0.50d0))
B(2,3) = 4.0d0*beta2 + 1.0d0
B(2,4) = -((8.0d0*beta2/3.0d0) + (1.0d0/6.0d0))
B(2,5) = 2.0d0*beta2/3.0d0
B(N-1,N) = -((2.0d0*betaN/3.0d0) - (1.0d0/3.0d0))
B(N-1,N-1) = (8.0d0*betaN/3.0d0) + (0.50d0)
B(N-1,N-2) = -(4.0d0*betaN + 1.0d0)
B(N-1,N-3) = (8.0d0*betaN/3.0d0) + (1.0d0/6.0d0)
B(N-1,N-4) = -(2.0d0*betaN/3.0d0)
DO i = 3,N-2
B(i,i-2) = qm2   ! -(F/4.0d0) + (eta/300)
B(i,i-1) = qm1   ! -(E/2.0d0) + (eta/30)
B(i,i) = q0    ! -11.0d0*eta/150.0d0
B(i,i+1) = qp1   ! (E/2.0d0) + (eta/30)
B(i,i+2) =  qp2   ! (F/4.0d0) + (eta/300)
END DO

CALL Inverse_A
iota = (0.0d0,1.0d0)
 C  = matmul(AI,B)
!--------------------------------------------------
alp = 0.457    ! 0.47702202202203
a1 = -4.0d0*alp - 2.0d0
b1 = 16.0d0*alp + 8.0d0
DO i = 1,N
DO j = 1,N
   A(i,j) = 0.0d0
   B(i,j) = 0.0d0
END DO
END DO
A(1,1) = 1.0d0
A(2,2) = 1.0d0
A(N-1,N-1) = 1.0d0
A(N,N) = 1.0d0
DO i = 3,N-2
 A(i,i) = 1.0d0
 A(i,i-1) = alp 
 A(i,i+1) = alp 
END DO

B(1,1) = -1.5d0
B(1,2) = 2.0d0
B(1,3) =  -0.5d0
B(N,N) = 1.5d0
B(N,N-1) = -2.0d0
B(N,N-2) = 0.5d0
B(2,1) = (2.0d0*beta2/3.0d0) - (1.0d0/3.0d0)
B(2,2) = -((8.0d0*beta2/3.0d0) + (0.50d0))
B(2,3) = 4.0d0*beta2 + 1.0d0
B(2,4) = -((8.0d0*beta2/3.0d0) + (1.0d0/6.0d0))
B(2,5) = 2.0d0*beta2/3.0d0
B(N-1,N) = -((2.0d0*betaN/3.0d0) - (1.0d0/3.0d0))
B(N-1,N-1) = (8.0d0*betaN/3.0d0) + (0.50d0)
B(N-1,N-2) = -(4.0d0*betaN + 1.0d0)
B(N-1,N-3) = (8.0d0*betaN/3.0d0) + (1.0d0/6.0d0)
B(N-1,N-4) = -(2.0d0*betaN/3.0d0)
DO i = 3,N-2
B(i,i-2) = -b1/16.0d0
B(i,i-1) = -a1/2.0d0
B(i,i) = 0.0d0
B(i,i+1) = a1/2.0d0
B(i,i+2) = b1/16.0d0
END DO

CALL Inverse_A
 C2 = matmul(AI,B)

End Subroutine Define_AB

subroutine Inverse_A
use Var
       IMPLICIT NONE
       INTEGER:: L3,LI,LK,KP1,JJ,LJ,M2,M
       REAL:: P1,T,AB,BIG                                
         M = N + N                               
         M2 = N + 1
!      GENERATING THE AUGMENTED MATRIX AA
         DO I = 1, N            !ENTERING THE COEFFICIENT MATRIX A INTO  
         DO  J = 1, N           !THE AUGMENTED MATRIX AA       
           AA(I, J) = A(I, J)
         END DO
         END DO          
        DO  I = 1, N
        DO  J = M2, M
           AA(I, J) = 0.0
        END DO
        END DO       
        DO I = 1, N       !ENTERING THE IDENTITY MATRIX I INTO
            J = I + N         !THE AUGMENTED MATRIX AA
       AA(I, J) = 1.0 
        END DO             
!      STARTING ROW TRANSFORMATION
       DO LJ = 1, N
          K = LJ
          IF(K .LT. N)THEN
            JJ = K
            BIG = ABS(AA(K, K))
            KP1 = K + 1
            DO I = KP1, N
              AB = ABS(AA(I, K))
              IF((BIG - AB) .LT. 0.0)THEN
                 BIG = AB               !PERFORMING PIVOTING AND
                 JJ = I                 !ROW TRANSFORMATION OPERATIONS
              ENDIF
          END DO
            IF ((JJ - K) .GT. 0.0)THEN
              DO J = K, M
                 T = AA(JJ, J)
                 AA(JJ, J) = AA(K, J)
                 AA(K, J) = T
            END DO
            ENDIF
         ENDIF
          P1 = AA(LJ, LJ)
          DO I = LJ, M
             AA(LJ, I) = AA(LJ, I)/P1
          END DO
          DO LK = 1, N
            T = AA(LK, LJ)
            DO LI = LJ, M
              IF((LK - LJ) .NE. 0)THEN                 
                AA(LK, LI) = AA(LK, LI) - AA(LJ, LI)*T 
              ENDIF
      END DO
      END DO
   END DO
         DO I = 1, N
         DO J = M2, M
            L3 = J - N
            AI(I, L3) = AA(I, J) !GENERATING THE INVERTED MATRIX
   END DO
   END DO         
End subroutine Inverse_A

Subroutine Properties
use Var

!DO r = 1,N,1
 r = 26

DO k = 1,50
Ne = 0.01d0*k
DO j = 1,50
Nc = 0.01d0*j
DO i = 1,314
kh = 0.01d0*i

sumC = (0.0d0,0.0d0)
sumC2 = (0.0d0,0.0d0)

DO m1 = 1,N
sumC = sumC + C(r,m1)*exp(iota*kh*(m1-r))    ! node = r
END DO  

DO m1 = 1,N
sumC2 = sumC2 + C2(r,m1)*exp(iota*kh*(m1-r))    ! node = r
END DO

Aj(i,j,k) =  Nc*sumC - Ne*sumC2   

Gj(i,j,k) = 1.0d0 - Aj(i,j,k) + ((Aj(i,j,k))**2)/2.0d0 - ((Aj(i,j,k))**3)/6.0d0 + ((Aj(i,j,k))**4)/24.0d0

Gr = real(Gj(i,j,k))
Gi = aimag(Gj(i,j,k))

G(i,j,k) = sqrt(Gr**2 + Gi**2)  

pi = 4.0d0*atan(1.0d0)

!betaj(i,j,k) = (atan(-Gi/Gr))

IF((Gr > 0.0d0) .and. (-Gi>0.0d0))THEN
 betaj(i,j,k) = (atan(-Gi/Gr))
ELSEIF((Gr < 0.0d0) .and. (-Gi > 0.0d0))THEN
 betaj(i,j,k) = pi + (atan(-Gi/Gr))        ! was pi +     
ELSEIF((Gr < 0.0d0) .and. (-Gi < 0.0d0))THEN
 betaj(i,j,k) = pi + atan(-Gi/Gr)
ELSEIF((Gr > 0.0d0) .and. (-Gi < 0.0d0))THEN
 betaj(i,j,k) = 2.0d0*pi +  atan(-Gi/Gr)   ! was 2pi +                 
ENDIF

   cnbyc(i,j,k) = betaj(i,j,k)/(kh*Nc + Ne*kh*kh*kh) 

   phaser(i,j,k) = (betaj(i,j,k)/kh) - Nc - Ne*kh*kh  ! this is dt*(cn-c)

   term(i,j,k) = Nc + Ne*kh*kh
END DO
END DO
END DO

DO k = 1,50
Ne = 0.01d0*k
DO j = 1,50
Nc = 0.01d0*j
DO i = 1,314
kh = 0.01d0*i
      IF(i.eq.1) THEN
      dphaser(i,j,k)= (phaser(i+1,j,k)-phaser(i,j,k))/0.01
      ELSEIF(i.eq.314) THEN
      dphaser(i,j,k)= (phaser(i,j,k)-phaser(i-1,j,k))/0.01
      ELSE
      dphaser(i,j,k)= (phaser(i+1,j,k)-phaser(i-1,j,k))/(2.0d0*0.01)
      ENDIF
     dcnbyc(i,j,k) = dphaser(i,j,k)    ! this is dt*[d(cn-c)/dk]
END DO
END DO
END DO

DO k = 1,50
Ne = 0.01d0*k
DO j = 1,50
Nc = 0.01d0*j
DO i = 1,314
kh = 0.01d0*i
      IF(i.eq.1) THEN
      dbetaj(i,j,k)= (betaj(i+1,j,k)-betaj(i,j,k))/0.01
      dterm(i,j,k)= (term(i+1,j,k)-term(i,j,k))/0.01
      ELSEIF(i.eq.314) THEN
      dbetaj(i,j,k)= (betaj(i,j,k)-betaj(i-1,j,k))/0.01
      dterm(i,j,k)= (term(i,j,k)-term(i-1,j,k))/0.01
      ELSE
      dbetaj(i,j,k)= (betaj(i+1,j,k)-betaj(i-1,j,k))/(2.0d0*0.01)
      dterm(i,j,k)= (term(i+1,j,k)-term(i-1,j,k))/(2.0d0*0.01)
      ENDIF
     vgnbyvg(i,j,k) = dbetaj(i,j,k)/ ((Nc+kh*kh*Ne) + kh*(dterm(i,j,k)))  
END DO
END DO
END DO

write(filename,'(a,I2,a)') trim("O3_rk4_j_"),r,'.dat'
open(5,file = filename)
write(5,*)'Variables = Ne, Nc, kh, G, cn/c, (dt/h)*dcn/dk, (dt/h)*phase_err, vgnby/vg'
Write(5,*)'zone I=50 J=50 K=314 F=POINT'

DO i = 1,314
kh = 0.01d0*i
DO j = 1,50
Nc = 0.01d0*j
DO k = 1,50
Ne = 0.01d0*k
write(5,55) Ne, Nc, kh, G(i,j,k), cnbyc(i,j,k), dcnbyc(i,j,k), phaser(i,j,k), vgnbyvg(i,j,k)
55 FORMAT(3F12.8,F24.4,4F16.8)
END DO
END DO
END DO

print*,"node complete = ",r
!END DO
End Subroutine Properties
!-----------------------------------------------------------------------------------------------------------------------------
PROGRAM CNO3
use Var

CALL Define_AB
CALL Inverse_A
CALL Properties

END PROGRAM
