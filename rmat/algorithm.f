ccccccc
        module algorithm
            contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!gauss-legendre integral, N:integral mesh number; 
!x1,x2: integral interval;    
!x,w: mesh point and weight
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N
        REAL*8 x1,x2,X(N),W(N)
        REAL*8 z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j

        pi=acos(-1.0)
        tol=1.E-12

        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)

         DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))

 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende

         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
        END SUBROUTINE gauleg 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        end module algorithm