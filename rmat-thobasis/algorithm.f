ccccccc
!algorithm and special functions
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
ccccccc
c *** Calculate du(r)/dr using five points derivative formula
c     f(ndim)=function to make derivative
c     h      =step
c     j      =point for derivative
      complex*16 function deriv1(f,h,ndim,j)
        implicit none
        integer ndim,j
        complex*16::f(ndim),h

        if ((j.eq.1).or.(j.eq.2)) then
           deriv1=(-f(j+2)+4d0*f(j+1)-3d0*f(j))/2d0/h
        else if (j.eq.ndim-1) then
           deriv1=(3d0*f(j)-4d0*f(j-1)+f(j-2))/2d0/h
        else if (j.eq.ndim) then
           deriv1=0 !!!CHECK
        else ! five points formula
           deriv1=(f(j-2)-8*f(j-1)+8*f(j+1)-f(j+2))/h/12.
         end if
      end function deriv1
ccccccc        
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! five points derivative formula for second derivative
! y: function value array
! d2y: second derivative array
! n:size of the array
! uniform grid 
!!complex type 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine second_derivative(y,d2y,n,dx)
       implicit none
       integer,intent(in)::n
       complex*16,intent(in)::dx
       complex*16,dimension(1:n),intent(in)::y
       complex*16,dimension(1:n),intent(out)::d2y
       integer::i 
       d2y(1)=(35.d0/12.d0*y(1)-26.d0/3.d0*y(2)+19.d0/2.d0*y(3)
     &   -14.d0/3.d0*y(4)+11.d0/12.d0*y(5))/(dx**2)
       d2y(2)=(11.d0/12.d0*y(1)-5.d0/3.d0*y(2)+1.d0/2.d0*y(3)
     &  +1.d0/3.d0*y(4)-1.d0/12.d0*y(5))/(dx**2)
       d2y(n-1)=(-1.d0/12.d0*y(N-4)+1.d0/3.d0*y(N-3)+1.d0/2.d0*y(N-2)
     &  -5.d0/3.d0*y(N-1)+11.d0/12.d0*y(N))/(dx**2)
       d2y(n)=(11.d0/12.d0*y(N-4)-14.d0/3.d0*y(N-3)+19.d0/2.d0*y(N-2)
     &  -26.d0/3.d0*y(N-1)+35.d0/12.d0*y(N))/(dx**2)
       do i=3,n-2
       d2y(i)=(-y(i-2)+16.d0*y(i-1)-30.d0*y(i)+ 
     & 16.d0*y(i+1)-y(i+2))/(12.d0*dx**2)
          !  write(*,*) d2y(i)
       end do
       end subroutine second_derivative
ccccccc
!complex interpolation function for uniform grids
      FUNCTION FFC(PP,F,N)
      COMPLEX*16 FFC,F(N)
      REAL*8 PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC=F(N)
      RETURN
      END function
ccccccc
        end module algorithm