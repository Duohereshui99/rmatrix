!!hobasis (complex)
        module basis
            use parameter
            implicit none
            integer::nbasis
            real(8)::alpha,gamma,m,b
!phi: basis, (1:Nr,1:Nbasis,1:Nch), each channel (L) in line with Nbasis basis functions
!d2phi: 2nd derivative of basis functions, size the same as that of phi
!phi1: the value of basis at uniformed mesh, size (1:ndiff)
!d2phi1:: the 2nd derivative of basis at uniformed mesh,size (1:ndiff)
!phia: the value of basis at the boundary r=a, phi(a), size: (1:Nbasis,1:Nch)
!phipa: the derivative of basis at the boundary r=a, phi'(a), size: (1:Nbasis,1:Nch)
            complex(16),allocatable::phi(:,:,:)
            complex(16),allocatable::d2phi(:,:,:)
            complex(16),allocatable::phi1(:)
            complex(16),allocatable::d2phi1(:)
            complex(16),allocatable::phia(:,:)
            complex(16),allocatable::phipa(:,:)
            contains
!!norm of ho basis
            real*8 function normho(nu,n,l)
            implicit none
            real*8::nu
            integer::n,l
            normho=sqrt(sqrt(2*nu**3/pi)*2**(n+2*l+3)*fact(n)
     &      *nu**l/doublefact(2*n+2*l+1))
      !   if(abs(normho)<1e-6) then
      !      write(*,*)'nu,n,l,normho',nu,n,l,normho
      !      write(*,*)'fact(n)',fact(n)
      !      write(*,*)'fact(n+l)=',fact(n+l)
      !      write(*,*)'fact(2n+2l+1)=',fact(2*n+2*l+1)
      !      stop
      !   endif
            end function
!!hobasis
            complex*16 function ho3d(n,l,nu,r)            !!3d hobasis
            implicit none                           !!nu=mu*omega/(2hbar)
            integer l,n
            real*8::norma,nu
            complex*16::r
            norma=normho(nu,n,l)
      !        if (norma<1e-6) then
      !           write(*,*)'ho3d: Norm=0!!!for  nu,n,l',nu,n,l
      !        endif
            ho3d=norma*r**l*exp(-nu*r**2)*
     &       generalized_laguerre(n,l+0.5d0,cmplx(2d0*nu*r**2,kind=16)) !!convert argument x into complex 16 type
            end function ho3d
!! fact            
            function fact(n)           
                  implicit none
                  integer n
                  real*8 fact, dgamma,x
                  x=dfloat(n+1)
                  fact=dgamma(x)
            end function fact
!!double fact
      real*8 function doublefact(n)       !double factorial
            implicit none
            integer::n,i
            real::s
            s=1.0
            if(mod(n,2)==0) then
                  do i=n,2,-2
                        s=s*i
                  end do
            else
                  do i=n,1,-2
                        s=s*i
                  end do
            end if
            doublefact=s
      end function
!!generalized_laguerre, non-recursive,alpha=l+1/2
      function generalized_laguerre(n, alpha, x) result(Ln_alpha_x)
            implicit none
            integer, intent(in) :: n
            real(8), intent(in) :: alpha
            complex(16), intent(in) :: x  
            complex(16) :: Ln_alpha_x          
            integer :: i
            complex(16) :: L0, L1, L2
!initialize the first two polynomials
            L0=1.0d0
            L1=1.0d0+alpha-x
ccccccc
            if (n==0) then
                Ln_alpha_x=L0
                return
            endif
ccccccc
            if (n==1) then
                Ln_alpha_x=L1
                return
            endif
ccccccc
            do i=1,n-1
                L2=((2.0d0*i+1.0d0+alpha-x)*L1-(i+alpha)*L0)/(i+1.0d0)
                L0=L1
                L1=L2
            end do
ccccccc
            Ln_alpha_x=L2
        end function generalized_laguerre
ccccccc
        complex*16 function D1Laguerre(n,alpha,x) !d/dx L(n,l+0.5,x),alpha=l+0.5
            implicit none
            integer::n
            real*8::alpha
            complex*16::x
            D1Laguerre=1d0/x*(n*generalized_laguerre(n,alpha,cmplx(x,kind=16))
     &       -(n+alpha)*generalized_laguerre(n-1,alpha,cmplx(x,kind=16)))
       end function D1Laguerre
ccccccc
        complex*16 function D2Laguerre(n,alpha,x) !d^2/dx^2 L(n,l+0.5,x),alpha=l+0.5
            implicit none
            integer::n
            real*8::alpha
            complex*16::x
            D2Laguerre=1d0/x**2*
     &       ((n**2d0-n)*generalized_laguerre(n,alpha,cmplx(x,kind=16))  
     &   -(2d0*n-2)*(n+alpha)*generalized_laguerre(n-1,alpha,cmplx(x,kind=16))
     &   +(n+alpha)*(n-1+alpha)*generalized_laguerre(n-2,alpha,cmplx(x,kind=16)))
      end function D2Laguerre
ccccccc       
      complex*16 function LSTFUN(gamma,m,r) !LST transformation s(r)
      implicit none
      real*8::gamma,m
      complex*16::r
      LSTFUN=(1/((1.d0/r)**m+
     & (1.d0/gamma/sqrt(r))**m))**(1.d0/m)
      end function
ccccccc
      complex*16 function D1LSTFUN(gamma,m,r)     !s'(r) analytic expression 
      implicit none
      real*8::gamma,m
      complex*16::x,y
      complex*16::r
      x=(1.d0/r)**(m+1.d0)+1.d0/(2.0*r)*(1.d0/(gamma*sqrt(r)))**m
      y=(1.d0/r)**m+(1.d0/(gamma*sqrt(r)))**m
      D1LSTFUN=LSTFUN(gamma,m,r)*x/y
      end function
ccccccc
!!alpha = nu in this module 
      complex*16 function THOFUNC(n,l,alpha,gamma,m,r)    !alpha: dimensionless parameter in hobasis
      implicit none
      integer::n,l
      real*8::alpha,gamma,m
      complex*16::r
      THOFUNC=sqrt(D1LSTFUN(gamma,m,r))
     & *ho3d(n,l,alpha,LSTFUN(gamma,m,r))
     & *LSTFUN(gamma,m,r)/r
      end function
ccccccc
      complex*16 function D2LSTFUN(gamma,m,r)     !s''(r) analytic expression
      implicit none
      real(8)::gamma,m
      complex(16)::r
      D2LSTFUN=-(1d0/((1d0/r)**m+(1d0/sqrt(r)/gamma)**m))**(2d0+(1d0/m))
     & *((2d0+m)*(1d0/r)**m+(1d0/sqrt(r)/gamma)**m)*(1d0/sqrt(r)/gamma)**m
     & /(4d0*r**2)
      end function
ccccccc
      complex*16 function D3LSTFUN(gamma,m,r)     !s'''(r) analytic expression
      implicit none
      real(8)::gamma,m
      complex(16)::r
      D3LSTFUN=1d0/8d0/r**3d0*(1d0/((1d0/r)**m+(1/sqrt(r)/gamma)**m))**(3d0+1d0/m)
     & *((m**2-4d0)*(1d0/r)**(2d0*m)+(8d0+3d0*m+m**2d0)*(1d0/r)**m
     &   *(1d0/sqrt(r)/gamma)**m+3d0*(1d0/sqrt(r)/gamma)**(2d0*m))
     & *(1d0/sqrt(r)/gamma)**m
      end function
ccccccc
      complex*16 function D1THOFUNC(n,l,alpha,gamma,m,r)    !d/dr R_{nl}^{THO}(r)
      implicit none
      integer::n,l
      real*8::alpha,gamma,m
      complex*16::r
      D1THOFUNC=(1d0/2d0/sqrt(D1LSTFUN(gamma,m,r)))*LSTFUN(gamma,m,r)**(l-1)*exp(-alpha*LSTFUN(gamma,m,r)**2)
     &  *normho(alpha,n,l)*
     & ( 
     &  8*alpha*LSTFUN(gamma,m,r)**2*D1LSTFUN(gamma,m,r)**2
     &  *D1Laguerre(n,l+0.5d0,2d0*alpha*LSTFUN(gamma,m,r)**2)
     &  +(2d0*D1LSTFUN(gamma,m,r)**2)*(l-2d0*alpha*LSTFUN(gamma,m,r)**2)   
     &  +LSTFUN(gamma,m,r)*D2LSTFUN(gamma,m,cmplx(r,kind=16))*
     *  generalized_laguerre(n,l+0.5d0,cmplx(2*alpha*LSTFUN(gamma,m,r)**2,kind=16))                         
     &  ) 
      end function
ccccccc
      complex*16 function D2THOFUNC(n,l,alpha,gamma,m,r)    !d^2/dr^2 R_{nl}^{THO}(r)
      implicit none
      integer::n,l
      real*8::alpha,gamma,m
      complex*16::r,sr,d1sr,d2sr,d3sr,LL,D1L,D2L
ccccccc
      sr=LSTFUN(gamma,m,r)
      d1sr=D1LSTFUN(gamma,m,r)
      d2sr=D2LSTFUN(gamma,m,cmplx(r,kind=16))
      d3sr=D3LSTFUN(gamma,m,cmplx(r,kind=16))
      LL=generalized_laguerre(n,l+0.5d0,cmplx(2d0*alpha*sr**2,kind=16))
      D1L=D1Laguerre(n,l+0.5d0,2d0*alpha*sr**2)
      D2L=D2Laguerre(n,l+0.5d0,2d0*alpha*sr**2)
ccccccc
      D2THOFUNC=1d0/4d0/d1sr**(3d0/2d0)*sr**(l-2d0)*exp(-alpha*sr**2)*normho(alpha,n,l)
     &  *(16d0*alpha*sr**2*d1sr**2*
     &  (4d0*alpha*sr**2*d1sr**2*D2L+
     &  (d1sr**2*(2d0*l-4d0*alpha*sr**2+1d0)+2d0*sr*d2sr)*D1L)
     &  +(4*d1sr**4  *( (2*alpha*sr**2) * (2*alpha*sr**2-2d0*l-1d0) + l*(l-1d0)
     &  +8d0*sr*d1sr**2*d2sr*(l-2*alpha*sr**2)-sr**2*d2sr*2+2d0*sr**2*d3sr*d1sr))   
     &  *LL) 
      end function

      end module