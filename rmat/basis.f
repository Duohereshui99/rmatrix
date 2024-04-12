ccccccc
!ho basis and tho basis and other possible basis
!
      module basis 
            use parameter
            contains

      function normho(nu,n,l)
        implicit none
        real*8:: nu,pi,normho
        integer:: n,l,m,p
        pi=acos(-1d0)

        normho=sqrt(sqrt(2*nu**3/pi)*2**(n+2*l+3)*fact(n)
     &  *nu**l/doublefact(2*n+2*l+1))
        if(normho<1e-6) then
           write(*,*)'nu,n,l,normho',nu,n,l,normho
           write(*,*)'fact(n)',fact(n)
           write(*,*)'fact(n+l)=',fact(n+l)
           write(*,*)'doublefact(2n+2l+1)=',doublefact(2*n+2*l+1)
           stop
        endif
      end function normho
cccccccccccccc factorial
          function fact(n)           
            implicit none
            integer n
            real*8 fact, dgamma,x
            x=dfloat(n+1)
            fact=dgamma(x)
          end function fact
ccccccc doublefact
      real*8 function doublefact(n)
            implicit none
            integer::n,i
            real::s
            if(n<0) then 
                write(*,*)'doublefact: n<0!',n
                stop
            else if(n==0.or.n==1) then 
                s=1.0
            else 
                s=n*1d0
                do i=n-2,1,-2
                    s=s*i
                end do
            end if 
            doublefact=s
      end function
ccccccc
      real*8 function ho3d(n,l,nu,r)          !3d hobasis of standard form
      implicit none                           !nu=mu*omega/(2hbar)
      integer l,n,p
      real*8::r,tp,norma,nu
      norma=normho(nu,n,l)
!        if (norma<1e-6) then
!           write(*,*)'ho3d: Norm=0!!!for  nu,n,l',nu,n,l
!        endif
      ho3d=norma*r**l*dexp(-nu*r**2)*
     & general_laguerre(2*nu*r**2,n,l+0.5d0) 
      end function ho3d
ccccccc   
c     Generalized Laguerre function L(n,l+1/2,x)
c *** -----------------------------------------------
      function laguerre(n,l,x)
        implicit none
        integer n,l,p
        real*8:: x,eps,tp,laguerre
        parameter(eps=1e-6)

        if(x<eps) x=eps
        tp=(-x)**n/fact(n)

        do p=1,n
!           write(*,*)'p=',p
           tp=tp-(n+l+1.5-p)*(n+1-p)/p/x*tp
        enddo
        laguerre=tp
        end function laguerre
ccccccc
      recursive function general_laguerre(x, k, alpha) result(L)
      real*8, intent(in) :: x
      integer, intent(in) :: k
      real*8, intent(in) :: alpha
      real*8 :: L

      if (k == 0) then
         L = 1.0
      else if (k == 1) then
         L = 1.0 + alpha - x
      else
         L = ((2*k - 1 + alpha - x) * general_laguerre(x, k-1, alpha) 
     &    - (k - 1 + alpha) * general_laguerre(x, k-2, alpha)) / k
      end if

      end function general_laguerre
ccccccc       
      real*8 function LSTFUN(gamma,m,r) !transformation function s(r)
      implicit none
      real*8::gamma,m,r
      LSTFUN=(1/((1.d0/r)**m+
     & (1.d0/gamma/sqrt(r))**m))**(1.d0/m)
      end function
ccccccc
      real*8 function D1LSTFUN(gamma,m,r)     !s'(r)
      implicit none
      real*8::gamma,m,r
      real*8::x,y
      x=(1.d0/r)**(m+1.d0)+1.d0/(2.0*r)*(1.d0/(gamma*sqrt(r)))**m
      y=(1.d0/r)**m+(1.d0/(gamma*sqrt(r)))**m
      D1LSTFUN=LSTFUN(gamma,m,r)*x/y
      end function
ccccccc
      real*8 function THOFUNC(n,l,alpha,gamma,m,r)    !alpha: dimensionless parameter of HOBASIS
      implicit none
      integer::n,l
      real*8::alpha,gamma,m,r
      THOFUNC=sqrt(D1LSTFUN(gamma,m,r))
     & *ho3d(n,l,alpha,LSTFUN(gamma,m,r))
     & *LSTFUN(gamma,m,r)/r 
      end function

        end module

