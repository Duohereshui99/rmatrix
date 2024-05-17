ccccccc
        module rmatmod
                use algorithm
                use basis
                use parameter
                use channels
                use whittaker
                use matinv
                use coulfunc
                use deltaf
                use coulvar
                use wtkvar
                use rmatvar
                use mesh
                use system
                use potential
                use potvar
                implicit none
        contains
!this subroutine gives initializion for rmat 
!before rmat calculation, the initialization involves allocating the arrays 
!and the basis functions, integral radius and weight, size: nr (also integral mesh number)
!E: total E
!mu: reduced mass, z1,z2: charge number
ccccccc
        subroutine rmat_int(nr,nbasis,rmax,alpha,gamma,m,E,mu,z1,z2)
                implicit none
!input variables
                integer,intent(in)::nr,nbasis
                real(8),intent(in)::rmax
                real(8),intent(in)::alpha,gamma,m
                real(8),intent(in)::E
                real(8),intent(in)::mu
                real(8),intent(in)::z1,z2
!sum variables
                integer::i,j,k                                 
ccccccc
                call getchannelalphaD()
ccccccc
                if(allocated(WTK)) deallocate(WTK)
                if(allocated(WTKP)) deallocate(WTKP)
                if(allocated(r)) deallocate(r)
                if(allocated(rr)) deallocate(rr)
                if(allocated(rw)) deallocate(rw)
                if(allocated(phi)) deallocate(phi)
                if(allocated(d2phi)) deallocate(d2phi)
                if(allocated(phi1)) deallocate(phi1)
                if(allocated(d2phi1)) deallocate(d2phi1)
                if(allocated(phia)) deallocate(phia)
                if(allocated(phipa)) deallocate(phipa)
                if(allocated(Cmat)) deallocate(Cmat)
                if(allocated(Vcouple)) deallocate(Vcouple)
                if(allocated(T)) deallocate(T)
                if(allocated(B_i)) deallocate(B_i)
                if(allocated(Ech)) deallocate(Ech)
                if(allocated(C)) deallocate(C)
                if(allocated(Rmat)) deallocate(Rmat)
                if(allocated(Z_O)) deallocate(Z_O)
                if(allocated(Z_I)) deallocate(Z_I)
                if(allocated(Smat)) deallocate(Smat)
ccccccc
                allocate(WTK(1:beta%nchmax+1),WTKP(1:beta%nchmax+1))
                allocate(r(1:nr),rw(1:nr))
                allocate(rr(1:ndiff))
                allocate(phi(1:nr,1:nbasis,1:beta%nchmax))
                allocate(d2phi(1:nr,1:nbasis,1:beta%nchmax))
                allocate(phi1(1:ndiff))
                allocate(d2phi1(1:ndiff))
                allocate(phia(1:nbasis,1:beta%nchmax))
                allocate(phipa(1:nbasis,1:beta%nchmax))
                allocate(Cmat(1:nbasis,1:nbasis,1:beta%nchmax,1:beta%nchmax))
                allocate(Vcouple(1:nbasis,1:nbasis,1:beta%nchmax,1:beta%nchmax))
                allocate(T(1:nbasis,1:nbasis,1:beta%nchmax))
                allocate(B_i(1:beta%nchmax))
                allocate(Ech(1:nbasis,1:nbasis,1:beta%nchmax))
                allocate(C(1:nbasis*beta%nchmax,1:nbasis*beta%nchmax))
                allocate(Rmat(1:beta%nchmax,1:beta%nchmax),Smat(1:beta%nchmax,1:beta%nchmax))
                allocate(Z_O(1:beta%nchmax,1:beta%nchmax),Z_I(1:beta%nchmax,1:beta%nchmax))
ccccccc
                call gauleg(nr,0d0,rmax,r,rw)
ccccccc
!quantum number n should be chosen from 0, in line with i-1
!and this is the basis function initialization
!initialization should also include the 2nd derivative of basis functions 
!and 1st boundary derivative of basis
!d/dr (rR)=R+rR'
!d^2/dr^2 (rR)=2R'+rR''
ccccccc
                hcm=rmax/ndiff
                hcm1=rmax/ndiff*cmplx(1d0,0d0)
ccccccc
                do i=1,ndiff
                        rr(i)=hcm*i
                end do
ccccccc
                do i=1,beta%nchmax
                        do j=1,nbasis
ccccccc
                                do k=1,nr
        phi(k,j,i)=r(k)*THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,r(k)*cmplx(1d0,0d0))    
                                end do
ccccccc
                                do k=1,ndiff
                                phi1(k)=rr(k)*THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,rr(k)*cmplx(1d0,0d0))                                          
                                end do  
ccccccc                              
                                call  second_derivative(phi1,d2phi1,ndiff,hcm*cmplx(1d0,0d0))
ccccccc
                                do k=1,nr
                                d2phi(k,j,i)=FFC(r(k)/hcm,d2phi1,ndiff)        
                                end do
ccccccc
                         phia(j,i)=THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,rmax*cmplx(1d0,0d0))
                         phipa(j,i)=deriv1(phi1,hcm1,ndiff,ndiff-1)
ccccccc
                        end do
                end do
ccccccc
!                 do i=1,beta%nchmax
!                         do j=1,nbasis

! !d/dr (rR)=R+rR'
                
!                 !THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,rmax*cmplx(1d0,0d0))+
!      !&           rmax*D1THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,rmax*cmplx(1d0,0d0))
!                         end do
!                 end do
ccccccc
!constants B_i in Bloch operator
ccccccc
                do i=1,beta%nchmax
                      if(E>Ec(i)) then
                                B_i(i)=0d0
                        else
                        ki=sqrt(2d0*mu*abs(E-Ec(i))/hbarc**2)
                        eta=z1*z2*e2*mu/hbarc**2/ki
                        call WHIT(eta,rmax,ki,E,int(lc(i),4),WTK,WTKP,0)
                                B_i(i)=2*ki*rmax*WTKP(i)/WTK(i)
                      end if
                end do
ccccccc
        end subroutine rmat_int
ccccccc


        end module