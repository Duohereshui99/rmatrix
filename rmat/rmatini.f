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
                if(allocated(rw)) deallocate(rw)
                if(allocated(phi)) deallocate(phi)
                if(allocated(d2phi)) deallocate(d2phi)
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
                allocate(phi(1:nr,1:nbasis,1:beta%nchmax))
                allocate(d2phi(1:nr,1:nbasis,1:beta%nchmax))
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
                do i=1,beta%nchmax
                        do j=1,nbasis
                                do k=1,nr
                        phi(k,j,i)=THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,r(k)*cmplx(1d0,0d0))
                        d2phi(k,j,i)=D2THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,r(k)*cmplx(1d0,0d0)) 
                                end do
                        end do
                end do
ccccccc
                do i=1,beta%nchmax
                        do j=1,nbasis
                phia(j,i)=THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,rmax*cmplx(1d0,0d0))
                phipa(j,i)=D1THOFUNC(j-1,int(lc(i),4),alpha,gamma,m,rmax*cmplx(1d0,0d0))
                        end do
                end do
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