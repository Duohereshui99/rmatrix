ccccccc
!rmax:cut-off radius rmax=a
!nr: number of integral mesh points
!r,rw: integral radius vector, integral weight, which should be (1:Nr)
!phi: basis, (1:Nr,1:Nbasis,1:Nch), each channel (L) in line with Nbasis basis functions
!Cmat: Cmatrix before diagonalization C_{βn,β'm}, (1:Nbasis,1:Nbasis,1:Nch,1:Nch)
!Vcouple: coupled potential, with the size (1:Nbasis,1:Nbasis,1:Nch,1:Nch)
!T: kinetic energy matrix: (1:Nbasis,1:Nbasis,1:Nch)
!Blo: Bloch operator matrix elements, size: (1:Nbasis,1:Nbasis,1:Nch)
!Ech: channel Energy, size: (1:Nbasis,1:Nbasis,1:Nch)
!C: Cmatrix reconstructed for diagonalization, size: (1:Nbasis*Nch,1:Nbasis*Nch)
!C_minus: inverse of Cmatrix C, size the same as that of C: (1:Nbasis*Nch,1:Nbasis*Nch)
!Rmat: Rmatrix (R_{ββ'}),size (1:Nch,1:Nch)
!Z_O,Z_I: (k_ja)^{-1/2}(O_{l_i}(k_ia)delta_{ij}-k_jaR_{ij}O'_{L_j}(k_ja)),size (1:Nch,1:Nch)
!Smat: S-matrix, size (1:Nch,1:Nch)
!where Nch=beta%Nchmax calculated before in subroutine getchannelalphaD()
ccccccc
        module rmatmod
                use algorithm
                use basis
                use channels
                implicit none
                real(8)::rmax
                integer::nr
                real(8),allocatable::r(:),rw(:)
                complex(16),allocatable::phi(:,:,:)
                complex(16),allocatable::Cmat(:,:,:,:)
                complex(16),allocatable::Vcouple(:,:,:,:)
                complex(16),allocatable::T(:,:,:)
                complex(16),allocatable::Blo(:,:,:)
                complex(16),allocatable::Ech(:,:,:)
                complex(16),allocatable::C(:,:)
                complex(16),allocatable::C_minus(:,:)
                complex(16),allocatable::Rmat(:,:)
                complex(16),allocatable::Z_O(:,:)
                complex(16),allocatable::Z_I(:,:)
                complex(16),allocatable::Smat(:,:)
        contains
!this subroutine gives initializion for rmat 
!before rmat calculation, the initialization involves allocating the arrays 
!and give the basis functions, integral radius and weight, size: nr (also integral mesh number)
ccccccc
        subroutine rmat_int(nr,nbasis,rmax,alpha,gamma,m)
                implicit none
!input variables
                integer,intent(in)::nr,nbasis
                real(8),intent(in)::rmax
                real(8),intent(in)::alpha,gamma,m
!sum variables
                integer::i,j,k                                 
ccccccc
                call getchannelalphaD()
ccccccc
                if(allocated(r)) deallocate(r)
                if(allocated(rw)) deallocate(rw)
                if(allocated(phi)) deallocate(phi)
                if(allocated(Cmat)) deallocate(Cmat)
                if(allocated(Vcouple)) deallocate(Vcouple)
                if(allocated(T)) deallocate(T)
                if(allocated(Blo)) deallocate(Blo)
                if(allocated(Ech)) deallocate(Ech)
                if(allocated(C)) deallocate(C)
                if(allocated(C_minus)) deallocate(C_minus)
                if(allocated(Rmat)) deallocate(Rmat)
                if(allocated(Z_O)) deallocate(Z_O)
                if(allocated(Z_I)) deallocate(Z_I)
                if(allocated(Smat)) deallocate(Smat)
ccccccc
                allocate(r(1:nr),rw(1:nr))
                allocate(phi(1:nr,1:nbasis,1:beta%nchmax))
                allocate(Cmat(1:nbasis,1:nbasis,1:beta%nchmax,1:beta%nchmax))
                allocate(Vcouple(1:nbasis,1:nbasis,1:beta%nchmax,1:beta%nchmax))
                allocate(T(1:nbasis,1:nbasis,1:beta%nchmax))
                allocate(Blo(1:nbasis,1:nbasis,1:beta%nchmax))
                allocate(Ech(1:nbasis,1:nbasis,1:beta%nchmax))
                allocate(C(1:nbasis*beta%nchmax,1:nbasis*beta%nchmax))
                allocate(C_minus(1:nbasis*beta%nchmax,1:nbasis*beta%nchmax))
                allocate(Rmat(1:beta%nchmax,1:beta%nchmax),Smat(1:beta%nchmax,1:beta%nchmax))
                allocate(Z_O(1:beta%nchmax,1:beta%nchmax),Z_I(1:beta%nchmax,1:beta%nchmax))
ccccccc
                call gauleg(nr,0d0,rmax,r,rw)
ccccccc
!quantum number n should be chosen from 0, in line with i-1
                do i=1,beta%nchmax
                        do j=1,nbasis
                                do k=1,nr
                                        phi(k,j,i)=THOFUNC(j-1,int(beta%l(1),4),alpha,gamma,m,r(k)) 
                                end do
                        end do
                end do
        end subroutine rmat_int





        end module