!rmatvar: variables for r-matrix calculation
!Cmat: Cmatrix before diagonalization C_{βn,β'm}, (1:Nbasis,1:Nbasis,1:Nch,1:Nch)
!Vcouple: coupled potential in Rmatrix calculation, with the size (1:Nbasis,1:Nbasis,1:Nch,1:Nch)
!T: kinetic energy matrix: (1:Nbasis,1:Nbasis,1:Nch), we put T and Bloch operator L(B) together in this array T
!B_i: Constants B_i in Bloch operator
!Ech: channel Energy, size: (1:Nbasis,1:Nbasis,1:Nch),here we set Ech matrix include E_i-E
!C: Cmatrix reconstructed for diagonalization, size: (1:Nbasis*Nch,1:Nbasis*Nch)
!Rmat: Rmatrix (R_{ββ'}),size (1:Nch,1:Nch)
!Z_O,Z_I: (k_ja)^{-1/2}(O_{l_i}(k_ia)delta_{ij}-k_jaR_{ij}O'_{L_j}(k_ja)),size (1:Nch,1:Nch)
!Smat: S-matrix, size (1:Nch,1:Nch)
!where Nch=beta%Nchmax calculated before in subroutine getchannelalphaD()
        module rmatvar
        implicit none
ccccccc
        complex(16),allocatable::Cmat(:,:,:,:)
        complex(16),allocatable::Vcouple(:,:,:,:)
        complex(16),allocatable::T(:,:,:)
        real(8),allocatable::B_i(:)
        complex(16),allocatable::Ech(:,:,:)
        complex(16),allocatable::C(:,:)
        complex(16),allocatable::Rmat(:,:)
        complex(16),allocatable::Z_O(:,:)
        complex(16),allocatable::Z_I(:,:)
        complex(16),allocatable::Smat(:,:)
        end module

ccccccc
!rmax:cut-off radius rmax=a
!nr: number of integral mesh points
!r,rw: integral radius vector, integral weight, which should be (1:Nr)
!rc: complex type r
        module mesh 
            implicit none
            real(8)::rmax
            integer::nr
            real(8),allocatable::r(:),rw(:)
            complex(16),allocatable::rc(:)
        end module
ccccccc
!2b sysmtem variables
        module system
            implicit none
            real(8)::mass1
            real(8)::mass2
            real(8)::mu
            real(8)::z1
            real(8)::z2
            real(8)::z12
        end module

!COUL90's variables
        module coulvar
        implicit none
ccccccc
        integer::KFN,IFAIL      !coul90 variables
        real*8::k_i,k_j !channel wave numbers
        real*8,allocatable::FC_i(:),GC_i(:),FCP_i(:),GCP_i(:) !coul90 variables for channel i
        real*8,allocatable::FC_j(:),GC_j(:),FCP_j(:),GCP_j(:) !coul90 variables for channel j
        complex*16::hlp_i,hln_i,dhln_i,dhlp_i   !hankel H^{\pm} and its derivative for channel i
        complex*16::hlp_j,hln_j,dhln_j,dhlp_j   !hankel H^{\pm} and its derivative for channel j
ccccccc
        end module
ccccccc
!whittaker function's variables
!ki: channel wave numbers k_i, eta=Z1Z2e^2mu/hbar^2/k_i,Sommerfeld parameter
!WTK,WTKP: Whittaker functions and their derivatives
        module wtkvar
        implicit none
        real(8)::ki,eta
        real(8),allocatable::WTK(:),WTKP(:)
        end module
ccccccc
!potvar: variables in the potential function
!str: type of potential, char type var
        module potvar
        implicit none
        character(len=7)::str
        real*8::v0,r0,a
        end module 
ccccccc
