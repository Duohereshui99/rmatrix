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
                real*8,intent(in)::rmax
                real*8,intent(in)::alpha,gamma,m
                real*8,intent(in)::E
                real*8,intent(in)::mu
                real*8,intent(in)::z1,z2
!sum variables
                integer::i,j,k                                 
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
                if(allocated(V1)) deallocate(V1)
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
                allocate(phi(1:nr,1:nbasis))
                allocate(d2phi(1:nr,1:nbasis))
                allocate(phi1(1:ndiff))
                allocate(d2phi1(1:ndiff))
                allocate(phia(1:nbasis))
                allocate(phipa(1:nbasis))
                allocate(Cmat(1:nbasis,1:nbasis,1:beta%nchmax,1:beta%nchmax))
                allocate(V1(1:ndiff))
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
                hcm1=hcm*cmplx(1d0,0d0,kind=8)
ccccccc
                do i=1,ndiff
                        rr(i)=hcm*i
                end do
ccccccc
                    do j=1,nbasis
                        do k=1,ndiff
                             phi1(k)=rr(k)*THOFUNC(j-1,0,alpha,gamma,m,rr(k)*cmplx(1d0,0d0,kind=8))                                          
                        end do  
ccccccc
                        do k=1,nr
                            phi(k,j)=FFC(r(k)/hcm,phi1,ndiff)
                        end do
ccccccc                              
                        call second_derivative(phi1,d2phi1,ndiff,hcm1)
ccccccc
                        do k=1,nr
                            d2phi(k,j)=FFC(r(k)/hcm,d2phi1,ndiff)        
                        end do
ccccccc
                        phia(j)=phi1(ndiff)
                        phipa(j)=(phi1(ndiff)-phi1(ndiff-1))/hcm
ccccccc
                    end do

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
!this subroutine gives the potential
        subroutine getpot(str)
            implicit none
            integer::i,j,k
            character(len=*)::str
ccccccc
            real*8::xx,zz,vcen,vtens,vls,hm,rmu
ccccccc
        if(allocated(Vc)) deallocate(Vc)
        allocate(Vc(1:nr,1:beta%nchmax,1:beta%nchmax))    
ccccccc
        select case(str)
ccccccc
            case('g')
                    do i=1,beta%nchmax
                        do j=1,beta%nchmax
!pot at uniformed mesh points
                            do k=1,ndiff
                                V1(k)=gausspot(rr(k),v0,r0,a)*cmplx(1d0,0d0,kind=8)
                            end do
!only channel diagonal elements pot                 
                        if (i==j) then                 
!then FFC -> interpolate   gauss mesh
                            do k=1,nr
                                Vc(k,i,j)=FFC(r(k)/hcm,V1,ndiff)
                            end do 
                        end if
ccccccc
                        end do
                    end do
ccccccc
!t: tensor force term included in the coupled pot for neutron-proton scattering
!!(only for 2 channels l=0,2)
            case('t') 
c Reid neutron-proton potential (T=1, soft core)
        rmu=mu/amu
        hm=20.736d0/rmu
        do i=1,nr
            xx=0.7d0*r(i)
            zz=exp(-xx)
            vcen=(-10.463d0*zz+105.468d0*zz**2-3187.8d0*zz**4+9924.3d0*zz**6)/xx
            vtens=-10.463d0*((1+3/xx+3/xx**2)*zz-(12/xx+3/xx**2)*zz**4)/xx+351.77d0*zz**4/xx-1673.5d0*zz**6/xx
            vls=708.91d0*zz**4/xx-2713.1d0*zz**6/xx
            Vc(i,1,1)=vcen-2*(beta%j_tot-1)*vtens/(2*beta%j_tot+1)+(beta%j_tot-1)*vls
            Vc(i,1,2)=6*vtens*sqrt(beta%j_tot*(beta%j_tot+1.0d0))/(2*beta%j_tot+1)
            Vc(i,2,1)=Vc(i,1,2)
            Vc(i,2,2)=vcen-2*(beta%j_tot+2)*vtens/(2*beta%j_tot+1)-(beta%j_tot+2)*vls
        end do
         !   Vc=Vc/hm
        end select
        end subroutine


ccccccc
        subroutine rmatrix()
!This subroutine implements the calculation of Rmatrix
! with input: 
! rmax: cut-off radius rmax=a
! Vc: coupled potential V_{ββ'}(r), size: (1:Nr,1:Nch,1:Nch),Vc is calculated additionlly       
! E:total energy
! mu: reduced mass
! z1,z2: charge number
ccccccc
        implicit none
        integer::i,j,k,mm,nn !sum variables
        integer::li,lj       !lc(i),lc(j)
ccccccc
!calculate C matrix
!C matrix is divided into 4 parts: kinetic T, Bloch operator B
!channel excitation Ei-E, and coupling potentials
ccccccc
ccccccc Ei-E matrix elements
        do mm=1,nbasis
            do nn=1,nbasis
                do i=1,beta%nchmax
                    Ech(mm,nn,i)=0d0
ccccccc
                    do k=1,nr
                        Ech(mm,nn,i)=Ech(mm,nn,i)+(Ec(i)-E)*conjg(phi(k,mm))*phi(k,nn)*rw(k)                     
                    end do     
ccccccc
                end do
            end do
        end do        
ccccccc coupled potential matrix elements Vcouple_{im,jn}
            do i=1,beta%nchmax
                do j=1,beta%nchmax
                    do mm=1,nbasis
                        do nn=1,nbasis
                            Vcouple(mm,nn,i,j)=0d0
                            do k=1,nr
                                Vcouple(mm,nn,i,j)=Vcouple(mm,nn,i,j)+conjg(phi(k,mm))*Vc(k,i,j)*phi(k,nn)*rw(k)   
                            end do
                    end do
                end do
            end do
        end do
ccccccc
!we put the Bloch operator L(B) and and kinetic energy T together into array T , since T+L(B) is Hermitian
!in real case, Hermitian means symmetric
        do mm=1,nbasis
            do nn=1,nbasis
                do i=1,beta%nchmax
                    do k=1,nr
                  T(mm,nn,i)=T(mm,nn,i)+(-hbarc**2/2/mu)*phi(k,mm)*d2phi(k,nn)*rw(k)  !rw(k)*conjg(phi(k,mm,i))*(-hbarc**2/2/mu*d2phi(k,nn,i))+
     &           +(hbarc**2/2/mu)*lc(i)*(lc(i)+1)/r(k)**2*phi(k,mm)*phi(k,nn)*rw(k)
                    end do
                T(mm,nn,i)=T(mm,nn,i)+hbarc**2/2/mu*conjg(phia(mm))*(phipa(nn)-B_i(i)/rmax*phia(nn))
                end do
            end do
        end do
!then we add the several matrix elements together to get Cmatrix before the reconstruction
        do mm=1,nbasis
            do nn=1,nbasis
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        if(i==j) then 
                        Cmat(mm,nn,i,j)=Ech(mm,nn,i)+T(mm,nn,i)
                        end if
                        Cmat(mm,nn,i,j)=Cmat(mm,nn,i,j)+Vcouple(mm,nn,i,j)
                    end do
                end do
            end do
        end do     
ccccccc
!reconstruct the Cmatrix for inversion, the size of C is (nbasis*beta%nchmax,nbasis*beta%nchmax)
ccccccc
        do mm=1,nbasis
            do nn=1,nbasis
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        C((i-1)*nbasis+mm,(j-1)*nbasis+nn)=Cmat(mm,nn,i,j)
                    end do
                end do
            end do
        end do
!get the inversion of Cmatrix, and the inversion is stored just in C.
        call mat_inv(C,nbasis*beta%nchmax,nbasis*beta%nchmax)
ccccccc
!Rmatrix , R_{ij}=hbar^2/(2mu a)*\sum_{mn}φ_n(a)(C^{-1})_{in,jm}φ_m(a)
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                Rmat(i,j)=0d0
                do mm=1,nbasis
                    do nn=1,nbasis
                        Rmat(i,j) = Rmat(i,j)+hbarc**2/2d0/mu/rmax*phia(mm)*C((i-1)*nbasis+mm,(j-1)*nbasis+nn)*phia(nn)                
                    end do
                end do
            end do
        end do
ccccccc
!Zmatrix: Z_O,Z_I, Smatrix: S=(Z_O)^{-1}Z_I
ccccccc
        KFN=0
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                k_i=sqrt(2d0*mu*abs(Ec(i)-E)/hbarc**2)
                k_j=sqrt(2d0*mu*abs(Ec(j)-E)/hbarc**2)
ccccccc
                li=int(lc(i),4)
                lj=int(lc(j),4)
ccccccc
                allocate(FC_i(0:li),GC_i(0:li),FCP_i(0:li),GCP_i(0:li))
                allocate(FC_j(0:lj),GC_j(0:lj),FCP_j(0:lj),GCP_j(0:lj))
ccccccc
                call COUL90(k_i*rmax,z1*z2*e2*mu/hbarc**2/k_i,0d0,li,FC_i,GC_i,FCP_i,GCP_i,KFN,IFAIL)
                call COUL90(k_j*rmax,z1*z2*e2*mu/hbarc**2/k_j,0d0,lj,FC_j,GC_j,FCP_j,GCP_j,KFN,IFAIL)
ccccccc H^{+}=G+iF, H^{-}=G-iF, and their derivatives
                hlp_i=cmplx(GC_i(li),FC_i(li),kind=8)
                hln_i=cmplx(GC_i(li),-FC_i(li),kind=8)
                dhlp_i=cmplx(GCP_i(li),FCP_i(li),kind=8)
                dhln_i=cmplx(GCP_i(li),-FCP_i(li),kind=8)
                hlp_j=cmplx(GC_j(lj),FC_j(lj),kind=8)
                hln_j=cmplx(GC_j(lj),-FC_j(lj),kind=8)
                dhlp_j=cmplx(GCP_j(lj),FCP_j(lj),kind=8)
                dhln_j=cmplx(GCP_j(lj),-FCP_j(lj),kind=8)
ccccccc
                Z_O(i,j)=(k_j*rmax)**(-0.5d0)*(hlp_i*delta(i,j)-k_j*rmax*Rmat(i,j)*dhlp_j)
                Z_I(i,j)=(k_j*rmax)**(-0.5d0)*(hln_i*delta(i,j)-k_j*rmax*Rmat(i,j)*dhln_j)
ccccccc
                deallocate(FC_i,GC_i,FCP_i,GCP_i)
                deallocate(FC_j,GC_j,FCP_j,GCP_j)
            end do
        end do
ccccccc
! Smatrix: S=(Z_O)^{-1}Z_I,first we get the inverse of Z_O
                call mat_inv(Z_O,beta%nchmax,beta%nchmax)
                Smat=matmul(Z_O,Z_I)
ccccccc
                write(*,*) 'Smatrix:'
                do i=1,beta%nchmax
                    write(*,*) Smat(i,:)
                end do
                write(*,*) 'the module of S'
                do i=1,beta%nchmax
                    write(*,*) abs(Smat(i,:))
                end do
                write(*,*) 'the phase of S'
                do i=1,beta%nchmax
                        write(*,*) atan2(aimag(Smat(i,:)),real(Smat(i,:)))/2
                end do

        end subroutine
ccccccc

        end module