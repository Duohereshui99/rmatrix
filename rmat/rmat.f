ccccccc        
        program main
            use rmatmod 
            implicit none
ccccccc
            integer::i,j,k
            complex(16),allocatable::Vc(:,:,:)
            complex(16)::sum
ccccccc
            namelist /channelbeta/ beta
            namelist /meshs/ nr,ndiff,rmax
            namelist /basisvar/ nbasis,b,gamma,m
            namelist /systems/  E,mass1,mass2,z1,z2
            namelist /potentials/ str,v0,r0,a
ccccccc
            call get_info()
ccccccc
            read(5,nml=channelbeta)
            read(5,nml=meshs)
            read(5,nml=basisvar)
            read(5,nml=systems)
            read(5,nml=potentials)
ccccccc
            write(*,*) 'testmessages'
ccccccc
            z12=z1*z2
            mu=amu*(mass1*mass2/(mass1+mass2))
            alpha=1d0/2d0/b**2
ccccccc
        call rmat_int(nr,nbasis,rmax,alpha,gamma,m,E,mu,z1,z2)
ccccccc
        if(allocated(Vc)) deallocate(Vc)
        allocate(Vc(1:nr,1:beta%nchmax,1:beta%nchmax))    
ccccccc
        select case(str)
ccccccc
            case('g')
                    do i=1,beta%nchmax
                            do j=1,beta%nchmax
                                    do k=1,nr
                                            Vc(k,i,j)=gausspot(r(k),v0,r0,a)*cmplx(1d0,0d0)
                                    end do 
                            end do
                    end do
        end select
ccccccc
        call rmatrix(rmax,Vc,E,mu,z1,z2)
ccccccc
                do i=1,beta%nchmax
                        write(27,*) Smat(i,:)
                end do
ccccccc
                do k=1,nr
                        write(29,*) r(k),real((d2phi(k,1,1)))
                        write(30,*) r(k),r(k)*real(THOFUNC(0,0,alpha,gamma,m,r(k)*cmplx(1d0,0d0)))
                end do
ccccccc
                sum=0
                do k=1,nr
                        sum=sum+rw(k)*conjg(phi(k,1,1))*phi(k,2,1)
                end do
                write(*,*) sum
ccccccc

                do k=1,nr
                        write(31,*)  r(k),-hbarc**2/2/mu*real(conjg(phi(k,1,1))*(d2phi(k,1,1)))
                end do



                deallocate(Vc)
ccccccc
!!
!!
!!----------------------------------------------------------------------------------
!!
!!
            contains
ccccccc
        subroutine rmatrix(rmax,Vc,E,mu,z1,z2)
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
ccccccc
        real(8),intent(in)::rmax
        complex(16),intent(in)::Vc(1:nr,1:beta%nchmax,1:beta%nchmax)
        complex(16)::CC(nbasis*beta%nchmax,nbasis*beta%nchmax),CI(nbasis*beta%nchmax,nbasis*beta%nchmax)
        real(8),intent(in)::E   
        real(8),intent(in)::mu
        real(8),intent(in)::z1,z2
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
                                do j=1,nr
                                        Ech(mm,nn,i)=Ech(mm,nn,i)+(Ec(i)-E)*conjg(phi(j,mm,i))*phi(j,nn,i)*rw(j)                             
                                end do
                        end do
                end do
        end do        
ccccccc coupled potential matrix elements Vcouple_{im,jn}
        do mm=1,nbasis
                do nn=1,nbasis
                        do i=1,beta%nchmax
                                do j=1,beta%nchmax
                                        Vcouple(mm,nn,i,j)=0d0
                                        do k=1,nr
                                             Vcouple(mm,nn,i,j)=Vcouple(mm,nn,i,j)+conjg(phi(k,mm,i))*Vc(k,i,j)*phi(k,nn,j)*rw(k)   
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
                                T(mm,nn,i)=0d0
                                do k=1,nr
                T(mm,nn,i)=T(mm,nn,i)+rw(k)*conjg(phi(k,mm,i))*(-hbarc**2/2/mu*d2phi(k,nn,i)+lc(i)*(lc(i)+1)/r(k)**2*phi(k,nn,i))
                                end do
                T(mm,nn,i)=T(mm,nn,i)+hbarc**2/2/mu*(conjg(phia(mm,i))*phipa(nn,i)-conjg(phia(mm,i))*B_i(i)/rmax*phia(nn,i))
                        end do
                end do
        end do
ccccccc
100      format(5F18.6)
        do i=1,nbasis
                write(789,100) real(T(i,:,1))
        end do
ccccccc
!then we add the several matrix elements together to get Cmatrix before the reconstruction
        do mm=1,nbasis
                do nn=1,nbasis
                        do i=1,beta%nchmax
                                do j=1,beta%nchmax
                                        Cmat(mm,nn,i,j)=Ech(mm,nn,i)+T(mm,nn,i)+Vcouple(mm,nn,i,j)
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
                                        C(mm*i,nn*j)=Cmat(mm,nn,i,j)
                                end do
                        end do
                end do
        end do
        CC=C

        do i=1,nbasis*beta%nchmax
                write(666,*) (CC(i,:))
        end do
ccccccc
!get the inversion of Cmatrix, and the inversion is stored just in C.
        call mat_inv(C,nbasis*beta%nchmax,nbasis*beta%nchmax)
ccccccc
        CI=matmul(C,CC)      
        do i=1,nbasis*beta%nchmax
                write(667,*) (C(i,:))
        end do
        do i=1,nbasis*beta%nchmax
               write(777,*) CI(i,:)
        end do
!Rmatrix , R_{ij}=hbar^2/(2mu a)*\sum_{mn}φ_n(a)(C^{-1})_{in,jm}φ_m(a)
        do i=1,beta%nchmax
                do j=1,beta%nchmax
                        Rmat(i,j)=0d0
                        do mm=1,nbasis
                                do nn=1,nbasis
                        Rmat(i,j)=Rmat(i,j)+hbarc**2/2d0/mu/rmax
     &          *phia(mm,i)*C(mm*i,nn*j)*phia(nn,j)                
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
                        allocate(FC_i(0:lc(i)),GC_i(0:lc(i)),FCP_i(0:lc(i)),GCP_i(0:lc(i)))
                        allocate(FC_j(0:lc(j)),GC_j(0:lc(j)),FCP_j(0:lc(j)),GCP_j(0:lc(j)))
ccccccc
                        call COUL90(k_i*rmax,z1*z2*e2*mu/hbarc**2/k_i,0d0,int(lc(i),4),FC_i,GC_i,FCP_i,GCP_i,KFN,IFAIL)
                        call COUL90(k_j*rmax,z1*z2*e2*mu/hbarc**2/k_j,0d0,int(lc(j),4),FC_j,GC_j,FCP_j,GCP_j,KFN,IFAIL)
ccccccc H^{+}=G+iF, H^{-}=G-iF, and their derivatives
                        hlp_i=cmplx(GC_i(lc(i)),FC_i(lc(i)))
                        hln_i=cmplx(GC_i(lc(i)),-FC_i(lc(i)))
                        dhlp_i=cmplx(GCP_i(lc(i)),FCP_i(lc(i)))
                        dhln_i=cmplx(GCP_i(lc(i)),-FCP_i(lc(i)))
                        hlp_j=cmplx(GC_j(lc(j)),FC_j(lc(j)))
                        hln_j=cmplx(GC_j(lc(j)),-FC_j(lc(j)))
                        dhlp_j=cmplx(GCP_j(lc(j)),FCP_j(lc(j)))
                        dhln_j=cmplx(GCP_j(lc(j)),-FCP_j(lc(j)))
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
        end subroutine
ccccccc
! 
! 
! 
! 
! 
!
ccccccc
            subroutine get_info()
#ifdef BASE
        print *, 'Base directory: ', BASE
#endif

#ifdef VERDATE
        print *, 'Version date: ', VERDATE
#endif

#ifdef VERREV
        print *, 'Version revision: ', VERREV
#endif

#ifdef COMPDATE
        print *, 'Compilation date: ', COMPDATE
#endif
            end subroutine            
        end program 