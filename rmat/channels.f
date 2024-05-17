ccccccc
        module channels
        implicit none
!channel index for decay type(alpha), including        
!daughter nucleus spin sd, daughter nucleus excitation Ed,
!orbit angular momenta L, total(d+alpha) parity jpi
!total j(parent nucleus jt)
!minimum and maximum orbit angular momenta lmin,lmax
!jd is the parity of the daughter nucleus 
!where lmin=|j-s|,lmax=j+s
!channel number variable nch
!number of channels nchmax  
!allocatable Sc stores the spin of daughter nucleus of different channels 
!jc stores the total angular momentum of the 2b system J 
!lc stores the orbit L between the 2 bodies
!with size being the number of the channel
!daughter excitation Ec
        type channel
        real(8)::Sd
        real(8)::Ed
        real(8)::jt
        real(8)::jpi
        real(8)::jd
        integer::lmin,lmax
        integer::nch
        integer::nchmax
        end type
ccccccc
        real(8),allocatable::Sc(:)
        real(8),allocatable::jc(:)
        real(8),allocatable::lc(:)
        real(8),allocatable::Ec(:)
!channel notation β=|l Sd j>,2b system only requires one index β (1 configuration)
        type(channel)::beta 
!total energy E
        real(8)::E 
ccccccc
        contains
!this subroutine calculates the number of channels for the 2b decay of alpha type
!and gives each channel index beta, then stores them in alpha type allocatable vector
!where J,Sd are already known
!output: nchmax,S,j,l,E
        subroutine getchannelalphaD()
        implicit none
        integer::l
        integer::nch
ccccccc         
        beta%nchmax=0
ccccccc       
        beta%lmin=nint(abs(beta%jt-beta%Sd))
        beta%lmax=nint(beta%jt+beta%Sd)
        
        do l=beta%lmin,beta%lmax
             if((-1d0)**l*beta%jd==beta%jpi) then 
                beta%nchmax=beta%nchmax+1 
             end if
        end do
!allocate different channel index vector
!with size being the number of channels nchmax
        if (allocated(Sc)) deallocate(Sc)
        if (allocated(jc)) deallocate(jc)
        if (allocated(lc)) deallocate(lc)
        if (allocated(Ec)) deallocate(Ec)
ccccccc
        allocate(Sc(1:beta%nchmax))
        allocate(jc(1:beta%nchmax))
        allocate(lc(1:beta%nchmax))
        allocate(Ec(1:beta%nchmax))
ccccccc
        beta%nch=1
        do l=beta%lmin,beta%lmax
             if((-1d0)**l*beta%jd==beta%jpi) then 
                Sc(beta%nch)=beta%Sd
                jc(beta%nch)=beta%jt
                lc(beta%nch)=l
                Ec(beta%nch)=beta%Ed
ccccccc
                beta%nch=beta%nch+1
             end if
        end do
ccccccc
        end subroutine getchannelalphaD


        end module



