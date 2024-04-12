ccccccc
        module channels
        implicit none
!channel index for decay type(alpha), including        
!daughter nucleus spin sd, daughter nucleus excitation Ed,
!orbit angular momenta L, total(d+alpha) parity jpi
!total j(parent nucleus j)
!minimum and maximum orbit angular momenta lmin,lmax
!jd is the parity of the daughter nucleus 
!where lmin=|j-s|,lmax=j+s
!channel number variable nch
!number of channels nchmax  
!allocatable S stores the spin of daughter nucleus of different channels 
!j stores the total angular momentum of the 2b system J 
!l stores the orbit L between the 2 bodies
!with size being the number of the channel
!daughter excitation E
        type channel
        real(8)::Sd
        real(8)::Ed
        real(8)::jt
        real(8)::jpi
        real(8)::jd
        integer::lmin,lmax
        integer::nch
        integer::nchmax
        real(8),allocatable::S(:)
        real(8),allocatable::j(:)
        real(8),allocatable::l(:)
        real(8),allocatable::E(:)
        end type
!channel notation β=|l Sd j>,2b system only requires one index β (1 configuration)
        type(channel)::beta 
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
        if (allocated(beta%S)) deallocate(beta%S)
        if (allocated(beta%j)) deallocate(beta%j)
        if (allocated(beta%l)) deallocate(beta%l)
        if (allocated(beta%E)) deallocate(beta%E)
ccccccc
        allocate(beta%S(1:beta%nchmax))
        allocate(beta%j(1:beta%nchmax))
        allocate(beta%l(1:beta%nchmax))
        allocate(beta%E(1:beta%nchmax))
ccccccc
        beta%nch=1
        do l=beta%lmin,beta%lmax
             if((-1d0)**l*beta%jd==beta%jpi) then 
                beta%S(beta%nch)=beta%Sd
                beta%j(beta%nch)=beta%jt
                beta%l(beta%nch)=l
                beta%E(beta%nch)=beta%Ed
ccccccc
                beta%nch=beta%nch+1
             end if
        end do
ccccccc
        end subroutine getchannelalphaD


        end module



