ccccccc        
        program main
            use channels
            use rmatmod
            use basis
            implicit none
ccccccc
            integer::i
ccccccc
            call get_info()
ccccccc
            write(*,*) 'testmessages'
            beta%jpi=1d0
            beta%jd=1d0
            beta%Sd=2d0
            beta%jt=2d0
            beta%Ed=0d0
        !     write(*,*) 'L','S','J','E'
            
        !     call getchannelalphaD()
        !     do i=1,beta%nchmax
        !         write(*,*) beta%L(i),beta%S(i),beta%J(i),beta%E(i)
        !     end do
        !     write(*,*) beta%nchmax
        !     write(*,*) beta%lmin
        !     write(*,*) beta%lmax

        !     rmax=10
        !     nr=300
        !     call rmat_int(nr,20,beta%nchmax,rmax)
        !     do i=1,nr
        !         write(*,*) 'mesh r',r(i),rw(i)
        !     end do
         !   rmat_int(nr,nbasis,rmax,alpha,gamma,m)
                call rmat_int(200,5,20d0,1d0,1.5d0,4d0)
                do i=1,200
                    write(14,*) r(i),real(phi(i,1,1))
                end do


            contains
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