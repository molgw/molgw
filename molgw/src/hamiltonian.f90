!=========================================================================
subroutine setup_overlap(print_matrix_,basis,s_matrix)
 use m_definitions
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(basis%nbf,basis%nbf)
!====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
!====

 write(stdout,*) 'Setup overlap matrix S'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call overlap_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo
     s_matrix(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                                                   MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) )


     deallocate(matrix_cart)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='=== Overlap matrix S ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix)


end subroutine setup_overlap


!=========================================================================
subroutine setup_kinetic(print_matrix_,basis,hamiltonian_kinetic)
 use m_definitions
 use m_timing
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(basis%nbf,basis%nbf)
!====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
!====

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call kinetic_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo
     hamiltonian_kinetic(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) )


     deallocate(matrix_cart)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='===  Kinetic energy contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic)

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic


!=========================================================================
subroutine setup_nucleus(print_matrix_,basis,hamiltonian_nucleus)
 use m_definitions
 use m_timing
 use m_basis_set
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
!====

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian'
 if( nproc > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank /= MODULO(iatom-1,nproc) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     matrix_cart(:,:) = 0.0_dp
     do iatom=1,natom
       if( rank /= MODULO(iatom-1,nproc) ) cycle
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           call nucleus_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),zatom(iatom),x(:,iatom),vnucleus_ij)
           matrix_cart(i_cart,j_cart) = matrix_cart(i_cart,j_cart) + vnucleus_ij
         enddo
       enddo
     enddo
     hamiltonian_nucleus(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) ) 


     deallocate(matrix_cart)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 !
 ! Reduce operation
 call xsum(hamiltonian_nucleus)

 title='===  Nucleus potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus


!=========================================================================
subroutine setup_effective_core(print_matrix_,basis,hamiltonian_nucleus)
 use m_definitions
 use m_timing
 use m_basis_set
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!====
 integer              :: natom_local
 integer              :: iproj,jbf
 integer              :: iproj_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vecp_proji_j
!FBFB
 logical,parameter    :: normalized=.FALSE.
 type(basis_set)      :: proj
 real(dp)             :: x0(3)
 real(dp)             :: alpha(1)
 real(dp)             :: coeff(1)
 integer              :: shell_index
 real(dp),allocatable :: projector_l(:,:)
 real(dp),allocatable :: projector_r(:,:)
 real(dp),allocatable :: proj_coef(:)
!====

 call start_clock(timing_hamiltonian_ecp)
 write(stdout,'(/,a)') ' Setup effective core potential part of the Hamiltonian'

 proj%nbf      =   1 ! + 3 + 6
 proj%nbf_cart =   1 ! + 3 + 6
 proj%gaussian_type = 'PURE'
 proj%ammax    = 1
 allocate(proj%bf(proj%nbf_cart))

 allocate(projector_l(basis%nbf,proj%nbf)) ! FBFB deallocate
 allocate(projector_r(proj%nbf,basis%nbf)) ! FBFB deallocate
 allocate(proj_coef(proj%nbf))             ! FBFB deallocate
 
 !
 ! Creating a temporary basis that contains the projectors
 shell_index = 0
 iproj = 0
 do iatom=1,natom
   x0(:) = 0.0_dp

#if 0
   shell_index = shell_index + 1
   alpha(1) = 2.653000000 * 0.50_dp
   coeff(1) = 1.0_dp 
   iproj = iproj + 1
   proj_coef(iproj) = 13.325000000 
   call init_basis_function(normalized,1,0,0,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
   proj%bf(iproj)%g(1)%norm_factor = 1.0_dp / SQRT( 4.0_dp * pi )
   write(*,*) proj%bf(iproj)%g(1)%norm_factor

   shell_index = shell_index + 1
   alpha(1) = 3.120000000 * 0.50_dp
   coeff(1) = 1.0_dp 
   iproj = iproj + 1
   proj_coef(iproj) = -1.574000000
   call init_basis_function(normalized,1,1,0,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
   proj%bf(iproj)%g(1)%norm_factor = SQRT(3.0_dp) / SQRT( 4.0_dp * pi )
   write(*,*) proj%bf(iproj)%g(1)%norm_factor

   iproj = iproj + 1
   proj_coef(iproj) =  -1.574000000
   call init_basis_function(normalized,1,0,1,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
   proj%bf(iproj)%g(1)%norm_factor = SQRT(3.0_dp) / SQRT( 4.0_dp * pi )
   write(*,*) proj%bf(iproj)%g(1)%norm_factor

   iproj = iproj + 1
   proj_coef(iproj) = -1.574000000
   call init_basis_function(normalized,1,0,0,1,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
   proj%bf(iproj)%g(1)%norm_factor = SQRT(3.0_dp) / SQRT( 4.0_dp * pi )
   write(*,*) proj%bf(iproj)%g(1)%norm_factor
#else
! s proj
   shell_index = shell_index + 1
   alpha(1) = 0.73277 ! 1.732000000 * 0.50_dp
   coeff(1) = 1.0_dp 
   iproj = iproj + 1
   proj_coef(iproj) = 1.0 ! 14.676000000 
   call init_basis_function(normalized,1,0,0,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
   proj%bf(iproj)%g(1)%norm_factor = 1.0_dp / SQRT( 4.0_dp * pi )
   write(*,*) proj%bf(iproj)%g(1)%norm_factor

! ! p proj
!    shell_index = shell_index + 1
!    alpha(1) = 1.115000000 * 0.50_dp
!    coeff(1) = 1.0_dp 
!    iproj = iproj + 1
!    proj_coef(iproj) = 5.175700000
!    call init_basis_function(normalized,1,1,0,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(3.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
! 
!    iproj = iproj + 1
!    proj_coef(iproj) = 5.175700000
!    call init_basis_function(normalized,1,0,1,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(3.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
! 
!    iproj = iproj + 1
!    proj_coef(iproj) = 5.175700000
!    call init_basis_function(normalized,1,0,0,1,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(3.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
! 
! ! d proj
!    shell_index = shell_index + 1
!    alpha(1) = 1.203000000 * 0.50_dp
!    coeff(1) = 1.0_dp 
!    iproj = iproj + 1
!    proj_coef(iproj) = -1.816000000
!    call init_basis_function(normalized,1,2,0,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(15.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
!    iproj = iproj + 1
!    proj_coef(iproj) = -1.816000000
!    call init_basis_function(normalized,1,1,1,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(15.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
!    iproj = iproj + 1
!    proj_coef(iproj) = -1.816000000
!    call init_basis_function(normalized,1,1,0,1,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(15.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
!    iproj = iproj + 1
!    proj_coef(iproj) = -1.816000000
!    call init_basis_function(normalized,1,0,2,0,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(15.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
!    iproj = iproj + 1
!    proj_coef(iproj) = -1.816000000
!    call init_basis_function(normalized,1,0,1,1,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(15.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor
!    iproj = iproj + 1
!    proj_coef(iproj) = -1.816000000
!    call init_basis_function(normalized,1,0,0,2,iatom,x0,alpha,coeff,shell_index,proj%bf(iproj))
!    proj%bf(iproj)%g(1)%norm_factor = SQRT(15.0_dp) / SQRT( 4.0_dp * pi )
!    write(*,*) proj%bf(iproj)%g(1)%norm_factor


#endif

 enddo


 !
 ! Index iproj is running over the projectors in the ECP
 ! It is always a 'PURE' Gaussian
 iproj_cart = 1
 jbf_cart = 1
 iproj    = 1
 jbf      = 1
 do while(iproj_cart<=proj%nbf_cart)
   li      = proj%bf(iproj_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am('PURE',li)
   write(*,*) 'li',li
   write(*,'(3(f12.6,2x))') cart_to_pure(li)%matrix(:,:)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     matrix_cart(:,:) = 0.0_dp
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call overlap_basis_function(proj%bf(iproj_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo

     projector_r(iproj:iproj+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                                                         MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) ) 


     deallocate(matrix_cart)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   iproj      = iproj      + ni
   iproj_cart = iproj_cart + ni_cart

 enddo

 forall(iproj=1:proj%nbf)
   projector_l(:,iproj) = proj_coef(iproj) * projector_r(iproj,:)
 end forall
 
 hamiltonian_nucleus(:,:) = hamiltonian_nucleus(:,:) + MATMUL(  projector_l(:,:) , projector_r(:,:) ) 

 title='===  Effective core potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_ecp)

end subroutine setup_effective_core


!=========================================================================
subroutine setup_hartree(print_matrix_,nbf,nspin,p_matrix,pot_hartree,ehartree)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_hartree(nbf,nbf)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Hartree term'
 call start_clock(timing_hartree)

 pot_hartree(:,:)=0.0_dp

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
 do jbf=1,nbf
   do ibf=1,nbf
     if( negligible_basispair(ibf,jbf) ) cycle
     do lbf=1,nbf
       !
       ! symmetry k <-> l
       do kbf=1,lbf-1 ! nbf
         if( negligible_basispair(kbf,lbf) ) cycle
         !
         ! symmetry (ij|kl) = (kl|ij) has been used to loop in the fast order
         pot_hartree(ibf,jbf) = pot_hartree(ibf,jbf) &
                    + eri(kbf,lbf,ibf,jbf) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
       enddo
       pot_hartree(ibf,jbf) = pot_hartree(ibf,jbf) &
                  + eri(lbf,lbf,ibf,jbf) * SUM( p_matrix(lbf,lbf,:) )
     enddo
   enddo
 enddo
!$OMP END DO
!$OMP END PARALLEL


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,pot_hartree)

 ehartree = 0.5_dp*SUM(pot_hartree(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(pot_hartree(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree


!=========================================================================
subroutine setup_hartree_ri(print_matrix_,nbf,nspin,p_matrix,pot_hartree,ehartree)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_hartree(nbf,nbf)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 integer              :: ibf_auxil,ipair
 integer              :: index_ij,index_kl
 real(dp),allocatable :: partial_sum(:)
 real(dp)             :: rtmp
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity'
 call start_clock(timing_hartree)


 allocate(partial_sum(nauxil_3center))
 partial_sum(:) = 0.0_dp
 do ipair=1,npair
   kbf = index_basis(1,ipair)
   lbf = index_basis(2,ipair)
   ! Factor 2 comes from the symmetry of p_matrix
   partial_sum(:) = partial_sum(:) + eri_3center(:,ipair) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
   ! Then diagonal terms have been counted twice and should be removed once.
   if( kbf == lbf ) &
     partial_sum(:) = partial_sum(:) - eri_3center(:,ipair) * SUM( p_matrix(kbf,kbf,:) )
 enddo

 ! Hartree potential is not sensitive to spin
 pot_hartree(:,:) = 0.0_dp
 do ipair=1,npair
   ibf = index_basis(1,ipair)
   jbf = index_basis(2,ipair)
   rtmp = DOT_PRODUCT( eri_3center(:,ipair) , partial_sum(:) )
   pot_hartree(ibf,jbf) = rtmp
   pot_hartree(jbf,ibf) = rtmp
 enddo

 deallocate(partial_sum)

 !
 ! Sum up the different contribution from different procs only if needed
 call xsum(pot_hartree)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,pot_hartree)

 ehartree = 0.5_dp*SUM(pot_hartree(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(pot_hartree(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree_ri


!=========================================================================
subroutine setup_exchange(print_matrix_,nbf,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 use m_inputparam,only: nspin,spin_fact
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Exchange term'
 call start_clock(timing_exchange)


 pot_exchange(:,:,:)=0.0_dp

 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO 
   do jbf=1,nbf
     do lbf=1,nbf
       if( negligible_basispair(lbf,jbf) ) cycle
       do kbf=1,nbf
!         if( ABS(p_matrix(kbf,lbf,ispin)) <  1.0e-12_dp ) cycle 
         do ibf=1,nbf
           if( negligible_basispair(ibf,kbf) ) cycle
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           pot_exchange(ibf,jbf,ispin) = pot_exchange(ibf,jbf,ispin) &
                      - eri(ibf,kbf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo

 title='=== Exchange contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange


!=========================================================================
subroutine setup_exchange_ri(print_matrix_,nbf,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 use m_inputparam,only: nspin,spin_fact
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: p_matrix_sqrt(nbf,nbf)
 real(dp),allocatable :: p_matrix_sqrt_occ(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
!=====

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)


 pot_exchange(:,:,:)=0.0_dp

 allocate(tmp(nauxil_3center,nbf))

 do ispin=1,nspin

   !
   ! Calculate the square-root of the density matrix
   ! This is obtained with a diagonalization
   !
   ! Maybe a Cholesky decomposition could work too in the future
   ! but right now the p_matrix has some noise and it is not exactly positive
   ! semi-definite...
   p_matrix_sqrt(:,:) = p_matrix(:,:,ispin)

#ifndef HAVE_SCALAPACK
   call diagonalize(nbf,p_matrix_sqrt,eigval)
#else
   call diagonalize_scalapack(nbf,p_matrix_sqrt,eigval)
#endif

   ! Denombrate the strictly positive eigenvalues
   nocc = COUNT( eigval(:) > completely_empty )

   allocate(p_matrix_sqrt_occ(nbf,nocc))
   ibf = 0
   do jbf=1,nbf
     if( eigval(jbf) > completely_empty ) then
       ibf = ibf + 1
       p_matrix_sqrt_occ(:,ibf) = p_matrix_sqrt(:,jbf) * SQRT( eigval(jbf) )
     endif
   enddo

   do istate=1,nocc
     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + p_matrix_sqrt_occ(jbf,istate) * eri_3center(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + p_matrix_sqrt_occ(ibf,istate) * eri_3center(:,ipair)
     enddo

     pot_exchange(:,:,ispin) = pot_exchange(:,:,ispin) &
                        - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
   enddo

   deallocate(p_matrix_sqrt_occ)

 enddo
 deallocate(tmp)

 call xsum(pot_exchange)

 call dump_out_matrix(print_matrix_,'=== Exchange contribution ===',nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri


!=========================================================================
subroutine setup_exchange_longrange_ri(print_matrix_,nbf,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 use m_inputparam,only: nspin,spin_fact
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: p_matrix_sqrt(nbf,nbf)
 real(dp),allocatable :: p_matrix_sqrt_occ(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
!=====

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)


 pot_exchange(:,:,:)=0.0_dp

 allocate(tmp(nauxil_3center_lr,nbf))

 do ispin=1,nspin

   !
   ! Calculate the square-root of the density matrix
   ! This is obtained with a diagonalization
   !
   ! Maybe a Cholesky decomposition could work too in the future
   ! but right now the p_matrix has some noise and it is not exactly positive
   ! semi-definite...
   p_matrix_sqrt(:,:) = p_matrix(:,:,ispin)
   call diagonalize(nbf,p_matrix_sqrt,eigval)

   ! Denombrate the strictly positive eigenvalues
   nocc = COUNT( eigval(:) > completely_empty )

   allocate(p_matrix_sqrt_occ(nbf,nocc))
   ibf = 0
   do jbf=1,nbf
     if( eigval(jbf) > completely_empty ) then
       ibf = ibf + 1
       p_matrix_sqrt_occ(:,ibf) = p_matrix_sqrt(:,jbf) * SQRT( eigval(jbf) )
     endif
   enddo

   do istate=1,nocc
     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + p_matrix_sqrt_occ(jbf,istate) * eri_3center_lr(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + p_matrix_sqrt_occ(ibf,istate) * eri_3center_lr(:,ipair)
     enddo

     pot_exchange(:,:,ispin) = pot_exchange(:,:,ispin) &
                        - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
   enddo

   deallocate(p_matrix_sqrt_occ)

 enddo

 deallocate(tmp)

 call xsum(pot_exchange)

 call dump_out_matrix(print_matrix_,'=== LR Exchange contribution ===',nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange_ri


!=========================================================================
subroutine setup_exchange_longrange(print_matrix_,nbf,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 use m_inputparam,only: nspin,spin_fact
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Long-Range Exchange term'
 call start_clock(timing_exchange)


 pot_exchange(:,:,:)=0.0_dp

 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           pot_exchange(ibf,jbf,ispin) = pot_exchange(ibf,jbf,ispin) &
                      - eri_lr(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo

 title='=== Exchange contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange


!=========================================================================
subroutine read_potential(print_matrix_,nbf,nspin,p_matrix,pot_read,eread)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_read(nbf,nbf,nspin)
 real(dp),intent(out) :: eread
!=====
 integer              :: ibf,jbf,ispin
 character(len=100)   :: title
 logical              :: file_exists
 integer              :: potfile
!=====

 pot_read(:,:,:)=0.0_dp

 inquire(file='manual_potential',exist=file_exists)
 if(file_exists) then
   open(newunit=potfile,file='manual_potential',status='old')
   do ispin=1,nspin
     do jbf=1,nbf
       do ibf=1,nbf
         read(potfile,*) pot_read(ibf,jbf,ispin)
       enddo
     enddo
   enddo
   close(potfile)

 else
   call die('file not found: manual_potential')
 endif


 title='=== Read potential contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,nspin,pot_read)

 eread = 0.5_dp*SUM(pot_read(:,:,:)*p_matrix(:,:,:))

end subroutine read_potential


!=========================================================================
subroutine setup_density_matrix(nbf,nspin,c_matrix,occupation,p_matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(in)  :: occupation(nbf,nspin)
 real(dp),intent(out) :: p_matrix(nbf,nbf,nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       p_matrix(ibf,jbf,ispin) = SUM( occupation(:,ispin) * c_matrix(ibf,:,ispin) * c_matrix(jbf,:,ispin) )
     enddo
   enddo
 enddo


end subroutine setup_density_matrix


!=========================================================================
subroutine test_density_matrix(nbf,nspin,p_matrix,s_matrix)
 use m_definitions
 use m_warning
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin),s_matrix(nbf,nbf)
!=====
 integer              :: ispin
 real(dp)             :: matrix(nbf,nbf)
 character(len=100)   :: title
!=====

 write(stdout,*) 'Check equality PSP = P'
 write(stdout,*) ' valid only for integer occupation numbers'
 do ispin=1,nspin

   !
   ! Calculate PSP
   matrix(:,:) = MATMUL( p_matrix(:,:,ispin), MATMUL( s_matrix(:,:) , p_matrix(:,:,ispin) ) )


   title='=== PSP ==='
   call dump_out_matrix(1,title,nbf,1,matrix)
   title='===  P  ==='
   call dump_out_matrix(1,title,nbf,1,p_matrix(:,:,ispin))

 enddo


end subroutine test_density_matrix


!=========================================================================
subroutine set_occupation(nbf,temperature,electrons,magnetization,energy,occupation)
 use m_definitions
 use m_mpi
 use m_warning
 use m_inputparam,only: nspin,spin_fact,print_matrix_
 implicit none
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: electrons,magnetization,temperature
 real(dp),intent(in)  :: energy(nbf,nspin)
 real(dp),intent(out) :: occupation(nbf,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin)
 integer              :: ibf,nlines,ilines
 logical              :: file_exists
 integer              :: occfile
!=====

 if( temperature < 1.0e-8_dp ) then

   occupation(:,:)=0.0_dp
  
   inquire(file='manual_occupations',exist=file_exists)
  
   if(.NOT. file_exists) then
     remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
     if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)
  
     do ibf=1,nbf
       occupation(ibf,:) = MIN(remaining_electrons(:), spin_fact)
       remaining_electrons(:)  = remaining_electrons(:) - occupation(ibf,:)
     end do
   else
     write(stdout,*)
     write(stdout,*) 'occupations are read from file: manual_occupations'
     msg='reading manual occupations from file'
     call issue_warning(msg)
     open(newunit=occfile,file='manual_occupations',status='old')
     !
     ! read nlines, all other occupations are set to zero
     read(occfile,*) nlines
     do ilines=1,nlines
       read(occfile,*) occupation(ilines,:)
     enddo
     close(occfile)
     write(stdout,*) 'occupations set, closing file'
   endif

 else
   !
   ! Finite temperature case
   !
   call die('Finite temperature not implemented yet')

 endif

 !
 ! final check
 if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0e-7_dp ) then
   write(stdout,*) 'occupation set up failed to give the right number of electrons'
   write(stdout,*) 'sum of occupations',SUM(occupation(:,:))
   write(stdout,*) 'electrons',electrons
   do ibf=1,nbf
     write(stdout,*) ibf,occupation(ibf,:)
   enddo
   call die('FAILURE in set_occupation')
 endif 

 call dump_out_occupation('=== Occupations ===',nbf,nspin,occupation)

end subroutine set_occupation


!=========================================================================
subroutine matrix_basis_to_eigen(nspin,nbf,c_matrix,matrix_inout)
 use m_definitions
 implicit none
 integer,intent(in)      :: nspin,nbf
 real(dp),intent(in)     :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(inout)  :: matrix_inout(nbf,nbf,nspin)
!====
 integer                 :: ispin
!====


 do ispin=1,nspin
   matrix_inout(:,:,ispin) = MATMUL( TRANSPOSE( c_matrix(:,:,ispin) ) , MATMUL( matrix_inout(:,:,ispin) , c_matrix(:,:,ispin) ) )
 enddo


end subroutine matrix_basis_to_eigen


!=========================================================================
subroutine evaluate_s2_operator(nspin,nbf,occupation,c_matrix,s_matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)      :: nspin,nbf
 real(dp),intent(in)     :: occupation(nbf,nspin)
 real(dp),intent(in)     :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(in)     :: s_matrix(nbf,nbf)
!====
 integer                 :: ispin,istate,jstate
 real(dp)                :: s2,s2_exact
 real(dp)                :: n1,n2,nmax,nmin
!====

 if(nspin /= 2) return

 n1 = SUM( occupation(:,1) )
 n2 = SUM( occupation(:,2) )
 nmax = MAX(n1,n2)
 nmin = MIN(n1,n2)

 s2_exact = (nmax-nmin)/2.0_dp * ( (nmax-nmin)/2.0_dp + 1.0_dp )
 s2       = s2_exact + nmin
 do istate=1,nbf
   if( occupation(istate,1) < completely_empty ) cycle
   do jstate=1,nbf
     if( occupation(jstate,2) < completely_empty ) cycle

     s2 = s2 - ABS( DOT_PRODUCT( c_matrix(:,istate,1) , MATMUL( s_matrix(:,:) , c_matrix(:,jstate,2) ) )  &
                      * occupation(istate,1) * occupation(jstate,2) )**2

   enddo
 enddo


 write(stdout,'(/,a,f8.4)') ' Total Spin S**2 = ',s2
 write(stdout,'(a,f8.4)')   ' Instead of        ',s2_exact


end subroutine evaluate_s2_operator


!=========================================================================
subroutine level_shifting(nbf,s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
 use m_definitions
 use m_warning,only: die
 use m_mpi
 use m_inputparam,only: nspin
 implicit none
 integer,intent(in)     :: nbf
 real(dp),intent(in)    :: s_matrix(nbf,nbf)
 real(dp),intent(in)    :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(in)    :: occupation(nbf,nspin)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(nbf,nbf,nspin)
!====
 integer  :: ispin
 integer  :: ibf
 real(dp) :: sqrt_level_shifting(nbf)
 real(dp) :: matrix_tmp(nbf,nbf)
!====

 write(stdout,'(/,a)')     ' Level shifting switched on'
 write(stdout,'(a,f12.6)') '   energy shift (eV):',level_shifting_energy * Ha_eV

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif


 do ispin=1,nspin
   !
   ! Shift up empty states only
   do ibf=1,nbf
     if( occupation(ibf,ispin) < completely_empty ) then
       sqrt_level_shifting(ibf) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(ibf) = 0.0_dp
     endif
   enddo
   forall(ibf=1:nbf)
     matrix_tmp(:,ibf) =  c_matrix(:,ibf,ispin) * sqrt_level_shifting(ibf)
   end forall

   ! 
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp , TRANSPOSE(matrix_tmp) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + matrix_tmp(:,:)

 enddo
 

end subroutine level_shifting


!=========================================================================
subroutine diagonalize_hamiltonian(nspin_local,nbf,nstate,hamiltonian,s_matrix_sqrt_inv,energy,c_matrix)
 use m_definitions
 use m_timing
 use m_tools
 implicit none

 integer,intent(in)   :: nspin_local,nbf,nstate
 real(dp),intent(in)  :: hamiltonian(nbf,nbf,nspin_local)
 real(dp),intent(in)  :: s_matrix_sqrt_inv(nbf,nstate)
 real(dp),intent(out) :: c_matrix(nbf,nbf,nspin_local)
 real(dp),intent(out) :: energy(nbf,nspin_local)
!=====
 integer  :: ispin,ibf,jbf,istate
 real(dp) :: h_small(nstate,nstate)
!=====

 energy(:,:) = 1.0e+10_dp
 c_matrix(:,:,:) = 0.0_dp
 do ibf=1,nbf
   c_matrix(ibf,ibf,:) = 1.0_dp
 enddo

 do ispin=1,nspin_local
   write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
   call start_clock(timing_diago_hamiltonian)

!   !
!   ! First symmetrize the hamiltonian from the lower triangle to the full matrix
!   do jbf=1,nbf
!     do ibf=1,jbf-1
!       hamiltonian(jbf,ibf,ispin) = hamiltonian(ibf,jbf,ispin)
!     enddo
!   enddo


   h_small(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
                            MATMUL( hamiltonian(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )


   call diagonalize(nstate,h_small,energy(1:nstate,ispin))

   c_matrix(:,1:nstate,ispin) = MATMUL( s_matrix_sqrt_inv(:,:) , h_small(:,:) )


   call stop_clock(timing_diago_hamiltonian)
 enddo

end subroutine diagonalize_hamiltonian


!=========================================================================
