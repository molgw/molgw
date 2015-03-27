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
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         do iatom=1,natom
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


 title='===  Nucleus potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus


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
 real(dp),intent(out) :: pot_hartree(nbf,nbf,nspin)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Hartree term'
 call start_clock(timing_hartree)

 pot_hartree(:,:,:)=0.0_dp

 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
   do jbf=1,nbf
     do ibf=1,nbf
       if( negligible_basispair(ibf,jbf) ) cycle
       if( .NOT. is_my_task(ibf,jbf) ) cycle
       do lbf=1,nbf
         !
         ! symmetry k <-> l
         do kbf=1,lbf-1 ! nbf
           if( negligible_basispair(kbf,lbf) ) cycle
           !
           ! symmetry (ij|kl) = (kl|ij) has been used to loop in the fast order
           pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                      + eri(kbf,lbf,ibf,jbf) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
         enddo
         pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                    + eri(lbf,lbf,ibf,jbf) * SUM( p_matrix(lbf,lbf,:) )
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo
 !
 ! Sum up the different contribution from different procs only if needed
 if( parallel_integral) call xsum(pot_hartree)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,nspin,pot_hartree)

 ehartree = 0.5_dp*SUM(pot_hartree(:,:,:)*p_matrix(:,:,:))

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
 real(dp),intent(out) :: pot_hartree(nbf,nbf,nspin)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 integer              :: nbf_auxil,ibf_auxil
 integer              :: index_ij,index_kl
 real(dp),allocatable :: partial_sum(:)
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity'
 call start_clock(timing_hartree)


 nbf_auxil = SIZE( eri_3center(:,:), DIM=1 )

 allocate(partial_sum(nbf_auxil))
 partial_sum(:) = 0.0_dp
 do lbf=1,nbf
   do kbf=1,nbf
     if( negligible_basispair(kbf,lbf) ) cycle
     index_kl = index_prod(kbf,lbf)
     partial_sum(:) = partial_sum(:) + eri_3center(:,index_kl) * SUM( p_matrix(kbf,lbf,:) )
   enddo
 enddo

 ! Hartree potential is not sensitive to spin
 pot_hartree(:,:,:)=0.0_dp
 do jbf=1,nbf
   do ibf=1,nbf
     if( negligible_basispair(ibf,jbf) ) cycle
     if( .NOT. is_my_task(ibf,jbf) ) cycle
     index_ij = index_prod(ibf,jbf)
     pot_hartree(ibf,jbf,1) = DOT_PRODUCT( eri_3center(:,index_ij) , partial_sum(:) )
   enddo
 enddo
 if( nspin==2) pot_hartree(:,:,nspin) = pot_hartree(:,:,1)

 deallocate(partial_sum)

 !
 ! Sum up the different contribution from different procs only if needed
 if( parallel_integral) call xsum(pot_hartree)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,nspin,pot_hartree)

 ehartree = 0.5_dp*SUM(pot_hartree(:,:,:)*p_matrix(:,:,:))

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
subroutine setup_exchange_ri(print_matrix_,nbf,c_matrix,occupation,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 use m_inputparam,only: nspin,spin_fact
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: c_matrix(nbf,nbf,nspin),occupation(nbf,nspin)
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nbf_auxil,nocc
 real(dp)             :: occ_sqrt_istate
 real(dp),allocatable :: tmp(:,:)
!=====

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)

 nbf_auxil = SIZE( eri_3center(:,:), DIM=1 )

 pot_exchange(:,:,:)=0.0_dp

 allocate(tmp(nbf_auxil,nbf))

 do ispin=1,nspin
   do istate=1,nbf
     if( occupation(istate,ispin) > completely_empty ) nocc = istate
   enddo

   do istate=1,nocc
     occ_sqrt_istate = SQRT( occupation(istate,ispin) )
     tmp(:,:) = 0.0_dp
     do jbf=1,nbf
       do ibf=1,nbf
         if( negligible_basispair(ibf,jbf) ) cycle
         index_ij = index_prod(ibf,jbf)
         tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center(:,index_ij) * occ_sqrt_istate
       enddo
     enddo

     pot_exchange(:,:,ispin) = pot_exchange(:,:,ispin) &
                        - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
   enddo

 enddo

 deallocate(tmp)



 call dump_out_matrix(print_matrix_,'=== Exchange contribution ===',nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri


!=========================================================================
subroutine setup_exchange_longrange_ri(print_matrix_,nbf,c_matrix,occupation,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 use m_inputparam,only: nspin,spin_fact
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: c_matrix(nbf,nbf,nspin),occupation(nbf,nspin)
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nbf_auxil,nocc
 real(dp)             :: occ_sqrt_istate
 real(dp),allocatable :: tmp(:,:)
!=====

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)

 nbf_auxil = SIZE( eri_3center_lr(:,:), DIM=1 )

 pot_exchange(:,:,:)=0.0_dp

 allocate(tmp(nbf_auxil,nbf))

 do ispin=1,nspin
   do istate=1,nbf
     if( occupation(istate,ispin) > completely_empty ) nocc = istate
   enddo

   do istate=1,nocc
     occ_sqrt_istate = SQRT( occupation(istate,ispin) )
     tmp(:,:) = 0.0_dp
     do jbf=1,nbf
       do ibf=1,nbf
         if( negligible_basispair(ibf,jbf) ) cycle
         index_ij = index_prod(ibf,jbf)
         tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center_lr(:,index_ij) * occ_sqrt_istate
       enddo
     enddo

     pot_exchange(:,:,ispin) = pot_exchange(:,:,ispin) &
                        - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
   enddo

 enddo

 deallocate(tmp)



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
   stop'file not found: manual_potential'
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
subroutine set_occupation(electrons,magnetization,nbf,occupation)
 use m_definitions
 use m_mpi
 use m_warning
 use m_inputparam,only: nspin,spin_fact,print_matrix_
 implicit none
 real(dp),intent(in)  :: electrons,magnetization
 integer,intent(in)   :: nbf
 real(dp),intent(out) :: occupation(nbf,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin)
 integer              :: ibf,nlines,ilines
 logical              :: file_exists
 integer              :: occfile
!=====


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
 
  !
  ! final check
  if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0d-7 ) then
    write(stdout,*) 'occupation set up failed to give the right number of electrons'
    write(stdout,*) 'sum of occupations',SUM(occupation(:,:))
    write(stdout,*) 'electrons',electrons
    do ibf=1,nbf
      write(stdout,*) ibf,occupation(ibf,:)
    enddo
    stop'FAILURE in set_occupation'
  endif 

 call dump_out_occupation('=== Occupations ===',nbf,nspin,occupation)

end subroutine set_occupation


!=========================================================================
subroutine guess_starting_c_matrix(nbf,nspin,c_matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(out) :: c_matrix(nbf,nbf,nspin)
!=====
 integer :: ibf
!=====

 !
 ! fill the c_matrix with the identity
 c_matrix(:,:,:)=0.0_dp
 do ibf=1,nbf
   c_matrix(ibf,ibf,:) = 1.0_dp
!   c_matrix(ibf,modulo(ibf,nbf)+1,:) = 1.0_dp
 enddo

end subroutine guess_starting_c_matrix


!=========================================================================
subroutine guess_starting_c_matrix_new(basis,c_matrix)
 use m_definitions
 use m_mpi
 use m_gaussian
 use m_basis_set
 use m_inputparam,only: nspin
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: c_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer  :: ibf,jbf,kbf,ig
 real(dp) :: alpha_max_bf(basis%nbf),alpha_max_remaining
!=====

 !
 ! find the sharpest gaussians
 alpha_max_bf(:)=0.0_dp
 do ibf=1,basis%nbf
   do ig=1,basis%bf(ibf)%ngaussian
     alpha_max_bf(ibf)=MAX(basis%bf(ibf)%g(ig)%alpha,alpha_max_bf(ibf))
   enddo
!   write(stdout,*) ibf,alpha_max_bf(ibf)
 enddo

 !
 ! fill the c_matrix 
 c_matrix(:,:,:)=0.0_dp
 do ibf=1,basis%nbf
   alpha_max_remaining=0.0_dp
   do jbf=1,basis%nbf
     if( alpha_max_bf(jbf) > alpha_max_remaining ) then
       alpha_max_remaining = alpha_max_bf(jbf)
       kbf = jbf
     endif
   enddo
   c_matrix(kbf,ibf,:) = 1.0_dp
!   write(stdout,*) 'chosen',ibf,kbf,alpha_max_bf(kbf)
   alpha_max_bf(kbf)   = -1.0_dp
   
 enddo

end subroutine guess_starting_c_matrix_new


!=========================================================================
subroutine setup_initial_c_matrix(print_matrix_,nbf,nspin,hamiltonian_nucleus,s_matrix,occupation,c_matrix)
 use m_definitions
 use m_mpi
 use m_tools
 use m_timing
 implicit none
 integer,intent(in)         :: print_matrix_,nspin,nbf
 real(dp),intent(in)        :: hamiltonian_nucleus(nbf,nbf),s_matrix(nbf,nbf)
 real(dp),intent(in)        :: occupation(nbf,nspin)
 real(dp),intent(out)       :: c_matrix(nbf,nbf,nspin)
!=====
 integer                    :: ibf,jbf,kbf
 real(dp)                   :: hamiltonian(nbf,nbf),matrix(nbf,nbf),energy(nbf)
 real(dp)                   :: coeff_max,bonding
 real(dp)                   :: line(nbf)
 character(len=100)         :: title
!=====


 !
 ! Diagonalize a spin independant hamiltonian
 ! to obtain a starting point for matrix C
 hamiltonian(:,:) = hamiltonian_nucleus(:,:)

 title='=== bare hamiltonian ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,hamiltonian)

 write(stdout,*) 'Diagonalization of an initial hamiltonian'
 call start_clock(timing_diago_hamiltonian)
 call diagonalize_generalized_sym(nbf,hamiltonian,s_matrix,energy,matrix)
 call stop_clock(timing_diago_hamiltonian)

 title='=== Energies ==='
 call dump_out_eigenenergy(title,nbf,1,occupation,energy)


 ibf=1
 do while( ibf < nbf )
!   if( ALL( occupation(ibf,:) < completely_empty ) cycle

   jbf = ibf + 1
   !
   ! Find degenerate energies
   if( ABS( energy(ibf) - energy(jbf) ) < 1.0e-3_dp ) then
     !
     ! Find the bonding and anti bonding states
     coeff_max = MAXVAL( ABS( matrix(:,ibf) ) )
     bonding = 1.0_dp
     do kbf=1,nbf
       if( ABS( matrix(kbf,ibf) ) > coeff_max - 1.0e-3_dp ) then
         bonding = SIGN( 1.0_dp, matrix(kbf,ibf) ) * bonding
       endif
     enddo
     if( bonding < 0.0_dp ) then
       line(:)       = matrix(:,ibf)
       matrix(:,ibf) = matrix(:,jbf)
       matrix(:,jbf) = line(:)
     endif

     ibf = ibf + 2
   else
     ibf = ibf + 1
   endif

 enddo

 c_matrix(:,:,1)     = matrix(:,:)
 c_matrix(:,:,nspin) = matrix(:,:)

 matrix(:,:) = TRANSPOSE( matrix(:,:) )
 title='=== C matrix ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,matrix)


end subroutine setup_initial_c_matrix


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
