!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min


contains


!=========================================================================
subroutine setup_overlap(print_matrix_,basis,s_matrix)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(basis%nbf,basis%nbf)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
!=====

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup overlap matrix S'

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

 call stop_clock(timing_overlap)


end subroutine setup_overlap


!=========================================================================
subroutine setup_overlap_mixedbasis(print_matrix_,basis1,basis2,s_matrix)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis1,basis2
 real(dp),intent(out)       :: s_matrix(basis1%nbf,basis2%nbf)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
!=====

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis1%nbf_cart)
   li      = basis1%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis1%gaussian_type,li)

   do while(jbf_cart<=basis2%nbf_cart)
     lj      = basis2%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis2%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call overlap_basis_function(basis1%bf(ibf_cart+i_cart-1),basis2%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
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


end subroutine setup_overlap_mixedbasis


!=========================================================================
subroutine setup_kinetic(print_matrix_,basis,hamiltonian_kinetic)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(basis%nbf,basis%nbf)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
!=====

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
 use m_basis_set
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
!=====

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
subroutine setup_hartree(print_matrix_,nbf,p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: hartree_ij(nbf,nbf)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Hartree term'
 call start_clock(timing_hartree)

 hartree_ij(:,:)=0.0_dp

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
         hartree_ij(ibf,jbf) = hartree_ij(ibf,jbf) &
                    + eri(kbf,lbf,ibf,jbf) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
       enddo
       hartree_ij(ibf,jbf) = hartree_ij(ibf,jbf) &
                  + eri(lbf,lbf,ibf,jbf) * SUM( p_matrix(lbf,lbf,:) )
     enddo
   enddo
 enddo


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,hartree_ij)

 ehartree = 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree


!=========================================================================
subroutine setup_hartree_ri(print_matrix_,nbf,p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: hartree_ij(nbf,nbf)
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
 hartree_ij(:,:) = 0.0_dp
 do ipair=1,npair
   ibf = index_basis(1,ipair)
   jbf = index_basis(2,ipair)
   rtmp = DOT_PRODUCT( eri_3center(:,ipair) , partial_sum(:) )
   hartree_ij(ibf,jbf) = rtmp
   hartree_ij(jbf,ibf) = rtmp
 enddo

 deallocate(partial_sum)

 !
 ! Sum up the different contribution from different procs only if needed
 call xsum(hartree_ij)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,hartree_ij)

 ehartree = 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree_ri


!=========================================================================
subroutine setup_exchange(print_matrix_,nbf,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Exchange term'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 do ispin=1,nspin
   do jbf=1,nbf
     do lbf=1,nbf
       if( negligible_basispair(lbf,jbf) ) cycle
       do kbf=1,nbf
!         if( ABS(p_matrix(kbf,lbf,ispin)) <  1.0e-12_dp ) cycle 
         do ibf=1,nbf
           if( negligible_basispair(ibf,kbf) ) cycle
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           exchange_ij(ibf,jbf,ispin) = exchange_ij(ibf,jbf,ispin) &
                      - eri(ibf,kbf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
 enddo


 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange


!=========================================================================
subroutine setup_exchange_ri(print_matrix_,nbf,p_matrix_occ,p_matrix_sqrt,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix_occ(nbf,nspin)
 real(dp),intent(in)  :: p_matrix_sqrt(nbf,nbf,nspin)
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
!=====

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)

 exchange_ij(:,:,:)=0.0_dp

 allocate(tmp(nauxil_3center,nbf))

 do ispin=1,nspin

   do istate=1,nbf
     if( p_matrix_occ(istate,ispin) < completely_empty)  cycle

     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + p_matrix_sqrt(jbf,istate,ispin) * eri_3center(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + p_matrix_sqrt(ibf,istate,ispin) * eri_3center(:,ipair)
     enddo

     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nauxil_3center,-1.0_dp/spin_fact,tmp,nauxil_3center,1.0_dp,exchange_ij(:,:,ispin),nbf)
   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo

 enddo
 deallocate(tmp)

 call xsum(exchange_ij)

 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri


!=========================================================================
subroutine setup_exchange_longrange_ri(print_matrix_,nbf,p_matrix_occ,p_matrix_sqrt,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix_occ(nbf,nspin)
 real(dp),intent(in)  :: p_matrix_sqrt(nbf,nbf,nspin)
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
!=====

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 allocate(tmp(nauxil_3center_lr,nbf))

 do ispin=1,nspin

   do istate=1,nbf
     if( p_matrix_occ(istate,ispin) < completely_empty)  cycle

     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + p_matrix_sqrt(jbf,istate,ispin) * eri_3center_lr(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + p_matrix_sqrt(ibf,istate,ispin) * eri_3center_lr(:,ipair)
     enddo

     !exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                   - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nauxil_3center_lr,-1.0_dp/spin_fact,tmp,nauxil_3center_lr,1.0_dp,exchange_ij(:,:,ispin),nbf)
   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo

 enddo
 deallocate(tmp)

 call xsum(exchange_ij)

 call dump_out_matrix(print_matrix_,'=== LR Exchange contribution ===',nbf,nspin,exchange_ij)

 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange_ri


!=========================================================================
subroutine setup_exchange_longrange(print_matrix_,nbf,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
!=====

 write(stdout,*) 'Calculate Long-Range Exchange term'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           exchange_ij(ibf,jbf,ispin) = exchange_ij(ibf,jbf,ispin) &
                      - eri_lr(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
 enddo


 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange


!=========================================================================
subroutine read_potential(print_matrix_,nbf,nspin,p_matrix,pot_read,eread)
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
subroutine setup_density_matrix(nbf,nstate,c_matrix,occupation,p_matrix)
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
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
subroutine set_occupation(nstate,temperature,electrons,magnetization,energy,occupation)
 use m_inputparam,only: print_matrix_
 implicit none
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: electrons,magnetization,temperature
 real(dp),intent(in)  :: energy(nstate,nspin)
 real(dp),intent(out) :: occupation(nstate,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin)
 integer              :: istate,nlines,ilines
 logical              :: file_exists
 integer              :: occfile
!=====

 if( temperature < 1.0e-8_dp ) then

   occupation(:,:)=0.0_dp
  
   inquire(file='manual_occupations',exist=file_exists)
  
   if(.NOT. file_exists) then
     remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
     if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)
  
     do istate=1,nstate
       occupation(istate,:) = MIN(remaining_electrons(:), spin_fact)
       remaining_electrons(:)  = remaining_electrons(:) - occupation(istate,:)
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
   do istate=1,nstate
     write(stdout,*) istate,occupation(istate,:)
   enddo
   call die('FAILURE in set_occupation')
 endif 

 call dump_out_occupation('=== Occupations ===',nstate,nspin,occupation)

end subroutine set_occupation


!=========================================================================
subroutine matrix_basis_to_eigen(nbf,nstate,c_matrix,matrix_inout)
 implicit none
 integer,intent(in)      :: nbf,nstate
 real(dp),intent(in)     :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(inout)  :: matrix_inout(nbf,nbf,nspin)
!=====
 integer                 :: ispin
!=====


 do ispin=1,nspin
   matrix_inout(1:nstate,1:nstate,ispin) = MATMUL( TRANSPOSE( c_matrix(:,:,ispin) ) , MATMUL( matrix_inout(:,:,ispin) , c_matrix(:,:,ispin) ) )
 enddo


end subroutine matrix_basis_to_eigen


!=========================================================================
subroutine evaluate_s2_operator(nbf,nstate,occupation,c_matrix,s_matrix)
 implicit none
 integer,intent(in)      :: nbf,nstate
 real(dp),intent(in)     :: occupation(nstate,nspin)
 real(dp),intent(in)     :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)     :: s_matrix(nbf,nbf)
!=====
 integer                 :: ispin,istate,jstate
 real(dp)                :: s2,s2_exact
 real(dp)                :: n1,n2,nmax,nmin
!=====

 if(nspin /= 2) return

 n1 = SUM( occupation(:,1) )
 n2 = SUM( occupation(:,2) )
 nmax = MAX(n1,n2)
 nmin = MIN(n1,n2)

 s2_exact = (nmax-nmin)/2.0_dp * ( (nmax-nmin)/2.0_dp + 1.0_dp )
 s2       = s2_exact + nmin
 do istate=1,nstate
   if( occupation(istate,1) < completely_empty ) cycle
   do jstate=1,nstate
     if( occupation(jstate,2) < completely_empty ) cycle

     s2 = s2 - ABS( DOT_PRODUCT( c_matrix(:,istate,1) , MATMUL( s_matrix(:,:) , c_matrix(:,jstate,2) ) )  &
                      * occupation(istate,1) * occupation(jstate,2) )**2

   enddo
 enddo


 write(stdout,'(/,a,f8.4)') ' Total Spin S**2: ',s2
 write(stdout,'(a,f8.4)')   ' Instead of:      ',s2_exact


end subroutine evaluate_s2_operator


!=========================================================================
subroutine level_shifting(nbf,nstate,s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
 implicit none
 integer,intent(in)     :: nbf,nstate
 real(dp),intent(in)    :: s_matrix(nbf,nbf)
 real(dp),intent(in)    :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)    :: occupation(nstate,nspin)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(nbf,nbf,nspin)
!=====
 integer  :: ispin
 integer  :: ibf,istate
 real(dp) :: sqrt_level_shifting(nstate)
 real(dp) :: matrix_tmp(nbf,nbf)
!=====

 write(stdout,'(/,a)')     ' Level shifting switched on'
 write(stdout,'(a,f12.6)') '   energy shift (eV):',level_shifting_energy * Ha_eV

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif


 do ispin=1,nspin
   !
   ! Shift up empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   ! 
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + matrix_tmp(:,:)

 enddo
 

end subroutine level_shifting


!=========================================================================
subroutine setup_sqrt_overlap(TOL_OVERLAP,nbf,s_matrix,nstate,s_matrix_sqrt_inv)
 use m_tools
 implicit none

 real(dp),intent(in)                :: TOL_OVERLAP
 integer,intent(in)                 :: nbf
 real(dp),intent(in)                :: s_matrix(nbf,nbf)
 integer,intent(out)                :: nstate
 real(dp),allocatable,intent(inout) :: s_matrix_sqrt_inv(:,:)
!=====
 integer  :: ibf,jbf
 real(dp) :: s_eigval(nbf)
 real(dp) :: matrix_tmp(nbf,nbf)
!=====

 write(stdout,'(/,a)') ' Calculate overlap matrix square-root S^{1/2}'

 matrix_tmp(:,:) = s_matrix(:,:)
 ! Diagonalization with or without SCALAPACK
 call diagonalize_scalapack(scalapack_block_min,nbf,matrix_tmp,s_eigval)

 nstate = COUNT( s_eigval(:) > TOL_OVERLAP )

 allocate(s_matrix_sqrt_inv(nbf,nstate))

 write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
 write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
 write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',TOL_OVERLAP
 write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',nbf

 ibf=0
 do jbf=1,nbf
   if( s_eigval(jbf) > TOL_OVERLAP ) then
     ibf = ibf + 1
     s_matrix_sqrt_inv(:,ibf) = matrix_tmp(:,jbf) / SQRT( s_eigval(jbf) )
   endif
 enddo


end subroutine setup_sqrt_overlap


!=========================================================================
subroutine setup_sqrt_density_matrix(nbf,p_matrix,p_matrix_sqrt,p_matrix_occ)
 use m_tools
 implicit none

 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: p_matrix_sqrt(nbf,nbf,nspin)
 real(dp),intent(out) :: p_matrix_occ(nbf,nspin)
!=====
 integer              :: ispin,ibf
!=====

 write(stdout,*) 'Calculate the square root of the density matrix'
 call start_clock(timing_sqrt_density_matrix)

 do ispin=1,nspin
   p_matrix_sqrt(:,:,ispin) = p_matrix(:,:,ispin)
   ! Diagonalization with or without SCALAPACK
   call diagonalize_scalapack(scalapack_block_min,nbf,p_matrix_sqrt(:,:,ispin),p_matrix_occ(:,ispin))
   do ibf=1,nbf
     ! this is to avoid instabilities
     if( p_matrix_occ(ibf,ispin) < 1.0e-8_dp ) then
       p_matrix_occ(ibf,ispin)    = 0.0_dp
       p_matrix_sqrt(:,ibf,ispin) = 0.0_dp
     else
       p_matrix_sqrt(:,ibf,ispin) = p_matrix_sqrt(:,ibf,ispin) * SQRT( p_matrix_occ(ibf,ispin) )
     endif
   enddo
 enddo

 call stop_clock(timing_sqrt_density_matrix)

end subroutine setup_sqrt_density_matrix


!=========================================================================
subroutine dft_approximate_vhxc(basis,vhxc_ij)
 use m_basis_set
 use m_dft_grid
 use m_eri_calculate
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: vhxc_ij(basis%nbf,basis%nbf)
!=====
 integer              :: idft_xc
 integer              :: igrid,ibf,jbf,ispin
 real(dp)             :: rr(3)
 real(dp)             :: normalization
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: rhor
 real(dp)             :: vxc,excr,exc
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhartree
 real(dp)             :: vhgau(basis%nbf,basis%nbf)
 integer              :: iatom,igau,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
!=====

 call start_clock(timing_approx_ham)

 vhxc_ij(:,:) = 0.0_dp

 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities'

 do iatom=1,natom
   if( rank /= MODULO(iatom,nproc) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),coeff,alpha)


   call calculate_eri_approximate_hartree(basis,basis%nbf,basis%nbf,x(:,iatom),ngau,coeff,alpha,vhgau)
   vhxc_ij(:,:) = vhxc_ij(:,:) + vhgau(:,:) 

   deallocate(alpha,coeff)
 enddo

 write(stdout,*) 'Simple LDA functional on a coarse grid'

 !
 ! Create a temporary grid with low quality
 ! This grid is to be destroyed at the end of the present subroutine
 call init_dft_grid(low)

 !
 ! If it is the first time, set up the stored arrays
 !
 if( .NOT. ALLOCATED(bfr) ) call prepare_basis_functions_r(basis)

 normalization = 0.0_dp
 exc           = 0.0_dp
 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)
   weight = w_grid(igrid)

   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r(basis,igrid,basis_function_r)

   !
   ! calculate the density at point r for spin up and spin down
   call setup_atomic_density(rr,rhor,vhartree)

   !
   ! Normalization
   normalization = normalization + rhor * weight

   call teter_lda_vxc_exc(rhor,vxc,excr)

   !
   ! XC energy
   exc = exc + excr * weight * rhor

   !
   ! HXC
   do jbf=1,basis%nbf
     do ibf=1,basis%nbf 
       vhxc_ij(ibf,jbf) =  vhxc_ij(ibf,jbf) + weight &
           *  vxc * basis_function_r(ibf) * basis_function_r(jbf)
     enddo
   enddo

 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum(normalization)
 call xsum(exc)
 call xsum(vhxc_ij)

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization
 write(stdout,  '(a,2(2x,f12.6))') '      XC energy (Ha):',exc

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()

 call stop_clock(timing_approx_ham)

end subroutine dft_approximate_vhxc


!=========================================================================
subroutine virtual_fno(basis,nstate,occupation,energy,c_matrix,energy_ref,c_matrix_ref)
 use m_memory
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(inout)                :: energy(nstate,nspin)
 real(dp),intent(inout)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),allocatable,intent(inout)    :: energy_ref(:,:)
 real(dp),allocatable,intent(inout)    :: c_matrix_ref(:,:,:)
!=====
 real(dp),parameter                    :: alpha_mp2=1.0_dp
 integer                               :: istate,jstate,astate,bstate,cstate
 integer                               :: ispin
 integer                               :: nocc,ncore,nvirtual,nvirtual_kept
 real(dp),allocatable                  :: p_matrix_mp2(:,:),ham_virtual(:,:),ham_virtual_kept(:,:)
 real(dp),allocatable                  :: occupation_mp2(:),energy_virtual_kept(:)
 real(dp)                              :: eri_ci_aj,eri_ci_bj
 real(dp)                              :: den_ca_ij,den_cb_ij
 real(dp)                              :: en_mp2
 logical                               :: file_found
 integer                               :: file_unit
 integer                               :: nvirtualmin
 real(dp),allocatable                  :: eri_ci_i(:)
!=====

 call start_clock(timing_fno)

 write(stdout,'(/,x,a)') 'Prepare optimized empty states with Frozen Natural Orbitals'

 nvirtualmin = MIN(nvirtualw,nvirtualg)

 if( nvirtualmin > nstate ) then 
   call issue_warning('virtual_fno is on, however nvirtualw and nvirtualg are not set. Skipping the Frozen Natural Orbitals generation.')
   return
 endif

 if( .NOT. has_auxil_basis ) then
   call issue_warning('virtual_fno not implemented when no auxiliary basis is provided')
   return
 endif

 if( nspin > 1 ) then
   call issue_warning('virtual_fno not implemented yet for spin polarized calculations')
   return
 endif

 ncore = ncoreg
 if(is_frozencore) then
   if( ncore == 0) ncore = atoms_core_states()
 endif


 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix)


 do ispin=1,nspin
   !
   ! First, set up the list of occupied states
   nocc = ncore
   do istate=ncore+1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     nocc = istate
   enddo

!   call issue_warning('FBFB cheating here')
!   nocc = nocc + 1
 


   nvirtual = nstate - nocc

   allocate(p_matrix_mp2(nvirtual,nvirtual))
   p_matrix_mp2(:,:) = 0.0_dp

   allocate(eri_ci_i(nocc+1:nstate))

#if 1
   ! Approximation by Aquilante et al.
   do istate=ncore+1,nocc
     do cstate=nocc+1,nstate

       do astate=nocc+1,nstate
         eri_ci_i(astate) = eri_eigen_ri(cstate,istate,ispin,astate,istate,ispin) &
                              / ( energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(cstate,ispin) )
       enddo
       call xsum(eri_ci_i)

       do bstate=nocc+1,nstate
         do astate=nocc+1,nstate

           p_matrix_mp2(astate-nocc,bstate-nocc) = &
                p_matrix_mp2(astate-nocc,bstate-nocc)  &
                   + 0.50_dp * eri_ci_i(astate) * eri_ci_i(bstate)

         enddo
       enddo

     enddo
   enddo
   deallocate(eri_ci_i)
#else

   do bstate=nocc+1,nstate
     do astate=nocc+1,nstate

#if 0
       ! Full calculation of the MP2 density matrix on virtual orbitals (See Taube and Bartlett)
       do cstate=nocc+1,nstate
         do istate=ncore+1,nocc
           do jstate=ncore+1,nocc

             den_cb_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(bstate,ispin) - energy(cstate,ispin)
             den_ca_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(astate,ispin) - energy(cstate,ispin)

             eri_ci_aj = eri_eigen_ri(cstate,istate,ispin,astate,jstate,ispin) &
                          - eri_eigen_ri(cstate,jstate,ispin,astate,istate,ispin)  * alpha_mp2
             eri_ci_bj = eri_eigen_ri(cstate,istate,ispin,bstate,jstate,ispin) &
                         - eri_eigen_ri(cstate,jstate,ispin,bstate,istate,ispin)   * alpha_mp2

             p_matrix_mp2(astate-nocc,bstate-nocc) = &
                  p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                     + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

           enddo
         enddo
       enddo
#elif 1
       ! Approximation by Aquilante et al.
       do cstate=nocc+1,nstate
         do istate=ncore+1,nocc
           den_cb_ij = energy(istate,ispin) + energy(istate,ispin) - energy(bstate,ispin) - energy(cstate,ispin)
           den_ca_ij = energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(cstate,ispin)

           eri_ci_aj = eri_eigen_ri(cstate,istate,ispin,astate,istate,ispin) 
           eri_ci_bj = eri_eigen_ri(cstate,istate,ispin,bstate,istate,ispin) 

           p_matrix_mp2(astate-nocc,bstate-nocc) = &
                p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                   + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

         enddo
       enddo
#elif 0
       do istate=ncore+1,nocc
         ! Approximation by Aquilante et al.
         den_cb_ij = energy(istate,ispin) + energy(istate,ispin) - energy(bstate,ispin) - energy(bstate,ispin)
         den_ca_ij = energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(astate,ispin)

         eri_ci_aj = eri_eigen_ri(astate,istate,ispin,astate,istate,ispin) 
         eri_ci_bj = eri_eigen_ri(bstate,istate,ispin,bstate,istate,ispin) 

         p_matrix_mp2(astate-nocc,bstate-nocc) = &
              p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                 + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )
       enddo
#else
         do istate=ncore+1,nocc
           do jstate=ncore+1,nocc

             den_cb_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(bstate,ispin) - energy(bstate,ispin)
             den_ca_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(astate,ispin) - energy(astate,ispin)

             eri_ci_aj = eri_eigen_ri(astate,istate,ispin,astate,jstate,ispin) 
             eri_ci_bj = eri_eigen_ri(bstate,istate,ispin,bstate,jstate,ispin) 

             p_matrix_mp2(astate-nocc,bstate-nocc) = &
                  p_matrix_mp2(astate-nocc,bstate-nocc)  &
                     + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

           enddo
         enddo

#endif

     enddo
   enddo


#endif

   allocate(occupation_mp2(nvirtual))
   call diagonalize(nvirtual,p_matrix_mp2,occupation_mp2)

!   write(stdout,*) 
!   do astate=1,nvirtual
!     write(stdout,'(i4,2x,es16.4)') astate,occupation_mp2(astate)
!   enddo
!   write(stdout,*) 

   allocate(ham_virtual(nvirtual,nvirtual))

   ham_virtual(:,:) = 0.0_dp
   do astate=1,nvirtual
     ham_virtual(astate,astate) = energy(astate+nocc,ispin)
   enddo

   nvirtual_kept = nvirtualmin - 1 - nocc

   write(stdout,'(/,x,a,i5)')    'Max state index included ',nvirtual_kept + nocc
   write(stdout,'(x,a,i5,a,i5)') 'Retain ',nvirtual_kept,' virtual orbitals out of ',nvirtual
   write(stdout,'(x,a,es14.6)')  'Occupation number of the last virtual state',occupation_mp2(nvirtual - (nocc + nvirtual_kept))
   write(stdout,'(x,a,es14.6,4x,es14.6)') '  to be compared to the first and last virtual states',occupation_mp2(nvirtual),occupation_mp2(1)
   write(stdout,*)

   allocate(ham_virtual_kept(nvirtual_kept,nvirtual_kept))
   allocate(energy_virtual_kept(nvirtual_kept))

   ham_virtual_kept(:,:)  = MATMUL( TRANSPOSE( p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) ) ,   &
                                        MATMUL( ham_virtual , p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) ) )

   deallocate(ham_virtual)

   call diagonalize(nvirtual_kept,ham_virtual_kept,energy_virtual_kept)

!   write(stdout,'(/,x,a)') ' virtual state    FNO energy (eV)   reference energy (eV)'
!   do astate=1,nvirtual_kept
!     write(stdout,'(i4,4x,f16.5,2x,f16.5)') astate,energy_virtual_kept(astate)*Ha_eV,energy(astate+nocc,ispin)*Ha_eV
!   enddo
!   write(stdout,*)
   

   !
   ! Save the original c_matrix and energies
   if( .NOT. ALLOCATED(energy_ref) ) then
     allocate(energy_ref(nocc+1:nocc+nvirtual_kept,nspin))
     allocate(c_matrix_ref(basis%nbf,nocc+1:nocc+nvirtual_kept,nspin))   ! TODO clean_allocate of this

     energy_ref(nocc+1:nocc+nvirtual_kept,:)     = energy(nocc+1:nocc+nvirtual_kept,:)
     c_matrix_ref(:,nocc+1:nocc+nvirtual_kept,:) = c_matrix(:,nocc+1:nocc+nvirtual_kept,:)
   endif

   !
   ! And then override the c_matrix and the energy with the fresh new ones
   energy(nocc+1:nocc+nvirtual_kept,ispin)     = energy_virtual_kept(:)
   c_matrix(:,nocc+1:nocc+nvirtual_kept,ispin) = MATMUL( c_matrix(:,nocc+1:,ispin) ,  &
                                                    MATMUL( p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) , &
                                                       ham_virtual_kept(:,:) ) )


   deallocate(energy_virtual_kept)
   deallocate(ham_virtual_kept)
   deallocate(p_matrix_mp2,occupation_mp2)
 enddo



 call destroy_eri_3center_eigen()

 write(stdout,'(x,a)') 'Optimized empty states with Frozen Natural Orbitals'

 call stop_clock(timing_fno)

end subroutine virtual_fno


!=========================================================================
subroutine destroy_fno(basis,nstate,energy,c_matrix,energy_ref,c_matrix_ref)
 use m_memory
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(inout)                :: energy(nstate,nspin)
 real(dp),intent(inout)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),allocatable,intent(inout)    :: energy_ref(:,:)
 real(dp),allocatable,intent(inout)    :: c_matrix_ref(:,:,:)
!=====
 integer                               :: lowerb,upperb
!=====

 write(stdout,'(/,x,a)') 'Deallocate the Frozen Natural Orbitals'

 lowerb = LBOUND(energy_ref,DIM=1)
 upperb = UBOUND(energy_ref,DIM=1)


 energy    (lowerb:upperb,:) = energy_ref    (lowerb:upperb,:)
 c_matrix(:,lowerb:upperb,:) = c_matrix_ref(:,lowerb:upperb,:)


 deallocate(energy_ref)
 deallocate(c_matrix_ref)   ! TODO clean_allocate of this


end subroutine destroy_fno




!=========================================================================
end module m_hamiltonian
!=========================================================================
