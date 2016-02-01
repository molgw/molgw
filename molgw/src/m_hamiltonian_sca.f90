!=========================================================================
module m_hamiltonian_sca
 use m_definitions
 use m_mpi
 use m_timing
 use m_warning
 use m_inputparam,only: nspin,spin_fact


contains


!=========================================================================
subroutine matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,matrix_local)
 use m_basis_set
 implicit none

 integer,intent(in)     :: ibf,jbf
 integer,intent(in)     :: li,lj
 integer,intent(in)     :: ni_cart,nj_cart,ni,nj
 integer,intent(in)     :: m_ham,n_ham
 real(dp),intent(in)    :: matrix_cart(ni_cart,nj_cart)
 real(dp),intent(inout) :: matrix_local(m_ham,n_ham)
!=====
 integer  :: iglobal,jglobal
 integer  :: ilocal,jlocal
 real(dp) :: matrix_final(ibf:ibf+ni-1,jbf:jbf+nj-1)
!=====

#ifdef HAVE_SCALAPACK

 matrix_final(:,:) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                             MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) )

 do jglobal=jbf,jbf+nj-1
   jlocal = colindex_global_to_local('H',jglobal)
   if( jlocal == 0 ) cycle

   do iglobal=ibf,ibf+ni-1
     ilocal = rowindex_global_to_local('H',iglobal)
     if( ilocal == 0 ) cycle

     matrix_local(ilocal,jlocal) = matrix_final(iglobal,jglobal)

   enddo
 enddo


#endif

end subroutine matrix_cart_to_local


!=========================================================================
subroutine setup_overlap_sca(print_matrix_,basis,m_ham,n_ham,s_matrix)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: s_matrix(m_ham,n_ham)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
!=====

#ifdef HAVE_SCALAPACK

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup overlap matrix S: SCALAPACK'

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

     call matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,s_matrix)


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


#endif

end subroutine setup_overlap_sca


!=========================================================================
subroutine setup_kinetic_sca(print_matrix_,basis,m_ham,n_ham,hamiltonian_kinetic)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: hamiltonian_kinetic(m_ham,n_ham)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart,iglobal,jglobal,ilocal,jlocal
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp),allocatable :: matrix_final(:,:)
!=====

#ifdef HAVE_SCALAPACK

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian: SCALAPACK'

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
     allocate(matrix_final(ni,nj))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call kinetic_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo

     call matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,hamiltonian_kinetic)

     deallocate(matrix_cart,matrix_final)
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

#endif

end subroutine setup_kinetic_sca


!=========================================================================
subroutine setup_nucleus_sca(print_matrix_,basis,m_ham,n_ham,hamiltonian_nucleus)
 use m_basis_set
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: hamiltonian_nucleus(m_ham,n_ham)
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

#ifdef HAVE_SCALAPACK

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian: SCALAPACK'
 if( nproc_local > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_local /= MODULO(iatom-1,nproc_local) ) cycle
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
       if( rank_local /= MODULO(iatom-1,nproc_local) ) cycle
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           call nucleus_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),zatom(iatom),x(:,iatom),vnucleus_ij)
           matrix_cart(i_cart,j_cart) = matrix_cart(i_cart,j_cart) + vnucleus_ij
         enddo
       enddo
     enddo

     call matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,hamiltonian_nucleus)


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
 call xlocal_sum(hamiltonian_nucleus)

 title='===  Nucleus potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

#endif

end subroutine setup_nucleus_sca


!=========================================================================
subroutine setup_hartree_ri_sca(print_matrix_,nbf,m_ham,n_ham,p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: hartree_ij(m_ham,n_ham)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ilocal,jlocal
 integer              :: iglobal,jglobal
 integer              :: ibf,jbf,kbf,lbf,ispin
 integer              :: ibf_auxil,ipair
 integer              :: index_ij,index_kl
 real(dp),allocatable :: partial_sum(:)
 real(dp)             :: rtmp
 character(len=100)   :: title
!=====

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity: SCALAPACK'
 call start_clock(timing_hartree)

 allocate(partial_sum(nauxil_3center))
 partial_sum(:) = 0.0_dp

 do jlocal=1,n_ham
   jglobal = colindex_local_to_global('H',jlocal)

   do ilocal=1,m_ham
     iglobal = rowindex_local_to_global('H',ilocal)
     if( negligible_basispair(iglobal,jglobal) ) cycle

     partial_sum(:) = partial_sum(:) + eri_3center(:,index_pair(iglobal,jglobal)) * SUM( p_matrix(ilocal,jlocal,:) )  

   enddo
 enddo

 call xtrans_sum(partial_sum)


 ! Hartree potential is not sensitive to spin
 hartree_ij(:,:) = 0.0_dp
 do jlocal=1,n_ham
   jglobal = colindex_local_to_global('H',jlocal)

   do ilocal=1,m_ham
     iglobal = rowindex_local_to_global('H',ilocal)
     if( negligible_basispair(iglobal,jglobal) ) cycle

     hartree_ij(ilocal,jlocal) = SUM( partial_sum(:) * eri_3center(:,index_pair(iglobal,jglobal)) )

   enddo
 enddo

 
 deallocate(partial_sum)

 !
 ! Sum up the different contribution from different procs only if needed
 call xlocal_sum(hartree_ij)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,hartree_ij)

 !
 ! Calculate the Hartree energy
 if( cntxt_ham > 0 ) then
   ehartree = 0.5_dp*SUM(hartree_ij(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
 else
   ehartree = 0.0_dp
 endif
 call xsum(ehartree)


 call stop_clock(timing_hartree)


#endif

end subroutine setup_hartree_ri_sca


!=========================================================================
subroutine setup_exchange_ri_sca(print_matrix_,nbf,m_ham,n_ham,p_matrix_occ,p_matrix_sqrt,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix_occ(nbf,nspin)
 real(dp),intent(in)  :: p_matrix_sqrt(m_ham,n_ham,nspin)
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: exchange_ij(m_ham,n_ham,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nocc
 real(dp),allocatable :: tmpa(:,:)
 real(dp),allocatable :: tmpb(:,:)
 real(dp),allocatable :: tmpc(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
 real(dp)             :: p_matrix_i(nbf)
 integer              :: nbf_trans
 integer              :: iglobal,jglobal,ilocal,jlocal
 integer              :: ii
 integer              :: ibf_global,ier,to,iprow_recv,ipcol_recv
 integer,external     :: INDXG2P
!=====

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity: SCALAPACK'
 call start_clock(timing_exchange)

 nbf_trans = 0
 do ibf=1,nbf
   if( MODULO(ibf-1,nproc_trans) == rank_trans ) then
     nbf_trans = nbf_trans + 1
   endif
 enddo


 allocate(tmpa(nauxil_3center,m_ham))
 allocate(tmpb(nauxil_3center,n_ham))


 exchange_ij(:,:,:) = 0.0_dp

 do ispin=1,nspin
   do istate=1,nbf

     if( p_matrix_occ(istate,ispin) < completely_empty ) cycle

     !
     ! First all processors must have the p_matrix for (istate, ispin)
     p_matrix_i(:) = 0.0_dp
     if( cntxt_ham > 0 ) then
       jlocal = colindex_global_to_local('H',istate)
       if( jlocal /= 0 ) then
         do ilocal=1,m_ham
           iglobal = rowindex_local_to_global('H',ilocal)
           p_matrix_i(iglobal) = p_matrix_sqrt(ilocal,jlocal,ispin)
         enddo
       endif
     endif
     call xsum(p_matrix_i)


     allocate(tmpc(nauxil_3center,nbf_trans))

     tmpc(:,:) = 0.0_dp
     do ibf=1,nbf_trans
       ibf_global = rank_trans + (ibf-1) * nproc_trans + 1
       do ii=1,nbf
         if( negligible_basispair(ii,ibf_global) ) cycle
         tmpc(:,ibf) = tmpc(:,ibf) + eri_3center(:,index_pair(ii,ibf_global)) * p_matrix_i(ii)
       enddo
     enddo


     tmpa(:,:) = 0.0_dp
     tmpb(:,:) = 0.0_dp
     
     do ibf=1,nbf_trans
       ibf_global = rank_trans + (ibf-1) * nproc_trans + 1

       iprow_recv = INDXG2P(ibf_global,block_row,0,first_row,nprow_ham)
       do ipcol_recv=0,npcol_ham-1
         to = rank_sca_to_mpi(iprow_recv,ipcol_recv)
         call MPI_SEND(tmpc(:,ibf),nauxil_3center,MPI_DOUBLE_PRECISION,to,ibf_global,comm_trans,ier) 
       enddo

       ipcol_recv = INDXG2P(ibf_global,block_col,0,first_col,npcol_ham)
       do iprow_recv=0,nprow_ham-1
         to = rank_sca_to_mpi(iprow_recv,ipcol_recv)
         call MPI_SEND(tmpc(:,ibf),nauxil_3center,MPI_DOUBLE_PRECISION,to,nbf+ibf_global,comm_trans,ier) 
       enddo

     enddo

     do ilocal=1,m_ham
       iglobal = rowindex_local_to_global('H',ilocal)
       call MPI_RECV(tmpa(:,ilocal),nauxil_3center,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,iglobal,comm_trans,MPI_STATUS_IGNORE,ier) 
     enddo

     do jlocal=1,n_ham
       jglobal = colindex_local_to_global('H',jlocal)
       call MPI_RECV(tmpb(:,jlocal),nauxil_3center,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,nbf+jglobal,comm_trans,MPI_STATUS_IGNORE,ier) 
     enddo


     deallocate(tmpc)


     exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin)  &
                        - MATMUL( TRANSPOSE(tmpa(:,:)) , tmpb(:,:) ) / spin_fact


   enddo
 enddo

 call xlocal_sum(exchange_ij)

 !
 ! Calculate the Hartree energy
 if( cntxt_ham > 0 ) then
   eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )
 else
   eexchange = 0.0_dp
 endif
 call xsum(eexchange)

 call stop_clock(timing_exchange)


#endif

end subroutine setup_exchange_ri_sca


!=========================================================================
subroutine setup_exchange_longrange_ri_sca(print_matrix_,nbf,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: occupation(nbf,nspin)
 real(dp),intent(in)  :: c_matrix(nbf,nbf,nspin)
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

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity: SCALAPACK'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 allocate(tmp(nauxil_3center_lr,nbf))

 do ispin=1,nspin

   ! Denombrate the strictly positive eigenvalues
   nocc = COUNT( occupation(:,ispin) > completely_empty )

   do istate=1,nocc
     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix(jbf,istate,ispin) * eri_3center_lr(:,ipair) * SQRT( occupation(istate,ispin) )
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center_lr(:,ipair) * SQRT( occupation(istate,ispin) )
     enddo

     exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
                        - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
   enddo

 enddo
 deallocate(tmp)

 call xsum(exchange_ij)

 call dump_out_matrix(print_matrix_,'=== LR Exchange contribution ===',nbf,nspin,exchange_ij)

 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

#endif

end subroutine setup_exchange_longrange_ri_sca


!=========================================================================
subroutine setup_density_matrix_sca(nbf,nstate,m_c,n_c,c_matrix,occupation,m_ham,n_ham,p_matrix)
 implicit none
 integer,intent(in)   :: nbf,nstate
 integer,intent(in)   :: m_ham,n_ham,m_c,n_c
 real(dp),intent(in)  :: c_matrix(m_c,n_c,nspin)
 real(dp),intent(in)  :: occupation(nbf,nspin)
 real(dp),intent(out) :: p_matrix(m_ham,n_ham,nspin)
!=====
 integer  :: ispin,jlocal,jglobal
 real(dp) :: matrix_tmp(m_ham,n_ham)
!=====

#ifdef HAVE_SCALAPACK

 if( cntxt_ham > 0 ) then
   do ispin=1,nspin
     do jlocal=1,n_c
       jglobal = colindex_local_to_global('H',jlocal)
       matrix_tmp(:,jlocal) = c_matrix(:,jlocal,ispin) * SQRT( occupation(jglobal,ispin) )
     enddo

     call PDGEMM('N','T',nbf,nbf,nstate,1.0_dp,matrix_tmp,1,1,desc_c,       &
                  matrix_tmp,1,1,desc_c,0.0_dp,                             &
                  p_matrix,1,1,desc_ham)


   enddo
 else
   p_matrix(:,:,:) = 0.0_dp
 endif

 ! Poor man distribution
 call xlocal_sum(p_matrix)


#endif

end subroutine setup_density_matrix_sca


!=========================================================================
subroutine diagonalize_hamiltonian_sca(nspin_local,nbf,nstate,m_ham,n_ham,  &
                                       hamiltonian,s_matrix_sqrt_inv,energy,m_c,n_c,c_matrix)
 implicit none

 integer,intent(in)   :: nspin_local,nbf,nstate
 integer,intent(in)   :: m_ham,n_ham
 integer,intent(in)   :: m_c,n_c
 real(dp),intent(in)  :: hamiltonian(m_ham,n_ham,nspin_local)
 real(dp),intent(in)  :: s_matrix_sqrt_inv(m_c,n_c)
 real(dp),intent(out) :: c_matrix(m_c,n_c,nspin_local)
 real(dp),intent(out) :: energy(nstate,nspin_local)
!=====
 integer  :: ispin,ibf,jbf,istate
 integer  :: ilocal,jlocal,jglobal
 integer  :: m_small,n_small
 real(dp),allocatable :: h_small(:,:)
!=====

#ifdef HAVE_SCALAPACK



 if(cntxt_ham > 0 ) then
   call init_desc('H',nstate,nstate,desc_small,m_small,n_small)
   allocate(h_small(m_small,n_small))

   do ispin=1,nspin_local
     write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
     call start_clock(timing_diago_hamiltonian)

!     h_small(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
!                              MATMUL( hamiltonian(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )

     !
     ! H_small = ^tS^{-1/2} H S^{-1/2}
     call PDGEMM('N','N',nbf,nstate,nbf,                          &
                  1.0_dp,hamiltonian(:,:,ispin),1,1,desc_ham,     &
                  s_matrix_sqrt_inv,1,1,desc_c,                  &
                  0.0_dp,c_matrix(:,:,ispin),1,1,desc_c)

     call PDGEMM('T','N',nstate,nstate,nbf,                       &
                  1.0_dp,s_matrix_sqrt_inv,1,1,desc_c,           &
                  c_matrix(:,:,ispin),1,1,desc_c,                         &
                  0.0_dp,h_small,1,1,desc_small)



     call diagonalize_sca(desc_small,nstate,m_small,n_small,h_small,energy(:,ispin))


!     c_matrix(:,1:nstate,ispin) = MATMUL( s_matrix_sqrt_inv(:,:) , h_small(:,:) )

     !
     ! C = S^{-1/2} C_small 
     call PDGEMM('N','N',nbf,nstate,nstate,                   &
                  1.0_dp,s_matrix_sqrt_inv,1,1,desc_c,       &
                  h_small,1,1,desc_small,                     &
                  0.0_dp,c_matrix(:,:,ispin),1,1,desc_c) 


     call stop_clock(timing_diago_hamiltonian)
   enddo

   deallocate(h_small)

 else
   energy(:,:) = 0.0_dp
   c_matrix(:,:,:) = 0.0_dp
 endif

 ! Poor man distribution
 call xlocal_sum(energy)
 call xlocal_sum(c_matrix)



#endif

end subroutine diagonalize_hamiltonian_sca


!=========================================================================
subroutine setup_sqrt_overlap_sca(TOL_OVERLAP,nbf,m_ham,n_ham,s_matrix,nstate,m_c,n_c,s_matrix_sqrt_inv)
 use m_tools
 implicit none

 real(dp),intent(in)                :: TOL_OVERLAP
 integer,intent(in)                 :: nbf,m_ham,n_ham
 real(dp),intent(in)                :: s_matrix(m_ham,n_ham)
 integer,intent(out)                :: nstate,m_c,n_c
 real(dp),allocatable,intent(inout) :: s_matrix_sqrt_inv(:,:)
!=====
 real(dp) :: TOL_OVERLAP_FAKE=-1.0_dp
 real(dp) :: matrix_tmp(m_ham,n_ham)
 integer  :: ibf,jbf
 integer  :: ilocal,jlocal
 integer  :: iglobal,jglobal
 real(dp) :: s_eigval(nbf)
 real(dp),allocatable :: diag(:,:)
!=====

#ifdef HAVE_SCALAPACK

 write(stdout,'(/,a)') ' Calculate overlap matrix square-root S^{1/2}: SCALAPACK'

 if( cntxt_ham > 0 ) then
   matrix_tmp(:,:) = s_matrix(:,:)
   call diagonalize_sca(desc_ham,nbf,m_ham,n_ham,matrix_tmp,s_eigval)

   nstate = COUNT( s_eigval(:) > TOL_OVERLAP )

   ! 
   ! Initialize the descriptor of the rectangular matric S^{-1/2}
   call init_desc('H',nbf,nstate,desc_c,m_c,n_c)
   
 else
   nstate = 0
   m_c    = 0
   n_c    = 0
 endif
 ! Propagate nstate
 call xlocal_max(nstate)
 call xlocal_max(m_c)
 call xlocal_max(n_c)


 allocate(s_matrix_sqrt_inv(m_c,n_c))

 write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
 write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
 write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',TOL_OVERLAP
 write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',nbf

 !
 ! Whether a filtering is necessary
 if( nstate == nbf ) then

   if( cntxt_ham > 0 ) then

     do jlocal=1,n_ham
       jglobal = colindex_local_to_global('H',jlocal)
       s_matrix_sqrt_inv(:,jlocal) = matrix_tmp(:,jlocal) / SQRT( s_eigval(jglobal) )
     enddo

   else
     s_matrix_sqrt_inv(:,:) = 0.0_dp
   endif

 else

   if( cntxt_ham > 0 ) then

     ! Create the diagonal matrix than transforms a square matrix into a
     ! rectangular one
     allocate(diag(m_c,n_c))
     diag(:,:) = 0.0_dp

     jglobal = 0
     do iglobal=1,nbf
       if( s_eigval(iglobal) > TOL_OVERLAP ) then
         jglobal = jglobal + 1
         ilocal = rowindex_global_to_local('H',iglobal)
         jlocal = colindex_global_to_local('H',jglobal)
         if( ilocal * jlocal /= 0 ) then
           diag(ilocal,jlocal) = 1.0_dp / SQRT( s_eigval(iglobal) )
         endif
       endif
     enddo


     call PDGEMM('N','N',nbf,nstate,nbf,               &
                  1.0_dp,matrix_tmp,1,1,desc_ham,      &
                  diag,1,1,desc_c,                    &
                  0.0_dp,s_matrix_sqrt_inv,1,1,desc_c)

     deallocate(diag)

   else
     s_matrix_sqrt_inv(:,:) = 0.0_dp
   endif


 endif
 call xlocal_sum(s_matrix_sqrt_inv)


#endif

end subroutine setup_sqrt_overlap_sca


!=========================================================================
subroutine setup_sqrt_density_matrix_sca(nbf,m_ham,n_ham,p_matrix,p_matrix_sqrt,p_matrix_occ)
 use m_tools
 implicit none

 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: p_matrix_sqrt(m_ham,n_ham,nspin)
 real(dp),intent(out) :: p_matrix_occ(nbf,nspin)
!=====
 integer              :: ispin,jlocal,jglobal
!=====

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate the square root of the density matrix: SCALAPACK'
 call start_clock(timing_sqrt_density_matrix)

 if( cntxt_ham > 0 ) then
   do ispin=1,nspin
     p_matrix_sqrt(:,:,ispin) = p_matrix(:,:,ispin)
     call diagonalize_sca(desc_ham,nbf,m_ham,n_ham,p_matrix_sqrt(:,:,ispin),p_matrix_occ(:,ispin))
   enddo
   do jlocal=1,n_ham
     jglobal = colindex_local_to_global('H',jlocal)
     p_matrix_sqrt(:,jlocal,ispin) = p_matrix_sqrt(:,jlocal,ispin) * SQRT( p_matrix_occ(jglobal,ispin) )
   enddo
 else
   p_matrix_sqrt(:,:,:) = 0.0_dp
   p_matrix_occ(:,:)    = 0.0_dp
 endif

 call xlocal_sum(p_matrix_sqrt)
 call xlocal_sum(p_matrix_occ)

 call stop_clock(timing_sqrt_density_matrix)


#endif

end subroutine setup_sqrt_density_matrix_sca


!=========================================================================
subroutine dft_approximate_vhxc_sca(basis,m_ham,n_ham,vhxc_ij)
 use m_basis_set
 use m_dft_grid
 use m_eri_calculate
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: vhxc_ij(m_ham,n_ham)
!=====
 integer              :: idft_xc
 integer              :: igrid,ibf,jbf,ispin
 real(dp)             :: rr(3)
 real(dp)             :: normalization
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: rhor
 real(dp)             :: vxc,exc
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhartree
 real(dp)             :: vhgau(m_ham,n_ham)
 integer              :: iatom,igau,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

#ifdef HAVE_SCALAPACK

 vhxc_ij(:,:) = 0.0_dp


 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities: SCALAPACK'

 do iatom=1,natom
   if( rank_local /= MODULO(iatom,nproc_local) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),coeff,alpha)


   do igau=1,ngau
     call calculate_eri_approximate_hartree(.FALSE.,basis,m_ham,n_ham,x(:,iatom),alpha(igau),vhgau)
     vhxc_ij(:,:) = vhxc_ij(:,:) + vhgau(:,:) * coeff(igau) / 2.0_dp**1.25_dp / pi**0.75_dp * alpha(igau)**1.5_dp
   enddo

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

   call teter_lda_vxc_exc(rhor,vxc,exc)

   !
   ! HXC
   do jlocal=1,n_ham
     jglobal = colindex_local_to_global('H',jlocal)
     do ilocal=1,m_ham
       iglobal = rowindex_local_to_global('H',ilocal)

       vhxc_ij(ilocal,jlocal) =  vhxc_ij(ilocal,jlocal) &
            + weight * vxc * basis_function_r(iglobal)  &
                           * basis_function_r(jglobal)

     enddo
   enddo

 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum(normalization)
 call xlocal_sum(vhxc_ij)

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()


#endif

end subroutine dft_approximate_vhxc_sca


!=========================================================================
end module m_hamiltonian_sca
!=========================================================================
