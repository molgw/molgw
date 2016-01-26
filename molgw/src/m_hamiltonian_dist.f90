!=========================================================================
module m_hamiltonian_dist
 use m_definitions
 use m_mpi
 use m_timing
 use m_warning
 use m_inputparam,only: nspin,spin_fact


contains


!=========================================================================
subroutine reduce_hamiltonian_sca(nbf,matrix_global,m_ham,n_ham,matrix_local)
 implicit none

 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: matrix_global(nbf,nbf)
 real(dp),intent(out) :: matrix_local(m_ham,n_ham)
!=====
 integer              :: ipcol,iprow,rank_dest
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: m_block,n_block
 real(dp),allocatable :: matrix_block(:,:)
!=====

#ifdef HAVE_SCALAPACK
 call start_clock(timing_sca_distr)

 ! Loops over the SCALAPACK grid
 do ipcol=0,npcol_ham-1
   do iprow=0,nprow_ham-1

     ! Identify the destination processor
     rank_dest = rank_ham_sca_to_mpi(iprow,ipcol)

     m_block = row_block_size(nbf,iprow,nprow_ham)
     n_block = col_block_size(nbf,ipcol,npcol_ham)
     allocate(matrix_block(m_block,n_block))

     do jlocal=1,n_block
       jglobal = colindex_local_to_global(ipcol,npcol_ham,jlocal)
       do ilocal=1,m_block
         iglobal = rowindex_local_to_global(iprow,nprow_ham,ilocal)

         matrix_block(ilocal,jlocal) = matrix_global(iglobal,jglobal)

       enddo
     enddo


     call xsum(rank_dest,matrix_block)

     if( rank == rank_dest ) then
       matrix_local(:,:) = matrix_block(:,:)
     endif
     deallocate(matrix_block)

   enddo
 enddo


 call stop_clock(timing_sca_distr)

#else
 call die('SCALAPACK is required for distribute_hamiltonian_sca')
#endif

end subroutine reduce_hamiltonian_sca


!=========================================================================
subroutine broadcast_hamiltonian_sca(nbf,matrix_global,m_ham,n_ham,matrix_local)
 implicit none

 integer,intent(in)     :: nbf,m_ham,n_ham
 real(dp),intent(inout) :: matrix_global(nbf,nbf)
 real(dp),intent(in)    :: matrix_local(m_ham,n_ham)
!=====
 integer              :: ipcol,iprow,rank_orig
 integer              :: ier
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: m_block,n_block
 real(dp),allocatable :: matrix_block(:,:)
!=====

#ifdef HAVE_SCALAPACK
 call start_clock(timing_sca_distr)

 ! Loops over the SCALAPACK grid
 do ipcol=0,npcol_ham-1
   do iprow=0,nprow_ham-1

     ! Identify the destination processor
     rank_orig = rank_ham_sca_to_mpi(iprow,ipcol)

     m_block = row_block_size(nbf,iprow,nprow_ham)
     n_block = col_block_size(nbf,ipcol,npcol_ham)
     allocate(matrix_block(m_block,n_block))

     if( rank == rank_orig ) then
       matrix_block(:,:) = matrix_local(:,:)
     endif

!     call MPI_BCAST(matrix_block,m_block * n_block,MPI_DOUBLE_PRECISION,rank_orig,comm_world,ier)

     call xbcast(rank_orig,matrix_block)


     do jlocal=1,n_block
       jglobal = colindex_local_to_global(ipcol,npcol_ham,jlocal)
       do ilocal=1,m_block
         iglobal = rowindex_local_to_global(iprow,nprow_ham,ilocal)

         matrix_global(iglobal,jglobal) = matrix_global(iglobal,jglobal) + matrix_block(ilocal,jlocal)

       enddo
     enddo

     deallocate(matrix_block)

   enddo
 enddo


 call stop_clock(timing_sca_distr)

#else
 call die('SCALAPACK is required for distribute_hamiltonian_sca')
#endif

end subroutine broadcast_hamiltonian_sca


!=========================================================================
subroutine setup_nucleus_buffer_sca(print_matrix_,basis,m_ham,n_ham,hamiltonian_nucleus)
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
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
 real(dp),allocatable :: buffer(:,:)
!=====

#ifdef HAVE_SCALAPACK

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian: SCALAPACK buffer'
 if( nproc > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank /= MODULO(iatom-1,nproc) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif

 allocate(buffer(basis%nbf,basis%nbf))

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
     buffer(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
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


 ! Sum up the buffers and store the result in the sub matrix hamiltonian_nucleus
 call reduce_hamiltonian_sca(basis%nbf,buffer,m_ham,n_ham,hamiltonian_nucleus)

 deallocate(buffer)


 call stop_clock(timing_hamiltonian_nuc)


#endif


end subroutine setup_nucleus_buffer_sca


!=========================================================================
subroutine setup_hartree_ri_buffer_sca(print_matrix_,nbf,m_ham,n_ham,p_matrix,pot_hartree,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: pot_hartree(m_ham,n_ham)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf
 integer              :: ipair
 real(dp),allocatable :: partial_sum(:)
 real(dp)             :: rtmp
 real(dp),allocatable :: buffer(:,:)
!=====


#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity: SCALAPACK buffer'
 call start_clock(timing_hartree)

 allocate(buffer(nbf,nbf))
 buffer(:,:) = 0.0_dp

 call broadcast_hamiltonian_sca(nbf,buffer,m_ham,n_ham,p_matrix(:,:,1))
 if( nspin == 2 ) then
   call broadcast_hamiltonian_sca(nbf,buffer,m_ham,n_ham,p_matrix(:,:,2))
 endif

 allocate(partial_sum(nauxil_3center))

 partial_sum(:) = 0.0_dp
 do ipair=1,npair
   kbf = index_basis(1,ipair)
   lbf = index_basis(2,ipair)
   ! Factor 2 comes from the symmetry of p_matrix
   partial_sum(:) = partial_sum(:) + eri_3center(:,ipair) * buffer(kbf,lbf) * 2.0_dp
   ! Then diagonal terms have been counted twice and should be removed once.
   if( kbf == lbf ) &
     partial_sum(:) = partial_sum(:) - eri_3center(:,ipair) * buffer(kbf,kbf)
 enddo


 ! Hartree potential is not sensitive to spin
 buffer(:,:) = 0.0_dp
 do ipair=1,npair
   ibf = index_basis(1,ipair)
   jbf = index_basis(2,ipair)
   rtmp = DOT_PRODUCT( eri_3center(:,ipair) , partial_sum(:) )
   buffer(ibf,jbf) = rtmp
   buffer(jbf,ibf) = rtmp
 enddo

 deallocate(partial_sum)

 ! Sum up the buffers and store the result in the sub matrix pot_exchange
 call reduce_hamiltonian_sca(nbf,buffer,m_ham,n_ham,pot_hartree)
 deallocate(buffer)

 !
 ! Calculate the Hartree energy
 if( cntxt_ham > 0 ) then
   ehartree = 0.5_dp*SUM(pot_hartree(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
 else
   ehartree = 0.0_dp
 endif
 call xsum(ehartree)


 call stop_clock(timing_hartree)

#endif


end subroutine setup_hartree_ri_buffer_sca


!=========================================================================
subroutine setup_exchange_ri_buffer_sca(print_matrix_,nbf,m_ham,n_ham,p_matrix_occ,p_matrix_sqrt,p_matrix,pot_exchange,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix_occ(nbf,nspin)
 real(dp),intent(in)  :: p_matrix_sqrt(m_ham,n_ham,nspin)
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: pot_exchange(m_ham,n_ham,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,ispin,istate
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
 real(dp)             :: p_matrix_i(nbf)
 integer              :: iglobal,ilocal,jlocal
 real(dp),allocatable :: buffer(:,:)
!=====


 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity: SCALAPACK buffer'
 call start_clock(timing_exchange)

 allocate(buffer(nbf,nbf))


 buffer(:,:) = 0.0_dp

 allocate(tmp(nauxil_3center,nbf))

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
           p_matrix_i(iglobal) = p_matrix_sqrt(ilocal,jlocal,ispin) * SQRT( p_matrix_occ(istate,ispin) )
         enddo
       endif
     endif
     call xsum(p_matrix_i)


     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)

       tmp(:,ibf) = tmp(:,ibf) + p_matrix_i(jbf) * eri_3center(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + p_matrix_i(ibf) * eri_3center(:,ipair)
     enddo

     buffer(:,:) = buffer(:,:) &
                        - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
   enddo

 enddo
 deallocate(tmp)


 ! Sum up the buffers and store the result in the sub matrix pot_exchange
 call reduce_hamiltonian_sca(nbf,buffer,m_ham,n_ham,pot_exchange)
 deallocate(buffer)

 !
 ! Calculate the exchange energy
 if( cntxt_ham > 0 ) then
   eexchange = 0.5_dp * SUM( pot_exchange(:,:,:) * p_matrix(:,:,:) )
 else
   eexchange = 0.0_dp
 endif
 call xsum(eexchange)

 call stop_clock(timing_exchange)



end subroutine setup_exchange_ri_buffer_sca


!=========================================================================
end module m_hamiltonian_dist
!=========================================================================
