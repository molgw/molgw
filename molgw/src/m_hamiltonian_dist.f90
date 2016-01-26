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


 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity with large buffer'
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
