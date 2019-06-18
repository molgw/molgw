!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham hamiltonian with SCALAPACK
! when a single buffer allocation is permitted.
!
!=========================================================================
module m_hamiltonian_buffer
 use m_definitions
 use m_mpi
 use m_timing
 use m_tddft_variables
 use m_warning
 use m_memory
 use m_scalapack
 use m_cart_to_pure
 use m_inputparam
 use m_basis_set
 use m_density_tools
 use m_hamiltonian_onebody
 use m_dft_grid
 use m_eri_calculate
 use m_linear_algebra,only: matrix_trace
 use m_dft_grid
 use m_libxc_tools



 real(dp),private,allocatable :: buffer(:,:)


contains


!=========================================================================
subroutine allocate_parallel_buffer(nbf)
 implicit none

 integer,intent(in) :: nbf
!=====

 write(stdout,'(/,1x,a)') 'For SCALAPACK buffer, only this buffer is not distributed'
 call clean_allocate('non distributed buffer',buffer,nbf,nbf)

end subroutine allocate_parallel_buffer


!=========================================================================
subroutine destroy_parallel_buffer()
 implicit none

!=====

 if( .NOT. ALLOCATED(buffer) ) call die('destroy_parallel_buffer: this should not happen')

 call clean_deallocate('non distributed buffer',buffer)

end subroutine destroy_parallel_buffer


!=========================================================================
subroutine reduce_hamiltonian_sca(matrix_local)
 implicit none

 real(dp),intent(out) :: matrix_local(:,:)
!=====
 integer              :: m_ham,n_ham
 integer              :: nbf
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 call start_clock(timing_sca_distr1)

 nbf   = SIZE(buffer(:,:),DIM=1)
 m_ham = SIZE(matrix_local,DIM=1)
 n_ham = SIZE(matrix_local,DIM=2)

#ifdef HAVE_SCALAPACK

 ! Dirty coding.. however surprisingly efficient
 call xsum_world(buffer)
 do jlocal=1,n_ham
   jglobal = colindex_local_to_global('H',jlocal)
   do ilocal=1,m_ham
     iglobal = rowindex_local_to_global('H',ilocal)

     matrix_local(ilocal,jlocal) = buffer(iglobal,jglobal)

   enddo
 enddo


#else

 call xsum_world(buffer)
 matrix_local(:,:) = buffer(:,:)

#endif

 call stop_clock(timing_sca_distr1)

end subroutine reduce_hamiltonian_sca


!=========================================================================
subroutine broadcast_hamiltonian_sca(matrix_local)
 implicit none

 real(dp),intent(in)  :: matrix_local(:,:)
!=====
 integer              :: m_ham,n_ham
 integer              :: nbf
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 call start_clock(timing_sca_distr2)

 nbf = SIZE(buffer(:,:),DIM=1)
 m_ham = SIZE(matrix_local,DIM=1)
 n_ham = SIZE(matrix_local,DIM=2)

#ifdef HAVE_SCALAPACK

 ! buffer is first divided by the number of procs,
 ! since there will be an allreduce operation at the end
 buffer(:,:) = buffer(:,:) / REAL( nproc_world , dp )
 if( cntxt_ham > 0 ) then
   do jlocal=1,n_ham
     jglobal = colindex_local_to_global('H',jlocal)
     do ilocal=1,m_ham
       iglobal = rowindex_local_to_global('H',ilocal)

       buffer(iglobal,jglobal) = buffer(iglobal,jglobal) + matrix_local(ilocal,jlocal)

     enddo
   enddo
 endif
 call xsum_world(buffer)

#else

 buffer(:,:) = buffer(:,:) + matrix_local(:,:)

#endif

 call stop_clock(timing_sca_distr2)

end subroutine broadcast_hamiltonian_sca


!=========================================================================
subroutine setup_nucleus_buffer_sca(basis,hamiltonian_nucleus,atom_list)
 implicit none
 type(basis_set),intent(in)  :: basis
 real(dp),intent(out)        :: hamiltonian_nucleus(:,:)
 integer,intent(in),optional :: atom_list(:)
!=====
!=====

 buffer(:,:) = 0.0_dp
 if( PRESENT(atom_list) ) then
   call setup_nucleus(basis,buffer,atom_list)
 else
   call setup_nucleus(basis,buffer)
 endif

 ! Sum up the buffers and store the result in the sub matrix hamiltonian_nucleus
 buffer(:,:) = buffer(:,:) / REAL(nproc_world,dp)
 call reduce_hamiltonian_sca(hamiltonian_nucleus)

end subroutine setup_nucleus_buffer_sca


!=========================================================================
subroutine setup_hartree_ri_buffer_sca(p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: hartree_ij(:,:)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf
 integer              :: ipair
 real(dp),allocatable :: partial_sum(:)
 real(dp)             :: rtmp
!=====

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity: SCALAPACK buffer'


 !
 ! First the buffer contains the density matrix p_matrix
 buffer(:,:) = 0.0_dp

 call broadcast_hamiltonian_sca(p_matrix(:,:,1))
 if( nspin == 2 ) then
   call broadcast_hamiltonian_sca(p_matrix(:,:,2))
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

 ! Sum up the buffers and store the result in the sub matrix hartree_ij
 call reduce_hamiltonian_sca(hartree_ij)

 !
 ! Calculate the Hartree energy
 if( cntxt_ham > 0 ) then
   ehartree = 0.5_dp * SUM( hartree_ij(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
 else
   ehartree = 0.0_dp
 endif
 call xsum_world(ehartree)



end subroutine setup_hartree_ri_buffer_sca


!=========================================================================
subroutine setup_exchange_ri_buffer_sca(occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ij(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nbf
 integer              :: nstate,m_c
 integer              :: ibf,jbf,ispin,istate
 integer              :: ipair
 real(dp),allocatable :: tmp(:,:)
 real(dp),allocatable :: c_matrix_i(:)
 integer              :: iglobal,ilocal,jlocal
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity: SCALAPACK buffer'

 nbf    = SIZE(buffer(:,:),DIM=1)
 m_c    = SIZE(c_matrix,DIM=1)
 nstate = SIZE(occupation,DIM=1)


 allocate(tmp(nauxil_3center,nbf))
 allocate(c_matrix_i(nbf))

 do ispin=1,nspin

   buffer(:,:) = 0.0_dp

   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle

     !
     ! First all processors must have the c_matrix for (istate, ispin)
     c_matrix_i(:) = 0.0_dp
     if( cntxt_ham > 0 ) then
       jlocal = colindex_global_to_local('H',istate)
       if( jlocal /= 0 ) then
         do ilocal=1,m_c
           iglobal = rowindex_local_to_global('H',ilocal)
           c_matrix_i(iglobal) = c_matrix(ilocal,jlocal,ispin)
         enddo
       endif
     endif
     call xsum_world(c_matrix_i)


     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix_i(jbf) * eri_3center(:,ipair)
       if( ibf /= jbf ) &
         tmp(:,jbf) = tmp(:,jbf) + c_matrix_i(ibf) * eri_3center(:,ipair)
     enddo


     ! buffer(:,:) = buffer(:,:) - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact * occ(i)
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nauxil_3center,-occupation(istate,ispin)/spin_fact,tmp,nauxil_3center,1.0_dp,buffer,nbf)

   enddo

   !
   ! Need to symmetrize buffer
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       buffer(ibf,jbf) = buffer(jbf,ibf)
     enddo
   enddo

   ! Sum up the buffers and store the result in the sub matrix exchange_ij
   call reduce_hamiltonian_sca(exchange_ij(:,:,ispin))

 enddo
 deallocate(tmp)


 !
 ! Calculate the exchange energy
 if( cntxt_ham > 0 ) then
   eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )
 else
   eexchange = 0.0_dp
 endif
 call xsum_world(eexchange)

 call stop_clock(timing_exchange)



end subroutine setup_exchange_ri_buffer_sca


!=========================================================================
subroutine setup_exchange_longrange_ri_buffer_sca(occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ij(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nbf
 integer              :: nstate,m_c
 integer              :: ibf,jbf,ispin,istate
 real(dp),allocatable :: tmp(:,:)
 integer              :: ipair
 real(dp)             :: c_matrix_i(desc_c(M_))
 integer              :: iglobal,ilocal,jlocal
!=====


 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity: SCALAPACK buffer'

 if( npcol_3center > 1 ) call die('setup_exchange_longrange_ri_buffer_sca: npcol_3center > 1 not allowed')

 nbf    = desc_c(M_)
 m_c    = SIZE(c_matrix,DIM=1)
 nstate = SIZE(occupation,DIM=1)


 allocate(tmp(nauxil_3center_lr,nbf))

 do ispin=1,nspin

   buffer(:,:) = 0.0_dp

   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle

     !
     ! First all processors must have the c_matrix for (istate, ispin)
     c_matrix_i(:) = 0.0_dp
     if( cntxt_ham > 0 ) then
       jlocal = colindex_global_to_local('H',istate)
       if( jlocal /= 0 ) then
         do ilocal=1,m_c
           iglobal = rowindex_local_to_global('H',ilocal)
           c_matrix_i(iglobal) = c_matrix(ilocal,jlocal,ispin)
         enddo
       endif
     endif
     call xsum_world(c_matrix_i)


     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix_i(jbf) * eri_3center_lr(:,ipair)
       if( ibf /= jbf ) &
          tmp(:,jbf) = tmp(:,jbf) + c_matrix_i(ibf) * eri_3center_lr(:,ipair)
     enddo


     ! buffer(:,:) = buffer(:,:) - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact * occ(i)
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nauxil_3center_lr,-occupation(istate,ispin)/spin_fact,tmp,nauxil_3center_lr,1.0_dp,buffer,nbf)

   enddo

   !
   ! Need to symmetrize buffer
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       buffer(ibf,jbf) = buffer(jbf,ibf)
     enddo
   enddo

   ! Sum up the buffers and store the result in the sub matrix exchange_ij
   call reduce_hamiltonian_sca(exchange_ij(:,:,ispin))

 enddo
 deallocate(tmp)


 !
 ! Calculate the exchange energy
 if( cntxt_ham > 0 ) then
   eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )
 else
   eexchange = 0.0_dp
 endif
 call xsum_world(eexchange)

 call stop_clock(timing_exchange)



end subroutine setup_exchange_longrange_ri_buffer_sca



!=========================================================================
end module m_hamiltonian_buffer
!=========================================================================
