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
 use m_warning
 use m_memory
 use m_scalapack
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact
 use m_basis_set
 use m_hamiltonian_onebody


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
subroutine setup_overlap_buffer_sca(basis,overlap)
 use m_atoms
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: overlap(:,:)
!=====
!=====

 buffer(:,:) = 0.0_dp
 call setup_overlap(basis,buffer)

 ! Sum up the buffers and store the result in the sub matrix overlap
 buffer(:,:) = buffer(:,:) / REAL(nproc_world,dp)
 call reduce_hamiltonian_sca(overlap)

end subroutine setup_overlap_buffer_sca


!=========================================================================
subroutine setup_kinetic_buffer_sca(basis,hamiltonian_kinetic)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(:,:)
!=====
!=====

 buffer(:,:) = 0.0_dp
 call setup_kinetic(basis,buffer)

 ! Sum up the buffers and store the result in the sub matrix hamiltonian_kinetic
 buffer(:,:) = buffer(:,:) / REAL(nproc_world,dp)
 call reduce_hamiltonian_sca(hamiltonian_kinetic)

end subroutine setup_kinetic_buffer_sca


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
 call start_clock(timing_hartree)


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


 call stop_clock(timing_hartree)


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
     do ipair=1,nbf
       ibf = index_basis(1,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix_i(ibf) * eri_3center(:,ipair)
     enddo
     do ipair=nbf+1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix_i(jbf) * eri_3center(:,ipair)
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
     do ipair=1,nbf
       ibf = index_basis(1,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix_i(ibf) * eri_3center_lr(:,ipair)
     enddo
     do ipair=nbf+1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix_i(jbf) * eri_3center_lr(:,ipair)
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
subroutine dft_exc_vxc_buffer_sca(batch_size,basis,occupation,c_matrix,vxc_ij,exc_xc)
 use m_inputparam
 use m_dft_grid
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 integer,intent(in)         :: batch_size
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(:,:)
 real(dp),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(out)       :: vxc_ij(:,:,:)
 real(dp),intent(out)       :: exc_xc
!=====
 real(dp),parameter   :: TOL_RHO=1.0e-9_dp
 integer              :: nstate
 integer              :: idft_xc
 integer              :: igrid_start,igrid_end,ir,nr
 integer              :: ibf,jbf,ispin
 real(dp)             :: normalization(nspin)
 real(dp)             :: rhor(nspin,ngrid)
 real(dp)             :: grad_rhor(nspin,ngrid,3)
 real(dp),allocatable :: exc_batch(:)
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: rhor_batch(:,:)
 real(dp),allocatable :: grad_rhor_batch(:,:,:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: basis_function_gradr_batch(:,:,:)
 real(dp),allocatable :: dedd_r_batch(:,:)
 real(dp),allocatable :: dedgd_r_batch(:,:,:)
 real(dp),allocatable :: sigma_batch(:,:)
 real(dp),allocatable :: vrho_batch(:,:)
 real(dp),allocatable :: vsigma_batch(:,:)
 real(dp),allocatable :: tmp_batch(:,:)
!=====

! if( nspin/=1 ) call die('DFT XC potential: SCALAPACK buffer not implemented for spin unrestricted')

 nstate = SIZE(occupation,DIM=1)
 exc_xc = 0.0_dp
 vxc_ij(:,:,:) = 0.0_dp
 if( ndft_xc == 0 ) return

 call start_clock(timing_dft)


#ifdef HAVE_LIBXC

 write(stdout,'(1x,a)') 'Calculate DFT XC potential: SCALAPACK buffer'
 if( batch_size /= 1 ) write(stdout,'(1x,a,1x,i4)') 'Using batches of size',batch_size


 normalization(:) = 0.0_dp

 do ispin=1,nspin
   !
   ! Buffer constains the c_matrix for a spin channel ispin
   buffer(:,:) = 0.0_dp
   call broadcast_hamiltonian_sca(c_matrix(:,:,ispin))


   !
   ! Loop over batches of grid points
   !
   do igrid_start=1,ngrid,batch_size
     igrid_end = MIN(ngrid,igrid_start+batch_size-1)
     nr = igrid_end - igrid_start + 1

     allocate(weight_batch(nr))
     allocate(basis_function_r_batch(basis%nbf,nr))
     allocate(rhor_batch(1,nr))

     if( dft_xc_needs_gradient ) then
       allocate(basis_function_gradr_batch(basis%nbf,nr,3))
       allocate(grad_rhor_batch(1,nr,3))
     endif

     weight_batch(:) = w_grid(igrid_start:igrid_end)


     !
     ! Get all the functions at point r
     call get_basis_functions_r_batch(basis,igrid_start,nr,basis_function_r_batch)

     if( dft_xc_needs_gradient ) call get_basis_functions_gradr_batch(basis,igrid_start,nr,basis_function_gradr_batch)

     !
     ! Calculate the density at points r for spin up and spin down
     ! Calculate grad rho at points r for spin up and spin down
     if( .NOT. dft_xc_needs_gradient ) then
       call calc_density_r_batch(1,basis%nbf,nstate,nr,occupation(:,ispin),buffer(:,1:nstate),basis_function_r_batch,rhor_batch)
     else
       call calc_density_gradr_batch(1,basis%nbf,nstate,nr,occupation(:,ispin),buffer(:,1:nstate), &
                                     basis_function_r_batch,basis_function_gradr_batch,rhor_batch,grad_rhor_batch)
     endif

     ! Save the whole rhor and gradr
     rhor(ispin,igrid_start:igrid_end) = rhor_batch(1,1:nr)
     if( dft_xc_needs_gradient ) then
       grad_rhor(ispin,igrid_start:igrid_end,:) = grad_rhor_batch(1,1:nr,:)
     endif

     !
     ! Normalization
     normalization(ispin) = normalization(ispin) + DOT_PRODUCT( rhor_batch(1,:) , weight_batch(:) )

     deallocate(weight_batch)
     deallocate(basis_function_r_batch)
     deallocate(rhor_batch)

     if( dft_xc_needs_gradient ) then
       deallocate(basis_function_gradr_batch)
       deallocate(grad_rhor_batch)
     endif


   enddo
 enddo


 do ispin=1,nspin

   !
   ! buffer now contains the vxc_ij for one spin
   buffer(:,:) = 0.0_dp

   !
   ! Loop over batches of grid points
   !
   do igrid_start=1,ngrid,batch_size
     igrid_end = MIN(ngrid,igrid_start+batch_size-1)
     nr = igrid_end - igrid_start + 1


     ! Skip if the density is too small for the whole batch
     if( ALL( rhor(:,igrid_start:igrid_end) < TOL_RHO ) ) cycle

     allocate(weight_batch(nr))
     allocate(dedd_r_batch(nspin,nr))
     allocate(vrho_batch(nspin,nr))
     allocate(basis_function_r_batch(basis%nbf,nr))
     allocate(exc_batch(nr))
     allocate(rhor_batch(nspin,nr))

     if( dft_xc_needs_gradient ) then
       allocate(grad_rhor_batch(nspin,nr,3))
       allocate(sigma_batch(2*nspin-1,nr))
       allocate(vsigma_batch(2*nspin-1,nr))
       allocate(dedgd_r_batch(3,nr,nspin))
       allocate(basis_function_gradr_batch(basis%nbf,nr,3))
     endif

     weight_batch(:) = w_grid(igrid_start:igrid_end)

     !
     ! Get the batch density and gradient from the saved rhor and gradr
     rhor_batch(:,1:nr) = rhor(:,igrid_start:igrid_end)
     if( dft_xc_needs_gradient ) then
        grad_rhor_batch(:,1:nr,:) = grad_rhor(:,igrid_start:igrid_end,:)
     endif

     if( dft_xc_needs_gradient ) then
       do ir=1,nr
         sigma_batch(1,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(1,ir,:) )
         if( nspin == 2 ) then
           sigma_batch(2,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(2,ir,:) )
           sigma_batch(3,ir) = DOT_PRODUCT( grad_rhor_batch(2,ir,:) , grad_rhor_batch(2,ir,:) )
         endif
       enddo
     endif

     !
     ! LIBXC calls
     !
     dedd_r_batch(:,:) = 0.0_dp
     if( dft_xc_needs_gradient ) dedgd_r_batch(:,:,:) = 0.0_dp

     do idft_xc=1,ndft_xc
       if( ABS(dft_xc_coef(idft_xc)) < 1.0e-6_dp ) cycle

       select case(xc_f90_info_family(calc_type%xc_info(idft_xc)))
       case(XC_FAMILY_LDA)
         call xc_f90_lda_exc_vxc(calc_type%xc_func(idft_xc),nr,rhor_batch(1,1),exc_batch(1),vrho_batch(1,1))

       case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         call xc_f90_gga_exc_vxc(calc_type%xc_func(idft_xc),nr,rhor_batch(1,1),sigma_batch(1,1),exc_batch(1),vrho_batch(1,1),vsigma_batch(1,1))

         ! Remove too small densities to stabilize the computation
         ! especially useful for Becke88
         do ir=1,nr
           if( ALL( rhor_batch(:,ir) < TOL_RHO ) ) then
             exc_batch(ir)      = 0.0_dp
             vrho_batch(:,ir)   = 0.0_dp
             vsigma_batch(:,ir) = 0.0_dp
           endif
         enddo

       case default
         call die('functional is not LDA nor GGA nor hybrid')
       end select


       exc_xc = exc_xc + SUM( weight_batch(:) * exc_batch(:) * rhor_batch(ispin,:) ) * dft_xc_coef(idft_xc)


       dedd_r_batch(:,:) = dedd_r_batch(:,:) + vrho_batch(:,:) * dft_xc_coef(idft_xc)

       !
       ! Set up divergence term if needed (GGA case)
       !
       if( dft_xc_needs_gradient ) then
         do ir=1,nr
           if( nspin == 1 ) then

             dedgd_r_batch(:,ir,1) = dedgd_r_batch(:,ir,1)  &
                        + 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) * dft_xc_coef(idft_xc)

           else

             dedgd_r_batch(:,ir,1) = dedgd_r_batch(:,ir,1) &
                       + ( 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) &
                                   + vsigma_batch(2,ir) * grad_rhor_batch(2,ir,:) ) * dft_xc_coef(idft_xc)

             dedgd_r_batch(:,ir,2) = dedgd_r_batch(:,ir,2) &
                       + ( 2.0_dp * vsigma_batch(3,ir) * grad_rhor_batch(2,ir,:) &
                                   + vsigma_batch(2,ir) * grad_rhor_batch(1,ir,:) ) * dft_xc_coef(idft_xc)
           endif

         enddo
       endif

     enddo ! loop on the XC functional

     !
     ! Get all the functions at point r
     call get_basis_functions_r_batch(basis,igrid_start,nr,basis_function_r_batch)

     if( dft_xc_needs_gradient ) call get_basis_functions_gradr_batch(basis,igrid_start,nr,basis_function_gradr_batch)

     !
     ! LDA and GGA
     allocate(tmp_batch(basis%nbf,nr))
     forall(ir=1:nr)
       tmp_batch(:,ir) = weight_batch(ir) * dedd_r_batch(ispin,ir) * basis_function_r_batch(:,ir)
     end forall

     call DGEMM('N','T',basis%nbf,basis%nbf,nr,1.0d0,tmp_batch,basis%nbf,basis_function_r_batch,basis%nbf,1.0d0,buffer,basis%nbf)

     !
     ! GGA-only
     if( dft_xc_needs_gradient ) then
       do ir=1,nr
         tmp_batch(:,ir) = MATMUL( basis_function_gradr_batch(:,ir,:) , dedgd_r_batch(:,ir,ispin) * weight_batch(ir) )
       enddo

       call DSYR2K('L','N',basis%nbf,nr,1.0d0,basis_function_r_batch,basis%nbf,tmp_batch,basis%nbf,1.0d0,buffer,basis%nbf)

     endif
     deallocate(tmp_batch)

     deallocate(weight_batch)
     deallocate(dedd_r_batch)
     deallocate(vrho_batch)
     deallocate(basis_function_r_batch)
     deallocate(exc_batch)
     deallocate(rhor_batch)

     if( dft_xc_needs_gradient ) then
       deallocate(sigma_batch)
       deallocate(vsigma_batch)
       deallocate(dedgd_r_batch)
       deallocate(basis_function_gradr_batch)
       deallocate(grad_rhor_batch)
     endif


   enddo ! loop on the batches

   ! Symmetrize now
   do jbf=1,basis%nbf
     do ibf=jbf+1,basis%nbf
       buffer(jbf,ibf) = buffer(ibf,jbf)
     enddo
   enddo

   call reduce_hamiltonian_sca(vxc_ij(:,:,ispin))

 enddo ! spin

 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(exc_xc)


#else
 write(stdout,*) 'XC energy and potential set to zero'
 write(stdout,*) 'LIBXC is not present'
#endif

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
 write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

 call stop_clock(timing_dft)

end subroutine dft_exc_vxc_buffer_sca


!=========================================================================
subroutine dft_approximate_vhxc_buffer_sca(basis,m_ham,n_ham,vhxc_ij)
 use m_dft_grid
 use m_eri_calculate
 use m_tools,only: matrix_trace
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: vhxc_ij(m_ham,n_ham)
!=====
 integer              :: igrid
 real(dp)             :: rr(3)
 real(dp)             :: normalization
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: rhor
 real(dp)             :: vxc,exc,excr
 integer              :: iatom,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 call start_clock(timing_approx_ham)

 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities: buffer SCALAPACK'

 vhxc_ij(:,:) = 0.0_dp

 buffer(:,:) = 0.0_dp
 do iatom=1,natom
   if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zvalence(iatom),zatom(iatom),coeff,alpha)

   !
   ! buffer +=  < i | v_h | j >
   ! buffer is not full. It needs to be symmetrized
   call calculate_eri_approximate_hartree(basis,basis%nbf,basis%nbf,xatom(:,iatom),ngau,coeff,alpha,buffer)

   deallocate(alpha,coeff)
 enddo

 buffer(:,:) = buffer(:,:) + TRANSPOSE(buffer(:,:))


 write(stdout,*) 'Simple LDA functional on a coarse grid'

 !
 ! Create a temporary grid with low quality
 ! This grid is to be destroyed at the end of the present subroutine
 call init_dft_grid(basis,low,.FALSE.,.FALSE.,1)

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
   call setup_atomic_density(rr,rhor)

   !
   ! Normalization
   normalization = normalization + rhor * weight

   call teter_lda_vxc_exc(1,rhor,vxc,excr)

   !
   ! XC energy
   exc = exc + excr * weight * rhor


   !
   ! HXC
   do jglobal=1,basis%nbf
     do iglobal=1,basis%nbf

       buffer(iglobal,jglobal) =  buffer(iglobal,jglobal) &
            + weight * vxc * basis_function_r(iglobal)  &
                           * basis_function_r(jglobal)

     enddo
   enddo
 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(exc)
 call xsum_grid(buffer)

 !
 ! Distribute the result
 if( cntxt_ham > 0 ) then
   do jlocal=1,n_ham
     jglobal = colindex_local_to_global('H',jlocal)
     do ilocal=1,m_ham
       iglobal = rowindex_local_to_global('H',ilocal)

       vhxc_ij(ilocal,jlocal) =  vhxc_ij(ilocal,jlocal) + buffer(iglobal,jglobal)

     enddo
   enddo
 endif


 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization
 write(stdout,  '(a,2(2x,f12.6))') '      XC energy (Ha):',exc

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()


 call stop_clock(timing_approx_ham)

end subroutine dft_approximate_vhxc_buffer_sca


!=========================================================================
end module m_hamiltonian_buffer
!=========================================================================
