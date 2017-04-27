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
subroutine reduce_hamiltonian_sca(m_ham,n_ham,matrix_local)
 implicit none

 integer,intent(in)   :: m_ham,n_ham
 real(dp),intent(out) :: matrix_local(m_ham,n_ham)
!=====
 integer              :: nbf
 integer              :: ipcol,iprow,rank_dest
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: m_block,n_block
 real(dp),allocatable :: matrix_block(:,:)
!=====

 call start_clock(timing_sca_distr1)

 nbf = SIZE(buffer(:,:),DIM=1)

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
subroutine broadcast_hamiltonian_sca(m_ham,n_ham,matrix_local)
 implicit none

 integer,intent(in)   :: m_ham,n_ham
 real(dp),intent(in)  :: matrix_local(m_ham,n_ham)
!=====
 integer              :: nbf
 integer              :: ipcol,iprow,rank_orig
 integer              :: ier
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: m_block,n_block
 real(dp),allocatable :: matrix_block(:,:)
!=====

 call start_clock(timing_sca_distr2)

 nbf = SIZE(buffer(:,:),DIM=1)

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
subroutine setup_nucleus_buffer_sca(print_matrix_,basis,m_ham,n_ham,hamiltonian_nucleus)
 use m_basis_set
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: hamiltonian_nucleus(m_ham,n_ham)
!=====
 integer              :: gt
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: natom_local
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
!=====

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian: SCALAPACK buffer'
 gt = get_gaussian_type_tag(basis%gaussian_type)

 buffer(:,:) = 0.0_dp

 if( nproc_world > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif


 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj        = number_basis_function_am(basis%gaussian_type,lj)
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni        = number_basis_function_am(basis%gaussian_type,li)
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(matrix_cart(ni_cart,nj_cart))
     matrix_cart(:,:) = 0.0_dp
     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           call nucleus_basis_function(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),zvalence(iatom),xatom(:,iatom),vnucleus_ij)
           matrix_cart(i_cart,j_cart) = matrix_cart(i_cart,j_cart) + vnucleus_ij
         enddo
       enddo
     enddo
     buffer(ibf1:ibf2,jbf1:jbf2) = MATMUL( TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)) , &
                                           MATMUL( matrix_cart(:,:) , cart_to_pure(lj,gt)%matrix(:,:) ) )


     deallocate(matrix_cart)
   enddo
 enddo


 ! Sum up the buffers and store the result in the sub matrix hamiltonian_nucleus
 call reduce_hamiltonian_sca(m_ham,n_ham,hamiltonian_nucleus)


 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus_buffer_sca


!=========================================================================
subroutine setup_hartree_ri_buffer_sca(print_matrix_,nbf,m_ham,n_ham,p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: hartree_ij(m_ham,n_ham)
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

 call broadcast_hamiltonian_sca(m_ham,n_ham,p_matrix(:,:,1))
 if( nspin == 2 ) then
   call broadcast_hamiltonian_sca(m_ham,n_ham,p_matrix(:,:,2))
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
 call reduce_hamiltonian_sca(m_ham,n_ham,hartree_ij)

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
subroutine setup_exchange_ri_buffer_sca(nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 integer,intent(in)   :: nbf,m_ham,n_ham
 integer,intent(in)   :: nstate,m_c,n_c
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(in)  :: c_matrix(m_c,n_c,nspin)
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: exchange_ij(m_ham,n_ham,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,ispin,istate
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
 real(dp)             :: c_matrix_i(nbf)
 integer              :: iglobal,ilocal,jlocal
!=====

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity: SCALAPACK buffer'
 call start_clock(timing_exchange)


 allocate(tmp(nauxil_3center,nbf))

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
   call reduce_hamiltonian_sca(m_ham,n_ham,exchange_ij(:,:,ispin))

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
subroutine setup_exchange_longrange_ri_buffer_sca(nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 integer,intent(in)   :: nbf,m_ham,n_ham
 integer,intent(in)   :: nstate,m_c,n_c
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(in)  :: c_matrix(m_c,n_c,nspin)
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: exchange_ij(m_ham,n_ham,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,ispin,istate
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
 real(dp)             :: c_matrix_i(nbf)
 integer              :: iglobal,ilocal,jlocal
!=====


 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity: SCALAPACK buffer'
 call start_clock(timing_exchange)


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
   call reduce_hamiltonian_sca(m_ham,n_ham,exchange_ij(:,:,ispin))

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
subroutine dft_exc_vxc_buffer_sca(basis,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,vxc_ij,exc_xc)
 use m_inputparam
 use m_basis_set
 use m_dft_grid
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 integer,intent(in)         :: m_c,n_c
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(m_c,n_c,nspin)
 real(dp),intent(in)        :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out)       :: vxc_ij(m_ham,n_ham,nspin)
 real(dp),intent(out)       :: exc_xc
!=====

 real(dp),parameter :: TOL_RHO=1.0e-10_dp
 integer  :: idft_xc
 integer  :: igrid,ibf,jbf,ispin
 real(dp) :: normalization(nspin)
 real(dp) :: weight

 real(dp) :: basis_function_r(basis%nbf)
 real(dp) :: basis_function_gradr(3,basis%nbf)
 real(dp) :: basis_function_laplr(3,basis%nbf)

 real(dp) :: rhor(nspin,ngrid)
 real(dp) :: grad_rhor(3,nspin,ngrid)
 real(dp) :: sigma(2*nspin-1)
 real(dp) :: tau(nspin),lapl_rhor(nspin)
 real(dp) :: vxc_libxc(nspin)
 real(dp) :: vxc_dummy(nspin)
 real(dp) :: exc_libxc(1)
 real(dp) :: vsigma(2*nspin-1)
 real(dp) :: vlapl_rho(nspin),vtau(nspin)
 real(dp) :: vxc_av(nspin)
 real(dp) :: dedd_r(nspin)
 real(dp) :: dedgd_r(3,nspin)
 real(dp) :: gradtmp(basis%nbf)
!=====

! if( nspin/=1 ) call die('DFT XC potential: SCALAPACK buffer not implemented for spin unrestricted')

 exc_xc = 0.0_dp
 vxc_ij(:,:,:) = 0.0_dp
 if( ndft_xc == 0 ) return

 call start_clock(timing_dft)


#ifdef HAVE_LIBXC

 write(stdout,*) 'Calculate DFT XC potential: SCALAPACK buffer'
 

 normalization(:)=0.0_dp


 do ispin=1,nspin
   !
   ! Buffer constains the c_matrix for a spin channel ispin
   buffer(:,:) = 0.0_dp
   call broadcast_hamiltonian_sca(m_c,n_c,c_matrix(:,:,ispin))


   do igrid=1,ngrid

     weight = w_grid(igrid)

     !
     ! Get all the functions at point r
     call get_basis_functions_r(basis,igrid,basis_function_r)
     !
     ! Calculate the density at point r for spin ispin
     call calc_density_r(1,basis%nbf,nstate,occupation(:,ispin),buffer,basis_function_r,rhor(ispin,igrid))

     ! Skip all the rest if the density is too small
     if( rhor(ispin,igrid) < TOL_RHO ) cycle

     if( dft_xc_needs_gradient ) then
       call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
     endif

     !
     ! Normalization
     normalization(ispin) = normalization(ispin) + rhor(ispin,igrid) * weight


     if( dft_xc_needs_gradient ) then 
       call calc_density_gradr(1,basis%nbf,nstate,occupation(:,ispin),buffer,basis_function_r,basis_function_gradr,grad_rhor(:,ispin,igrid))
     endif

   enddo
 enddo


 do ispin=1,nspin

   !
   ! buffer now contains the vxc_ij
   buffer(:,:) = 0.0_dp

   do igrid=1,ngrid

     ! Skip if the density is too small
     if( rhor(ispin,igrid) < TOL_RHO ) cycle

     weight = w_grid(igrid)

     if( dft_xc_needs_gradient ) then
       sigma(1) = SUM( grad_rhor(:,1,igrid)**2 )
       if(nspin==2) then
         sigma(2) = SUM( grad_rhor(:,1,igrid) * grad_rhor(:,2,igrid) )
         sigma(3) = SUM( grad_rhor(:,2,igrid)**2 )
       endif
     endif

     !
     ! LIBXC calls
     !
     dedd_r(:)    = 0.0_dp
     dedgd_r(:,:) = 0.0_dp

     do idft_xc=1,ndft_xc

       select case(xc_f90_info_family(calc_type%xc_info(idft_xc)))

       case(XC_FAMILY_LDA)
         if( dft_xc_type(idft_xc) < 1000 ) then 
           call xc_f90_lda_exc_vxc(calc_type%xc_func(idft_xc),1,rhor(1,igrid),exc_libxc(1),vxc_libxc(1))
         else
           call my_lda_exc_vxc(nspin,dft_xc_type(idft_xc),rhor(:,igrid),exc_libxc(1),vxc_libxc)
         endif

       case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         if( dft_xc_type(idft_xc) < 2000 ) then 
           !
           ! Remove too small densities to stabilize the computation
           ! especially useful for Becke88
           if( ANY( rhor(:,igrid) > 1.0e-9_dp ) ) then
             call xc_f90_gga_exc_vxc(calc_type%xc_func(idft_xc),1,rhor(1,igrid),sigma(1),exc_libxc(1),vxc_libxc(1),vsigma(1))
           else
             exc_libxc(:)     = 0.0_dp
             vxc_libxc(:)     = 0.0_dp
             vsigma(:)        = 0.0_dp
           endif
         else
           call my_gga_exc_vxc_hjs(gamma_hybrid,rhor(1,igrid),sigma(1),exc_libxc(1),vxc_libxc(1),vsigma(1))
         endif

       case default
         call die('functional is not LDA nor GGA nor hybrid nor meta-GGA')
       end select

       exc_xc = exc_xc + weight * exc_libxc(1) * rhor(ispin,igrid) * dft_xc_coef(idft_xc)

       dedd_r(:) = dedd_r(:) + vxc_libxc(:) * dft_xc_coef(idft_xc)

       !
       ! Set up divergence term if needed (GGA case)
       !
       if( xc_f90_info_family(calc_type%xc_info(idft_xc)) == XC_FAMILY_GGA &
          .OR. xc_f90_info_family(calc_type%xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) then
         if(nspin==1) then

           dedgd_r(:,1) = dedgd_r(:,1) + 2.0_dp * vsigma(1) * grad_rhor(:,1,igrid) * dft_xc_coef(idft_xc)

         else

           dedgd_r(:,1) = dedgd_r(:,1) + ( 2.0_dp * vsigma(1) * grad_rhor(:,1,igrid) & 
                                                  + vsigma(2) * grad_rhor(:,2,igrid) ) * dft_xc_coef(idft_xc)

           dedgd_r(:,2) = dedgd_r(:,2) + ( 2.0_dp * vsigma(3) * grad_rhor(:,2,igrid) &
                                                  + vsigma(2) * grad_rhor(:,1,igrid) ) * dft_xc_coef(idft_xc)
         endif

       endif


     enddo ! loop on the XC functional


     !
     ! Get all the functions at point r
     call get_basis_functions_r(basis,igrid,basis_function_r)
     if( dft_xc_needs_gradient ) then
       call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
     endif

     !
     ! Eventually set up the vxc term
     !
     if( .NOT. dft_xc_needs_gradient ) then 
       ! LDA
       call DSYR('L',basis%nbf,weight*dedd_r(ispin),basis_function_r,1,buffer,basis%nbf)

     else 
       ! GGA
       gradtmp(:) = MATMUL( dedgd_r(:,ispin) , basis_function_gradr(:,:) )
       call DSYR('L',basis%nbf,weight*dedd_r(ispin),basis_function_r,1,buffer,basis%nbf)
       call DSYR2('L',basis%nbf,weight,basis_function_r,1,gradtmp,1,buffer,basis%nbf)
     endif

   enddo ! loop on the grid point

   ! Symmetrize now
   do jbf=1,basis%nbf
     do ibf=jbf+1,basis%nbf
       buffer(jbf,ibf) = buffer(ibf,jbf)
     enddo
   enddo

   call reduce_hamiltonian_sca(m_ham,n_ham,vxc_ij(:,:,ispin))


 enddo

 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(exc_xc)

! !
! ! Destroy operations
! do idft_xc=1,ndft_xc
!   call xc_f90_func_end(calc_type%xc_func(idft_xc))
! enddo

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
 use m_basis_set
 use m_dft_grid
 use m_eri_calculate
 use m_tools,only: matrix_trace
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
 real(dp)             :: vxc,exc,excr
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhgau(m_ham,n_ham)
 integer              :: iatom,igau,ngau
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
