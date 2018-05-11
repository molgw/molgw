!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian_cmplx
 use m_definitions
 use m_timing
 use m_tddft_variables
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min


contains


!=========================================================================
subroutine setup_exchange_ri_cmplx(nbf,nstate,nocc,occupation,c_matrix_cmplx,p_matrix_cmplx, &
                                   exchange_ij_cmplx,eexchange)
 use m_eri
 implicit none
 integer,intent(in)         :: nbf,nstate,nocc
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(out)       :: eexchange
 complex(dp),intent(in)    :: c_matrix_cmplx(nbf,nocc,nspin)
 complex(dp),intent(in)    :: p_matrix_cmplx(nbf,nbf,nspin)
 complex(dp),intent(out)   :: exchange_ij_cmplx(nbf,nbf,nspin)
!=====
 integer                    :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer                    :: index_ij
 real(dp)                   :: eigval(nbf)
 integer                    :: ipair
 complex(dp),allocatable   :: tmp_cmplx(:,:)
!=====

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'

 call start_clock(timing_tddft_exchange)

 exchange_ij_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )
 allocate(tmp_cmplx(nauxil_3center,nbf))
 do ispin=1,nspin

   do istate=1,nocc
     if( MODULO( istate-1 , nproc_ortho ) /= rank_ortho ) cycle

     tmp_cmplx(:,:) = ( 0.0_dp, 0.0_dp )
     !$OMP PARALLEL PRIVATE(ibf,jbf) 
     !$OMP DO REDUCTION(+:tmp_cmplx)
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp_cmplx(:,ibf) = tmp_cmplx(:,ibf) + c_matrix_cmplx(jbf,istate,ispin) * eri_3center(:,ipair)
       if( ibf /= jbf ) then
          tmp_cmplx(:,jbf) = tmp_cmplx(:,jbf) + c_matrix_cmplx(ibf,istate,ispin) * eri_3center(:,ipair)
       end if
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^H * A + C ; C - exchange_ij(:,:,ispin); A - tmp
     !    exchange_ij_cmplx(:,:,ispin) = exchange_ij_cmplx(:,:,ispin) - & 
     !            MATMUL( TRANSPOSE(tmp_cmplx(:,:)) , CONJG(tmp_cmplx(:,:)) ) * occupation(istate,ispin)/ spin_fact
     call ZHERK('L','C',nbf,nauxil_3center,-occupation(istate,ispin)/spin_fact,tmp_cmplx,nauxil_3center,1.0d0,exchange_ij_cmplx(:,:,ispin),nbf)
       
   enddo
   exchange_ij_cmplx(:,:,ispin)=CONJG(exchange_ij_cmplx(:,:,ispin))
   !
   ! Need to make exchange_ij Hermitian (not just symmetric)
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij_cmplx(ibf,jbf,ispin) = conjg( exchange_ij_cmplx(jbf,ibf,ispin) )
     enddo
   enddo

 enddo ! end of loop do ispin=1,nspin
 deallocate(tmp_cmplx)
 ! This interface should work also for complex exchange_ij_cmplx 
 call xsum_world(exchange_ij_cmplx)
 eexchange = REAL( 0.5_dp * SUM( exchange_ij_cmplx(:,:,:) * CONJG(p_matrix_cmplx(:,:,:)) ),dp)

 call stop_clock(timing_tddft_exchange)

end subroutine setup_exchange_ri_cmplx


!=========================================================================
subroutine setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
 implicit none

 complex(dp),intent(in)  :: c_matrix_cmplx(:,:,:)
 real(dp),intent(in)     :: occupation(:,:)
 complex(dp),intent(out) :: p_matrix_cmplx(:,:,:)
!=====
 integer :: nbf,nstate,nocc
 integer :: ispin,ibf,jbf
 integer :: istate
 complex(dp),allocatable :: c_matrix_tmp(:,:)
!=====

 call start_clock(timing_density_matrix_cmplx)

 nbf    = SIZE(c_matrix_cmplx(:,:,:),DIM=1)
 nocc   = SIZE(c_matrix_cmplx(:,:,:),DIM=2)
 nstate = SIZE(occupation(:,:),DIM=1)

 allocate(c_matrix_tmp(nbf,nstate))

 p_matrix_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )
 do ispin=1,nspin

   do istate=1,nocc
     c_matrix_tmp(:,istate) = c_matrix_cmplx(:,istate,ispin) * SQRT(occupation(istate,ispin))
   enddo
   call ZHERK('L','N',nbf,nocc,1.0d0,c_matrix_tmp,nbf,0.0d0,p_matrix_cmplx(:,:,ispin),nbf)
 
   !do istate=1,nocc
   !  call ZHER('L',nbf,occupation(istate,ispin),c_matrix_cmplx(:,istate,ispin),1,p_matrix_cmplx(:,:,ispin),nbf)
   !enddo
  

   ! Hermitianize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix_cmplx(jbf,ibf,ispin) = CONJG( p_matrix_cmplx(ibf,jbf,ispin) )
     enddo
   enddo

 enddo

 deallocate(c_matrix_tmp)

 call stop_clock(timing_density_matrix_cmplx)


end subroutine setup_density_matrix_cmplx


!=========================================================================
subroutine dft_exc_vxc_batch_cmplx(batch_size,basis,nstate,nocc,occupation,c_matrix_cmplx,vxc_ij,exc_xc)
 use m_inputparam
 use m_basis_set
 use m_dft_grid
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 integer,intent(in)         :: batch_size
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 integer,intent(in)         :: nocc
 real(dp),intent(in)        :: occupation(nstate,nspin)
 complex(dp),intent(in)     :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 real(dp),intent(out)       :: vxc_ij(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: exc_xc
!=====
 real(dp),parameter   :: TOL_RHO=1.0e-9_dp
 integer              :: idft_xc
 integer              :: ibf,jbf,ispin
 integer              :: igrid_start,igrid_end,ir,nr
 real(dp)             :: normalization(nspin)
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: tmp_batch(:,:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: basis_function_gradr_batch(:,:,:)
 real(dp),allocatable :: exc_batch(:)
 real(dp),allocatable :: rhor_batch(:,:)
 real(dp),allocatable :: vrho_batch(:,:)
 real(dp),allocatable :: dedd_r_batch(:,:)
 real(dp),allocatable :: grad_rhor_batch(:,:,:)
 real(dp),allocatable :: sigma_batch(:,:)
 real(dp),allocatable :: dedgd_r_batch(:,:,:)
 real(dp),allocatable :: vsigma_batch(:,:)
!=====

 exc_xc = 0.0_dp
 vxc_ij(:,:,:) = 0.0_dp
 if( ndft_xc == 0 ) return

 call start_clock(timing_tddft_xc)


#ifdef HAVE_LIBXC

! write(stdout,*) 'Calculate DFT XC potential'
! if( batch_size /= 1 ) write(stdout,*) 'Using batches of size',batch_size
 

 normalization(:) = 0.0_dp

 !
 ! Loop over batches of grid points
 !
 do igrid_start=1,ngrid,batch_size
   igrid_end = MIN(ngrid,igrid_start+batch_size-1)
   nr = igrid_end - igrid_start + 1

   allocate(weight_batch(nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(exc_batch(nr))
   allocate(rhor_batch(nspin,nr))
   allocate(vrho_batch(nspin,nr))
   allocate(dedd_r_batch(nspin,nr))

   if( dft_xc_needs_gradient ) then 
     allocate(basis_function_gradr_batch(basis%nbf,nr,3))
     allocate(grad_rhor_batch(nspin,nr,3))
     allocate(dedgd_r_batch(nspin,nr,3))
     allocate(sigma_batch(2*nspin-1,nr))
     allocate(vsigma_batch(2*nspin-1,nr))
   endif

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   call get_basis_functions_r_batch(basis,igrid_start,nr,basis_function_r_batch)
   !
   ! Get the gradient at points r
   if( dft_xc_needs_gradient ) call get_basis_functions_gradr_batch(basis,igrid_start,nr,basis_function_gradr_batch)

   !
   ! Calculate the density at points r for spin up and spin down
   ! Calculate grad rho at points r for spin up and spin down
   if( .NOT. dft_xc_needs_gradient ) then 
     call calc_density_r_batch_cmplx(nspin,basis%nbf,nstate,nocc,nr,occupation,c_matrix_cmplx,basis_function_r_batch,rhor_batch)

   else
     call calc_density_gradr_batch_cmplx(nspin,basis%nbf,nstate,nocc,nr,occupation,c_matrix_cmplx, &
                                   basis_function_r_batch,basis_function_gradr_batch,rhor_batch,grad_rhor_batch)
     do ir=1,nr
       sigma_batch(1,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(1,ir,:) )
       if( nspin == 2 ) then
         sigma_batch(2,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(2,ir,:) )
         sigma_batch(3,ir) = DOT_PRODUCT( grad_rhor_batch(2,ir,:) , grad_rhor_batch(2,ir,:) )
       endif
     enddo
   endif

   ! Normalization
   normalization(:) = normalization(:) + MATMUL( rhor_batch(:,:) , weight_batch(:) )

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

     ! XC energy
     exc_xc = exc_xc + SUM( weight_batch(:) * exc_batch(:) * SUM(rhor_batch(:,:),DIM=1) ) * dft_xc_coef(idft_xc)

     dedd_r_batch(:,:) = dedd_r_batch(:,:) + vrho_batch(:,:) * dft_xc_coef(idft_xc)

     !
     ! Set up divergence term if needed (GGA case)
     !
     if( dft_xc_needs_gradient ) then
       do ir=1,nr
         if(nspin==1) then

           dedgd_r_batch(1,ir,:) = dedgd_r_batch(1,ir,:)  &
                      + 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) * dft_xc_coef(idft_xc)

         else

           dedgd_r_batch(1,ir,:) = dedgd_r_batch(1,ir,:) &
                     + ( 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) &
                                 + vsigma_batch(2,ir) * grad_rhor_batch(2,ir,:) ) * dft_xc_coef(idft_xc) 

           dedgd_r_batch(2,ir,:) = dedgd_r_batch(2,ir,:) &
                     + ( 2.0_dp * vsigma_batch(3,ir) * grad_rhor_batch(2,ir,:) &
                                 + vsigma_batch(2,ir) * grad_rhor_batch(1,ir,:) ) * dft_xc_coef(idft_xc)
         endif

       enddo
     endif

   enddo ! loop on the XC functional



   !
   ! Eventually set up the vxc term
   !

   !
   ! LDA and GGA
   allocate(tmp_batch(basis%nbf,nr))
   do ispin=1,nspin
     forall(ir=1:nr)
       tmp_batch(:,ir) = weight_batch(ir) * dedd_r_batch(ispin,ir) * basis_function_r_batch(:,ir)
     end forall

     call DGEMM('N','T',basis%nbf,basis%nbf,nr,1.0d0,tmp_batch,basis%nbf,basis_function_r_batch,basis%nbf,1.0d0,vxc_ij(:,:,ispin),basis%nbf)
   enddo

   !
   ! GGA-only
   if( dft_xc_needs_gradient ) then 

     do ispin=1,nspin

       do ir=1,nr
         tmp_batch(:,ir) = MATMUL( basis_function_gradr_batch(:,ir,:) , dedgd_r_batch(ispin,ir,:) * weight_batch(ir) )
       enddo

       call DSYR2K('L','N',basis%nbf,nr,1.0d0,basis_function_r_batch,basis%nbf,tmp_batch,basis%nbf,1.0d0,vxc_ij(:,:,ispin),basis%nbf)

     enddo
   endif
   deallocate(tmp_batch)

    

   deallocate(weight_batch)
   deallocate(basis_function_r_batch)
   deallocate(exc_batch)
   deallocate(rhor_batch)
   deallocate(vrho_batch)
   deallocate(dedd_r_batch)
   if( dft_xc_needs_gradient ) then 
     deallocate(basis_function_gradr_batch)
     deallocate(grad_rhor_batch)
     deallocate(sigma_batch)
     deallocate(dedgd_r_batch)
     deallocate(vsigma_batch)
   endif

 enddo ! loop on the batches




 ! Symmetrize now
 do ispin=1,nspin
   do jbf=1,basis%nbf
     do ibf=jbf+1,basis%nbf 
       vxc_ij(jbf,ibf,ispin) = vxc_ij(ibf,jbf,ispin)
     enddo
   enddo
 enddo

 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(vxc_ij)
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

! write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
! write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

 call stop_clock(timing_tddft_xc)

end subroutine dft_exc_vxc_batch_cmplx

!=========================================================================
subroutine static_dipole_fast_cmplx(basis,p_matrix_cmplx,dipole_basis,dipole)
 use m_basis_set
 use m_atoms
 use m_timing
 implicit none

!=====
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: dipole_basis(basis%nbf,basis%nbf,3)
 complex(dp),intent(in)      :: p_matrix_cmplx(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)        :: dipole(3)
!=====
 integer                     :: iatom,idir

 ! Minus sign for electrons
 do idir=1,3
   dipole(idir) = real( -SUM( dipole_basis(:,:,idir) * SUM( p_matrix_cmplx(:,:,:) , DIM=3 ) ),dp)
 enddo

 do iatom=1,natom
   dipole(:) = dipole(:) + zatom(iatom) * xatom(:,iatom)
 enddo

end subroutine static_dipole_fast_cmplx

end module m_hamiltonian_cmplx
!=========================================================================


