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
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min
 use m_basis_set
 use m_eri_calculate

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
       exchange_ij_cmplx(ibf,jbf,ispin) = CONJG( exchange_ij_cmplx(jbf,ibf,ispin) )
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
subroutine setup_exchange_versatile_ri_cmplx(occupation,c_matrix_cmplx,p_matrix_cmplx,exchange_ij_cmplx,eexchange)
 implicit none
 real(dp),intent(in)     :: occupation(:,:)
 complex(dp),intent(in)  :: c_matrix_cmplx(:,:,:)
 complex(dp),intent(in)  :: p_matrix_cmplx(:,:,:)
 complex(dp),intent(out) :: exchange_ij_cmplx(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer                 :: nauxil_local,npair_local
 integer                 :: nbf,nstate
 integer                 :: ibf,jbf,ispin,istate
 integer                 :: nocc
 complex(dp),allocatable :: tmp_cmplx(:,:)
 integer                 :: ipair,ipair_local
 integer                 :: ibf_auxil_first,nbf_auxil_core
!=====

 call start_clock(timing_tddft_exchange)

 write(stdout,*) 'Calculate Complex Exchange term with Resolution-of-Identity: versatile version'

 exchange_ij_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )

 nbf         = SIZE(exchange_ij_cmplx,DIM=1)
 nstate      = SIZE(occupation(:,:),DIM=1)

 ! Find highest occupied state
 nocc = 0
 do ispin=1,nspin
   do istate=1,nstate
     if( ABS(occupation(istate,ispin)) < completely_empty )  cycle
     nocc = MAX(nocc,istate)
   enddo
 enddo

 nauxil_local = SIZE(eri_3center,DIM=1)
 npair_local  = SIZE(eri_3center,DIM=2)

 if( npcol_3center > 1 ) then
   ! Prepare a sub-division of tasks after the reduction DGSUM2D
   nbf_auxil_core  = CEILING( REAL(nauxil_local,dp) / REAL(npcol_3center,dp) )
   ibf_auxil_first = ipcol_3center * nbf_auxil_core + 1
   ! Be careful when a core has fewer data or no data at all
   if( ibf_auxil_first + nbf_auxil_core - 1 > nauxil_local ) nbf_auxil_core = nauxil_local - ibf_auxil_first + 1
   if( ibf_auxil_first > nauxil_local ) nbf_auxil_core = 0

   allocate(tmp_cmplx(nauxil_local,nbf))

   do ispin=1,nspin

     do istate=1,nocc

       if( ABS(occupation(istate,ispin)) < completely_empty ) cycle

       tmp_cmplx(:,:) = ( 0.0_dp , 0.0_dp )
       !$OMP PARALLEL PRIVATE(ibf,jbf,ipair)
       !$OMP DO REDUCTION(+:tmp_cmplx)
       do ipair_local=1,npair_local
         ipair = INDXL2G(ipair_local,block_col,ipcol_3center,first_col,npcol_3center)
         ibf = index_basis(1,ipair)
         jbf = index_basis(2,ipair)
         tmp_cmplx(:,ibf) = tmp_cmplx(:,ibf) + c_matrix_cmplx(jbf,istate,ispin) * eri_3center(:,ipair_local)
         if( ibf /= jbf) &
           tmp_cmplx(:,jbf) = tmp_cmplx(:,jbf) + c_matrix_cmplx(ibf,istate,ispin) * eri_3center(:,ipair_local)
       enddo
       !$OMP END DO
       !$OMP END PARALLEL
       call ZGSUM2D(cntxt_3center,'R',' ',nauxil_local,nbf,tmp_cmplx,nauxil_local,-1,-1)

       ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
       !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
       ! C = A^T * A + C
       if( nbf_auxil_core > 0 ) &
         call ZHERK('L','C',nbf,nbf_auxil_core,-occupation(istate,ispin)/spin_fact,tmp_cmplx(ibf_auxil_first,1),nauxil_local,1.0_dp,exchange_ij_cmplx(:,:,ispin),nbf)

     enddo

   enddo
   deallocate(tmp_cmplx)

 else
   ! npcol_row is equal to 1 => no parallelization over basis pairs

   allocate(tmp_cmplx(nauxil_local,nbf))

   do ispin=1,nspin

     do istate=1,nocc
       if( MODULO( istate - 1 , nproc_ortho ) /= rank_ortho ) cycle

       if( ABS(occupation(istate,ispin)) < completely_empty ) cycle

       tmp_cmplx(:,:) = ( 0.0_dp, 0.0_dp )
       !$OMP PARALLEL PRIVATE(ibf,jbf,ipair)
       !$OMP DO REDUCTION(+:tmp_cmplx)
       do ipair=1,npair
         ibf = index_basis(1,ipair)
         jbf = index_basis(2,ipair)
         tmp_cmplx(:,ibf) = tmp_cmplx(:,ibf) + c_matrix_cmplx(jbf,istate,ispin) * eri_3center(:,ipair)
         if( ibf /= jbf) &
           tmp_cmplx(:,jbf) = tmp_cmplx(:,jbf) + c_matrix_cmplx(ibf,istate,ispin) * eri_3center(:,ipair)
       enddo
       !$OMP END DO
       !$OMP END PARALLEL

       ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
       !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
       ! C = A^T * A + C
       call ZHERK('L','C',nbf,nauxil_local,-occupation(istate,ispin)/spin_fact,tmp_cmplx,nauxil_local,1.0_dp,exchange_ij_cmplx(:,:,ispin),nbf)

     enddo
   enddo
   deallocate(tmp_cmplx)

 endif

 !
 ! Need to make exchange_ij Hermitian
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij_cmplx(ibf,jbf,ispin) = CONJG(exchange_ij_cmplx(jbf,ibf,ispin))
     enddo
   enddo
   exchange_ij_cmplx(:,:,ispin)=CONJG(exchange_ij_cmplx(:,:,ispin))
 enddo
 call xsum_world(exchange_ij_cmplx)

 eexchange = REAL( 0.5_dp * SUM( exchange_ij_cmplx(:,:,:) * CONJG(p_matrix_cmplx(:,:,:)) ) ,dp)

 call stop_clock(timing_tddft_exchange)

end subroutine setup_exchange_versatile_ri_cmplx

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
 complex(dp),allocatable :: c_matrix_sqrtocc(:,:)
!=====

 call start_clock(timing_density_matrix_cmplx)

 nbf    = SIZE(c_matrix_cmplx(:,:,:),DIM=1)
 nocc   = SIZE(c_matrix_cmplx(:,:,:),DIM=2)
 nstate = SIZE(occupation(:,:),DIM=1)

 if( ANY( occupation(:,:) < 0.0_dp ) ) call die('setup_density_matrix_cmplx: negative occupation number should not happen here.')

 allocate(c_matrix_sqrtocc(nbf,nocc))

 p_matrix_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )
 do ispin=1,nspin

   do istate=1,nocc
     c_matrix_sqrtocc(:,istate) = c_matrix_cmplx(:,istate,ispin) * SQRT(occupation(istate,ispin))
   enddo
   call ZHERK('L','N',nbf,nocc,1.0d0,c_matrix_sqrtocc,nbf,0.0d0,p_matrix_cmplx(1,1,ispin),nbf)
 

   ! Hermitianize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix_cmplx(jbf,ibf,ispin) = CONJG( p_matrix_cmplx(ibf,jbf,ispin) )
     enddo
   enddo

 enddo

 deallocate(c_matrix_sqrtocc)

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

 write(stdout,*) 'Calculate DFT XC potential'
 if( batch_size /= 1 ) write(stdout,*) 'Using batches of size',batch_size
 

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
     allocate(dedgd_r_batch(3,nr,nspin))
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
   ! Eventually set up the vxc term
   !

   allocate(tmp_batch(basis%nbf,nr))
   do ispin=1,nspin
     !
     ! LDA and GGA
     do ir=1,nr
       tmp_batch(:,ir) = weight_batch(ir) * dedd_r_batch(ispin,ir) * basis_function_r_batch(:,ir) * 0.50_dp
     enddo
     !
     ! GGA-only
     if( dft_xc_needs_gradient ) then
       do ir=1,nr
         tmp_batch(:,ir) = tmp_batch(:,ir) &
                          +  MATMUL( basis_function_gradr_batch(:,ir,:) , dedgd_r_batch(:,ir,ispin) * weight_batch(ir) )
       enddo
     endif

     call DSYR2K('L','N',basis%nbf,nr,1.0d0,basis_function_r_batch,basis%nbf,tmp_batch,basis%nbf,1.0d0,vxc_ij(:,:,ispin),basis%nbf)
   enddo
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

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
 write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

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

!=========================================================================
subroutine calc_density_in_disc_cmplx_dft_grid(batch_size,basis,occupation,c_matrix_cmplx,num,time_cur)
 use m_inputparam
 use m_dft_grid
 use m_atoms
 implicit none

 integer,intent(in)         :: batch_size
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(:,:)
 complex(dp),intent(in)     :: c_matrix_cmplx(:,:,:)
 integer,intent(in)         :: num
 real(dp),intent(in)        :: time_cur
!=====
 real(dp)             :: length
 integer              :: nstate,ndisc
 integer              :: ibf,jbf,ispin
 integer              :: idft_xc
 integer              :: igrid_start,igrid_end,ir,nr
 character(len=200)   :: file_name(2)
 real(dp)             :: z_min,z_max
 integer              :: file_out(2),igrid
 !vectors in the plane
 real(dp)             :: vec_r(3) 
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: tmp_batch(:,:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: rhor_batch(:,:)
 real(dp)             :: dz_disc
 real(dp),allocatable :: charge_disc(:,:)
 real(dp),allocatable :: charge_out(:,:)
 integer(dp),allocatable :: count_z_section(:,:)
 integer              :: idisc,i_max_atom,nocc
 logical              :: file_exists
 integer              :: imanual
!=====

 call start_clock(timing_calc_dens_disc)

 nocc = SIZE(c_matrix_cmplx(:,:,:),DIM=2)

 if( excit_type%form==EXCIT_PROJECTILE ) then
   i_max_atom=natom-nprojectile
 else
   i_max_atom=natom
 endif

 inquire(file='manual_disc_dft_grid',exist=file_exists)
 if(file_exists) then
   open(newunit=imanual,file='manual_disc_dft_grid',status='old')
   read(imanual,*) ndisc
   read(imanual,*) length
   close(imanual)
 else
   ndisc=100
   length=10.0_dp
   call issue_warning('calc_density_in_disc_cmplx_dft_grid: manual file was not found')
 endif

 z_min =MIN(MINVAL( xatom(3,1:i_max_atom) ),MINVAL( xbasis(3,:) )) - length
 z_max =MAX(MAXVAL( xatom(3,1:i_max_atom) ),MAXVAL( xbasis(3,:) )) + length

 nstate = SIZE(occupation,DIM=1)

 write(stdout,*) 'Calculate electronic density in discs'

 !
 ! Loop over batches of grid points
 !

 allocate(charge_disc(ndisc,nspin))
 allocate(charge_out(2,nspin))
 allocate(count_z_section(ndisc,nspin))
 charge_disc(:,:) = 0.0_dp
 charge_out(:,:) = 0.0_dp
 count_z_section(:,:) = 0

 dz_disc=(z_max-z_min)/ndisc

 do igrid_start=1,ngrid,batch_size
   igrid_end = MIN(ngrid,igrid_start+batch_size-1)
   nr = igrid_end - igrid_start + 1

   allocate(weight_batch(nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(rhor_batch(nspin,nr))

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   call get_basis_functions_r_batch(basis,igrid_start,nr,basis_function_r_batch)

   call calc_density_r_batch_cmplx(nspin,basis%nbf,nstate,nocc,nr,occupation,c_matrix_cmplx,basis_function_r_batch,rhor_batch)

   do ir=1,nr
     igrid = igrid_start + ir - 1
     vec_r=rr_grid(1:3,igrid)
     do ispin=1,nspin
       idisc = INT((vec_r(3)-z_min)/dz_disc) + 1
       if( idisc > 0 .AND. idisc <= ndisc .AND. (vec_r(1)**2+vec_r(2)**2)**0.5_dp <= r_disc ) then
         charge_disc(idisc,ispin)=charge_disc(idisc,ispin)+rhor_batch(ispin,ir) * weight_batch(ir)
         count_z_section(idisc,ispin)=count_z_section(idisc,ispin)+1
       end if
       if( idisc <= 0 ) then
         charge_out(1,ispin)=charge_out(1,ispin)+rhor_batch(ispin,ir) * weight_batch(ir)
       end if
       if( idisc > ndisc ) then
         charge_out(2,ispin)=charge_out(2,ispin)+rhor_batch(ispin,ir) * weight_batch(ir)
       end if
     enddo ! loop on the ispin
   enddo ! loop on ir

   deallocate(weight_batch)
   deallocate(basis_function_r_batch)
   deallocate(rhor_batch)

 enddo ! loop on the batches
 call xsum_grid(charge_disc(:,:))

 if( is_iomaster ) then

   do ispin=1,nspin
     write(file_name(ispin),'(a,i4.4,a,i1,a,i3.3,f0.3,a)') 'disc_dens_',num, "_s_",ispin,"_r_",INT(r_disc),r_disc-INT(r_disc),".dat"
     open(newunit=file_out(ispin),file=file_name(ispin))
   enddo

   do ispin=1,nspin
     write(file_out(ispin),'(a,F12.6,a,3F12.6)') '# Time: ',time_cur, '  Projectile position (A): ',xatom(:,natom+nghost)*bohr_A
     do idisc=1,ndisc  
       write(file_out(ispin),'(F16.4,F19.10,i6)') (z_min+idisc*dz_disc)*bohr_A,charge_disc(idisc,ispin),count_z_section(idisc,ispin)
     end do
     close(file_out(ispin))
   end do
   write(stdout,'(A,2F12.6)') "Charge out of region, left and right:",charge_out(:,1)

 end if
 !
 ! Sum up the contributions from all procs only if needed

 deallocate(charge_disc)

 call stop_clock(timing_calc_dens_disc)

end subroutine calc_density_in_disc_cmplx_dft_grid

end module m_hamiltonian_cmplx
!=========================================================================


