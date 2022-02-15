!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the routines to calculate the polarizability within RPA, TDDFT or BSE
!
!=========================================================================
#include "molgw.h"
subroutine polarizability(enforce_rpa,calculate_w,basis,nstate,occupation,energy,c_matrix,en_rpa,en_gw,wpol_out)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_scalapack
 use m_cart_to_pure
 use m_block_diago
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_spectra
 implicit none

 logical,intent(in)                    :: enforce_rpa,calculate_w
 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)                  :: en_rpa,en_gw
 type(spectral_function),intent(inout) :: wpol_out
!=====
 type(spectral_function)   :: wpol_static
 logical                   :: is_bse
 integer                   :: nmat,nexc
 real(dp)                  :: alpha_local
 real(dp),allocatable      :: amb_diag_rpa(:)
 real(dp),allocatable      :: amb_matrix(:,:),apb_matrix(:,:)
 real(dp),allocatable      :: xpy_matrix(:,:),xmy_matrix(:,:)
 real(dp),allocatable      :: eigenvalue(:)
 real(dp)                  :: energy_qp(nstate,nspin)
 logical                   :: is_tddft,is_rpa
 logical                   :: has_manual_tdhf
 integer                   :: reading_status
 integer                   :: tdhffile
 integer                   :: m_apb,n_apb,m_x,n_x
! SCALAPACK variables
 integer                   :: desc_apb(NDEL),desc_x(NDEL)
 integer                   :: info
!=====

 call start_clock(timing_pola)
 en_rpa = 0.0_dp
 en_gw  = 0.0_dp

 write(stdout,'(/,a)') ' Calculating the polarizability'
 if(is_triplet) then
   write(stdout,'(a)') ' Triplet final state'
 else
   write(stdout,'(a)') ' Singlet final state'
 endif

 if( has_auxil_basis ) &
   call calculate_eri_3center_eigen(c_matrix,ncore_W+1,nvirtual_W-1,ncore_W+1,nvirtual_W-1,timing=timing_aomo_pola)

 ! Set up all the switches to be able to treat
 ! GW, BSE, TDHF, TDDFT (semilocal or hybrid)

 !
 ! Set up flag is_tddft and is_bse
 is_tddft = calc_type%is_td .AND. calc_type%is_dft .AND. .NOT. enforce_rpa
 is_bse   = calc_type%is_bse .AND. .NOT. enforce_rpa

 !
 ! Set up exchange content alpha_local
 ! manual_tdhf can override anything
 inquire(file='manual_tdhf',exist=has_manual_tdhf)
 if(has_manual_tdhf) then
   open(newunit=tdhffile,file='manual_tdhf',status='old')
   read(tdhffile,*) alpha_local
   close(tdhffile)
   write(msg,'(a,f12.6,3x,f12.6)') 'calculating the TDHF polarizability with alpha ',alpha_local
   call issue_warning(msg)
 else
   if(calc_type%is_td) then        ! TDDFT or TDHF
     alpha_local = alpha_hybrid
   else if(is_bse .AND. .NOT. calc_type%no_bse_kernel) then  ! BSE
     alpha_local = 1.0_dp
   else                  ! RPA or no_bse_kernel
     alpha_local = 0.0_dp
   endif
 endif

 is_rpa = .NOT.(is_tddft) .AND. .NOT.(is_bse) .AND. (ABS(alpha_local)<1.0e-5_dp)

 call start_clock(timing_build_h2p)

 !
 ! Prepare the QP energies
 !
 if( is_bse ) then
   ! Get energy_qp
   call get_energy_qp(nstate,energy,occupation,energy_qp)
 else
   ! For any other type of calculation, just fill energy_qp array with energy
   energy_qp(:,:) = energy(:,:)
 endif

 !
 ! BSE needs the static screening from a previous calculation
 ! It is stored in object wpol_static
 !
 if( is_bse ) then
   call init_spectral_function(nstate,occupation,1,wpol_static)
   call read_spectral_function(wpol_static,reading_status)

   ! If a SCREENED_COULOMB file cannot be found,
   ! then recalculate it from scratch
   if( reading_status /= 0 ) then
     if( .NOT. has_auxil_basis ) then
       call die('polarizability: BSE calculation without having a precalculated SCREENED_COULOMB file is impossible ' &
                // 'unless when using an auxiliary basis')
     endif
     wpol_static%nprodbasis = nauxil_3center
     call static_polarizability(nstate,occupation,energy,wpol_static)
   endif

 endif

 !
 ! Prepare the big matrices (A+B) and (A-B)
 !
 nmat = wpol_out%npole_reso
 !
 ! The distribution of the two matrices have to be the same for A-B and A+B
 ! This is valid also when SCALAPACK is not used!
 m_apb = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
 n_apb = NUMROC(nmat,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desc_apb,nmat,nmat,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,m_apb),info)
 call clean_allocate('A+B',apb_matrix,m_apb,n_apb)
 call clean_allocate('A-B',amb_matrix,m_apb,n_apb)
 allocate(amb_diag_rpa(nmat))


 !
 ! Build the (A+B) and (A-B) matrices in 3 steps
 ! to span all the possible approximations
 ! Only the lower triangle is calculated
 ! the upper part will be filled later by symmetry
 !

 apb_matrix(:,:) = 0.0_dp
 amb_matrix(:,:) = 0.0_dp
 write(stdout,'(/,a)') ' Build the electron-hole hamiltonian'

 if( has_auxil_basis) then

   !
   ! Step 1
   call build_amb_apb_diag_auxil(nmat,nstate,energy_qp,wpol_out,m_apb,n_apb,amb_matrix,apb_matrix,amb_diag_rpa)

#if defined(HAVE_SCALAPACK)
   call build_apb_hartree_auxil_scalapack(desc_apb,wpol_out,m_apb,n_apb,apb_matrix)
#else
   call build_apb_hartree_auxil(desc_apb,wpol_out,m_apb,n_apb,apb_matrix)
#endif

   call get_rpa_correlation(nmat,m_apb,n_apb,amb_matrix,apb_matrix,en_rpa)



   !
   ! Step 2
   if(is_tddft) call build_apb_tddft(nmat,nstate,basis,c_matrix,occupation,wpol_out,m_apb,n_apb,apb_matrix)

   !
   ! Step 3
   if(alpha_local > 1.0e-6_dp) then
     call build_amb_apb_screened_exchange_auxil(alpha_local,desc_apb,wpol_out,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
   endif

   if( is_bse ) then
     call destroy_spectral_function(wpol_static)
   endif


 else

   !
   ! Step 1
   call build_amb_apb_common(nmat,basis%nbf,nstate,c_matrix,energy_qp,wpol_out,alpha_local, &
                             m_apb,n_apb,amb_matrix,apb_matrix,amb_diag_rpa,en_rpa)

   !
   ! Step 2
   if(is_tddft) call build_apb_tddft(nmat,nstate,basis,c_matrix,occupation,wpol_out,m_apb,n_apb,apb_matrix)


   !
   ! Step 3
   if( is_bse ) then
     call build_amb_apb_bse(wpol_out,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
     call destroy_spectral_function(wpol_static)
   endif

 endif



 ! Warning if Tamm-Dancoff flag is on
 if(is_tda) then
   msg='Tamm-Dancoff approximation is switched on'
   call issue_warning(msg)
   ! Tamm-Dancoff approximation consists in setting B matrix to zero
   ! Then A+B = A-B = A
   apb_matrix(:,:) = 0.5_dp * ( apb_matrix(:,:) + amb_matrix(:,:) )
   amb_matrix(:,:) = apb_matrix(:,:)
 endif
 ! Construction done!
 !if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_build_h2p)

 if( is_rpa .AND. .NOT. is_tda ) call clean_deallocate('A-B',amb_matrix)


 !
 ! Prepare the second dimension of xpy_matrix and xmy_matrix
 nexc = nexcitation
 if( nexc == 0 ) nexc = nmat

 allocate(eigenvalue(nexc))

 ! Allocate (X+Y)
 ! Allocate (X-Y) only if actually needed
 m_x = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
 n_x = NUMROC(nexc,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desc_x,nmat,nexc,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,m_x),info)

 call clean_allocate('X+Y',xpy_matrix,m_x,n_x)
 if( .NOT. is_rpa .OR. is_tda ) &
   call clean_allocate('X-Y',xmy_matrix,m_x,n_x)

 !
 ! Diago using the 4 block structure and the symmetry of each block
 ! With or Without SCALAPACK
 !
 if( .NOT. is_rpa .OR. is_tda ) then
   if( nexcitation == 0 ) then
     ! The following call works with AND without SCALAPACK
     call diago_4blocks_chol(amb_matrix,apb_matrix,desc_apb,eigenvalue,xpy_matrix,xmy_matrix,desc_x)

   else ! Partial diagonalization with Davidson
     ! The following call works with AND without SCALAPACK
     call diago_4blocks_davidson(toldav,nstep_dav,amb_diag_rpa,amb_matrix,apb_matrix,desc_apb, &
                                 eigenvalue,xpy_matrix,xmy_matrix,desc_x)
   endif
 else
   ! The following call works with AND without SCALAPACK
   call diago_4blocks_rpa_sca(amb_diag_rpa,apb_matrix,desc_apb,eigenvalue,xpy_matrix,desc_x)
 endif


 ! Deallocate the non-necessary matrices
 deallocate(amb_diag_rpa)
 write(stdout,*) 'Deallocate (A+B) and possibly (A-B)'
 call clean_deallocate('A+B',apb_matrix)
 !
 ! (A-B) may have been already deallocated earlier in the case of RPA
 ! Relax: this is indeed tolerated by clean_deallocate
 call clean_deallocate('A-B',amb_matrix)


 !
 ! Second part of the RPA correlation energy: sum over positive eigenvalues
 en_rpa = en_rpa + 0.50_dp * SUM( ABS(eigenvalue(:)) )
 if( is_rpa ) then
   write(stdout,'(/,a)') ' Calculate the RPA energy using the Tamm-Dancoff decomposition'
   write(stdout,'(a)')   ' Eq. (9) from J. Chem. Phys. 132, 234114 (2010)'
   write(stdout,'(/,a,f16.10)') ' RPA correlation energy (Ha): ',en_rpa
 endif

 write(stdout,'(/,a,f12.6)') ' Lowest neutral excitation energy (eV):',MINVAL(ABS(eigenvalue(1:nexc)))*Ha_eV

 !if( has_auxil_basis ) call calculate_eri_3center_eigen(c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1,timing=timing_aomo_pola)

 !
 ! Calculate the optical sprectrum
 ! and the dynamic dipole tensor
 !
 if( calc_type%is_td .OR. is_bse ) then
   call optical_spectrum(basis,occupation,c_matrix,wpol_out,xpy_matrix,xmy_matrix,eigenvalue)
   select case(TRIM(lower(stopping)))
   case('spherical')
     call stopping_power(basis,c_matrix,wpol_out,xpy_matrix,eigenvalue)
   case('3d')
     call stopping_power_3d(basis,c_matrix,wpol_out,xpy_matrix,desc_x,eigenvalue)
   end select
 endif

 !
 ! Now only the sum ( X + Y ) is needed in fact.
 ! Free ( X - Y ) at once.
 call clean_deallocate('X-Y',xmy_matrix)

 !
 ! Calculate Wp = v * chi * v    if necessary
 ! and then write it down on file
 !
 if( print_w_ .OR. calculate_w ) then
   if( has_auxil_basis) then
     call chi_to_sqrtvchisqrtv_auxil(desc_x,m_x,n_x,xpy_matrix,eigenvalue,wpol_out,en_gw)
     ! This following coding of the Galitskii-Migdal correlation energy is only working with
     ! an auxiliary basis
     if( is_rpa ) then
       write(stdout,'(/,1x,a)')        'Correlation energy in the Galitskii-Migdal formula'
       write(stdout,'(1x,a,f19.10,/)') '                        1/2 Tr[ Sig_c * G ] (Ha): ',en_gw
     endif

   else
     call chi_to_vchiv(basis%nbf,nstate,c_matrix,xpy_matrix,eigenvalue,wpol_out)
   endif


   ! If requested write the spectral function on file
   if( print_w_ ) call write_spectral_function(wpol_out)

 else
   call destroy_spectral_function(wpol_out)
 endif


 write(stdout,*) 'Deallocate eigenvector array'
 call clean_deallocate('X+Y',xpy_matrix)

 if(has_auxil_basis) call destroy_eri_3center_eigen()

 if(ALLOCATED(eigenvalue)) deallocate(eigenvalue)

 call stop_clock(timing_pola)


end subroutine polarizability


!=========================================================================
subroutine polarizability_onering(basis,nstate,energy,c_matrix,vchi0v)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_scalapack
 use m_cart_to_pure
 use m_block_diago
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: energy(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(inout) :: vchi0v
!=====
 integer :: t_jb
 integer :: jstate,bstate,jbspin
!=====

 call allocate_spectral_function(nauxil_3center,vchi0v)

 call calculate_eri_3center_eigen(c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1,timing=timing_aomo_pola)


 do t_jb=1,vchi0v%npole_reso
   jstate = vchi0v%transition_table(1,t_jb)
   bstate = vchi0v%transition_table(2,t_jb)
   jbspin = vchi0v%transition_table(3,t_jb)

   vchi0v%residue_left(:,t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin) * SQRT(spin_fact)
   vchi0v%pole(t_jb)           = energy(bstate,jbspin) - energy(jstate,jbspin)

 enddo

 call destroy_eri_3center_eigen()

end subroutine polarizability_onering


!=========================================================================
subroutine get_energy_qp(nstate,energy,occupation,energy_qp)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_io
 implicit none

 integer,intent(in)                  :: nstate
 real(dp),intent(in)                 :: energy(nstate,nspin)
 real(dp),intent(in)                 :: occupation(nstate,nspin)
 real(dp),intent(out)                :: energy_qp(nstate,nspin)
!=====
 integer  :: reading_status
 integer  :: mspin,mstate
!=====

 ! If the keyword scissor is used in the input file,
 ! then use it and ignore the ENERGY_QP file
 if( ABS(scissor) > 1.0e-5_dp ) then

   call issue_warning('Using a manual scissor to open up the fundamental gap')

   write(stdout,'(a,2(1x,f12.6))') ' Scissor operator with value (eV):',scissor*Ha_eV
   do mspin=1,nspin
     do mstate=1,nstate
       if( occupation(mstate,mspin) > completely_empty/spin_fact ) then
         energy_qp(mstate,mspin) = energy(mstate,mspin)
       else
         energy_qp(mstate,mspin) = energy(mstate,mspin) + scissor
       endif
     enddo
   enddo
   write(stdout,'(/,a)') ' Scissor updated energies'
   do mstate=1,nstate
     write(stdout,'(i5,4(2x,f16.6))') mstate,energy(mstate,:)*Ha_eV,energy_qp(mstate,:)*Ha_eV
   enddo
   write(stdout,*)

 else

   call read_energy_qp(nstate,energy_qp,reading_status)

   select case(reading_status)
   case(0)
     write(stdout,'(a)') ' Reading OK'
   case(1,2)
     write(stdout,'(a,/,a)') ' Something happened during the reading of energy_qp file',' Fill up the QP energies with KS energies'
     energy_qp(:,:) = energy(:,:)
   case default
     call die('reading_status BUG')
   end select

 endif

end subroutine get_energy_qp


!=========================================================================
subroutine chi_to_vchiv(nbf,nstate,c_matrix,xpy_matrix,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none

 integer,intent(in)                    :: nbf,nstate
 real(dp),intent(in)                   :: c_matrix(nbf,nstate,nspin)
 type(spectral_function),intent(inout) :: wpol
 real(dp),intent(in)                   :: xpy_matrix(wpol%npole_reso,wpol%npole_reso)
 real(dp),intent(in)                   :: eigenvalue(wpol%npole_reso)
!=====
 integer                               :: t_jb,jbspin,mpspin
 integer                               :: mstate,pstate,jstate,bstate,mpstate_spin
 integer                               :: kbstate_min
 integer                               :: kbstate_max
 integer                               :: nmat,nprodbasis
 real(dp)                              :: eri_eigen_klij
 real(dp),allocatable                  :: eri_eigenstate_klmin(:,:,:,:)
!=====

 call start_clock(timing_vchiv)

 write(stdout,'(/,a)') ' Build Wp = v * chi * v'
 if(has_auxil_basis) then
   call die('you should not be here')
 endif

 allocate(eri_eigenstate_klmin(nbf,nbf,nbf,nspin))
 ! Set this to zero and then enforce the calculation of the first array of Coulomb integrals
 ! If removed, calculate_eri_4center_eigen might not calculate the first term!
 eri_eigenstate_klmin(:,:,:,:) = 0.0_dp

 nprodbasis = index_prodstate(nvirtual_W-1,nvirtual_W-1) * nspin
 call allocate_spectral_function(nprodbasis,wpol)

 wpol%pole(1:wpol%npole_reso) = eigenvalue(:)

 nmat = wpol%npole_reso

 wpol%residue_left(:,:) = 0.0_dp


 do t_jb=1,nmat
   jstate = wpol%transition_table(1,t_jb)
   bstate = wpol%transition_table(2,t_jb)
   jbspin = wpol%transition_table(3,t_jb)

   kbstate_min = MIN(jstate,bstate)
   kbstate_max = MAX(jstate,bstate)
   call calculate_eri_4center_eigen(c_matrix,kbstate_min,jbspin,eri_eigenstate_klmin)


   ! COLLAPSE is used because nspin is much smaller than # of threads.
   !$OMP PARALLEL
   !$OMP DO PRIVATE(eri_eigen_klij,mpstate_spin) COLLAPSE(2)
   do mpspin=1,nspin
     do pstate=1,nstate
       do mstate = 1,pstate

         ! Unique ordering for mpstate_spin so to please OPENMP
         mpstate_spin = ( mpspin - 1 ) * ( nstate * ( nstate + 1 ) ) / 2 + ( ( pstate - 1 ) * pstate ) / 2 + mstate

         eri_eigen_klij = eri_eigenstate_klmin(kbstate_max,mstate,pstate,mpspin)

         ! Use the symmetry ( k l | i j ) to regroup (kl) and (lk) contributions
         ! and the block structure of eigenvector | X  Y |
         !                                        | Y  X |
         wpol%residue_left(mpstate_spin,:) = wpol%residue_left(mpstate_spin,:) &
                              + eri_eigen_klij * xpy_matrix(t_jb,:)

       enddo
     enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

 enddo

 !$OMP PARALLEL WORKSHARE
 wpol%residue_left(:,:) = wpol%residue_left(:,:) * SQRT(spin_fact)
 !$OMP END PARALLEL WORKSHARE


 if(ALLOCATED(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)

 call stop_clock(timing_vchiv)

end subroutine chi_to_vchiv


!=========================================================================
subroutine chi_to_sqrtvchisqrtv_auxil(desc_x,m_x,n_x,xpy_matrix,eigenvalue,wpol,energy_gm)
 use m_definitions
 use m_warning
 use m_scalapack
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none

 integer,intent(in)                    :: m_x,n_x
 integer,intent(in)                    :: desc_x(NDEL)
 real(dp),intent(inout)                :: xpy_matrix(m_x,n_x)
 type(spectral_function),intent(inout) :: wpol
 real(dp),intent(in)                   :: eigenvalue(wpol%npole_reso)
 real(dp),intent(out)                  :: energy_gm
!=====
 integer                               :: t_jb,t_jb_global,jbspin
 integer                               :: nmat
 integer                               :: jstate,bstate
 real(dp),allocatable                  :: eri_3tmp(:,:)
 real(dp),allocatable                  :: eri_3tmp_sd(:,:)
 real(dp),allocatable                  :: vsqrt_xpy(:,:)
 integer                               :: desc_auxil(NDEL),desc_sd(NDEL)
 integer                               :: mlocal,nlocal
 integer                               :: info
!=====

 call start_clock(timing_vchiv)

 write(stdout,'(/,a)') ' Build v^{1/2} * chi * v^{1/2}'

 call allocate_spectral_function(nauxil_3center,wpol)
 wpol%pole(1:wpol%npole_reso) = eigenvalue(:)

 nmat = wpol%npole_reso

#if !defined(HAVE_SCALAPACK)

 allocate(eri_3tmp(nauxil_3center,nmat))
 do t_jb=1,nmat
   jstate = wpol%transition_table(1,t_jb)
   bstate = wpol%transition_table(2,t_jb)
   jbspin = wpol%transition_table(3,t_jb)
   eri_3tmp(:,t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin)
 enddo

 ! Use the symmetry ( I | k l ) to regroup (kl) and (lk) contributions
 ! and the block structure of eigenvector | X  Y |
 !                                        | Y  X |
 ! => only needs (X+Y)
 !wpol%residue_left(:,:) = MATMUL( eri_3tmp , xpy_matrix(:,:) ) * SQRT(spin_fact)
 call DGEMM('N','N',nauxil_3center,nmat,nmat,DSQRT(spin_fact),eri_3tmp,nauxil_3center, &
                                                              xpy_matrix,nmat, &
                                                        0.0d0,wpol%residue_left,nauxil_3center)

 energy_gm = 0.5_dp * ( SUM( wpol%residue_left(:,:)**2 ) - spin_fact * SUM( eri_3tmp(:,:)**2 ) )
 !
 ! Since wpol%residue_left and eri_3tmp are distributed, we have to sum up
 call auxil%sum(energy_gm)

 deallocate(eri_3tmp)

#else

 call clean_allocate('TMP 3-center integrals',eri_3tmp,nauxil_3center,nmat)
 do t_jb=1,nmat
   jstate = wpol%transition_table(1,t_jb)
   bstate = wpol%transition_table(2,t_jb)
   jbspin = wpol%transition_table(3,t_jb)
   eri_3tmp(:,t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin)
 enddo

 !
 ! Descriptors
 mlocal = NUMROC(nauxil_2center,MB_eri3_mo,iprow_eri3_mo,first_row,nprow_eri3_mo)
 call DESCINIT(desc_auxil,nauxil_2center,nmat,MB_eri3_mo,NB_eri3_mo,first_row,first_col,cntxt_eri3_mo,MAX(1,mlocal),info)

 mlocal = NUMROC(nauxil_2center,block_row,iprow_sd,first_row,nprow_sd)
 nlocal = NUMROC(nmat          ,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desc_sd,nauxil_2center,nmat,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mlocal),info)

 call clean_allocate('TMP 3-center integrals',eri_3tmp_sd,mlocal,nlocal)

 call PDGEMR2D(nauxil_2center,nmat,eri_3tmp,1,1,desc_auxil, &
                                eri_3tmp_sd,1,1,desc_sd,cntxt_sd)

 call clean_deallocate('TMP 3-center integrals',eri_3tmp)

 !
 !   SQRT(spin_fact) * v**1/2 * ( X + Y )
 call clean_allocate('TMP v**1/2 * (X+Y)',vsqrt_xpy,mlocal,nlocal)
 call PDGEMM('N','N',nauxil_2center,nmat,nmat, &
               DSQRT(spin_fact),eri_3tmp_sd,1,1,desc_sd,  &
                                 xpy_matrix,1,1,desc_x,   &
                      0.0_dp,     vsqrt_xpy,1,1,desc_sd)

 call clean_deallocate('TMP 3-center integrals',eri_3tmp_sd)


 call PDGEMR2D(nauxil_2center,nmat,vsqrt_xpy,1,1,desc_sd, &
                            wpol%residue_left,1,1,desc_auxil,cntxt_sd)
 !
 ! Do not forget ortho parallelization direction
 if( ortho%nproc > 1 ) then
   call ortho%bcast(0,wpol%residue_left)
 endif

 call clean_deallocate('TMP v**1/2 * (X+Y)',vsqrt_xpy)


 energy_gm = 0.0_dp
 do t_jb_global=1,nmat
   jstate = wpol%transition_table(1,t_jb_global)
   bstate = wpol%transition_table(2,t_jb_global)
   jbspin = wpol%transition_table(3,t_jb_global)
   energy_gm = energy_gm - SUM( eri_3center_eigen(:,jstate,bstate,jbspin)**2 ) * spin_fact * 0.5_dp
 enddo

 energy_gm = energy_gm + 0.5_dp * ( SUM( wpol%residue_left(:,:)**2 ) )
 call auxil%sum(energy_gm)


#endif



 call stop_clock(timing_vchiv)

end subroutine chi_to_sqrtvchisqrtv_auxil


!=========================================================================
