!=========================================================================
! This module contains
! the routines to calculate the polarizability within RPA, TDDFT or BSE
! and the corresponding optical spectra
!=========================================================================
module m_timedependent
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi

contains


!=========================================================================
subroutine polarizability(basis,prod_basis,auxil_basis,nstate,occupation,energy,c_matrix,rpa_correlation,wpol_out)
 use m_tools
 use m_basis_set
 use m_spectral_function
 implicit none

 type(basis_set),intent(in)            :: basis,prod_basis,auxil_basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(basis%nbf,nspin)
 real(dp),intent(in)                   :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)                  :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol_out
!=====
 integer                   :: nstate0
 integer                   :: t_ij,t_kl
 type(spectral_function)   :: wpol_static
 integer                   :: nmat
 real(dp)                  :: energy_gm
 real(dp)                  :: alpha_local
 real(prec_td),allocatable :: amb_matrix(:,:),apb_matrix(:,:)
 real(prec_td),allocatable :: a_diag(:)
 real(prec_td),allocatable :: bigx(:,:),bigy(:,:)
 real(dp),allocatable      :: eigenvalue(:)
 real(dp)                  :: energy_qp(basis%nbf,nspin)
 logical                   :: is_tddft
 logical                   :: is_ij
 logical                   :: is_rpa
 logical                   :: has_manual_tdhf
 integer                   :: reading_status
 integer                   :: tdhffile
 integer                   :: m_apb,n_apb,m_x,n_x
! Scalapack variables
 integer                   :: desc_apb(ndel),desc_x(ndel)
!=====

 call start_clock(timing_pola)
 nstate0 = basis%nbf

 write(stdout,'(/,a)') ' Calculating the polarizability'
 if(is_triplet) then
   write(stdout,'(a)') ' Triplet state'
 else
   write(stdout,'(a)') ' Singlet state'
 endif
 
 ! Set up all the switches to be able to treat
 ! GW, BSE, TDHF, TDDFT (semilocal or hybrid)

 !
 ! Set up flag is_rpa
 inquire(file='manual_rpa',exist=is_rpa)
 if(is_rpa) then
   msg='RPA calculation is enforced'
   call issue_warning(msg)
 endif
 is_rpa   = calc_type%is_gw .OR. is_rpa

 ! 
 ! Set up flag is_tddft
 is_tddft = calc_type%is_td .AND. calc_type%is_dft .AND. .NOT. is_rpa

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
   if(is_rpa) then
     alpha_local = 0.0_dp
   else if(is_tddft) then
     alpha_local = alpha_hybrid
   else ! TDHF or BSE case
     alpha_local = 1.0_dp
   endif
 endif


 call start_clock(timing_build_h2p)

 !
 ! Prepare the QP energies
 !
 if( calc_type%is_bse .OR. calc_type%gwmethod==GnWn ) then
   ! Get energy_qp 
   call get_energy_qp(basis%nbf,energy,occupation,energy_qp)
 else
   ! For any other type of calculation, just fill energy_qp array with energy
   energy_qp(:,:) = energy(:,:)
 endif

 ! 
 ! BSE needs the static screening from a previous calculation
 ! It is stored in object wpol_static
 !
 if( calc_type%is_bse ) then
   call read_spectral_function(wpol_static,reading_status)

   ! If a SCREENED_COULOMB file cannot be found,
   ! then recalculate it from scratch
   if( reading_status /= 0 ) then
     call init_spectral_function(basis%nbf,nstate,occupation,wpol_static)
     wpol_static%nprodbasis = auxil_basis%nbf_local
     call static_polarizability(basis,auxil_basis,occupation,energy,wpol_static)
   endif

 endif


 !
 ! Prepare the big matrices (A+B) and (A-B)
 ! 
 nmat = wpol_out%npole_reso_apb
 !
 ! The distribution of the two matrices have to be the same for A-B and A+B
 ! This is valid also when SCALAPACK is not used!
 call init_desc('S',nmat,nmat,desc_apb,m_apb,n_apb)
 call clean_allocate('A+B',apb_matrix,m_apb,n_apb)
 call clean_allocate('A-B',amb_matrix,m_apb,n_apb)

 ! A diagonal is owned by all procs (= no distribution)
 ! wpol_out%npole_reso_spa are the pole not explictely counted in wpol_out%npole_reso_apb
 allocate(a_diag(wpol_out%npole_reso_spa))

 !
 ! Build the (A+B) and (A-B) matrices in 3 steps
 ! to span all the possible approximations
 ! Only the lower triangle is calculated
 ! the upper part will be filled later by symmetry
 !
 apb_matrix(:,:) = 0.0_dp
 amb_matrix(:,:) = 0.0_dp
 write(stdout,'(/,a)') ' Build the electron-hole hamiltonian'
 ! Step 1
 call build_amb_apb_common(nmat,basis%nbf,nstate0,c_matrix,energy_qp,wpol_out,alpha_local,m_apb,n_apb,amb_matrix,apb_matrix,rpa_correlation)

 ! Calculate the diagonal separately: it's needed for the single pole approximation
 if( nvirtual_SPA < nvirtual_W .AND. is_rpa ) & 
     call build_a_diag_common(nmat,basis%nbf,nstate0,c_matrix,energy_qp,wpol_out,a_diag)

 ! Step 2
 if(is_tddft) call build_apb_tddft(nmat,basis,c_matrix,occupation,wpol_out,m_apb,n_apb,apb_matrix)

 ! Step 3
 if(calc_type%is_bse .AND. .NOT. is_rpa) then
   if(.NOT. has_auxil_basis ) then
     call build_amb_apb_bse(basis%nbf,prod_basis,c_matrix,wpol_out,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
   else
     call build_amb_apb_bse_auxil(nmat,basis%nbf,c_matrix,wpol_out,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
   endif
   call destroy_spectral_function(wpol_static)
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

 call stop_clock(timing_build_h2p)


 allocate(eigenvalue(nmat))
 ! bigX, bigY
 call init_desc('S',nmat,nmat,desc_x,m_x,n_x)
 write(stdout,*) 'Allocate eigenvector array'
 call clean_allocate('X',bigx,m_x,n_x)
 call clean_allocate('Y',bigy,m_x,n_x)

 !
 ! Diago using the 4 block structure and the symmetry of each block
 ! With or Without SCALAPACK
 !
#ifndef HAVE_SCALAPACK
 call diago_4blocks_sqrt(nmat,amb_matrix,apb_matrix,eigenvalue,bigx,bigy)
#else
 call diago_4blocks_chol(nmat,desc_apb,m_apb,n_apb,amb_matrix,apb_matrix,eigenvalue,&
                              desc_x,m_x,n_x,bigx,bigy)
#endif

 ! Deallocate the non-necessary matrices
 write(stdout,*) 'Deallocate (A+B) and (A-B) matrices'
 call clean_deallocate('A+B',apb_matrix)
 call clean_deallocate('A-B',amb_matrix)


 !
 ! Second part of the RPA correlation energy: sum over positive eigenvalues
 rpa_correlation = rpa_correlation + 0.50_dp * SUM( ABS(eigenvalue(:)) )
 if(is_rpa) then
  write(stdout,'(/,a)') ' Calculate the RPA energy using the Tamm-Dancoff decomposition'
  write(stdout,'(a)')   ' Eq. (9) from J. Chem. Phys. 132, 234114 (2010)'
  write(stdout,'(/,a,f16.10)') ' RPA energy (Ha): ',rpa_correlation
 endif

 write(stdout,'(/,a,f12.6)') ' Lowest neutral excitation energy (eV):',MINVAL(ABS(eigenvalue(:)))*Ha_eV

 !
 ! Calculate the optical sprectrum
 ! and the dynamic dipole tensor
 !
 if( calc_type%is_td .OR. calc_type%is_bse ) then
   call optical_spectrum(basis,occupation,c_matrix,wpol_out,m_x,n_x,bigx,bigy,eigenvalue)
   call stopping_power(basis,occupation,c_matrix,wpol_out,m_x,n_x,bigx,bigy,eigenvalue)
 endif

 !
 ! Calculate Wp= v * chi * v    if necessary
 ! and then write it down on file
 !
 if( print_w_ .OR. calc_type%is_gw ) then
   if( has_auxil_basis) then
     call chi_to_sqrtvchisqrtv_auxil(basis%nbf,auxil_basis%nbf_local,desc_x,m_x,n_x,bigx,bigy,eigenvalue,wpol_out,energy_gm)
     ! This following coding of the Galitskii-Migdal correlation energy is only working with
     ! an auxiliary basis
     if(is_rpa) write(stdout,'(a,f16.10,/)') ' Correlation energy in the Galitskii-Migdal formula (Ha): ',energy_gm
     
     ! Add the single pole approximation for the poles that have been neglected
     ! in the diagonalization
     if( nvirtual_SPA < nvirtual_W .AND. is_rpa ) & 
        call chi_to_sqrtvchisqrtv_auxil_spa(basis%nbf,auxil_basis%nbf_local,a_diag,wpol_out)

   else
     call chi_to_vchiv(basis%nbf,nstate0,prod_basis,c_matrix,bigx,bigy,eigenvalue,wpol_out)
   endif
  
 
   ! If requested write the spectral function on file
   if( print_w_ ) call write_spectral_function(wpol_out)

 endif

 if( .NOT. calc_type%is_gw ) call destroy_spectral_function(wpol_out)

 write(stdout,*) 'Deallocate eigenvector arrays'
 call clean_deallocate('X',bigx)
 call clean_deallocate('Y',bigy)


 if(ALLOCATED(eigenvalue)) deallocate(eigenvalue)
 if(ALLOCATED(a_diag))     deallocate(a_diag)

 call stop_clock(timing_pola)


end subroutine polarizability


!=========================================================================
subroutine build_amb_apb_common(nmat,nbf,nstate,c_matrix,energy,wpol,alpha_local,m_apb,n_apb,amb_matrix,apb_matrix,rpa_correlation)
 use m_tools 
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                 :: nmat,nbf,nstate
 real(dp),intent(in)                :: energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(in)                :: alpha_local
 integer,intent(in)                 :: m_apb,n_apb
 real(prec_td),intent(out)          :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
 real(dp),intent(out)               :: rpa_correlation
!=====
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
 integer              :: klmin
 real(dp),allocatable :: eri_eigenstate_klmin(:,:,:,:)
 real(dp)             :: eri_eigen_ijkl
 real(dp)             :: eri_eigen_ikjl,eri_eigen_iljk
 logical              :: k_is_klmin
 integer              :: iprow,ipcol
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: amb_block(:,:)
 real(dp),allocatable :: apb_block(:,:)
!=====

 call start_clock(timing_build_common)

 write(stdout,'(a)') ' Build Common part: Energies + Hartree + possibly Exchange'
 write(stdout,'(a,f8.3)') ' Content of Exchange: ',alpha_local

 if( .NOT. has_auxil_basis) then
   allocate(eri_eigenstate_klmin(nstate,nstate,nstate,nspin))
   ! Set this to zero and then enforce the calculation of the first series of
   ! Coulomb integrals
   eri_eigenstate_klmin(:,:,:,:) = 0.0_dp
 endif

 if( nprow_sd * npcol_sd > 1 ) &
    write(stdout,'(a,i4,a,i4)') ' SCALAPACK grid    :',nprow_sd,' x ',npcol_sd

 rpa_correlation = 0.0_dp
 !
 ! Set up energy+hartree+exchange contributions to matrices (A+B) and (A-B) 
 !
 
 ! First loops over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_apb_block = row_block_size(nmat,iprow,nprow_sd)
     n_apb_block = col_block_size(nmat,ipcol,npcol_sd)

     allocate(amb_block(m_apb_block,n_apb_block))
     allocate(apb_block(m_apb_block,n_apb_block))

     ! Then loop inside each of the SCALAPACK blocks
     do t_kl=1,n_apb_block
       t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)
       kstate = wpol%transition_table_apb(1,t_kl_global)
       lstate = wpol%transition_table_apb(2,t_kl_global)
       klspin = wpol%transition_table_apb(3,t_kl_global)
    
       if( .NOT. has_auxil_basis ) then
         klmin = MIN(kstate,lstate)
         k_is_klmin = (klmin == kstate)
         call calculate_eri_4center_eigen(nbf,nstate,c_matrix,klmin,klspin,eri_eigenstate_klmin)
       endif
    
       do t_ij=1,m_apb_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)
    
         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ij_global < t_kl_global ) cycle
    
         if(has_auxil_basis) then
           eri_eigen_ijkl = eri_eigen_ri_paral(istate,jstate,ijspin,kstate,lstate,klspin)
         else
           if(k_is_klmin) then ! treating (k,l)
             eri_eigen_ijkl = eri_eigenstate_klmin(lstate,istate,jstate,ijspin)
           else                ! treating (l,k)
             eri_eigen_ijkl = eri_eigenstate_klmin(kstate,istate,jstate,ijspin)
           endif
         endif
    
         if( .NOT. is_triplet) then
           apb_block(t_ij,t_kl) = 2.0_dp * eri_eigen_ijkl * spin_fact
         else
           apb_block(t_ij,t_kl) = 0.0_dp
         endif
         amb_block(t_ij,t_kl) = 0.0_dp
    
         if( alpha_local > 1.0e-6_dp ) then
           if(ijspin==klspin) then
             if(has_auxil_basis) then
               eri_eigen_ikjl = eri_eigen_ri_paral(istate,kstate,ijspin,jstate,lstate,klspin)
               eri_eigen_iljk = eri_eigen_ri_paral(istate,lstate,ijspin,jstate,kstate,klspin)
             else
               if(k_is_klmin) then
                 eri_eigen_ikjl = eri_eigenstate_klmin(istate,jstate,lstate,klspin)
                 eri_eigen_iljk = eri_eigenstate_klmin(jstate,istate,lstate,klspin)
               else
                 eri_eigen_ikjl = eri_eigenstate_klmin(jstate,istate,kstate,klspin)
                 eri_eigen_iljk = eri_eigenstate_klmin(istate,jstate,kstate,klspin)
               endif
             endif
             apb_block(t_ij,t_kl) = apb_block(t_ij,t_kl) - eri_eigen_ikjl * alpha_local - eri_eigen_iljk * alpha_local
             amb_block(t_ij,t_kl) = amb_block(t_ij,t_kl) - eri_eigen_ikjl * alpha_local + eri_eigen_iljk * alpha_local
           endif
         endif
    
         if( t_ij_global == t_kl_global ) then
           !
           ! Only one proc should add the diagonal
           if( rank == 0 ) then
             apb_block(t_ij,t_kl) =  apb_block(t_ij,t_kl) + ( energy(lstate,klspin) - energy(kstate,klspin) )
             amb_block(t_ij,t_kl) =  amb_block(t_ij,t_kl) + ( energy(lstate,klspin) - energy(kstate,klspin) )
           endif
           ! First part of the RPA correlation energy: sum over diagonal terms
           rpa_correlation = rpa_correlation - 0.25_dp * ( apb_block(t_ij,t_kl) + amb_block(t_ij,t_kl) )
         endif
    
       enddo 
    
     enddo 

     call xsum(amb_block)
     call xsum(apb_block)

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       amb_matrix(:,:) = amb_block(:,:)
       apb_matrix(:,:) = apb_block(:,:)
     endif
     deallocate(amb_block)
     deallocate(apb_block)
   enddo 
 enddo 

#ifdef HAVE_SCALAPACK
 call xsum(rpa_correlation)
#endif

 if(ALLOCATED(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)


 call stop_clock(timing_build_common)

end subroutine build_amb_apb_common


!=========================================================================
subroutine build_a_diag_common(nmat,nbf,nstate,c_matrix,energy,wpol,a_diag)
 use m_spectral_function
 use m_eri_ao_mo
 use m_tools 
 implicit none

 integer,intent(in)                 :: nmat,nbf,nstate
 real(dp),intent(in)                :: energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 real(prec_td),intent(out)          :: a_diag(wpol%npole_reso_spa)
!=====
 integer              :: t_kl,t_kl_global
 integer              :: kstate,lstate
 integer              :: klspin
 integer              :: klmin
 real(dp),allocatable :: eri_eigenstate_klmin(:,:,:,:)
 real(dp)             :: eri_eigen_klkl
 logical              :: k_is_klmin
 real(dp),parameter   :: empirical_fact=1.5_dp
 character(len=100)   :: ctmp
!=====

 call start_clock(timing_build_common)

 write(stdout,'(a)') ' Build diagonal of the RPA part: Energies + Hartree'
 if( ABS( empirical_fact - 1.0_dp ) > 1.0e-6_dp ) then
   write(ctmp,'(a,1x,f6.2)') 'Empirical parameter',empirical_fact
   call issue_warning(ctmp)
 endif

 if( .NOT. has_auxil_basis) then
   allocate(eri_eigenstate_klmin(nstate,nstate,nstate,nspin))
   ! Set this to zero and then enforce the calculation of the first series of
   ! Coulomb integrals
   eri_eigenstate_klmin(:,:,:,:) = 0.0_dp
 endif

 !
 ! Set up energy+hartree+exchange contributions to matrices (A+B) and (A-B) 
 !


 ! Then loop inside each of the SCALAPACK blocks
 do t_kl=1,wpol%npole_reso_spa
   kstate = wpol%transition_table_spa(1,t_kl)
   lstate = wpol%transition_table_spa(2,t_kl)
   klspin = wpol%transition_table_spa(3,t_kl)

   if( .NOT. has_auxil_basis ) then
     klmin = MIN(kstate,lstate)
     k_is_klmin = (klmin == kstate)
     call calculate_eri_4center_eigen(nbf,nstate,c_matrix,klmin,klspin,eri_eigenstate_klmin)
   endif

   if(has_auxil_basis) then
     eri_eigen_klkl = eri_eigen_ri(kstate,lstate,klspin,kstate,lstate,klspin)
   else
     if(k_is_klmin) then ! treating (k,l)
       eri_eigen_klkl = eri_eigenstate_klmin(lstate,kstate,lstate,klspin)
     else                ! treating (l,k)
       eri_eigen_klkl = eri_eigenstate_klmin(kstate,kstate,lstate,klspin)
     endif
   endif

   a_diag(t_kl) = eri_eigen_klkl * spin_fact + energy(lstate,klspin) - energy(kstate,klspin)
   a_diag(t_kl) = a_diag(t_kl) * empirical_fact

 enddo 

 if(ALLOCATED(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)


 call stop_clock(timing_build_common)

end subroutine build_a_diag_common


!=========================================================================
subroutine build_apb_tddft(nmat,basis,c_matrix,occupation,wpol,m_apb,n_apb,apb_matrix)
 use m_spectral_function
 use m_basis_set
 use m_dft_grid
 implicit none

 integer,intent(in)                 :: nmat
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)                :: occupation(basis%nbf,nspin)
 type(spectral_function),intent(in) :: wpol
 integer,intent(in)                 :: m_apb,n_apb
 real(prec_td),intent(inout)        :: apb_matrix(m_apb,n_apb)
!=====
 integer              :: nspin_tddft
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
 logical              :: require_gradient
 real(dp),allocatable :: grad_ij(:,:,:),grad_kl(:,:,:)
 real(dp),allocatable :: dot_ij_kl(:,:),dot_rho_ij(:,:),dot_rho_kl(:,:)
 real(dp),allocatable :: v2rho2(:,:)
 real(dp),allocatable :: vsigma(:,:)
 real(dp),allocatable :: v2rhosigma(:,:)
 real(dp),allocatable :: v2sigma2(:,:)
 real(dp),allocatable :: wf_r(:,:,:)
 real(dp),allocatable :: wf_gradr(:,:,:,:)
 real(dp),allocatable :: rho_gradr(:,:,:)
 real(dp)             :: xctmp
 integer              :: ipcol,iprow
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: apb_block(:,:)
!=====

 call start_clock(timing_build_tddft)

 write(stdout,'(a)') ' Build fxc part'

 if( is_triplet ) then
   nspin_tddft = 2
 else
   nspin_tddft = nspin
 endif
 !
 ! Prepare TDDFT calculations
 call prepare_tddft(nspin_tddft,basis,c_matrix,occupation,v2rho2,vsigma,v2rhosigma,v2sigma2,wf_r,wf_gradr,rho_gradr)
 require_gradient = .FALSE.
 if(ALLOCATED(v2sigma2)) then ! GGA case
   require_gradient = .TRUE.
   allocate(grad_ij(3,ngrid,nspin))
   allocate(grad_kl(3,ngrid,nspin))
   allocate(dot_ij_kl(ngrid,nspin))
   allocate(dot_rho_ij(ngrid,nspin))
   allocate(dot_rho_kl(ngrid,nspin))
 endif


 if( nprow_sd * npcol_sd > 1 ) &
    write(stdout,'(a,i4,a,i4)') ' SCALAPACK grid    :',nprow_sd,' x ',npcol_sd


 ! First loops over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_apb_block = row_block_size(nmat,iprow,nprow_sd)
     n_apb_block = col_block_size(nmat,ipcol,npcol_sd)

     allocate(apb_block(m_apb_block,n_apb_block))
     apb_block(:,:) = 0.0_dp


     !
     ! Set up fxc contributions to matrices (A+B)
     !
     do t_kl=1,n_apb_block
       t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)
       kstate = wpol%transition_table_apb(1,t_kl_global)
       lstate = wpol%transition_table_apb(2,t_kl_global)
       klspin = wpol%transition_table_apb(3,t_kl_global)
    
       do t_ij=1,m_apb_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)
    
         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ij_global < t_kl_global ) cycle
    
    
         if( nspin_tddft == 1 ) then
    
           xctmp = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                              * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                              * v2rho2(:,ijspin) * 2.0_dp )
    
    
           if(require_gradient) then
    
             grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
             dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                                  + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
             dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                                   + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
             dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                                   + rho_gradr(3,:,1) * grad_kl(3,:,klspin)
    
             xctmp = xctmp   &
                   +  SUM( dot_ij_kl(:,1) * 4.0_dp * vsigma(:,1) ) &
                   +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * 8.0_dp * v2sigma2(:,1) ) &
                   +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                           + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                             * 4.0_dp * v2rhosigma(:,1) )
    
           endif
    
         else if( .NOT. is_triplet ) then
    
           xctmp = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                              * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                              * ( v2rho2(:,1) + v2rho2(:,2) ) )
    
    
           if(require_gradient) then
    
             grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
             dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                                  + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
             dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                                   + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
             dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                                   + rho_gradr(3,:,1) * grad_kl(3,:,klspin)
    
             xctmp = xctmp   &
                   +  SUM( dot_ij_kl(:,1) * ( 2.0_dp * vsigma(:,1) + vsigma(:,2) ) ) &
                   +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * ( v2sigma2(:,1) + 0.5_dp * v2sigma2(:,2) +  2.0_dp * v2sigma2(:,3) + v2sigma2(:,4) ) ) &
                   +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                          + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                             * ( v2rhosigma(:,1) + v2rhosigma(:,2) + v2rhosigma(:,3) )  )
           endif
    
    
           
         else ! triplet case
    
           xctmp = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                              * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                              * ( v2rho2(:,1) - v2rho2(:,2) ) )
    
           if(require_gradient) then
    
             grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
             dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                                  + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
             dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                                   + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
             dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                                   + rho_gradr(3,:,1) * grad_kl(3,:,klspin)
    
             xctmp = xctmp   &
                   +  SUM( dot_ij_kl(:,1) * ( 2.0_dp * vsigma(:,1) - vsigma(:,2) ) ) &
                   +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * ( v2sigma2(:,1) - v2sigma2(:,3) ) ) &   !FIXME 3 or 5 are working, but only one is correct in principle
                   +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                          + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                             * ( v2rhosigma(:,1) - v2rhosigma(:,4) )  )   !FIXME 3 and 4 are working, but only one is correct in principle
           endif
    
    
    
         endif
    
         ! The factor two accounts for (A+B), and not A or B.
         apb_block(t_ij,t_kl) = apb_block(t_ij,t_kl) + 2.0_dp * xctmp
    
       enddo
     enddo

     call xsum(apb_block)

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       apb_matrix(:,:) = apb_matrix(:,:) + apb_block(:,:)
     endif
     deallocate(apb_block)

   enddo
 enddo


 deallocate(wf_r)

 if(require_gradient) then
   deallocate(grad_ij,grad_kl)
   deallocate(dot_ij_kl,dot_rho_ij,dot_rho_kl)
   deallocate(v2rho2)
   deallocate(vsigma)
   deallocate(v2rhosigma)
   deallocate(v2sigma2)
   deallocate(wf_gradr)
   deallocate(rho_gradr)
 endif


 call stop_clock(timing_build_tddft)

end subroutine build_apb_tddft


!=========================================================================
subroutine build_amb_apb_bse(nbf,prod_basis,c_matrix,wpol,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
 use m_spectral_function
 use m_basis_set
 use m_eri_ao_mo
 use m_tools 
 implicit none

 integer,intent(in)                 :: nbf
 type(basis_set),intent(in)         :: prod_basis
 real(dp),intent(in)                :: c_matrix(nbf,nbf,nspin)
 type(spectral_function),intent(in) :: wpol,wpol_static
 integer,intent(in)                 :: m_apb,n_apb
 real(prec_td),intent(out)          :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
!=====
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
 integer              :: kbf
 real(dp),allocatable :: bra(:),ket(:)
 real(dp)             :: wtmp
!=====

 if( has_auxil_basis ) call die('Have an auxil basis. This should not happen here')

 call start_clock(timing_build_bse)

 write(stdout,'(a)') ' Build W part'

 !
 ! Prepare the bra and ket for BSE
 allocate(bra(wpol_static%npole_reso))
 allocate(ket(wpol_static%npole_reso))

 !
 ! Set up -W contributions to matrices (A+B) and (A-B)
 !
 do t_kl=1,n_apb
   t_kl_global = colindex_local_to_global('C',t_kl)
   kstate = wpol%transition_table_apb(1,t_kl_global)
   lstate = wpol%transition_table_apb(2,t_kl_global)
   klspin = wpol%transition_table_apb(3,t_kl_global)

   do t_ij=1,m_apb
     t_ij_global = rowindex_local_to_global('C',t_ij)
     istate = wpol%transition_table_apb(1,t_ij_global)
     jstate = wpol%transition_table_apb(2,t_ij_global)
     ijspin = wpol%transition_table_apb(3,t_ij_global)

     !
     ! Only calculate the lower triangle
     ! Symmetrization will be performed later (in the diago subroutines)
     if( t_ij_global < t_kl_global ) cycle

     if(ijspin/=klspin) cycle

     kbf = prod_basis%index_prodbasis(istate,kstate)+prod_basis%nbf*(ijspin-1)
     bra(:) = wpol_static%residu_left(kbf,:)
     kbf = prod_basis%index_prodbasis(jstate,lstate)+prod_basis%nbf*(klspin-1)
     ket(:) = wpol_static%residu_left(kbf,:)

     wtmp =  SUM( 2.0_dp * bra(:)*ket(:)/(-wpol_static%pole(:)) )   ! Factor two comes from Resonant and Anti-resonant transitions
     apb_matrix(t_ij,t_kl) =  apb_matrix(t_ij,t_kl) - wtmp
     amb_matrix(t_ij,t_kl) =  amb_matrix(t_ij,t_kl) - wtmp

     kbf = prod_basis%index_prodbasis(istate,lstate)+prod_basis%nbf*(ijspin-1)
     bra(:) = wpol_static%residu_left(kbf,:)
     kbf = prod_basis%index_prodbasis(jstate,kstate)+prod_basis%nbf*(klspin-1)
     ket(:) = wpol_static%residu_left(kbf,:)  

     wtmp =  SUM( 2.0_dp * bra(:)*ket(:)/(-wpol_static%pole(:)) )   ! Factor two comes from Resonant and Anti-resonant transitions
     apb_matrix(t_ij,t_kl) =  apb_matrix(t_ij,t_kl) - wtmp
     amb_matrix(t_ij,t_kl) =  amb_matrix(t_ij,t_kl) + wtmp


   enddo

 enddo

 deallocate(bra,ket)


 call stop_clock(timing_build_bse)


end subroutine build_amb_apb_bse


!=========================================================================
subroutine build_amb_apb_bse_auxil(nmat,nbf,c_matrix,wpol,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
 use m_spectral_function
 use m_basis_set
 use m_eri_ao_mo
 use m_tools 
 implicit none

 integer,intent(in)                 :: nmat,nbf
 real(dp),intent(in)                :: c_matrix(nbf,nbf,nspin)
 type(spectral_function),intent(in) :: wpol,wpol_static
 integer,intent(in)                 :: m_apb,n_apb
 real(prec_td),intent(out)          :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
!=====
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: kstate_prev
 integer              :: ijspin,klspin
 integer              :: kbf
 real(dp)             :: wtmp
 integer              :: kstate_max
 integer              :: ipole,ibf_auxil,jbf_auxil,ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable :: vsqrt_chi_vsqrt(:,:)
 real(dp),allocatable :: vsqrt_chi_vsqrt_i(:),residu_i(:),wp0_i(:,:)
 real(dp),allocatable :: wp0(:,:,:,:),w0_local(:)
 integer              :: iprow,ipcol
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: amb_block(:,:)
 real(dp),allocatable :: apb_block(:,:)
!=====

 call start_clock(timing_build_bse)
 if( .NOT. has_auxil_basis ) call die('Does not have auxil basis. This should not happen')

 write(stdout,'(a)') ' Build W part Auxil' 

 kstate_max = MAXVAL( wpol%transition_table_apb(1,1:wpol%npole_reso_apb) )

 call clean_allocate('Temporary array for W',wp0,1,nauxil_3center,ncore_W+1,nvirtual_W-1,ncore_W+1,kstate_max,1,nspin)


 do ijspin=1,nspin
#ifndef HAVE_SCALAPACK

   allocate(vsqrt_chi_vsqrt(nauxil_2center,nauxil_2center))

   !
   ! Test if w0 is already available or if we need to calculate it first
   if( ALLOCATED(wpol_static%w0) ) then

     vsqrt_chi_vsqrt(:,:) = wpol_static%w0(:,:)

   else ! w0 is not available

    
     vsqrt_chi_vsqrt(:,:) = 0.0_dp
     do ipole=1,wpol_static%npole_reso
       do jbf_auxil=1,nauxil_2center
         vsqrt_chi_vsqrt(:,jbf_auxil) = vsqrt_chi_vsqrt(:,jbf_auxil) &
                  - wpol_static%residu_left(:,ipole)                 &
                    * wpol_static%residu_left(jbf_auxil,ipole) * 2.0_dp / wpol_static%pole(ipole)
       enddo
     enddo
   endif
    
   !
   ! The last index of wp0 only runs on occupied states (to save memory and CPU time)
   ! Be careful in the following not to forget it
   do kstate=ncore_W+1,kstate_max
     wp0(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) = MATMUL( vsqrt_chi_vsqrt(:,:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) )
   enddo
  
   deallocate(vsqrt_chi_vsqrt)


#else

   !
   ! Test if w0 is already available or if we need to calculate it first
   if( ALLOCATED(wpol_static%w0) ) then

     allocate(wp0_i(ncore_W+1:nvirtual_W-1,ncore_W+1:kstate_max))
     allocate(w0_local(nauxil_3center))

     do ibf_auxil_global=1,nauxil_2center

       do jbf_auxil=1,nauxil_3center
         jbf_auxil_global = ibf_auxil_g(jbf_auxil)
         w0_local(jbf_auxil) = wpol_static%w0(ibf_auxil_global,jbf_auxil_global)
       enddo

       do kstate=ncore_W+1,kstate_max
         wp0_i(ncore_W+1:nvirtual_W-1,kstate) = MATMUL( w0_local(:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) )
       enddo
       call xsum(wp0_i)

       if( iproc_ibf_auxil(ibf_auxil_global) == rank ) then
         wp0(ibf_auxil_l(ibf_auxil_global),:,:,ijspin) = wp0_i(:,:)
       endif

     enddo
     deallocate(w0_local)

   else ! w0 is not available

     allocate(vsqrt_chi_vsqrt_i(nauxil_3center))
     allocate(residu_i(wpol_static%npole_reso))
     allocate(wp0_i(ncore_W+1:nvirtual_W-1,ncore_W+1:kstate_max))
    
     do ibf_auxil=1,nauxil_2center
    
       if( iproc_ibf_auxil(ibf_auxil) == rank ) then
         residu_i(:) = wpol_static%residu_left(ibf_auxil_l(ibf_auxil),:)
       else
         residu_i(:) = 0.0_dp
       endif
       call xsum(residu_i)
    
       vsqrt_chi_vsqrt_i(:) = 0.0_dp
       do ipole=1,wpol_static%npole_reso
         vsqrt_chi_vsqrt_i(:) = vsqrt_chi_vsqrt_i(:) &
                  - residu_i(ipole) * wpol_static%residu_left(:,ipole) * 2.0_dp / wpol_static%pole(ipole)
       enddo
       !
       ! The last index of wp0 only runs on occupied states (to save memory and CPU time)
       ! Be careful in the following not to forget it
       do kstate=ncore_W+1,kstate_max
         wp0_i(ncore_W+1:nvirtual_W-1,kstate) = MATMUL( vsqrt_chi_vsqrt_i(:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) )
       enddo
       call xsum(wp0_i)
    
       if( iproc_ibf_auxil(ibf_auxil) == rank ) then
         wp0(ibf_auxil_l(ibf_auxil),:,:,ijspin) = wp0_i(:,:)
       endif
    
     enddo
    
     deallocate(vsqrt_chi_vsqrt_i,residu_i,wp0_i)

   endif

#endif
 enddo
 

 if( nprow_sd * npcol_sd > 1 ) &
    write(stdout,'(a,i4,a,i4)') ' SCALAPACK grid    :',nprow_sd,' x ',npcol_sd
 
 ! First loops over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_apb_block = row_block_size(nmat,iprow,nprow_sd)
     n_apb_block = col_block_size(nmat,ipcol,npcol_sd)

     allocate(amb_block(m_apb_block,n_apb_block))
     allocate(apb_block(m_apb_block,n_apb_block))


     ! Set up -W contributions to matrices (A+B) and (A-B)
     !
     do t_kl=1,n_apb_block
       t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)
       kstate = wpol%transition_table_apb(1,t_kl_global)
       lstate = wpol%transition_table_apb(2,t_kl_global)
       klspin = wpol%transition_table_apb(3,t_kl_global)

       do t_ij=1,m_apb_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)

         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ij_global < t_kl_global ) cycle

         if(ijspin/=klspin) cycle

         wtmp = DOT_PRODUCT( eri_3center_eigen(:,jstate,lstate,ijspin) , wp0(:,istate,kstate,ijspin) )

         apb_block(t_ij,t_kl) = -wtmp
         amb_block(t_ij,t_kl) = -wtmp


         wtmp = DOT_PRODUCT( eri_3center_eigen(:,istate,lstate,ijspin) , wp0(:,jstate,kstate,ijspin) )

         apb_block(t_ij,t_kl) =  apb_block(t_ij,t_kl) - wtmp
         amb_block(t_ij,t_kl) =  amb_block(t_ij,t_kl) + wtmp


       enddo

     enddo

     call xsum(amb_block)
     call xsum(apb_block)
     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       amb_matrix(:,:) = amb_matrix(:,:) + amb_block(:,:)
       apb_matrix(:,:) = apb_matrix(:,:) + apb_block(:,:)
     endif
     deallocate(amb_block)
     deallocate(apb_block)


   enddo
 enddo

 call clean_deallocate('Temporary array for W',wp0)


 call stop_clock(timing_build_bse)


end subroutine build_amb_apb_bse_auxil


!=========================================================================
subroutine diago_4blocks_sqrt(nmat,amb_matrix,apb_matrix,eigenvalue,bigx,bigy)
 use m_spectral_function
! use m_eri
 use m_tools 
 implicit none

 integer,intent(in)           :: nmat
 real(prec_td),intent(inout)  :: amb_matrix(nmat,nmat)
 real(prec_td),intent(inout)  :: apb_matrix(nmat,nmat)  ! apb_matrix constains (A+B) in the input, however it used a matrix buffer after
 real(dp),intent(out)         :: eigenvalue(nmat)
 real(prec_td),intent(out)    :: bigx(nmat,nmat),bigy(nmat,nmat)
!=====
 integer                   :: t_ij,t_kl
 real(prec_td),allocatable :: amb_eigval(:),bigomega(:)
!=====

 call start_clock(timing_diago_h2p)

 ! First symmetrize the matrices since only the lower triangle was calculated
 do t_kl=1,nmat
   do t_ij=t_kl+1,nmat
     amb_matrix(t_kl,t_ij) = amb_matrix(t_ij,t_kl)
     apb_matrix(t_kl,t_ij) = apb_matrix(t_ij,t_kl)
   enddo
 enddo



 write(stdout,'(/,a)') ' Performing the block diago with square root of matrices'

 !
 ! Calculate (A-B)^{1/2}
 ! First diagonalize (A-B):
 ! (A-B) R = R D
 ! (A-B) is real symmetric, hence R is orthogonal R^{-1} = tR
 ! (A-B)       = R D tR 
 ! (A-B)^{1/2} = R D^{1/2} tR 
 write(stdout,'(a,i8,a,i8)') ' Diago to get (A - B)^{1/2}                   ',nmat,' x ',nmat
 allocate(amb_eigval(nmat))
 call diagonalize(nmat,amb_matrix,amb_eigval)


 ! bigx contains the (A-B)**1/2
 ! bigy contains the (A-B)**-1/2
 forall(t_kl=1:nmat)
   bigx(:,t_kl) = amb_matrix(:,t_kl)*SQRT(amb_eigval(t_kl))
   bigy(:,t_kl) = amb_matrix(:,t_kl)/SQRT(amb_eigval(t_kl))
 end forall
 deallocate(amb_eigval)

 amb_matrix = TRANSPOSE( amb_matrix )
 bigx(:,:) = MATMUL( bigx(:,:) , amb_matrix(:,:) )
 bigy(:,:) = MATMUL( bigy(:,:) , amb_matrix(:,:) )
 
 ! Use amb_matrix as a temporary matrix here:
 amb_matrix(:,:) = MATMUL( apb_matrix , bigx )
 apb_matrix(:,:)  = MATMUL( bigx, amb_matrix )

! write(stdout,*) 'CC ',matrix_is_symmetric(nmat,apb_matrix)


 write(stdout,'(a,i8,a,i8)') ' Diago (A - B)^{1/2} * (A + B) * (A - B)^{1/2}',nmat,' x ',nmat
 allocate(bigomega(nmat))
 call diagonalize(nmat,apb_matrix,bigomega)

 bigomega(:) = SQRT(bigomega(:))

 forall(t_kl=1:nmat)
   apb_matrix(:,t_kl) = apb_matrix(:,t_kl) / SQRT(bigomega(t_kl))
   eigenvalue(t_kl) = bigomega(t_kl)
 end forall

 ! Save (A-B)**-1/2 in amb_matrix 
 amb_matrix(:,:) = bigy(:,:)

 bigx(:,:) = 0.5_dp * MATMUL( bigx(:,:)   , apb_matrix(:,:) )
 bigy(:,:) = bigx(:,:)

 apb_matrix(:,:) = 0.5_dp * MATMUL( amb_matrix(:,:) , apb_matrix(:,:) )
 forall(t_kl=1:nmat)
   apb_matrix(:,t_kl) = apb_matrix(:,t_kl) * bigomega(t_kl)
 end forall
 deallocate(bigomega)

 ! Finalize Resonant (positive excitations second index from 1 to nmat)
 bigx(:,:) = bigx(:,:) + apb_matrix(:,:)
 bigy(:,:) = bigy(:,:) - apb_matrix(:,:)

 call stop_clock(timing_diago_h2p)

end subroutine diago_4blocks_sqrt


!=========================================================================
subroutine diago_4blocks_chol(nmat,desc_apb,m_apb,n_apb,amb_matrix,apb_matrix,&
                              eigenvalue,desc_x,m_x,n_x,bigx,bigy)
 use m_spectral_function
! use m_eri
 use m_tools 
 implicit none

 integer,intent(in)     :: nmat,m_apb,n_apb,m_x,n_x
 integer,intent(in)     :: desc_apb(ndel),desc_x(ndel)
 real(dp),intent(inout) :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
 real(dp),intent(out)   :: eigenvalue(nmat)
 real(dp),intent(out)   :: bigx(m_x,n_x)
 real(dp),intent(out)   :: bigy(m_x,n_x)
!=====
 integer  :: info
 integer  :: lwork,liwork
 real(dp),allocatable :: work(:)
 integer,allocatable :: iwork(:)
!=====

#ifdef HAVE_SCALAPACK
 call start_clock(timing_diago_h2p)

 ! First allocate eigenvector
 write(stdout,'(/,a)') ' Performing the block diago with Cholesky'

 allocate(work(1))
 allocate(iwork(1))
 lwork=-1
 liwork=-1
 call pdbssolver1(nmat,apb_matrix,1,1,desc_apb,amb_matrix,1,1,desc_apb,    &
                  eigenvalue,bigx,1,1,desc_x,bigy,                         &
                  work,lwork,iwork,liwork,info)
 if(info/=0) call die('SCALAPACK failed')

 lwork  = NINT(work(1))
 deallocate(work)
 call clean_allocate('Buffer array for SCALAPACK diago',work,lwork)

 liwork = iwork(1)
 deallocate(iwork)
 allocate(iwork(liwork))

 call pdbssolver1(nmat,apb_matrix,1,1,desc_apb,amb_matrix,1,1,desc_apb,    &
                  eigenvalue,bigx,1,1,desc_x,bigy,                         &
                  work,lwork,iwork,liwork,info)
 if(info/=0) call die('SCALAPACK failed')

 call clean_deallocate('Buffer array for SCALAPACK diago',work)
 deallocate(iwork)



 call stop_clock(timing_diago_h2p)

#else
 call die('Cholesky diago cannot run without SCALAPACK')
#endif

end subroutine diago_4blocks_chol


!=========================================================================
subroutine optical_spectrum(basis,occupation,c_matrix,chi,m_x,n_x,bigx,bigy,eigenvalue)
 use m_tools
 use m_basis_set
 use m_eri
 use m_dft_grid
 use m_spectral_function
 use m_atoms
 implicit none

 integer,intent(in)                 :: m_x,n_x
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 type(spectral_function),intent(in) :: chi
 real(prec_td),intent(in)           :: bigx(m_x,n_x)
 real(prec_td),intent(in)           :: bigy(m_x,n_x)
 real(dp),intent(in)                :: eigenvalue(chi%npole_reso_apb)
!=====
 integer                            :: t_ij,t_kl
 integer                            :: t_ij_global,t_kl_global
 integer                            :: nmat
 integer                            :: istate,jstate,ijspin
 integer                            :: ibf,jbf
 integer                            :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart,ibf_cart,jbf_cart
 integer                            :: iomega,idir,jdir
 integer,parameter                  :: nomega=600
 complex(dp)                        :: omega(nomega)
 real(dp)                           :: coeff,trace
 real(dp)                           :: dynamical_pol(nomega,3,3),photoabsorp_cross(nomega,3,3)
 real(dp)                           :: static_polarizability(3,3)
 real(dp)                           :: oscillator_strength,trk_sumrule,mean_excitation
 real(dp),allocatable               :: dipole_basis(:,:,:),dipole_tmp(:,:,:),dipole_state(:,:,:,:)
 real(dp),allocatable               :: dipole_cart(:,:,:)
 real(dp),allocatable               :: residu_left(:,:)
 integer                            :: dynpolfile
 integer                            :: photocrossfile
 integer                            :: parityi,parityj,reflectioni,reflectionj
 integer,external                   :: wfn_parity
 integer,external                   :: wfn_reflection
 character(len=32)                  :: symsymbol
!=====


 call start_clock(timing_spectrum)
 !
 ! Calculate the spectrum now
 !

 write(stdout,'(/,a)') ' Calculate the optical spectrum'

 if (nspin/=1) then
   msg='no nspin/=1 allowed'
   call issue_warning(msg)
   return
 endif

 !
 ! First precalculate all the needed dipole in the basis set
 !
 allocate(dipole_basis(3,basis%nbf,basis%nbf))
 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   jbf_cart = 1
   jbf      = 1
   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(dipole_cart(3,ni_cart,nj_cart))


     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_dipole(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),dipole_cart(:,i_cart,j_cart))
       enddo
     enddo

     do idir=1,3
       dipole_basis(idir,ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE( cart_to_pure(li)%matrix(:,:) ) , &
             MATMUL(  dipole_cart(idir,:,:) , cart_to_pure(lj)%matrix(:,:) ) )
     enddo

     deallocate(dipole_cart)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo

 !
 ! Get the dipole oscillator strength on states
 allocate(dipole_state(3,basis%nbf,basis%nbf,nspin))
 allocate(dipole_tmp(3,basis%nbf,basis%nbf))

 do ijspin=1,nspin
   do jbf=1,basis%nbf
     do istate=1,basis%nbf
       dipole_tmp(:,istate,jbf) = MATMUL( dipole_basis(:,:,jbf) , c_matrix(:,istate,ijspin) )
     enddo
   enddo

   do jstate=1,basis%nbf
     do istate=1,basis%nbf
       dipole_state(:,istate,jstate,ijspin) = MATMUL( dipole_tmp(:,istate,:) , c_matrix(:,jstate,ijspin) )
     enddo
   enddo

 enddo
 deallocate(dipole_basis,dipole_tmp)


 allocate(residu_left(3,chi%npole_reso_apb))

 nmat=chi%npole_reso_apb
 residu_left(:,:) = 0.0_dp
 do t_ij=1,m_x
   t_ij_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ij)
   istate = chi%transition_table_apb(1,t_ij_global)
   jstate = chi%transition_table_apb(2,t_ij_global)
   ijspin = chi%transition_table_apb(3,t_ij_global)

   ! Let use (i <-> j) symmetry to halve the loop
   do t_kl=1,n_x
     t_kl_global = colindex_local_to_global(ipcol_sd,npcol_sd,t_kl)

     residu_left(:,t_kl_global) = residu_left(:,t_kl_global) &
                  + dipole_state(:,istate,jstate,ijspin) * ( bigx(t_ij,t_kl) + bigy(t_ij,t_kl) ) * SQRT(spin_fact)
   enddo

 enddo
 call xsum(residu_left)

 deallocate(dipole_state)


 !
 ! Calculate the dynamical dipole polarizability
 ! and the static dipole polarizability
 !
 ! Set the frequency mesh
 omega(1)     =MAX( 0.0_dp      ,MINVAL(ABS(eigenvalue(:)))-10.00/Ha_eV)
 omega(nomega)=MIN(50.0_dp/Ha_eV,MAXVAL(ABS(eigenvalue(:)))+10.00/Ha_eV)
 do iomega=2,nomega-1
   omega(iomega) = omega(1) + ( omega(nomega)-omega(1) ) /REAL(nomega-1,dp) * (iomega-1) 
 enddo
 ! Add the broadening
 omega(:) = omega(:) + im * 0.10/Ha_eV

 dynamical_pol(:,:,:) = 0.0_dp
 static_polarizability(:,:) = 0.0_dp
 do t_ij=1,nmat
   do idir=1,3
     do jdir=1,3
       dynamical_pol(:,idir,jdir) = dynamical_pol(:,idir,jdir) &
                            + residu_left(idir,t_ij) * residu_left(jdir,t_ij) &
                              * ( AIMAG( -1.0_dp  / ( omega(:) - eigenvalue(t_ij) ) ) - AIMAG( -1.0_dp  / ( omega(:) + eigenvalue(t_ij) ) ) )
       static_polarizability(idir,jdir) = static_polarizability(idir,jdir) &
                      + 2.0_dp * residu_left(idir,t_ij) * residu_left(jdir,t_ij) / eigenvalue(t_ij)
     enddo
   enddo
 enddo
 !
 ! Get the photoabsorption cross section
 do iomega=1,nomega
   photoabsorp_cross(iomega,:,:) = 4.0_dp * pi * REAL(omega(iomega),dp) / c_speedlight * dynamical_pol(iomega,:,:)
 enddo


 write(stdout,'(/,a)') ' Excitation energies (eV)     Oscil. strengths   [Symmetry] '  
 trk_sumrule=0.0_dp
 mean_excitation=0.0_dp
 do t_kl_global=1,nmat
   t_kl = colindex_global_to_local('S',t_kl_global)

   if( is_triplet ) then 
     oscillator_strength = 0.0_dp
   else
     oscillator_strength = 2.0_dp/3.0_dp * DOT_PRODUCT(residu_left(:,t_kl_global),residu_left(:,t_kl_global)) * eigenvalue(t_kl_global)
   endif
   trk_sumrule = trk_sumrule + oscillator_strength
   mean_excitation = mean_excitation + oscillator_strength * LOG( eigenvalue(t_kl_global) )

   if(t_kl_global<=30) then

     if( is_triplet ) then
       symsymbol='3'
     else
       symsymbol='1'
     endif

     !
     ! Test the parity in case of molecule with inversion symmetry
    
     t_ij_global = 0
     do t_ij=1,m_x
       ! t_kl is zero if the proc is not in charge of this process
       if( t_kl /=0 ) then 
         if( ABS(bigx(t_ij,t_kl)) > 0.1_dp ) then
           t_ij_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ij)
           exit
         endif
       endif
     enddo
     call xmax(t_ij_global)

     istate = chi%transition_table_apb(1,t_ij_global)
     jstate = chi%transition_table_apb(2,t_ij_global)
     ijspin = chi%transition_table_apb(3,t_ij_global)
     if(planar) then
       reflectioni = wfn_reflection(basis,c_matrix,istate,ijspin)
       reflectionj = wfn_reflection(basis,c_matrix,jstate,ijspin)
       select case(reflectioni*reflectionj)
       case( 1)
         symsymbol=TRIM(symsymbol)//'(A1, B2 or Ap )'
       case(-1)
         symsymbol=TRIM(symsymbol)//'(A2, B1 or App)'
       end select
     endif
     if(inversion) then
       parityi = wfn_parity(basis,c_matrix,istate,ijspin)
       parityj = wfn_parity(basis,c_matrix,jstate,ijspin)
       select case(parityi*parityj)
       case( 1)
         symsymbol=TRIM(symsymbol)//'g'
       case(-1)
         symsymbol=TRIM(symsymbol)//'u'
       end select
     endif

     write(stdout,'(1x,i4.4,a3,2(f18.8,2x),5x,a32)') t_kl_global,' : ', &
                  eigenvalue(t_kl_global)*Ha_eV,oscillator_strength,symsymbol

     !
     ! Output the transition coefficients
     do t_ij_global=1,nmat
       t_ij = rowindex_global_to_local('S',t_ij_global)
       istate = chi%transition_table_apb(1,t_ij_global)
       jstate = chi%transition_table_apb(2,t_ij_global)
     
       coeff = 0.0_dp
       if( t_ij /= 0 .AND. t_kl /=0 ) then 
         if( ABS(bigx(t_ij,t_kl)) / SQRT(2.0_dp) > 0.1_dp ) then
           coeff = bigx(t_ij,t_kl) / SQRT(2.0_dp)
         endif
       endif
       call xsum(coeff)
       if( ABS(coeff) > 0.1_dp ) write(stdout,'(8x,i4,a,i4,1x,f12.5)') istate,' -> ',jstate,coeff

       coeff = 0.0_dp
       if( t_ij /= 0 .AND. t_kl /=0 ) then
         if( ABS(bigy(t_ij,t_kl)) / SQRT(2.0_dp) > 1.0e-1_dp ) then
           coeff = bigy(t_ij,t_kl) / SQRT(2.0_dp)
         endif
       endif
       call xsum(coeff)
       if( ABS(coeff) > 0.1_dp ) write(stdout,'(8x,i4,a,i4,1x,f12.5)') istate,' <- ',jstate,coeff
     enddo




     write(stdout,*)
   endif
 enddo

 if( is_triplet ) return

 write(stdout,'(/,a)')     ' TRK sum rule: the two following numbers should compare well'
 write(stdout,'(a,f12.6)') ' Sum over oscillator strengths: ',trk_sumrule
 write(stdout,'(a,f12.6)') '   Number of valence electrons: ',SUM( occupation(ncore_W+1:,:) )

 write(stdout,'(/,a,f12.6)') ' Mean excitation energy (eV): ',EXP( mean_excitation / trk_sumrule ) * Ha_eV

 write(stdout,'(/,a)') ' Static dipole polarizability:'
 trace = 0.0_dp
 do idir=1,3
   write(stdout,'(3(4x,f12.6))') static_polarizability(idir,:)
   trace = trace + static_polarizability(idir,idir) / 3.0_dp
 enddo
 write(stdout,'(a,f12.6)') ' Static dipole polarizability trace: ',trace

 if( is_iomaster ) then

   open(newunit=dynpolfile,file='dynamical_dipole_polarizability.dat',form='formatted')
   open(newunit=photocrossfile,file='photoabsorption_cross_section.dat',form='formatted')
   write(dynpolfile,'(a)') '#  Imaginary part of dynamical dipole polarizability'
   write(dynpolfile,'(a)') '#  omega (eV)   Average     xx    yx    zx    xy    yy    zy    xz    yz    zz'
   write(photocrossfile,'(a)') '#  Imaginary part of dynamical dipole polarizability'
   write(photocrossfile,'(a)') '#  omega (eV)   Average     xx    yx    zx    xy    yy    zy    xz    yz    zz'
   do iomega=1,nomega
     write(dynpolfile,'(11(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,                                      &
                                          (dynamical_pol(iomega,1,1)+dynamical_pol(iomega,2,2)+dynamical_pol(iomega,3,3))/3.0_dp, &
                                          dynamical_pol(iomega,:,:)
     write(photocrossfile,'(11(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,                                      &
                                              (photoabsorp_cross(iomega,1,1)+photoabsorp_cross(iomega,2,2)+photoabsorp_cross(iomega,3,3))/3.0_dp, &
                                              photoabsorp_cross(iomega,:,:)
   enddo 

   close(dynpolfile)
   close(photocrossfile)

 endif


 deallocate(residu_left)

 call stop_clock(timing_spectrum)

end subroutine optical_spectrum


!=========================================================================
subroutine stopping_power(basis,occupation,c_matrix,chi,m_x,n_x,bigx,bigy,eigenvalue)
 use m_tools
 use m_basis_set
 use m_eri
 use m_dft_grid
 use m_spectral_function
 use m_atoms
 implicit none

 integer,intent(in)                 :: m_x,n_x
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 type(spectral_function),intent(in) :: chi
 real(prec_td),intent(in)           :: bigx(m_x,n_x)
 real(prec_td),intent(in)           :: bigy(m_x,n_x)
 real(dp),intent(in)                :: eigenvalue(chi%npole_reso_apb)
!=====
 integer                            :: t_ij,t_kl
 integer                            :: t_ij_global,t_kl_global
 integer                            :: nmat
 integer                            :: istate,jstate,ijspin
 integer                            :: ibf,jbf
 integer                            :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart,ibf_cart,jbf_cart
 integer                            :: iomega,idir,jdir
 integer,parameter                  :: nomega=600
 complex(dp)                        :: omega(nomega)
 real(dp)                           :: coeff
 real(dp)                           :: dynamical_pol(nomega),structure_factor(nomega)
 complex(dpc)                       :: bethe_sumrule
 complex(dpc),allocatable           :: gos_basis(:,:),gos_tmp(:,:),gos_state(:,:,:)
 complex(dpc),allocatable           :: gos_cart(:,:)
 complex(dpc),allocatable           :: residu_left(:)
 real(dp)                           :: qvec(3)
 integer,parameter                  :: nq=0 ! 1000
 integer                            :: iq
 real(dp)                           :: fnq(chi%npole_reso_apb)
 integer,parameter                  :: nv=20
 integer                            :: iv
 real(dp)                           :: stopping(nv)
 real(dp)                           :: vv
!=====


 call start_clock(timing_spectrum)
 !
 ! Calculate the spectrum now
 !

 write(stdout,'(/,a)') ' Calculate the stopping power'

 if (nspin/=1) then
   msg='no nspin/=1 allowed'
   call issue_warning(msg)
   return
 endif

 !
 ! Prepare the precalculated table of coefficients
 call setup_gos_llp()


 !
 ! Calculate the dynamical dipole polarizability
 ! and the static dipole polarizability
 !
 ! Set the frequency mesh
 omega(1)     =0.1_dp ! MAX( 0.0_dp      ,MINVAL(ABS(eigenvalue(:)))-3.00/Ha_eV)
 omega(nomega)=4.0_dp ! MIN(20.0_dp/Ha_eV,MAXVAL(ABS(eigenvalue(:)))+3.00/Ha_eV)
 do iomega=2,nomega-1
   omega(iomega) = omega(1) + ( omega(nomega)-omega(1) ) /REAL(nomega-1,dp) * (iomega-1)
 enddo
 ! Add the broadening
 omega(:) = omega(:) + im * 0.10/Ha_eV
  

 
!TESTINGOK§ call basis_function_dipole(basis%bf(1),basis%bf(14),qvec)
!TESTINGOK§ write(*,*) 'dipole',qvec(:)
!TESTINGOK§ call overlap_basis_function(basis%bf(1),basis%bf(14),qvec(1))
!TESTINGOK§ write(*,*) 'overlap',qvec(1)
!TESTINGOK§ qvec(1)=0.00001
!TESTINGOK§ qvec(2)=0.00000
!TESTINGOK§ qvec(3)=0.00000
!TESTINGOK§ call gos_basis_function(basis%bf(1),basis%bf(14),qvec,bethe_sumrule)
!TESTINGOK§ write(*,*) 'bethe_sumrule',bethe_sumrule / qvec(1)/im
!TESTINGOK§ call die('ENOUGH')

 do iq=1,nq
   qvec(1) = 0.0_dp
   qvec(2) = 0.0_dp
   qvec(3) = iq*0.01_dp

   !
   ! First precalculate all the needed GOS in the basis set
   !
   allocate(gos_basis(basis%nbf,basis%nbf))
   ibf_cart = 1
   ibf      = 1
   do while(ibf_cart<=basis%nbf_cart)
     li      = basis%bf(ibf_cart)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
  
     jbf_cart = 1
     jbf      = 1
     do while(jbf_cart<=basis%nbf_cart)
       lj      = basis%bf(jbf_cart)%am
       nj_cart = number_basis_function_am('CART',lj)
       nj      = number_basis_function_am(basis%gaussian_type,lj)
  
       allocate(gos_cart(ni_cart,nj_cart))
  
  
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           call gos_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),qvec,gos_cart(i_cart,j_cart))
         enddo
       enddo
  
       gos_basis(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE( cart_to_pure(li)%matrix(:,:) ) , &
             MATMUL(  gos_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) )
  
       deallocate(gos_cart)
  
       jbf      = jbf      + nj
       jbf_cart = jbf_cart + nj_cart
     enddo
  
     ibf      = ibf      + ni
     ibf_cart = ibf_cart + ni_cart
   enddo
  
   !
   ! Get the gos oscillator strength on states
   allocate(gos_state(basis%nbf,basis%nbf,nspin))
   allocate(gos_tmp(basis%nbf,basis%nbf))
  
   do ijspin=1,nspin
     gos_tmp(:,:) = MATMUL( TRANSPOSE( c_matrix(:,:,ijspin) ) , gos_basis(:,:) )
  
     gos_state(:,:,ijspin) = MATMUL( gos_tmp(:,:) , c_matrix(:,:,ijspin) )
  
   enddo
   deallocate(gos_basis,gos_tmp)
  
  
   nmat=chi%npole_reso_apb
   allocate(residu_left(chi%npole_reso_apb))
  
   residu_left(:) = 0.0_dp
   do t_ij=1,m_x
     t_ij_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ij)
     istate = chi%transition_table_apb(1,t_ij_global)
     jstate = chi%transition_table_apb(2,t_ij_global)
     ijspin = chi%transition_table_apb(3,t_ij_global)
  
     ! Let use (i <-> j) symmetry to halve the loop
     do t_kl=1,n_x
       t_kl_global = colindex_local_to_global(ipcol_sd,npcol_sd,t_kl)
  
       residu_left(t_kl_global) = residu_left(t_kl_global) &
                    + gos_state(istate,jstate,ijspin) * ( bigx(t_ij,t_kl) + bigy(t_ij,t_kl) ) * SQRT(spin_fact)
     enddo
  
   enddo
   call xsum(residu_left)
  
   deallocate(gos_state)

   fnq(:) = 2.0_dp * ABS( residu_left(:) )**2 * eigenvalue(:) / SUM( qvec(:)**2 )

   write(stdout,*) 'bethe_sumrule',NORM2(qvec(:)),SUM(fnq(:))


  
   dynamical_pol(:) = 0.0_dp
   do t_ij=1,nmat
     dynamical_pol(:) = dynamical_pol(:) &
                       + ABS(residu_left(t_ij))**2 &
                        * ( AIMAG( -1.0_dp  / ( omega(:) - eigenvalue(t_ij) ) ) - AIMAG( -1.0_dp  / ( omega(:) + eigenvalue(t_ij) ) ) )
   enddo
!   !
!   ! Get the structure factor
!   write(999,*) '# qvec',qvec(:)
!   do iomega=1,nomega
!     structure_factor(iomega) = 4.0_dp * pi * REAL(omega(iomega),dp) / c_speedlight * dynamical_pol(iomega) * SUM( qvec(:)**2 )
!     write(999,*) REAL(omega(iomega),dp)*Ha_eV,structure_factor(iomega)
!   enddo
!   write(999,*)


   deallocate(residu_left)

!   write(998,*) SUM( qvec(:)**2 ), fnq(6)

!   do iv=1,nv
!     vv = iv * 0.1_dp
!     do t_ij=1,nmat
!       if( NORM2(qvec) < eigenvalue(t_ij) / vv )   &
!          stopping(iv) = stopping(iv) + 1.0_dp / ( pi * vv**2 )  * fnq(t_ij)  * NORM2(qvec)**2
!     enddo
!
!   enddo


 enddo 

! do iv=1,nv
!   vv = iv * 0.1_dp
!   write(997,*) vv,stopping(iv)
! enddo



end subroutine stopping_power


!=========================================================================
subroutine prepare_tddft(nspin_tddft,basis,c_matrix,occupation,v2rho2,vsigma,v2rhosigma,v2sigma2,wf_r,wf_gradr,rho_gradr)
 use,intrinsic ::  iso_c_binding, only: C_INT,C_DOUBLE
 use m_dft_grid
 use m_basis_set
 use m_hamiltonian
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 integer,intent(in)               :: nspin_tddft
 type(basis_set),intent(in)       :: basis
 real(dp),intent(in)              :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)              :: occupation(basis%nbf,nspin)
 real(dp),allocatable,intent(out) :: v2rho2(:,:)
 real(dp),allocatable,intent(out) :: vsigma(:,:)
 real(dp),allocatable,intent(out) :: v2rhosigma(:,:)
 real(dp),allocatable,intent(out) :: v2sigma2(:,:)
 real(dp),allocatable,intent(out) :: wf_r(:,:,:)
 real(dp),allocatable,intent(out) :: wf_gradr(:,:,:,:)
 real(dp),allocatable,intent(out) :: rho_gradr(:,:,:)
!=====
 real(dp),parameter :: kernel_capping=1.0e14_dp
 integer  :: idft_xc,igrid
 integer  :: ispin
 real(dp) :: rr(3)
 real(dp) :: basis_function_r(basis%nbf)
 real(dp) :: basis_function_gradr(3,basis%nbf)
 real(dp) :: rhor_r(nspin)
 real(dp) :: grad_rhor(3,nspin)
 real(dp) :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp) :: max_v2sigma2
 logical  :: require_gradient
 character(len=256) :: string
#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func(ndft_xc),xc_functest
 type(xc_f90_pointer_t) :: xc_info(ndft_xc),xc_infotest
#endif
 real(C_DOUBLE) :: rho_c(nspin_tddft)
 real(C_DOUBLE) :: v2rho2_c(2*nspin_tddft-1)
 real(C_DOUBLE) :: sigma_c(2*nspin_tddft-1)
 real(C_DOUBLE) :: vrho_c(nspin_tddft)
 real(C_DOUBLE) :: vsigma_c(2*nspin_tddft-1)
 real(C_DOUBLE) :: v2rhosigma_c(5*nspin_tddft-4)
 real(C_DOUBLE) :: v2sigma2_c(5*nspin_tddft-4)
!=====

 !
 ! Prepare DFT kernel calculation with Libxc
 !
 require_gradient  =.FALSE.
 do idft_xc=1,ndft_xc
   if(nspin_tddft==1) then
     call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_UNPOLARIZED)
   else
     call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_POLARIZED)
   endif
   call xc_f90_info_name(xc_info(idft_xc),string)
   write(stdout,'(a,i4,a,i6,5x,a)') '   XC functional ',idft_xc,' :  ',xc_f90_info_number(xc_info(idft_xc)),&
         TRIM(string)
   if( MODULO(xc_f90_info_flags( xc_info(idft_xc)),XC_FLAGS_HAVE_FXC*2) < XC_FLAGS_HAVE_FXC ) then
     call die('This functional does not have the kernel implemented in Libxc')
   endif
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_GGA     ) require_gradient  =.TRUE.
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) require_gradient  =.TRUE.
 enddo

 !
 ! calculate rho, grad rho and the kernel
 ! 
 ! Get the density matrix P from C
 call setup_density_matrix(basis%nbf,c_matrix,occupation,p_matrix)

 allocate(v2rho2(ngrid,2*nspin_tddft-1),wf_r(ngrid,basis%nbf,nspin))
 v2rho2(:,:) = 0.0_dp

 if(require_gradient) then
   allocate(vsigma(ngrid,2*nspin_tddft-1))
   allocate(v2rhosigma(ngrid,5*nspin_tddft-4))
   allocate(v2sigma2(ngrid,5*nspin_tddft-4))
   allocate(wf_gradr(3,ngrid,basis%nbf,nspin))
   allocate(rho_gradr(3,ngrid,nspin))
   vsigma(:,:)     = 0.0_dp
   v2rhosigma(:,:) = 0.0_dp
   v2sigma2(:,:)   = 0.0_dp
 endif


 max_v2sigma2 = -1.0_dp
 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)

   if( .NOT. ALLOCATED(bfr) ) call prepare_basis_functions_r(basis)
   if( require_gradient .AND. .NOT. ALLOCATED(bfgr) ) call prepare_basis_functions_gradr(basis)
   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r(basis,igrid,basis_function_r)
   !
   ! store the wavefunction in r
   do ispin=1,nspin
     wf_r(igrid,:,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,:,ispin) )
   enddo

   if( require_gradient ) then
     call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
     !
     ! store the wavefunction in r
     do ispin=1,nspin
       wf_gradr(:,igrid,:,ispin) = MATMUL( basis_function_gradr(:,:) , c_matrix(:,:,ispin) )
     enddo
   endif


   call calc_density_pmatrix(nspin,basis,p_matrix,rr,basis_function_r,rhor_r)
   if( require_gradient ) then
     call calc_density_gradr_pmatrix(nspin,basis%nbf,p_matrix,basis_function_r,basis_function_gradr,grad_rhor)

     rho_gradr(:,igrid,:) = grad_rhor(:,:)

     if( nspin_tddft==1 ) then
       sigma_c(1) = SUM( grad_rhor(:,1)**2 )
     else if( nspin==2 ) then
       sigma_c(2) = SUM( grad_rhor(:,1) * grad_rhor(:,2) )
       sigma_c(3) = SUM( grad_rhor(:,2)**2 )
     else ! triplet excitations from singlet ground-state
       sigma_c(:) = SUM( grad_rhor(:,1)**2 ) * 0.25_dp
     endif

   endif


   if( nspin_tddft==1 ) then
     rho_c(1) = rhor_r(1)
   else if( nspin==2 ) then
     rho_c(:) = rhor_r(:)
   else ! triplet excitations from singlet ground-state
     rho_c(:) = rhor_r(1) * 0.5_dp
   endif

   !
   ! Calculate the kernel
   ! 
   do idft_xc=1,ndft_xc
     select case(xc_f90_info_family(xc_info(idft_xc)))
     case(XC_FAMILY_LDA)
       call xc_f90_lda_fxc(xc_func(idft_xc),1_C_INT,rho_c(1),v2rho2_c(1))
     case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
       call xc_f90_gga_vxc(xc_func(idft_xc),1_C_INT,rho_c(1),sigma_c(1),vrho_c(1),vsigma_c(1))
       call xc_f90_gga_fxc(xc_func(idft_xc),1_C_INT,rho_c(1),sigma_c(1),v2rho2_c(1),v2rhosigma_c(1),v2sigma2_c(1))
     case default
       call die('Other kernels not yet implemented')
     end select
     !
     ! Remove the too large values for stability
     v2rho2_c(:) = MIN( v2rho2_c(:), kernel_capping )
     if(require_gradient) then
       max_v2sigma2 = MAX(ABS(v2sigma2_c(1)),max_v2sigma2)
       vsigma_c(:)     = MIN( vsigma_c(:), kernel_capping )
       v2rhosigma_c(:) = MIN( v2rhosigma_c(:), kernel_capping )
       v2sigma2_c(:)   = MIN( v2sigma2_c(:), kernel_capping )
     endif

     ! Store the result with the weight
     ! Remove too large values for stability
     v2rho2(igrid,:)     = v2rho2(igrid,:) + v2rho2_c(:) * w_grid(igrid) * dft_xc_coef(idft_xc)
     if(require_gradient) then
       vsigma(igrid,:)     = vsigma(igrid,:)     + vsigma_c(:)     * w_grid(igrid) * dft_xc_coef(idft_xc)
       v2rhosigma(igrid,:) = v2rhosigma(igrid,:) + v2rhosigma_c(:) * w_grid(igrid) * dft_xc_coef(idft_xc)
       v2sigma2(igrid,:)   = v2sigma2(igrid,:)   + v2sigma2_c(:)   * w_grid(igrid) * dft_xc_coef(idft_xc)
     endif

   enddo
 enddo
 if(require_gradient) then
   call xmax(max_v2sigma2)
   write(stdout,'(a,e18.6)') ' Maximum numerical value for fxc: ',max_v2sigma2
 endif


end subroutine prepare_tddft


!=========================================================================
subroutine get_energy_qp(nbf,energy,occupation,energy_qp)
 implicit none

 integer,intent(in)                  :: nbf
 real(dp),intent(in)                 :: energy(nbf,nspin)
 real(dp),intent(in)                 :: occupation(nbf,nspin)
 real(dp),intent(out)                :: energy_qp(nbf,nspin)
!=====
 integer  :: reading_status
 real(dp) :: scissor_energy(nspin)
 integer  :: ispin,istate
!=====

 ! If the keyword scissor is used in the input file,
 ! then use it and ignore the ENERGY_QP file
 if( ABS(scissor) > 1.0e-5_dp ) then

   call issue_warning('Using a manual scissor to open up the fundamental gap')

   write(stdout,'(a,2(1x,f12.6))') ' Scissor operator with value (eV):',scissor*Ha_eV
   do ispin=1,nspin
     do istate=1,nbf
       if( occupation(istate,ispin) > completely_empty/spin_fact ) then
         energy_qp(istate,ispin) = energy(istate,ispin)
       else
         energy_qp(istate,ispin) = energy(istate,ispin) + scissor
       endif
     enddo
   enddo
   write(stdout,'(/,a)') ' Scissor updated energies'
   do istate=1,nbf
     write(stdout,'(i5,4(2x,f16.6))') istate,energy(istate,:)*Ha_eV,energy_qp(istate,:)*Ha_eV
   enddo
   write(stdout,*)

 else

   call read_energy_qp(nspin,nbf,energy_qp,reading_status)

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
subroutine chi_to_vchiv(nbf,nstate,prod_basis,c_matrix,bigx,bigy,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,nstate
 type(basis_set),intent(in)            :: prod_basis
 real(dp),intent(in)                   :: c_matrix(nbf,nstate,nspin)
 type(spectral_function),intent(inout) :: wpol
 real(prec_td),intent(in)              :: bigx(wpol%npole_reso_apb,wpol%npole_reso_apb)
 real(prec_td),intent(in)              :: bigy(wpol%npole_reso_apb,wpol%npole_reso_apb)
 real(dp),intent(in)                   :: eigenvalue(wpol%npole_reso_apb)
!=====
 integer                               :: t_kl,klspin,ijspin
 integer                               :: istate,jstate,kstate,lstate,ijstate,ijstate_spin
 integer                               :: klstate_min
 integer                               :: klstate_max
 integer                               :: nmat,nprodspin
 real(dp)                              :: eri_eigen_klij
 real(dp),allocatable                  :: eri_eigenstate_klmin(:,:,:,:)
 real(dp)                              :: rtmp
!=====

 call start_clock(timing_buildw)

 write(stdout,'(/,a)') ' Build W = v * chi * v'
 if(has_auxil_basis) then
   call die('you should not be here')
 endif

 allocate(eri_eigenstate_klmin(nbf,nbf,nbf,nspin))
 ! Set this to zero and then enforce the calculation of the first array of Coulomb integrals
 eri_eigenstate_klmin(:,:,:,:) = 0.0_dp

 call allocate_spectral_function(prod_basis%nbf*nspin,wpol)

 wpol%pole(1:wpol%npole_reso_apb) = eigenvalue(:)

 nmat = wpol%npole_reso_apb

 wpol%residu_left(:,:) = 0.0_dp
 do t_kl=1,nmat 
   kstate = wpol%transition_table_apb(1,t_kl)
   lstate = wpol%transition_table_apb(2,t_kl)
   klspin = wpol%transition_table_apb(3,t_kl)

   klstate_min = MIN(kstate,lstate)
   klstate_max = MAX(kstate,lstate)
   call calculate_eri_4center_eigen(nbf,nstate,c_matrix,klstate_min,klspin,eri_eigenstate_klmin)

   do ijspin=1,nspin
     do ijstate=1,prod_basis%nbf
       istate = prod_basis%index_ij(1,ijstate)
       jstate = prod_basis%index_ij(2,ijstate)

       ijstate_spin = ijstate+prod_basis%nbf*(ijspin-1)

       eri_eigen_klij = eri_eigenstate_klmin(klstate_max,istate,jstate,ijspin)

       ! Use the symmetry ( k l | i j ) to regroup (kl) and (lk) contributions
       ! and the block structure of eigenvector | X  Y |
       !                                        | Y  X |
       wpol%residu_left(ijstate_spin,:) = wpol%residu_left(ijstate_spin,:) &
                              + eri_eigen_klij * ( bigx(t_kl,:) + bigy(t_kl,:) )

     enddo
   enddo
 enddo

 wpol%residu_left(:,:) = wpol%residu_left(:,:) * SQRT(spin_fact)



 if(ALLOCATED(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)

 call stop_clock(timing_buildw)

end subroutine chi_to_vchiv


!=========================================================================
subroutine chi_to_sqrtvchisqrtv_auxil(nbf,nbf_auxil,desc_x,m_x,n_x,bigx,bigy,eigenvalue,wpol,energy_gm)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,nbf_auxil,m_x,n_x
 integer,intent(in)                    :: desc_x(ndel)
 real(prec_td),intent(inout)           :: bigx(m_x,n_x)
 real(prec_td),intent(in)              :: bigy(m_x,n_x)
 type(spectral_function),intent(inout) :: wpol
 real(dp),intent(in)                   :: eigenvalue(wpol%npole_reso_apb)
 real(dp),intent(out)                  :: energy_gm
!=====
 integer                               :: t_kl,t_kl_global,klspin
 integer                               :: t_ij,t_ij_global,ijspin
 integer                               :: ibf_auxil,ibf_auxil_global
 integer                               :: nmat
 integer                               :: kstate,lstate
 integer                               :: istate,jstate
 real(dp),allocatable                  :: eri_3center_mat(:,:),residu_local(:,:)
 integer                               :: desc_3center_eigen(ndel)
 integer                               :: desc_residu(ndel)
 integer                               :: m_3center,n_3center
 real(dp)                              :: rtmp
 integer                               :: iprow,ipcol
 integer                               :: m_bigx_block,n_bigx_block
 real(dp),allocatable                  :: bigx_block(:,:)
!=====

 call start_clock(timing_buildw)

 write(stdout,'(/,a)') ' Build v^{1/2} * chi * v^{1/2}'

 call allocate_spectral_function(nbf_auxil,wpol)
 wpol%pole(1:wpol%npole_reso_apb) = eigenvalue(:)

 nmat = wpol%npole_reso_apb

#ifndef HAVE_SCALAPACK

 allocate(eri_3center_mat(nbf_auxil,nmat))
 do t_kl=1,nmat
   kstate = wpol%transition_table_apb(1,t_kl)
   lstate = wpol%transition_table_apb(2,t_kl)
   klspin = wpol%transition_table_apb(3,t_kl)
   eri_3center_mat(:,t_kl) = eri_3center_eigen(:,kstate,lstate,klspin)
 end do

 ! Use the symmetry ( I | k l ) to regroup (kl) and (lk) contributions
 ! and the block structure of eigenvector | X  Y |
 !                                        | Y  X |
 wpol%residu_left(:,:) = MATMUL( eri_3center_mat , bigx(:,:) + bigy(:,:) ) * SQRT(spin_fact)

 energy_gm = 0.5_dp * ( SUM( wpol%residu_left(:,:)**2 ) - spin_fact * SUM( eri_3center_mat(:,:)**2 ) )
 !
 ! Since wpol%residu_left and eri_3center_mat are distributed, we have to sum up
 call xsum(energy_gm)

 deallocate(eri_3center_mat)

#else 


 ! bigx = bigx + bigy 
 call PDGEADD('N',nmat,nmat,1.d0,bigy,1,1,desc_x,1.d0,bigx,1,1,desc_x)

 bigx(:,:) = bigx(:,:) * SQRT(spin_fact)

 wpol%residu_left(:,:) = 0.0_dp
 ! First loop over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_bigx_block = row_block_size(nmat,iprow,nprow_sd)
     n_bigx_block = col_block_size(nmat,ipcol,npcol_sd)

     allocate(bigx_block(m_bigx_block,n_bigx_block))
     if( ipcol == ipcol_sd .AND. iprow == iprow_sd ) then
       bigx_block(:,:) = bigx(:,:)
     else
       bigx_block(:,:) = 0.0_dp
     endif
     call xsum(bigx_block)

     do t_kl=1,n_bigx_block
       t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)

       do t_ij=1,m_bigx_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)

         wpol%residu_left(:,t_kl_global) = wpol%residu_left(:,t_kl_global) &
                  + eri_3center_eigen(:,istate,jstate,ijspin) * bigx_block(t_ij,t_kl)

       enddo
     enddo


     deallocate(bigx_block)
   enddo
 enddo

 energy_gm = 0.0_dp
 do t_ij_global=1,nmat
   istate = wpol%transition_table_apb(1,t_ij_global)
   jstate = wpol%transition_table_apb(2,t_ij_global)
   ijspin = wpol%transition_table_apb(3,t_ij_global)
   energy_gm = energy_gm - SUM( eri_3center_eigen(:,istate,jstate,ijspin)**2 ) * spin_fact * 0.5_dp
 enddo

 energy_gm = energy_gm + 0.5_dp * ( SUM( wpol%residu_left(:,:)**2 ) )
 call xsum(energy_gm)


#endif



 call stop_clock(timing_buildw)

end subroutine chi_to_sqrtvchisqrtv_auxil


!=========================================================================
subroutine chi_to_sqrtvchisqrtv_auxil_spa(nbf,nbf_auxil,a_diag,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,nbf_auxil
 type(spectral_function),intent(inout) :: wpol
 real(dp),intent(in)                   :: a_diag(wpol%npole_reso_spa)
!=====
 integer                               :: t_kl,klspin
 integer                               :: kstate,lstate
!=====

 call start_clock(timing_buildw)

 write(stdout,'(/,a)') ' Build v^{1/2} * chi * v^{1/2} part from single pole approximation'

 wpol%pole(wpol%npole_reso_apb+1:wpol%npole_reso) = a_diag(:)

 do t_kl=1,wpol%npole_reso_spa
   kstate = wpol%transition_table_spa(1,t_kl)
   lstate = wpol%transition_table_spa(2,t_kl)
   klspin = wpol%transition_table_spa(3,t_kl)


   wpol%residu_left(:,wpol%npole_reso_apb+t_kl) = eri_3center_eigen(:,kstate,lstate,klspin) * SQRT(spin_fact) 

 end do


 call stop_clock(timing_buildw)

end subroutine chi_to_sqrtvchisqrtv_auxil_spa


end module m_timedependent
!=========================================================================
