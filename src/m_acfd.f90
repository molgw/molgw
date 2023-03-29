!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval, Mauricio Rodriguez-Mayorga
!
! This module contains
! the routines to calculate the ACFD total energy
!
!=========================================================================
#include "molgw.h"
module m_acfd
  use m_definitions
  use m_timing
  use m_warning
  use m_memory
  use m_inputparam
  use m_mpi
  use m_scalapack
  use m_cart_to_pure
  use m_basis_set
  use m_dft_grid
  use m_hamiltonian_twobodies
  use m_spectral_function
  use m_spectra
  use m_scf
  use m_gw_selfenergy_grid
  use m_linear_response
  use m_numerical_tools,only: coeffs_gausslegint

#if !defined(NO_LIBXC)
#include <xc_funcs.h>
#endif

contains


!=========================================================================
subroutine acfd_total_energy(basis,nstate,occupation,energy,c_matrix,en_mbpt)
  type(basis_set),intent(in)  :: basis
  integer,intent(in)          :: nstate
  real(dp),intent(in)         :: occupation(:,:),energy(:,:)
  real(dp),intent(in)         :: c_matrix(:,:,:)
  real(dp),allocatable        :: matrix_tmp(:,:,:)
  type(energy_contributions),intent(inout) :: en_mbpt
  !=====
  type(spectral_function)    :: wpol
  real(dp)                   :: erpa_singlet,erpa_triplet,egw_tmp,erpa_sie_KP,erpa_state
  real(dp)                   :: wlambda(acfd_nlambda),lambda(acfd_nlambda)
  real(dp),allocatable       :: x_matrix(:,:),y_matrix(:,:)
  real(dp),allocatable       :: a_matrix(:,:),b_matrix(:,:)
  integer                    :: nmat,desc_x(NDEL)
  !=====

  !
  ! First the "+" part of RPA+ 
  !
  if( TRIM(postscf) == 'RPAP' .OR. TRIM(postscf) == 'RPAP_IM' ) then
    if( ALLOCATED(dft_xc) ) then
       write(stdout,'(/,a,/)') ' Deallocate dft_xc object before re-allocating it for RPA+ correction.'
       call destroy_libxc_info(dft_xc)
    endif
    allocate(dft_xc(2))
    dft_xc(:)%id = 0
    dft_xc(:)%nspin = nspin
    dft_xc(1)%id = XC_LDA_C_PW      ! HEG
    dft_xc(2)%id = XC_LDA_C_PW_RPA  ! RPA@HEG
    dft_xc(1)%coeff = one
    dft_xc(2)%coeff = -one
    call init_libxc_info(dft_xc)
    call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
    call clean_allocate('XC operator RPA+',matrix_tmp,basis%nbf,basis%nbf,nspin)
    call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix,matrix_tmp,erpa_sie_KP)
    call destroy_dft_grid()
    call clean_deallocate('XC operator RPA+',matrix_tmp)
    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' RPA+ Energy scaled by :',kappa_hybrid
      erpa_sie_KP = kappa_hybrid * erpa_sie_KP
    else
      erpa_sie_KP = 0.0_dp
    endif
    write(stdout,'(/,a,f19.10)') ' E(LHEG) - E(LRPA) correlation energy (Ha):',erpa_sie_KP
  endif
  



  select case(TRIM(postscf))
  case('RPA','RPAP')
    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol,enforce_spin_multiplicity=1)
    call destroy_spectral_function(wpol)
    en_mbpt%rpa = erpa_singlet
    write(stdout,'(a,2x,f19.10)') ' RPA Energy      (Ha):',en_mbpt%rpa

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' RPA Energy scaled by :',kappa_hybrid
      en_mbpt%rpa = kappa_hybrid * en_mbpt%rpa
    endif

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa + erpa_sie_KP 
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus &
                    + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa + erpa_sie_KP
    endif

    if( TRIM(postscf) == 'RPAP' ) then
      write(stdout,'(/,a,f19.10)') ' RPA+ correlation energy (Ha):',en_mbpt%rpa + erpa_sie_KP
      write(stdout,'(a,2x,f19.10)') ' RPA+ Total Energy (Ha):',en_mbpt%total
    else
      write(stdout,'(a,2x,f19.10)') ' RPA Total Energy (Ha):',en_mbpt%total
    endif

  case('RPA_IM','RPAP_IM')
    call init_spectral_function(nstate,occupation,nomega_chi_imag,wpol)
    call polarizability_grid_scalapack(basis,occupation,energy,c_matrix,erpa_state,egw_tmp,wpol)
    call destroy_spectral_function(wpol)
    en_mbpt%rpa = erpa_state

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' RPA Energy scaled by :',kappa_hybrid
      en_mbpt%rpa = kappa_hybrid * en_mbpt%rpa
    endif

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa + erpa_sie_KP 
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus &
                    + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa + erpa_sie_KP
    endif

    if( TRIM(postscf) == 'RPAP_IM' ) then
      write(stdout,'(/,a,f19.10)') ' RPA+ correlation energy (Ha):',en_mbpt%rpa + erpa_sie_KP
    endif
    write(stdout,'(a,2x,f19.10)') ' RPA Total Energy (Ha):',en_mbpt%total

  case('RPAX','RPAX-II')
    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol,enforce_spin_multiplicity=1)
    call destroy_spectral_function(wpol)
    en_mbpt%rpa = 0.50_dp * erpa_singlet

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' Singlet RPAx Energy scaled by :',kappa_hybrid
      en_mbpt%rpa=kappa_hybrid*en_mbpt%rpa
    endif
    write(stdout,'(a,2x,f19.10)') ' Singlet RPAx Energy contribution      (Ha):',en_mbpt%rpa

    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_triplet,egw_tmp,wpol,enforce_spin_multiplicity=3)
    call destroy_spectral_function(wpol)
    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' Triplet RPAx Energy scaled by :',kappa_hybrid
      erpa_triplet=kappa_hybrid*erpa_triplet
    endif
    write(stdout,'(a,2x,f19.10)') ' Triplet RPAx Energy contribution      (Ha):',1.50_dp * erpa_triplet
    en_mbpt%rpa = en_mbpt%rpa + 1.50_dp * erpa_triplet
    write(stdout,'(a,2x,f19.10)') ' RPAx Energy      (Ha):',en_mbpt%rpa

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    endif

    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPAx Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  case('RPA-I')
    call die('acfd_total_energy: RPA-I formula not yet tested')
    write(stdout,'(/,1x,a,i4)') 'RPA with integration over lambda with acfd_nlambda: ',acfd_nlambda
    call coeffs_gausslegint(0.0_dp,1.0_dp,lambda,wlambda,acfd_nlambda)

    call init_spectral_function(nstate,occupation,0,wpol)
    nmat = wpol%npole_reso
    m_x = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
    n_x = NUMROC(nmat,block_col,ipcol_sd,first_col,npcol_sd)
    call DESCINIT(desc_x,nmat,nmat,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,m_x),info)

    call clean_allocate('X matrix',x_matrix,m_x,n_x)
    call clean_allocate('Y matrix',y_matrix,m_x,n_x)

    call clean_allocate('A matrix',a_matrix,m_x,n_x)
    call clean_allocate('B matrix',b_matrix,m_x,n_x)
    ! Get A and B
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                        enforce_spin_multiplicity=1,lambda=1.0_dp,a_matrix=a_matrix,b_matrix=b_matrix)
    call destroy_spectral_function(wpol)

    en_mbpt%rpa = 0.0_dp
    do ilambda=1,acfd_nlambda
      write(stdout,'(1x,a,i4,a,i4)') '=== Lambda',ilambda,' / ',acfd_nlambda
      call init_spectral_function(nstate,occupation,0,wpol)
      call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                          enforce_spin_multiplicity=1,lambda=lambda(ilambda),x_matrix=x_matrix,y_matrix=y_matrix)
      call destroy_spectral_function(wpol)

      call calculate_ec_acft(desc_x,a_matrix,b_matrix,x_matrix,y_matrix,erpa_singlet)

      en_mbpt%rpa = en_mbpt%rpa + erpa_singlet * wlambda(ilambda)
    enddo
    if(has_auxil_basis) call destroy_eri_3center_eigen()

    call clean_deallocate('X matrix',x_matrix)
    call clean_deallocate('Y matrix',y_matrix)
    call clean_deallocate('A matrix',a_matrix)
    call clean_deallocate('B matrix',b_matrix)

    write(stdout,'(a,2x,f19.10)') ' RPA Energy      (Ha):',en_mbpt%rpa

    en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPA Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  case('RPAX-I')
#if defined(HAVE_SCALAPACK)
    call die('acfd_total_energy: RPAX-I formula only available in sequential. Please recompile the code.')
#endif
    write(stdout,'(/,1x,a,i4)') 'RPAx with integration over lambda with acfd_nlambda: ',acfd_nlambda
    call coeffs_gausslegint(0.0_dp,1.0_dp,lambda,wlambda,acfd_nlambda)

    call init_spectral_function(nstate,occupation,0,wpol)
    nmat = wpol%npole_reso
    m_x = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
    n_x = NUMROC(nmat,block_col,ipcol_sd,first_col,npcol_sd)
    call DESCINIT(desc_x,nmat,nmat,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,m_x),info)

    call clean_allocate('X matrix',x_matrix,m_x,n_x)
    call clean_allocate('Y matrix',y_matrix,m_x,n_x)

    call clean_allocate('A matrix',a_matrix,m_x,n_x)
    call clean_allocate('B matrix',b_matrix,m_x,n_x)
    !
    ! Only singlet needed for RPAx-I
    !
    ! Get A and B
    call polarizability(.TRUE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                        enforce_spin_multiplicity=1,lambda=1.0_dp,a_matrix=a_matrix,b_matrix=b_matrix)
    call destroy_spectral_function(wpol)

    en_mbpt%rpa = 0.0_dp
    do ilambda=1,acfd_nlambda
      write(stdout,'(1x,a,i4,a,i4)') '=== Lambda',ilambda,' / ',acfd_nlambda
      write(stdout,'(1x,a,f12.8)')   'lambda: ',lambda(ilambda)
      call init_spectral_function(nstate,occupation,0,wpol)
      call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                          enforce_spin_multiplicity=1,lambda=lambda(ilambda),x_matrix=x_matrix,y_matrix=y_matrix)
      call destroy_spectral_function(wpol)

      call calculate_ec_acft(desc_x,a_matrix,b_matrix,x_matrix,y_matrix,erpa_singlet)

      en_mbpt%rpa = en_mbpt%rpa + erpa_singlet * wlambda(ilambda)
    enddo
    call destroy_spectral_function(wpol)

    if(has_auxil_basis) call destroy_eri_3center_eigen()

    call clean_deallocate('X matrix',x_matrix)
    call clean_deallocate('Y matrix',y_matrix)
    call clean_deallocate('A matrix',a_matrix)
    call clean_deallocate('B matrix',b_matrix)

    write(stdout,'(a,2x,f19.10)') ' RPAx-I Energy      (Ha):',en_mbpt%rpa

    en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPAx-I Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  case default
    call die('acfd_total_energy: postscf option not recognized')
  end select

  call print_energy_yaml('mbpt energy',en_mbpt)

end subroutine acfd_total_energy


!=========================================================================
subroutine calculate_ec_acft(desc_x,a_matrix,b_matrix,x_matrix,y_matrix,erpa)
  implicit none

  integer,intent(in)   :: desc_x(NDEL)
  real(dp),intent(in)  :: x_matrix(:,:),y_matrix(:,:),a_matrix(:,:),b_matrix(:,:)
  real(dp),intent(out) :: erpa
  !=====
  integer :: imat,nmat
  real(dp),allocatable :: m_matrix(:,:)
  !=====

  allocate(m_matrix,MOLD=a_matrix)
  nmat=SIZE(a_matrix,DIM=1)
  m_matrix(:,:) = -a_matrix(:,:)
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(x_matrix) , MATMUL(a_matrix,x_matrix) )
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(y_matrix) , MATMUL(a_matrix,y_matrix) )
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(x_matrix) , MATMUL(b_matrix,y_matrix) )
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(y_matrix) , MATMUL(b_matrix,x_matrix) )

  erpa = 0.0_dp
  do imat=1,nmat
    erpa = erpa + 0.5_dp * m_matrix(imat,imat)
  enddo
  deallocate(m_matrix)


end subroutine calculate_ec_acft


end module m_acfd
!=========================================================================
