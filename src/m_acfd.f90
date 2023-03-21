!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval, Mauricio Rodriguez-Mayorga
!
! This module contains
! the routines to calculate the ACFD total energy
!
!=========================================================================
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
  use m_spectral_function
  use m_spectra
  use m_scf
  use m_linear_response




contains


subroutine acfd_total_energy(basis,nstate,occupation,energy,c_matrix,en_mbpt)
  type(basis_set),intent(in)               :: basis
  integer,intent(in)                       :: nstate
  real(dp),intent(in)                      :: occupation(:,:),energy(:,:)
  real(dp),intent(in)                      :: c_matrix(:,:,:)
  type(energy_contributions),intent(inout) :: en_mbpt
  !=====
  type(spectral_function)    :: wpol
  real(dp)                   :: erpa_singlet,erpa_triplet,egw_tmp
  !=====

  call init_spectral_function(nstate,occupation,0,wpol)
  call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol,enforce_spin_multiplicity=1)
  call destroy_spectral_function(wpol)
  if( TRIM(postscf) == 'RPA' ) then
    en_mbpt%rpa = erpa_singlet

    if(kappa_hybrid/=zero) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' RPA Energy scaled by :',kappa_hybrid
      en_mbpt%rpa=kappa_hybrid*en_mbpt%rpa
    endif
    write(stdout,'(a,2x,f19.10)') ' RPA Energy      (Ha):',en_mbpt%rpa

    if(kappa_hybrid/=zero) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    endif

    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPA Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  else
    en_mbpt%rpa = 0.50_dp * erpa_singlet
    if(kappa_hybrid/=zero) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' Singlet RPAx Energy scaled by :',kappa_hybrid
      en_mbpt%rpa=kappa_hybrid*en_mbpt%rpa
    endif
    write(stdout,'(a,2x,f19.10)') ' Singlet RPAx Energy contribution      (Ha):',en_mbpt%rpa

    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_triplet,egw_tmp,wpol,enforce_spin_multiplicity=3)
    call destroy_spectral_function(wpol)
    if(kappa_hybrid/=zero) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout,'(/,a,f16.10)') ' Triplet RPAx Energy scaled by :',kappa_hybrid
      erpa_triplet=kappa_hybrid*erpa_triplet
    endif
    write(stdout,'(a,2x,f19.10)') ' Triplet RPAx Energy contribution      (Ha):',1.50_dp * erpa_triplet
    en_mbpt%rpa = en_mbpt%rpa + 1.50_dp * erpa_triplet
    write(stdout,'(a,2x,f19.10)') ' RPAx Energy      (Ha):',en_mbpt%rpa

    if(kappa_hybrid/=zero) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    endif

    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPAx Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  endif
  call print_energy_yaml('mbpt energy',en_mbpt)

end subroutine acfd_total_energy


end module m_acfd
!=========================================================================
