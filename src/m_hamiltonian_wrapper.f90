!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
#include "molgw.h"
module m_hamiltonian_wrapper
  use m_definitions
  use m_timing
  use m_mpi
  use m_scalapack
  use m_warning
  use m_inputparam,only: has_auxil_basis,incore_
  use m_basis_set
  use m_hamiltonian_twobodies
  use m_hamiltonian_tools
  use m_scf

  interface calculate_exchange
    module procedure calculate_exchange_real
  end interface

contains


!=========================================================================
subroutine calculate_hartree(basis,p_matrix,hhartree,eh)
  implicit none
  type(basis_set),intent(inout) :: basis
  class(*),intent(in)           :: p_matrix(:,:,:)
  real(dp),intent(out)          :: hhartree(:,:)
  real(dp),intent(out),optional :: eh
  !=====
  real(dp) :: ehartree
  real(dp),allocatable :: rho_coeff(:,:)
  !=====


  !
  if( .NOT. has_auxil_basis ) then
    select type(p_matrix)
    type is(real(dp))
      if( incore_ ) then
        call setup_hartree(p_matrix,hhartree,ehartree)
      else
        call setup_hartree_oneshell(basis,p_matrix,hhartree,ehartree)
      endif
    class default
      call die('calculate_hartree: should not happen')
    end select
  else
    if( .NOT. eri3_genuine_ ) then
      call setup_hartree_ri(p_matrix,hhartree,ehartree)
    else
      call calculate_density_auxilbasis(p_matrix,rho_coeff)
      call setup_hartree_genuine_ri(p_matrix,rho_coeff,hhartree,ehartree)
      deallocate(rho_coeff)
    endif
  endif

  if( PRESENT(eh) ) eh = ehartree

end subroutine calculate_hartree


!=========================================================================
subroutine calculate_exchange_real(basis,p_matrix,hexx,ex,occupation,c_matrix)
  implicit none
  type(basis_set),intent(in)    :: basis
  real(dp),intent(in)           :: p_matrix(:,:,:)
  real(dp),intent(out)          :: hexx(:,:,:)
  real(dp),intent(out),optional :: ex
  real(dp),intent(in),optional  :: occupation(:,:)
  real(dp),intent(in),optional  :: c_matrix(:,:,:)
  !=====
  real(dp),allocatable :: c_matrix_tmp(:,:,:)
  real(dp),allocatable :: occupation_tmp(:,:)
  real(dp)             :: eexx
  !=====


  if( .NOT. has_auxil_basis ) then
    if( incore_ ) then
      call setup_exchange(p_matrix,hexx,eexx)
    else
      call issue_warning('no out-of-core exchange implemented yet')
      hexx(:,:,:) = 0.0_dp
      eexx = 0.0_dp
    endif
  else
    if( PRESENT(occupation) .AND. PRESENT(c_matrix) ) then
      call setup_exchange_ri(occupation,c_matrix,p_matrix,hexx,eexx)
    else
      !
      ! c_matrix is not provided, then calculate it from the square-root of P
      call get_c_matrix_from_p_matrix(p_matrix,c_matrix_tmp,occupation_tmp)
      call setup_exchange_ri(occupation_tmp,c_matrix_tmp,p_matrix,hexx,eexx)
      deallocate(c_matrix_tmp)
      deallocate(occupation_tmp)

    endif
  endif

  if( PRESENT(ex) ) ex = eexx

end subroutine calculate_exchange_real


!=========================================================================
subroutine calculate_exchange_lr(basis,p_matrix,hexx,ex,occupation,c_matrix)
  implicit none
  type(basis_set),intent(in)    :: basis
  real(dp),intent(in)           :: p_matrix(:,:,:)
  real(dp),intent(out)          :: hexx(:,:,:)
  real(dp),intent(out),optional :: ex
  real(dp),intent(in),optional  :: occupation(:,:)
  real(dp),intent(in),optional  :: c_matrix(:,:,:)
  !=====
  real(dp),allocatable :: c_matrix_tmp(:,:,:)
  real(dp),allocatable :: occupation_tmp(:,:)
  real(dp)             :: eexx
  !=====


  if( .NOT. has_auxil_basis ) then
    call setup_exchange_longrange(p_matrix,hexx,eexx)
  else
    if( PRESENT(occupation) .AND. PRESENT(c_matrix) ) then
      call setup_exchange_longrange_ri(occupation,c_matrix,p_matrix,hexx,eexx)
    else
      !
      ! c_matrix is not provided, then calculate it from the square-root of P
      call get_c_matrix_from_p_matrix(p_matrix,c_matrix_tmp,occupation_tmp)
      call setup_exchange_longrange_ri(occupation_tmp,c_matrix_tmp,p_matrix,hexx,eexx)
      deallocate(c_matrix_tmp)
      deallocate(occupation_tmp)

    endif
  endif

  if( PRESENT(ex) ) ex = eexx

end subroutine calculate_exchange_lr


!=========================================================================
subroutine calculate_hamiltonian_hxc(basis,nstate,occupation,c_matrix,p_matrix,hamiltonian_hxc,en_inout)
  implicit none

  type(basis_set),intent(inout) :: basis
  integer,intent(in)            :: nstate
  real(dp),intent(in)           :: occupation(nstate,nspin)
  real(dp),intent(in)           :: c_matrix(:,:,:)
  real(dp),intent(in)           :: p_matrix(:,:,:)
  real(dp),intent(out)          :: hamiltonian_hxc(:,:,:)
  type(energy_contributions),intent(inout) :: en_inout
  !=====
  integer              :: ispin
  real(dp),allocatable :: hamiltonian_nospin_real(:,:)
  real(dp),allocatable :: hamiltonian_spin_real(:,:,:)
  real(dp)             :: exx_lr
  !=====

  en_inout%hartree = 0.0_dp
  en_inout%xc      = 0.0_dp
  en_inout%exx     = 0.0_dp
  en_inout%exx_hyb = 0.0_dp
  hamiltonian_hxc(:,:,:) = 0.0_dp

  !
  ! For a core only calculation, no need to go any further
  ! no hartree, no exchange-correlation
  if( calc_type%is_core ) return


  allocate(hamiltonian_nospin_real,MOLD=hamiltonian_hxc(:,:,1))
  allocate(hamiltonian_spin_real,MOLD=hamiltonian_hxc(:,:,:))

  !
  ! Hartree contribution to the Hamiltonian
  !
  call calculate_hartree(basis,p_matrix,hamiltonian_nospin_real,eh=en_inout%hartree)

  do ispin=1,nspin
    hamiltonian_hxc(:,:,ispin) = hamiltonian_nospin_real(:,:)
  enddo


  !
  !  XC part of the Hamiltonian
  !

  !
  ! DFT XC potential is added here
  !
  if( calc_type%is_dft ) then
    hamiltonian_spin_real(:,:,:) = 0.0_dp

    call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix,hamiltonian_spin_real,en_inout%xc)

    hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_real(:,:,:)
  endif


  !
  ! LR Exchange contribution to the Hamiltonian
  !
  if(calc_type%need_exchange_lr) then
    hamiltonian_spin_real(:,:,:) = 0.0_dp

    call calculate_exchange_lr(basis,p_matrix,hamiltonian_spin_real,ex=exx_lr,occupation=occupation,c_matrix=c_matrix)
    ! Rescale with beta_hybrid for range-separated hybrid functionals
    en_inout%exx_hyb = beta_hybrid * exx_lr
    hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_real(:,:,:) * beta_hybrid

  endif


  !
  ! Exchange contribution to the Hamiltonian
  !
  if( calc_type%need_exchange ) then
    hamiltonian_spin_real(:,:,:) = 0.0_dp

    call calculate_exchange(basis,p_matrix,hamiltonian_spin_real,ex=en_inout%exx,occupation=occupation,c_matrix=c_matrix)
    ! Rescale with alpha_hybrid for hybrid functionals
    en_inout%exx_hyb = en_inout%exx_hyb + alpha_hybrid * en_inout%exx
    hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_real(:,:,:) * alpha_hybrid

  endif

  deallocate(hamiltonian_spin_real)
  deallocate(hamiltonian_nospin_real)


end subroutine calculate_hamiltonian_hxc


!=========================================================================
subroutine calculate_hamiltonian_hxc_ri_cmplx(basis,                  &
                                              occupation,             &
                                              c_matrix_cmplx,         &
                                              p_matrix_cmplx,         &
                                              hamiltonian_hxc_cmplx,  &
                                              en_inout)
  implicit none

  type(basis_set),intent(inout) :: basis
  real(dp),intent(in)           :: occupation(:,:)
  complex(dp),intent(in)        :: c_matrix_cmplx(:,:,:)
  complex(dp),intent(in)        :: p_matrix_cmplx(:,:,:)
  complex(dp),intent(out)       :: hamiltonian_hxc_cmplx(:,:,:)
  type(energy_contributions),intent(inout) :: en_inout
  !=====
  integer                    :: nstate
  integer                    :: ispin
  real(dp),allocatable       :: hamiltonian_nospin_real(:,:)
  real(dp),allocatable       :: hamiltonian_spin_real(:,:,:)
  complex(dp),allocatable    :: hamiltonian_spin_cmplx(:,:,:)
  !=====

  en_inout%hartree = 0.0_dp
  en_inout%xc      = 0.0_dp
  en_inout%exx     = 0.0_dp
  en_inout%exx_hyb = 0.0_dp

  nstate = SIZE(occupation,DIM=1)

  ! Initialize real arrays
  hamiltonian_hxc_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )

  !
  ! For a core only calculation, no need to go any further
  ! no hartree, no exchange-correlation
  if( calc_type%is_core ) return

  !
  ! Exchange contribution to the Hamiltonian
  !
  if( calc_type%need_exchange ) then
    call setup_exchange_ri_cmplx(occupation,c_matrix_cmplx,p_matrix_cmplx,hamiltonian_hxc_cmplx,en_inout%exx)

    hamiltonian_hxc_cmplx(:,:,:) = hamiltonian_hxc_cmplx(:,:,:) * alpha_hybrid
  endif

  allocate(hamiltonian_nospin_real(basis%nbf,basis%nbf))
  !
  ! Hartree contribution to the Hamiltonian
  ! Hartree contribution is real and depends only on real(p_matrix) but we pass the full p_matrix_cmplx any way
  !
  call calculate_hartree(basis,p_matrix_cmplx,hamiltonian_nospin_real,eh=en_inout%hartree)

  do ispin=1,nspin
    hamiltonian_hxc_cmplx(:,:,ispin) = hamiltonian_hxc_cmplx(:,:,ispin) + hamiltonian_nospin_real(:,:)
  enddo
  deallocate(hamiltonian_nospin_real)

  !
  !  XC part of the Hamiltonian
  !

  !
  ! DFT XC potential is added here
  !
  if( calc_type%is_dft ) then
    allocate(hamiltonian_spin_real(basis%nbf,basis%nbf,nspin))
    call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix_cmplx,hamiltonian_spin_real,en_inout%xc)

    hamiltonian_hxc_cmplx(:,:,:) = hamiltonian_hxc_cmplx(:,:,:) + hamiltonian_spin_real(:,:,:)
    deallocate(hamiltonian_spin_real)
  endif

  !
  ! LR Exchange contribution to the Hamiltonian
  !
  if(calc_type%need_exchange_lr) then
    allocate(hamiltonian_spin_cmplx(basis%nbf,basis%nbf,nspin))
    hamiltonian_spin_cmplx(:,:,:) = 0.0_dp
  
    call setup_exchange_longrange_ri_cmplx(occupation,c_matrix_cmplx,p_matrix_cmplx,hamiltonian_spin_cmplx,en_inout%exx_hyb)
  
    hamiltonian_hxc_cmplx(:,:,:) = hamiltonian_hxc_cmplx(:,:,:) + hamiltonian_spin_cmplx(:,:,:) * beta_hybrid
    deallocate(hamiltonian_spin_cmplx)
  else
    en_inout%exx_hyb = 0.0_dp
  endif

  ! Rescale exchange energy with alpha_hybrid and beta_hybrid for hybrid functionals
  en_inout%exx_hyb = alpha_hybrid * en_inout%exx + beta_hybrid * en_inout%exx_hyb



end subroutine calculate_hamiltonian_hxc_ri_cmplx


!=========================================================================
end module m_hamiltonian_wrapper
!=========================================================================
