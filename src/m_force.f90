!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the calculation of the forces (requires LIBCINT or LIBINT with gradients)
!
!=========================================================================
#include "molgw.h"
module m_force
  use m_definitions
  use m_warning
  use m_timing
  use m_atoms
  use m_inputparam
  use m_basis_set
  use m_eri
  use m_eri_calculate
  use m_hamiltonian_tools
  use m_hamiltonian_onebody
  use m_linear_algebra, only: diagonalize
  use m_dft_grid
  use m_hamiltonian_twobodies, only: dft_exc_vxc_batch
  use m_libcint_tools


  implicit none

contains


!=========================================================================
subroutine calculate_force(basis, auxil_basis, occupation, energy, c_matrix)

  type(basis_set), intent(inout) :: basis
  type(basis_set), intent(inout) :: auxil_basis
  real(dp), intent(in)           :: occupation(:, :)
  real(dp), intent(in)           :: energy(:, :)
  real(dp), intent(in)           :: c_matrix(:, :, :)
  !=====
  integer                 :: nstate
  integer                 :: iatom
  integer                 :: ibf
  real(dp), allocatable    :: r_matrix(:, :)
  real(dp), allocatable    :: p_matrix(:, :, :)
  real(dp), allocatable    :: grad_onebody(:, :, :)
  real(dp), allocatable    :: grad_nucleus(:, :, :, :)
  real(dp), allocatable    :: vxc_ao(:, :, :)
  real(dp)                 :: exc_xc
  !=====

  if( .NOT. MOLGW_has_gradient) then
    call issue_warning('calculate_force: impossible to calculate gradient without LIBCINT or LIBINT with gradients')
    return
  end if

  call timer_force%start()

  nstate = SIZE(energy, DIM=1)

  write(stdout, '(/,1x,a)') 'Calculate the forces'

  allocate(p_matrix(basis%nbf, basis%nbf, nspin))
  allocate(r_matrix(basis%nbf, basis%nbf))
  call setup_density_matrix(c_matrix, occupation, p_matrix)
  call setup_energy_density_matrix(c_matrix, occupation, energy, r_matrix)

  !
  ! Nucleus-nucleus repulsion forces
  !
  call nucleus_nucleus_force()


  !
  ! Energy density matrix forces
  !
  allocate(grad_onebody(basis%nbf, basis%nbf, 3))
  call setup_overlap_grad(basis, grad_onebody)

  ! F_A = + 2 \sum_{i \in A, j} W_ij d S_ij / d R_i
  force_ovp(:, :) = 0.0_dp
  do ibf=1, basis%nbf
    iatom = basis%bff(ibf)%icenter
    force_ovp(:, iatom) = force_ovp(:, iatom) + 2.0_dp * MATMUL( r_matrix(ibf, :) , grad_onebody(ibf, :, :) )
  end do
  deallocate(grad_onebody)


  !
  ! Kinetic energy forces
  !
  allocate(grad_onebody(basis%nbf, basis%nbf, 3))
  call setup_kinetic_grad(basis, grad_onebody)

  ! F_A = - 2 \sum_{i \in A, j} P_ij d T_ij / d R_i
  force_kin(:, :) = 0.0_dp
  do ibf=1, basis%nbf
    iatom = basis%bff(ibf)%icenter
    force_kin(:, iatom) = force_kin(:, iatom) - 2.0_dp * MATMUL( SUM( p_matrix(ibf, :, :), DIM=2) , grad_onebody(ibf, :, :) )
  end do
  deallocate(grad_onebody)


  !
  ! Nucleus energy forces
  !
  allocate(grad_nucleus(basis%nbf, basis%nbf, natom+1, 3))
  call setup_nucleus_grad(basis, grad_nucleus)

  ! Pulay: F_A = - 2 \sum_{i \in A, j} P_ij d V_ij / d R_i
  force_nuc(:, :) = 0.0_dp
  do ibf=1, basis%nbf
    iatom = basis%bff(ibf)%icenter
    force_nuc(:, iatom) = force_nuc(:, iatom) &
        - 2.0_dp * MATMUL( SUM( p_matrix(ibf, :, :), DIM=2) , grad_nucleus(ibf, :, natom+1, :) )
  end do

  ! Hellmann-Feynman: F_C = - \sum_ij P_ij ( i | d/dR_C -Z_C/|r-R_C| | j )
  force_hellfeyn(:, :) = 0.0_dp
  do iatom=1, natom
    force_hellfeyn(1, iatom) = force_hellfeyn(1, iatom) - SUM( SUM( p_matrix(:, :, :), DIM=3 ) * grad_nucleus(:, :, iatom, 1) )
    force_hellfeyn(2, iatom) = force_hellfeyn(2, iatom) - SUM( SUM( p_matrix(:, :, :), DIM=3 ) * grad_nucleus(:, :, iatom, 2) )
    force_hellfeyn(3, iatom) = force_hellfeyn(3, iatom) - SUM( SUM( p_matrix(:, :, :), DIM=3 ) * grad_nucleus(:, :, iatom, 3) )
  end do
  force_nuc(:, :) = force_nuc(:, :) + force_hellfeyn(:, :)


  deallocate(grad_nucleus)


  !
  ! Two-body (Hartree and exchange) forces
  !
  if( has_auxil_basis ) then
    call force_twobody_ri(basis, auxil_basis, occupation, c_matrix, p_matrix)
  else
    call force_twobody_4center(basis, p_matrix)
  end if

  !
  ! Exchange-correlation forces (including grid response)
  !
  force_exc(:, :) = 0.0_dp
  if( calc_type%is_dft ) then
    if( calc_type%need_exchange_lr ) call die('calculate_force: range-separated hybrids are not implemented')
    call init_dft_grid(basis, grid_level, dft_xc(1)%needs_gradient, .FALSE., BATCH_SIZE)
    allocate(vxc_ao(basis%nbf, basis%nbf, nspin))
    call dft_exc_vxc_batch(BATCH_SIZE, basis, occupation, c_matrix, vxc_ao, exc_xc, xc_force=force_exc)
    deallocate(vxc_ao)
    call destroy_dft_grid()
  end if

  ! is_core is an inefficient way to get the Kinetic+Nucleus hamiltonian
  if( calc_type%is_core ) force_har(:, :) = 0.0_dp



  !
  ! Total forces
  !
  force(:, :) = force_nuc_nuc(:, :) + force_ovp(:, :) &
               + force_kin(:, :) + force_nuc(:, :) + force_har(:, :) + force_exx(:, :) * alpha_hybrid + force_exc(:, :)


  write(stdout, '(/,1x,a)') ' ====== Hellman-Feynman Forces ====== '
  write(*, '(1x,a,22x,a,19x,a,19x,a)') 'Atoms', 'Fx', 'Fy','Fz'
  do iatom=1, natom
    write(*, '(3x,a,i4,a,2x,3(2x,f19.10))') 'H-F force atom   ', iatom, ':', force_hellfeyn(:, iatom)
  end do
  write(stdout, '(/,1x,a)') ' ====== Pulay Forces ====== '
  write(*, '(1x,a,22x,a,19x,a,19x,a)') 'Atoms', 'Fx', 'Fy','Fz'
  do iatom=1, natom
    write(*, '(3x,a,i4,a,2x,3(2x,f19.10))') 'Pulay force atom ', iatom, ':', force(:, iatom) - force_hellfeyn(:, iatom)
  end do
  write(stdout, '(/,1x,a)') ' ====== Total Forces ====== '
  write(*, '(1x,a,22x,a,19x,a,19x,a)') 'Atoms', 'Fx', 'Fy','Fz'
  do iatom=1, natom
    write(*, '(3x,a,i4,a,2x,3(2x,f19.10))') 'Total force atom ', iatom, ':', force(:, iatom)
  end do
  write(stdout, '(1x,a,/)') ' ==================== '


  deallocate(p_matrix)
  deallocate(r_matrix)

  call timer_force%stop()


end subroutine calculate_force


!=========================================================================
! Write the last calculated total force in the yaml file
! (before the removal of the net force and torque done for the relaxation)
subroutine print_force_yaml()
  !=====
  integer  :: iatom
  real(dp) :: force_total(3, ncenter_nuclei)
  !=====

  if( .NOT. ( print_yaml_ .AND. is_iomaster ) ) return

  force_total(:, :) = force_nuc_nuc(:, :) + force_ovp(:, :) + force_kin(:, :) + force_nuc(:, :) + force_har(:, :) &
                      + alpha_hybrid * force_exx(:, :) + force_exc(:, :)

  write(unit_yaml, '(/,a)') 'forces:'
  write(unit_yaml, '(4x,a)') 'unit: Ha/bohr'
  write(unit_yaml, '(4x,a)') 'total:'
  do iatom=1, ncenter_nuclei
    write(unit_yaml, '(8x,a,"[ ",es18.10,", ",es18.10,", ",es18.10,"]")') '- ', force_total(:, iatom)
  end do

end subroutine print_force_yaml


!=========================================================================
! Hartree and exchange forces with 4-center integral gradients
!
subroutine force_twobody_4center(basis, p_matrix)

  type(basis_set), intent(inout) :: basis
  real(dp), intent(in)           :: p_matrix(:, :, :)
  !=====
  real(dp), parameter      :: TOL_DENSITY_MATRIX=1.0e-6
  integer                 :: ijshellpair, klshellpair
  integer                 :: iatom
  integer                 :: ibf, jbf, kbf, lbf
  integer                 :: ishell, jshell, kshell, lshell
  integer                 :: ni, nj, nk, nl
  real(dp)                :: fact
  real(dp), allocatable    :: shell_gradA(:, :, :, :, :)
  real(dp), allocatable    :: shell_gradB(:, :, :, :, :)
  real(dp), allocatable    :: shell_gradC(:, :, :, :, :)
  real(dp), allocatable    :: shell_gradD(:, :, :, :, :)
  logical, allocatable     :: skip_shellpair(:, :)
  !=====

  !
  ! Filter out the low density matrix shells
  allocate(skip_shellpair(basis%nshell, basis%nshell))
  skip_shellpair(:, :) = .TRUE.
  do ijshellpair=1, nshellpair
    ishell = index_shellpair(1, ijshellpair)
    jshell = index_shellpair(2, ijshellpair)
    do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
      do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
        if( ANY( ABS(p_matrix(ibf, jbf, :) / spin_fact) > TOL_DENSITY_MATRIX ) ) then
          skip_shellpair(ishell, jshell) = .FALSE.
          skip_shellpair(jshell, ishell) = .FALSE.
        end if
      end do
    end do
  end do
  write(stdout, '(1x,a,i6,a,i6)') 'Shell pair skipped due to low density matrix screening:', &
                                 COUNT( skip_shellpair(:, :) ), ' / ', basis%nshell**2



  !
  ! Hartree-Fock energy forces
  !
  write(stdout, '(1x,a)') 'Calculate the Hartree-Fock part with 4-center integrals gradient'


  force_har(:, :) = 0.0_dp
  force_exx(:, :) = 0.0_dp
  do klshellpair=1, nshellpair
    kshell = index_shellpair(1, klshellpair)
    lshell = index_shellpair(2, klshellpair)
    nk = number_basis_function_am( basis%gaussian_type , basis%shell(kshell)%am )
    nl = number_basis_function_am( basis%gaussian_type , basis%shell(lshell)%am )

    do ijshellpair=1, nshellpair
      ishell = index_shellpair(1, ijshellpair)
      jshell = index_shellpair(2, ijshellpair)
      ni = number_basis_function_am( basis%gaussian_type , basis%shell(ishell)%am )
      nj = number_basis_function_am( basis%gaussian_type , basis%shell(jshell)%am )

      if( skip_shellpair(ishell, jshell) .AND. skip_shellpair(kshell, lshell) &
              .AND.  skip_shellpair(ishell, kshell) .AND. skip_shellpair(jshell, lshell)  &
              .AND.  skip_shellpair(ishell, lshell) .AND. skip_shellpair(jshell, kshell) ) then
        cycle
      end if

      ! Libint ordering is strict!
      if( basis%shell(ishell)%am + basis%shell(jshell)%am > basis%shell(kshell)%am + basis%shell(lshell)%am ) cycle
      allocate(shell_gradA(ni, nj, nk, nl, 3))
      allocate(shell_gradB(ni, nj, nk, nl, 3))
      allocate(shell_gradC(ni, nj, nk, nl, 3))
      allocate(shell_gradD(ni, nj, nk, nl, 3))
      call calculate_eri_4center_shell_grad(basis, 0.0_dp, ijshellpair, klshellpair, &
                                            shell_gradA, shell_gradB, shell_gradC, shell_gradD)

      !
      ! Hartree
      !
      if( .NOT. ( skip_shellpair(ishell, jshell) .AND. skip_shellpair(kshell, lshell) ) ) then
        if( kshell /= lshell ) then
          fact = -4.0_dp
        else
          fact = -2.0_dp
        end if

        iatom = basis%shell(ishell)%icenter
        do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
          do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
            do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
              do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                force_har(:, iatom) = force_har(:, iatom)  &
                                      + fact * SUM( p_matrix(ibf, jbf, :) * p_matrix(kbf, lbf, :) ) &
                                         * shell_gradA(ibf-basis%shell(ishell)%istart+1,       &
                                                       jbf-basis%shell(jshell)%istart+1,       &
                                                       kbf-basis%shell(kshell)%istart+1,       &
                                                       lbf-basis%shell(lshell)%istart+1, :)
              end do
            end do
          end do
        end do

        if( ishell /= jshell ) then
          iatom = basis%shell(jshell)%icenter
          do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
            do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
              do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                  force_har(:, iatom) = force_har(:, iatom)  &
                                        + fact * SUM( p_matrix(ibf, jbf, :) * p_matrix(kbf, lbf, :) ) &
                                           * shell_gradB(ibf-basis%shell(ishell)%istart+1,       &
                                                         jbf-basis%shell(jshell)%istart+1,       &
                                                         kbf-basis%shell(kshell)%istart+1,       &
                                                         lbf-basis%shell(lshell)%istart+1, :)
                end do
              end do
            end do
          end do
        end if

        !
        ! When the opposite is not calculated by LIBINT:
        if( basis%shell(ishell)%am + basis%shell(jshell)%am < basis%shell(kshell)%am + basis%shell(lshell)%am ) then
          if( ishell /= jshell ) then
            fact = -4.0_dp
          else
            fact = -2.0_dp
          end if

          iatom = basis%shell(kshell)%icenter
          do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
            do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
              do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                  force_har(:, iatom) = force_har(:, iatom)  &
                                        + fact * SUM( p_matrix(kbf, lbf, :) * p_matrix(ibf, jbf, :) ) &
                                           * shell_gradC(ibf-basis%shell(ishell)%istart+1,       &
                                                         jbf-basis%shell(jshell)%istart+1,       &
                                                         kbf-basis%shell(kshell)%istart+1,       &
                                                         lbf-basis%shell(lshell)%istart+1, :)
                end do
              end do
            end do
          end do

          if( kshell /= lshell ) then
            iatom = basis%shell(lshell)%icenter
            do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
              do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
                do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                  do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                    force_har(:, iatom) = force_har(:, iatom)  &
                                          + fact * SUM( p_matrix(lbf, kbf, :) * p_matrix(ibf, jbf, :) ) &
                                             * shell_gradD(ibf-basis%shell(ishell)%istart+1,       &
                                                           jbf-basis%shell(jshell)%istart+1,       &
                                                           kbf-basis%shell(kshell)%istart+1,       &
                                                           lbf-basis%shell(lshell)%istart+1, :)
                  end do
                end do
              end do
            end do
          end if


        end if
      end if

      !
      ! Exchange
      !
      if( .NOT. ( skip_shellpair(ishell, kshell) .AND. skip_shellpair(jshell, lshell) ) &
          .OR. .NOT. ( skip_shellpair(ishell, lshell) .AND. skip_shellpair(jshell, kshell) ) ) then
        iatom = basis%shell(ishell)%icenter
        do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
          do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
            do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
              do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                force_exx(:, iatom) = force_exx(:, iatom)  &
                                      + 2.0_dp / spin_fact * SUM( p_matrix(ibf, kbf, :) * p_matrix(jbf, lbf, :) ) &
                                         * shell_gradA(ibf-basis%shell(ishell)%istart+1,       &
                                                       jbf-basis%shell(jshell)%istart+1,       &
                                                       kbf-basis%shell(kshell)%istart+1,       &
                                                       lbf-basis%shell(lshell)%istart+1, :)
              end do
            end do
          end do
        end do

        if( kshell /= lshell ) then
          iatom = basis%shell(ishell)%icenter
          do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
            do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
              do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                  force_exx(:, iatom) = force_exx(:, iatom)  &
                                        + 2.0_dp / spin_fact * SUM( p_matrix(ibf, lbf, :) * p_matrix(jbf, kbf, :) ) &
                                           * shell_gradA(ibf-basis%shell(ishell)%istart+1,       &
                                                         jbf-basis%shell(jshell)%istart+1,       &
                                                         kbf-basis%shell(kshell)%istart+1,       &
                                                         lbf-basis%shell(lshell)%istart+1, :)
                end do
              end do
            end do
          end do
        end if

        if( ishell /= jshell ) then
          iatom = basis%shell(jshell)%icenter
          do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
            do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
              do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                  force_exx(:, iatom) = force_exx(:, iatom)  &
                                        + 2.0_dp / spin_fact * SUM( p_matrix(jbf, kbf, :) * p_matrix(ibf, lbf, :) ) &
                                           * shell_gradB(ibf-basis%shell(ishell)%istart+1,       &
                                                         jbf-basis%shell(jshell)%istart+1,       &
                                                         kbf-basis%shell(kshell)%istart+1,       &
                                                         lbf-basis%shell(lshell)%istart+1, :)
                end do
              end do
            end do
          end do

          if( kshell /= lshell ) then
            iatom = basis%shell(jshell)%icenter
            do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
              do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
                do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                  do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                    force_exx(:, iatom) = force_exx(:, iatom)  &
                                          + 2.0_dp / spin_fact * SUM( p_matrix(jbf, lbf, :) * p_matrix(ibf, kbf, :) ) &
                                             * shell_gradB(ibf-basis%shell(ishell)%istart+1,       &
                                                           jbf-basis%shell(jshell)%istart+1,       &
                                                           kbf-basis%shell(kshell)%istart+1,       &
                                                           lbf-basis%shell(lshell)%istart+1, :)
                  end do
                end do
              end do
            end do
          end if
        end if

        !
        ! When the opposite is not calculated by LIBINT:
        if( basis%shell(ishell)%am + basis%shell(jshell)%am /= basis%shell(kshell)%am + basis%shell(lshell)%am ) then
          iatom = basis%shell(kshell)%icenter
          do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
            do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
              do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                  force_exx(:, iatom) = force_exx(:, iatom)  &
                                        + 2.0_dp / spin_fact * SUM( p_matrix(jbf, lbf, :) * p_matrix(ibf, kbf, :) ) &
                                           * shell_gradC(ibf-basis%shell(ishell)%istart+1,       &
                                                         jbf-basis%shell(jshell)%istart+1,       &
                                                         kbf-basis%shell(kshell)%istart+1,       &
                                                         lbf-basis%shell(lshell)%istart+1, :)
                end do
              end do
            end do
          end do

          if( ishell /= jshell ) then
            iatom = basis%shell(kshell)%icenter
            do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
              do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
                do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                  do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                    force_exx(:, iatom) = force_exx(:, iatom)  &
                                          + 2.0_dp / spin_fact * SUM( p_matrix(ibf, lbf, :) * p_matrix(jbf, kbf, :) ) &
                                             * shell_gradC(ibf-basis%shell(ishell)%istart+1,       &
                                                           jbf-basis%shell(jshell)%istart+1,       &
                                                           kbf-basis%shell(kshell)%istart+1,       &
                                                           lbf-basis%shell(lshell)%istart+1, :)
                  end do
                end do
              end do
            end do
          end if

          if( kshell /= lshell ) then
            iatom = basis%shell(lshell)%icenter
            do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
              do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
                do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                  do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                    force_exx(:, iatom) = force_exx(:, iatom)  &
                                          + 2.0_dp / spin_fact * SUM( p_matrix(jbf, kbf, :) * p_matrix(ibf, lbf, :) ) &
                                             * shell_gradD(ibf-basis%shell(ishell)%istart+1,       &
                                                           jbf-basis%shell(jshell)%istart+1,       &
                                                           kbf-basis%shell(kshell)%istart+1,       &
                                                           lbf-basis%shell(lshell)%istart+1, :)
                  end do
                end do
              end do
            end do
            if( ishell /= jshell ) then
              iatom = basis%shell(lshell)%icenter
              do lbf=basis%shell(lshell)%istart, basis%shell(lshell)%iend
                do kbf=basis%shell(kshell)%istart, basis%shell(kshell)%iend
                  do jbf=basis%shell(jshell)%istart, basis%shell(jshell)%iend
                    do ibf=basis%shell(ishell)%istart, basis%shell(ishell)%iend
                      force_exx(:, iatom) = force_exx(:, iatom)  &
                                            + 2.0_dp / spin_fact * SUM( p_matrix(ibf, kbf, :) * p_matrix(jbf, lbf, :) ) &
                                               * shell_gradD(ibf-basis%shell(ishell)%istart+1,       &
                                                             jbf-basis%shell(jshell)%istart+1,       &
                                                             kbf-basis%shell(kshell)%istart+1,       &
                                                             lbf-basis%shell(lshell)%istart+1, :)
                    end do
                  end do
                end do
              end do

            end if

          end if

        end if
      end if


      deallocate(shell_gradA, shell_gradB, shell_gradC, shell_gradD)
    end do
  end do

  deallocate(skip_shellpair)

end subroutine force_twobody_4center


!=========================================================================
! Hartree and exchange forces with RI (Coulomb metric V_PQ = (P|Q))
!  E_J = 1/2 t^T V^{-1} t                            t_P   = \sum_ik P_ik (ik|P)
!  E_K = -sf/2 \sum_s \sum_ab y_ab^T V^{-1} y_ab     y_abP = \sum_ik ct_ai ct_bk (ik|P),  ct_ai = sqrt(f_a/sf) C_ia
!  dE  = \sum_ikP G_ikP d(ik|P)/dR + \sum_PQ W_PQ d(P|Q)/dR
!  G_ikP = P_ik d_P - sf \sum_s (ct^T Z_P ct)_ik     d = V^{-1} t,  Z_ab = V^{-1} y_ab
!  W_PQ  = -1/2 d_P d_Q + sf/2 \sum_s \sum_ab Z_abP Z_abQ
!
subroutine force_twobody_ri(basis, auxil_basis, occupation, c_matrix, p_matrix)

  type(basis_set), intent(in) :: basis
  type(basis_set), intent(in) :: auxil_basis
  real(dp), intent(in)        :: occupation(:, :)
  real(dp), intent(in)        :: c_matrix(:, :, :)
  real(dp), intent(in)        :: p_matrix(:, :, :)
  !=====
  integer               :: nbf, nauxil, nocc, nkept
  integer               :: ispin, istate, ipair, iswap, idir, iauxil
  integer               :: ishell, kshell, pshell, qshell
  integer               :: ibf1, ibf2, kbf1, kbf2, pbf1, pbf2, qbf1, qbf2, pbf
  integer               :: iatom, patom
  real(dp)              :: val_j(3), val_k(3)
  real(dp), allocatable :: p_total(:, :)
  real(dp), allocatable :: ct(:, :, :)
  real(dp), allocatable :: t_vector(:), d_vector(:)
  real(dp), allocatable :: y_matrix(:, :, :), z_matrix(:, :, :)
  real(dp), allocatable :: v_matrix(:, :), v_inv(:, :), eigvec(:, :), eigval(:)
  real(dp), allocatable :: w_j(:, :), w_k(:, :)
  real(dp), allocatable :: eri3(:, :, :), g_j(:, :, :), g_k(:, :, :)
  real(dp), allocatable :: integrals(:, :, :), integrals2(:, :)
  real(dp), allocatable :: grad3(:, :, :, :), grad2(:, :, :)
  real(dp), allocatable :: grad_j(:, :), grad_k(:, :)
  !=====

  write(stdout, '(1x,a)') 'Calculate the Hartree and exchange forces with RI'

  nbf    = basis%nbf
  nauxil = auxil_basis%nbf
  nocc   = get_number_occupied_states(occupation)

  allocate(p_total(nbf, nbf))
  p_total(:, :) = SUM(p_matrix(:, :, :), DIM=3)

  ! ct_ai = sqrt(f_a / sf) C_ia
  allocate(ct(nocc, nbf, nspin))
  do ispin=1, nspin
    do istate=1, nocc
      ct(istate, :, ispin) = c_matrix(:, istate, ispin) * SQRT( MAX(occupation(istate, ispin), 0.0_dp) / spin_fact )
    end do
  end do

  !
  ! t_P and y_abP from the 3-center integrals
  allocate(t_vector(nauxil), y_matrix(nocc*nocc, nauxil, nspin))
  t_vector(:) = 0.0_dp
  y_matrix(:, :, :) = 0.0_dp
  do pshell=1, auxil_basis%nshell
    if( MODULO(pshell - 1, world%nproc) /= world%rank ) cycle
    pbf1 = auxil_basis%shell(pshell)%istart
    pbf2 = auxil_basis%shell(pshell)%iend
    allocate(eri3(nbf, nbf, pbf2 - pbf1 + 1))
    eri3(:, :, :) = 0.0_dp
    do ipair=1, nshellpair
      ishell = index_shellpair(1, ipair)
      kshell = index_shellpair(2, ipair)
      ibf1 = basis%shell(ishell)%istart
      ibf2 = basis%shell(ishell)%iend
      kbf1 = basis%shell(kshell)%istart
      kbf2 = basis%shell(kshell)%iend
      call libcint_3center_shell(basis, auxil_basis, ishell, kshell, pshell, integrals)
      eri3(ibf1:ibf2, kbf1:kbf2, :) = integrals(:, :, :)
      eri3(kbf1:kbf2, ibf1:ibf2, :) = RESHAPE(integrals, [kbf2-kbf1+1, ibf2-ibf1+1, pbf2-pbf1+1], ORDER=[2, 1, 3])
    end do
    do pbf=pbf1, pbf2
      t_vector(pbf) = SUM( p_total(:, :) * eri3(:, :, pbf-pbf1+1) )
      do ispin=1, nspin
        y_matrix(:, pbf, ispin) = RESHAPE( MATMUL( ct(:, :, ispin), &
                                           MATMUL( eri3(:, :, pbf-pbf1+1), TRANSPOSE(ct(:, :, ispin)) ) ), [nocc*nocc] )
      end do
    end do
    deallocate(eri3)
  end do
  call world%sum(t_vector)
  call world%sum(y_matrix)

  !
  ! V^{-1} with the same eigenvalue filtering as in the SCF
  allocate(v_matrix(nauxil, nauxil))
  do pshell=1, auxil_basis%nshell
    pbf1 = auxil_basis%shell(pshell)%istart
    pbf2 = auxil_basis%shell(pshell)%iend
    do qshell=1, auxil_basis%nshell
      qbf1 = auxil_basis%shell(qshell)%istart
      qbf2 = auxil_basis%shell(qshell)%iend
      call libcint_2center_shell(auxil_basis, pshell, qshell, integrals2)
      v_matrix(pbf1:pbf2, qbf1:qbf2) = integrals2(:, :)
    end do
  end do
  allocate(eigval(nauxil), eigvec(nauxil, nauxil), v_inv(nauxil, nauxil))
  call diagonalize(' ', v_matrix, eigval, eigvec)
  nkept = 0
  v_inv(:, :) = 0.0_dp
  do iauxil=1, nauxil
    if( eigval(iauxil) <= TOO_LOW_EIGENVAL ) cycle
    nkept = nkept + 1
    call DSYR('L', nauxil, 1.0_dp / eigval(iauxil), eigvec(:, iauxil), 1, v_inv, nauxil)
  end do
  do iauxil=1, nauxil
    v_inv(iauxil, iauxil+1:) = v_inv(iauxil+1:, iauxil)
  end do
  if( nkept < nauxil ) then
    call issue_warning('force_twobody_ri: auxiliary eigenvalues were filtered out, RI forces are not exact')
  end if
  deallocate(v_matrix, eigval, eigvec)

  ! d = V^{-1} t,  Z_ab = V^{-1} y_ab
  allocate(d_vector(nauxil), z_matrix(nocc*nocc, nauxil, nspin))
  d_vector(:) = MATMUL(v_inv, t_vector)
  do ispin=1, nspin
    z_matrix(:, :, ispin) = MATMUL(y_matrix(:, :, ispin), v_inv)
  end do
  deallocate(v_inv, y_matrix, t_vector)

  ! W_J = -1/2 d d^T,  W_K = sf/2 \sum_s Z^T Z
  allocate(w_j(nauxil, nauxil), w_k(nauxil, nauxil))
  do iauxil=1, nauxil
    w_j(:, iauxil) = -0.5_dp * d_vector(:) * d_vector(iauxil)
  end do
  w_k(:, :) = 0.0_dp
  do ispin=1, nspin
    w_k(:, :) = w_k(:, :) + 0.5_dp * spin_fact * MATMUL( TRANSPOSE(z_matrix(:, :, ispin)), z_matrix(:, :, ispin) )
  end do

  allocate(grad_j(3, natom), grad_k(3, natom))
  grad_j(:, :) = 0.0_dp
  grad_k(:, :) = 0.0_dp

  !
  ! 3-center derivatives: dE/dR_A = 2 \sum_{i \in A} \sum_kP G_ikP d(ik|P)/dR_i - 2 \sum_{P \in A} \sum_ik G_ikP d(ik|P)/dR_i
  do pshell=1, auxil_basis%nshell
    if( MODULO(pshell - 1, world%nproc) /= world%rank ) cycle
    pbf1  = auxil_basis%shell(pshell)%istart
    pbf2  = auxil_basis%shell(pshell)%iend
    patom = auxil_basis%shell(pshell)%icenter
    allocate(g_j(nbf, nbf, pbf2 - pbf1 + 1), g_k(nbf, nbf, pbf2 - pbf1 + 1))
    do pbf=pbf1, pbf2
      g_j(:, :, pbf-pbf1+1) = p_total(:, :) * d_vector(pbf)
      g_k(:, :, pbf-pbf1+1) = 0.0_dp
      do ispin=1, nspin
        g_k(:, :, pbf-pbf1+1) = g_k(:, :, pbf-pbf1+1) - spin_fact &
                   * MATMUL( TRANSPOSE(ct(:, :, ispin)), MATMUL( RESHAPE(z_matrix(:, pbf, ispin), [nocc, nocc]), ct(:, :, ispin) ) )
      end do
    end do

    do ipair=1, nshellpair
      do iswap=1, 2
        if( iswap == 1 ) then
          ishell = index_shellpair(1, ipair)
          kshell = index_shellpair(2, ipair)
        else
          if( index_shellpair(1, ipair) == index_shellpair(2, ipair) ) cycle
          ishell = index_shellpair(2, ipair)
          kshell = index_shellpair(1, ipair)
        end if
        ibf1  = basis%shell(ishell)%istart
        ibf2  = basis%shell(ishell)%iend
        kbf1  = basis%shell(kshell)%istart
        kbf2  = basis%shell(kshell)%iend
        iatom = basis%shell(ishell)%icenter
        ! grad3 = ( \nabla_r i k | P ) = - d(ik|P)/dR_i
        call libcint_3center_ip1_shell(basis, auxil_basis, ishell, kshell, pshell, grad3)
        do idir=1, 3
          val_j(idir) = -SUM( g_j(ibf1:ibf2, kbf1:kbf2, :) * grad3(:, :, :, idir) )
          val_k(idir) = -SUM( g_k(ibf1:ibf2, kbf1:kbf2, :) * grad3(:, :, :, idir) )
        end do
        grad_j(:, iatom) = grad_j(:, iatom) + 2.0_dp * val_j(:)
        grad_j(:, patom) = grad_j(:, patom) - 2.0_dp * val_j(:)
        grad_k(:, iatom) = grad_k(:, iatom) + 2.0_dp * val_k(:)
        grad_k(:, patom) = grad_k(:, patom) - 2.0_dp * val_k(:)
      end do
    end do
    deallocate(g_j, g_k)
  end do

  !
  ! 2-center derivatives: dE/dR_A = 2 \sum_{P \in A} \sum_Q W_PQ d(P|Q)/dR_P
  do pshell=1, auxil_basis%nshell
    if( MODULO(pshell - 1, world%nproc) /= world%rank ) cycle
    pbf1  = auxil_basis%shell(pshell)%istart
    pbf2  = auxil_basis%shell(pshell)%iend
    patom = auxil_basis%shell(pshell)%icenter
    do qshell=1, auxil_basis%nshell
      qbf1 = auxil_basis%shell(qshell)%istart
      qbf2 = auxil_basis%shell(qshell)%iend
      ! grad2 = ( \nabla_r P | Q ) = - d(P|Q)/dR_P
      call libcint_2center_ip1_shell(auxil_basis, pshell, qshell, grad2)
      do idir=1, 3
        grad_j(idir, patom) = grad_j(idir, patom) - 2.0_dp * SUM( w_j(pbf1:pbf2, qbf1:qbf2) * grad2(:, :, idir) )
        grad_k(idir, patom) = grad_k(idir, patom) - 2.0_dp * SUM( w_k(pbf1:pbf2, qbf1:qbf2) * grad2(:, :, idir) )
      end do
    end do
  end do

  call world%sum(grad_j)
  call world%sum(grad_k)

  ! F = - dE/dR
  force_har(:, :) = -grad_j(:, :)
  force_exx(:, :) = -grad_k(:, :)

  deallocate(grad_j, grad_k, w_j, w_k, d_vector, z_matrix, ct, p_total)

end subroutine force_twobody_ri


end module m_force
!=========================================================================
