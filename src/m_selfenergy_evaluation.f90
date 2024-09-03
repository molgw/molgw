!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the driver for the different self-energy methods:
! PT2, PT3, GW, evGW, COHSEX, GWGamma, etc.
!
!=========================================================================
#include "molgw.h"
module m_selfenergy_evaluation
  use m_definitions
  use m_timing
  use m_warning
  use m_memory
  use m_inputparam
  use m_eri
  use m_eri_calculate
  use m_eri_ao_mo
  use m_scf,only: energy_contributions
  use m_spectral_function
  use m_selfenergy_tools
  use m_virtual_orbital_space
  use m_io
  use m_gw_selfenergy_grid
  use m_linear_response
  use m_g3w2_selfenergy
  use m_tdhf_selfenergy



contains


!=========================================================================
subroutine selfenergy_evaluation(basis,occupation,energy,c_matrix,exchange_m_vxc,en_mbpt)
  implicit none

  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: occupation(:,:)
  real(dp),intent(inout)     :: energy(:,:)
  real(dp),intent(inout)     :: c_matrix(:,:,:)
  real(dp),intent(in)        :: exchange_m_vxc(:,:,:)
  type(energy_contributions),intent(inout) :: en_mbpt
  !=====
  integer                 :: nstate
  type(selfenergy_grid)   :: se,se2,se3,se_sox,se_gwpt3,se_gw,se_sosex,se_g3w2
  logical                 :: enforce_rpa
  character(len=36)       :: selfenergy_tag
  integer                 :: reading_status
  integer                 :: nstate_small
  type(spectral_function) :: wpol
  real(dp),allocatable    :: zz(:,:)
  real(dp),allocatable    :: energy_qp_new(:,:),energy_qp_z(:,:)
  integer                 :: iomega
  integer                 :: istep_gw,pstate
  real(dp),allocatable    :: exchange_m_vxc_diag(:,:)
  real(dp),allocatable    :: energy_g(:,:)
  real(dp),allocatable    :: energy_w(:,:)
  !=====

  write(stdout,'(/,/,1x,a)') '=================================================='
  write(stdout,'(1x,a)')     'Self-energy evaluation starts here'
  write(stdout,'(1x,a,/)')   '=================================================='

  nstate = SIZE(occupation(:,:),DIM=1)
  allocate(energy_g(nstate,nspin),energy_w(nstate,nspin))
  allocate(exchange_m_vxc_diag(nstate,nspin))
  do pstate=1,nstate
    exchange_m_vxc_diag(pstate,:) = exchange_m_vxc(pstate,pstate,:)
  enddo


  do istep_gw=1,nstep_gw

    !
    ! Set the character string for the calculation we are currently doing
    select case(calc_type%selfenergy_approx)
    case(GW)
      selfenergy_tag='GW'
    case(GnW0)
      write(selfenergy_tag,'(i3)') istep_gw-1
      selfenergy_tag='G'//TRIM(ADJUSTL(selfenergy_tag))//'W0'
    case(GnWn)
      write(selfenergy_tag,'(i3)') istep_gw-1
      selfenergy_tag='G'//TRIM(ADJUSTL(selfenergy_tag))//'W'//TRIM(ADJUSTL(selfenergy_tag))
    case(PT2)
      selfenergy_tag='PT2'
    case(TWO_RINGS)
      selfenergy_tag='TWO_RINGS'
    case(PT3)
      selfenergy_tag='PT3'
    case(ONE_RING)
      selfenergy_tag='ONE_RING'
    case(SOX)
      selfenergy_tag='SOX'
    case(GWSOX)
      selfenergy_tag='GW+SOX'
    case(GWPT3)
      selfenergy_tag='GW+PT3'
    case(GWSOSEX)
      selfenergy_tag='GW+SOSEX'
    case(COHSEX)
      selfenergy_tag='COHSEX'
    case(G3W2,G3W2_NUMERICAL)
      selfenergy_tag='GW+G3W2'
    case(GW0GW0G)
      selfenergy_tag='GW+GW0GW0G'
    case(GWGW0G)
      selfenergy_tag='GW+GWGW0G'
    case(GWGW0RPAG)
      selfenergy_tag='GW+GWGW0RPAG'
    case(SIGMA_TDHF)
      selfenergy_tag='SIGMA_TDHF'
    case default
      write(stdout,*) 'selfenergy approx not listed:',calc_type%selfenergy_approx
      call die('selfenergy_evaluation: bug')
    end select


    !
    ! Set the range of states on which to evaluate the self-energy
    call selfenergy_set_state_range(nstate,occupation)

    !
    ! If requested,
    ! prepare an optmized virtual subspace based on
    ! Frozen Natural Orbitals technique
    if( is_virtual_fno ) then
      !
      ! Be aware that the energies and the c_matrix for virtual orbitals are altered after this point
      ! and until they are restored in destroy_fno
      !
      call calculate_virtual_fno(basis,nstate,nsemax,occupation,energy,c_matrix)
    endif
    !
    ! Or alternatively use the small basis technique
    if( has_small_basis ) then
      call setup_virtual_smallbasis(basis,nstate,occupation,nsemax,energy,c_matrix,nstate_small)
      !
      ! Set the range again after the change of the virtual space
      ! to nstate
      call selfenergy_set_state_range(nstate_small,occupation)
    else
      nstate_small = nstate
    endif


    !
    ! Choose which one-electron energies to use in G and in W
    !
    if( calc_type%selfenergy_technique == EVSC .OR. force_energy_qp_ ) then
      call read_energy_qp(nstate,energy_g,reading_status)
      if( reading_status /=0 ) then
        call issue_warning('File energy_qp not found: assuming 1st iteration')
        energy_g(:,:) = energy(:,:)
      endif

      !
      ! For GnWn, update both the energy in G and in W
      if( calc_type%selfenergy_approx == GnWn ) then
        energy_w(:,:) = energy_g(:,:)
      else
        energy_w(:,:) = energy(:,:)
      endif

    else
      energy_g(:,:) = energy(:,:)
      energy_w(:,:) = energy(:,:)
    endif


    call se%init(calc_type%selfenergy_technique,energy_g)


    !
    ! selfenergy = GW or COHSEX
    !
    if(    calc_type%selfenergy_approx == GW          &
      .OR. calc_type%selfenergy_approx == COHSEX      &
      .OR. calc_type%selfenergy_approx == GnW0        &
      .OR. calc_type%selfenergy_approx == GnWn   ) then

      !
      ! First calculate W except if performing GnW0 for the second and following times
      !
      if( calc_type%selfenergy_approx == GnW0 .AND. istep_gw > 1 ) then
        write(stdout,'(/,1x,a,/)') 'GnW0 calculations skip the re-calculation of W'
      else

        ! Try to read a spectral function file in order to skip the polarizability calculation
        ! Skip the reading if GnWn (=evGW) is requested
        if( calc_type%selfenergy_approx /= GnWn ) then
          call read_spectral_function(wpol,reading_status)
        else
          write(stdout,'(/,1x,a)') 'For GnWn calculations, never try to read file SCREENED_COULOMB'
          reading_status = 1
        endif
        ! If reading has failed, then do the calculation
        if( reading_status /= 0 ) then
          select case(calc_type%selfenergy_technique)
          case(imaginary_axis_pade,imaginary_axis_homolumo)
            call wpol%init(nstate_small,occupation,nomega_chi_imag,grid_type=IMAGINARY_QUAD)
            call polarizability_grid_scalapack(occupation,energy_w,c_matrix,en_mbpt%rpa,en_mbpt%gw,wpol)
          case(contour_deformation)
            ! no need for chi, it will be calculated directly inside
          case default
            ! in case of BSE calculation, enforce RPA here
            enforce_rpa = calc_type%is_bse
            call wpol%init(nstate_small,occupation,0)
            call polarizability(enforce_rpa,.TRUE.,basis,occupation,energy_w,c_matrix,en_mbpt%rpa,en_mbpt%gw,wpol)
          end select
        endif

        en_mbpt%total = en_mbpt%total + en_mbpt%rpa
        en_mbpt%total = en_mbpt%total - en_mbpt%xc - en_mbpt%exx_hyb + en_mbpt%exx

        if( ABS(en_mbpt%rpa) > 1.e-6_dp ) then
          write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en_mbpt%total
        endif

      endif

      select case(calc_type%selfenergy_technique)
      case(contour_deformation)
        call gw_selfenergy_contour(energy_g,occupation,c_matrix,se)
      case(imaginary_axis_pade)
        call gw_selfenergy_imag_scalapack(energy_g,c_matrix,wpol,se)
        call se%pade_fit()
      case(imaginary_axis_homolumo)
        call gw_selfenergy_imag_scalapack(energy_g,c_matrix,wpol,se)
        call self_energy_polynomial(se)
      case(exact_dyson)
        call gw_selfenergy_analytic(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,wpol,exchange_m_vxc)
      case default
        ! The SCALAPACK implementation only works for plain vanilla GW
#if defined(HAVE_SCALAPACK)
        if( has_auxil_basis &
           .AND. (calc_type%selfenergy_approx == GW .OR. calc_type%selfenergy_approx == GnW0  &
             .OR. calc_type%selfenergy_approx == GnWn) ) then
          call gw_selfenergy_scalapack(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,wpol,se)
        else
#endif
          call gw_selfenergy(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,wpol,se)
#if defined(HAVE_SCALAPACK)
        endif
#endif
      end select

      if( ABS(en_mbpt%gw) > 1.0e-5_dp ) then
        write(stdout,'(/,a,f19.10)') ' Galitskii-Migdal Total energy (Ha): ',en_mbpt%total - en_mbpt%rpa + en_mbpt%gw
      endif

      if( .NOT. ( calc_type%selfenergy_approx == GnW0 .AND. istep_gw < nstep_gw ) ) then
        call wpol%destroy()
      endif

      if( has_small_basis ) then
        !
        ! Output the G0W0 results in the small basis first
        allocate(energy_qp_z(nstate,nspin))
        allocate(energy_qp_new(nstate,nspin))
        allocate(zz(nstate,nspin))
        call find_qp_energy_linearization(se,exchange_m_vxc_diag,energy,energy_qp_z,zz)
        call find_qp_energy_graphical(se,exchange_m_vxc_diag,energy,energy_qp_new,zz)
        call output_qp_energy('GW small basis',energy,exchange_m_vxc_diag,1,se,energy_qp_z,energy_qp_new,zz)
        deallocate(zz)
        deallocate(energy_qp_z)
        call output_homolumo('GW small basis',occupation,energy_qp_new,nsemin,nsemax)
        deallocate(energy_qp_new)

        call se2%init(static_selfenergy,energy)
        call se3%init(static_selfenergy,energy)

        ! Sigma^2 = Sigma^{1-ring}_small
        call onering_selfenergy(nstate_small,basis,occupation(1:nstate_small,:), &
                                 energy_g(1:nstate_small,:),c_matrix(:,1:nstate_small,:),se2,en_mbpt%mp2)

        ! Reset wavefunctions, eigenvalues and number of virtual orbitals in G
        call destroy_fno(basis,nstate,energy,c_matrix)
        energy_g(:,:) = energy(:,:)
        call selfenergy_set_state_range(nstate,occupation)

        ! Sigma^3 = Sigma^{1-ring}_big
        call onering_selfenergy(nstate,basis,occupation,energy_g,c_matrix,se3,en_mbpt%mp2)

        if( print_sigma_ ) then
          call write_selfenergy_omega('selfenergy_GW_small'   ,exchange_m_vxc_diag,occupation,energy_g,se)
          call write_selfenergy_omega('selfenergy_1ring_big'  ,exchange_m_vxc_diag,occupation,energy_g,se3)
          call write_selfenergy_omega('selfenergy_1ring_small',exchange_m_vxc_diag,occupation,energy_g,se2)
        endif

        !
        ! Extrapolated Sigma(omega) = Sigma^{GW}_small(omega) + Sigma^{1-ring}_big(0) - Sigma^{1-ring}_small(0)
        do iomega=-se%nomega,se%nomega
          se%sigma(iomega,:,:) = se%sigma(iomega,:,:) + se3%sigma(0,:,:) - se2%sigma(0,:,:)
        enddo

        call se2%destroy()
        call se3%destroy()

      endif

    endif

    !
    !  sigma_TDHF self-energy (See Vacondio et al., Forster-Bruneval)
    !
    if( calc_type%selfenergy_approx == SIGMA_TDHF ) then
      call tdhf_selfenergy(basis,occupation,energy_g,c_matrix,se)
    endif

    !
    ! GW+SOX or 
    ! GW+SOSEX or GW+2SOSEX or GW+G3W2
    !
    if( calc_type%selfenergy_approx == GWSOX &
        .OR. calc_type%selfenergy_approx == GWSOSEX &
        .OR. calc_type%selfenergy_approx == G3W2 &
        .OR. calc_type%selfenergy_approx == G3W2_NUMERICAL &
        .OR. calc_type%selfenergy_approx == GW0GW0G &
        .OR. calc_type%selfenergy_approx == GWGW0G &
        .OR. calc_type%selfenergy_approx == GWGW0RPAG &
      ) then

      !
      ! Two implementations available:
      ! analytic or imaginary frequencies
      !
      if( calc_type%selfenergy_technique /= imaginary_axis_pade ) then
        call wpol%init(nstate,occupation,0)
        call read_spectral_function(wpol,reading_status)
        ! If reading has failed, then do the calculation
        if( reading_status /= 0 ) then
          call polarizability(.FALSE.,.TRUE.,basis,occupation,energy_w,c_matrix,en_mbpt%rpa,en_mbpt%gw,wpol)
        endif

        call se_gw%init(calc_type%selfenergy_technique,energy_g)
        call gw_selfenergy(GW,nstate,basis,occupation,energy_g,c_matrix,wpol,se_gw)

        !
        ! Output the G0W0 results first
        allocate(energy_qp_z(nstate,nspin))
        allocate(energy_qp_new(nstate,nspin))
        allocate(zz(nstate,nspin))
        call find_qp_energy_linearization(se_gw,exchange_m_vxc_diag,energy,energy_qp_z,zz)
        call find_qp_energy_graphical(se_gw,exchange_m_vxc_diag,energy,energy_qp_new,zz)
        call output_qp_energy('GW',energy,exchange_m_vxc_diag,1,se_gw,energy_qp_z,energy_qp_new,zz)
        call output_qp_energy_yaml('GW',energy,exchange_m_vxc_diag,se_gw,energy_qp_z,energy_qp_new,zz)
        call output_homolumo('GW',occupation,energy_qp_new,nsemin,nsemax)
        call dump_out_energy_yaml('gw energies',energy_qp_new,nsemin,nsemax)
        deallocate(zz)
        deallocate(energy_qp_z)

        if( g3w2_static_approximation_ ) then
          call issue_warning("selfenergy_evaluation: use Arno's approximation for G3W2")
          ! enforce a single frequency located at the GW qp energy
          call se_g3w2%init(static_selfenergy,energy_qp_new)
        else
          call se_g3w2%init(calc_type%selfenergy_technique,energy_g)
        endif

        !
        ! selfenergy = GWSOX
        !
        select case(calc_type%selfenergy_approx)
        case(GWSOX)
          !
          ! Perform a standard SOX calculation
          !
          call se_sox%init(calc_type%selfenergy_technique,energy_g)
          call pt2_selfenergy(SOX,nstate,basis,occupation,energy_g,c_matrix,se_sox,en_mbpt%mp2)

          !
          ! Finally add up the contributions and then destroy the se_sox object
          !
          call se%add(se_gw)
          call se%add(se_sox)
          call se_sox%destroy()

        case(GW0GW0G,GWGW0G,GWGW0RPAG)
          call gwgw0g_selfenergy(occupation,energy_g,c_matrix,wpol,se_g3w2)

          call se%add(se_gw)
          call se%add(se_g3w2)
          call se_g3w2%destroy()

        case(GWSOSEX)
          call sosex_selfenergy(basis,occupation,energy_g,c_matrix,wpol,se_g3w2)

          call se%add(se_gw)
          call se%add(se_g3w2)

        case(G3W2)
          call sosex_selfenergy(basis,occupation,energy_g,c_matrix,wpol,se_g3w2)
          !call se%reset()
          call se%add(se_gw)
          call se%add(se_g3w2)

          ! Output the GW+SOSEX qp energies before over-riding them
          allocate(energy_qp_z(nstate,nspin))
          allocate(zz(nstate,nspin))
          call find_qp_energy_linearization(se,exchange_m_vxc_diag,energy,energy_qp_z,zz)
          call find_qp_energy_graphical(se,exchange_m_vxc_diag,energy,energy_qp_new,zz)
          call output_qp_energy('GW+SOSEX2',energy,exchange_m_vxc_diag,1,se,energy_qp_z,energy_qp_new,zz)
          call output_qp_energy_yaml('GW+SOSEX2',energy,exchange_m_vxc_diag,se,energy_qp_z,energy_qp_new,zz)
          call output_homolumo('GW+SOSEX2',occupation,energy_qp_new,nsemin,nsemax)
          call dump_out_energy_yaml('gw+sosex2 energies',energy_qp_new,nsemin,nsemax)
          deallocate(zz)
          deallocate(energy_qp_z)

          call g3w2_selfenergy(occupation,energy_g,c_matrix,wpol,se_g3w2)

          ! se_g3w2 contains SOSEX2+G3W2
          call se%reset()
          call se%add(se_gw)
          call se%add(se_g3w2)
          call se_g3w2%destroy()

        case(G3W2_NUMERICAL)
          call sosex_selfenergy(basis,occupation,energy_g,c_matrix,wpol,se_g3w2)

          call g3w2_selfenergy_real_grid(basis,occupation,energy_g,c_matrix,se_g3w2)
          call se%add(se_gw)
          call se%add(se_g3w2)
          call se_g3w2%destroy()
        case default
          call die('selfenergy_evaluation: selfenergy type not recognized')
        end select

        deallocate(energy_qp_new)
        call wpol%destroy()
        call se_gw%destroy()

      else
        !
        ! Imaginary frequencies implementation
        call gw_selfenergy_grid(basis,occupation,energy_g,c_matrix,se)
        call sox_selfenergy_imag_grid(occupation,energy_g,c_matrix,se)
        if( calc_type%selfenergy_approx == GWSOSEX .OR. calc_type%selfenergy_approx == G3W2 ) then
          call sosex_selfenergy_imag_grid(basis,occupation,energy_g,c_matrix,se)
        endif
        if( calc_type%selfenergy_approx == G3W2 ) then
          call g3w2_selfenergy_imag_grid(basis,occupation,energy_g,c_matrix,se)
        endif
        call se%pade_fit()
      endif
    endif


    !
    ! Selfenergy = PT2
    !
    if(   calc_type%selfenergy_approx == PT2       &
    .OR. calc_type%selfenergy_approx == ONE_RING  &
    .OR. calc_type%selfenergy_approx == SOX ) then

      call pt2_selfenergy(calc_type%selfenergy_approx,nstate,basis,occupation,energy_g,c_matrix,se,en_mbpt%mp2)

      if( ABS( en_mbpt%mp2 ) > 1.0e-8 ) then
        write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en_mbpt%mp2
        write(stdout,*)
        en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%mp2

        write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en_mbpt%total
        write(stdout,*)
      endif

    endif

    !
    ! Selfenergy = PT3 or 2-rings
    !
    if( calc_type%selfenergy_approx == PT3 .OR. calc_type%selfenergy_approx == TWO_RINGS ) then
      call pt3_selfenergy(calc_type%selfenergy_approx,calc_type%selfenergy_technique, &
                         nstate,basis,occupation,energy_g,c_matrix,se,en_mbpt%mp2)
    endif

    !
    ! selfenergy = GWPT3
    !
    if( calc_type%selfenergy_approx == GWPT3 ) then
      !
      ! First perform a standard GW calculation
      !
      call wpol%init(nstate,occupation,0)
      call read_spectral_function(wpol,reading_status)
      ! If reading has failed, then do the calculation
      if( reading_status /= 0 ) then
        call polarizability(.FALSE.,.TRUE.,basis,occupation,energy_w,c_matrix,en_mbpt%rpa,en_mbpt%gw,wpol)
      endif
      call gw_selfenergy(GW,nstate,basis,occupation,energy_g,c_matrix,wpol,se)

      !
      ! Second perform a PT3 calculation minus the ring diagrams
      !
      call se_gwpt3%init(calc_type%selfenergy_technique,energy_g)
      call pt3_selfenergy(GWPT3,calc_type%selfenergy_technique,nstate,basis,occupation,energy_g,c_matrix,se_gwpt3,en_mbpt%mp2)

      !
      ! Finally add up the contributions and then destroy the se_sox object
      !
      se%sigma(:,:,:) = se%sigma(:,:,:) + se_gwpt3%sigma(:,:,:)

      call se_gwpt3%destroy()

    endif


    !
    ! Final output the quasiparticle energies, the self-energy etc.
    !

    if( print_sigma_ ) then
      call write_selfenergy_omega('selfenergy_'//TRIM(selfenergy_tag),exchange_m_vxc_diag,occupation,energy_g,se)
    endif

    allocate(energy_qp_new(nstate,nspin))

    select case(calc_type%selfenergy_technique)
    case(EVSC)
      call find_qp_energy_linearization(se,exchange_m_vxc_diag,energy,energy_qp_new)
      call output_qp_energy(TRIM(selfenergy_tag),energy,exchange_m_vxc_diag,1,se,energy_qp_new)
    case(exact_dyson)
      ! Fake new QP energies in this case
      ! because it is not obvious to find which are the QP and which are not.
      energy_qp_new(:,:) = energy(:,:)
    case default
      allocate(energy_qp_z(nstate,nspin))
      allocate(zz(nstate,nspin))
      call find_qp_energy_linearization(se,exchange_m_vxc_diag,energy,energy_qp_z,zz)
      call find_qp_energy_graphical(se,exchange_m_vxc_diag,energy,energy_qp_new,zz)

      call output_qp_energy(TRIM(selfenergy_tag),energy,exchange_m_vxc_diag,1,se,energy_qp_z,energy_qp_new,zz)
      call output_qp_energy_yaml(TRIM(selfenergy_tag),energy,exchange_m_vxc_diag,se,energy_qp_z,energy_qp_new,zz)
      deallocate(zz)
      deallocate(energy_qp_z)
    end select

    if( calc_type%selfenergy_approx == GW ) then
      call selfenergy_convergence_prediction(basis,c_matrix,energy_qp_new)
    endif

    !
    ! Write the QP energies on disk: ENERGY_QP file
    !
    call write_energy_qp(energy_qp_new)
    call dump_out_energy_yaml(TRIM(selfenergy_tag)//' energies',energy_qp_new,nsemin,nsemax)

    !
    ! Output the new HOMO and LUMO energies
    !
    call output_homolumo(TRIM(selfenergy_tag),occupation,energy_qp_new,nsemin,nsemax)

    deallocate(energy_qp_new)

    !
    ! Deallocations
    !
    call se%destroy()

    ! Synchronization of all CPUs before going on
    call world%barrier()
  enddo ! nstep_gw

  deallocate(exchange_m_vxc_diag)

end subroutine selfenergy_evaluation


end module m_selfenergy_evaluation
!=========================================================================
