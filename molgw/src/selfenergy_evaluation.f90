!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the driver for the different self-energy methods:
! PT2, GW, evGW, COHSEX, GWGamma
!
!=========================================================================
subroutine selfenergy_evaluation(basis,auxil_basis,nstate,m_ham,n_ham,occupation,energy,c_matrix, &
                                 s_matrix,hamiltonian_exx,hamiltonian_xc)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_eri
 use m_eri_calculate
 use m_eri_lr_calculate
 use m_eri_ao_mo
 use m_dft_grid
 use m_scf,only: en
 use m_hamiltonian
 use m_spectral_function
 use m_selfenergy_tools
 use m_virtual_orbital_space
 implicit none

 type(basis_set),intent(in) :: basis
 type(basis_set),intent(in) :: auxil_basis
 integer,intent(in)         :: m_ham,n_ham
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(inout)     :: energy(nstate,nspin)
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: s_matrix(m_ham,n_ham)
 real(dp),intent(in)        :: hamiltonian_exx(m_ham,n_ham,nspin)
 real(dp),intent(in)        :: hamiltonian_xc(m_ham,n_ham,nspin)
!=====
 integer                 :: reading_status
 integer                 :: ispin
 type(spectral_function) :: wpol
 real(dp),allocatable    :: exchange_m_vxc_diag(:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: sigc(:,:)
 real(dp),allocatable    :: selfenergy_omega(:,:,:)
 real(dp)                :: energy_g(nstate,nspin)
 real(dp)                :: energy_w(nstate,nspin)
#ifdef ACTIVATE_EXPERIMENTAL
 real(dp),allocatable    :: p_matrix(:,:,:),p_matrix_sqrt(:,:,:),p_matrix_occ(:,:)
 integer                 :: istate
 real(dp)                :: exc
#endif
!=====

 !
 ! Prepare the diagonal of the matrix Sigma_x - Vxc
 ! for the forthcoming GW or PT2 corrections
 allocate(exchange_m_vxc_diag(nstate,nspin))
 call setup_exchange_m_vxc_diag(basis,nstate,m_ham,n_ham,occupation,c_matrix,hamiltonian_exx,hamiltonian_xc,exchange_m_vxc_diag)
 !
 ! Set the range of states on which to evaluate the self-energy
 call selfenergy_set_omega_grid()

 !
 ! If requested,
 ! prepare an optmized virtual subspace based on 
 ! Frozen Natural Orbitals technique
 if( is_virtual_fno ) then
   !
   ! Be aware that the energies and the c_matrix for virtual orbitals are altered after this point
   ! and until they are restored in destroy_fno
   !
   call virtual_fno(basis,nstate,occupation,energy,c_matrix)
 endif
 !
 ! Or alternatively use the small basis technique
 if( has_small_basis ) then
   if( .NOT. parallel_ham ) then  
     call setup_virtual_subspace(basis,nstate,occupation,energy,c_matrix)
   endif
 endif

 !
 ! Set the range after the change of the virtual space
 call selfenergy_set_state_range(nstate,occupation)


 allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nspin))


 !
 ! Choose the one-electron energies to use
 ! one_shot or EVGW

 if( calc_type%selfenergy_technique == EVSC ) then
   call read_energy_qp(nstate,energy_g,reading_status)
   if(reading_status/=0) then
     call issue_warning('File energy_qp not found: assuming 1st iteration')
     energy_g(:,:) = energy(:,:)
   endif

   if( calc_type%selfenergy_approx == GnWn ) then
     energy_w(:,:) = energy_g(:,:)
   else
     energy_w(:,:) = energy(:,:)
   endif

 else
   energy_g(:,:) = energy(:,:)
   energy_w(:,:) = energy(:,:)
 endif






 !
 ! selfenergy = GW or COHSEX
 !
 if( calc_type%is_gw .AND. &
       ( calc_type%selfenergy_approx == GV .OR. calc_type%selfenergy_approx == GSIGMA .OR.  calc_type%selfenergy_approx == LW &
    .OR. calc_type%selfenergy_approx == LW2 &
    .OR. calc_type%selfenergy_approx == G0W0_IOMEGA .OR. calc_type%selfenergy_approx == GWTILDE &
    .OR. calc_type%selfenergy_approx == G0W0 .OR. calc_type%selfenergy_approx == COHSEX   &
    .OR. calc_type%selfenergy_approx == GnW0 .OR. calc_type%selfenergy_approx == GnWn ) ) then

   call init_spectral_function(nstate,occupation,wpol)

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
     call polarizability(basis,auxil_basis,nstate,occupation,energy_w,c_matrix,en%rpa,wpol)
   endif

   en%tot = en%tot + en%rpa
   if( calc_type%is_dft ) en%tot = en%tot - en%xc - en%exx_hyb + en%exx 
   write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en%tot


   call gw_selfenergy(nstate,calc_type%selfenergy_approx,basis,occupation,energy_g,exchange_m_vxc_diag,c_matrix,wpol,selfenergy_omega,en%gw)

   if( ABS(en%gw) > 1.0e-5_dp ) then
     write(stdout,'(/,a,f19.10)') ' Galitskii-Migdal Total energy (Ha): ',en%tot - en%rpa + en%gw
   endif

   call destroy_spectral_function(wpol)

   if( is_virtual_fno .OR. has_small_basis ) then
     call destroy_fno(basis,nstate,energy,c_matrix)
   endif

 endif

 !
 ! GWGamma
 !
 if( calc_type%selfenergy_approx == G0W0GAMMA0 .OR. calc_type%selfenergy_approx == G0W0SOX0 ) then
   call init_spectral_function(nstate,occupation,wpol)
   call read_spectral_function(wpol,reading_status)
   ! If reading has failed, then do the calculation
   if( reading_status /= 0 ) then
     call polarizability(basis,auxil_basis,nstate,occupation,energy_w,c_matrix,en%rpa,wpol)
   endif

   call gw_selfenergy(nstate,G0W0,basis,occupation,energy_g,exchange_m_vxc_diag,c_matrix,wpol,selfenergy_omega,en%gw)
   call gwgamma_selfenergy(nstate,basis,occupation,energy_g,exchange_m_vxc_diag,c_matrix,wpol,selfenergy_omega)
   call destroy_spectral_function(wpol)
 endif


 !
 ! Selfenergy = PT2
 ! 
 if( calc_type%selfenergy_approx == PT2 ) then

   allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
   call mp2_selfenergy(nstate,basis,occupation,energy_g,exchange_m_vxc_diag,c_matrix,s_matrix,matrix_tmp,en%mp2)
   deallocate(matrix_tmp)
   if( ABS( en%mp2 ) > 1.0e-8 ) then
     write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en%mp2
     write(stdout,*)
     en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%mp2

     write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en%tot
     write(stdout,'(a,2x,f19.10)') ' SE+MP2  Total En (Ha):',en%tot+en%se
     write(stdout,*)
   endif

 endif


 !
 ! EXPERIMENTAL COHSEX implementation
 ! final evaluation for perturbative COHSEX
 !
 if( calc_type%selfenergy_approx == COHSEX_DEVEL .OR. calc_type%selfenergy_approx == TUNED_COHSEX ) then

   if( .NOT. has_auxil_basis ) call die('cohsex needs an auxiliary basis')
   call init_spectral_function(nstate,occupation,wpol)
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1)
   !
   ! Calculate v^{1/2} \chi v^{1/2}
   call static_polarizability(nstate,occupation,energy_w,wpol)

   call destroy_eri_3center_eigen()

   !
   allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
   allocate(sigc(nstate,nspin))

#ifdef ACTIVATE_EXPERIMENTAL
   ! Calculate the DFT potential part
   if( ABS( delta_cohsex ) > 1.0e-6_dp ) then

     allocate(p_matrix(basis%nbf,basis%nbf,nspin))
     allocate(p_matrix_sqrt(basis%nbf,basis%nbf,nspin))
     allocate(p_matrix_occ(basis%nbf,nspin))
     call init_dft_grid(grid_level)
     call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
     call setup_sqrt_density_matrix(basis%nbf,p_matrix,p_matrix_sqrt,p_matrix_occ)

     ! Override the DFT XC correlation settings
     call init_dft_type('HJSx',calc_type)
#ifdef HAVE_LIBXC
     call xc_f90_gga_x_hjs_set_par(calc_type%xc_func(1),1.0_dp/rcut_mbpt)
#endif
     call dft_exc_vxc(basis,p_matrix_occ,p_matrix_sqrt,p_matrix,matrix_tmp,exc)
 
     write(stdout,*) '===== SigX SR ======'
     do ispin=1,nspin
       do istate=1,nstate
         sigc(istate,ispin) = DOT_PRODUCT( c_matrix(:,istate,ispin) , &
                                   MATMUL( matrix_tmp(:,:,ispin) , c_matrix(:,istate,ispin ) ) )
         write(stdout,*) istate,ispin,sigc(istate,ispin) * Ha_eV
       enddo
     enddo
     write(stdout,*) '===================='
     sigc(istate,ispin) = sigc(istate,ispin) * delta_cohsex 

     deallocate(p_matrix)
     deallocate(p_matrix_sqrt)
     deallocate(p_matrix_occ)
     call destroy_dft_grid()

   else

     sigc(:,:) = 0.0_dp

   endif

#endif

   call cohsex_selfenergy(nstate,basis,occupation,energy_g,exchange_m_vxc_diag, & 
                          c_matrix,s_matrix,wpol,matrix_tmp,sigc,en%gw)


   !
   ! A section under development for the range-separated RPA
   if( calc_type%is_lr_mbpt ) then

     ! 2-center integrals
     call calculate_eri_2center_lr(auxil_basis,rcut_mbpt)
     ! Prepare the distribution of the 3-center integrals
     call distribute_auxil_basis_lr(nauxil_2center_lr,nauxil_3center_lr)
     ! 3-center integrals
     call calculate_eri_3center_lr(basis,auxil_basis,rcut_mbpt)

     call cohsex_selfenergy_lr(nstate,basis,occupation,energy_g,exchange_m_vxc_diag, &
                               c_matrix,s_matrix,wpol,matrix_tmp,sigc,en%gw)
   endif

   deallocate(matrix_tmp)
   deallocate(sigc)

 endif ! COHSEX
 !
 ! end of EXPERIMENTAL COHSEX implementation
 !





 !
 ! Deallocations
 !
 call selfenergy_destroy_omega_grid()
 if( ALLOCATED(exchange_m_vxc_diag) ) deallocate(exchange_m_vxc_diag)


end subroutine selfenergy_evaluation


!=========================================================================
