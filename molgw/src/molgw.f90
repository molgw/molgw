!=========================================================================
!
! MOLGW
! is a free software
! 
! All legal information can be found elsewhere
!
!=========================================================================
program molgw
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_tools
 use m_scf
 use m_atoms
 use m_gaussian
 use m_basis_set
 use m_eri
 use m_dft_grid
 use m_spectral_function
 use m_timedependent
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none

!=====
 type(basis_set)         :: basis
 type(basis_set)         :: auxil_basis
 type(basis_set)         :: prod_basis
 type(spectral_function) :: wpol
 integer                 :: reading_status
 integer                 :: ibf,jbf
 integer                 :: ispin,istate
 logical                 :: is_restart,is_big_restart
 character(len=100)      :: title
 real(dp)                :: energy_tmp
 real(dp),allocatable    :: hamiltonian_tmp(:,:)
 real(dp),allocatable    :: hamiltonian_kinetic(:,:)
 real(dp),allocatable    :: hamiltonian_nucleus(:,:)
 real(dp),allocatable    :: hamiltonian_exx(:,:,:)
 real(dp),allocatable    :: hamiltonian_xc(:,:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: s_matrix(:,:)
 real(dp),allocatable    :: c_matrix(:,:,:)
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: energy(:,:)
 real(dp),allocatable    :: occupation(:,:)
 real(dp),allocatable    :: exchange_m_vxc_diag(:,:)
!FBFB
 real(dp) :: rtmp(3)
 complex(dpc) :: ctmp(3)
 real(dp) :: qvec(3)
!=============================

 call init_mpi()

 call init_scalapack()
 !
 ! initialize the warning counters
 call init_warning()

 !
 ! start counting time here
 call init_timing()
 call start_clock(timing_total)
 call start_clock(timing_prescf)

 !
 ! Output some welcome message and compilation options
 call header()

 !
 ! Reading input file: the input parameters are stored in the module m_inputparam
 call read_inputfile_namelist()

 !
 ! Nucleus-nucleus repulsion contribution to the energy
 call nucleus_nucleus_energy(en%nuc_nuc)

 !
 ! Build up the basis set 
 !
 write(stdout,*) 'Setting up the basis set for wavefunctions'
 call init_basis_set(basis_path,basis_name,gaussian_type,basis)
 call setup_cart_to_pure_transforms(gaussian_type)

 !
 ! First attempt to distribute the work load among procs
 ! TODO: to be removed
 call init_distribution(basis%nbf)
 
 !
 ! Allocate the main arrays
 allocate(occupation(basis%nbf,nspin))
 allocate(energy(basis%nbf,nspin))
 allocate(c_matrix(basis%nbf,basis%nbf,nspin))
 allocate(s_matrix(basis%nbf,basis%nbf))
 allocate(p_matrix(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_kinetic(basis%nbf,basis%nbf))
 allocate(hamiltonian_nucleus(basis%nbf,basis%nbf))
 allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_exx(basis%nbf,basis%nbf,nspin) )
 allocate(hamiltonian_xc(basis%nbf,basis%nbf,nspin) )
 allocate(exchange_m_vxc_diag(basis%nbf,nspin))

 !
 ! Some required initializations
 energy(:,:) = 0.0_dp

 !
 ! Build up the overlap matrix S
 ! S only depends onto the basis set
 call setup_overlap(print_matrix_,basis,s_matrix)

 !
 ! Set up the electron repulsion integrals
 !
 ! ERI are stored "privately" in the module m_eri
 call prepare_eri(basis,0.0_dp,BUFFER1)
 if( .NOT. is_full_auxil) then
   call calculate_eri(print_eri_,basis,0.0_dp,BUFFER1)
 endif

! call refine_negligible_basispair()

 !
 ! for Range-separated hybrids, calculate the long-range ERI
 if(calc_type%need_exchange_lr) then
   call prepare_eri(basis,rcut,BUFFER2)
   if( .NOT. is_full_auxil) then
     call calculate_eri(print_eri_,basis,rcut,BUFFER2)
   endif
 endif
! call negligible_eri(1.0e-10_dp)

 !
 ! In case of GW or BSE run, set up the product basis 
 if( calc_type%is_gw .OR. calc_type%is_td .OR. calc_type%is_bse) call init_product_basis_set(basis,prod_basis)

 !
 ! Build the occupation array
 ! with zero temperature since we do not have the energy yet
 call set_occupation(basis%nbf,0.0_dp,electrons,magnetization,energy,occupation)

 !
 ! Try to read a RESTART file if it exists
 call read_any_restart(basis%nbf,occupation,c_matrix,energy,hamiltonian_exx,hamiltonian_xc,is_restart,is_big_restart)

 !
 ! Setup the grids for the quadrature of DFT potential/energy
 if( calc_type%is_dft .AND. .NOT. is_big_restart) then
   call setup_dft_grid()
 endif

 !
 ! Calculate the parts of the hamiltonian that does not change along
 ! with the SCF cycles
 if( .NOT. is_big_restart ) then
   !
   ! Kinetic energy contribution
   call setup_kinetic(print_matrix_,basis,hamiltonian_kinetic)
  
   !
   ! Nucleus-electron interaction
   call setup_nucleus(print_matrix_,basis,hamiltonian_nucleus)
#ifdef TODAY
   call setup_effective_core(print_matrix_,basis,hamiltonian_nucleus)
#endif
 endif

 if( .NOT. is_restart) then
   !
   ! Setup the initial c_matrix by diagonalizing the bare hamiltonian
   allocate(hamiltonian_tmp(basis%nbf,basis%nbf))
   !
   ! Calculate a very approximate vhxc based on simple gaussians placed on atoms
   call dft_approximate_vhxc(basis,hamiltonian_tmp)
   hamiltonian_tmp(:,:) = hamiltonian_tmp(:,:) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)

   write(stdout,'(/,a)') ' Diagonalization of an approximate hamiltonian'
   call diagonalize_generalized_sym(basis%nbf,hamiltonian_tmp,s_matrix,&
                                    energy(:,1),c_matrix(:,:,1))
   deallocate(hamiltonian_tmp)
   ! The hamiltonian is still spin-independent:
   c_matrix(:,:,nspin) = c_matrix(:,:,1)
  
   if( print_matrix_ ) then
     do ispin=1,nspin
       matrix_tmp(:,:,ispin) = TRANSPOSE( c_matrix(:,:,ispin) )
     enddo
     title='=== Initial C matrix ==='
     call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,matrix_tmp)
   endif

 endif

 !
 ! Setup the density matrix: p_matrix
 call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
!!
!! Test PSP = P
! call test_density_matrix(basis%nbf,nspin,p_matrix,s_matrix)

 title='=== 1st density matrix P ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,p_matrix)

 !
 ! Initialize the SCF mixing procedure
 call init_scf(basis%nbf)

 !
 ! If an auxiliary basis is given,
 ! then set it up now and calculate the required ERI: 2- and 3-center integrals
 !
 if( has_auxil_basis ) then
   write(stdout,*) 'Setting up the auxiliary basis set for Coulomb integrals'
   call init_basis_set(basis_path,auxil_basis_name,gaussian_type,auxil_basis)

   ! 2-center integrals
   call calculate_eri_2center(print_eri_,auxil_basis)
   ! Prepare the distribution of the 3-center integrals
   call distribute_auxil_basis(auxil_basis)
   ! 3-center integrals
   call calculate_eri_3center(print_eri_,basis,auxil_basis)

   ! If Range-Separated Hybrid are requested
   ! If is_big_restart, these integrals are NOT needed
   if(calc_type%need_exchange_lr .AND. .NOT. is_big_restart) then
     ! 2-center integrals
     call calculate_eri_2center_lr(print_eri_,auxil_basis,rcut)
     ! Prepare the distribution of the 3-center integrals
     call distribute_auxil_basis_lr(auxil_basis)
     ! 3-center integrals
     call calculate_eri_3center_lr(print_eri_,basis,auxil_basis,rcut)
   endif

 endif


 call stop_clock(timing_prescf)

 !
 ! Big SCF loop is in there
 ! Only do it if the calculation is NOT a big restart
 !
 if( .NOT. is_big_restart) then
   call scf_loop(basis,prod_basis,auxil_basis,                                           &
                 s_matrix,c_matrix,p_matrix,                                             &
                 hamiltonian_kinetic,hamiltonian_nucleus,hamiltonian_exx,hamiltonian_xc, &
                 occupation,energy)
 endif
 
 call start_clock(timing_postscf)

 if( print_wfn_ ) call plot_wfn(basis,c_matrix)
 if( print_wfn_ ) call plot_rho(basis,occupation,c_matrix)
 if( print_cube_ ) call plot_cube_wfn(basis,c_matrix)
 if( print_pdos_ ) call mulliken_pdos(basis,s_matrix,c_matrix,occupation,energy)


 !
 ! Some deallocations here
 !
 ! If an auxiliary basis is given, the 4-center integrals are not needed anymore
 if( has_auxil_basis ) call deallocate_eri_buffer()
 ! If RSH calculations were performed, then deallocate the LR integrals which
 ! are not needed anymore
 if( calc_type%need_exchange_lr .AND. .NOT. is_big_restart ) call deallocate_eri_buffer_lr()
 if( has_auxil_basis .AND. calc_type%need_exchange_lr ) call destroy_eri_3center_lr()

 !
 ! CI calculation is done here
 ! implemented for 2 electrons only!
 if(calc_type%is_ci) then
   if(nspin/=1) call die('for CI, nspin should be 1')
   if( ABS( electrons - 2.0_dp ) > 1.e-5_dp ) call die('CI is implemented for 2 electrons only')
   call full_ci_2electrons_spin(print_wfn_,0,basis,hamiltonian_kinetic+hamiltonian_nucleus,c_matrix,en%nuc_nuc)
 endif

 !
 ! Time Dependent calculations
 ! works for DFT, HF, and hybrid
 if(calc_type%is_td .OR. calc_type%is_bse) then

   if(calc_type%is_td .AND. calc_type%is_dft) call setup_dft_grid()
   if(has_auxil_basis) then
     call prepare_eri_3center_eigen(c_matrix)
     call destroy_eri_3center()
   endif
   call init_spectral_function(basis%nbf,occupation,wpol)
   call polarizability(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,en%rpa,wpol)
   call destroy_spectral_function(wpol)

   if(calc_type%is_td .AND. calc_type%is_dft) call destroy_dft_grid()
   if(has_auxil_basis) call destroy_eri_3center_eigen()

 endif
  



 exchange_m_vxc_diag(:,:) = 0.0_dp
 do ispin=1,nspin
   do istate=1,basis%nbf
     do ibf=1,basis%nbf
       do jbf=1,basis%nbf
         exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                 + c_matrix(ibf,istate,ispin) * ( hamiltonian_exx(ibf,jbf,ispin) - hamiltonian_xc(ibf,jbf,ispin) )&
                  * c_matrix(jbf,istate,ispin)
       enddo
     enddo
   enddo
 enddo

 !
 ! final evaluation for G0W0
 if( calc_type%is_gw .AND. &
       ( calc_type%gwmethod == GV .OR. calc_type%gwmethod == GSIGMA .OR.  calc_type%gwmethod == LW &
    .OR. calc_type%gwmethod == LW2 &
    .OR. calc_type%gwmethod == GSIGMA3 & ! FBFB LW testing purposes to be removed
    .OR. calc_type%gwmethod == G0W0 .OR. calc_type%gwmethod == COHSEX   &
    .OR. calc_type%gwmethod == GnW0 .OR. calc_type%gwmethod == GnWn ) ) then

   !
   ! A section under development for the range-separated RPA
   if( calc_type%is_lr_mbpt ) call die('lr_mbpt code removed. Does not exist anymore')

   if(has_auxil_basis) then
     call prepare_eri_3center_eigen(c_matrix)
     if( calc_type%gwmethod == LW .OR. calc_type%gwmethod == LW2 .OR. calc_type%gwmethod == GSIGMA ) &
         call prepare_eri_3center_eigen_mixed(c_matrix) ! FBFB LW
     call destroy_eri_3center()
   endif

   call init_spectral_function(basis%nbf,occupation,wpol)
   ! Try to read a spectral function file in order to skip the calculation
   call read_spectral_function(wpol,reading_status)
   ! If reading has failed, then do the calculation
   if( reading_status /= 0 ) then
     call polarizability(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,en%rpa,wpol)
   endif

   en%tot = en%tot + en%rpa
   if( calc_type%is_dft ) en%tot = en%tot - en%xc - en%exx_hyb + en%exx 
   write(stdout,'(/,a,f19.10)') ' RPA Total energy [Ha]: ',en%tot

   call gw_selfenergy(calc_type%gwmethod,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix_tmp,en%gw)
   if(has_auxil_basis) call destroy_eri_3center_eigen()

   if( ABS(en%gw) > 1.0e-5_dp ) then
     write(stdout,'(/,a,f19.10)') ' Galitskii-Migdal Total energy [Ha]: ',en%tot - en%rpa + en%gw
   endif

   title='=== Self-energy === (in the eigenstate basis)'
   call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,matrix_tmp)
   call destroy_spectral_function(wpol)

 endif ! G0W0

 !
 ! final evaluation for MP2
 if( calc_type%is_mp2 .AND. calc_type%gwmethod == perturbative ) then

! This routine is faster but only gives the correlation energy
   if(has_auxil_basis) then
     call prepare_eri_3center_eigen(c_matrix)
     call mp2_energy_ri(basis,occupation,energy,en%mp2)
     call destroy_eri_3center_eigen()
   else
     call mp2_energy_fast(basis,occupation,c_matrix,energy,en%mp2)
   endif

! This routine is slower but gives both the correlation energy and the self-energy
!   call mp2_selfenergy(calc_type%gwmethod,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,hamiltonian_exx,en%mp2)
   write(stdout,'(a,2x,f19.10)') ' MP2 Energy       [Ha]:',en%mp2
   write(stdout,*) 
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%mp2

   write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy [Ha]:',en%tot
   write(stdout,'(a,2x,f19.10)') ' SE+MP2  Total En [Ha]:',en%tot+en%se

   title='=== Self-energy === (in the eigenstate basis)'
   call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,matrix_tmp)

 endif


 !
 ! Cleanly exiting the code
 !
 deallocate(s_matrix,c_matrix,p_matrix)
 deallocate(hamiltonian_kinetic,hamiltonian_nucleus)
 deallocate(hamiltonian_exx,hamiltonian_xc)
 deallocate(energy,occupation)

 deallocate(matrix_tmp)
 deallocate(exchange_m_vxc_diag)

 call deallocate_eri()
 if(has_auxil_basis) call destroy_eri_3center()

 call destroy_basis_set(basis)
 if(has_auxil_basis) call destroy_basis_set(auxil_basis)
 if(calc_type%is_gw .OR. calc_type%is_td .OR. calc_type%is_bse ) call destroy_basis_set(prod_basis)
 call destroy_atoms()

 call total_memory_statement()

 call stop_clock(timing_postscf)
 call stop_clock(timing_total)

 call output_timing()

 call output_all_warnings()

 write(stdout,'(/,a,/)') ' This is the end'

 call finish_mpi()


end program molgw


!=========================================================================


