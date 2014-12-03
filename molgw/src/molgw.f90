#include "macros.h"
!=========================================================================
program molgw
 use m_definitions
 use m_mpi
 use m_timing
 use m_warning
 use m_calculation_type
 use m_inputparam
 use m_tools
 use m_scf
 use m_atoms
 use m_gaussian
 use m_basis_set
 use m_eri
 use m_dft_grid
 use m_spectral_function
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
 integer                 :: ispin,istate,ncore
 logical                 :: file_exists
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
!=============================

 call init_mpi()

 call init_scalapack()
 !
 ! initialize the warning counter
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
 call read_inputparameter_molecule()

 !
 ! Nucleus-nucleus repulsion contribution to the energy
 call nucleus_nucleus_energy(en%nuc_nuc)

 !
 ! Build up the basis set and the auxiliary basis set if needed
 !
 call init_basis_set(print_basis,basis_name      ,gaussian_type,basis)
 call init_basis_set(print_basis,auxil_basis_name,gaussian_type,auxil_basis)
 call setup_cart_to_pure_transforms(MAX(basis%ammax,auxil_basis%ammax),gaussian_type)

 !
 ! First attempt to distribute the work load among procs
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
 energy(:,:)            = 0.0_dp

 !
 ! Build up the overlap matrix S
 ! S only depends onto the basis set
 call setup_overlap(print_matrix,basis,s_matrix)

 !
 ! Set up the electron repulsion integrals
 !
 ! ERI are stored "privately" in the module m_eri
 call start_clock(timing_eri)
 call allocate_eri(basis,0.0_dp,BUFFER1)
 call calculate_eri(print_eri,basis,0.0_dp,BUFFER1)
 call stop_clock(timing_eri)

 call refine_negligible_basispair()

 !
 ! for HSE functionals, calculate the long-range ERI
 if(calc_type%is_screened_hybrid) then
   call allocate_eri(basis,rcut,BUFFER2)
   call start_clock(timing_eri)
   call calculate_eri(print_eri,basis,rcut,BUFFER2)
   call stop_clock(timing_eri)
 endif
! call negligible_eri(1.0e-10_dp)

 !
 ! In case of GW or BSE run, set up the product basis 
 if( calc_type%is_gw .OR. calc_type%is_td .OR. calc_type%is_bse) call init_product_basis_set(basis,prod_basis)

 !
 ! Build the occupation array
 call set_occupation(electrons,magnetization,basis%nbf,nspin,occupation)
 title='=== Occupations ==='
 call dump_out_occupation(title,basis%nbf,nspin,occupation)

 !
 ! Calculate the parts of the hamiltonian that does not change along
 ! with the SCF cycles
 !
 ! Kinetic energy contribution
 call setup_kinetic(print_matrix,basis,hamiltonian_kinetic)

 !
 ! Nucleus-electron interaction
 call setup_nucleus(print_matrix,basis,hamiltonian_nucleus)

 !
 ! Setup the initial c_matrix by diagonalizing the bare hamiltonian
 allocate(hamiltonian_tmp(basis%nbf,basis%nbf))
 hamiltonian_tmp(:,:) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
 WRITE_MASTER(*,*) 'Diagonalization of the bare hamiltonian'
 call diagonalize_generalized_sym(basis%nbf,hamiltonian_tmp,s_matrix,&
                                  energy(:,1),c_matrix(:,:,1))
 deallocate(hamiltonian_tmp)
 ! The hamiltonian is still spin-independent:
 c_matrix(:,:,nspin) = c_matrix(:,:,1)

! call guess_starting_c_matrix(basis%nbf,nspin,c_matrix)
! call guess_starting_c_matrix_new(basis,nspin,c_matrix)
 if( print_matrix ) then
   do ispin=1,nspin
     matrix_tmp(:,:,ispin) = transpose( c_matrix(:,:,ispin) )
   enddo
   title='=== Initial C matrix ==='
   call dump_out_matrix(print_matrix,title,basis%nbf,nspin,matrix_tmp)
 endif


 !
 ! Setup the density matrix: p_matrix
 call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
 !
 ! Read the density matrix if asked and override the previously guessed matrix
 call read_density_matrix(basis%nbf,nspin,p_matrix)
!!
!! Test PSP = P
! call test_density_matrix(basis%nbf,nspin,p_matrix,s_matrix)

 title='=== 1st density matrix P ==='
 call dump_out_matrix(print_matrix,title,basis%nbf,nspin,p_matrix)

 !
 ! Initialize the SCF mixing procedure
 call init_scf(nscf,basis%nbf,nspin,alpha_mixing)

 !
 ! Setup the grids for the quadrature of DFT potential/energy
 if( calc_type%is_dft ) then
   call setup_dft_grid()
 endif

 call stop_clock(timing_prescf)
 !
 ! Big SCF loop is in there
 !
 call scf_loop(basis,auxil_basis,prod_basis,s_matrix,c_matrix,p_matrix,                &
               hamiltonian_kinetic,hamiltonian_nucleus,hamiltonian_exx,hamiltonian_xc, &
               occupation,energy)
 
 call start_clock(timing_postscf)


 if( print_densitymatrix )   call write_density_matrix(nspin,basis%nbf,p_matrix)
 if( print_wfn ) call plot_wfn(nspin,basis,c_matrix)
! if( print_wfn ) call plot_cube_wfn(nspin,basis,c_matrix)


 !
 ! If an auxiliary basis is set up, 
 ! calculate the required ERI: 2- and 3-center integrals
 if( is_auxil_basis ) then
   call deallocate_eri_buffer()
   call allocate_eri_auxil(auxil_basis)
   ! 2-center integrals
   call calculate_eri_2center(print_eri,auxil_basis)
   ! 3-center integrals
   call calculate_eri_3center(print_eri,basis,auxil_basis)
 endif

 !
 ! CI calculation is done here
 ! implemented for 2 electrons only!
 if(calc_type%is_ci) then
   if(nspin/=1) stop'for CI, nspin should be 1'
   if( ABS( electrons - 2.0_dp ) > 1.e-5_dp ) stop'CI is implemented for 2 electrons only'
   call full_ci_2electrons_spin(print_wfn,0,basis,hamiltonian_kinetic+hamiltonian_nucleus,c_matrix,en%nuc_nuc)
 endif

 !
 ! Time Dependent calculations
 ! works for DFT, HF, and hybrid
 if(calc_type%is_td .OR. calc_type%is_bse) then
   call init_spectral_function(basis%nbf,prod_basis%nbf,nspin,occupation,wpol)

   ! For BSE calculation, obtain the wpol object from a previous calculation
   if(calc_type%is_bse) then
     call read_spectral_function(wpol,reading_status)
     if(reading_status/=0) then 
       stop'BSE requires a previous GW calculation stored in a spectral_file'
     endif
   endif
   call polarizability_td(basis,prod_basis,occupation,energy,c_matrix,wpol)
   call destroy_spectral_function(wpol)
 endif
  
!%!   inquire(file='manual_coresplitting',exist=file_exists)
!%!   if(file_exists) then
!%!     WRITE_MASTER(*,*) 'TESTING CORE-VALENCE SPLITTING'
!%!     open(13,file='manual_coresplitting')
!%!     read(13,*) ncore
!%!     close(13)
!%!     WRITE_MASTER(msg,'(a,i4,2x,i4)') 'core-valence splitting switched on up to state = ',ncore
!%!     call issue_warning(msg)
!%!     do istate=1,ncore
!%!       occupation(istate,:) = 0.0_dp
!%!     enddo
!%!     call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
!%!     call dft_exc_vxc(nspin,basis,ndft_xc,dft_xc_type,dft_xc_coef,p_matrix,ehomo,hamiltonian_xc,en%xc)
!%!   endif

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
 if( calc_type%is_gw .AND. ( calc_type%gwmethod == perturbative .OR. calc_type%gwmethod == COHSEX ) ) then

   !
   ! A section under development for the range-separated RPA
   if( calc_type%is_lr_mbpt ) call testing_tobe_removed()


   call init_spectral_function(basis%nbf,prod_basis%nbf,nspin,occupation,wpol)
   call start_clock(timing_pola)
#ifdef CASIDA
   call polarizability_casida(nspin,basis,prod_basis,occupation,energy,c_matrix,en%rpa,wpol)
#else
   call polarizability_rpa(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,en%rpa,wpol)
#endif
   call stop_clock(timing_pola)
   en%tot = en%tot + en%rpa
   if( calc_type%is_dft ) en%tot = en%tot - en%xc + en%exx * ( 1.0_dp - alpha_hybrid )
   WRITE_MASTER(*,'(/,a,f16.10)') ' RPA Total energy [Ha]: ',en%tot

   call start_clock(timing_self)
#ifndef CASIDA
   call gw_selfenergy(calc_type%gwmethod,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix_tmp)
#endif
   call stop_clock(timing_self)

   title='=== Self-energy === (in the orbital basis)'
   call dump_out_matrix(print_matrix,title,basis%nbf,nspin,matrix_tmp)
   call destroy_spectral_function(wpol)

 endif ! G0W0

 !
 ! final evaluation for MP2
 if( calc_type%is_mp2 .AND. calc_type%gwmethod == perturbative ) then

   call setup_exchange(print_matrix,basis%nbf,nspin,p_matrix,matrix_tmp,en%exx)
   WRITE_MASTER(*,*) 'EXX     [Ha]:',en%exx
#ifdef CASIDA
   call start_clock(timing_mp2_energy)
!   call mp2_energy(nspin,basis,occupation,c_matrix,energy,en%mp2)
   call mp2_energy_fast(nspin,basis,occupation,c_matrix,energy,en%mp2)
   call stop_clock(timing_mp2_energy)
#else
   call start_clock(timing_mp2_self)
   call mp2_selfenergy(calc_type%gwmethod,nspin,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,matrix_tmp,en%mp2)
   call stop_clock(timing_mp2_self)
#endif
   WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Energy       [Ha]:',en%mp2
   WRITE_MASTER(*,*) 
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%mp2

   WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Total Energy [Ha]:',en%tot
   WRITE_MASTER(*,'(a,2x,f16.10)') ' SE+MP2  Total En [Ha]:',en%tot+en%se

   title='=== Self-energy === (in the orbital basis)'
   call dump_out_matrix(print_matrix,title,basis%nbf,nspin,matrix_tmp)

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
 if( calc_type%is_dft ) call destroy_dft_grid()

 call destroy_basis_set(basis)
 if(is_auxil_basis) call destroy_basis_set(auxil_basis)
 if(calc_type%is_gw .OR. calc_type%is_td .OR. calc_type%is_bse ) call destroy_basis_set(prod_basis)

 call stop_clock(timing_postscf)
 call stop_clock(timing_total)

 call output_timing()

 call output_all_warnings()

 WRITE_MASTER(*,'(/,a,/)') ' This is the end'

 call finish_mpi()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BEGIN OF SECTION TO BE REMOVED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testing_tobe_removed()
 real(dp),allocatable    :: ehomo(:)
!=====
 allocate(ehomo(nspin))
   write(*,*) '========== FABIEN ============'
   call setup_exchange(print_matrix,basis%nbf,nspin,p_matrix,matrix_tmp,en%exx)
   WRITE_MASTER(*,*) 'EXX [Ha]:',en%exx
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * matrix_tmp(ibf,jbf,ispin) &
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
       write(*,*) istate,ispin,exchange_m_vxc_diag(istate,ispin)*Ha_eV
     enddo
   enddo
   write(*,*) '========== END FABIEN ========'
   !
   ! Hard-coded cutoff radius
   rcut= 9.0909_dp ! 0.5000_dp
   write(msg,'(a,f10.4)') 'hard coded cutoff radius  ',rcut
   call issue_warning(msg)
   call deallocate_eri()
   call allocate_eri(basis,rcut,BUFFER1)
   call calculate_eri(print_matrix,basis,rcut,BUFFER1)

   write(*,*) '========== FABIEN ============'
   call setup_exchange(print_matrix,basis%nbf,nspin,p_matrix,matrix_tmp,en%exx)
   WRITE_MASTER(*,*) 'EXX [Ha]:',en%exx
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * matrix_tmp(ibf,jbf,ispin) &
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
       write(*,*) istate,ispin,exchange_m_vxc_diag(istate,ispin)*Ha_eV
     enddo
   enddo
   write(*,*) '========== END FABIEN ========'
   write(*,*) '========== FABIEN ============'
   if( .NOT. ALLOCATED( hamiltonian_xc ) ) allocate( hamiltonian_xc(basis%nbf,basis%nbf,nspin) )
   if( ndft_xc == 0 ) call setup_dft_grid()
    call dft_exc_vxc(nspin,basis,1,(/XC_GGA_X_PBE/),(/1.0_dp/),p_matrix,ehomo,hamiltonian_xc,energy_tmp)
#ifdef HAVE_LIBXC
   call dft_exc_vxc(nspin,basis,1,(/XC_HYB_GGA_XC_HSE06/),(/1.0_dp/),p_matrix,ehomo,hamiltonian_xc,energy_tmp)
#endif
   WRITE_MASTER(*,'(/,a,f16.10)') '    PBEx energy [Ha]: ',energy_tmp
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * ( hamiltonian_xc(ibf,jbf,ispin) )&
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
       write(*,*) istate,ispin,exchange_m_vxc_diag(istate,ispin)*Ha_eV
     enddo
   enddo
   write(*,*) '========== END FABIEN ========'

   if( .NOT. ALLOCATED( hamiltonian_xc ) ) allocate( hamiltonian_xc(basis%nbf,basis%nbf,nspin) )
   call dft_exc_vxc(nspin,basis,1,(/1000/),(/1.0_dp/),p_matrix,ehomo,hamiltonian_xc,energy_tmp)
   WRITE_MASTER(*,'(/,a,f16.10)') '    RPA LDA energy [Ha]: ',energy_tmp
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * ( hamiltonian_xc(ibf,jbf,ispin) )&
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
     enddo
   enddo
   WRITE_MASTER(*,*) 'test'
   WRITE_MASTER(*,*) '<vxc RPA FR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
   WRITE_MASTER(*,*) 'test'

   hamiltonian_xc(:,:,:) = 0.0_dp
   call dft_exc_vxc(nspin,basis,1,(/1005/),(/1.0_dp/),p_matrix,ehomo,hamiltonian_xc,energy_tmp)
   WRITE_MASTER(*,'(/,a,f16.10)') ' LR-RPA LDA energy [Ha]: ',energy_tmp
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * ( hamiltonian_xc(ibf,jbf,ispin) )&
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
     enddo
   enddo
   WRITE_MASTER(*,*) 'test'
   WRITE_MASTER(*,*) '<vxc RPA LR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
   WRITE_MASTER(*,*) 'test'

   call setup_exchange(print_matrix,basis%nbf,nspin,p_matrix,matrix_tmp,en%exx)  
   WRITE_MASTER(*,*) 'LR EXX [Ha]:',en%exx
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * matrix_tmp(ibf,jbf,ispin)&
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
     enddo
   enddo
   WRITE_MASTER(*,*) 'test'
   WRITE_MASTER(*,*) '<sigx LR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
   WRITE_MASTER(*,*) 'test'

   hamiltonian_xc(:,:,:) = 0.0_dp
   call dft_exc_vxc(nspin,basis,1,(/2001/),(/1.0_dp/),p_matrix,ehomo,hamiltonian_xc,energy_tmp)
   WRITE_MASTER(*,'(/,a,f16.10)') ' SR-EXX GGA energy [Ha]: ',energy_tmp
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * ( hamiltonian_xc(ibf,jbf,ispin))&
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
     enddo
   enddo
   WRITE_MASTER(*,*) 'test'
   WRITE_MASTER(*,*) '<vxc EXX SR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
   WRITE_MASTER(*,*) 'test'



   exchange_m_vxc_diag(:,:) = 0.0_dp

 deallocate(ehomo)
end subroutine testing_tobe_removed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END OF SECTION TO BE REMOVED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program molgw
!=========================================================================
