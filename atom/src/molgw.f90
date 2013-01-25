#include "macros.h"
!=========================================================================
program molgw
 use m_definitions
 use m_mpi
 use m_timing
 use m_warning
 use m_calculation_type
 use m_tools
 use m_scf
 use m_atoms
 use m_gaussian
 use m_basis_set
 use m_eri
 use m_gw
#ifdef OPENMP
 use omp_lib
#endif
 implicit none

 !
 ! Input parameters will be set in read_inputparameters
 type(calculation_type)       :: calc_type
 integer                      :: nspin,nscf
 real(dp)                     :: alpha_mixing
 integer                      :: print_volume
 character(len=100)           :: basis_name
 integer                      :: gaussian_type
 real(dp)                     :: electrons
 real(dp)                     :: magnetization
!===== variables for testing
 type(gaussian) :: gatmp,gbtmp
 type(basis_function) :: bftmp1,bftmp2
 real(dp) :: rtmp
!=====
 type(basis_set)         :: basis
 type(basis_set)         :: prod_basis
 type(spectral_function) :: wpol
 integer                 :: ibf,jbf,kbf,lbf,ijbf,klbf
 integer                 :: ispin,iscf,istate,jstate,astate,iatom,ncore
 logical                 :: scf_loop_convergence
 logical                 :: file_exists
 character(len=100)      :: title
 real(dp)                :: spin_fact
 real(dp)                :: energy_tmp,overlap_tmp
 real(dp)                :: dipole(3)
 real(dp),allocatable    :: ehomo(:),elumo(:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_xc(:,:,:)
 real(dp),allocatable    :: hamiltonian_kinetic(:,:)
 real(dp),allocatable    :: hamiltonian_nucleus(:,:)
 real(dp),allocatable    :: matrix(:,:,:)
 real(dp),allocatable    :: matrix3(:,:,:)
 real(dp),allocatable    :: vxc_matrix(:,:,:)
 real(dp),allocatable    :: s_matrix(:,:)
 real(dp),allocatable    :: c_matrix(:,:,:)
 real(dp),allocatable    :: p_matrix(:,:,:),p_matrix_old(:,:,:)
 real(dp),allocatable    :: energy(:,:)
 real(dp),allocatable    :: occupation(:,:)
 real(dp),allocatable    :: s_filtered_basis(:,:)
 real(dp),allocatable    :: sinv_filtered_basis(:,:)
 real(dp),allocatable    :: v_filtered_basis(:,:)
 real(dp),allocatable    :: sinv_v_sinv_filtered(:,:),sinv_v_sinv(:,:)
 real(dp),allocatable    :: exchange_m_vxc_diag(:,:)
 real(dp),allocatable    :: eri2(:,:),matrix2(:,:)
 real(dp),allocatable    :: self_energy_old(:,:,:)
! tmp elements to be removed eventually
 real(dp),allocatable    :: eigval(:),eigvec(:,:),matrix_tmp(:,:) 
!=====
 type energy_contributions
   real(dp) :: nuc_nuc= 0.0_dp
   real(dp) :: kin    = 0.0_dp
   real(dp) :: nuc    = 0.0_dp
   real(dp) :: hart   = 0.0_dp
   real(dp) :: exx    = 0.0_dp
   real(dp) :: xc     = 0.0_dp
   real(dp) :: se     = 0.0_dp      ! single-excitation contribution
   real(dp) :: mp2    = 0.0_dp
   real(dp) :: rpa    = 0.0_dp
   real(dp) :: tot    = 0.0_dp
 end type
!=====
 type(energy_contributions) :: en
!=============================

 call init_mpi()

 call header()

 !
 ! start counting time here
 call init_timing()
 call start_clock(timing_total)
 !
 ! initialize the warning counter
 call init_warning()
#ifdef CHI0
 msg='CHI0 preprocessing option is set'
 call issue_warning(msg)
#endif
#ifdef OPENMP
 write(msg,'(i6)') OMP_get_max_threads()
 msg='OPENMP option is activated with threads number'//msg
 call issue_warning(msg)
#endif
#ifdef LOW_MEMORY1
 msg='LOW_MEMORY version 1 option is swichted on'
 call issue_warning(msg)
#endif
#ifdef LOW_MEMORY2
 msg='LOW_MEMORY version 2 option is swichted on'
 call issue_warning(msg)
#endif
#ifdef MPI
#ifndef LOW_MEMORY1
 stop'MPI and LOW_MEMORY1 should be set to on'
#endif
#endif


 !
 ! Reading input file
 call read_inputparameter_molecule(calc_type,nspin,nscf,alpha_mixing,print_volume,&
                                   basis_name,gaussian_type,electrons,magnetization)

 !
 ! Nucleus-nucleus repulsion contribution to the energy
 call nucleus_nucleus_energy(en%nuc_nuc)

 !
 ! Build up the basis set
 !
 call init_basis_set(print_volume,basis_name,gaussian_type,basis)
 call setup_cart_to_pure_transforms(basis)

 !
 ! First attempt to distribute the work load among procs
 call init_distribution(basis%nbf)
 
 !
 ! Allocate the main arrays
 allocate(hamiltonian(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_xc(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_kinetic(basis%nbf,basis%nbf))
 allocate(hamiltonian_nucleus(basis%nbf,basis%nbf))
 allocate(matrix(basis%nbf,basis%nbf,nspin))
 allocate(c_matrix(basis%nbf,basis%nbf,nspin))
 allocate(s_matrix(basis%nbf,basis%nbf))
 allocate(p_matrix(basis%nbf,basis%nbf,nspin))
 allocate(p_matrix_old(basis%nbf,basis%nbf,nspin))
 allocate(energy(basis%nbf,nspin))
 allocate(occupation(basis%nbf,nspin))
 allocate(exchange_m_vxc_diag(basis%nbf,nspin))
 allocate(self_energy_old(basis%nbf,basis%nbf,nspin))
 allocate(ehomo(nspin))
 allocate(elumo(nspin))
 if( ndft_xc /= 0 ) allocate( vxc_matrix(basis%nbf,basis%nbf,nspin) )

 !
 ! Some required initializations
 energy(:,:)            = 0.0_dp
 self_energy_old(:,:,:) = 0.0_dp

 !
 ! Build up the overlap matrix S
 ! S only depends onto the basis set
 call setup_overlap(print_volume,basis,s_matrix)

 !
 ! Set up the electron repulsion integrals
 !
 ! ERI are stored "privately" in the module m_eri
 call start_clock(timing_integrals)
 call allocate_eri(basis%nbf)
 call calculate_eri(print_volume,basis,0.0_dp,BUFFER1)
 !
 ! for HSE functionals, calculate the long-range ERI
 if(calc_type%is_screened_hybrid) then
   call allocate_eri_lr(basis%nbf)
   call calculate_eri(print_volume,basis,rcut,BUFFER2)
 endif
 call stop_clock(timing_integrals)
! call negligible_eri(1.0e-10_dp)

 !
 ! In case of GW run, set up the product basis 
 if( calc_type%is_gw ) call init_product_basis_set(basis,prod_basis)

!========================================================
! AUXILIARY basis set GW
!========================================================
! not supported any longer...
! should be removed in the future (or not)
!
#ifdef AUXIL_BASIS
 if( calc_type%is_gw ) then
   call start_clock(timing_prodbasis)
   allocate(v_filtered_basis(prod_basis%nbf_filtered,prod_basis%nbf_filtered))
   allocate(s_filtered_basis(prod_basis%nbf_filtered,prod_basis%nbf_filtered))
   allocate(sinv_filtered_basis(prod_basis%nbf_filtered,prod_basis%nbf_filtered))

   allocate(eri2(prod_basis%nbf,prod_basis%nbf))
   do klbf=1,prod_basis%nbf
     do ijbf=1,prod_basis%nbf
       ibf = prod_basis%index_ij(1,ijbf)
       jbf = prod_basis%index_ij(2,ijbf)
       kbf = prod_basis%index_ij(1,klbf)
       lbf = prod_basis%index_ij(2,klbf)
       eri2(ijbf,klbf) = eri(ibf,jbf,kbf,lbf)
     enddo
   enddo

   v_filtered_basis = MATMUL( prod_basis%rotation , MATMUL( eri2 , TRANSPOSE(prod_basis%rotation) ) )
   deallocate(eri2)

   title='=== Coulomb matrix in product basis ==='
   call dump_out_matrix(print_volume,title,prod_basis%nbf_filtered,1,v_filtered_basis)

   allocate(matrix2(prod_basis%nbf,prod_basis%nbf))
   do klbf=1,prod_basis%nbf
     do ijbf=1,klbf ! prod_basis%nbf
       call overlap_basis_function(prod_basis%bf(ijbf),prod_basis%bf(klbf),overlap_tmp)
       matrix2(ijbf,klbf) = overlap_tmp
       matrix2(klbf,ijbf) = overlap_tmp
     enddo
   enddo
   s_filtered_basis = MATMUL( prod_basis%rotation , MATMUL( matrix2 , TRANSPOSE(prod_basis%rotation) ) )
   deallocate(matrix2)

   ! Remember that the product basis is unnormalized !
   title='=== product overlap matrix S ==='
   call dump_out_matrix(print_volume,title,prod_basis%nbf_filtered,1,s_filtered_basis)

   !
   ! calculate S^{-1}
   call invert(prod_basis%nbf_filtered,s_filtered_basis,sinv_filtered_basis)

   !
   ! S^-1 V S^-1 is first calculated on the filtered product basis
   allocate(sinv_v_sinv_filtered(prod_basis%nbf_filtered,prod_basis%nbf_filtered))
   sinv_v_sinv_filtered = matmul( v_filtered_basis , sinv_filtered_basis )
   sinv_v_sinv_filtered = matmul( sinv_filtered_basis , sinv_v_sinv_filtered )
   title='=== S-1 V S-1 ==='
   call dump_out_matrix(print_volume,title,prod_basis%nbf_filtered,1,sinv_v_sinv_filtered)

   !
   ! S^-1 V S^-1 is then transposed to the full product basis
   deallocate(v_filtered_basis,s_filtered_basis,sinv_filtered_basis)
   allocate(sinv_v_sinv(prod_basis%nbf,prod_basis%nbf))
   sinv_v_sinv = MATMUL( TRANSPOSE(prod_basis%rotation), MATMUL( sinv_v_sinv_filtered , prod_basis%rotation ) )

   deallocate(sinv_v_sinv_filtered)

   call stop_clock(timing_prodbasis)
 endif
#endif
!========================================================
! end of AUXILIARY basis set GW
!========================================================
   

 !
 ! Build occupation array, initial c_matrix
 call set_occupation(electrons,magnetization,basis%nbf,nspin,occupation)
 title='=== Occupations ==='
 call dump_out_occupation(title,basis%nbf,nspin,occupation)
! call guess_starting_c_matrix(basis%nbf,nspin,c_matrix)
 call guess_starting_c_matrix_new(basis,nspin,c_matrix)
 do ispin=1,nspin
   matrix(:,:,ispin) = transpose( c_matrix(:,:,ispin) )
 enddo
 title='=== Initial C matrix ==='
 call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)

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
 call dump_out_matrix(print_volume,title,basis%nbf,nspin,p_matrix)

 !
 ! Initialize the SCF mixing procedure
 call init_scf(nscf,basis%nbf,nspin,alpha_mixing)


 !
 ! Calculate the parts of the hamiltonian that does not change along
 ! with the SCF cycles
 !
 ! Kinetic energy contribution
 call setup_kinetic(print_volume,basis,hamiltonian_kinetic)

 !
 ! Nucleus-electron interaction
 call setup_nucleus(print_volume,basis,hamiltonian_nucleus)

! call setup_initial_c_matrix(print_volume,basis%nbf,nspin,hamiltonian_nucleus,s_matrix,occupation,c_matrix)
!
! call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
! title='=== 1st density matrix P ==='
! call dump_out_matrix(print_volume,title,basis%nbf,nspin,p_matrix)


 !
 ! start the big scf loop
 !
 do iscf=1,nscf
   WRITE_MASTER(*,'(/,a)') '-------------------------------------------'
   WRITE_MASTER(*,'(a,x,i4,/)') ' *** SCF cycle No:',iscf

   call output_homolumo(basis%nbf,nspin,occupation,energy,ehomo,elumo)

   en%kin  = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
   en%nuc  = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix(:,:,:),DIM=3) )

   !
   ! Setup kinetic and nucleus contributions (that are independent of the
   ! density matrix and therefore of spin polarization)
   !
   hamiltonian(:,:,1)        = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:) 
   hamiltonian(:,:,nspin)    = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:) 

   !
   ! Reset XC part of the Hamiltonian
   hamiltonian_xc(:,:,:) = 0.0_dp

   if( calc_type%read_potential ) then
     call read_potential(print_volume,basis%nbf,nspin,p_matrix,matrix,en%hart)
   else
     !
     ! Hartree contribution to the Hamiltonian
     !
     call setup_hartree(print_volume,basis%nbf,nspin,p_matrix,matrix,en%hart)
   endif

   hamiltonian(:,:,:)    = hamiltonian(:,:,:) + matrix(:,:,:)
  
   !
   ! Exchange contribution to the Hamiltonian
   if( calc_type%need_exchange ) then

     call setup_exchange(print_volume,basis%nbf,nspin,p_matrix,matrix,energy_tmp)
     en%exx = alpha_hybrid * energy_tmp
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix(:,:,:) * alpha_hybrid

     if(calc_type%is_screened_hybrid) then
       call setup_exchange_longrange(print_volume,basis%nbf,nspin,p_matrix,matrix,energy_tmp)
       en%exx = en%exx + alpha_hybrid_lr * energy_tmp
       hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix(:,:,:) * alpha_hybrid_lr
     endif


   endif

   !
   ! DFT XC potential is added here
   if( ndft_xc /= 0 ) then
     call start_clock(timing_dft)
     call dft_exc_vxc(nspin,basis,ndft_xc,dft_xc_type,dft_xc_coef,p_matrix,ehomo,vxc_matrix,en%xc)
     call stop_clock(timing_dft)

     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + vxc_matrix(:,:,:)

     title='=== DFT XC contribution ==='
     call dump_out_matrix(print_volume,title,basis%nbf,nspin,vxc_matrix)
   endif

   !
   ! QPscGW self energy
   if( calc_type%is_gw .AND. calc_type%method == QS .AND. iscf > 5 ) then

     call init_spectral_function(basis%nbf,prod_basis%nbf,nspin,occupation,wpol)
     call start_clock(timing_pola)
#ifdef AUXIL_BASIS
     call polarizability_rpa(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,en%rpa,wpol)
#else
     call polarizability_rpa_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,en%rpa,wpol)
#endif
     call stop_clock(timing_pola)
     en%tot = en%tot + en%rpa
     WRITE_MASTER(*,'(/,a,f16.10)') ' RPA Total energy [Ha]: ',en%tot

     call start_clock(timing_self)
     exchange_m_vxc_diag(:,:)=0.0_dp
#ifdef AUXIL_BASIS
     call gw_selfenergy(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#else
     call gw_selfenergy_noaux(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#endif
     call stop_clock(timing_self)

     matrix = alpha_mixing * matrix + (1.0_dp-alpha_mixing) * self_energy_old
     self_energy_old = matrix
     title='=== Self-energy ==='
     call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)
     call destroy_spectral_function(wpol)

     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix(:,:,:)

   endif

   !
   ! QPscMP2
   if( calc_type%is_mp2 .AND. calc_type%method == QS .AND. iscf > 5 ) then

!     call start_clock(timing_mp2_energy)
!     call mp2_energy(nspin,basis,occupation,c_matrix,energy,en%mp2)
!     call stop_clock(timing_mp2_energy)

     exchange_m_vxc_diag(:,:)=0.0_dp
     call start_clock(timing_mp2_self)
     call mp2_selfenergy(calc_type%method,nspin,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,matrix,en%mp2)
     call stop_clock(timing_mp2_self)
     WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Energy       [Ha]:',en%mp2
     WRITE_MASTER(*,*) 
     en%tot = en%tot + en%mp2
     WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Total Energy [Ha]:',en%tot

     matrix = alpha_mixing * matrix + (1.0_dp-alpha_mixing) * self_energy_old
     self_energy_old = matrix
     title='=== Self-energy ==='
     call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)
  
     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix(:,:,:)

   endif

   !
   ! Add the XC part of the hamiltonian to the total hamiltonian
   hamiltonian(:,:,:) = hamiltonian(:,:,:) + hamiltonian_xc(:,:,:)
   
   title='=== Total Hamiltonian ==='
   call dump_out_matrix(print_volume,title,basis%nbf,nspin,hamiltonian)
  
   !
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \phi = E S \phi
   ! save the old eigenvalues
   do ispin=1,nspin
     WRITE_MASTER(*,*) 'Diagonalization for spin polarization',ispin
     call diagonalize_generalized_sym(basis%nbf,&
                                      hamiltonian(:,:,ispin),s_matrix(:,:),&
                                      energy(:,ispin),c_matrix(:,:,ispin))
   enddo
  
   title='=== Energies ==='
   call dump_out_eigenenergy(title,basis%nbf,nspin,occupation,energy)
  
   !
   ! REMEMBER:
   ! \phi_a = \sum_i C_{ia} \varphi_i
   ! 
   ! hence transpose the c_matrix for a correct output by dump_out_matrix
   do ispin=1,nspin
     matrix(:,:,ispin) = transpose( c_matrix(:,:,ispin) )
   enddo
   title='=== C coefficients ==='
   call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)
!   matrix(:,:,1) = matmul( c_matrix(:,:,1), matmul( s_matrix(:,:), transpose(c_matrix(:,:,1)) ) )
!   title='=== C S C^T = identity ? ==='
!   call dump_out_matrix(print_volume,title,basis%nbf,1,matrix)
   matrix(:,:,1) = matmul( transpose(c_matrix(:,:,1)), matmul( s_matrix(:,:), c_matrix(:,:,1) ) )
   title='=== C^T S C = identity ? ==='
   call dump_out_matrix(print_volume,title,basis%nbf,1,matrix)

  
   !
   ! Setup the new density matrix: p_matrix
   ! Save the old one for the convergence criterium
   p_matrix_old(:,:,:) = p_matrix(:,:,:)
   call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
   title='=== density matrix P ==='
   call dump_out_matrix(print_volume,title,basis%nbf,nspin,p_matrix)
  
   !
   ! Output the total energy and its components
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Nucleus-Nucleus [Ha]:',en%nuc_nuc
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Kinetic Energy  [Ha]:',en%kin
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Nucleus Energy  [Ha]:',en%nuc
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Hartree Energy  [Ha]:',en%hart
   if(calc_type%need_exchange) then
     WRITE_MASTER(*,'(a25,x,f16.10)') 'Exchange Energy [Ha]:',en%exx
   endif
   if( ndft_xc /= 0 ) then
     WRITE_MASTER(*,'(a25,x,f16.10)') 'XC Energy       [Ha]:',en%xc
   endif
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%xc
   WRITE_MASTER(*,'(/,a25,x,f16.10,/)') 'Total Energy    [Ha]:',en%tot

   !
   ! Store the history of residuals
   call store_residual(p_matrix_old,p_matrix)
   !
   ! Then check the convergence
   call check_convergence(scf_loop_convergence)
   !
   ! Produce the next density matrix
   call new_p_matrix(p_matrix)
   if(scf_loop_convergence) exit
   
 !
 ! end of the big SCF loop
 enddo

 call destroy_scf()


 WRITE_MASTER(*,*) '=================================================='
 WRITE_MASTER(*,*) 'The SCF loop ends here'
 WRITE_MASTER(*,*) '=================================================='
 WRITE_MASTER(*,'(/,/,a25,x,f16.10,/,/)') 'SCF Total Energy [Ha]:',en%tot

 if(MODULO(print_volume/1000 ,2)>0) call write_density_matrix(nspin,basis%nbf,p_matrix)
 if(MODULO(print_volume/10000,2)>0) call plot_wfn(nspin,basis,c_matrix)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TESTING SECTION TO BE REMOVED IN THE FUTURE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
 WRITE_MASTER(*,*) '==================== TESTS ==================='
 allocate(matrix3(basis%nbf,basis%nbf,nspin))

 do iatom=1,3
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,*) 'CHECK TRK sum-rule along axis',iatom 
   WRITE_MASTER(*,*)

   do jbf=1,basis%nbf
     do ibf=1,basis%nbf
       call basis_function_dipole(basis%bf(ibf),basis%bf(jbf),dipole)
       matrix(ibf,jbf,1) = dipole(iatom)
     enddo
   enddo
   title='=== Dipole matrix === in gaussian basis'
   call dump_out_matrix(print_volume,title,basis%nbf,1,matrix(:,:,1))
  
   matrix(:,:,1) = MATMUL( TRANSPOSE( c_matrix(:,:,1) ) , MATMUL( matrix(:,:,1) , c_matrix(:,:,1) ) )

   title='=== Dipole matrix === in orbital basis'
   call dump_out_matrix(print_volume,title,basis%nbf,1,matrix(:,:,1))

   do jbf=1,basis%nbf
     do ibf=1,basis%nbf
       call basis_function_dipole_sq(basis%bf(ibf),basis%bf(jbf),dipole)
       matrix3(ibf,jbf,1) = dipole(iatom)
     enddo
   enddo
   title='=== Dipole squared matrix === in gaussian basis'
   call dump_out_matrix(print_volume,title,basis%nbf,1,matrix3(:,:,1))
  
   matrix3(:,:,1) = MATMUL( TRANSPOSE( c_matrix(:,:,1) ) , MATMUL( matrix3(:,:,1) , c_matrix(:,:,1) ) )

   title='=== Dipole squared matrix === in orbital basis'
   call dump_out_matrix(print_volume,title,basis%nbf,1,matrix3(:,:,1))
   

   
   do istate=1,MIN(basis%nbf,2)
     WRITE_MASTER(*,*) '______ state _____',istate
     energy_tmp=0.0_dp
     do jstate=1,basis%nbf
       energy_tmp = energy_tmp + matrix(istate,jstate,1)**2
     enddo
     WRITE_MASTER(*,*) 'test completeness \sum_j <i|r|j><j|r|i>'
     WRITE_MASTER(*,*) 'result 1:',istate,energy_tmp
     WRITE_MASTER(*,*) 'test completeness <i|r^2|i>'
     WRITE_MASTER(*,*) 'result 2:',istate,matrix3(istate,istate,1)
     WRITE_MASTER(*,*)
   enddo

   rtmp = -0.2
   WRITE_MASTER(*,*) 'high energy [Ha]',rtmp


   do jstate=1,MIN(basis%nbf,2)

     energy_tmp=0.0_dp
     do istate=1,basis%nbf
       energy_tmp = energy_tmp + ( energy(istate,1) - energy(jstate,1) ) * matrix(istate,jstate,1)**2
     enddo
     WRITE_MASTER(*,*) 'TRK result:',jstate,energy_tmp
!     energy_tmp=0.0_dp

     do istate=1,basis%nbf
       energy_tmp = energy_tmp - ( rtmp             - energy(jstate,1) ) * matrix(istate,jstate,1)**2
     enddo
     WRITE_MASTER(*,*) 'TRK result:',jstate,energy_tmp
!     energy_tmp=0.0_dp

     energy_tmp = energy_tmp + ( rtmp             - energy(jstate,1) )  * matrix3(jstate,jstate,1)
     WRITE_MASTER(*,*) 'TRK result:',jstate,energy_tmp

   enddo


   energy_tmp=0.0_dp
   do jstate=1,basis%nbf
     do istate=1,basis%nbf
       if( occupation(jstate,1) - occupation(istate,1)  < 1.0e-6_dp ) cycle
!       WRITE_MASTER(*,*) '----------',istate,jstate
!       WRITE_MASTER(*,*) occupation(jstate,1) - occupation(istate,1),energy(istate,1) - energy(jstate,1), matrix(istate,jstate,1)**2 
       energy_tmp = energy_tmp + matrix(istate,jstate,1)**2 / ( energy(istate,1) - energy(jstate,1) ) * ( occupation(jstate,1) - occupation(istate,1) )
     enddo
   enddo
   !
   ! factor 2 from resonant + anti resonant
   WRITE_MASTER(*,*) 'polarizability',2.0_dp*energy_tmp


   WRITE_MASTER(*,*) '---------with completeness'

!   energy_tmp=0.0_dp
   do jstate=1,basis%nbf
     do istate=1,basis%nbf
       energy_tmp = energy_tmp - matrix(istate,jstate,1)**2 / ( rtmp             - energy(jstate,1) ) *  occupation(jstate,1)
     enddo
   enddo
   !
   ! factor 2 from resonant + anti resonant
   WRITE_MASTER(*,*) 'polarizability without delta',2.0_dp*energy_tmp

!   energy_tmp=0.0_dp
   do jstate=1,basis%nbf
     energy_tmp = energy_tmp + occupation(jstate,1) * matrix3(jstate,jstate,1) /  ( rtmp             - energy(jstate,1) )  
   enddo
   WRITE_MASTER(*,*) 'polarizability with    delta',2.0_dp*energy_tmp


 enddo
 WRITE_MASTER(*,*) '========= END   OF   TESTS ==================='
 stop'ENOUGH'
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END OF SECTION TO BE REMOVED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !
 ! CI calculation is done here
 ! implemented for 2 electrons only!
 if(calc_type%is_ci) then
   if(nspin/=1) stop'for CI, nspin should be 1'
   if( ABS( electrons - 2.0_dp ) > 1.e-5_dp ) stop'CI is implemented for 2 electrons only'
   call full_ci_2electrons_spin(print_volume,0,basis,hamiltonian_kinetic+hamiltonian_nucleus,c_matrix,en%nuc_nuc)
 endif
  

 !
 ! Single excitation term
 if( .TRUE.) then

   !
   ! Obtain the Fock matrix
   call setup_exchange(print_volume,basis%nbf,nspin,p_matrix,matrix,en%exx)
   matrix(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + matrix(:,:,:)

   !
   ! Rotate the Fock matrix to the eigenstate basis
   call matrix_basis_to_eigen(nspin,basis%nbf,c_matrix,matrix)

   title='=== Fock matrix ==='
   call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)

   spin_fact = REAL(-nspin+3,dp)
   en%se = 0.0_dp
   do ispin=1,nspin
     ! loop on occupied states
     do istate=1,basis%nbf
       if( occupation(istate,ispin) < completely_empty ) cycle
       ! loop on virtual states
       do astate=1,basis%nbf
         if( occupation(astate,ispin) > spin_fact - completely_empty ) cycle
         en%se = en%se + matrix(istate,astate,ispin)**2 / ( energy(istate,ispin) - energy(astate,ispin) ) * spin_fact
       enddo
     enddo
   enddo

   WRITE_MASTER(*,'(a,2x,f16.10)') ' Etotal EXX       [Ha]:',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx 
   WRITE_MASTER(*,'(a,2x,f16.10)') ' Single-Excitation[Ha]:',en%se

 endif



 !
 ! in case of DFT + GW
 if( calc_type%need_final_exchange ) then

   inquire(file='manual_coresplitting',exist=file_exists)
   if(file_exists) then
     WRITE_MASTER(*,*) 'TESTING CORE-VALENCE SPLITTING'
     open(13,file='manual_coresplitting')
     read(13,*) ncore
     close(13)
     WRITE_MASTER(msg,'(a,i4,2x,i4)') 'core-valence splitting switched on up to state = ',ncore
     call issue_warning(msg)
     do istate=1,ncore
       occupation(istate,:) = 0.0_dp
     enddo
     call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
     call dft_exc_vxc(nspin,basis,ndft_xc,dft_xc_type,dft_xc_coef,p_matrix,ehomo,vxc_matrix,en%xc)
     hamiltonian_xc(:,:,:) = vxc_matrix
   endif


   call setup_exchange(print_volume,basis%nbf,nspin,p_matrix,matrix,en%exx)
   WRITE_MASTER(*,*) 'EXX [Ha]:',en%exx

   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * ( matrix(ibf,jbf,ispin) - hamiltonian_xc(ibf,jbf,ispin) )&
                    * c_matrix(jbf,istate,ispin)
         enddo
       enddo
     enddo
   enddo

 else
   exchange_m_vxc_diag(:,:) = 0.0_dp
 endif

 !
 ! final evaluation for G0W0
 if( calc_type%is_gw .AND. calc_type%method == perturbative ) then

   if( calc_type%is_lr_mbpt ) then
     !
     ! Hard-coded cutoff radius
     rcut= 9.0909_dp
     write(msg,'(a,f10.4)') 'hard coded cutoff radius  ',rcut
     call issue_warning(msg)
     call deallocate_eri()
     call allocate_eri(basis%nbf)
     call calculate_eri(print_volume,basis,rcut,BUFFER1)

     if( .NOT. ALLOCATED( vxc_matrix ) ) allocate( vxc_matrix(basis%nbf,basis%nbf,nspin) )
     call dft_exc_vxc(nspin,basis,1,(/1000/),(/1.0_dp/),p_matrix,ehomo,vxc_matrix,energy_tmp)
     WRITE_MASTER(*,'(/,a,f16.10)') '    RPA LDA energy [Ha]: ',energy_tmp
     exchange_m_vxc_diag(:,:) = 0.0_dp
     do ispin=1,nspin
       do istate=1,basis%nbf
         do ibf=1,basis%nbf
           do jbf=1,basis%nbf
             exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                     + c_matrix(ibf,istate,ispin) * ( vxc_matrix(ibf,jbf,ispin) )&
                      * c_matrix(jbf,istate,ispin)
           enddo
         enddo
       enddo
     enddo
     WRITE_MASTER(*,*) 'test'
     WRITE_MASTER(*,*) '<vxc RPA FR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
     WRITE_MASTER(*,*) 'test'

     vxc_matrix(:,:,:) = 0.0_dp
     call dft_exc_vxc(nspin,basis,1,(/1005/),(/1.0_dp/),p_matrix,ehomo,vxc_matrix,energy_tmp)
     WRITE_MASTER(*,'(/,a,f16.10)') ' LR-RPA LDA energy [Ha]: ',energy_tmp
     exchange_m_vxc_diag(:,:) = 0.0_dp
     do ispin=1,nspin
       do istate=1,basis%nbf
         do ibf=1,basis%nbf
           do jbf=1,basis%nbf
             exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                     + c_matrix(ibf,istate,ispin) * ( vxc_matrix(ibf,jbf,ispin) )&
                      * c_matrix(jbf,istate,ispin)
           enddo
         enddo
       enddo
     enddo
     WRITE_MASTER(*,*) 'test'
     WRITE_MASTER(*,*) '<vxc RPA LR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
     WRITE_MASTER(*,*) 'test'

     call setup_exchange(print_volume,basis%nbf,nspin,p_matrix,matrix,en%exx)  
     WRITE_MASTER(*,*) 'LR EXX [Ha]:',en%exx
     exchange_m_vxc_diag(:,:) = 0.0_dp
     do ispin=1,nspin
       do istate=1,basis%nbf
         do ibf=1,basis%nbf
           do jbf=1,basis%nbf
             exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                     + c_matrix(ibf,istate,ispin) * matrix(ibf,jbf,ispin)&
                      * c_matrix(jbf,istate,ispin)
           enddo
         enddo
       enddo
     enddo
     WRITE_MASTER(*,*) 'test'
     WRITE_MASTER(*,*) '<sigx LR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
     WRITE_MASTER(*,*) 'test'

     vxc_matrix(:,:,:) = 0.0_dp
     call dft_exc_vxc(nspin,basis,1,(/2001/),(/1.0_dp/),p_matrix,ehomo,vxc_matrix,energy_tmp)
     WRITE_MASTER(*,'(/,a,f16.10)') ' SR-EXX GGA energy [Ha]: ',energy_tmp
     exchange_m_vxc_diag(:,:) = 0.0_dp
     do ispin=1,nspin
       do istate=1,basis%nbf
         do ibf=1,basis%nbf
           do jbf=1,basis%nbf
             exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                     + c_matrix(ibf,istate,ispin) * ( vxc_matrix(ibf,jbf,ispin))&
                      * c_matrix(jbf,istate,ispin)
           enddo
         enddo
       enddo
     enddo
     WRITE_MASTER(*,*) 'test'
     WRITE_MASTER(*,*) '<vxc EXX SR>', exchange_m_vxc_diag(1:MIN(basis%nbf,15),:)*27.211
     WRITE_MASTER(*,*) 'test'



     exchange_m_vxc_diag(:,:) = 0.0_dp

   endif

   call init_spectral_function(basis%nbf,prod_basis%nbf,nspin,occupation,wpol)
   call start_clock(timing_pola)
#ifdef AUXIL_BASIS
   call polarizability_rpa(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,en%rpa,wpol)
#else
#ifdef CASIDA
   call polarizability_casida_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,en%rpa,wpol)
#else
   call polarizability_rpa_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,en%rpa,wpol)
#endif
#endif
   call stop_clock(timing_pola)
   en%tot = en%tot + en%rpa
   if( ndft_xc /= 0 ) en%tot = en%tot - en%xc + en%exx * ( 1.0_dp - alpha_hybrid )
   WRITE_MASTER(*,'(/,a,f16.10)') ' RPA Total energy [Ha]: ',en%tot

!   call write_spectral_function(wpol)
!   msg='write spectral function'
!   call read_spectral_function(wpol)
!   msg='read spectral function'
!   call issue_warning(msg)

   call start_clock(timing_self)
#ifndef CASIDA
#ifdef AUXIL_BASIS
   call gw_selfenergy(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#else
   call gw_selfenergy_noaux(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#endif
#endif
   call stop_clock(timing_self)

   title='=== Self-energy === (in the orbital basis)'
   call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)
   call destroy_spectral_function(wpol)

   if(allocated(sinv_v_sinv)) deallocate(sinv_v_sinv)
 endif ! G0W0

 !
 ! final evaluation for MP2
 if( calc_type%is_mp2 .AND. calc_type%method == perturbative ) then
!   call start_clock(timing_mp2_energy)
!   call mp2_energy(nspin,basis,occupation,c_matrix,energy,en%mp2)
!   call stop_clock(timing_mp2_energy)

!   call setup_hartree(print_volume,basis%nbf,nspin,p_matrix,matrix,en%hart)
!   WRITE_MASTER(*,*) 'Hartree [Ha]:',en%hart
   call setup_exchange(print_volume,basis%nbf,nspin,p_matrix,matrix,en%exx)
   WRITE_MASTER(*,*) 'EXX     [Ha]:',en%exx

   call start_clock(timing_mp2_self)
   call mp2_selfenergy(calc_type%method,nspin,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,matrix,en%mp2)
   call stop_clock(timing_mp2_self)
   WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Energy       [Ha]:',en%mp2
   WRITE_MASTER(*,*) 
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%mp2

   WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Total Energy [Ha]:',en%tot

   title='=== Self-energy === (in the orbital basis)'
   call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)

 endif



 !
 ! Cleanly exiting the code
 !
 deallocate(hamiltonian,hamiltonian_xc,hamiltonian_kinetic,hamiltonian_nucleus)
 deallocate(matrix,s_matrix,c_matrix,p_matrix,p_matrix_old)
 deallocate(energy,occupation,exchange_m_vxc_diag)
 deallocate(self_energy_old)
 call deallocate_eri()
 call deallocate_eri_lr()
 if( ndft_xc /= 0 ) deallocate( vxc_matrix )

 call destroy_basis_set(basis)
 if(calc_type%is_gw) call destroy_basis_set(prod_basis)

 call stop_clock(timing_total)
 call output_timing()

 call output_all_warnings()

 WRITE_MASTER(*,'(/,a,/)') ' This is the end'

 call finish_mpi()

end program molgw
!=========================================================================
