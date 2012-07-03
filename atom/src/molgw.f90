!=========================================================================
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

 integer,parameter            :: gaussian_type=CARTESIAN
 !
 ! Input parameters will be set in read_inputparameters
 type(calculation_type)       :: calc_type
 integer                      :: nspin,nscf
 real(dp)                     :: alpha_mixing
 integer                      :: print_volume
 character(len=100)           :: basis_name
 real(dp)                     :: electrons
 real(dp)                     :: magnetization
!===== variables for testing
 type(gaussian) :: gatmp,gbtmp
 type(basis_function) :: bftmp1,bftmp2
 real(dp) :: rtmp
!=====
 type(basis_set)         :: basis
 type(basis_set)         :: basis_cart
 type(basis_set)         :: prod_basis
 type(spectral_function) :: wpol
 integer                 :: ibf,jbf,kbf,lbf,ijbf,klbf
 integer                 :: ispin,iscf,istate,jstate,iatom
 logical                 :: scf_loop_convergence
 logical                 :: file_exists
 character(len=100)      :: title
 real(dp)                :: energy_tmp,overlap_tmp
 real(dp)                :: dipole(3)
 real(dp),allocatable    :: ehomo(:),elumo(:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_xc(:,:,:)
 real(dp),allocatable    :: hamiltonian_kinetic(:,:,:)           !TODO remove spin
 real(dp),allocatable    :: hamiltonian_nucleus(:,:,:)           !TODO remove spin
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
   real(dp) :: nuc_nuc=0.0_dp
   real(dp) :: kin    =0.0_dp
   real(dp) :: nuc    =0.0_dp
   real(dp) :: hart   =0.0_dp
   real(dp) :: exx    =0.0_dp
   real(dp) :: xc     =0.0_dp
   real(dp) :: mp2    =0.0_dp
   real(dp) :: rpa    =0.0_dp
   real(dp) :: tot    =0.0_dp
 end type
!=====
 type(energy_contributions) :: en
!=============================

 call init_mpi()

 call header()


 !
 ! Development tests are commented below
#if 0

 call init_gaussian_general(1,2,3,0.4361_dp,(/0.0_dp,0.0_dp,-1.0_dp/),gatmp)
 call print_gaussian(gatmp)

 call init_gaussian_general(1,4,1,0.8120_dp,(/0.0_dp,0.0_dp,1.2_dp/),gbtmp)
 call print_gaussian(gbtmp)

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) ' === CHECK OVERLAP === '

! call overlap_normalized(gatmp,gbtmp,rtmp)
! WRITE_MASTER(*,*) 'normalized S_ab',rtmp

 call overlap_recurrence(gatmp,gbtmp,rtmp)
 WRITE_MASTER(*,*) 'normalized S_ab from recurrence',rtmp
 call overlap_recurrence(gbtmp,gatmp,rtmp)
 WRITE_MASTER(*,*) 'normalized S_ba from recurrence',rtmp

! call numerical_overlap(gatmp,gbtmp)


 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) ' === CHECK KINETIC === '

! call kinetic_gaussian(gatmp,gbtmp,rtmp)
! WRITE_MASTER(*,*) 'kinetic matrix element [Ha]',rtmp

 call kinetic_recurrence(gatmp,gbtmp,rtmp)
 WRITE_MASTER(*,*) 'new kinetic matrix element K_ab',rtmp
 call kinetic_recurrence(gbtmp,gatmp,rtmp)
 WRITE_MASTER(*,*) 'new kinetic matrix element K_ba',rtmp

! call numerical_kinetic(gatmp,gbtmp)

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) ' === CHECK NUCLEUS === '

! call nucleus_pot_gaussian(gatmp,gbtmp,1.0_dp,rtmp)
! WRITE_MASTER(*,*) 'nucleus pot [Ha]',rtmp

 call nucleus_recurrence(1.0_dp,(/0.0_dp,0.0_dp,0.0_dp/),gatmp,gbtmp,rtmp)
 WRITE_MASTER(*,*) 'new nucleus matrix element V_ba',rtmp
 call nucleus_recurrence(1.0_dp,(/0.0_dp,0.0_dp,0.0_dp/),gbtmp,gatmp,rtmp)
 WRITE_MASTER(*,*) 'new nucleus matrix element V_ba',rtmp

 call numerical_nucleus(gatmp,gbtmp)


 stop'ENOUGH FOR TODAY'

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) '                   END OF THE TESTS'
 WRITE_MASTER(*,*) '==========================================================='
 WRITE_MASTER(*,*) '==========================================================='
 WRITE_MASTER(*,*) '==========================================================='
 WRITE_MASTER(*,*) '                   START REAL LIFE CALCULATION'
 WRITE_MASTER(*,*)
#endif 

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
 ! reading input file
 call read_inputparameter_molecule(calc_type,nspin,nscf,alpha_mixing,print_volume,&
                                   basis_name,electrons,magnetization)

 !
 ! Nucleus-nucleus repulsion contribution to the energy
 call nucleus_nucleus_energy(en%nuc_nuc)

 !
 ! Build up the basis set
 !
 select case(gaussian_type)
 case(CARTESIAN)
   call init_basis_set(print_volume,basis_name,CARTESIAN,basis)
 case(PURE)
   call init_basis_set(print_volume,basis_name,CARTESIAN,basis_cart)
   call init_basis_set(print_volume,basis_name,PURE     ,basis)
   call set_cart_to_pure_transform(basis_cart,basis)
 end select

 !
 ! first attempt to distribute the work load among procs
 call init_distribution(basis%nbf)
 
 !
 ! allocate everything
 allocate(hamiltonian(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_xc(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_kinetic(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_nucleus(basis%nbf,basis%nbf,nspin))
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
 do jbf=1,basis%nbf
   do ibf=1,jbf ! basis%nbf ! use the symmetry of the matrix
     call overlap_basis_function(basis%bf(ibf),basis%bf(jbf),overlap_tmp)
     s_matrix(ibf,jbf) = overlap_tmp
     s_matrix(jbf,ibf) = overlap_tmp
   enddo
   if( ABS( s_matrix(jbf,jbf) - 1.0_dp ) > 1.0d-4 ) then
     WRITE_MASTER(*,*) 'A diagonal term of the overlap matrix is not equal to 1.0'
     WRITE_MASTER(*,*) jbf,jbf,s_matrix(jbf,jbf)
     WRITE_MASTER(*,*) 'check the basis set definition'
     stop'ERROR'
   endif
 enddo
 title='=== overlap matrix S ==='
 call dump_out_matrix(print_volume,title,basis%nbf,1,s_matrix)

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
 ! in case of GW run, set up product basis 
 !                    set up coulomb integrals
 if( calc_type%is_gw ) then
   call start_clock(timing_prodbasis)

   call init_product_basis_set(basis,prod_basis)

#ifdef AUXIL_BASIS
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
#endif
   
   call stop_clock(timing_prodbasis)
 endif

 !
 ! build occupation array, initial c_matrix
 call set_occupation(electrons,magnetization,basis%nbf,nspin,occupation)
 title='=== occupations ==='
 call dump_out_array(.FALSE.,title,basis%nbf,nspin,occupation)
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
 ! Possibility offered to override the automatic setup of the initial density
 ! matrix
 inquire(file='p_matrix_diag.in',exist=file_exists)
 if(file_exists) then
   WRITE_MASTER(*,*) 'reading input density matrix from file'
   open(unit=11,file='p_matrix_diag.in',status='old')
   p_matrix(:,:,:) = 0.0_dp
   do ispin=1,nspin
     do ibf=1,basis%nbf
       read(11,*) p_matrix(ibf,ibf,ispin) 
     enddo
   enddo
   close(11)
   if( ABS( SUM(p_matrix(:,:,:)) - electrons ) > 1.d-8 ) &
     stop'input density matrix does not contain the right number of electrons'
   msg='manual input of the initial density matrix diagonal'
   call issue_warning(msg)
 endif
 title='=== 1st density matrix P ==='
 call dump_out_matrix(print_volume,title,basis%nbf,nspin,p_matrix)
! matrix(:,:,1) = matmul( p_matrix(:,:,1) ,p_matrix(:,:,1) )
! matrix(:,:,nspin) = matmul( p_matrix(:,:,nspin) ,p_matrix(:,:,nspin) )
! title='=== 1st density matrix P ==='
! call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)

 !
 ! Initialize the SCF mixing procedure
 call init_scf(nscf,basis%nbf,nspin,alpha_mixing)


 !
 ! Calculate the parts of the hamiltonian that does not change along
 ! with the SCF cycles
 !
 ! Kinetic energy contribution
 do jbf=1,basis%nbf
   do ibf=1,basis%nbf
     call kinetic_basis_function(basis%bf(ibf),basis%bf(jbf),energy_tmp)
     hamiltonian_kinetic(ibf,jbf,:) = energy_tmp
   enddo
 enddo
 title='=== kinetic energy contribution ==='
 call dump_out_matrix(print_volume,title,basis%nbf,nspin,hamiltonian_kinetic)

 !
 ! nucleus-electron interaction
 hamiltonian_nucleus(:,:,:) =  0.0_dp
 do iatom=1,natom
   do ibf=1,basis%nbf
     do jbf=1,basis%nbf
       call nucleus_pot_basis_function(basis%bf(ibf),basis%bf(jbf),zatom(iatom),x(:,iatom),energy_tmp)
       hamiltonian_nucleus(ibf,jbf,:) =  hamiltonian_nucleus(ibf,jbf,:) + energy_tmp
     enddo
   enddo
 enddo
 title='=== nucleus contribution ==='
 call dump_out_matrix(print_volume,title,basis%nbf,nspin,hamiltonian_nucleus)


 !
 ! start the big scf loop
 !
 do iscf=1,nscf
   WRITE_MASTER(*,'(/,a)') '-------------------------------------------'
   WRITE_MASTER(*,'(a,x,i4,/)') ' *** SCF cycle No:',iscf

   call output_homolumo(basis%nbf,nspin,occupation,energy,ehomo,elumo)

   !
   ! Skip the first iteration
   if(iscf>1) call new_p_matrix(p_matrix)
   

   en%kin  = SUM(hamiltonian_kinetic(:,:,:)*p_matrix(:,:,:))
   en%nuc  = SUM(hamiltonian_nucleus(:,:,:)*p_matrix(:,:,:))

   !
   ! Setup kinetic and nucleus contributions (that are independant of the
   ! density matrix)
   !
   hamiltonian(:,:,:)    = hamiltonian_kinetic(:,:,:) + hamiltonian_nucleus(:,:,:) 

   !
   ! Reset XC part of the Hamiltonian
   hamiltonian_xc(:,:,:) = 0.0_dp

   !
   ! Hartree contribution to the Hamiltonian
   !
   call setup_hartree(print_volume,basis%nbf,nspin,p_matrix,matrix,en%hart)

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
     call polarizability_casida(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,en%rpa,wpol)
#else
     call polarizability_casida_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,en%rpa,wpol)
#endif
     call stop_clock(timing_pola)
     en%tot = en%tot + en%rpa
     WRITE_MASTER(*,'(/,a,f14.8)') ' RPA Total energy [Ha]: ',en%tot

     call start_clock(timing_self)
     exchange_m_vxc_diag(:,:)=0.0_dp
#ifdef AUXIL_BASIS
     call gw_selfenergy_casida(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#else
     call gw_selfenergy_casida_noaux(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
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
     WRITE_MASTER(*,'(a,2x,f14.8)') ' MP2 Energy       [Ha]:',en%mp2
     WRITE_MASTER(*,*) 
     en%tot = en%tot + en%mp2
     WRITE_MASTER(*,'(a,2x,f14.8)') ' MP2 Total Energy [Ha]:',en%tot

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
   call dump_out_array(.TRUE.,title,basis%nbf,nspin,energy)
  
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
   if(scf_loop_convergence) exit
   
 !
 ! end of the big SCF loop
 enddo

 call destroy_scf()


 WRITE_MASTER(*,*) '=================================================='
 WRITE_MASTER(*,*) 'The SCF loop ends here'
 WRITE_MASTER(*,*) '=================================================='
 WRITE_MASTER(*,*)

 if(MODULO(print_volume,100)>4) call plot_wfn(nspin,basis,c_matrix)

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
   call full_ci_2electrons_spin(0,basis,hamiltonian_kinetic+hamiltonian_nucleus,c_matrix,en%nuc_nuc)
 endif
  
 !
 ! in case of DFT + GW
 if( calc_type%need_final_exchange ) then

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

!  WRITE_MASTER(*,*) 'test'
!  WRITE_MASTER(*,*) '<Sigx>', exchange_m_vxc_diag(1:2,:)*27.211
!  WRITE_MASTER(*,*) 'test'

 else
   exchange_m_vxc_diag(:,:) = 0.0_dp
 endif

 !
 ! final evaluation for G0W0
 if( calc_type%is_gw .AND. calc_type%method == perturbative ) then

   if( calc_type%is_lr_mbpt ) then
     !
     ! Hard-coded cutoff radius
     rcut= 0.5_dp
     write(msg,'(a,f10.4)') 'hard coded cutoff radius  ',rcut
     call issue_warning(msg)
     call deallocate_eri()
     call allocate_eri(basis%nbf)
     call calculate_eri(print_volume,basis,rcut,BUFFER1)

     if( .NOT. ALLOCATED( vxc_matrix ) ) allocate( vxc_matrix(basis%nbf,basis%nbf,nspin) )
     call dft_exc_vxc(nspin,basis,1,(/1000/),(/1.0_dp/),p_matrix,ehomo,vxc_matrix,energy_tmp)
     WRITE_MASTER(*,'(/,a,f14.8)') '    RPA LDA energy [Ha]: ',energy_tmp
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
     WRITE_MASTER(*,'(/,a,f14.8)') ' LR-RPA LDA energy [Ha]: ',energy_tmp
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
     exchange_m_vxc_diag(:,:) = 0.0_dp

   endif

   call init_spectral_function(basis%nbf,prod_basis%nbf,nspin,occupation,wpol)
   call start_clock(timing_pola)
#ifdef AUXIL_BASIS
   call polarizability_casida(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,en%rpa,wpol)
#else
   call polarizability_casida_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,en%rpa,wpol)
#endif
   call stop_clock(timing_pola)
   en%tot = en%tot + en%rpa
   if( ndft_xc /= 0 ) en%tot = en%tot - en%xc + en%exx * ( 1.0_dp - alpha_hybrid )
   WRITE_MASTER(*,'(/,a,f14.8)') ' RPA Total energy [Ha]: ',en%tot

   call start_clock(timing_self)
#ifdef AUXIL_BASIS
   call gw_selfenergy_casida(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#else
   call gw_selfenergy_casida_noaux(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
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
   call start_clock(timing_mp2_self)
   call mp2_selfenergy(calc_type%method,nspin,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,matrix,en%mp2)
   call stop_clock(timing_mp2_self)
   WRITE_MASTER(*,'(a,2x,f12.6)') ' MP2 Energy       [Ha]:',en%mp2
   WRITE_MASTER(*,*) 
   en%tot = en%tot + en%mp2
   if( ndft_xc /= 0 ) en%tot = en%tot - en%xc + en%exx * ( 1.0_dp - alpha_hybrid )
   WRITE_MASTER(*,'(a,2x,f12.6)') ' MP2 Total Energy [Ha]:',en%tot

   title='=== Self-energy === (in the orbital basis)'
   call dump_out_matrix(print_volume,title,basis%nbf,nspin,matrix)

 endif

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
subroutine setup_density_matrix(nbf,nspin,c_matrix,occupation,p_matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(in)  :: occupation(nbf,nspin)
 real(dp),intent(out) :: p_matrix(nbf,nbf,nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       p_matrix(ibf,jbf,ispin) = SUM( occupation(:,ispin) * c_matrix(ibf,:,ispin) * c_matrix(jbf,:,ispin) )
     enddo
   enddo
 enddo


end subroutine setup_density_matrix

!=========================================================================
subroutine  set_occupation(electrons,magnetization,nbf,nspin,occupation)
 use m_definitions
 use m_mpi
 use m_warning
 implicit none
 real(dp),intent(in)  :: electrons,magnetization
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(out) :: occupation(nbf,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin),spin_fact
 integer              :: ibf,nlines,ilines
 logical              :: file_exists
!=====


  occupation(:,:)=0.0_dp
  spin_fact = REAL(-nspin+3,dp)

  inquire(file='manual_occupations',exist=file_exists)

  if(.NOT. file_exists) then
    remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
    if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)

    do ibf=1,nbf
      occupation(ibf,:) = MIN(remaining_electrons(:), spin_fact)
      remaining_electrons(:)  = remaining_electrons(:) - occupation(ibf,:)
    end do
  else
    WRITE_MASTER(*,*)
    WRITE_MASTER(*,*) 'occupations are read from file: manual_occupations'
    msg='reading manual occupations from file'
    call issue_warning(msg)
    open(unit=12,file='manual_occupations',status='old')
    !
    ! read nlines, all other occupations are set to zeroalpha_max_bf(jbf)
    read(12,*) nlines
    do ilines=1,nlines
      read(12,*) occupation(ilines,:)
    enddo
    close(12)
    WRITE_MASTER(*,*) 'occupations set, closing file'
  endif
 
  !
  ! final check
  if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0d-7 ) then
    WRITE_MASTER(*,*) 'occupation set up failed to give the right number of electrons'
    WRITE_MASTER(*,*) 'sum of occupations',SUM(occupation(:,:))
    WRITE_MASTER(*,*) 'electrons',electrons
    do ibf=1,nbf
      WRITE_MASTER(*,*) ibf,occupation(ibf,:)
    enddo
    stop'FAILURE in set_occupations'
  endif 

end subroutine set_occupation

!=========================================================================
subroutine guess_starting_c_matrix(nbf,nspin,c_matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(out) :: c_matrix(nbf,nbf,nspin)
!=====
 integer :: ibf
!=====

 !
 ! fill the c_matrix with the identity
 c_matrix(:,:,:)=0.0_dp
 do ibf=1,nbf
   c_matrix(ibf,ibf,:) = 1.0_dp
!   c_matrix(ibf,modulo(ibf,nbf)+1,:) = 1.0_dp
 enddo

end subroutine guess_starting_c_matrix

!=========================================================================
subroutine guess_starting_c_matrix_new(basis,nspin,c_matrix)
 use m_definitions
 use m_mpi
 use m_gaussian
 use m_basis_set
 implicit none
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nspin
 real(dp),intent(out)       :: c_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer  :: ibf,jbf,kbf,ig
 real(dp) :: alpha_max_bf(basis%nbf),alpha_max_remaining
!=====

 !
 ! find the sharpest gaussians
 alpha_max_bf(:)=0.0_dp
 do ibf=1,basis%nbf
   do ig=1,basis%bf(ibf)%ngaussian
     alpha_max_bf(ibf)=MAX(basis%bf(ibf)%g(ig)%alpha,alpha_max_bf(ibf))
   enddo
!   WRITE_MASTER(*,*) ibf,alpha_max_bf(ibf)
 enddo

 !
 ! fill the c_matrix 
 c_matrix(:,:,:)=0.0_dp
 do ibf=1,basis%nbf
   alpha_max_remaining=0.0_dp
   do jbf=1,basis%nbf
     if( alpha_max_bf(jbf) > alpha_max_remaining ) then
       alpha_max_remaining = alpha_max_bf(jbf)
       kbf = jbf
     endif
   enddo
   c_matrix(kbf,ibf,:) = 1.0_dp
!   WRITE_MASTER(*,*) 'chosen',ibf,kbf,alpha_max_bf(kbf)
   alpha_max_bf(kbf)   = -1.0_dp
   
 enddo

end subroutine guess_starting_c_matrix_new
