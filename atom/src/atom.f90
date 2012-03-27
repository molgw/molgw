!=========================================================================
program atom
 use m_definitions
 use m_timing
 use m_warning
 use m_calculation_type
 use m_tools
 use m_gaussian
 use m_basis_set
 use m_eri
 use m_gw
 use m_dft
#ifdef OPENMP
 use omp_lib
#endif
 implicit none

 !
 ! Input parameters will be set in read_inputparameters
 type(calculation_type)       :: calc_type
 integer                      :: nspin,nscf
 real(dp)                     :: alpha_mixing
 integer                      :: PRINT_VOLUME
 integer                      :: basis_element
 character(len=100)           :: basis_name
 real(dp)                     :: zatom
 real(dp)                     :: electrons
 real(dp)                     :: magnetization
!=====
 real(dp),parameter           :: alpha_hybrid=0.25_dp
!=====
 type(gaussian) :: gatmp,gbtmp
 type(basis_function) :: bftmp
 real(dp) :: rtmp,rtmp2,x(3),dx,dh,dhx(3),dhy(3),dhz(3)
 integer :: ix,iy,iz,ntmp,ng
 real(dp),allocatable :: alpha(:),coeff(:)
!=====
 type(basis_set)         :: basis
 type(basis_set)         :: prod_basis
 type(spectral_function) :: wpol
 integer                 :: ibf,jbf,kbf,lbf,ijbf,klbf
 integer                 :: ispin,iscf,istate
 integer                 :: iprodbf,jprodbf,nprodbf_max
 character(len=100)      :: title
 real(dp)                :: energy_tmp,overlap_tmp,spin_fact,rms
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_kinetic(:,:,:)           !TODO remove spin
 real(dp),allocatable    :: hamiltonian_nucleus(:,:,:)           !TODO remove spin
 real(dp),allocatable    :: matrix(:,:,:)
 real(dp),allocatable    :: vxc_matrix(:,:,:)
 real(dp),allocatable    :: s_matrix(:,:)                      !TODO remove spin
 real(dp),allocatable    :: c_matrix(:,:,:)
 real(dp),allocatable    :: p_matrix(:,:,:),p_matrix_old(:,:,:)
 real(dp),allocatable    :: energy(:,:),energy_old(:,:)
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
   real(dp) :: kin
   real(dp) :: nuc
   real(dp) :: hart
   real(dp) :: exx
   real(dp) :: xc
   real(dp) :: mp2
   real(dp) :: tot
 end type
!=====
 type(energy_contributions) :: en
 real(dp),external       :: hartree_fock_energy_band
!=============================

 write(*,*)
 write(*,*) 'Welcome to the fascinating world of ATOM'
 write(*,*)
 write(*,*)


 !
 ! Development tests are commented below
#if 0

 call init_gaussian_general(1,2,3,0.4361_dp,(/0.0_dp,0.0_dp,-1.0_dp/),gatmp)
 call print_gaussian(gatmp)

 call init_gaussian_general(1,4,1,0.8120_dp,(/0.0_dp,0.0_dp,1.2_dp/),gbtmp)
 call print_gaussian(gbtmp)

 write(*,*)
 write(*,*) ' === CHECK OVERLAP === '

! call overlap_normalized(gatmp,gbtmp,rtmp)
! write(*,*) 'normalized S_ab',rtmp

 call overlap_recurrence(gatmp,gbtmp,rtmp)
 write(*,*) 'normalized S_ab from recurrence',rtmp
 call overlap_recurrence(gbtmp,gatmp,rtmp)
 write(*,*) 'normalized S_ba from recurrence',rtmp

! call numerical_overlap(gatmp,gbtmp)


 write(*,*)
 write(*,*) ' === CHECK KINETIC === '

! call kinetic_gaussian(gatmp,gbtmp,rtmp)
! write(*,*) 'kinetic matrix element [Ha]',rtmp

 call kinetic_recurrence(gatmp,gbtmp,rtmp)
 write(*,*) 'new kinetic matrix element K_ab',rtmp
 call kinetic_recurrence(gbtmp,gatmp,rtmp)
 write(*,*) 'new kinetic matrix element K_ba',rtmp

! call numerical_kinetic(gatmp,gbtmp)

 write(*,*)
 write(*,*) ' === CHECK NUCLEUS === '

! call nucleus_pot_gaussian(gatmp,gbtmp,1.0_dp,rtmp)
! write(*,*) 'nucleus pot [Ha]',rtmp

 call nucleus_recurrence(1.0_dp,(/0.0_dp,0.0_dp,0.0_dp/),gatmp,gbtmp,rtmp)
 write(*,*) 'new nucleus matrix element V_ba',rtmp
 call nucleus_recurrence(1.0_dp,(/0.0_dp,0.0_dp,0.0_dp/),gbtmp,gatmp,rtmp)
 write(*,*) 'new nucleus matrix element V_ba',rtmp

 call numerical_nucleus(gatmp,gbtmp)

 stop'ENOUGH FOR TODAY'
#endif

#if 0
! TESTING CONTRACTION
!H STO-3G
 ng=3
 allocate(alpha(ng),coeff(ng))
 alpha(1) = 0.1098180_dp
 alpha(2) = 0.4057710_dp
 alpha(3) = 2.2276600_dp
 coeff(1) = 0.4446350_dp
 coeff(2) = 0.5353280_dp
 coeff(3) = 0.1543290_dp
 call init_basis_function(ng,0,0,0,alpha,coeff,bftmp)
 deallocate(alpha,coeff)

 call print_basis_function(bftmp)
 call overlap_basis_function(bftmp,bftmp,rtmp)
 write(*,*) 'overlap_basis_function',rtmp

 call kinetic_basis_function(bftmp,bftmp,rtmp)
 rtmp2 = rtmp
 write(*,*) 'kinetic',rtmp
 call nucleus_pot_basis_function(bftmp,bftmp,zatom,rtmp)
 rtmp2 = rtmp2 + rtmp
 write(*,*) 'nucleus',rtmp
 write(*,*) 'matrix element basis function',rtmp2


 call nucleus_pot_gaussian(gatmp,gbtmp,zatom,rtmp)
 write(*,*) 'nucleus pot [Ha]',rtmp


 write(*,*)
 write(*,*) '                   END OF THE TESTS'
 write(*,*) '==========================================================='
 write(*,*) '==========================================================='
 write(*,*) '==========================================================='
 write(*,*) '                   START REAL LIFE CALCULATION'
 write(*,*)
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
 !
 ! 
 call read_inputparameter(calc_type,nspin,nscf,alpha_mixing,print_volume,basis_name,zatom,electrons,magnetization,basis_element)

 !
 ! start build up the basis set
 call init_basis_set(PRINT_VOLUME,zatom,basis_name,basis_element,basis)
 
 !
 ! allocate everything
 allocate(hamiltonian(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_kinetic(basis%nbf,basis%nbf,nspin))
 allocate(hamiltonian_nucleus(basis%nbf,basis%nbf,nspin))
 allocate(matrix(basis%nbf,basis%nbf,nspin))
 allocate(c_matrix(basis%nbf,basis%nbf,nspin))
 allocate(s_matrix(basis%nbf,basis%nbf))
 allocate(p_matrix(basis%nbf,basis%nbf,nspin))
 allocate(p_matrix_old(basis%nbf,basis%nbf,nspin))
 allocate(energy(basis%nbf,nspin))
 allocate(energy_old(basis%nbf,nspin))
 allocate(occupation(basis%nbf,nspin))
 allocate(exchange_m_vxc_diag(basis%nbf,nspin))
 allocate(self_energy_old(basis%nbf,basis%nbf,nspin))
 self_energy_old(:,:,:) = 0.0_dp
 if( calc_type%need_dft_xc ) allocate( vxc_matrix(basis%nbf,basis%nbf,nspin) )

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
     write(*,*) 'A diagonal term of the overlap matrix is not equal to 1.0'
     write(*,*) jbf,jbf,s_matrix(jbf,jbf)
     write(*,*) 'check the basis set definition'
     stop'ERROR'
   endif
 enddo
 title='=== overlap matrix S ==='
 call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,1,s_matrix)

 !
 ! Set up the electron repulsion integrals
 !
 ! ERI are stored privately in the module m_eri
 call start_clock(timing_integrals)
 call allocate_eri(basis%nbf)
 call calculate_eri2(basis)
 call stop_clock(timing_integrals)
 !
 ! In all the following cases, one needs the Coulomb integral
 ! in the eigenvector basis
 if( calc_type%is_gw .OR. calc_type%is_mp2 .OR. calc_type%type==CI) then
   call allocate_eri_eigen(nspin)
 endif


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

   call start_clock(timing_tmp3)
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
   call dump_out_matrix(PRINT_VOLUME,title,prod_basis%nbf_filtered,1,v_filtered_basis)
   call stop_clock(timing_tmp3)

   call start_clock(timing_tmp4)
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
   call dump_out_matrix(PRINT_VOLUME,title,prod_basis%nbf_filtered,1,s_filtered_basis)
   call stop_clock(timing_tmp4)

   !
   ! calculate S^{-1}
   call invert(prod_basis%nbf_filtered,s_filtered_basis,sinv_filtered_basis)

   !
   ! S^-1 V S^-1 is first calculated on the filtered product basis
   allocate(sinv_v_sinv_filtered(prod_basis%nbf_filtered,prod_basis%nbf_filtered))
   sinv_v_sinv_filtered = matmul( v_filtered_basis , sinv_filtered_basis )
   sinv_v_sinv_filtered = matmul( sinv_filtered_basis , sinv_v_sinv_filtered )
   title='=== S-1 V S-1 ==='
   call dump_out_matrix(PRINT_VOLUME,title,prod_basis%nbf_filtered,1,sinv_v_sinv_filtered)

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
 call guess_starting_c_matrix(basis%nbf,nspin,c_matrix)

 !
 ! Setup the density matrix: p_matrix
 call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
 title='=== 1st density matrix P ==='
 call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,p_matrix)

 !
 ! Kinetic energy contribution
 do jbf=1,basis%nbf
   do ibf=1,basis%nbf
     call kinetic_basis_function(basis%bf(ibf),basis%bf(jbf),energy_tmp)
     hamiltonian_kinetic(ibf,jbf,:) = energy_tmp
   enddo
 enddo
 title='=== kinetic energy contribution ==='
 call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,hamiltonian_kinetic)

 !
 ! nucleus-electron interaction
 do ibf=1,basis%nbf
   do jbf=1,basis%nbf
     call nucleus_pot_basis_function(basis%bf(ibf),basis%bf(jbf),zatom,energy_tmp)
     hamiltonian_nucleus(ibf,jbf,:) = energy_tmp
   enddo
 enddo
 title='=== nucleus contribution ==='
 call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,hamiltonian_nucleus)

 !
 ! start the big scf loop
 !
 do iscf=1,nscf
   write(*,'(/,a)') '-------------------------------------------'
   write(*,'(a,x,i4,/)') ' *** SCF cycle No:',iscf

   !
   ! Hartree contribution
   call start_clock(timing_hartree)
   matrix(:,:,:)=0.0_dp
!!! !$OMP PARALLEL DEFAULT(SHARED)
!!! !$OMP DO SCHEDULE(STATIC) REDUCTION(+:matrix)
   do lbf=1,basis%nbf
     do kbf=1,basis%nbf
       do jbf=1,basis%nbf
         do ibf=1,basis%nbf
           matrix(ibf,jbf,:) = matrix(ibf,jbf,:) &
                      + eri(ibf,jbf,kbf,lbf) * SUM( p_matrix(kbf,lbf,:) )
         enddo
       enddo
     enddo
   enddo
!!! !$OMP END DO
!!! !$OMP END PARALLEL
   en%kin  = SUM(hamiltonian_kinetic(:,:,:)*p_matrix(:,:,:))
   en%nuc  = SUM(hamiltonian_nucleus(:,:,:)*p_matrix(:,:,:))
   en%hart = 0.5_dp*SUM(matrix(:,:,:)*p_matrix(:,:,:))

   call stop_clock(timing_hartree)
   title='=== Hartree contribution ==='
   call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,matrix)
   hamiltonian(:,:,:) = hamiltonian_kinetic(:,:,:) + hamiltonian_nucleus(:,:,:) + matrix(:,:,:)
  
   !
   ! Exchange contribution
   if( calc_type%need_exchange ) then
     call start_clock(timing_exchange)

     spin_fact = REAL(-nspin+3,dp)
     matrix(:,:,:)=0.0_dp
     do ispin=1,nspin
       do lbf=1,basis%nbf
         do jbf=1,basis%nbf
           do kbf=1,basis%nbf
             do ibf=1,basis%nbf
               matrix(ibf,jbf,ispin) = matrix(ibf,jbf,ispin) &
                          - eri(ibf,kbf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact
             enddo
           enddo
         enddo
       enddo
     enddo
     en%exx = 0.5_dp*SUM(matrix(:,:,:)*p_matrix(:,:,:))
     call stop_clock(timing_exchange)

     title='=== Exchange contribution ==='
     call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,matrix)

     if( .NOT. calc_type%need_dft_xc ) then
       hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix(:,:,:) 
     else ! this is an hybrid functional
       write(*,*) 'this is an hybrid functional ',alpha_hybrid
       hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix(:,:,:) * alpha_hybrid
     endif
   endif

   !
   ! DFT XC potential is added here
   if( calc_type%need_dft_xc ) then
     call start_clock(timing_dft)
     if( .NOT. calc_type%need_exchange ) then
       call dft_exc_vxc(nspin,basis,(/calc_type%dft_x ,calc_type%dft_c/),p_matrix,vxc_matrix,en%xc)
     else 
       write(*,*) 'this is an hybrid functional ',alpha_hybrid
       matrix(:,:,:) = 0.0_dp
       call dft_exc_vxc(nspin,basis,(/calc_type%dft_x,0/),p_matrix,matrix,en%xc)
       call dft_exc_vxc(nspin,basis,(/calc_type%dft_c,0/),p_matrix,vxc_matrix,en%xc)
       vxc_matrix = vxc_matrix + (1.0_dp - alpha_hybrid) * matrix
     endif
     call stop_clock(timing_dft)
     hamiltonian(:,:,:) = hamiltonian(:,:,:) + vxc_matrix(:,:,:)

     title='=== DFT XC contribution ==='
     call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,vxc_matrix)
   endif

   !
   ! QPscGW self energy
   if( calc_type%is_gw .AND. calc_type%method == QS .AND. iscf > nscf/3 ) then

     call init_spectral_function(basis%nbf,prod_basis%nbf,nspin,occupation,wpol)
     call start_clock(timing_pola)
#ifdef AUXIL_BASIS
     call polarizability_casida(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,energy_tmp,wpol)
#else
     call polarizability_casida_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,energy_tmp,wpol)
#endif
     call stop_clock(timing_pola)
     write(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',energy_tmp

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
     call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,matrix)
     call destroy_spectral_function(wpol)

     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix(:,:,:)

   endif

   !
   ! QPscMP2
   if( calc_type%is_mp2 .AND. calc_type%method == QS .AND. iscf > nscf/3 ) then

!     call start_clock(timing_mp2_energy)
!     call mp2_energy(nspin,basis,occupation,c_matrix,energy,en%mp2)
!     call stop_clock(timing_mp2_energy)

     exchange_m_vxc_diag(:,:)=0.0_dp
     call start_clock(timing_mp2_self)
     call mp2_selfenergy(calc_type%method,nspin,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,matrix,en%mp2)
     call stop_clock(timing_mp2_self)
     write(*,'(a,2x,f12.6)') ' MP2 Energy       [Ha]:',en%mp2
     write(*,*) 
     en%tot = en%tot + en%mp2
     write(*,'(a,2x,f12.6)') ' MP2 Total Energy [Ha]:',en%tot

     matrix = alpha_mixing * matrix + (1.0_dp-alpha_mixing) * self_energy_old
     self_energy_old = matrix
     title='=== Self-energy ==='
     call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,matrix)
  
     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix(:,:,:)

   endif

  
   title='=== Total Hamiltonian ==='
   call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,hamiltonian)
  
   !
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \phi = E S \phi
   ! save the old eigenvalues
   energy_old = energy
   do ispin=1,nspin
    write(*,*) 'Diagonalization for spin polarization',ispin
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
   call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,c_matrix)
   ! test matrix equality C^T S C = I
   ! matrix(:,:,1) = matmul( transpose(c_matrix(:,:,1)), matmul( s_matrix(:,:), c_matrix(:,:,1) ) )
   ! title='=== C^T S C = identity ? ==='
   ! call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,1,matrix)

  
   !
   ! Setup the new density matrix: p_matrix
   ! Save the old one
   p_matrix_old(:,:,:) = p_matrix(:,:,:)
   call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
   title='=== density matrix P ==='
   call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,p_matrix)
  
   call check_convergence(basis%nbf,nspin,p_matrix_old,p_matrix,rms)
   write(*,*) 'convergence criterium on the density matrix',rms
   !
   ! Simple mixing
   p_matrix = alpha_mixing * p_matrix + ( 1.0_dp - alpha_mixing ) * p_matrix_old
   
!   energy = alpha_mixing * energy  + ( 1.0_dp - alpha_mixing ) * energy_old

 
   write(*,*)
   write(*,*) 'Kinetic Energy  [Ha]:',en%kin
   write(*,*) 'Nucleus Energy  [Ha]:',en%nuc
   write(*,*) 'Hartree Energy  [Ha]:',en%hart
   if(calc_type%need_exchange) write(*,*) 'Exchange Energy [Ha]:',en%exx
   if(calc_type%need_dft_xc)   write(*,*) 'XC Energy       [Ha]:',en%xc
   en%tot = en%kin + en%nuc + en%hart + en%exx + en%xc
   write(*,*)
   write(*,*) 'Total Energy    [Ha]:',en%tot
   write(*,*)
   
!   write(*,'(/,a,f15.8,/)') ' HF Total Energy [Ha]: ',&
!              hartree_fock_energy_band(basis%nbf,nspin,occupation,p_matrix,hamiltonian_kinetic+hamiltonian_nucleus,energy)

 !
 ! end of the big SCF loop
 enddo
 if(rms>1.d-4) then
   msg='SCF convergence is poor'
   call issue_warning(msg)
 endif

 write(*,*) '=================================================='
 write(*,*) 'The SCF loop ends here'
 write(*,*) '=================================================='
 write(*,*)

#if 0
 call plot_wfn(nspin,basis,c_matrix)
#endif

 !
 ! CI calculation is done here
 ! implemented for 2 electrons only!
 if(calc_type%type==CI) then
   if (nspin/=1) stop'for CI, nspin should be 1'
   call full_ci_2electrons_spin(0,basis%nbf,hamiltonian_kinetic+hamiltonian_nucleus,c_matrix)
 endif
  
 !
 ! in case of DFT + GW
 if( calc_type%need_final_exchange ) then

   spin_fact = REAL(-nspin+3,dp)
   matrix(:,:,:)=0.0_dp
   do ispin=1,nspin
     do ibf=1,basis%nbf
       do jbf=1,basis%nbf
         do kbf=1,basis%nbf
           do lbf=1,basis%nbf
             matrix(ibf,jbf,ispin) = matrix(ibf,jbf,ispin) &
                        - eri(ibf,kbf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact
           enddo
         enddo
       enddo
     enddo
   enddo

   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,basis%nbf
       do ibf=1,basis%nbf
         do jbf=1,basis%nbf
           exchange_m_vxc_diag(istate,ispin) = exchange_m_vxc_diag(istate,ispin) &
                   + c_matrix(ibf,istate,ispin) * ( matrix(ibf,jbf,ispin) - vxc_matrix(ibf,jbf,ispin) ) * c_matrix(jbf,istate,ispin)
         enddo
       enddo
     enddo
   enddo

#if 0
   msg='hacking here'
   call issue_warning(msg)
   occupation(1:5,:)=0.0_dp
   call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
#endif

   en%exx = 0.5_dp*SUM(matrix(:,:,:)*p_matrix(:,:,:))
   write(*,*) 'EXX [Ha]:',en%exx

 else
   exchange_m_vxc_diag(:,:) = 0.0_dp
 endif

 !
 ! final evaluation for G0W0
 if( calc_type%is_gw .AND. calc_type%method == perturbative ) then

   call init_spectral_function(basis%nbf,prod_basis%nbf,nspin,occupation,wpol)
   call start_clock(timing_pola)
#ifdef AUXIL_BASIS
   call polarizability_casida(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,energy_tmp,wpol)
#else
   call polarizability_casida_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,energy_tmp,wpol)
#endif
   call stop_clock(timing_pola)
   write(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',energy_tmp
   call start_clock(timing_self)
#ifdef AUXIL_BASIS
   call gw_selfenergy_casida(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#else
   call gw_selfenergy_casida_noaux(calc_type%method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix)
#endif
   call stop_clock(timing_self)

   title='=== Self-energy === (in the orbital basis)'
   call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,matrix)
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
   write(*,'(a,2x,f12.6)') ' MP2 Energy       [Ha]:',en%mp2
   write(*,*) 
   en%tot = en%tot + en%mp2
   write(*,'(a,2x,f12.6)') ' MP2 Total Energy [Ha]:',en%tot

   title='=== Self-energy === (in the orbital basis)'
   call dump_out_matrix(PRINT_VOLUME,title,basis%nbf,nspin,matrix)

 endif

 deallocate(hamiltonian,hamiltonian_kinetic,hamiltonian_nucleus)
 deallocate(matrix,s_matrix,c_matrix,p_matrix,p_matrix_old)
 deallocate(energy,occupation,exchange_m_vxc_diag)
 deallocate(self_energy_old)
!ERI deallocate(eri)
 call deallocate_eri()
 if( calc_type%need_dft_xc ) deallocate( vxc_matrix )

 call destroy_basis_set(basis)
 if(calc_type%is_gw) call destroy_basis_set(prod_basis)

 call stop_clock(timing_total)
 call output_timing()

 call output_all_warnings()

 write(*,'(/,a,/)') ' This is the end'

end program atom

!=========================================================================
subroutine setup_density_matrix(nbf,nspin,c_matrix,occupation,p_matrix)
 use m_definitions
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
subroutine check_convergence(nbf,nspin,p_matrix_old,p_matrix,rms)
 use m_definitions
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix_old(nbf,nbf,nspin),p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: rms
!=====

 rms = SQRT( SUM( ( p_matrix_old(:,:,:) - p_matrix(:,:,:) )**2 ) )

end subroutine check_convergence

!=========================================================================
subroutine  set_occupation(electrons,magnetization,nbf,nspin,occupation)
 use m_definitions
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
    write(*,*)
    write(*,*) 'occupations are read from file: manual_occupations'
    msg='reading manual occupations from file'
    call issue_warning(msg)
    open(unit=12,file='manual_occupations',status='old')
    !
    ! read nlines, all other occupations are set to zero
    read(12,*) nlines
    do ilines=1,nlines
      read(12,*) occupation(ilines,:)
    enddo
    close(12)
    write(*,*) 'occupations set, closing file'
  endif
 
  !
  ! final check
  if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0d-7 ) then
    write(*,*) 'occupation set up failed to give the right number of electrons'
    write(*,*) 'sum of occupations',SUM(occupation(:,:))
    write(*,*) 'electrons',electrons
    do ibf=1,nbf
      write(*,*) ibf,occupation(ibf,:)
    enddo
    stop'FAILURE in set_occupations'
  endif 

end subroutine set_occupation

!=========================================================================
function hartree_fock_energy_band(nbf,nspin,occupation,p_matrix,hamiltonian_constant,energy)
 use m_definitions
 implicit none

 integer,intent(in)  :: nbf,nspin
 real(dp),intent(in) :: hamiltonian_constant(nbf,nbf,nspin)
 real(dp),intent(in) :: energy(nbf,nspin),p_matrix(nbf,nbf,nspin)
 real(dp),intent(in) :: occupation(nbf,nspin)
 real(dp)            :: hartree_fock_energy_band
 integer             :: ibf,jbf,ispin
!=====

 hartree_fock_energy_band = 0.0_dp
 do ispin=1,nspin
   do ibf=1,nbf
     hartree_fock_energy_band = hartree_fock_energy_band &
            + 0.5_dp * occupation(ibf,ispin) * energy(ibf,ispin)
     do jbf=1,nbf
       hartree_fock_energy_band =  hartree_fock_energy_band &
              + 0.5_dp * p_matrix(ibf,jbf,ispin) * hamiltonian_constant(ibf,jbf,ispin)
     enddo
   enddo
 enddo

end function hartree_fock_energy_band

!=========================================================================
subroutine guess_starting_c_matrix(nbf,nspin,c_matrix)
 use m_definitions
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
 enddo

end subroutine guess_starting_c_matrix

