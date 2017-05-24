!=========================================================================
! This file is part of MOLGW.
! Author: Ivan Maliyov
!
! This module contains
! propagation of the wavefunction in time
!
!=========================================================================
module m_tddft_propagator
 use m_atoms
 use m_definitions
 use m_basis_set 
 use m_scf_loop
 use m_memory
 use m_hamiltonian
 use m_hamiltonian_sca
 use m_hamiltonian_onebody
 use m_hamiltonian_buffer
 use m_hamiltonian_cmplx
 use m_inputparam
 use m_dft_grid
 use m_tools
 use m_scf
 use m_warning

 interface propagate_orth
  module procedure propagate_orth_ham_1
  module procedure propagate_orth_ham_2
 end interface propagate_orth

 interface print_2d_matrix
  module procedure print_2d_matrix_real
  module procedure print_2d_matrix_cmplx
 end interface print_2d_matrix

! real(dp),allocatable       :: dipole_basis(:,:,:)
 real(dp)                   :: time_read
 real(dp)                   :: excit_dir_norm(3)
 real(dp),allocatable       :: s_matrix_inv(:,:)
 real(dp),allocatable       :: xatom_start(:,:)
 complex(dp)                :: m_excit_field_dir
 integer                    :: nocc
contains


!=========================================================================
subroutine calculate_propagation(nstate,              &  
                                 basis,               &  
                                 occupation,          &   
                                 s_matrix,            & 
                                 s_matrix_sqrt_inv,   &
                                 c_matrix,            &
                                 hamiltonian_kinetic, & 
                                 hamiltonian_nucleus)
 use ISO_C_BINDING
 use m_timing
 implicit none

!=====
 type(basis_set),intent(in)      :: basis
 integer,intent(in)              :: nstate 
 real(dp),intent(in)             :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: c_matrix(basis%nbf,nstate,nspin) 
 real(dp),intent(in)             :: occupation(nstate,nspin)
 real(dp),intent(in)             :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(inout)          :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: s_matrix_sqrt_inv(basis%nbf,nstate)
!=====
 integer,parameter          :: BATCH_SIZE=64
 integer                    :: ntau, itau,idir, info, ispin, ibf,nomega,iomega
 integer                    :: istate
 integer                    :: nwrite_step, iwrite_step
 integer                    :: file_dipolar_spectra
 integer                    :: file_real_dipolar_spectra
 integer                    :: file_aimag_dipolar_spectra, file_transforms,file_eexcit
 real(dp)                   :: time_cur, time_min=0.0_dp
 real(dp)                   :: omega_factor
 real(dp),allocatable       :: dipole_basis(:,:,:)
 complex(dp),allocatable    :: dipole_time_ref(:,:)
 type(C_PTR)                :: plan
!===== variables for the calc_p_matrix_error
 complex(dp),allocatable    :: p_matrix_time_test_cmplx(:,:,:,:)
 complex(dp),allocatable    :: p_matrix_time_ref_cmplx(:,:,:,:)
 complex(dp),allocatable    :: dipole_time_test(:,:)
 real(dp),allocatable       :: m_time_steps(:)
 integer,allocatable        :: m_n_hists(:)
 integer,allocatable        :: m_n_iters(:)
 integer                    :: istep, iprop_type,ipred_corr,ihist,iiter
 character(len=4),allocatable  :: m_prop_types(:), m_pred_corrs(:)
 integer                    :: file_p_matrix_error, file_dipole_error
 character(len=100)         :: name_p_matrix_error, name_dipole_error
 integer                    :: iwrite, nwrite, mod_write, mod_write_cur 
 integer                    :: n_time_steps,n_prop_types,n_pred_corrs,n_n_hists,n_n_iters
 
 nocc=0
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     if( istate > nocc ) nocc = istate
   enddo
 end do

 if(write_step / time_step - NINT( write_step / time_step ) > 0.0E-10_dp) then
   call die("Tddft error: write_step is not multiple of time_step.")
 end if
 mod_write = NINT( write_step / time_step ) ! write each write_step, assumed that time_step <= write_step

 if( calc_type%is_dft ) then
    call init_dft_grid(basis,grid_level,dft_xc_needs_gradient,.TRUE.,BATCH_SIZE)
 endif

 excit_dir_norm(:)=excit_type%dir(:)/NORM2(excit_type%dir(:))

 time_read=0.0_dp
 time_min=0.0_dp
 allocate(s_matrix_inv(basis%nbf,basis%nbf))
 if(excit_type%is_light) then
   allocate(dipole_basis(basis%nbf,basis%nbf,3))
   call calculate_dipole_basis(basis,dipole_basis)
 end if
 call invert(basis%nbf,s_matrix,s_matrix_inv)  

 ntau=NINT((time_sim-time_min)/time_step)
 nwrite_step=NINT((time_sim - time_min)/write_step) 
 if(excit_type%is_light) then
   allocate(dipole_time_ref(nwrite_step,3))
   dipole_time_ref(:,:)= ( 0.0_dp , 0.0_dp )
 end if
 if(calc_p_matrix_error_) then
   allocate(p_matrix_time_ref_cmplx(basis%nbf,basis%nbf,nspin,nwrite_step))
   allocate(p_matrix_time_test_cmplx(basis%nbf,basis%nbf,nspin,nwrite_step))
   if(excit_type%is_light) allocate(dipole_time_test(nwrite_step,3))
 end if

 write(stdout,"(x,A,F8.2,A,F9.5,A,I8)") "Calculate the reference tddft loop for time_sim=", time_sim, " time_step=", time_step, " number of iterations", ntau
 call start_clock(timing_tddft_loop)
 call tddft_time_loop(nstate,                           &
                      basis,                            &
                      occupation,                       &
                      dipole_basis,                     &
                      s_matrix,                         &
                      s_matrix_sqrt_inv,                &
                      c_matrix,                         &
                      hamiltonian_kinetic,              &
                      hamiltonian_nucleus,              &
                      pred_corr,                        &
                      prop_type,                        &
                      time_step,                        &
                      n_hist,                           &
                      n_iter,                           &
                      p_matrix_time_ref_cmplx,          &
                      ntau,                             &
                      dipole_time_ref,                  &
                      .true.)
 call stop_clock(timing_tddft_loop)
 write(stdout,*) "End of the reference tddft loop"
 write(stdout,*)
 
 if( calc_p_matrix_error_ ) then

   call get_number_of_elements(error_pred_corrs,n_pred_corrs)
   allocate(m_pred_corrs(n_pred_corrs))
   read(error_pred_corrs,*)m_pred_corrs(:)

   call get_number_of_elements(error_prop_types,n_prop_types)
   allocate(m_prop_types(n_prop_types))
   read(error_prop_types,*)m_prop_types(:) 
  
   call get_number_of_elements(error_time_steps,n_time_steps)
   allocate(m_time_steps(n_time_steps))
   read(error_time_steps,*)m_time_steps(:)

   call get_number_of_elements(error_n_hists,n_n_hists)
   allocate(m_n_hists(n_n_hists))
   read(error_n_hists,*)m_n_hists(:)

   call get_number_of_elements(error_n_iters,n_n_iters)
   allocate(m_n_iters(n_n_iters))
   read(error_n_iters,*)m_n_iters(:)

   do istep=1, size(m_time_steps)
     ntau=NINT((time_sim-time_min)/m_time_steps(istep))
     call flush(stdout)
     mod_write_cur = NINT( write_step / m_time_steps(istep) )
     do ipred_corr=1, size(m_pred_corrs)
       do iprop_type=1, size(m_prop_types)
         do ihist=1, size(m_n_hists)
           do iiter=1, size(m_n_iters)
             if(write_step / m_time_steps(istep) - NINT( write_step / m_time_steps(istep) ) > 0.0E-10_dp) then
               write(stdout,*) "The time_step ", m_time_steps(istep), " is not a multiple of write_step: ", write_step
               call die("Tddft error: write_step is not multiple of one of time steps for the p_matrix error calculation.")
             end if
             if( m_time_steps(istep)>write_step) then
               write(stdout,*) "The time_step ", m_time_steps(istep), " is larger than write_step: ", write_step
               call die("Tddft error: one of time steps for the p_matrix error calculation is larger than write_step.")
             end if
             write(stdout,"(x,5A,F9.5,A,I1,A,I1)") "Start of tddft loop for the p_matrix error for pred_corr: ", & 
                      m_pred_corrs(ipred_corr), "; prop_type: ", m_prop_types(iprop_type), " time_step: ", m_time_steps(istep), &
                      " n_hist: ",m_n_hists(ihist), " and n_iter: ", m_n_iters(iiter)
             call tddft_time_loop(nstate,                           &
                                  basis,                            &
                                  occupation,                       &
                                  dipole_basis,                     &
                                  s_matrix,                         &
                                  s_matrix_sqrt_inv,                &
                                  c_matrix,                         &
                                  hamiltonian_kinetic,              &
                                  hamiltonian_nucleus,              &
                                  m_pred_corrs(ipred_corr),         &
                                  m_prop_types(iprop_type),         &
                                  m_time_steps(istep),              &
                                  m_n_hists(ihist),                 &
                                  m_n_iters(iiter),                 &
                                  p_matrix_time_test_cmplx,         &
                                  ntau,                             &
                                  dipole_time_test,                 &
                                  .false.)
             write(stdout,"(x,A)") "End of the tddft loop"
             if( is_iomaster ) then
               write(name_p_matrix_error,'(5A,F5.3,A,I1,A,I1,A)') "p_matrix_error_", TRIM(m_pred_corrs(ipred_corr)), "_", TRIM(m_prop_types(iprop_type)), &
                     "_dt_", m_time_steps(istep),"_hist_",m_n_hists(ihist), "_iter_",m_n_iters(iiter),".dat" 
               open(newunit=file_p_matrix_error,file=name_p_matrix_error)
               if(excit_type%is_light) then
                 write(name_dipole_error,'(5A,F5.3,A,I1,A,I1,A)') "dipole_error_", TRIM(m_pred_corrs(ipred_corr)), "_", TRIM(m_prop_types(iprop_type)), & 
                       "_dt_", m_time_steps(istep),"_hist_",m_n_hists(ihist),"_iter_",m_n_iters(iiter), ".dat"
                 open(newunit=file_dipole_error,file=name_dipole_error)
               end if
               time_min = time_read + write_step
               time_cur = time_min
               do iwrite_step=1, nwrite_step 
                 write(file_p_matrix_error,*) time_cur, SUM(ABS(p_matrix_time_ref_cmplx(:,:,:,iwrite_step) - p_matrix_time_test_cmplx(:,:,:,iwrite_step))) / (basis%nbf)**2
                 if(excit_type%is_light) then
                   write(file_dipole_error,*) time_cur, SUM(ABS(dipole_time_ref(iwrite_step,:) - dipole_time_test(iwrite_step,:))) / (basis%nbf)**2
                 end if
                 time_cur=time_min+iwrite_step*write_step
               end do
               close(file_p_matrix_error)     
               close(file_dipole_error)
             end if
           end do ! n_iter
         end do ! n_hist
       end do ! prop_type
     end do ! pred_corr
   end do !time_step

 deallocate(p_matrix_time_test_cmplx)
 end if

 if( calc_type%is_dft ) call destroy_dft_grid()

end subroutine calculate_propagation



!=======================================
subroutine tddft_time_loop(nstate,                           &
                           basis,                            &
                           occupation,                       &
                           dipole_basis,                     &
                           s_matrix,                         &
                           s_matrix_sqrt_inv,                &
                           c_matrix,                         &
                           hamiltonian_kinetic,              &
                           hamiltonian_nucleus,              &
                           pred_corr_cur,                    &
                           prop_type_cur,                    &
                           time_step_cur,                    &
                           n_hist_cur,                       &
                           n_iter_cur,                       &
                           p_matrix_time_cmplx,              &
                           ntau,                             &
                           dipole_time,                      &
                           ref_)

 implicit none

!=====
 type(basis_set),intent(in)             :: basis
 integer,intent(in)                     :: nstate
 integer,intent(in)                     :: n_hist_cur
 integer,intent(in)                     :: n_iter_cur
 real(dp),intent(in)                    :: s_matrix(basis%nbf,basis%nbf)
 real(dp),allocatable,intent(in)        :: dipole_basis(:,:,:)
 real(dp),intent(in)                    :: c_matrix(basis%nbf,nstate,nspin) 
 real(dp),intent(in)                    :: occupation(nstate,nspin)
 real(dp),intent(in)                    :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(inout)                 :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),intent(in)                    :: s_matrix_sqrt_inv(basis%nbf,nstate)
 real(dp),intent(in)                    :: time_step_cur
 character(len=4),intent(in)            :: prop_type_cur,pred_corr_cur
 complex(dp),allocatable,intent(inout)  :: p_matrix_time_cmplx(:,:,:,:)
 complex(dp),allocatable,intent(inout)  :: dipole_time(:,:)
 logical,intent(in)                     :: ref_
!=====
 integer                    :: itau, idir, info, ispin, ibf, ntau, mod_write
 integer                    :: i_iter, iwrite_step
 integer                    :: file_time_data, file_excit_field
 integer                    :: file_dipole_time,file_iter_norm 
 integer                    :: file_q_matrix_ii(2),file_tmp(10)
 integer                    :: min_elem_q_mat,max_elem_q_mat,i_elem_q_mat
 real(dp)                   :: time_cur, excit_field(3),time_min, time_one_iter
 real(dp)                   :: dipole(3)
 real(dp)                   :: energies_inst(nstate)
 complex(dp),allocatable    :: c_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: c_matrix_0_cmplx(:,:,:)
 complex(dp),allocatable    :: c_matrix_orth_cmplx(:,:,:)
 complex(dp),allocatable    :: q_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: hamiltonian_fock_cmplx(:,:,:)
 complex(dp),allocatable    :: p_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: h_small_hist_cmplx(:,:,:,:)
 complex(dp),allocatable    :: c_matrix_orth_hist_cmplx(:,:,:,:)
 complex(dp),allocatable    :: h_small_cmplx(:,:,:)
 complex(dp),allocatable    :: c_matrix_buf_cmplx(:,:,:)
 character(len=50)          :: name_time_data,name_dipole_time
 character(len=50)          :: name_iter_norm
 character(len=50)          :: name_file_q_matrix_ii
 logical                    :: calc_excit_ 
 logical                    :: is_identity_
!==variables for extrapolation
 integer                    :: iextr,ham_dim_cur
 real(dp),allocatable       :: m_nods(:)
 real(dp),allocatable       :: extrap_coefs(:)
 real(dp)                   :: x_pred 
 logical                    :: restart_is_correct
!=====

 mod_write = NINT( write_step / time_step_cur )

 min_elem_q_mat=1
 max_elem_q_mat=10
 max_elem_q_mat=MIN(nocc,max_elem_q_mat)

 allocate(c_matrix_cmplx(basis%nbf,nocc,nspin))
 allocate(c_matrix_0_cmplx(basis%nbf,nstate,nspin))
 allocate(c_matrix_orth_cmplx(nstate,nocc,nspin))
! allocate(q_matrix_cmplx(nstate,nocc,nspin))
 allocate(hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(p_matrix_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(h_small_cmplx(nstate,nstate,nspin))
 allocate(c_matrix_buf_cmplx(basis%nbf,nstate,nspin))

!OPENING FILES
 if( is_iomaster ) then
   if(ref_) then
     open(newunit=file_time_data,file="time_data.dat")
     do ispin=1,nspin
       write(name_file_q_matrix_ii,"(a,i1,a)") "q_matrix_ii_spin_", ispin, ".dat" 
!       open(newunit=file_q_matrix_ii(ispin),file=name_file_q_matrix_ii)
     end do
     if(excit_type%is_light) then
       open(newunit=file_dipole_time,file="dipole_time.dat")
       open(newunit=file_excit_field,file="excitation_time.dat")
       write(file_excit_field,*) "# time(au)                      E_field_excit_dir(au)"
     end if
   else
     write(name_time_data,'(5A,F5.3,A,I1,A,I1,A)') "time_data_", TRIM(pred_corr_cur), "_", TRIM(prop_type_cur), "_dt_", time_step_cur, &
                   "_hist_",n_hist_cur, "_iter_",n_iter_cur,".dat"
     if(excit_type%is_light) then
       write(name_dipole_time,'(5A,F5.3,A,I1,A,I1,A)') "dipole_time_", TRIM(pred_corr_cur), "_", TRIM(prop_type_cur), "_dt_", time_step_cur, &
                     "_hist_",n_hist_cur,"_iter_",n_iter_cur,".dat"
       open(newunit=file_dipole_time,file=name_dipole_time)
     end if
     open(newunit=file_time_data,file=name_time_data)
   end if  
   write(file_time_data,"(A)") " # time(au)     e_total             enuc_nuc             enuc            ekin               ehart            &
                               eexx_hyb            exc             eexcit            P*S trace (# of els)" 
   if(excit_type%is_light) write(file_dipole_time,*) "# time(au)                      Dipole_x(D)               Dipole_y(D)               Dipole_z(D)"
 end if

 time_min=0.0_dp
!===INITIAL CONDITIONS===
 if(ref_) then
   allocate(xatom_start(3,natom))
   xatom_start=xatom
 end if
 ! Getting c_matrix_cmplx(t=0) whether using RESTART_TDDFT file, whether using c_matrix
 if(.NOT. ignore_tddft_restart_) then
   call read_restart_tddft(nstate,time_read,c_matrix_orth_cmplx,restart_is_correct)
   time_min = time_read
   do ispin=1,nspin
     c_matrix_cmplx(:,:,ispin) = MATMUL( s_matrix_sqrt_inv(:,:) , c_matrix_orth_cmplx(:,:,ispin) )
   end do
 end if

! call print_2d_matrix("c_matrix_orth_cmplx of beginning",c_matrix_orth_cmplx(:,:,1),basis%nbf,nocc,stdout,4)
! call print_2d_matrix("c_matrix_cmplx of beginning",c_matrix_cmplx(:,:,1),basis%nbf,nocc,stdout,4)

 if(ignore_tddft_restart_ .OR. (.NOT. restart_is_correct)) then
   c_matrix_cmplx(1:basis%nbf,1:nocc,1:nspin)=c_matrix(1:basis%nbf,1:nocc,1:nspin)
 end if

 c_matrix_0_cmplx(:,:,:)=c_matrix(:,:,:)
! do ispin=1,nspin
!   q_matrix_cmplx(:,:,ispin)=MATMUL(MATMUL(CONJG(TRANSPOSE(c_matrix_0_cmplx(:,:,ispin))),s_matrix),c_matrix_cmplx(:,:,ispin))
!!   call print_square_2d_matrix_cmplx("q_matrix_cmplx",q_matrix_cmplx,nstate,stdout,4)
! end do

 ! Getting starting value of the Hamiltonian
 ! itau=0 to avoid excitation calculation
 call setup_hamiltonian_fock_cmplx( basis,                   &
                                    nstate,                  &
                                    0,                       &
                                    time_min,                &
                                    time_step_cur,           &
                                    occupation,              &
                                    c_matrix_cmplx,          &
                                    hamiltonian_kinetic,     &
                                    hamiltonian_nucleus,     &
                                    h_small_cmplx,           &
                                    s_matrix_sqrt_inv,       &
                                    dipole_basis,            &
                                    hamiltonian_fock_cmplx,  &
                                    ref_)
 
 ! In case of no restart, find the c_matrix_orth_cmplx by diagonalizing h_small
 if(ignore_tddft_restart_ .OR. (.NOT. restart_is_correct)) then
   do ispin=1, nspin
     call diagonalize(nstate,h_small_cmplx(:,:,ispin),energies_inst(:),c_matrix_buf_cmplx(1:basis%nbf,1:nstate,ispin))
   end do
   c_matrix_orth_cmplx(1:nstate,1:nocc,1:nspin)=c_matrix_buf_cmplx(1:nstate,1:nocc,1:nspin)
 end if

 !Printing starting values of energy and dipole taken from SCF of RESTART_TDDFT
 call setup_density_matrix_cmplx(basis%nbf,nstate,nocc,c_matrix_cmplx,occupation,p_matrix_cmplx)

 en%kin = real(SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
 en%nuc = real(SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
 en%tot = en%nuc + en%kin + en%nuc_nuc + en%hart + en%exx_hyb + en%xc + en%excit

 if(excit_type%is_light) call static_dipole_fast_cmplx(basis,p_matrix_cmplx,dipole_basis,dipole)

 if( is_iomaster ) then
 ! Here time_min point coresponds to the end of calculation written in the RESTART_TDDFT or to 0 a.u.
   if(excit_type%is_projectile) call output_projectile_position()
   if( print_cube_rho_tddft_ ) call plot_cube_wfn_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0)
   if(excit_type%is_light) write(file_dipole_time,*) time_min, dipole(:) * au_debye
   write(file_time_data,"(F9.4,8(2x,es16.8E3),2x,2(2x,F7.2))") &
      time_min, en%tot, en%nuc_nuc, en%nuc, en%kin, en%hart, en%exx_hyb, en%xc, en%excit, matrix_trace_cmplx(MATMUL(p_matrix_cmplx(:,:,1),s_matrix(:,:)))
   write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Simulation time (au):', time_min
   write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Total Energy    (Ha):',en%tot
   if(excit_type%is_light) write(stdout,'(a31,1x,3f19.10)') 'RT-TDDFT Dipole Moment    (D):', dipole(:) * au_debye

 end if
! if( .NOT. ignore_tddft_restart_ .AND. restart_is_correct ) then
 time_min=time_min+time_step_cur
 
! call print_square_2d_matrix_cmplx("c_matrix_cmplx of the beginning",c_matrix_cmplx,nstate,stdout,4)

! call print_2d_matrix("check of beginning",MATMUL(MATMUL(s_matrix(:,:) ,c_matrix_cmplx(:,:,nspin) ) ,TRANSPOSE(CONJG(c_matrix_cmplx(:,:,nspin)))  ),basis%nbf,nocc,stdout,4)

 if(pred_corr_cur /= 'PC0' ) then

   allocate(m_nods(n_hist_cur),extrap_coefs(n_hist_cur))

   if(pred_corr_cur=='PC1') then
     ham_dim_cur=2
   end if

   if(pred_corr_cur=='PC2') then
     ham_dim_cur=n_hist_cur
     do iextr=1,n_hist_cur
       m_nods(iextr)=0.5_dp+(iextr-1.0_dp)
     end do
     x_pred=0.5_dp+(n_hist_cur-1.0_dp)+3.0_dp/4.0_dp
     call get_extrap_coefs_lagr(m_nods,x_pred,extrap_coefs,n_hist_cur)
   end if

   if(pred_corr_cur=='PC2B' .OR. pred_corr_cur=='PC2C') then
     ham_dim_cur=n_hist_cur
     do iextr=1,n_hist_cur
       m_nods(iextr)=(iextr-1.0_dp)*0.5_dp
     end do
     x_pred=(n_hist_cur-1.0_dp)*0.5_dp+0.25_dp
     if(pred_corr_cur=='PC2B') call get_extrap_coefs_lagr(m_nods,x_pred,extrap_coefs,n_hist_cur)
     if(pred_corr_cur=='PC2C') call get_extrap_coefs_aspc(extrap_coefs,n_hist_cur)
   end if

   if(pred_corr_cur=='PC2D' .OR. pred_corr_cur=='PC2E') then
     ham_dim_cur=n_hist_cur+1
     do iextr=1,n_hist_cur
       m_nods(iextr)=(iextr-1.0_dp)*0.5_dp
     end do
     x_pred=(n_hist_cur-1.0_dp)*0.5_dp+0.5_dp
     if(pred_corr_cur=='PC2D') call get_extrap_coefs_lagr(m_nods,x_pred,extrap_coefs,n_hist_cur)
     if(pred_corr_cur=='PC2E') call get_extrap_coefs_aspc(extrap_coefs,n_hist_cur)
   end if

   if(pred_corr_cur=='PC3' .OR. pred_corr_cur=='PC4') then
     ham_dim_cur=n_hist_cur+2
     do iextr=1,n_hist_cur
       m_nods(iextr)=iextr-1.0_dp
     end do
     x_pred=n_hist_cur
     if(pred_corr_cur=='PC3') call get_extrap_coefs_lagr(m_nods,x_pred,extrap_coefs,n_hist_cur)
     if(pred_corr_cur=='PC4') call get_extrap_coefs_aspc(extrap_coefs,n_hist_cur)
   end if

   if(pred_corr_cur=='PC5' .OR. pred_corr_cur=='PC6') then
     ham_dim_cur=n_hist_cur+1
     do iextr=1,n_hist_cur
       m_nods(iextr)=iextr-1.0_dp
     end do
     x_pred=n_hist_cur
     if(pred_corr_cur=='PC5') call get_extrap_coefs_lagr(m_nods,x_pred,extrap_coefs,n_hist_cur)
     if(pred_corr_cur=='PC6') call get_extrap_coefs_aspc(extrap_coefs,n_hist_cur)
   end if

   if(pred_corr_cur=='PC7' ) then
     ham_dim_cur=2
   end if

   allocate(h_small_hist_cmplx(nstate,nstate,nspin,ham_dim_cur))
   allocate(c_matrix_orth_hist_cmplx(nstate,nocc,nspin,1))
   do iextr=1,ham_dim_cur
     h_small_hist_cmplx(:,:,:,iextr)=h_small_cmplx(:,:,:)
   end do
   c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
   m_nods(:)=0.0_dp
 end if
!===END INITIAL CONDITIONS===

!********start time loop*************
 time_cur=time_min
 iwrite_step = 1
 itau = 1
 in_tddft_loop=.true.
 do while ( (time_cur - time_sim) < 1.0e-10 )
   if(itau==3) call start_clock(timing_tddft_one_iter)
   if(pred_corr_cur=='PC0') then
     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)
     call setup_hamiltonian_fock_cmplx( basis,                   &
                                        nstate,                  &
                                        itau,                    &
                                        time_cur,                &
                                        time_step_cur,           &
                                        occupation,              &
                                        c_matrix_cmplx,          &
                                        hamiltonian_kinetic,     &
                                        hamiltonian_nucleus,     &
                                        h_small_cmplx,           &
                                        s_matrix_sqrt_inv,       &
                                        dipole_basis,            &
                                        hamiltonian_fock_cmplx,  &
                                        ref_)
   end if
   ! ///////////////////////////////////
   ! Following Van Voorhis, Phys Rev B 74 (2006)
     if(pred_corr_cur=='PC1') then
     !--1--PREDICTOR----| H(2/4),H(6/4)-->H(9/4)
         h_small_cmplx= -3.0_dp/4.0_dp*h_small_hist_cmplx(:,:,:,1)+7.0_dp/4.0_dp*h_small_hist_cmplx(:,:,:,2)  
     !--2--PREDICTOR----| C(8/4)---U[H(9/4)]--->C(10/4)
     call propagate_orth(nstate,basis,time_step_cur/2.0_dp,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)  

     !--3--CORRECTOR----| C(10/4)-->H(10/4)
     call setup_hamiltonian_fock_cmplx( basis,                   &
                                        nstate,                  &
                                        itau,                    &
                                        time_cur,                &
                                        time_step_cur,           &
                                        occupation,              &
                                        c_matrix_cmplx,          &
                                        hamiltonian_kinetic,     &
                                        hamiltonian_nucleus,     &
                                        h_small_cmplx,           &
                                        s_matrix_sqrt_inv,       &
                                        dipole_basis,            &
                                        hamiltonian_fock_cmplx,  &
                                        ref_)
       
     !--4--PROPAGATION----| C(8/4)---U[H(10/4)]--->C(12/4)
     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)

     !--5--UPDATE----| C(12/4)-->C(8/4); H(6/4)-->H(2/4); H(10/4)-->H(6/4)
     c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
     h_small_hist_cmplx(:,:,:,1)=h_small_hist_cmplx(:,:,:,2)
     h_small_hist_cmplx(:,:,:,2)=h_small_cmplx(:,:,:)
   end if

   ! ///////////////////////////////////
   ! Lagrange interpolation, the rest is the same as for PC1
   if(pred_corr_cur=='PC2') then 
     hamiltonian_fock_cmplx=(0.0_dp,0.0_dp)
     h_small_cmplx=(0.0_dp,0.0_dp)
     !--1--PREDICTOR---- 
     do iextr=1,n_hist_cur
       h_small_cmplx=h_small_cmplx+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
     end do
     !--2--PREDICTOR----| hamiltonian_fock_cmplx in the 1/4 of the current time interval

     call propagate_orth(nstate,basis,time_step_cur/2.0_dp,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)  
     !--3--CORRECTOR----| C(1/2)-->H(1/2)

     call setup_hamiltonian_fock_cmplx( basis,                   &
                                        nstate,                  &
                                        itau,                    &
                                        time_cur,                &
                                        time_step_cur,           &
                                        occupation,              &
                                        c_matrix_cmplx,          &
                                        hamiltonian_kinetic,     &
                                        hamiltonian_nucleus,     &
                                        h_small_cmplx,           &
                                        s_matrix_sqrt_inv,       &
                                        dipole_basis,            &
                                        hamiltonian_fock_cmplx,  &
                                        ref_)
       
     !--4--PROPAGATION----| C(0)---U[H(1/2)]--->C(1)
     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)

     !--5--UPDATE----| C(12/4)-->C(8/4); H(6/4)-->H(2/4); H(10/4)-->H(6/4)
     c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
  
     do iextr=1,n_hist_cur-1
       h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+1)
     end do
     h_small_hist_cmplx(:,:,:,n_hist_cur)=h_small_cmplx(:,:,:)
   end if

   ! ///////////////////////////////////

   if(pred_corr_cur=='PC2B' .OR. pred_corr_cur=='PC2C') then 
     h_small_cmplx=(0.0_dp,0.0_dp)
     !--0--EXTRAPOLATE---- 
     do iextr=1,n_hist_cur
       h_small_cmplx=h_small_cmplx+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
     end do
     ! n_hist_cur-2  <---  n_hist_cur; n_hist_cur-3 <-- n_hist_cur-1
     if(n_hist_cur > 2) then
       do iextr=1,n_hist_cur-2
         h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+2)
       end do
     end if

     !--1--PROPAGATE----| C(t)--U[H(1/4dt)]-->C(t+dt/2)
     call propagate_orth(nstate,basis,time_step_cur/2.0_dp,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)  
     !--2--CALCULATE- H(t+dt/4) 
     call setup_hamiltonian_fock_cmplx( basis,                   &
                                        nstate,                  &
                                        itau,                    &
                                        time_cur,                &
                                        time_step_cur,           &
                                        occupation,              &
                                        c_matrix_cmplx,          &
                                        hamiltonian_kinetic,     &
                                        hamiltonian_nucleus,     &
                                        h_small_cmplx,           &
                                        s_matrix_sqrt_inv,       &
                                        dipole_basis,            &
                                        hamiltonian_fock_cmplx,  &
                                        ref_)
       
    if (n_hist_cur > 1) h_small_hist_cmplx(:,:,:,n_hist_cur-1)=h_small_cmplx 
    !--3--PROPAGATION----| 
     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)

     call setup_hamiltonian_fock_cmplx( basis,                   &
                                        nstate,                  &
                                        itau,                    &
                                        time_cur,                &
                                        time_step_cur,           &
                                        occupation,              &
                                        c_matrix_cmplx,          &
                                        hamiltonian_kinetic,     &
                                        hamiltonian_nucleus,     &
                                        h_small_cmplx,           &
                                        s_matrix_sqrt_inv,       &
                                        dipole_basis,            &
                                        hamiltonian_fock_cmplx,  &
                                        ref_)
       
     !--5--UPDATE----| 
     c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
     h_small_hist_cmplx(:,:,:,n_hist_cur)=h_small_cmplx
   end if

   ! ///////////////////////////////////

   if(pred_corr_cur=='PC2D' .OR. pred_corr_cur=='PC2E') then 
     h_small_hist_cmplx(:,:,:,n_hist_cur+1)=(0.0_dp,0.0_dp)
     !--0--EXTRAPOLATE---- 
     do iextr=1,n_hist_cur
       h_small_hist_cmplx(:,:,:,n_hist_cur+1)=h_small_hist_cmplx(:,:,:,n_hist_cur+1)+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
     end do
     !--1--PROPAGATE----| C(t)--U[H(1/4dt)]-->C(t+dt/2)
     call propagate_orth(nstate,basis,time_step_cur/2.0_dp,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_hist_cmplx(:,:,:,n_hist_cur:n_hist_cur+1),s_matrix_sqrt_inv,"ETRS")  

     ! n_hist_cur-2  <---  n_hist_cur; n_hist_cur-3 <-- n_hist_cur-1
     if(n_hist_cur > 2) then
       do iextr=1,n_hist_cur-2
         h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+2)
       end do
     end if

     !--2--CALCULATE- H(t+dt/4) 
     call setup_hamiltonian_fock_cmplx( basis,                   &
                                        nstate,                  &
                                        itau,                    &
                                        time_cur,                &
                                        time_step_cur,           &
                                        occupation,              &
                                        c_matrix_cmplx,          &
                                        hamiltonian_kinetic,     &
                                        hamiltonian_nucleus,     &
                                        h_small_cmplx,           &
                                        s_matrix_sqrt_inv,       &
                                        dipole_basis,            &
                                        hamiltonian_fock_cmplx,  &
                                        ref_)
       
    if (n_hist_cur > 1) h_small_hist_cmplx(:,:,:,n_hist_cur-1)=h_small_cmplx 
    !--3--PROPAGATION----| 
     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,"MAG2")

     call setup_hamiltonian_fock_cmplx( basis,                   &
                                        nstate,                  &
                                        itau,                    &
                                        time_cur,                &
                                        time_step_cur,           &
                                        occupation,              &
                                        c_matrix_cmplx,          &
                                        hamiltonian_kinetic,     &
                                        hamiltonian_nucleus,     &
                                        h_small_cmplx,           &
                                        s_matrix_sqrt_inv,       &
                                        dipole_basis,            &
                                        hamiltonian_fock_cmplx,  &
                                        ref_)
       
     !--5--UPDATE----| 
     c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
     h_small_hist_cmplx(:,:,:,n_hist_cur)=h_small_cmplx
   end if

   ! ///////////////////////////////////

   ! Iterative propagation with Lagrange interpolation
   ! ---------------|t-dt|-----------------|t|----------|t+dt/2)|----------|t+dt|
   !............|H(n_hist-1)|..........|H(n_hist)|....|H(n_hist+1)|.....|H(n_hist+2)|
   if(pred_corr_cur=='PC3' .OR. pred_corr_cur=='PC4') then
     hamiltonian_fock_cmplx=(0.0_dp,0.0_dp)
     h_small_hist_cmplx(:,:,:,n_hist_cur+2)=(0.0_dp,0.0_dp)
     !--1--EXTRAPOLATION WITH PREVIOUS STEPS---- 
     do iextr=1,n_hist_cur
       h_small_hist_cmplx(:,:,:,n_hist_cur+2)=h_small_hist_cmplx(:,:,:,n_hist_cur+2)+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
     end do

     !--2--LOCAL LINEAR INTERPOLATION----| hamiltonian_fock_cmplx in the 1/2 of the current time interval
     h_small_hist_cmplx(:,:,:,n_hist_cur+1)=0.5_dp*(h_small_hist_cmplx(:,:,:,n_hist_cur)+h_small_hist_cmplx(:,:,:,n_hist_cur+2))      

!     if ( is_iomaster .AND. mod(itau-1,mod_write)==0 ) then
!       write(name_iter_norm,"(3A,I4.4,A)") "./iter_norm/", TRIM(pred_corr), "_norm_itau_",itau,".dat"
!       open(newunit=file_iter_norm,file=name_iter_norm)
!     end if
     do i_iter=1,n_iter_cur
       h_small_cmplx(:,:,:)=h_small_hist_cmplx(:,:,:,n_hist_cur+1)
       !--3--PREDICTOR (propagation of C(0)-->C(1))
       call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)  
  
       !--4--CORRECTOR----| C(1)-->H(1)
       call setup_hamiltonian_fock_cmplx( basis,                   &
                                          nstate,                  &
                                          itau,                    &
                                          time_cur,                &
                                          time_step_cur,           &
                                          occupation,              &
                                          c_matrix_cmplx,          &
                                          hamiltonian_kinetic,     &
                                          hamiltonian_nucleus,     &
                                          h_small_hist_cmplx(:,:,:,n_hist_cur+2),           &
                                          s_matrix_sqrt_inv,       &
                                          dipole_basis,            &
                                          hamiltonian_fock_cmplx,  &
                                          ref_)

       c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
       !--2B--LOCAL LINEAR INTERPOLATION----| hamiltonian_fock_cmplx in the 1/2 of the current time interval
       h_small_hist_cmplx(:,:,:,n_hist_cur+1)=0.5_dp*(h_small_hist_cmplx(:,:,:,n_hist_cur)+h_small_hist_cmplx(:,:,:,n_hist_cur+2))      

       !**COMPARISON**
!       if( is_iomaster ) then      
!         if(mod(itau-1,mod_write)==0 ) then
!           write(file_iter_norm,*) i_iter, NORM2(ABS(h_small_hist_cmplx(:,:,:,n_hist_cur+1)-h_small_cmplx(:,:,:)))      
!         end if
!       end if
     end do ! i_iter
!     close(file_iter_norm)

     !--5--PROPAGATION----| C(0)---U[H(1/2)]--->C(1)
     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_hist_cmplx(:,:,:,n_hist_cur+1),s_matrix_sqrt_inv,prop_type_cur)
 
     !--6--UPDATE----|C(1)-->C(0)
     c_matrix_orth_cmplx(:,:,:)=c_matrix_orth_hist_cmplx(:,:,:,1) 
     do iextr=1,n_hist_cur-1
       h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+1)
     end do
     h_small_hist_cmplx(:,:,:,n_hist_cur)=h_small_hist_cmplx(:,:,:,n_hist_cur+2)
   end if

   ! ///////////////////////////////////

   ! Iterative ETRS - enforced time-reversal symmetry
   ! ---------------|t-dt|-----------------|t|--------------------|t+dt|
   !............|H(n_hist-1)|..........|H(n_hist)|.............|H(n_hist+1)|
   if(pred_corr_cur=='PC5' .OR. pred_corr_cur=='PC6') then
     hamiltonian_fock_cmplx=(0.0_dp,0.0_dp)
     h_small_hist_cmplx(:,:,:,n_hist_cur+1)=(0.0_dp,0.0_dp)
     !--1--EXTRAPOLATION WITH PREVIOUS STEPS---- 
     do iextr=1,n_hist_cur
       h_small_hist_cmplx(:,:,:,n_hist_cur+1)=h_small_hist_cmplx(:,:,:,n_hist_cur+1)+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
     end do

!     if( is_iomaster ) then
!       if(mod(itau-1,mod_write)==0 ) then
!         write(name_iter_norm,"(3A,I4.4,A)") "./iter_norm/", TRIM(pred_corr), "_norm_itau_",itau,".dat"
!         open(newunit=file_iter_norm,file=name_iter_norm)
!       end if
!     end if
     do i_iter=1,n_iter_cur
       c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
       h_small_cmplx(:,:,:)=h_small_hist_cmplx(:,:,:,n_hist_cur+1)
       !--3--PREDICTOR (propagation of C(0)-->C(1))
       call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_hist_cmplx(:,:,:,n_hist_cur:n_hist_cur+1),s_matrix_sqrt_inv,"ETRS")  
       !--4--CORRECTOR----| C(1)-->H(1)
       call setup_hamiltonian_fock_cmplx( basis,                   &
                                          nstate,                  &
                                          itau,                    &
                                          time_cur,                &
                                          time_step_cur,           &
                                          occupation,              &
                                          c_matrix_cmplx,          &
                                          hamiltonian_kinetic,     &
                                          hamiltonian_nucleus,     &
                                          h_small_hist_cmplx(:,:,:,n_hist_cur+1) ,    &
                                          s_matrix_sqrt_inv,       &
                                          dipole_basis,            &
                                          hamiltonian_fock_cmplx,  &
                                          ref_)

!       c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)

       !**COMPARISON**
!       if( is_iomaster ) then
!         if(mod(itau-1,mod_write)==0 ) then
!           write(file_iter_norm,*) i_iter, NORM2(ABS(h_small_hist_cmplx(:,:,:,n_hist_cur+1)-h_small_cmplx(:,:,:)))      
!         end if
!       end if
     end do ! i_iter
!     close(file_iter_norm)

!     !--5--PROPAGATION----| C(0)---U[H(1/2)]--->C(!1)
!     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_hist_cmplx(:,:,:,n_hist_cur:n_hist_cur+1),s_matrix_sqrt_inv,"ETRS")
 
     !--6--UPDATE----|
!     c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
     c_matrix_orth_cmplx(:,:,:)=c_matrix_orth_hist_cmplx(:,:,:,1)
     do iextr=1,n_hist_cur
       h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+1)
     end do
   end if

   ! ///////////////////////////////////

   ! ETRS proposed by Xavier Andrade
   ! ---------------|t-dt|-----------------|t|--------------------|t+dt|
   !............|---------|...............|H(1)|..................|H(2)|
   if(pred_corr_cur=='PC7' ) then
     hamiltonian_fock_cmplx=(0.0_dp,0.0_dp)
     h_small_hist_cmplx(:,:,:,2)=(0.0_dp,0.0_dp)

     call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_hist_cmplx(:,:,:,1),s_matrix_sqrt_inv,"MAG2")  

     do i_iter=1,n_iter_cur
       call setup_hamiltonian_fock_cmplx( basis,                   &
                                          nstate,                  &
                                          itau,                    &
                                          time_cur,                &
                                          time_step_cur,           &
                                          occupation,              &
                                          c_matrix_cmplx,          &
                                          hamiltonian_kinetic,     &
                                          hamiltonian_nucleus,     &
                                          h_small_hist_cmplx(:,:,:,2),   &
                                          s_matrix_sqrt_inv,       &
                                          dipole_basis,            &
                                          hamiltonian_fock_cmplx,  &
                                          ref_)
  
       if(i_iter/=n_iter_cur) then
         c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
         call propagate_orth(nstate,basis,time_step_cur,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_hist_cmplx(:,:,:,1:2),s_matrix_sqrt_inv,"ETRS")  
       end if
     end do
 
     !--6--UPDATE----|C(1)-->C(0)
     c_matrix_orth_cmplx(:,:,:)=c_matrix_orth_hist_cmplx(:,:,:,1)
     h_small_hist_cmplx(:,:,:,1)=h_small_hist_cmplx(:,:,:,2)
   end if

    
!   call check_identity_cmplx(basis%nbf,basis%nbf,MATMUL(MATMUL(TRANSPOSE(CONJG(c_matrix_cmplx(:,:,nspin))),s_matrix(:,:)), c_matrix_cmplx(:,:,nspin) ),is_identity_) 
!   if(.NOT. is_identity_) then
!     write(stdout,*) "C**H*S*C is not identity at itau= ", itau
!   end if

!   call print_2d_matrix("check",MATMUL(MATMUL(TRANSPOSE(CONJG(c_matrix_cmplx(:,:,nspin))),s_matrix(:,:)), c_matrix_cmplx(:,:,nspin) )   ,nocc,nocc,stdout,4)
!   call print_2d_matrix("check",s_matrix,basis%nbf,basis%nbf,stdout,4)
!   call print_2d_matrix("c_matrix_orth_cmplx",c_matrix_orth_cmplx(:,:,1),basis%nbf,nocc,stdout,4)
!   call print_2d_matrix("c_matrix_cmplx",c_matrix_cmplx(:,:,1),basis%nbf,nocc,stdout,4)

   if( is_iomaster .AND. ABS(time_cur / (write_step)- NINT(time_cur / (write_step))) < 1.0e-7 ) then

     if(excit_type%is_projectile) call output_projectile_position()

     call setup_density_matrix_cmplx(basis%nbf,nstate,nocc,c_matrix_cmplx,occupation,p_matrix_cmplx)
     en%kin = real(SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
     en%nuc = real(SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
     en%tot = en%nuc + en%kin + en%nuc_nuc + en%hart + en%exx_hyb + en%xc + en%excit

     if(ref_) then
       if( print_cube_rho_tddft_ ) call plot_cube_wfn_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,iwrite_step)
       if(excit_type%is_light) write(file_excit_field,*) time_cur, REAL(m_excit_field_dir)
     end if
     if(calc_p_matrix_error_) then
       p_matrix_time_cmplx(:,:,:,iwrite_step)=p_matrix_cmplx(:,:,:)
     end if
     if(excit_type%is_light) then
       call static_dipole_fast_cmplx(basis,p_matrix_cmplx,dipole_basis,dipole)
       dipole_time(iwrite_step,:)=dipole(:)
       write(file_dipole_time,*) time_cur, dipole(:) * au_debye
     end if
     write(file_time_data,"(F9.4,8(2x,es16.8E3),2x,2(2x,F7.2))") &
        time_cur, en%tot, en%nuc_nuc, en%nuc, en%kin, en%hart, en%exx_hyb, en%xc, en%excit, matrix_trace_cmplx(MATMUL(p_matrix_cmplx(:,:,1),s_matrix(:,:)))
     write(stdout,*)
     write(stdout,'(1x,a31,1x,f19.10)')  'RT-TDDFT Simulation time  (au):', time_cur
     write(stdout,'(1x,a31,1x,f19.10)')  'RT-TDDFT Total Energy     (Ha):', en%tot
     if(excit_type%is_light) write(stdout,'(1x,a31,1x,3f19.10)') 'RT-TDDFT Dipole Moment      (D):', dipole(:) * au_debye

     iwrite_step = iwrite_step + 1

   end if

!Q matrix (occupations) calculation
!   if(ref_) then
!     do ispin=1,nspin
!       q_matrix_cmplx(:,:,ispin)=MATMUL(MATMUL(CONJG(TRANSPOSE(c_matrix_0_cmplx(:,:,ispin))),s_matrix),c_matrix_cmplx(:,:,ispin))
!       if( is_iomaster ) then
!         if( time_cur / write_step - NINT(time_cur / write_step) < 1.0e-7_dp  ) then 
!           write(file_q_matrix_ii(ispin),"(F9.4)",advance='no') time_cur
!           do i_elem_q_mat=min_elem_q_mat, max_elem_q_mat
!             write(file_q_matrix_ii(ispin),"(2x,F9.4)",advance='no') ABS((q_matrix_cmplx(i_elem_q_mat,i_elem_q_mat,ispin)))**2
!           end do
!           write(file_q_matrix_ii(ispin),*)
!         end if
!       end if
!     end do
!   end do
!   write(stdout,"(a,F9.4)") "time_cur  ", time_cur
!   call print_square_2d_matrix_cmplx("q_matrix_cmplx",q_matrix_cmplx,nstate,stdout,4)
!--TIMING
   if(itau==3) then
     call stop_clock(timing_tddft_one_iter)
     time_one_iter=timing(timing_tddft_one_iter)
     write(stdout,"(1x,a30,2x,es14.6,1x,a)") "Time of one iteration is", time_one_iter,"s"
     write(stdout,"(1x,a30,2x,3(f12.2,1x,a))") "Estimated calculation time is", time_one_iter*ntau, "s = ", time_one_iter*ntau/60, "min = ", time_one_iter*ntau/3600, "hrs"
     call flush(stdout)
   end if

   if( print_tddft_restart_ .AND. mod(itau,100)==0 ) then
     call write_restart_tddft(nstate,time_cur,c_matrix_orth_cmplx)
   end if

  time_cur = time_min + itau*time_step_cur
  itau = itau + 1 
!---
 end do
 in_tddft_loop=.false.
!********end time loop*******************
 
 if(print_tddft_restart_) then
   !time_cur-time_step_cur to be consistent with the actual last moment of the simulation
   call write_restart_tddft(nstate,time_cur-time_step_cur,c_matrix_orth_cmplx)
 end if

 if( is_iomaster) then 
   close(file_time_data)
   if( excit_type%is_light) close(file_dipole_time)
   if( ref_ .AND. excit_type%is_light ) close(file_excit_field)
 end if
! if( ref_ ) then
!   do ispin=1,nspin
!     close(file_q_matrix_ii(ispin))
!   end do
! end if
 deallocate(c_matrix_cmplx)
 deallocate(c_matrix_0_cmplx)
 deallocate(c_matrix_orth_cmplx)
! deallocate(q_matrix_cmplx)
 deallocate(hamiltonian_fock_cmplx)
 deallocate(p_matrix_cmplx)
 deallocate(h_small_cmplx)

 if(pred_corr_cur /= 'PC0') then
   deallocate(h_small_hist_cmplx)
   deallocate(c_matrix_orth_hist_cmplx)
   deallocate(m_nods,extrap_coefs)
 end if

end subroutine tddft_time_loop

!==========================================
subroutine check_identity_cmplx(n,m,mat_cmplx,ans)
 implicit none
 integer,intent(in)     :: n, m
 complex(dp),intent(in) :: mat_cmplx(n,m)
 logical,intent(inout)  :: ans               
!=====
 integer   ::  i,j
 real(dp),parameter :: tol=1.0e-9_dp
!=====

 ans=.true.
 do i=1,n
   do j=1,m
     if(i==j) then
       if(ABS(mat_cmplx(i,j)-1.0_dp)>tol) then
         write(stdout,*) "M(i,i)/=1 for: ",i,j,mat_cmplx(i,j)
         ans=.false.
       end if
     else
       if(ABS(mat_cmplx(i,j))>tol) then
         write(stdout,*) "M(i,j)/=0 for: ",i,j,mat_cmplx(i,j)
         ans=.false.
       end if
     end if
   end do 
 end do

end subroutine check_identity_cmplx
!==========================================


!==========================================
subroutine write_restart_tddft(nstate,time_cur,c_matrix_orth_cmplx)
 use m_definitions
 implicit none
 integer,intent(in)         :: nstate
 complex(dp),intent(in)     :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 real(dp),intent(in)        :: time_cur
!===
 integer                    :: restartfile
 integer                    :: ibf, istate,ispin
!=====

 if( .NOT. is_iomaster) return
 
 write(stdout,'(/,a)') ' Writing a RESTART_TDDFT file'

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',action='write')
 ! Atomic structure
 write(restartfile) natom
 write(restartfile) zatom(1:natom)
 write(restartfile) xatom(:,1:natom)
 ! nocc
 write(restartfile) nocc
 ! Nstate
 write(restartfile) nstate
 ! nspin
 write(restartfile) nspin
 ! current time
 write(restartfile) time_cur
 ! Complex wavefunction coefficients C
 do ispin=1,nspin
   do istate=1,nocc
     write(restartfile) c_matrix_orth_cmplx(:,istate,ispin)
   enddo
 enddo

 close(restartfile)

end subroutine write_restart_tddft

!==========================================
subroutine read_restart_tddft(nstate,time_min,c_matrix_orth_cmplx,restart_is_correct)
 use m_definitions
 implicit none
 integer,intent(in)         :: nstate
 complex(dp),intent(inout)  :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 logical,intent(out)        :: restart_is_correct
 real(dp),intent(inout)     :: time_min
!===
 logical                    :: file_exists,same_geometry
 integer                    :: restartfile
 integer                    :: ibf, istate,ispin
 integer                    :: natom_read
 real(dp),allocatable       :: zatom_read(:),x_read(:,:)
 real(dp)                   :: nbf_read
 integer                    :: nstate_read, nspin_read,nocc_read
!=====

 if( .NOT. is_iomaster) return

 write(stdout,'(/,a)') ' Reading a RESTART_TDDFT file'

 restart_is_correct=.true.

 inquire(file='RESTART_TDDFT',exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,a)') ' No RESTART file found'
   restart_is_correct=.false.
   return
 endif

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',status='old',action='read')

 ! Atomic structure
 read(restartfile) natom_read
 allocate(zatom_read(natom_read),x_read(3,natom_read))
 read(restartfile) zatom_read(1:natom_read)
 read(restartfile) x_read(:,1:natom_read)
 if( natom_read /= natom  &
  .OR. ANY( ABS( zatom_read(1:MIN(natom_read,natom)) - zatom(1:MIN(natom_read,natom)) ) > 1.0e-5_dp ) &
  .OR. ANY( ABS(   x_read(:,1:MIN(natom_read,natom)) - xatom(:,1:MIN(natom_read,natom))   ) > 1.0e-5_dp ) ) then
   same_geometry = .FALSE.
   call issue_warning('RESTART_TDDFT file: Geometry has changed')
 else
   same_geometry = .TRUE.
 endif
 deallocate(zatom_read,x_read)

 ! nocc
 read(restartfile) nocc_read
 if(nocc /= nocc_read) then
   call issue_warning('RESTART_TDDFT file: nocc is not the same, restart file will not be used')
   restart_is_correct=.false.
   return
 end if

 ! Nstate
 read(restartfile) nstate_read
 if(nstate /= nstate_read) then
   call issue_warning('RESTART_TDDFT file: nstate is not the same, restart file will not be used')
   restart_is_correct=.false.
   return
 end if

 ! nspin
 read(restartfile) nspin_read
 if(nspin /= nspin_read) then
   call issue_warning('RESTART_TDDFT file: nspin is not the same, restart file will not be used')
   restart_is_correct=.false.
   return
 end if
 
 ! current time
 read(restartfile) time_min
 write(stdout,"(1x,a,f7.3)") "time_min= ", time_min

 ! Complex wavefunction coefficients C
 do ispin=1,nspin
   do istate=1,nocc
     read(restartfile) c_matrix_orth_cmplx(:,istate,ispin)
   enddo
 enddo

 close(restartfile)

end subroutine read_restart_tddft



!==========================================
subroutine get_extrap_coefs_lagr(m_nods,x_pred,extrap_coefs,n_hist_cur)
 use m_definitions 
 implicit none
 integer, intent(in)          :: n_hist_cur
 real(dp),intent(in)          :: m_nods(n_hist_cur)
 real(dp),intent(in)          :: x_pred
 real(dp),intent(inout)       :: extrap_coefs(n_hist_cur)
!=====
 integer      :: i,j
!=====

 extrap_coefs(:)=1.0_dp
 do i=1,n_hist_cur
   do j=1,n_hist_cur
     if(i==j) cycle
     extrap_coefs(i)=extrap_coefs(i)*(x_pred-m_nods(j))/(m_nods(i)-m_nods(j))
   end do
 end do
end subroutine get_extrap_coefs_lagr

!==========================================
subroutine get_extrap_coefs_aspc(extrap_coefs,n_hist_cur)
 use m_definitions
 implicit none
 integer, intent(in)          :: n_hist_cur
 real(dp),intent(inout)       :: extrap_coefs(n_hist_cur)
!=====

 if(n_hist_cur==1) then
   extrap_coefs(1)=1.0_dp
 end if


 if(n_hist_cur==2) then
   extrap_coefs(1)=-1.0_dp
   extrap_coefs(2)=2.0_dp
 end if

 if(n_hist_cur==3) then
   extrap_coefs(1)=0.5_dp
   extrap_coefs(2)=-2.0_dp
   extrap_coefs(3)=2.5_dp
 end if

 if(n_hist_cur==4) then
   extrap_coefs(1)=-0.2_dp
   extrap_coefs(2)=1.2_dp
   extrap_coefs(3)=-2.8_dp
   extrap_coefs(4)=2.8_dp
 end if

 if(n_hist_cur==5) then
   extrap_coefs(1)=1.0_dp/14.0_dp
   extrap_coefs(2)=-4.0_dp/7.0_dp
   extrap_coefs(3)=27.0_dp/14.0_dp
   extrap_coefs(4)=-24.0_dp/7.0_dp
   extrap_coefs(5)=3.0_dp
 end if

 if(n_hist_cur==6) then
   extrap_coefs(1)=-1.0_dp/42.0_dp
   extrap_coefs(2)=5.0_dp/21.0_dp
   extrap_coefs(3)=-22.0_dp/21.0_dp
   extrap_coefs(4)=55.0_dp/21.0_dp
   extrap_coefs(5)=-55.0_dp/14.0_dp
   extrap_coefs(6)=22.0_dp/7.0_dp
 end if


end subroutine get_extrap_coefs_aspc

!==========================================
subroutine propagate_orth_ham_1(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type_cur)
 use m_timing
 implicit none
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: time_step_cur
 complex(dp),intent(inout)   :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 complex(dp), intent(inout)  :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(in)      :: h_small_cmplx(nstate,nstate,nspin)
 real(dp),intent(in)         :: s_matrix_sqrt_inv(basis%nbf,nstate)
 character(len=4),intent(in) :: prop_type_cur
!=====
 integer                    :: ispin
 integer                    :: ibf
!==variables for the MAG2 propagator
 complex(dp),allocatable    :: a_matrix_orth_cmplx(:,:)
 real(dp),allocatable       :: energies_inst(:)
 complex(dp),allocatable    :: propagator_eigen(:,:)
!==variables for the CN propagator
 complex(dp),allocatable    :: l_matrix_cmplx(:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al, 
 complex(dp),allocatable    :: b_matrix_cmplx(:,:) ! TDDFT Book, Springer (2006), !p205
!=====

 call start_clock(timing_tddft_propagation)

! a_matrix_cmplx(:,1:nstate) = MATMUL( s_matrix_sqrt_inv(:,:) , a_matrix_cmplx(:,:) )

 do ispin =1, nspin
!   h_small_cmplx(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
 !                   MATMUL( hamiltonian_fock_cmplx(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )
   select case (prop_type_cur)
   case('CN')
     allocate(l_matrix_cmplx(nstate,nstate))
     allocate(b_matrix_cmplx(nstate,nstate))
     l_matrix_cmplx(:,:)= im * time_step_cur / 2.0_dp * h_small_cmplx(:,:,ispin)
     b_matrix_cmplx(:,:)=-l_matrix_cmplx(:,:)
     do ibf=1,nstate
       b_matrix_cmplx(ibf,ibf)=b_matrix_cmplx(ibf,ibf)+1.0_dp
       l_matrix_cmplx(ibf,ibf)=l_matrix_cmplx(ibf,ibf)+1.0_dp
     end do
     call invert(nstate , l_matrix_cmplx(:,:))
  
     c_matrix_orth_cmplx(:,:,ispin) = MATMUL( l_matrix_cmplx(:,:),MATMUL( b_matrix_cmplx(:,:),c_matrix_orth_cmplx(:,:,ispin) ) )
     deallocate(l_matrix_cmplx)
     deallocate(b_matrix_cmplx)
   case('MAG2')
     allocate(propagator_eigen(nstate,nstate))
     allocate(a_matrix_orth_cmplx(nstate,nstate))
     allocate(energies_inst(nstate))
     call diagonalize(nstate,h_small_cmplx(:,:,ispin),energies_inst(:),a_matrix_orth_cmplx)
     propagator_eigen(:,:) = ( 0.0_dp , 0.0_dp ) 
    
     do ibf=1,nstate
       propagator_eigen(ibf,ibf) = exp(-im*time_step_cur*energies_inst(ibf))
     end do

     c_matrix_orth_cmplx(:,:,ispin) = MATMUL( MATMUL( MATMUL( a_matrix_orth_cmplx(:,:), propagator_eigen(:,:)  ) , &
             conjg(transpose(a_matrix_orth_cmplx(:,:)))  ), c_matrix_orth_cmplx(:,:,ispin) )
     deallocate(propagator_eigen)
     deallocate(a_matrix_orth_cmplx)
     deallocate(energies_inst)
   case default
     call die('Invalid choice for the propagation algorithm. Change prop_type_cur value in the input file')
   end select
   c_matrix_cmplx(:,:,ispin) = MATMUL( s_matrix_sqrt_inv(:,:) , c_matrix_orth_cmplx(:,:,ispin) )
 end do

 call stop_clock(timing_tddft_propagation)


end subroutine propagate_orth_ham_1

!==========================================
subroutine propagate_orth_ham_2(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_hist2_cmplx,s_matrix_sqrt_inv,prop_type_cur)
 use m_timing
 implicit none
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: time_step_cur
 complex(dp),intent(inout)   :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 complex(dp), intent(inout)  :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(in)      :: h_small_hist2_cmplx(nstate,nstate,nspin,2)
 real(dp),intent(in)         :: s_matrix_sqrt_inv(basis%nbf,nstate)
 character(len=4),intent(in) :: prop_type_cur
!=====
 integer             :: ispin,iham
 integer             :: ibf
 complex(dp)         :: a_matrix_orth_cmplx(nstate,nstate,2)
 real(dp)            :: energies_inst(nstate)
 complex(dp)         :: propagator_eigen(nstate,nstate,2)
!=====

 call start_clock(timing_tddft_propagation)
! a_matrix_cmplx(:,1:nstate) = MATMUL( s_matrix_sqrt_inv(:,:) , a_matrix_cmplx(:,:) )

 do ispin =1, nspin
   select case (prop_type_cur)
   case('ETRS')
     do iham=1,2
       call diagonalize(nstate,h_small_hist2_cmplx(:,:,ispin,iham),energies_inst(:),a_matrix_orth_cmplx(:,:,iham))
       propagator_eigen(:,:,iham) = ( 0.0_dp , 0.0_dp ) 
       do ibf=1,nstate
         propagator_eigen(ibf,ibf,iham) = exp(-im*time_step_cur/2.d0*energies_inst(ibf))
       end do
     end do
     c_matrix_orth_cmplx(:,:,ispin) = & 
         MATMUL(MATMUL(MATMUL(MATMUL( MATMUL( MATMUL( a_matrix_orth_cmplx(:,:,2), propagator_eigen(:,:,2)  ) , conjg(transpose(a_matrix_orth_cmplx(:,:,2)))  ),  & 
                            a_matrix_orth_cmplx(:,:,1)), propagator_eigen(:,:,1)), conjg(transpose(a_matrix_orth_cmplx(:,:,1))) ), c_matrix_orth_cmplx(:,:,ispin) )
   case default
     call die('Invalid choice of the propagation algorithm for the given PC scheme. Change prop_type value in the input file')
   end select
    c_matrix_cmplx(:,:,ispin) = MATMUL( s_matrix_sqrt_inv(:,:) , c_matrix_orth_cmplx(:,:,ispin) )
 end do

call stop_clock(timing_tddft_propagation)

end subroutine propagate_orth_ham_2

!==========================================
subroutine setup_hamiltonian_fock_cmplx( basis,                   &
                                         nstate,                  &   
                                         itau,                    & 
                                         time_cur,                &     
                                         time_step_cur,           &
                                         occupation,              &       
                                         c_matrix_cmplx,          &           
                                         hamiltonian_kinetic,     &                
                                         hamiltonian_nucleus,     &                
                                         h_small_cmplx,           &
                                         s_matrix_sqrt_inv,       &
                                         dipole_basis,            &         
                                         hamiltonian_fock_cmplx,  &
                                         ref_)                      

 use m_timing
 implicit none
!=====
 type(basis_set),intent(in)       :: basis
 logical,intent(in)               :: ref_
 integer,intent(in)               :: nstate
 integer,intent(in)               :: itau
 real(dp),intent(in)              :: time_cur
 real(dp),intent(in)              :: time_step_cur
 real(dp),intent(in)              :: occupation(nstate,nspin)
 real(dp),intent(in)              :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(inout)           :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),allocatable,intent(in)  :: dipole_basis(:,:,:)
 real(dp),intent(in)              :: s_matrix_sqrt_inv(basis%nbf,nstate)
 complex(dp),intent(in)           :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(inout)        :: hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin)
 complex(dp),intent(inout)        :: h_small_cmplx(nstate,nstate,nspin)
!=====
 logical        :: calc_excit_
 integer        :: ispin, idir
 real(dp)       :: excit_field(3)
 complex(dp)    :: p_matrix_cmplx(basis%nbf,basis%nbf,nspin)
!=====

 call start_clock(timing_tddft_hamiltonian_fock)

 call setup_density_matrix_cmplx(basis%nbf,nstate,nocc,c_matrix_cmplx,occupation,p_matrix_cmplx)
!--Hamiltonian - Hartree Exchange Correlation---
 call calculate_hamiltonian_hxc_ri_cmplx(basis,                    &
                                         nstate,                   &
                                         nocc,                     &
                                         basis%nbf,                &
                                         basis%nbf,                &
                                         basis%nbf,                &
                                         nocc,                     &      
                                         occupation,               &       
                                         c_matrix_cmplx,           &          
                                         p_matrix_cmplx,           &          
                                         hamiltonian_fock_cmplx)                   
 do ispin=1, nspin                                                  
   en%excit=0.0_dp
   if(excit_type%is_light) then
     !--Hamiltonian - Excitation--
     excit_field=0.0_dp
     calc_excit_ = .false.
     calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'GAU' )
     calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'HSW'  .AND. abs(time_cur - excit_type%time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )  
     calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'STEP' .AND. abs(time_cur - excit_type%time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )
     calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'DEL'  .AND. abs(time_cur - excit_type%time0)<=time_step_cur ) 
     if(itau==0) calc_excit_=.false.
     if(excit_type%name == 'NO') calc_excit_=.false.   
     if ( calc_excit_ ) then
       call calculate_excit_field(time_cur,excit_field)
       if(ref_)  m_excit_field_dir=NORM2(excit_field(:))
       do idir=1,3
         hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) - dipole_basis(:,:,idir) * excit_field(idir)
         en%excit=en%excit+real(SUM(dipole_basis(:,:,idir)*excit_field(idir)*p_matrix_cmplx(:,:,ispin)),dp)
       end do     
     end if
   end if ! light excitation
   if(excit_type%is_projectile) then

     xatom(:,natom)=xatom_start(:,natom) + vel(:,natom)*time_cur
     call nucleus_nucleus_energy(en%nuc_nuc)
  
     !
     ! Nucleus-electron interaction
     if( parallel_ham ) then
       if( parallel_buffer ) then
         call setup_nucleus_buffer_sca(.false.,basis,basis%nbf,basis%nbf,hamiltonian_nucleus)
       else
         call setup_nucleus_sca(.false.,basis,basis%nbf,basis%nbf,hamiltonian_nucleus)
       endif
     else
       call setup_nucleus(.false.,basis,hamiltonian_nucleus)

       if( nelement_ecp > 0 ) then
         call setup_nucleus_ecp(.false.,basis,hamiltonian_nucleus)
       endif
     endif
     !-------------------------------
   end if
   hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
   h_small_cmplx(:,:,ispin) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
                   MATMUL( hamiltonian_fock_cmplx(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )
 end do ! spin loop

 call stop_clock(timing_tddft_hamiltonian_fock)

end subroutine setup_hamiltonian_fock_cmplx


!==========================================
subroutine calculate_excit_field(time_cur,excit_field)
 implicit none
 real(dp),intent(in)      :: time_cur ! time in au
 real(dp),intent(inout)   :: excit_field(3) ! electric field in 3 dimensions
!=====


 select case(excit_type%name)
 case('GAU') !Gaussian electic field
   excit_field(:) = excit_type%kappa * exp( -( time_cur-excit_type%time0 )**2 / 2.0_dp / excit_omega**2 ) * &
                  & excit_dir_norm(:)
 case('HSW') !Hann sine window
   excit_field(:) = excit_type%kappa * sin( pi / excit_omega * ( time_cur - excit_type%time0  ) )**2 * excit_dir_norm(:) 
 case('DEL') ! Delta excitation
   excit_field(:) = excit_type%kappa * excit_dir_norm(:) 
 case('STEP') ! Step excitation
   excit_field(:) = excit_type%kappa * excit_dir_norm(:) 
 case default
    call die('Invalid choice for the excitation type. Change excit_type value in the input file')
 end select
 
end subroutine calculate_excit_field

!=======================================
subroutine fill_unity(unity_matrix_cmplx,M)
 implicit none
 integer, intent(in) :: M
 complex(dp) :: unity_matrix_cmplx(M,M)
 integer :: i,j
 do i=1,M
   do j=1,M
     if (i == j) then
       unity_matrix_cmplx(i,j)=1.0_dp
     else
      unity_matrix_cmplx(i,j)=0.0_dp
     end if
   end do
 end do
end subroutine

!=======================================
subroutine print_2d_matrix_cmplx(desc,matrix_cmplx,size_n,size_m,write_unit,prec)
 implicit none
 integer, intent(in)      :: prec ! precision
 integer, intent(in)      :: size_n,size_m, write_unit
 complex(dp),intent(in)  :: matrix_cmplx(size_m,size_m)
 character(*),intent(in)  :: desc
!=====
 character(100)  :: write_format1, write_format2, write_form
 integer            :: ivar,beg

! beg=4
 beg=3
 write(write_format1,*) '(',size_m," ('( ',F", prec+beg, ".", prec,"' ,',F", prec+beg, ".",prec,",' )  ') " ,')' ! (  1.01 ,  -0.03)  (  0.04 ,  0.10) 
 write(write_format2,*) '(',size_m," (F", prec+beg, ".", prec,"' +  i',F", prec+beg, ".",prec,",'  ') " ,')'   ! 1.01 +  i  -0.03    0.03 +  i  0.10
 write(write_unit,*) desc
 do ivar=1,size_n
   write(write_unit,write_format1) matrix_cmplx(ivar,:)
 end do

end subroutine print_2d_matrix_cmplx

!=======================================
subroutine print_2d_matrix_real(desc,matrix_real,size_n,size_m,write_unit,prec)
 implicit none
 integer, intent(in)      :: prec ! precision
 integer, intent(in)      :: size_n,size_m, write_unit
 real(dp),intent(in)      :: matrix_real(size_m,size_m)
 character(*),intent(in)  :: desc
!=====
 character(100)  :: write_format1
 integer            :: ivar

 write(write_format1,*) '(',size_m," (F", prec+4, ".", prec,') ' ,')' 
 write(write_unit,*) desc
 do ivar=1,size_n
   write(write_unit,write_format1) matrix_real(ivar,:)          
 end do
end subroutine print_2d_matrix_real


!==========================================
subroutine propagate_non_orth(nstate,basis,time_step_cur,c_matrix_cmplx,hamiltonian_fock_cmplx,s_matrix,s_matrix_sqrt_inv,prop_type_cur)
 implicit none
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: time_step_cur
 complex(dp),intent(inout)   :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(in)      :: hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)         :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)         :: s_matrix_sqrt_inv(basis%nbf,nstate)
 character(len=4),intent(in)            :: prop_type_cur
!=====
 integer                    :: ispin
!==variables for the MAG2 propagator
 integer        :: ibf
 complex(dp),allocatable    :: h_small_cmplx(:,:)
 complex(dp),allocatable    :: a_matrix_cmplx(:,:)
 complex(dp),allocatable    :: propagator_eigen(:,:)
 real(dp),allocatable       :: energies_inst(:)
!==variables for the CN propagator
 complex(dp),allocatable    :: l_matrix_cmplx(:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al, 
 complex(dp),allocatable    :: b_matrix_cmplx(:,:) ! TDDFT Book, Springer (2006), !p205
!=====

 allocate(h_small_cmplx(nstate,nstate))
 h_small_cmplx(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
                       MATMUL( hamiltonian_fock_cmplx(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )
 

 do ispin =1, nspin
   select case (prop_type_cur)
   case('CN')
     allocate(l_matrix_cmplx(basis%nbf,basis%nbf))
     allocate(b_matrix_cmplx(basis%nbf,basis%nbf))
     l_matrix_cmplx(:,:)= im * time_step_cur / 2.0_dp * MATMUL( s_matrix_inv,hamiltonian_fock_cmplx(:,:,ispin) )
     b_matrix_cmplx(:,:)=-l_matrix_cmplx(:,:)
     l_matrix_cmplx(:,:)= im * time_step_cur / 2.0_dp * h_small_cmplx(:,:)
     do ibf=1,basis%nbf
       b_matrix_cmplx(ibf,ibf)=b_matrix_cmplx(ibf,ibf)+1.0_dp
       l_matrix_cmplx(ibf,ibf)=l_matrix_cmplx(ibf,ibf)+1.0_dp
     end do
     call invert(basis%nbf , l_matrix_cmplx(:,:))
  
     c_matrix_cmplx(:,:,ispin) = MATMUL( l_matrix_cmplx(:,:),MATMUL( b_matrix_cmplx(:,:),c_matrix_cmplx(:,:,ispin) ) )
    deallocate(l_matrix_cmplx)
    deallocate(b_matrix_cmplx)
   case('MAG2')
     allocate(energies_inst(nstate))
     allocate(propagator_eigen(basis%nbf,basis%nbf))
     allocate(a_matrix_cmplx(basis%nbf,nstate))
     propagator_eigen(:,:) = ( 0.0_dp , 0.0_dp ) 
     call diagonalize(nstate,h_small_cmplx(:,:),energies_inst(:),a_matrix_cmplx)
     a_matrix_cmplx(:,1:nstate) = MATMUL( s_matrix_sqrt_inv(:,:) , a_matrix_cmplx(:,:) )
    
     do ibf=1,basis%nbf
       propagator_eigen(ibf,ibf) = exp(-im*time_step_cur*energies_inst(ibf))
     end do
     !remove traspose
     c_matrix_cmplx(:,:,ispin) = MATMUL( MATMUL( MATMUL( MATMUL( a_matrix_cmplx(:,:), propagator_eigen(:,:)  ) , &
             conjg(transpose(a_matrix_cmplx(:,:)))  ), s_matrix(:,:) ), c_matrix_cmplx(:,:,ispin) )
     deallocate(energies_inst)
     deallocate(propagator_eigen)
     deallocate(a_matrix_cmplx)
   case default
     call die('Invalid choice of the propagation algorithm. Change prop_type value in the input file')
   end select
 end do


end subroutine propagate_non_orth

!=========================================================================
end module m_tddft_propagator
!=========================================================================


