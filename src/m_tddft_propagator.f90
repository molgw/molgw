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
 use m_tddft_variables
 use m_timing

 interface propagate_orth
  module procedure propagate_orth_ham_1
  module procedure propagate_orth_ham_2
 end interface propagate_orth

 integer,private                    :: nocc
 real(dp),private                   :: dipole(3)
 real(dp),private                   :: time_read
 real(dp),allocatable,private       :: xatom_start(:,:)
 complex(dp),private                :: excit_field_norm
 ! hamiltonian extrapolation variables
 real(dp),allocatable,private       :: extrap_coefs(:)
 complex(dp),allocatable,private    :: h_small_hist_cmplx(:,:,:,:)
 complex(dp),allocatable,private    :: c_matrix_orth_hist_cmplx(:,:,:,:)
 !-----------------------------------
 complex(dp),allocatable    :: q_matrix_cmplx(:,:,:)
 integer,private            :: ntau


contains


!=========================================================================
subroutine calculate_propagation(basis,occupation,c_matrix)
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(in)        :: occupation(:,:)
!=====
 integer,parameter          :: BATCH_SIZE = 128
 integer                    :: fixed_atom_list(natom-nprojectile)
 integer                    :: ispin 
 integer                    :: istate,nstate_tmp
 integer                    :: nwrite_step 
 real(dp)                   :: time_min
 real(dp),allocatable       :: dipole_basis(:,:,:)
 real(dp),allocatable       :: s_matrix(:,:)
 real(dp),allocatable       :: s_matrix_sqrt_inv(:,:)
 real(dp),allocatable       :: hamiltonian_kinetic(:,:)
 real(dp),allocatable       :: hamiltonian_nucleus(:,:)
!=====initial values
 integer                    :: nstate
 real(dp),allocatable       :: energies_inst(:)
 complex(dp),allocatable    :: c_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: c_matrix_orth_cmplx(:,:,:)
 complex(dp),allocatable    :: hamiltonian_fock_cmplx(:,:,:)
 complex(dp),allocatable    :: h_small_cmplx(:,:,:)
!=====TDDFT loop variables=============================
 integer                    :: iatom
 integer                    :: itau
 integer                    :: iwrite_step
 integer                    :: file_time_data,file_excit_field
 integer                    :: file_dipole_time
 real(dp)                   :: time_cur,time_one_iter
 complex(dp),allocatable    :: p_matrix_cmplx(:,:,:)
 logical                    :: is_identity_ ! keep this varibale
!==cube_diff varibales===
 real(dp),allocatable       :: cube_density_start(:,:,:,:)
 integer                    :: nx,ny,nz,unit_cube_diff
 logical                    :: file_exists
!==qmatrix==
 integer,allocatable        :: istate_cut(:)
 integer                    :: file_q_matrix(2)
 integer                    :: iocc
 complex(dp),allocatable    :: c_matrix_orth_start_complete_cmplx(:,:,:)
!=====

 call issue_warning("calculate_propagation: not available in this version. Contact the developing team.")

end subroutine calculate_propagation


!==========================================
subroutine echo_tddft_variables()
 implicit none

 write(stdout,'(/,1x,a)') 'The most important variables of this section:'
 write(stdout,'(2x,a32,2x,es16.8)') 'Simulation time: time_sim',time_sim
 write(stdout,'(2x,a32,2x,es16.8)') 'Time step: time_step',time_step
 write(stdout,'(2x,a32,2x,i8)') 'Number of iterations: ntau',NINT((time_sim)/time_step)
 write(stdout,'(2x,a32,6x,a)')      'Predictor-corrector: pred_corr',pred_corr
 write(stdout,'(2x,a32,6x,a)')      'Propagator: prop_type',prop_type
 write(stdout,'(2x,a32,2x,i8)')     'Number of occupied states: nocc',nocc
 write(stdout,'(2x,a32,2x,i8)')     'Historic of Hamiltonian: n_hist',n_hist
 write(stdout,*) 

end subroutine echo_tddft_variables

!==========================================
subroutine output_timing_one_iter()
 implicit none
 real(dp)           :: time_one_iter
!=====
 
  time_one_iter=timing(timing_tddft_one_iter)
  write(stdout,'(/,1x,a)') '**********************************'
  write(stdout,"(1x,a30,2x,es14.6,1x,a)") "Time of one iteration is", time_one_iter,"s"
  write(stdout,"(1x,a30,2x,3(f12.2,1x,a))") "Estimated calculation time is", time_one_iter*ntau, "s  = ", &
                                            time_one_iter*ntau/60.0_dp, &
                                            "min  = ", time_one_iter*ntau/3600.0_dp, "hrs"
  write(stdout,'(1x,a)') '**********************************'
  flush(stdout)

end subroutine output_timing_one_iter


!==========================================
subroutine initialize_extrap_coefs(c_matrix_orth_cmplx,h_small_cmplx)
 implicit none
 complex(dp),intent(in)    :: c_matrix_orth_cmplx(:,:,:)
 complex(dp),intent(in)    :: h_small_cmplx(:,:,:)
!=====
 integer               :: iextr,ham_hist_dim,nstate
 real(dp)              :: x_pred
 real(dp),allocatable  :: m_nodes(:)
!=====

 nstate = SIZE(c_matrix_orth_cmplx,DIM=1)

 allocate(m_nodes(n_hist),extrap_coefs(n_hist))

 select case (pred_corr)
 case('PC1') 
   ham_hist_dim=2

 case('PC2B') 
   ham_hist_dim=n_hist
   do iextr=1,n_hist
     m_nodes(iextr)=(iextr-1.0_dp)*0.5_dp
   end do
   x_pred=(n_hist-1.0_dp)*0.5_dp+0.25_dp
   call get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist)

 case('PC3','PC4') 
   ham_hist_dim=n_hist+2
   do iextr=1,n_hist
     m_nodes(iextr)=iextr-1.0_dp
   end do
   x_pred=n_hist
   if(pred_corr=='PC3') call get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist)
   if(pred_corr=='PC4') call get_extrap_coefs_aspc(extrap_coefs,n_hist)

 case('PC5','PC6') 
   ham_hist_dim=n_hist+1
   do iextr=1,n_hist
     m_nodes(iextr)=iextr-1.0_dp
   end do
   x_pred=n_hist
   if(pred_corr=='PC5') call get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist)
   if(pred_corr=='PC6') call get_extrap_coefs_aspc(extrap_coefs,n_hist)

 case('PC7' ) 
   ham_hist_dim=2

 end select

 if(pred_corr /= 'PC0' ) then
   call clean_allocate('h_small_hist_cmplx for TDDFT',h_small_hist_cmplx,nstate,nstate,nspin,ham_hist_dim)
   call clean_allocate('c_matrix_orth_hist_cmplx for TDDFT',c_matrix_orth_hist_cmplx,nstate,nocc,nspin,1)
   do iextr=1,ham_hist_dim
     h_small_hist_cmplx(:,:,:,iextr)=h_small_cmplx(:,:,:)
   end do
   c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
 end if

 deallocate(m_nodes)

end subroutine initialize_extrap_coefs

!==========================================
subroutine print_tddft_values(time_cur,file_time_data,file_dipole_time,file_excit_field,itau)
 implicit none
 integer,intent(in)    :: file_time_data,file_dipole_time,file_excit_field
 real(dp),intent(in)   :: time_cur
 integer,intent(in)    :: itau 
!=====

 if( .NOT. is_iomaster ) return

 write(stdout,'(/,1x,a)')    '==================================================================================================='
 write(stdout,'(1x,a,i8,a)') '===================== RT-TDDFT values for the iteration  ',itau,' ================================='
 write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Simulation time (au):', time_cur
 write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Total Energy    (Ha):', en%tot

 select case(excit_type%form)
 case(EXCIT_PROJECTILE)
   write(file_time_data,"(F9.4,8(2x,es16.8E3))") &
      time_cur, en%tot, xatom(3,natom), en%nuc_nuc, en%nuc, en%kin, en%hart, en%exx_hyb, en%xc
   call output_projectile_position()

 case(EXCIT_LIGHT)
   write(file_time_data,"(F9.4,8(2x,es16.8E3))") &
    time_cur, en%tot, en%nuc_nuc, en%nuc, en%kin, en%hart, en%exx_hyb, en%xc, en%excit
   write(file_dipole_time,'(4f19.10)') time_cur, dipole(:) * au_debye
   write(file_excit_field,'(2f19.10)') time_cur, REAL(excit_field_norm)
   write(stdout,'(a31,1x,3f19.10)') 'RT-TDDFT Dipole Moment    (D):', dipole(:) * au_debye
 end select

 write(stdout,'(1x,a,/)')      '==================================================================================================='

end subroutine print_tddft_values


!==========================================
subroutine initialize_files(file_time_data,file_dipole_time,file_excit_field)
 implicit none
 integer,intent(inout)    :: file_time_data,file_excit_field,file_dipole_time
!=====

 if( .NOT. is_iomaster ) return

 open(newunit=file_time_data,file="time_data.dat")

 if(excit_type%form==EXCIT_LIGHT) then

   open(newunit=file_dipole_time,file="dipole_time.dat")
   open(newunit=file_excit_field,file="excitation_time.dat")

   write(file_excit_field,"(A)") "# time(au)                      E_field_excit_dir(au)"

 end if

!---------------------------------
 select case(excit_type%form) 
 case(EXCIT_PROJECTILE)
   write(file_time_data,"(A)") "  # time(au)     e_total        z_projectile        enuc_nuc            enuc             ekin              ehart&
                             &           eexx_hyb            exc"

 case(EXCIT_LIGHT) 
   write(file_time_data,"(A)") " # time(au)     e_total             enuc_nuc             enuc            ekin               ehart            &
                             &eexx_hyb            exc             eexcit"

   write(file_dipole_time,"(A)") "# time(au)                      Dipole_x(D)               Dipole_y(D)               Dipole_z(D)" 
 end select

end subroutine initialize_files


!==========================================
subroutine initialize_q(nstate,nocc,nspin,c_matrix_orth_start_complete_cmplx,h_small_cmplx,istate_cut,file_q_matrix)
 implicit none
 integer,intent(in)                    :: nstate, nocc, nspin
 integer,allocatable,intent(out)       :: istate_cut(:)
 complex(dp),allocatable,intent(out)   :: c_matrix_orth_start_complete_cmplx(:,:,:)
 complex(dp),allocatable,intent(in)    :: h_small_cmplx(:,:,:)
 integer,intent(out)                   :: file_q_matrix(2)
!=====
 character(len=50)           :: name_file_q_matrix
 integer                    :: ispin
 real(dp),allocatable       :: energies_inst(:)
 logical                    :: file_exists
 integer                    :: file_q_matrix_param
!=====

 allocate(istate_cut(10))
 call clean_allocate('q_matrix for TDDFT',q_matrix_cmplx,nstate,nocc,nspin)
 call clean_allocate('c_matrix_orth_start for TDDFT',c_matrix_orth_start_complete_cmplx,nstate,nstate,nspin)
 allocate(energies_inst(nstate))
 do ispin=1, nspin
   call diagonalize(postscf_diago_flavor,h_small_cmplx(:,:,ispin),energies_inst,c_matrix_orth_start_complete_cmplx(:,:,ispin))
 end do
 deallocate(energies_inst)
 istate_cut(4)=nstate
 inquire(file='manual_q_matrix_param',exist=file_exists)
 if(file_exists) then
   open(newunit=file_q_matrix_param,file='manual_q_matrix_param',status='old')
   read(file_q_matrix_param,*) istate_cut(1), istate_cut(2), istate_cut(3)
   close(file_q_matrix_param)
 else
   istate_cut(1)=1
   istate_cut(2)=natom-1
   istate_cut(3)=natom+INT((natom-1)/2)
   call issue_warning('plot_rho_traj_bunch_contrib: manual_q_matrix_param file was not found')
 endif
 
 if( is_iomaster ) then
   do ispin=1,nspin
     write(name_file_q_matrix,"(a,i1,a)") "q_matrix_", ispin, ".dat"
     open(newunit=file_q_matrix(ispin),file=name_file_q_matrix)
   end do
 end if

end subroutine initialize_q

!==========================================
subroutine calc_q_matrix(occupation,c_matrix_orth_start_complete_cmplx,c_matrix_orth_cmplx,istate_cut,file_q_matrix,time_cur)
 implicit none
 real(dp),intent(in)      :: occupation(:,:)
 complex(dp),intent(in)   :: c_matrix_orth_start_complete_cmplx(:,:,:)
 complex(dp),intent(in)   :: c_matrix_orth_cmplx(:,:,:)
 integer,intent(in)       :: istate_cut(:)
 integer,intent(in)       :: file_q_matrix(:)
 real(dp),intent(in)      :: time_cur
!=====
 integer                  :: istate,iocc,ispin
 real(dp)                 :: q_occ(10)
!=====

 q_occ=0.0_dp

 do ispin=1,nspin
   q_matrix_cmplx(:,:,ispin)=MATMUL(CONJG(TRANSPOSE(c_matrix_orth_start_complete_cmplx(:,:,ispin))),c_matrix_orth_cmplx(:,:,ispin))

   do istate=istate_cut(1),istate_cut(2)
     do iocc=1,nocc
       q_occ(1)=q_occ(1)+ABS(q_matrix_cmplx(istate,iocc,ispin))**2*occupation(iocc,ispin)
     end do
   end do

   do istate=istate_cut(2)+1,istate_cut(3)
     do iocc=1,nocc
       q_occ(2)=q_occ(2)+ABS(q_matrix_cmplx(istate,iocc,ispin)**2)*occupation(iocc,ispin)
     end do
   end do

   do istate=istate_cut(3)+1,istate_cut(4)
     do iocc=1,nocc
       q_occ(3)=q_occ(3)+ABS(q_matrix_cmplx(istate,iocc,ispin))**2*occupation(iocc,ispin)
     end do
   end do

   write(file_q_matrix(ispin),"(F9.4,10(2x,es16.8E3))") time_cur, q_occ(:)
 end do

end subroutine calc_q_matrix

!==========================================
subroutine check_identity_cmplx(n,m,mat_cmplx,ans)
 implicit none
 integer,intent(in)     :: n, m
 complex(dp),intent(in) :: mat_cmplx(n,m)
 logical,intent(inout)  :: ans
!=====
 integer   ::  imat,jmat
 real(dp),parameter :: tol=1.0e-9_dp
!=====

 ans=.TRUE.
 do imat=1,n
   do jmat=1,m
     if(imat==jmat) then
       if(ABS(mat_cmplx(imat,jmat)-1.0_dp)>tol) then
         write(stdout,*) "M(imat,imat)/=1 for: ",imat,jmat,mat_cmplx(imat,jmat)
         ans=.FALSE.
       end if
     else
       if(ABS(mat_cmplx(imat,jmat))>tol) then
         write(stdout,*) "M(imat,jmat)/=0 for: ",imat,jmat,mat_cmplx(imat,jmat)
         ans=.FALSE.
       end if
     end if
   end do
 end do

end subroutine check_identity_cmplx


!==========================================
subroutine write_restart_tddft(nstate,time_cur,c_matrix_orth_cmplx)
 use m_definitions
 implicit none
 integer,intent(in)         :: nstate
 complex(dp),intent(in)     :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 real(dp),intent(in)        :: time_cur
!===
 integer                    :: restartfile
 integer                    :: istate,ispin
!=====

 if( .NOT. is_iomaster) return

 call start_clock(timing_restart_tddft_file)

 if (.NOT. in_tddft_loop) then
   write(stdout,'(/,a,f19.10)') ' Writing a RESTART_TDDFT file, time_cur= ', time_cur
 endif
 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',action='write')
 ! current time
 write(restartfile) time_cur
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
 ! Complex wavefunction coefficients C
 do ispin=1,nspin
   do istate=1,nocc
     write(restartfile) c_matrix_orth_cmplx(:,istate,ispin)
   enddo
 enddo

 close(restartfile)

 call stop_clock(timing_restart_tddft_file)

end subroutine write_restart_tddft


!==========================================
subroutine check_restart_tddft(nstate,occupation,restart_is_correct)
 use m_definitions
 use m_density_tools
 implicit none
 logical,intent(out)        :: restart_is_correct
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
!===
 logical                    :: file_exists
 integer                    :: nocc_check
 integer                    :: restartfile
 integer                    :: istate,ispin
 integer                    :: natom_read
 real(dp)                   :: time_cur_read
 real(dp),allocatable       :: zatom_read(:),x_read(:,:)
 integer                    :: nstate_read, nspin_read,nocc_read
!=====

 write(stdout,'(/,a)') ' Checking RESTART_TDDFT file'

 restart_is_correct=.TRUE.

 ! Find highest occupied state
 nocc_check = get_number_occupied_states(occupation)

 inquire(file='RESTART_TDDFT',exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,a)') ' No RESTART file found'
   restart_is_correct=.FALSE.
   return
 endif

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',status='old',action='read')

 ! current time
 read(restartfile) time_cur_read

 !Different number of atoms in restart and input files is not provided for tddft restart
 !natom
 read(restartfile) natom_read
 if( natom_read /= natom ) then
   call issue_warning('RESTART_TDDFT file: natom is not the same.')
   restart_is_correct=.FALSE.
   close(restartfile)
   return 
 end if

 allocate(zatom_read(natom_read),x_read(3,natom_read))
 read(restartfile) zatom_read(1:natom_read)
 read(restartfile) x_read(:,1:natom_read)
 if( natom_read /= natom  &
  .OR. ANY( ABS( zatom_read(1:MIN(natom_read,natom)) - zatom(1:MIN(natom_read,natom)) ) > 1.0e-5_dp ) &
  .OR. ANY( ABS(   x_read(:,1:MIN(natom_read,natom-nprojectile)) - xatom(:,1:MIN(natom_read,natom-nprojectile))   ) > 1.0e-5_dp ) ) then
   call issue_warning('RESTART_TDDFT file: Geometry has changed')
 endif
 deallocate(zatom_read,x_read)

 ! nocc
 read(restartfile) nocc_read
 if(nocc_check /= nocc_read) then
   call issue_warning('RESTART_TDDFT file: nocc is not the same, restart file will not be used')
   restart_is_correct=.FALSE.
   close(restartfile)
   return
 end if

 ! nstate
 read(restartfile) nstate_read
 if(nstate /= nstate_read) then
   call issue_warning('RESTART_TDDFT file: nstate is not the same, restart file will not be used')
   restart_is_correct=.FALSE.
   close(restartfile)
   return
 end if

 ! nspin
 read(restartfile) nspin_read
 if(nspin /= nspin_read) then
   call issue_warning('RESTART_TDDFT file: nspin is not the same, restart file will not be used')
   restart_is_correct=.FALSE.
   close(restartfile)
   return
 end if

 close(restartfile)

end subroutine check_restart_tddft


!==========================================
subroutine get_time_min_restart(time_min)
 use m_definitions
 implicit none
 real(dp),intent(inout)     :: time_min
!===
 logical                    :: file_exists
 integer                    :: restartfile
!=====

 write(stdout,'(/,a)') ' Getting time_min from RESTART_TDDFT file'

 inquire(file='RESTART_TDDFT',exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,a)') ' No RESTART file found'
   call die("No RESTART_TDDFT file found for the second read")
 endif

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',status='old',action='read')

 ! current time
 read(restartfile) time_min
 write(stdout,"(1x,a,f7.3)") "time_min= ", time_min

 close(restartfile)

end subroutine get_time_min_restart


!==========================================
subroutine read_restart_tddft(nstate,time_min,c_matrix_orth_cmplx)
 use m_definitions
 implicit none
 complex(dp),intent(inout)  :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 real(dp),intent(inout)     :: time_min
 integer,intent(in)         :: nstate
!===
 logical                    :: file_exists
 integer                    :: restartfile
 integer                    :: istate,ispin
 integer                    :: natom_read
 real(dp),allocatable       :: zatom_read(:),x_read(:,:)
 integer                    :: nstate_read, nspin_read,nocc_read
!=====

 write(stdout,'(/,a)') ' Reading a RESTART_TDDFT file'

 inquire(file='RESTART_TDDFT',exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,a)') ' No RESTART file found'
   return
 endif

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',status='old',action='read')

 ! current time
 read(restartfile) time_min
 write(stdout,"(1x,a,f7.3)") "time_min= ", time_min

 ! Atomic structure
 read(restartfile) natom_read
 allocate(zatom_read(natom_read),x_read(3,natom_read))
 read(restartfile) zatom_read(1:natom_read)
 read(restartfile) xatom_start(:,1:natom_read)
 deallocate(zatom_read)

 ! nocc
 read(restartfile) nocc_read

 ! Nstate
 read(restartfile) nstate_read

 ! nspin
 read(restartfile) nspin_read

 ! Complex wavefunction coefficients C
 do ispin=1,nspin
   do istate=1,nocc
     read(restartfile) c_matrix_orth_cmplx(:,istate,ispin)
   enddo
 enddo

 close(restartfile)

end subroutine read_restart_tddft


!==========================================
subroutine get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist_cur)
 use m_definitions
 implicit none
 integer, intent(in)          :: n_hist_cur
 real(dp),intent(in)          :: m_nodes(n_hist_cur)
 real(dp),intent(in)          :: x_pred
 real(dp),intent(inout)       :: extrap_coefs(n_hist_cur)
!=====
 integer      :: inode,jnode
!=====

 extrap_coefs(:)=1.0_dp
 do inode=1,n_hist_cur
   do jnode=1,n_hist_cur
     if(inode==jnode) cycle
     extrap_coefs(inode)=extrap_coefs(inode)*(x_pred-m_nodes(jnode))/(m_nodes(inode)-m_nodes(jnode))
   end do
 end do
end subroutine get_extrap_coefs_lagr


!==========================================
subroutine get_extrap_coefs_aspc(extrap_coefs,n_hist_cur)
 use m_definitions
 implicit none
 integer, intent(in)          :: n_hist_cur!
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
subroutine propagate_orth_ham_1(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,s_matrix_sqrt_inv,prop_type)
 implicit none
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: time_step_cur
 complex(dp),intent(inout)   :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 complex(dp), intent(inout)  :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(in)      :: h_small_cmplx(nstate,nstate,nspin)
 real(dp),intent(in)         :: s_matrix_sqrt_inv(basis%nbf,nstate)
 character(len=4),intent(in) :: prop_type
!=====
 integer                    :: ispin
 integer                    :: istate,jstate
 complex(dp),allocatable    :: m_tmp_1(:,:)
 complex(dp),allocatable    :: m_tmp_3(:,:)
!==variables for the MAG2 propagator
 complex(dp),allocatable    :: a_matrix_orth_cmplx(:,:)
 real(dp),allocatable       :: energies_inst(:)
!==variables for the CN propagator
 complex(dp),allocatable    :: l_matrix_cmplx(:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al,
 complex(dp),allocatable    :: b_matrix_cmplx(:,:) ! TDDFT Book, Springer (2006), !p205
!=====
 complex(dp)                :: s_matrix_sqrt_inv_cmplx(basis%nbf,nstate)
!=====

 call start_clock(timing_tddft_propagation)

 s_matrix_sqrt_inv_cmplx = s_matrix_sqrt_inv

 do ispin=1, nspin
   select case(prop_type)
   case('CN')
     allocate(l_matrix_cmplx(nstate,nstate))
     allocate(b_matrix_cmplx(nstate,nstate))
     l_matrix_cmplx(:,:)= im * time_step_cur / 2.0_dp * h_small_cmplx(:,:,ispin)
     b_matrix_cmplx(:,:)=-l_matrix_cmplx(:,:)
     do istate=1,nstate
       b_matrix_cmplx(istate,istate)=b_matrix_cmplx(istate,istate)+1.0_dp
       l_matrix_cmplx(istate,istate)=l_matrix_cmplx(istate,istate)+1.0_dp
     end do
     call invert(l_matrix_cmplx)

     b_matrix_cmplx(:,:)            = MATMUL( l_matrix_cmplx(:,:),b_matrix_cmplx(:,:))
     c_matrix_orth_cmplx(:,:,ispin) = MATMUL( b_matrix_cmplx(:,:),c_matrix_orth_cmplx(:,:,ispin))

!    c_matrix_orth_cmplx(:,:,ispin) = MATMUL( l_matrix_cmplx(:,:),MATMUL( b_matrix_cmplx(:,:),c_matrix_orth_cmplx(:,:,ispin) ) )
     deallocate(l_matrix_cmplx)
     deallocate(b_matrix_cmplx)
   case('MAG2')
     allocate(a_matrix_orth_cmplx(nstate,nstate))
     allocate(energies_inst(nstate))

     !
     ! First part, diagonalize
     call start_clock(timing_propagate_diago)
     a_matrix_orth_cmplx(:,:) = h_small_cmplx(:,:,ispin)
     call diagonalize_scalapack(postscf_diago_flavor,scalapack_block_min,a_matrix_orth_cmplx,energies_inst)
     call stop_clock(timing_propagate_diago)

     !
     ! Second part, multiply matrices
     call start_clock(timing_propagate_matmul)

     allocate(m_tmp_1(nstate,nstate))
     forall (jstate=1:nstate)
       m_tmp_1(:,jstate) = a_matrix_orth_cmplx(:,jstate) * EXP(-im*time_step_cur*energies_inst(jstate) )
     end forall

     allocate(m_tmp_3(nstate,nocc))
     call matmul_abc_scalapack(scalapack_block_min,m_tmp_1,CONJG(TRANSPOSE(a_matrix_orth_cmplx(:,:))),c_matrix_orth_cmplx(:,:,ispin),m_tmp_3  )
     c_matrix_orth_cmplx(:,:,ispin) = m_tmp_3

     deallocate(m_tmp_3)
     deallocate(m_tmp_1)
     deallocate(energies_inst)
     deallocate(a_matrix_orth_cmplx)
     call stop_clock(timing_propagate_matmul)

     case default
       call die('Invalid choice for the propagation algorithm. Change prop_type or error_prop_types value in the input file')
   end select

   call matmul_ab_scalapack(scalapack_block_min,s_matrix_sqrt_inv_cmplx,c_matrix_orth_cmplx(:,:,ispin),c_matrix_cmplx(:,:,ispin))

 end do

 call stop_clock(timing_tddft_propagation)


end subroutine propagate_orth_ham_1


!==========================================
subroutine propagate_orth_ham_2(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_hist2_cmplx,s_matrix_sqrt_inv,prop_type)
 implicit none
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: time_step_cur
 complex(dp),intent(inout)   :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 complex(dp), intent(inout)  :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(in)      :: h_small_hist2_cmplx(nstate,nstate,nspin,2)
 real(dp),intent(in)         :: s_matrix_sqrt_inv(basis%nbf,nstate)
 character(len=4),intent(in) :: prop_type
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
   select case (prop_type)
   case('ETRS')
     do iham=1,2
       call diagonalize(postscf_diago_flavor,h_small_hist2_cmplx(:,:,ispin,iham),energies_inst,a_matrix_orth_cmplx(:,:,iham))
       propagator_eigen(:,:,iham) = ( 0.0_dp , 0.0_dp )
       do ibf=1,nstate
         propagator_eigen(ibf,ibf,iham) = EXP(-im*time_step_cur/2.d0*energies_inst(ibf))
       end do
     end do
     c_matrix_orth_cmplx(:,:,ispin) = &
         MATMUL(MATMUL(MATMUL(MATMUL( MATMUL( MATMUL( a_matrix_orth_cmplx(:,:,2), propagator_eigen(:,:,2)  ) , CONJG(TRANSPOSE(a_matrix_orth_cmplx(:,:,2)))  ),  &
                            a_matrix_orth_cmplx(:,:,1)), propagator_eigen(:,:,1)), CONJG(TRANSPOSE(a_matrix_orth_cmplx(:,:,1))) ), c_matrix_orth_cmplx(:,:,ispin) )
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
                                         hamiltonian_fock_cmplx)
                                         
 implicit none
!=====
 type(basis_set),intent(in)      :: basis
 integer,intent(in)              :: nstate
 integer,intent(in)              :: itau
 real(dp),intent(in)             :: time_cur
 real(dp),intent(in)             :: time_step_cur
 real(dp),intent(in)             :: occupation(nstate,nspin)
 real(dp),intent(in)             :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),allocatable,intent(in) :: dipole_basis(:,:,:)
 real(dp),intent(in)             :: s_matrix_sqrt_inv(basis%nbf,nstate)
 complex(dp),intent(in)          :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(out)         :: hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin)
 complex(dp),intent(out)         :: h_small_cmplx(nstate,nstate,nspin)
!=====
 logical              :: calc_excit_
 integer              :: ispin, idir
 real(dp)             :: excit_field(3)
 complex(dp)          :: p_matrix_cmplx(basis%nbf,basis%nbf,nspin)
 complex(dp)          :: s_matrix_sqrt_inv_cmplx(basis%nbf,nstate)
 integer              :: projectile_list(1)
 real(dp),allocatable :: hamiltonian_projectile(:,:)
!=====

 call start_clock(timing_tddft_hamiltonian)

 s_matrix_sqrt_inv_cmplx = s_matrix_sqrt_inv

 call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)

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



 !
 ! Excitation part of the Hamiltonian
 !
 en%excit = 0.0_dp

 select case(excit_type%form)
 !
 ! Light excitation
 case(EXCIT_LIGHT)
   excit_field=0.0_dp
   calc_excit_ = .FALSE.
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'GAU' )
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'HSW'  .AND. abs(time_cur - excit_type%time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'STEP' .AND. abs(time_cur - excit_type%time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'DEL'  .AND. abs(time_cur - excit_type%time0)<=time_step_cur )
   if(itau==0) calc_excit_=.FALSE.
   if ( calc_excit_ ) then
     call calculate_excit_field(time_cur,excit_field)
     excit_field_norm=NORM2(excit_field(:))
     do idir=1,3
       do ispin=1, nspin
         hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) - dipole_basis(:,:,idir) * excit_field(idir)
         en%excit = en%excit + REAL(SUM(dipole_basis(:,:,idir) * excit_field(idir) * p_matrix_cmplx(:,:,ispin)),dp)
       enddo
     end do
   end if

 !
 ! Projectile excitation
 case(EXCIT_PROJECTILE)

   !
   ! Move the projectile 
   call change_position_one_atom(natom,xatom_start(:,natom) + vel(:,natom) * ( time_cur - time_read ))

   call nucleus_nucleus_energy(en%nuc_nuc)

   !
   ! Nucleus-electron interaction due to the projectile only
   projectile_list(1) = natom
   allocate(hamiltonian_projectile(basis%nbf,basis%nbf))
   call setup_nucleus(basis,hamiltonian_projectile,projectile_list)
   
   do ispin=1,nspin
     hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) + hamiltonian_projectile(:,:)
   enddo
   en%excit = REAL( SUM( hamiltonian_projectile(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
   deallocate(hamiltonian_projectile)

 end select

 do ispin=1,nspin
   hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
 enddo



 !   h_small_cmplx(:,:,ispin) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
 !                   MATMUL( hamiltonian_fock_cmplx(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )
 call start_clock(timing_tddft_ham_orthobasis)
 do ispin=1,nspin
   call matmul_transaba_scalapack(scalapack_block_min,s_matrix_sqrt_inv_cmplx,hamiltonian_fock_cmplx(:,:,ispin),h_small_cmplx(:,:,ispin))
 end do ! spin loop
 call stop_clock(timing_tddft_ham_orthobasis)

 ! kinetic and nuclei-electrons energy contributions
 en%kin = REAL( SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
 en%nuc = REAL( SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)

 call stop_clock(timing_tddft_hamiltonian)

end subroutine setup_hamiltonian_fock_cmplx


!==========================================
subroutine calculate_excit_field(time_cur,excit_field)
 implicit none
 real(dp),intent(in)      :: time_cur       ! time in au
 real(dp),intent(inout)   :: excit_field(3) ! electric field in 3 dimensions
!=====
 real(dp)                 :: excit_dir_norm(3)
!=====

 excit_dir_norm(:)=excit_type%dir(:)/NORM2(excit_type%dir(:))

 select case(excit_type%name)
 case('GAU') !Gaussian electic field
   excit_field(:) = excit_type%kappa * EXP( -( time_cur-excit_type%time0 )**2 / 2.0_dp / excit_omega**2 ) * &
                  & excit_dir_norm(:)
 case('HSW') !Hann sine window
   excit_field(:) = excit_type%kappa * SIN( pi / excit_omega * ( time_cur - excit_type%time0  ) )**2 * excit_dir_norm(:)
 case('DEL') ! Delta excitation
   excit_field(:) = excit_type%kappa * excit_dir_norm(:)
 case('STEP') ! Step excitation
   excit_field(:) = excit_type%kappa * excit_dir_norm(:)
 case default
    call die('Invalid choice for the excitation type. Change excit_type value in the input file')
 end select

end subroutine calculate_excit_field

!=========================================================================
end module m_tddft_propagator
!=========================================================================
