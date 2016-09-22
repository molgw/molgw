!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! subroutine to generate clever representations of the virtual orbital space
!
!=========================================================================
module m_virtual_orbital_space
 use m_definitions
 use m_timing
 use m_mpi
 use m_warning
 use m_memory

 real(dp),allocatable,private :: energy_ref(:,:)
 real(dp),allocatable,private :: c_matrix_ref(:,:,:)

 real(dp),allocatable,private :: c_local(:,:,:)              ! C coefficients:  basis%nbf x nstate
 integer,private              :: desc_bb_sb(ndel)            ! and the corresponding descriptor


contains

!=========================================================================
subroutine setup_virtual_smallbasis(basis,nstate,occupation,nsemax,energy,c_matrix,nstate_small)
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_hamiltonian
 use m_hamiltonian_sca
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 integer,intent(in)                    :: nsemax
 real(dp),intent(inout)                :: energy(nstate,nspin)
 real(dp),intent(inout)                :: c_matrix(basis%nbf,nstate,nspin)
 integer,intent(out)                   :: nstate_small
!=====
 integer                               :: ispin
 integer                               :: ibf
 integer                               :: istate,jstate
 type(basis_set)                       :: basis_small
 real(dp),allocatable                  :: s_bigsmall(:,:)
 real(dp),allocatable                  :: s_small(:,:)
 real(dp),allocatable                  :: s_small_sqrt_inv(:,:)
 real(dp),allocatable                  :: h_small(:,:,:)
 real(dp),allocatable                  :: energy_small(:,:)
 real(dp),allocatable                  :: c_small(:,:,:)
 real(dp),allocatable                  :: c_big(:,:,:)
 real(dp),allocatable                  :: s_matrix(:,:)
 real(dp),allocatable                  :: h_big(:,:,:)
 real(dp),allocatable                  :: s_matrix_inv(:,:)
 real(dp),allocatable                  :: matrix_tmp(:,:)
 real(dp),allocatable                  :: s_bar(:,:),h_bar(:,:,:),s_bar_sqrt_inv(:,:)
 real(dp),allocatable                  :: energy_bar(:,:),c_bar(:,:,:)
 integer                               :: nstate_bar
 integer                               :: nocc,nfrozen
!=====

 call start_clock(timing_fno)

 write(stdout,'(/,1x,a)') 'Prepare optimized empty states using a smaller basis set'


 ! Remember how to go from the small basis set to the big one
 ! 
 ! | \phi^small_a > = \sum_{BC} | \phi^big_B > * S^-1_CB * Sbs_aC
 !
 ! This is key to everything else!


 ! Initialize the small basis set
 write(stdout,*) 'Set up a smaller basis set to optimize the virtual orbital space'
 call init_basis_set(basis_path,small_basis_name,gaussian_type,basis_small)
 call issue_warning('Reduce the virtual orbital subspace by using a smaller basis set: '//TRIM(small_basis_name))

 ! Get the overlap matrix of the wavefunction basis set S: s_matrix
 call clean_allocate('Overlap matrix S',s_matrix,basis%nbf,basis%nbf)
 call clean_allocate('Overlap inverse S^{-1}',s_matrix_inv,basis%nbf,basis%nbf)
 call setup_overlap(.FALSE.,basis,s_matrix)
 call invert(basis%nbf,s_matrix,s_matrix_inv)

 ! Calculate the mixed overlap matrix Sbs: s_bigsmall
 call clean_allocate('Big-Small overlap Sbs',s_bigsmall,basis%nbf,basis_small%nbf)
 call setup_overlap_mixedbasis(.FALSE.,basis,basis_small,s_bigsmall)

 ! Calculate the overlap matrix in the small basis:
 !  tilde S = Sbs**T *  S**-1 * Sbs
 call clean_allocate('Overlap matrix Ssmall',s_small,basis_small%nbf,basis_small%nbf)
 s_small(:,:) = MATMUL( TRANSPOSE(s_bigsmall) , MATMUL( s_matrix_inv , s_bigsmall ) )

 ! Calculate ( tilde S )^{-1/2}
 call setup_sqrt_overlap(min_overlap,basis_small%nbf,s_small,nstate_small,s_small_sqrt_inv)
 call clean_deallocate('Overlap matrix Ssmall',s_small)


 ! Obtain the Hamiltonian in the big basis and in the small basis
 !
 ! H = S * C * E * C**T * S
 ! and 
 ! tilde H = Sbs**T * S**-1 * H * S**-1 * Sbs
 !         = Sbs**T *     C * E * C**T  * Sbs
 !
 call clean_allocate('Hamiltonian H',h_big,basis%nbf,basis%nbf,nspin)
 call clean_allocate('Hamiltonian small basis',h_small,basis_small%nbf,basis_small%nbf,nspin)

 call clean_allocate('Tmp matrix',matrix_tmp,basis%nbf,basis%nbf)
 do ispin=1,nspin
   
   ! M = E * C**T
   do istate=1,nstate
     matrix_tmp(istate,:) = energy(istate,ispin) * c_matrix(:,istate,ispin)
   enddo
   ! M = C * E * C**T
   matrix_tmp(:,:) = MATMUL( c_matrix(:,1:nstate,ispin) , matrix_tmp(1:nstate,:) )
   
   ! H = S * M * S
   h_big(:,:,ispin) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )
   ! Hsmall = Sbs**T * M * Sbs
   h_small(:,:,ispin) = MATMUL( TRANSPOSE(s_bigsmall) , MATMUL( matrix_tmp , s_bigsmall ) )

 enddo
 call clean_deallocate('Tmp matrix',matrix_tmp)

 ! Diagonalize the small Hamiltonian in the small basis
 !
 ! tilde H * tilde C = tilde S * tilde C * tilde E
 !
 allocate(energy_small(nstate_small,nspin))
 call clean_allocate('Coefficients small basis',c_small,basis_small%nbf,nstate_small,nspin)
 call diagonalize_hamiltonian_scalapack(nspin,basis_small%nbf,nstate_small,h_small,s_small_sqrt_inv,energy_small,c_small)
 call dump_out_energy('=== Energies in the initial small basis ===',&
              nstate_small,nspin,occupation(1:nstate_small,:),energy_small)

 call clean_deallocate('Hamiltonian small basis',h_small)
 call clean_deallocate('Overlap sqrt S^{-1/2}',s_small_sqrt_inv)
 deallocate(energy_small)


 !
 ! Transform the wavefunction coefficients from the small basis to the big basis
 ! tilde C -> Cbig
 call clean_allocate('Small wavefunction coeff C',c_big,basis%nbf,nstate_small,nspin)
 ! M = S^-1

 ! Cbig = S**-1 * Sbs * tilde C
 do ispin=1,nspin
   c_big(:,:,ispin) = MATMUL( s_matrix_inv(:,:) , MATMUL( s_bigsmall(:,:) , c_small(:,:,ispin) ) )
 enddo
 call clean_deallocate('Coefficients small basis',c_small)
 call clean_deallocate('Overlap inverse S^{-1}',s_matrix_inv)
 call clean_deallocate('Big-Small overlap Sbs',s_bigsmall)

 ! 
 ! Frozen orbitals for occupied state plus the selfenergy braket
 !
 ! Find the highest occupied state
 do istate=1,nstate
   if( ALL(occupation(istate,:) < completely_empty) ) cycle
   nocc = istate
 enddo

 ! Override the Cbig coefficients with the original C coefficients up to max(nocc,nsemax)
 nfrozen = MAX(nocc,nsemax)

 ! Avoid separating degenerate states
 do while( ANY( ABS(energy(nfrozen+1,:)-energy(nfrozen,:)) < 1.0e-4_dp ) )
   nfrozen = nfrozen + 1
   if( nfrozen == nstate_small ) exit
 end do
 
 write(stdout,'(1x,a,i6)') 'Leave the first states frozen up to: ',nfrozen
 c_big(:,1:nfrozen,:) = c_matrix(:,1:nfrozen,:)

 !
 ! Final diagonalization of in the composite basis
 ! with frozen orbitals complemented with the small basis
 !
 ! Calculate the corresponding overlap matrix Sbar and hamiltonian Hbar
 call clean_allocate('Overlap selected states',s_bar,nstate_small,nstate_small)
 call clean_allocate('Hamiltonian selected states',h_bar,nstate_small,nstate_small,nspin)
 s_bar(:,:) = MATMUL( TRANSPOSE(c_big(:,:,1)) , MATMUL( s_matrix , c_big(:,:,1) ) )

 call clean_deallocate('Overlap matrix S',s_matrix)

 call setup_sqrt_overlap(min_overlap,nstate_small,s_bar,nstate_bar,s_bar_sqrt_inv)
 if( nstate_small /= nstate_bar ) call die('virtual_smallbasis: this usually never happens')
 call clean_deallocate('Overlap selected states',s_bar)

 do ispin=1,nspin
   h_bar(:,:,ispin) = MATMUL( TRANSPOSE(c_big(:,:,ispin)) , MATMUL( h_big(:,:,ispin) , c_big(:,:,ispin) ) )
 enddo
 call clean_deallocate('Hamiltonian H',h_big)

 allocate(energy_bar(nstate_bar,nspin))
 call clean_allocate('Selected states coeffs C',c_bar,nstate_small,nstate_bar,nspin)
 call diagonalize_hamiltonian_scalapack(nspin,nstate_small,nstate_bar,h_bar,s_bar_sqrt_inv,energy_bar,c_bar)


 do ispin=1,nspin
   c_big(:,1:nstate_bar,ispin) = MATMUL( c_big(:,:,ispin) , c_bar(:,:,ispin) )
 enddo

 call dump_out_energy('=== Energies in the final small basis ===',&
              nstate_bar,nspin,occupation(1:nstate_bar,:),energy_bar)

 call clean_deallocate('Overlap sqrt S^{-1/2}',s_bar_sqrt_inv)
 call clean_deallocate('Hamiltonian selected states',h_bar)
 call clean_deallocate('Selected states coeffs C',c_bar)

 
 !
 ! Save the original c_matrix and energies
 if( .NOT. ALLOCATED(energy_ref) ) then
   allocate(energy_ref(nstate_bar,nspin))
   allocate(c_matrix_ref(basis%nbf,nstate_bar,nspin))   ! TODO clean_allocate of this

   energy_ref(:,:)     = energy(1:nstate_bar,:)
   c_matrix_ref(:,:,:) = c_matrix(:,1:nstate_bar,:)
 endif

 !
 ! And then override the c_matrix and the energy with the fresh new ones
 energy(1:nstate_bar,:)     = energy_bar(:,:)
 c_matrix(:,1:nstate_bar,:) = c_big(:,:,:)

 nstate_small = nstate_bar

 deallocate(energy_bar)
 call clean_deallocate('Small wavefunction coeff C',c_big)

 call destroy_basis_set(basis_small)

 write(stdout,'(1x,a)') 'Optimized empty states with a smaller basis set'

 call stop_clock(timing_fno)

end subroutine setup_virtual_smallbasis


!=========================================================================
subroutine setup_virtual_smallbasis_sca(basis,nstate,occupation,nsemax,energy,c_matrix,nstate_small)
 use m_inputparam
 use m_basis_set
 use m_hamiltonian
 use m_hamiltonian_sca
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 integer,intent(in)         :: nsemax
 real(dp),intent(inout)     :: energy(nstate,nspin)
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
 integer,intent(out)        :: nstate_small
!=====
 integer                    :: desc_bb_bb(ndel)            ! basis%nbf x basis%nbf
 integer                    :: desc_bb_bs(ndel)            ! basis%nbf x basis_small%nbf
 integer                    :: desc_bs_bs(ndel)            ! basis_small%nbf x basis_small%nbf
 integer                    :: desc_bs_ss(ndel)            ! basis_small%nbf x nstate_small
 integer                    :: desc_bb_ss(ndel)            ! basis%nbf x nstate_small
 integer                    :: desc_ss_ss(ndel)            ! nstate_small x nstate_small
 integer                    :: desc_tmp(ndel)              ! Dummy
 integer                    :: nprow,npcol,iprow,ipcol
 integer                    :: cntxt
 integer                    :: rank_master
 integer                    :: ma,na
 integer                    :: mb,nb
 integer                    :: mc,nc
 integer                    :: md,nd
 integer                    :: me,ne
 integer                    :: mf,nf
 integer                    :: mg,ng
 integer                    :: ispin
 integer                    :: ibf
 integer                    :: istate,jstate
 type(basis_set)            :: basis_small
 real(dp),allocatable       :: s_bigsmall_global(:,:)   !TODO: remove this in the future
 real(dp),allocatable       :: s_bigsmall(:,:)
 real(dp),allocatable       :: s_small(:,:)
 real(dp),allocatable       :: s_small_sqrt_inv(:,:)
 real(dp),allocatable       :: h_small(:,:,:)
 real(dp),allocatable       :: energy_small(:,:)
 real(dp),allocatable       :: c_small(:,:,:)
 real(dp),allocatable       :: c_big(:,:,:)
 real(dp),allocatable       :: s_matrix(:,:)
 real(dp),allocatable       :: h_big(:,:,:)
 real(dp),allocatable       :: s_matrix_inv(:,:)
 real(dp),allocatable       :: matrix_tmp1(:,:)
 real(dp),allocatable       :: matrix_tmp2(:,:)
 real(dp),allocatable       :: s_bar(:,:),h_bar(:,:,:),s_bar_sqrt_inv(:,:)
 real(dp),allocatable       :: energy_bar(:,:),c_bar(:,:,:)
 integer                    :: nocc,nfrozen
 integer                    :: info
 integer                    :: jglobal,jlocal
!=====

 call start_clock(timing_fno)

 write(stdout,'(/,1x,a)') 'Prepare optimized empty states using a smaller basis set'

 !
 ! All procs save the original energies
 allocate(energy_ref(nstate,nspin))
 energy_ref(:,:)     = energy(:,:)

 ! Remember how to go from the small basis set to the big one
 ! 
 ! | \phi^small_a > = \sum_{AB} | \phi^big_A > * (S**-1)_AB * Sbs_Ba
 !
 ! This is key to everything else!


 ! First set up the SCALAPACK grid
 cntxt = desc_ham(CTXT_A)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK with a grid',nprow,' x ',npcol

 ! Find the master
 if( iprow == 0 .AND. ipcol == 0 ) then
   rank_master = rank_world
 else
   rank_master = -1
 endif
 call xmax_world(rank_master)

 ! Initialize the small basis set
 write(stdout,*) 'Set up a smaller basis set to optimize the virtual orbital space'
 call init_basis_set(basis_path,small_basis_name,gaussian_type,basis_small)
 call issue_warning('Reduce the virtual orbital subspace by using a smaller basis set: '//TRIM(small_basis_name))

 if( cntxt > 0 ) then

   !
   ! Descriptor desc_bb_sb   mf,nf
   mf = NUMROC(basis%nbf,block_row,iprow,first_row,nprow)
   nf = NUMROC(nstate   ,block_col,ipcol,first_col,npcol)
   call DESCINIT(desc_bb_sb,basis%nbf,nstate,block_row,block_col,first_row,first_col,cntxt,MAX(1,mf),info) 
   call clean_allocate('Distributed C',c_local,mf,nf,nspin)
   !
   ! Distribute C -> C_local
   call create_distributed_copy(c_matrix,desc_bb_sb,c_local)

   !
   ! Descriptor desc_bb_bb   ma,na
   ma = NUMROC(basis%nbf,block_row,iprow,first_row,nprow)
   na = NUMROC(basis%nbf,block_col,ipcol,first_col,npcol)
   call DESCINIT(desc_bb_bb,basis%nbf,basis%nbf,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info) 
   call clean_allocate('Overlap matrix S',s_matrix,ma,na)
   call clean_allocate('Overlap inverse S^{-1}',s_matrix_inv,ma,na)
   call setup_overlap_sca(.FALSE.,basis,ma,na,s_matrix)
   call invert_sca(desc_bb_bb,s_matrix,s_matrix_inv)

   ! Calculate the mixed overlap matrix Sbs: s_bigsmall
   !TODO: Distribute this from the beginning
   call clean_allocate('Big-Small overlap Sbs',s_bigsmall_global,basis%nbf,basis_small%nbf)
   call setup_overlap_mixedbasis(.FALSE.,basis,basis_small,s_bigsmall_global)
   !
   ! Descriptor desc_bb_bs  mb,nb
   mb = NUMROC(basis%nbf      ,block_row,iprow,first_row,nprow)
   nb = NUMROC(basis_small%nbf,block_col,ipcol,first_col,npcol)
   call DESCINIT(desc_bb_bs,basis%nbf,basis_small%nbf,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info) 
   call clean_allocate('Big-Small overlap Sbs',s_bigsmall,mb,nb)
   call create_distributed_copy(s_bigsmall_global,desc_bb_sb,s_bigsmall)
   call clean_deallocate('Big-Small overlap Sbs',s_bigsmall_global)
      
   !
   ! Descriptor desc_bs_bs  mc,nc
   mc = NUMROC(basis_small%nbf,block_row,iprow,first_row,nprow)
   nc = NUMROC(basis_small%nbf,block_col,ipcol,first_col,npcol)
   call DESCINIT(desc_bs_bs,basis_small%nbf,basis_small%nbf,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info) 

   ! Calculate the overlap matrix in the small basis:
   ! tilde S = Sbs**T *  S**-1 * Sbs
   call clean_allocate('Overlap matrix Ssmall',s_small,mc,nc)
!   s_small(:,:) = MATMUL( TRANSPOSE(s_bigsmall) , MATMUL( s_matrix_inv , s_bigsmall ) )
   call product_transaba_sca(desc_bb_bs,s_bigsmall,desc_bb_bb,s_matrix_inv,desc_bs_bs,s_small)

   ! Calculate ( tilde S )^{-1/2}

   call setup_sqrt_overlap_sca(min_overlap,desc_bs_bs,s_small,desc_bs_ss,s_small_sqrt_inv)
   nstate_small = desc_bs_ss(N_A)
   md = NUMROC(basis_small%nbf,block_row,iprow,first_row,nprow)
   nd = NUMROC(nstate_small   ,block_col,ipcol,first_col,npcol)
   !
   ! Descriptor desc_bs_ss  md,nd
!   TODO remove this line
!   call DESCINIT(desc_bs_ss,basis_small%nbf,nstate_small,block_row,block_col,first_row,first_col,cntxt,MAX(1,md),info) 

   call clean_deallocate('Overlap matrix Ssmall',s_small)


   ! Obtain the Hamiltonian in the big basis and in the small basis
   !
   ! H = S * C * E * C**T * S
   ! and 
   ! tilde H = Sbs**T * S**-1 * H * S**-1 * Sbs
   !         = Sbs**T *     C * E * C**T  * Sbs
   !
   call clean_allocate('Hamiltonian H',h_big,ma,na,nspin)
   call clean_allocate('Hamiltonian small basis',h_small,mc,nc,nspin)


   do ispin=1,nspin
     ! M1 = C * E
     call clean_allocate('Tmp matrix',matrix_tmp1,mf,nf)
     matrix_tmp1(:,:) = c_local(:,:,ispin)
     call matmul_diag_sca('R',energy(:,ispin),desc_bb_sb,matrix_tmp1)
     ! M2 = M1 * C**T = C * E * C**T
     call clean_allocate('Tmp matrix',matrix_tmp2,ma,na)
     call PDGEMM('N','T',basis%nbf,basis%nbf,nstate,  &
                  1.0_dp,matrix_tmp1,1,1,desc_bb_sb,  &
                  c_local(1,1,ispin),1,1,desc_bb_sb,  &
                  0.0_dp,matrix_tmp2,1,1,desc_bb_bb)
     ! H = S**T * M2 * S
     call product_transaba_sca(desc_bb_bb,s_matrix,desc_bb_bb,matrix_tmp2,desc_bb_bb,h_big(:,:,ispin))
     
     ! tilde H = Sbs**T * M2 * Sbs
     call product_transaba_sca(desc_bb_bs,s_bigsmall,desc_bb_bb,matrix_tmp2,desc_bs_bs,h_small(:,:,ispin))

     call clean_deallocate('Tmp matrix',matrix_tmp1)
     call clean_deallocate('Tmp matrix',matrix_tmp2)
   enddo

   ! Diagonalize the small Hamiltonian in the small basis
   !
   ! tilde H * tilde C = tilde S * tilde C * tilde E
   !
   allocate(energy_small(nstate_small,nspin))
   call clean_allocate('Coefficients small basis',c_small,md,nd,nspin)
   call diagonalize_hamiltonian_sca(1,nspin,desc_bs_bs,h_small,desc_bs_ss,s_small_sqrt_inv,energy_small,c_small)

   call dump_out_energy('=== Energies in the initial small basis ===',&
                        nstate_small,nspin,occupation(1:nstate_small,:),energy_small)

   call clean_deallocate('Hamiltonian small basis',h_small)
   call clean_deallocate('Overlap sqrt S^{-1/2}',s_small_sqrt_inv)
   deallocate(energy_small)

   !
   ! Transform the wavefunction coefficients from the small basis to the big basis
   ! tilde C -> Cbig
   me = NUMROC(basis%nbf   ,block_row,iprow,first_row,nprow)
   ne = NUMROC(nstate_small,block_col,ipcol,first_col,npcol)
   call DESCINIT(desc_bb_ss,basis%nbf,nstate_small,block_row,block_col,first_row,first_col,cntxt,MAX(1,me),info) 
   call clean_allocate('Small wavefunction coeff C',c_big,me,ne,nspin)


   ! Cbig = S**-1 * Sbs * tilde C
   do ispin=1,nspin
     call product_abc_sca(desc_bb_bb,s_matrix_inv,desc_bb_bs,s_bigsmall,desc_bs_ss,c_small(:,:,ispin),desc_bb_ss,c_big(:,:,ispin))
   enddo


   call clean_deallocate('Coefficients small basis',c_small)
   call clean_deallocate('Overlap inverse S^{-1}',s_matrix_inv)
   call clean_deallocate('Big-Small overlap Sbs',s_bigsmall)


   ! 
   ! Frozen orbitals for occupied state plus the selfenergy braket
   !
   ! Find the highest occupied state
   do istate=1,nstate
     if( ALL(occupation(istate,:) < completely_empty) ) cycle
     nocc = istate
   enddo
  
   ! Override the Cbig coefficients with the original C coefficients up to max(nocc,nsemax)
   nfrozen = MAX(nocc,nsemax)
  
   ! Avoid separating degenerate states
   do while( ANY( ABS(energy(nfrozen+1,:)-energy(nfrozen,:)) < 1.0e-4_dp ) )
     nfrozen = nfrozen + 1
     if( nfrozen == nstate_small ) exit
   end do
   
   write(stdout,'(1x,a,i6)') 'Leave the first states frozen up to: ',nfrozen

   do jlocal=1,ne
     jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
     if( jglobal > nfrozen ) cycle
     c_big(:,jlocal,:) = c_local(:,jlocal,:)
   enddo

   


   !
   ! Final diagonalization of in the composite basis
   ! with frozen orbitals complemented with the small basis
   !
   ! Calculate the corresponding overlap matrix Sbar and hamiltonian Hbar
   mg = NUMROC(nstate_small,block_row,iprow,first_row,nprow)
   ng = NUMROC(nstate_small,block_col,ipcol,first_col,npcol)
   call DESCINIT(desc_ss_ss,nstate_small,nstate_small,block_row,block_col,first_row,first_col,cntxt,MAX(1,mg),info) 
   call clean_allocate('Overlap selected states',s_bar,mg,ng)
   call clean_allocate('Hamiltonian selected states',h_bar,mg,ng,nspin)

   ! Sbar = C'**T * S * C'
   ! s_bar(:,:) = MATMUL( TRANSPOSE(c_big(:,:,1)) , MATMUL( s_matrix , c_big(:,:,1) ) )
   call product_transaba_sca(desc_bb_ss,c_big(:,:,1),desc_bb_bb,s_matrix,desc_ss_ss,s_bar)

   call clean_deallocate('Overlap matrix S',s_matrix)

   call setup_sqrt_overlap_sca(min_overlap,desc_ss_ss,s_bar,desc_tmp,s_bar_sqrt_inv)
   if( nstate_small /= desc_tmp(N_A) ) call die('virtual_smallbasis: this usually never happens')
   call clean_deallocate('Overlap selected states',s_bar)

   do ispin=1,nspin
     ! Sbar = C'**T * H * C'
     ! h_bar(:,:,ispin) = MATMUL( TRANSPOSE(c_big(:,:,ispin)) , MATMUL( h_big(:,:,ispin) , c_big(:,:,ispin) ) )
     call product_transaba_sca(desc_bb_ss,c_big(:,:,1),desc_bb_bb,h_big(:,:,ispin),desc_ss_ss,h_bar(:,:,ispin))
   enddo

   call clean_deallocate('Hamiltonian H',h_big)

   allocate(energy_bar(nstate_small,nspin))
   call clean_allocate('Selected states coeffs C',c_bar,mg,ng,nspin)
   call diagonalize_hamiltonian_sca(1,nspin,desc_ss_ss,h_bar,desc_ss_ss,s_bar_sqrt_inv,energy_bar,c_bar)


   call clean_allocate('Tmp matrix',matrix_tmp1,me,ne)
   do ispin=1,nspin
     ! C' <- C' * Cbar
     ! c_big(:,:,ispin) = MATMUL( c_big(:,:,ispin) , c_bar(:,:,ispin) )
     matrix_tmp1(:,:) = c_big(:,:,ispin)
     call PDGEMM('N','N',basis%nbf,nstate_small,nstate_small,  &
                  1.0_dp,     matrix_tmp1,1,1,desc_bb_ss,  &
                         c_bar(1,1,ispin),1,1,desc_ss_ss,  &
                  0.0_dp,c_big(1,1,ispin),1,1,desc_bb_ss)
   enddo
   call clean_deallocate('Tmp matrix',matrix_tmp1)


   call dump_out_energy('=== Energies in the final small basis ===',&
                        nstate_small,nspin,occupation(1:nstate_small,:),energy_bar)

   call clean_deallocate('Overlap sqrt S^{-1/2}',s_bar_sqrt_inv)
   call clean_deallocate('Hamiltonian selected states',h_bar)
   call clean_deallocate('Selected states coeffs C',c_bar)

   ! Override the energies now
   energy(1:nstate_small,:) = energy_bar(:,:)
 endif
 !
 ! And then override the c_matrix and the energy with the fresh new ones
 call xbcast_world(rank_master,energy)

 call gather_distributed_copy(desc_bb_ss,c_big,c_matrix(:,1:nstate_small,:))
 if( cntxt > 0 ) then
   call clean_deallocate('Small wavefunction coeff C',c_big)
 endif

 deallocate(energy_bar)

 call destroy_basis_set(basis_small)

 write(stdout,'(1x,a)') 'Optimized empty states with a smaller basis set'

 call stop_clock(timing_fno)

end subroutine setup_virtual_smallbasis_sca


!=========================================================================
subroutine virtual_fno(basis,nstate,occupation,energy,c_matrix)
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(inout)                :: energy(nstate,nspin)
 real(dp),intent(inout)                :: c_matrix(basis%nbf,nstate,nspin)
!=====
 real(dp),parameter                    :: alpha_mp2=1.0_dp
 integer                               :: istate,jstate,astate,bstate,cstate
 integer                               :: ispin
 integer                               :: nocc,ncore,nvirtual,nvirtual_kept
 real(dp),allocatable                  :: p_matrix_mp2(:,:),ham_virtual(:,:),ham_virtual_kept(:,:)
 real(dp),allocatable                  :: occupation_mp2(:),energy_virtual_kept(:)
 real(dp)                              :: eri_ci_aj,eri_ci_bj
 real(dp)                              :: den_ca_ij,den_cb_ij
 real(dp)                              :: en_mp2
 integer                               :: nvirtualmin
 real(dp),allocatable                  :: eri_ci_i(:)
!=====

 call start_clock(timing_fno)

 write(stdout,'(/,1x,a)') 'Prepare optimized empty states with Frozen Natural Orbitals'

 call assert_experimental()

 nvirtualmin = MIN(nvirtualw,nvirtualg)

 if( nvirtualmin > nstate ) then 
   call issue_warning('virtual_fno is on, however nvirtualw and nvirtualg are not set. Skipping the Frozen Natural Orbitals generation.')
   return
 endif

 if( .NOT. has_auxil_basis ) then
   call issue_warning('virtual_fno not implemented when no auxiliary basis is provided')
   return
 endif

 if( nspin > 1 ) then
   call issue_warning('virtual_fno not implemented yet for spin polarized calculations')
   return
 endif

 ncore = ncoreg
 if(is_frozencore) then
   if( ncore == 0) ncore = atoms_core_states()
 endif


 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,1,nstate,1,nstate)


 do ispin=1,nspin
   !
   ! First, set up the list of occupied states
   nocc = ncore
   do istate=ncore+1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     nocc = istate
   enddo



   nvirtual = nstate - nocc

   allocate(p_matrix_mp2(nvirtual,nvirtual))
   p_matrix_mp2(:,:) = 0.0_dp

   allocate(eri_ci_i(nocc+1:nstate))

#if 1
   ! Approximation by Aquilante et al.
   do istate=ncore+1,nocc
     do cstate=nocc+1,nstate

       do astate=nocc+1,nstate
         eri_ci_i(astate) = eri_eigen_ri_paral(cstate,istate,ispin,astate,istate,ispin) &
                              / ( energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(cstate,ispin) )
       enddo
       call xsum_auxil(eri_ci_i)

       do bstate=nocc+1,nstate
         do astate=nocc+1,nstate

           p_matrix_mp2(astate-nocc,bstate-nocc) = &
                p_matrix_mp2(astate-nocc,bstate-nocc)  &
                   + 0.50_dp * eri_ci_i(astate) * eri_ci_i(bstate)

         enddo
       enddo

     enddo
   enddo
   deallocate(eri_ci_i)
#else

   do bstate=nocc+1,nstate
     do astate=nocc+1,nstate

#if 0
       ! Full calculation of the MP2 density matrix on virtual orbitals (See Taube and Bartlett)
       do cstate=nocc+1,nstate
         do istate=ncore+1,nocc
           do jstate=ncore+1,nocc

             den_cb_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(bstate,ispin) - energy(cstate,ispin)
             den_ca_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(astate,ispin) - energy(cstate,ispin)

             eri_ci_aj = eri_eigen_ri(cstate,istate,ispin,astate,jstate,ispin) &
                          - eri_eigen_ri(cstate,jstate,ispin,astate,istate,ispin)  * alpha_mp2
             eri_ci_bj = eri_eigen_ri(cstate,istate,ispin,bstate,jstate,ispin) &
                         - eri_eigen_ri(cstate,jstate,ispin,bstate,istate,ispin)   * alpha_mp2

             p_matrix_mp2(astate-nocc,bstate-nocc) = &
                  p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                     + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

           enddo
         enddo
       enddo
#elif 1
       ! Approximation by Aquilante et al.
       do cstate=nocc+1,nstate
         do istate=ncore+1,nocc
           den_cb_ij = energy(istate,ispin) + energy(istate,ispin) - energy(bstate,ispin) - energy(cstate,ispin)
           den_ca_ij = energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(cstate,ispin)

           eri_ci_aj = eri_eigen_ri(cstate,istate,ispin,astate,istate,ispin) 
           eri_ci_bj = eri_eigen_ri(cstate,istate,ispin,bstate,istate,ispin) 

           p_matrix_mp2(astate-nocc,bstate-nocc) = &
                p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                   + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

         enddo
       enddo
#elif 0
       do istate=ncore+1,nocc
         ! Approximation by Aquilante et al.
         den_cb_ij = energy(istate,ispin) + energy(istate,ispin) - energy(bstate,ispin) - energy(bstate,ispin)
         den_ca_ij = energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(astate,ispin)

         eri_ci_aj = eri_eigen_ri(astate,istate,ispin,astate,istate,ispin) 
         eri_ci_bj = eri_eigen_ri(bstate,istate,ispin,bstate,istate,ispin) 

         p_matrix_mp2(astate-nocc,bstate-nocc) = &
              p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                 + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )
       enddo
#else
         do istate=ncore+1,nocc
           do jstate=ncore+1,nocc

             den_cb_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(bstate,ispin) - energy(bstate,ispin)
             den_ca_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(astate,ispin) - energy(astate,ispin)

             eri_ci_aj = eri_eigen_ri(astate,istate,ispin,astate,jstate,ispin) 
             eri_ci_bj = eri_eigen_ri(bstate,istate,ispin,bstate,jstate,ispin) 

             p_matrix_mp2(astate-nocc,bstate-nocc) = &
                  p_matrix_mp2(astate-nocc,bstate-nocc)  &
                     + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

           enddo
         enddo

#endif

     enddo
   enddo


#endif

   allocate(occupation_mp2(nvirtual))
   call diagonalize(nvirtual,p_matrix_mp2,occupation_mp2)

!   write(stdout,*) 
!   do astate=1,nvirtual
!     write(stdout,'(i4,2x,es16.4)') astate,occupation_mp2(astate)
!   enddo
!   write(stdout,*) 

   allocate(ham_virtual(nvirtual,nvirtual))

   ham_virtual(:,:) = 0.0_dp
   do astate=1,nvirtual
     ham_virtual(astate,astate) = energy(astate+nocc,ispin)
   enddo

   nvirtual_kept = nvirtualmin - 1 - nocc

   write(stdout,'(/,1x,a,i5)')    'Max state index included ',nvirtual_kept + nocc
   write(stdout,'(1x,a,i5,a,i5)') 'Retain ',nvirtual_kept,' virtual orbitals out of ',nvirtual
   write(stdout,'(1x,a,es14.6)')  'Occupation number of the last virtual state',occupation_mp2(nvirtual - (nocc + nvirtual_kept))
   write(stdout,'(1x,a,es14.6,4x,es14.6)') '  to be compared to the first and last virtual states',occupation_mp2(nvirtual),occupation_mp2(1)
   write(stdout,*)

   allocate(ham_virtual_kept(nvirtual_kept,nvirtual_kept))
   allocate(energy_virtual_kept(nvirtual_kept))

   ham_virtual_kept(:,:)  = MATMUL( TRANSPOSE( p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) ) ,   &
                                        MATMUL( ham_virtual , p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) ) )

   deallocate(ham_virtual)

   call diagonalize(nvirtual_kept,ham_virtual_kept,energy_virtual_kept)

!   write(stdout,'(/,1x,a)') ' virtual state    FNO energy (eV)   reference energy (eV)'
!   do astate=1,nvirtual_kept
!     write(stdout,'(i4,4x,f16.5,2x,f16.5)') astate,energy_virtual_kept(astate)*Ha_eV,energy(astate+nocc,ispin)*Ha_eV
!   enddo
!   write(stdout,*)
   

   !
   ! Save the original c_matrix and energies
   if( .NOT. ALLOCATED(energy_ref) ) then
     allocate(energy_ref(nocc+1:nocc+nvirtual_kept,nspin))
     allocate(c_matrix_ref(basis%nbf,nocc+1:nocc+nvirtual_kept,nspin))   ! TODO clean_allocate of this

     energy_ref(nocc+1:nocc+nvirtual_kept,:)     = energy(nocc+1:nocc+nvirtual_kept,:)
     c_matrix_ref(:,nocc+1:nocc+nvirtual_kept,:) = c_matrix(:,nocc+1:nocc+nvirtual_kept,:)
   endif

   !
   ! And then override the c_matrix and the energy with the fresh new ones
   energy(nocc+1:nocc+nvirtual_kept,ispin)     = energy_virtual_kept(:)
   c_matrix(:,nocc+1:nocc+nvirtual_kept,ispin) = MATMUL( c_matrix(:,nocc+1:,ispin) ,  &
                                                    MATMUL( p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) , &
                                                       ham_virtual_kept(:,:) ) )


   deallocate(energy_virtual_kept)
   deallocate(ham_virtual_kept)
   deallocate(p_matrix_mp2,occupation_mp2)
 enddo



 call destroy_eri_3center_eigen()

 write(stdout,'(1x,a)') 'Optimized empty states with Frozen Natural Orbitals'

 call stop_clock(timing_fno)

end subroutine virtual_fno


!=========================================================================
subroutine destroy_fno(basis,nstate,energy,c_matrix)
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_eri_ao_mo
 use m_scalapack
 implicit none

 type(basis_set),intent(in)           :: basis
 integer,intent(in)                   :: nstate
 real(dp),intent(inout)               :: energy(nstate,nspin)
 real(dp),intent(inout)               :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                              :: lowerb,upperb
!=====

 write(stdout,'(/,1x,a)') 'Deallocate the Frozen Natural Orbitals'

 lowerb = LBOUND(energy_ref,DIM=1)
 upperb = UBOUND(energy_ref,DIM=1)


 energy(lowerb:upperb,:) = energy_ref(lowerb:upperb,:)
 deallocate(energy_ref)

 if( ALLOCATED(c_matrix_ref) ) then
   c_matrix(:,lowerb:upperb,:) = c_matrix_ref(:,lowerb:upperb,:)
   deallocate(c_matrix_ref)   ! TODO clean_allocate of this
 else

   call gather_distributed_copy(desc_bb_sb,c_local,c_matrix)
   if( desc_bb_sb(CTXT_A) > 0 ) then
     call clean_deallocate('Distributed C',c_local)
   endif

 endif


end subroutine destroy_fno

!=========================================================================
end module m_virtual_orbital_space
