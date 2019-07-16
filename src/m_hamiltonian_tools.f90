!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to manipulate the Kohn-Sham Hamiltonian and wavefunctions
!
!=========================================================================
module m_hamiltonian_tools
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_basis_set
 use m_linear_algebra
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min,scf_diago_flavor


 interface setup_density_matrix
   module procedure setup_density_matrix_real
   module procedure setup_density_matrix_cmplx
 end interface setup_density_matrix

contains


!=========================================================================
pure function get_number_occupied_states(occupation) result(nocc)
 implicit none

 real(dp),intent(in) :: occupation(:,:)
 integer             :: nocc
!=====
 integer :: nstate,istate,ispin,nspin_local
!=====

 nstate      = SIZE(occupation(:,:),DIM=1)
 nspin_local = SIZE(occupation(:,:),DIM=2)

 ! Find highest occupied state
 ! Take care of negative occupations, this can happen if C comes from P^{1/2}
 nocc = 0
 do ispin=1,nspin_local
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty )  cycle
     nocc = MAX(nocc,istate)
   enddo
 enddo


end function get_number_occupied_states


!=========================================================================
subroutine setup_density_matrix_real(c_matrix,occupation,p_matrix)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(out) :: p_matrix(:,:,:)
!=====
 integer :: nbf,nstate,nocc
 integer :: ispin,ibf,jbf
 integer :: istate
 real(dp),allocatable :: c_matrix_sqrtocc(:,:)
!=====

 call start_clock(timing_density_matrix)

 write(stdout,'(1x,a)') 'Build density matrix'

 nbf    = SIZE(c_matrix(:,:,:),DIM=1)
 nstate = SIZE(c_matrix(:,:,:),DIM=2)

 if( ANY( occupation(:,:) < -1.0e-5_dp ) ) then
   write(stdout,*) 'Min occupation:',MINVAL(occupation)
   call die('setup_density_matrix: negative occupation number should not happen here.')
 endif
 ! Find the number of occupatied states
 nocc = get_number_occupied_states(occupation)

 allocate(c_matrix_sqrtocc(nbf,nocc))

 p_matrix(:,:,:) = 0.0_dp
 do ispin=1,nspin

   do istate=1,nocc
     c_matrix_sqrtocc(:,istate) = c_matrix(:,istate,ispin) * SQRT(occupation(istate,ispin))
   enddo

   call DSYRK('L','N',nbf,nocc,1.0d0,c_matrix_sqrtocc,nbf,0.0d0,p_matrix(1,1,ispin),nbf)

   ! Symmetrize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix(jbf,ibf,ispin) = p_matrix(ibf,jbf,ispin)
     enddo
   enddo
 enddo

 deallocate(c_matrix_sqrtocc)

 call stop_clock(timing_density_matrix)


end subroutine setup_density_matrix_real


!=========================================================================
subroutine setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
 implicit none

 complex(dp),intent(in)  :: c_matrix_cmplx(:,:,:)
 real(dp),intent(in)     :: occupation(:,:)
 complex(dp),intent(out) :: p_matrix_cmplx(:,:,:)
!=====
 integer :: nbf,nstate,nocc
 integer :: ispin,ibf,jbf
 integer :: istate
 complex(dp),allocatable :: c_matrix_sqrtocc(:,:)
!=====

 call start_clock(timing_density_matrix_cmplx)

 nbf    = SIZE(c_matrix_cmplx(:,:,:),DIM=1)
 nocc   = SIZE(c_matrix_cmplx(:,:,:),DIM=2)
 nstate = SIZE(occupation(:,:),DIM=1)

 if( ANY( occupation(:,:) < 0.0_dp ) ) call die('setup_density_matrix_cmplx: negative occupation number should not happen here.')

 allocate(c_matrix_sqrtocc(nbf,nocc))

 p_matrix_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )
 do ispin=1,nspin

   do istate=1,nocc
     c_matrix_sqrtocc(:,istate) = c_matrix_cmplx(:,istate,ispin) * SQRT(occupation(istate,ispin))
   enddo
   call ZHERK('L','N',nbf,nocc,1.0d0,c_matrix_sqrtocc,nbf,0.0d0,p_matrix_cmplx(1,1,ispin),nbf)


   ! Hermitianize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix_cmplx(jbf,ibf,ispin) = CONJG( p_matrix_cmplx(ibf,jbf,ispin) )
     enddo
   enddo

 enddo

 deallocate(c_matrix_sqrtocc)

 call stop_clock(timing_density_matrix_cmplx)


end subroutine setup_density_matrix_cmplx


!=========================================================================
subroutine setup_energy_density_matrix(c_matrix,occupation,energy,q_matrix)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: energy(:,:)
 real(dp),intent(out) :: q_matrix(:,:)
!=====
 integer :: nbf,nstate
 integer :: ispin,ibf,jbf
 integer :: istate
!=====

 call start_clock(timing_density_matrix)

 write(stdout,'(1x,a)') 'Build energy-density matrix'

 nbf    = SIZE(c_matrix(:,:,:),DIM=1)
 nstate = SIZE(c_matrix(:,:,:),DIM=2)

 q_matrix(:,:) = 0.0_dp
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     call DSYR('L',nbf,occupation(istate,ispin)*energy(istate,ispin),c_matrix(:,istate,ispin),1,q_matrix(:,:),nbf)
   enddo
 enddo


 ! Symmetrize
 do jbf=1,nbf
   do ibf=jbf+1,nbf
     q_matrix(jbf,ibf) = q_matrix(ibf,jbf)
   enddo
 enddo

 call stop_clock(timing_density_matrix)


end subroutine setup_energy_density_matrix


!=========================================================================
subroutine test_density_matrix(p_matrix,s_matrix)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(in)  :: s_matrix(:,:)
!=====
 integer              :: nbf
 integer              :: ispin
 real(dp),allocatable :: matrix(:,:)
 character(len=100)   :: title
!=====

 nbf = SIZE(p_matrix(:,:,:),DIM=1)
 allocate(matrix(nbf,nbf))

 write(stdout,*) 'Check equality PSP = P'
 write(stdout,*) ' valid only for integer occupation numbers'

 do ispin=1,nspin

   !
   ! Calculate PSP
   matrix(:,:) = MATMUL( p_matrix(:,:,ispin), MATMUL( s_matrix(:,:) , p_matrix(:,:,ispin) ) )


   title='=== PSP ==='
   !call dump_out_matrix(1,title,nbf,1,matrix)
   title='===  P  ==='
   !call dump_out_matrix(1,title,nbf,1,p_matrix(:,:,ispin))

 enddo

 deallocate(matrix)

end subroutine test_density_matrix


!=========================================================================
subroutine set_occupation(temperature,electrons,magnetization,energy,occupation)
 implicit none
 real(dp),intent(in)  :: electrons,magnetization,temperature
 real(dp),intent(in)  :: energy(:,:)
 real(dp),intent(out) :: occupation(:,:)
!=====
 integer              :: nstate
 real(dp)             :: remaining_electrons(nspin)
 real(dp)             :: electrons_mu,mu,delta_mu,grad_electrons
 real(dp)             :: mu_change
 integer              :: istate,nlines,ilines
 logical              :: file_exists
 integer              :: occfile
 integer              :: iter
!=====

 nstate = SIZE(occupation,DIM=1)

 if( temperature < 1.0e-8_dp ) then

   occupation(:,:)=0.0_dp

   inquire(file='manual_occupations',exist=file_exists)

   if(.NOT. file_exists) then
     remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
     if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)

     do istate=1,nstate
       occupation(istate,:) = MIN(remaining_electrons(:), spin_fact)
       remaining_electrons(:)  = remaining_electrons(:) - occupation(istate,:)
     enddo
   else
     write(stdout,*)
     write(stdout,*) 'occupations are read from file: manual_occupations'
     msg='reading manual occupations from file'
     call issue_warning(msg)
     open(newunit=occfile,file='manual_occupations',status='old')
     !
     ! read nlines, all other occupations are set to zero
     read(occfile,*) nlines
     do ilines=1,nlines
       read(occfile,*) occupation(ilines,:)
     enddo
     close(occfile)
     write(stdout,*) 'occupations set, closing file'
   endif

 else

   !
   ! Finite temperature case
   !
   write(stdout,'(1x,a,f12.6,3x,f15.3)') 'Find new the occupations and Fermi level for temperature (Ha) (K): ', &
                                         temperature,temperature * Ha_K

   ! First, set mu half way between the HOMO and the LUMO
   mu = 0.50_dp * ( energy(NINT(electrons/2.0_dp)+1,1) + energy(NINT(electrons/2.0_dp),1) )

   delta_mu = 1.0e-5_dp
   electrons_mu = -1.0_dp
   iter = 0
   mu_change = 0.0_dp

   do while( ABS( electrons_mu - electrons ) > 1.0e-8_dp .AND. iter <= 100 )

     iter = iter + 1
     mu = mu + mu_change

     occupation(:,:) = fermi_dirac(energy,mu)
     electrons_mu    = SUM( occupation(:,:) )

     grad_electrons = ( SUM( fermi_dirac(energy,mu+delta_mu) ) - SUM( fermi_dirac(energy,mu-delta_mu) ) ) / ( 2.0_dp* delta_mu )

     ! Maximum change is made bounded within +/- 0.10 Hartree
     mu_change = -( electrons_mu - electrons ) / grad_electrons
     mu_change = MAX( MIN( mu_change , 0.10_dp / REAL(iter) ), -0.10_dp / REAL(iter) )

!     write(*,*) iter,mu,mu_change,0.10_dp / REAL(iter),electrons_mu

   enddo

   write(stdout,'(1x,a,f12.6)') 'Fermi level (eV): ', mu * Ha_eV

 endif

 !
 ! final check
 if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0e-4_dp ) then
   write(stdout,*) 'occupation set up failed to give the right number of electrons'
   write(stdout,*) 'sum of occupations',SUM(occupation(:,:))
   write(stdout,*) 'electrons',electrons
   do istate=1,nstate
     write(stdout,*) istate,occupation(istate,:)
   enddo
   call die('FAILURE in set_occupation')
 endif

 call dump_out_occupation('=== Occupations ===',occupation)

contains

function fermi_dirac(energy,mu)
 implicit none
 real(dp),intent(in) :: energy(nstate,nspin)
 real(dp),intent(in) :: mu
 real(dp)            :: fermi_dirac(nstate,nspin)
!=====

 fermi_dirac(:,:) = spin_fact / ( 1.0_dp + EXP( ( energy(:,:) - mu ) / temperature ) )

end function fermi_dirac

end subroutine set_occupation


!=========================================================================
subroutine dump_out_occupation(title,occupation)
 implicit none
 character(len=*),intent(in) :: title
 real(dp),intent(in)         :: occupation(:,:)
!=====
 integer :: ihomo
 integer :: istate,nstate
!=====

 nstate = SIZE(occupation,DIM=1)

 write(stdout,'(/,1x,a)') TRIM(title)

 if( nspin == 2 ) then
   write(stdout,'(a)') '           spin 1       spin 2 '
 endif
 do istate=1,nstate
   if( ANY(occupation(istate,:) > 0.001_dp) ) ihomo = istate
 enddo

 select case(nspin)
 case(1)
   do istate=MAX(1,ihomo-5),MIN(ihomo+5,nstate)
     write(stdout,'(1x,i3,2(2(1x,f12.5)),2x)') istate,occupation(istate,1)
   enddo
 case(2)
   do istate=MAX(1,ihomo-5),MIN(ihomo+5,nstate)
     write(stdout,'(1x,i3,2(2(1x,f12.5)),2x)') istate,occupation(istate,1),occupation(istate,2)
   enddo
 end select
 write(stdout,*)

end subroutine dump_out_occupation


!=========================================================================
subroutine dump_out_energy(title,occupation,energy)
 implicit none
 character(len=*),intent(in) :: title
 real(dp),intent(in)         :: occupation(:,:),energy(:,:)
!=====
 integer,parameter :: MAXSIZE=300
 integer  :: istate,nocc,nstate
!=====

 nocc   = get_number_occupied_states(occupation)
 nstate = SIZE(occupation,DIM=1)

 write(stdout,'(/,1x,a)') TRIM(title)

 if(nspin==1) then
   write(stdout,'(a)') '   #       (Ha)         (eV)      '
 else
   write(stdout,'(a)') '   #              (Ha)                      (eV)      '
   write(stdout,'(a)') '           spin 1       spin 2       spin 1       spin 2'
 endif
 do istate=MAX(1,nocc-MAXSIZE/2),MIN(nstate,nocc+MAXSIZE/2)
   select case(nspin)
   case(1)
     write(stdout,'(1x,i3,2(1x,f12.5),4x,f8.4)') istate,energy(istate,1),energy(istate,1)*Ha_eV,occupation(istate,1)
   case(2)
     write(stdout,'(1x,i3,2(2(1x,f12.5)),4x,2(f8.4,2x))') istate,energy(istate,1),energy(istate,2), &
                                                          energy(istate,1)*Ha_eV,energy(istate,2)*Ha_eV, &
                                                          occupation(istate,1),occupation(istate,2)
   end select
   if(istate < nstate) then
     if( ANY( occupation(istate+1,:) < spin_fact/2.0_dp .AND. occupation(istate,:) > spin_fact/2.0_dp ) ) then
        if(nspin==1) then
          write(stdout,'(a)') '  -----------------------------'
        else
          write(stdout,'(a)') '  -------------------------------------------------------'
        endif
     endif
   endif
 enddo

 write(stdout,*)

end subroutine dump_out_energy


!=========================================================================
subroutine output_homolumo(calculation_name,occupation,energy,istate_min,istate_max)
 implicit none

 character(len=*),intent(in) :: calculation_name
 integer,intent(in)          :: istate_min,istate_max
 real(dp),intent(in)         :: occupation(:,:),energy(:,:)
!=====
 real(dp) :: ehomo_tmp,elumo_tmp
 real(dp) :: ehomo(nspin),elumo(nspin)
 integer  :: ispin,istate,nstate
!=====

 nstate = SIZE(occupation(:,:),DIM=1)

 do ispin=1,nspin
   ehomo_tmp=-HUGE(1.0_dp)
   elumo_tmp= HUGE(1.0_dp)

   do istate=istate_min,istate_max

     if( occupation(istate,ispin)/spin_fact > completely_empty ) then
       ehomo_tmp = MAX( ehomo_tmp , energy(istate,ispin) )
     endif

     if( occupation(istate,ispin)/spin_fact < 1.0_dp - completely_empty ) then
       elumo_tmp = MIN( elumo_tmp , energy(istate,ispin) )
     endif

   enddo

   ehomo(ispin) = ehomo_tmp
   elumo(ispin) = elumo_tmp

 enddo


 write(stdout,*)
 if( ALL( ehomo(:) > -1.0e6_dp ) ) then
   write(stdout,'(1x,a,1x,a,2(3x,f12.6))') TRIM(calculation_name),'HOMO energy    (eV):',ehomo(:) * Ha_eV
 endif
 if( ALL( elumo(:) <  1.0e6_dp ) ) then
   write(stdout,'(1x,a,1x,a,2(3x,f12.6))') TRIM(calculation_name),'LUMO energy    (eV):',elumo(:) * Ha_eV
 endif
 if( ALL( ehomo(:) > -1.0e6_dp ) .AND. ALL( elumo(:) <  1.0e6_dp ) ) then
   write(stdout,'(1x,a,1x,a,2(3x,f12.6))') TRIM(calculation_name),'HOMO-LUMO gap  (eV):',( elumo(:)-ehomo(:) ) * Ha_eV
 endif
 write(stdout,*)


end subroutine output_homolumo


!=========================================================================
subroutine matrix_ao_to_mo_diag(c_matrix,matrix_in,diag_out)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: matrix_in(:,:,:)
 real(dp),intent(out) :: diag_out(:,:)
!=====
 integer              :: nbf,nstate,nspin_ham
 integer              :: ispin,istate,ispin_ham
 real(dp),allocatable :: vector_tmp(:)
!=====

 nbf       = SIZE(c_matrix(:,:,:),DIM=1)
 nstate    = SIZE(c_matrix(:,:,:),DIM=2)
 nspin_ham = SIZE(matrix_in,DIM=3)


 allocate(vector_tmp(nbf))

 do ispin=1,nspin

   ispin_ham = MIN(ispin,nspin_ham)

   !matrix_inout(1:nstate,1:nstate,ispin) = MATMUL( TRANSPOSE( c_matrix(:,:,ispin) ) , MATMUL( matrix_inout(:,:,ispin) , c_matrix(:,:,ispin) ) )
   !diag_i =  DOT_PRODUCT( c_matrix_restart(:,istate,ispin) , MATMUL( hamiltonian_hartree(:,:) , c_matrix_restart(:,istate,ispin) ) )
   do istate=1,nstate

     ! H * C_i
     call DSYMV('L',nbf,1.0d0,matrix_in(1,1,ispin_ham),nbf, &
                              c_matrix(1,istate,ispin),1,  &
                        0.0d0,vector_tmp,1)
     ! C_i**T * (H * C_i)
     diag_out(istate,ispin) = DOT_PRODUCT( c_matrix(:,istate,ispin) , vector_tmp(:) )

   enddo
 enddo
 deallocate(vector_tmp)

end subroutine matrix_ao_to_mo_diag


!=========================================================================
subroutine matrix_ao_to_mo(c_matrix,matrix_in,matrix_out)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: matrix_in(:,:,:)
 real(dp),intent(out) :: matrix_out(:,:,:)
!=====
 integer                 :: nbf,nstate
 integer                 :: ispin
 real(dp),allocatable    :: matrix_tmp(:,:)
!=====

 nbf    = SIZE(c_matrix(:,:,:),DIM=1)
 nstate = SIZE(c_matrix(:,:,:),DIM=2)


 allocate(matrix_tmp(nbf,nstate))

 do ispin=1,nspin
   !matrix_inout(1:nstate,1:nstate,ispin) = MATMUL( TRANSPOSE( c_matrix(:,:,ispin) ) , MATMUL( matrix_inout(:,:,ispin) , c_matrix(:,:,ispin) ) )

   ! H * C
   call DSYMM('L','L',nbf,nstate,1.0d0,matrix_in(1,1,ispin),nbf, &
                                       c_matrix(1,1,ispin),nbf,  &
                                 0.0d0,matrix_tmp,nbf)

   ! C**T * (H * C)
#if defined(HAVE_MKL)
   call DGEMMT('L','T','N',nstate,nbf,1.0d0,c_matrix(1,1,ispin),nbf, &
                                            matrix_tmp,nbf,          &
                                      0.0d0,matrix_out(1,1,ispin),nstate)
   call matrix_lower_to_full(matrix_out(:,:,1))
#else
   call DGEMM('T','N',nstate,nstate,nbf,1.0d0,c_matrix(1,1,ispin),nbf, &
                                              matrix_tmp,nbf,          &
                                        0.0d0,matrix_out(1,1,ispin),nstate)
#endif

 enddo
 deallocate(matrix_tmp)

end subroutine matrix_ao_to_mo


!=========================================================================
subroutine matrix_mo_to_ao(c_matrix,matrix_in,matrix_out)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: matrix_in(:,:,:)
 real(dp),intent(out) :: matrix_out(:,:,:)
!=====
 integer              :: nbf,nstate
 integer              :: ispin
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 nbf    = SIZE(c_matrix(:,:,:),DIM=1)
 nstate = SIZE(c_matrix(:,:,:),DIM=2)


 allocate(matrix_tmp(nbf,nstate))

 do ispin=1,nspin
   !matrix_out(1:nbf,1:nbf,ispin) = MATMUL( c_matrix(:,:,ispin) , MATMUL( matrix_in(:,:,ispin) , TRANSPOSE( c_matrix(:,:,ispin) ) ) )

   !  C * H
   call DSYMM('R','L',nbf,nstate,1.0d0,matrix_in(1,1,ispin),nbf, &
                                       c_matrix(1,1,ispin),nbf,  &
                                 0.0d0,matrix_tmp,nbf)

   ! (C * H) * C**T
#if defined(HAVE_MKL)
   call DGEMMT('L','N','T',nbf,nstate,1.0d0,matrix_tmp,nbf,          &
                                            c_matrix(1,1,ispin),nbf, &
                                      0.0d0,matrix_out(1,1,ispin),nbf)
   call matrix_lower_to_full(matrix_out(:,:,1))
#else
   call DGEMM('N','T',nbf,nbf,nstate,1.0d0,matrix_tmp,nbf,        &
                                           c_matrix(1,1,ispin),nbf,  &
                                     0.0d0,matrix_out(1,1,ispin),nbf)
#endif

 enddo
 deallocate(matrix_tmp)

end subroutine matrix_mo_to_ao


!=========================================================================
subroutine evaluate_s2_operator(occupation,c_matrix,s_matrix)
 implicit none
 real(dp),intent(in)     :: occupation(:,:)
 real(dp),intent(in)     :: c_matrix(:,:,:)
 real(dp),intent(in)     :: s_matrix(:,:)
!=====
 integer                 :: nstate
 integer                 :: istate,jstate
 real(dp)                :: s2,s2_exact
 real(dp)                :: n1,n2,nmax,nmin
!=====

 if(nspin /= 2) return

 nstate = SIZE(occupation,DIM=1)

 n1 = SUM( occupation(:,1) )  ! Number of spin up   electrons
 n2 = SUM( occupation(:,2) )  ! Number of spin down electrons
 nmax = MAX(n1,n2)
 nmin = MIN(n1,n2)

 s2_exact = (nmax-nmin)/2.0_dp * ( (nmax-nmin)/2.0_dp + 1.0_dp )
 s2       = s2_exact + nmin
 do istate=1,nstate
   if( occupation(istate,1) < completely_empty ) cycle
   do jstate=1,nstate
     if( occupation(jstate,2) < completely_empty ) cycle

     s2 = s2 - ABS( DOT_PRODUCT( c_matrix(:,istate,1) , MATMUL( s_matrix(:,:) , c_matrix(:,jstate,2) ) )  &
                      * occupation(istate,1) * occupation(jstate,2) )**2

   enddo
 enddo


 write(stdout,'(/,a,f8.4)') ' Total Spin S**2: ',s2
 write(stdout,'(a,f8.4)')   ' Instead of:      ',s2_exact


end subroutine evaluate_s2_operator


!=========================================================================
subroutine level_shifting_up(s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
 implicit none
 real(dp),intent(in)    :: s_matrix(:,:)
 real(dp),intent(in)    :: c_matrix(:,:,:)
 real(dp),intent(in)    :: occupation(:,:)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(:,:,:)
!=====
 integer              :: nstate
 integer              :: ispin,istate
 real(dp),allocatable :: sqrt_level_shifting(:)
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 write(stdout,'(/,a)')     ' Level shifting switched on'
 write(stdout,'(a,f12.6)') '   energy shift (eV):',level_shifting_energy * Ha_eV

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif

 nstate = SIZE(occupation,DIM=1)
 allocate(matrix_tmp,MOLD=s_matrix)
 allocate(sqrt_level_shifting(nstate))

 do ispin=1,nspin
   !
   ! Shift up empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   !
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + matrix_tmp(:,:)

 enddo

 deallocate(matrix_tmp)
 deallocate(sqrt_level_shifting)


end subroutine level_shifting_up


!=========================================================================
subroutine level_shifting_down(s_matrix,c_matrix,occupation,level_shifting_energy,energy,hamiltonian)
 implicit none
 real(dp),intent(in)    :: s_matrix(:,:)
 real(dp),intent(in)    :: c_matrix(:,:,:)
 real(dp),intent(in)    :: occupation(:,:)
 real(dp),intent(inout) :: energy(:,:)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(:,:,:)
!=====
 integer              :: nstate
 integer              :: ispin,istate
 real(dp),allocatable :: sqrt_level_shifting(:)
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif

 nstate = SIZE(occupation,DIM=1)
 allocate(matrix_tmp,MOLD=s_matrix)
 allocate(sqrt_level_shifting(nstate))

 !
 ! Shift down the energies of the virtual orbitals
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       energy(istate,ispin) = energy(istate,ispin) - level_shifting_energy
     endif
   enddo
 enddo


 do ispin=1,nspin
   !
   ! Shift down empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   !
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) - matrix_tmp(:,:)

 enddo

 deallocate(matrix_tmp)
 deallocate(sqrt_level_shifting)

end subroutine level_shifting_down


!=========================================================================
subroutine setup_sqrt_overlap(TOL_OVERLAP,s_matrix,nstate,x_matrix)
 implicit none

 real(dp),intent(in)                :: TOL_OVERLAP
 real(dp),intent(in)                :: s_matrix(:,:)
 integer,intent(out)                :: nstate
 real(dp),allocatable,intent(inout) :: x_matrix(:,:)
!=====
 integer  :: nbf
 integer  :: istate,jbf
 real(dp),allocatable :: s_eigval(:)
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 write(stdout,'(/,a)') ' Calculate overlap matrix square-root S^{1/2}'

 nbf = SIZE(s_matrix,DIM=1)

 allocate(matrix_tmp(nbf,nbf))
 allocate(s_eigval(nbf))

 matrix_tmp(:,:) = s_matrix(:,:)
 ! Diagonalization with or without SCALAPACK
 call diagonalize_scalapack(scf_diago_flavor,scalapack_block_min,matrix_tmp,s_eigval)

 nstate = COUNT( s_eigval(:) > TOL_OVERLAP )

 !
 ! S**{-1} = X * X**T
 !
 call clean_allocate('Overlap X * X**H = S**-1',x_matrix,nbf,nstate)

 write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
 write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
 write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',TOL_OVERLAP
 write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',nbf

 istate = 0
 do jbf=1,nbf
   if( s_eigval(jbf) > TOL_OVERLAP ) then
     istate = istate + 1
     x_matrix(:,istate) = matrix_tmp(:,jbf) / SQRT( s_eigval(jbf) )
   endif
 enddo

 deallocate(matrix_tmp,s_eigval)

end subroutine setup_sqrt_overlap


!=========================================================================
subroutine setup_sqrt_density_matrix(p_matrix,p_matrix_sqrt,p_matrix_occ)
 implicit none

 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: p_matrix_sqrt(:,:,:)
 real(dp),intent(out) :: p_matrix_occ(:,:)
!=====
 integer              :: nbf
 integer              :: ispin,ibf
!=====

 nbf = SIZE( p_matrix(:,:,:), DIM=1 )

 write(stdout,*) 'Calculate the square root of the density matrix'
 call start_clock(timing_sqrt_density_matrix)

 do ispin=1,nspin
   p_matrix_sqrt(:,:,ispin) = p_matrix(:,:,ispin)
   ! Diagonalization with or without SCALAPACK
   call diagonalize_scalapack(scf_diago_flavor,scalapack_block_min,p_matrix_sqrt(:,:,ispin),p_matrix_occ(:,ispin))
   do ibf=1,nbf
     ! this is to avoid instabilities
     if( p_matrix_occ(ibf,ispin) < 1.0e-8_dp ) then
       p_matrix_occ(ibf,ispin)    = 0.0_dp
       p_matrix_sqrt(:,ibf,ispin) = 0.0_dp
     else
       p_matrix_sqrt(:,ibf,ispin) = p_matrix_sqrt(:,ibf,ispin) * SQRT( p_matrix_occ(ibf,ispin) )
     endif
   enddo
 enddo

 call stop_clock(timing_sqrt_density_matrix)

end subroutine setup_sqrt_density_matrix


!=========================================================================
subroutine get_c_matrix_from_p_matrix(p_matrix,c_matrix,occupation)
 implicit none

 real(dp),intent(in)              :: p_matrix(:,:,:)
 real(dp),allocatable,intent(out) :: c_matrix(:,:,:)
 real(dp),allocatable,intent(out) :: occupation(:,:)
!=====
 real(dp),allocatable :: p_matrix_sqrt(:,:,:),occupation_tmp(:,:)
 integer              :: nbf,nstate
 integer              :: ispin
!=====

 nbf    = SIZE( p_matrix(:,:,:), DIM=1 )
 allocate(p_matrix_sqrt(nbf,nbf,nspin))
 allocate(occupation_tmp(nbf,nspin))

 write(stdout,*) 'Calculate the square root of the density matrix to obtain the C matrix'
 call start_clock(timing_sqrt_density_matrix)

 ! Minus the p_matrix so that the eigenvalues are ordered from the largest to the lowest
 p_matrix_sqrt(:,:,:) = -p_matrix(:,:,:)

 do ispin=1,nspin

   ! Diagonalization with or without SCALAPACK
   call diagonalize_scalapack(scf_diago_flavor,scalapack_block_min,p_matrix_sqrt(:,:,ispin),occupation_tmp(:,ispin))

 enddo

 occupation_tmp(:,:) = -occupation_tmp(:,:)

 nstate = COUNT( ALL( occupation_tmp(:,:) > -0.0001_dp , DIM=2 ) )

 if( nstate /= nbf ) then
   call issue_warning('get_c_matrix_from_p_matrix: negative occupation numbers')
   write(stdout,'(1x,a,i4)')     'Number of negative eigenvalues: ',nbf - nstate
   write(stdout,'(1x,a,*(1x,es18.8))') 'Most negative eigenvalue: ',occupation_tmp(nbf,:)
 endif


 if( ALLOCATED(c_matrix) ) deallocate(c_matrix)
 if( ALLOCATED(occupation) ) deallocate(occupation)
 allocate(c_matrix(nbf,nstate,nspin))
 allocate(occupation(nstate,nspin))
 c_matrix(:,:,:) = p_matrix_sqrt(:,1:nstate,:)
 occupation(:,:) = occupation_tmp(:,:)


 deallocate(p_matrix_sqrt,occupation_tmp)


 call stop_clock(timing_sqrt_density_matrix)

end subroutine get_c_matrix_from_p_matrix


!=========================================================================
subroutine diagonalize_hamiltonian_scalapack(hamiltonian,x_matrix,energy,c_matrix)
 implicit none

 real(dp),intent(in)  :: hamiltonian(:,:,:)
 real(dp),intent(in)  :: x_matrix(:,:)
 real(dp),intent(out) :: c_matrix(:,:,:)
 real(dp),intent(out) :: energy(:,:)
!=====
 integer :: nspin_local,nbf,nstate
 integer :: mh,nh,mc,nc,ms,ns
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
#if defined(HAVE_SCALAPACK)
 integer :: cntxt
 integer :: rank_sca,nprocs_sca
 integer :: desch(NDEL),descc(NDEL),descs(NDEL)
#endif
 integer  :: ispin
 integer  :: ilocal,jlocal,iglobal,jglobal
 integer  :: m_small,n_small
 real(dp),allocatable :: h_small(:,:),h_small2(:,:)
 real(dp),allocatable :: ham_local(:,:),c_matrix_local(:,:),s_matrix_local(:,:)
!=====

 nbf         = SIZE(c_matrix,DIM=1)
 nstate      = SIZE(c_matrix,DIM=2)
 nspin_local = SIZE(c_matrix,DIM=3)

#if defined(HAVE_SCALAPACK)

 nprow = MIN(nprow_sd,nbf/scalapack_block_min)
 npcol = MIN(npcol_sd,nbf/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 if( nprow * npcol > 1 ) then
   write(stdout,'(1x,a,i4,a,i4)') 'Generalized diagonalization using SCALAPACK with a grid',nprow,' x ',npcol
   call BLACS_PINFO( rank_sca, nprocs_sca )

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

   c_matrix(:,:,:) = 0.0_dp

   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if(cntxt > 0 ) then

     mh = NUMROC(nbf   ,block_row,iprow,first_row,nprow)
     nh = NUMROC(nbf   ,block_col,ipcol,first_col,npcol)
     mc = NUMROC(nbf   ,block_row,iprow,first_row,nprow)
     nc = NUMROC(nstate,block_col,ipcol,first_col,npcol)
     ms = NUMROC(nstate,block_row,iprow,first_row,nprow)
     ns = NUMROC(nstate,block_col,ipcol,first_col,npcol)


     call DESCINIT(desch,nbf   ,nbf   ,block_row,block_col,first_row,first_col,cntxt,MAX(1,mh),info)
     call DESCINIT(descc,nbf   ,nstate,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)
     call DESCINIT(descs,nstate,nstate,block_row,block_col,first_row,first_col,cntxt,MAX(1,ms),info)


     allocate(ham_local(mh,nh))
     allocate(c_matrix_local(mc,nc))
     allocate(s_matrix_local(mc,nc))
     allocate(h_small(ms,ns))

     !
     ! Set up the local copy of x_matrix
     do jlocal=1,nc
       jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
       do ilocal=1,mc
         iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
         s_matrix_local(ilocal,jlocal) = x_matrix(iglobal,jglobal)
       enddo
     enddo



     do ispin=1,nspin_local
       write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
       call start_clock(timing_diago_hamiltonian)

       !
       ! Set up the local copy of hamiltonian
       do jlocal=1,nh
         jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
         do ilocal=1,mh
           iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
           ham_local(ilocal,jlocal) = hamiltonian(iglobal,jglobal,ispin)
         enddo
       enddo

!       h_small(:,:) = MATMUL( TRANSPOSE(x_matrix(:,:)) , &
!                                MATMUL( hamiltonian(:,:,ispin) , x_matrix(:,:) ) )

       !
       ! H_small = ^tS^{-1/2} H S^{-1/2}
       call PDGEMM('N','N',nbf,nstate,nbf,                &
                    1.0_dp,ham_local,1,1,desch,           &
                    s_matrix_local,1,1,descc,             &
                    0.0_dp,c_matrix_local,1,1,descc)

       call PDGEMM('T','N',nstate,nstate,nbf,             &
                    1.0_dp,s_matrix_local,1,1,descc,      &
                    c_matrix_local,1,1,descc,             &
                    0.0_dp,h_small,1,1,descs)



       call diagonalize_sca(scf_diago_flavor,h_small,descs,energy(:,ispin))


!       c_matrix(:,:,ispin) = MATMUL( x_matrix(:,:) , h_small(:,:) )

       !
       ! C = S^{-1/2} C_small
       call PDGEMM('N','N',nbf,nstate,nstate,             &
                    1.0_dp,s_matrix_local,1,1,descc,      &
                    h_small,1,1,descs,                    &
                    0.0_dp,c_matrix_local,1,1,descc)


       do jlocal=1,nc
         jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
         do ilocal=1,mc
           iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
           c_matrix(iglobal,jglobal,ispin) = c_matrix_local(ilocal,jlocal)
         enddo
       enddo


      ! Nullify the eigval array for all CPUs but one, so that the all_reduce
      ! operation in the end yields the correct value
      ! Of course, using a broadcast instead would be a better solution, but I'm so lazy
       if( rank_sca /= 0 ) energy(:,ispin) = 0.0_dp


       call stop_clock(timing_diago_hamiltonian)
     enddo

     deallocate(ham_local,c_matrix_local,s_matrix_local,h_small)

     call BLACS_GRIDEXIT( cntxt )

   else
     energy(:,:) = 0.0_dp
   endif


   ! Poor man distribution TODO replace by a broadcast
   call xsum_world(energy)
   call xsum_world(c_matrix)

 else ! only one proc selected
#endif

   allocate(h_small2(nstate,nstate))

   do ispin=1,nspin_local
     write(stdout,'(1x,a,i3)') 'Generalized diagonalization for spin: ',ispin
     call start_clock(timing_diago_hamiltonian)

     allocate(h_small(nbf,nstate))
     ! h_small(:,:) = MATMUL( TRANSPOSE(x_matrix(:,:)) , &
     !                          MATMUL( hamiltonian(:,:,ispin) , x_matrix(:,:) ) )

     ! H * U
     call DGEMM('N','N',nbf,nstate,nbf,1.0d0,hamiltonian(:,:,ispin),nbf, &
                                             x_matrix,nbf,      &
                                       0.0d0,h_small,nbf)
     ! U**T * (H * U)
     call DGEMM('T','N',nstate,nstate,nbf,1.0d0,x_matrix,nbf,  &
                                                h_small,nbf,            &
                                          0.0d0,h_small2,nstate)
     deallocate(h_small)

     ! H * C' = C' * E
     call diagonalize(scf_diago_flavor,h_small2,energy(:,ispin))

     !c_matrix(:,1:nstate,ispin) = MATMUL( x_matrix(:,:) , h_small2(:,:) )
     ! C = U * C'
     call DGEMM('N','N',nbf,nstate,nstate,1.0d0,x_matrix,nbf, &
                                                h_small2,nstate,       &
                                          0.0d0,c_matrix(:,:,ispin),nbf)


     call stop_clock(timing_diago_hamiltonian)
   enddo

   deallocate(h_small2)

#if defined(HAVE_SCALAPACK)
 endif
#endif

end subroutine diagonalize_hamiltonian_scalapack


!=========================================================================
end module m_hamiltonian_tools
!=========================================================================
