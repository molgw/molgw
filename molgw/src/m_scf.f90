!=========================================================================
! This file is part of MOLGW.
!=========================================================================
module m_scf
 use m_definitions
 use m_warning
 use m_memory
 use m_mpi
 use m_inputparam


 integer,private              :: nhistmax
 integer,private              :: nhist_current
 integer,private              :: m_ham_scf,n_ham_scf            ! nbf    x nbf
 integer,private              :: m_c_scf,n_c_scf                ! nbf    x nstate
 integer,private              :: m_r_scf,n_r_scf                ! nstate x nstate

 real(dp),allocatable,private :: ham_hist(:,:,:,:)
 real(dp),allocatable,private :: res_hist(:,:,:,:)

 integer,private              :: iscf
 integer,private              :: desc_r(ndel)

 type energy_contributions
   real(dp) :: nuc_nuc = 0.0_dp
   real(dp) :: kin     = 0.0_dp
   real(dp) :: nuc     = 0.0_dp
   real(dp) :: hart    = 0.0_dp
   real(dp) :: exx_hyb = 0.0_dp
   real(dp) :: exx     = 0.0_dp
   real(dp) :: xc      = 0.0_dp
   real(dp) :: se      = 0.0_dp      ! single-excitation contribution
   real(dp) :: mp2     = 0.0_dp
   real(dp) :: rpa     = 0.0_dp
   real(dp) :: gw      = 0.0_dp
   real(dp) :: tot     = 0.0_dp
 end type
 type(energy_contributions) :: en


contains


!=========================================================================
subroutine init_scf(m_ham,n_ham,m_c,n_c,nstate)
 implicit none
 integer,intent(in)  :: m_ham,n_ham,m_c,n_c,nstate
!=====

 m_ham_scf     = m_ham
 n_ham_scf     = n_ham
 m_c_scf       = m_c
 n_c_scf       = n_c
 iscf          = 0
 nhist_current = 0

 select case(mixing_scheme)
 case('SIMPLE')
   nhistmax = 0
 case('PULAY','DIIS')
   nhistmax = npulay_hist
 case default
   call die('mixing scheme not implemented')
 end select

 if( cntxt_ham > 0 ) then
   call init_desc('H',nstate,nstate,desc_r,m_r_scf,n_r_scf)
 else
   m_r_scf = 0
   n_r_scf = 0
 endif
 call xlocal_max(m_r_scf)
 call xlocal_max(n_r_scf)


 call clean_allocate('Hamiltonian history',ham_hist,m_ham,n_ham,nspin,nhistmax)
 call clean_allocate('Residual history',res_hist,n_r_scf,n_r_scf,nspin,nhistmax)
 
end subroutine init_scf


!=========================================================================
subroutine destroy_scf()
 implicit none
!=====

 if(ALLOCATED(ham_hist))         call clean_deallocate('Hamiltonian history',ham_hist)
 if(ALLOCATED(res_hist))         call clean_deallocate('Residual history',res_hist)

end subroutine destroy_scf


!=========================================================================
subroutine hamiltonian_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)
 implicit none
 real(dp),intent(in)    :: s_matrix(m_ham_scf,n_ham_scf)
 real(dp),intent(in)    :: s_matrix_sqrt_inv(m_c_scf,n_c_scf)
 real(dp),intent(in)    :: p_matrix(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(inout) :: ham(m_ham_scf,n_ham_scf,nspin)
!=====
 integer                :: ihist
!=====

 iscf = iscf + 1
 nhist_current  = MIN(nhist_current+1,nhistmax) 

 if( mixing_scheme=='SIMPLE') return

 !
 ! Shift the old matrices and then store the new ones
 ! the newest is 1
 ! the oldest is nhistmax
 do ihist=nhistmax-1,1,-1
   res_hist(:,:,:,ihist+1) = res_hist(:,:,:,ihist)
   ham_hist(:,:,:,ihist+1) = ham_hist(:,:,:,ihist)
 enddo
 ham_hist(:,:,:,1) = ham(:,:,:)


 ! New DIIS prediction here !
 call diis_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)

end subroutine hamiltonian_prediction


!=========================================================================
subroutine simple_mixing_p_matrix(p_matrix_old,p_matrix_new)
 implicit none
 real(dp),intent(out)   :: p_matrix_old(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(inout) :: p_matrix_new(m_ham_scf,n_ham_scf,nspin)
!=====

 if( ABS(alpha_mixing-1.0_dp) < 1.0e-5_dp ) return

 if( mixing_scheme/='SIMPLE' .AND. iscf > mixing_first_nscf ) return

 write(stdout,'(/,x,a,x,f8.4)') 'Simple mixing of the density matrix with alpha_mixing:',alpha_mixing

 p_matrix_new(:,:,:) = alpha_mixing * p_matrix_new(:,:,:) + (1.0_dp - alpha_mixing) * p_matrix_old(:,:,:)
 
end subroutine simple_mixing_p_matrix


!=========================================================================
subroutine diis_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)
 implicit none
 real(dp),intent(in)    :: s_matrix(m_ham_scf,n_ham_scf)
 real(dp),intent(in)    :: s_matrix_sqrt_inv(m_c_scf,n_c_scf)
 real(dp),intent(in)    :: p_matrix(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(inout) :: ham(m_ham_scf,n_ham_scf,nspin)
!=====
 integer                :: ispin
 integer                :: ihist,jhist
 real(dp)               :: matrix_tmp(m_ham_scf,n_ham_scf)
 real(dp),allocatable   :: a_matrix(:,:)
 real(dp),allocatable   :: a_matrix_inv(:,:)
 real(dp),allocatable   :: alpha_diis(:)
 real(dp)               :: residual_pred(m_r_scf,n_r_scf,nspin)
!=====

 write(stdout,'(/,x,a)') 'Pulay DIIS mixing'

 allocate(a_matrix(nhist_current+1,nhist_current+1))
 allocate(a_matrix_inv(nhist_current+1,nhist_current+1))
 allocate(alpha_diis(nhist_current))

 !
 ! Calculate the residuals as proposed in
 ! P. Pulay, J. Comput. Chem. 3, 554 (1982).
 do ispin=1,nspin

   matrix_tmp(:,:) = MATMUL( MATMUL( ham(:,:,ispin) , p_matrix(:,:,ispin) ) , s_matrix )  &
                      - MATMUL( s_matrix , MATMUL( p_matrix(:,:,ispin) , ham(:,:,ispin) ) )

   res_hist(:,:,ispin,1) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv) , MATMUL( matrix_tmp , s_matrix_sqrt_inv ) )

 enddo


 !
 ! a_matrix contains the scalar product of residuals
 do jhist=1,nhist_current
   do ihist=1,nhist_current
     a_matrix(ihist,jhist) = SUM( res_hist(:,:,:,ihist) * res_hist(:,:,:,jhist) )
   enddo
 enddo

 write(stdout,'(a,4x,30(2x,es12.5))') '  Residuals:',( SQRT(a_matrix(ihist,ihist)) , ihist=1,nhist_current )

 a_matrix(1:nhist_current,nhist_current+1) = -1.0_dp
 a_matrix(nhist_current+1,1:nhist_current) = -1.0_dp
 a_matrix(nhist_current+1,nhist_current+1) =  0.0_dp

 call invert(nhist_current+1,a_matrix,a_matrix_inv)

 alpha_diis(1:nhist_current) = -a_matrix_inv(1:nhist_current,nhist_current+1)

 ! Renormalize the coefficients
 ! It should not be needed in principle, but sometimes it is
 if( ABS( SUM(alpha_diis(1:nhist_current)) -1.0_dp ) > 1.0e-4_dp ) then
   call issue_warning('DIIS coefficients rescaled')
   alpha_diis(1:nhist_current) = alpha_diis(1:nhist_current) / SUM( alpha_diis(1:nhist_current) )
 endif

 write(stdout,'(a,30(2x,f12.6))') ' Alpha DIIS: ',alpha_diis(1:nhist_current)

 residual_pred(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   residual_pred(:,:,:) = residual_pred(:,:,:) + alpha_diis(ihist) * res_hist(:,:,:,ihist)
 enddo
 write(stdout,'(a,2x,es12.5)') ' DIIS predicted residual:',SQRT( SUM( residual_pred(:,:,:)**2 ))
 write(stdout,*)


 ham(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   ham(:,:,:) = ham(:,:,:) + alpha_diis(ihist) * ham_hist(:,:,:,ihist) 
 enddo


 deallocate(a_matrix)
 deallocate(a_matrix_inv)
 deallocate(alpha_diis)



end subroutine diis_prediction


!=========================================================================
function check_converged(p_matrix_old,p_matrix_new)
 implicit none

 logical               :: check_converged
 real(dp),intent(in)   :: p_matrix_old(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(in)   :: p_matrix_new(m_ham_scf,n_ham_scf,nspin)
!=====
 real(dp)              :: rms
!=====

 rms = SQRT( SUM( ( p_matrix_new(:,:,:) - p_matrix_old(:,:,:) )**2 ) )

 call xtrans_sum(rms)

 write(stdout,'(x,a,2x,es12.5)') 'Convergence criterium on the density matrix',rms
 if( rms < tolscf ) then 
   check_converged = .TRUE.
   write(stdout,*) ' ===> convergence has been reached'
   write(stdout,*)
 else
   check_converged = .FALSE.
   write(stdout,*) ' ===> convergence not reached yet'
   write(stdout,*)

   if(iscf == nscf) then
     if(rms>1.d-3) then
       call issue_warning('SCF convergence is very poor')
     else if(rms>1.d-5) then
       call issue_warning('SCF convergence is poor')
     endif
   endif

 endif

end function check_converged


!=========================================================================
end module m_scf
