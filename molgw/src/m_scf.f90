!=========================================================================
module m_scf
 use m_definitions
 use m_warning
 use m_memory
 use m_inputparam


 integer,private              :: nhistmax
 integer,private              :: nhist_current
 integer,private              :: m_ham_scf,n_ham_scf
 integer,private              :: m_c_scf,n_c_scf

 real(dp),allocatable,private :: p_matrix_in_hist(:,:,:,:)
 real(dp),allocatable,private :: residual_hist(:,:,:,:)

 real(dp),allocatable,private :: ham_hist(:,:,:,:)
 real(dp),allocatable,private :: res_hist(:,:,:,:)

 integer,private              :: iscf

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
subroutine init_scf(m_ham,n_ham,m_c,n_c)
 implicit none
 integer,intent(in)  :: m_ham,n_ham,m_c,n_c
!=====

 m_ham_scf     = m_ham
 n_ham_scf     = n_ham
 m_c_scf       = m_c
 n_c_scf       = n_c
 iscf          = 1                 ! initialize with 1, since the new_p_matrix is not called for the first scf cycle
 nhist_current = 0

 select case(mixing_scheme)
 case('SIMPLE')
   nhistmax = 1
 case('PULAY')
   nhistmax = npulay_hist
 case('DIIS')
   nhistmax = npulay_hist
 case default
   call die('mixing scheme not implemented')
 end select

 call clean_allocate('P matrix history',p_matrix_in_hist,m_ham,n_ham,nspin,nhistmax)
 call clean_allocate('Residual history',residual_hist,m_ham,n_ham,nspin,nhistmax)
 call clean_allocate('Hamiltonian history',ham_hist,m_ham,n_ham,nspin,nhistmax)
 call clean_allocate('Residual history',res_hist,m_c,n_c,nspin,nhistmax)
 
end subroutine init_scf


!=========================================================================
subroutine destroy_scf()
 implicit none
!=====

 if(ALLOCATED(p_matrix_in_hist)) call clean_deallocate('P matrix history',p_matrix_in_hist)
 if(ALLOCATED(residual_hist))    call clean_deallocate('Residual history',residual_hist)
 if(ALLOCATED(ham_hist))         call clean_deallocate('Hamiltonian history',ham_hist)

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

 !
 ! Shift the old matrices and then store the new ones
 ! the newest is 1
 ! the oldest is nhistmax
 do ihist=nhistmax-1,1,-1
   res_hist(:,:,:,ihist+1) = res_hist(:,:,:,ihist)
   ham_hist(:,:,:,ihist+1) = ham_hist(:,:,:,ihist)
 enddo
 ham_hist(:,:,:,1) = ham(:,:,:)


 select case(mixing_scheme)
 case('SIMPLE')
   ! New simple mixing prediction here !
   call simple_mixing_prediction(ham)
 case('PULAY')
   ! New DIIS prediction here !
   call diis_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)
 end select

end subroutine hamiltonian_prediction


!=========================================================================
subroutine simple_mixing_prediction(ham)
 implicit none
 real(dp),intent(out) :: ham(m_ham_scf,n_ham_scf,nspin)
!=====

 write(stdout,*) 'A simple mixing of the Hamiltonian is used'

 ham(:,:,:) = alpha_mixing * ham_hist(:,:,:,1) + (1.0_dp - alpha_mixing) * ham_hist(:,:,:,MIN(2,nhist_current))

end subroutine simple_mixing_prediction


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
 real(dp)               :: residual_pred(m_ham_scf,n_ham_scf,nspin)
!=====


 allocate(a_matrix(nhist_current+1,nhist_current+1))
 allocate(a_matrix_inv(nhist_current+1,nhist_current+1))
 allocate(alpha_diis(nhist_current))

 !
 ! Calcualte the residuals as proposed in
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

 write(stdout,'(/,a,30(2x,e12.6))') '  Residuals:',( SQRT(a_matrix(ihist,ihist)) , ihist=1,nhist_current )

 a_matrix(1:nhist_current,nhist_current+1) = -1.0_dp
 a_matrix(nhist_current+1,1:nhist_current) = -1.0_dp
 a_matrix(nhist_current+1,nhist_current+1) =  0.0_dp

 call invert(nhist_current+1,a_matrix,a_matrix_inv)

 alpha_diis(1:nhist_current) = -a_matrix_inv(1:nhist_current,nhist_current+1)

 write(stdout,'(/,a,30(2x,f12.6))') ' Alpha DIIS:',alpha_diis(1:nhist_current)

 residual_pred(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   residual_pred(:,:,:) = residual_pred(:,:,:) + alpha_diis(ihist) * res_hist(:,:,:,ihist)
 enddo
 write(stdout,*) 'DIIS predicted residual',SQRT( SUM( residual_pred(:,:,:)**2 ))
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
subroutine store_residual(p_matrix_in,p_matrix_out)
 implicit none
 real(dp),intent(in)  :: p_matrix_in(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(in)  :: p_matrix_out(m_ham_scf,n_ham_scf,nspin)
!=====
 integer              :: ihist
!=====

 !
 ! Shift the old matrices and then store the new ones
 ! the newest is 1
 ! the oldest is nhistmax
 do ihist=nhistmax-1,1,-1
   residual_hist   (:,:,:,ihist+1) = residual_hist   (:,:,:,ihist)
   p_matrix_in_hist(:,:,:,ihist+1) = p_matrix_in_hist(:,:,:,ihist)
 enddo
 residual_hist   (:,:,:,1) = p_matrix_out(:,:,:) - p_matrix_in(:,:,:)
 p_matrix_in_hist(:,:,:,1) = p_matrix_in(:,:,:)

end subroutine store_residual

 
!=========================================================================
subroutine new_p_matrix(p_matrix_in)
 implicit none
 real(dp),intent(out) :: p_matrix_in (m_ham_scf,n_ham_scf,nspin)
!=====
 real(dp),allocatable :: alpha_diis(:)
!=====

 iscf = iscf + 1
 nhist_current  = MIN(nhist_current+1,nhistmax) 

 select case(mixing_scheme)
 case('SIMPLE')
   call do_simple_mixing(p_matrix_in)
 case('PULAY')
   if( iscf <= 3 ) then ! for safety, just do simple mixing at the begining
     call do_simple_mixing(p_matrix_in)
   else
     allocate(alpha_diis(nhist_current))
     call do_pulay_mixing(p_matrix_in,alpha_diis)
     do while( ANY( ABS(alpha_diis(:)) > 4.0_dp ) .AND. nhist_current > 3 )
       nhist_current = nhist_current - 1
       deallocate(alpha_diis)
       allocate(alpha_diis(nhist_current))
       call do_pulay_mixing(p_matrix_in,alpha_diis)
     enddo
     deallocate(alpha_diis)
   endif
 case default
   call die(TRIM(mixing_scheme)//' mixing scheme not implemented')
 end select

end subroutine new_p_matrix


!=========================================================================
subroutine do_simple_mixing(p_matrix_in)
 implicit none
 real(dp),intent(out) :: p_matrix_in(m_ham_scf,n_ham_scf,nspin)
!=====

 write(stdout,*) 'A simple mixing of the density matrix is used'

 p_matrix_in(:,:,:) = alpha_mixing * residual_hist(:,:,:,1) + p_matrix_in_hist(:,:,:,1)

end subroutine do_simple_mixing


!=========================================================================
subroutine do_pulay_mixing(p_matrix_in,alpha_diis)
 use m_tools,only:       invert
 implicit none
 real(dp),intent(out) :: p_matrix_in(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(out) :: alpha_diis(nhist_current)
!=====
 integer  :: ihist,jhist
 real(dp) :: a_matrix(nhist_current+1,nhist_current+1)
 real(dp) :: a_matrix_inv(nhist_current+1,nhist_current+1)
 real(dp) :: residual_pred(m_ham_scf,n_ham_scf,nspin)
!=====

 write(stdout,*) 'A Pulay mixing of the density matrix is used'

 !
 ! a_matrix contains the scalar product of residuals
 do jhist=1,nhist_current
   do ihist=1,nhist_current
     a_matrix(ihist,jhist) = SUM( residual_hist(:,:,:,ihist) * residual_hist(:,:,:,jhist) )
   enddo
 enddo
! call xtrans_sum(a_matrix)

 a_matrix(1:nhist_current,nhist_current+1) = -1.0_dp
 a_matrix(nhist_current+1,1:nhist_current) = -1.0_dp
 a_matrix(nhist_current+1,nhist_current+1) =  0.0_dp

 call invert(nhist_current+1,a_matrix,a_matrix_inv)

 alpha_diis(1:nhist_current) = -a_matrix_inv(1:nhist_current,nhist_current+1)


 write(stdout,'(/,a,30(2x,f12.6))') ' alpha DIIS:',alpha_diis(1:nhist_current)
 
 residual_pred(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   residual_pred(:,:,:) = residual_pred(:,:,:) + alpha_diis(ihist) * residual_hist(:,:,:,ihist)
 enddo
 write(stdout,*) 'DIIS predicted residual',SQRT( SUM( residual_pred(:,:,:)**2 ) )
 write(stdout,*)

 p_matrix_in(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   p_matrix_in(:,:,:) = p_matrix_in(:,:,:) + alpha_diis(ihist) &
        * ( p_matrix_in_hist(:,:,:,ihist) + alpha_mixing * residual_hist(:,:,:,ihist) )
 enddo
 

end subroutine do_pulay_mixing


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

 write(stdout,*) 'convergence criterium on the density matrix',rms
 if( rms < tolscf ) then 
   check_converged = .TRUE.
   write(stdout,*) ' ===> convergence has been reached'
   write(stdout,*)
 else
   check_converged = .FALSE.
   write(stdout,*) ' ===> convergence not reached yet'
   write(stdout,*)

   if(iscf > nscf) then
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
