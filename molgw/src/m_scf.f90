!=========================================================================
#include "macros.h"
!=========================================================================
module m_scf
 use m_definitions
 use m_mpi
 use m_inputparam


 integer,private              :: nhistmax
 integer,private              :: nhist_current
 integer,private              :: nbf_scf

 real(dp),allocatable,private :: p_matrix_in_hist(:,:,:,:)
 real(dp),allocatable,private :: residual_hist(:,:,:,:)

 integer,private              :: iscf

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
 type(energy_contributions) :: en


contains


!=========================================================================
subroutine init_scf(nbf)
 implicit none
 integer,intent(in)  :: nbf
!=====

 nbf_scf       = nbf
 iscf          = 1                 ! initialize with 1, since the new_p_matrix is not called for the first scf cycle
 nhist_current = 0

 select case(mixing_scheme)
 case('SIMPLE')
   nhistmax=1
 case('PULAY')
   nhistmax=8
 case default
   stop'mixing scheme not implemented'
 end select

 allocate(p_matrix_in_hist(nbf_scf,nbf_scf,nspin,nhistmax))
 allocate(residual_hist(nbf_scf,nbf_scf,nspin,nhistmax))
 
end subroutine init_scf


!=========================================================================
subroutine destroy_scf()
 implicit none
!=====

 deallocate(p_matrix_in_hist)
 if(allocated(residual_hist)) deallocate(residual_hist)

end subroutine destroy_scf


!=========================================================================
subroutine store_residual(p_matrix_in,p_matrix_out)
 implicit none
 real(dp),intent(in)  :: p_matrix_in(nbf_scf,nbf_scf,nspin)
 real(dp),intent(in)  :: p_matrix_out(nbf_scf,nbf_scf,nspin)
!=====
 integer              :: ihist
!=====

 !
 ! shift the old matrices and then store the new ones
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
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin)
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
     do while( ANY( ABS(alpha_diis(:)) > 4.0_dp ) .AND. nhist_current>3 )
       nhist_current = nhist_current - 1
       deallocate(alpha_diis)
       allocate(alpha_diis(nhist_current))
       call do_pulay_mixing(p_matrix_in,alpha_diis)
     enddo
     deallocate(alpha_diis)
   endif
 case default
   WRITE_MASTER(*,*) mixing_scheme
   stop'mixing scheme not implemented'
 end select

end subroutine new_p_matrix


!=========================================================================
subroutine do_simple_mixing(p_matrix_in)
 implicit none
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin)
!=====

 WRITE_MASTER(*,*) 'A simple mixing of the density matrix is used'

 p_matrix_in(:,:,:) = alpha_mixing * residual_hist(:,:,:,1) + p_matrix_in_hist(:,:,:,1)

end subroutine do_simple_mixing


!=========================================================================
subroutine do_pulay_mixing(p_matrix_in,alpha_diis)
 use m_tools,only:       invert
 implicit none
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin)
 real(dp),intent(out) :: alpha_diis(nhist_current)
!=====
 integer              :: ihist,jhist
 real(dp),allocatable :: amat(:,:),amat_inv(:,:)
 real(dp)             :: residual_pred(nbf_scf,nbf_scf,nspin)
!=====

 WRITE_MASTER(*,*) 'A Pulay mixing of the density matrix is used'

 allocate(amat    (nhist_current+1,nhist_current+1))
 allocate(amat_inv(nhist_current+1,nhist_current+1))

 !
 ! amat contains the scalar product of residuals
 do jhist=1,nhist_current
   do ihist=1,nhist_current
     amat(ihist,jhist) = SUM( residual_hist(:,:,:,ihist) * residual_hist(:,:,:,jhist) )
   enddo
 enddo
 amat(1:nhist_current,nhist_current+1) = -1.0_dp
 amat(nhist_current+1,1:nhist_current) = -1.0_dp
 amat(nhist_current+1,nhist_current+1) =  0.0_dp

 call invert(nhist_current+1,amat,amat_inv)

 alpha_diis(1:nhist_current) = -amat_inv(1:nhist_current,nhist_current+1)

 deallocate(amat,amat_inv)

 WRITE_MASTER(*,'(/,a,30(2x,f12.6))') ' alpha DIIS:',alpha_diis(1:nhist_current)
 
 residual_pred(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   residual_pred(:,:,:) = residual_pred(:,:,:) + alpha_diis(ihist) * residual_hist(:,:,:,ihist)
 enddo
 WRITE_MASTER(*,*) 'DIIS predicted residual',SQRT( SUM( residual_pred(:,:,:)**2 ) )
 WRITE_MASTER(*,*)

 p_matrix_in(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   p_matrix_in(:,:,:) = p_matrix_in(:,:,:) + alpha_diis(ihist) &
        * ( p_matrix_in_hist(:,:,:,ihist) + alpha_mixing * residual_hist(:,:,:,ihist) )
 enddo
 

end subroutine do_pulay_mixing


!=========================================================================
function check_converged()
 use m_definitions
 use m_warning
 implicit none
 logical               :: check_converged
!=====
 real(dp)              :: rms
!=====

 rms = SQRT( SUM( residual_hist(:,:,:,1)**2 ) )

 WRITE_MASTER(*,*) 'convergence criterium on the density matrix',rms
 if( rms < tolscf ) then 
   check_converged = .TRUE.
   WRITE_MASTER(*,*) ' ===> convergence has been reached'
 else
   check_converged = .FALSE.
   WRITE_MASTER(*,*) ' ===> convergence not reached yet'
 endif

 WRITE_MASTER(*,*)

 if(iscf == nscf) then
   if(rms>1.d-4) then
     msg='SCF convergence is very poor'
     call issue_warning(msg)
   else if(rms>1.d-6) then
     msg='SCF convergence is poor'
     call issue_warning(msg)
   endif
 endif

end function check_converged


!=========================================================================
end module m_scf
