!=========================================================================
#include "macros.h"
!=========================================================================
module m_scf
 use m_definitions

 private
 
 public :: simple_mixing,pulay_mixing,mixing_scheme,&
           init_scf,destroy_scf,store_residual,new_p_matrix,check_convergence

 integer,parameter    :: simple_mixing = 1
 integer,parameter    :: pulay_mixing  = 2

 integer              :: mixing_scheme

 integer              :: nhist
 integer              :: nhist_current
 integer              :: nbf_scf,nspin_scf                 ! internally saved data
 real(dp)             :: alpha_scf

 real(dp),allocatable :: p_matrix_in_hist(:,:,:,:)
 real(dp),allocatable :: residual_hist(:,:,:,:)

 integer              :: n_scf,n_scf_max


contains

!=========================================================================
subroutine init_scf(nscf,nbf,nspin,alpha_mixing)
 implicit none
 integer,intent(in)  :: nscf,nbf,nspin
 real(dp),intent(in) :: alpha_mixing
!=====

 n_scf_max     = nscf
 nbf_scf       = nbf
 nspin_scf     = nspin
 alpha_scf     = alpha_mixing
 n_scf         = 1                 ! initialize with 1, since the new_p_matrix is not called for the first scf cycle

 select case(mixing_scheme)
 case(simple_mixing)
   nhist=1
 case(pulay_mixing)
   nhist=6
 case default
   stop'mixing scheme not implemented'
 end select

 allocate(p_matrix_in_hist(nbf_scf,nbf_scf,nspin_scf,nhist))
 allocate(residual_hist(nbf_scf,nbf_scf,nspin_scf,nhist))
 
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
 real(dp),intent(in)  :: p_matrix_in(nbf_scf,nbf_scf,nspin_scf)
 real(dp),intent(in)  :: p_matrix_out(nbf_scf,nbf_scf,nspin_scf)
!=====
 integer              :: ihist
!=====

 !
 ! shift the old matrices and then store the new ones
 ! the newest is 1
 ! the oldest is nhist
 do ihist=nhist-1,1,-1
   residual_hist   (:,:,:,ihist+1) = residual_hist   (:,:,:,ihist)
   p_matrix_in_hist(:,:,:,ihist+1) = p_matrix_in_hist(:,:,:,ihist)
 enddo
 residual_hist   (:,:,:,1) = p_matrix_out(:,:,:) - p_matrix_in(:,:,:)
 p_matrix_in_hist(:,:,:,1) = p_matrix_in(:,:,:)

end subroutine store_residual
 
!=========================================================================
subroutine new_p_matrix(p_matrix_in)
 implicit none
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin_scf)
!=====

 n_scf          = n_scf + 1
 nhist_current  = MIN(nhist,n_scf-1)

 select case(mixing_scheme)
 case(simple_mixing)
   call do_simple_mixing(p_matrix_in)
 case(pulay_mixing)
   if(n_scf<=3) then ! for safety, just do simple mixing at the begining
     call do_simple_mixing(p_matrix_in)
   else
     call do_pulay_mixing(p_matrix_in)
   endif
 case default
   stop'mixing scheme not implemented'
 end select

end subroutine new_p_matrix

!=========================================================================
subroutine do_simple_mixing(p_matrix_in)
 implicit none
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin_scf)
!=====

 WRITE_MASTER(*,*) 'A simple mixing of the density matrix is used'

 p_matrix_in(:,:,:) = alpha_scf * residual_hist(:,:,:,1) + p_matrix_in_hist(:,:,:,1)

end subroutine do_simple_mixing

!=========================================================================
subroutine do_pulay_mixing(p_matrix_in)
 use m_tools,only:       invert
 implicit none
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin_scf)
!=====
 integer              :: ihist,jhist
 real(dp),allocatable :: amat(:,:),amat_inv(:,:)
 real(dp),allocatable :: alpha_diis(:)
 real(dp)             :: residual_pred(nbf_scf,nbf_scf,nspin_scf)
!=====

 WRITE_MASTER(*,*) 'A Pulay mixing of the density matrix is used'

 allocate(amat    (nhist_current+1,nhist_current+1))
 allocate(amat_inv(nhist_current+1,nhist_current+1))
 allocate(alpha_diis(nhist_current))

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
        * ( p_matrix_in_hist(:,:,:,ihist) + alpha_scf * residual_hist(:,:,:,ihist) )
 enddo
 
 deallocate(alpha_diis)

end subroutine do_pulay_mixing

!=========================================================================
subroutine check_convergence(scf_loop_converged)
 use m_definitions
 use m_warning
 implicit none
 logical,intent(out)   :: scf_loop_converged
!=====
 real(dp)              :: rms
!=====

 rms = SQRT( SUM( residual_hist(:,:,:,1)**2 ) )

 WRITE_MASTER(*,*) 'convergence criterium on the density matrix',rms
 if( rms < 1.0e-8_dp ) then 
   scf_loop_converged= .TRUE.
   WRITE_MASTER(*,*) ' ===> convergence has been reached'
 else
   scf_loop_converged= .FALSE.
   WRITE_MASTER(*,*) ' ===> convergence not reached yet'
 endif

 WRITE_MASTER(*,*)

 if(n_scf == n_scf_max) then
   if(rms>1.d-4) then
     msg='SCF convergence is very poor'
     call issue_warning(msg)
   else if(rms>1.d-6) then
     msg='SCF convergence is poor'
     call issue_warning(msg)
   endif
 endif

end subroutine check_convergence

!=========================================================================
end module m_scf
