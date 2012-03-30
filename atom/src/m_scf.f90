!=========================================================================
module m_scf
 use m_definitions

 integer,parameter    :: simple_mixing = 1
 integer,parameter    :: rmdiis        = 2

 integer              :: mixing_scheme
 integer              :: nhist
 integer              :: nbf_scf,nspin_scf                 ! internally saved data
 real(dp)             :: alpha_scf

 real(dp),allocatable :: p_matrix_hist(:,:,:,:)
 real(dp),allocatable :: residual_hist(:,:,:,:)

 integer              :: n_scf


contains

!=========================================================================
subroutine init_scf(nbf,nspin,mixing_scheme_set,alpha_mixing)
 implicit none
 integer,intent(in)  :: nbf,nspin,mixing_scheme_set
 real(dp),intent(in) :: alpha_mixing
!=====

 nbf_scf       = nbf
 nspin_scf     = nspin
 mixing_scheme = mixing_scheme_set
 alpha_scf     = alpha_mixing
 n_scf         = 0

 select case(mixing_scheme)
 case(simple_mixing)
   nhist=1
 case(rmdiis)
   nhist=4
 case default
   stop'mixing scheme not implemented'
 end select

 allocate(p_matrix_hist(nbf_scf,nbf_scf,nspin_scf,nhist))
 
end subroutine init_scf

!=========================================================================
subroutine destroy_scf()
 implicit none
!=====

 deallocate(p_matrix_hist)
 if(allocated(residual_hist)) deallocate(residual_hist)

end subroutine destroy_scf

!=========================================================================
subroutine new_p_matrix(p_matrix_out,p_matrix_in)
 implicit none
 real(dp),intent(in)  :: p_matrix_out(nbf_scf,nbf_scf,nspin_scf)
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin_scf)
!=====

 n_scf = n_scf + 1

 !
 ! special treatment for the first iteration
 !
 if( n_scf==1 ) then
   p_matrix_hist(:,:,:,1) = p_matrix_out(:,:,:)
 endif

 select case(mixing_scheme)
 case(simple_mixing)
   call do_simple_mixing(p_matrix_out,p_matrix_in)
 case(rmdiis)
   call do_rmdiis(p_matrix_out,p_matrix_in)
 case default
   stop'mixing scheme not implemented'
 end select

end subroutine new_p_matrix

!=========================================================================
subroutine do_simple_mixing(p_matrix_out,p_matrix_in)
 implicit none
 real(dp),intent(in)  :: p_matrix_out(nbf_scf,nbf_scf,nspin_scf)
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin_scf)
!=====

 p_matrix_in(:,:,:) = alpha_scf * p_matrix_out(:,:,:) + ( 1.0_dp - alpha_scf ) * p_matrix_hist(:,:,:,1)

 p_matrix_hist(:,:,:,1) = p_matrix_in(:,:,:)

end subroutine do_simple_mixing

!=========================================================================
subroutine do_rmdiis(p_matrix_out,p_matrix_in)
 implicit none
 real(dp),intent(in)  :: p_matrix_out(nbf_scf,nbf_scf,nspin_scf)
 real(dp),intent(out) :: p_matrix_in (nbf_scf,nbf_scf,nspin_scf)
!=====

 p_matrix_in(:,:,:) = alpha_scf * p_matrix_out(:,:,:) + ( 1.0_dp - alpha_scf ) * p_matrix_hist(:,:,:,1)

 p_matrix_hist(:,:,:,1) = p_matrix_in(:,:,:)

end subroutine do_rmdiis

!=========================================================================
subroutine check_convergence(p_matrix_old,p_matrix,rms)
 use m_definitions
 implicit none
 real(dp),intent(in)  :: p_matrix_old(nbf_scf,nbf_scf,nspin_scf),p_matrix(nbf_scf,nbf_scf,nspin_scf)
 real(dp),intent(out) :: rms
!=====

 rms = SQRT( SUM( ( p_matrix_old(:,:,:) - p_matrix(:,:,:) )**2 ) )

end subroutine check_convergence

!=========================================================================
end module m_scf
