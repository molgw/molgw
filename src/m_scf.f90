!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the self-consistent field cycle methods (DIIS, simple mixing)
!
!=========================================================================
module m_scf
 use m_definitions
 use m_warning
 use m_memory
 use m_timing
 use m_mpi
 use m_inputparam
 use m_scalapack
 use m_tools,only: invert


 integer,private              :: nhistmax
 integer,private              :: nhist_current

 integer,private              :: nstate_scf,nbf_scf             ! Physical dimensions
                                                                ! Storage  dimensions
 integer,private              :: m_ham_scf,n_ham_scf            ! nbf    x nbf
 integer,private              :: m_c_scf,n_c_scf                ! nbf    x nstate
 integer,private              :: m_r_scf,n_r_scf                ! nstate x nstate

 real(dp),allocatable,private :: ham_hist(:,:,:,:)
 real(dp),allocatable,private :: res_hist(:,:,:,:)
 real(dp),allocatable,private :: p_matrix_hist(:,:,:,:)
 real(dp),allocatable,private :: a_matrix_hist(:,:)
 real(dp),allocatable,private :: p_dot_h_hist(:,:)

 logical,private              :: adiis_regime

 integer,private              :: iscf
 integer,private              :: desc_r(NDEL)

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
subroutine init_scf(m_ham,n_ham,m_c,n_c,nbf,nstate)
 implicit none
 integer,intent(in)  :: m_ham,n_ham,m_c,n_c,nbf,nstate
!=====

 adiis_regime = .FALSE.

 nbf_scf    = nbf
 nstate_scf = nstate
 
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
   allocate(a_matrix_hist(nhistmax,nhistmax))
 case('ADIIS')
   adiis_regime = .TRUE.
   nhistmax = npulay_hist
   allocate(a_matrix_hist(nhistmax,nhistmax))
   allocate(p_dot_h_hist(nhistmax,nhistmax))
 case default
   call die('mixing scheme not implemented')
 end select

 if( cntxt_ham > 0 ) then
   call init_desc('H',nstate,nstate,desc_r,m_r_scf,n_r_scf)
 else
   m_r_scf = 0
   n_r_scf = 0
 endif
 call xmax_local(m_r_scf)
 call xmax_local(n_r_scf)


 call clean_allocate('Hamiltonian history',ham_hist,m_ham,n_ham,nspin,nhistmax)
 call clean_allocate('Residual history',res_hist,m_r_scf,n_r_scf,nspin,nhistmax)
 if( mixing_scheme == 'ADIIS' ) then
   call clean_allocate('Density matrix history',p_matrix_hist,m_ham,n_ham,nspin,nhistmax)
 endif
 
end subroutine init_scf


!=========================================================================
subroutine destroy_scf()
 implicit none
!=====

 if(ALLOCATED(ham_hist))         call clean_deallocate('Hamiltonian history',ham_hist)
 if(ALLOCATED(res_hist))         call clean_deallocate('Residual history',res_hist)
 if(ALLOCATED(p_matrix_hist))    call clean_deallocate('Density matrix history',p_matrix_hist)
 if(ALLOCATED(a_matrix_hist))    deallocate(a_matrix_hist)
 if(ALLOCATED(p_dot_h_hist))     deallocate(p_dot_h_hist)

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
   ham_hist(:,:,:,ihist+1) = ham_hist(:,:,:,ihist)
   res_hist(:,:,:,ihist+1) = res_hist(:,:,:,ihist)
   a_matrix_hist(:,ihist+1) = a_matrix_hist(:,ihist)
   a_matrix_hist(ihist+1,:) = a_matrix_hist(ihist,:)

   if( mixing_scheme == 'ADIIS' ) then
     p_matrix_hist(:,:,:,ihist+1) = p_matrix_hist(:,:,:,ihist)
     p_dot_h_hist(:,ihist+1) = p_dot_h_hist(:,ihist)
     p_dot_h_hist(ihist+1,:) = p_dot_h_hist(ihist,:)
   endif

 enddo

 if( cntxt_ham > 0 ) then
   ham_hist(:,:,:,1) = ham(:,:,:)
   if( mixing_scheme == 'ADIIS' ) then
     p_matrix_hist(:,:,:,1) = p_matrix(:,:,:)
   endif
 endif

 ! Standard Pulay DIIS prediction here !
 call diis_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)

 ! If ADIIS prediction, overwrite the hamiltonian with a new one
 if( mixing_scheme == 'ADIIS' .AND. adiis_regime ) then
   call adiis_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)
 endif

end subroutine hamiltonian_prediction


!=========================================================================
subroutine simple_mixing_p_matrix(p_matrix_old,p_matrix_new)
 implicit none
 real(dp),intent(in)    :: p_matrix_old(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(inout) :: p_matrix_new(m_ham_scf,n_ham_scf,nspin)
!=====

 if( ABS(alpha_mixing-1.0_dp) < 1.0e-5_dp ) return

 if( mixing_scheme/='SIMPLE' .AND. iscf > mixing_first_nscf ) return

 write(stdout,'(/,1x,a,1x,f8.4)') 'Simple mixing of the density matrix with alpha_mixing:',alpha_mixing

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
 integer                :: ihist
 real(dp),allocatable   :: matrix_tmp1(:,:)
 real(dp),allocatable   :: matrix_tmp2(:,:)
 real(dp),allocatable   :: a_matrix(:,:)
 real(dp),allocatable   :: a_matrix_inv(:,:)
 real(dp),allocatable   :: alpha_diis(:)
 real(dp)               :: residual_pred(m_r_scf,n_r_scf,nspin)
 real(dp)               :: residual,work(1)
#ifdef HAVE_SCALAPACK
 real(dp),external      :: PDLANGE
#endif
!=====

 call start_clock(timing_diis)


 write(stdout,'(/,1x,a)') 'Pulay DIIS mixing'


 allocate(a_matrix(nhist_current+1,nhist_current+1))
 allocate(a_matrix_inv(nhist_current+1,nhist_current+1))
 allocate(alpha_diis(nhist_current))

 !
 ! Calculate the new residual as proposed in
 ! P. Pulay, J. Comput. Chem. 3, 554 (1982).
 !
 !  R =  U^T * [ H * P * S - S * P * H ] * U

 if( parallel_ham ) then

#ifdef HAVE_SCALAPACK
   if( cntxt_ham > 0 ) then

     do ispin=1,nspin

       allocate(matrix_tmp1(m_ham_scf,n_ham_scf))
       allocate(matrix_tmp2(m_ham_scf,n_ham_scf))

       ! M1 = H * P
       call PDGEMM('N','N',nbf_scf,nbf_scf,nbf_scf,1.0_dp,ham(:,:,ispin),1,1,desc_ham,          &
                   p_matrix(:,:,ispin),1,1,desc_ham,0.0_dp,matrix_tmp1,1,1,desc_ham)

       ! M2 = ( H * P ) * S
       call PDGEMM('N','N',nbf_scf,nbf_scf,nbf_scf,1.0_dp,matrix_tmp1,1,1,desc_ham,          &
                   s_matrix,1,1,desc_ham,0.0_dp,matrix_tmp2,1,1,desc_ham)

       ! M1 = S * P
       call PDGEMM('N','N',nbf_scf,nbf_scf,nbf_scf,1.0_dp,s_matrix,1,1,desc_ham,          &
                   p_matrix(:,:,ispin),1,1,desc_ham,0.0_dp,matrix_tmp1,1,1,desc_ham)

       ! M2 = M2 - ( S * P ) * H 
       call PDGEMM('N','N',nbf_scf,nbf_scf,nbf_scf,1.0_dp,matrix_tmp1,1,1,desc_ham,          &
                   ham(:,:,ispin),1,1,desc_ham,-1.0_dp,matrix_tmp2,1,1,desc_ham)

     
       deallocate(matrix_tmp1)
       allocate(matrix_tmp1(m_c_scf,n_c_scf))

       ! M1 = M2 * U
       call PDGEMM('N','N',nbf_scf,nstate_scf,nbf_scf,1.0_dp,matrix_tmp2,1,1,desc_ham,      &
                   s_matrix_sqrt_inv,1,1,desc_c,0.0_dp,matrix_tmp1,1,1,desc_c)

       ! R = U^T * M1
       call PDGEMM('T','N',nstate_scf,nstate_scf,nbf_scf,1.0_dp,s_matrix_sqrt_inv,1,1,desc_c,      &
                   matrix_tmp1,1,1,desc_c,0.0_dp,res_hist(:,:,ispin,1),1,1,desc_r)


       deallocate(matrix_tmp1,matrix_tmp2)
     enddo

   endif
#endif

 else

   allocate(matrix_tmp1(m_ham_scf,n_ham_scf))
   allocate(matrix_tmp2(m_ham_scf,n_ham_scf))

   do ispin=1,nspin

     !
     ! M1 = H * P * S
     call product_abc_scalapack(scalapack_block_min,ham(:,:,ispin),p_matrix(:,:,ispin),s_matrix,matrix_tmp1)
     !
     ! M2 = S * P * H
     call product_abc_scalapack(scalapack_block_min,s_matrix,p_matrix(:,:,ispin),ham(:,:,ispin),matrix_tmp2)

     ! M1 = M1 + M2
     matrix_tmp1(:,:) = matrix_tmp1(:,:) - matrix_tmp2(:,:)

     !
     ! R = U^T * M1 * U
     ! Remeber that S = U * U^T
     call product_transaba_scalapack(scalapack_block_min,s_matrix_sqrt_inv,matrix_tmp1,res_hist(:,:,ispin,1))

   enddo

   deallocate(matrix_tmp1,matrix_tmp2)

 endif

 !
 ! Build up a_matrix that contains the scalar product of residuals
 !

 !
 ! The older parts of a_matrix are saved in a_matrix_hist
 ! Just calculate the new ones
 if( parallel_ham ) then

   if( cntxt_ham > 0 ) then
     a_matrix_hist(1,1) = SUM( res_hist(:,:,:,1)**2 ) * nspin
     do ihist=2,nhist_current
       a_matrix_hist(ihist,1) = SUM( res_hist(:,:,:,ihist) * res_hist(:,:,:,1) ) * nspin
       a_matrix_hist(1,ihist) = a_matrix_hist(ihist,1)
     enddo
   else
     a_matrix_hist(1,1:nhist_current) = 0.0_dp
     a_matrix_hist(1:nhist_current,1) = 0.0_dp
   endif
   call xsum_world(a_matrix_hist(1,1))
   call xsum_world(a_matrix_hist(1,2:nhist_current))
   call xsum_world(a_matrix_hist(2:nhist_current,1))

 else

   a_matrix_hist(1,1) = SUM( res_hist(:,:,:,1)**2 ) * nspin

   do ihist=2,nhist_current
     a_matrix_hist(ihist,1) = SUM( res_hist(:,:,:,ihist) * res_hist(:,:,:,1) ) * nspin
     a_matrix_hist(1,ihist) = a_matrix_hist(ihist,1) 
   enddo

 endif

 a_matrix(1:nhist_current,1:nhist_current) = a_matrix_hist(1:nhist_current,1:nhist_current)


 !
 ! DIIS algorithm from Pulay (1980)
 ! 
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

 !
 ! Output the residual history and coefficients
 !
 write(stdout,'(a,4x,30(2x,es12.5))') '  Residuals:',( SQRT(a_matrix(ihist,ihist)) , ihist=1,nhist_current )
 write(stdout,'(a,30(2x,f12.6))') ' Alpha DIIS: ',alpha_diis(1:nhist_current)


 !
 ! Calculate the predicted hamiltonian
 !
 residual_pred(:,:,:) = 0.0_dp
 ham(:,:,:) = 0.0_dp
 if( parallel_ham ) then

#ifdef HAVE_SCALAPACK
   if( cntxt_ham > 0 ) then

     residual = 0.0_dp
     do ispin=1,nspin

       do ihist=1,nhist_current
         call PDGEADD('N',nstate_scf,nstate_scf,alpha_diis(ihist),res_hist(:,:,ispin,ihist),1,1,desc_r,1.0_dp,residual_pred(:,:,ispin),1,1,desc_r)
         call PDGEADD('N',nbf_scf,nbf_scf,alpha_diis(ihist),ham_hist(:,:,ispin,ihist),1,1,desc_ham,1.0_dp,ham(:,:,ispin),1,1,desc_ham)
       enddo

       residual = residual + PDLANGE('F',nstate_scf,nstate_scf,residual_pred(:,:,ispin),1,1,desc_r,work)**2
     enddo

   endif
   write(stdout,'(a,2x,es12.5,/)') ' DIIS predicted residual:',SQRT( residual * nspin )
#endif

 else

   do ihist=1,nhist_current
     residual_pred(:,:,:) = residual_pred(:,:,:) + alpha_diis(ihist) * res_hist(:,:,:,ihist)
     ham(:,:,:)           = ham(:,:,:)           + alpha_diis(ihist) * ham_hist(:,:,:,ihist) 
   enddo
   write(stdout,'(a,2x,es12.5,/)') ' DIIS predicted residual:',NORM2( residual_pred(:,:,:) ) * SQRT(REAL(nspin,dp))

 endif




 deallocate(a_matrix)
 deallocate(a_matrix_inv)
 deallocate(alpha_diis)


 call stop_clock(timing_diis)

end subroutine diis_prediction


!=========================================================================
subroutine adiis_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)
 use m_tools,only: matrix_trace,random
 use m_lbfgs
 implicit none
 real(dp),intent(in)    :: s_matrix(m_ham_scf,n_ham_scf)
 real(dp),intent(in)    :: s_matrix_sqrt_inv(m_c_scf,n_c_scf)
 real(dp),intent(in)    :: p_matrix(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(out)   :: ham(m_ham_scf,n_ham_scf,nspin)
!=====
 integer                :: ispin
 integer                :: ihist,khist
 real(dp),allocatable   :: matrix_tmp1(:,:)
 real(dp),allocatable   :: matrix_tmp2(:,:)
 real(dp),allocatable   :: a_matrix(:,:)
 real(dp),allocatable   :: a_matrix_inv(:,:)
 real(dp),allocatable   :: alpha_diis(:),alpha_diis_min(:)
 real(dp)               :: residual_pred(m_r_scf,n_r_scf,nspin)
 real(dp)               :: residual,work(1)
 real(dp)               :: ph_trace
 real(dp),allocatable   :: ph_matrix(:,:)
 real(dp),allocatable   :: ti(:),gradf(:),ci(:)
 real(dp)               :: sum_ti2
#ifdef HAVE_SCALAPACK
 real(dp),external      :: PDLANGE
#endif
 integer :: info
 integer :: iter
 integer,parameter :: niter=50 ! 000000
 real(dp) :: f_adiis,f_adiis_min
!=====

 call start_clock(timing_diis)


 write(stdout,'(/,1x,a)') 'ADIIS mixing'

 p_dot_h_hist(1,:) = 0.0_dp
 p_dot_h_hist(:,1) = 0.0_dp

 !
 ! Fill the p_dot_h matrix with the new row and column
 do ispin=1,nspin
   call trace_transab_scalapack(scalapack_block_min,p_matrix_hist(:,:,ispin,1),ham_hist(:,:,ispin,1),ph_trace)
   p_dot_h_hist(1,1) =  p_dot_h_hist(1,1) + ph_trace
 enddo
 write(stdout,*) 'FBFB < P | H >',p_dot_h_hist(1,1)

 do ihist=2,nhist_current
   do ispin=1,nspin
     call trace_transab_scalapack(scalapack_block_min,p_matrix_hist(:,:,ispin,ihist),ham_hist(:,:,ispin,1),ph_trace)
     p_dot_h_hist(ihist,1) =  p_dot_h_hist(ihist,1) + ph_trace
   enddo
 enddo

 do ihist=2,nhist_current
   do ispin=1,nspin
     call trace_transab_scalapack(scalapack_block_min,p_matrix_hist(:,:,ispin,1),ham_hist(:,:,ispin,ihist),ph_trace)
     p_dot_h_hist(1,ihist) =  p_dot_h_hist(1,ihist) + ph_trace
   enddo
 enddo


 allocate(alpha_diis(nhist_current))
 allocate(alpha_diis_min(nhist_current))
 allocate(ph_matrix(nhist_current,nhist_current))
 ph_matrix(:,:) = p_dot_h_hist(1:nhist_current,1:nhist_current)

#if 0
 f_adiis_min = HUGE(1.0_dp)
 do iter=1,niter

   do ihist=1,nhist_current
     alpha_diis(ihist)=random()
   enddo
   alpha_diis(:) = alpha_diis(:) / SUM( alpha_diis(:) )

   f_adiis =  DOT_PRODUCT( alpha_diis , ph_matrix(:,1) ) &
              + DOT_PRODUCT( alpha_diis , MATMUL( ph_matrix(:,:) , alpha_diis ) ) &
              - DOT_PRODUCT( ph_matrix(1,:) , alpha_diis )  &
              - ph_matrix(1,1)
   
   if( f_adiis < f_adiis_min ) then
     write(stdout,'(1x,i6,12(2x,f12.6))') iter,alpha_diis(1:nhist_current),f_adiis
     f_adiis_min = f_adiis
     alpha_diis_min(:) = alpha_diis(:)
   endif

 enddo

#else

 allocate(ti(nhist_current),ci(nhist_current))
 allocate(gradf(nhist_current))

 do ihist=1,nhist_current
   ti(ihist)=EXP( -REAL(ihist,dp) )
 enddo
 
 if( nhist_current > 1 ) then

   call setup_lbfgs(nhist_current)

   do iter=1,niter

     sum_ti2 = SUM( ti(:)**2 )
     ci(:) = ti(:)**2 / sum_ti2

     ! Evaluate function
     f_adiis =  DOT_PRODUCT( ci , ph_matrix(:,1) ) &
                + DOT_PRODUCT( ci , MATMUL( ph_matrix , ci ) ) &
                - DOT_PRODUCT( ph_matrix(1,:) , ci )  &
                - ph_matrix(1,1)


     do khist=1,nhist_current
       gradf(khist) = ph_matrix(khist,1) - DOT_PRODUCT( ti(:)**2 , ph_matrix(:,1) ) / sum_ti2  &
                     - ph_matrix(1,khist) + DOT_PRODUCT( ph_matrix(1,:) , ti(:)**2 ) / sum_ti2  &
                     + DOT_PRODUCT( ph_matrix(khist,:) , ti(:)**2 ) / sum_ti2  &
                     + DOT_PRODUCT( ti(:)**2 , ph_matrix(:,khist) ) / sum_ti2  &
                     - 2.0_dp * DOT_PRODUCT( ti(:)**2 , MATMUL( ph_matrix(:,:) , ti(:)**2 ) ) / sum_ti2**2
     enddo

     gradf(:) = gradf(:) * 2.0_dp * ti(:) / sum_ti2

     write(stdout,'(1x,i6,12(2x,f12.6))') iter,ci(:),f_adiis

     info = lbfgs_wrapper(ti,f_adiis,gradf)

     if( info <= 0 ) exit

   enddo

   call destroy_lbfgs()

 endif


 alpha_diis_min(:) = ti(:)**2 / SUM( ti(:)**2 )

 write(stdout,'(1x,a,12(2x,es14.6))') 'Final coefficients',ci(:)

 deallocate(ti,ci,gradf)

#endif

 
 ham(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   ham(:,:,:) = ham(:,:,:) + alpha_diis_min(ihist) * ham_hist(:,:,:,ihist) 
 enddo

 deallocate(alpha_diis)
 deallocate(alpha_diis_min)
 deallocate(ph_matrix)

 call stop_clock(timing_diis)

end subroutine adiis_prediction


!=========================================================================
function check_converged(p_matrix_old,p_matrix_new)
 implicit none

 logical               :: check_converged
 real(dp),intent(in)   :: p_matrix_old(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(in)   :: p_matrix_new(m_ham_scf,n_ham_scf,nspin)
!=====
 real(dp)              :: rms
!=====

 if( parallel_ham ) then
   if( cntxt_ham > 0 ) then
     rms = NORM2( p_matrix_new(:,:,:) - p_matrix_old(:,:,:) )**2
   else
     rms = 0.0_dp
   endif
   call xsum_world(rms)
 else
   rms = NORM2( p_matrix_new(:,:,:) - p_matrix_old(:,:,:) )**2
 endif


 rms = SQRT( rms * nspin )

 write(stdout,'(1x,a,2x,es12.5)') 'Convergence criterium on the density matrix',rms

 if( mixing_scheme == 'ADIIS' .AND. adiis_regime ) then
   if( rms < diis_switch ) then
     write(stdout,'(1x,a,es12.5)') 'Fair convergence has been reached: lower than ',diis_switch
     write(stdout,*) 'Now switch on regular DIIS'
     adiis_regime = .FALSE.
   endif
 endif

 if( rms < tolscf ) then 
   check_converged = .TRUE.
   write(stdout,*) ' ===> convergence has been reached'
   write(stdout,*)
 else
   check_converged = .FALSE.
   write(stdout,*) ' ===> convergence not reached yet'
   write(stdout,*)

   if(iscf == nscf) then
     if(rms>1.0e-2_dp) then
       call issue_warning('SCF convergence is very poor')
     else if(rms>1.0e-4_dp) then
       call issue_warning('SCF convergence is poor')
     endif
   endif

 endif


end function check_converged


!=========================================================================
end module m_scf
