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
 integer,private              :: nhistmax_pmatrix
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
 real(dp),allocatable,private :: en_hist(:)
 real(dp),allocatable,private :: alpha_diis(:)

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
   real(dp) :: mp3     = 0.0_dp
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
   nhistmax         = 2
   nhistmax_pmatrix = 1

 case('PULAY','DIIS')
   nhistmax         = npulay_hist
   nhistmax_pmatrix = 1
   allocate(a_matrix_hist(nhistmax,nhistmax))

 case('ADIIS','EDIIS')
   adiis_regime = .TRUE.
   nhistmax         = npulay_hist
   nhistmax_pmatrix = npulay_hist
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
 call clean_allocate('Density matrix history',p_matrix_hist,m_ham,n_ham,nspin,nhistmax_pmatrix)
 if( mixing_scheme /= 'SIMPLE' ) then 
   call clean_allocate('Residual history',res_hist,m_r_scf,n_r_scf,nspin,nhistmax)
 endif

 allocate(en_hist(nhistmax))

 
end subroutine init_scf


!=========================================================================
subroutine destroy_scf()
 implicit none
!=====

 if(ALLOCATED(en_hist))          deallocate(en_hist)
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
 real(dp),intent(inout) :: p_matrix(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(inout) :: ham(m_ham_scf,n_ham_scf,nspin)
!=====
!=====

 iscf = iscf + 1
 nhist_current  = MIN(nhist_current+1,nhistmax) 

 allocate(alpha_diis(nhist_current))

 !
 ! Shift the old matrices and then store the new ones
 ! the newest is 1
 ! the oldest is nhistmax
 !
 en_hist(2:)          = en_hist(1:nhistmax-1)
 ham_hist(:,:,:,2:)   = ham_hist(:,:,:,1:nhistmax-1)
 if( ALLOCATED(res_hist) ) then
   res_hist(:,:,:,2:)   = res_hist(:,:,:,1:nhistmax-1)
 endif
 if( ALLOCATED(a_matrix_hist) ) then
   a_matrix_hist(2:,2:) = a_matrix_hist(1:nhistmax-1,1:nhistmax-1)
   a_matrix_hist(1,:)   = 0.0_dp
   a_matrix_hist(:,1)   = 0.0_dp
 endif

 if( mixing_scheme == 'ADIIS' .OR. mixing_scheme == 'EDIIS' ) then
   p_matrix_hist(:,:,:,2:) = p_matrix_hist(:,:,:,1:nhistmax_pmatrix-1)

   p_dot_h_hist(2:,2:) = p_dot_h_hist(1:nhistmax_pmatrix-1,1:nhistmax_pmatrix-1)
   p_dot_h_hist(1,:)   = 0.0_dp
   p_dot_h_hist(:,1)   = 0.0_dp
 endif

 !
 ! Set the newest values in history
 ! if already available here
 en_hist(1) = en%tot
 if( cntxt_ham > 0 ) then
   ham_hist(:,:,:,1)      = ham(:,:,:)
   p_matrix_hist(:,:,:,1) = p_matrix(:,:,:)
 endif


 ! Standard Pulay DIIS prediction here !
 if( mixing_scheme /= 'SIMPLE' ) then
   call diis_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,ham)
 else
   call simple_prediction(ham)
 endif


 ! If ADIIS or EDIIS prediction, overwrite the hamiltonian with a new one
 if( ( mixing_scheme == 'ADIIS' .OR. mixing_scheme == 'EDIIS' ) .AND. adiis_regime ) then
   call xdiis_prediction(p_matrix,ham)
 endif

 deallocate(alpha_diis)

end subroutine hamiltonian_prediction


!=========================================================================
subroutine simple_prediction(ham)
 implicit none
 real(dp),intent(inout) :: ham(m_ham_scf,n_ham_scf,nspin)
!=====
!=====

 if( nhist_current >= 2 ) then
   alpha_diis(1) = alpha_mixing
   alpha_diis(2) = 1.0_dp - alpha_mixing

   write(stdout,'(/,1x,a,f8.4)') 'Simple mixing with parameter:',alpha_mixing
   ham(:,:,:) = alpha_mixing * ham_hist(:,:,:,1) + (1.0_dp - alpha_mixing) * ham_hist(:,:,:,2)

 endif
 ham_hist(:,:,:,1) = ham(:,:,:)

end subroutine simple_prediction


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
     call matmul_abc_scalapack(scalapack_block_min,ham(:,:,ispin),p_matrix(:,:,ispin),s_matrix,matrix_tmp1)
     !
     ! M2 = S * P * H
     call matmul_abc_scalapack(scalapack_block_min,s_matrix,p_matrix(:,:,ispin),ham(:,:,ispin),matrix_tmp2)

     ! M1 = M1 + M2
     matrix_tmp1(:,:) = matrix_tmp1(:,:) - matrix_tmp2(:,:)

     !
     ! R = U^T * M1 * U
     ! Remeber that S = U * U^T
     call matmul_transaba_scalapack(scalapack_block_min,s_matrix_sqrt_inv,matrix_tmp1,res_hist(:,:,ispin,1))

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

 call invert(a_matrix,a_matrix_inv)

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

   else
     residual = -1.0_dp
   endif
   call xmax_world(residual)
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


 call stop_clock(timing_diis)

end subroutine diis_prediction


!=========================================================================
subroutine xdiis_prediction(p_matrix,ham)
 use m_lbfgs
 implicit none
 real(dp),intent(out)   :: p_matrix(m_ham_scf,n_ham_scf,nspin)
 real(dp),intent(out)   :: ham(m_ham_scf,n_ham_scf,nspin)
!=====
 type(lbfgs_state)      :: lbfgs_plan
 integer                :: ispin
 integer                :: ihist,jhist,khist
 real(dp),allocatable   :: alpha_diis_mc(:)
 real(dp)               :: ph_trace
 real(dp),allocatable   :: half_ph(:,:)
 real(dp),allocatable   :: ti(:),gradf(:),ci(:),dcdt(:,:),diag(:)
 real(dp)               :: sum_ti2
 integer :: info,iproc
 integer :: imc,ibfgs
 integer :: nseed,iseed
 integer,allocatable :: seed(:)
 integer,parameter :: nmc = 1000000
 integer,parameter :: nbfgs = 20
 real(dp) :: f_xdiis,f_xdiis_min
 real(dp),parameter :: alpha_max=0.60_dp
#ifdef HAVE_SCALAPACK
 real(dp),external      :: PDLANGE
#endif
!=====

 call start_clock(timing_diis)


 write(stdout,'(/,1x,a)') TRIM(mixing_scheme)//' mixing'

 p_dot_h_hist(1,:) = 0.0_dp
 p_dot_h_hist(:,1) = 0.0_dp

 do ihist=1,nhist_current
   do ispin=1,nspin
     call trace_transab_scalapack(scalapack_block_min,p_matrix_hist(:,:,ispin,ihist),ham_hist(:,:,ispin,1),ph_trace)
     p_dot_h_hist(ihist,1) =  p_dot_h_hist(ihist,1) + ph_trace
   enddo
 enddo

 do jhist=2,nhist_current
   do ispin=1,nspin
     call trace_transab_scalapack(scalapack_block_min,p_matrix_hist(:,:,ispin,1),ham_hist(:,:,ispin,jhist),ph_trace)
     p_dot_h_hist(1,jhist) =  p_dot_h_hist(1,jhist) + ph_trace
   enddo
 enddo


 allocate(alpha_diis_mc(nhist_current))
 allocate(half_ph(nhist_current,nhist_current))
 half_ph(:,:) = p_dot_h_hist(1:nhist_current,1:nhist_current) * 0.5_dp


 allocate(diag(nhist_current))

 do ihist=1,nhist_current
   diag(ihist) = en_hist(ihist) - half_ph(ihist,ihist)
 enddo


 allocate(ti(nhist_current),ci(nhist_current))
 allocate(dcdt(nhist_current,nhist_current))
 allocate(gradf(nhist_current))

 ci(1)               = 1.0_dp 
 ci(2:nhist_current) = 0.0_dp
 alpha_diis(:)  = ci(:)
 f_xdiis_min = eval_f_xdiis(ci)
   
 if( nhist_current > 1 ) then

   do ihist=1,nhist_current
     ti(1)  = 1.0_dp
     ti(2:) = 0.2_dp
   enddo

   call lbfgs_init(lbfgs_plan,nhist_current,5)

   do ibfgs=1,nbfgs

     sum_ti2 = SUM( ti(:)**2 )
     ci(:) = ti(:)**2 / sum_ti2

     do jhist=1,nhist_current
       do ihist=1,nhist_current
         dcdt(ihist,jhist) = - 2.0_dp * ti(ihist)**2 * ti(jhist) / sum_ti2**2
       enddo
       dcdt(jhist,jhist) = dcdt(jhist,jhist) + 2.0_dp * ti(jhist) / sum_ti2 
     enddo
  

     ! Evaluate XDIIS function
     f_xdiis =  eval_f_xdiis(ci)

     gradf(:) = eval_gradf_xdiis(ci,dcdt) 

     ! Perform a LBGS step
     info = lbfgs_execute(lbfgs_plan,ti,f_xdiis,gradf)

     !
     ! If the coefficient ci are identical within 1.0e-4, then consider they are converged
     if( ALL( ABS(ci(:) - ti(:)**2 / SUM( ti(:)**2 ) ) < 1.0e-4_dp ) ) then
        exit
     endif

     if( info <= 0 ) exit

   enddo

   call lbfgs_destroy(lbfgs_plan)

   sum_ti2 = SUM( ti(:)**2 )
   ci(:) = ti(:)**2 / sum_ti2
   f_xdiis = eval_f_xdiis(ci)

   if( f_xdiis < f_xdiis_min ) then
     alpha_diis(:) = ci(:)
     f_xdiis_min = f_xdiis
   endif


   ! If a coefficient is too large, start again the minimization 
   if( ANY( alpha_diis(2:) > alpha_max ) ) then

     !
     ! Find the offender
     khist = MAXLOC(alpha_diis(2:),DIM=1) + 1
     write(stdout,'(1x,a,i4,1x,f12.6)') 'Performing a sub-optimal XDIIS because one too-old coefficient is too large: ', &
                                        khist,alpha_diis(khist)

     call RANDOM_SEED(SIZE=nseed)
     allocate(seed(nseed))
     do iseed=1,nseed
       seed(iseed) = NINT( rank_world * iseed * pi * 27.21 )
     enddo
     call RANDOM_SEED(PUT=seed)
     deallocate(seed)

     alpha_diis(1)   = alpha_max
     alpha_diis(2:)  = (1.0_dp - alpha_max) / REAL(nhist_current-1,dp)
     f_xdiis_min = eval_f_xdiis(alpha_diis)
    
     do imc=1,nmc
       if( MODULO( imc - 1 , nproc_world ) /= rank_world ) cycle
     
       ! Find random coefficients that keep alpha_k = alpha_max and that sum up to 1
       call RANDOM_NUMBER(alpha_diis_mc)
       alpha_diis_mc(khist) = alpha_max
       sum_ti2 = ( SUM( alpha_diis_mc(:khist-1) ) + SUM( alpha_diis_mc(khist+1:) ) )
       alpha_diis_mc(:khist-1) = alpha_diis_mc(:khist-1) / sum_ti2 * (1.0_dp - alpha_max)
       alpha_diis_mc(khist+1:) = alpha_diis_mc(khist+1:) / sum_ti2 * (1.0_dp - alpha_max)
    
       f_xdiis = eval_f_xdiis(alpha_diis_mc)
    
       if( f_xdiis < f_xdiis_min ) then
         f_xdiis_min = f_xdiis
         alpha_diis(:) = alpha_diis_mc(:)
       endif
    
     enddo

     ! Propage f_xdiis_min and alpha_diis to all procs
     f_xdiis = f_xdiis_min
     call xmin_world(f_xdiis_min)

     if( ABS( f_xdiis_min - f_xdiis ) < 1.0e-14_dp ) then
       iproc = rank_world
     else
       iproc = -1
     endif
     call xmax_world(iproc)
     call xbcast_world(iproc,alpha_diis)


   endif


 endif


 deallocate(ti,ci,gradf)
 deallocate(dcdt)
 deallocate(diag)

 write(stdout,'(1x,a,12(2x,f16.6))') TRIM(mixing_scheme)//' final coefficients:',alpha_diis(:)
 write(stdout,'(1x,a,12(2x,f16.8))') 'Total energy history:    ',en_hist(1:nhist_current)
 write(stdout,'(1x,a,12(2x,f16.8))') TRIM(mixing_scheme)//' final energy:      ',f_xdiis_min

 ham(:,:,:)      = 0.0_dp
 p_matrix(:,:,:) = 0.0_dp
 do ihist=1,nhist_current
   ham(:,:,:)      = ham(:,:,:)      + alpha_diis(ihist) * ham_hist(:,:,:,ihist) 
   p_matrix(:,:,:) = p_matrix(:,:,:) + alpha_diis(ihist) * p_matrix_hist(:,:,:,ihist) 
 enddo

 deallocate(alpha_diis_mc)
 deallocate(half_ph)

 call stop_clock(timing_diis)


contains


function eval_f_xdiis(xi)
 real(dp),intent(in) :: xi(nhist_current)
 real(dp) :: eval_f_xdiis
!=====
 
 select case(mixing_scheme)
 case('EDIIS')
   eval_f_xdiis = DOT_PRODUCT( xi , MATMUL( half_ph , xi ) ) + DOT_PRODUCT( xi , diag )
 case('ADIIS')
   eval_f_xdiis =  en_hist(1) - 2.0_dp * half_ph(1,1)  &
                + 2.0_dp * DOT_PRODUCT( xi , half_ph(:,1) ) &
                - 2.0_dp * DOT_PRODUCT( half_ph(1,:) , xi )  &
                + 2.0_dp * DOT_PRODUCT( xi , MATMUL( half_ph , xi ) )  
 case default
   call die('eval_f_xdiis: sheme not allowed')
 end select

end function eval_f_xdiis


function eval_gradf_xdiis(xi,dxdt)
 real(dp),intent(in) :: xi(nhist_current)
 real(dp),intent(in) :: dxdt(nhist_current,nhist_current)
 real(dp) :: eval_gradf_xdiis(nhist_current)
!=====
 
 select case(mixing_scheme)
 case('EDIIS')
   eval_gradf_xdiis(:) = MATMUL( diag , dxdt ) + MATMUL( TRANSPOSE(dxdt) , MATMUL( half_ph , xi ) )  &
                                               + MATMUL( xi , MATMUL( half_ph , dxdt ) )
 case('ADIIS')
   eval_gradf_xdiis(:) =  &
                 2.0_dp * MATMUL( TRANSPOSE(dxdt) , half_ph(:,1) ) &
               - 2.0_dp * MATMUL( half_ph(1,:) , dxdt )  &
               + 2.0_dp * MATMUL( TRANSPOSE(dxdt) , MATMUL( half_ph , xi ) )  &
               + 2.0_dp * MATMUL( xi , MATMUL( half_ph , dxdt ) )  
 case default
   call die('eval_gradf_xdiis: sheme not allowed')
 end select

end function eval_gradf_xdiis

end subroutine xdiis_prediction


!=========================================================================
function check_converged(p_matrix_new)
 implicit none

 logical               :: check_converged
 real(dp),intent(in)   :: p_matrix_new(m_ham_scf,n_ham_scf,nspin)
!=====
 real(dp)              :: rms
!=====

 if( parallel_ham ) then
   if( cntxt_ham > 0 ) then
     rms = NORM2( p_matrix_new(:,:,:) - p_matrix_hist(:,:,:,1) )**2
   else
     rms = 0.0_dp
   endif
   call xsum_world(rms)
 else
   rms = NORM2( p_matrix_new(:,:,:) - p_matrix_hist(:,:,:,1) )**2
 endif


 rms = SQRT( rms * nspin )

 write(stdout,'(1x,a,2x,es12.5)') 'Convergence criterium on the density matrix',rms

 if( ( mixing_scheme == 'ADIIS' .OR. mixing_scheme == 'EDIIS' ) .AND. adiis_regime ) then
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

   if( iscf == nscf ) then
     if( rms > 1.0e-2_dp ) then
       call issue_warning('SCF convergence is very poor')
     else if( rms > 1.0e-4_dp ) then
       call issue_warning('SCF convergence is poor')
     endif
   endif

 endif


end function check_converged


!=========================================================================
end module m_scf
