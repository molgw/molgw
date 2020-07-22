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
 use m_linear_algebra,only: invert


 integer,private              :: nhistmax
 integer,private              :: nhist_current

 integer,private              :: nstate,nbf    ! Internal copu of the array dimensions

 real(dp),allocatable,private :: p_matrix_in(:,:,:)
 real(dp),allocatable,private :: ham_hist(:,:,:,:)
 real(dp),allocatable,private :: res_hist(:,:,:,:)
 real(dp),allocatable,private :: p_matrix_hist(:,:,:,:)
 real(dp),allocatable,private :: a_matrix_hist(:,:)
 real(dp),allocatable,private :: p_dot_h_hist(:,:)
 real(dp),allocatable,private :: en_hist(:)
 real(dp),allocatable,private :: alpha_diis(:)

 logical,private              :: adiis_regime

 integer,private              :: iscf

 type energy_contributions
   real(dp) :: nuc_nuc = 0.0_dp
   real(dp) :: kinetic = 0.0_dp
   real(dp) :: nucleus = 0.0_dp
   real(dp) :: hartree = 0.0_dp
   real(dp) :: exx_hyb = 0.0_dp
   real(dp) :: exx     = 0.0_dp
   real(dp) :: xc      = 0.0_dp
   real(dp) :: mp2     = 0.0_dp
   real(dp) :: mp3     = 0.0_dp
   real(dp) :: rpa     = 0.0_dp
   real(dp) :: gw      = 0.0_dp
   real(dp) :: total   = 0.0_dp
   real(dp) :: excit   = 0.0_dp      ! TDDFT excitation energy
   real(dp) :: id      = 0.0_dp      ! H-iD correction
 end type

#if defined(HAVE_SCALAPACK)
 real(dp),external :: PDLANGE
#endif


contains


!=========================================================================
subroutine init_scf(nbf_in,nstate_in)
 implicit none
 integer,intent(in)  :: nbf_in,nstate_in
!=====

 adiis_regime = .FALSE.

 nbf    = nbf_in
 nstate = nstate_in

 iscf          = 0
 nhist_current = 0

 select case(mixing_scheme)
 case('SIMPLE')
   nhistmax         = 2

 case('PULAY','DIIS')
   nhistmax         = npulay_hist
   allocate(a_matrix_hist(nhistmax,nhistmax))

 case('ADIIS','EDIIS')
   adiis_regime = .TRUE.
   nhistmax         = npulay_hist
   allocate(a_matrix_hist(nhistmax,nhistmax))
   allocate(p_dot_h_hist(nhistmax,nhistmax))

 case default
   call die('mixing scheme not implemented')
 end select


 call clean_allocate('Input density matrix',p_matrix_in,nbf,nbf,nspin)
 call clean_allocate('Hamiltonian history',ham_hist,nbf,nbf,nspin,nhistmax)
 call clean_allocate('Density matrix history',p_matrix_hist,nbf,nbf,nspin,nhistmax)
 if( mixing_scheme /= 'SIMPLE' ) then
   call clean_allocate('Residual history',res_hist,nstate,nstate,nspin,nhistmax)
 endif

 allocate(en_hist(nhistmax))


end subroutine init_scf


!=========================================================================
subroutine destroy_scf()
 implicit none
!=====

 if(ALLOCATED(en_hist))          deallocate(en_hist)
 if(ALLOCATED(p_matrix_in))      call clean_deallocate('Input density matrix',p_matrix_in)
 if(ALLOCATED(ham_hist))         call clean_deallocate('Hamiltonian history',ham_hist)
 if(ALLOCATED(res_hist))         call clean_deallocate('Residual history',res_hist)
 if(ALLOCATED(p_matrix_hist))    call clean_deallocate('Density matrix history',p_matrix_hist)
 if(ALLOCATED(a_matrix_hist))    deallocate(a_matrix_hist)
 if(ALLOCATED(p_dot_h_hist))     deallocate(p_dot_h_hist)

end subroutine destroy_scf


!=========================================================================
subroutine hamiltonian_prediction(s_matrix,x_matrix,p_matrix,ham,etot)
 implicit none
 real(dp),intent(in)    :: s_matrix(:,:)
 real(dp),intent(in)    :: x_matrix(:,:)
 real(dp),intent(inout) :: p_matrix(:,:,:)
 real(dp),intent(inout) :: ham(:,:,:)
 real(dp),intent(in)    :: etot
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
 en_hist(2:)             = en_hist(1:nhistmax-1)
 ham_hist(:,:,:,2:)      = ham_hist(:,:,:,1:nhistmax-1)
 p_matrix_hist(:,:,:,2:) = p_matrix_hist(:,:,:,1:nhistmax-1)
 if( ALLOCATED(res_hist) ) then
   res_hist(:,:,:,2:)   = res_hist(:,:,:,1:nhistmax-1)
 endif
 if( ALLOCATED(a_matrix_hist) ) then
   a_matrix_hist(2:,2:) = a_matrix_hist(1:nhistmax-1,1:nhistmax-1)
   a_matrix_hist(1,:)   = 0.0_dp
   a_matrix_hist(:,1)   = 0.0_dp
 endif


 if( mixing_scheme == 'ADIIS' .OR. mixing_scheme == 'EDIIS' ) then
   p_dot_h_hist(2:,2:) = p_dot_h_hist(1:nhistmax-1,1:nhistmax-1)
   p_dot_h_hist(1,:)   = 0.0_dp
   p_dot_h_hist(:,1)   = 0.0_dp
 endif

 !
 ! Set the newest values in history
 ! if already available here
 en_hist(1) = etot
 ham_hist(:,:,:,1)      = ham(:,:,:)
 p_matrix_hist(:,:,:,1) = p_matrix(:,:,:)


 ! Standard Pulay DIIS prediction here !
 if( mixing_scheme /= 'SIMPLE' ) then
   call diis_prediction(s_matrix,x_matrix,p_matrix,ham)
 else
   call simple_prediction(p_matrix,ham)
 endif


 ! If ADIIS or EDIIS prediction, overwrite the hamiltonian with a new one
 if( ( mixing_scheme == 'ADIIS' .OR. mixing_scheme == 'EDIIS' ) .AND. adiis_regime ) then
   call xdiis_prediction(p_matrix,ham)
 endif

 p_matrix_in(:,:,:) = p_matrix(:,:,:)
 deallocate(alpha_diis)

end subroutine hamiltonian_prediction


!=========================================================================
subroutine simple_prediction(p_matrix,ham)
 implicit none
 real(dp),intent(inout) :: p_matrix(:,:,:)
 real(dp),intent(inout) :: ham(:,:,:)
!=====
!=====

 if( nhist_current >= 2 ) then
   alpha_diis(1) = alpha_mixing
   alpha_diis(2) = 1.0_dp - alpha_mixing

   write(stdout,'(/,1x,a,f8.4)') 'Simple mixing with parameter:',alpha_mixing
   ham(:,:,:)      = alpha_mixing * ham_hist(:,:,:,1) + (1.0_dp - alpha_mixing) * ham_hist(:,:,:,2)
   p_matrix(:,:,:) = alpha_mixing * p_matrix_hist(:,:,:,1) + (1.0_dp - alpha_mixing) * p_matrix_hist(:,:,:,2)

 endif
 ham_hist(:,:,:,1) = ham(:,:,:)

end subroutine simple_prediction


!=========================================================================
subroutine diis_prediction(s_matrix,x_matrix,p_matrix,ham)
 implicit none
 real(dp),intent(in)    :: s_matrix(:,:)
 real(dp),intent(in)    :: x_matrix(:,:)
 real(dp),intent(inout) :: p_matrix(:,:,:)
 real(dp),intent(inout) :: ham(:,:,:)
!=====
 integer                :: ispin
 integer                :: ihist
 real(dp),allocatable   :: matrix_tmp1(:,:)
 real(dp),allocatable   :: matrix_tmp2(:,:)
 real(dp),allocatable   :: a_matrix(:,:)
 real(dp),allocatable   :: a_matrix_inv(:,:)
 real(dp),allocatable   :: residual_pred(:,:,:)
 real(dp)               :: residual,work(1)
!=====

 call start_clock(timing_diis)


 write(stdout,'(/,1x,a)') 'Pulay DIIS mixing'


 allocate(a_matrix(nhist_current+1,nhist_current+1))
 allocate(a_matrix_inv(nhist_current+1,nhist_current+1))

 !
 ! Calculate the new residual as proposed in
 ! P. Pulay, J. Comput. Chem. 3, 554 (1982).
 !
 !  R =  X^T * [ H * P * S - S * P * H ] * X

 allocate(matrix_tmp1(nbf,nbf))
 allocate(matrix_tmp2(nbf,nbf))

 do ispin=1,nspin

   ! M2 = P * S
   call DSYMM('L','L',nbf,nbf,1.0d0,p_matrix(1,1,ispin),nbf,s_matrix,nbf,0.0d0,matrix_tmp2,nbf)

   ! M1 = H * M2 = H * P * S
   call DSYMM('L','L',nbf,nbf,1.0d0,ham(1,1,ispin),nbf,matrix_tmp2,nbf,0.0d0,matrix_tmp1,nbf)

   ! M2 = S * P * H = ( H * P * S )**T = M1**T
   matrix_tmp2(:,:) = TRANSPOSE( matrix_tmp1(:,:) )

   ! M1 = H * P * S - S * P * H = M1 - M2
   matrix_tmp1(:,:) = matrix_tmp1(:,:) - matrix_tmp2(:,:)

   !
   ! R = X^T * M1 * X
   ! Remeber that S**-1 = X * X^T
   call matmul_transaba_scalapack(scalapack_block_min,x_matrix,matrix_tmp1,res_hist(:,:,ispin,1))

 enddo

 deallocate(matrix_tmp1)
 deallocate(matrix_tmp2)


 !
 ! Build up a_matrix that contains the scalar product of residuals
 !

 !
 ! The older parts of a_matrix are saved in a_matrix_hist
 ! Just calculate the new ones
 a_matrix_hist(1,1) = SUM( res_hist(:,:,:,1)**2 ) * nspin

 do ihist=2,nhist_current
   a_matrix_hist(ihist,1) = SUM( res_hist(:,:,:,ihist) * res_hist(:,:,:,1) ) * nspin
   a_matrix_hist(1,ihist) = a_matrix_hist(ihist,1)
 enddo


 a_matrix(1:nhist_current,1:nhist_current) = a_matrix_hist(1:nhist_current,1:nhist_current)


 !
 ! DIIS algorithm from Pulay (1980)
 !
 a_matrix(1:nhist_current,nhist_current+1) = -1.0_dp
 a_matrix(nhist_current+1,1:nhist_current) = -1.0_dp
 a_matrix(nhist_current+1,nhist_current+1) =  0.0_dp

 a_matrix_inv(:,:) = a_matrix(:,:)
 call invert(a_matrix_inv)

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
 write(stdout,'(a,30(2x,f12.6))')     ' Alpha DIIS: ',alpha_diis(1:nhist_current)


 !
 ! Calculate the predicted hamiltonian
 !
 !residual_pred(:,:,:) = 0.0_dp
 !ham(:,:,:)           = 0.0_dp
 !p_matrix(:,:,:)      = 0.0_dp
 !do ihist=1,nhist_current
 !  residual_pred(:,:,:) = residual_pred(:,:,:) + alpha_diis(ihist) * res_hist(:,:,:,ihist)
 !  ham(:,:,:)           = ham(:,:,:)           + alpha_diis(ihist) * ham_hist(:,:,:,ihist)
 !  p_matrix(:,:,:)      = p_matrix(:,:,:)      + alpha_diis(ihist) * p_matrix_hist(:,:,:,ihist)
 !enddo

 allocate(residual_pred(SIZE(res_hist,DIM=1),SIZE(res_hist,DIM=2),nspin))

 call DGEMV('N',nstate*nstate*nspin,nhist_current,1.0d0,res_hist     ,nstate*nstate*nspin, &
            alpha_diis,1,0.0d0,residual_pred,1)
 call DGEMV('N',nbf*nbf*nspin      ,nhist_current,1.0d0,ham_hist     ,nbf*nbf*nspin      , &
            alpha_diis,1,0.0d0,ham,1)
 call DGEMV('N',nbf*nbf*nspin      ,nhist_current,1.0d0,p_matrix_hist,nbf*nbf*nspin      , &
            alpha_diis,1,0.0d0,p_matrix,1)

 write(stdout,'(a,2x,es12.5,/)') ' DIIS predicted residual:',NORM2( residual_pred(:,:,:) ) * SQRT(REAL(nspin,dp))

 deallocate(residual_pred)

 deallocate(a_matrix)
 deallocate(a_matrix_inv)


 call stop_clock(timing_diis)

end subroutine diis_prediction


!=========================================================================
subroutine xdiis_prediction(p_matrix,ham)
 use m_lbfgs
 implicit none
 real(dp),intent(out)   :: p_matrix(:,:,:)
 real(dp),intent(out)   :: ham(:,:,:)
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
#if defined(HAVE_SCALAPACK)
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
subroutine density_matrix_preconditioning(hkin,s_matrix,p_matrix_new)
 implicit none

 real(dp),intent(in)      :: hkin(:,:)
 real(dp),intent(in)      :: s_matrix(:,:)
 real(dp),intent(inout)   :: p_matrix_new(:,:,:)
!=====
 real(dp),allocatable     :: hkin_tmp(:,:)
 real(dp),allocatable     :: hkin_inv(:,:)
 real(dp),allocatable     :: delta_p_matrix(:,:)
 integer :: mlocal,nlocal
 integer :: ilocal,jlocal
 integer :: iglobal,jglobal
 integer :: nbf,ispin,ihist
 real(dp),allocatable :: matrix(:,:)
 real(dp) :: trace_ref,trace_current
!=====


 if( kerker_k0 > 1.0e-6_dp ) then

   !
   ! EXPERIMENTAL CODING which is not functional as of today
   call assert_experimental()

   nbf = SIZE(hkin,DIM=1)

   write(stdout,'(1x,a,f8.3)') 'Preconditioning a la Kerker for the density matrix with k0: ',kerker_k0

   allocate(hkin_tmp,SOURCE=hkin)
   allocate(hkin_inv,MOLD=hkin)
   allocate(delta_p_matrix,MOLD=hkin)

   do iglobal=1,nbf
     hkin_tmp(iglobal,iglobal) = hkin_tmp(iglobal,iglobal) + 0.5_dp * kerker_k0**2
   enddo
   hkin_inv(:,:) = hkin_tmp(:,:)
   call invert(hkin_inv)
   hkin_inv(:,:) = -hkin_inv(:,:) * 0.5_dp * kerker_k0**2
   do iglobal=1,nbf
     hkin_inv(iglobal,iglobal) = hkin_inv(iglobal,iglobal) + 1.0_dp
   enddo

   write(stdout,*) '================================'
   do iglobal=1,20
     write(stdout,'(*(2x,f6.2))') hkin(iglobal,1:20)
   enddo
   write(stdout,*) '================================'
   do iglobal=1,20
     write(stdout,'(*(2x,f6.2))') hkin_inv(iglobal,1:20)
   enddo
   write(stdout,*) '================================'

   write(stdout,*) '=============P_MATRIX BEFORE===='
   do iglobal=1,20
     write(stdout,'(*(2x,f6.2))') p_matrix_new(iglobal,1:20,1)
   enddo
   write(stdout,*) '================================'

   allocate(matrix(nbf,nbf))
   matrix(:,:) = MATMUL( p_matrix_new(:,:,1) , s_matrix )
   write(stdout,*) 'TRACE PS',matrix_trace(matrix)



   do ispin=1,nspin
     matrix(:,:) = MATMUL( p_matrix_in(:,:,ispin) , s_matrix )
     trace_ref = matrix_trace(matrix)

     p_matrix_new(:,:,ispin) = p_matrix_in(:,:,ispin) + MATMUL( hkin_inv , p_matrix_new(:,:,ispin) - p_matrix_in(:,:,ispin) )

     matrix(:,:) = MATMUL( p_matrix_new(:,:,ispin) , s_matrix )
     trace_current = matrix_trace(matrix)

     p_matrix_new(:,:,ispin) = p_matrix_new(:,:,ispin) + p_matrix_in(:,:,ispin) * (trace_ref-trace_current) / trace_ref

   enddo

   write(stdout,*) '=============P_MATRIX AFTER ===='
   do iglobal=1,20
     write(stdout,'(*(2x,f6.2))') p_matrix_new(iglobal,1:20,1)
   enddo
   write(stdout,*) '================================'

   matrix(:,:) = MATMUL( p_matrix_new(:,:,1) , s_matrix )
   write(stdout,*) 'TRACE PS after',matrix_trace(matrix)
   deallocate(matrix)

   deallocate(hkin_tmp)
   deallocate(hkin_inv)
   deallocate(delta_p_matrix)

 endif

 if( density_matrix_damping > 1.0e-6 ) then

   write(stdout,'(1x,a,f8.4)') 'Apply a density damping with mixing: ',density_matrix_damping
   do ispin=1,nspin
     p_matrix_new(:,:,:) = p_matrix_in(:,:,:) +  ( p_matrix_new(:,:,:) - p_matrix_in(:,:,:) ) * ( 1.0_dp - density_matrix_damping)
   enddo

 endif


end subroutine density_matrix_preconditioning


!=========================================================================
function check_converged(p_matrix_new)
 implicit none

 logical               :: check_converged
 real(dp),intent(in)   :: p_matrix_new(:,:,:)
!=====
 real(dp)              :: rms
!=====

 rms = NORM2( p_matrix_new(:,:,:) - p_matrix_hist(:,:,:,1) ) * SQRT( REAL(nspin,dp) )


 write(stdout,'(1x,a,es12.5)') 'Convergence criterium on the density matrix: ',rms

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
subroutine print_energy_yaml(name,en)
 implicit none
 character(len=*),intent(in)           :: name
 type(energy_contributions),intent(in) :: en
!=====
!=====

 if( .NOT. ( print_yaml_ .AND. is_iomaster ) ) return

 write(unit_yaml,'(/,a,a)') TRIM(name),':'
 write(unit_yaml,'(4x,a,a)')           'unit:                ','    Ha'
 write(unit_yaml,'(4x,a,1x,es18.8)')   'total:               ',en%total
 write(unit_yaml,'(4x,a,1x,es18.8)')   'nucleus-nucleus:     ',en%nuc_nuc
 write(unit_yaml,'(4x,a,1x,es18.8)')   'kinetic:             ',en%kinetic
 write(unit_yaml,'(4x,a,1x,es18.8)')   'electron-nucleus:    ',en%nucleus
 write(unit_yaml,'(4x,a,1x,es18.8)')   'hartree:             ',en%hartree
 if( ABS(en%exx) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'exchange:            ',en%exx
 if( ABS(en%exx_hyb) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'hybrid exchange:     ',en%exx_hyb
 if( ABS(en%xc) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'exchange-correlation:',en%xc
 if( ABS(en%rpa) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'rpa correlation:     ',en%rpa
 if( ABS(en%gw) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'gw correlation:      ',en%gw
 if( ABS(en%mp2) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'mp2 correlation:     ',en%mp2
 if( ABS(en%mp3) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'mp3 correlation:     ',en%mp3
 if( ABS(en%excit) > 1.0e-10_dp ) &
   write(unit_yaml,'(4x,a,1x,es18.8)') 'excitation:          ',en%excit


end subroutine print_energy_yaml

!=========================================================================


!=========================================================================
end module m_scf
