!!****m* DoNOF/m_optorb
!! NAME
!!  m_optorb
!!
!! FUNCTION
!!  Module prepared to compute orb optimization for a fixed set of occs and 2-RDM
!!
!! PARENTS
!!  m_noft_driver
!!
!! CHILDREN
!!  m_elag
!!  m_diagf
!!  m_e_grad_occ
!!
!! SOURCE
module m_optorb

 use m_nofoutput
 use m_rdmd
 use m_integd
 use m_elag
 use m_hessian
 use m_diagf
 use m_e_grad_occ
 use m_e_grad_occ_cpx
 use m_anti2unit

 implicit none

 private :: lambda_conv
!!***

 public :: opt_orb,s2_calc,num_grad_hess_orb,dm2_JK_3d  
!!***

contains

!!***
!!****f* DoNOF/opt_orb
!! NAME
!!  opt_orb
!!
!! FUNCTION
!!  Call the F-matrix method or QC (with Hessian) for orb optimization 
!!
!! INPUTS
!!  iter=Number of global iteration 
!!  imethod=Method to use 1 -> F diag. 
!!  Vnn=Nuclear-nuclear rep. energy
!!  mo_ints=External subroutine that computes the hCORE and ERI integrals
!!
!! OUTPUT
!!  Energy=Sum of nuclear and electronic energy
!!  NO_COEF=Nat. orb. coefficients
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine opt_orb(iter,imethod,ELAGd,RDMd,INTEGd,HESSIANd,Vnn,Energy,maxdiff,mo_ints,Phases,NO_COEF,NO_COEF_cmplx)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::iter,imethod
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy,maxdiff
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 type(hessian_t),intent(inout)::HESSIANd
 interface 
  subroutine mo_ints(nbf,nstate_occ,nstate_kji,Occ,DM2_JK,NO_COEF,hCORE,ERImol,ERImolJsr,ERImolLsr,&
     &             NO_COEF_cmplx,hCORE_cmplx,ERImol_cmplx,ERImolJsr_cmplx,ERImolLsr_cmplx,all_ERIs,&
     &             Edft_xc,do_xc_dft)
  use m_definitions
  implicit none

  logical,optional,intent(in)     :: all_ERIs,do_xc_dft
  integer,intent(in)              :: nbf,nstate_occ,nstate_kji
  real(dp),intent(in)             :: Occ(nstate_occ)
  real(dp),optional,intent(inout) :: Edft_xc
  real(dp),optional,intent(in)    :: DM2_JK(2,nstate_occ,nstate_occ)
  real(dp),optional,intent(in)    :: NO_COEF(nbf,nbf)
  real(dp),optional,intent(inout) :: hCORE(nbf,nbf)
  real(dp),optional,intent(inout) :: ERImol(nbf,nstate_kji,nstate_kji,nstate_kji)
  real(dp),optional,intent(inout) :: ERImolJsr(nbf,nstate_kji,nstate_kji)
  real(dp),optional,intent(inout) :: ERImolLsr(nbf,nstate_kji,nstate_kji)
  complex(dp),optional,intent(in)    :: NO_COEF_cmplx(nbf,nbf)
  complex(dp),optional,intent(inout) :: hCORE_cmplx(nbf,nbf)
  complex(dp),optional,intent(inout) :: ERImol_cmplx(nbf,nstate_kji,nstate_kji,nstate_kji)
  complex(dp),optional,intent(inout) :: ERImolJsr_cmplx(nbf,nstate_kji,nstate_kji)
  complex(dp),optional,intent(inout) :: ERImolLsr_cmplx(nbf,nstate_kji,nstate_kji)
  end subroutine mo_ints
 end interface
!arrays
 real(dp),intent(inout) :: Phases(RDMd%NBF_occ,RDMd%NBF_occ)
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF_cmplx
!Local variables ------------------------------
!scalars
 logical::convLambda,nogamma,diddiis,allocated_DMNs,all_ERIs
 logical::F_meth_printed,NR_meth_printed
 integer::icall,icall_method,iorbmax1,iorbmax2,imethod_in
 real(dp)::sumdiff,maxdiff_all,Ediff,Energy_old,Edft_xc
!arrays
 real(dp),allocatable,dimension(:,:,:)::DM2_JK
 real(dp),allocatable,dimension(:)::DM2_L_saved
 real(dp),allocatable,dimension(:,:)::U_mat,kappa_mat
 complex(dp),allocatable,dimension(:,:)::U_mat_cmplx,kappa_mat_cmplx
 character(len=200)::msg
!************************************************************************

 icall_method=30
 imethod_in=imethod

 Ediff=zero; Energy=zero; Energy_old=zero; Edft_xc=zero; convLambda=.false.;nogamma=.true.;
 allocated_DMNs=.false.;F_meth_printed=.false.;NR_meth_printed=.false.;all_ERIs=.false.;

 ! If RS-NOFT
 if(INTEGd%irange_sep/=0) then
  allocate(DM2_JK(2,RDMd%NBF_occ,RDMd%NBF_occ))
 endif 

 if((imethod_in==1).and.(iter==0)) then
  ELAGd%sumdiff_old=zero
 endif

 ! Select the method
 if(imethod_in/=1 .and. iter>5 .and. abs(maxdiff)<tol3) then ! TODO Fix the size of the step in the QC method.
  !ELAGd%itolLambda=8; ELAGd%tolE=1e-12;
  write(msg,'(a)') 'Performing QC method for orbital optimization'
  call write_output(msg)
  allocated_DMNs=.true.;all_ERIs=.true.;
  if(INTEGd%complex_ints) then
   allocate(U_mat_cmplx(RDMd%NBF_tot,RDMd%NBF_tot),kappa_mat_cmplx(RDMd%NBF_tot,RDMd%NBF_tot))
  else
   allocate(U_mat(RDMd%NBF_tot,RDMd%NBF_tot),kappa_mat(RDMd%NBF_tot,RDMd%NBF_tot))
  endif
  ! Allocate density matrices for QC method
  allocate(DM2_L_saved(RDMd%NBF_occ*RDMd%NBF_occ))
  DM2_L_saved=RDMd%DM2_L
  RDMd%DM2_K=RDMd%DM2_K+RDMd%DM2_L; RDMd%DM2_L=zero; ! Time-rev. sym DM2_L - added to -> DM2_K
 else
  imethod_in=1
  write(msg,'(a)') 'Building F matrix for orbital optimization'
  call write_output(msg)
 endif
 
 if(INTEGd%complex_ints) then
  if(INTEGd%irange_sep/=0) then
   call dm2_JK_3d(RDMd%NBF_occ,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,RDMd%DM2_iiii,DM2_JK)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,DM2_JK=DM2_JK,NO_COEF_cmplx=NO_COEF_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx,ERImolJsr_cmplx=INTEGd%ERImolJsr_cmplx,  &
  & ERImolLsr_cmplx=INTEGd%ERImolLsr_cmplx,all_ERIs=all_ERIs,Edft_xc=Edft_xc)
  else
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx,all_ERIs=all_ERIs)
  endif
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_old,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
 else
  if(INTEGd%irange_sep/=0) then
   call dm2_JK_3d(RDMd%NBF_occ,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,RDMd%DM2_iiii,DM2_JK)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,DM2_JK=DM2_JK,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, &
   & ERImol=INTEGd%ERImol,ERImolJsr=INTEGd%ERImolJsr,ERImolLsr=INTEGd%ERImolLsr,all_ERIs=all_ERIs,Edft_xc=Edft_xc)
  else
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, & 
   & ERImol=INTEGd%ERImol,all_ERIs=all_ERIs)
  endif
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_old,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
  & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
 endif

 ! Optimization loop
 icall=0
 do
  ! If we used a DIIS step, do not stop after DIIS for small Energy dif.
  diddiis=.false.

  ! Build Lambda matrix
  call ELAGd%build(RDMd,INTEGd,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,RDMd%DM2_Jsr,RDMd%DM2_Lsr)

  ! Check if these NO_COEF with the RDMs are already the solution
  call lambda_conv(ELAGd,RDMd,convLambda,sumdiff,maxdiff,maxdiff_all,iorbmax1,iorbmax2)
  if(convLambda) then
   write(msg,'(a)') 'Lambda_qp - Lambda_pq* converged for the Hemiticty of Lambda'
   call write_output(msg)
   exit
  else
   if(imethod_in==1.and.icall==0) then                                        ! F method: adjust MaxScaling for the rest of orb. icall iterations
    if(iter>2.and.iter>ELAGd%itscale.and.(sumdiff>ELAGd%sumdiff_old)) then ! Parameters chosen from experience to
     ELAGd%itscale=iter+10                                                 ! ensure convergence. Maybe we can set them as input variables?
     ELAGd%MaxScaling=ELAGd%MaxScaling+1
     if(ELAGd%MaxScaling>ELAGd%itolLambda) then
      ELAGd%MaxScaling=2                                                   ! One more empirical/experience parameter used for convergence =(
     endif
    endif
    ELAGd%sumdiff_old=sumdiff
   endif
  endif

  ! Update NO_COEF
  if(imethod_in==1) then ! Build F matrix for iterative diagonalization
   if(INTEGd%complex_ints) then
    call diagF_to_coef(iter,icall,maxdiff,diddiis,ELAGd,RDMd,NO_COEF_cmplx=NO_COEF_cmplx) ! Build new NO_COEF and set icall=icall+1
    ! Build all integrals in the new NO_COEF basis (including arrays for ERI_J and ERI_K)
    if(INTEGd%irange_sep/=0) then
     call dm2_JK_3d(RDMd%NBF_occ,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,RDMd%DM2_iiii,DM2_JK)
     call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,DM2_JK=DM2_JK,NO_COEF_cmplx=NO_COEF_cmplx, &
     & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx,ERImolJsr_cmplx=INTEGd%ERImolJsr_cmplx, &
     & ERImolLsr_cmplx=INTEGd%ERImolLsr_cmplx,all_ERIs=all_ERIs,Edft_xc=Edft_xc)
    else
     call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_cmplx, &
     & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx,all_ERIs=all_ERIs)
    endif
    call INTEGd%eritoeriJKL(RDMd%NBF_occ)
    call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
    & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   else
    call diagF_to_coef(iter,icall,maxdiff,diddiis,ELAGd,RDMd,NO_COEF=NO_COEF) ! Build new NO_COEF and set icall=icall+1
    ! Build all integrals in the new NO_COEF basis (including arrays for ERI_J and ERI_K)
    if(INTEGd%irange_sep/=0) then
     call dm2_JK_3d(RDMd%NBF_occ,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,RDMd%DM2_iiii,DM2_JK)
     call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,DM2_JK=DM2_JK,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, &
     & ERImol=INTEGd%ERImol,ERImolJsr=INTEGd%ERImolJsr,ERImolLsr=INTEGd%ERImolLsr,all_ERIs=all_ERIs,Edft_xc=Edft_xc)
    else
     call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, &
     & ERImol=INTEGd%ERImol,all_ERIs=all_ERIs)
    endif
    call INTEGd%eritoeriJKL(RDMd%NBF_occ)
    call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
    & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
    endif
  else                ! Use QC method to produce new COEFs
   call HESSIANd%build(ELAGd,RDMd,INTEGd,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L)
   if(INTEGd%complex_ints) then
    call HESSIANd%quadratic_conver(icall,RDMd%NBF_tot,kappa_mat_cmplx=kappa_mat_cmplx) ! kappa = - H^-1 g -> norm(kappa)
    call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=kappa_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
    NO_COEF_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
    ! Build all integrals in the new NO_COEF basis (including arrays for ERI_J and ERI_K)
    call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_cmplx, &
    & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx,all_ERIs=all_ERIs)
    call INTEGd%eritoeriJKL(RDMd%NBF_occ)
    call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
    & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   else
    if(INTEGd%irange_sep/=0) then
     write(msg,'(a)') 'Warning! The Hessian of the range-sep is not available'
     call write_output(msg)
     stop
    endif
    call HESSIANd%quadratic_conver(icall,RDMd%NBF_tot,kappa_mat=kappa_mat)             ! kappa = - H^-1 g -> norm(kappa)
    call anti_2_unitary(RDMd%NBF_tot,X_mat=kappa_mat,U_mat=U_mat)                             
    NO_COEF=matmul(NO_COEF,U_mat)
    ! Build all integrals in the new NO_COEF basis (including arrays for ERI_J and ERI_K)
    call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, &
    & ERImol=INTEGd%ERImol,all_ERIs=all_ERIs)
    call INTEGd%eritoeriJKL(RDMd%NBF_occ)
    call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
    & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
   endif
  endif

  ! Check if we did Diag[(Lambda_pq + Lambda_qp*)/2] for F method (first iteration)
  if((imethod_in==1.and.iter==0).and.ELAGd%diagLpL_done) then ! For F method if we did Diag[(Lambda_pq + Lambda_qp*)/2].
   exit                                                    ! -> Do only one icall iteration before the occ. opt.
  endif

  ! For this icall using the new NO_COEF (and fixed RDMs). Is the Energy still changing?
  Ediff=Energy-Energy_old
  if((icall>1).and.(dabs(Ediff)<ELAGd%tolE).and.(.not.diddiis)) then ! The energy is not changing anymore (not stopping for DIIS itertation)
   write(msg,'(a)') 'Lambda_qp - Lambda_pq* converged for small energy differences'
   call write_output(msg)
   exit
  endif
  Energy_old=Energy

  ! We allow at most 30 generations of new NO_COEF updates (and integrals) in Piris Ugalde Algoritms
  if(icall==icall_method) exit
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 enddo

 ! Restore density matrices used in QC method  
 if(allocated_DMNs) then
  RDMd%DM2_L=DM2_L_saved
  RDMd%DM2_K=RDMd%DM2_K-RDMd%DM2_L
  deallocate(DM2_L_saved)
  if(INTEGd%complex_ints) then
   deallocate(U_mat_cmplx,kappa_mat_cmplx)
  else
   deallocate(U_mat,kappa_mat)
  endif
 endif

 ! If RS-NOFT
 if(INTEGd%irange_sep/=0) then
  deallocate(DM2_JK)
 endif 

 ! Calc. the final Energy using fixed RDMs and the new NO_COEF (before going back to occ. optimization)
 if(INTEGd%complex_ints) then
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
 else
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
  & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
 endif
 write(msg,'(a,f15.6,a,i6,a)') 'Orb. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
 call write_output(msg)
 if(abs(Edft_xc)>tol8) then
  write(msg,'(a,f15.6)') 'Exc sr-DFT energy    = ',Edft_xc
  call write_output(msg)
 endif
 if(iter>0) then
  write(msg,'(a,f15.6,a,i5,a,i5,a)') 'Max. [Lambda_qp - Lambda_pq*]= ',maxdiff_all,' pair (',iorbmax1,',',iorbmax2,')'
  call write_output(msg)
  write(msg,'(a,f19.10)') 'Energy difference orb. opt.=',Ediff
  call write_output(msg)
 endif
 
 !if(icall==30) write(*,'(a)') 'Warning! Max. number of iterations (30) reached in orb. optimization'

end subroutine opt_orb
!!***

!!***
!!****f* DoNOF/lambda_conv
!! NAME
!!  lambda_conv
!!
!! FUNCTION
!!  Check if the Lambda matrix already fulfils the condition Lambda_qp-Lambda_pq^* <= tol_dif_lambda.
!!  (i.e. the gradient of the exp^-kappa parametrization)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine lambda_conv(ELAGd,RDMd,converg_lamb,sumdiff,maxdiff,maxdiff_all,iorbmax1,iorbmax2)
!Arguments ------------------------------------
!scalars
 logical,intent(inout)::converg_lamb
 integer,intent(inout)::iorbmax1,iorbmax2
 real(dp),intent(inout)::sumdiff,maxdiff,maxdiff_all
 type(elag_t),intent(in)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 real(dp)::diff,tol_dif_Lambda
!arrays
!************************************************************************

 tol_dif_Lambda=ten**(-ELAGd%itolLambda)
 converg_lamb=.true.; sumdiff=zero; maxdiff=zero; maxdiff_all=zero;
 
 do iorb=1,RDMd%NBF_tot
  do iorb1=1,RDMd%NBF_tot
   diff=two*dabs(ELAGd%Lambdas(iorb1,iorb)-ELAGd%Lambdas(iorb,iorb1))
   if(ELAGd%cpx_lambdas .and. iorb/=iorb1) then
    diff=diff*diff+(two*(ELAGd%Lambdas_im(iorb1,iorb)+ELAGd%Lambdas_im(iorb,iorb1)))**two
    diff=dsqrt(diff)
   endif
   if(ELAGd%cpx_lambdas .and. iorb1==iorb) then
    diff=two*dabs(ELAGd%Lambdas_im(iorb1,iorb)+ELAGd%Lambdas_im(iorb,iorb1))
   endif
   sumdiff=sumdiff+diff
   if((diff>=tol_dif_Lambda) .and. converg_lamb) then
    converg_lamb=.false.
   endif
   if(ELAGd%cpx_lambdas) then
    if(diff>maxdiff .and.iorb/=iorb1) then
     maxdiff=diff
    endif
   else
    if(diff>maxdiff) then
     maxdiff=diff
    endif
   endif
   if(diff>maxdiff_all) then
    maxdiff_all=diff
    iorbmax1=iorb
    iorbmax2=iorb1
   endif
  enddo
 enddo

end subroutine lambda_conv
!!***

!!***
!!****f* DoNOF/s2_calc
!! NAME
!!  s2_calc
!!
!! FUNCTION
!!  Compute <S^2>
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine s2_calc(RDMd,INTEGd,NO_COEF,NO_COEF_cmplx) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF_cmplx
!Local variables ------------------------------
!scalars
 logical::intrinsic_cmplx=.false.
 integer::NBF2,NsdORBs,NsdVIRT
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5,iorb6,iorb7,iorbmin,iorbmax
 real(dp)::Delem,Nelect
 complex(dp)::S2_val
!arrays
 integer,allocatable,dimension(:)::coup
 real(dp),allocatable,dimension(:)::occ,occ_dyn
 real(dp),allocatable,dimension(:,:)::S_pq_NO
 complex(dp),allocatable,dimension(:,:)::S_pq_NO_cmplx
 character(len=200)::msg
!************************************************************************
!************************************************************************

 NBF2=2*RDMd%NBF_occ
 NsdORBs=RDMd%Nfrozen+RDMd%Npairs
 NsdVIRT=RDMd%NBF_occ-NsdORBs
 allocate(coup(NBF2),occ(NBF2),occ_dyn(NBF2))
 occ=zero;occ_dyn=zero;Nelect=zero;

 ! Get sw 1-RDM
 do iorb=1,RDMd%NBF_occ
  Delem=RDMd%occ(iorb)
  Nelect=Nelect+two*Delem
  if(dabs(Delem)>tol8) then
   occ(2*iorb-1)=Delem
   occ(2*iorb)  =occ(2*iorb-1)
  endif
  Delem=RDMd%occ_dyn(iorb)
  if(dabs(Delem)>tol8) then
   occ_dyn(2*iorb-1)=Delem
   occ_dyn(2*iorb)  =occ_dyn(2*iorb-1)
  endif
 enddo
 do iorb=1,2*RDMd%Nfrozen
  occ_dyn(iorb)=zero
 enddo

 ! Prepare coupling map
 coup=-1
 do iorb=1,NsdORBs
  iorbmin=NsdORBs+RDMd%Ncoupled*(NsdORBs-iorb)+1
  iorbmax=NsdORBs+RDMd%Ncoupled*(NsdORBs-iorb)+RDMd%Ncoupled
  do iorb1=1,NsdVIRT
   iorb2=iorb1+NsdORBs
   if((iorbmin<=iorb2.and.iorb2<=iorbmax).and.(iorbmax<=RDMd%NBF_occ)) then
    coup(2*iorb-1) =iorb
    coup(2*iorb)   =iorb
    coup(2*iorb2-1)=iorb
    coup(2*iorb2)  =iorb
   endif
  enddo
 enddo
 
 S2_val=-Nelect*(Nelect-four)/four
 if(INTEGd%complex_ints) then
  allocate(S_pq_NO_cmplx(RDMd%NBF_tot,RDMd%NBF_tot))
  S_pq_NO_cmplx=matmul(transpose(NO_COEF_cmplx),matmul(INTEGd%Overlap,NO_COEF_cmplx)) ! S_pq ^alpha beta not
                                                                                      ! using conjugate
  do iorb=1,RDMd%NBF_tot
   do iorb1=1,iorb-1
    if(abs(S_pq_NO_cmplx(iorb,iorb1))>tol8) intrinsic_cmplx=.true.
   enddo
   if(abs(abs(S_pq_NO_cmplx(iorb,iorb))-ONE)>tol8) intrinsic_cmplx=.true.
  enddo
  if(intrinsic_cmplx) then
   write(msg,'(a,f19.10)') " Orbitals are intrinsically complex ( S_pq ^ss' /= delta_pq Exp[i theta_p] )"
   call write_output(msg)
  else
   write(msg,'(a,f19.10)') " Orbitals are intrinsically real ( S_pq ^ss' = delta_pq Exp[i theta_p] )"
   call write_output(msg)
  endif
  ! Calc. <S^2> with the sw 2-RDM
  do iorb=1,NBF2
   do iorb1=1,NBF2
    do iorb2=1,NBF2
     do iorb3=1,NBF2
      ! Calculate the ^2D_ij,kl element
      if((mod(iorb,2)==mod(iorb2,2)).and.(mod(iorb1,2)==mod(iorb3,2))) then
       if(mod(iorb,2)==0) then
        iorb4=iorb/2
        iorb6=iorb2/2
       else
        iorb4=(iorb+1)/2
        iorb6=(iorb2+1)/2
       endif
       if(mod(iorb1,2)==0) then
        iorb5=iorb1/2
        iorb7=iorb3/2
       else
        iorb5=(iorb1+1)/2
        iorb7=(iorb3+1)/2
       endif
       if((mod(iorb,2)==mod(iorb1,2)).and.(iorb==iorb2.and.iorb1==iorb3)) then ! alpha alpha alpha alpha or beta beta beta beta
        call compu_swdm2(RDMd,occ,occ_dyn,coup,iorb,iorb1,iorb2,iorb3,NBF2,NsdORBs,Delem)
        S2_val=S2_val+Delem
       endif
       if((mod(iorb,2)/=mod(iorb1,2))) then ! alpha beta alpha beta or beta alpha beta alpha
        call compu_swdm2(RDMd,occ,occ_dyn,coup,iorb,iorb1,iorb2,iorb3,NBF2,NsdORBs,Delem)
        if(mod(iorb,2)/=0) then ! alpha beta alpha beta
         S2_val=S2_val-Delem*S_pq_NO_cmplx(iorb4,iorb7)*S_pq_NO_cmplx(iorb5,iorb6)
        else ! beta alpha beta alpha
         S2_val=S2_val-Delem*conjg(S_pq_NO_cmplx(iorb4,iorb7)*S_pq_NO_cmplx(iorb5,iorb6))
        endif
       endif
      endif
     enddo
    enddo
   enddo
  enddo
  deallocate(S_pq_NO_cmplx)
 else
  allocate(S_pq_NO(RDMd%NBF_tot,RDMd%NBF_tot))
  S_pq_NO=matmul((transpose(NO_COEF)),matmul(INTEGd%Overlap,NO_COEF))
  ! Calc. <S^2> with the sw 2-RDM
  do iorb=1,NBF2
   do iorb1=1,NBF2
    do iorb2=1,NBF2
     do iorb3=1,NBF2
      ! Calculate the ^2D_ij,kl element
      if((mod(iorb,2)==mod(iorb2,2)).and.(mod(iorb1,2)==mod(iorb3,2))) then
       if(mod(iorb,2)==0) then
        iorb4=iorb/2
        iorb6=iorb2/2
       else
        iorb4=(iorb+1)/2
        iorb6=(iorb2+1)/2
       endif
       if(mod(iorb1,2)==0) then
        iorb5=iorb1/2
        iorb7=iorb3/2
       else
        iorb5=(iorb1+1)/2
        iorb7=(iorb3+1)/2
       endif
       if((mod(iorb,2)==mod(iorb1,2)).and.(iorb==iorb2.and.iorb1==iorb3)) then ! alpha alpha alpha alpha or beta beta beta beta
        call compu_swdm2(RDMd,occ,occ_dyn,coup,iorb,iorb1,iorb2,iorb3,NBF2,NsdORBs,Delem)
        S2_val=S2_val+Delem
       endif
       if((mod(iorb,2)/=mod(iorb1,2))) then ! alpha beta alpha beta or beta alpha beta alpha
        call compu_swdm2(RDMd,occ,occ_dyn,coup,iorb,iorb1,iorb2,iorb3,NBF2,NsdORBs,Delem)
        if(mod(iorb,2)/=0) then ! alpha beta alpha beta
         S2_val=S2_val-Delem*S_pq_NO(iorb4,iorb7)*S_pq_NO(iorb5,iorb6)
        else ! beta alpha beta alpha
         S2_val=S2_val-Delem*S_pq_NO(iorb4,iorb7)*S_pq_NO(iorb5,iorb6)
        endif
       endif
      endif
     enddo
    enddo
   enddo
  enddo
  deallocate(S_pq_NO)
 endif
 
 write(msg,'(a,f10.5,a)') ' Computed <S^2> = ',abs(real(S2_val))
 call write_output(msg)
 if(abs(aimag(S2_val))>tol8) then
  write(msg,'(a,f10.3,a)') ' Warning! Imaginary contribution found evaluating <S^2> = ',aimag(S2_val)
  call write_output(msg)
 endif
 deallocate(coup,occ,occ_dyn) 

end subroutine s2_calc
!!***

!!***
!!****f* DoNOF/num_grad_hess_orb
!! NAME
!!  num_grad_hess_orb
!!
!! FUNCTION
!!  Calculate the numerical orbital rotation gradient and Hessian for a k_pq k_rs rot.
!!  Note: This subroutine calls mo_ints NOT building all ERIs. Only builds the ones 
!!        that are needed! all_ERIs=.false.
!!
!! INPUTS
!!  iorbx=orbitals to use
!!  NO_COEF and NO_COEF_cmplx are the current Nat. orb. coefs. 
!!  mo_ints=External subroutine that computes the hCORE and ERI integrals
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine num_grad_hess_orb(iorbp,iorbq,iorbr,iorbs,RDMd,INTEGd,Vnn,mo_ints,Phases,NO_COEF,NO_COEF_cmplx) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::iorbp,iorbq,iorbr,iorbs
 real(dp),intent(in)::Vnn
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 interface 
  subroutine mo_ints(NBF_tot,NBF_occ,NBF_jkl,Occ,NO_COEF,hCORE,ERImol,ERImolJsr,ERImolLsr,&
  & NO_COEF_cmplx,hCORE_cmplx,ERImol_cmplx,all_ERIs)
  use m_definitions
  implicit none
  logical,optional,intent(in)::all_ERIs
  integer,intent(in)::NBF_tot,NBF_occ,NBF_jkl
  real(dp),intent(in)::Occ(NBF_occ)
  real(dp),optional,intent(in)::NO_COEF(NBF_tot,NBF_tot)
  real(dp),optional,intent(inout)::hCORE(NBF_tot,NBF_tot)
  real(dp),optional,intent(inout)::ERImol(NBF_tot,NBF_jkl,NBF_jkl,NBF_jkl)
  real(dp),optional,intent(inout)::ERImolJsr(NBF_tot,NBF_jkl,NBF_jkl)
  real(dp),optional,intent(inout)::ERImolLsr(NBF_tot,NBF_jkl,NBF_jkl)
  complex(dp),optional,intent(in)::NO_COEF_cmplx(NBF_tot,NBF_tot)
  complex(dp),optional,intent(inout)::hCORE_cmplx(NBF_tot,NBF_tot)
  complex(dp),optional,intent(inout)::ERImol_cmplx(NBF_tot,NBF_jkl,NBF_jkl,NBF_jkl)
  end subroutine mo_ints
 end interface
!arrays
 real(dp),intent(inout)::Phases(RDMd%NBF_occ,RDMd%NBF_occ)
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF_cmplx
!Local variables ------------------------------
!scalars
 logical::nogamma=.true.
 real(dp)::Energy_in,Energy_00,Energy_h0,Energy_mh0,Energy_0h,Energy_0mh,Energy_hh,Energy_mhmh
 real(dp)::Gradient_pq,Gradient_rs,Hessian_pqrs,step
 real(dp)::Energy_00_00,Energy_h0_00,Energy_0h_00,Energy_00_h0,Energy_00_0h
 real(dp)::Energy_mh0_00,Energy_0mh_00,Energy_00_mh0,Energy_00_0mh
 real(dp)::Energy_h0_h0,Energy_mh0_mh0,Energy_0h_0h,Energy_0mh_0mh
 real(dp)::Energy_h0_0h,Energy_mh0_0mh,Energy_0h_h0,Energy_0mh_mh0
 complex(dp)::Gradient_pq_cmplx,Gradient_rs_cmplx,Hessian_pqrs_cmplx
 complex(dp)::Hessian_pqrs_cmplx_RR,Hessian_pqrs_cmplx_II,Hessian_pqrs_cmplx_RI,Hessian_pqrs_cmplx_IR
!arrays
 character(len=200)::msg
 real(dp),allocatable,dimension(:,:)::U_mat,X_mat,NO_COEF_in
 complex(dp),allocatable,dimension(:,:)::U_mat_cmplx,X_mat_cmplx,NO_COEF_in_cmplx

!************************************************************************
 
 step=tol3;

 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') 'Computing orb. rot. numerical Gradient and Hessian'
 call write_output(msg)

 if(INTEGd%complex_ints) then ! Complex

  allocate(U_mat_cmplx(RDMd%NBF_tot,RDMd%NBF_tot),X_mat_cmplx(RDMd%NBF_tot,RDMd%NBF_tot))
  allocate(NO_COEF_in_cmplx(RDMd%NBF_tot,RDMd%NBF_tot))
  U_mat_cmplx=complex_zero;
  Gradient_pq_cmplx=complex_zero;Gradient_rs_cmplx=complex_zero;Hessian_pqrs_cmplx=complex_zero;
  Energy_00_00=zero;Energy_h0_00=zero;Energy_0h_00=zero;Energy_00_h0=zero;Energy_00_0h=zero;
  Energy_mh0_00=zero;Energy_0mh_00=zero;Energy_00_mh0=zero;Energy_00_0mh=zero; 
  Energy_h0_h0=zero;Energy_mh0_mh0=zero;Energy_0h_0h=zero;Energy_0mh_0mh=zero; 
  Energy_h0_0h=zero;Energy_mh0_0mh=zero;Energy_0h_h0=zero;Energy_0mh_mh0=zero;
  ! kappa_pq_R=0,kappa_pq_I=0,kappa_rs_R=0,kappa_rs_I=0
  X_mat_cmplx=complex_zero;
  call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
  NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
  Energy_00_00=Energy_in+Vnn
  ! kappa_pq_R=step,kappa_pq_I=0,kappa_rs_R=0,kappa_rs_I=0
  if(iorbp/=iorbq) then
   X_mat_cmplx=complex_zero;
   X_mat_cmplx(iorbp,iorbq)=step;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_h0_00=Energy_in+Vnn
  endif
  ! kappa_pq_R=0,kappa_pq_I=step,kappa_rs_R=0,kappa_rs_I=0
  X_mat_cmplx=complex_zero;
  X_mat_cmplx(iorbp,iorbq)=step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
  call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
  NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
  Energy_0h_00=Energy_in+Vnn
  ! kappa_pq_R=0,kappa_pq_I=0,kappa_rs_R=step,kappa_rs_I=0
  if(iorbr/=iorbs) then
   X_mat_cmplx=complex_zero;
   X_mat_cmplx(iorbr,iorbs)=step;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_00_h0=Energy_in+Vnn
  endif
  ! kappa_pq_R=0,kappa_pq_I=0,kappa_rs_R=0,kappa_rs_I=step
  X_mat_cmplx=complex_zero;
  X_mat_cmplx(iorbr,iorbs)=step*im;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
  call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
  NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
  Energy_00_0h=Energy_in+Vnn
  ! kappa_pq_R=-step,kappa_pq_I=0,kappa_rs_R=0,kappa_rs_I=0
  if(iorbp/=iorbq) then
   X_mat_cmplx=complex_zero;
   X_mat_cmplx(iorbp,iorbq)=-step;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_mh0_00=Energy_in+Vnn
  endif
  ! kappa_pq_R=0,kappa_pq_I=-step,kappa_rs_R=0,kappa_rs_I=0
  X_mat_cmplx=complex_zero;
  X_mat_cmplx(iorbp,iorbq)=-step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
  call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
  NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
  Energy_0mh_00=Energy_in+Vnn
  ! kappa_pq_R=0,kappa_pq_I=0,kappa_rs_R=-step,kappa_rs_I=0
  if(iorbr/=iorbs) then
   X_mat_cmplx=complex_zero;
   X_mat_cmplx(iorbr,iorbs)=-step;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_00_mh0=Energy_in+Vnn
  endif
  ! kappa_pq_R=0,kappa_pq_I=0,kappa_rs_R=0,kappa_rs_I=-step
  X_mat_cmplx=complex_zero;
  X_mat_cmplx(iorbr,iorbs)=-step*im;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
  call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
  NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
  Energy_00_0mh=Energy_in+Vnn
  ! kappa_pq_R=step,kappa_pq_I=0,kappa_rs_R=step,kappa_rs_I=0
  if(iorbp/=iorbq .and. iorbr/=iorbs) then
   X_mat_cmplx=complex_zero;
   X_mat_cmplx(iorbp,iorbq)=step;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   X_mat_cmplx(iorbr,iorbs)=step;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_h0_h0=Energy_in+Vnn
  endif
  ! kappa_pq_R=-step,kappa_pq_I=0,kappa_rs_R=-step,kappa_rs_I=0
  if(iorbp/=iorbq .and. iorbr/=iorbs) then
   X_mat_cmplx=complex_zero;
   X_mat_cmplx(iorbp,iorbq)=-step;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   X_mat_cmplx(iorbr,iorbs)=-step;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_mh0_mh0=Energy_in+Vnn
  endif
  ! kappa_pq_R=0,kappa_pq_I=step,kappa_rs_R=0,kappa_rs_I=step
  X_mat_cmplx=complex_zero;
  X_mat_cmplx(iorbp,iorbq)=step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
  X_mat_cmplx(iorbr,iorbs)=step*im;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
  call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
  NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
  Energy_0h_0h=Energy_in+Vnn
  ! kappa_pq_R=0,kappa_pq_I=-step,kappa_rs_R=0,kappa_rs_I=-step
  X_mat_cmplx=complex_zero;
  X_mat_cmplx(iorbp,iorbq)=-step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
  X_mat_cmplx(iorbr,iorbs)=-step*im;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
  call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
  NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
  Energy_0mh_0mh=Energy_in+Vnn
  ! kappa_pq_R=step,kappa_pq_I=0,kappa_rs_R=0,kappa_rs_I=step
  if(iorbp/=iorbq) then
   X_mat_cmplx=complex_zero;
   if(iorbp==iorbr .and. iorbq==iorbs) then
    X_mat_cmplx(iorbp,iorbq)=step+step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   else
    X_mat_cmplx(iorbp,iorbq)=step;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
    X_mat_cmplx(iorbr,iorbs)=step*im;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   endif
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_h0_0h=Energy_in+Vnn
  endif
  ! kappa_pq_R=-step,kappa_pq_I=0,kappa_rs_R=0,kappa_rs_I=-step
  if(iorbp/=iorbq) then
   X_mat_cmplx=complex_zero;
   if(iorbp==iorbr .and. iorbq==iorbs) then
    X_mat_cmplx(iorbp,iorbq)=-step-step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   else
    X_mat_cmplx(iorbp,iorbq)=-step;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
    X_mat_cmplx(iorbr,iorbs)=-step*im;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   endif
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_mh0_0mh=Energy_in+Vnn
  endif
  ! kappa_pq_R=0,kappa_pq_I=step,kappa_rs_R=step,kappa_rs_I=0
  if(iorbr/=iorbs) then
   X_mat_cmplx=complex_zero;
   if(iorbp==iorbr .and. iorbq==iorbs) then
    X_mat_cmplx(iorbp,iorbq)=step+step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   else
    X_mat_cmplx(iorbp,iorbq)=step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
    X_mat_cmplx(iorbr,iorbs)=step;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   endif
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_0h_h0=Energy_in+Vnn
  endif
  ! kappa_pq_R=0,kappa_pq_I=-step,kappa_rs_R=-step,kappa_rs_I=0
  if(iorbr/=iorbs) then
   X_mat_cmplx=complex_zero;
   if(iorbp==iorbr .and. iorbq==iorbs) then
    X_mat_cmplx(iorbp,iorbq)=-step-step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
   else
    X_mat_cmplx(iorbp,iorbq)=-step*im;X_mat_cmplx(iorbq,iorbp)=-conjg(X_mat_cmplx(iorbp,iorbq));
    X_mat_cmplx(iorbr,iorbs)=-step;X_mat_cmplx(iorbs,iorbr)=-conjg(X_mat_cmplx(iorbr,iorbs));
   endif
   call anti_2_unitary(RDMd%NBF_tot,X_mat_cmplx=X_mat_cmplx,U_mat_cmplx=U_mat_cmplx)
   NO_COEF_in_cmplx=matmul(NO_COEF_cmplx,U_mat_cmplx)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_in_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,nogamma=nogamma)
   Energy_0mh_mh0=Energy_in+Vnn
  endif
  ! Gradient and Hessian
  Gradient_pq_cmplx=(Energy_h0_00-Energy_mh0_00)/(two*step)+(Energy_0h_00-Energy_0mh_00)*im/(two*step)
  Gradient_rs_cmplx=(Energy_00_h0-Energy_00_mh0)/(two*step)+(Energy_00_0h-Energy_00_0mh)*im/(two*step)
  if(iorbp==iorbq .or. iorbr==iorbs) then
   Hessian_pqrs_cmplx_RR=complex_zero
  else
   if((iorbp==iorbr .and. iorbq==iorbs).or.(iorbp==iorbs .and. iorbq==iorbr)) then
    Hessian_pqrs_cmplx_RR=(Energy_h0_00-two*Energy_00_00+Energy_mh0_00)/(step*step)
    if(iorbp==iorbs .and. iorbq==iorbr) Hessian_pqrs_cmplx_RR=-Hessian_pqrs_cmplx_RR
   else
    Hessian_pqrs_cmplx_RR=(Energy_h0_h0-Energy_h0_00-Energy_00_h0+two*Energy_00_00 &
    &                     -Energy_mh0_00-Energy_00_mh0+Energy_mh0_mh0)/(two*step*step)
   endif
  endif
  if((iorbp==iorbr .and. iorbq==iorbs).or.(iorbp==iorbs .and. iorbq==iorbr)) then
   Hessian_pqrs_cmplx_II=(Energy_0h_00-two*Energy_00_00+Energy_0mh_00)/(step*step)
   if(iorbp==iorbs .and. iorbq==iorbr) Hessian_pqrs_cmplx_II=-Hessian_pqrs_cmplx_II
  else
   Hessian_pqrs_cmplx_II=(Energy_0h_0h-Energy_0h_00-Energy_00_0h+two*Energy_00_00 &
   &                     -Energy_0mh_00-Energy_00_0mh+Energy_0mh_0mh)/(two*step*step)
  endif
  if(iorbp/=iorbq) then
   Hessian_pqrs_cmplx_RI=-(Energy_h0_0h-Energy_h0_00-Energy_00_0h+two*Energy_00_00 &
  &                     -Energy_mh0_00-Energy_00_0mh+Energy_mh0_0mh)/(two*step*step)
  endif
  if(iorbr/=iorbs) then
   Hessian_pqrs_cmplx_IR=(Energy_0h_h0-Energy_0h_00-Energy_00_h0+two*Energy_00_00 &
  &                     -Energy_0mh_00-Energy_00_mh0+Energy_0mh_mh0)/(two*step*step)
  endif
  Hessian_pqrs_cmplx= Hessian_pqrs_cmplx_RR+Hessian_pqrs_cmplx_II &
  &                 +(Hessian_pqrs_cmplx_RI+Hessian_pqrs_cmplx_IR)*im
  if(iorbp==iorbq .or. iorbr==iorbs) Hessian_pqrs_cmplx=two*Hessian_pqrs_cmplx
  write(msg,'(a,f15.6)') 'Energy = ',Energy_00_00
  call write_output(msg)
  write(msg,'(a,i5,a,i5,a,f15.6,a,f15.6,a)') 'Grad(',iorbp,',',iorbq,') = (', &
  & real(Gradient_pq_cmplx),',',aimag(Gradient_pq_cmplx),')'
  call write_output(msg)
  write(msg,'(a,i5,a,i5,a,f15.6,a,f15.6,a)') 'Grad(',iorbr,',',iorbs,') = (', &
  & real(Gradient_rs_cmplx),',',aimag(Gradient_rs_cmplx),')'
  call write_output(msg)
  write(msg,'(a,i5,a,i5,a,i5,a,i5,a,f15.6,a,f15.6,a)') 'Hessian(',iorbp,',',iorbq,',',iorbr,',',iorbs,') = (', &
  & real(Hessian_pqrs_cmplx),',',aimag(Hessian_pqrs_cmplx),')'
  call write_output(msg)
  write(msg,'(a)') ' '
  call write_output(msg)
  ! Recover initial integrals and deallocate arrays
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  deallocate(U_mat_cmplx,X_mat_cmplx,NO_COEF_in_cmplx)

 else ! Real

  allocate(U_mat(RDMd%NBF_tot,RDMd%NBF_tot),X_mat(RDMd%NBF_tot,RDMd%NBF_tot))
  allocate(NO_COEF_in(RDMd%NBF_tot,RDMd%NBF_tot))
  U_mat=zero;Gradient_pq=zero;Gradient_rs=zero;Hessian_pqrs=zero;
  Energy_00=zero;Energy_h0=zero;Energy_mh0=zero;Energy_0h=zero;Energy_0mh=zero;Energy_hh=zero;Energy_mhmh=zero;
  ! kappa_pq=0,kappa_rs=0
  X_mat=zero;
  call anti_2_unitary(RDMd%NBF_tot,X_mat=X_mat,U_mat=U_mat)
  NO_COEF_in=matmul(NO_COEF,U_mat)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF_in,hCORE=INTEGd%hCORE, &
  & ERImol=INTEGd%ERImol)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
  & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
  Energy_00=Energy_in+Vnn
  if(iorbp/=iorbq) then ! avoiding k_pp
   ! kappa_pq=step,kappa_rs=0
   X_mat=zero;
   X_mat(iorbp,iorbq)=step;X_mat(iorbq,iorbp)=-X_mat(iorbp,iorbq);
   call anti_2_unitary(RDMd%NBF_tot,X_mat=X_mat,U_mat=U_mat)
   NO_COEF_in=matmul(NO_COEF,U_mat)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF_in,hCORE=INTEGd%hCORE, &
   & ERImol=INTEGd%ERImol)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
   Energy_h0=Energy_in+Vnn
   ! kappa_pq=-step,kappa_rs=0
   X_mat=zero;
   X_mat(iorbp,iorbq)=-step;X_mat(iorbq,iorbp)=-X_mat(iorbp,iorbq);
   call anti_2_unitary(RDMd%NBF_tot,X_mat=X_mat,U_mat=U_mat)
   NO_COEF_in=matmul(NO_COEF,U_mat)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF_in,hCORE=INTEGd%hCORE, &
   & ERImol=INTEGd%ERImol)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
   Energy_mh0=Energy_in+Vnn
  endif
  if(iorbr/=iorbs) then ! avoiding k_rr
   ! kappa_pq=0,kappa_rs=step
   X_mat=zero;
   X_mat(iorbr,iorbs)=step;X_mat(iorbs,iorbr)=-X_mat(iorbr,iorbs);
   call anti_2_unitary(RDMd%NBF_tot,X_mat=X_mat,U_mat=U_mat)
   NO_COEF_in=matmul(NO_COEF,U_mat)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF_in,hCORE=INTEGd%hCORE, &
   & ERImol=INTEGd%ERImol)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
   Energy_0h=Energy_in+Vnn
   ! kappa_pq=0,kappa_rs=-step
   X_mat=zero;
   X_mat(iorbr,iorbs)=-step;X_mat(iorbs,iorbr)=-X_mat(iorbr,iorbs);
   call anti_2_unitary(RDMd%NBF_tot,X_mat=X_mat,U_mat=U_mat)
   NO_COEF_in=matmul(NO_COEF,U_mat)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF_in,hCORE=INTEGd%hCORE, &
   & ERImol=INTEGd%ERImol)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
   Energy_0mh=Energy_in+Vnn
  endif
  if(iorbp/=iorbq .and. iorbr/=iorbs) then ! avoiding k_pp and k_rr
   ! kappa_pq=step,kappa_rs=step
   X_mat=zero;
   X_mat(iorbp,iorbq)=step;X_mat(iorbq,iorbp)=-X_mat(iorbp,iorbq);
   X_mat(iorbr,iorbs)=step;X_mat(iorbs,iorbr)=-X_mat(iorbr,iorbs);
   call anti_2_unitary(RDMd%NBF_tot,X_mat=X_mat,U_mat=U_mat)
   NO_COEF_in=matmul(NO_COEF,U_mat)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF_in,hCORE=INTEGd%hCORE, &
   & ERImol=INTEGd%ERImol)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
   Energy_hh=Energy_in+Vnn
   ! kappa_pq=-step,kappa_rs=-step
   X_mat=zero;
   X_mat(iorbp,iorbq)=-step;X_mat(iorbq,iorbp)=-X_mat(iorbp,iorbq);
   X_mat(iorbr,iorbs)=-step;X_mat(iorbs,iorbr)=-X_mat(iorbr,iorbs);
   call anti_2_unitary(RDMd%NBF_tot,X_mat=X_mat,U_mat=U_mat)
   NO_COEF_in=matmul(NO_COEF,U_mat)
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF_in,hCORE=INTEGd%hCORE, &
   & ERImol=INTEGd%ERImol)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_in,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
   Energy_mhmh=Energy_in+Vnn
  endif
  ! Gradient and Hessian
  Gradient_pq=(Energy_h0-Energy_mh0)/(two*step)
  Gradient_rs=(Energy_0h-Energy_0mh)/(two*step)
  if(iorbp/=iorbq .and. iorbr/=iorbs) then ! avoiding k_pp and k_rr
   if((iorbp==iorbr .and. iorbq==iorbs).or.(iorbp==iorbs .and. iorbq==iorbr)) then
    Hessian_pqrs=(Energy_h0-two*Energy_00+Energy_mh0)/(step*step)
    if(iorbp==iorbs .and. iorbq==iorbr) Hessian_pqrs=-Hessian_pqrs
   else
    Hessian_pqrs=(Energy_hh-Energy_h0-Energy_0h+two*Energy_00-Energy_mh0-Energy_0mh+Energy_mhmh)/(two*step*step)
   endif
  endif
  write(msg,'(a,f15.6)') 'Energy = ',Energy_00
  call write_output(msg)
  write(msg,'(a,i5,a,i5,a,f15.6)') 'Grad(',iorbp,',',iorbq,') = ',Gradient_pq
  call write_output(msg)
  write(msg,'(a,i5,a,i5,a,f15.6)') 'Grad(',iorbr,',',iorbs,') = ',Gradient_rs
  call write_output(msg)
  write(msg,'(a,i5,a,i5,a,i5,a,i5,a,f15.6)') 'Hessian(',iorbp,',',iorbq,';',iorbr,',',iorbs,') = ',Hessian_pqrs
  call write_output(msg)
  write(msg,'(a)') ' '
  call write_output(msg)
  ! Recover initial integrals and deallocate arrays
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, &
  & ERImol=INTEGd%ERImol)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  deallocate(U_mat,X_mat,NO_COEF_in)

 endif

end subroutine num_grad_hess_orb
!!***

!!***
!!****f* DoNOF/dm2_JK_3d
!! NAME
!! dm2_JK_3d
!!
!! FUNCTION
!!  Build the DM2_JK to be send to the main code for computing PI(r)
!!  (i.e. the on-top pair density)
!!
!! INPUTS
!! DM2_iiii=DM2 same orb elements
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DM2_L=DM2 elements that use L integrals 
!!
!! OUTPUT
!! DM2_JK
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_JK_3d(NBF_occ,DM2_J,DM2_K,DM2_L,DM2_iiii,DM2_JK)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_occ
!arrays
 real(dp),dimension(NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(NBF_occ,NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(2,NBF_occ,NBF_occ),intent(inout)::DM2_JK
!Local variables ------------------------------
!scalars
 integer::iorb
!arrays
!************************************************************************

 DM2_JK=zero

 DM2_JK(1,:,:)=DM2_J(:,:)
 DM2_JK(2,:,:)=DM2_K(:,:)+DM2_L(:,:) ! Time-rev. symmetry

 do iorb=1,NBF_occ
  DM2_JK(1,iorb,iorb)=DM2_iiii(iorb)
  DM2_JK(2,iorb,iorb)=zero
 enddo
 
!-----------------------------------------------------------------------
end subroutine dm2_JK_3d
!!***

end module m_optorb
!!***
