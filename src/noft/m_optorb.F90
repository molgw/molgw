!!****m* DoNOF/m_optorb
!! NAME
!!  m_optorb
!!
!! FUNCTION
!!  Module prepared to compute orb optimization for a fixed set of OCCs and 2-RDM
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
 use m_diagf
 use m_e_grad_occ

 implicit none

 private :: lambda_conv
!!***

 public :: opt_orb 
!!***

contains

!!***
!!****f* DoNOF/opt_orb
!! NAME
!!  opt_orb
!!
!! FUNCTION
!!  Call the F-matrix method or Newton (Hessian) for orb optimization 
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

subroutine opt_orb(iter,imethod,ELAGd,RDMd,INTEGd,Vnn,Energy,NO_COEF,mo_ints) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::iter,imethod
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 logical::convLambda,nogamma,diddiis
 integer::icall
 real(dp)::sumdiff,maxdiff,Ediff,Energy_old
!arrays
 character(len=200)::msg
!************************************************************************

 Energy=0.0d0; Energy_old=0.0d0; convLambda=.false.;nogamma=.true.;
 if((imethod==1).and.(iter==0)) then
  ELAGd%sumdiff_old=0.0d0
 endif
 
 icall=0
 call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,NO_COEF,INTEGd%hCORE,INTEGd%ERImol)
 call INTEGd%eritoeriJK(RDMd%NBF_occ)
 do
  ! If we used a DIIS step, do not stop after DIIS for small Energy dif.
  diddiis=.false.

  ! Build Lambda matrix
  call ELAGd%build(RDMd,INTEGd,RDMd%DM2_J,RDMd%DM2_K)

  ! Check if these NO_COEF with the RDMs are already the solution =)
  call lambda_conv(ELAGd,RDMd,convLambda,sumdiff,maxdiff)
  if(convLambda) then
   write(msg,'(a)') 'Lambda_pq - Lambda_qp* converged for the Hemiticty of Lambda'
   call write_output(msg)
   exit
  else
   if(imethod==1.and.icall==0) then                                        ! F method: adjust MaxScaling for the rest of orb. icall iterations
    if(iter>2.and.iter>ELAGd%itscale.and.(sumdiff>ELAGd%sumdiff_old)) then ! Parameters chosen empirically (i.e. experience) to
     ELAGd%itscale=iter+10                                                 ! ensure convergence. Maybe we can set them as input variables?
     ELAGd%MaxScaling=ELAGd%MaxScaling+1
     if(ELAGd%MaxScaling>ELAGd%itolLambda) then
      ELAGd%MaxScaling=2                                                   ! One more empirical/experience parameter =(
     endif
    endif
    ELAGd%sumdiff_old=sumdiff
   endif
  endif

  ! Update NO_COEF
  if(imethod==1) then ! Build F matrix for iterative diagonalization
   call diagF_to_coef(iter,icall,maxdiff,diddiis,ELAGd,RDMd,NO_COEF) ! Build new NO_COEF and set icall=icall+1
  else                ! Use Newton method to produce new COEFs
   ! TODO
  endif

  ! Build all integrals in the new NO_COEF basis (including arrays for ERI_J and ERI_K)
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,NO_COEF,INTEGd%hCORE,INTEGd%ERImol)
  call INTEGd%eritoeriJK(RDMd%NBF_occ)
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K,nogamma=nogamma)
 
  ! Check if we did Diag[(Lambda_pq + Lambda_qp*)/2] for F method (first iteration)
  if((imethod==1.and.iter==0).and.ELAGd%diagLpL_done) then ! For F method if we did Diag[(Lambda_pq + Lambda_qp*)/2].
   exit                                                    ! -> Do only one icall iteration before the occ. opt.
  endif

  ! For this icall using the new NO_COEF (and fixed RDMs). Is the Energy still changing?
  Ediff=Energy_old-Energy
  if((icall>1).and.(dabs(Ediff)<ELAGd%tolE).and.(.not.diddiis)) then ! The energy is not changing anymore (not stopping for DIIS)
   write(msg,'(a)') 'Lambda_pq - Lambda_qp* converged for small energy differences'
   call write_output(msg)
   exit
  endif
  Energy_old=Energy

  ! We allow at most 30 generations of new NO_COEF updates (and integrals)
  if(icall==30) exit
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 enddo
 
 ! Calc. the final Energy using fixed RDMs and the new NO_COEF (before going back to occ. optimization)
 call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K,nogamma=nogamma)
 write(msg,'(a,f15.6,a,i6,a)') 'Orb. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
 call write_output(msg)
 if(imethod==1.and.iter>0) then
  write(msg,'(a,f15.6)') 'Max. [Lambda_pq - Lambda_qp*]= ',maxdiff
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
!!  Check if the Lambda matrix already fulfils the condition Lambda_pq-Lambda_qp^* <= tol_dif_lambda.
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

subroutine lambda_conv(ELAGd,RDMd,converg_lamb,sumdiff,maxdiff)
!Arguments ------------------------------------
!scalars
 logical,intent(inout)::converg_lamb
 real(dp),intent(inout)::sumdiff,maxdiff
 type(elag_t),intent(in)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 real(dp)::diff,tol_dif_Lambda
!arrays
!************************************************************************

 tol_dif_Lambda=1.0d1**(-ELAGd%itolLambda)
 converg_lamb=.true.; sumdiff=0.0d0; maxdiff=0.0d0;
 
 do iorb=1,RDMd%NBF_tot
  do iorb1=1,RDMd%NBF_tot
   diff=dabs(ELAGd%Lambdas(iorb1,iorb)-ELAGd%Lambdas(iorb,iorb1))
   sumdiff=sumdiff+diff
   if((diff>=tol_dif_Lambda) .and. converg_lamb) then
    converg_lamb=.false.
   endif
   if(diff>maxdiff) then
    maxdiff=diff
   endif
  enddo
 enddo

end subroutine lambda_conv
!!***

end module m_optorb
!!***
