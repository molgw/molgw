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
 use m_diagf
 use m_e_grad_occ
 use m_e_grad_occ_cpx

 implicit none

 private :: lambda_conv
!!***

 public :: opt_orb,s2_calc 
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

subroutine opt_orb(iter,imethod,ELAGd,RDMd,INTEGd,Vnn,Energy,mo_ints,NO_COEF,NO_COEF_cmplx) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::iter,imethod
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 interface 
  subroutine mo_ints(NBF_tot,NBF_occ,NBF_jkl,Occ,NO_COEF,hCORE,ERImol,ERImolJsr,ERImolLsr,&
  & NO_COEF_cmplx,hCORE_cmplx,ERImol_cmplx)
  use m_definitions
  implicit none
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
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF_cmplx
!Local variables ------------------------------
!scalars
 logical::convLambda,nogamma,diddiis
 integer::icall,iorbmax1,iorbmax2
 real(dp)::sumdiff,maxdiff,maxdiff_all,Ediff,Energy_old
!arrays
 character(len=200)::msg
!************************************************************************

 Ediff=zero; Energy=zero; Energy_old=zero; convLambda=.false.;nogamma=.true.;
 if((imethod==1).and.(iter==0)) then
  ELAGd%sumdiff_old=zero
 endif
 
 icall=0
 if(INTEGd%complex_ints) then
  call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_cmplx, &
  & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_old,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,nogamma=nogamma)
 else
  if(INTEGd%irange_sep/=0) then
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, & 
   & ERImol=INTEGd%ERImol,ERImolJsr=INTEGd%ERImolJsr,ERImolLsr=INTEGd%ERImolLsr)
  else
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, & 
   & ERImol=INTEGd%ERImol)
  endif
  call INTEGd%eritoeriJKL(RDMd%NBF_occ)
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_old,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
  & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
 endif
 do
  ! If we used a DIIS step, do not stop after DIIS for small Energy dif.
  diddiis=.false.

  ! Build Lambda matrix
  call ELAGd%build(RDMd,INTEGd,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,RDMd%DM2_Jsr,RDMd%DM2_Lsr)

  ! Check if these NO_COEF with the RDMs are already the solution =)
  call lambda_conv(ELAGd,RDMd,convLambda,sumdiff,maxdiff,maxdiff_all,iorbmax1,iorbmax2)
  if(convLambda) then
   write(msg,'(a)') 'Lambda_qp - Lambda_pq* converged for the Hemiticty of Lambda'
   call write_output(msg)
   exit
  else
   if(imethod==1.and.icall==0) then                                        ! F method: adjust MaxScaling for the rest of orb. icall iterations
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
  if(imethod==1) then ! Build F matrix for iterative diagonalization
   if(INTEGd%complex_ints) then
    call diagF_to_coef(iter,icall,maxdiff,diddiis,ELAGd,RDMd,NO_COEF_cmplx=NO_COEF_cmplx) ! Build new NO_COEF and set icall=icall+1
   else
    call diagF_to_coef(iter,icall,maxdiff,diddiis,ELAGd,RDMd,NO_COEF=NO_COEF) ! Build new NO_COEF and set icall=icall+1
   endif
  else                ! Use Newton method to produce new COEFs
   ! TODO
  endif

  ! Build all integrals in the new NO_COEF basis (including arrays for ERI_J and ERI_K)
  if(INTEGd%complex_ints) then
   call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF_cmplx=NO_COEF_cmplx, &
   & hCORE_cmplx=INTEGd%hCORE_cmplx,ERImol_cmplx=INTEGd%ERImol_cmplx)
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,nogamma=nogamma)
  else
   if(INTEGd%irange_sep/=0) then
    call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, &
    & ERImol=INTEGd%ERImol,ERImolJsr=INTEGd%ERImolJsr,ERImolLsr=INTEGd%ERImolLsr)
   else
    call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,RDMd%occ,NO_COEF=NO_COEF,hCORE=INTEGd%hCORE, &
    & ERImol=INTEGd%ERImol)
   endif
   call INTEGd%eritoeriJKL(RDMd%NBF_occ)
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
  endif
 
  ! Check if we did Diag[(Lambda_pq + Lambda_qp*)/2] for F method (first iteration)
  if((imethod==1.and.iter==0).and.ELAGd%diagLpL_done) then ! For F method if we did Diag[(Lambda_pq + Lambda_qp*)/2].
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

  ! We allow at most 30 generations of new NO_COEF updates (and integrals)
  if(icall==30) exit
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 enddo
 
 ! Calc. the final Energy using fixed RDMs and the new NO_COEF (before going back to occ. optimization)
 if(INTEGd%complex_ints) then
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,nogamma=nogamma)
 else
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
  & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
 endif
 if(icall>0) then
  write(msg,'(a,f15.6,a,i6,a)') 'Orb. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
  call write_output(msg)
  if(imethod==1.and.iter>0) then
   write(msg,'(a,f15.6,a,i5,a,i5,a)') 'Max. [Lambda_qp - Lambda_pq*]= ',maxdiff_all,' pair (',iorbmax1,',',iorbmax2,')'
   call write_output(msg)
   write(msg,'(a,f19.10)') 'Energy difference orb. opt.=',Ediff
   call write_output(msg)
  endif
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
   diff=dabs(ELAGd%Lambdas(iorb1,iorb)-ELAGd%Lambdas(iorb,iorb1))
   if(ELAGd%cpx_lambdas .and. iorb/=iorb1) then
    diff=diff*diff+(ELAGd%Lambdas_im(iorb1,iorb)+ELAGd%Lambdas_im(iorb,iorb1))**two
    diff=dsqrt(diff)
   endif
   if(ELAGd%cpx_lambdas .and. iorb1==iorb) then
    diff=dabs(ELAGd%Lambdas_im(iorb1,iorb)+ELAGd%Lambdas_im(iorb,iorb1))
   endif
   sumdiff=sumdiff+diff
   if((diff>=tol_dif_Lambda) .and. converg_lamb) then
    converg_lamb=.false.
   endif
   if(ELAGd%cpx_lambdas) then
    if(diff>maxdiff .and.iorb/=iorb1) then ! TODO : check for complex if we really need this
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

end module m_optorb
!!***
