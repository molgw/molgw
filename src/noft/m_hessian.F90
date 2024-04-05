!!****m* DoNOF/m_hessian
!! NAME
!!  m_hessian
!!
!! FUNCTION
!!  Module to build the Hessian 
!!
!! COPYRIGHT
!! This file is distributed under the terms of the
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!!
!! PARENTS
!!  m_optorb
!!  m_noft_driver
!!
!! CHILDREN
!!  m_rdmd
!!  m_elag
!!
!! SOURCE

module m_hessian

 use m_nofoutput
 use m_definitions
 use m_rdmd
 use m_elag
 use m_integd

 implicit none


!!***
!!****t* m_hessian/hessian_t
!! NAME
!! hessian_t
!!
!! FUNCTION
!! Datatype building the (orbitals) Hessian matrix and gradients
!!
!! SOURCE

 type,public :: hessian_t

  logical::cpx_hessian=.false.  ! True for complex Hessian (i.e. complex orbitals)
  integer::NDIM_hess            ! Size of the HESSIAN
! arrays 
  real(dp),allocatable,dimension(:)::Gradient_vec             ! F_pq - F_qp (Gradient matrix)
  complex(dp),allocatable,dimension(:)::Gradient_vec_cmplx    ! F_pq - F_qp* (Gradient matrix, complex)
  real(dp),allocatable,dimension(:,:)::Hessian_mat            ! Hessian matrix 
  complex(dp),allocatable,dimension(:,:)::Hessian_mat_cmplx   ! Hessian matrix (complex)


 contains 
   procedure :: free => hessian_free
   ! Destructor.

   procedure :: build => build_hessian
   ! Use integrals and the 1,2-RDM to build the Hessian and Gradient.
   
   procedure :: build_brut => build_hessian_brut
   ! Use integrals and the 1,2-RDM to build the Hessian and Gradient.

   procedure :: diag => diag_hessian
   ! Diagonalize the Hessian matrix and analyze the eigenvalues

 end type hessian_t

 public :: hessian_init 
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/hessian_init
!! NAME
!! hessian_init
!!
!! FUNCTION
!!  Initialize the data type hessian_t 
!!
!! INPUTS
!! NBF_tot=Number of total orbitals
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine hessian_init(HESSIANd,NBF_tot,cpx_mos)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::cpx_mos
 integer,intent(in)::NBF_tot
 type(hessian_t),intent(inout)::HESSIANd
!Local variables ------------------------------
!scalars
 real(dp)::totMEM
!arrays
 character(len=200)::msg
!************************************************************************

 HESSIANd%cpx_hessian=cpx_mos
 HESSIANd%NDIM_hess=NBF_tot*NBF_tot
 ! Calculate memory needed
 if(cpx_mos) then
  totMEM=2*(NBF_tot*NBF_tot*NBF_tot*NBF_tot+NBF_tot*NBF_tot)
 else
  totMEM=NBF_tot*NBF_tot*NBF_tot*NBF_tot+NBF_tot*NBF_tot
 endif
 totMEM=8*totMEM       ! Bytes
 totMEM=totMEM*tol6    ! Bytes to Mb  
 if(totMEM>thousand) then  ! Mb to Gb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing HESSIANd object ',totMEM*tol3,' Gb'
 elseif(totMEM<one) then   ! Mb to Kb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing HESSIANd object ',totMEM*thousand,' Kb'
 else                      ! Mb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing HESSIANd object ',totMEM,' Mb'
 endif
 call write_output(msg)
 ! Allocate arrays
 if(cpx_mos) then
  allocate(HESSIANd%Gradient_vec_cmplx(HESSIANd%NDIM_hess))
  allocate(HESSIANd%Hessian_mat_cmplx(HESSIANd%NDIM_hess,HESSIANd%NDIM_hess)) 
  HESSIANd%Gradient_vec_cmplx=complex_zero; HESSIANd%Hessian_mat_cmplx=complex_zero;
 else 
  allocate(HESSIANd%Gradient_vec(HESSIANd%NDIM_hess))
  allocate(HESSIANd%Hessian_mat(HESSIANd%NDIM_hess,HESSIANd%NDIM_hess)) 
  HESSIANd%Gradient_vec=zero; HESSIANd%Hessian_mat=zero;
 endif
 
end subroutine hessian_init
!!***

!!***
!!****f* DoNOF/hessian_free
!! NAME
!! hessian_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type hessian_t 
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

subroutine hessian_free(HESSIANd)
!Arguments ------------------------------------
!scalars
 class(hessian_t),intent(inout)::HESSIANd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 if(HESSIANd%cpx_hessian) then
  deallocate(HESSIANd%Gradient_vec_cmplx) 
  deallocate(HESSIANd%Hessian_mat_cmplx) 
 else
  deallocate(HESSIANd%Gradient_vec) 
  deallocate(HESSIANd%Hessian_mat) 
 endif

end subroutine hessian_free
!!***

!!****f* DoNOF/build_hessian
!! NAME
!! build_hessian
!!
!! FUNCTION
!!  Build the Hessian matrix and Gradient vector. Note that the electron rep. integrals are given in DIRAC's format
!!
!!  For complex orbs. with time-reversal symmetry [i.e. states p_alpha = (p_beta)* ]: 
!!      < p_alpha q_beta | r_alpha s_beta > =  < p_alpha s_alpha | r_alpha q_alpha >
!!
!! INPUTS
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  INTEGd=Object containg all integrals
!!  ELAGd=Object containg (orbital) Lagrange multipliers matrix (Lambda_pq)
!!
!! OUTPUT
!!  HESSIANd%Hessian_mat=Matrix build with the Hessian
!!  HESSIANd%Gradient_vec=Vector build with the Gradient
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine build_hessian(HESSIANd,ELAGd,RDMd,INTEGd,DM2_J,DM2_K,DM2_L)
!Arguments ------------------------------------
!scalars
 class(hessian_t),intent(inout)::HESSIANd
 class(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
!Local variables ------------------------------
!scalars
 integer::iorbp,iorbq,iorbr,iorbs,iorbt
 real(dp)::G_pqrs,grad_pq
 complex(dp)::G_pqrs_cmplx,grad_pq_cmplx
!arrays
 character(len=200)::msg
!************************************************************************
 
 if(HESSIANd%cpx_hessian) then
  write(msg,'(a)') 'Computing the complex Hessian Matrix'
  call write_output(msg)
  do iorbp=1,RDMd%NBF_tot ! p
   do iorbq=1,RDMd%NBF_tot ! q
    grad_pq_cmplx=two*(ELAGd%Lambdas(iorbq,iorbp)-ELAGd%Lambdas(iorbp,iorbq)) &
  &              +two*im*(ELAGd%Lambdas_im(iorbq,iorbp)+ELAGd%Lambdas_im(iorbp,iorbq))
    ! TODO
    HESSIANd%Gradient_vec_cmplx(iorbp+RDMd%NBF_tot*(iorbq-1))=grad_pq_cmplx
   enddo
  enddo
 else
  write(msg,'(a)') 'Computing the real Hessian matrix '
  call write_output(msg)
  do iorbp=1,RDMd%NBF_tot ! p
   do iorbq=1,RDMd%NBF_tot ! q
    grad_pq=two*(ELAGd%Lambdas(iorbq,iorbp)-ELAGd%Lambdas(iorbp,iorbq))
    do iorbr=1,RDMd%NBF_tot ! r
     do iorbs=1,RDMd%NBF_tot ! s
      G_pqrs=zero
      if(iorbq==iorbs) then ! q=s
       G_pqrs=G_pqrs-half*(ELAGd%Lambdas(iorbr,iorbp)+ELAGd%Lambdas(iorbp,iorbr))
       if(iorbq<=RDMd%NBF_occ) then ! q is occ
        G_pqrs=G_pqrs+two*RDMd%occ(iorbq)*INTEGd%hCORE(iorbp,iorbr)
        do iorbt=1,RDMd%NBF_occ ! t
         G_pqrs=G_pqrs+two*DM2_J(iorbq,iorbt) &
        & *(INTEGd%ERImol(iorbp,iorbt,iorbr,iorbt)-INTEGd%ERImol(iorbp,iorbt,iorbt,iorbr))
        enddo
       endif
      endif
      if(iorbp==iorbr) then ! p=r
       G_pqrs=G_pqrs-half*(ELAGd%Lambdas(iorbq,iorbs)+ELAGd%Lambdas(iorbs,iorbq))
       if(iorbp<=RDMd%NBF_occ) then ! p is occ
        G_pqrs=G_pqrs+two*RDMd%occ(iorbp)*INTEGd%hCORE(iorbs,iorbq)
        do iorbt=1,RDMd%NBF_occ ! t
         G_pqrs=G_pqrs+two*DM2_J(iorbp,iorbt) &
        & *(INTEGd%ERImol(iorbs,iorbt,iorbq,iorbt)-INTEGd%ERImol(iorbs,iorbt,iorbt,iorbq))
        enddo
       endif
      endif
      if(iorbp==iorbs .and. iorbp<=RDMd%NBF_occ) then ! p=s and p is occ
       do iorbt=1,RDMd%NBF_occ ! t
        G_pqrs=G_pqrs-two*DM2_L(iorbt,iorbp)*INTEGd%ERImol(iorbt,iorbt,iorbq,iorbr)
       enddo
       G_pqrs=G_pqrs-two*RDMd%DM2_iiii(iorbp)*INTEGd%ERImol(iorbp,iorbp,iorbq,iorbr)
      endif
      if(iorbq==iorbr .and. iorbq<=RDMd%NBF_occ) then ! q=r and q is occ
       do iorbt=1,RDMd%NBF_occ ! t
        G_pqrs=G_pqrs-two*DM2_L(iorbq,iorbt)*INTEGd%ERImol(iorbp,iorbs,iorbt,iorbt)
       enddo
       G_pqrs=G_pqrs-two*RDMd%DM2_iiii(iorbq)*INTEGd%ERImol(iorbp,iorbs,iorbq,iorbq)
      endif
      if(iorbp<=RDMd%NBF_occ .and. iorbs<=RDMd%NBF_occ) then  ! p & s are occ
       G_pqrs=G_pqrs-two*DM2_J(iorbp,iorbs)*INTEGd%ERImol(iorbp,iorbs,iorbq,iorbr) &
      &             -two*DM2_K(iorbp,iorbs)*INTEGd%ERImol(iorbs,iorbp,iorbq,iorbr) 
      endif
      if(iorbq<=RDMd%NBF_occ .and. iorbr<=RDMd%NBF_occ) then ! q & r are occ
       G_pqrs=G_pqrs-two*DM2_J(iorbq,iorbr)*INTEGd%ERImol(iorbp,iorbs,iorbq,iorbr) &
      &             -two*DM2_K(iorbq,iorbr)*INTEGD%ERImol(iorbp,iorbs,iorbr,iorbq)
      endif
      if(iorbp<=RDMd%NBF_occ .and. iorbr<=RDMd%NBF_occ) then  ! p & r are occ
       G_pqrs=G_pqrs+two*DM2_K(iorbp,iorbr)    &
      & *(INTEGd%ERImol(iorbs,iorbp,iorbq,iorbr)-INTEGd%ERImol(iorbs,iorbp,iorbr,iorbq))
       G_pqrs=G_pqrs+two*DM2_L(iorbr,iorbp)    &
      & *(INTEGd%ERImol(iorbs,iorbr,iorbq,iorbp)-INTEGd%ERImol(iorbs,iorbr,iorbp,iorbq))
       if(iorbp==iorbr) then ! p=r
        G_pqrs=G_pqrs+two*RDMd%DM2_iiii(iorbp) &
      & *(INTEGd%ERImol(iorbs,iorbp,iorbq,iorbp)-INTEGd%ERImol(iorbs,iorbp,iorbp,iorbq))
       endif
      endif
      if(iorbq<=RDMd%NBF_occ .and. iorbs<=RDMd%NBF_occ) then ! q & s are occ
       G_pqrs=G_pqrs+two*DM2_K(iorbq,iorbs) &
      & *(INTEGd%ERImol(iorbp,iorbs,iorbr,iorbq)-INTEGd%ERImol(iorbp,iorbs,iorbq,iorbr))
       G_pqrs=G_pqrs+two*DM2_L(iorbq,iorbs) &
      & *(INTEGd%ERImol(iorbp,iorbq,iorbr,iorbs)-INTEGd%ERImol(iorbp,iorbq,iorbs,iorbr))
       if(iorbq==iorbs) then ! q=s
        G_pqrs=G_pqrs+two*RDMd%DM2_iiii(iorbq) &
      & *(INTEGd%ERImol(iorbp,iorbq,iorbr,iorbq)-INTEGd%ERImol(iorbp,iorbq,iorbq,iorbr))
       endif
      endif
      HESSIANd%Hessian_mat(iorbp+RDMd%NBF_tot*(iorbq-1),iorbr+RDMd%NBF_tot*(iorbs-1))=G_pqrs
      ! write(*,*) iorbp,iorbq,iorbr,iorbs,G_pqrs
     enddo
    enddo
    HESSIANd%Gradient_vec(iorbp+RDMd%NBF_tot*(iorbq-1))=grad_pq
    ! write(*,*) iorbp,iorbq,grad_pq
   enddo
  enddo
 endif

! Check if the Hessian is a symmetric matrix
! 
! do iorbp=1,HESSIANd%NDIM_hess
!  do iorbq=1,HESSIANd%NDIM_hess
!   if(abs(HESSIANd%Hessian_mat(iorbp,iorbq)-HESSIANd%Hessian_mat(iorbq,iorbp))>tol8) then
!    write(*,*) iorbp,iorbq,HESSIANd%Hessian_mat(iorbp,iorbq),HESSIANd%Hessian_mat(iorbq,iorbp)
!   endif
!  enddo
! enddo

end subroutine build_hessian
!!***

!!****f* DoNOF/diag_hessian
!! NAME
!! diag_hessian
!!
!! FUNCTION
!!  Diagonalize the Hessian matrix and evaluate the eigenvalues
!!
!! INPUTS
!!
!! OUTPUT
!!  HESSIANd%Hessian_mat=Eigenvectors 
!!  HESSIANd%Hessian_mat_cmplx=Eigenvectors 
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine diag_hessian(HESSIANd)
!Arguments ------------------------------------
!scalars
 class(hessian_t),intent(inout)::HESSIANd
!arrays
!Local variables ------------------------------
 integer::iindex,lwork,info,nneg
 real(dp)::sum_neg,max_neg
!scalars
 real(dp),allocatable,dimension(:)::Work,RWork,Eigeval
 complex(dp),allocatable,dimension(:)::Work_cmplx
!arrays
 character(len=200)::msg
!************************************************************************

 max_neg=zero; sum_neg=zero; nneg=0;

 allocate(Eigeval(HESSIANd%NDIM_hess))
 Eigeval=zero
 
 if(HESSIANd%cpx_hessian) then
  write(*,*) 'banana'
 else
  allocate(Work(1))
  lwork=-1
  call DSYEV('V','L',HESSIANd%NDIM_hess,HESSIANd%Hessian_mat,HESSIANd%NDIM_hess,Eigeval,Work,lwork,info)
  lwork=nint(Work(1))
  if(info==0) then
   deallocate(Work)
   allocate(Work(lwork))
   call DSYEV('V','L',HESSIANd%NDIM_hess,HESSIANd%Hessian_mat,HESSIANd%NDIM_hess,Eigeval,Work,lwork,info)
  endif
  deallocate(Work)
 endif

 do iindex=1,HESSIANd%NDIM_hess
  if(abs(Eigeval(iindex))<tol6) then
   Eigeval(iindex)=zero
  endif
  if(Eigeval(iindex)<zero) then
   nneg=nneg+1
   sum_neg=sum_neg+Eigeval(iindex)
   if(Eigeval(iindex)<max_neg) max_neg=Eigeval(iindex)
  endif
 enddo
 
 write(msg,'(a,I10)') 'Number of Negative eigenvalues',nneg
 call write_output(msg)
 write(msg,'(a,F10.5)') 'Max Negative eigenvalue       ',max_neg
 call write_output(msg)
 write(msg,'(a,F10.5)') 'Sum Negative eigenvalues      ',sum_neg
 call write_output(msg)

 deallocate(Eigeval)

end subroutine diag_hessian
!!***

!!****f* DoNOF/build_hessian_brut
!! NAME
!! build_hessian_brut
!!
!! FUNCTION
!!  Build the Hessian matrix and Gradient vector. Note that the electron rep. integrals are given in DIRAC's format
!!
!!  For complex orbs. with time-reversal symmetry [i.e. states p_alpha = (p_beta)* ]: 
!!      < p_alpha q_beta | r_alpha s_beta > =  < p_alpha s_alpha | r_alpha q_alpha >
!!
!! INPUTS
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  INTEGd=Object containg all integrals
!!
!! OUTPUT
!!  HESSIANd%Hessian_mat=Matrix build with the Hessian
!!  HESSIANd%Gradient_vec=Vector build with the Gradient
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine build_hessian_brut(HESSIANd,RDMd,INTEGd,DM2_J,DM2_K,DM2_L)
!Arguments ------------------------------------
!scalars
 class(hessian_t),intent(inout)::HESSIANd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
!Local variables ------------------------------
!scalars
 integer::iorbp,iorbq,iorbr,iorbs,iorbt,iorbu,iorbv
 real(dp)::G_pqrs,G_qprs,G_pqsr,G_qpsr,grad_pq!,h_core,vee
!arrays
 real(dp),allocatable,dimension(:,:)::DM1
 real(dp),allocatable,dimension(:,:,:,:)::DM2
 character(len=200)::msg
!************************************************************************

 !h_core=zero; vee=zero; 

 allocate(DM1(RDMd%NBF_tot,RDMd%NBF_tot),DM2(RDMd%NBF_tot,RDMd%NBF_tot,RDMd%NBF_tot,RDMd%NBF_tot))
 DM1=zero; DM2=zero;
 ! Build full 1,2-RDMs (this could be done outside this subroutine)
 ! Tr[1D] = N/2 , Tr[2D]= N(N-1)/2
 do iorbp=1,RDMd%NBF_occ ! p
  DM1(iorbp,iorbp)=RDMd%occ(iorbp)
  do iorbq=1,RDMd%NBF_occ ! q
   DM2(iorbp,iorbq,iorbp,iorbq)=DM2_J(iorbp,iorbq)
   DM2(iorbp,iorbq,iorbq,iorbp)=DM2_K(iorbp,iorbq)
   DM2(iorbp,iorbp,iorbq,iorbq)=DM2_L(iorbp,iorbq)
  enddo
  DM2(iorbp,iorbp,iorbp,iorbp)=RDMd%DM2_iiii(iorbp)
 enddo

 ! Build Hessian
 if(.not.HESSIANd%cpx_hessian) then
  do iorbp=1,RDMd%NBF_tot ! p
   do iorbq=1,RDMd%NBF_tot ! q
!    h_core=h_core+DM1(iorbp,iorbq)*INTEGd%Hcore(iorbp,iorbq)
    !
    ! Calc. gradient
    !
    grad_pq=zero
    do iorbt=1,RDMd%NBF_tot !t
     grad_pq=grad_pq+two*DM1(iorbt,iorbq)*INTEGd%Hcore(iorbt,iorbp) &
    &               -two*DM1(iorbp,iorbt)*INTEGd%Hcore(iorbq,iorbt)
     do iorbu=1,RDMd%NBF_tot ! u
      do iorbv=1,RDMd%NBF_tot ! v
       grad_pq=grad_pq+two*DM2(iorbu,iorbt,iorbv,iorbq)*INTEGd%ERImol(iorbu,iorbt,iorbv,iorbp) &
       &              -two*DM2(iorbp,iorbu,iorbt,iorbv)*INTEGd%ERImol(iorbq,iorbu,iorbt,iorbv) 
      enddo
     enddo
    enddo
    ! 
    !
    do iorbr=1,RDMd%NBF_tot ! r
     do iorbs=1,RDMd%NBF_tot ! s
!      vee=vee+half*DM2(iorbp,iorbq,iorbr,iorbs)*INTEGd%ERImol(iorbp,iorbq,iorbr,iorbs)
      ! G_pqrs
      G_pqrs=-two*(DM1(iorbr,iorbq)*INTEGd%Hcore(iorbs,iorbp)+DM1(iorbp,iorbs)*INTEGd%Hcore(iorbq,iorbr))
      if(iorbr==iorbq) then ! q=r
       do iorbt=1,RDMd%NBF_tot !t
        G_pqrs=G_pqrs+DM1(iorbt,iorbs)*INTEGd%Hcore(iorbt,iorbp) &
       &             +DM1(iorbp,iorbt)*INTEGd%Hcore(iorbs,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_pqrs=G_pqrs+DM2(iorbu,iorbv,iorbs,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbp,iorbt) &
          &            +DM2(iorbp,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbs,iorbt,iorbu,iorbv) 
         enddo
        enddo
       enddo 
      endif
      if(iorbp==iorbs) then ! p=s
       do iorbt=1,RDMd%NBF_tot !t
        G_pqrs=G_pqrs+DM1(iorbt,iorbq)*INTEGd%Hcore(iorbt,iorbr) &
       &             +DM1(iorbr,iorbt)*INTEGd%Hcore(iorbq,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_pqrs=G_pqrs+DM2(iorbr,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbq,iorbt,iorbu,iorbv) &
          &            +DM2(iorbu,iorbv,iorbq,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbr,iorbt) 
         enddo
        enddo
       enddo 
      endif
      do iorbt=1,RDMd%NBF_tot !t
       do iorbu=1,RDMd%NBF_tot !u
        G_pqrs=G_pqrs-two*DM2(iorbr,iorbt,iorbq,iorbu)*INTEGd%ERImol(iorbs,iorbt,iorbp,iorbu) &
        &            -two*DM2(iorbt,iorbr,iorbq,iorbu)*INTEGd%ERImol(iorbt,iorbs,iorbp,iorbu) &
        &            -two*DM2(iorbp,iorbu,iorbs,iorbt)*INTEGd%ERImol(iorbq,iorbu,iorbr,iorbt) &
        &            -two*DM2(iorbp,iorbu,iorbt,iorbs)*INTEGd%ERImol(iorbq,iorbu,iorbt,iorbr) &
        &            +two*DM2(iorbt,iorbu,iorbq,iorbs)*INTEGd%ERImol(iorbt,iorbu,iorbp,iorbr) &
        &            +two*DM2(iorbp,iorbr,iorbt,iorbu)*INTEGd%ERImol(iorbq,iorbs,iorbt,iorbu) 
       enddo
      enddo
      ! G_qprs
      G_qprs=two*(DM1(iorbr,iorbp)*INTEGd%Hcore(iorbs,iorbq)+DM1(iorbq,iorbs)*INTEGd%Hcore(iorbp,iorbr))
      if(iorbr==iorbp) then ! p=r
       do iorbt=1,RDMd%NBF_tot !t
        G_qprs=G_qprs-DM1(iorbt,iorbs)*INTEGd%Hcore(iorbt,iorbq) &
       &             -DM1(iorbq,iorbt)*INTEGd%Hcore(iorbs,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_qprs=G_qprs-DM2(iorbu,iorbv,iorbs,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbq,iorbt) &
          &            -DM2(iorbq,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbs,iorbt,iorbu,iorbv)
         enddo
        enddo
       enddo
      endif
      if(iorbq==iorbs) then ! q=s
       do iorbt=1,RDMd%NBF_tot !t
        G_qprs=G_qprs-DM1(iorbt,iorbp)*INTEGd%Hcore(iorbt,iorbr) &
       &             -DM1(iorbr,iorbt)*INTEGd%Hcore(iorbp,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_qprs=G_qprs-DM2(iorbr,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbp,iorbt,iorbu,iorbv) &
          &            -DM2(iorbu,iorbv,iorbp,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbr,iorbt)
         enddo
        enddo
       enddo
      endif
      do iorbt=1,RDMd%NBF_tot !t
       do iorbu=1,RDMd%NBF_tot !u
        G_qprs=G_qprs+two*DM2(iorbr,iorbt,iorbp,iorbu)*INTEGd%ERImol(iorbs,iorbt,iorbq,iorbu) &
        &            +two*DM2(iorbt,iorbr,iorbp,iorbu)*INTEGd%ERImol(iorbt,iorbs,iorbq,iorbu) &
        &            +two*DM2(iorbq,iorbu,iorbs,iorbt)*INTEGd%ERImol(iorbp,iorbu,iorbr,iorbt) &
        &            +two*DM2(iorbq,iorbu,iorbt,iorbs)*INTEGd%ERImol(iorbp,iorbu,iorbt,iorbr) &
        &            -two*DM2(iorbt,iorbu,iorbp,iorbs)*INTEGd%ERImol(iorbt,iorbu,iorbq,iorbr) &
        &            -two*DM2(iorbq,iorbr,iorbt,iorbu)*INTEGd%ERImol(iorbp,iorbs,iorbt,iorbu)
       enddo
      enddo
      ! G_pqsr
      G_pqsr=two*(DM1(iorbs,iorbq)*INTEGd%Hcore(iorbr,iorbp)+DM1(iorbp,iorbr)*INTEGd%Hcore(iorbq,iorbs))
      if(iorbs==iorbq) then ! q=s
       do iorbt=1,RDMd%NBF_tot !t
        G_pqsr=G_pqsr-DM1(iorbt,iorbr)*INTEGd%Hcore(iorbt,iorbp) &
       &             -DM1(iorbp,iorbt)*INTEGd%Hcore(iorbr,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_pqsr=G_pqsr-DM2(iorbu,iorbv,iorbr,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbp,iorbt) &
          &            -DM2(iorbp,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbr,iorbt,iorbu,iorbv)
         enddo
        enddo
       enddo
      endif
      if(iorbp==iorbr) then ! p=r
       do iorbt=1,RDMd%NBF_tot !t
        G_pqsr=G_pqsr-DM1(iorbt,iorbq)*INTEGd%Hcore(iorbt,iorbs) &
       &             -DM1(iorbs,iorbt)*INTEGd%Hcore(iorbq,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_pqsr=G_pqsr-DM2(iorbs,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbq,iorbt,iorbu,iorbv) &
          &            -DM2(iorbu,iorbv,iorbq,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbs,iorbt)
         enddo
        enddo
       enddo
      endif
      do iorbt=1,RDMd%NBF_tot !t
       do iorbu=1,RDMd%NBF_tot !u
        G_pqsr=G_pqsr+two*DM2(iorbs,iorbt,iorbq,iorbu)*INTEGd%ERImol(iorbr,iorbt,iorbp,iorbu) &
        &            +two*DM2(iorbt,iorbs,iorbq,iorbu)*INTEGd%ERImol(iorbt,iorbr,iorbp,iorbu) &
        &            +two*DM2(iorbp,iorbu,iorbr,iorbt)*INTEGd%ERImol(iorbq,iorbu,iorbs,iorbt) &
        &            +two*DM2(iorbp,iorbu,iorbt,iorbr)*INTEGd%ERImol(iorbq,iorbu,iorbt,iorbs) &
        &            -two*DM2(iorbt,iorbu,iorbq,iorbr)*INTEGd%ERImol(iorbt,iorbu,iorbp,iorbs) &
        &            -two*DM2(iorbp,iorbs,iorbt,iorbu)*INTEGd%ERImol(iorbq,iorbr,iorbt,iorbu)
       enddo
      enddo
      ! G_qpsr
      G_qpsr=-two*(DM1(iorbs,iorbp)*INTEGd%Hcore(iorbr,iorbq)+DM1(iorbq,iorbr)*INTEGd%Hcore(iorbp,iorbs))
      if(iorbs==iorbp) then ! p=s
       do iorbt=1,RDMd%NBF_tot !t
        G_qpsr=G_qpsr+DM1(iorbt,iorbr)*INTEGd%Hcore(iorbt,iorbq) &
       &             +DM1(iorbq,iorbt)*INTEGd%Hcore(iorbr,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_qpsr=G_qpsr+DM2(iorbu,iorbv,iorbr,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbq,iorbt) &
          &            +DM2(iorbq,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbr,iorbt,iorbu,iorbv)
         enddo
        enddo
       enddo
      endif
      if(iorbq==iorbr) then ! q=r
       do iorbt=1,RDMd%NBF_tot !t
        G_qpsr=G_qpsr+DM1(iorbt,iorbp)*INTEGd%Hcore(iorbt,iorbs) &
       &             +DM1(iorbs,iorbt)*INTEGd%Hcore(iorbp,iorbt)
        do iorbu=1,RDMd%NBF_tot ! u
         do iorbv=1,RDMd%NBF_tot ! v
          G_qpsr=G_qpsr+DM2(iorbs,iorbt,iorbu,iorbv)*INTEGd%ERImol(iorbp,iorbt,iorbu,iorbv) &
          &            +DM2(iorbu,iorbv,iorbp,iorbt)*INTEGd%ERImol(iorbu,iorbv,iorbs,iorbt)
         enddo
        enddo
       enddo
      endif
      do iorbt=1,RDMd%NBF_tot !t
       do iorbu=1,RDMd%NBF_tot !u
        G_qpsr=G_qpsr-two*DM2(iorbs,iorbt,iorbp,iorbu)*INTEGd%ERImol(iorbr,iorbt,iorbq,iorbu) &
        &            -two*DM2(iorbt,iorbs,iorbp,iorbu)*INTEGd%ERImol(iorbt,iorbr,iorbq,iorbu) &
        &            -two*DM2(iorbq,iorbu,iorbr,iorbt)*INTEGd%ERImol(iorbp,iorbu,iorbs,iorbt) &
        &            -two*DM2(iorbq,iorbu,iorbt,iorbr)*INTEGd%ERImol(iorbp,iorbu,iorbt,iorbs) &
        &            +two*DM2(iorbt,iorbu,iorbp,iorbr)*INTEGd%ERImol(iorbt,iorbu,iorbq,iorbs) &
        &            +two*DM2(iorbq,iorbs,iorbt,iorbu)*INTEGd%ERImol(iorbp,iorbr,iorbt,iorbu)
       enddo
      enddo
      HESSIANd%Hessian_mat(iorbp+RDMd%NBF_tot*(iorbq-1),iorbr+RDMd%NBF_tot*(iorbs-1))=G_pqrs+G_qprs+G_pqsr+G_qpsr
      ! write(*,*) iorbp,iorbq,iorbr,iorbs,G_pqrs+G_qprs+G_pqsr+G_qpsr
     enddo
    enddo
    HESSIANd%Gradient_vec(iorbp+RDMd%NBF_tot*(iorbq-1))=two*grad_pq
    ! write(*,*) iorbp,iorbq,two*grad_pq
   enddo
  enddo
 endif

 deallocate(DM1,DM2)

! Check that energy contributions are fine
!write(*,*) h_core,vee

! Check if the Hessian is a symmetric matrix
 
! do iorbp=1,HESSIANd%NDIM_hess
!  do iorbq=1,HESSIANd%NDIM_hess
!   if(abs(HESSIANd%Hessian_mat(iorbp,iorbq)-HESSIANd%Hessian_mat(iorbq,iorbp))>tol8) then
!    write(*,*) iorbp,iorbq,HESSIANd%Hessian_mat(iorbp,iorbq),HESSIANd%Hessian_mat(iorbq,iorbp)
!   endif
!  enddo
! enddo

end subroutine build_hessian_brut
!!***

end module m_hessian
!!***
