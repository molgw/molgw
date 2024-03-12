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
 integer::iorb,iorb1,iorb2,iorb3,iorb4
 real(dp)::hess_pqrs,grad_pq
 complex(dp)::hess_pqrs_cmplx,grad_pq_cmplx
!arrays
 character(len=200)::msg
!************************************************************************
 
 if(HESSIANd%cpx_hessian) then
  write(msg,'(a)') 'Computing the complex Hessian Matrix'
  call write_output(msg)
  do iorb=1,RDMd%NBF_tot ! p
   do iorb1=1,RDMd%NBF_tot ! q
    grad_pq_cmplx=ELAGd%Lambdas(iorb1,iorb)-ELAGd%Lambdas(iorb,iorb1) &
  &              +im*(ELAGd%Lambdas_im(iorb1,iorb)+ELAGd%Lambdas_im(iorb,iorb1))
    ! TODO
    HESSIANd%Gradient_vec_cmplx(iorb+RDMd%NBF_tot*(iorb1-1))=grad_pq_cmplx
   enddo
  enddo
 else
  write(msg,'(a)') 'Computing the real Hessian matrix '
  call write_output(msg)
  do iorb=1,RDMd%NBF_tot ! p
   do iorb1=1,RDMd%NBF_tot ! q
    grad_pq=ELAGd%Lambdas(iorb1,iorb)-ELAGd%Lambdas(iorb,iorb1)
    do iorb2=1,RDMd%NBF_tot ! r
     do iorb3=1,RDMd%NBF_tot ! s
      hess_pqrs=zero
      if(iorb1==iorb3) then ! q=s
       hess_pqrs=hess_pqrs-(ELAGd%Lambdas(iorb2,iorb)+ELAGd%Lambdas(iorb,iorb2))
       if(iorb1<=RDMd%NBF_occ) then ! q is occ
        hess_pqrs=hess_pqrs+two*RDMd%occ(iorb1)*INTEGd%hCORE(iorb,iorb2)
        do iorb4=1,RDMd%NBF_occ ! t
         hess_pqrs=hess_pqrs+four*DM2_J(iorb1,iorb4) &
        & *(INTEGd%ERImol(iorb,iorb4,iorb2,iorb4)-INTEGd%ERImol(iorb,iorb4,iorb4,iorb2))
        enddo
       endif
      endif
      if(iorb==iorb2) then ! p=r
       hess_pqrs=hess_pqrs-(ELAGd%Lambdas(iorb1,iorb3)+ELAGd%Lambdas(iorb3,iorb1))
       if(iorb<=RDMd%NBF_occ) then ! p is occ
        hess_pqrs=hess_pqrs+two*RDMd%occ(iorb)*INTEGd%hCORE(iorb3,iorb1)
        do iorb4=1,RDMd%NBF_occ ! t
         hess_pqrs=hess_pqrs+four*DM2_J(iorb,iorb4) &
        & *(INTEGd%ERImol(iorb3,iorb4,iorb1,iorb4)-INTEGd%ERImol(iorb3,iorb4,iorb4,iorb1))
        enddo
       endif
      endif
      if(iorb==iorb3 .and. iorb<=RDMd%NBF_occ) then ! p=s and p is occ
       do iorb4=1,RDMd%NBF_occ ! t
        hess_pqrs=hess_pqrs-four*DM2_L(iorb4,iorb)*INTEGd%ERImol(iorb4,iorb4,iorb1,iorb2)
       enddo
       hess_pqrs=hess_pqrs-four*RDMd%DM2_iiii(iorb)*INTEGd%ERImol(iorb,iorb,iorb1,iorb2)
      endif
      if(iorb1==iorb2 .and. iorb1<=RDMd%NBF_occ) then ! q=r and q is occ
       do iorb4=1,RDMd%NBF_occ ! t
        hess_pqrs=hess_pqrs-four*DM2_L(iorb1,iorb4)*INTEGd%ERImol(iorb,iorb3,iorb4,iorb4)
       enddo
       hess_pqrs=hess_pqrs-four*RDMd%DM2_iiii(iorb1)*INTEGd%ERImol(iorb,iorb3,iorb1,iorb1)
      endif
      if(iorb<=RDMd%NBF_occ .and. iorb3<=RDMd%NBF_occ) then  ! p & s are occ
       hess_pqrs=hess_pqrs-four*DM2_J(iorb,iorb3)*INTEGd%ERImol(iorb,iorb3,iorb1,iorb2) &
      &                   -four*DM2_K(iorb,iorb3)*INTEGd%ERImol(iorb3,iorb,iorb1,iorb2) 
      endif
      if(iorb1<=RDMd%NBF_occ .and. iorb2<=RDMd%NBF_occ) then ! q & r are occ
       hess_pqrs=hess_pqrs-four*DM2_J(iorb1,iorb2)*INTEGd%ERImol(iorb,iorb3,iorb1,iorb2) &
      &                   -four*DM2_K(iorb1,iorb2)*INTEGD%ERImol(iorb,iorb3,iorb2,iorb1)
      endif
      if(iorb<=RDMd%NBF_occ .and. iorb2<=RDMd%NBF_occ) then  ! p & r are occ
       hess_pqrs=hess_pqrs+four*DM2_K(iorb,iorb2)    &
      & *(INTEGd%ERImol(iorb3,iorb,iorb1,iorb2)-INTEGd%ERImol(iorb3,iorb,iorb2,iorb1))
       hess_pqrs=hess_pqrs+four*DM2_L(iorb2,iorb)    &
      & *(INTEGd%ERImol(iorb3,iorb2,iorb1,iorb)-INTEGd%ERImol(iorb3,iorb2,iorb,iorb1))
       if(iorb==iorb2) then ! p=r
        hess_pqrs=hess_pqrs+four*RDMd%DM2_iiii(iorb) &
      & *(INTEGd%ERImol(iorb3,iorb,iorb1,iorb)-INTEGd%ERImol(iorb3,iorb,iorb,iorb1))
       endif
      endif
      if(iorb1<=RDMd%NBF_occ .and. iorb3<=RDMd%NBF_occ) then ! q & s are occ
       hess_pqrs=hess_pqrs+four*DM2_K(iorb1,iorb3) &
      & *(INTEGd%ERImol(iorb,iorb3,iorb2,iorb1)-INTEGd%ERImol(iorb,iorb3,iorb1,iorb2))
       hess_pqrs=hess_pqrs+four*DM2_L(iorb1,iorb3) &
      & *(INTEGd%ERImol(iorb,iorb1,iorb2,iorb3)-INTEGd%ERImol(iorb,iorb1,iorb3,iorb2))
       if(iorb1==iorb3) then ! q=s
        hess_pqrs=hess_pqrs+four*RDMd%DM2_iiii(iorb1) &
      & *(INTEGd%ERImol(iorb,iorb1,iorb2,iorb1)-INTEGd%ERImol(iorb,iorb1,iorb1,iorb2))
       endif
      endif
      HESSIANd%Hessian_mat(iorb+RDMd%NBF_tot*(iorb1-1),iorb2+RDMd%NBF_tot*(iorb3-1))=hess_pqrs
     enddo
    enddo
    HESSIANd%Gradient_vec(iorb+RDMd%NBF_tot*(iorb1-1))=grad_pq
   enddo
  enddo
 endif

! Check if the Hessian is a symmetric matrix
! 
! do iorb=1,HESSIANd%NDIM_hess
!  do iorb1=1,HESSIANd%NDIM_hess
!   if(abs(HESSIANd%Hessian_mat(iorb,iorb1)-HESSIANd%Hessian_mat(iorb1,iorb))>tol8) then
!    write(*,*) iorb,iorb1,HESSIANd%Hessian_mat(iorb,iorb1),HESSIANd%Hessian_mat(iorb1,iorb)
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
  if(abs(Eigeval(iindex))<tol8) then
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

end module m_hessian
!!***
