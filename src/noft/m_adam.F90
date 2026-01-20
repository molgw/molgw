!!****m* DoNOF/m_adam
!! NAME
!!  m_adam
!!
!! FUNCTION
!!  Module to build the ADAM 
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

module m_adam

 use m_nofoutput
 use m_definitions
 use m_rdmd
 use m_elag

 implicit none


!!***
!!****t* m_adam/adam_t
!! NAME
!! adam_t
!!
!! FUNCTION
!! Datatype building the (orbitals) ADAM matrix and gradients
!!
!! SOURCE

 type,public :: adam_t

  logical::restart=.false.      ! Comunicate that ADAM requires restart orb opt
  logical::cpx_adam=.false.     ! True for complex ADAM (i.e. complex orb)
  integer::icall_max_add=0      ! Increase icall_max
  integer::NDIM_adam            ! Size of the ADAM
  real(dp)::l_rate=0.01d0       ! Learning rate
  real(dp)::fact_rate=0.2d0     ! Learning rate factor
  real(dp)::beta1=0.7d0         ! beta1 parameter
  real(dp)::beta2=0.9d0         ! beta2 parameter
! arrays 
  real(dp),allocatable,dimension(:)::mom1_vec                 ! first_moment
  real(dp),allocatable,dimension(:)::mom2_vec                 ! second_moment
  real(dp),allocatable,dimension(:)::mom2_vec_max             ! second_moment_max
  complex(dp),allocatable,dimension(:)::mom1_vec_cmplx
  complex(dp),allocatable,dimension(:)::mom2_vec_cmplx
  complex(dp),allocatable,dimension(:)::mom2_vec_max_cmplx

 contains 
   procedure :: free => adam_free
   ! Destructor.

   procedure :: build => build_adam
   ! Use ELAGd to build the ADAM approx to kappa.

   procedure :: clean => clean_adam
   ! Use ELAGd to build the ADAM approx to kappa.
   
 end type adam_t

 public :: adam_init 
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/adam_init
!! NAME
!! adam_init
!!
!! FUNCTION
!!  Initialize the data type adam_t 
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

subroutine adam_init(ADAMd,NBF_tot,cpx_mos)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::cpx_mos
 integer,intent(in)::NBF_tot
 type(adam_t),intent(inout)::ADAMd
!Local variables ------------------------------
!scalars
 real(dp)::totMEM
!arrays
 character(len=200)::msg
!************************************************************************

 ADAMd%cpx_adam=cpx_mos
 ADAMd%NDIM_adam=NBF_tot*(NBF_tot-1)/2
 ! Calculate memory needed
 if(cpx_mos) then
  ADAMd%NDIM_adam=ADAMd%NDIM_adam+NBF_tot
  totMEM=8*ADAMd%NDIM_adam
 else
  totMEM=4*ADAMd%NDIM_adam
 endif
 totMEM=8*totMEM       ! Bytes
 totMEM=totMEM*tol6    ! Bytes to Mb  
 if(totMEM>thousand) then  ! Mb to Gb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing ADAMd object    ',totMEM*tol3,' Gb'
 elseif(totMEM<one) then   ! Mb to Kb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing ADAMd object    ',totMEM*thousand,' Kb'
 else                      ! Mb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing ADAMd object    ',totMEM,' Mb'
 endif
 call write_output(msg)
 ! Allocate arrays
 if(cpx_mos) then
  allocate(ADAMd%mom1_vec_cmplx(ADAMd%NDIM_adam))
  allocate(ADAMd%mom2_vec_cmplx(ADAMd%NDIM_adam))
  allocate(ADAMd%mom2_vec_max_cmplx(ADAMd%NDIM_adam))
  ADAMd%mom1_vec_cmplx=complex_zero;
  ADAMd%mom2_vec_cmplx=complex_zero;
  ADAMd%mom2_vec_max_cmplx=complex_zero;
 else 
  allocate(ADAMd%mom1_vec(ADAMd%NDIM_adam))
  allocate(ADAMd%mom2_vec(ADAMd%NDIM_adam))
  allocate(ADAMd%mom2_vec_max(ADAMd%NDIM_adam))
  ADAMd%mom1_vec=zero;
  ADAMd%mom2_vec=zero;
  ADAMd%mom2_vec_max=zero;
 endif
 
end subroutine adam_init
!!***

!!***
!!****f* DoNOF/adam_free
!! NAME
!! adam_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type adam_t 
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

subroutine adam_free(ADAMd)
!Arguments ------------------------------------
!scalars
 class(adam_t),intent(inout)::ADAMd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 if(ADAMd%cpx_adam) then
  deallocate(ADAMd%mom1_vec_cmplx) 
  deallocate(ADAMd%mom2_vec_cmplx) 
  deallocate(ADAMd%mom2_vec_max_cmplx) 
 else
  deallocate(ADAMd%mom1_vec) 
  deallocate(ADAMd%mom2_vec) 
  deallocate(ADAMd%mom2_vec_max) 
 endif

end subroutine adam_free
!!***

!!****f* DoNOF/build_adam
!! NAME
!! build_adam
!!
!! FUNCTION
!!  Build kappa with the ADAM.
!!
!! INPUTS
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  ELAGd=Object containg (orbital) Lagrange multipliers matrix (Lambda_pq)
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine build_adam(ADAMd,ELAGd,RDMd,icall,kappa_mat,kappa_mat_cmplx)
!Arguments ------------------------------------
!scalars
 integer,intent(inout)::icall
 class(adam_t),intent(inout)::ADAMd
 class(elag_t),intent(in)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),optional,intent(out)::kappa_mat
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),optional,intent(out)::kappa_mat_cmplx
!Local variables ------------------------------
!scalars
 integer::ipair
 integer::iorbp,iorbq
 real(dp)::grad_pq,mom1_hat,mom2_hat
 complex(dp)::grad_pq_cmplx,mom1_hat_cmplx,mom2_hat_cmplx
!arrays
! character(len=200)::msg
!************************************************************************

 icall=icall+1
 
 ! Build ADAM
 if(ADAMd%cpx_adam) then ! Complex

  kappa_mat_cmplx=complex_zero;

  ipair=0;
  do iorbp=1,RDMd%NBF_tot ! p
   do iorbq=iorbp,RDMd%NBF_tot ! q
    ipair=ipair+1
    !
    ! Calc. gradient
    !
    grad_pq_cmplx=(ELAGd%Lambdas(iorbq,iorbp)-ELAGd%Lambdas(iorbp,iorbq)) &
    &            +(ELAGd%Lambdas_im(iorbq,iorbp)+ELAGd%Lambdas_im(iorbp,iorbq))*im
    grad_pq_cmplx=two*grad_pq_cmplx
    ADAMd%mom1_vec_cmplx(ipair)=ADAMd%beta1*ADAMd%mom1_vec_cmplx(ipair)+(one-ADAMd%beta1)*grad_pq_cmplx
    ADAMd%mom2_vec_cmplx(ipair)=ADAMd%beta2*ADAMd%mom2_vec_cmplx(ipair)+(one-ADAMd%beta2)*grad_pq_cmplx*conjg(grad_pq_cmplx)
    mom1_hat_cmplx=ADAMd%mom1_vec_cmplx(ipair)/(one-ADAMd%beta1**icall)
    mom2_hat_cmplx=ADAMd%mom2_vec_cmplx(ipair)/(one-ADAMd%beta2**icall)
    ADAMd%mom2_vec_max_cmplx(ipair)=max(real(ADAMd%mom2_vec_max_cmplx(ipair)),real(mom2_hat_cmplx))
    kappa_mat_cmplx(iorbp,iorbq)=-ADAMd%l_rate*mom1_hat_cmplx/(sqrt(ADAMd%mom2_vec_max_cmplx(ipair))+1d-15)
    kappa_mat_cmplx(iorbq,iorbp)=-conjg(kappa_mat_cmplx(iorbp,iorbq))
   enddo
  enddo

 else ! Real

  kappa_mat=zero;

  ipair=0;
  do iorbp=1,RDMd%NBF_tot ! p
   do iorbq=iorbp+1,RDMd%NBF_tot ! q
    ipair=ipair+1
    !
    ! Calc. gradient
    !   Note: The k_pp does not have a real part 
    !
    grad_pq=ELAGd%Lambdas(iorbq,iorbp)-ELAGd%Lambdas(iorbp,iorbq)
    grad_pq=two*grad_pq
    ADAMd%mom1_vec(ipair)=ADAMd%beta1*ADAMd%mom1_vec(ipair)+(one-ADAMd%beta1)*grad_pq
    ADAMd%mom2_vec(ipair)=ADAMd%beta2*ADAMd%mom2_vec(ipair)+(one-ADAMd%beta2)*grad_pq**two
    mom1_hat=ADAMd%mom1_vec(ipair)/(one-ADAMd%beta1**icall)
    mom2_hat=ADAMd%mom2_vec(ipair)/(one-ADAMd%beta2**icall)
    ADAMd%mom2_vec_max(ipair)=max(ADAMd%mom2_vec_max(ipair),mom2_hat)
    kappa_mat(iorbp,iorbq)=-ADAMd%l_rate*mom1_hat/(sqrt(ADAMd%mom2_vec_max(ipair))+1d-15)
    kappa_mat(iorbq,iorbp)=-kappa_mat(iorbp,iorbq)
   enddo
  enddo

 endif

end subroutine build_adam
!!***

!!****f* DoNOF/clean_adam
!! NAME
!! clean_adam
!!
!! FUNCTION
!!  Set to 0 the ADAM arrays.
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

subroutine clean_adam(ADAMd)
!Arguments ------------------------------------
!scalars
 class(adam_t),intent(inout)::ADAMd
!arrays
!Local variables ------------------------------
!scalars
!arrays
! character(len=200)::msg
!************************************************************************

 ADAMd%restart=.false. ! ADAM will allow convergece for small energy differences
  
 ! Clean ADAM arrays
 if(ADAMd%cpx_adam) then ! Complex

  ADAMd%mom1_vec_cmplx=complex_zero;
  ADAMd%mom2_vec_cmplx=complex_zero;
  ADAMd%mom2_vec_max_cmplx=complex_zero;

 else ! Real

  ADAMd%mom1_vec=zero;
  ADAMd%mom2_vec=zero;
  ADAMd%mom2_vec_max=zero;
 endif

end subroutine clean_adam
!!***


end module m_adam
!!***
