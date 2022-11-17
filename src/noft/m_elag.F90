!!****m* DoNOF/m_elag
!! NAME
!!  m_elag
!!
!! FUNCTION
!!  Module to build the Lagrange multipliers Lambda 
!!
!! COPYRIGHT
!! This file is distributed under the terms of the
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!!
!! PARENTS
!!  m_optorb
!!
!! CHILDREN
!!  m_rdmd
!!
!! SOURCE

module m_elag

 use m_nofoutput
 use m_definitions
 use m_rdmd
 use m_integd

 implicit none

 private :: dyson_orbs
!!***
!!****t* m_rdmd/rdm_t
!! NAME
!! rdm_t
!!
!! FUNCTION
!! Datatype storing noft quantities and arrays needed
!!
!! SOURCE

 type,public :: elag_t

  logical::cpx_lambdas=.false.  ! True for complex Lambdas (i.e. complex orbitals)
  logical::diagLpL=.true.       ! Do the diag. using (lambda+lambda)/2?
  logical::diagLpL_done=.false. ! Did we use use (lambda+lambda)/2?
  integer::imethod=1            ! Method used for optimization (1-> Diag F matrix)
  integer::MaxScaling=0         ! Max scaling reductions employed to avoid divergence of diag[F]
  integer::itscale=1            ! Above this number of iterations we do MaxScaling=MaxScaling+1
  integer::itolLambda=5         ! Integer used to define 10**-itolLambda as threshold of Lambda_qp-Lambda_pq* convergence
  integer::itoldiis=3           ! Integer used to define 10**-itoldiis as threshold of DIIS trigger
  integer::idiis=0              ! Current DIIS iteration
  integer::ndiis=5              ! The number of iterations required to apply DIIS is ndiis+1
  integer::ndiis_array          ! Size of the arrays used in DIIS (ndiis+2)
  real(dp)::sumdiff_old         ! Old value of sum_pq |F_pq|  for p/=q 
  real(dp)::tolE                ! Tolerance that will be imposed in Energy convergence
! arrays 
  real(dp),allocatable,dimension(:)::F_diag         ! F_pp (Diag. part of the F matrix)
  real(dp),allocatable,dimension(:)::Coef_DIIS      ! DIIS coefs. used to build linear comb. of F matrices
  real(dp),allocatable,dimension(:,:)::Lambdas      ! Lambda_pq (Lagrange multipliers matrix)
  real(dp),allocatable,dimension(:,:)::Lambdas_im   ! Lambda_pq (Lagrange multipliers matrix imag part)
  real(dp),allocatable,dimension(:,:)::DIIS_mat     ! DIIS matrix used to solve the system of eqs. DIIS_MAT*Coef_DIIS = (0 0 ... 0 1) 
  real(dp),allocatable,dimension(:,:,:)::F_DIIS     ! F matrices used by DIIS
  complex(dp),allocatable,dimension(:)::Coef_DIIS_cmplx  ! DIIS coefs. used to build linear comb. of F matrices (complex)
  complex(dp),allocatable,dimension(:,:)::DIIS_mat_cmplx ! DIIS matrix used to solve the system of eqs. DIIS_MAT*Coef_DIIS = (0 0 ... 0 1) (complex)
  complex(dp),allocatable,dimension(:,:,:)::F_DIIS_cmplx ! F matrices used by DIIS (complex)


 contains 
   procedure :: free => elag_free
   ! Destructor.

   procedure :: build => build_elag
   ! Use integrals and the 1,2-RDM to build Lambdas matrix.

   procedure :: diag_lag => diag_lambda_ekt
   ! Diagonalize the matrix Lambdas (or divided by occ. numbers) to compute canonical orbs. or EKT.

   procedure :: clean_diis => wipeout_diis
   ! Set to ZERO all arrays employed by DIIS.

   procedure :: print_Fdiag => print_F_diag
   ! Print the F_diag vector to un unformated file.

 end type elag_t

 public :: elag_init 
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/elag_init
!! NAME
!! elag_init
!!
!! FUNCTION
!!  Initialize the data type elag_t 
!!
!! INPUTS
!! NBF_tot=Number of total orbitals
!! diagLpL_in=Diagonalize 0.5 (Lambda+Lambda) for the first iteration?
!! itolLambda=Used as 10**-itolLambda to check for Lambda_qp-Lambda_pq* convergence
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine elag_init(ELAGd,NBF_tot,diagLpL_in,itolLambda_in,ndiis_in,imethod_in,tolE_in,cpx_mos)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::diagLpL_in,cpx_mos
 integer,intent(in)::NBF_tot,itolLambda_in,ndiis_in,imethod_in
 real(dp),intent(in)::tolE_in
 type(elag_t),intent(inout)::ELAGd
!Local variables ------------------------------
!scalars
 real(dp)::totMEM
!arrays
 character(len=200)::msg
!************************************************************************

 ELAGd%cpx_lambdas=cpx_mos
 ELAGd%imethod=imethod_in
 ELAGd%itolLambda=itolLambda_in
 ELAGd%diagLpL=diagLpL_in
 ELAGd%ndiis=ndiis_in
 ELAGd%ndiis_array=ELAGd%ndiis+2
 ELAGd%tolE=tolE_in
 ! Calculate memory needed
 if(cpx_mos) then
  totMEM=NBF_tot+2*NBF_tot*NBF_tot
  if(ELAGd%ndiis>0.and.ELAGd%imethod==1) then
   totMEM=totMEM+2*ELAGd%ndiis_array+8*ELAGd%ndiis_array*NBF_tot*NBF_tot+4*ELAGd%ndiis_array*ELAGd%ndiis_array
  endif
 else
  totMEM=NBF_tot+NBF_tot*NBF_tot
  if(ELAGd%ndiis>0.and.ELAGd%imethod==1) then
   totMEM=totMEM+ELAGd%ndiis_array+ELAGd%ndiis_array*NBF_tot*NBF_tot+ELAGd%ndiis_array*ELAGd%ndiis_array
  endif
 endif
 totMEM=8*totMEM       ! Bytes
 totMEM=totMEM*tol6    ! Bytes to Mb  
 if(totMEM>thousand) then  ! Mb to Gb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing ELAGd object  ',totMEM*tol3,' Gb'
 elseif(totMEM<one) then   ! Mb to Kb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing ELAGd object  ',totMEM*thousand,' Kb'
 else                      ! Mb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing ELAGd object  ',totMEM,' Mb'
 endif
 call write_output(msg)
 ! Allocate arrays
 allocate(ELAGd%F_diag(NBF_tot))
 allocate(ELAGd%Lambdas(NBF_tot,NBF_tot)) 
 if(cpx_mos) then
  allocate(ELAGd%Lambdas_im(NBF_tot,NBF_tot)) 
  if(ELAGd%ndiis>0.and.ELAGd%imethod==1) then
   allocate(ELAGd%Coef_DIIS_cmplx(ELAGd%ndiis_array))
   allocate(ELAGd%F_DIIS_cmplx(ELAGd%ndiis_array,NBF_tot,NBF_tot))
   allocate(ELAGd%DIIS_mat_cmplx(ELAGd%ndiis_array,ELAGd%ndiis_array)) 
  endif
 else 
  if(ELAGd%ndiis>0.and.ELAGd%imethod==1) then
   allocate(ELAGd%Coef_DIIS(ELAGd%ndiis_array))
   allocate(ELAGd%F_DIIS(ELAGd%ndiis_array,NBF_tot,NBF_tot))
   allocate(ELAGd%DIIS_mat(ELAGd%ndiis_array,ELAGd%ndiis_array)) 
  endif
 endif
 
end subroutine elag_init
!!***

!!***
!!****f* DoNOF/elag_free
!! NAME
!! elag_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type elag_t 
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

subroutine elag_free(ELAGd)
!Arguments ------------------------------------
!scalars
 class(elag_t),intent(inout)::ELAGd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 deallocate(ELAGd%F_diag) 
 deallocate(ELAGd%Lambdas)
 
 if(ELAGd%cpx_lambdas) then
  deallocate(ELAGd%Lambdas_im)
  if(ELAGd%ndiis>0.and.ELAGd%imethod==1) then
   deallocate(ELAGd%Coef_DIIS_cmplx)
   deallocate(ELAGd%F_DIIS_cmplx)
   deallocate(ELAGd%DIIS_mat_cmplx) 
  endif 
 else
  if(ELAGd%ndiis>0.and.ELAGd%imethod==1) then
   deallocate(ELAGd%Coef_DIIS)
   deallocate(ELAGd%F_DIIS)
   deallocate(ELAGd%DIIS_mat) 
  endif 
 endif

end subroutine elag_free
!!***

!!****f* DoNOF/build_elag
!! NAME
!! build_elag
!!
!! FUNCTION
!!  Build the Lagrange multipliers Lambda matrix. Nothe that the electron rep. integrals are given in DoNOF format
!!
!!  For complex orbs. with time-reversal symmetry [i.e. states p_alpha = (p_beta)* ]: 
!!      < p_alpha q_beta | r_alpha s_beta > =  < p_alpha s_alpha | r_alpha q_beta >
!!             and      Lij integral (alpha beta) = Kij integral (alpha alpha)
!!
!! INPUTS
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  INTEGd=Object containg all integrals
!!
!! OUTPUT
!!  ELAGd%Lambdas=Matrix build with the Lagrange multipliers Lambda_pq
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine build_elag(ELAGd,RDMd,INTEGd,DM2_J,DM2_K,DM2_L,DM2_H)
!Arguments ------------------------------------
!scalars
 class(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L,DM2_H
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorbv,iorbv1
!arrays
 character(len=200)::msg
!************************************************************************

 ELAGd%Lambdas=zero
 if(ELAGd%cpx_lambdas) then
  ELAGd%Lambdas_im=zero
  do iorb=1,RDMd%NBF_occ
   ELAGd%Lambdas(iorb,:)=RDMd%occ(iorb)*real(INTEGd%hCORE_cmplx(:,iorb))                                        ! Init: Lambda_pq = n_p hCORE_qp
   ELAGd%Lambdas_im(iorb,:)=RDMd%occ(iorb)*aimag(INTEGd%hCORE_cmplx(:,iorb))                                    ! Init: Lambda_pq = n_p hCORE_qp
   if(INTEGd%iERItyp/=-1) then
    ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+RDMd%DM2_iiii(iorb)*real(INTEGd%ERImol_cmplx(:,iorb,iorb,iorb))          ! any->iorb,iorb->iorb
    ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+RDMd%DM2_iiii(iorb)*aimag(INTEGd%ERImol_cmplx(:,iorb,iorb,iorb))   ! any->iorb,iorb->iorb
    do iorb1=1,RDMd%NBF_occ
     if(iorb/=iorb1) then
      if(INTEGd%iERItyp==0) then ! DoNOF notation {ij|lk}
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb,iorb1)) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb,iorb1)) ! any->iorb1,iorb->iorb1
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_J(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_K(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb,iorb1)) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_L(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb,iorb1)) ! any->iorb1,iorb->iorb1
      elseif(INTEGd%iERItyp==1) then ! <ij|kl>
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb,iorb1)) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb->iorb1
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_J(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb,iorb1)) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_K(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_L(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb->iorb1
      elseif(INTEGd%iERItyp==2) then ! (ik|jl)
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb,iorb1,iorb1)) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*real(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb->iorb1
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_J(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb,iorb1,iorb1)) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_K(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_L(iorb,iorb1)*aimag(INTEGd%ERImol_cmplx(:,iorb1,iorb1,iorb)) ! any->iorb1,iorb->iorb1
      else
       ! Nth
      endif
     endif
    enddo
   else
    iorbv= (iorb-1)*(INTEGd%NBF2+INTEGd%NBF3+INTEGd%NBF4)+1
    iorbv1=(iorb-1)*(INTEGd%NBF2+INTEGd%NBF3+INTEGd%NBF4)+INTEGd%NBF2
    ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+RDMd%DM2_iiii(iorb)*real(INTEGd%ERImolv_cmplx(iorbv:iorbv1))   ! any->iorb,iorb->iorb
    ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+RDMd%DM2_iiii(iorb)*aimag(INTEGd%ERImolv_cmplx(iorbv:iorbv1))   ! any->iorb,iorb->iorb
    do iorb1=1,RDMd%NBF_occ
     iorbv= (iorb1-1)*(INTEGd%NBF2+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF3+1
     iorbv1=(iorb1-1)*(INTEGd%NBF2+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF3+INTEGd%NBF2
     ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*real(INTEGd%ERImolv_cmplx(iorbv:iorbv1)) ! any->iorb,iorb1->iorb1
     ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_J(iorb,iorb1)*aimag(INTEGd%ERImolv_cmplx(iorbv:iorbv1)) ! any->iorb,iorb1->iorb1
     iorbv= (iorb1-1)*(INTEGd%NBF2+INTEGd%NBF3)+(iorb-1)*INTEGd%NBF4+1
     iorbv1=(iorb1-1)*(INTEGd%NBF2+INTEGd%NBF3)+(iorb-1)*INTEGd%NBF4+INTEGd%NBF2
     ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*real(INTEGd%ERImolv_cmplx(iorbv:iorbv1)) ! any->iorb1,iorb1->iorb
     ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_K(iorb,iorb1)*aimag(INTEGd%ERImolv_cmplx(iorbv:iorbv1)) ! any->iorb1,iorb1->iorb
     !iorbv= (iorb1-1)*(INTEGd%NBF3+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF2+1
     !iorbv1=(iorb1-1)*(INTEGd%NBF3+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF2+INTEGd%NBF2
     ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*real(INTEGd%ERImolv_cmplx(iorbv:iorbv1)) ! any->iorb1,iorb->iorb1
     ELAGd%Lambdas_im(iorb,:)=ELAGd%Lambdas_im(iorb,:)+DM2_L(iorb,iorb1)*aimag(INTEGd%ERImolv_cmplx(iorbv:iorbv1)) ! any->iorb1,iorb->iorb1
    enddo
   endif
  enddo
 else
  do iorb=1,RDMd%NBF_occ
   ELAGd%Lambdas(iorb,:)=RDMd%occ(iorb)*INTEGd%hCORE(:,iorb)                                          ! Init: Lambda_pq = n_p hCORE_qp
   if(INTEGd%iERItyp/=-1) then
    ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+RDMd%DM2_iiii(iorb)*INTEGd%ERImol(:,iorb,iorb,iorb)   ! any->iorb,iorb->iorb
    do iorb1=1,RDMd%NBF_occ
     if(INTEGd%range_sep) ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_H(iorb,iorb1)*INTEGd%ERImolH(:,iorb1,iorb) ! any->iorb,iorb1->iorb1 | rs-NOFT
     if(iorb/=iorb1) then
      if(INTEGd%iERItyp==0) then ! DoNOF notation {ij|lk}
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb1,iorb) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb,iorb1) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*INTEGd%ERImol(:,iorb,iorb1,iorb1) ! any->iorb1,iorb->iorb1
      elseif(INTEGd%iERItyp==1) then ! <ij|kl>
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb,iorb1) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb1,iorb) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*INTEGd%ERImol(:,iorb,iorb1,iorb1) ! any->iorb1,iorb->iorb1
      elseif(INTEGd%iERItyp==2) then ! (ik|jl)
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*INTEGd%ERImol(:,iorb,iorb1,iorb1) ! any->iorb,iorb1->iorb1
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb1,iorb) ! any->iorb1,iorb1->iorb
       ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb,iorb1) ! any->iorb1,iorb->iorb1
      else
       ! Nth
      endif
     endif
    enddo
   else
    iorbv= (iorb-1)*(INTEGd%NBF2+INTEGd%NBF3+INTEGd%NBF4)+1
    iorbv1=(iorb-1)*(INTEGd%NBF2+INTEGd%NBF3+INTEGd%NBF4)+INTEGd%NBF2
    ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+RDMd%DM2_iiii(iorb)*INTEGd%ERImolv(iorbv:iorbv1)   ! any->iorb,iorb->iorb
    do iorb1=1,RDMd%NBF_occ
     if(INTEGd%range_sep) ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_H(iorb,iorb1)*INTEGd%ERImolH(:,iorb1,iorb) ! any->iorb,iorb1->iorb1 | rs-NOFT
     iorbv= (iorb1-1)*(INTEGd%NBF2+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF3+1
     iorbv1=(iorb1-1)*(INTEGd%NBF2+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF3+INTEGd%NBF2
     ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*INTEGd%ERImolv(iorbv:iorbv1) ! any->iorb,iorb1->iorb1
     iorbv= (iorb1-1)*(INTEGd%NBF2+INTEGd%NBF3)+(iorb-1)*INTEGd%NBF4+1
     iorbv1=(iorb1-1)*(INTEGd%NBF2+INTEGd%NBF3)+(iorb-1)*INTEGd%NBF4+INTEGd%NBF2
     ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_K(iorb,iorb1)*INTEGd%ERImolv(iorbv:iorbv1) ! any->iorb1,iorb1->iorb
     iorbv= (iorb1-1)*(INTEGd%NBF3+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF2+1
     iorbv1=(iorb1-1)*(INTEGd%NBF3+INTEGd%NBF4)+(iorb-1)*INTEGd%NBF2+INTEGd%NBF2
     ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_L(iorb,iorb1)*INTEGd%ERImolv(iorbv:iorbv1) ! any->iorb1,iorb->iorb1
    enddo
   endif
  enddo 
 endif
 !ELAGd%Lambdas=two*ELAGd%Lambdas ! We only need half for 'alpha' orbs to define gradients

 ! TODO 
 if(RDMd%Nsingleocc>0) then
  write(msg,'(a)') 'Error! The Lambda_pq matrix construction is not implemented for Nsingleocc>0'
  call write_output(msg)
 endif

end subroutine build_elag
!!***

!!****f* DoNOF/diag_lambda_ekt
!! NAME
!! diag_lambda_ekt
!!
!! FUNCTION
!!  Diagonalize the Lagrange multipliers Lambda matrix (produce either the 'canonical orbitals' or EKT). 
!!
!! INPUTS
!!  ELAGd%Lambdas=Matrix containing the Lagrange multipliers Lambda_pq
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  INTEGd=Object containg all integrals
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine diag_lambda_ekt(ELAGd,RDMd,INTEGd,NO_COEF,NO_COEF_cmplx,ekt)
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::ekt
 class(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(in)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF_cmplx
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,lwork,info
 real(dp)::sqrt_occ_iorb,sqrt_occ_iorb1
!arrays
 character(len=10)::coef_file
 character(len=200)::msg
 real(dp),allocatable,dimension(:)::Eigval,Eigval_nocc,Work,RWork
 real(dp),allocatable,dimension(:,:)::Eigvec,CANON_COEF
 complex(dp),allocatable,dimension(:)::Work_cmplx
 complex(dp),allocatable,dimension(:,:)::Eigvec_cmplx,CANON_COEF_cmplx 
!************************************************************************

 allocate(Eigval_nocc(RDMd%NBF_occ),Eigval(RDMd%NBF_tot),RWork(3*RDMd%NBF_tot-2))
 RWork=zero
 
 allocate(Work(1),Work_cmplx(1))
 if(ELAGd%cpx_lambdas) then
  allocate(Eigvec(1,1),Eigvec_cmplx(RDMd%NBF_tot,RDMd%NBF_tot))
  Eigvec_cmplx=ELAGd%Lambdas+ELAGd%Lambdas_im*im
 else
  allocate(Eigvec(RDMd%NBF_tot,RDMd%NBF_tot),Eigvec_cmplx(1,1))
  Eigvec=ELAGd%Lambdas
 endif

 if(present(ekt)) then
  do iorb=1,RDMd%NBF_tot
   do iorb1=1,RDMd%NBF_tot
    if(iorb<=RDMd%NBF_occ.and.iorb1<=RDMd%NBF_occ) then
     sqrt_occ_iorb =dsqrt(RDMd%occ(iorb))
     sqrt_occ_iorb1=dsqrt(RDMd%occ(iorb1))
     if((dabs(sqrt_occ_iorb)>tol6).and.(dabs(sqrt_occ_iorb1)>tol6)) then
      if(ELAGd%cpx_lambdas) then
       Eigvec_cmplx(iorb,iorb1)=Eigvec_cmplx(iorb,iorb1)/(sqrt_occ_iorb*sqrt_occ_iorb1)
      else
       Eigvec(iorb,iorb1)=Eigvec(iorb,iorb1)/(sqrt_occ_iorb*sqrt_occ_iorb1)
      endif
     else
      if(ELAGd%cpx_lambdas) then
       Eigvec_cmplx(iorb,iorb1)=complex_zero
      else
       Eigvec(iorb,iorb1)=zero
      endif
     endif
    else
     if(ELAGd%cpx_lambdas) then
      Eigvec_cmplx(iorb,iorb1)=complex_zero
     else
      Eigvec(iorb,iorb1)=zero
     endif
    endif
   enddo
  enddo
 endif

 ! Diagonalize
 if(ELAGd%cpx_lambdas) then
  lwork=-1
  call ZHEEV('V','L',RDMd%NBF_tot,Eigvec_cmplx,RDMd%NBF_tot,Eigval,Work_cmplx,lwork,RWork,info)
  lwork=nint(real(Work_cmplx(1)))
  if(info==0) then
   deallocate(Work_cmplx)
   allocate(Work_cmplx(lwork))
   call ZHEEV('V','L',RDMd%NBF_tot,Eigvec_cmplx,RDMd%NBF_tot,Eigval,Work_cmplx,lwork,RWork,info)
  endif
 else
  lwork=-1
  call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,Eigval,Work,lwork,info)
  lwork=nint(Work(1))
  if(info==0) then
   deallocate(Work)
   allocate(Work(lwork)) 
   call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,Eigval,Work,lwork,info)
  endif
 endif

 ! Print final eigenvalues and orbs.
 write(msg,'(a)') ' '
 call write_output(msg)
 if(present(ekt)) then
  Eigval=-Eigval
  if(ELAGd%cpx_lambdas) then
   call dyson_orbs(RDMd,INTEGd,Eigvec_cmplx=Eigvec_cmplx,NO_COEF_cmplx=NO_COEF_cmplx)
  else
   call dyson_orbs(RDMd,INTEGd,Eigvec=Eigvec,NO_COEF=NO_COEF)
  endif
  write(msg,'(a)') 'EKT ionization potentials (a.u.)'
  call write_output(msg)
 else
  coef_file='CANON_COEF'
  if(ELAGd%cpx_lambdas) then
   allocate(CANON_COEF_cmplx(RDMd%NBF_tot,RDMd%NBF_tot))
   CANON_COEF_cmplx=matmul(NO_COEF_cmplx,Eigvec_cmplx)
   call RDMd%print_orbs(coef_file,COEF_cmplx=CANON_COEF_cmplx)
   deallocate(CANON_COEF_cmplx)
  else
   allocate(CANON_COEF(RDMd%NBF_tot,RDMd%NBF_tot))
   CANON_COEF=matmul(NO_COEF,Eigvec)
   call RDMd%print_orbs(coef_file,COEF=CANON_COEF)
   deallocate(CANON_COEF)
  endif
  write(msg,'(a)') 'Canonical orbital eigenvalues (a.u.)'
  call write_output(msg)
 endif

 Eigval_nocc(1:RDMd%NBF_occ)=Eigval(1:RDMd%NBF_occ)
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(msg,'(f12.6,9f11.6)') Eigval_nocc(iorb:iorb+9)
  call write_output(msg)
 enddo
 iorb=(RDMd%NBF_occ/10)*10+1
 write(msg,'(f12.6,*(f11.6))') Eigval_nocc(iorb:)
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
  
 deallocate(Eigvec,Eigvec_cmplx,Work)
 deallocate(Eigval,Eigval_nocc,RWork,Work_cmplx)

end subroutine diag_lambda_ekt
!!***

!!****f* DoNOF/dyson_orbs
!! NAME
!!  dyson_orbs   
!!
!! FUNCTION
!!  Compute Dyson orbitals and corrected occ numbers. 
!!
!! INPUTS
!!  RDMd=Object containg all required variables about sizes
!!  INTEGd=Object containg all integrals
!!  Eigvec=Eigenvectos obtaing from EKT diag.
!!  NO_COEF=Nat. orb. coefs
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine dyson_orbs(RDMd,INTEGd,Eigvec,Eigvec_cmplx,NO_COEF,NO_COEF_cmplx)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(in)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::Eigvec
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::Eigvec_cmplx
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF_cmplx
!Local variables ------------------------------
!scalars
 logical::cpx_mos=.false.
 integer::iorb,iorb1,iorb2
!arrays
 character(len=10)::coef_file
 character(len=200)::msg
 real(dp),allocatable,dimension(:)::DYSON_occ
 real(dp),allocatable,dimension(:,:)::DYSON_COEF
 complex(dp),allocatable,dimension(:)::DYSON_occ_cmplx
 complex(dp),allocatable,dimension(:,:)::DYSON_COEF_cmplx
!************************************************************************

 allocate(DYSON_occ(RDMd%NBF_occ))
 if(present(NO_COEF_cmplx).and.present(Eigvec_cmplx)) cpx_mos=.true.
 ! Compute DYSON_COEFs
 if(cpx_mos) then
  allocate(DYSON_occ_cmplx(RDMd%NBF_occ))
  allocate(DYSON_COEF_cmplx(RDMd%NBF_tot,RDMd%NBF_tot))
  DYSON_COEF_cmplx=complex_zero; DYSON_occ_cmplx=complex_zero;
  ! Unnormalized DYSON_COEF
  do iorb=1,RDMd%NBF_tot
   do iorb1=1,RDMd%NBF_occ
    DYSON_COEF_cmplx(iorb,iorb1)=complex_zero
    do iorb2=1,RDMd%NBF_occ
     DYSON_COEF_cmplx(iorb,iorb1)=DYSON_COEF_cmplx(iorb,iorb1)+&
     & dsqrt(RDMd%occ(iorb2))*NO_COEF_cmplx(iorb,iorb2)*Eigvec_cmplx(iorb2,iorb1)
    enddo
   enddo
  enddo
  ! Occ numbers of DYSON_COEF
  do iorb=1,RDMd%NBF_occ
   DYSON_occ_cmplx(iorb)=complex_zero
   do iorb1=1,RDMd%NBF_tot
    DYSON_occ_cmplx(iorb)=DYSON_occ_cmplx(iorb)+conjg(DYSON_COEF_cmplx(iorb1,iorb))&
   &     *sum(INTEGd%Overlap(iorb1,:)*DYSON_COEF_cmplx(:,iorb))
   enddo
  enddo
  ! Normalized DYSON_COEF
  do iorb=1,RDMd%NBF_occ
   do iorb1=1,RDMd%NBF_tot
    DYSON_COEF_cmplx(iorb1,iorb)=DYSON_COEF_cmplx(iorb1,iorb)/cdsqrt(DYSON_occ_cmplx(iorb))
   enddo
  enddo
  ! Print DYSON_COEF
  coef_file='DYSON_COEF'
  call RDMd%print_orbs(coef_file,COEF_cmplx=DYSON_COEF_cmplx)
  DYSON_occ(:)=real(DYSON_occ_cmplx(:))
  deallocate(DYSON_COEF_cmplx,DYSON_occ_cmplx)
 else
  allocate(DYSON_COEF(RDMd%NBF_tot,RDMd%NBF_tot))
  DYSON_COEF=zero; DYSON_occ=zero;
  ! Unnormalized DYSON_COEF
  do iorb=1,RDMd%NBF_tot
   do iorb1=1,RDMd%NBF_occ
    DYSON_COEF(iorb,iorb1)=zero
    do iorb2=1,RDMd%NBF_occ
     DYSON_COEF(iorb,iorb1)=DYSON_COEF(iorb,iorb1)+&
     & dsqrt(RDMd%occ(iorb2))*NO_COEF(iorb,iorb2)*Eigvec(iorb2,iorb1)
    enddo
   enddo
  enddo 
  ! Occ numbers of DYSON_COEF
  do iorb=1,RDMd%NBF_occ
   DYSON_occ(iorb)=zero
   do iorb1=1,RDMd%NBF_tot
    DYSON_occ(iorb)=DYSON_occ(iorb)+DYSON_COEF(iorb1,iorb)&
   &     *sum(INTEGd%Overlap(iorb1,:)*DYSON_COEF(:,iorb))
   enddo
  enddo
  ! Normalized DYSON_COEF
  do iorb=1,RDMd%NBF_occ
   do iorb1=1,RDMd%NBF_tot
    DYSON_COEF(iorb1,iorb)=DYSON_COEF(iorb1,iorb)/dsqrt(DYSON_occ(iorb))
   enddo
  enddo
  ! Print DYSON_COEF
  coef_file='DYSON_COEF'
  call RDMd%print_orbs(coef_file,COEF=DYSON_COEF)
  deallocate(DYSON_COEF)
 endif
 ! Print DYSON occ. numbers
 DYSON_occ(:)=two*DYSON_occ(:)
 write(msg,'(a,f10.5,a)') 'Dyson occ ',sum(DYSON_occ(:)),'. Dyson occ. numbers '
 call write_output(msg)
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(msg,'(f12.6,9f11.6)') DYSON_occ(iorb:iorb+9)
  call write_output(msg)
 enddo
 iorb=(RDMd%NBF_occ/10)*10+1
 write(msg,'(f12.6,*(f11.6))') DYSON_occ(iorb:)
 call write_output(msg)

 ! Deallocate array
 deallocate(DYSON_occ)

end subroutine dyson_orbs
!!***

!!****f* DoNOF/wipeout_diis
!! NAME
!! wipeout_diis
!!
!! FUNCTION
!!  Build the Lagrange multipliers Lambda matrix. Nothe that the electron rep. integrals are given in DoNOF format
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

subroutine wipeout_diis(ELAGd)
!Arguments ------------------------------------
!scalars
 class(elag_t),intent(inout)::ELAGd
!arrays
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 ELAGd%idiis=0
 if(ELAGd%cpx_lambdas) then
  if(ELAGd%ndiis>0) then
   ELAGd%Coef_DIIS_cmplx=complex_zero
   ELAGd%F_DIIS_cmplx=complex_zero
   ELAGd%DIIS_mat_cmplx=complex_zero
  endif
 else
  if(ELAGd%ndiis>0) then
   ELAGd%Coef_DIIS=zero
   ELAGd%F_DIIS=zero
   ELAGd%DIIS_mat=zero
  endif
 endif

end subroutine wipeout_diis
!!***

!!***
!!****f* DoNOF/print_F_diag
!! NAME
!! print_F_diag
!!
!! FUNCTION
!!  Print the diagonal elements of the F matrix to a file (used by restart)
!!
!! INPUTS
!!  NBF_tot=Size of the F diag. array
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine print_F_diag(ELAGd,NBF_tot)
!Arguments ------------------------------------
!scalars
 class(elag_t),intent(in)::ELAGd
 integer,intent(in)::NBF_tot
!arrays
!Local variables ------------------------------
!scalars
integer::iorb,iunit=312
!arrays

!************************************************************************

 ! Print F_diag vector
 open(unit=iunit,form='unformatted',file='F_DIAG')
 do iorb=1,NBF_tot
  write(iunit) iorb,ELAGd%F_diag(iorb)
 enddo
 write(iunit) 0,zero
 close(iunit)

end subroutine print_F_diag
!!***

end module m_elag
!!***
