!!****m* DoNOF/m_anti2unit
!! NAME
!!  m_anti2unit
!!
!! FUNCTION
!! Compute unitary matrix from X (real or complex) anti-sym/herm. matrix
!!
!! PARENTS
!!  m_noft_driver
!!
!! CHILDREN
!!
!! SOURCE
module m_anti2unit
 
 use m_definitions

 implicit none

!! private ::
!!***

 public :: anti_2_unitary 
!!***

contains

!!***
!!****f* DoNOF/anti_2_unitary
!! NAME
!!  anti_2_unitary
!!
!! FUNCTION
!!  Compute form X_complex anti-Hermitian the (complex) Unitary U_mat_complx matrix
!!  Compute form X anti-symmetric the (real) Unitary U_mat matrix
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

subroutine anti_2_unitary(NBF_tot,X_mat,U_mat,X_mat_cmplx,U_mat_cmplx) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_tot
!arrays
 real(dp),optional,dimension(NBF_tot,NBF_tot),intent(inout)::U_mat,X_mat
 complex(dp),optional,dimension(NBF_tot,NBF_tot),intent(inout)::U_mat_cmplx,X_mat_cmplx
!Local variables ------------------------------
!scalars
 integer::iorb,lwork,info
!arrays
 real(dp),allocatable,dimension(:)::delta_diag,RWork
 complex(dp),allocatable,dimension(:)::Work_cmplx
 complex(dp),allocatable,dimension(:,:)::Eigvec_cmplx,Tmp_mat,Tmp_mat2
!************************************************************************
!************************************************************************

 allocate(delta_diag(NBF_tot),Work_cmplx(1),RWork(3*NBF_tot-2))
 allocate(Tmp_mat(NBF_tot,NBF_tot),Tmp_mat2(NBF_tot,NBF_tot),Eigvec_cmplx(NBF_tot,NBF_tot))

 ! Build iX (Hermitian)
 if(present(X_mat_cmplx).and.present(U_mat_cmplx)) then
  Eigvec_cmplx(:,:)=im*X_mat_cmplx(:,:)
  U_mat_cmplx=complex_zero 
 else
  Eigvec_cmplx(:,:)=im*X_mat(:,:) 
  U_mat=zero
 endif
  
 ! Diag iX ->  V delta V^dagger
 lwork=-1
 call ZHEEV('V','L',NBF_tot,Eigvec_cmplx,NBF_tot,delta_diag,Work_cmplx,lwork,RWork,info)
 lwork=nint(real(Work_cmplx(1)))
 if(info==0) then
  deallocate(Work_cmplx)
  allocate(Work_cmplx(lwork))
  delta_diag=zero
  call ZHEEV('V','L',NBF_tot,Eigvec_cmplx,NBF_tot,delta_diag,Work_cmplx,lwork,RWork,info)
 endif

 ! iX = V delta V^dagger --times i--> -X = i V delta V^dagger -> X = V (-i delta) V^dagger
 ! exp(X) = exp(V (-i delta) V^dagger) = V exp(-i delta) V^dagger
 Tmp_mat=complex_zero
 do iorb=1,NBF_tot
  Tmp_mat(iorb,iorb)=exp(-im*delta_diag(iorb))
 enddo
 Tmp_mat2=matmul(matmul(Eigvec_cmplx,Tmp_mat),conjg(transpose(Eigvec_cmplx)))

 ! Set U = V exp(-i delta) V^dagger 
 if(present(X_mat_cmplx).and.present(U_mat_cmplx)) then
  U_mat_cmplx(:,:)=Tmp_mat2(:,:)
 else
  U_mat(:,:)=real(Tmp_mat2(:,:))
 endif

 deallocate(Work_cmplx,RWork,delta_diag,Tmp_mat,Tmp_mat2,Eigvec_cmplx)

end subroutine anti_2_unitary
!!***

end module m_anti2unit
!!***
