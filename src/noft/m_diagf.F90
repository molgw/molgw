!!****m* DoNOF/m_diagf
!! NAME
!!  m_diagf
!!
!! FUNCTION
!!  Module prepared to compute orb optimization using the F diag method as in DoNOF
!!
!! PARENTS
!!  m_optorb
!!
!! CHILDREN
!!
!! SOURCE
module m_diagf

 use m_rdmd
 use m_elag

 implicit none

 private :: scale_F,diis_F,traceF
!!***

 public :: diagF_to_coef
!!***

contains

!!***
!!****f* DoNOF/diagF_to_coef
!! NAME
!!  diagF_to_coef
!!
!! FUNCTION
!!  Build the F-matrix, diagonalize it and update the NO_COEF
!!
!! INPUTS
!!  iter=Number of global iteration 
!!  icall=Number of call from opt_orb subroutine
!!
!! OUTPUT
!!  NO_COEF=Updated Nat. orb. coefs.
!!  ELAGd%F_diag=Update the diag elements of the F_pq matrix 
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine diagF_to_coef(iter,icall,maxdiff,diddiis,ELAGd,RDMd,NO_COEF) 
!Arguments ------------------------------------
!scalars
 logical,intent(inout)::diddiis
 integer,intent(in)::iter
 integer,intent(inout)::icall
 real(dp),intent(in)::maxdiff
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,lwork,info
 real(dp)::thresholddiis
!arrays
 real(dp),allocatable,dimension(:)::Work
 real(dp),allocatable,dimension(:,:)::Eigvec,New_NO_COEF ! Eigvec is initailly the F matrix
!************************************************************************
 
 thresholddiis=1.0d1**(-ELAGd%itoldiis)
 allocate(New_NO_COEF(RDMd%NBF_tot,RDMd%NBF_tot),Eigvec(RDMd%NBF_tot,RDMd%NBF_tot),Work(1))

 if((icall==0.and.iter==0).and.(ELAGd%diagLpL.and.(.not.ELAGd%diagLpL_done))) then
  ELAGd%diagLpL_done=.true. 
  do iorb=1,RDMd%NBF_tot 
   do iorb1=1,iorb-1
    Eigvec(iorb,iorb1)=0.5d0*(ELAGd%Lambdas(iorb,iorb1)+ELAGd%Lambdas(iorb1,iorb))
    Eigvec(iorb1,iorb)=Eigvec(iorb,iorb1)
   enddo
   Eigvec(iorb,iorb)=ELAGd%Lambdas(iorb,iorb)
  enddo
 else
  do iorb=1,RDMd%NBF_tot 
   do iorb1=1,iorb-1
    Eigvec(iorb,iorb1)=ELAGd%Lambdas(iorb1,iorb)-ELAGd%Lambdas(iorb,iorb1)
    call scale_F(ELAGd%MaxScaling+9,Eigvec(iorb,iorb1)) ! Scale the Fpq element to avoid divergence
    Eigvec(iorb1,iorb)=Eigvec(iorb,iorb1)               ! Fpq=Fqp
   enddo
   Eigvec(iorb,iorb)=ELAGd%F_diag(iorb)
  enddo  
 endif

 ! Decide whether to do DIIS before diagonalizing
 if(maxdiff<thresholddiis.and.ELAGd%ndiis>0) then
  call diis_F(diddiis,RDMd,ELAGd,Eigvec)
 endif 

 ! Prepare F_pq diagonalization (stored as Eigvec) and diagonalize it to produce the rot. matrix
 lwork=-1
 call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,ELAGd%F_diag,Work,lwork,info)
 lwork=nint(Work(1))
 if(info==0) then
  deallocate(Work)
  allocate(Work(lwork))
  ELAGd%F_diag=0.0d0
  call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,ELAGd%F_diag,Work,lwork,info)
 endif

 ! Update the NO_COEF
 New_NO_COEF=matmul(NO_COEF,Eigvec)
 NO_COEF=New_NO_COEF

 ! Increase icall (iterator accumulator)
 icall=icall+1

 deallocate(New_NO_COEF,Eigvec,Work)

end subroutine diagF_to_coef
!!***

!!***
!!****f* DoNOF/scale_F
!! NAME
!!  scale_F
!!
!! FUNCTION
!!  Scale the F-matrix element
!!
!! INPUTS
!!  MaxScaling=Number of zeros used to scale F
!!
!! OUTPUT
!!  Fpq=Scaled F_pq matrix element
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine scale_F(MaxScaling,Fpq) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::MaxScaling
 real(dp),intent(inout)::Fpq
!arrays
!Local variables ------------------------------
!scalars
 integer::iscale
!arrays
 real(dp)::Abs_Fpq
!************************************************************************
 do iscale=1,MaxScaling
  Abs_Fpq=dabs(Fpq)
  if(Abs_Fpq>1.0d1**(9-iscale).and.Abs_Fpq<1.0d1**(10-iscale)) then
   Fpq=0.1d0*Fpq
  endif
 enddo 
end subroutine scale_F
!!***

!!***
!!****f* DoNOF/diis_F
!! NAME
!!  diis_F
!!
!! FUNCTION
!!  Use DIIS on the F-matrix before diagonalizing it
!!
!! INPUTS
!!  RDMd=Density matrix descriptor
!!  ELAGd=Langragian descriptor
!!
!! OUTPUT
!!  Eigvec=F matrix on input, DIIS F matrix on output
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine diis_F(diddiis,RDMd,ELAGd,Eigvec) 
!Arguments ------------------------------------
!scalars
 logical,intent(inout)::diddiis
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::Eigvec
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,idiis1,idiisp1,info
!arrays
 integer,allocatable,dimension(:)::IPIV
!************************************************************************
 ELAGd%idiis=ELAGd%idiis+1 
 ELAGd%F_DIIS(ELAGd%idiis,:,:)=Eigvec(:,:)
 idiisp1=ELAGd%idiis+1
 do idiis1=1,ELAGd%idiis
  ELAGd%DIIS_mat(idiis1,ELAGd%idiis) = traceF(RDMd,ELAGd,idiis1)
  ELAGd%DIIS_mat(ELAGd%idiis,idiis1) = ELAGd%DIIS_mat(idiis1,ELAGd%idiis)
  ELAGd%DIIS_mat(idiis1,idiisp1) = -1.0d0
  ELAGd%DIIS_mat(idiisp1,idiis1) = -1.0d0
 enddo
 ELAGd%DIIS_mat(idiisp1,idiisp1) = 0.0d0
 if(ELAGd%idiis>ELAGd%ndiis) then
  diddiis=.true.
  allocate(IPIV(ELAGd%ndiis_array))
  IPIV=0
  ELAGd%Coef_DIIS=0.0d0
  ELAGd%Coef_DIIS(ELAGd%ndiis_array) = -1.0d0
  call DGESV(ELAGd%ndiis_array,1,ELAGd%DIIS_mat,ELAGd%ndiis_array,IPIV,ELAGd%Coef_DIIS,ELAGd%ndiis_array,info)
  deallocate(IPIV)
  do iorb=1,RDMd%NBF_tot
   do iorb1=1,iorb-1
    Eigvec(iorb,iorb1)=sum(ELAGd%Coef_DIIS(1:ELAGd%ndiis_array-1)*ELAGd%F_DIIS(1:ELAGd%ndiis_array-1,iorb,iorb1))
    Eigvec(iorb1,iorb)=Eigvec(iorb,iorb1)
   enddo
  enddo
  call ELAGd%clean_diis()
 endif

end subroutine diis_F
!!***

!!***
!!****f* DoNOF/traceF
!! NAME
!!  traceF
!!
!! FUNCTION
!!  Multiply F matrices to produce DIIS matrix elements.
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

function traceF(RDMd,ELAGd,idiis_in) result(traceFF)
!Arguments ------------------------------------
!scalars
 integer::idiis_in
 real(dp)::traceFF
 type(elag_t),intent(in)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
!arrays
!************************************************************************
traceFF = 0.0d0
 do iorb=1,RDMd%NBF_tot
  do iorb1=1,iorb-1
   traceFF=traceFF+ELAGd%F_DIIS(idiis_in,iorb,iorb1)*ELAGd%F_DIIS(ELAGd%idiis,iorb1,iorb)
  enddo
 enddo 
end function traceF
!!***

end module m_diagf
!!***
