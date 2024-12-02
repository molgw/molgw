!!****m* DoNOF/m_optocc
!! NAME
!!  m_optocc
!!
!! FUNCTION
!!  Module prepared to compute occ optimization for a fixed hCORE and ERIs
!!
!! PARENTS
!!  m_noft_driver
!!
!! CHILDREN
!!  m_e_grad_occ
!!  m_lbfgs
!!
!! SOURCE
module m_optocc

 use m_nofoutput
 use m_rdmd
 use m_lbfgs_intern
 use m_e_grad_occ
 use m_e_grad_occ_cpx

 implicit none

!!private :: 
!!***

 public :: opt_occ,occ_chempot 
!!***

contains

!!***
!!****f* DoNOF/opt_occ
!! NAME
!!  opt_occ
!!
!! FUNCTION
!!  Call the L-BFGS subroutine for occ optimization 
!!
!! INPUTS
!!  imethocc=Method to use for the occ. optimization (default = 0, i.e. conjugate gradients)
!!  hCORE=One-body integrals (h_pq) 
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  ERI_L=Lower triangular part of the L_pq matrix
!!  Vnn=Nuclear-nuclear rep. energy
!!
!! OUTPUT
!!  Energy=Sum of nuclear and electronic energy
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine opt_occ(iter,imethod,keep_occs,RDMd,Vnn,Energy,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr,&
& hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::keep_occs
 integer,intent(inout)::iter
 integer,intent(in)::imethod
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
 real(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L
 real(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr,ERI_Lsr
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE_cmplx
 complex(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx 
 complex(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr_cmplx,ERI_Lsr_cmplx
!Local variables ------------------------------
!scalars
 logical::diagco,conveg=.false.,debug=.false.,cpx_mos=.false.
 integer,parameter::msave=7
 integer::iorb,igamma,iflag,icall,Mtosave,Nwork
!arrays
 character(len=200)::msg
 integer,dimension(2)::info_print
 real(dp),allocatable,dimension(:)::GAMMAs,Grad_GAMMAs,diag,Work
!************************************************************************

 if(present(hCORE_cmplx).and.present(ERI_J_cmplx).and.present(ERI_K_cmplx).and.present(ERI_L_cmplx)) then
  cpx_mos=.true.
 endif
 Energy=zero
 allocate(GAMMAs(RDMd%Ngammas),Grad_GAMMAs(RDMd%Ngammas))
 Grad_GAMMAs=zero
 if((iter==-1).and.RDMd%GAMMAs_nread) then 
  GAMMAs=pi/four               ! Perturbed occ. numbers (i.e pi/4) -> occ(i<Fermi level) = 0.75
  if(RDMd%INOF==0) GAMMAs=zero ! HF set occ. to 0 or 1. open-shell TODO
 else
  if(RDMd%INOF==0) then        ! HF set occ. to 0 or 1. open-shell TODO
   GAMMAs=zero
  else
   GAMMAs=RDMd%GAMMAs_old  ! Read from previous run
  endif
 endif

 if(RDMd%INOF<0) then
  write(msg,'(a)') 'Error: for pCCD we should not enter m_optocc module.'
  call write_output(msg)
  error stop
 endif

 ! Check if the current GAMMAs already solve the problem. Is it converged? 
 if(cpx_mos) then
  call calc_E_occ_cmplx(RDMd,GAMMAs,Energy,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,&
  &    ERI_Jsr_cmplx,ERI_Lsr_cmplx)
  call calc_Grad_occ_cmplx(RDMd,Grad_GAMMAs,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,&
  &    ERI_Jsr_cmplx,ERI_Lsr_cmplx)
  !call num_calc_Grad_occ_cmplx(RDMd,GAMMAs,Grad_GAMMAs,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx)
 else
  call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
  call calc_Grad_occ(RDMd,Grad_GAMMAs,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
  !call num_calc_Grad_occ(RDMd,GAMMAs,Grad_GAMMAs,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
 endif
 conveg=.true.
 do igamma=1,RDMd%Ngammas
  if(dabs(Grad_GAMMAs(igamma))>tol6.and.conveg) then
   conveg=.false.
  endif 
  Grad_GAMMAs(igamma)=zero
 enddo

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 ! Do iterations if the current GAMMAs do not produce small gradients
 icall=0
 if((.not.conveg).and.(.not.keep_occs)) then 
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
  if(imethod==1) then ! L-BFGS
   write(msg,'(a)') 'Calling L-BFGS to optimize occ. numbers'
   call write_output(msg)
   Nwork=RDMd%Ngammas*(2*msave+1)+2*msave
   Mtosave=5; info_print(1)= -1; info_print(2)= 0; diagco= .false.;
   icall=0; iflag=0;
   allocate(Work(Nwork),diag(RDMd%Ngammas))
   do
    if(cpx_mos) then
     call calc_E_occ_cmplx(RDMd,GAMMAs,Energy,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,&
     &    ERI_Jsr_cmplx,ERI_Lsr_cmplx)
     call calc_Grad_occ_cmplx(RDMd,Grad_GAMMAs,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,&
     &    ERI_Jsr_cmplx,ERI_Lsr_cmplx)
     !call num_calc_Grad_occ_cmplx(RDMd,GAMMAs,Grad_GAMMAs,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx)
    else
     call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
     call calc_Grad_occ(RDMd,Grad_GAMMAs,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
     !call num_calc_Grad_occ(RDMd,GAMMAs,Grad_GAMMAs,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
    endif
    call LBFGS_INTERN(RDMd%Ngammas,Mtosave,GAMMAs,Energy,Grad_GAMMAs,diagco,diag,info_print,tol5,tol16,Work,iflag)
    if(iflag<=0) exit
    icall=icall+1
    !  We allow at most 2000 evaluations of Energy and Gradient
    if(icall==2000) exit
   enddo
   deallocate(Work,diag)
  endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 endif 
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 
 iter=iter+1
 if(RDMd%INOF==0) then ! Ensure that for HF we keep occ. 0 or 1
  RDMd%GAMMAs_old=zero
 else
  RDMd%GAMMAs_old=GAMMAs
 endif
 if(cpx_mos) then
  call calc_E_occ_cmplx(RDMd,GAMMAs,Energy,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,&
  &    ERI_Jsr_cmplx,ERI_Lsr_cmplx)
 else
  call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
 endif
 write(msg,'(a,f15.6,a,i6,a)') 'Occ. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
 call write_output(msg)
 Grad_GAMMAs(:)=dabs(Grad_GAMMAs(:))
 write(msg,'(a,f15.6)') 'Max. [|Grad Energy w.r.t. GAMMAS|]= ',maxval(Grad_GAMMAs(:))
 call write_output(msg)
 write(msg,'(a)') 'Optimized occ. numbers '
 call write_output(msg)
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(msg,'(f12.6,9f11.6)') two*RDMd%occ(iorb:iorb+9)
  call write_output(msg)
 enddo
 iorb=(RDMd%NBF_occ/10)*10+1
 write(msg,'(f12.6,*(f11.6))') two*RDMd%occ(iorb:)
 call write_output(msg)
 write(msg,'(a,i6)') 'Number of global iter. ',iter
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 
 if(icall==2000) then
  write(msg,'(a)') 'Warning! Max. number of iterations (2000) reached in occ. optimization'
  call write_output(msg)
 endif

 deallocate(GAMMAs,Grad_GAMMAs)

end subroutine opt_occ
!!***

!!***
!!****f* DoNOF/occ_chempot
!! NAME
!!  occ_chempot
!!
!! FUNCTION
!!  Compute the 2RDM from GAMMAs and the gradients.
!!  Then, compute the chem. pot. for each orbital. 
!!
!! INPUTS
!!  hCORE=One-body integrals (h_pq) 
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  ERI_L=Lower triangular part of the L_pq matrix
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine occ_chempot(RDMd,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr,hCORE_cmplx,&
 &     ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
 real(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L
 real(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr,ERI_Lsr
 complex(dp),optional,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE_cmplx
 complex(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx 
 complex(dp),optional,dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr_cmplx,ERI_Lsr_cmplx 
!Local variables ------------------------------
!scalars
 logical::chempot,cpx_mos=.false.
 real(dp)::Energy
!arrays
 real(dp),allocatable,dimension(:)::GAMMAs,Grad_GAMMAs
!************************************************************************

 if(present(hCORE_cmplx).and.present(ERI_J_cmplx).and.present(ERI_K_cmplx).and.present(ERI_L_cmplx)) then
  cpx_mos=.true.
 endif
 Energy=zero
 allocate(GAMMAs(RDMd%Ngammas),Grad_GAMMAs(RDMd%Ngammas))
 GAMMAs=RDMd%GAMMAs_old  ! Read from previous run

 ! Calc. the 2RDM and derivatives in RDMd
 if(cpx_mos) then
  call calc_E_occ_cmplx(RDMd,GAMMAs,Energy,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,&
  &    ERI_Jsr_cmplx,ERI_Lsr_cmplx,chempot=chempot)
  call calc_Chem_pot_cmplx(RDMd,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx)
 else
  ! MRM: Warning! In rs-NOFT, the chemical potential could be computed without the d Exc / dn contribution...
  call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr,chempot=chempot)
  call calc_Chem_pot(RDMd,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,ERI_Lsr)
 endif

 deallocate(GAMMAs,Grad_GAMMAs)

end subroutine occ_chempot
!!***

end module m_optocc 
!!***
