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
!!  m_cg
!!  m_lbfgs
!!
!! SOURCE
module m_optocc

 use m_nofoutput
 use m_rdmd
 use m_cg
 use m_lbfgs_intern
 use m_e_grad_occ

 implicit none

!!private :: 
!!***

 public :: opt_occ 
!!***

contains

!!***
!!****f* DoNOF/opt_occ
!! NAME
!!  opt_occ
!!
!! FUNCTION
!!  Call the CG of LBFGS subroutines for occ optimization 
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

subroutine opt_occ(iter,imethod,RDMd,Vnn,Energy,hCORE,ERI_J,ERI_K,ERI_L)
!Arguments ------------------------------------
!scalars
 integer,intent(inout)::iter
 integer,intent(in)::imethod
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L 
!Local variables ------------------------------
!scalars
 logical::diagco,conveg=.false.
 integer,parameter::msave=7,nextv=47,nfcall=6,nfgcal=7,g=28,toobig=2,vneed=4
 integer::igamma,iflag,ig,icall,icall1,Mtosave,Nwork,Nwork2
 real(dp)::eps,xtol,tolgamma=1.0d-6
!arrays
 character(len=200)::msg
 integer,dimension(2)::info_print
 integer,allocatable,dimension(:)::iWork
 real(dp),allocatable,dimension(:)::GAMMAs,Grad_GAMMAs,diag,Work,Work2
!************************************************************************

 Energy=0.0d0
 allocate(GAMMAs(RDMd%Ngammas),Grad_GAMMAs(RDMd%Ngammas))
 Grad_GAMMAs=0.0d0
 if((iter==-1).and.RDMd%GAMMAs_nread) then 
  GAMMAs=0.785398163       ! Perturbed occ. numbers (i.e pi/4) -> occ(i<Fermi level) = 0.75
 else
  GAMMAs=RDMd%GAMMAs_old   ! Read from previous run
 endif

 ! Check if the current GAMMAs already solve the problem. Is it converged? 
 call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L)
 call calc_Grad_occ(RDMd,Grad_GAMMAs,hCORE,ERI_J,ERI_K,ERI_L)
 conveg=.true.
 do igamma=1,RDMd%Ngammas
  if(dabs(Grad_GAMMAs(igamma))>tolgamma.and.conveg) then
   conveg=.false.
  endif 
  Grad_GAMMAs(igamma)=0.0d0
 enddo

 ! Do iterations if the current GAMMAs do not produce small gradients
 icall=0
 if(.not.conveg) then 
  if(imethod==1) then ! Conjugate gradients. (The subroutine uses goto. It is not clean but needed)
   write(msg,'(a)') 'Calling CG to optimize occ. numbers'
   call write_output(msg)
   Nwork=60; Nwork2=71+RDMd%Ngammas*(RDMd%Ngammas+15)/2; 
   allocate(iWork(Nwork),Work(RDMd%Ngammas),Work2(Nwork2))
   iWork=0; Work=0.1d0;   
   if (iWork(1)==0) call deflt(2,iWork, Nwork, Nwork2, Work2)
   iflag = iWork(1)
   if (iflag == 12 .or. iflag == 13) iWork(vneed) = iWork(vneed) + RDMd%Ngammas
   if (iflag == 14) goto 10
   if (iflag > 2 .and. iflag < 12) goto 10
   ig = 1
   if (iflag == 12) iWork(1) = 13
   goto 20

10  ig = iWork(g)

20  call sumit(Work, Energy, Work2(ig), iWork, Nwork, Nwork2, RDMd%Ngammas, Work2, GAMMAs)
   if(iWork(1)-2<0) then 
    goto 30
   elseif(iWork(1)-2==0) then
    goto 40
   else
    goto 50
   endif

30  icall1 = iWork(nfcall)
   call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L)
   icall=icall+1
   if(icall==2000) goto 60
   if(icall1<=0) iWork(toobig) = 1
   goto 20

40  call calc_Grad_occ(RDMd,Grad_GAMMAs,hCORE,ERI_J,ERI_K,ERI_L)
   Work2(ig:ig+RDMd%Ngammas)=Grad_GAMMAs(1:RDMd%Ngammas)
   goto 20

50  if(iWork(1) /= 14) then
       goto 60
    endif
!
!  Storage allocation
!
   iWork(g) = iWork(nextv)
   iWork(nextv) = iWork(g) + RDMd%Ngammas
   if(iflag /= 13) goto 10

60  deallocate(iWork,Work,Work2)

  else ! LBFGS
   write(msg,'(a)') 'Calling LBFGS to optimize occ. numbers'
   call write_output(msg)
   Nwork=RDMd%Ngammas*(2*msave+1)+2*msave
   Mtosave=5; info_print(1)= -1; info_print(2)= 0; diagco= .false.;
   eps= 1.0d-5; xtol= 1.0d-16; icall=0; iflag=0;
   allocate(Work(Nwork),diag(RDMd%Ngammas))
   do
    call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L)
    call calc_Grad_occ(RDMd,Grad_GAMMAs,hCORE,ERI_J,ERI_K,ERI_L)
    call LBFGS_INTERN(RDMd%Ngammas,Mtosave,GAMMAs,Energy,Grad_GAMMAs,diagco,diag,info_print,eps,xtol,Work,iflag)
    if(iflag<=0) exit
    icall=icall+1
!  We allow at most 2000 evaluations of Energy and Gradient
    if(icall==2000) exit
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
   enddo
   deallocate(Work,diag)
  endif
 endif 
 
 iter=iter+1
 RDMd%GAMMAs_old=GAMMAs
 call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L)
 write(msg,'(a,f15.6,a,i6,a)') 'Occ. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
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

end module m_optocc 
!!***
