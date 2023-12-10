!!****m* DoNOF/m_tz_pCCD_amplitudes
!! NAME
!!  m_tz_pCCD_amplitudes
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
module m_tz_pCCD_amplitudes

 use m_nofoutput
 use m_gammatodm2 
 use m_rdmd
 use m_integd
 use m_elag
 use m_lbfgs_intern
 use m_diagf
 use m_e_grad_occ
 use m_e_grad_occ_cpx

 implicit none

 private :: calc_t_residues,calc_z_residues,num_calc_Grad_t_amp,num_calc_Grad_z_amp
!!***

 public :: calc_tz_pCCD_amplitudes
!!***

contains

!!***
!!****f* DoNOF/calc_tz_pCCD_amplitudes
!! NAME
!!  calc_tz_pCCD_amplitudes
!!
!! FUNCTION
!!
!! INPUTS
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

subroutine calc_tz_pCCD_amplitudes(ELAGd,RDMd,INTEGd,Vnn,Energy,iter_global,imethod,keep_occs) 
!Arguments ------------------------------------
!scalars
 logical,intent(in)::keep_occs
 integer,intent(in)::imethod
 integer,intent(inout)::iter_global
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
!Local variables ------------------------------
!scalars
 logical::diagco
 integer,parameter::msave=7
 integer::iter_t,iter_z,iorb,iorb1,iorb2,iorb3,iorb4,iorb5
 integer::iflag,Mtosave,Nwork
 real(dp)::tol10=1e-10
 real(dp)::sumdiff_t,sumdiff_z,maxdiff_t,maxdiff_z
 real(dp)::Ecorr_new,Ecorr_old,Ecorr_diff,Ediff
!arrays
 character(len=200)::msg
 integer,dimension(2)::info_print
 real(dp),allocatable,dimension(:)::diag,Work,diag_tz,Grad_residue
 real(dp),allocatable,dimension(:,:)::y_ij,y_ab
!************************************************************************

 Ecorr_new=zero; Ecorr_old=zero; Ecorr_diff=zero;
 maxdiff_t=zero; maxdiff_z=zero; Ediff=Energy;
 iter_t=0; iter_z=0;

 ! Build diag elements of the Lambda matrix (with HF 1-RDM and 2-RDM)
 call ELAGd%build_sd_diag(RDMd,INTEGd)
 allocate(y_ij(RDMd%Npairs,RDMd%Npairs))
 allocate(y_ab(RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs),RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)))

 ! RDMd%t_pccd(RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
 ! t-amplitudes

 ! Init. t_pccd
 if(iter_global==-1 .and. RDMd%t_amplitudes_nread) then
  do iorb=1,RDMd%Npairs ! Occ
   iorb1=iorb+RDMd%Nfrozen
   do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
    if(INTEGd%complex_ints) then
     RDMd%t_pccd_old(iorb,iorb2)=real(INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1))  &
  &  /(two*(ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1))+tol10)
    else
     RDMd%t_pccd_old(iorb,iorb2)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)  &
  &  /(two*(ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1))+tol10)
    endif
   enddo
  enddo
 endif

 if(.not.keep_occs) then
         
  if(imethod/=1) then

   iter_t=0

   do

    ! Build t_ia (using the intermediate y_ij)
    call calc_t_residues(ELAGd,RDMd,INTEGd,y_ij) 
    do iorb=1,RDMd%Npairs ! Occ
     iorb1=iorb+RDMd%Nfrozen
     do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
      iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
      RDMd%t_pccd(iorb,iorb2)=RDMd%t_pccd_old(iorb,iorb2)-RDMd%tz_residue(iorb,iorb2) &
      & /(two*(ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1))+tol10)
     enddo
    enddo

    ! Increment iter. accumulator
    iter_t=iter_t+1
   
    ! Check convergence
    Ecorr_new=zero
    sumdiff_t=zero
    maxdiff_t=-one
    do iorb=1,RDMd%Npairs ! Occ
     iorb1=iorb+RDMd%Nfrozen
     do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
      iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
      if(INTEGd%complex_ints) then
       Ecorr_new=Ecorr_new+RDMd%t_pccd(iorb,iorb2)*real(INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1))
      else        
       Ecorr_new=Ecorr_new+RDMd%t_pccd(iorb,iorb2)*INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
      endif 
      sumdiff_t=sumdiff_t+abs(RDMd%t_pccd(iorb,iorb2)-RDMd%t_pccd_old(iorb,iorb2))
      if(abs(RDMd%t_pccd(iorb,iorb2)-RDMd%t_pccd_old(iorb,iorb2))>maxdiff_t) then
       maxdiff_t=abs(RDMd%t_pccd(iorb,iorb2)-RDMd%t_pccd_old(iorb,iorb2))
      endif   
     enddo
    enddo
    Ecorr_diff=Ecorr_new-Ecorr_old
    Ecorr_old=Ecorr_new

    ! Update old t_ia
    RDMd%t_pccd_old=RDMd%t_pccd

    ! Exit if converged
    if(maxdiff_t<tol6 .and. sumdiff_t<tol5) exit

    if(abs(Ecorr_diff)<tol8) then
     write(msg,'(a,f15.6)') 't-amplitudes converged for small Ecorr dif. ',Ecorr_diff
     call write_output(msg)
     exit
    endif

    ! Exit max iter
    if(iter_t==2000) exit

   enddo

  else

   ! L-BFGS
   write(msg,'(a)') 'Calling L-BFGS to optimize t-amplitudes'
   call write_output(msg)
   Nwork=RDMd%Namplitudes*(2*msave+1)+2*msave
   Mtosave=5; info_print(1)= -1; info_print(2)= 0; diagco= .false.;
   iter_t=0; iflag=0;
   allocate(Work(Nwork),diag(RDMd%Namplitudes),diag_tz(RDMd%Namplitudes),Grad_residue(RDMd%Namplitudes))
   diag_tz=reshape(RDMd%t_pccd_old,(/RDMd%Namplitudes/))
   do
    RDMd%t_pccd_old=reshape(diag_tz,(/RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)/))
    call calc_t_residues(ELAGd,RDMd,INTEGd,y_ij)
    sumdiff_t=dsqrt(sum(RDMd%tz_residue(:,:)**two)) ! The function we are minimizing is sqrt(Sum_ia residue_ia^2)
    call num_calc_Grad_t_amp(ELAGd,RDMd,INTEGd,y_ij,Grad_residue)
    call LBFGS_INTERN(RDMd%Namplitudes,Mtosave,diag_tz,sumdiff_t,Grad_residue,diagco,diag,info_print,tol6,tol16,Work,iflag)
    if(iflag<=0) exit
    iter_t=iter_t+1
    !  We allow at most 2000 evaluations of Energy and Gradient
    if(iter_t==2000) exit
   enddo
   deallocate(Work,diag,diag_tz,Grad_residue)

   ! Update final t_ia and error
   RDMd%t_pccd=RDMd%t_pccd_old
   sumdiff_t=dsqrt(sum(RDMd%tz_residue(:,:)**two)) ! The function we are minimizing is sqrt(Sum_ia residue_ia^2)

  endif

 endif

 ! RDMd%z_pccd(RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
 ! z-amplitudes (with fixed t-amplitudes)

 ! Init. z_pccd
 if(iter_global==-1 .and. RDMd%z_amplitudes_nread) then
  RDMd%z_pccd_old=RDMd%t_pccd
 endif

 if(.not.keep_occs) then

  ! Build intermediate y_ij = sum_b v_bb^ii t_j^b 
  y_ij=zero 
  do iorb=1,RDMd%Npairs ! Occ
   iorb1=iorb+RDMd%Nfrozen
   do iorb2=1,RDMd%Npairs ! Occ
    iorb3=iorb2+RDMd%Nfrozen
    do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
     iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
     if(INTEGd%complex_ints) then
      y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd(iorb2,iorb4)*real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
     else
      y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd(iorb2,iorb4)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     endif
    enddo
   enddo
  enddo

  ! Build intermediate y_ab = sum_j v_aa^jj t_j^b 
  y_ab=zero 
  do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
   iorb1=iorb+RDMd%Nfrozen+RDMd%Npairs
   do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
    do iorb4=1,RDMd%Npairs ! Occ
     iorb5=iorb4+RDMd%Nfrozen
     if(INTEGd%complex_ints) then
      y_ab(iorb,iorb2)=y_ab(iorb,iorb2)+RDMd%t_pccd(iorb4,iorb2)*real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
     else
      y_ab(iorb,iorb2)=y_ab(iorb,iorb2)+RDMd%t_pccd(iorb4,iorb2)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     endif
    enddo
   enddo
  enddo

  if(imethod/=1) then

   iter_z=0

   do

    ! Build z_ia
    call calc_z_residues(ELAGd,RDMd,INTEGd,y_ij,y_ab) 
    do iorb=1,RDMd%Npairs ! Occ
     iorb1=iorb+RDMd%Nfrozen
     do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
      iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
      RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd_old(iorb,iorb2)-RDMd%tz_residue(iorb,iorb2) &
      & /(two*(ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1))+tol10)
     enddo
    enddo

    ! Increment iter. accumulator
    iter_z=iter_z+1
   
    ! Check convergence
    sumdiff_z=zero
    maxdiff_z=-one
    do iorb=1,RDMd%Npairs ! Occ
     iorb1=iorb+RDMd%Nfrozen
     do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
      iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
      sumdiff_z=sumdiff_z+abs(RDMd%z_pccd(iorb,iorb2)-RDMd%z_pccd_old(iorb,iorb2))
      if(abs(RDMd%z_pccd(iorb,iorb2)-RDMd%z_pccd_old(iorb,iorb2))>maxdiff_z) then
       maxdiff_z=abs(RDMd%z_pccd(iorb,iorb2)-RDMd%z_pccd_old(iorb,iorb2))
      endif   
     enddo
    enddo

    ! Update old z_ia
    RDMd%z_pccd_old=RDMd%z_pccd

    ! Exit if converged
    if(maxdiff_z<tol6 .and. sumdiff_z<tol5) exit

    ! Exit max iter
    if(iter_z==2000) exit

   enddo

  else

   ! L-BFGS
   write(msg,'(a)') 'Calling L-BFGS to optimize z-amplitudes'
   call write_output(msg)
   Nwork=RDMd%Namplitudes*(2*msave+1)+2*msave
   Mtosave=5; info_print(1)= -1; info_print(2)= 0; diagco= .false.;
   iter_z=0; iflag=0;
   allocate(Work(Nwork),diag(RDMd%Namplitudes),diag_tz(RDMd%Namplitudes),Grad_residue(RDMd%Namplitudes))
   diag_tz=reshape(RDMd%z_pccd_old,(/RDMd%Namplitudes/))
   do
    RDMd%z_pccd_old=reshape(diag_tz,(/RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)/))
    call calc_z_residues(ELAGd,RDMd,INTEGd,y_ij,y_ab)
    sumdiff_z=dsqrt(sum(RDMd%tz_residue(:,:)**two)) ! The function we are minimizing is sqrt(Sum_ia residue_ia^2)
    call num_calc_Grad_z_amp(ELAGd,RDMd,INTEGd,y_ij,y_ab,Grad_residue)
    call LBFGS_INTERN(RDMd%Namplitudes,Mtosave,diag_tz,sumdiff_z,Grad_residue,diagco,diag,info_print,tol6,tol16,Work,iflag)
    if(iflag<=0) exit
    iter_z=iter_z+1
    !  We allow at most 2000 evaluations of Energy and Gradient
    if(iter_z==2000) exit
   enddo
   deallocate(Work,diag,diag_tz,Grad_residue)

   ! Update final z_ia
   RDMd%z_pccd=RDMd%z_pccd_old
   sumdiff_z=dsqrt(sum(RDMd%tz_residue(:,:)**two)) ! The function we are minimizing is sqrt(Sum_ia residue_ia^2)

  endif

 endif

 ! Calc. the final Energy using new RDMs
 iter_global=iter_global+1
 if(INTEGd%complex_ints) then
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx)
 else
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
  & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr)
 endif
 write(msg,'(a,f15.6,a,i6,a,i6,a)') 'T-,Z-amp. opt. energy= ',Energy+Vnn,' after ',iter_t,' t-iter. and',&
  & iter_z,' z-iter.'
 call write_output(msg)
 if(iter_global>0) then
  Ediff=Energy-Ediff
  if(imethod/=1) then
   write(msg,'(a,f15.6)') 'Max. [t_pq^i+1 - t_pq^i]=      ',maxdiff_t
   call write_output(msg)
   write(msg,'(a,f15.6)') 'Max. [z_pq^i+1 - z_pq^i]=      ',maxdiff_z
   call write_output(msg)
  endif  
  write(msg,'(a,f15.6)') 'Error t-residues        =      ',sumdiff_t
  call write_output(msg)
  write(msg,'(a,f15.6)') 'Error z-residues        =      ',sumdiff_z
  call write_output(msg)
  ! This is w.r.t. |0>  (use it only for debug)
  ! This is w.r.t. |0>  (use it only for debug)
  !write(msg,'(a,f15.6)') 'Correlation Energy   =',Ecorr_new
  !call write_output(msg)
  write(msg,'(a,f19.10)') 'Energy difference amp. opt.=',Ediff
  call write_output(msg)
 endif
 write(msg,'(a,i6)') 'Number of global iter. ',iter_global
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)

 deallocate(y_ij,y_ab)

end subroutine calc_tz_pCCD_amplitudes
!!***

!!***
!!****f* DoNOF/num_calc_Grad_t_amp
!! NAME
!!  num_calc_Grad_t_amp
!!
!! FUNCTION
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

subroutine num_calc_Grad_t_amp(ELAGd,RDMd,INTEGd,y_ij,Grad_residue) 
!Arguments ------------------------------------
!scalars
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),dimension(RDMd%Namplitudes)::Grad_residue
 real(dp),dimension(RDMd%Npairs,RDMd%Npairs)::y_ij
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 real(dp)::saved_t,sum_tmp,step=tol6
!arrays
 real(dp),allocatable,dimension(:,:)::Grad_res_tmp
!************************************************************************
 
 allocate(Grad_res_tmp(RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)))
 Grad_res_tmp=zero

 do iorb=1,RDMd%Npairs
  do iorb1=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)
   saved_t=RDMd%t_pccd_old(iorb,iorb1)
   ! 2*step
   RDMd%t_pccd_old(iorb,iorb1)=RDMd%t_pccd_old(iorb,iorb1)+two*step
   call calc_t_residues(ELAGd,RDMd,INTEGd,y_ij)
   sum_tmp=-dsqrt(sum(RDMd%tz_residue(:,:)**two)) 
   ! step
   RDMd%t_pccd_old(iorb,iorb1)=RDMd%t_pccd_old(iorb,iorb1)-step
   call calc_t_residues(ELAGd,RDMd,INTEGd,y_ij)
   sum_tmp=sum_tmp+eight*dsqrt(sum(RDMd%tz_residue(:,:)**two))
   ! -step
   RDMd%t_pccd_old(iorb,iorb1)=RDMd%t_pccd_old(iorb,iorb1)-two*step
   call calc_t_residues(ELAGd,RDMd,INTEGd,y_ij)
   sum_tmp=sum_tmp-eight*dsqrt(sum(RDMd%tz_residue(:,:)**two))
   ! -2*step
   RDMd%t_pccd_old(iorb,iorb1)=RDMd%t_pccd_old(iorb,iorb1)-step
   call calc_t_residues(ELAGd,RDMd,INTEGd,y_ij)
   sum_tmp=sum_tmp+dsqrt(sum(RDMd%tz_residue(:,:)**two)) 
   ! Save the gradient
   Grad_res_tmp(iorb,iorb1)=sum_tmp/(twelve*step) 
   ! Recover t_pccd_old
   RDMd%t_pccd_old(iorb,iorb1)=saved_t
  enddo
 enddo
 
 Grad_residue=reshape(Grad_res_tmp,(/RDMd%Namplitudes/))

 deallocate(Grad_res_tmp)

end subroutine num_calc_Grad_t_amp
!!***

!!***
!!****f* DoNOF/num_calc_Grad_z_amp
!! NAME
!!  num_calc_Grad_z_amp
!!
!! FUNCTION
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

subroutine num_calc_Grad_z_amp(ELAGd,RDMd,INTEGd,y_ij,y_ab,Grad_residue) 
!Arguments ------------------------------------
!scalars
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),dimension(RDMd%Namplitudes)::Grad_residue
 real(dp),dimension(RDMd%Npairs,RDMd%Npairs)::y_ij
 real(dp),dimension(RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs),RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))::y_ab
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 real(dp)::saved_z,sum_tmp,step=tol6
!arrays
 real(dp),allocatable,dimension(:,:)::Grad_res_tmp
!************************************************************************
 
 allocate(Grad_res_tmp(RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)))
 Grad_res_tmp=zero

 do iorb=1,RDMd%Npairs
  do iorb1=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)
   saved_z=RDMd%z_pccd_old(iorb,iorb1)
   ! 2*step
   RDMd%z_pccd_old(iorb,iorb1)=RDMd%z_pccd_old(iorb,iorb1)+two*step
   call calc_z_residues(ELAGd,RDMd,INTEGd,y_ij,y_ab)
   sum_tmp=-dsqrt(sum(RDMd%tz_residue(:,:)**two))
   ! step
   RDMd%z_pccd_old(iorb,iorb1)=RDMd%z_pccd_old(iorb,iorb1)-step
   call calc_z_residues(ELAGd,RDMd,INTEGd,y_ij,y_ab)
   sum_tmp=sum_tmp+eight*dsqrt(sum(RDMd%tz_residue(:,:)**two)) 
   ! -step
   RDMd%z_pccd_old(iorb,iorb1)=RDMd%z_pccd_old(iorb,iorb1)-two*step
   call calc_z_residues(ELAGd,RDMd,INTEGd,y_ij,y_ab)
   sum_tmp=sum_tmp-eight*dsqrt(sum(RDMd%tz_residue(:,:)**two))
   ! -2*step
   RDMd%z_pccd_old(iorb,iorb1)=RDMd%z_pccd_old(iorb,iorb1)-step
   call calc_z_residues(ELAGd,RDMd,INTEGd,y_ij,y_ab)
   sum_tmp=sum_tmp+dsqrt(sum(RDMd%tz_residue(:,:)**two))
   ! Save the gradient
   Grad_res_tmp(iorb,iorb1)=sum_tmp/(twelve*step) 
   ! Recover z_pccd_old
   RDMd%z_pccd_old(iorb,iorb1)=saved_z
  enddo
 enddo
 
 Grad_residue=reshape(Grad_res_tmp,(/RDMd%Namplitudes/))

 deallocate(Grad_res_tmp)

end subroutine num_calc_Grad_z_amp
!!***

!!***
!!****f* DoNOF/calc_t_residues
!! NAME
!!  calc_t_residues
!!
!! FUNCTION
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

subroutine calc_t_residues(ELAGd,RDMd,INTEGd,y_ij) 
!Arguments ------------------------------------
!scalars
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),dimension(RDMd%Npairs,RDMd%Npairs)::y_ij
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
 real(dp)::sum_tmp
!arrays
!************************************************************************

 ! Build intermediate y_ij = sum_b v_bb^jj t_i^b 
 y_ij=zero 
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  do iorb2=1,RDMd%Npairs ! Occ
   iorb3=iorb2+RDMd%Nfrozen
   do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
    if(INTEGd%complex_ints) then
     y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)*real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
    else
     y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
    endif 
   enddo
  enddo
 enddo

 ! Build the t residue
 RDMd%tz_residue=zero
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
   iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
   if(INTEGd%complex_ints) then
    RDMd%tz_residue(iorb,iorb2)=real(INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1))
   else
    RDMd%tz_residue(iorb,iorb2)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
   endif
   sum_tmp=ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1)
   do iorb4=1,RDMd%Npairs ! Occ
    iorb5=iorb4+RDMd%Nfrozen
    if(INTEGd%complex_ints) then
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb4,iorb2)* &
     &        real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb4,iorb2)*real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
    else       
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb4,iorb2)* &
     &        INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb4,iorb2)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
    endif 
    RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+y_ij(iorb,iorb4)*RDMd%t_pccd_old(iorb4,iorb2)
   enddo
   do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
    if(INTEGd%complex_ints) then
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)* &
     &        real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb,iorb4)*real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
    else
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)* &
     &        INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
    endif
   enddo
   if(INTEGd%complex_ints) then
    sum_tmp=sum_tmp-real(two*INTEGd%ERImol_cmplx(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1) &
   &         -INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd_old(iorb,iorb2))
   else
    sum_tmp=sum_tmp-(two*INTEGd%ERImol(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1) &
   &         -INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd_old(iorb,iorb2))
   endif
   sum_tmp=two*sum_tmp*RDMd%t_pccd_old(iorb,iorb2)
   RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+sum_tmp
  enddo
 enddo

end subroutine calc_t_residues
!!***

!!***
!!****f* DoNOF/calc_z_residues
!! NAME
!!  calc_z_residues
!!
!! FUNCTION
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

subroutine calc_z_residues(ELAGd,RDMd,INTEGd,y_ij,y_ab) 
!Arguments ------------------------------------
!scalars
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),dimension(RDMd%Npairs,RDMd%Npairs)::y_ij
 real(dp),dimension(RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs),RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))::y_ab
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
 real(dp)::sum_tmp,zt_ii,zt_aa
!arrays
!************************************************************************

 ! Build the z residue
 RDMd%tz_residue=zero
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
   iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
   if(INTEGd%complex_ints) then
    RDMd%tz_residue(iorb,iorb2)=real(INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1))
   else        
    RDMd%tz_residue(iorb,iorb2)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
   endif 
   sum_tmp=ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1)
   zt_aa=zero; zt_ii=zero;
   do iorb4=1,RDMd%Npairs ! Occ
    iorb5=iorb4+RDMd%Nfrozen
    if(INTEGd%complex_ints) then
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%z_pccd_old(iorb4,iorb2)* & 
     &        real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
     sum_tmp=sum_tmp-RDMd%t_pccd(iorb4,iorb2)*real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
    else
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%z_pccd_old(iorb4,iorb2)* &
     &        INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     sum_tmp=sum_tmp-RDMd%t_pccd(iorb4,iorb2)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
    endif 
    RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+y_ij(iorb,iorb4)*RDMd%z_pccd_old(iorb4,iorb2)
    zt_aa=zt_aa+RDMd%t_pccd(iorb4,iorb2)*RDMd%z_pccd_old(iorb4,iorb2)
   enddo
   do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
    if(INTEGd%complex_ints) then
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%z_pccd_old(iorb,iorb4)* &
     &   real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
     sum_tmp=sum_tmp-RDMd%t_pccd(iorb,iorb4)*real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
    else
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%z_pccd_old(iorb,iorb4)* &
     &       INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
     sum_tmp=sum_tmp-RDMd%t_pccd(iorb,iorb4)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
    endif 
    RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+y_ab(iorb,iorb4)*RDMd%z_pccd_old(iorb,iorb4)
    zt_ii=zt_ii+RDMd%t_pccd(iorb,iorb4)*RDMd%z_pccd_old(iorb,iorb4)
   enddo
   if(INTEGd%complex_ints) then
    sum_tmp=sum_tmp-real(two*INTEGd%ERImol_cmplx(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1) &
    &       -two*INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd(iorb,iorb2))
    RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)-two*(zt_aa+zt_ii)* & 
    &        real(INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1))
   else
    sum_tmp=sum_tmp-(two*INTEGd%ERImol(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1) &
    &       -two*INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd(iorb,iorb2))
    RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)-two*(zt_aa+zt_ii)* &
    &        INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
   endif
   sum_tmp=two*sum_tmp*RDMd%z_pccd_old(iorb,iorb2)
   RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+sum_tmp
  enddo
 enddo

end subroutine calc_z_residues
!!***

end module m_tz_pCCD_amplitudes
!!***
