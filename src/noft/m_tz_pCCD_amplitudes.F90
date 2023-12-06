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
 use m_diagf
 use m_e_grad_occ
 use m_e_grad_occ_cpx

 implicit none

!! private :: lambda_conv
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

subroutine calc_tz_pCCD_amplitudes(ELAGd,RDMd,INTEGd,Vnn,Energy,iter_global) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::iter_global
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
!Local variables ------------------------------
!scalars
 logical::dodiis
 integer::iter_t,iter_z,iorb,iorb1,iorb2,iorb3,iorb4,iorb5
 real(dp)::tol10=1e-10
 real(dp)::sumdiff_t,sumdiff_z,maxdiff_t,maxdiff_z,sum_val,zt_aa,zt_ii
 real(dp)::Ecorr_new,Ecorr_old,Ecorr_diff,Ediff
!arrays
 character(len=200)::msg
 real(dp),allocatable,dimension(:,:)::y_ij,y_ab
!************************************************************************

 Ecorr_new=zero; Ecorr_old=zero; Ecorr_diff=zero;
 maxdiff_t=zero; maxdiff_z=zero; Ediff=Energy;

 ! Build diag elements of the Lambda matrix (with HF 1-RDM and 2-RDM)
 call ELAGd%build_sd_diag(RDMd,INTEGd)
 allocate(y_ij(RDMd%Npairs,RDMd%Npairs))
 allocate(y_ab(RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs),RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)))

 ! RDMd%t_pccd(RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
 ! t-amplitudes

 ! Init. t_pccd
 if(iter_global==-1) then
  do iorb=1,RDMd%Npairs ! Occ
   iorb1=iorb+RDMd%Nfrozen
   do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
    RDMd%t_pccd_old(iorb,iorb2)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)  &
  & /(two*(ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1))+tol10)
   enddo
  enddo
 endif

 iter_t=0
 do

  ! Build intermediate y_ij = sum_b v_bb^jj t_i^b 
  y_ij=zero 
  do iorb=1,RDMd%Npairs ! Occ
   iorb1=iorb+RDMd%Nfrozen
   do iorb2=1,RDMd%Npairs ! Occ
    iorb3=iorb2+RDMd%Nfrozen
    do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
     iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
     y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
    enddo
   enddo
  enddo

  ! Build t_ia
  RDMd%t_pccd=zero
  do iorb=1,RDMd%Npairs ! Occ
   iorb1=iorb+RDMd%Nfrozen
   do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
    RDMd%t_pccd(iorb,iorb2)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
    sum_val=ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1)
    do iorb4=1,RDMd%Npairs ! Occ
     iorb5=iorb4+RDMd%Nfrozen
     RDMd%t_pccd(iorb,iorb2)=RDMd%t_pccd(iorb,iorb2)+RDMd%t_pccd_old(iorb4,iorb2)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     RDMd%t_pccd(iorb,iorb2)=RDMd%t_pccd(iorb,iorb2)+y_ij(iorb,iorb4)*RDMd%t_pccd_old(iorb4,iorb2)
     sum_val=sum_val-RDMd%t_pccd_old(iorb4,iorb2)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
    enddo
    do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
     iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
     RDMd%t_pccd(iorb,iorb2)=RDMd%t_pccd(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
     sum_val=sum_val-RDMd%t_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
    enddo
    sum_val=sum_val-(two*INTEGd%ERImol(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1) &
    &        -INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd_old(iorb,iorb2))
    sum_val=two*sum_val*RDMd%t_pccd_old(iorb,iorb2)
    RDMd%t_pccd(iorb,iorb2)=RDMd%t_pccd(iorb,iorb2)+sum_val
    RDMd%t_pccd(iorb,iorb2)=RDMd%t_pccd_old(iorb,iorb2)-RDMd%t_pccd(iorb,iorb2) &
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
    Ecorr_new=Ecorr_new+RDMd%t_pccd(iorb,iorb2)*INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
    sumdiff_t=sumdiff_t+abs(RDMd%t_pccd(iorb,iorb2)-RDMd%t_pccd_old(iorb,iorb2))
    if(abs(RDMd%t_pccd(iorb,iorb2)-RDMd%t_pccd_old(iorb,iorb2))>maxdiff_t) then
     maxdiff_t=abs(RDMd%t_pccd(iorb,iorb2)-RDMd%t_pccd_old(iorb,iorb2))
    endif   
   enddo
  enddo

  ! Update old quantities
  RDMd%t_pccd_old=RDMd%t_pccd
  Ecorr_diff=Ecorr_new-Ecorr_old
  Ecorr_old=Ecorr_new

  ! Exit if converged
  if(maxdiff_t<tol6 .and. sumdiff_t<tol5) exit

  if(abs(Ecorr_diff)<tol8) then
   write(msg,'(a,f15.6)') 't-amplitudes converged for small Ecorr dif. ',Ecorr_diff
   call write_output(msg)
   exit
  endif

  ! Exit max iter
  if(iter_t==1000) exit

 enddo

 ! RDMd%z_pccd(RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
 ! z-amplitudes (with fixed t-amplitudes)

 ! Init. t_pccd
 if(iter_global==-1) then
  RDMd%z_pccd_old=RDMd%t_pccd
 endif
 iter_z=0

 ! Build intermediate y_ij = sum_b v_bb^ii t_j^b 
 y_ij=zero 
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  do iorb2=1,RDMd%Npairs ! Occ
   iorb3=iorb2+RDMd%Nfrozen
   do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
    y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd(iorb2,iorb4)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
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
    y_ab(iorb,iorb2)=y_ab(iorb,iorb2)+RDMd%t_pccd(iorb4,iorb2)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
   enddo
  enddo
 enddo

 do

  ! Build z_ia
  RDMd%z_pccd=zero
  do iorb=1,RDMd%Npairs ! Occ
   iorb1=iorb+RDMd%Nfrozen
   do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
    RDMd%z_pccd(iorb,iorb2)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
    sum_val=ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1)
    zt_aa=zero; zt_ii=zero;
    do iorb4=1,RDMd%Npairs ! Occ
     iorb5=iorb4+RDMd%Nfrozen
     RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd(iorb,iorb2)+RDMd%z_pccd_old(iorb4,iorb2)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd(iorb,iorb2)+y_ij(iorb,iorb4)*RDMd%z_pccd_old(iorb4,iorb2)
     sum_val=sum_val-RDMd%t_pccd(iorb4,iorb2)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
     zt_aa=zt_aa+RDMd%t_pccd(iorb4,iorb2)*RDMd%z_pccd_old(iorb4,iorb2)
    enddo
    do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
     iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
     RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd(iorb,iorb2)+RDMd%z_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
     RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd(iorb,iorb2)+y_ab(iorb,iorb4)*RDMd%z_pccd_old(iorb,iorb4)
     sum_val=sum_val-RDMd%t_pccd(iorb,iorb4)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     zt_ii=zt_ii+RDMd%t_pccd(iorb,iorb4)*RDMd%z_pccd_old(iorb,iorb4)
    enddo
    sum_val=sum_val-(two*INTEGd%ERImol(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1) &
    &        -two*INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd(iorb,iorb2))
    sum_val=two*sum_val*RDMd%z_pccd_old(iorb,iorb2)
    RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd(iorb,iorb2)-two*(zt_aa+zt_ii)*INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
    RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd(iorb,iorb2)+sum_val
    RDMd%z_pccd(iorb,iorb2)=RDMd%z_pccd_old(iorb,iorb2)-RDMd%z_pccd(iorb,iorb2) &
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

  ! Exit if converged
  RDMd%z_pccd_old=RDMd%z_pccd

  if(maxdiff_z<tol6 .and. sumdiff_z<tol5) exit

  ! Exit max iter
  if(iter_z==1000) exit

 enddo

 ! Calc. the final Energy using new RDMs
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
  Ediff=Ediff-Energy
  write(msg,'(a,f15.6)') 'Max. [t_pq^i+1 - t_pq^i]=      ',maxdiff_t
  call write_output(msg)
  write(msg,'(a,f15.6)') 'Max. [z_pq^i+1 - z_pq^i]=      ',maxdiff_z
  call write_output(msg)
  write(msg,'(a,f15.6)') 'Correlation Energy   =',Ecorr_new
  call write_output(msg)
  write(msg,'(a,f19.10)') 'Energy difference amp. opt.=',Ediff
  call write_output(msg)
 endif

 deallocate(y_ij,y_ab)

end subroutine calc_tz_pCCD_amplitudes
!!***

end module m_tz_pCCD_amplitudes
!!***
