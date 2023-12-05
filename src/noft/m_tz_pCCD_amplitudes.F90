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

subroutine calc_tz_pCCD_amplitudes(ELAGd,RDMd,INTEGd,Vnn,Energy) 
!Arguments ------------------------------------
!scalars
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
!Local variables ------------------------------
!scalars
 logical::dodiis,nogamma
 integer::iter,iorb
 real(dp)::sumdiff,maxdiff_t,maxdiff_z,Ediff,Energy_old
!arrays
 character(len=200)::msg
!************************************************************************

 Ediff=zero; Energy=zero; Energy_old=zero; nogamma=.true.;
 maxdiff_t=zero; maxdiff_z=zero;
 iter=-1;

 !! Building single-det. ELAGd!!
 ! Build Lambda matrix (with HF 1-RDM and 2-RDM)
 call ELAGd%build_sd_diag(RDMd,INTEGd)

 write(*,*)' '
 write(*,*) INTEGd%NBF_jkl
 do iorb=1,RDMd%NBF_tot
  write(*,*) ELAGd%Lambdas_pp(iorb)
 enddo
 write(*,*) ' '




 iter=10





 ! Calc. the final Energy using new RDMs
 if(INTEGd%complex_ints) then
  call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
  & INTEGd%ERI_L_cmplx,nogamma=nogamma)
 else
  call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
  & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,nogamma=nogamma)
 endif
 write(msg,'(a,f15.6,a,i6,a)') 'T-,Z-amp. opt. energy= ',Energy+Vnn,' after ',iter,' iter.'
 call write_output(msg)
 if(iter>0) then
  write(msg,'(a,f15.6)') 'Max. [t_pq^i+1 - t_pq^i]= ',maxdiff_t
  call write_output(msg)
  write(msg,'(a,f15.6)') 'Max. [z_pq^i+1 - z_pq^i]= ',maxdiff_z
  call write_output(msg)
  write(msg,'(a,f19.10)') 'Energy difference amp. opt. =',Ediff
  call write_output(msg)
 endif

end subroutine calc_tz_pCCD_amplitudes
!!***

end module m_tz_pCCD_amplitudes
!!***
