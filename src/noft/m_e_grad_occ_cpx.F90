!!****m* DoNOF/m_e_grad_occ_cpx
!! NAME
!!  m_e_grad_occ_cpx
!!
!! FUNCTION
!!  Module prepared to compute energies from Gammas and gradients of the Energy w.r.t Gammas
!!    				 using complex integrals
!!
!!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  Nfrozen |            Npairs_p_sing     |              Nvirt                               = NBF               
!!  Nfrozen |         Npairs + Nsingleocc  |     Ncoupled*Npairs                   + Nempty   = NBF               
!!                           | Nsingleocc  |   NBF_occ - Npairs_p_sing - Nfrozen   | Nempty   = NBF
!!                           Nbeta         Nalpha                            NBF_occ
!!- - - - - - - - - - - - - - -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!!
!! PARENTS
!!  m_optocc
!!  m_optorb
!!
!! CHILDREN
!!  m_nofoutput
!!  m_rdmd
!!  m_gammatodm2
!!
!! SOURCE
module m_e_grad_occ_cpx

 use m_nofoutput
 use m_rdmd
 use m_gammatodm2

 implicit none

 private :: dm2_x_eri_cmplx,Ddm2_gamma_x_ERI_cmplx
!!***

 public :: calc_E_occ_cmplx,calc_Grad_occ_cmplx,num_calc_Grad_occ_cmplx,calc_Chem_pot_cmplx
!!***

contains

!!***
!!****f* DoNOF/calc_E_occ_cmplx
!! NAME
!!  calc_E_occ_cmplx
!!
!! FUNCTION
!!  Calculate the Energy from gamma independent parameters 
!!  Note: In the HF case, we play the trick that occ=occ^2 as occ=0,1
!!
!! INPUTS
!!  GAMMAs=Indep. variables used in the occ. optimization 
!!  hCORE=One-body integrals (h_pq) we only use the diag part (h_pp)
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  ERI_L=Lower triangular part of the L_pq matrix
!!  nogamma=Do not build Occ, DM2_J, DM2_K, etc from GAMMAs (use the stored ones).
!!  chempot=Create the DM2 and the DDM2_w.r.t_occs matrices.
!!
!! OUTPUT
!!  Energy = Energy computed from the occs (actually from gammas)
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_E_occ_cmplx(RDMd,GAMMAs,Energy,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx,&
 &                          nogamma,chempot) 
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::nogamma,chempot
 real(dp),intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%Ngammas),intent(in)::GAMMAs
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr_cmplx,ERI_Lsr_cmplx 
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE_cmplx
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,ipair
 complex(dp)::Energy_cmplx
!arrays
 character(len=200)::msg
!************************************************************************

 if(.not.present(nogamma)) then 
  if(present(chempot)) then
   call gamma_to_2rdm(RDMd,GAMMAs,chempot=chempot)
  else
   call gamma_to_2rdm(RDMd,GAMMAs)
  endif
 endif 

 Energy=zero; Energy_cmplx=complex_zero;

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if(RDMd%Ncoupled==1) then     ! Perfect Pairing (Ncoupled=1)
   do iorb=1,RDMd%Nfrozen
    Energy_cmplx = Energy_cmplx + RDMd%occ(iorb) * two*hCORE_cmplx(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                                                &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_J,ERI_J_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_K,ERI_K_cmplx)    &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOFT 
     Energy_cmplx = Energy_cmplx + RDMd%DM2_iiii(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                        &
    &       + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOFT 
     Energy_cmplx = Energy_cmplx + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) &
    &       + dm2_x_eri_cmplx(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx)
    endif
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy_cmplx = Energy_cmplx + RDMd%occ(iorb) * two*hCORE_cmplx(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                                                &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_J,ERI_J_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_K,ERI_K_cmplx)    &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOFT 
     Energy_cmplx = Energy_cmplx + RDMd%DM2_iiii(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                                   &
    &       + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOFT 
     Energy_cmplx = Energy_cmplx + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) &
    &       + dm2_x_eri_cmplx(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx)
    endif
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
    Energy_cmplx = Energy_cmplx + RDMd%occ(iorb) * two*hCORE_cmplx(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                                                &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_J,ERI_J_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_K,ERI_K_cmplx)    &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOFT 
     Energy_cmplx = Energy_cmplx + RDMd%DM2_iiii(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                                   &
    &       + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOFT 
     Energy = Energy + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) &
    &       + dm2_x_eri_cmplx(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx)
    endif
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                  ! Extended PNOF (Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do iorb=1,RDMd%Nfrozen
    Energy_cmplx = Energy_cmplx + RDMd%occ(iorb) * two*hCORE_cmplx(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                                                &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_J,ERI_J_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_K,ERI_K_cmplx)    &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOFT 
     Energy_cmplx = Energy_cmplx + RDMd%DM2_iiii(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                        &
    &       + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOFT 
     Energy_cmplx = Energy_cmplx + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) &
    &        + dm2_x_eri_cmplx(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx)
    endif
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy_cmplx = Energy_cmplx + RDMd%occ(iorb) * two*hCORE_cmplx(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                                                &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_J,ERI_J_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_K,ERI_K_cmplx)    &
    &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOFT 
     Energy_cmplx = Energy_cmplx + RDMd%DM2_iiii(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                        &
    &       + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOFT 
     Energy_cmplx = Energy_cmplx + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) &
    &        + dm2_x_eri_cmplx(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx)
    endif
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
     Energy_cmplx = Energy_cmplx + RDMd%occ(iorb) * two*hCORE_cmplx(iorb,iorb)                                 &
     &      + RDMd%DM2_iiii(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                                               &
     &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_J,ERI_J_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_K,ERI_K_cmplx)   &
     &      + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_L,ERI_L_cmplx)
     if(RDMd%irange_sep==1) then ! Intra rs-NOFT 
      Energy_cmplx = Energy_cmplx + RDMd%DM2_iiii(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                       &
     &       + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Lsr,ERI_Lsr_cmplx)
     endif
     if(RDMd%irange_sep==2) then ! Hartree rs-NOFT 
      Energy_cmplx = Energy_cmplx + dm2_x_eri_cmplx(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx) &
     &        + dm2_x_eri_cmplx(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr_cmplx)
     endif

    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 else ! TODO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 Energy=real(Energy_cmplx)
 if(dabs(aimag(Energy_cmplx))>tol8) then
  write(msg,'(a,f10.5,a)') 'Warning! Imaginary[Energy] = ',aimag(Energy_cmplx),' (a.u.)'
  call write_output(msg)
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_E_occ_cmplx
!!***

!!***
!!****f* DoNOF/calc_Grad_occ_cmplx
!! NAME
!!  calc_Grad_occ_cmplx
!!
!! FUNCTION
!!  Calculate the Gradient of the energy from gamma independent parameters 
!!
!! INPUTS
!!  hCORE=One-body integrals (h_pq) we only use the diag part (h_pp)
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  ERI_L=Lower triangular part of the L_pq matrix
!!
!! OUTPUT
!!  Energy=Energy computed from the occs (actually from gammas)
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_Grad_occ_cmplx(RDMd,Grad,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 real(dp),dimension(RDMd%Ngammas),intent(inout)::Grad
!arrays
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr_cmplx,ERI_Lsr_cmplx
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE_cmplx
!Local variables ------------------------------
!scalars
 integer::igamma,iorb,iorb1,ipair
 real(dp)::max_imag
!arrays
 character(len=200)::msg
 complex(dp),allocatable,dimension(:)::Grad_cmplx
!************************************************************************
 
 allocate(Grad_cmplx(RDMd%Ngammas))

 Grad = zero; Grad_cmplx=complex_zero;

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE_cmplx(iorb,iorb) &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)         &
    &         + two * ( Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)                          &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)                                  &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx) )
     if(RDMd%irange_sep==1) then ! Intra rs-NOFT
     Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb)        &
    &         * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                                                                           &
    &         + two * ( Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)                     &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx) )
     endif
     if(RDMd%irange_sep==2) then ! Hartree rs-NOFT
      Grad_cmplx(igamma) = Grad_cmplx(igamma) + two * Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
     endif      
     iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
     Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE_cmplx(iorb,iorb) &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)         &
    &         + two * ( Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)                          &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)                                  &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx) )
     if(RDMd%irange_sep==1) then ! Intra rs-NOFT
     Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb)        &
    &         * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                                                                           &
    &         + two * ( Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)                     &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx) )
     endif
     if(RDMd%irange_sep==2) then ! Hartree rs-NOFT
      Grad_cmplx(igamma) = Grad_cmplx(igamma) + two * Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
     endif      
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                 ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE_cmplx(iorb,iorb) &
     &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)         &
     &        + two * ( Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)                          &
     &        + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)                                  &
     &        + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx) )
     if(RDMd%irange_sep==1) then ! Intra rs-NOFT
     Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb)        &
    &         * ERI_Jsr_cmplx(iorb*(iorb+1)/2)                                                                           &
    &         + two * ( Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)                     &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx) )
     endif
     if(RDMd%irange_sep==2) then ! Hartree rs-NOFT
      Grad_cmplx(igamma) = Grad_cmplx(igamma) + two * Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
     endif      
     do iorb1=1,RDMd%Ncoupled
      iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
      Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE_cmplx(iorb,iorb) &
     &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)         &
     &         + two * (Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)                           &
     &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)                                  &
     &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx) )
      if(RDMd%irange_sep==1) then ! Intra rs-NOFT
       Grad_cmplx(igamma) = Grad_cmplx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb)      &
    &          *  ERI_Jsr_cmplx(iorb*(iorb+1)/2)                                                                         &
    &          + two * ( Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)                    &
    &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx) )
      endif
      if(RDMd%irange_sep==2) then ! Hartree rs-NOFT
       Grad_cmplx(igamma) = Grad_cmplx(igamma)+two*Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
    &          + Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
      endif      
     enddo
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
 else ! TODO
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 Grad(:)=real(Grad_cmplx(:))
 max_imag=maxval(aimag(Grad_cmplx(:)))
 if(dabs(max_imag)>tol8) then
  write(msg,'(a,f10.5,a)') 'Warning! Imaginary[Grad_Energy] = ',max_imag
  call write_output(msg)
 endif
  
 deallocate(Grad_cmplx)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_Grad_occ_cmplx
!!***

!!***
!!****f* DoNOF/num_calc_Grad_occ_cmplx
!! NAME
!!  num_calc_Grad_occ_cmplx
!!
!! FUNCTION
!!  Calculate the (numerical) Gradient of the energy from gamma independent parameters 
!!
!! INPUTS
!!  GAMMAs=Indep. variables used in the occ. optimization 
!!  hCORE=One-body integrals (h_pq) we only use the diag part (h_pp)
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  ERI_L=Lower triangular part of the L_pq matrix
!!
!! OUTPUT
!!  Energy=Energy computed from the occs (actually from gammas)
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine num_calc_Grad_occ_cmplx(RDMd,GAMMAs,Grad,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 real(dp),dimension(RDMd%Ngammas),intent(inout)::Grad
!arrays
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr_cmplx,ERI_Lsr_cmplx
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE_cmplx
 real(dp),dimension(RDMd%Ngammas),intent(in)::GAMMAs
!Local variables ------------------------------
!scalars
 integer::igamma
 real(dp)::Energy_num,grad_igamma,step=tol3
!arrays
 real(dp),allocatable,dimension(:)::GAMMAs_num
!************************************************************************
 
 allocate(GAMMAs_num(RDMd%Ngammas)) 
 Grad = zero
 do igamma=1,RDMd%Ngammas
  GAMMAs_num=GAMMAs;
  ! 2*step
  GAMMAS_num(igamma)=GAMMAS_num(igamma)+two*step
  call calc_E_occ_cmplx(RDMd,GAMMAs_num,Energy_num,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx)
  grad_igamma=-Energy_num
  ! step
  GAMMAS_num(igamma)=GAMMAS_num(igamma)-step
  call calc_E_occ_cmplx(RDMd,GAMMAs_num,Energy_num,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx)
  grad_igamma=grad_igamma+eight*Energy_num
  ! -step 
  GAMMAS_num(igamma)=GAMMAS_num(igamma)-two*step
  call calc_E_occ_cmplx(RDMd,GAMMAs_num,Energy_num,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx)
  grad_igamma=grad_igamma-eight*Energy_num
  ! -2step 
  GAMMAS_num(igamma)=GAMMAS_num(igamma)-step
  call calc_E_occ_cmplx(RDMd,GAMMAs_num,Energy_num,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx)
  grad_igamma=grad_igamma+Energy_num
  ! Save the gradient
  Grad(igamma)=grad_igamma/(twelve*step)
 enddo
 deallocate(GAMMAs_num) 

end subroutine num_calc_Grad_occ_cmplx
!!***

!!***
!!****f* DoNOF/calc_Chem_pot_cmplx
!! NAME
!!  calc_Chem_pot_cmplx
!!
!! FUNCTION
!!  Calculate the chemical potential for each orbital
!!
!! INPUTS
!!  hCORE=One-body integrals (h_pq) we only use the diag part (h_pp)
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  ERI_L=Lower triangular part of the L_pq matrix
!!
!! OUTPUT
!!  RDMd%chempot_orb=conatins chem. pot. for each orbital
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_Chem_pot_cmplx(RDMd,hCORE_cmplx,ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx,ERI_Jsr_cmplx,ERI_Lsr_cmplx) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_Jsr_cmplx,ERI_Lsr_cmplx
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE_cmplx
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,ipair,igamma=1
 real(dp)::max_imag
!arrays
 complex(dp),allocatable,dimension(:)::chempot_orb
 character(len=200)::msg
!************************************************************************
 
 allocate(chempot_orb(RDMd%NBF_occ))

 RDMd%chempot_orb=zero; chempot_orb=complex_zero; 

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
         chempot_orb(iorb) = hCORE_cmplx(iorb,iorb)                                       &
   &         + half * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                   &
   &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)   &
   &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)   &
   &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx) 
    if(RDMd%irange_sep==1) then ! Intra rs-NOF 
     RDMd%chempot_orb(iorb) = half * RDMd%Dfni_ni(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)      &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)  &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOF 
     RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                             &
   &      + half * Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
   &      + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
    endif
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
         chempot_orb(iorb) = hCORE_cmplx(iorb,iorb)                                       &
   &         + half * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                   &
   &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)   &
   &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)   &
   &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOF 
     RDMd%chempot_orb(iorb) = half * RDMd%Dfni_ni(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)      &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)  &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOF 
     RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                             &
   &      + half * Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
   &      + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
    endif
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                 ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
         chempot_orb(iorb) = hCORE_cmplx(iorb,iorb)                                       &
    &        + half * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                   &
    &        + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)   &
    &        + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)   &
    &        + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOF 
     RDMd%chempot_orb(iorb) = half * RDMd%Dfni_ni(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)        &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)    &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx)
    endif
    if(RDMd%irange_sep==2) then ! Hartree rs-NOF 
     RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                             &
   &      + half * Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
   &      + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
    endif
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
          chempot_orb(iorb) = hCORE_cmplx(iorb,iorb)                                      &
    &         + half * RDMd%Dfni_ni(iorb) * ERI_J_cmplx(iorb*(iorb+1)/2)                  &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J_cmplx)  &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K_cmplx)  &
    &         + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L_cmplx)
    if(RDMd%irange_sep==1) then ! Intra rs-NOF 
      RDMd%chempot_orb(iorb) = half * RDMd%Dfni_ni(iorb) * ERI_Jsr_cmplx(iorb*(iorb+1)/2)    &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
   &          + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Lsr,ERI_Lsr_cmplx)
     endif
     if(RDMd%irange_sep==2) then ! Hartee rs-NOF 
      RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                             &
   &       + half * Ddm2_gamma_x_ERI_cmplx(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx) &
   &       + Ddm2_gamma_x_ERI_cmplx(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr_cmplx)
     endif
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
 endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
 RDMd%chempot_orb(:)=real(chempot_orb(:))
 max_imag=maxval(aimag(chempot_orb(:)))
 if(dabs(max_imag)>tol8) then
  write(msg,'(a,f10.5,a)') 'Warning! Imaginary[Chem. pot.] = ',max_imag
  call write_output(msg)
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 deallocate(chempot_orb)

end subroutine calc_Chem_pot_cmplx
!!***

!!***
!!****f* DoNOF/dm2_x_eri_cmplx
!! NAME
!!  dm2_x_eri_cmplx
!!
!! FUNCTION
!!  Multiply the 2RDM elements by the two electron repulsion integrals (ERIs) to produce energy contributions.
!!  Note: Term with iorb=iorb1 is not included
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

function dm2_x_eri_cmplx(RDMd,icase,iorb,DM2_JKL,ERI) result(E_dm2ERI_iorb)
!Arguments ------------------------------------
!scalars
 complex(dp)::E_dm2ERI_iorb
 integer,intent(in)::icase,iorb
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(in)::DM2_JKL 
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI 
!Local variables ------------------------------
!scalars
 integer::iorb1
!arrays
!************************************************************************
E_dm2ERI_iorb = complex_zero
select case(icase)
!-----------------------------------------------------------------------
 case(-1)
!---------------------------------
!    DM2_JKL*ERI.
!---------------------------------
  E_dm2ERI_iorb = DM2_JKL(iorb,iorb)*ERI(iorb*(iorb+1)/2)
!-----------------------------------------------------------------------
 case(0)
!---------------------------------
!    DM2_JKL*ERI.
!---------------------------------
  do iorb1=1,iorb-1
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%NBF_occ
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*conjg(ERI(iorb+iorb1*(iorb1-1)/2))
  enddo
!-----------------------------------------------------------------------
end select

end function dm2_x_eri_cmplx
!!***

!!***
!!****f* DoNOF/Ddm2_gamma_x_ERI_cmplx
!! NAME
!!  Ddm2_gamma_x_ERI_cmplx
!!
!! FUNCTION
!!  Multiply the Derivative of DM2 w.r.t. gamma by the two electron repulsion integrals (ERIs) to produce gradients.
!!  Note: Term with iorb=iorb1 is not included
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

function Ddm2_gamma_x_ERI_cmplx(RDMd,icase,iorb,igamma,DDM2_JorKorL,ERI) result(Grad_Ddm2_ERI_iorb)
!Arguments ------------------------------------
!scalars
 complex(dp)::Grad_Ddm2_ERI_iorb
 integer,intent(in)::icase,iorb,igamma
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(in)::DDM2_JorKorL
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI 
!Local variables ------------------------------
!scalars
 integer::iorb1
!arrays
!************************************************************************
Grad_Ddm2_ERI_iorb = complex_zero
select case(icase)
!-----------------------------------------------------------------------
 case(-1)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI. 
!-----------------------------------------------------------------------
  Grad_Ddm2_ERI_iorb = DDM2_JorKorL(iorb,iorb,igamma)*real(ERI(iorb*(iorb+1)/2))
!-----------------------------------------------------------------------
 case(0)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI. 
!-----------------------------------------------------------------------
  do iorb1=1,iorb-1
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*real(ERI(iorb1+iorb*(iorb-1)/2))
  enddo
  do iorb1=iorb+1,RDMd%NBF_occ
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*real(ERI(iorb+iorb1*(iorb1-1)/2))
  enddo
!-----------------------------------------------------------------------
end select

end function Ddm2_gamma_x_ERI_cmplx
!!***

end module m_e_grad_occ_cpx 
!!***
