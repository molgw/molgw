!!****m* DoNOF/m_e_grad_occ
!! NAME
!!  m_e_grad_occ
!!
!! FUNCTION
!!  Module prepared to compute energies from Gammas and gradients of the Energy w.r.t Gammas
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
!!  m_gammatodm2
!!
!! SOURCE
module m_e_grad_occ

 use m_rdmd
 use m_gammatodm2

 implicit none

 private :: dm2_x_eri,occ_x_eri,Ddm2_gamma_x_ERI,Docc_gamma_x_ERI
!!***

 public :: calc_E_occ,calc_Grad_occ
!!***

contains

!!***
!!****f* DoNOF/calc_E_occ
!! NAME
!!  calc_E_occ
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
!!  nogamma=Do not build OCC, DM2_J, DM2_K, etc from GAMMAs (use the stored ones).
!!
!! OUTPUT
!!  Energy=Energy computed from the occs (actually from gammas)
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L,nogamma) 
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::nogamma
 real(dp),intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%Ngammas),intent(in)::GAMMAs
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L 
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,ipair
!arrays
!************************************************************************

 if(.not.present(nogamma)) then 
  call gamma_to_2rdm(RDMd,GAMMAs)
 endif
 Energy=0.0d0
 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if(RDMd%Ncoupled==1) then     ! Perfect Pairing (Ncoupled=1)
   do iorb=1,RDMd%Nfrozen
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                         &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                         &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                         &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                  ! Extended PNOF (Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do iorb=1,RDMd%Nfrozen
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                         &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                         &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
     Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                        &
     &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                            &
     &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)       &
     &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 else 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do iorb=1,RDMd%Nfrozen
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                      &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                          &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_K,ERI_K)     &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_L,ERI_L)                                               &
    &      + 2.0d0*occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_K)
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                      &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                          &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_K,ERI_K)     &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_L,ERI_L)                                               &
    &      + 2.0d0*occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_K)
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                      &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                          &
    &      + dm2_x_eri(RDMd,2,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,2,iorb,RDMd%DM2_K,ERI_K)     &
    &      + dm2_x_eri(RDMd,2,iorb,RDMd%DM2_L,ERI_L)                                               &
    &      + 2.0d0*occ_x_eri(RDMd,2,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,2,iorb,RDMd%occ,ERI_K)
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   do ipair=RDMd%Npairs+1,RDMd%Npairs_p_sing
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb)*hCORE(iorb,iorb)                                              &
    &      + 0.5d0*(occ_x_eri(RDMd,0,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,0,iorb,RDMd%occ,ERI_K))
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                  ! PNOFi(Nc): Extended PNOF (Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do iorb=1,RDMd%Nfrozen
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                      &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                          &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_K,ERI_K)     &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_L,ERI_L)                                               &
    &      + 2.0d0*occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_K)
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                      &
    &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                          &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_K,ERI_K)     &
    &      + dm2_x_eri(RDMd,1,iorb,RDMd%DM2_L,ERI_L)                                               &
    &      + 2.0d0*occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,1,iorb,RDMd%occ,ERI_K)
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
     Energy = Energy + RDMd%occ(iorb) * 2.0d0*hCORE(iorb,iorb)                                     &
     &      + RDMd%DM2_IIII(iorb) * ERI_J(iorb*(iorb+1)/2)                                         &
     &      + dm2_x_eri(RDMd,2,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,2,iorb,RDMd%DM2_K,ERI_K)    &
     &      + dm2_x_eri(RDMd,2,iorb,RDMd%DM2_L,ERI_L)                                              &
     &      + 2.0d0*occ_x_eri(RDMd,2,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,2,iorb,RDMd%occ,ERI_K)
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   do ipair=RDMd%Npairs+1,RDMd%Npairs_p_sing
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb)*hCORE(iorb,iorb)                                              &
    &      + 0.5d0*(occ_x_eri(RDMd,0,iorb,RDMd%occ,ERI_J) - occ_x_eri(RDMd,0,iorb,RDMd%occ,ERI_K))
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_E_occ
!!***

!!***
!!****f* DoNOF/calc_Grad_occ
!! NAME
!!  calc_Grad_occ
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

subroutine calc_Grad_occ(RDMd,Grad,hCORE,ERI_J,ERI_K,ERI_L) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 real(dp),dimension(RDMd%Ngammas),intent(inout)::Grad
!arrays
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L 
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
!Local variables ------------------------------
!scalars
 integer::igamma,iorb,iorb1,ipair
!arrays
!************************************************************************
 
 Grad = 0.0d0
 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)    &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)  &
    &         + 2.0d0 * ( Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                       &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                 &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
     iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
     Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)    &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)  &
    &         + 2.0d0 * ( Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                       &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                 &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                 ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)    &
     &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)  &
     &        + 2.0d0 * ( Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                       &
     &        + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                 &
     &        + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
     do iorb1=1,RDMd%Ncoupled
      iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
      Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)   &
     &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2) &
     &         + 2.0d0 * (Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                       &
     &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                &
     &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
     enddo
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
 else
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     GRAD(igamma) = GRAD(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)   &
     &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2) &
     &        + 2.0d0 * ( Ddm2_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                      &
     &        +           Ddm2_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                      &
     &        +           Ddm2_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L))                     &
     &        + 2.0d0 *   Docc_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%Docc_gamma,ERI_J)                        &
     &        +           Docc_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%Docc_gamma,ERI_K)
     iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
     GRAD(igamma) = GRAD(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)   &
     &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2) &
     &        + 2.0d0 * ( Ddm2_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                      &
     &        +           Ddm2_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                      &
     &        +           Ddm2_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L))                     &
     &        + 2.0d0 *   Docc_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%Docc_gamma,ERI_J)                        &
     &        +           Docc_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%Docc_gamma,ERI_K)
    enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  else                  ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     GRAD(igamma) = GRAD(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)    &
     &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)  &
     &        + 2.0d0 * ( Ddm2_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                       &
     &        +          Ddm2_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                        &
     &        +          Ddm2_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L))                       &
     &        + 2.0d0 *   Docc_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%Docc_gamma,ERI_J)                         &
     &        +           Docc_gamma_x_ERI(RDMd,1,iorb,igamma,RDMd%Docc_gamma,ERI_K)
     do iorb1=1,RDMd%Ncoupled
      iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
      GRAD(igamma) = GRAD(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * 2.0d0*hCORE(iorb,iorb)   &
      &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2) &
      &        + 2.0d0 * (Ddm2_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                       &
      &        +          Ddm2_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                       &
      &        +          Ddm2_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L))                      &
      &        + 2.0d0 *  Docc_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%Docc_gamma,ERI_J)                         &
      &        +          Docc_gamma_x_ERI(RDMd,2,iorb,igamma,RDMd%Docc_gamma,ERI_K)
     enddo
    enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -          
   enddo
  endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_Grad_occ
!!***

!!***
!!****f* DoNOF/dm2_x_eri
!! NAME
!!  dm2_x_eri
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

function dm2_x_eri(RDMd,icase,iorb,DM2_JKL,ERI) result(E_dm2ERI_iorb)
!Arguments ------------------------------------
!scalars
 real(dp)::E_dm2ERI_iorb
 integer,intent(in)::icase,iorb
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(in)::DM2_JKL 
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI 
!Local variables ------------------------------
!scalars
 integer::iorb1
!arrays
!************************************************************************
E_dm2ERI_iorb = 0.0d0
select case(icase)
!-----------------------------------------------------------------------
 case(0)
!---------------------------------
!    DM2_JKL*ERI.
!---------------------------------
  do iorb1=1,iorb-1
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%NBF_occ
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
 case(1)
!-----------------------------------------------------------------------
!    DM2_JKL*ERI, iorb<Nbeta_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
  do iorb1=1,iorb-1
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%Nbeta_elect
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
  do iorb1=RDMd%Nalpha_elect+1,RDMd%NBF_occ
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
 case(2)
!-----------------------------------------------------------------------
!    DM2_JKL*ERI, iorb>Nalpha_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
  do iorb1=1,RDMd%Nbeta_elect
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=RDMd%Nalpha_elect+1,iorb-1
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%NBF_occ
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
end select

end function dm2_x_eri
!!***

!!***
!!****f* DoNOF/occ_x_eri
!! NAME
!!  occ_x_eri
!!
!! FUNCTION
!!  Multiply the OCCs by the two electron repulsion integrals (ERIs) to produce energy contributions.
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

function occ_x_eri(RDMd,icase,iorb,OCC,ERI) result(E_occERI_iorb)
!Arguments ------------------------------------
!scalars
 real(dp)::E_occERI_iorb
 integer,intent(in)::icase,iorb
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::OCC
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI 
!Local variables ------------------------------
!scalars
 integer::iorb1
!arrays
!************************************************************************
E_occERI_iorb = 0.0d0
select case(icase)
!-----------------------------------------------------------------------
 case(0)
!-----------------------------------------------------------------------
!     OCC*ERI, Sum only for Nsingleocc
!-----------------------------------------------------------------------
  do iorb1=RDMd%Nbeta_elect+1,iorb-1
   E_occERI_iorb = E_occERI_iorb + OCC(iorb1)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%Nalpha_elect
   E_occERI_iorb = E_occERI_iorb + OCC(iorb1)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!--------------------------------------------------------------------      
 case(1)
!-----------------------------------------------------------------------
!     OCC*ERI, iorb<Nbeta_elect<iorb1 (Sum only for Nsingleocc)
!-----------------------------------------------------------------------
  do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   E_occERI_iorb = E_occERI_iorb + OCC(iorb)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------      
 case(2)
!-----------------------------------------------------------------------
!     OCC*ERI, iorb>Nalpha_elect>iorb1 (Sum only for Nsingleocc)
!-----------------------------------------------------------------------
  do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   E_occERI_iorb = E_occERI_iorb + OCC(iorb)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
!-----------------------------------------------------------------------
end select

end function occ_x_eri
!!***

!!***
!!****f* DoNOF/Ddm2_gamma_x_ERI
!! NAME
!!  Ddm2_gamma_x_ERI
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

function Ddm2_gamma_x_ERI(RDMd,icase,iorb,igamma,DDM2_JorKorL,ERI) result(Grad_Ddm2_ERI_iorb)
!Arguments ------------------------------------
!scalars
 real(dp)::Grad_Ddm2_ERI_iorb
 integer,intent(in)::icase,iorb,igamma
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(in)::DDM2_JorKorL
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI 
!Local variables ------------------------------
!scalars
 integer::iorb1
!arrays
!************************************************************************
Grad_Ddm2_ERI_iorb = 0.0d0
select case(icase)
!-----------------------------------------------------------------------
 case(0)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI. 
!-----------------------------------------------------------------------
  do iorb1=1,iorb-1
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%NBF_occ
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
 case(1)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI, iorb<Nbeta_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
  do iorb1=1,iorb-1
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%Nbeta_elect
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
  do iorb1=RDMd%Nalpha_elect+1,RDMd%NBF_occ
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
 case(2)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI, iorb>Nalpha_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
  do iorb1=1,RDMd%Nbeta_elect
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=RDMd%Nalpha_elect+1,iorb-1
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
  do iorb1=iorb+1,RDMd%NBF_occ
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
end select

end function Ddm2_gamma_x_ERI
!!***

!!***
!!****f* DoNOF/Docc_gamma_x_ERI
!! NAME
!!  Docc_gamma_x_ERI
!!
!! FUNCTION
!!  Multiply the Derivative of occ w.r.t. gamma by the two electron repulsion integrals (ERIs) to produce gradients.
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

function Docc_gamma_x_ERI(RDMd,icase,iorb,igamma,Docc_gamma,ERI) result(Grad_Docc_ERI_iorb)
!Arguments ------------------------------------
!scalars
 real(dp)::Grad_Docc_ERI_iorb
 integer,intent(in)::icase,iorb,igamma
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Docc_gamma
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI 
!Local variables ------------------------------
!scalars
 integer::iorb1
!arrays
!************************************************************************
Grad_Docc_ERI_iorb = 0.0d0
select case(icase)
!-----------------------------------------------------------------------
 case(1)
!-----------------------------------------------------------------------
!     Docc_gamma*ERI, iorb<Nbeta_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
  do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   Grad_Docc_ERI_iorb = Grad_Docc_ERI_iorb + Docc_gamma(iorb,igamma)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
 case(2)
!-----------------------------------------------------------------------
!     Docc_gamma*ERI, iorb>Nalpha_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
  do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   Grad_Docc_ERI_iorb = Grad_Docc_ERI_iorb + Docc_gamma(iorb,igamma)*ERI(iorb1+iorb*(iorb-1)/2)
  enddo
!-----------------------------------------------------------------------
end select

end function Docc_gamma_x_ERI
!!***

end module m_e_grad_occ 
!!***
