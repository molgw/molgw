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
!!  m_rdmd
!!  m_gammatodm2
!!
!! SOURCE
module m_e_grad_occ

 use m_rdmd
 use m_gammatodm2

 implicit none

 private :: dm2_x_eri,Ddm2_gamma_x_ERI
!!***

 public :: calc_E_occ,calc_Grad_occ,num_calc_Grad_occ,calc_Chem_pot
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
!!  ERI_Jsr=Lower triangular part of the Jsr_pq matrix (zero when rs-NOFT=false)
!!  nogamma=Do not build occ, DM2_J, DM2_K, etc from GAMMAs (use the stored ones).
!!  chempot=Create the DM2 and the DDM2_w.r.t_occs matrices.
!!
!! OUTPUT
!!  Energy=Energy computed from the occs (actually from gammas)
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr,nogamma,chempot) 
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::nogamma,chempot
 real(dp),intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%Ngammas),intent(in)::GAMMAs
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L,ERI_Jsr 
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,ipair
!arrays
!************************************************************************

 Energy=zero
 if(.not.present(nogamma)) then 
  if(present(chempot)) then
   call gamma_to_2rdm(RDMd,GAMMAs,chempot=chempot)
  else
   call gamma_to_2rdm(RDMd,GAMMAs)
  endif
 endif 

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if(RDMd%Ncoupled==1) then     ! Perfect Pairing (Ncoupled=1)
   do iorb=1,RDMd%Nfrozen
    Energy = Energy + RDMd%occ(iorb) * two*hCORE(iorb,iorb)                                           &
    &      + RDMd%DM2_iiii(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    Energy = Energy + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr) + dm2_x_eri(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr) ! rs-NOFT
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb) * two*hCORE(iorb,iorb)                                           &
    &      + RDMd%DM2_iiii(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    Energy=Energy+dm2_x_eri(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr)+dm2_x_eri(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr)       ! rs-NOFT
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
    Energy = Energy + RDMd%occ(iorb) * two*hCORE(iorb,iorb)                                           &
    &      + RDMd%DM2_iiii(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    Energy = Energy + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr) + dm2_x_eri(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr) ! rs-NOFT
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                  ! Extended PNOF (Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do iorb=1,RDMd%Nfrozen
    Energy = Energy + RDMd%occ(iorb) * two*hCORE(iorb,iorb)                                           &
    &      + RDMd%DM2_iiii(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    Energy = Energy + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr) + dm2_x_eri(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr) ! rs-NOFT
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy = Energy + RDMd%occ(iorb) * two*hCORE(iorb,iorb)                                           &
    &      + RDMd%DM2_iiii(iorb) * ERI_J(iorb*(iorb+1)/2)                                             &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)        &
    &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
    Energy = Energy + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr) + dm2_x_eri(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr) ! rs-NOFT
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
     Energy = Energy + RDMd%occ(iorb) * two*hCORE(iorb,iorb)                                          &
     &      + RDMd%DM2_iiii(iorb) * ERI_J(iorb*(iorb+1)/2)                                            &
     &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_J,ERI_J) + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_K,ERI_K)       &
     &      + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_L,ERI_L)
     Energy = Energy + dm2_x_eri(RDMd,0,iorb,RDMd%DM2_Jsr,ERI_Jsr) + dm2_x_eri(RDMd,-1,iorb,RDMd%DM2_Jsr,ERI_Jsr) ! rs-NOFT
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 else ! TODO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
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
!!  ERI_Jsr=Lower triangular part of the Jsr_pq matrix (zero when rs-NOFT=false)
!!
!! OUTPUT
!!  Energy=Energy computed from the occs (actually from gammas)
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_Grad_occ(RDMd,Grad,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 real(dp),dimension(RDMd%Ngammas),intent(inout)::Grad
!arrays
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L,ERI_Jsr 
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
!Local variables ------------------------------
!scalars
 integer::igamma,iorb,iorb1,ipair
!arrays
!************************************************************************
 
 Grad = zero
 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE(iorb,iorb)      &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)  &
    &         + two * ( Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                         &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                 &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
     Grad(igamma) = Grad(igamma) + two * Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)    & ! rs-NOFT 
    &         + Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)       
     iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
     Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE(iorb,iorb)      &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)  &
    &         + two * ( Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                         &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                 &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
     Grad(igamma) = Grad(igamma) + two * Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)    & ! rs-NOFT 
    &         + Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)       
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                 ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE(iorb,iorb)      &
     &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)  &
     &        + two * ( Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                         &
     &        + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                 &
     &        + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
     Grad(igamma) = Grad(igamma) + two * Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)    & ! rs-NOFT
     &        + Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)       
     do iorb1=1,RDMd%Ncoupled
      iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
      Grad(igamma) = Grad(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCORE(iorb,iorb)     &
     &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2) &
     &         + two * (Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)                         &
     &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)                                &
     &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) )
      Grad(igamma) = Grad(igamma) + two * Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)   & ! rs-NOFT    
     &         + Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)       
     enddo
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
 else ! TODO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      High-Spin Multiplet State (S>0,Ms=S)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_Grad_occ
!!***

!!***
!!****f* DoNOF/num_calc_Grad_occ
!! NAME
!!  num_calc_Grad_occ
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
!!  ERI_Jsr=Lower triangular part of the Jsr_pq matrix
!!
!! OUTPUT
!!  Energy=Energy computed from the occs (actually from gammas)
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine num_calc_Grad_occ(RDMd,GAMMAs,Grad,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 real(dp),dimension(RDMd%Ngammas),intent(inout)::Grad
!arrays
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L,ERI_Jsr
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
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
  call calc_E_occ(RDMd,GAMMAs_num,Energy_num,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr)
  grad_igamma=-Energy_num
  ! step
  GAMMAS_num(igamma)=GAMMAS_num(igamma)-step
  call calc_E_occ(RDMd,GAMMAs_num,Energy_num,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr)
  grad_igamma=grad_igamma+eight*Energy_num
  ! -step 
  GAMMAS_num(igamma)=GAMMAS_num(igamma)-two*step
  call calc_E_occ(RDMd,GAMMAs_num,Energy_num,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr)
  grad_igamma=grad_igamma-eight*Energy_num
  ! -2step 
  GAMMAS_num(igamma)=GAMMAS_num(igamma)-step
  call calc_E_occ(RDMd,GAMMAs_num,Energy_num,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr)
  grad_igamma=grad_igamma+Energy_num
  ! Save the gradient
  Grad(igamma)=grad_igamma/(twelve*step)
 enddo
 deallocate(GAMMAs_num) 

end subroutine num_calc_Grad_occ
!!***

!!***
!!****f* DoNOF/calc_Chem_pot
!! NAME
!!  calc_Chem_pot
!!
!! FUNCTION
!!  Calculate the chemical potential for each orbital
!!
!! INPUTS
!!  hCORE=One-body integrals (h_pq) we only use the diag part (h_pp)
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  ERI_L=Lower triangular part of the L_pq matrix
!!  ERI_Jsr=Lower triangular part of the Jsr_pq matrix (zero for rs-NOFT=false)
!!
!! OUTPUT
!!  RDMd%chempot_orb=conatins chem. pot. for each orbital
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_Chem_pot(RDMd,hCORE,ERI_J,ERI_K,ERI_L,ERI_Jsr) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K,ERI_L,ERI_Jsr 
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCORE
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,ipair,igamma=1
!arrays
!************************************************************************
 
 RDMd%chempot_orb=zero

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    RDMd%chempot_orb(iorb) = hCORE(iorb,iorb)                                        &
   &         + half * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)                    &
   &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)          &
   &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)          &
   &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) 
    RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                  & ! rs-NOFT
   &     + half * Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)  &
   &     + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr) 
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
    RDMd%chempot_orb(iorb) = hCORE(iorb,iorb)                                        &
   &         + half * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)                    &
   &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)          &
   &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)          &
   &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) 
    RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                  & ! rs-NOFT
   &     + half * Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr)  &
   &     + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr) 
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                 ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    RDMd%chempot_orb(iorb) = hCORE(iorb,iorb)                                        &
    &        + half * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)                    &
    &        + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)          &
    &        + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)          &
    &        + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) 
    RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                  & ! rs-NOFT
    &     + half * Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr) &
    &     + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr) 
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
     RDMd%chempot_orb(iorb) = hCORE(iorb,iorb)                                       &
    &         + half * RDMd%Dfni_ni(iorb) * ERI_J(iorb*(iorb+1)/2)                   &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERI_J)         &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERI_K)         &
    &         + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERI_L) 
     RDMd%chempot_orb(iorb) = RDMd%chempot_orb(iorb)                                 & ! rs-NOFT
    &     + half * Ddm2_gamma_x_ERI(RDMd,-1,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr) &
    &     + Ddm2_gamma_x_ERI(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_Jsr,ERI_Jsr) 
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --     
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_Chem_pot
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
E_dm2ERI_iorb = zero
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
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*ERI(iorb+iorb1*(iorb1-1)/2)
  enddo
!-----------------------------------------------------------------------
end select

end function dm2_x_eri
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
Grad_Ddm2_ERI_iorb = zero
select case(icase)
!-----------------------------------------------------------------------
 case(-1)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI. 
!-----------------------------------------------------------------------
  Grad_Ddm2_ERI_iorb = DDM2_JorKorL(iorb,iorb,igamma)*ERI(iorb*(iorb+1)/2)
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
end select

end function Ddm2_gamma_x_ERI
!!***

end module m_e_grad_occ 
!!***
