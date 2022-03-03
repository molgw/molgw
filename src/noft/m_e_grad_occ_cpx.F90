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

 private :: dm2_x_eriC,Ddm2_gamma_x_ERIc
!!***

 public :: calc_E_occ_cpx,calc_Grad_occ_cpx,calc_Chem_pot_cpx
!!***

contains

!!***
!!****f* DoNOF/calc_E_occ_cpx
!! NAME
!!  calc_E_occ_cpx
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

subroutine calc_E_occ_cpx(RDMd,GAMMAs,Energy,hCOREc,ERIc_J,ERIc_K,ERIc_L,nogamma,chempot) 
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::nogamma,chempot
 real(dp),intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%Ngammas),intent(in)::GAMMAs
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERIc_J,ERIc_K,ERIc_L 
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCOREc
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,ipair
 complex(dp)::Energy_cpx
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

 Energy=zero; Energy_cpx=COMPLEX_ZERO;

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!      Singlet State (S=0,Ms=0) and Multiplet States (S>0,Ms=0)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if(RDMd%Ncoupled==1) then     ! Perfect Pairing (Ncoupled=1)
   do iorb=1,RDMd%Nfrozen
    Energy_cpx = Energy_cpx + RDMd%occ(iorb) * two*hCOREc(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERIc_J(iorb*(iorb+1)/2)                                            &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_J,ERIc_J) + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_K,ERIc_K)    &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_L,ERIc_L)
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy_cpx = Energy_cpx + RDMd%occ(iorb) * two*hCOREc(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERIc_J(iorb*(iorb+1)/2)                                            &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_J,ERIc_J) + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_K,ERIc_K)    &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_L,ERIc_L)
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
    Energy_cpx = Energy_cpx + RDMd%occ(iorb) * two*hCOREc(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERIc_J(iorb*(iorb+1)/2)                                            &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_J,ERIc_J) + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_K,ERIc_K)    &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_L,ERIc_L)
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                  ! Extended PNOF (Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do iorb=1,RDMd%Nfrozen
    Energy_cpx = Energy_cpx + RDMd%occ(iorb) * two*hCOREc(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERIc_J(iorb*(iorb+1)/2)                                            &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_J,ERIc_J) + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_K,ERIc_K)    &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_L,ERIc_L)
   enddo
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
    Energy_cpx = Energy_cpx + RDMd%occ(iorb) * two*hCOREc(iorb,iorb)                                  &
    &      + RDMd%DM2_iiii(iorb) * ERIc_J(iorb*(iorb+1)/2)                                            &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_J,ERIc_J) + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_K,ERIc_K)    &
    &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_L,ERIc_L)
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
     Energy_cpx = Energy_cpx + RDMd%occ(iorb) * two*hCOREc(iorb,iorb)                                 &
     &      + RDMd%DM2_iiii(iorb) * ERIc_J(iorb*(iorb+1)/2)                                           &
     &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_J,ERIc_J) + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_K,ERIc_K)   &
     &      + dm2_x_eriC(RDMd,0,iorb,RDMd%DM2_L,ERIc_L)
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 else ! TODO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 Energy=real(Energy_cpx)
 if(dabs(aimag(Energy_cpx))>tol8) then
  write(msg,'(a,f10.5,a)') 'Warning! Imaginary[Energy] = ',aimag(Energy_cpx),' (a.u.)'
  call write_output(msg)
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_E_occ_cpx
!!***

!!***
!!****f* DoNOF/calc_Grad_occ_cpx
!! NAME
!!  calc_Grad_occ_cpx
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

subroutine calc_Grad_occ_cpx(RDMd,Grad,hCOREc,ERIc_J,ERIc_K,ERIc_L) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 real(dp),dimension(RDMd%Ngammas),intent(inout)::Grad
!arrays
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERIc_J,ERIc_K,ERIc_L 
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCOREc
!Local variables ------------------------------
!scalars
 integer::igamma,iorb,iorb1,ipair
 real(dp)::max_imag
!arrays
 character(len=200)::msg
 complex(dp),allocatable,dimension(:)::Grad_cpx
!************************************************************************
 
 allocate(Grad_cpx(RDMd%Ngammas))

 Grad = zero; Grad_cpx=COMPLEX_ZERO;

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad_cpx(igamma) = Grad_cpx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCOREc(iorb,iorb) &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)     &
    &         + two * ( Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)                           &
    &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)                                   &
    &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) )
     iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
     Grad_cpx(igamma) = Grad_cpx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCOREc(iorb,iorb) &
    &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)     &
    &         + two * ( Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)                           &
    &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)                                   &
    &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) )
    enddo
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                 ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do igamma=1,RDMd%Ngammas
    do ipair=1,RDMd%Npairs
     iorb = RDMd%Nfrozen+ipair
     Grad_cpx(igamma) = Grad_cpx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCOREc(iorb,iorb) &
     &        + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)     &
     &        + two * ( Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)                           &
     &        + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)                                   &
     &        + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) )
     do iorb1=1,RDMd%Ncoupled
      iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
      Grad_cpx(igamma) = Grad_cpx(igamma) + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * two*hCOREc(iorb,iorb) &
     &         + RDMd%Docc_gamma(iorb+(igamma-1)*RDMd%NBF_occ) * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)     &
     &         + two * (Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)                            &
     &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)                                   &
     &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) )
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
 Grad(:)=real(Grad_cpx(:))
 max_imag=maxval(aimag(Grad_cpx(:)))
 if(dabs(max_imag)>tol8) then
  write(msg,'(a,f10.5,a)') 'Warning! Imaginary[Grad_Energy] = ',max_imag
  call write_output(msg)
 endif
  
 deallocate(Grad_cpx)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine calc_Grad_occ_cpx
!!***

!!***
!!****f* DoNOF/calc_Chem_pot_cpx
!! NAME
!!  calc_Chem_pot_cpx
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

subroutine calc_Chem_pot_cpx(RDMd,hCOREc,ERIc_J,ERIc_K,ERIc_L) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 complex(dp),dimension(RDMd%NBF_ldiag),intent(in)::ERIc_J,ERIc_K,ERIc_L 
 complex(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::hCOREc
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,ipair,igamma=1
 real(dp)::max_imag
!arrays
 complex(dp),allocatable,dimension(:)::chempot_orb
 character(len=200)::msg
!************************************************************************
 
 allocate(chempot_orb(RDMd%NBF_occ))

 RDMd%chempot_orb=zero; chempot_orb=COMPLEX_ZERO; 

 if(RDMd%Nsingleocc==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(RDMd%Ncoupled==1) then       ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
         chempot_orb(iorb) = hCOREc(iorb,iorb)                                       &
   &         + half * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)                   &
   &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)        &
   &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)        &
   &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) 
    iorb = RDMd%Nalpha_elect+RDMd%Npairs-ipair+1
         chempot_orb(iorb) = hCOREc(iorb,iorb)                                       &
   &         + half * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)                   &
   &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)        &
   &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)        &
   &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) 
   enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else                 ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   do ipair=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+ipair
         chempot_orb(iorb) = hCOREc(iorb,iorb)                                       &
    &        + half * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)                   &
    &        + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)        &
    &        + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)        &
    &        + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) 
    do iorb1=1,RDMd%Ncoupled
     iorb = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-ipair)+iorb1
          chempot_orb(iorb) = hCOREc(iorb,iorb)                                      &
    &         + half * RDMd%Dfni_ni(iorb) * ERIc_J(iorb*(iorb+1)/2)                  &
    &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_J,ERIc_J)       &
    &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_K,ERIc_K)       &
    &         + Ddm2_gamma_x_ERIc(RDMd,0,iorb,igamma,RDMd%DDM2_gamma_L,ERIc_L) 
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

end subroutine calc_Chem_pot_cpx
!!***

!!***
!!****f* DoNOF/dm2_x_eriC
!! NAME
!!  dm2_x_eriC
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

function dm2_x_eriC(RDMd,icase,iorb,DM2_JKL,ERI) result(E_dm2ERI_iorb)
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
E_dm2ERI_iorb = COMPLEX_ZERO
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
   E_dm2ERI_iorb = E_dm2ERI_iorb + DM2_JKL(iorb,iorb1)*conjg(ERI(iorb+iorb1*(iorb1-1)/2))
  enddo
!-----------------------------------------------------------------------
 case(1) ! TODO  
!-----------------------------------------------------------------------
!    DM2_JKL*ERI, iorb<Nbeta_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 case(2) ! TODO 
!-----------------------------------------------------------------------
!    DM2_JKL*ERI, iorb>Nalpha_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end select

end function dm2_x_eriC
!!***

!!***
!!****f* DoNOF/Ddm2_gamma_x_ERIc
!! NAME
!!  Ddm2_gamma_x_ERIc
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

function Ddm2_gamma_x_ERIc(RDMd,icase,iorb,igamma,DDM2_JorKorL,ERI) result(Grad_Ddm2_ERI_iorb)
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
Grad_Ddm2_ERI_iorb = COMPLEX_ZERO
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
   Grad_Ddm2_ERI_iorb = Grad_Ddm2_ERI_iorb + DDM2_JorKorL(iorb,iorb1,igamma)*conjg(ERI(iorb+iorb1*(iorb1-1)/2))
  enddo
!-----------------------------------------------------------------------
 case(1)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI, iorb<Nbeta_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 case(2)
!-----------------------------------------------------------------------
!     DDM2_JorK*ERI, iorb>Nalpha_elect (Nsingleocc is excluded from the Sum)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end select

end function Ddm2_gamma_x_ERIc
!!***

end module m_e_grad_occ_cpx 
!!***
