!!****m* DoNOF/m_gammatodm2
!! NAME
!!  m_gammatodm2
!!
!! FUNCTION
!!  Module prepared to perform all procedures required for building the DM2 and their derivatives w.r.t. using gamma
!!
!!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  Nfrozen |            Npairs_p_sing     |              Nvirt                               = NBF               
!!  Nfrozen |         Npairs + Nsingleocc  |     Ncoupled*Npairs                   + Nempty   = NBF               
!!                           | Nsingleocc  |   NBF_occ - Npairs_p_sing - Nfrozen   | Nempty   = NBF
!!                           Nbeta         Nalpha                                  NBF_occ
!!- - - - - - - - - - - - - - -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!!
!! PARENTS
!!  m_e_grad_occ
!!
!! CHILDREN
!!
!! SOURCE

module m_gammatodm2

 use m_rdmd

 implicit none

 private :: dm2_hf,dm2_mbb,dm2_power,dm2_pnof5,dm2_pnof7
!!***

 public :: gamma_to_2rdm
!!***

contains
!!***

!!****f* DoNOF/gamma_to_2rdm
!! NAME
!! gamma_to_2rdm
!!
!! FUNCTION
!!  Build the occ numbers and its derivatives from gamma_i (i.e. the independent variables used in the unconstrained optmization)
!!  Then, call the construction of the 2-RDM element in DM2_J and DM2_K & its derivatives w.r.t occ numbers
!!
!! INPUTS
!!
!! OUTPUT
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  GAMMAs=Independent variables used in the unconstrained optimization as cos^2 (gamma) + sin^2 (gamma) = 1
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine gamma_to_2rdm(RDMd,GAMMAs)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%Ngammas),intent(in)::GAMMAs
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5,iorb6,iorb7,iorb8
 integer::igamma,igamma1,igamma2
 integer::mult
 real(dp)::occ_orb,hole_orb,sqrt_occ_orb,sqrt_hole_orb,sqrthole_orb
!arrays
 real(dp),allocatable,dimension(:)::Docc_gamma0,sqrt_occ,Dsqrt_occ_gamma0,hole
 real(dp),allocatable,dimension(:,:)::Dsqrt_occ_gamma,Dhole_gamma,Docc_gamma 
!************************************************************************

!-----------------------------------------------------------------------
!                 Occupancies and their Derivatives
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(Docc_gamma0(RDMd%NBF_occ),sqrt_occ(RDMd%NBF_occ),Dsqrt_occ_gamma0(RDMd%NBF_occ))
 allocate(Dsqrt_occ_gamma(RDMd%NBF_occ,RDMd%Ngammas),Docc_gamma(RDMd%NBF_occ,RDMd%Ngammas))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 RDMd%occ = 0.0d0; sqrt_occ = 0.0d0; Docc_gamma0 = 0.0d0; Dsqrt_occ_gamma0 = 0.0d0;
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (1,RDMd%Nfrozen)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(RDMd%Nfrozen>0) then
  RDMd%occ(1:RDMd%Nfrozen) = 1.0d0
  sqrt_occ(1:RDMd%Nfrozen) = 1.0d0
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (RDMd%Nfrozen+1,RDMd%Nfrozen+RDMd%Npairs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 do igamma=1,RDMd%Npairs
  iorb1 = RDMd%Nfrozen+igamma
  RDMd%occ(iorb1) = dcos(GAMMAs(igamma))
  RDMd%occ(iorb1) = 0.5d0 + 0.5d0*RDMd%occ(iorb1)*RDMd%occ(iorb1)
  Docc_gamma0(iorb1) = -0.5d0*dsin(2.0d0*GAMMAs(igamma))
  sqrt_occ(iorb1) = dsqrt(RDMd%occ(iorb1))
  Dsqrt_occ_gamma0(iorb1) = 0.5d0*Docc_gamma0(iorb1)/sqrt_occ(iorb1)
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (RDMd%Nfrozen+RDMd%Npairs+1,Nalpha=RDMd%Nfrozen+RDMd%Npairs_p_sing)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(RDMd%Nsingleocc>0) then
  do iorb=RDMd%Npairs+1,RDMd%Npairs_p_sing
   iorb1 = RDMd%Nfrozen+iorb
   mult = RDMd%Nalpha_elect - RDMd%Nbeta_elect
   RDMd%occ(iorb1)  = 0.5d0*mult/RDMd%Nsingleocc
   Docc_gamma0(iorb1) = 0.0d0
   sqrt_occ(iorb1) = dsqrt(RDMd%occ(iorb1))
   Dsqrt_occ_gamma0(iorb1) = 0.0d0
  enddo
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (Nalpha+1,RDMd%NBF_occ)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(RDMd%Ncoupled==1) then            ! PNOFi(1): Perfect Pairing (RDMd%Ncoupled=1)
  Docc_gamma = 0.0d0; Dsqrt_occ_gamma = 0.0d0;
  do igamma=1,RDMd%Npairs
   iorb = RDMd%Nfrozen+igamma         ! iorb=RDMd%Nfrozen+1,RDMd%Nbeta_elect 
   Docc_gamma(iorb,igamma) = Docc_gamma0(iorb)
   Dsqrt_occ_gamma(iorb,igamma) = Dsqrt_occ_gamma0(iorb)
   iorb1 = RDMd%Nalpha_elect+RDMd%Npairs-igamma+1  ! iorb1=RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-igamma)+RDMd%Ncoupled with RDMd%Ncoupled=1
   sqrt_occ(iorb1) = dsin(GAMMAs(igamma))*dsqrt(0.5d0)
   Dsqrt_occ_gamma0(iorb1) = dcos(GAMMAs(igamma))*dsqrt(0.5d0)
   Dsqrt_occ_gamma(iorb1,igamma) = Dsqrt_occ_gamma0(iorb1)
   RDMd%occ(iorb1) = sqrt_occ(iorb1)*sqrt_occ(iorb1)
   Docc_gamma0(iorb1) = 0.5d0*dsin(2.0d0*GAMMAs(igamma))
   Docc_gamma(iorb1,igamma) = Docc_gamma0(iorb1)
  enddo
 else                            ! PNOFi(Nc): Extended PNOF (RDMd%Ncoupled>1)
  allocate(hole(RDMd%Ngammas-RDMd%Npairs),Dhole_gamma(RDMd%Ngammas-RDMd%Npairs,RDMd%Ngammas))
  Docc_gamma = 0.0d0;  Dsqrt_occ_gamma = 0.0d0;  hole = 0.0d0;  Dhole_gamma = 0.0d0;
  do igamma1=1,RDMd%Npairs                      ! igamma=igamma1=1,RDMd%Npairs
   iorb = RDMd%Nfrozen+igamma1                  ! iorb=RDMd%Nfrozen+1,RDMd%Nbeta_elect
   Docc_gamma(iorb,igamma1) = Docc_gamma0(iorb)
   Dsqrt_occ_gamma(iorb,igamma1) = Dsqrt_occ_gamma0(iorb)
   iorb1 = (RDMd%Ncoupled-1)*(igamma1-1)+1
   iorb2 = (RDMd%Ncoupled-1)*igamma1
   hole(iorb1:iorb2)  = 1.0d0 - RDMd%occ(iorb)
   Dhole_gamma(iorb1:iorb2,igamma1) = - Docc_gamma0(iorb)
!- - - -- - - - - - - - - - (igamma1,iorb3) <-> iorb4,igamma,iorb  - - - - - - - - - - - -
   do iorb3=1,RDMd%Ncoupled-1
    iorb4 = (RDMd%Ncoupled-1)*(igamma1-1)+iorb3                           ! iorb4=1,RDMd%Npairs*(RDMd%Ncoupled-1)
    igamma = RDMd%Npairs+iorb4                                            ! igamma=RDMd%Npairs+1,RDMd%Npairs*RDMd%Ncoupled
    iorb5 = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-igamma1)+iorb3   ! iorb5=RDMd%Nalpha_elect+1,RDMd%Nalpha_elect+RDMd%Ncoupled*RDMd%Npairs-1
    sqrt_occ_orb = dsin(GAMMAs(igamma))
    occ_orb = sqrt_occ_orb*sqrt_occ_orb
    Docc_gamma0(iorb5) = dsin(2.0d0*GAMMAs(igamma))
    Dsqrt_occ_gamma0(iorb5) = dcos(GAMMAs(igamma))
    RDMd%occ(iorb5) =  hole(iorb4)*occ_orb
    sqrt_hole_orb = dsqrt(hole(iorb4))
    sqrt_occ(iorb5) = sqrt_hole_orb*sqrt_occ_orb
    Docc_gamma(iorb5,igamma1) = Dhole_gamma(iorb4,igamma1)*occ_orb
    if(sqrt_hole_orb>0.0d0) then
     Dsqrt_occ_gamma(iorb5,igamma1) = 0.5d0*Dhole_gamma(iorb4,igamma1)*sqrt_occ_orb/sqrt_hole_orb
    else
     Dsqrt_occ_gamma(iorb5,igamma1) = 0.0d0
    endif
    do iorb6=iorb1,iorb4-1                   !     iorb1 < iorb6   < iorb4-1
     igamma2 = RDMd%Npairs+iorb6             !   igamma1 < igamma2 < igamma
     Docc_gamma(iorb5,igamma2) =  Dhole_gamma(iorb4,igamma2)*occ_orb
     if(sqrt_hole_orb>0.0d0) then
      Dsqrt_occ_gamma(iorb5,igamma2) = 0.5d0*Dhole_gamma(iorb4,igamma2)*sqrt_occ_orb/sqrt_hole_orb
     else
      Dsqrt_occ_gamma(iorb5,igamma2) = 0.0d0
     endif
    enddo
    Docc_gamma(iorb5,igamma) = hole(iorb4)*Docc_gamma0(iorb5)
    Dsqrt_occ_gamma(iorb5,igamma) = sqrt_hole_orb*Dsqrt_occ_gamma0(iorb5)
!- - - - hole(iorb4+1) - - - - - - - - - - - - - - - -
    if(iorb3<RDMd%Ncoupled-1) then
     do iorb7=1,RDMd%Ncoupled-1-iorb3
      iorb6 = iorb4+iorb7                          ! iorb4 < iorb6 < igamma1*(RDMd%Ncoupled-1)
      hole(iorb6) = hole(iorb6) - RDMd%occ(iorb5)
      Dhole_gamma(iorb6,igamma1) = Dhole_gamma(iorb6,igamma1) - Docc_gamma(iorb5,igamma1)
      do iorb8=iorb1,iorb4-1
       igamma2 = RDMd%Npairs+iorb8
       Dhole_gamma(iorb6,igamma2) = Dhole_gamma(iorb6,igamma2) - Docc_gamma(iorb5,igamma2)
      enddo
      Dhole_gamma(iorb6,igamma) = Dhole_gamma(iorb6,igamma) - Docc_gamma(iorb5,igamma)
     enddo
    endif
!- - - - hole(iorb4+1) - - - - - - - - - - - - - - - -
   enddo
!- - - -- - - - - - - - - - (igamma1,iorb3) <-> iorb4,igamma,iorb  - - - - - - - - - - - -
!- - - - iorb4 = iorb2 - last occ  - - - - - - - - - - - - - -
   igamma = RDMd%Npairs+iorb2               ! igamma=RDMd%Npairs+igamma1*(RDMd%Ncoupled-1)
   iorb5 = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-igamma1)+RDMd%Ncoupled
   sqrt_occ_orb = dcos(GAMMAs(igamma))
   hole_orb = sqrt_occ_orb*sqrt_occ_orb 
   Docc_gamma0(iorb5)  = -dsin(2.0d0*GAMMAs(igamma))
   Dsqrt_occ_gamma0(iorb5) = -dsin(GAMMAs(igamma))
   RDMd%occ(iorb5) = hole(iorb2)*hole_orb
   sqrthole_orb = dsqrt(hole(iorb2))
   sqrt_occ(iorb5)= sqrthole_orb*sqrt_occ_orb
   Docc_gamma(iorb5,igamma1) = Dhole_gamma(iorb2,igamma1)*hole_orb
   if(sqrthole_orb>0.0d0) then
    Dsqrt_occ_gamma(iorb5,igamma1) = 0.5d0*Dhole_gamma(iorb2,igamma1)*sqrt_occ_orb/sqrthole_orb
   else
    Dsqrt_occ_gamma(iorb5,igamma1) = 0.0d0
   endif
   do iorb6=iorb1,iorb2-1            !     iorb1 < iorb6   < iorb2-1
    igamma2 = RDMd%Npairs+iorb6      !   igamma1 < igamma2 < igamma
    Docc_gamma(iorb5,igamma2) = Dhole_gamma(iorb2,igamma2)*hole_orb
    if(sqrthole_orb>0.0d0) then
     Dsqrt_occ_gamma(iorb5,igamma2) = 0.5d0*Dhole_gamma(iorb2,igamma2)*sqrt_occ_orb/sqrthole_orb
    else
     Dsqrt_occ_gamma(iorb5,igamma2) = 0.0d0
    endif
   enddo
   Docc_gamma(iorb5,igamma) = hole(iorb2) *Docc_gamma0(iorb5)
   Dsqrt_occ_gamma(iorb5,igamma) = sqrthole_orb*Dsqrt_occ_gamma0(iorb5)
!- - - - iorb4 = iorb2 - last occ  - - - - - - - - - - - - - -
  enddo
  deallocate(hole,Dhole_gamma)
 endif
 deallocate(Docc_gamma0,Dsqrt_occ_gamma0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Sum of the Holes below the Fermi Level (RDMd%Sums)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 RDMd%Sums = DFLOAT(RDMd%Nbeta_elect)
 do iorb=1,RDMd%Nbeta_elect
  RDMd%Sums = RDMd%Sums - RDMd%occ(iorb)
 enddo
!-----------------------------------------------------------------------
!                   DM2_J, DM2_K, DDM2_gamma_J, DDM2_gamma_K
!Comment: This is not the cleanest way to call them but it is fastest way 
!         to program them without requiring extra memory allocations or pointers.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 RDMd%Docc_gamma=reshape(Docc_gamma,(/RDMd%NBF_occ*RDMd%Ngammas/))
 if(RDMd%INOF==0) then
  call dm2_hf(RDMd,RDMd%Docc_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==-1) then
  call dm2_mbb(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==-2) then
  call dm2_power(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==5) then
  call dm2_pnof5(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==7) then
  call dm2_pnof7(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 else
  ! Nth
 endif
!-----------------------------------------------------------------------
 deallocate(sqrt_occ,Dsqrt_occ_gamma,Docc_gamma)
 
end subroutine gamma_to_2rdm

!!***
!!****f* DoNOF/dm2_hf
!! NAME
!! dm2_hf
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for HF
!!  Note: For spin-uncompensated this is a weird HF. Warning!
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
!! Dsqrt_occ_gamma=Matrix with the derivative of sqrt(occ numbers) vs gamma
!!
!! OUTPUT
!! DM2_iiii=DM2 same orb elements
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DM2_K=DM2 elements that use L integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals 
!! DDM2_gamma_K=Derivative of the DM2 elements w.r.t. gamma that use K integrals
!! DDM2_gamma_L=Derivative of the DM2 elements w.r.t. gamma that use L integrals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_hf(RDMd,Docc_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Docc_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,igamma
!arrays
!************************************************************************

 DM2_L=0.0d0; DDM2_gamma_L=0.0d0; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -NpNq [ DDM2_Jpqk = 2DNpk*Nq, DDM2_Kpq = -DNpk*Nq ]
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = 2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1)
   do igamma=1,RDMd%Ngammas
    DDM2_gamma_J(iorb,iorb1,igamma) = 2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    DDM2_gamma_K(iorb,iorb1,igamma) = -Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
   enddo
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
 if(RDMd%Nsingleocc>1) then
  do iorb=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
    DM2_K(iorb,iorb1) = -2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
    do igamma=1,RDMd%Ngammas
     DDM2_gamma_K(iorb,iorb1,igamma) = -2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    enddo
   enddo
  enddo
 end if
!-----------------------------------------------------------------------
!                 DM2(iorb,iorb,iorb,iorb)=OCC(iorb)*OCC(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)*RDMd%occ(iorb)
  DM2_J(iorb,iorb)=0.0d0
  DM2_K(iorb,iorb)=0.0d0
  RDMd%Dfni_ni(iorb)=2.0d0*RDMd%occ(iorb)
  DDM2_gamma_J(iorb,iorb,:)=0.0d0
  DDM2_gamma_K(iorb,iorb,:)=0.0d0
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_hf
!!***

!!***
!!****f* DoNOF/dm2_mbb
!! NAME
!! dm2_mbb
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for MBB
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
!! Dsqrt_occ_gamma=Matrix with the derivative of sqrt(occ numbers) vs gamma
!!
!! OUTPUT
!! DM2_iiii=DM2 same orb elements
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DM2_L=DM2 elements that use L integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals 
!! DDM2_gamma_K=Derivative of the DM2 elements w.r.t. gamma that use K integrals
!! DDM2_gamma_L=Derivative of the DM2 elements w.r.t. gamma that use L integrals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_mbb(RDMd,Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
& DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Dsqrt_occ_gamma,Docc_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,igamma
!arrays
!************************************************************************

 DM2_L=0.0d0; DDM2_gamma_L=0.0d0; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -NpNq [ DDM2_Jpqk = 2DNpk*Nq, DDM2_Kpq = -DNpk*Nq ] 
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = 2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -sqrt_occ(iorb)*sqrt_occ(iorb1)
   do igamma=1,RDMd%Ngammas
    DDM2_gamma_J(iorb,iorb1,igamma) = 2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    DDM2_gamma_K(iorb,iorb1,igamma) = -Dsqrt_occ_gamma(iorb,igamma)*sqrt_occ(iorb1)
   enddo
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
 if(RDMd%Nsingleocc>1) then
  do iorb=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
    DM2_K(iorb,iorb1) = -2.0d0*sqrt_occ(iorb)*sqrt_occ(iorb1)
    do igamma=1,RDMd%Ngammas
     DDM2_gamma_K(iorb,iorb1,igamma) = -2.0d0*Dsqrt_occ_gamma(iorb,igamma)*sqrt_occ(iorb1)
    enddo
   enddo
  enddo
 end if
!-----------------------------------------------------------------------
!          DM2(iorb,iorb,iorb,iorb)=2*OCC(iorb)*OCC(iorb)-OCC(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb)-RDMd%occ(iorb)
  DM2_J(iorb,iorb)=0.0d0
  DM2_K(iorb,iorb)=0.0d0
  RDMd%Dfni_ni(iorb)=4.0d0*RDMd%occ(iorb)-1.0d0
  DDM2_gamma_J(iorb,iorb,:)=0.0d0
  DDM2_gamma_K(iorb,iorb,:)=0.0d0
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_mbb
!!***

!!***
!!****f* DoNOF/dm2_power
!! NAME
!! dm2_power
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for POWER 
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
!! Dsqrt_occ_gamma=Matrix with the derivative of sqrt(occ numbers) vs gamma
!!
!! OUTPUT
!! DM2_iiii=DM2 same orb elements
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DM2_L=DM2 elements that use L integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals 
!! DDM2_gamma_K=Derivative of the DM2 elements w.r.t. gamma that use K integrals
!! DDM2_gamma_L=Derivative of the DM2 elements w.r.t. gamma that use L integrals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_power(RDMd,Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
& DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Dsqrt_occ_gamma,Docc_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,igamma
 real(dp)::ONEmLpower
!arrays
!************************************************************************

 ONEmLpower=1.0d0-RDMd%Lpower
 DM2_L=0.0d0; DDM2_gamma_L=0.0d0; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -NpNq [ DDM2_Jpqk = 2DNpk*Nq, DDM2_Kpq = -DNpk*Nq ] 
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = 2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -((RDMd%occ(iorb)*RDMd%occ(iorb1))**RDMd%Lpower)
   do igamma=1,RDMd%Ngammas
    DDM2_gamma_J(iorb,iorb1,igamma) = 2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    DDM2_gamma_K(iorb,iorb1,igamma) = -RDMd%Lpower*RDMd%occ(iorb1)*Docc_gamma(iorb,igamma)&
&                                   *((RDMd%occ(iorb)*RDMd%occ(iorb1))**ONEmLpower)
   enddo
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
 if(RDMd%Nsingleocc>1) then
  do iorb=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
    DM2_K(iorb,iorb1) = -2.0d0*sqrt_occ(iorb)*sqrt_occ(iorb1)
    do igamma=1,RDMd%Ngammas
     DDM2_gamma_K(iorb,iorb1,igamma) = -2.0d0*Dsqrt_occ_gamma(iorb,igamma)*sqrt_occ(iorb1)
    enddo
   enddo
  enddo
 end if
!-----------------------------------------------------------------------
!          DM2(iorb,iorb,iorb,iorb)=2*OCC(iorb)*OCC(iorb)-OCC(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb)-(RDMd%occ(iorb)**(2.0d0*RDMd%Lpower))
  DM2_J(iorb,iorb)=0.0d0
  DM2_K(iorb,iorb)=0.0d0
  RDMd%Dfni_ni(iorb)=4.0d0*RDMd%occ(iorb)-2.0d0*RDMd%Lpower*(RDMd%occ(iorb)**(2.0d0*RDMd%Lpower-1.0d0))
  DDM2_gamma_J(iorb,iorb,:)=0.0d0
  DDM2_gamma_K(iorb,iorb,:)=0.0d0
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_power
!!***

!!****f* DoNOF/dm2_pnof5
!! NAME
!! dm2_pnof5
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for PNOF5
!!  JCP 134, 164102, 2011; JCP 139, 234109, 2013
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
!! Dsqrt_occ_gamma=Matrix with the derivative of sqrt(occ numbers) vs gamma
!!
!! OUTPUT
!! DM2_iiii=DM2 same orb elements
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DM2_L=DM2 elements that use L integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals 
!! DDM2_gamma_K=Derivative of the DM2 elements w.r.t. gamma that use K integrals
!! DDM2_gamma_L=Derivative of the DM2 elements w.r.t. gamma that use L integrals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_pnof5(RDMd,Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
& DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Dsqrt_occ_gamma,Docc_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5,igamma
!arrays
!************************************************************************

!-----------------------------------------------------------------------
!                Inter-pair interactions for PNOF5 (Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     DM2_Jpq = 2NpNq, DM2_Kpq = -NpNq [ DDM2_Jpqk = 2DNpk*Nq, DDM2_Kpq = -DNpk*Nq ] 
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = 2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_L(iorb,iorb1) = 0.0d0
   do igamma=1,RDMd%Ngammas
    DDM2_gamma_J(iorb,iorb1,igamma) = 2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    DDM2_gamma_K(iorb,iorb1,igamma) = -Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    DDM2_gamma_L(iorb,iorb1,igamma) = 0.0d0 
   enddo
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
 if(RDMd%Nsingleocc>1) then ! TODO
  do iorb=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
    DM2_K(iorb,iorb1) = -2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
    do igamma=1,RDMd%Ngammas
     DDM2_gamma_K(iorb,iorb1,igamma) = -2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    enddo
   enddo
  enddo
 end if
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for PNOF5(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb2=1,RDMd%Npairs
  iorb3 = RDMd%Nfrozen+iorb2
  do iorb1=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
   iorb4 = RDMd%Nfrozen+iorb1
   DM2_J(iorb3,iorb4) = 0.0d0
   DM2_J(iorb4,iorb3) = 0.0d0
   DM2_K(iorb3,iorb4) = 0.0d0
   DM2_K(iorb4,iorb3) = 0.0d0
   DM2_L(iorb3,iorb4) = -sqrt_occ(iorb3)*sqrt_occ(iorb4)
   DM2_L(iorb4,iorb3) = -sqrt_occ(iorb4)*sqrt_occ(iorb3)
   do igamma=1,RDMd%Ngammas
    DDM2_gamma_J(iorb3,iorb4,igamma) = 0.0d0
    DDM2_gamma_J(iorb4,iorb3,igamma) = 0.0d0
    DDM2_gamma_K(iorb3,iorb4,igamma) = 0.0d0
    DDM2_gamma_K(iorb4,iorb3,igamma) = 0.0d0
    DDM2_gamma_L(iorb3,iorb4,igamma) = -Dsqrt_occ_gamma(iorb3,igamma)*sqrt_occ(iorb4)
    DDM2_gamma_L(iorb4,iorb3,igamma) = -Dsqrt_occ_gamma(iorb4,igamma)*sqrt_occ(iorb3)
   enddo
   do iorb=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
    iorb5 = RDMd%Nfrozen+iorb
    DM2_J(iorb5,iorb4) = 0.0d0
    DM2_K(iorb5,iorb4) = 0.0d0
    DM2_L(iorb5,iorb4) = sqrt_occ(iorb5)*sqrt_occ(iorb4)
    do igamma=1,RDMd%Ngammas
     DDM2_gamma_J(iorb5,iorb4,igamma) = 0.0d0
     DDM2_gamma_K(iorb5,iorb4,igamma) = 0.0d0
     DDM2_gamma_L(iorb5,iorb4,igamma) = Dsqrt_occ_gamma(iorb5,igamma)*sqrt_occ(iorb4)
    enddo
   enddo
  enddo
 enddo
!-----------------------------------------------------------------------
!                 DM2(iorb,iorb,iorb,iorb)=OCC(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=0.0d0
  DM2_K(iorb,iorb)=0.0d0
  DM2_L(iorb,iorb)=0.0d0
  RDMd%Dfni_ni(iorb)=1.0d0
  DDM2_gamma_J(iorb,iorb,:)=0.0d0
  DDM2_gamma_K(iorb,iorb,:)=0.0d0
  DDM2_gamma_L(iorb,iorb,:)=0.0d0
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_pnof5
!!***

!!****f* DoNOF/dm2_pnof7
!! NAME
!! dm2_pnof7
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for PNOF7
!!  PRL 119, 063002, 2017; PRA 100, 032508, 2019
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
!! Dsqrt_occ_gamma=Matrix with the derivative of sqrt(occ numbers) vs gamma
!!
!! OUTPUT
!! DM2_iiii=DM2 same orb elements
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DM2_L=DM2 elements that use L integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals 
!! DDM2_gamma_K=Derivative of the DM2 elements w.r.t. gamma that use K integrals
!! DDM2_gamma_L=Derivative of the DM2 elements w.r.t. gamma that use L integrals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_pnof7(RDMd,Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
& DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Dsqrt_occ_gamma,Docc_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5,igamma
!arrays
 real(dp),allocatable,dimension(:)::FIs
 real(dp),allocatable,dimension(:,:)::DFIs
!************************************************************************

!-----------------------------------------------------------------------
!     FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(FIs(RDMd%NBF_occ),DFIs(RDMd%NBF_occ,RDMd%Ngammas))
 FIs = 0.0d0; DFIs = 0.0d0;
 if(RDMd%Ista==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (Np*Hp)^1/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = dsqrt( RDMd%occ(iorb)*(1.0d0-RDMd%occ(iorb)) )
   if(FIs(iorb)>1.0d-20) then
    do igamma=1,RDMd%Ngammas
     DFIs(iorb,igamma) = (0.5d0-RDMd%occ(iorb))*Docc_gamma(iorb,igamma)/FIs(iorb)
    enddo
   endif
  enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
 else if(RDMd%Ista==1) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (4*Np*Hp)^0.5*(Np*Hp)^0.5 = 2*Np*Hp
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = 2.0d0*RDMd%occ(iorb)*(1.0d0-RDMd%occ(iorb))
   do igamma=1,RDMd%Ngammas
    DFIs(iorb,igamma) = 2.0d0*(1.0d0-2.0d0*RDMd%occ(iorb))*Docc_gamma(iorb,igamma)
   enddo
  enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - -       
 end if
!-----------------------------------------------------------------------
!                Inter-pair interactions for PNOF7 (Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = 2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1) 
   DM2_L(iorb,iorb1) = -FIs(iorb)*FIs(iorb1)
   do igamma=1,RDMd%Ngammas
    DDM2_gamma_J(iorb,iorb1,igamma) = 2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    DDM2_gamma_K(iorb,iorb1,igamma) = -Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    DDM2_gamma_L(iorb,iorb1,igamma) = -DFIs(iorb,igamma)*FIs(iorb1)
   enddo
  enddo
 enddo
 deallocate(FIs,DFIs)
!- - - - - - - - - - - - - - - - - - - - - - - -              
 if(RDMd%Nsingleocc>1) then ! TODO
  do iorb=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
   do iorb1=RDMd%Nbeta_elect+1,RDMd%Nalpha_elect
    DM2_K(iorb,iorb1) = -2.0d0*RDMd%occ(iorb)*RDMd%occ(iorb1)
    do igamma=1,RDMd%Ngammas
     DDM2_gamma_K(iorb,iorb1,igamma) = -2.0d0*Docc_gamma(iorb,igamma)*RDMd%occ(iorb1)
    enddo
   enddo
  enddo
 end if
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for PNOF7(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb2=1,RDMd%Npairs
  iorb3 = RDMd%Nfrozen+iorb2
  do iorb1=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
   iorb4 = RDMd%Nfrozen+iorb1
   DM2_J(iorb3,iorb4) = 0.0d0
   DM2_J(iorb4,iorb3) = 0.0d0
   DM2_K(iorb3,iorb4) = 0.0d0
   DM2_K(iorb4,iorb3) = 0.0d0
   DM2_L(iorb3,iorb4) = -sqrt_occ(iorb3)*sqrt_occ(iorb4)
   DM2_L(iorb4,iorb3) = -sqrt_occ(iorb4)*sqrt_occ(iorb3)
   do igamma=1,RDMd%Ngammas
    DDM2_gamma_J(iorb3,iorb4,igamma) = 0.0d0
    DDM2_gamma_J(iorb4,iorb3,igamma) = 0.0d0
    DDM2_gamma_K(iorb3,iorb4,igamma) = 0.0d0
    DDM2_gamma_K(iorb4,iorb3,igamma) = 0.0d0
    DDM2_gamma_L(iorb3,iorb4,igamma) = -Dsqrt_occ_gamma(iorb3,igamma)*sqrt_occ(iorb4)
    DDM2_gamma_L(iorb4,iorb3,igamma) = -Dsqrt_occ_gamma(iorb4,igamma)*sqrt_occ(iorb3)
   enddo
   do iorb=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
    iorb5 = RDMd%Nfrozen+iorb
    DM2_J(iorb5,iorb4) = 0.0d0
    DM2_K(iorb5,iorb4) = 0.0d0
    DM2_L(iorb5,iorb4) = sqrt_occ(iorb5)*sqrt_occ(iorb4)
    do igamma=1,RDMd%Ngammas
     DDM2_gamma_J(iorb5,iorb4,igamma) = 0.0d0
     DDM2_gamma_K(iorb5,iorb4,igamma) = 0.0d0
     DDM2_gamma_L(iorb5,iorb4,igamma) = Dsqrt_occ_gamma(iorb5,igamma)*sqrt_occ(iorb4)
    enddo
   enddo
  enddo
 enddo
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                 DM2(iorb,iorb,iorb,iorb)=OCC(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=0.0d0
  DM2_K(iorb,iorb)=0.0d0
  DM2_L(iorb,iorb)=0.0d0
  RDMd%Dfni_ni(iorb)=1.0d0
  DDM2_gamma_J(iorb,iorb,:)=0.0d0
  DDM2_gamma_K(iorb,iorb,:)=0.0d0
  DDM2_gamma_L(iorb,iorb,:)=0.0d0
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_pnof7
!!***

end module m_gammatodm2
!!***
