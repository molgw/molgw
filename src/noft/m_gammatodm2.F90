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

 use m_nofoutput
 use m_rdmd

 implicit none

 private :: dm2_hartree,dm2_hf,dm2_mbb,dm2_ca,dm2_cga,dm2_gu,dm2_power,dm2_pnof5,dm2_pnof7,dm2_gnof
 private :: dm2_intra,dm2_pccd,dm2_pnof7_sup
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
!!  chempot=When looking for chemical potential calc. the d occ(iorb)/d GAMMAs_igamma =1 (i.e. we only want dDM2/docc(iorb))
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine gamma_to_2rdm(RDMd,GAMMAs,chempot)
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::chempot
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%Ngammas),intent(in)::GAMMAs
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5,iorb6,iorb7,iorb8
 integer::igamma,igamma1,igamma2
 real(dp)::occ_orb,hole_orb,sqrt_occ_orb,sqrt_hole_orb,sqrthole_orb
 real(dp)::exponential,dexponent,hole_dyn
!arrays
 real(dp),allocatable,dimension(:)::Docc_gamma0,sqrt_occ,Dsqrt_occ_gamma0,hole
 real(dp),allocatable,dimension(:,:)::Dsqrt_occ_gamma,Dhole_gamma,Docc_gamma,Docc_dyn 
!************************************************************************
!-----------------------------------------------------------------------
!                 Occupancies and their Derivatives
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(Docc_gamma0(RDMd%NBF_occ),sqrt_occ(RDMd%NBF_occ),Dsqrt_occ_gamma0(RDMd%NBF_occ))
 allocate(Dsqrt_occ_gamma(RDMd%NBF_occ,RDMd%Ngammas),Docc_gamma(RDMd%NBF_occ,RDMd%Ngammas))
 allocate(Docc_dyn(RDMd%NBF_occ,RDMd%Ngammas))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 RDMd%occ = zero; RDMd%occ_dyn = zero; Docc_dyn = zero; sqrt_occ = zero; 
 Docc_gamma0 = zero; Dsqrt_occ_gamma0 = zero;
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (1,RDMd%Nfrozen)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(RDMd%Nfrozen>0) then
  RDMd%occ(1:RDMd%Nfrozen) = one
  sqrt_occ(1:RDMd%Nfrozen) = one
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (RDMd%Nfrozen+1,RDMd%Nfrozen+RDMd%Npairs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 do igamma=1,RDMd%Npairs
  iorb = RDMd%Nfrozen+igamma
  RDMd%occ(iorb) = half + half*dcos(GAMMAs(igamma))*dcos(GAMMAs(igamma)) 
  Docc_gamma0(iorb) = -half*dsin(two*GAMMAs(igamma))
  sqrt_occ(iorb) = dsqrt(RDMd%occ(iorb))
  Dsqrt_occ_gamma0(iorb) = half*Docc_gamma0(iorb)/sqrt_occ(iorb)
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (RDMd%Nfrozen+RDMd%Npairs+1,Nalpha=RDMd%Nfrozen+RDMd%Npairs_p_sing)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !if(RDMd%Nsingleocc>0) then
 ! TODO
 !endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Occupancies (Nalpha+1,RDMd%NBF_occ)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(RDMd%Ncoupled==1) then            ! NOF: Perfect Pairing (RDMd%Ncoupled=1)
  Docc_gamma = zero; Dsqrt_occ_gamma = zero;
  do igamma=1,RDMd%Npairs
   ! iorb=RDMd%Nfrozen+1,RDMd%Nbeta_elect 
   iorb = RDMd%Nfrozen+igamma
   Docc_gamma(iorb,igamma) = Docc_gamma0(iorb)
   Dsqrt_occ_gamma(iorb,igamma) = Dsqrt_occ_gamma0(iorb)
   ! iorb=RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-igamma)+RDMd%Ncoupled with RDMd%Ncoupled=1
   iorb2 = RDMd%Nalpha_elect+RDMd%Npairs-igamma+1
   RDMd%occ(iorb2) = one-RDMd%occ(iorb)
   Docc_gamma0(iorb2) = -Docc_gamma0(iorb)
   Docc_gamma(iorb2,igamma)=Docc_gamma0(iorb2)
   sqrt_occ(iorb2)=dsqrt(RDMd%occ(iorb2)) 
   if(sqrt_occ(iorb2)>tol20) then
    Dsqrt_occ_gamma0(iorb2)=half*Docc_gamma0(iorb2)/sqrt_occ(iorb2)
   else
    Dsqrt_occ_gamma0(iorb2)=zero
   endif
   Dsqrt_occ_gamma(iorb2,igamma)=Dsqrt_occ_gamma0(iorb2) 
  enddo
 else                                 ! NOF: Extended NOF (RDMd%Ncoupled>1)
  allocate(hole(RDMd%Ngammas-RDMd%Npairs),Dhole_gamma(RDMd%Ngammas-RDMd%Npairs,RDMd%Ngammas))
  Docc_gamma = zero;  Dsqrt_occ_gamma = zero;  hole = zero;  Dhole_gamma = zero;
   do igamma=1,RDMd%Npairs                           ! igamma1=igamma=1,RDMd%Npairs
    iorb = RDMd%Nfrozen+igamma                                ! iorb=no1+1,nb
    Docc_gamma(iorb,igamma) = Docc_gamma0(iorb)
    Dsqrt_occ_gamma(iorb,igamma) = Dsqrt_occ_gamma0(iorb)
    iorb7 = (RDMd%Ncoupled-1)*(igamma-1)+1
    iorb6 = (RDMd%Ncoupled-1)*igamma
    hole(iorb7:iorb6)  = one - RDMd%occ(iorb)
    Dhole_gamma(iorb7:iorb6,igamma)= - Docc_gamma0(iorb)
!- -- - - - - - - - - - (igamma,iorb1) <-> iorb2,igamma1,iorb5  - - - - - - - - - - - -
    do iorb1=1,RDMd%Ncoupled-1
     iorb2 = (RDMd%Ncoupled-1)*(igamma-1)+iorb1             ! iorb2=1,RDMd%Npairs*(RDMd%Ncoupled-1)
! igamma1=RDMd%Npairs+1,RDMd%Npairs*RDMd%Ncoupled
     igamma1 = RDMd%Npairs+iorb2                       
! iorb5=RDMd%Nalpha_elect+1,RDMd%Nalpha_elect+RDMd%Ncoupled*RDMd%Npairs-1
     iorb5 = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-igamma)+iorb1           
     occ_orb = dsin(GAMMAs(igamma1))*dsin(GAMMAs(igamma1))
     Docc_gamma0(iorb5) = dsin(two*GAMMAs(igamma1))
     sqrthole_orb = dsqrt(occ_orb)
     if(sqrthole_orb>tol20) then
      Dsqrt_occ_gamma0(iorb5) = half*Docc_gamma0(iorb5)/sqrthole_orb
     else
      Dsqrt_occ_gamma0(iorb5) = zero
     end if
     RDMd%occ(iorb5) = hole(iorb2)*occ_orb
     sqrt_occ_orb = dsqrt(hole(iorb2))
     sqrt_occ(iorb5) = sqrt_occ_orb*sqrthole_orb
     Docc_gamma(iorb5,igamma) = Dhole_gamma(iorb2,igamma)*occ_orb
     if(sqrt_occ_orb>tol20) then
      Dsqrt_occ_gamma(iorb5,igamma) = half*Dhole_gamma(iorb2,igamma)*sqrthole_orb/sqrt_occ_orb
     else
      Dsqrt_occ_gamma(iorb5,igamma) = zero
     endif
     do iorb8=iorb7,iorb2-1                           ! iorb7 < iorb8 < iorb2-1
      igamma2 = RDMd%Npairs+iorb8                     !   igamma < igamma2 < igamma1
      Docc_gamma(iorb5,igamma2) = Dhole_gamma(iorb2,igamma2)*occ_orb
      if(sqrt_occ_orb>tol20) then
       Dsqrt_occ_gamma(iorb5,igamma2) = half*Dhole_gamma(iorb2,igamma2)*sqrthole_orb/sqrt_occ_orb
      else
       Dsqrt_occ_gamma(iorb5,igamma2) = zero
      endif
     enddo
     Docc_gamma(iorb5,igamma1) = hole(iorb2)*Docc_gamma0(iorb5)
     Dsqrt_occ_gamma(iorb5,igamma1) = sqrt_occ_orb*Dsqrt_occ_gamma0(iorb5)
!- - hole(iorb2+1) - - - - - - - - - - - - - - - -
     if(iorb1<RDMd%Ncoupled-1) then
      do iorb3=1,RDMd%Ncoupled-1-iorb1
       iorb8 = iorb2+iorb3                          ! iorb2 < iorb8 < igamma*(RDMd%Ncoupled-1)
       hole(iorb8)  =  hole(iorb8) - RDMd%occ(iorb5)
       Dhole_gamma(iorb8,igamma)= Dhole_gamma(iorb8,igamma) - Docc_gamma(iorb5,igamma)
       do iorb4=iorb7,iorb2-1
        igamma2 = RDMd%Npairs+iorb4
        Dhole_gamma(iorb8,igamma2)= Dhole_gamma(iorb8,igamma2) - Docc_gamma(iorb5,igamma2)
       enddo
       Dhole_gamma(iorb8,igamma1)= Dhole_gamma(iorb8,igamma1) - Docc_gamma(iorb5,igamma1)
      enddo
     endif
!- - hole(iorb2+1) - - - - - - - - - - - - - - - -
    enddo
!- -- - - - - - - - - - (igamma,iorb1) <-> iorb2,igamma1,iorb5  - - - - - - - - - - - -
!- - iorb2 = iorb6 - last RDMd%occ  - - - - - - - - - - - - - -
    igamma1 = RDMd%Npairs+iorb6               ! igamma1=RDMd%Npairs+igamma*(RDMd%Ncoupled-1)
    iorb5 = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-igamma)+RDMd%Ncoupled
    hole_orb = dcos(GAMMAs(igamma1))*dcos(GAMMAs(igamma1))
    Docc_gamma0(iorb5) = -dsin(two*GAMMAs(igamma1))
    sqrthole_orb = dsqrt(hole_orb)
    if(sqrthole_orb>tol20) then
     Dsqrt_occ_gamma0(iorb5) = half*Docc_gamma0(iorb5)/sqrthole_orb
    else
     Dsqrt_occ_gamma0(iorb5) = zero
    end if
    RDMd%occ(iorb5)  = hole(iorb6)*hole_orb
    sqrt_hole_orb = dsqrt(hole(iorb6))
    sqrt_occ(iorb5)= sqrt_hole_orb*sqrthole_orb
    Docc_gamma(iorb5,igamma) = Dhole_gamma(iorb6,igamma)*hole_orb
    if(sqrt_hole_orb>tol20) then
     Dsqrt_occ_gamma(iorb5,igamma) = half*Dhole_gamma(iorb6,igamma)*sqrthole_orb/sqrt_hole_orb
    else
     Dsqrt_occ_gamma(iorb5,igamma) = zero
    endif
    do iorb8=iorb7,iorb6-1            ! iorb7 < iorb8 < iorb6-1
     igamma2 = RDMd%Npairs+iorb8      !   igamma < igamma2 < igamma1
     Docc_gamma(iorb5,igamma2) = Dhole_gamma(iorb6,igamma2)*hole_orb
     if(sqrt_hole_orb>tol20) then
      Dsqrt_occ_gamma(iorb5,igamma2) = half*Dhole_gamma(iorb6,igamma2)*sqrthole_orb/sqrt_hole_orb
     else
      Dsqrt_occ_gamma(iorb5,igamma2) = zero
     endif
    enddo
    Docc_gamma(iorb5,igamma1) = hole(iorb6) *Docc_gamma0(iorb5)
    Dsqrt_occ_gamma(iorb5,igamma1) = sqrt_hole_orb*Dsqrt_occ_gamma0(iorb5)
!- - iorb2 = iorb6 - last occ  - - - - - - - - - - - - - -
   enddo
  deallocate(hole,Dhole_gamma)
 endif
 deallocate(Docc_gamma0,Dsqrt_occ_gamma0)
!-----------------------------------------------------------------------
!  Compute dynamic occ and its derivative
!-----------------------------------------------------------------------
 RDMd%occ_dyn(1:RDMd%Nfrozen)=one
 do iorb=1,RDMd%Npairs
  iorb1=RDMd%Nfrozen+iorb
  hole_dyn=(one-RDMd%occ(iorb1))/RDMd%Hcut 
  exponential=dexp(-(hole_dyn**two)) 
  dexponent=-two*hole_dyn/RDMd%Hcut 
  RDMd%occ_dyn(iorb1)=RDMd%occ(iorb1)*exponential
  Docc_dyn(iorb1,:)=exponential*Docc_gamma(iorb1,:)*(one-RDMd%occ(iorb1)*dexponent)
  if(RDMd%Ncoupled>1) then  ! Extended
   do iorb2=1,RDMd%Ncoupled  
    iorb3=RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-iorb)+iorb2
    RDMd%occ_dyn(iorb3)=RDMd%occ(iorb3)*exponential
    Docc_dyn(iorb3,:)=exponential*(Docc_gamma(iorb3,:)-RDMd%occ(iorb3)*Docc_gamma(iorb1,:)*dexponent)
   enddo 
  else                      ! Perfect-pairing
   iorb2=RDMd%Nalpha_elect+(RDMd%Npairs-iorb)+1
   RDMd%occ_dyn(iorb2)=RDMd%occ(iorb2)*exponential
   Docc_dyn(iorb2,:)=exponential*(Docc_gamma(iorb2,:)-RDMd%occ(iorb2)*Docc_gamma(iorb1,:)*dexponent)
  endif 
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   If we are computing the chemical potential(s) mu = h_ii + d Vee/dn_i
!(We pretend that occ(iorb) = GAMMAs(iorb) and d occ(iorb) = d GAMMAs(iorb))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(present(chempot)) then
  Docc_gamma = one; Dsqrt_occ_gamma = zero;
  do iorb=1,RDMd%NBF_occ
   if(dabs(RDMd%occ(iorb))>tol20) then
    Dsqrt_occ_gamma(iorb,:)=half/dsqrt(RDMd%occ(iorb))
   endif
  enddo
 endif 
!-----------------------------------------------------------------------
!       DM2_J, DM2_K, DM2_L, DDM2_gamma_J, DDM2_gamma_K, DDM2_gamma_L
!Comment: This is not the cleanest way to call them but it is fastest way 
!         to program them without requiring extra memory allocations or pointers.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 RDMd%Docc_gamma=reshape(Docc_gamma,(/RDMd%NBF_occ*RDMd%Ngammas/))
 if(RDMd%INOF==0) then
  call dm2_hf(RDMd,RDMd%Docc_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==-1) then
  call dm2_pccd(RDMd,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L)
 elseif(RDMd%INOF==100) then
  call dm2_mbb(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==101) then
  call dm2_power(RDMd,RDMd%Docc_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==102) then
  call dm2_ca(RDMd,RDMd%Docc_gamma,sqrt_occ,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==103) then
  call dm2_cga(RDMd,RDMd%Docc_gamma,sqrt_occ,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==104) then
  call dm2_gu(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==5) then
  call dm2_pnof5(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==7) then
  call dm2_pnof7(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==8) then
  call dm2_gnof(RDMd,RDMd%Docc_gamma,Docc_dyn,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 elseif(RDMd%INOF==70) then
  call dm2_pnof7_sup(RDMd,RDMd%Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_J,RDMd%DM2_K,RDMd%DM2_L,&
  & RDMd%DDM2_gamma_J,RDMd%DDM2_gamma_K,RDMd%DDM2_gamma_L)
 else
  ! Nth
 endif
!-----------------------------------------------------------------------
!       DM2_Jsr and DDM2_gamma_Jsr
!                 &
!       DM2_Lsr and DDM2_gamma_Lsr
!          (except for pCDD)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(RDMd%irange_sep/=0 .and. RDMd%INOF>-1) then
  call dm2_hartree(RDMd,RDMd%Docc_gamma,RDMd%DM2_Jsr,RDMd%DDM2_gamma_Jsr)
  if(RDMd%irange_sep==1) then
   call dm2_intra(RDMd,sqrt_occ,Dsqrt_occ_gamma,RDMd%DM2_iiii,RDMd%DM2_Jsr,RDMd%DM2_Lsr,&
  &  RDMd%DDM2_gamma_Jsr,RDMd%DDM2_gamma_Lsr)
  endif
 endif
!-----------------------------------------------------------------------
 deallocate(sqrt_occ,Dsqrt_occ_gamma,Docc_gamma,Docc_dyn)
 
end subroutine gamma_to_2rdm

!!***
!!****f* DoNOF/dm2_hartree
!! NAME
!! dm2_hartree
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for Hartree
!!
!! INPUTS
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
!!
!! OUTPUT
!! DM2_J=DM2 elements that use J integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_hartree(RDMd,Docc_gamma,DM2_J,DDM2_gamma_J)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Docc_gamma
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
!arrays
!************************************************************************

 DM2_J=zero; DDM2_gamma_J=zero; 
!     DM2_Jpq = 2NpNq
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!                 DM2(iorb,iorb,iorb,iorb)=2*occ(iorb)*occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_J(iorb,iorb)=two*RDMd%occ(iorb)*RDMd%occ(iorb)
  DDM2_gamma_J(iorb,iorb,:)=four*Docc_gamma(iorb,:)*RDMd%occ(iorb)
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_hartree
!!***

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
 integer::iorb,iorb1
!arrays
!************************************************************************

 DM2_L=zero; DDM2_gamma_L=zero; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -NpNq [ DDM2_Jpqk = 2DNpk*Nq, DDM2_Kpq = -DNpk*Nq ]
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Docc_gamma(iorb,:)*RDMd%occ(iorb1)
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!                 DM2(iorb,iorb,iorb,iorb)=occ(iorb)*occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)*RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=two*RDMd%occ(iorb)
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
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
 integer::iorb,iorb1
!arrays
!************************************************************************

 DM2_L=zero; DDM2_gamma_L=zero; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -sqrt(NpNq) 
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -sqrt_occ(iorb)*sqrt_occ(iorb1)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Dsqrt_occ_gamma(iorb,:)*sqrt_occ(iorb1)
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!          DM2(iorb,iorb,iorb,iorb)=2*occ(iorb)*occ(iorb)-occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=two*RDMd%occ(iorb)*RDMd%occ(iorb)-RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=four*RDMd%occ(iorb)-one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_mbb
!!***

!!***
!!****f* DoNOF/dm2_ca
!! NAME
!! dm2_ca
!! Equivalent to Constrained-Pairing Mean-Fielf Theory (CPMFT) if used with perfect-pairing
!! (see T. Tsuchimochi et al. J. Chem. Phys., 133, 13, 2010)
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for CA
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
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

subroutine dm2_ca(RDMd,Docc_gamma,sqrt_occ,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
& DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Docc_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 real(dp)::sqrt_hole1,sqrt_hole2
!arrays
!************************************************************************

!     DM2_Jpq = 2NpNq, DM2_Kpq = -NpNq, DM2_Lpq = -sqrt[Np (1-Np) Nq (1-Nq)]
 do iorb=1,RDMd%NBF_occ
  sqrt_hole1=dsqrt(one-RDMd%occ(iorb))
  do iorb1=1,RDMd%NBF_occ
   sqrt_hole2=dsqrt(one-RDMd%occ(iorb1))
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1)  
   DM2_L(iorb,iorb1) = -sqrt_occ(iorb)*sqrt_occ(iorb1)*sqrt_hole1*sqrt_hole2
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_L(iorb,iorb1,:) = -Docc_gamma(iorb,:)*(one-two*RDMd%occ(iorb))*  &
   & ( sqrt_occ(iorb1)*sqrt_hole2/(two*sqrt_occ(iorb)*sqrt_hole1+tol20) )
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!    DM2(iorb,iorb,iorb,iorb)=2*occ(iorb)*occ(iorb)-occ(iorb)*occ(iorb)-occ(iorb)*(1-occ(iorb))
!                            =occ(iorb)*(2*occ(iorb)-1)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)*(two*RDMd%occ(iorb)-one)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  DM2_L(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=four*RDMd%occ(iorb)-one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
  DDM2_gamma_L(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_ca
!!***

!!***
!!****f* DoNOF/dm2_cga
!! NAME
!! dm2_cga
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for CGA
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Docc_gamma=Matrix with the derivative of occ numbers vs gamma
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

subroutine dm2_cga(RDMd,Docc_gamma,sqrt_occ,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
& DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Docc_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 real(dp)::sqrt_hole1,sqrt_hole2
!arrays
!************************************************************************

 DM2_L=zero; DDM2_gamma_L=zero; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -f(Np,Nq)
 do iorb=1,RDMd%NBF_occ
  sqrt_hole1=dsqrt(two-RDMd%occ(iorb))
  do iorb1=1,RDMd%NBF_occ
   sqrt_hole2=dsqrt(two-RDMd%occ(iorb1))
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -half*sqrt_occ(iorb)*sqrt_occ(iorb1)*sqrt_hole1*sqrt_hole2 &
                     & -half*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -half*Docc_gamma(iorb,:)*(one-RDMd%occ(iorb))*  &
   & ( sqrt_occ(iorb1)*sqrt_hole2/(sqrt_occ(iorb)*sqrt_hole1+tol20) ) &
   & -half*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!    DM2(iorb,iorb,iorb,iorb)=2*occ(iorb)*occ(iorb)-[occ(iorb)*occ(iorb)+occ(iorb)*(2-occ(iorb))]/2
!                            =occ(iorb)*(2*occ(iorb)-1)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)*(two*RDMd%occ(iorb)-one)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=four*RDMd%occ(iorb)-one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_cga
!!***

!!***
!!****f* DoNOF/dm2_gu
!! NAME
!! dm2_gu
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for GU
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

subroutine dm2_gu(RDMd,Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
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
 integer::iorb,iorb1
!arrays
!************************************************************************

 DM2_L=zero; DDM2_gamma_L=zero; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -f(Np,Nq)
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -sqrt_occ(iorb)*sqrt_occ(iorb1)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Dsqrt_occ_gamma(iorb,:)*sqrt_occ(iorb1)
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!     DM2(iorb,iorb,iorb,iorb)=2*occ(iorb)*occ(iorb)-occ(iorb)*occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)*RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=two*RDMd%occ(iorb)
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_gu
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

subroutine dm2_power(RDMd,Docc_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
& DDM2_gamma_L)
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
 integer::iorb,iorb1
 real(dp)::ONEmLpower
!arrays
!************************************************************************

 ONEmLpower=one-RDMd%Lpower
 DM2_L=zero; DDM2_gamma_L=zero; 
!     DM2_Jpq = 2NpNq, DM2_Kpq = -f(Np,Nq)
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -((RDMd%occ(iorb)*RDMd%occ(iorb1))**RDMd%Lpower)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -RDMd%Lpower*RDMd%occ(iorb1)*&
   &                            Docc_gamma(iorb,:)*((RDMd%occ(iorb)*RDMd%occ(iorb1))**ONEmLpower)
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!          DM2(iorb,iorb,iorb,iorb)=2*occ(iorb)*occ(iorb)- (occ(iorb)**(2*power)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=two*RDMd%occ(iorb)*RDMd%occ(iorb)-(RDMd%occ(iorb)**(two*RDMd%Lpower))
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=four*RDMd%occ(iorb)-two*RDMd%Lpower*(RDMd%occ(iorb)**(two*RDMd%Lpower-one))
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_power
!!***

!!****f* DoNOF/dm2_intra
!! NAME
!! dm2_intra
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for PNOF5-intra part
!!  JCP 134, 164102, 2011; JCP 139, 234109, 2013
!!
!! INPUTS
!! sqrt_occ=Square root of the occupancies of the frozen + active orbitals
!! Dsqrt_occ_gamma=Matrix with the derivative of sqrt(occ numbers) vs gamma
!!
!! OUTPUT
!! DM2_J=DM2 elements that use J integrals 
!! DM2_L=DM2 elements that use L integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals 
!! DDM2_gamma_L=Derivative of the DM2 elements w.r.t. gamma that use L integrals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_intra(RDMd,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_L,DDM2_gamma_J,DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Dsqrt_occ_gamma
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
!arrays
!************************************************************************

!- - - - - - - - - - - - - - - - - - - - - - - -              
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for PNOF5(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb2=1,RDMd%Npairs
  iorb3 = RDMd%Nfrozen+iorb2
  do iorb1=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
   iorb4 = RDMd%Nfrozen+iorb1
   DM2_J(iorb3,iorb4) = zero
   DM2_J(iorb4,iorb3) = zero
   DM2_L(iorb3,iorb4) = -sqrt_occ(iorb3)*sqrt_occ(iorb4)
   DM2_L(iorb4,iorb3) = -sqrt_occ(iorb4)*sqrt_occ(iorb3)
   DDM2_gamma_J(iorb3,iorb4,:) = zero
   DDM2_gamma_J(iorb4,iorb3,:) = zero
   DDM2_gamma_L(iorb3,iorb4,:) = -Dsqrt_occ_gamma(iorb3,:)*sqrt_occ(iorb4)
   DDM2_gamma_L(iorb4,iorb3,:) = -Dsqrt_occ_gamma(iorb4,:)*sqrt_occ(iorb3)
   do iorb=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
    iorb5 = RDMd%Nfrozen+iorb
    DM2_J(iorb5,iorb4) = zero
    DM2_L(iorb5,iorb4) = sqrt_occ(iorb5)*sqrt_occ(iorb4)
    DDM2_gamma_J(iorb5,iorb4,:) = zero
    DDM2_gamma_L(iorb5,iorb4,:) = Dsqrt_occ_gamma(iorb5,:)*sqrt_occ(iorb4)
   enddo
  enddo
 enddo
!-----------------------------------------------------------------------
!                 DM2(iorb,iorb,iorb,iorb)=occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_L(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_L(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_intra
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
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
!arrays
!************************************************************************

!-----------------------------------------------------------------------
!                Inter-pair interactions for PNOF5 (Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!     DM2_Jpq = 2NpNq, DM2_Kpq = -NpNq [ DDM2_Jpqk = 2DNpk*Nq, DDM2_Kpq = -DNpk*Nq ] 
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_L(iorb,iorb1) = zero
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_L(iorb,iorb1,:) = zero 
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for PNOF5(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb2=1,RDMd%Npairs
  iorb3 = RDMd%Nfrozen+iorb2
  do iorb1=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
   iorb4 = RDMd%Nfrozen+iorb1
   DM2_J(iorb3,iorb4) = zero
   DM2_J(iorb4,iorb3) = zero
   DM2_K(iorb3,iorb4) = zero
   DM2_K(iorb4,iorb3) = zero
   DM2_L(iorb3,iorb4) = -sqrt_occ(iorb3)*sqrt_occ(iorb4)
   DM2_L(iorb4,iorb3) = -sqrt_occ(iorb4)*sqrt_occ(iorb3)
   DDM2_gamma_J(iorb3,iorb4,:) = zero
   DDM2_gamma_J(iorb4,iorb3,:) = zero
   DDM2_gamma_K(iorb3,iorb4,:) = zero
   DDM2_gamma_K(iorb4,iorb3,:) = zero
   DDM2_gamma_L(iorb3,iorb4,:) = -Dsqrt_occ_gamma(iorb3,:)*sqrt_occ(iorb4)
   DDM2_gamma_L(iorb4,iorb3,:) = -Dsqrt_occ_gamma(iorb4,:)*sqrt_occ(iorb3)
   do iorb=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
    iorb5 = RDMd%Nfrozen+iorb
    DM2_J(iorb5,iorb4) = zero
    DM2_K(iorb5,iorb4) = zero
    DM2_L(iorb5,iorb4) = sqrt_occ(iorb5)*sqrt_occ(iorb4)
    DDM2_gamma_J(iorb5,iorb4,:) = zero
    DDM2_gamma_K(iorb5,iorb4,:) = zero
    DDM2_gamma_L(iorb5,iorb4,:) = Dsqrt_occ_gamma(iorb5,:)*sqrt_occ(iorb4)
   enddo
  enddo
 enddo
!-----------------------------------------------------------------------
!                 DM2(iorb,iorb,iorb,iorb)=occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  DM2_L(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
  DDM2_gamma_L(iorb,iorb,:)=zero
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
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
!arrays
 real(dp),allocatable,dimension(:)::FIs
 real(dp),allocatable,dimension(:,:)::DFIs
!************************************************************************

!-----------------------------------------------------------------------
!     FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(FIs(RDMd%NBF_occ),DFIs(RDMd%NBF_occ,RDMd%Ngammas))
 FIs = zero; DFIs = zero;
 if(RDMd%Ista==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (Np*Hp)^1/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = dsqrt( RDMd%occ(iorb)*(one-RDMd%occ(iorb)) )
   if(FIs(iorb)>tol20) then
    DFIs(iorb,:) = (half-RDMd%occ(iorb))*Docc_gamma(iorb,:)/FIs(iorb)
   endif
  enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
 else if(RDMd%Ista==1) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (4*Np*Hp)^0.5*(Np*Hp)^0.5 = 2*Np*Hp
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = two*RDMd%occ(iorb)*(one-RDMd%occ(iorb))
   DFIs(iorb,:) = two*(one-two*RDMd%occ(iorb))*Docc_gamma(iorb,:)
  enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - -       
 end if
!-----------------------------------------------------------------------
!                Inter-pair interactions for PNOF7 (Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1) 
   DM2_L(iorb,iorb1) = -FIs(iorb)*FIs(iorb1)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_L(iorb,iorb1,:) = -DFIs(iorb,:)*FIs(iorb1)
  enddo
 enddo
 deallocate(FIs,DFIs)
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for PNOF7(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb2=1,RDMd%Npairs
  iorb3 = RDMd%Nfrozen+iorb2
  do iorb1=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
   iorb4 = RDMd%Nfrozen+iorb1
   DM2_J(iorb3,iorb4) = zero
   DM2_J(iorb4,iorb3) = zero
   DM2_K(iorb3,iorb4) = zero
   DM2_K(iorb4,iorb3) = zero
   DM2_L(iorb3,iorb4) = -sqrt_occ(iorb3)*sqrt_occ(iorb4)
   DM2_L(iorb4,iorb3) = -sqrt_occ(iorb4)*sqrt_occ(iorb3)
   DDM2_gamma_J(iorb3,iorb4,:) = zero
   DDM2_gamma_J(iorb4,iorb3,:) = zero
   DDM2_gamma_K(iorb3,iorb4,:) = zero
   DDM2_gamma_K(iorb4,iorb3,:) = zero
   DDM2_gamma_L(iorb3,iorb4,:) = -Dsqrt_occ_gamma(iorb3,:)*sqrt_occ(iorb4)
   DDM2_gamma_L(iorb4,iorb3,:) = -Dsqrt_occ_gamma(iorb4,:)*sqrt_occ(iorb3)
   do iorb=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
    iorb5 = RDMd%Nfrozen+iorb
    DM2_J(iorb5,iorb4) = zero
    DM2_K(iorb5,iorb4) = zero
    DM2_L(iorb5,iorb4) = sqrt_occ(iorb5)*sqrt_occ(iorb4)
    DDM2_gamma_J(iorb5,iorb4,:) = zero
    DDM2_gamma_K(iorb5,iorb4,:) = zero
    DDM2_gamma_L(iorb5,iorb4,:) = Dsqrt_occ_gamma(iorb5,:)*sqrt_occ(iorb4)
   enddo
  enddo
 enddo
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                 DM2(iorb,iorb,iorb,iorb)=occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  DM2_L(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
  DDM2_gamma_L(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_pnof7
!!***

!!****f* DoNOF/dm2_gnof
!! NAME
!! dm2_gnof
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for GNOF
!!  PRL 127, 233001, 2021
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

subroutine dm2_gnof(RDMd,Docc_gamma,Docc_dyn,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,&
& DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(in)::sqrt_occ
 real(dp),dimension(RDMd%NBF_occ,RDMd%Ngammas),intent(in)::Dsqrt_occ_gamma,Docc_gamma,Docc_dyn
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ,RDMd%Ngammas),intent(inout)::DDM2_gamma_J,DDM2_gamma_K,DDM2_gamma_L
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
!arrays
 real(dp),allocatable,dimension(:)::FIs,sqrt_occ_dyn
 real(dp),allocatable,dimension(:,:)::DFIs,Dsqrt_occ_dyn
!************************************************************************

!-----------------------------------------------------------------------
!     Dynamic occupation numbers and derivatives
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(sqrt_occ_dyn(RDMd%NBF_occ),Dsqrt_occ_dyn(RDMd%NBF_occ,RDMd%Ngammas))
 sqrt_occ_dyn = zero; Dsqrt_occ_dyn = zero;
 sqrt_occ_dyn(:)=dsqrt(RDMd%occ_dyn(:))
 do iorb=1,RDMd%NBF_occ 
  if(sqrt_occ_dyn(iorb)>tol20) then
   Dsqrt_occ_dyn(iorb,:)=half*Docc_dyn(iorb,:)/sqrt_occ_dyn(iorb)
  endif
 enddo
!-----------------------------------------------------------------------
!     FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(FIs(RDMd%NBF_occ),DFIs(RDMd%NBF_occ,RDMd%Ngammas))
 FIs = zero; DFIs = zero;
 if(RDMd%Ista==0 .or. RDMd%Ista==3) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (Np*Hp)^1/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = dsqrt( RDMd%occ(iorb)*(one-RDMd%occ(iorb)) )
   if(FIs(iorb)>tol20) then
    DFIs(iorb,:) = (half-RDMd%occ(iorb))*Docc_gamma(iorb,:)/FIs(iorb)
   endif
  enddo
 else if(RDMd%Ista==1) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (4*Np*Hp)^0.5*(Np*Hp)^0.5 = 2*Np*Hp
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = two*RDMd%occ(iorb)*(one-RDMd%occ(iorb))
   DFIs(iorb,:) = two*(one-two*RDMd%occ(iorb))*Docc_gamma(iorb,:)
  enddo
 endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
!-----------------------------------------------------------------------
!                Inter-pair interactions for GNOF (Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
!   HF-like (keeping only this term is PNOF5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1) 
   DM2_L(iorb,iorb1) = zero
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_L(iorb,iorb1,:) = zero
  enddo
 enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
!   Non-dynamic
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
 do iorb=RDMd%Nfrozen+1,RDMd%Nbeta_elect
  do iorb1=RDMd%Nalpha_elect+1,RDMd%NBF_occ
   DM2_L(iorb,iorb1) = -FIs(iorb)*FIs(iorb1)
   DDM2_gamma_L(iorb,iorb1,:) = -DFIs(iorb,:)*FIs(iorb1)
  enddo
!
!  ! Including this and removing the dynamic contrib., it is PNOF7
!  do iorb1=RDMd%Nfrozen+1,RDMd%Nbeta_elect
!   DM2_L(iorb,iorb1) = -FIs(iorb)*FIs(iorb1)
!   DDM2_gamma_L(iorb,iorb1,:) = -DFIs(iorb,:)*FIs(iorb1)
!  enddo
!
 enddo
 do iorb=RDMd%Nalpha_elect+1,RDMd%NBF_occ
  do iorb1=RDMd%Nfrozen+1,RDMd%Nbeta_elect
   DM2_L(iorb,iorb1) = -FIs(iorb)*FIs(iorb1)  
   DDM2_gamma_L(iorb,iorb1,:) = -DFIs(iorb,:)*FIs(iorb1)
  enddo
  do iorb1=RDMd%Nalpha_elect+1,RDMd%NBF_occ
   DM2_L(iorb,iorb1) = -FIs(iorb)*FIs(iorb1)  
   DDM2_gamma_L(iorb,iorb1,:) = -DFIs(iorb,:)*FIs(iorb1)
  enddo
 enddo
 deallocate(FIs,DFIs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   Dynamic
! (if Ista = 1 or 3, it is ignored and we only retain the non-dyn/static contrib)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(RDMd%Ista==0) then
  do iorb=RDMd%Nfrozen+1,RDMd%Nbeta_elect
   do iorb1=RDMd%Nalpha_elect+1,RDMd%NBF_occ
    DM2_L(iorb,iorb1) = DM2_L(iorb,iorb1)-sqrt_occ_dyn(iorb)*sqrt_occ_dyn(iorb1)+RDMd%occ_dyn(iorb)*RDMd%occ_dyn(iorb1)
    DDM2_gamma_L(iorb,iorb1,:) = DDM2_gamma_L(iorb,iorb1,:)-Dsqrt_occ_dyn(iorb,:)*sqrt_occ_dyn(iorb1) &
   &  +Docc_dyn(iorb,:)*RDMd%occ_dyn(iorb1)      
   enddo
  enddo
  do iorb=RDMd%Nalpha_elect+1,RDMd%NBF_occ
   do iorb1=RDMd%Nfrozen+1,RDMd%Nbeta_elect
    DM2_L(iorb,iorb1) = DM2_L(iorb,iorb1)-sqrt_occ_dyn(iorb)*sqrt_occ_dyn(iorb1)+RDMd%occ_dyn(iorb)*RDMd%occ_dyn(iorb1)
    DDM2_gamma_L(iorb,iorb1,:) = DDM2_gamma_L(iorb,iorb1,:)-Dsqrt_occ_dyn(iorb,:)*sqrt_occ_dyn(iorb1) &
   &  +Docc_dyn(iorb,:)*RDMd%occ_dyn(iorb1)      
   enddo
   do iorb1=RDMd%Nalpha_elect+1,RDMd%NBF_occ
    DM2_L(iorb,iorb1) = DM2_L(iorb,iorb1)+sqrt_occ_dyn(iorb)*sqrt_occ_dyn(iorb1)+RDMd%occ_dyn(iorb)*RDMd%occ_dyn(iorb1)
    DDM2_gamma_L(iorb,iorb1,:) = DDM2_gamma_L(iorb,iorb1,:)+Dsqrt_occ_dyn(iorb,:)*sqrt_occ_dyn(iorb1) &
   &  +Docc_dyn(iorb,:)*RDMd%occ_dyn(iorb1)      
   enddo
  enddo
 endif
 deallocate(sqrt_occ_dyn,Dsqrt_occ_dyn)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for GNOF(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb2=1,RDMd%Npairs
  iorb3 = RDMd%Nfrozen+iorb2
  do iorb1=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
   iorb4 = RDMd%Nfrozen+iorb1
   DM2_J(iorb3,iorb4) = zero
   DM2_J(iorb4,iorb3) = zero
   DM2_K(iorb3,iorb4) = zero
   DM2_K(iorb4,iorb3) = zero
   DM2_L(iorb3,iorb4) = -sqrt_occ(iorb3)*sqrt_occ(iorb4)
   DM2_L(iorb4,iorb3) = -sqrt_occ(iorb4)*sqrt_occ(iorb3)
   DDM2_gamma_J(iorb3,iorb4,:) = zero
   DDM2_gamma_J(iorb4,iorb3,:) = zero
   DDM2_gamma_K(iorb3,iorb4,:) = zero
   DDM2_gamma_K(iorb4,iorb3,:) = zero
   DDM2_gamma_L(iorb3,iorb4,:) = -Dsqrt_occ_gamma(iorb3,:)*sqrt_occ(iorb4)
   DDM2_gamma_L(iorb4,iorb3,:) = -Dsqrt_occ_gamma(iorb4,:)*sqrt_occ(iorb3)
   do iorb=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
    iorb5 = RDMd%Nfrozen+iorb
    DM2_J(iorb5,iorb4) = zero
    DM2_K(iorb5,iorb4) = zero
    DM2_L(iorb5,iorb4) = sqrt_occ(iorb5)*sqrt_occ(iorb4)
    DDM2_gamma_J(iorb5,iorb4,:) = zero
    DDM2_gamma_K(iorb5,iorb4,:) = zero
    DDM2_gamma_L(iorb5,iorb4,:) = Dsqrt_occ_gamma(iorb5,:)*sqrt_occ(iorb4)
   enddo
  enddo
 enddo
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                 DM2(iorb,iorb,iorb,iorb)=occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  DM2_L(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
  DDM2_gamma_L(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_gnof
!!***

!!***
!!****f* DoNOF/dm2_pccd
!! NAME
!! dm2_pccd
!!
!! FUNCTION
!!  Build from the occ numbers and the 2-RDM elements for pCCD
!!  We may need to build RDMd%DM2_Jsr in here for range-sep
!!
!! INPUTS
!!
!! OUTPUT
!! DM2_iiii=DM2 same orb elements
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DM2_L=DM2 elements that use L integrals 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dm2_pccd(RDMd,DM2_iiii,DM2_J,DM2_K,DM2_L)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_occ),intent(inout)::DM2_iiii
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K,DM2_L
!Local variables ------------------------------
!scalars
 logical::file_exists
 integer::iorb,iorb1,iorb2,iorb3
 real(dp)::occ,err_HermJ,err_HermK,err_HermL
!arrays
 real(dp),allocatable,dimension(:,:)::xij,xab,xia
 character(len=200)::msg
!************************************************************************

!- - - - - - - - - - - - - - - - - - - - - - - -              
!  Build intermediates
!- - - - - - - - - - - - - - - - - - - - - - - -              
 allocate(xij(RDMd%Npairs,RDMd%Npairs));xij=zero;
 allocate(xab(RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs),RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)));xab=zero;
 allocate(xia(RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)));xia=zero;
 
 ! xij
 do iorb=1,RDMd%Npairs ! Occ
  do iorb1=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
   xij(iorb,:)=xij(iorb,:)+RDMd%t_pccd(iorb,iorb1)*RDMd%z_pccd(:,iorb1)
  enddo
 enddo

 ! xab
 do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
  do iorb1=1,RDMd%Npairs ! Occ
   xab(iorb,:)=xab(iorb,:)+RDMd%t_pccd(iorb1,:)*RDMd%z_pccd(iorb1,iorb)
  enddo
 enddo

 ! Fix the xij and the xab elements that lead to N-rep. violations
 occ=zero
 do iorb=1,RDMd%Npairs ! Occ
  occ=occ+(one-xij(iorb,iorb))
 enddo
 do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
  occ=occ+xab(iorb,iorb)
 enddo
 do iorb=1,RDMd%Npairs ! Occ
  if(xij(iorb,iorb)>one .or. xij(iorb,iorb)<zero) then
   write(msg,'(a,i5,f10.5)') 'Fixing the xii element ',iorb,xij(iorb,iorb)
   xij(iorb,iorb)=0.5d0
  endif
  occ=occ-(one-xij(iorb,iorb))
 enddo
 do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
  if(.not.(xab(iorb,iorb)>one .or. xab(iorb,iorb)<zero)) then
   occ=occ-xab(iorb,iorb) 
  endif
 enddo
 do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
  if((xab(iorb,iorb)>one .or. xab(iorb,iorb)<zero) .and. occ>zero) then
   write(msg,'(a,i5,f10.5)') 'Fixing the xaa element ',iorb,xab(iorb,iorb)
   if(occ>=0.5d0) then
    xab(iorb,iorb)=0.5d0
    occ=occ-0.5d0
   else
    xab(iorb,iorb)=occ
    occ=zero
   endif
  endif
 enddo

 ! xia
 do iorb=1,RDMd%Npairs ! Occ
  do iorb1=1,RDMd%Npairs ! Occ
   xia(iorb,:)=xia(iorb,:)+xij(iorb,iorb1)*RDMd%t_pccd(iorb1,:) 
  enddo
 enddo

 ! Build Occ
 RDMd%occ=zero
 RDMd%occ(1:RDMd%Nfrozen)=one
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  RDMd%occ(iorb1)=one-xij(iorb,iorb)
 enddo
 do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
  iorb1=iorb+RDMd%Nfrozen+RDMd%Npairs
  RDMd%occ(iorb1)=xab(iorb,iorb)
 enddo 

 ! Hartree terms DM2_Jpq = 2NpNq
 DM2_L=zero; DM2_K=zero;
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
  enddo
 enddo

 ! Build D_ii^ii
 do iorb=1,RDMd%Nfrozen
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
 enddo
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  DM2_iiii(iorb1)=one-xij(iorb,iorb)
  DM2_J(iorb1,iorb1)=zero
 enddo
 do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
  iorb1=iorb+RDMd%Nfrozen+RDMd%Npairs
  DM2_iiii(iorb1)=xab(iorb,iorb)
  DM2_J(iorb1,iorb1)=zero
 enddo

 ! Build D_ii^pp, D_pp^ii, D_ip,ip, and D_pi,pi 
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  ! D_ii^jj and D_ij^ij for i/=j
  do iorb2=1,RDMd%Npairs ! Occ
   iorb3=iorb2+RDMd%Nfrozen
   if(iorb/=iorb2) then
    DM2_J(iorb1,iorb3)=two*(one-xij(iorb,iorb)-xij(iorb2,iorb2))       
    DM2_L(iorb1,iorb3)=xij(iorb,iorb2)  
   endif
  enddo
  ! D_ii^aa, D_aa^ii D_ia^ia
  do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
   iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
   DM2_J(iorb1,iorb3)=two*(xab(iorb2,iorb2)-RDMd%t_pccd(iorb,iorb2)*RDMd%z_pccd(iorb,iorb2))
   DM2_J(iorb3,iorb1)=DM2_J(iorb1,iorb3)
   DM2_L(iorb1,iorb3)=RDMd%t_pccd(iorb,iorb2)+xia(iorb,iorb2)-two*RDMd%t_pccd(iorb,iorb2) &
  &                  *(xab(iorb2,iorb2)+xij(iorb,iorb)-RDMd%t_pccd(iorb,iorb2)*RDMd%z_pccd(iorb,iorb2))
   DM2_L(iorb3,iorb1)=RDMd%z_pccd(iorb,iorb2)
  enddo
 enddo  
 
 ! Build D_aa^bb for a/=b
 do iorb=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
  iorb1=iorb+RDMd%Nfrozen+RDMd%Npairs
  do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
   iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
   if(iorb/=iorb2) then
    DM2_L(iorb1,iorb3)=xab(iorb,iorb2)
    DM2_J(iorb1,iorb3)=zero
   endif
  enddo
 enddo

 ! Exchange elements D_pq^qp = -1/2 D_pq^pq for p/=q
 DM2_K=-half*DM2_J

 ! Check for the Hermiticity of the 2-RDM 
 err_HermJ=zero;err_HermK=zero;err_HermL=zero;
 do iorb=1,RDMd%NBF_occ
  do iorb1=iorb,RDMd%NBF_occ
   err_HermJ=err_HermJ+abs(DM2_J(iorb,iorb1)-DM2_J(iorb1,iorb))
   err_HermK=err_HermK+abs(DM2_K(iorb,iorb1)-DM2_K(iorb1,iorb))
   err_HermL=err_HermL+abs(DM2_L(iorb,iorb1)-DM2_L(iorb1,iorb))
  enddo
 enddo

 file_exists=.false.
 inquire(file='no_pccd_hermiticity', exist=file_exists)
 if(.not.file_exists) then
  if(abs(err_HermJ)>tol6) then
   write(msg,'(a,f15.6)') 'Hermiticity error DM2_J       =',err_HermJ
   call write_output(msg)
   write(msg,'(a)') 'Enforcing Hermiticity in DM2_J'
   call write_output(msg)
   DM2_J=HALF*(DM2_J+transpose(DM2_J))
  endif
  if(abs(err_HermK)>tol6) then
   write(msg,'(a,f15.6)') 'Hermiticity error DM2_K       =',err_HermK
   call write_output(msg)
   write(msg,'(a)') 'Enforcing Hermiticity in DM2_K'
   call write_output(msg)
   DM2_K=HALF*(DM2_K+transpose(DM2_K))
  endif
  if(abs(err_HermL)>tol6) then
   write(msg,'(a,f15.6)') 'Hermiticity error DM2_L       =',err_HermL
   call write_output(msg)
   write(msg,'(a)') 'Enforcing Hermiticity in DM2_L'
   call write_output(msg)
   DM2_L=HALF*(DM2_L+transpose(DM2_L))
  endif
 endif

!-----------------------------------------------------------------------
 deallocate(xij,xab,xia)
!-----------------------------------------------------------------------
end subroutine dm2_pccd
!!***


!!****f* DoNOF/dm2_pnof7_sup
!! NAME
!! dm2_pnof7_sup
!!
!! FUNCTION
!!  Build from the occ numbers and its derivatives the 2-RDM elements and its derivatives w.r.t. gamma for PNOF7 for superconductors
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

subroutine dm2_pnof7_sup(RDMd,Docc_gamma,sqrt_occ,Dsqrt_occ_gamma,DM2_iiii,DM2_J,DM2_K,DM2_L,DDM2_gamma_J,DDM2_gamma_K,&
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
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
!arrays
 real(dp),allocatable,dimension(:)::FIs
 real(dp),allocatable,dimension(:,:)::DFIs
!************************************************************************

!-----------------------------------------------------------------------
!     FIs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(FIs(RDMd%NBF_occ),DFIs(RDMd%NBF_occ,RDMd%Ngammas))
 FIs = zero; DFIs = zero;
 if(RDMd%Ista==0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (Np*Hp)^1/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = dsqrt( RDMd%occ(iorb)*(one-RDMd%occ(iorb)) )
   if(FIs(iorb)>tol20) then
    DFIs(iorb,:) = (half-RDMd%occ(iorb))*Docc_gamma(iorb,:)/FIs(iorb)
   endif
  enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - -              
 else if(RDMd%Ista==1) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  
!      FIs = (4*Np*Hp)^0.5*(Np*Hp)^0.5 = 2*Np*Hp
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
  do iorb=RDMd%Nfrozen+1,RDMd%NBF_occ
   FIs(iorb) = two*RDMd%occ(iorb)*(one-RDMd%occ(iorb))
   DFIs(iorb,:) = two*(one-two*RDMd%occ(iorb))*Docc_gamma(iorb,:)
  enddo
!- - - - - - - - - - - - - - - - - - - - - - - - - - -       
 end if
!-----------------------------------------------------------------------
!                Inter-pair interactions for PNOF7 (Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb=1,RDMd%NBF_occ
  do iorb1=1,RDMd%NBF_occ
   DM2_J(iorb,iorb1) = two*RDMd%occ(iorb)*RDMd%occ(iorb1)
   DM2_K(iorb,iorb1) = -RDMd%occ(iorb)*RDMd%occ(iorb1) 
   DM2_L(iorb,iorb1) = -RDMd%PNOF7sup*FIs(iorb)*FIs(iorb1)
   DDM2_gamma_J(iorb,iorb1,:) = two*Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_K(iorb,iorb1,:) = -Docc_gamma(iorb,:)*RDMd%occ(iorb1)
   DDM2_gamma_L(iorb,iorb1,:) = -RDMd%PNOF7sup*DFIs(iorb,:)*FIs(iorb1)
  enddo
 enddo
 deallocate(FIs,DFIs)
!- - - - - - - - - - - - - - - - - - - - - - - -              
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                Intra-pair interactions for PNOF7(Nc)
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 do iorb2=1,RDMd%Npairs
  iorb3 = RDMd%Nfrozen+iorb2
  do iorb1=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
   iorb4 = RDMd%Nfrozen+iorb1
   DM2_J(iorb3,iorb4) = zero
   DM2_J(iorb4,iorb3) = zero
   DM2_K(iorb3,iorb4) = zero
   DM2_K(iorb4,iorb3) = zero
   DM2_L(iorb3,iorb4) = -RDMd%PNOF7sup*sqrt_occ(iorb3)*sqrt_occ(iorb4)
   DM2_L(iorb4,iorb3) = -RDMd%PNOF7sup*sqrt_occ(iorb4)*sqrt_occ(iorb3)
   DDM2_gamma_J(iorb3,iorb4,:) = zero
   DDM2_gamma_J(iorb4,iorb3,:) = zero
   DDM2_gamma_K(iorb3,iorb4,:) = zero
   DDM2_gamma_K(iorb4,iorb3,:) = zero
   DDM2_gamma_L(iorb3,iorb4,:) = -RDMd%PNOF7sup*Dsqrt_occ_gamma(iorb3,:)*sqrt_occ(iorb4)
   DDM2_gamma_L(iorb4,iorb3,:) = -RDMd%PNOF7sup*Dsqrt_occ_gamma(iorb4,:)*sqrt_occ(iorb3)
   do iorb=RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2)+1,RDMd%Npairs_p_sing+RDMd%Ncoupled*(RDMd%Npairs-iorb2+1)
    iorb5 = RDMd%Nfrozen+iorb
    DM2_J(iorb5,iorb4) = zero
    DM2_K(iorb5,iorb4) = zero
    DM2_L(iorb5,iorb4) = RDMd%PNOF7sup*sqrt_occ(iorb5)*sqrt_occ(iorb4)
    DDM2_gamma_J(iorb5,iorb4,:) = zero
    DDM2_gamma_K(iorb5,iorb4,:) = zero
    DDM2_gamma_L(iorb5,iorb4,:) = RDMd%PNOF7sup*Dsqrt_occ_gamma(iorb5,:)*sqrt_occ(iorb4)
   enddo
  enddo
 enddo
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!                 DM2(iorb,iorb,iorb,iorb)=occ(iorb)
!-----------------------------------------------------------------------
 do iorb=1,RDMd%NBF_occ
  DM2_iiii(iorb)=RDMd%occ(iorb)
  DM2_J(iorb,iorb)=zero
  DM2_K(iorb,iorb)=zero
  DM2_L(iorb,iorb)=zero
  RDMd%Dfni_ni(iorb)=one
  DDM2_gamma_J(iorb,iorb,:)=zero
  DDM2_gamma_K(iorb,iorb,:)=zero
  DDM2_gamma_L(iorb,iorb,:)=zero
 enddo
!-----------------------------------------------------------------------
end subroutine dm2_pnof7_sup
!!***

end module m_gammatodm2
!!***
