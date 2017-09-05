!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the perturbation theory to 2nd order evaluation of the self-energy
!
!=========================================================================
subroutine pt2_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,se,emp2)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_warning
 use m_timing
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: selfenergy_approx,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
 real(dp),intent(out)       :: emp2
!=====
 integer                 :: pstate,qstate
 complex(dp),allocatable :: selfenergy_ring(:,:,:)
 complex(dp),allocatable :: selfenergy_sox(:,:,:)
 integer                 :: iomega
 integer                 :: istate,jstate,kstate
 integer                 :: pqispin,jkspin
 real(dp)                :: fact_occ1,fact_occ2
 real(dp)                :: fi,fj,fk,ei,ej,ek
 complex(dp)             :: omega
 complex(dp)             :: fact_comp
 real(dp)                :: fact_energy
 real(dp)                :: emp2_sox,emp2_ring
 real(dp),allocatable    :: eri_eigenstate_i(:,:,:,:)
 real(dp)                :: coul_iqjk,coul_ijkq,coul_ipkj
!=====

 call start_clock(timing_mp2_self)

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp


 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 write(stdout,*) 'with the perturbative approach'

 

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
 endif



 allocate(selfenergy_ring(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(selfenergy_sox (-se%nomega:se%nomega,nsemin:nsemax,nspin))


 selfenergy_ring(:,:,:) = 0.0_dp
 selfenergy_sox(:,:,:)  = 0.0_dp

 do pqispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !LOOP of the first Green's function
     if( MODULO( istate - (ncore_G+1) , nproc_ortho ) /= rank_ortho ) cycle

     if( .NOT. has_auxil_basis ) then
       call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,istate,pqispin,eri_eigenstate_i)
     endif

     fi = occupation(istate,pqispin)
     ei = energy(istate,pqispin)

     do pstate=nsemin,nsemax ! external loop ( bra )
       qstate=pstate         ! external loop ( ket )


       do jkspin=1,nspin
         do jstate=ncore_G+1,nvirtual_G-1  !LOOP of the second Green's function
           fj = occupation(jstate,jkspin)
           ej = energy(jstate,jkspin)

           do kstate=ncore_G+1,nvirtual_G-1 !LOOP of the third Green's function
             fk = occupation(kstate,jkspin)
             ek = energy(kstate,jkspin)

             fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**3
             fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**3

             if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle

             if( has_auxil_basis ) then
               coul_ipkj = eri_eigen_ri(istate,pstate,pqispin,kstate,jstate,jkspin)
               coul_iqjk = eri_eigen_ri(istate,qstate,pqispin,jstate,kstate,jkspin)
               if( pqispin == jkspin ) then
                 coul_ijkq = eri_eigen_ri(istate,jstate,pqispin,kstate,qstate,pqispin) 
               endif
             else
               coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
               coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
               if( pqispin == jkspin ) then
                 coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin) 
               endif
             endif

             do iomega=-se%nomega,se%nomega
               omega = energy(qstate,pqispin) + se%omega(iomega)

               fact_comp   = fact_occ1 / ( omega - ei + ej - ek + ieta ) &
                           + fact_occ2 / ( omega - ei + ej - ek - ieta )
               fact_energy = REAL( fact_occ1 / (energy(pstate,pqispin) - ei + ej - ek + ieta) , dp )

               selfenergy_ring(iomega,pstate,pqispin) = selfenergy_ring(iomega,pstate,pqispin) &
                        + fact_comp * coul_ipkj * coul_iqjk * spin_fact

               if(iomega==0 .AND. occupation(pstate,pqispin)>completely_empty) then
                 emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                       * fact_energy * coul_ipkj * coul_iqjk * spin_fact
               endif

               if( pqispin == jkspin ) then

                 selfenergy_sox(iomega,pstate,pqispin) = selfenergy_sox(iomega,pstate,pqispin) &
                          - fact_comp * coul_ipkj * coul_ijkq

                 if(iomega==0 .AND. occupation(pstate,pqispin)>completely_empty) then
                   emp2_sox = emp2_sox - occupation(pstate,pqispin) &
                             * fact_energy * coul_ipkj * coul_ijkq
                 endif

               endif
  
  
             enddo ! iomega
  
           enddo
         enddo
       enddo 
     enddo
   enddo 
 enddo ! pqispin

 call xsum_ortho(selfenergy_ring)
 call xsum_ortho(selfenergy_sox)
 call xsum_ortho(emp2_ring)
 call xsum_ortho(emp2_sox)

 emp2_ring = 0.5_dp * emp2_ring
 emp2_sox  = 0.5_dp * emp2_sox

 if( selfenergy_approx == ONE_RING ) then
   emp2_sox = 0.0_dp
   selfenergy_sox(:,:,:) = 0.0_dp
 endif
 if( selfenergy_approx == SOX ) then
   emp2_ring = 0.0_dp
   selfenergy_ring(:,:,:) = 0.0_dp
 endif

 if( nsemin <= ncore_G+1 .AND. nsemax >= nhomo_G ) then
   emp2 = emp2_ring + emp2_sox
   write(stdout,'(/,a)')       ' MP2 Energy'
   write(stdout,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
   write(stdout,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
   write(stdout,'(a,f14.8,/)') ' MP2 correlation :',emp2
 else
   emp2 = 0.0_dp
 endif

 se%sigma(:,:,:) = selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:)


 if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)
 if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_mp2_self)

end subroutine pt2_selfenergy


!=========================================================================
subroutine onering_selfenergy(nstate,basis,occupation,energy,c_matrix,se,emp2)
 use m_definitions
 use m_mpi
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_spectral_function
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
 real(dp),intent(out)       :: emp2
!=====
 type(spectral_function) :: vchi0v
!=====

 call start_clock(timing_mp2_self)

 if( .NOT. has_auxil_basis ) &
   call die('onering_selfenergy: only implemented when an auxiliary basis is available')

 emp2 = 0.0_dp


 write(stdout,'(/,a)') ' Perform the one-ring self-energy calculation'
 write(stdout,*) 'with the perturbative approach'

 call init_spectral_function(nstate,occupation,0,vchi0v)

 call polarizability_onering(basis,nstate,energy,c_matrix,vchi0v)

#ifdef HAVE_SCALAPACK
 call gw_selfenergy_scalapack(ONE_RING,nstate,basis,occupation,energy,c_matrix,vchi0v,se)
#else
 call gw_selfenergy(ONE_RING,nstate,basis,occupation,energy,c_matrix,vchi0v,se,emp2)
#endif
 
 call destroy_spectral_function(vchi0v)

 call stop_clock(timing_mp2_self)


end subroutine onering_selfenergy


!=========================================================================
subroutine pt2_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,selfenergy,emp2)
 use m_definitions
 use m_mpi
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(out)       :: selfenergy(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: emp2
!=====
 integer                 :: pstate,qstate
 complex(dp),allocatable :: selfenergy_ring(:,:,:)
 complex(dp),allocatable :: selfenergy_sox(:,:,:)
 integer                 :: istate,jstate,kstate
 integer                 :: pqispin,jkspin
 real(dp)                :: fact_occ1,fact_occ2
 real(dp)                :: fi,fj,fk,ei,ej,ek,ep,eq
 complex(dp)             :: fact_comp
 real(dp)                :: fact_energy
 real(dp)                :: emp2_sox,emp2_ring
 real(dp),allocatable    :: eri_eigenstate_i(:,:,:,:)
 real(dp)                :: coul_iqjk,coul_ijkq,coul_ipkj
!=====

 call start_clock(timing_mp2_self)

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp


 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 write(stdout,*) 'with the QP self-consistent approach'

 

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
 endif



 allocate(selfenergy_ring (nsemin:nsemax,nsemin:nsemax,nspin))
 allocate(selfenergy_sox  (nsemin:nsemax,nsemin:nsemax,nspin))


 selfenergy_ring(:,:,:) = 0.0_dp
 selfenergy_sox(:,:,:)  = 0.0_dp

 do pqispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !LOOP of the first Green's function
     if( MODULO( istate - (ncore_G+1) , nproc_ortho ) /= rank_ortho ) cycle

     if( .NOT. has_auxil_basis ) then
       call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,istate,pqispin,eri_eigenstate_i)
     endif

     fi = occupation(istate,pqispin)
     ei = energy(istate,pqispin)

     do pstate=nsemin,nsemax ! external loop ( bra )
       do qstate=nsemin,nsemax   ! external loop ( ket )


         do jkspin=1,nspin
           do jstate=ncore_G+1,nvirtual_G-1  !LOOP of the second Green's function
             fj = occupation(jstate,jkspin)
             ej = energy(jstate,jkspin)
  
             do kstate=ncore_G+1,nvirtual_G-1 !LOOP of the third Green's function
               fk = occupation(kstate,jkspin)
               ek = energy(kstate,jkspin)
  
               fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**3
               fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**3
  
               if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle

               if( has_auxil_basis ) then
                 coul_ipkj = eri_eigen_ri(istate,pstate,pqispin,kstate,jstate,jkspin)
                 coul_iqjk = eri_eigen_ri(istate,qstate,pqispin,jstate,kstate,jkspin)
                 if( pqispin == jkspin ) then
                   coul_ijkq = eri_eigen_ri(istate,jstate,pqispin,kstate,qstate,pqispin) 
                 endif
               else
                 coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
                 coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
                 if( pqispin == jkspin ) then
                   coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin) 
                 endif
               endif

               ep = energy(pstate,pqispin) 
               eq = energy(qstate,pqispin) 
  
               fact_comp   = fact_occ1 / ( eq - ei + ej - ek + ieta) &
                           + fact_occ2 / ( eq - ei + ej - ek - ieta)
               fact_energy = REAL( fact_occ1 / ( ep - ei + ej - ek + ieta) , dp )
  
               selfenergy_ring(pstate,qstate,pqispin) = selfenergy_ring(pstate,qstate,pqispin) &
                        + fact_comp * coul_ipkj * coul_iqjk * spin_fact

               if(pstate==qstate .AND. occupation(pstate,pqispin)>completely_empty) then
                 emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                       * fact_energy * coul_ipkj * coul_iqjk * spin_fact
               endif
 
               if( pqispin == jkspin ) then

                 selfenergy_sox(pstate,qstate,pqispin) = selfenergy_sox(pstate,qstate,pqispin) &
                          - fact_comp * coul_ipkj * coul_ijkq

                 if(pstate==qstate .AND. occupation(pstate,pqispin)>completely_empty) then
                   emp2_sox = emp2_sox - occupation(pstate,pqispin) &
                             * fact_energy * coul_ipkj * coul_ijkq
                 endif

               endif
    
    
    
             enddo
           enddo
         enddo
       enddo 
     enddo
   enddo 
 enddo ! pqispin

 call xsum_ortho(selfenergy_ring)
 call xsum_ortho(selfenergy_sox)
 call xsum_ortho(emp2_ring)
 call xsum_ortho(emp2_sox)

 emp2_ring = 0.5_dp * emp2_ring
 emp2_sox  = 0.5_dp * emp2_sox
 if( nsemin <= ncore_G+1 .AND. nsemax >= nhomo_G ) then
   emp2 = emp2_ring + emp2_sox
   write(stdout,'(/,a)')       ' MP2 Energy'
   write(stdout,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
   write(stdout,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
   write(stdout,'(a,f14.8,/)') ' MP2 correlation :',emp2
 else
   emp2 = 0.0_dp
 endif


 selfenergy(:,:,:) = REAL( selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:) ,dp)

 call apply_qs_approximation(basis%nbf,nstate,s_matrix,c_matrix,selfenergy)


 if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)
 if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_mp2_self)

end subroutine pt2_selfenergy_qs


!=========================================================================
subroutine pt3_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,se,emp3)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_warning
 use m_timing
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: selfenergy_approx,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
 real(dp),intent(out)       :: emp3
!=====
 integer                 :: pstate,qstate
 complex(dp),allocatable :: selfenergy(:,:,:)
 complex(dp),allocatable :: selfenergy2(:,:,:)
 complex(dp),allocatable :: selfenergygw(:,:,:)
 complex(dp),allocatable :: selfenergy_1rings(:,:,:)
 complex(dp),allocatable :: selfenergy_2rings(:,:,:)
 integer                 :: iomega
 integer                 :: istate,jstate,kstate,lstate
 integer                 :: astate,bstate,cstate,dstate
 integer                 :: jkspin,pqspin
 complex(dp)             :: omega
 complex(dp)             :: denom1,denom2,denom3
 real(dp)                :: num1,num2,num3
 real(dp)                :: num1a,num1b,num2a,num2b,num3a,num3b
 real(dp)                :: numgw
 complex(dp)             :: ieta_tmp
!=====

 ieta_tmp = 0.0_dp  ! ieta
 emp3 = 0.0_dp
 write(stdout,'(/,a)') ' Perform the third-order self-energy calculation'
 write(stdout,*) 'with the perturbative approach'

 if( nspin /= 1 ) call die('pt3_selfenergy: only implemented for spin restricted calculations')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif

 allocate(selfenergy(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(selfenergy2(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(selfenergygw(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(selfenergy_1rings(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(selfenergy_2rings(-se%nomega:se%nomega,nsemin:nsemax,nspin))

 write(stdout,'(/,1x,a)') ' state    3rd order diagrams    2nd order diagrams   Ring diagrams (contained in PT2 and PT3)'

 selfenergy(:,:,:) = 0.0_dp
 selfenergy2(:,:,:) = 0.0_dp
 selfenergygw(:,:,:) = 0.0_dp
 selfenergy_1rings(:,:,:) = 0.0_dp
 selfenergy_2rings(:,:,:) = 0.0_dp

 pqspin = 1
 do pstate=nsemin,nsemax
   qstate = pstate


   ! B1 i,j    a
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
       num1 = 2.0_dp * eri_eigen(pstate,istate,pqspin,astate,jstate,pqspin) - eri_eigen(pstate,jstate,pqspin,astate,istate,pqspin)
       num2 = eri_eigen(qstate,istate,pqspin,astate,jstate,pqspin)
       numgw = 2.0_dp * eri_eigen(pstate,istate,pqspin,astate,jstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(jstate,pqspin) - ieta
         selfenergy2(iomega,pstate,pqspin) = selfenergy2(iomega,pstate,pqspin) + num1 * num2 / denom1
         selfenergy_1rings(iomega,pstate,pqspin) = selfenergy_1rings(iomega,pstate,pqspin) + numgw * num2 / denom1
         selfenergygw(iomega,pstate,pqspin) = selfenergygw(iomega,pstate,pqspin) + numgw * num2 / denom1
       enddo
     enddo
   enddo
   enddo

   ! B2 i    a,b
   do istate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       num1 = 2.0_dp * eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)
       num2 = eri_eigen(qstate,astate,pqspin,istate,bstate,pqspin)
       numgw = 2.0_dp * eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
         selfenergy2(iomega,pstate,pqspin) = selfenergy2(iomega,pstate,pqspin) + num1 * num2 / denom1 
         selfenergy_1rings(iomega,pstate,pqspin) = selfenergy_1rings(iomega,pstate,pqspin) + numgw * num2 / denom1 
         selfenergygw(iomega,pstate,pqspin) = selfenergygw(iomega,pstate,pqspin) + numgw * num2 / denom1 
       enddo
     enddo
     enddo
   enddo

   ! A1   i,j,k   a,b
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
       denom2 = energy(kstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
       num1 = 2.0_dp * eri_eigen(pstate,qstate,pqspin,kstate,jstate,pqspin) - eri_eigen(pstate,jstate,pqspin,kstate,qstate,pqspin)
       num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
       num3 = eri_eigen(astate,kstate,pqspin,bstate,istate,pqspin)
       selfenergy(:,pstate,pqspin) = selfenergy(:,pstate,pqspin) - num1 * num2 * num3 / ( denom1 * denom2 ) 
     enddo
     enddo
   enddo
   enddo
   enddo

   ! A2   i,j   a,b,c
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
     do cstate=nhomo_G+1,nvirtual_G-1
       denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
       denom2 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(cstate,pqspin)
       num1 = 2.0_dp * eri_eigen(pstate,qstate,pqspin,cstate,bstate,pqspin) - eri_eigen(pstate,bstate,pqspin,cstate,qstate,pqspin)
       num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
       num3 = eri_eigen(istate,cstate,pqspin,jstate,astate,pqspin)
       selfenergy(:,pstate,pqspin) = selfenergy(:,pstate,pqspin) + num1 * num2 * num3 / ( denom1 * denom2 ) 
     enddo
     enddo
     enddo
   enddo
   enddo

   ! A3,A4   i,j   a,b,c
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
     do cstate=nhomo_G+1,nvirtual_G-1
       denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
       denom2 = energy(jstate,pqspin) - energy(cstate,pqspin)
       num1 = 2.0_dp * eri_eigen(pstate,qstate,pqspin,cstate,jstate,pqspin) - eri_eigen(pstate,jstate,pqspin,cstate,qstate,pqspin)
       num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
       num3 = eri_eigen(astate,cstate,pqspin,bstate,istate,pqspin)
       selfenergy(:,pstate,pqspin) = selfenergy(:,pstate,pqspin) + 2.0_dp * num1 * num2 * num3 / ( denom1 * denom2 )
     enddo
     enddo
     enddo
   enddo
   enddo

   ! A5,A6   i,j,k   a,b
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
       denom2 = energy(kstate,pqspin) - energy(bstate,pqspin)
       num1 = 2.0_dp * eri_eigen(pstate,qstate,pqspin,bstate,kstate,pqspin) - eri_eigen(pstate,kstate,pqspin,bstate,qstate,pqspin)
       num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
       num3 = eri_eigen(istate,kstate,pqspin,jstate,astate,pqspin)
       selfenergy(:,pstate,pqspin) = selfenergy(:,pstate,pqspin) - 2.0_dp * num1 * num2 * num3 / ( denom1 * denom2 )
     enddo
     enddo
   enddo
   enddo
   enddo

   ! C1   i   a,b,c,d
   do istate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
     do cstate=nhomo_G+1,nvirtual_G-1
     do dstate=nhomo_G+1,nvirtual_G-1
       num1 = 2.0_dp * eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)
       num2 = eri_eigen(astate,cstate,pqspin,bstate,dstate,pqspin)
       num3 = eri_eigen(qstate,cstate,pqspin,istate,dstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
         denom2 = omega + energy(istate,pqspin) - energy(cstate,pqspin) - energy(dstate,pqspin) - ieta
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) + num1 * num2 * num3 / ( denom1 * denom2 )
       enddo
     enddo
     enddo
     enddo
   enddo
   enddo

   ! C6   i,j,k,l   a
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
   do lstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
       num1 = 2.0_dp * eri_eigen(pstate,kstate,pqspin,astate,lstate,pqspin) - eri_eigen(pstate,lstate,pqspin,astate,kstate,pqspin)
       num2 = eri_eigen(kstate,istate,pqspin,lstate,jstate,pqspin)
       num3 = eri_eigen(qstate,istate,pqspin,astate,jstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(jstate,pqspin) - ieta_tmp
         denom2 = omega + energy(astate,pqspin) - energy(kstate,pqspin) - energy(lstate,pqspin) - ieta_tmp
         ! Minus sign from Domcke-Cederbaum book chapter 1977 (forgotten in von niessen review in 1983)
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) - num1 * num2 * num3 / ( denom1 * denom2 ) 
       enddo
     enddo
   enddo
   enddo
   enddo
   enddo

   ! C2+C3   i,j,k   a,b 
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       num1 = 2.0_dp * eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)
       num2 = eri_eigen(astate,jstate,pqspin,bstate,kstate,pqspin)
       num3 = eri_eigen(qstate,jstate,pqspin,istate,kstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) - ieta
         denom2 = energy(jstate,pqspin) + energy(kstate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) + 2.0_dp * num1 * num2 * num3 / ( denom1 * denom2 )
       enddo
     enddo
     enddo
   enddo
   enddo
   enddo

   ! C4+C5   i,j   a,b,c
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
     do cstate=nhomo_G+1,nvirtual_G-1
       num1 = 2.0_dp * eri_eigen(pstate,istate,pqspin,astate,jstate,pqspin) - eri_eigen(pstate,jstate,pqspin,astate,istate,pqspin)
       num2 = eri_eigen(istate,bstate,pqspin,jstate,cstate,pqspin)
       num3 = eri_eigen(qstate,bstate,pqspin,astate,cstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(jstate,pqspin) + ieta
         denom2 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(bstate,pqspin) - energy(cstate,pqspin)
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) + num1 * num2 * num3 / ( denom1 * denom2 ) 
       enddo
     enddo
     enddo
     enddo
   enddo
   enddo

   ! D1   i,j   a,b,c
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
     do cstate=nhomo_G+1,nvirtual_G-1
       num1a = eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin)
       num1b = eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)
       num2a = eri_eigen(astate,istate,pqspin,jstate,cstate,pqspin)
       num2b = eri_eigen(astate,cstate,pqspin,jstate,istate,pqspin)
       num3a = eri_eigen(qstate,cstate,pqspin,jstate,bstate,pqspin) - 2.0_dp * eri_eigen(qstate,bstate,pqspin,jstate,cstate,pqspin)
       num3b = eri_eigen(qstate,bstate,pqspin,jstate,cstate,pqspin) - 2.0_dp * eri_eigen(qstate,cstate,pqspin,jstate,bstate,pqspin)
       numgw = 4.0_dp * eri_eigen(qstate,bstate,pqspin,jstate,cstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta_tmp
         denom2 = omega + energy(jstate,pqspin) - energy(bstate,pqspin) - energy(cstate,pqspin) + ieta_tmp
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) &
                    + ( num1a * ( num2a * num3a + num2b * num3b ) + num1b * ( num2a * (-2.0_dp)*num3a + num2b * num3a ) )   / ( denom1 * denom2 )
         selfenergy_2rings(iomega,pstate,pqspin) = selfenergy_2rings(iomega,pstate,pqspin) &
                    + ( num1b * num2a * numgw ) / ( denom1 * denom2 )
       enddo
     enddo
     enddo
     enddo
   enddo
   enddo


   ! D2+D3   i,j   a,b,c
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
     do cstate=nhomo_G+1,nvirtual_G-1
       num1a = eri_eigen(pstate,cstate,pqspin,istate,astate,pqspin)
       num1b = eri_eigen(pstate,astate,pqspin,istate,cstate,pqspin)
       num2a = eri_eigen(astate,istate,pqspin,bstate,jstate,pqspin)
       num2b = eri_eigen(astate,jstate,pqspin,bstate,istate,pqspin)
       num3a = eri_eigen(qstate,jstate,pqspin,bstate,cstate,pqspin) - 2.0_dp * eri_eigen(qstate,cstate,pqspin,bstate,jstate,pqspin)
       num3b = eri_eigen(qstate,cstate,pqspin,bstate,jstate,pqspin) - 2.0_dp * eri_eigen(qstate,jstate,pqspin,bstate,cstate,pqspin)
       numgw = 4.0_dp * eri_eigen(qstate,cstate,pqspin,bstate,jstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(cstate,pqspin) + ieta
         denom2 = energy(jstate,pqspin) + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) &
                    + 2.0_dp * ( num1a * ( num2a * -(2.0_dp)*num3a + num2b * num3a ) + num1b * ( num2a * num3b + num2b * num3a ) )   / ( denom1 * denom2 )
         selfenergygw(iomega,pstate,pqspin) = selfenergygw(iomega,pstate,pqspin) &
                    + 2.0_dp * num1a * num2a * numgw  / ( denom1 * denom2 )
       enddo
     enddo
     enddo
     enddo
   enddo
   enddo


   ! D4+D5   i,j,k   a,b
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       num1a = eri_eigen(pstate,kstate,pqspin,astate,jstate,pqspin)
       num1b = eri_eigen(pstate,jstate,pqspin,astate,kstate,pqspin)
       num2a = eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin)
       num2b = eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
       num3a = eri_eigen(qstate,bstate,pqspin,istate,kstate,pqspin) - 2.0_dp * eri_eigen(qstate,kstate,pqspin,istate,bstate,pqspin)
       num3b = eri_eigen(qstate,kstate,pqspin,istate,bstate,pqspin) - 2.0_dp * eri_eigen(qstate,bstate,pqspin,istate,kstate,pqspin)
       numgw = 4.0_dp * eri_eigen(qstate,kstate,pqspin,istate,bstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(astate,pqspin) - energy(jstate,pqspin) - energy(kstate,pqspin) - ieta
         denom2 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) &
                    + 2.0_dp * ( num1a * ( num2a * -(2.0_dp)*num3a + num2b * num3a ) + num1b * ( num2a * num3a + num2b * num3b ) )   / ( denom1 * denom2 )
         selfenergygw(iomega,pstate,pqspin) = selfenergygw(iomega,pstate,pqspin) &
                    + 2.0_dp * num1a * num2a * numgw  / ( denom1 * denom2 )
       enddo
     enddo
     enddo
   enddo
   enddo
   enddo

   ! D6   i,j,k   a,b
   do istate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       num1a = eri_eigen(pstate,kstate,pqspin,astate,istate,pqspin)
       num1b = eri_eigen(pstate,istate,pqspin,astate,kstate,pqspin)
       num2a = eri_eigen(istate,astate,pqspin,bstate,jstate,pqspin)
       num2b = eri_eigen(istate,jstate,pqspin,bstate,astate,pqspin)
       num3a = eri_eigen(qstate,jstate,pqspin,bstate,kstate,pqspin) - 2.0_dp * eri_eigen(qstate,kstate,pqspin,bstate,jstate,pqspin)
       num3b = eri_eigen(qstate,kstate,pqspin,bstate,jstate,pqspin) - 2.0_dp * eri_eigen(qstate,jstate,pqspin,bstate,kstate,pqspin)
       numgw = 4.0_dp * eri_eigen(qstate,kstate,pqspin,bstate,jstate,pqspin)
       do iomega=-se%nomega,se%nomega
         omega = energy(qstate,pqspin) + se%omega(iomega)
         denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(kstate,pqspin) - ieta_tmp
         denom2 = omega + energy(bstate,pqspin) - energy(jstate,pqspin) - energy(kstate,pqspin) - ieta_tmp
         selfenergy(iomega,pstate,pqspin) = selfenergy(iomega,pstate,pqspin) &
                    -( num1a * ( num2a * -(2.0_dp)*num3a + num2b * num3a ) + num1b * ( num2a * num3a + num2b * num3b ) )   / ( denom1 * denom2 )
         selfenergy_2rings(iomega,pstate,pqspin) = selfenergy_2rings(iomega,pstate,pqspin) &
                    -( num1a * num2a * numgw )   / ( denom1 * denom2 )
       enddo
     enddo
     enddo
   enddo
   enddo
   enddo


   write(*,'(i4,*(7x,f14.6))') pstate,REAL(selfenergy(0,pstate,:),dp) * Ha_eV,REAL(selfenergy2(0,pstate,pqspin),dp) * Ha_eV,REAL(selfenergygw(0,pstate,pqspin),dp) * Ha_eV
 enddo

 se%sigma(:,:,:) = selfenergy2(:,:,:) + selfenergy(:,:,:) ! selfenergy_2rings(:,:,:) ! selfenergy2(:,:,:) + selfenergy(:,:,:) ! FBFB

 deallocate(selfenergy)
 deallocate(selfenergy2)
 deallocate(selfenergygw)
 deallocate(selfenergy_2rings)


end subroutine pt3_selfenergy


!=========================================================================
