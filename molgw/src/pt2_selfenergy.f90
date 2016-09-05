!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the perturbation theory to 2nd order evaluation of the self-energy
!
!=========================================================================
subroutine pt2_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,selfenergy_omega,emp2)
 use m_definitions
 use m_mpi
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: selfenergy_approx,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)       :: selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nspin)
 real(dp),intent(out)       :: emp2
!=====
 integer               :: pstate,qstate
 real(dp),allocatable  :: selfenergy_ring(:,:,:)
 real(dp),allocatable  :: selfenergy_sox(:,:,:)
 integer               :: iomegai
 integer               :: istate,jstate,kstate
 integer               :: pqispin,jkspin
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek
 real(dp)              :: omega
 real(dp)              :: fact_real,fact_energy
 real(dp)              :: emp2_sox,emp2_ring
 real(dp),allocatable  :: eri_eigenstate_i(:,:,:,:)
 real(dp)              :: coul_iqjk,coul_ijkq,coul_ipkj
!=====

 call start_clock(timing_mp2_self)

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp


 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 write(stdout,*) 'with the perturbative approach'

 

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
 endif



 allocate(selfenergy_ring(-nomegai:nomegai,nsemin:nsemax,nspin))
 allocate(selfenergy_sox (-nomegai:nomegai,nsemin:nsemax,nspin))


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

             do iomegai=-nomegai,nomegai
               omega = energy(qstate,pqispin) + omegai(iomegai)

               fact_real   = REAL( fact_occ1 / (omega-ei+ej-ek+ieta) + fact_occ2 / (omega-ei+ej-ek-ieta) , dp)
               fact_energy = REAL( fact_occ1 / (energy(pstate,pqispin)-ei+ej-ek+ieta) , dp )

               selfenergy_ring(iomegai,pstate,pqispin) = selfenergy_ring(iomegai,pstate,pqispin) &
                        + fact_real * coul_ipkj * coul_iqjk * spin_fact

               if(iomegai==0 .AND. occupation(pstate,pqispin)>completely_empty) then
                 emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                       * fact_energy * coul_ipkj * coul_iqjk * spin_fact
               endif

               if( pqispin == jkspin ) then

                 selfenergy_sox(iomegai,pstate,pqispin) = selfenergy_sox(iomegai,pstate,pqispin) &
                          - fact_real * coul_ipkj * coul_ijkq

                 if(iomegai==0 .AND. occupation(pstate,pqispin)>completely_empty) then
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

 selfenergy_omega(:,:,:) = selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:)


 if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)
 if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_mp2_self)

end subroutine pt2_selfenergy


!=========================================================================
subroutine onering_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,selfenergy_omega,emp2)
 use m_definitions
 use m_mpi
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_spectral_function
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: selfenergy_approx,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)       :: selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nspin)
 real(dp),intent(out)       :: emp2
!=====
 integer               :: pstate,qstate
 integer               :: iomegai
 integer               :: istate,jstate,kstate
 integer               :: pqispin,jkspin
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek
 real(dp)              :: omega
 real(dp)              :: fact_real,fact_energy
 real(dp)              :: emp2_ring
 real(dp)              :: coul_iqjk,coul_ijkq,coul_ipkj

 type(spectral_function) :: vsqrtchi0vsqrt
 integer                 :: bstate,jbspin,t_jb
!=====

 call start_clock(timing_mp2_self)

 if( .NOT. has_auxil_basis ) &
   call die('onering_selfenergy: only implemented when an auxiliary basis is available')

 emp2_ring = 0.0_dp


 write(stdout,'(/,a)') ' Perform the one-ring self-energy calculation'
 write(stdout,*) 'with the perturbative approach'

 
 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)

 call init_spectral_function(nstate,occupation,vsqrtchi0vsqrt)
 call allocate_spectral_function(nauxil_3center,vsqrtchi0vsqrt)

 do t_jb=1,vsqrtchi0vsqrt%npole_reso
   jstate = vsqrtchi0vsqrt%transition_table_apb(1,t_jb)
   bstate = vsqrtchi0vsqrt%transition_table_apb(2,t_jb)
   jbspin = vsqrtchi0vsqrt%transition_table_apb(3,t_jb)

   vsqrtchi0vsqrt%residu_left(:,t_jb) = eri_3center_eigen(:,jstate,bstate,jbspin) * SQRT(spin_fact)
   vsqrtchi0vsqrt%pole(t_jb)          = energy(bstate,jbspin) - energy(jstate,jbspin)

 end do

 call destroy_eri_3center_eigen()

 call gw_selfenergy(GW,nstate,basis,occupation,energy,c_matrix,vsqrtchi0vsqrt,selfenergy_omega,emp2_ring)
 
 call destroy_spectral_function(vsqrtchi0vsqrt)

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
 integer               :: pstate,qstate,qstate2
 real(dp),allocatable  :: selfenergy_ring(:,:,:)
 real(dp),allocatable  :: selfenergy_sox(:,:,:)
 integer               :: istate,jstate,kstate
 integer               :: pqispin,jkspin
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek,ep,eq
 real(dp)              :: fact_real,fact_energy
 real(dp)              :: emp2_sox,emp2_ring
 real(dp),allocatable  :: eri_eigenstate_i(:,:,:,:)
 real(dp)              :: coul_iqjk,coul_ijkq,coul_ipkj
!=====

 call start_clock(timing_mp2_self)

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp


 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 write(stdout,*) 'with the QP self-consistent approach'

 

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
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
  
               fact_real   = REAL( fact_occ1 / ( eq - ei + ej - ek + ieta) &
                                 + fact_occ2 / ( eq - ei + ej - ek - ieta) , dp)
               fact_energy = REAL( fact_occ1 / ( ep - ei + ej - ek + ieta) , dp )
  
               selfenergy_ring(pstate,qstate,pqispin) = selfenergy_ring(pstate,qstate,pqispin) &
                        + fact_real * coul_ipkj * coul_iqjk * spin_fact

               if(pstate==qstate .AND. occupation(pstate,pqispin)>completely_empty) then
                 emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                       * fact_energy * coul_ipkj * coul_iqjk * spin_fact
               endif
 
               if( pqispin == jkspin ) then

                 selfenergy_sox(pstate,qstate,pqispin) = selfenergy_sox(pstate,qstate,pqispin) &
                          - fact_real * coul_ipkj * coul_ijkq

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


 selfenergy(:,:,:) = selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:)

 call apply_qs_approximation(basis%nbf,nstate,s_matrix,c_matrix,selfenergy)


 if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)
 if(has_auxil_basis) call destroy_eri_3center_eigen()

 call stop_clock(timing_mp2_self)

end subroutine pt2_selfenergy_qs


!=========================================================================
