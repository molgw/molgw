!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the GW self-energy with vertex function
! within different flavors: G0W0GAMMA0 G0W0SOX0
!
!=========================================================================
subroutine gwgamma_selfenergy(nstate,basis,occupation,energy,c_matrix,wpol,se)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_timing
 use m_inputparam
 use m_warning
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_selfenergy_tools
 use m_tddft_fxc
 implicit none

 integer,intent(in)                 :: nstate
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer                 :: iomega
 complex(dp),allocatable :: sigma_gw(:,:,:)
 complex(dp),allocatable :: sigma_gamma(:,:,:)
 complex(dp),allocatable :: sigma_sox(:,:,:)
 integer                 :: astate,bstate,cstate
 integer                 :: istate,jstate,kstate,ispin,spole
 integer                 :: mstate
 real(dp),allocatable    :: bra(:,:)
 real(dp)                :: vcoul,vcoul1,vcoul2
 real(dp)                :: pole_s
 real(dp)                :: fxc
!=====

 call start_clock(timing_gwgamma_self)

 write(stdout,*)
 select case(calc_type%selfenergy_approx)
 case(G0W0SOX0)
   write(stdout,*) 'Perform a one-shot G0W0SOX0 calculation'
 case(G0W0GAMMA0)
   if( calc_type%selfenergy_technique == EVSC ) then
     write(stdout,*) 'Perform an eigenvalue-self-consistent GWGAMMA calculation'
   else
     write(stdout,*) 'Perform a one-shot G0W0GAMMA0 calculation'
   endif
 case default
   call die('gwgamma_selfenergy: calculation type unknown')
 end select


 if( gwgamma_tddft_ ) then
   write(stdout,*) 'Include a TDDFT kernel contribution to the vertex'
   write(stdout,'(1x,a,f12.4)') 'Exact-exchange amount: ',alpha_hybrid
   call prepare_tddft(nstate,basis,c_matrix,occupation)
 endif

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif


 call clean_allocate('Temporary array',bra,ncore_G+1,nvirtual_G-1,ncore_G+1,MAX(nhomo_G,nsemax))



 !
 !
 allocate(sigma_gamma(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(sigma_sox(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 allocate(sigma_gw(-se%nomega:se%nomega,nsemin:nsemax,nspin))

 sigma_gamma(:,:,:)  = 0.0_dp
 sigma_sox(:,:,:)  = 0.0_dp


 write(stdout,*) 'Calculate SOX'

 do ispin=1,nspin

   !==========================
   do bstate=ncore_G+1,nvirtual_G-1
     if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty) cycle
     if( MODULO( bstate-(ncore_G+1) , nproc_ortho ) /= rank_ortho ) cycle

     do istate=ncore_G+1,nvirtual_G-1
       if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
       do kstate=ncore_G+1,nvirtual_G-1
         if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle

         do mstate=nsemin,nsemax

           vcoul1 = eri_eigen(mstate,istate,ispin,bstate,kstate,ispin)
           vcoul2 = eri_eigen(istate,bstate,ispin,kstate,mstate,ispin)
           if( gwgamma_tddft_ ) then
             fxc = eval_fxc_rks_singlet(istate,bstate,ispin,kstate,mstate,ispin)
             call xsum_grid(fxc)
             vcoul2 = alpha_hybrid * vcoul2 - fxc

!             if( ABS( eri_eigen(istate,bstate,ispin,kstate,mstate,ispin) -vcoul2)> 0.10 ) then
!               write(*,'(4(i4,1x),4(1x,f12.6))') istate,bstate,kstate,mstate, &
!                  eri_eigen(istate,bstate,ispin,kstate,mstate,ispin), &
!                  vcoul2
!               write(*,*) 'Hack'
!               vcoul2 = eri_eigen(istate,bstate,ispin,kstate,mstate,ispin)
!             endif

           endif
           !
           ! calculate only the diagonal !
           do iomega=-se%nomega,se%nomega
             sigma_sox(iomega,mstate,ispin) = sigma_sox(iomega,mstate,ispin) &
                 - vcoul1 * vcoul2            &
                   / ( energy(mstate,ispin) + se%omega(iomega) - energy(istate,ispin) - energy(kstate,ispin) + energy(bstate,ispin) - ieta )
           enddo
         enddo

       enddo
     enddo
   enddo

   !==========================
   do cstate=ncore_G+1,nvirtual_G-1
     if( (spin_fact - occupation(cstate,ispin)) / spin_fact < completely_empty) cycle
     if( MODULO( cstate-(ncore_G+1) , nproc_ortho ) /= rank_ortho ) cycle

     do jstate=ncore_G+1,nvirtual_G-1
       if( occupation(jstate,ispin) / spin_fact < completely_empty ) cycle
       do astate=ncore_G+1,nvirtual_G-1
         if( (spin_fact - occupation(astate,ispin)) / spin_fact < completely_empty) cycle

         do mstate=nsemin,nsemax

           vcoul1 = eri_eigen(mstate,astate,ispin,jstate,cstate,ispin)
           vcoul2 = eri_eigen(astate,jstate,ispin,cstate,mstate,ispin)
           if( gwgamma_tddft_ ) then
             fxc = eval_fxc_rks_singlet(astate,jstate,ispin,cstate,mstate,ispin)
             call xsum_grid(fxc)
             vcoul2 = alpha_hybrid * vcoul2 - fxc

!             if( ABS( eri_eigen(astate,jstate,ispin,cstate,mstate,ispin) -vcoul2 )> 0.10 ) then
!               write(*,'(4(i4,1x),4(1x,f12.6))') astate,jstate,cstate,mstate, &
!                  eri_eigen(astate,jstate,ispin,cstate,mstate,ispin), &
!                  vcoul2
!!               write(*,*) 'Hack'
!!               vcoul2 =  eri_eigen(astate,jstate,ispin,cstate,mstate,ispin)
!             endif

           endif
           !
           ! calculate only the diagonal !
           do iomega=-se%nomega,se%nomega
             sigma_sox(iomega,mstate,ispin) = sigma_sox(iomega,mstate,ispin) &
                 - vcoul1 * vcoul2            &
                   / ( energy(mstate,ispin) + se%omega(iomega) - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) + ieta )
           enddo
         enddo

       enddo
     enddo
   enddo


 enddo

 call xsum_ortho(sigma_sox)


 if( calc_type%selfenergy_approx == G0W0GAMMA0 ) then

   write(stdout,*) 'Calculate dynamical SOSEX'


   do ispin=1,nspin

     do spole=1,wpol%npole_reso

       if( MODULO( spole - 1 , nproc_ortho ) /= rank_ortho ) cycle
       write(stdout,*) 'SOSEX W poles:',spole,' / ',wpol%npole_reso

       pole_s = wpol%pole(spole)

       if(has_auxil_basis) then
         do mstate=ncore_G+1,MAX(nhomo_G,nsemax)
           ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
           bra(:,mstate)     = MATMUL( wpol%residue_left(:,spole) , eri_3center_eigen(:,:,mstate,ispin) )
         enddo
         call xsum_auxil(bra)
       else
         ! Here just grab the precalculated value
         forall(istate=ncore_G+1:nvirtual_G-1, mstate=ncore_G+1:MAX(nhomo_G,nsemax))
           bra(istate,mstate) = wpol%residue_left(index_prodstate(istate,mstate) + (ispin-1) * index_prodstate(nvirtual_W-1,nvirtual_W-1), &
                                                  spole)
         end forall
       endif


       !==========================
       do istate=ncore_G+1,nvirtual_G-1
         if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
         do bstate=ncore_G+1,nvirtual_G-1
           if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty) cycle
           do kstate=ncore_G+1,nvirtual_G-1
             if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle

             !
             ! calculate only the diagonal !
             do mstate=nsemin,nsemax

               vcoul = eri_eigen(istate,kstate,ispin,bstate,mstate,ispin)
               if( gwgamma_tddft_ ) then
                 fxc = eval_fxc_rks_singlet(istate,kstate,ispin,bstate,mstate,ispin)
                 call xsum_grid(fxc)
                 vcoul = alpha_hybrid * vcoul - fxc
               endif

               do iomega=-se%nomega,se%nomega
                 sigma_gamma(iomega,mstate,ispin) = sigma_gamma(iomega,mstate,ispin) &
                          - bra(kstate,mstate) * bra(bstate,istate) * vcoul                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + pole_s - ieta )  &
                             / ( -pole_s + energy(istate,ispin) - energy(bstate,ispin) + ieta )
               enddo
             enddo

           enddo
         enddo
       enddo

       !==========================
       do istate=ncore_G+1,nvirtual_G-1
         if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
         do bstate=ncore_G+1,nvirtual_G-1
           if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty ) cycle
           do cstate=ncore_G+1,nvirtual_G-1
             if( (spin_fact - occupation(cstate,ispin)) / spin_fact < completely_empty ) cycle

             !
             ! calculate only the diagonal !
             do mstate=nsemin,nsemax

               vcoul = eri_eigen(istate,cstate,ispin,bstate,mstate,ispin)
               if( gwgamma_tddft_ ) then
                 fxc = eval_fxc_rks_singlet(istate,cstate,ispin,bstate,mstate,ispin)
                 call xsum_grid(fxc)
                 vcoul = alpha_hybrid * vcoul - fxc
               endif

               do iomega=-se%nomega,se%nomega
                 sigma_gamma(iomega,mstate,ispin) = sigma_gamma(iomega,mstate,ispin) &
                          - bra(cstate,mstate) * bra(bstate,istate) * vcoul                          &
                            / ( energy(mstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )    &
                            / ( energy(mstate,ispin) + se%omega(iomega) - energy(cstate,ispin) + energy(istate,ispin) - energy(bstate,ispin) + ieta )


                 sigma_gamma(iomega,mstate,ispin) = sigma_gamma(iomega,mstate,ispin) &
                          + bra(cstate,mstate) * bra(bstate,istate) * vcoul                          &
                            / ( energy(mstate,ispin) + se%omega(iomega) - energy(bstate,ispin) - energy(cstate,ispin) + energy(istate,ispin) + ieta )  &
                            / ( energy(bstate,ispin) - energy(istate,ispin) + pole_s - ieta )

               enddo
             enddo

           enddo
         enddo
       enddo

       !==========================
       do astate=ncore_G+1,nvirtual_G-1
         if( (spin_fact - occupation(astate,ispin)) / spin_fact < completely_empty  ) cycle
         do jstate=ncore_G+1,nvirtual_G-1
           if( occupation(jstate,ispin) / spin_fact < completely_empty ) cycle
           do kstate=ncore_G+1,nvirtual_G-1
             if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle

             !
             ! calculate only the diagonal !
             do mstate=nsemin,nsemax

               vcoul = eri_eigen(astate,kstate,ispin,jstate,mstate,ispin)
               if( gwgamma_tddft_ ) then
                 fxc = eval_fxc_rks_singlet(astate,kstate,ispin,jstate,mstate,ispin)
                 call xsum_grid(fxc)
                 vcoul = alpha_hybrid * vcoul - fxc
               endif

               do iomega=-se%nomega,se%nomega
                 sigma_gamma(iomega,mstate,ispin) = sigma_gamma(iomega,mstate,ispin) &
                          - bra(kstate,mstate) * bra(astate,jstate) * vcoul                          &
                            / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                            / ( energy(jstate,ispin) - energy(astate,ispin) - pole_s + ieta )

                 sigma_gamma(iomega,mstate,ispin) = sigma_gamma(iomega,mstate,ispin) &
                          + bra(kstate,mstate) * bra(astate,jstate) * vcoul                          &
                            / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                            / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + pole_s - ieta )


               enddo
             enddo

           enddo
         enddo
       enddo

       !==========================
       do astate=ncore_G+1,nvirtual_G-1
         if( (spin_fact - occupation(astate,ispin)) / spin_fact < completely_empty  ) cycle
         do jstate=ncore_G+1,nvirtual_G-1
           if( occupation(jstate,ispin) / spin_fact < completely_empty ) cycle
           do cstate=ncore_G+1,nvirtual_G-1
             if( (spin_fact - occupation(cstate,ispin)) / spin_fact < completely_empty ) cycle

             !
             ! calculate only the diagonal !
             do mstate=nsemin,nsemax

               vcoul = eri_eigen(astate,cstate,ispin,jstate,mstate,ispin)
               if( gwgamma_tddft_ ) then
                 fxc = eval_fxc_rks_singlet(astate,cstate,ispin,jstate,mstate,ispin)
                 call xsum_grid(fxc)
                 vcoul = alpha_hybrid * vcoul - fxc
               endif

               do iomega=-se%nomega,se%nomega
                 sigma_gamma(iomega,mstate,ispin) = sigma_gamma(iomega,mstate,ispin) &
                          + bra(cstate,mstate) * bra(astate,jstate) * vcoul                          &
                            / ( energy(mstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )  &
                            / ( pole_s + energy(astate,ispin) - energy(jstate,ispin) - ieta )

               enddo
             enddo

           enddo
         enddo
       enddo



     enddo !spole
   enddo !ispin

   call xsum_ortho(sigma_gamma)

 endif


 write(stdout,'(a)') ' Sigma_c(omega) is calculated'


 !
 ! The input sigma contains the GW selfenergy
 sigma_gw(:,:,:) = se%sigma(:,:,:)


 forall(astate=nsemin:nsemax)
   se%sigma(:,astate,:) = sigma_gw(:,astate,:) + sigma_sox(:,astate,:) + sigma_gamma(:,astate,:)
 end forall


! if( print_sigma_) then
!   call write_selfenergy_omega('selfenergy_sox'    ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sox)
!   call write_selfenergy_omega('selfenergy_gamma'  ,energy,exchange_m_vxc_diag,occupation,energy,sigma_gamma)
!   call write_selfenergy_omega('selfenergy_gwgamma',energy,exchange_m_vxc_diag,occupation,energy,sigma)
! endif


 write(stdout,'(/,a)') ' G0W0Gamma0 self-energy contributions at E0 (eV)'
 if(nspin==1) then
   write(stdout,'(a)') '   #          E0             SigC_G0W0                 SigC_SOX                SigC_Gamma0                SigC_TOT'
 else
   write(stdout,'(a)') &
     '   #                E0                              SigC_G0W0            SigC_SOX             SigC_Gamma0                SigC_TOT'
 endif

 do astate=nsemin,nsemax
   write(stdout,'(i4,1x,20(1x,f12.6))') astate,energy(astate,:)*Ha_eV,          &
                                        sigma_gw(0,astate,:)*Ha_eV,   &
                                        sigma_sox(0,astate,:)*Ha_eV,  &
                                        sigma_gamma(0,astate,:)*Ha_eV,&
                                        se%sigma(0,astate,:)*Ha_eV
 enddo




 call clean_deallocate('Temporary array',bra)

 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 else
   call destroy_eri_4center_eigen_uks()
 endif

 if( gwgamma_tddft_ ) then
   call destroy_tddft()
 endif

 call stop_clock(timing_gwgamma_self)


end subroutine gwgamma_selfenergy


!=========================================================================
