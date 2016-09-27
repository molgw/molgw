!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the COHSEX self-energy
!
!=========================================================================
subroutine cohsex_selfenergy(nstate,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,wpol,se)
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam
 use m_warning,only: issue_warning,msg
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_selfenergy_tools
 implicit none

 integer,intent(in)                  :: nstate
 type(basis_set)                     :: basis
 real(dp),intent(in)                 :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)                 :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in)  :: wpol
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer               :: pstate
 integer               :: istate,ipspin
 real(dp)              :: fact_full_i,fact_empty_i
 integer               :: jbf_auxil
 integer               :: ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable  :: wp0(:,:),wp0_i(:),w0_local(:)
 real(dp)              :: sigx
!=====

 call start_clock(timing_self)

 if( .NOT. has_auxil_basis )    call die('cohsex: no RI is not coded')
 if( .NOT. ALLOCATED(wpol%w0) ) call die('cohsex: static W should be available here')


 select case(calc_type%selfenergy_approx)
 case(COHSEX_DEVEL)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case(TUNED_COHSEX)
   write(stdout,*) 'Perform a tuned COHSEX calculation'
 case default
   call die('cohsex_selfenergy: calculation type unknown')
 end select


 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,nsemin,nsemax,ncore_G+1,nvirtual_G-1)

 write(stdout,*) '=============='
 write(stdout,*) 'FBFB exchange'
 do ipspin=1,nspin
   do pstate=nsemin,nsemax
     sigx = 0.0_dp
     do istate=ncore_G+1,nvirtual_G-1
       fact_full_i   = occupation(istate,ipspin) / spin_fact
       if( fact_full_i < completely_empty ) cycle

       sigx = sigx - SUM( eri_3center_eigen(:,pstate,istate,ipspin)**2 ) * fact_full_i

     enddo
     call xsum_auxil(sigx)

     write(stdout,'(i4,4x,f16.6)') pstate,sigx * Ha_eV

     ! Store the result in se%sigma
     se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) + sigx * epsilon_cohsex

   enddo
 enddo
 write(stdout,*) '=============='




 se%sigma(:,:,:)  = 0.0_dp


 allocate(wp0(nauxil_3center,nsemin:nsemax))

 do ipspin=1,nspin
   !
   ! Apply the frozen core and frozen virtual approximation to G
   do istate=ncore_G+1,nvirtual_G-1

     if( MODULO( istate -(ncore_G+1) , nproc_ortho) /= rank_ortho ) cycle

     !
     ! Prepare the right hand side
     allocate(wp0_i(nsemin:nsemax))
     allocate(w0_local(nauxil_3center))

     do ibf_auxil_global=1,nauxil_2center

       do jbf_auxil=1,nauxil_3center
         jbf_auxil_global = ibf_auxil_g(jbf_auxil)
         w0_local(jbf_auxil) = wpol%w0(ibf_auxil_global,jbf_auxil_global)
       enddo

       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (sqrt(v) * chi * v)
       wp0_i(nsemin:nsemax) = MATMUL( w0_local(:) , eri_3center_eigen(:,nsemin:nsemax,istate,ipspin) )
       call xsum_auxil(wp0_i)

       if( iproc_ibf_auxil(ibf_auxil_global) == rank_auxil ) then
         wp0(ibf_auxil_l(ibf_auxil_global),:) = wp0_i(:)
       endif

     enddo
     deallocate(w0_local,wp0_i)





     ! The application of residue theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ipspin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ipspin)) / spin_fact



     select case(calc_type%selfenergy_approx)

     case(COHSEX_DEVEL)

       do pstate=nsemin,nsemax
         !
         ! SEX
         !
         se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) &
                    -  DOT_PRODUCT( eri_3center_eigen(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * fact_full_i * 1.0_dp  &
                          * beta_cohsex

         !
         ! COH
         !
         se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) &
                    +  DOT_PRODUCT( eri_3center_eigen(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * alpha_cohsex * 0.5_dp

       enddo

     case(TUNED_COHSEX)

       do pstate=nsemin,nsemax
         !
         ! SEX
         !
         se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) &
                    -  DOT_PRODUCT( eri_3center_eigen(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * fact_full_i * 1.0_dp  &
                          * ( beta_cohsex + gamma_cohsex )
         !
         ! No full range COH in tuned-COHSEX
         !
       enddo



     case default 
       call die('BUG')
     end select

   enddo !istate
 enddo !ipspin

 deallocate(wp0)

 ! Sum up the contribution from different procs
 call xsum_ortho(se%sigma)


 write(stdout,'(a)') ' Sigma_c is calculated'


 select case(calc_type%selfenergy_approx)
 case(TUNED_COHSEX) !==========================================================

   ! This coding has been corrupted !
   ! FIXME
   call die('the following line is obviously wrong')
   se%sigma(0,nsemin:nsemax,:) = se%sigma(0,nsemin:nsemax,:) + se%sigma(0,nsemin:nsemax,:)


 end select


 call destroy_eri_3center_eigen()

 call stop_clock(timing_self)


end subroutine cohsex_selfenergy


!=========================================================================
subroutine cohsex_selfenergy_lr(nstate,basis,occupation,energy,exchange_m_vxc_diag, &
                                c_matrix,wpol,se)
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam
 use m_warning,only: issue_warning,msg
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_eri
 use m_eri_calculate
 use m_eri_calculate_lr
 use m_selfenergy_tools
 implicit none

 integer,intent(in)                 :: nstate
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer               :: homo
 integer               :: pstate
 integer               :: istate,ipspin
 real(dp)              :: fact_full_i,fact_empty_i
 integer               :: jbf_auxil
 integer               :: ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable  :: wp0(:,:),wp0_i(:),w0_local(:)
 real(dp),allocatable  :: wp0_tmp(:,:),wp0_rotation(:,:)
 integer               :: nbf_auxil
 real(dp)              :: sigx
!=====

 call start_clock(timing_self)

 if( .NOT. has_auxil_basis )    call die('cohsex: no RI is not coded')
 if( .NOT. ALLOCATED(wpol%w0) ) call die('cohsex: static W should be available here')

 call assert_experimental()


 call calculate_eri_3center_eigen_lr(basis%nbf,nstate,c_matrix)


 write(stdout,*)
 select case(calc_type%selfenergy_approx)
 case(COHSEX)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case(TUNED_COHSEX)
   write(stdout,*) 'Perform a tuned COHSEX calculation'
 end select


 ! Rotation of W0

#ifdef COHSEX_DEVEL
 nbf_auxil    = SIZE( eri_2center_rotation   (:,:) , DIM=1)

 allocate( wp0_tmp(nbf_auxil,nbf_auxil) )
 allocate( wp0_rotation(nauxil_2center_lr,nauxil_2center_lr) )
 wp0_tmp(:,:)      = MATMUL( eri_2center_rotation(:,:) , MATMUL( wpol%w0(:,:) , TRANSPOSE( eri_2center_rotation(:,:) ) ) )
 wp0_rotation(:,:) = MATMUL( TRANSPOSE(eri_2center_rotation_lr(:,:)) , MATMUL( wp0_tmp(:,:) , eri_2center_rotation_lr(:,:) ) )
 deallocate( wp0_tmp )
#endif



 write(stdout,*) '=============='
 write(stdout,*) 'FBFB exchange'
 do ipspin=1,nspin
   do pstate=nsemin,nsemax
     sigx = 0.0_dp
     do istate=ncore_G+1,nstate
       fact_full_i   = occupation(istate,ipspin) / spin_fact
       if( fact_full_i < completely_empty ) cycle

       sigx = sigx - SUM( eri_3center_eigen_lr(:,pstate,istate,ipspin)**2 ) * fact_full_i

     enddo
     call xsum_auxil(sigx)

     write(stdout,'(i4,4x,f16.6)') pstate,sigx * Ha_eV

     ! Store the result in se%sigma
     se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) - sigx * epsilon_cohsex

   enddo
 enddo
 write(stdout,*) '=============='





 se%sigma(:,:,:)  = 0.0_dp


 allocate(wp0(nauxil_3center_lr,nsemin:nsemax))

 do ipspin=1,nspin
   !
   ! Apply the frozen core and frozen virtual approximation to G
   do istate=ncore_G+1,nvirtual_G-1

     if( MODULO( istate -(ncore_G+1) , nproc_ortho) /= rank_ortho ) cycle

     !
     ! Prepare the right hand side
     allocate(wp0_i(nsemin:nsemax))
     allocate(w0_local(nauxil_3center_lr))

     do ibf_auxil_global=1,nauxil_2center_lr

       do jbf_auxil=1,nauxil_3center_lr
         jbf_auxil_global = ibf_auxil_g(jbf_auxil)
         w0_local(jbf_auxil) =  wp0_rotation(ibf_auxil_global,jbf_auxil_global) ! wpol%w0(ibf_auxil_global,jbf_auxil_global)
       enddo

       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (sqrt(v) * chi * v)
       wp0_i(nsemin:nsemax) = MATMUL( w0_local(:) , eri_3center_eigen_lr(:,nsemin:nsemax,istate,ipspin) )
       call xsum_auxil(wp0_i)

       if( iproc_ibf_auxil(ibf_auxil_global) == rank_auxil ) then
         wp0(ibf_auxil_l(ibf_auxil_global),:) = wp0_i(:)
       endif

     enddo
     deallocate(w0_local,wp0_i)





     ! The application of residue theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ipspin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ipspin)) / spin_fact



     select case(calc_type%selfenergy_approx)

     case(COHSEX_DEVEL)

       do pstate=nsemin,nsemax
         !
         ! SEX
         !
         se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) &
                    -  DOT_PRODUCT( eri_3center_eigen_lr(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * fact_full_i * 1.0_dp  &
                          * beta_cohsex

         !
         ! COH
         !
         se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) &
                    +  DOT_PRODUCT( eri_3center_eigen_lr(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * alpha_cohsex * 0.5_dp

       enddo

     case(TUNED_COHSEX)

       do pstate=nsemin,nsemax
         !
         ! LR-SEX 
         !
         se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) &
                    -  DOT_PRODUCT( eri_3center_eigen_lr(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * fact_full_i * 1.0_dp  &
                          * (-gamma_cohsex)

         !
         ! LR-COH
         !
         se%sigma(0,pstate,ipspin) = se%sigma(0,pstate,ipspin) &
                    +  DOT_PRODUCT( eri_3center_eigen_lr(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * alpha_cohsex * 0.5_dp

       enddo


     case default 
       call die('BUG')
     end select

   enddo !istate
 enddo !ipspin

 deallocate(wp0)

 ! Sum up the contribution from different procs
 call xsum_ortho(se%sigma)


 write(stdout,'(a)') ' Sigma_c is calculated'



 select case(calc_type%selfenergy_approx)
 case(TUNED_COHSEX) !==========================================================

  call die('coding is corrupted')
! this is corrupted ! FIXME
!
!   sigc(nsemin:nsemax,:) = sigc(nsemin:nsemax,:) + se%sigma(0,nsemin:nsemax,:)
!
!   write(stdout,'(/,a)') ' COHSEX Eigenvalues (eV)'
!   if(nspin==1) then
!     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
!   else
!     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
!   endif
!   do pstate=nsemin,nsemax
!     energy_qp_new(pstate,:) = energy_qp(pstate,:) + sigc(pstate,:) + exchange_m_vxc_diag(pstate,:)
!
!     write(stdout,'(i4,1x,20(1x,f12.6))') pstate,energy_qp(pstate,:)*Ha_eV,               &
!                                                 exchange_m_vxc_diag(pstate,:)*Ha_eV,     &
!                                                 sigc(pstate,:)*Ha_eV, &
!                                           1.0_dp ,energy_qp_new(pstate,:)*Ha_eV
!   enddo
!
!   call write_energy_qp(nstate,energy_qp_new)



 end select



 call destroy_eri_3center_eigen_lr()


 call stop_clock(timing_self)


end subroutine cohsex_selfenergy_lr


!=========================================================================
