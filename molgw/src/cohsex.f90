!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the COHSEX self-energy
!
!=========================================================================
subroutine cohsex_selfenergy(nstate,gwmethod,basis,occupation,energy,exchange_m_vxc_diag, &
                             c_matrix,s_matrix,wpol,selfenergy,sigc,energy_gw)
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

 integer,intent(in)                 :: nstate,gwmethod
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
 real(dp),intent(inout)             :: sigc(nstate,nspin)
 real(dp),intent(out)               :: energy_gw
!=====
 real(dp),allocatable  :: selfenergy_omega(:,:,:,:)
 integer               :: pstate
 integer               :: istate,ipspin
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin)
 integer               :: jbf_auxil
 integer               :: ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable  :: wp0(:,:),wp0_i(:),w0_local(:)
 real(dp)              :: sigx
!=====

 call start_clock(timing_self)

 if( .NOT. has_auxil_basis )    call die('cohsex: no RI is not coded')
 if( .NOT. ALLOCATED(wpol%w0) ) call die('cohsex: static W should be available here')


 select case(gwmethod)
 case(COHSEX_DEVEL)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case(TUNED_COHSEX)
   write(stdout,*) 'Perform a tuned COHSEX calculation'
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

     ! Store the result in sigc
     sigc(pstate,ipspin) = sigc(pstate,ipspin) + sigx * epsilon_cohsex

   enddo
 enddo
 write(stdout,*) '=============='


 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 call issue_warning('small complex number is '//msg)


 !
 ! Which calculation type needs to update energy_qp
 !
!FBFB
! select case(gwmethod)
! case(GnW0,GnWn,GSIGMA3)
!   call read_energy_qp(nstate,energy_qp,reading_status)
!   if(reading_status/=0) then
!     call issue_warning('File energy_qp not found: assuming 1st iteration')
!     energy_qp(:,:) = energy(:,:)
!   endif
! case default
   energy_qp(:,:) = energy(:,:)

! end select

 !
 ! Which calculation type needs the diagonal only?
 ! Which calculation type needs a complex sigma?
 !
!FBFB
! select case(gwmethod)
! case(COHSEX,GV,GSIGMA,G0W0,GnW0,GnWn)   ! diagonal real
   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,1,nspin))

! case(QS,QSCOHSEX,GSIGMA3) ! matrix real
!   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nsemin:nsemax,nspin))
!
! case default
!   call die('GW case does not exist. Should not happen')
! end select
 if( ALLOCATED(selfenergy_omega) )  selfenergy_omega(:,:,:,:)  = 0.0_dp


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





     ! The application of residu theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ipspin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ipspin)) / spin_fact



     select case(gwmethod)

     case(COHSEX_DEVEL)

       do pstate=nsemin,nsemax
         !
         ! SEX
         !
         selfenergy_omega(0,pstate,1,ipspin) = selfenergy_omega(0,pstate,1,ipspin) &
                    -  DOT_PRODUCT( eri_3center_eigen(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * fact_full_i * 1.0_dp  &
                          * beta_cohsex

         !
         ! COH
         !
         selfenergy_omega(0,pstate,1,ipspin) = selfenergy_omega(0,pstate,1,ipspin) &
                    +  DOT_PRODUCT( eri_3center_eigen(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * alpha_cohsex * 0.5_dp

       enddo

     case(TUNED_COHSEX)

       do pstate=nsemin,nsemax
         !
         ! SEX
         !
         selfenergy_omega(0,pstate,1,ipspin) = selfenergy_omega(0,pstate,1,ipspin) &
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
 if( ALLOCATED(selfenergy_omega) ) then
   call xsum_ortho(selfenergy_omega)
 endif


 write(stdout,'(a)') ' Sigma_c is calculated'


 select case(gwmethod)

 case(QSCOHSEX)     !==========================================================
   write(stdout,*) 
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ipspin=1,nspin
     selfenergy(:,:,ipspin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ipspin) ) , &
                                MATMUL( selfenergy_omega(0,:,:,ipspin), &
                                  MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ipspin)), s_matrix(:,:) ) ) )
   enddo


 case(COHSEX_DEVEL) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(pstate=nsemin:nsemax)
     selfenergy(pstate,pstate,:) = selfenergy_omega(0,pstate,1,:)
   end forall
   
   call find_qp_energy_linearization(nomegai,omegai,nsemin,nsemax,selfenergy_omega(:,:,1,:),nstate,exchange_m_vxc_diag,energy_qp,energy_qp_new)
   call output_qp_energy('COHSEX',nstate,nsemin,nsemax,energy_qp,exchange_m_vxc_diag,1,selfenergy_omega(0,:,1,:),energy_qp_new)

   call write_energy_qp(nstate,energy_qp_new)


 case(TUNED_COHSEX) !==========================================================

   sigc(nsemin:nsemax,:) = sigc(nsemin:nsemax,:) + selfenergy_omega(0,nsemin:nsemax,1,:)


 end select


 !
 ! Output the new HOMO and LUMO energies
 !
 if( gwmethod == COHSEX_DEVEL ) then
   call output_new_homolumo('COHSEX',nstate,occupation,energy_qp_new,nsemin,nsemax)
 endif

 if(ALLOCATED(selfenergy_omega)) deallocate(selfenergy_omega)

 call destroy_eri_3center_eigen()

 call stop_clock(timing_self)


end subroutine cohsex_selfenergy


!=========================================================================
subroutine cohsex_selfenergy_lr(nstate,gwmethod,basis,occupation,energy,exchange_m_vxc_diag, &
                                c_matrix,s_matrix,wpol,selfenergy,sigc,energy_gw)
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
 use m_eri_lr_calculate
 use m_selfenergy_tools
 implicit none

 integer,intent(in)                 :: nstate,gwmethod
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
 real(dp),intent(inout)             :: sigc(nstate,nspin)
 real(dp),intent(out)               :: energy_gw
!=====
 integer               :: homo
 real(dp),allocatable  :: selfenergy_omega(:,:,:,:)
 integer               :: pstate
 integer               :: istate,ipspin
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin)
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
 select case(gwmethod)
 case(COHSEX)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case(TUNED_COHSEX)
   write(stdout,*) 'Perform a tuned COHSEX calculation'
 end select


 ! Rotation of W0

 nbf_auxil    = SIZE( eri_2center_rotation   (:,:) , DIM=1)

 allocate( wp0_tmp(nbf_auxil,nbf_auxil) )
 allocate( wp0_rotation(nauxil_2center_lr,nauxil_2center_lr) )
 wp0_tmp(:,:)      = MATMUL( eri_2center_rotation(:,:) , MATMUL( wpol%w0(:,:) , TRANSPOSE( eri_2center_rotation(:,:) ) ) )
 wp0_rotation(:,:) = MATMUL( TRANSPOSE(eri_2center_rotation_lr(:,:)) , MATMUL( wp0_tmp(:,:) , eri_2center_rotation_lr(:,:) ) )
 deallocate( wp0_tmp )



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

     ! Store the result in sigc
     sigc(pstate,ipspin) = sigc(pstate,ipspin) - sigx * epsilon_cohsex

   enddo
 enddo
 write(stdout,*) '=============='



 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 call issue_warning('small complex number is '//msg)


 !
 ! Which calculation type needs to update energy_qp
 !
!FBFB
! select case(gwmethod)
! case(GnW0,GnWn,GSIGMA3)
!   call read_energy_qp(nstate,energy_qp,reading_status)
!   if(reading_status/=0) then
!     call issue_warning('File energy_qp not found: assuming 1st iteration')
!     energy_qp(:,:) = energy(:,:)
!   endif
! case default
   energy_qp(:,:) = energy(:,:)

! end select

 !
 ! Which calculation type needs the diagonal only?
 ! Which calculation type needs a complex sigma?
 !
!FBFB
! select case(gwmethod)
! case(COHSEX,GV,GSIGMA,G0W0,GnW0,GnWn)   ! diagonal real
   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,1,nspin))

! case(QS,QSCOHSEX,GSIGMA3) ! matrix real
!   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nsemin:nsemax,nspin))
!
! case default
!   call die('GW case does not exist. Should not happen')
! end select
 if( ALLOCATED(selfenergy_omega) )  selfenergy_omega(:,:,:,:)  = 0.0_dp


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





     ! The application of residu theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ipspin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ipspin)) / spin_fact



     select case(gwmethod)

     case(COHSEX_DEVEL)

       do pstate=nsemin,nsemax
         !
         ! SEX
         !
         selfenergy_omega(0,pstate,1,ipspin) = selfenergy_omega(0,pstate,1,ipspin) &
                    -  DOT_PRODUCT( eri_3center_eigen_lr(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * fact_full_i * 1.0_dp  &
                          * beta_cohsex

         !
         ! COH
         !
         selfenergy_omega(0,pstate,1,ipspin) = selfenergy_omega(0,pstate,1,ipspin) &
                    +  DOT_PRODUCT( eri_3center_eigen_lr(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * alpha_cohsex * 0.5_dp

       enddo

     case(TUNED_COHSEX)

       do pstate=nsemin,nsemax
         !
         ! LR-SEX 
         !
         selfenergy_omega(0,pstate,1,ipspin) = selfenergy_omega(0,pstate,1,ipspin) &
                    -  DOT_PRODUCT( eri_3center_eigen_lr(:,pstate,istate,ipspin) , wp0(:,pstate) ) &
                          * fact_full_i * 1.0_dp  &
                          * (-gamma_cohsex)

         !
         ! LR-COH
         !
         selfenergy_omega(0,pstate,1,ipspin) = selfenergy_omega(0,pstate,1,ipspin) &
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
 if( ALLOCATED(selfenergy_omega) ) then
   call xsum_ortho(selfenergy_omega)
 endif


 write(stdout,'(a)') ' Sigma_c is calculated'



 select case(gwmethod)
 case(QSCOHSEX)     !==========================================================
   write(stdout,*) 
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ipspin=1,nspin
     selfenergy(:,:,ipspin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ipspin) ) , &
                                MATMUL( selfenergy_omega(0,:,:,ipspin), &
                                  MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ipspin)), s_matrix(:,:) ) ) )
   enddo


 case(COHSEX_DEVEL) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(pstate=nsemin:nsemax)
     selfenergy(pstate,pstate,:) = selfenergy_omega(0,pstate,1,:)
   end forall
   
   write(stdout,'(/,a)') ' COHSEX Eigenvalues (eV)'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do pstate=nsemin,nsemax
     energy_qp_new(pstate,:) = energy_qp(pstate,:) + selfenergy_omega(0,pstate,1,:) + exchange_m_vxc_diag(pstate,:)

     write(stdout,'(i4,1x,20(1x,f12.6))') pstate,energy_qp(pstate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(pstate,:)*Ha_eV,     &
                                                 selfenergy_omega(0,pstate,1,:)*Ha_eV, &
                                           1.0_dp ,energy_qp_new(pstate,:)*Ha_eV
   enddo

   call write_energy_qp(nstate,energy_qp_new)

 case(TUNED_COHSEX) !==========================================================


   sigc(nsemin:nsemax,:) = sigc(nsemin:nsemax,:) + selfenergy_omega(0,nsemin:nsemax,1,:)

   write(stdout,'(/,a)') ' COHSEX Eigenvalues (eV)'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do pstate=nsemin,nsemax
     energy_qp_new(pstate,:) = energy_qp(pstate,:) + sigc(pstate,:) + exchange_m_vxc_diag(pstate,:)

     write(stdout,'(i4,1x,20(1x,f12.6))') pstate,energy_qp(pstate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(pstate,:)*Ha_eV,     &
                                                 sigc(pstate,:)*Ha_eV, &
                                           1.0_dp ,energy_qp_new(pstate,:)*Ha_eV
   enddo

   call write_energy_qp(nstate,energy_qp_new)



 end select


 !
 ! Output the new HOMO and LUMO energies
 !
 do istate=1,nstate
   if( ANY(occupation(istate,:) > completely_empty) ) homo = istate
 enddo
 write(stdout,*)
 if( homo >= nsemin .AND. homo <= nsemax ) then
   write(stdout,'(a,2(2x,f12.6))') ' GW HOMO (eV):',energy_qp_new(homo,:)*Ha_eV
 endif
 if( homo+1 >= nsemin .AND. homo+1 <= nsemax ) then
   write(stdout,'(a,2(2x,f12.6))') ' GW LUMO (eV):',energy_qp_new(homo+1,:)*Ha_eV
 endif

 if(ALLOCATED(selfenergy_omega)) deallocate(selfenergy_omega)

 call destroy_eri_3center_eigen_lr()


 call stop_clock(timing_self)


end subroutine cohsex_selfenergy_lr


!=========================================================================
