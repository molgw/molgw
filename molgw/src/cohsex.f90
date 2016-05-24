!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the GW self-energy
! within different flavors: G0W0, GnW0, GnWn, COHSEX
!
!=========================================================================
subroutine cohsex_selfenergy(nstate,gwmethod,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy,energy_gw)
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam
 use m_warning,only: issue_warning,msg
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                 :: nstate,gwmethod
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)               :: energy_gw
!=====
 logical               :: file_exists=.FALSE.
 integer               :: homo
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 real(dp),allocatable  :: selfenergy_omega(:,:,:,:)
 real(dp),allocatable  :: sigma_xc_m_vxc_diag(:)
 integer               :: astate,bstate
 integer               :: istate,ispin
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin)
 integer               :: reading_status
 integer               :: selfenergyfile
 integer               :: nsemin,nsemax
 integer               :: ibf_auxil,jbf_auxil
 integer               :: ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable  :: wp0(:,:),wp0_i(:),w0_local(:)
!=====

 call start_clock(timing_self)

 if( .NOT. has_auxil_basis )    stop'Not coded'
 if( .NOT. ALLOCATED(wpol%w0) ) stop'static W should be available here'


 write(stdout,*)
 select case(gwmethod)
 case(COHSEX)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 end select

 !
 ! Set the range of states on which to evaluate the self-energy
 nsemin = MAX(ncore_G+1   ,selfenergy_state_min,1)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax


 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)


 nomegai = 0
 allocate(omegai(-nomegai:nomegai))
 omegai(0) = 0.0_dp

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

 do ispin=1,nspin
   !
   ! Apply the frozen core and frozen virtual approximation to G
   do istate=ncore_G+1,nvirtual_G-1

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
       wp0_i(nsemin:nsemax) = MATMUL( w0_local(:) , eri_3center_eigen(:,nsemin:nsemax,istate,ispin) )
       call xsum(wp0_i)

       if( iproc_ibf_auxil(ibf_auxil_global) == rank ) then
         wp0(ibf_auxil_l(ibf_auxil_global),:) = wp0_i(:)
       endif

     enddo
     deallocate(w0_local,wp0_i)





     ! The application of residu theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ispin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ispin)) / spin_fact



     select case(gwmethod)

     case(COHSEX_DEVEL)

       do astate=nsemin,nsemax
         !
         ! SEX
         !
         selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                    -  DOT_PRODUCT( eri_3center_eigen(:,astate,istate,ispin) , wp0(:,astate) ) &
                          * fact_full_i * 1.0_dp  &
                          * beta_cohsex

         !
         ! COH
         !
         selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                    +  DOT_PRODUCT( eri_3center_eigen(:,astate,istate,ispin) , wp0(:,astate) ) &
                          * alpha_cohsex * 0.5_dp

       enddo

!FBFB
!     case(QSCOHSEX) 
!       do bstate=nsemin,nsemax
!         do astate=nsemin,nsemax
!           !
!           ! SEX
!           !
!           selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
!                      + bra(ipole,astate) * bra(ipole,bstate)                                & 
!                            * fact_full_i / wpol%pole(ipole) * 2.0_dp                        &
!                            * beta_cohsex
!
!           !
!           ! COH
!           !
!           selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
!                      - bra(ipole,astate) * bra(ipole,bstate) & 
!                            / wpol%pole(ipole)                &
!                            * alpha_cohsex
!         enddo
!       enddo


     case default 
       call die('BUG')
     end select

   enddo !istate
 enddo !ispin

 deallocate(wp0)

 ! Sum up the contribution from different procs
 if( ALLOCATED(selfenergy_omega) ) then
   call xsum(selfenergy_omega)
 endif


 write(stdout,'(a)') ' Sigma_c is calculated'

 !
 ! Kotani's Hermitianization trick
 !
 if(gwmethod==QS) then
   do ispin=1,nspin
     selfenergy_omega(0,:,:,ispin) = 0.5_dp * ( selfenergy_omega(0,:,:,ispin) + TRANSPOSE(selfenergy_omega(0,:,:,ispin)) )
   enddo
 endif

 select case(gwmethod)
 case(QS)
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(QSCOHSEX)
   write(stdout,*) 
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(COHSEX_DEVEL) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall
   
   write(stdout,'(/,a)') ' COHSEX Eigenvalues (eV)'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do astate=nsemin,nsemax
     energy_qp_new(astate,:) = energy_qp(astate,:) + selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,     &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV, &
                                           1.0_dp ,energy_qp_new(astate,:)*Ha_eV
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

 if(ALLOCATED(omegai)) deallocate(omegai)
 if(ALLOCATED(selfenergy_omega)) deallocate(selfenergy_omega)

 call stop_clock(timing_self)


end subroutine cohsex_selfenergy


!=========================================================================
subroutine cohsex_selfenergy_lr(nstate,gwmethod,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy,energy_gw)
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
 implicit none

 integer,intent(in)                 :: nstate,gwmethod
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)               :: energy_gw
!=====
 logical               :: file_exists=.FALSE.
 integer               :: homo
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 real(dp),allocatable  :: selfenergy_omega(:,:,:,:)
 real(dp),allocatable  :: sigma_xc_m_vxc_diag(:)
 integer               :: astate,bstate
 integer               :: istate,ispin
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin)
 integer               :: reading_status
 integer               :: selfenergyfile
 integer               :: nsemin,nsemax
 integer               :: ibf_auxil,jbf_auxil
 integer               :: ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable  :: wp0(:,:),wp0_i(:),w0_local(:)
 real(dp),allocatable  :: wp0_tmp(:,:),wp0_rotation(:,:)
 integer               :: nbf_auxil
!=====

 call start_clock(timing_self)

 if( .NOT. has_auxil_basis )    stop'Not coded'
 if( .NOT. ALLOCATED(wpol%w0) ) stop'static W should be available here'

#ifndef TODAY
 stop '-DTODAY is required'
#endif

 write(stdout,*)
 select case(gwmethod)
 case(COHSEX)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 end select

 ! Rotation of W0

 nbf_auxil    = SIZE( eri_2center_rotation   (:,:) , DIM=1)

 allocate( wp0_tmp(nbf_auxil,nbf_auxil) )
 allocate( wp0_rotation(nauxil_2center_lr,nauxil_2center_lr) )
 wp0_tmp(:,:)      = MATMUL( eri_2center_rotation(:,:) , MATMUL( wpol%w0(:,:) , TRANSPOSE( eri_2center_rotation(:,:) ) ) )
 wp0_rotation(:,:) = MATMUL( TRANSPOSE(eri_2center_rotation_lr(:,:)) , MATMUL( wp0_tmp(:,:) , eri_2center_rotation_lr(:,:) ) )
 deallocate( wp0_tmp )


 !
 ! Set the range of states on which to evaluate the self-energy
 nsemin = MAX(ncore_G+1   ,selfenergy_state_min,1)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax


 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)


 nomegai = 0
 allocate(omegai(-nomegai:nomegai))
 omegai(0) = 0.0_dp

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

 do ispin=1,nspin
   !
   ! Apply the frozen core and frozen virtual approximation to G
   do istate=ncore_G+1,nvirtual_G-1

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
       wp0_i(nsemin:nsemax) = MATMUL( w0_local(:) , eri_3center_eigen_lr(:,nsemin:nsemax,istate,ispin) )
       call xsum(wp0_i)

       if( iproc_ibf_auxil(ibf_auxil_global) == rank ) then
         wp0(ibf_auxil_l(ibf_auxil_global),:) = wp0_i(:)
       endif

     enddo
     deallocate(w0_local,wp0_i)





     ! The application of residu theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ispin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ispin)) / spin_fact



     select case(gwmethod)

     case(COHSEX_DEVEL)

       do astate=nsemin,nsemax
         !
         ! SEX
         !
         selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                    -  DOT_PRODUCT( eri_3center_eigen_lr(:,astate,istate,ispin) , wp0(:,astate) ) &
                          * fact_full_i * 1.0_dp  &
                          * beta_cohsex

         !
         ! COH
         !
         selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                    +  DOT_PRODUCT( eri_3center_eigen_lr(:,astate,istate,ispin) , wp0(:,astate) ) &
                          * alpha_cohsex * 0.5_dp

       enddo

!FBFB
!     case(QSCOHSEX) 
!       do bstate=nsemin,nsemax
!         do astate=nsemin,nsemax
!           !
!           ! SEX
!           !
!           selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
!                      + bra(ipole,astate) * bra(ipole,bstate)                                & 
!                            * fact_full_i / wpol%pole(ipole) * 2.0_dp                        &
!                            * beta_cohsex
!
!           !
!           ! COH
!           !
!           selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
!                      - bra(ipole,astate) * bra(ipole,bstate) & 
!                            / wpol%pole(ipole)                &
!                            * alpha_cohsex
!         enddo
!       enddo


     case default 
       call die('BUG')
     end select

   enddo !istate
 enddo !ispin

 deallocate(wp0)

 ! Sum up the contribution from different procs
 if( ALLOCATED(selfenergy_omega) ) then
   call xsum(selfenergy_omega)
 endif


 write(stdout,'(a)') ' Sigma_c is calculated'

 !
 ! Kotani's Hermitianization trick
 !
 if(gwmethod==QS) then
   do ispin=1,nspin
     selfenergy_omega(0,:,:,ispin) = 0.5_dp * ( selfenergy_omega(0,:,:,ispin) + TRANSPOSE(selfenergy_omega(0,:,:,ispin)) )
   enddo
 endif

 select case(gwmethod)
 case(QS)
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(QSCOHSEX)
   write(stdout,*) 
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(COHSEX_DEVEL) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall
   
   write(stdout,'(/,a)') ' COHSEX Eigenvalues (eV)'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do astate=nsemin,nsemax
     energy_qp_new(astate,:) = energy_qp(astate,:) + selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,     &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV, &
                                           1.0_dp ,energy_qp_new(astate,:)*Ha_eV
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

 if(ALLOCATED(omegai)) deallocate(omegai)
 if(ALLOCATED(selfenergy_omega)) deallocate(selfenergy_omega)

 call stop_clock(timing_self)


end subroutine cohsex_selfenergy_lr


!=========================================================================
