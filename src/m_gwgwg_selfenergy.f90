!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains the calculation of the vertex function
! within different flavors: SOX, SOSEX, GWGWG, static GW0GW0G
!
!=========================================================================
#include "molgw.h"
module m_gwgwg_selfenergy
  use m_definitions
  use m_mpi
  use m_timing
  use m_inputparam
  use m_warning
  use m_basis_set
  use m_spectral_function
  use m_eri_ao_mo
  use m_selfenergy_tools
  use m_tddft_fxc
  use m_linear_response


contains


!=========================================================================
subroutine sosex_selfenergy(basis,occupation,energy,c_matrix,wpol,se)
  implicit none

  type(basis_set)                    :: basis
  real(dp),intent(in)                :: occupation(:,:),energy(:,:)
  real(dp),intent(in)                :: c_matrix(:,:,:)
  type(spectral_function),intent(inout) :: wpol
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer                 :: nstate
  integer                 :: iomega
  complex(dp),allocatable :: sigma_gw(:,:,:)
  complex(dp),allocatable :: sigma_sosex(:,:,:)
  complex(dp),allocatable :: sigma_sox(:,:,:)
  integer                 :: astate,bstate,cstate
  integer                 :: istate,jstate,kstate,ispin,spole
  integer                 :: pstate
  real(dp),allocatable    :: bra_s(:,:)
  real(dp)                :: vcoul,vcoul1,vcoul2
  real(dp)                :: pole_s
  real(dp)                :: fxc
  !=====

  call start_clock(timing_gwgamma_self)
  nstate = SIZE(energy,DIM=1)

  ! Turn dynamic into static
  !call issue_warning('FBFB HACK in SOSEX')
  !wpol%pole(:) = wpol%pole(:) * 1000.0d0
  !wpol%residue_left(:,:) = wpol%residue_left(:,:) * SQRT( 1000.0d0 )

  write(stdout,*)
  select case(calc_type%selfenergy_approx)
  case(GWSOX)
    write(stdout,*) 'Perform a one-shot GW+SOX calculation'
  case(GWSOSEX)
    if( calc_type%selfenergy_technique == EVSC ) then
      write(stdout,*) 'Perform an eigenvalue-self-consistent GW+SOSEX calculation'
    else
      write(stdout,*) 'Perform a one-shot GW+SOSEX calculation'
    endif
  case(GWGWG,GWGWG_numerical)
    write(stdout,*) 'Perform a one-shot GW+SOSEX calculation to prepare full GWGWG'
  case default
    call die('sosex_selfenergy: calculation type unknown')
  end select


  if( gwgamma_tddft_ ) then
    write(stdout,*) 'Include a TDDFT kernel contribution to the vertex'
    write(stdout,'(1x,a,f12.4)') 'Exact-exchange amount: ',alpha_hybrid
    call prepare_tddft(.FALSE.,nstate,basis,c_matrix,occupation)
  endif

  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
  else
    call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
  endif


  call clean_allocate('Temporary array',bra_s,ncore_G+1,nvirtual_G-1,ncore_G+1,MAX(nhomo_G,nsemax))


  !
  !
  allocate(sigma_sosex(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_sox(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gw(-se%nomega:se%nomega,nsemin:nsemax,nspin))

  sigma_sosex(:,:,:)  = 0.0_dp
  sigma_sox(:,:,:)  = 0.0_dp


  write(stdout,*) 'Calculate SOX'

  do ispin=1,nspin

    !==========================
    do bstate=ncore_G+1,nvirtual_G-1
      if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty) cycle
      if( MODULO( bstate-(ncore_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      do istate=ncore_G+1,nvirtual_G-1
        if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
        do kstate=ncore_G+1,nvirtual_G-1
          if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle

          do pstate=nsemin,nsemax

            vcoul1 = eri_eigen(pstate,istate,ispin,bstate,kstate,ispin)
            vcoul2 = eri_eigen(istate,bstate,ispin,kstate,pstate,ispin)
            if( gwgamma_tddft_ ) then
              fxc = eval_fxc_rks_singlet(istate,bstate,ispin,kstate,pstate,ispin)
              call grid%sum(fxc)
              vcoul2 = alpha_hybrid * vcoul2 - fxc

              !             if( ABS( eri_eigen(istate,bstate,ispin,kstate,pstate,ispin) -vcoul2)> 0.10 ) then
              !               write(*,'(4(i4,1x),4(1x,f12.6))') istate,bstate,kstate,pstate, &
              !                  eri_eigen(istate,bstate,ispin,kstate,pstate,ispin), &
              !                  vcoul2
              !               write(*,*) 'Hack'
              !               vcoul2 = eri_eigen(istate,bstate,ispin,kstate,pstate,ispin)
              !             endif

            endif
            !
            ! calculate only the diagonal !
            do iomega=-se%nomega,se%nomega
              sigma_sox(iomega,pstate,ispin) = sigma_sox(iomega,pstate,ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                        - energy(istate,ispin) - energy(kstate,ispin) + energy(bstate,ispin) - ieta )

            enddo
          enddo

        enddo
      enddo
    enddo

    !==========================
    do cstate=ncore_G+1,nvirtual_G-1
      if( (spin_fact - occupation(cstate,ispin)) / spin_fact < completely_empty) cycle
      if( MODULO( cstate-(ncore_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      do jstate=ncore_G+1,nvirtual_G-1
        if( occupation(jstate,ispin) / spin_fact < completely_empty ) cycle
        do astate=ncore_G+1,nvirtual_G-1
          if( (spin_fact - occupation(astate,ispin)) / spin_fact < completely_empty) cycle

          do pstate=nsemin,nsemax

            vcoul1 = eri_eigen(pstate,astate,ispin,jstate,cstate,ispin)
            vcoul2 = eri_eigen(astate,jstate,ispin,cstate,pstate,ispin)
            if( gwgamma_tddft_ ) then
              fxc = eval_fxc_rks_singlet(astate,jstate,ispin,cstate,pstate,ispin)
              call grid%sum(fxc)
              vcoul2 = alpha_hybrid * vcoul2 - fxc

              !             if( ABS( eri_eigen(astate,jstate,ispin,cstate,pstate,ispin) -vcoul2 )> 0.10 ) then
              !               write(*,'(4(i4,1x),4(1x,f12.6))') astate,jstate,cstate,pstate, &
              !                  eri_eigen(astate,jstate,ispin,cstate,pstate,ispin), &
              !                  vcoul2
              !!               write(*,*) 'Hack'
              !!               vcoul2 =  eri_eigen(astate,jstate,ispin,cstate,pstate,ispin)
              !             endif

            endif
            !
            ! calculate only the diagonal !
            do iomega=-se%nomega,se%nomega
              sigma_sox(iomega,pstate,ispin) = sigma_sox(iomega,pstate,ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                        - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) + ieta )
            enddo
          enddo

        enddo
      enddo
    enddo


  enddo

  call ortho%sum(sigma_sox)


  if( calc_type%selfenergy_approx == GWSOSEX .OR. calc_type%selfenergy_approx == GWGWG &
      .OR. calc_type%selfenergy_approx == GWGWG_NUMERICAL ) then

    write(stdout,*) 'Calculate dynamical SOSEX'


    do ispin=1,nspin

      do spole=1,wpol%npole_reso

        if( MODULO( spole - 1 , ortho%nproc ) /= ortho%rank ) cycle
        write(stdout,*) 'SOSEX W poles:',spole,' / ',wpol%npole_reso

        pole_s = wpol%pole(spole)

        if(has_auxil_basis) then
          do pstate=ncore_G+1,MAX(nhomo_G,nsemax)
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
            bra_s(:,pstate)     = MATMUL( wpol%residue_left(:,spole) , eri_3center_eigen(:,:,pstate,ispin) )
          enddo
          call auxil%sum(bra_s)
        else
          ! Here just grab the precalculated value
          forall(istate=ncore_G+1:nvirtual_G-1, pstate=ncore_G+1:MAX(nhomo_G,nsemax))
            bra_s(istate,pstate) = wpol%residue_left(index_prodstate(istate,pstate) &
                                                    + (ispin-1) * index_prodstate(nvirtual_W-1,nvirtual_W-1), &
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
              do pstate=nsemin,nsemax

                vcoul = eri_eigen(istate,kstate,ispin,bstate,pstate,ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(istate,kstate,ispin,bstate,pstate,ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,pstate,ispin) = sigma_sosex(iomega,pstate,ispin) &
                           - bra_s(kstate,pstate) * bra_s(bstate,istate) * vcoul                          &
                              / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + pole_s - ieta )  &
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
              do pstate=nsemin,nsemax

                vcoul = eri_eigen(istate,cstate,ispin,bstate,pstate,ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(istate,cstate,ispin,bstate,pstate,ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,pstate,ispin) = sigma_sosex(iomega,pstate,ispin) &
                           - bra_s(cstate,pstate) * bra_s(bstate,istate) * vcoul                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )    &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                                 - energy(cstate,ispin) + energy(istate,ispin) - energy(bstate,ispin) + ieta )


                  sigma_sosex(iomega,pstate,ispin) = sigma_sosex(iomega,pstate,ispin) &
                           + bra_s(cstate,pstate) * bra_s(bstate,istate) * vcoul                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(bstate,ispin) &
                                 - energy(cstate,ispin) + energy(istate,ispin) + ieta )  &
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
              do pstate=nsemin,nsemax

                vcoul = eri_eigen(astate,kstate,ispin,jstate,pstate,ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(astate,kstate,ispin,jstate,pstate,ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,pstate,ispin) = sigma_sosex(iomega,pstate,ispin) &
                           - bra_s(kstate,pstate) * bra_s(astate,jstate) * vcoul                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                             / ( energy(jstate,ispin) - energy(astate,ispin) - pole_s + ieta )

                  sigma_sosex(iomega,pstate,ispin) = sigma_sosex(iomega,pstate,ispin) &
                           + bra_s(kstate,pstate) * bra_s(astate,jstate) * vcoul                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                 + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + pole_s - ieta )


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
              do pstate=nsemin,nsemax

                vcoul = eri_eigen(astate,cstate,ispin,jstate,pstate,ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(astate,cstate,ispin,jstate,pstate,ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,pstate,ispin) = sigma_sosex(iomega,pstate,ispin) &
                           + bra_s(cstate,pstate) * bra_s(astate,jstate) * vcoul                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )  &
                             / ( pole_s + energy(astate,ispin) - energy(jstate,ispin) - ieta )

                enddo
              enddo

            enddo
          enddo
        enddo



      enddo !spole
    enddo !ispin

    call ortho%sum(sigma_sosex)

  endif


  write(stdout,'(a)') ' Sigma_c(omega) is calculated'


  !
  ! The input sigma contains the GW selfenergy
  sigma_gw(:,:,:) = se%sigma(:,:,:)


  forall(pstate=nsemin:nsemax)
    se%sigma(:,pstate,:) = sigma_gw(:,pstate,:) + sigma_sox(:,pstate,:) + factor_sosex * sigma_sosex(:,pstate,:)
  end forall


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_sox'    ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sox)
  !   call write_selfenergy_omega('selfenergy_sosex'  ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sosex)
  ! endif


  write(stdout,'(/,a)') ' GW+SOSEX self-energy contributions at E0 (eV)'
  write(stdout,'(a)') &
     '   #          E0       SigC_GW     SigC_SOX   SigC_G(W(w)-v)GvG SigC_TOT'

  do pstate=nsemin,nsemax
    write(stdout,'(i4,1x,20(1x,f12.6))') pstate,se%energy0(pstate,:)*Ha_eV,          &
                                         sigma_gw(0,pstate,:)%re*Ha_eV,   &
                                         sigma_sox(0,pstate,:)%re*Ha_eV,  &
                                         sigma_sosex(0,pstate,:)%re*Ha_eV,&
                                         se%sigma(0,pstate,:)%re*Ha_eV
  enddo




  call clean_deallocate('Temporary array',bra_s)

  if(has_auxil_basis) then
    call destroy_eri_3center_eigen()
  else
    call destroy_eri_4center_eigen_uks()
  endif

  if( gwgamma_tddft_ ) then
    call destroy_tddft()
  endif

  call stop_clock(timing_gwgamma_self)


end subroutine sosex_selfenergy


!=========================================================================
subroutine gwgw0g_selfenergy(nstate,basis,occupation,energy,c_matrix,wpol,se)
  implicit none

  integer,intent(in)                 :: nstate
  type(basis_set)                    :: basis
  real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
  type(spectral_function),intent(in) :: wpol
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer                 :: iomega,ibf_auxil
  complex(dp),allocatable :: sigma_gwgw0g(:,:,:)
  complex(dp),allocatable :: sigma_gw0gw0g(:,:,:)
  complex(dp),allocatable :: sigma_gw0gw0g_occ(:,:,:)
  complex(dp),allocatable :: sigma_gw0gw0g_vir(:,:,:)
  complex(dp),allocatable :: sigma_gvgw0g(:,:,:)
  integer                 :: astate,bstate,cstate
  integer                 :: istate,jstate,kstate,ispin,spole
  integer                 :: pstate
  real(dp),allocatable    :: bra_s(:,:)
  real(dp)                :: v_1,w0_1,w0_2
  real(dp)                :: pole_s
  real(dp)                :: fxc
  real(dp),allocatable    :: chi_static(:,:)
  real(dp),allocatable    :: ip(:),bk(:),ib(:),kp(:)
  real(dp),allocatable    :: pa(:),jc(:),aj(:),cp(:)
  real(dp),allocatable    :: ik(:),ac(:),bp(:),ak(:),jp(:),ic(:)
  !=====

  call start_clock(timing_gwgamma_self)
  if( .NOT. has_auxil_basis ) call die('gwgw0g_selfenergy: not implemented without an auxiliary basis')

  write(stdout,*)
  select case(calc_type%selfenergy_approx)
  case(GW0GW0G)
    write(stdout,*) 'Perform a one-shot GW+GW0GW0G calculation'
  case(GWGW0G)
      write(stdout,*) 'Perform a one-shot GW+GWGW0G calculation'
  case default
    call die('gwgw0g_selfenergy: calculation type unknown')
  end select


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)

  call clean_allocate('Temporary array',bra_s,ncore_G+1,nvirtual_G-1,ncore_G+1,MAX(nhomo_G,nsemax))


  call clean_allocate('chi static',chi_static,nauxil_global,nauxil_global)
  call wpol%evaluate((0.0_dp,0.0_dp),chi_static)

  ! Turn dynamic into static for debug purposes
  !call issue_warning('hack to recover GW+SOX for debug purposes')
  !chi_static(:,:) = 0.0d0   ! to recover GW+SOX
  do ibf_auxil=1,nauxil_global
    chi_static(ibf_auxil,ibf_auxil) = chi_static(ibf_auxil,ibf_auxil) + 1.0_dp
  enddo
  allocate(ib(nauxil_global))
  allocate(kp(nauxil_global))
  allocate(ip(nauxil_global))
  allocate(bk(nauxil_global))
  allocate(aj(nauxil_global))
  allocate(cp(nauxil_global))
  allocate(pa(nauxil_global))
  allocate(jc(nauxil_global))
  allocate(ac(nauxil_global))
  allocate(ik(nauxil_global))
  allocate(bp(nauxil_global))
  allocate(ak(nauxil_global))
  allocate(jp(nauxil_global))
  allocate(ic(nauxil_global))

  !
  !
  allocate(sigma_gwgw0g(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gvgw0g(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gw0gw0g(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gw0gw0g_occ(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gw0gw0g_vir(-se%nomega:se%nomega,nsemin:nsemax,nspin))

  sigma_gwgw0g(:,:,:)  = (0.0_dp, 0.0_dp)
  sigma_gvgw0g(:,:,:)  = (0.0_dp, 0.0_dp)
  sigma_gw0gw0g(:,:,:) = (0.0_dp, 0.0_dp)
  sigma_gw0gw0g_occ(:,:,:) = (0.0_dp, 0.0_dp)
  sigma_gw0gw0g_vir(:,:,:) = (0.0_dp, 0.0_dp)


  write(stdout,*) 'Calculate two static terms analog to SOX'

  do ispin=1,nspin

    !==========================
    do bstate=nhomo_G+1,nvirtual_G-1
      if( MODULO( bstate-(ncore_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      do istate=ncore_G+1,nhomo_G
        do kstate=ncore_G+1,nhomo_G

          do pstate=nsemin,nsemax

            !v_1 = eri_eigen(pstate,istate,ispin,bstate,kstate,ispin)
            !v_2 = eri_eigen(istate,bstate,ispin,kstate,pstate,ispin)

            ip(:) = eri_3center_eigen(:,pstate,istate,ispin)
            bk(:) = eri_3center_eigen(:,bstate,kstate,ispin)
            v_1  = DOT_PRODUCT( ip(:) , bk(:) )
            w0_1 = DOT_PRODUCT( ip(:) , MATMUL( chi_static(:,:) , bk(:) ) )

            ib(:) = eri_3center_eigen(:,istate,bstate,ispin)
            kp(:) = eri_3center_eigen(:,kstate,pstate,ispin)
            w0_2 = DOT_PRODUCT( ib(:) , MATMUL( chi_static(:,:) , kp(:) ) )

            !
            ! calculate only the diagonal !
            do iomega=-se%nomega,se%nomega
              sigma_gvgw0g(iomega,pstate,ispin) = sigma_gvgw0g(iomega,pstate,ispin) &
                  - v_1 * w0_2            &
                    / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                        - energy(istate,ispin) - energy(kstate,ispin) + energy(bstate,ispin) - ieta )
              sigma_gw0gw0g_occ(iomega,pstate,ispin) = sigma_gw0gw0g_occ(iomega,pstate,ispin) &
                  - w0_1 * w0_2            &
                    / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                        - energy(istate,ispin) - energy(kstate,ispin) + energy(bstate,ispin) - ieta )
            enddo
          enddo

        enddo
      enddo
    enddo

    !==========================
    do cstate=nhomo_G+1,nvirtual_G-1
      if( MODULO( cstate-(ncore_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      do jstate=ncore_G+1,nhomo_G
        do astate=nhomo_G+1,nvirtual_G-1

          do pstate=nsemin,nsemax

            !v_1 = eri_eigen(pstate,astate,ispin,jstate,cstate,ispin)
            !v_2 = eri_eigen(astate,jstate,ispin,cstate,pstate,ispin)

            pa(:) = eri_3center_eigen(:,pstate,astate,ispin)
            jc(:) = eri_3center_eigen(:,jstate,cstate,ispin)
            v_1  = DOT_PRODUCT( pa(:) , jc(:) )
            w0_1 = DOT_PRODUCT( pa(:) , MATMUL( chi_static(:,:) , jc(:) ) )

            aj(:) = eri_3center_eigen(:,astate,jstate,ispin)
            cp(:) = eri_3center_eigen(:,cstate,pstate,ispin)
            w0_2 = DOT_PRODUCT( aj(:) , MATMUL( chi_static(:,:) , cp(:) ) )

            !
            ! calculate only the diagonal !
            do iomega=-se%nomega,se%nomega
              sigma_gvgw0g(iomega,pstate,ispin) = sigma_gvgw0g(iomega,pstate,ispin) &
                  - v_1 * w0_2            &
                    / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                        - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) + ieta )
              sigma_gw0gw0g_vir(iomega,pstate,ispin) = sigma_gw0gw0g_vir(iomega,pstate,ispin) &
                  - w0_1 * w0_2            &
                    / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                        - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) + ieta )
            enddo
          enddo

        enddo
      enddo
    enddo


  enddo


  call ortho%sum(sigma_gvgw0g)
  call ortho%sum(sigma_gw0gw0g_occ)
  call ortho%sum(sigma_gw0gw0g_vir)
  sigma_gw0gw0g(:,:,:) = sigma_gw0gw0g_occ(:,:,:) + sigma_gw0gw0g_vir(:,:,:)


  if( calc_type%selfenergy_approx == GWGW0G ) then

    write(stdout,*) 'Calculate dynamical term analog to SOSEX'


    do ispin=1,nspin

      do spole=1,wpol%npole_reso

        if( MODULO( spole - 1 , ortho%nproc ) /= ortho%rank ) cycle
        write(stdout,*) 'GWGW0G W poles:',spole,' / ',wpol%npole_reso

        pole_s = wpol%pole(spole)

        if(has_auxil_basis) then
          do pstate=ncore_G+1,MAX(nhomo_G,nsemax)
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
            bra_s(:,pstate)     = MATMUL( wpol%residue_left(:,spole) , eri_3center_eigen(:,:,pstate,ispin) )
          enddo
          call auxil%sum(bra_s)
        else
          ! Here just grab the precalculated value
          forall(istate=ncore_G+1:nvirtual_G-1, pstate=ncore_G+1:MAX(nhomo_G,nsemax))
            bra_s(istate,pstate) = wpol%residue_left(index_prodstate(istate,pstate) &
                                     + (ispin-1) * index_prodstate(nvirtual_W-1,nvirtual_W-1),spole)
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
              do pstate=nsemin,nsemax

                !v_2 = eri_eigen(istate,kstate,ispin,bstate,pstate,ispin)
                ik(:) = eri_3center_eigen(:,istate,kstate,ispin)
                bp(:) = eri_3center_eigen(:,bstate,pstate,ispin)
                w0_2 = DOT_PRODUCT( ik(:) , MATMUL( chi_static(:,:) , bp(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,pstate,ispin) = sigma_gwgw0g(iomega,pstate,ispin) &
                           - bra_s(kstate,pstate) * bra_s(bstate,istate) * w0_2                          &
                              / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + pole_s - ieta )  &
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
              do pstate=nsemin,nsemax

                !v_2 = eri_eigen(istate,cstate,ispin,bstate,pstate,ispin)
                ic(:) = eri_3center_eigen(:,istate,cstate,ispin)
                bp(:) = eri_3center_eigen(:,bstate,pstate,ispin)
                w0_2 = DOT_PRODUCT( ic(:) , MATMUL( chi_static(:,:) , bp(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,pstate,ispin) = sigma_gwgw0g(iomega,pstate,ispin) &
                           - bra_s(cstate,pstate) * bra_s(bstate,istate) * w0_2                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )    &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) &
                                 - energy(cstate,ispin) + energy(istate,ispin) - energy(bstate,ispin) + ieta )


                  sigma_gwgw0g(iomega,pstate,ispin) = sigma_gwgw0g(iomega,pstate,ispin) &
                           + bra_s(cstate,pstate) * bra_s(bstate,istate) * w0_2                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(bstate,ispin) &
                                 - energy(cstate,ispin) + energy(istate,ispin) + ieta )  &
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
              do pstate=nsemin,nsemax

                !v_2 = eri_eigen(astate,kstate,ispin,jstate,pstate,ispin)
                ak(:) = eri_3center_eigen(:,astate,kstate,ispin)
                jp(:) = eri_3center_eigen(:,jstate,pstate,ispin)
                w0_2 = DOT_PRODUCT( ak(:) , MATMUL( chi_static(:,:) , jp(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,pstate,ispin) = sigma_gwgw0g(iomega,pstate,ispin) &
                           - bra_s(kstate,pstate) * bra_s(astate,jstate) * w0_2                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                             / ( energy(jstate,ispin) - energy(astate,ispin) - pole_s + ieta )

                  sigma_gwgw0g(iomega,pstate,ispin) = sigma_gwgw0g(iomega,pstate,ispin) &
                           + bra_s(kstate,pstate) * bra_s(astate,jstate) * w0_2                          &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                 + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(kstate,ispin) + pole_s - ieta )


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
              do pstate=nsemin,nsemax

                !v_2 = eri_eigen(astate,cstate,ispin,jstate,pstate,ispin)
                ac(:) = eri_3center_eigen(:,astate,cstate,ispin)
                jp(:) = eri_3center_eigen(:,jstate,pstate,ispin)
                w0_2 = DOT_PRODUCT( ac(:) , MATMUL( chi_static(:,:) , jp(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,pstate,ispin) = sigma_gwgw0g(iomega,pstate,ispin) &
                           + bra_s(cstate,pstate) * bra_s(astate,jstate) * w0_2             &
                             / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )  &
                             / ( pole_s + energy(astate,ispin) - energy(jstate,ispin) - ieta )

                enddo
              enddo

            enddo
          enddo
        enddo



      enddo !spole
    enddo !ispin

    call ortho%sum(sigma_gwgw0g)

  endif


  write(stdout,'(a)') ' Sigma_c(omega) is calculated'



  select case(calc_type%selfenergy_approx)
  case(GW0GW0G)
    sigma_gvgw0g(:,:,:) = (0.0_dp, 0.0_dp)
    forall(pstate=nsemin:nsemax)
      se%sigma(:,pstate,:) = sigma_gw0gw0g(:,pstate,:)
    end forall
    write(stdout,'(/,a)') ' GW0GW0G self-energy contributions at E0 (eV)'
  case(GWGW0G)
    forall(pstate=nsemin:nsemax)
      se%sigma(:,pstate,:) = 2.0_dp * sigma_gvgw0g(:,pstate,:) &
                           - sigma_gw0gw0g(:,pstate,:) + 2.0_dp * sigma_gwgw0g(:,pstate,:)
    end forall
    write(stdout,'(/,a)') ' GWGW0G self-energy contributions at E0 (eV)'
  case default
    call die('gwgw0g_selfenergy: calculation type unknown')
  end select


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_gvgw0g'    ,energy,exchange_m_vxc_diag,occupation,energy,sigma_gvgw0g)
  !   call write_selfenergy_omega('selfenergy_gwgw0g'  ,energy,exchange_m_vxc_diag,occupation,energy,sigma_gwgw0g)
  ! endif


  write(stdout,'(a)') &
   '   #       E0          GvGW0G    GW0GW0G_occ  GW0GW0G_vir     GW0GW0G     G(W-v)GW0G    G(W-W0)GW0G    Total'

  do pstate=nsemin,nsemax
    write(stdout,'(i4,1x,*(1x,f12.6))') pstate,se%energy0(pstate,:)*Ha_eV,          &
                                         sigma_gvgw0g(0,pstate,:)%re*Ha_eV,  &
                                         sigma_gw0gw0g_occ(0,pstate,:)%re*Ha_eV,   &
                                         sigma_gw0gw0g_vir(0,pstate,:)%re*Ha_eV,   &
                                         sigma_gw0gw0g(0,pstate,:)%re*Ha_eV,  &
                                         sigma_gwgw0g(0,pstate,:)%re*Ha_eV,&
     MERGE((sigma_gwgw0g(0,pstate,1)%re+sigma_gvgw0g(0,pstate,:)%re-sigma_gw0gw0g(0,pstate,1)%re)*Ha_eV, &
           0.0_dp,calc_type%selfenergy_approx==GWGW0G),&
                                         se%sigma(0,pstate,:)%re*Ha_eV
  enddo

  !DEBUG
  !do iomega=-se%nomega,se%nomega
  !  write(stdout,'(2(2x,f14.6))') (energy(nsemin,1)+se%omega(iomega)%re)*Ha_eV,sigma_gw0gw0g(iomega,nsemin,1)%re*Ha_eV
  !enddo



  call clean_deallocate('Temporary array',bra_s)

  call destroy_eri_3center_eigen()

  call clean_deallocate('chi static',chi_static)

  call stop_clock(timing_gwgamma_self)


end subroutine gwgw0g_selfenergy


!=========================================================================
subroutine gwgwg_selfenergy(nstate,basis,occupation,energy,c_matrix,wpol,se)
  implicit none

  integer,intent(in)                 :: nstate
  type(basis_set)                    :: basis
  real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
  type(spectral_function),intent(in) :: wpol
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer                 :: iomega
  complex(dp),allocatable :: sigma_rest(:,:,:)
  complex(dp),allocatable :: sigma_gwgwg(:,:,:,:)
  integer                 :: astate,bstate,cstate
  integer                 :: istate,jstate,kstate,pqspin,spole,tpole
  integer                 :: pstate,qstate
  real(dp),allocatable    :: bra_t(:,:),bra_s(:,:)
  real(dp)                :: vcoul,vcoul1,vcoul2
  real(dp)                :: Omega_s,Omega_t
  complex(dp)             :: denom1,denom2,denom3,denom4,num3
  real(dp)                :: omega,num1,num2,ei,ej,ek,ea,eb,ec
  !=====

  call start_clock(timing_gwgamma_self)

  write(stdout,*)
  write(stdout,*) 'Perform a one-shot GWGWG calculation'

  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
  else
    call die('not implemented')
  endif


  call clean_allocate('Temporary array',bra_s,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
  call clean_allocate('Temporary array',bra_t,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)



  !
  !
  allocate(sigma_gwgwg(-se%nomega:se%nomega,nsemin:nsemax,nspin,6))
  allocate(sigma_rest(-se%nomega:se%nomega,nsemin:nsemax,nspin))

  sigma_gwgwg(:,:,:,:)  = 0.0_dp


  do pqspin=1,nspin

    do spole=1,wpol%npole_reso
      if( MODULO( spole - 1 , ortho%nproc ) /= ortho%rank ) cycle
      write(stdout,*) 'W poles for GWGWG:',spole,' / ',wpol%npole_reso

      Omega_s = wpol%pole(spole)

      if(has_auxil_basis) then
        do pstate=ncore_G+1,nvirtual_G-1
          ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v) in MO
          bra_s(:,pstate)     = MATMUL( wpol%residue_left(:,spole) , eri_3center_eigen(:,:,pstate,pqspin) )
        enddo
      endif

      do tpole=1,wpol%npole_reso

        Omega_t = wpol%pole(tpole)
        if(has_auxil_basis) then
          do pstate=ncore_G+1,nvirtual_G-1
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v) in MO
            bra_t(:,pstate)     = MATMUL( wpol%residue_left(:,tpole) , eri_3center_eigen(:,:,pstate,pqspin) )
          enddo
        endif



        do pstate=nsemin,nsemax
          qstate=pstate
          do iomega=-se%nomega,se%nomega
            omega = se%energy0(pstate,pqspin) + se%omega(iomega)

            !
            ! 000
            !

            ! occ R occ R occ
            !  i  t  j  s  k
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(jstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(istate,jstate)

                  denom1 = omega + Omega_t - ei -2.0_dp*ieta
                  denom2 = omega + Omega_t + Omega_s - ej -3.0_dp*ieta
                  denom3 = omega + Omega_s - ek -2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,1) = sigma_gwgwg(iomega,pstate,pqspin,1) &
                            + num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            !
            ! 001
            !

            ! occ R occ R emp + occ R occ AR emp
            ! + emp AR occ R occ + emp AR occ R occ
            !  i  t  j  s  c
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(jstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(astate,jstate)

                  ! emp R occ R occ
                  denom1 = omega - ek + Omega_s - 2.0_dp*ieta
                  denom2 = omega - ej + Omega_s + Omega_t - 3.0_dp*ieta
                  denom3 = Omega_s + ea  - ej - 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,2) = sigma_gwgwg(iomega,pstate,pqspin,2) &
                            - 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                  ! emp AR occ AR emp
                  !denom1 = omega - ek + Omega_s - 2.0_dp*ieta
                  denom2 = omega - ea - Omega_t + 2.0_dp*ieta
                  !denom3 = Omega_s + ea  - ej - 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,2) = sigma_gwgwg(iomega,pstate,pqspin,2) &
                            + 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo


            !
            ! 010
            !

            ! occ R emp R occ + occ AR emp R occ + occ AR emp AR occ + occ R emp AR occ
            !  i  t  b  s  k
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(istate,bstate)

                  ! occ R emp R occ
                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = omega - ei - ek + eb - 3.0_dp*ieta
                  denom3 = omega + Omega_s - ek - 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,3) = sigma_gwgwg(iomega,pstate,pqspin,3) &
                            - num1 * num2 / denom1 / denom2 / denom3

                  ! occ AR emp R occ
                  ! occ R emp AR occ
                  denom1 = Omega_t + eb - ek - 3.0_dp*ieta
                  !denom2 = omega - ei - ek + eb - 3.0_dp*ieta
                  !denom3 = omega + Omega_s - ek - 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,3) = sigma_gwgwg(iomega,pstate,pqspin,3) &
                            - 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                  ! occ AR emp AR occ
                  num3   = (2.0_dp * eb - ei - ek + Omega_s + Omega_t - 6.0_dp * ieta )
                  denom1 = eb - ei + Omega_s - 3.0_dp*ieta
                  !denom2 = omega - ei - ek + eb - 3.0_dp*ieta
                  denom3 = eb - ek + Omega_t - 3.0_dp*ieta
                  denom4 = omega - eb - Omega_s - Omega_t + 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,3) = sigma_gwgwg(iomega,pstate,pqspin,3) &
                            + num1 * num2 * num3 / denom1 / denom2 / denom3 / denom4

                enddo
              enddo
            enddo

            !
            ! 011
            !

            ! occ AR emp AR emp + occ R  emp AR emp
            ! emp AR emp R  occ + emp AR emp AR occ
            !  i  t  b  s  c
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(istate,bstate)

                  ! occ AR emp AR emp
                  denom1 = omega - eb - Omega_s - Omega_t + 3.0_dp*ieta
                  denom2 = Omega_s - ei + eb - 3.0_dp*ieta
                  denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,4) = sigma_gwgwg(iomega,pstate,pqspin,4) &
                            + 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                  ! occ R emp AR emp
                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  !denom2 = Omega_s - ei + eb - 3.0_dp*ieta
                  !denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,4) = sigma_gwgwg(iomega,pstate,pqspin,4) &
                            - 2.0_dp * num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            !
            ! 101
            !

            ! emp R occ R emp
            ! emp AR occ R emp
            ! emp AR occ AR emp
            ! emp R occ AR emp
            !  a  t  j  s  c
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(jstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(astate,jstate)

                  ! emp R occ R emp
                  num3   = 2.0_dp * ej  - ea - ec - Omega_s - Omega_t + 6.0_dp * ieta 
                  denom1 = Omega_s - ej + ea - 3.0_dp*ieta
                  denom2 = Omega_t - ej + ec - 3.0_dp*ieta
                  denom3 = omega - ea - ec + ej + 3.0_dp*ieta
                  denom4 = omega - ej + Omega_s + Omega_t - 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,5) = sigma_gwgwg(iomega,pstate,pqspin,5) &
                            + num1 * num2 * num3 / denom1 / denom2 / denom3 / denom4

                  ! emp AR occ R  emp
                  ! emp R  occ AR emp
                  denom1 = omega - ea - Omega_t + 2.0_dp*ieta
                  denom2 = Omega_s + ea - ej - 3.0_dp*ieta
                  !denom3 = omega - ea  - ec + ej + 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,5) = sigma_gwgwg(iomega,pstate,pqspin,5) &
                            + 2.0_dp * num1 * num2 / denom1 / denom2 / denom3
                  ! emp AR occ AR emp
                  !denom1 = omega - ea - Omega_t + 2.0_dp*ieta
                  denom2 = omega - ec - Omega_s + 2.0_dp*ieta
                  !denom3 = omega - ea - ec + ej + 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,5) = sigma_gwgwg(iomega,pstate,pqspin,5) &
                            - num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo

            !
            ! 111
            !

            ! emp AR emp AR emp
            !  a  t  b  s  c
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(bstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(astate,bstate)

                  denom1 = omega - Omega_t - ea + 2.0_dp*ieta
                  denom2 = omega - eb - Omega_s - Omega_t + 3.0_dp*ieta
                  denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin,6) = sigma_gwgwg(iomega,pstate,pqspin,6) &
                            + num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

          enddo ! iomega
        enddo ! pstate


      enddo !tpole
    enddo !spole
  enddo !pqspin

  call ortho%sum(sigma_gwgwg)


  write(stdout,'(a)') ' Sigma_c(omega) is calculated'


  !
  ! The input sigma contains the GW+SOSEX2 selfenergy
  sigma_rest(:,:,:) = se%sigma(:,:,:)


  forall(pstate=nsemin:nsemax)
    se%sigma(:,pstate,:) = sigma_rest(:,pstate,:) + SUM(sigma_gwgwg(:,pstate,:,:),DIM=3)
  end forall


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_gwgwg',energy,exchange_m_vxc_diag,occupation,energy,sigma_gwgwg)
  ! endif


  write(stdout,'(/,a)') ' GWGWG self-energy contributions at E0 (eV)'
  write(stdout,'(a)') &
     '   #  SigC_G(W(w)-v)G(W(w)-v)G'
  write(stdout,'(a3,7(a13))') ' ','ooo','oov+voo','ovo','ovv+vvo','vov','vvv','tot'

  do pstate=nsemin,nsemax
    write(stdout,'(i4,1x,*(1x,f12.6))') pstate, &
                                        sigma_gwgwg(0,pstate,:,:)%re*Ha_eV,&
                                        SUM(sigma_gwgwg(0,pstate,:,:)%re,DIM=2)*Ha_eV
  enddo
  write(stdout,'(a)') &
     '   #          E0       SigC_GW+SOSEX2 SigC_G(W-v)G(W-v)G SigC_TOT'
  do pstate=nsemin,nsemax
    write(stdout,'(i4,1x,*(1x,f12.6))') pstate,se%energy0(pstate,:)*Ha_eV,          &
                                        sigma_rest(0,pstate,:)%re*Ha_eV,   &
                                        SUM(sigma_gwgwg(0,pstate,:,:)%re,DIM=2)*Ha_eV,&
                                        se%sigma(0,pstate,:)%re*Ha_eV
  enddo

  deallocate(sigma_rest)
  deallocate(sigma_gwgwg)


  call clean_deallocate('Temporary array',bra_s)
  call clean_deallocate('Temporary array',bra_t)

  call destroy_eri_3center_eigen()


  call stop_clock(timing_gwgamma_self)


end subroutine gwgwg_selfenergy


!=========================================================================
subroutine gwgwg_selfenergy_real_grid(basis,energy,occupation,c_matrix,se)
  implicit none

  type(basis_set),intent(in)          :: basis
  real(dp),intent(in)                 :: energy(:,:),occupation(:,:)
  real(dp),intent(in)                 :: c_matrix(:,:,:)
  type(selfenergy_grid),intent(inout) :: se
  !=====
  real(dp)             :: domega
  logical,parameter    :: static_fsos  = .FALSE.
  real(dp)             :: mu
  integer              :: nstate
  integer              :: iomega_sigma,iomega,iomegap,iomegapp
  real(dp)             :: df,braket1,braket2
  integer              :: pstate,qstate,rstate,pqspin
  integer              :: ustate,vstate,wstate
  integer              :: first_omega,last_omega
  complex(dp),allocatable :: sigmagwgwg(:,:,:)
  complex(dp)          :: num1,num2
  type(spectral_function) :: wpol_analytic
  real(dp) :: vw(nauxil_global),up(nauxil_global),uv(nauxil_global),qw(nauxil_global)
  real(dp) :: erpa,egw,omegap,omegapp
  complex(dp) :: ieta_u,ieta_v,ieta_w,g_u,g_v,g_w
  complex(dp),allocatable :: chi_omegap(:,:),chi_omegapp(:,:)
  !=====
  domega = eta * 0.5_dp

  call start_clock(timing_gw_self)

  if( nspin > 1 ) then
    call die('gwgwg_selfenergy_real_grid only for spin restricted')
  endif
  if( .NOT. has_auxil_basis ) then
    call die('gwgwg_selfenergy_real_grid requires an auxiliary basis')
  endif
  first_omega = -se%nomega
  last_omega  = se%nomega
  write(stdout,*) first_omega,last_omega


  write(stdout,'(/,1x,a)') 'GWGWG self-energy on a grid of real frequencies with double integration'
  write(stdout,'(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega,last_omega
    write(stdout,'(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma,se%omega(iomega_sigma)*Ha_eV
  enddo
  write(stdout,'(1x,a)') '========================================================'

  nstate = SIZE(energy,DIM=1)

  mu = 0.5_dp * ( MAXVAL(energy(nhomo_G,:)) + MINVAL(energy(nhomo_G+1,:)) )
  write(stdout,'(1x,a,f12.3)') 'Fermi energy mu (eV): ',mu*Ha_eV


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1,timing=timing_aomo_gw)

  if( analytic_chi_ ) then
    call wpol_analytic%init(nstate,occupation,0,grid=NO_GRID)
    call polarizability(.TRUE.,.TRUE.,basis,occupation,energy,c_matrix,erpa,egw,wpol_analytic)
  else
    call die('gwgwg_selfenergy_real_grid: analytic_chi is compulsory')
  endif

  allocate(chi_omegap(nauxil_global,nauxil_global))
  allocate(chi_omegapp(nauxil_global,nauxil_global))

  allocate(sigmagwgwg(first_omega:last_omega,nsemin:nsemax,nspin))
  sigmagwgwg(:,:,:) = 0.0_dp
  pqspin=1

  do iomegap=-nomega_chi_real,nomega_chi_real
    omegap = domega * iomegap
    write(stdout,'(1x,a,i4,es12.4)') 'External omega loop (eV): ',iomegap,omegap*Ha_eV
    call wpol_analytic%evaluate(omegap,chi_omegap)

    do iomegapp=-nomega_chi_real,nomega_chi_real
      omegapp = domega * (iomegapp-0.5_dp)
      call wpol_analytic%evaluate(omegapp,chi_omegapp)

      do pstate=nsemin,nsemax
        qstate=pstate ! diagonal only
        do ustate=ncore_G+1,nvirtual_G-1
        !do ustate=ncore_G+1,nhomo_G
        !do ustate=nhomo_G+1,nvirtual_G-1
          ieta_u = (0.0,1.0_dp) * SIGN(eta, energy(ustate,pqspin) - mu)
          do vstate=ncore_G+1,nvirtual_G-1
          !do vstate=ncore_G+1,nhomo_G
          !do vstate=nhomo_G+1,nvirtual_G-1
            ieta_v = (0.0,1.0_dp) * SIGN(eta, energy(vstate,pqspin) - mu)
            do wstate=ncore_G+1,nvirtual_G-1
            !do wstate=ncore_G+1,nhomo_G
            !do wstate=nhomo_G+1,nvirtual_G-1
              ieta_w = (0.0,1.0_dp) * SIGN(eta, energy(wstate,pqspin) - mu)

              qw(:) = eri_3center_eigen(:,qstate,wstate,pqspin)
              uv(:) = eri_3center_eigen(:,ustate,vstate,pqspin)
              up(:) = eri_3center_eigen(:,ustate,pstate,pqspin)
              vw(:) = eri_3center_eigen(:,vstate,wstate,pqspin)
              num1 = DOT_PRODUCT( vw(:) , MATMUL( chi_omegap(:,:), up(:) ) )
              num2 = DOT_PRODUCT( uv(:) , MATMUL( chi_omegapp(:,:) , qw(:) ) )

              do iomega_sigma=first_omega,last_omega
                g_u = 1.0_dp / ( se%energy0(pstate,pqspin) + se%omega(iomega_sigma) + omegap &
                                 -energy(ustate,pqspin) + ieta_u )
                g_v = 1.0_dp / ( se%energy0(pstate,pqspin) + se%omega(iomega_sigma) + omegap &
                                 + omegapp -energy(vstate,pqspin) + ieta_v )
                g_w = 1.0_dp / ( se%energy0(pstate,pqspin) + se%omega(iomega_sigma) + omegapp &
                                 -energy(wstate,pqspin) + ieta_w )

                sigmagwgwg(iomega_sigma,pstate,pqspin) = sigmagwgwg(iomega_sigma,pstate,pqspin) &
                                                    -1.0_dp / (2.0_dp * pi)**2 * g_u * num1 * g_v * num2 * g_w &
                                                       * domega**2
              enddo

            enddo
          enddo
        enddo
      enddo


    enddo
  enddo
  call world%sum(sigmagwgwg)
  write(stdout,*) 'Self-energy correction (eV): '
  do pstate=nsemin,nsemax
    do iomega_sigma=first_omega,last_omega
      write(stdout,'(4x,i4,*(4x,f14.6))') pstate,(se%energy0(pstate,1) + se%omega(iomega_sigma)%re)*Ha_eV, &
                                          sigmagwgwg(iomega_sigma,pstate,1)%re*Ha_eV
    enddo
  enddo

  se%sigma(:,:,:) = se%sigma(:,:,:) + sigmagwgwg(:,:,:)
  deallocate(sigmagwgwg)


  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)



end subroutine gwgwg_selfenergy_real_grid


!=========================================================================
subroutine gwgwg_selfenergy_imag_grid(basis,energy,occupation,c_matrix,se)
  implicit none

  type(basis_set),intent(in)          :: basis
  real(dp),intent(in)                 :: energy(:,:),occupation(:,:)
  real(dp),intent(in)                 :: c_matrix(:,:,:)
  type(selfenergy_grid),intent(inout) :: se
  !=====
  logical,parameter    :: static_fsos  = .FALSE.
  integer              :: nstate
  integer              :: iomega_sigma,iomegap,iomegapp
  real(dp)             :: df,braket1,braket2
  integer              :: pstate,qstate,rstate,pqspin
  integer              :: ustate,vstate,wstate
  integer              :: first_omega,last_omega
  complex(dp),allocatable :: sigmagwgwg(:,:,:)
  complex(dp)          :: denom1,denom2,denom3,denoms
  type(spectral_function) :: wpol_imag,wpol_analytic
  !type(spectral_function) :: wpol_one
  real(dp) :: vw(nauxil_global),up(nauxil_global),uv(nauxil_global),qw(nauxil_global)
  real(dp) :: erpa,egw
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('gwgwg_selfenergy_imag_grid: it requires an auxiliary basis')
  endif
  first_omega = LBOUND(se%omega_calc(:),DIM=1)
  last_omega  = UBOUND(se%omega_calc(:),DIM=1)

  call start_clock(timing_gw_self)

  write(stdout,'(/,1x,a)') 'GWGWG self-energy on a grid of imaginary frequencies centered on the HOMO-LUMO gap'
  write(stdout,'(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega,last_omega
    write(stdout,'(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma,se%omega_calc(iomega_sigma)*Ha_eV
  enddo
  write(stdout,'(1x,a)') '========================================================'

  nstate = SIZE(energy,DIM=1)


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1,timing=timing_aomo_gw)

  call wpol_imag%init(nstate,occupation,nomega_chi_imag,grid=IMAGINARY_QUAD)
  if( analytic_chi_ ) then
    call wpol_analytic%init(nstate,occupation,0,grid=NO_GRID)
    call polarizability(.TRUE.,.TRUE.,basis,occupation,energy,c_matrix,erpa,egw,wpol_analytic)
    call clean_allocate('Chi',wpol_imag%chi,nauxil_global,nauxil_global,wpol_imag%nomega,verbose=.FALSE.)
    call wpol_analytic%evaluate(wpol_imag%omega,wpol_imag%chi)
  else
    call wpol_imag%vsqrt_chi_vsqrt_rpa(occupation,energy,c_matrix,low_rank=.FALSE.)
  endif



  allocate(sigmagwgwg(first_omega:last_omega,nsemin:nsemax,nspin))
  sigmagwgwg(:,:,:) = 0.0_dp

  do pqspin=1,nspin
    do pstate=nsemin,nsemax
      qstate = pstate ! only the diagonal


      !
      ! First imaginary axis integral: i omega'
      !
      do iomegap=1,wpol_imag%nomega
        if( MODULO( iomegap - 1 , ortho%nproc) /= ortho%rank ) cycle
        write(stdout,'(1x,a,i4,a,i4)') 'Quadrature on omega prime: ',iomegap,' / ',wpol_imag%nomega


        do ustate=ncore_G+1,nvirtual_G-1
          do vstate=ncore_G+1,nvirtual_G-1
            do wstate=ncore_G+1,nvirtual_G-1



              qw(:) = eri_3center_eigen(:,qstate,wstate,pqspin)
              uv(:) = eri_3center_eigen(:,ustate,vstate,pqspin)
              up(:) = eri_3center_eigen(:,ustate,pstate,pqspin)
              vw(:) = eri_3center_eigen(:,vstate,wstate,pqspin)

              ! v * chi( +/- iw') * v
              braket1 = DOT_PRODUCT( vw(:), MATMUL( wpol_imag%chi(:,:,iomegap), up(:) ) )

              !
              ! Second imaginary axis integral: i omega''
              !
              do iomegapp=1,wpol_imag%nomega


                ! v * chi( +/- iw'') * v
                braket2 = DOT_PRODUCT( uv(:) , MATMUL( wpol_imag%chi(:,:,iomegapp) , qw(:) ) )

                do iomega_sigma=first_omega,last_omega
                  ! +w' +w''
                  denom1 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  - energy(ustate,pqspin)
                  denom2 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  &
                                                       + wpol_imag%omega(iomegapp) - energy(vstate,pqspin)
                  denom3 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegapp) - energy(wstate,pqspin)

                  denoms = 1.0_dp / ( denom1 * denom2 * denom3 )

                  ! -w' +w''
                  denom1 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  - energy(ustate,pqspin)
                  denom2 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  &
                                                       + wpol_imag%omega(iomegapp) - energy(vstate,pqspin)
                  denom3 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegapp) - energy(wstate,pqspin)

                  denoms = denoms + 1.0_dp / ( denom1 * denom2 * denom3 )

                  ! +w' -w''
                  denom1 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  - energy(ustate,pqspin)
                  denom2 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  &
                                                       - wpol_imag%omega(iomegapp) - energy(vstate,pqspin)
                  denom3 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegapp) - energy(wstate,pqspin)

                  denoms = denoms + 1.0_dp / ( denom1 * denom2 * denom3 )

                  ! -w' -w''
                  denom1 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  - energy(ustate,pqspin)
                  denom2 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  &
                                                       - wpol_imag%omega(iomegapp) - energy(vstate,pqspin)
                  denom3 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegapp) - energy(wstate,pqspin)

                  denoms = denoms + 1.0_dp / ( denom1 * denom2 * denom3 )



                  sigmagwgwg(iomega_sigma,pstate,pqspin) = sigmagwgwg(iomega_sigma,pstate,pqspin) &
                                + wpol_imag%weight_quad(iomegap) * wpol_imag%weight_quad(iomegapp)  &
                                   * ( braket1 * braket2 ) * denoms * ( -1.0_dp / (2.0_dp * pi) ) **2

                enddo
              enddo
            enddo

          enddo 

        enddo 
      enddo !iomegap

    enddo
  enddo
  call ortho%sum(sigmagwgwg)

  write(stdout,*) 'G (W-v) G (W-v) G'
  do pstate=nsemin,nsemax
    do iomega_sigma=first_omega,last_omega
      write(stdout,'(2(4x,i4),4x,f16.6,1x,f16.6)')  pstate,iomega_sigma, &
                                                    se%sigma_calc(iomega_sigma,pstate,1)%re*Ha_eV, &
                                                    sigmagwgwg(iomega_sigma,pstate,1)%re*Ha_eV
    enddo
  enddo
  se%sigma_calc(:,:,:) = se%sigma_calc(:,:,:) + sigmagwgwg(:,:,:)

  deallocate(sigmagwgwg)
  call wpol_imag%destroy()

  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)

end subroutine gwgwg_selfenergy_imag_grid


!=========================================================================
subroutine sosex_selfenergy_imag_grid(basis,energy,occupation,c_matrix,se)
  implicit none

  type(basis_set),intent(in)          :: basis
  real(dp),intent(in)                 :: energy(:,:),occupation(:,:)
  real(dp),intent(in)                 :: c_matrix(:,:,:)
  type(selfenergy_grid),intent(inout) :: se
  !=====
  logical,parameter    :: SOSEX_INCLUDING_SOX=.FALSE.
  integer              :: nstate
  integer              :: iomega_sigma,iomegap
  real(dp)             :: braket1,braket2
  integer              :: mstate,rstate,mpspin,iauxil
  integer              :: prange,isignp
  integer              :: first_omega,last_omega
  complex(dp),allocatable :: sigma_sosex(:,:,:)
  complex(dp)          :: denom1,denom2,omega_cmplx
  type(spectral_function) :: wpol_imag,wpol_one,wpol_analytic
  real(dp) :: erpa,egw
  ! DGEMM
  integer :: astate,istate
  real(dp),allocatable :: eri3_i_m(:,:),eri3_i_a(:,:),eri3_r_m(:,:),eri3_r_a(:,:),tmp(:,:),braket1_ri(:,:),braket2_ri(:,:)
  real(dp),allocatable :: eri3_a_m(:,:),eri3_a_i(:,:),eri3_r_i(:,:),braket1_ra(:,:),braket2_ra(:,:)
  real(dp) :: chi_wp(nauxil_global,nauxil_global),chi_wwp(nauxil_global,nauxil_global)
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('sosex_selfenergy_imag_grid: requires an auxiliary basis')
  endif
  first_omega = LBOUND(se%omega_calc(:),DIM=1)
  last_omega  = UBOUND(se%omega_calc(:),DIM=1)

  call start_clock(timing_gw_self)

  write(stdout,'(/,1x,a)') 'SOSEX self-energy on a grid of imaginary frequencies centered on the HOMO-LUMO gap'
  write(stdout,'(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega,last_omega
    write(stdout,'(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma,se%omega_calc(iomega_sigma)*Ha_eV
  enddo
  write(stdout,'(1x,a)') '========================================================'

  if( SOSEX_INCLUDING_SOX ) then
    write(stdout,'(1x,a)') 'Calculate SOSEX including SOX: G*W*G*v*G'
  else
    write(stdout,'(1x,a)') 'Calculate SOSEX excluding SOX: G*(W-v)*G*v*G (SOX is calculated elsewhere)'
  endif

  nstate = SIZE(energy,DIM=1)


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1,timing=timing_aomo_gw)

  !
  ! Initialize wpol_imag any way to obtain the quadrature grid points and weights
  call wpol_imag%init(nstate,occupation,nomega_chi_imag,grid=IMAGINARY_QUAD)
  if( analytic_chi_ ) then
    call wpol_analytic%init(nstate,occupation,0,grid=NO_GRID)
    call polarizability(.TRUE.,.TRUE.,basis,occupation,energy,c_matrix,erpa,egw,wpol_analytic)
  else
    call wpol_imag%vsqrt_chi_vsqrt_rpa(occupation,energy,c_matrix,low_rank=.FALSE.)
  endif


  prange = nvirtual_G - ncore_G - 1


  allocate(sigma_sosex(first_omega:last_omega,nsemin:nsemax,nspin))
  sigma_sosex(:,:,:) = 0.0_dp

  do mpspin=1,nspin
    do mstate=nsemin,nsemax


      !
      ! Imaginary axis integral
      !
      do iomegap=1,wpol_imag%nomega
        if( analytic_chi_ ) then
          call wpol_analytic%evaluate(wpol_imag%omega(iomegap),chi_wp)
        else
          chi_wp(:,:) = wpol_imag%chi(:,:,iomegap)
        endif
        if( SOSEX_INCLUDING_SOX ) then
          do iauxil=1,nauxil_global
            chi_wp(iauxil,iauxil) = chi_wp(iauxil,iauxil) + 1.0_dp
          enddo
        endif
        write(stdout,'(1x,a,i4,a,i4)') 'Quadrature on omega prime: ',iomegap,' / ',wpol_imag%nomega

        ! TODO eliminate the isignp loop and move iomega_sigma to become the most inner one
        !
        ! positive and negative omega'
        !
        do isignp=1,-1,-2

          !
          ! loop on Sigma(omega)
          !
          do iomega_sigma=first_omega,last_omega

            omega_cmplx = ABS(se%omega_calc(iomega_sigma)%im+isignp*wpol_imag%omega(iomegap)%im)*im
            chi_wwp(:,:) = 0.0_dp
            do iauxil=1,nauxil_global
              chi_wwp(iauxil,iauxil) = chi_wwp(iauxil,iauxil) + 1.0_dp
            enddo

            allocate(eri3_i_m(nauxil_global,ncore_G+1:nhomo_G))
            allocate(eri3_r_m(nauxil_global,ncore_G+1:nvirtual_G-1))
            allocate(eri3_i_a(nauxil_global,ncore_G+1:nhomo_G))
            allocate(eri3_r_a(nauxil_global,ncore_G+1:nvirtual_G-1))
            allocate(tmp(nauxil_global,ncore_G+1:nhomo_G))
            allocate(braket1_ri(ncore_G+1:nvirtual_G-1,ncore_G+1:nhomo_G))
            allocate(braket2_ri(ncore_G+1:nvirtual_G-1,ncore_G+1:nhomo_G))
            eri3_i_m(:,:) = eri_3center_eigen(:,ncore_G+1:nhomo_G,mstate,mpspin)
            eri3_r_m(:,:) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,mstate,mpspin)

            do astate=nhomo_G+1,nvirtual_G-1
              if( MODULO( astate - (nhomo_G+1) , ortho%nproc) /= ortho%rank ) cycle

              eri3_i_a(:,:) = eri_3center_eigen(:,ncore_G+1:nhomo_G,astate,mpspin)
              eri3_r_a(:,:) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,astate,mpspin)

              !
              ! Fix mstate and astate and then use BLAS level 3 for rstate, istate
              !

              !
              ! Chemist's notation:
              ! ( phi_r phi_m | W( +/- iw') | phi_i phi_a )
              call DGEMM('N','N',nauxil_global,nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp,chi_wp(:,:),nauxil_global, &
                                eri3_i_a,nauxil_global, &
                         0.0_dp,tmp,nauxil_global)
              call DGEMM('T','N',nvirtual_G-ncore_G-1,nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp,eri3_r_m(:,:),nauxil_global, &
                                tmp(:,:),nauxil_global, &
                         0.0_dp,braket1_ri(:,:),nvirtual_G-ncore_G-1)

              !
              ! Chemist's notation:
              ! ( phi_r phi_a | W( iw +/- iw') | phi_i phi_m )
              call DGEMM('N','N',nauxil_global,nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp,chi_wwp(:,:),nauxil_global, &
                                eri3_i_m,nauxil_global, &
                         0.0_dp,tmp,nauxil_global)
              call DGEMM('T','N',nvirtual_G-ncore_G-1,nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp,eri3_r_a(:,:),nauxil_global, &
                                tmp(:,:),nauxil_global, &
                         0.0_dp,braket2_ri(:,:),nvirtual_G-ncore_G-1)


              do rstate=ncore_G+1,nvirtual_G-1
                do istate=ncore_G+1,nhomo_G

                  braket1 = braket1_ri(rstate,istate)
                  braket2 = braket2_ri(rstate,istate)

                  denom1 = se%omega_calc(iomega_sigma) + isignp * wpol_imag%omega(iomegap) - energy(rstate,mpspin)
                  denom2 = isignp * wpol_imag%omega(iomegap) + energy(istate,mpspin) - energy(astate,mpspin)

                  sigma_sosex(iomega_sigma,mstate,mpspin) = sigma_sosex(iomega_sigma,mstate,mpspin) &
                                + wpol_imag%weight_quad(iomegap) &
                                    * braket1 / denom1 * braket2 / denom2 / (2.0_dp * pi)

                enddo
              enddo

            enddo

            deallocate(eri3_i_m)
            deallocate(eri3_i_a)
            deallocate(eri3_r_a)
            deallocate(tmp)
            deallocate(braket1_ri)
            deallocate(braket2_ri)


            allocate(eri3_a_m(nauxil_global,nhomo_G+1:nvirtual_G-1))
            allocate(eri3_a_i(nauxil_global,nhomo_G+1:nvirtual_G-1))
            allocate(eri3_r_i(nauxil_global,ncore_G+1:nvirtual_G-1))
            allocate(braket1_ra(ncore_G+1:nvirtual_G-1,nhomo_G+1:nvirtual_G-1))
            allocate(braket2_ra(ncore_G+1:nvirtual_G-1,nhomo_G+1:nvirtual_G-1))
            allocate(tmp(nauxil_global,nhomo_G+1:nvirtual_G-1))
            eri3_a_m(:,:) = eri_3center_eigen(:,nhomo_G+1:nvirtual_G-1,mstate,mpspin)

            do istate=ncore_G+1,nhomo_G
              if( MODULO( istate - (ncore_G+1) , ortho%nproc) /= ortho%rank ) cycle

              eri3_a_i(:,:) = eri_3center_eigen(:,nhomo_G+1:nvirtual_G-1,istate,mpspin)
              eri3_r_i(:,:) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,istate,mpspin)

              !
              ! Fix mstate and istate and then use BLAS level 3 for rstate, astate
              !

              !
              ! Chemist's notation:
              ! ( phi_r phi_m | W( +/-iw') | phi_i phi_a )
              call DGEMM('N','N',nauxil_global,nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp,chi_wp(:,:),nauxil_global, &
                                eri3_a_i(:,:),nauxil_global, &
                         0.0_dp,tmp(:,:),nauxil_global)
              call DGEMM('T','N',nvirtual_G-ncore_G-1,nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp,eri3_r_m(:,:),nauxil_global, &
                                tmp(:,:),nauxil_global, &
                         0.0_dp,braket1_ra(:,:),nvirtual_G-ncore_G-1)
              !braket1_ra(:,:) = MATMUL( TRANSPOSE(eri3_r_m(:,ncore_G+1:nvirtual_G-1)), tmp(:,nhomo_G+1:nvirtual_G-1) )

              ! v + v * chi(iw +/- iw') * v
              !
              ! Chemist's notation:
              ! ( phi_r phi_i | W( iw +/- iw') | phi_a phi_m )
              call DGEMM('N','N',nauxil_global,nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp,chi_wwp(:,:),nauxil_global, &
                                eri3_a_m,nauxil_global, &
                         0.0_dp,tmp,nauxil_global)
              call DGEMM('T','N',nvirtual_G-ncore_G-1,nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp,eri3_r_i(:,:),nauxil_global, &
                                tmp(:,:),nauxil_global, &
                         0.0_dp,braket2_ra(:,:),nvirtual_G-ncore_G-1)
              !braket2_ra(:,:) = MATMUL( TRANSPOSE(eri3_r_i(:,ncore_G+1:nvirtual_G-1)), tmp(:,nhomo_G+1:nvirtual_G-1) )


              do astate=nhomo_G+1,nvirtual_G-1
                do rstate=ncore_G+1,nvirtual_G-1


                  ! v + v * chi( +/- iw') * v
                  !braket1 = DOT_PRODUCT( mr(:), MATMUL( wpol_imag%chi(:,:,iomegap), ai(:) ) ) + DOT_PRODUCT( mr, ai )
                  braket1 = braket1_ra(rstate,astate)
                  braket2 = braket2_ra(rstate,astate)
                  !braket1 = DOT_PRODUCT( mr, ai )  !SOX

                  ! v + v * chi(iw +/- iw') * v
                  !braket2 = DOT_PRODUCT( ir(:) , MATMUL( wpol_one%chi(:,:,1) , ma(:) ) ) + DOT_PRODUCT( ir, ma)
                  !braket2 = DOT_PRODUCT( ir, ma)  !SOX or SOSEX

                  denom1 = se%omega_calc(iomega_sigma) + isignp * wpol_imag%omega(iomegap) - energy(rstate,mpspin)
                  denom2 = isignp * wpol_imag%omega(iomegap) + energy(astate,mpspin) - energy(istate,mpspin)

                  sigma_sosex(iomega_sigma,mstate,mpspin) = sigma_sosex(iomega_sigma,mstate,mpspin) &
                                - wpol_imag%weight_quad(iomegap) &
                                    * braket1 / denom1 * braket2 / denom2 / (2.0_dp * pi)

                enddo
              enddo

            enddo

            deallocate(eri3_a_m)
            deallocate(eri3_a_i)
            deallocate(eri3_r_i)
            deallocate(braket1_ra)
            deallocate(braket2_ra)
            deallocate(tmp)

            deallocate(eri3_r_m)

            if( .NOT. analytic_chi_ ) call wpol_one%destroy(verbose=.FALSE.)
          enddo !iomega_sigma
        enddo ! isignp
      enddo !iomegap

    enddo
  enddo
  call ortho%sum(sigma_sosex)

  se%sigma_calc(:,:,:) = se%sigma_calc(:,:,:) + factor_sosex * sigma_sosex(:,:,:)

  deallocate(sigma_sosex)
  call wpol_imag%destroy()
  if( analytic_chi_ ) call wpol_analytic%destroy()

  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)

end subroutine sosex_selfenergy_imag_grid


!=========================================================================
subroutine sox_selfenergy_imag_grid(basis,energy,occupation,c_matrix,se)
  implicit none

  type(basis_set)                    :: basis
  real(dp),intent(in)                :: occupation(:,:),energy(:,:)
  real(dp),intent(in)                :: c_matrix(:,:,:)
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer                 :: iomega_sigma
  integer                 :: first_omega,last_omega
  complex(dp),allocatable :: sigma_gw(:,:,:)
  complex(dp),allocatable :: sigma_sox(:,:,:)
  integer                 :: astate,bstate,cstate
  integer                 :: istate,jstate,kstate,ispin
  integer                 :: pstate
  real(dp)                :: vcoul1,vcoul2
  !=====

  call start_clock(timing_gwgamma_self)

  if( .NOT. has_auxil_basis ) then
    call die('sex_selfenergy_imag_grid: requires an auxiliary basis')
  endif
  first_omega = LBOUND(se%omega_calc(:),DIM=1)
  last_omega  = UBOUND(se%omega_calc(:),DIM=1)

  call start_clock(timing_gw_self)

  write(stdout,'(/,1x,a)') 'SOX self-energy on a grid of imaginary frequencies centered on the HOMO-LUMO gap'
  write(stdout,'(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega,last_omega
    write(stdout,'(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma,se%omega_calc(iomega_sigma)*Ha_eV
  enddo
  write(stdout,'(1x,a)') '========================================================'


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)



  !
  !
  allocate(sigma_sox(first_omega:last_omega,nsemin:nsemax,nspin))
  allocate(sigma_gw(first_omega:last_omega,nsemin:nsemax,nspin))

  sigma_sox(:,:,:)  = 0.0_dp


  write(stdout,*) 'Calculate SOX'

  do ispin=1,nspin

    !==========================
    do bstate=nhomo_G+1,nvirtual_G-1
      if( MODULO( bstate-(nhomo_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      do istate=ncore_G+1,nhomo_G
        do kstate=ncore_G+1,nhomo_G

          do pstate=nsemin,nsemax

            vcoul1 = eri_eigen(pstate,istate,ispin,bstate,kstate,ispin)
            vcoul2 = eri_eigen(istate,bstate,ispin,kstate,pstate,ispin)

            do iomega_sigma=first_omega,last_omega
              sigma_sox(iomega_sigma,pstate,ispin) = sigma_sox(iomega_sigma,pstate,ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%omega_calc(iomega_sigma) - energy(istate,ispin) - energy(kstate,ispin) + energy(bstate,ispin) )
            enddo
          enddo

        enddo
      enddo
    enddo

    !==========================
    do cstate=nhomo_G+1,nvirtual_G-1
      if( MODULO( cstate-(nhomo_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      do jstate=ncore_G+1,nhomo_G
        do astate=nhomo_G+1,nvirtual_G-1

          do pstate=nsemin,nsemax

            vcoul1 = eri_eigen(pstate,astate,ispin,jstate,cstate,ispin)
            vcoul2 = eri_eigen(astate,jstate,ispin,cstate,pstate,ispin)

            do iomega_sigma=first_omega,last_omega
              sigma_sox(iomega_sigma,pstate,ispin) = sigma_sox(iomega_sigma,pstate,ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%omega_calc(iomega_sigma) - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) )
            enddo
          enddo

        enddo
      enddo
    enddo

  enddo

  call ortho%sum(sigma_sox)
  !
  ! The input sigma contains the GW selfenergy
  sigma_gw(:,:,:) = se%sigma_calc(:,:,:)
  se%sigma_calc(:,:,:) = sigma_gw(:,:,:) + sigma_sox(:,:,:)


  deallocate(sigma_gw,sigma_sox)

  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)

end subroutine sox_selfenergy_imag_grid


end module m_gwgwg_selfenergy
!=========================================================================
