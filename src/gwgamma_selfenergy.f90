!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the GW self-energy with vertex function
! within different flavors: GWSOSEX GWSOX
!
!=========================================================================
#include "molgw.h"
#define ILL_ON
subroutine gwgamma_selfenergy(nstate,basis,occupation,energy,c_matrix,wpol,se)
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
  complex(dp),allocatable :: sigma_sosex(:,:,:)
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
  case(GWSOX)
    write(stdout,*) 'Perform a one-shot GW+SOX calculation'
  case(GWSOSEX)
    if( calc_type%selfenergy_technique == EVSC ) then
      write(stdout,*) 'Perform an eigenvalue-self-consistent GW+SOSEX calculation'
    else
      write(stdout,*) 'Perform a one-shot GW+SOSEX calculation'
    endif
  case(GWGWG)
    write(stdout,*) 'Perform a one-shot GW+SOSEX calculation to prepare full GWGWG'
  case default
    call die('gwgamma_selfenergy: calculation type unknown')
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


  call clean_allocate('Temporary array',bra,ncore_G+1,nvirtual_G-1,ncore_G+1,MAX(nhomo_G,nsemax))


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

          do mstate=nsemin,nsemax

            vcoul1 = eri_eigen(mstate,istate,ispin,bstate,kstate,ispin)
            vcoul2 = eri_eigen(istate,bstate,ispin,kstate,mstate,ispin)
            if( gwgamma_tddft_ ) then
              fxc = eval_fxc_rks_singlet(istate,bstate,ispin,kstate,mstate,ispin)
              call grid%sum(fxc)
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
                    / ( energy(mstate,ispin) + se%omega(iomega) &
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

          do mstate=nsemin,nsemax

            vcoul1 = eri_eigen(mstate,astate,ispin,jstate,cstate,ispin)
            vcoul2 = eri_eigen(astate,jstate,ispin,cstate,mstate,ispin)
            if( gwgamma_tddft_ ) then
              fxc = eval_fxc_rks_singlet(astate,jstate,ispin,cstate,mstate,ispin)
              call grid%sum(fxc)
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
                    / ( energy(mstate,ispin) + se%omega(iomega) &
                        - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) + ieta )
            enddo
          enddo

        enddo
      enddo
    enddo


  enddo

  call ortho%sum(sigma_sox)


  if( calc_type%selfenergy_approx == GWSOSEX .OR. calc_type%selfenergy_approx == GWGWG ) then

    write(stdout,*) 'Calculate dynamical SOSEX'


    do ispin=1,nspin

      do spole=1,wpol%npole_reso

        if( MODULO( spole - 1 , ortho%nproc ) /= ortho%rank ) cycle
        write(stdout,*) 'SOSEX W poles:',spole,' / ',wpol%npole_reso

        pole_s = wpol%pole(spole)

        if(has_auxil_basis) then
          do mstate=ncore_G+1,MAX(nhomo_G,nsemax)
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
            bra(:,mstate)     = MATMUL( wpol%residue_left(:,spole) , eri_3center_eigen(:,:,mstate,ispin) )
          enddo
          call auxil%sum(bra)
        else
          ! Here just grab the precalculated value
          forall(istate=ncore_G+1:nvirtual_G-1, mstate=ncore_G+1:MAX(nhomo_G,nsemax))
            bra(istate,mstate) = wpol%residue_left(index_prodstate(istate,mstate) &
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
              do mstate=nsemin,nsemax

                vcoul = eri_eigen(istate,kstate,ispin,bstate,mstate,ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(istate,kstate,ispin,bstate,mstate,ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,mstate,ispin) = sigma_sosex(iomega,mstate,ispin) &
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
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,mstate,ispin) = sigma_sosex(iomega,mstate,ispin) &
                           - bra(cstate,mstate) * bra(bstate,istate) * vcoul                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )    &
                             / ( energy(mstate,ispin) + se%omega(iomega) &
                                 - energy(cstate,ispin) + energy(istate,ispin) - energy(bstate,ispin) + ieta )


                  sigma_sosex(iomega,mstate,ispin) = sigma_sosex(iomega,mstate,ispin) &
                           + bra(cstate,mstate) * bra(bstate,istate) * vcoul                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(bstate,ispin) &
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
              do mstate=nsemin,nsemax

                vcoul = eri_eigen(astate,kstate,ispin,jstate,mstate,ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(astate,kstate,ispin,jstate,mstate,ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,mstate,ispin) = sigma_sosex(iomega,mstate,ispin) &
                           - bra(kstate,mstate) * bra(astate,jstate) * vcoul                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                             / ( energy(jstate,ispin) - energy(astate,ispin) - pole_s + ieta )

                  sigma_sosex(iomega,mstate,ispin) = sigma_sosex(iomega,mstate,ispin) &
                           + bra(kstate,mstate) * bra(astate,jstate) * vcoul                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                 + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
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
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega,se%nomega
                  sigma_sosex(iomega,mstate,ispin) = sigma_sosex(iomega,mstate,ispin) &
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

    call ortho%sum(sigma_sosex)
    sigma_sosex(:,:,:) = factor_sosex * sigma_sosex(:,:,:)  !  1.0 in the original paper
                                                            !  but should be 2.0 in reality

  endif


  write(stdout,'(a)') ' Sigma_c(omega) is calculated'


  !
  ! The input sigma contains the GW selfenergy
  sigma_gw(:,:,:) = se%sigma(:,:,:)


  forall(astate=nsemin:nsemax)
    se%sigma(:,astate,:) = sigma_gw(:,astate,:) + sigma_sox(:,astate,:) + sigma_sosex(:,astate,:)
  end forall


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_sox'    ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sox)
  !   call write_selfenergy_omega('selfenergy_sosex'  ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sosex)
  ! endif


  write(stdout,'(/,a)') ' GWSOSEX self-energy contributions at E0 (eV)'
  if(nspin==1) then
    write(stdout,'(a)') &
     '   #          E0             SigC_G0W0                 SigC_SOX                SigC_SOSEX                 SigC_TOT'
  else
    write(stdout,'(a)') &
      '   #                E0                              SigC_G0W0            SigC_SOX' &
       // '               SigC_SOSEX                 SigC_TOT'
  endif

  do astate=nsemin,nsemax
    write(stdout,'(i4,1x,20(1x,f12.6))') astate,energy(astate,:)*Ha_eV,          &
                                         sigma_gw(0,astate,:)*Ha_eV,   &
                                         sigma_sox(0,astate,:)*Ha_eV,  &
                                         sigma_sosex(0,astate,:)*Ha_eV,&
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
subroutine gwgw0g_selfenergy(nstate,basis,occupation,energy,c_matrix,wpol,se)
  use m_definitions
  use m_mpi
  use m_timing
  use m_inputparam
  use m_warning
  use m_basis_set
  use m_spectral_function
  use m_eri_ao_mo
  use m_selfenergy_tools
  implicit none

  integer,intent(in)                 :: nstate
  type(basis_set)                    :: basis
  real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
  type(spectral_function),intent(in) :: wpol
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer                 :: iomega,ibf_auxil
  complex(dp),allocatable :: sigma_gw(:,:,:)
  complex(dp),allocatable :: sigma_gwgw0g(:,:,:)
  complex(dp),allocatable :: sigma_gw0gw0g(:,:,:)
  complex(dp),allocatable :: sigma_gvgw0g(:,:,:)
  integer                 :: astate,bstate,cstate
  integer                 :: istate,jstate,kstate,ispin,spole
  integer                 :: mstate
  real(dp),allocatable    :: bra_s(:,:)
  real(dp)                :: v_1,w0_1,w0_2
  real(dp)                :: pole_s
  real(dp)                :: fxc
  real(dp),allocatable    :: chi_static(:,:)
  real(dp),allocatable    :: mi(:),bk(:),ib(:),km(:)
  real(dp),allocatable    :: ma(:),jc(:),aj(:),cm(:)
  real(dp),allocatable    :: ik(:),ac(:),bm(:),ak(:),jm(:),ic(:)
  !=====

  call start_clock(timing_gwgamma_self)
  if( .NOT. has_auxil_basis ) call die('gwgw0g_selfenergy: not implemented with auxiliary basi')

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
  case(GWGWG)
    write(stdout,*) 'Perform a one-shot GW+SOSEX calculation to prepare full GWGWG'
  case default
    call die('gwgamma_selfenergy: calculation type unknown')
  end select


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)

  call clean_allocate('Temporary array',bra_s,ncore_G+1,nvirtual_G-1,ncore_G+1,MAX(nhomo_G,nsemax))


  call clean_allocate('chi static',chi_static,nauxil_global,nauxil_global)
  call wpol%evaluate((0.0_dp,0.0_dp),chi_static)
  do ibf_auxil=1,nauxil_global
    chi_static(ibf_auxil,ibf_auxil) = chi_static(ibf_auxil,ibf_auxil) + 1.0_dp
  enddo
  allocate(ib(nauxil_global))
  allocate(km(nauxil_global))
  allocate(mi(nauxil_global))
  allocate(bk(nauxil_global))
  allocate(aj(nauxil_global))
  allocate(cm(nauxil_global))
  allocate(ma(nauxil_global))
  allocate(jc(nauxil_global))
  allocate(ac(nauxil_global))
  allocate(ik(nauxil_global))
  allocate(bm(nauxil_global))
  allocate(ak(nauxil_global))
  allocate(jm(nauxil_global))
  allocate(ic(nauxil_global))

  !
  !
  allocate(sigma_gwgw0g(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gvgw0g(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gw0gw0g(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gw(-se%nomega:se%nomega,nsemin:nsemax,nspin))

  sigma_gwgw0g(:,:,:)  = 0.0_dp
  sigma_gvgw0g(:,:,:)  = 0.0_dp
  sigma_gw0gw0g(:,:,:) = 0.0_dp


  write(stdout,*) 'Calculate two static terms analog to SOX'

  do ispin=1,nspin

    !==========================
    do bstate=ncore_G+1,nvirtual_G-1
      if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty) cycle
      if( MODULO( bstate-(ncore_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      do istate=ncore_G+1,nvirtual_G-1
        if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
        do kstate=ncore_G+1,nvirtual_G-1
          if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle

          do mstate=nsemin,nsemax

            !v_1 = eri_eigen(mstate,istate,ispin,bstate,kstate,ispin)
            !v_2 = eri_eigen(istate,bstate,ispin,kstate,mstate,ispin)

            mi(:) = eri_3center_eigen(:,mstate,istate,ispin)
            bk(:) = eri_3center_eigen(:,bstate,kstate,ispin)
            v_1 = DOT_PRODUCT( mi , bk )
            w0_1 = DOT_PRODUCT( mi(:) , MATMUL( chi_static(:,:) , bk(:) ) )

            ib(:) = eri_3center_eigen(:,istate,bstate,ispin)
            km(:) = eri_3center_eigen(:,kstate,mstate,ispin)
            w0_2 = DOT_PRODUCT( ib(:) , MATMUL( chi_static(:,:) , km(:) ) )

            !
            ! calculate only the diagonal !
            do iomega=-se%nomega,se%nomega
              sigma_gvgw0g(iomega,mstate,ispin) = sigma_gvgw0g(iomega,mstate,ispin) &
                  - v_1 * w0_2            &
                    / ( energy(mstate,ispin) + se%omega(iomega) &
                        - energy(istate,ispin) - energy(kstate,ispin) + energy(bstate,ispin) - ieta )
              sigma_gw0gw0g(iomega,mstate,ispin) = sigma_gw0gw0g(iomega,mstate,ispin) &
                  - w0_1 * w0_2            &
                    / ( energy(mstate,ispin) + se%omega(iomega) &
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

          do mstate=nsemin,nsemax

            !v_1 = eri_eigen(mstate,astate,ispin,jstate,cstate,ispin)
            !v_2 = eri_eigen(astate,jstate,ispin,cstate,mstate,ispin)

            ma(:) = eri_3center_eigen(:,mstate,astate,ispin)
            jc(:) = eri_3center_eigen(:,jstate,cstate,ispin)
            v_1 = DOT_PRODUCT( ma , jc )
            w0_1 = DOT_PRODUCT( ma(:) , MATMUL( chi_static(:,:) , jc(:) ) )

            aj(:) = eri_3center_eigen(:,astate,jstate,ispin)
            cm(:) = eri_3center_eigen(:,cstate,mstate,ispin)
            w0_2 = DOT_PRODUCT( aj(:) , MATMUL( chi_static(:,:) , cm(:) ) )

            !
            ! calculate only the diagonal !
            do iomega=-se%nomega,se%nomega
              sigma_gvgw0g(iomega,mstate,ispin) = sigma_gvgw0g(iomega,mstate,ispin) &
                  - v_1 * w0_2            &
                    / ( energy(mstate,ispin) + se%omega(iomega) &
                        - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) + ieta )
              sigma_gw0gw0g(iomega,mstate,ispin) = sigma_gvgw0g(iomega,mstate,ispin) &
                  - w0_1 * w0_2            &
                    / ( energy(mstate,ispin) + se%omega(iomega) &
                        - energy(astate,ispin) - energy(cstate,ispin) + energy(jstate,ispin) + ieta )
            enddo
          enddo

        enddo
      enddo
    enddo


  enddo

  call ortho%sum(sigma_gvgw0g)


  if( calc_type%selfenergy_approx == GWSOSEX .OR. calc_type%selfenergy_approx == GWGWG ) then

    write(stdout,*) 'Calculate dynamical term analog to SOSEX'


    do ispin=1,nspin

      do spole=1,wpol%npole_reso

        if( MODULO( spole - 1 , ortho%nproc ) /= ortho%rank ) cycle
        write(stdout,*) 'SOSEX W poles:',spole,' / ',wpol%npole_reso

        pole_s = wpol%pole(spole)

        if(has_auxil_basis) then
          do mstate=ncore_G+1,MAX(nhomo_G,nsemax)
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
            bra_s(:,mstate)     = MATMUL( wpol%residue_left(:,spole) , eri_3center_eigen(:,:,mstate,ispin) )
          enddo
          call auxil%sum(bra_s)
        else
          ! Here just grab the precalculated value
          forall(istate=ncore_G+1:nvirtual_G-1, mstate=ncore_G+1:MAX(nhomo_G,nsemax))
            bra_s(istate,mstate) = wpol%residue_left(index_prodstate(istate,mstate) &
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
              do mstate=nsemin,nsemax

                !v_2 = eri_eigen(istate,kstate,ispin,bstate,mstate,ispin)
                ik(:) = eri_3center_eigen(:,istate,kstate,ispin)
                bm(:) = eri_3center_eigen(:,bstate,mstate,ispin)
                w0_2 = DOT_PRODUCT( ik(:) , MATMUL( chi_static(:,:) , bm(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,mstate,ispin) = sigma_gwgw0g(iomega,mstate,ispin) &
                           - bra_s(kstate,mstate) * bra_s(bstate,istate) * w0_2                          &
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

                !v_2 = eri_eigen(istate,cstate,ispin,bstate,mstate,ispin)
                ic(:) = eri_3center_eigen(:,istate,cstate,ispin)
                bm(:) = eri_3center_eigen(:,bstate,mstate,ispin)
                w0_2 = DOT_PRODUCT( ic(:) , MATMUL( chi_static(:,:) , bm(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,mstate,ispin) = sigma_gwgw0g(iomega,mstate,ispin) &
                           - bra_s(cstate,mstate) * bra_s(bstate,istate) * w0_2                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )    &
                             / ( energy(mstate,ispin) + se%omega(iomega) &
                                 - energy(cstate,ispin) + energy(istate,ispin) - energy(bstate,ispin) + ieta )


                  sigma_gwgw0g(iomega,mstate,ispin) = sigma_gwgw0g(iomega,mstate,ispin) &
                           + bra_s(cstate,mstate) * bra_s(bstate,istate) * w0_2                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(bstate,ispin) &
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
              do mstate=nsemin,nsemax

                !v_2 = eri_eigen(astate,kstate,ispin,jstate,mstate,ispin)
                ak(:) = eri_3center_eigen(:,astate,kstate,ispin)
                jm(:) = eri_3center_eigen(:,jstate,mstate,ispin)
                w0_2 = DOT_PRODUCT( ak(:) , MATMUL( chi_static(:,:) , jm(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,mstate,ispin) = sigma_gwgw0g(iomega,mstate,ispin) &
                           - bra_s(kstate,mstate) * bra_s(astate,jstate) * w0_2                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
                             / ( energy(jstate,ispin) - energy(astate,ispin) - pole_s + ieta )

                  sigma_gwgw0g(iomega,mstate,ispin) = sigma_gwgw0g(iomega,mstate,ispin) &
                           + bra_s(kstate,mstate) * bra_s(astate,jstate) * w0_2                          &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(kstate,ispin) &
                                 + energy(astate,ispin) - energy(jstate,ispin)  - ieta )  &
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

                !v_2 = eri_eigen(astate,cstate,ispin,jstate,mstate,ispin)
                ac(:) = eri_3center_eigen(:,astate,cstate,ispin)
                jm(:) = eri_3center_eigen(:,jstate,mstate,ispin)
                w0_2 = DOT_PRODUCT( ac(:) , MATMUL( chi_static(:,:) , jm(:) ) )

                do iomega=-se%nomega,se%nomega
                  sigma_gwgw0g(iomega,mstate,ispin) = sigma_gwgw0g(iomega,mstate,ispin) &
                           + bra_s(cstate,mstate) * bra_s(astate,jstate) * w0_2             &
                             / ( energy(mstate,ispin) + se%omega(iomega) - energy(cstate,ispin) - pole_s + ieta )  &
                             / ( pole_s + energy(astate,ispin) - energy(jstate,ispin) - ieta )

                enddo
              enddo

            enddo
          enddo
        enddo



      enddo !spole
    enddo !ispin

    call ortho%sum(sigma_gwgw0g)
    sigma_gwgw0g(:,:,:) = 2.0_dp * sigma_gwgw0g(:,:,:)

  endif


  write(stdout,'(a)') ' Sigma_c(omega) is calculated'


  !
  ! The input sigma contains the GW selfenergy
  sigma_gw(:,:,:) = se%sigma(:,:,:)


  forall(astate=nsemin:nsemax)
    se%sigma(:,astate,:) = sigma_gw(:,astate,:) + sigma_gvgw0g(:,astate,:) + sigma_gwgw0g(:,astate,:)
  end forall


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_gvgw0g'    ,energy,exchange_m_vxc_diag,occupation,energy,sigma_gvgw0g)
  !   call write_selfenergy_omega('selfenergy_gwgw0g'  ,energy,exchange_m_vxc_diag,occupation,energy,sigma_gwgw0g)
  ! endif


  write(stdout,'(/,a)') ' GWSOSEX self-energy contributions at E0 (eV)'
  if(nspin==1) then
    write(stdout,'(a)') &
     '   #          E0             SigC_G0W0                 SigC_SOX                SigC_SOSEX                 SigC_TOT'
  else
    write(stdout,'(a)') &
      '   #                E0                              SigC_G0W0            SigC_SOX' &
       // '               SigC_SOSEX                 SigC_TOT'
  endif

  do astate=nsemin,nsemax
    write(stdout,'(i4,1x,20(1x,f12.6))') astate,energy(astate,:)*Ha_eV,          &
                                         sigma_gw(0,astate,:)*Ha_eV,   &
                                         sigma_gvgw0g(0,astate,:)*Ha_eV,  &
                                         sigma_gwgw0g(0,astate,:)*Ha_eV,&
                                         se%sigma(0,astate,:)*Ha_eV
  enddo




  call clean_deallocate('Temporary array',bra_s)

  call destroy_eri_3center_eigen()

  call clean_deallocate('chi static',chi_static)

  call stop_clock(timing_gwgamma_self)


end subroutine gwgw0g_selfenergy


!=========================================================================
subroutine gwgwg_selfenergy(nstate,basis,occupation,energy,c_matrix,wpol,se)
  use m_definitions
  use m_mpi
  use m_timing
  use m_inputparam
  use m_warning
  use m_basis_set
  use m_spectral_function
  use m_eri_ao_mo
  use m_selfenergy_tools
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
  complex(dp),allocatable :: sigma_gwgwg(:,:,:)
  integer                 :: astate,bstate,cstate
  integer                 :: istate,jstate,kstate,pqspin,spole,tpole
  integer                 :: pstate,qstate
  real(dp),allocatable    :: bra_t(:,:),bra_s(:,:)
  real(dp)                :: vcoul,vcoul1,vcoul2
  real(dp)                :: Omega_s,Omega_t
  complex(dp)             :: denom1,denom2,denom3
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
  allocate(sigma_gwgwg(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(sigma_gw(-se%nomega:se%nomega,nsemin:nsemax,nspin))

  sigma_gwgwg(:,:,:)  = 0.0_dp


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
            omega = energy(pstate,pqspin) + se%omega(iomega)

            !==========================================================
            ! t: RESONANT
            ! s: RESONANT
            !==========================================================

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

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            ! any R emp R emp
            !  u  t  j  s  k
            ! 0

            ! occ R emp R occ
            !  i  t  b  s  k
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(istate,bstate)

                  denom1 = omega + Omega_t - ei - 2.0_dp*ieta
                  denom2 = omega - ei - ek + eb - 3.0_dp*ieta
                  denom3 = omega + Omega_s - ek - 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            - num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            ! occ R occ R emp
            !  i  t  j  s  c
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(jstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(istate,jstate)

#if defined(ILL_ON)
                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = omega - ec + Omega_s               ! ill term
                  denom3 = ej - ec - Omega_t + 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = omega - ec + Omega_s               ! ill term
                  denom3 = omega - ej + Omega_t - 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
#endif

                enddo
              enddo
            enddo

            ! emp R occ R emp
            !  a  t  j  s  c
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(jstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(astate,jstate)

#if defined(ILL_ON)
                  denom1 = omega - ea + ej - ec + 3.0_dp*ieta
                  denom2 = omega - ec + Omega_s               ! ill term
                  denom3 = ej - ec - Omega_t + 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                  denom1 = ej - Omega_s - ea + 3.0_dp*ieta
                  denom2 = omega - ec + Omega_s               ! ill term
                  denom3 = omega - ej + Omega_s + Omega_t - 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            - num1 * num2 / denom1 / denom2 / denom3
#endif

                enddo
              enddo
            enddo


            ! emp R occ R occ
            !  a  t  j  s  k
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(jstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(astate,jstate)

                  denom1 = ej - Omega_s - ea + 3.0_dp*ieta
                  denom2 = omega - ej + Omega_s + Omega_t -3.0_dp*ieta
                  denom3 = omega - ek + Omega_s -2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo

            !==========================================================
            ! t: ANTIRESONANT
            ! s: RESONANT
            !==========================================================

            ! emp AR occ R occ
            !  a   t  j  s  k
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(jstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(astate,jstate)

                  denom1 = omega - Omega_t - ea + 2.0_dp*ieta
                  denom2 = ea - ej + Omega_s -3.0_dp*ieta
                  denom3 = omega - ek + Omega_s -2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo

            ! emp AR emp R occ
            !  a   t  b  s  k
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(bstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(astate,bstate)

                  denom1 = omega - Omega_t - ea + 2.0_dp*ieta
                  denom2 = ek - eb - Omega_t + 3.0_dp*ieta
                  denom3 = omega + Omega_s - ek - 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            ! occ AR emp R occ
            !  i   t  b  s  k
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(istate,bstate)

                  denom1 = omega + eb -ek - ei - 3.0_dp*ieta
                  denom2 = eb - ek + Omega_t - 3.0_dp*ieta
                  denom3 = omega + Omega_s - ek - 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            - num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo


            ! emp AR occ R emp
            !  a   t  j  s  c
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(jstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(astate,jstate)

#if defined(ILL_ON)
                  denom1 = omega - ea - Omega_t - 2.0_dp*ieta
                  denom2 = omega - ea  -ec + ej + 3.0_dp*ieta
                  denom3 = omega - ec + Omega_s  !ill term

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                  denom1 = omega - ea - Omega_t - 2.0_dp*ieta
                  denom2 = ea  - ej + Omega_s - 3.0_dp*ieta
                  denom3 = omega - ec + Omega_s  !ill term

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
#endif

                enddo
              enddo
            enddo

            !==========================================================
            ! t: RESONANT
            ! s: ANTIRESONANT
            !==========================================================

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

                  denom1 = omega - ea + ej - ec + 3.0_dp*ieta
                  denom2 = ej - ec - Omega_t + 3.0_dp*ieta
                  denom3 = omega - ec - Omega_s + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            - num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo

            ! occ R occ AR emp
            !  i  t  j   s  c
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(jstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(istate,jstate)

                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = Omega_t - ej  + ec - 3.0_dp*ieta
                  denom3 = omega - ec - Omega_s + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo

            ! any R occ AR occ
            !  u  t  j   s  j
            !  -> 0

            ! occ R emp AR occ
            !  i  t  b   s  k
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(istate,bstate)

#if defined(ILL_ON)
                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = omega - Omega_s - ek    ! ill term
                  denom3 = ei - eb - Omega_s + 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = omega - Omega_s - ek    ! ill term
                  denom3 = omega - ei - ek + eb - 3.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
#endif

                enddo
              enddo
            enddo

            ! occ R emp AR emp
            !  i  t  b   s  c
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(istate,bstate)

                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = ei - eb - Omega_s + 3.0_dp*ieta
                  denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            !==========================================================
            ! t: ANTIRESONANT
            ! s: ANTIRESONANT
            !==========================================================

            ! emp AR occ AR emp
            !  a   t  j   s  c
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do jstate=ncore_G+1,nhomo_G
                ej = energy(jstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(jstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(astate,jstate)

                  denom1 = omega - ea - Omega_t + 2.0_dp*ieta
                  denom2 = omega - ea - ec + ej + 3.0_dp*ieta
                  denom3 = omega - ec - Omega_s + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            - num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo

            ! occ AR emp AR emp
            !  i  t  b  s  c
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do cstate=nhomo_G+1,nvirtual_G-1
                  ec = energy(cstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,cstate)
                  num2 = bra_s(qstate,cstate) * bra_s(istate,bstate)

                  denom1 = Omega_s - ei + eb - 3.0_dp*ieta
                  denom2 = omega - eb - Omega_s - Omega_t - 3.0_dp*ieta
                  denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

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

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            - num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            ! occ AR emp AR occ
            !  i   t  b   s  k
            do istate=ncore_G+1,nhomo_G
              ei = energy(istate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,istate) * bra_t(bstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(istate,bstate)

#if defined(ILL_ON)
                  denom1 = Omega_s - ei + eb - 3.0_dp*ieta
                  denom2 = omega - eb - Omega_s - Omega_t + 3.0_dp*ieta
                  denom3 = omega - Omega_s - ek    ! ill term

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                  denom1 = omega - ei - ek + eb - 3.0_dp*ieta
                  denom2 = eb - ek + Omega_t - 3.0_dp*ieta
                  denom3 = omega - Omega_s - ek    ! ill term

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3
#endif

                enddo
              enddo
            enddo

            ! emp AR emp R occ
            !  a   t  b  s  k
            do astate=nhomo_G+1,nvirtual_G-1
              ea = energy(astate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eb = energy(bstate,pqspin)
                do kstate=ncore_G+1,nhomo_G
                  ek = energy(kstate,pqspin)
                  num1 = bra_t(pstate,astate) * bra_t(bstate,kstate)
                  num2 = bra_s(qstate,kstate) * bra_s(astate,bstate)

                  denom1 = omega - Omega_t - ea + 2.0_dp*ieta
                  denom2 = omega - eb - Omega_t + 3.0_dp*ieta
                  denom3 = omega - Omega_s - ek    ! ill term

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
                            + num1 * num2 / denom1 / denom2 / denom3

                  denom1 = omega - Omega_t - ea + 2.0_dp*ieta
                  denom2 = Omega_t - ek + eb + 3.0_dp*ieta
                  denom3 = omega - Omega_s - ek    ! ill term

                  sigma_gwgwg(iomega,pstate,pqspin) = sigma_gwgwg(iomega,pstate,pqspin) &
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
  sigma_gw(:,:,:) = se%sigma(:,:,:)


  forall(pstate=nsemin:nsemax)
    se%sigma(:,pstate,:) = sigma_gw(:,pstate,:) + sigma_gwgwg(:,pstate,:)
  end forall


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_gwgwg',energy,exchange_m_vxc_diag,occupation,energy,sigma_gwgwg)
  ! endif


  write(stdout,'(/,a)') ' GWGWG self-energy contributions at E0 (eV)'
  if(nspin==1) then
    write(stdout,'(a)') &
     '   #          E0             SigC_GW+SOSEX2         SigC_GWGWG                 SigC_TOT'
  else
    write(stdout,'(a)') &
      '   #                E0                              SigC_GW+SOSEX2      ' &
       // '               SigC_GWGWG                 SigC_TOT'
  endif

  do pstate=nsemin,nsemax
    write(stdout,'(i4,1x,20(1x,f12.6))') pstate,energy(pstate,:)*Ha_eV,          &
                                         sigma_gw(0,pstate,:)*Ha_eV,   &
                                         sigma_gwgwg(0,pstate,:)*Ha_eV,&
                                         se%sigma(0,pstate,:)*Ha_eV
  enddo




  call clean_deallocate('Temporary array',bra_s)
  call clean_deallocate('Temporary array',bra_t)

  call destroy_eri_3center_eigen()


  call stop_clock(timing_gwgamma_self)


end subroutine gwgwg_selfenergy


!=========================================================================
