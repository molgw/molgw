!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the calculation of the GW self-energy and RPA polarizability
! on a grid of imaginary frequencies
!
!=========================================================================
#include "molgw.h"
module m_gw_selfenergy_grid
  use m_definitions
  use m_timing
  use m_warning
  use m_memory
  use m_scalapack
  use m_inputparam
  use m_mpi
  use m_linear_algebra
  use m_basis_set
  use m_spectral_function
  use m_eri_ao_mo
  use m_selfenergy_tools


contains


!=========================================================================
subroutine polarizability_grid_scalapack(basis,occupation,energy,c_matrix,erpa,egw,wpol)
  implicit none

  type(basis_set),intent(in)            :: basis
  real(dp),intent(in)                   :: occupation(:,:)
  real(dp),intent(in)                   :: energy(:,:)
  real(dp),intent(in)                   :: c_matrix(:,:,:)
  real(dp),intent(out)                  :: erpa,egw
  type(spectral_function),intent(inout) :: wpol
  !=====
  integer              :: nstate
  integer              :: iomega
  integer              :: ilocal,jlocal
  integer              :: iglobal,jglobal
  integer              :: t_ia
  integer              :: istate,astate,iaspin
  integer              :: info
  real(dp)             :: docc,de,factor_sqrt
  real(dp),allocatable :: eri3_t(:,:)
  real(dp),allocatable :: chi0(:,:)
  real(dp),allocatable :: one_m_chi0(:,:)
  real(dp),allocatable :: one_m_chi0m1(:,:)
  real(dp)             :: eigval(nauxil_global)
  integer              :: desc_eri3_t(NDEL)
  integer              :: desc_eri3_final(NDEL)
  integer              :: meri3,neri3
#if defined(HAVE_SCALAPACK)
  real(dp),allocatable :: eri3_sca(:,:)
#endif
  !=====

  call start_clock(timing_rpa_dynamic)

  write(stdout,'(/,1x,a)') 'Calculation of RPA polarizability on imaginary axis grid'
#if defined(HAVE_SCALAPACK)
  write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow_sd,' x ',npcol_sd
#endif

  nstate = SIZE(occupation,DIM=1)


  if( wpol%nomega_quad < 1 ) call die('polarizability_grid_sca: nomega_chi_imag input variable should be greater than 1')

  if( .NOT. has_auxil_basis ) then
    call die('dynamical_polarizability_sca requires an auxiliary basis')
  endif



  wpol%nprodbasis = nauxil_local
  wpol%mchi = NUMROC(nauxil_global,block_row,iprow_sd,first_row,nprow_sd)
  wpol%nchi = NUMROC(nauxil_global,block_col,ipcol_sd,first_col,npcol_sd)
  call DESCINIT(wpol%desc_chi,nauxil_global,nauxil_global,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,wpol%mchi),info)
  call clean_allocate('Chi',wpol%chi,wpol%mchi,wpol%nchi,wpol%nomega_quad)

  write(stdout,'(1x,a,i7,a,i7)') 'Matrix sizes   ',nauxil_global,' x ',nauxil_global
  write(stdout,'(1x,a,i7,a,i7)') 'Distributed in ',wpol%mchi,' x ',wpol%nchi

  if( has_auxil_basis ) call calculate_eri_3center_eigen(c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1,timing=timing_aomo_pola)



  !
  ! Get the processor grid included in the input wpol%desc_chi
  meri3 = NUMROC(nauxil_global ,wpol%desc_chi(MB_),iprow_sd,wpol%desc_chi(RSRC_),nprow_sd)
  neri3 = NUMROC(wpol%npole_reso,wpol%desc_chi(NB_),ipcol_sd,wpol%desc_chi(CSRC_),npcol_sd)
  call DESCINIT(desc_eri3_final,nauxil_global,wpol%npole_reso,wpol%desc_chi(MB_),wpol%desc_chi(NB_), &
                wpol%desc_chi(RSRC_),wpol%desc_chi(CSRC_),wpol%desc_chi(CTXT_),MAX(1,meri3),info)

#if defined(HAVE_SCALAPACK)
  call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)
#endif
  call clean_allocate('TMP 3-center MO integrals',eri3_t,nauxil_local,wpol%npole_reso)
  call clean_allocate('Chi0',chi0,wpol%mchi,wpol%nchi)
  call clean_allocate('1-Chi0',one_m_chi0,wpol%mchi,wpol%nchi)
  call clean_allocate('(1-Chi0)**-1',one_m_chi0m1,wpol%mchi,wpol%nchi)

  call DESCINIT(desc_eri3_t,nauxil_global,wpol%npole_reso,MB_eri3_mo,NB_eri3_mo, &
                first_row,first_col,cntxt_eri3_mo,MAX(1,nauxil_local),info)


  erpa = 0.0_dp
  egw  = 0.0_dp
  do iomega=1,wpol%nomega_quad

    write(stdout,'(1x,a,i4,a,i4)') 'Loop on frequencies: ',iomega,' / ',wpol%nomega_quad

    !
    ! First evaluate v^{1/2} \chi_0 v^{1/2}
    !
    ! Loop over resonant transitions
    do t_ia=1,wpol%npole_reso
      istate = wpol%transition_table(1,t_ia)
      astate = wpol%transition_table(2,t_ia)
      iaspin = wpol%transition_table(3,t_ia)

      docc = occupation(istate,iaspin) - occupation(astate,iaspin)
      de   = energy(astate,iaspin)     - energy(istate,iaspin)
      factor_sqrt = SQRT( 2.0_dp * docc * de / ( wpol%omega_quad(iomega)**2 + de**2 ) )

      eri3_t(:,t_ia) = eri_3center_eigen(:,istate,astate,iaspin) * factor_sqrt

    enddo

#if defined(HAVE_SCALAPACK)
    call PDGEMR2D(nauxil_global,wpol%npole_reso,eri3_t,1,1,desc_eri3_t, &
                                                 eri3_sca,1,1,desc_eri3_final,wpol%desc_chi(CTXT_))
#endif

#if defined(HAVE_SCALAPACK)
    call PDSYRK('L','N',nauxil_global,wpol%npole_reso,1.0_dp,eri3_sca,1,1,desc_eri3_final,0.0_dp,chi0,1,1,wpol%desc_chi)
#else
    call DSYRK('L','N',nauxil_global,wpol%npole_reso,1.0_dp,eri3_t,nauxil_global,0.0_dp,chi0,nauxil_global)
#endif
    chi0(:,:) = -chi0(:,:)



    ! Symmetrize chi0
    call symmetrize_matrix_sca('L',nauxil_global,wpol%desc_chi,chi0,wpol%desc_chi,one_m_chi0)


    one_m_chi0(:,:) = -chi0(:,:)
    do jlocal=1,wpol%nchi
      jglobal = colindex_local_to_global_descriptor(wpol%desc_chi,jlocal)
      do ilocal=1,wpol%mchi
        iglobal = rowindex_local_to_global_descriptor(wpol%desc_chi,ilocal)
        if( iglobal == jglobal ) one_m_chi0(ilocal,jlocal) = one_m_chi0(ilocal,jlocal) + 1.0_dp
      enddo
    enddo


    one_m_chi0m1(:,:) = one_m_chi0(:,:)

    ! Diagonalize (1-chi0) in order to have RPA total energy.
    ! might be a bit time-consuming but we only calculate the eigenvalues
    call diagonalize_eigval_sca(postscf_diago_flavor,one_m_chi0m1,wpol%desc_chi,eigval)
    erpa = erpa + SUM( LOG(eigval(:)) + 1.0_dp - eigval(:) ) / (2.0_dp * pi) * wpol%weight_quad(iomega)
    egw  = egw + SUM( -( 1.0_dp - eigval(:) ) / eigval(:) + 1.0_dp - eigval(:) ) / (2.0_dp * pi) * wpol%weight_quad(iomega)

    call invert_sca(wpol%desc_chi,one_m_chi0,one_m_chi0m1)


#if defined(HAVE_SCALAPACK)
    call PDGEMM('N','N',nauxil_global,nauxil_global,nauxil_global, &
                1.0_dp,one_m_chi0m1        ,1,1,wpol%desc_chi,    &
                       chi0                ,1,1,wpol%desc_chi,    &
                0.0_dp,wpol%chi(:,:,iomega),1,1,wpol%desc_chi)
#else
    call DGEMM('N','N',nauxil_global,nauxil_global,nauxil_global, &
               1.0_dp,one_m_chi0m1,nauxil_global, &
                      chi0        ,nauxil_global, &
               0.0_dp,wpol%chi(:,:,iomega),nauxil_global)
#endif


  enddo

#if defined(HAVE_SCALAPACK)
  call clean_deallocate('TMP 3-center MO integrals',eri3_sca)
#endif
  call clean_deallocate('TMP 3-center MO integrals',eri3_t)
  call clean_deallocate('1-Chi0',one_m_chi0)
  call clean_deallocate('(1-Chi0)**-1',one_m_chi0m1)
  call clean_deallocate('Chi0',chi0)

  call destroy_eri_3center_eigen()

  write(stdout,'(/,1x,a,f16.10)') 'RPA correlation energy (Ha): ',erpa
  write(stdout,'(1x,a,f16.10)')   'GW  correlation energy (Ha): ',egw

  call stop_clock(timing_rpa_dynamic)


end subroutine polarizability_grid_scalapack


!=========================================================================
subroutine gw_selfenergy_imag_scalapack(basis,energy,c_matrix,wpol,se)
  implicit none

  type(basis_set),intent(in)          :: basis
  real(dp),intent(in)                 :: energy(:,:)
  real(dp),intent(in)                 :: c_matrix(:,:,:)
  type(spectral_function),intent(in)  :: wpol
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer              :: nstate
  integer              :: iomega_calc,iomega
  integer              :: info
  real(dp),allocatable :: eri3_sca(:,:)
  real(dp),allocatable :: chi_eri3_sca(:,:)
  real(dp)             :: v_chi_v_p
  integer              :: desc_eri3_t(NDEL)
  integer              :: iprow,ipcol,nprow,npcol
  integer              :: desc_eri3_final(NDEL)
  integer              :: meri3,neri3
  integer              :: mstate,pstate,mpspin
  integer              :: prange,plocal
  complex(dp),allocatable :: sigmaigw(:,:,:)
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('gw_selfenergy_imag_scalapack requires an auxiliary basis')
  endif

  call start_clock(timing_gw_self)

  write(stdout,'(/,1x,a)') 'GW self-energy on a grid of imaginary frequencies'
  write(stdout,'(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_calc=1,se%nomega_calc
    write(stdout,'(1x,i4,1x,f14.4,1x,f14.4)') iomega_calc,se%omega_calc(iomega_calc)*Ha_eV
  enddo
  write(stdout,'(1x,a)') '========================================================'

  nstate = SIZE(energy,DIM=1)

  nprow = 1
  npcol = 1
#if defined(HAVE_SCALAPACK)
  ! Get the processor grid included in the input wpol%desc_chi
  call BLACS_GRIDINFO(wpol%desc_chi(CTXT_),nprow,npcol,iprow,ipcol)
  write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow,' x ',npcol
#endif


  if( has_auxil_basis ) call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,nsemin,nsemax,timing=timing_aomo_gw)


  prange = nvirtual_G - ncore_G - 1

  meri3 = NUMROC(nauxil_global,wpol%desc_chi(MB_),iprow,wpol%desc_chi(RSRC_),nprow)
  neri3 = NUMROC(prange        ,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol)
  call DESCINIT(desc_eri3_final,nauxil_global,prange,wpol%desc_chi(MB_),wpol%desc_chi(NB_), &
                wpol%desc_chi(RSRC_),wpol%desc_chi(CSRC_),wpol%desc_chi(CTXT_),MAX(1,meri3),info)

  call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)
  call clean_allocate('TMP 3-center MO integrals',chi_eri3_sca,meri3,neri3)

  call DESCINIT(desc_eri3_t,nauxil_global,prange,MB_eri3_mo,NB_eri3_mo,first_row,first_col,cntxt_eri3_mo, &
                MAX(1,nauxil_local),info)

  allocate(sigmaigw(se%nomega_calc,nsemin:nsemax,nspin))
  sigmaigw(:,:,:) = 0.0_dp

  do mpspin=1,nspin
    do mstate=nsemin,nsemax

#if defined(HAVE_SCALAPACK)
      call PDGEMR2D(nauxil_global,prange,eri_3center_eigen(:,:,mstate,mpspin),1,1,desc_eri3_t, &
                                                                      eri3_sca,1,1,desc_eri3_final,wpol%desc_chi(CTXT_))
#else
      eri3_sca(:,1:prange) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,mstate,mpspin)
#endif


      do iomega=1,wpol%nomega_quad
#if defined(HAVE_SCALAPACK)
        call PDGEMM('N','N',nauxil_global,prange,nauxil_global,     &
                    1.0_dp,wpol%chi(:,:,iomega),1,1,wpol%desc_chi,    &
                           eri3_sca            ,1,1,desc_eri3_final,  &
                    0.0_dp,chi_eri3_sca        ,1,1,desc_eri3_final)
#else
        call DGEMM('N','N',nauxil_global,prange,nauxil_global,  &
                   1.0_dp,wpol%chi(:,:,iomega),nauxil_global,    &
                          eri3_sca            ,nauxil_global,    &
                   0.0_dp,chi_eri3_sca        ,nauxil_global)
#endif

        !$OMP PARALLEL PRIVATE(pstate, v_chi_v_p)
        !$OMP DO REDUCTION(+:sigmaigw)
        do plocal=1,neri3
          pstate = INDXL2G(plocal,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol) + ncore_G

          v_chi_v_p = DOT_PRODUCT( eri3_sca(:,plocal) , chi_eri3_sca(:,plocal) )

          sigmaigw(:,mstate,mpspin) = sigmaigw(:,mstate,mpspin) &
                        - wpol%weight_quad(iomega) &
                            * (  1.0_dp / ( ( se%omega_calc(:) - energy(pstate,mpspin) ) &
                                              + im * wpol%omega_quad(iomega) )   &
                               + 1.0_dp / ( ( se%omega_calc(:) - energy(pstate,mpspin) )  &
                                              - im * wpol%omega_quad(iomega) )  ) &
                           * v_chi_v_p /  (2.0_dp * pi)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

      enddo

    enddo
  enddo
  call world%sum(sigmaigw)

  se%sigma_calc(:,:,:) = sigmaigw(:,:,:)


  deallocate(sigmaigw)
  call clean_deallocate('TMP 3-center MO integrals',eri3_sca)
  call clean_deallocate('TMP 3-center MO integrals',chi_eri3_sca)

  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)

end subroutine gw_selfenergy_imag_scalapack


!=========================================================================
subroutine gw_selfenergy_contour(basis,energy,occupation,c_matrix,se)
  implicit none

  type(basis_set),intent(in)          :: basis
  real(dp),intent(in)                 :: energy(:,:),occupation(:,:)
  real(dp),intent(in)                 :: c_matrix(:,:,:)
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer              :: nstate
  integer              :: iomega,iomega_sigma
  integer              :: info
  real(dp),allocatable :: eri3_sca(:,:)
  real(dp)             :: v_chi_v_p,factor,de,de_max
  integer              :: desc_eri3_t(NDEL)
  integer              :: iprow,ipcol,nprow,npcol
  integer              :: desc_eri3_final(NDEL)
  integer              :: meri3,neri3
  integer              :: mstate,pstate,mpspin
  integer              :: prange,plocal
  integer              :: neig
  real(dp),allocatable :: tmp(:),tmp2(:,:)
  type(chi_type)       :: vchiv_sqrt_tmp
  complex(dp),allocatable :: sigmagw(:,:,:)
  type(spectral_function) :: wpol_imag,wpol_real
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('gw_selfenergy_coutour_scalapack requires an auxiliary basis')
  endif

  call start_clock(timing_gw_self)

  write(stdout,'(/,1x,a)') 'GW self-energy on a grid of real frequencies centered on the gKS energies'
  write(stdout,'(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega=-se%nomega,se%nomega
    write(stdout,'(1x,i4,1x,f14.4,1x,f14.4)') iomega,se%omega(iomega)*Ha_eV
  enddo
  write(stdout,'(1x,a)') '========================================================'

  nstate = SIZE(energy,DIM=1)

  nprow = 1
  npcol = 1
  iprow = 0
  ipcol = 0


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1,timing=timing_aomo_gw)

  call wpol_imag%init(nstate,occupation,nomega_chi_imag,grid=IMAGINARY_QUAD)
  call wpol_imag%vsqrt_chi_vsqrt_rpa(basis,occupation,energy,c_matrix,low_rank=.TRUE.)

  !
  ! Find largest real omega needed so to specify the frequency grid
  de_max = 0.0_dp
  do mpspin=1,nspin
    do mstate=nsemin,nsemax
      do pstate=ncore_G+1,nhomo_G
        do iomega_sigma=-se%nomega,se%nomega
          ! only poles in the first quadrant  \theta( e_p - omega)
          de = energy(pstate,mpspin) - (se%energy0(mstate,mpspin) + se%omega(iomega_sigma)%re)
          if( de < -eta ) cycle
          de_max = MAX(de_max,de)
        enddo
      enddo
      do pstate=nlumo_G,nvirtual_G-1
        do iomega_sigma=-se%nomega,se%nomega
          ! only poles in the third quadrant  \theta( omega - e_p)
          de = se%energy0(mstate,mpspin) + se%omega(iomega_sigma)%re - energy(pstate,mpspin)
          if( de < -eta ) cycle
          de_max = MAX(de_max,de)
        enddo
      enddo
    enddo
  enddo

  write(stdout,'(1x,a,f12.6)') 'Maximum real frequency needed (eV): ',de_max * Ha_eV
  call wpol_real%init(nstate,occupation,nomega_chi_real,grid=REAL_LINEAR,omega_max=de_max)
  do iomega=1,nomega_chi_real
    wpol_real%omega(iomega) = REAL(iomega-1,dp)/REAL(nomega_chi_real-1,dp) * de_max
  enddo
  call wpol_real%vsqrt_chi_vsqrt_rpa(basis,occupation,energy,c_matrix,low_rank=.TRUE.)


  prange = nvirtual_G - ncore_G - 1

  meri3 = nauxil_global
  neri3 = prange

  call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)

  call DESCINIT(desc_eri3_t,nauxil_global,prange,MB_eri3_mo,NB_eri3_mo,first_row,first_col,cntxt_eri3_mo, &
                MAX(1,nauxil_local),info)

  allocate(sigmagw(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  sigmagw(:,:,:) = 0.0_dp

  do mpspin=1,nspin
    do mstate=nsemin,nsemax

      eri3_sca(:,1:prange) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,mstate,mpspin)


      !
      ! Imaginary axis integral
      !
      do iomega=1,wpol_imag%nomega_quad

        neig = SIZE(wpol_imag%vchiv_sqrt(iomega)%matrix(:,:),DIM=2)
        if( neig == 0 ) cycle
        allocate(tmp2(neig,prange))

        call DGEMM('T','N',neig,prange,nauxil_global,  &
                   1.0_dp,wpol_imag%vchiv_sqrt(iomega)%matrix(:,:),nauxil_global,    &
                          eri3_sca            ,nauxil_global,    &
                   0.0_dp,tmp2                ,neig)

        !$OMP PARALLEL PRIVATE(pstate, v_chi_v_p)
        !$OMP DO REDUCTION(+:sigmagw)
        do plocal=1,neri3
          pstate = plocal + ncore_G

          v_chi_v_p = -SUM(tmp2(:,plocal)**2)

          ! Avoid the poles that are exactly at the origin of the contour
          where( ABS( se%energy0(mstate,mpspin) + se%omega(:) - energy(pstate,mpspin) ) > eta )

            sigmagw(:,mstate,mpspin) = sigmagw(:,mstate,mpspin) &
                          - wpol_imag%weight_quad(iomega) &
                              * (  1.0_dp / ( ( se%energy0(mstate,mpspin) + se%omega(:) - energy(pstate,mpspin) ) &
                                                + im * wpol_imag%omega_quad(iomega) )   &
                                 + 1.0_dp / ( ( se%energy0(mstate,mpspin) + se%omega(:) - energy(pstate,mpspin) )  &
                                                - im * wpol_imag%omega_quad(iomega) )  ) &
                             * v_chi_v_p /  (2.0_dp * pi)
          end where

        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        deallocate(tmp2)
      enddo


      !
      ! Residues for occupied states
      !
      do plocal=1,neri3
        pstate = INDXL2G(plocal,wpol_imag%desc_chi(NB_),ipcol,wpol_imag%desc_chi(CSRC_),npcol) + ncore_G
        ! only occupied states  \theta( \mu - e_p)
        if( pstate >= nlumo_G ) cycle

        do iomega_sigma=-se%nomega,se%nomega
          ! only poles in the first quadrant  \theta( e_p - omega)
          de = energy(pstate,mpspin) - (se%energy0(mstate,mpspin) + se%omega(iomega_sigma)%re)

          if( de < -eta ) cycle

          call wpol_real%interpolate_vsqrt_chi_vsqrt(ABS(de),vchiv_sqrt_tmp)


          allocate(tmp(SIZE(vchiv_sqrt_tmp%matrix(:,:),DIM=2)))
          tmp(:) = MATMUL( eri3_sca(:,plocal) , vchiv_sqrt_tmp%matrix(:,:) )

          factor = MERGE( 0.5_dp, 1.0_dp, ABS(de) < eta )

          sigmagw(iomega_sigma,mstate,mpspin) = sigmagw(iomega_sigma,mstate,mpspin) + SUM(tmp(:)**2) * factor

          deallocate(tmp)
          call vchiv_sqrt_tmp%destroy()

        enddo
      enddo

      !
      ! Residues for empty states
      !
      do plocal=1,neri3
        pstate = INDXL2G(plocal,wpol_imag%desc_chi(NB_),ipcol,wpol_imag%desc_chi(CSRC_),npcol) + ncore_G
        ! only empty states  \theta( e_p - \mu)
        if( pstate <= nhomo_G ) cycle

        do iomega_sigma=-se%nomega,se%nomega
          ! only poles in the third quadrant  \theta( \omega - e_p)
          de = (se%energy0(mstate,mpspin) + se%omega(iomega_sigma)%re) - energy(pstate,mpspin)
          if( de < -eta ) cycle

          call wpol_real%interpolate_vsqrt_chi_vsqrt(ABS(de),vchiv_sqrt_tmp)

          allocate(tmp(SIZE(vchiv_sqrt_tmp%matrix(:,:),DIM=2)))
          tmp(:) = MATMUL( eri3_sca(:,plocal) , vchiv_sqrt_tmp%matrix(:,:) )

          factor = MERGE( 0.5_dp, 1.0_dp, ABS(de) < eta )

          sigmagw(iomega_sigma,mstate,mpspin) = sigmagw(iomega_sigma,mstate,mpspin) - SUM(tmp(:)**2) * factor

          deallocate(tmp)
          call vchiv_sqrt_tmp%destroy()

        enddo
      enddo


    enddo
  enddo
  call world%sum(sigmagw)

  se%sigma(:,:,:) = sigmagw(:,:,:)

  deallocate(sigmagw)
  call clean_deallocate('TMP 3-center MO integrals',eri3_sca)
  call wpol_real%destroy()
  call wpol_imag%destroy()

  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)

end subroutine gw_selfenergy_contour


!=========================================================================
subroutine fsos_selfenergy_imag_scalapack(basis,energy,occupation,c_matrix,wpol,se)
  implicit none

  type(basis_set),intent(in)          :: basis
  real(dp),intent(in)                 :: energy(:,:),occupation(:,:)
  real(dp),intent(in)                 :: c_matrix(:,:,:)
  type(spectral_function),intent(in)  :: wpol
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer              :: nstate
  integer              :: iomega,iomegap
  integer              :: info
  real(dp),allocatable :: eri3_sca(:,:)
  real(dp),allocatable :: chi_eri3_sca(:,:)
  real(dp)             :: v_chi_v_p
  integer              :: desc_eri3_t(NDEL)
  integer              :: iprow,ipcol,nprow,npcol
  integer              :: desc_eri3_final(NDEL)
  integer              :: meri3,neri3
  integer              :: mstate,pstate,mpspin,qstate,rstate,jauxil
  integer              :: prange,plocal
  complex(dp),allocatable :: sigmaigw(:,:,:)
  real(dp) :: term1_iwp,term2_iwiwp,delta_occ_fp
  real(dp),allocatable :: ww(:,:,:),weig(:,:)
  integer :: lwork,nsing
  real(dp),allocatable :: a_matrix(:,:),work(:),u(:,:),vt(:,:),sigma(:)
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('fsos_selfenergy_imag_scalapack requires an auxiliary basis')
  endif

  call start_clock(timing_gw_self)

  write(stdout,'(/,1x,a)') 'FSOS self-energy on a grid of imaginary frequencies'
  write(stdout,'(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega=1,se%nomega_calc
    write(stdout,'(1x,i4,1x,f14.4,1x,f14.4)') iomega,se%omega_calc(iomega)*Ha_eV
  enddo
  write(stdout,'(1x,a)') '========================================================'

  nstate = SIZE(energy,DIM=1)

  nprow = 1
  npcol = 1
#if defined(HAVE_SCALAPACK)
  ! Get the processor grid included in the input wpol%desc_chi
  call BLACS_GRIDINFO(wpol%desc_chi(CTXT_),nprow,npcol,iprow,ipcol)
  write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow,' x ',npcol
#endif


  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1,timing=timing_aomo_gw)


  prange = nvirtual_G - ncore_G - 1

  meri3 = NUMROC(nauxil_global,wpol%desc_chi(MB_),iprow,wpol%desc_chi(RSRC_),nprow)
  neri3 = NUMROC(prange        ,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol)
  call DESCINIT(desc_eri3_final,nauxil_global,prange,wpol%desc_chi(MB_),wpol%desc_chi(NB_), &
                wpol%desc_chi(RSRC_),wpol%desc_chi(CSRC_),wpol%desc_chi(CTXT_),MAX(1,meri3),info)

  call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)
  call clean_allocate('TMP 3-center MO integrals',chi_eri3_sca,meri3,neri3)

  call DESCINIT(desc_eri3_t,nauxil_global,prange,MB_eri3_mo,NB_eri3_mo,first_row,first_col,cntxt_eri3_mo, &
                MAX(1,nauxil_local),info)

  allocate(sigmaigw(se%nomega_calc,nsemin:nsemax,nspin))
  sigmaigw(:,:,:) = 0.0_dp


  ! transform v^1/2* chi * v^1/2 -> 1 + v^1/2* chi * v^1/2
  allocate(ww,SOURCE=wpol%chi)
  allocate(weig(nauxil_global,wpol%nomega_quad))
  !do jauxil=1,SIZE(wpol%chi(:,:,:),DIM=2)
  !  ww(jauxil,jauxil,:) = ww(jauxil,jauxil,:) + 1.0_dp
  !enddo
  do iomegap=1,wpol%nomega_quad
    call diagonalize(' ',ww(:,:,iomegap),weig(:,iomegap))
    write(*,*) 'iomegap',iomegap
    write(*,*) COUNT( ABS(weig(:,iomegap)) < 1.0e-6_dp ),' / ',nauxil_local
    write(*,*) weig(1:4,iomegap)
    write(*,*) weig(SIZE(weig,DIM=1)-4:SIZE(weig,DIM=1),iomegap)
  enddo

  write(*,*) 'Dimensions:'
  write(*,*) nauxil_global,prange**2
  write(*,*)
  allocate(a_matrix(nauxil_global,prange**2))
  allocate(sigma(nauxil_global))
  allocate(u(nauxil_global,nauxil_global),vt(nauxil_global,prange**2))
  a_matrix(:,:) = RESHAPE( eri_3center_eigen(:,:,:,1) , [nauxil_global,prange**2] )
  write(*,*) '=========A'
  write(*,*) MAXVAL(ABS(a_matrix(:,:)))
  write(*,*) a_matrix(1:6,1:6)
  lwork = -1
  allocate(work(1))
  call DGESVD('S', 'S', nauxil_global, prange**2, a_matrix, nauxil_global, sigma, u, nauxil_global, vt, nauxil_global, &
              work, lwork, info )
  lwork = NINT(work(1))
  deallocate(work)
  allocate(work(lwork))
  call DGESVD('S', 'S', nauxil_global, prange**2, a_matrix, nauxil_global, sigma, u, nauxil_global, vt, nauxil_global, &
              work, lwork, info )
  deallocate(work)
  write(*,*) info
  write(*,*) sigma(:)
  nsing = COUNT(sigma(:) > 5.0e-1_dp)
  write(*,*) nsing, ' / ',nauxil_global

  do jauxil=1,nauxil_global
    if( jauxil <= nsing ) then
      u(:,jauxil) = u(:,jauxil) * sigma(jauxil)
    else
      u(:,jauxil) = 0.0_dp
    endif
  enddo
  write(*,*) '=========~~~A'
  a_matrix(:,:) = MATMUL( u , vt )
  write(*,*) MAXVAL(ABS(a_matrix(:,:)))
  write(*,*) a_matrix(1:6,1:6)
  stop 'ENOUGFH'

  do mpspin=1,nspin
    do mstate=nsemin,nsemax
      do iomegap=1,wpol%nomega_quad
        write(*,*) iomegap,wpol%omega_quad(iomegap)
        do iomega=1,se%nomega_calc

          write(*,*) iomega,se%omega_calc(iomega)
          do pstate=ncore_G+1,nvirtual_G-1
            do qstate=ncore_G+1,nvirtual_G-1
              do rstate=ncore_G+1,nvirtual_G-1

                delta_occ_fp = occupation(pstate,mpspin) - occupation(qstate,mpspin)
                term1_iwp = DOT_PRODUCT( eri_3center_eigen(:,qstate,mstate,mpspin) , &
                                     MATMUL( wpol%chi(:,:,iomegap) , eri_3center_eigen(:,rstate,pstate,mpspin) ) )
       !         term2_iwwp = DOT_PRODUCT( eri_3center_eigen(:,qstate,mstate,mpspin) , &



!        !$OMP PARALLEL PRIVATE(pstate, v_chi_v_p)
!        !$OMP DO REDUCTION(+:sigmaigw)
!        do plocal=1,neri3
!          pstate = INDXL2G(plocal,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol) + ncore_G
!
!          v_chi_v_p = DOT_PRODUCT( eri3_sca(:,plocal) , chi_eri3_sca(:,plocal) )
!
!          sigmaigw(:,mstate,mpspin) = sigmaigw(:,mstate,mpspin) &
!                        - wpol%weight_quad(iomega) &
!                            * (  1.0_dp / ( ( se%omega_calc(:) - energy(pstate,mpspin) ) &
!                                              + im * wpol%omega_quad(iomega) )   &
!                               + 1.0_dp / ( ( se%omega_calc(:) - energy(pstate,mpspin) )  &
!                                              - im * wpol%omega_quad(iomega) )  ) &
!                           * v_chi_v_p /  (2.0_dp * pi)
!        enddo
!        !$OMP END DO
!        !$OMP END PARALLEL

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  call world%sum(sigmaigw)

  se%sigma_calc(:,:,:) = sigmaigw(:,:,:)


  deallocate(sigmaigw)
  call clean_deallocate('TMP 3-center MO integrals',eri3_sca)
  call clean_deallocate('TMP 3-center MO integrals',chi_eri3_sca)

  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)

end subroutine fsos_selfenergy_imag_scalapack


!=========================================================================
end module m_gw_selfenergy_grid
!=========================================================================
