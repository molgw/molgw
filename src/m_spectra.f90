!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the routines to extract the optical spectra, stopping power in linear-response
!
!=========================================================================
#include "molgw.h"
module m_spectra
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_mpi
 use m_scalapack
 use m_cart_to_pure
 use m_inputparam
 use m_basis_set
 !use m_dft_grid
 use m_spectral_function
 use m_hamiltonian_onebody







contains


!=========================================================================
subroutine optical_spectrum(basis,occupation,c_matrix,chi,xpy_matrix,xmy_matrix,eigenvalue)
  implicit none

  type(basis_set),intent(in)         :: basis
  real(dp),intent(in)                :: occupation(:,:),c_matrix(:,:,:)
  type(spectral_function),intent(in) :: chi
  real(dp),intent(in)                :: xpy_matrix(:,:)
  real(dp),intent(in)                :: xmy_matrix(:,:)
  real(dp),intent(in)                :: eigenvalue(:)
  !=====
  integer                            :: nstate,m_x,n_x
  integer                            :: gt
  integer                            :: nexc,iexc
  integer                            :: t_ia,t_jb
  integer                            :: t_ia_global,t_jb_global
  integer                            :: istate,astate,iaspin
  integer                            :: mpspin
  integer                            :: iomega,idir,jdir
  integer,parameter                  :: nomega=600
  complex(dp)                        :: omega(nomega)
  real(dp)                           :: coeff(2*chi%npole_reso),trace
  real(dp)                           :: dynamical_pol(nomega,3,3),photoabsorp_cross(nomega,3,3)
  real(dp)                           :: static_polarizability(3,3)
  real(dp)                           :: oscillator_strength,trk_sumrule,mean_excitation
  real(dp),allocatable               :: dipole_ao(:,:,:),dipole_mo(:,:,:,:)
  real(dp),allocatable               :: residue(:,:)
  integer                            :: dynpolfile
  integer                            :: photocrossfile
  integer                            :: parityi,parityj,reflectioni,reflectionj
  character(len=32)                  :: symsymbol
  character(len=6)                   :: char6
  !=====


  call start_clock(timing_spectrum)
  !
  ! Calculate the spectrum now
  !

  write(stdout,'(/,a)') ' Calculate the optical spectrum'

  gt = get_gaussian_type_tag(basis%gaussian_type)
  nstate = SIZE(c_matrix,DIM=2)
  m_x = SIZE(xpy_matrix,DIM=1)
  n_x = SIZE(xpy_matrix,DIM=2)

  nexc = nexcitation
  if( nexc == 0 ) nexc = chi%npole_reso

  !
  ! First precalculate all the needed dipole in the basis set
  !
  call setup_dipole_ao(basis,dipole_ao)

  !
  ! Get the dipole oscillator strength on states
  allocate(dipole_mo(nstate,nstate,nspin,3))

  do idir=1,3
    do mpspin=1,nspin
      dipole_mo(:,:,mpspin,idir) = MATMUL( TRANSPOSE( c_matrix(:,:,mpspin) ) , &
                                              MATMUL( dipole_ao(:,:,idir) , c_matrix(:,:,mpspin) ) )
    enddo
  enddo

  deallocate(dipole_ao)


  allocate(residue(3,nexc))

  residue(:,:) = 0.0_dp
  do t_ia=1,m_x
    t_ia_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ia)
    istate = chi%transition_table(1,t_ia_global)
    astate = chi%transition_table(2,t_ia_global)
    iaspin = chi%transition_table(3,t_ia_global)

    ! Let use (i <-> j) symmetry to halve the loop
    do t_jb=1,n_x
      t_jb_global = colindex_local_to_global(ipcol_sd,npcol_sd,t_jb)

      if( t_jb_global <= nexc) then
        residue(:,t_jb_global) = residue(:,t_jb_global) &
                     + dipole_mo(istate,astate,iaspin,:) * xpy_matrix(t_ia,t_jb) * SQRT(spin_fact)
      endif
    enddo

  enddo
  call world%sum(residue)

  deallocate(dipole_mo)


  if( is_iomaster .AND. print_yaml_ ) then
    write(unit_yaml,'(/,a)')  'optical spectrum:'
    write(unit_yaml,'(4x,a)') 'excitations:'
    if( is_triplet ) then
      write(unit_yaml,'(8x,a)') 'spin multiplicity: triplet'
    else
      write(unit_yaml,'(8x,a)') 'spin multiplicity: singlet'
    endif
    write(unit_yaml,'(8x,a)') 'energies:'
    write(unit_yaml,'(12x,a)') 'units: eV'
    do iexc=1,nexc
      write(char6,'(i6)') iexc
      write(unit_yaml,'(12x,a6,a,1x,es18.8)') ADJUSTL(char6),':',eigenvalue(iexc) * Ha_eV
    enddo
    write(unit_yaml,'(8x,a)') 'oscillator strengths:'
  endif

  write(stdout,'(/,5x,a)') 'Excitation energies (eV)     Oscil. strengths   [Symmetry] '

  trk_sumrule=0.0_dp
  mean_excitation=0.0_dp
  do t_jb_global=1,nexc
    t_jb = colindex_global_to_local('S',t_jb_global)

    if( is_triplet ) then
      oscillator_strength = 0.0_dp
    else
      oscillator_strength = 2.0_dp/3.0_dp * DOT_PRODUCT(residue(:,t_jb_global),residue(:,t_jb_global)) * eigenvalue(t_jb_global)
    endif
    trk_sumrule = trk_sumrule + oscillator_strength
    mean_excitation = mean_excitation + oscillator_strength * LOG( eigenvalue(t_jb_global) )

    if( is_iomaster .AND. print_yaml_ ) then
      write(char6,'(i6)') t_jb_global
      write(unit_yaml,'(12x,a6,a,1x,es18.8)') ADJUSTL(char6),':',oscillator_strength
    endif

    if(t_jb_global<=30) then

      if( is_triplet ) then
        symsymbol='3'
      else
        symsymbol='1'
      endif

      !
      ! Test the parity in case of molecule with inversion symmetry

      t_ia_global = 0
      do t_ia=1,m_x
        ! t_jb is zero if the proc is not in charge of this process
        if( t_jb /=0 ) then
          if( 0.5_dp * ABS( xpy_matrix(t_ia,t_jb) + xmy_matrix(t_ia,t_jb) ) > 0.1_dp ) then
            t_ia_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ia)
            exit
          endif
        endif
      enddo
      call world%max(t_ia_global)
      if( t_ia_global == 0 ) cycle

      istate = chi%transition_table(1,t_ia_global)
      astate = chi%transition_table(2,t_ia_global)
      iaspin = chi%transition_table(3,t_ia_global)
      if(planar) then
        reflectioni = wfn_reflection(nstate,basis,c_matrix,istate,iaspin)
        reflectionj = wfn_reflection(nstate,basis,c_matrix,astate,iaspin)
        select case(reflectioni*reflectionj)
        case( 1)
          symsymbol=TRIM(symsymbol)//'(A1, B2 or Ap )'
        case(-1)
          symsymbol=TRIM(symsymbol)//'(A2, B1 or App)'
        end select
      endif
      if(inversion) then
        parityi = wfn_parity(nstate,basis,c_matrix,istate,iaspin)
        parityj = wfn_parity(nstate,basis,c_matrix,astate,iaspin)
        select case(parityi*parityj)
        case( 1)
          symsymbol=TRIM(symsymbol)//'g'
        case(-1)
          symsymbol=TRIM(symsymbol)//'u'
        end select
      endif

      write(stdout,'(1x,a,1x,i4.4,a3,2(f18.8,2x),5x,a32)') 'Exc.',t_jb_global,' : ', &
                   eigenvalue(t_jb_global)*Ha_eV,oscillator_strength,symsymbol

      !
      ! Output the transition coefficients
      coeff(:) = 0.0_dp
      do t_ia=1,m_x
        t_ia_global = rowindex_local_to_global('S',t_ia)
        istate = chi%transition_table(1,t_ia_global)
        astate = chi%transition_table(2,t_ia_global)
        if( t_jb /= 0 ) then
          ! Resonant
          coeff(                 t_ia_global) = 0.5_dp * ( xpy_matrix(t_ia,t_jb) + xmy_matrix(t_ia,t_jb) ) / SQRT(2.0_dp)
          ! Anti-Resonant
          coeff(chi%npole_reso + t_ia_global) = 0.5_dp * ( xpy_matrix(t_ia,t_jb) - xmy_matrix(t_ia,t_jb) ) / SQRT(2.0_dp)
        endif
      enddo
      call world%sum(coeff)


      do t_ia_global=1,chi%npole_reso
        istate = chi%transition_table(1,t_ia_global)
        astate = chi%transition_table(2,t_ia_global)
        ! Resonant
        if( ABS(coeff(                   t_ia_global)) > 0.1_dp )  &
          write(stdout,'(8x,i4,a,i4,1x,f12.5)') istate,' -> ',astate,coeff(t_ia_global)
        ! Anti-Resonant
        if( ABS(coeff(chi%npole_reso+t_ia_global)) > 0.1_dp )  &
          write(stdout,'(8x,i4,a,i4,1x,f12.5)') istate,' <- ',astate,coeff(chi%npole_reso+t_ia_global)
      enddo

      write(stdout,*)

    endif
  enddo

  !
  ! For some calculation conditions, the rest of the subroutine is irrelevant
  ! So skip it! Skip it!
  if( is_triplet ) then
    deallocate(residue)
    return
  endif


  !
  ! Calculate the dynamical dipole polarizability
  ! and the static dipole polarizability
  !
  ! Set the frequency mesh
  omega(1)     =MAX( 0.0_dp      ,MINVAL(ABS(eigenvalue(:)))-10.00/Ha_eV)
  omega(nomega)=MIN(50.0_dp/Ha_eV,MAXVAL(ABS(eigenvalue(:)))+10.00/Ha_eV)
  do iomega=2,nomega-1
    omega(iomega) = omega(1) + ( omega(nomega)-omega(1) ) /REAL(nomega-1,dp) * (iomega-1)
  enddo
  ! Add the broadening
  omega(:) = omega(:) + ieta

  dynamical_pol(:,:,:) = 0.0_dp
  static_polarizability(:,:) = 0.0_dp
  do t_ia=1,nexc
    forall(idir=1:3, jdir=1:3)
      dynamical_pol(:,idir,jdir) = dynamical_pol(:,idir,jdir) &
                           + residue(idir,t_ia) * residue(jdir,t_ia) &
                             * ( AIMAG( -1.0_dp  / ( omega(:) - eigenvalue(t_ia) ) ) &
                                 - AIMAG( -1.0_dp  / ( omega(:) + eigenvalue(t_ia) ) ) )
      static_polarizability(idir,jdir) = static_polarizability(idir,jdir) &
                     + 2.0_dp * residue(idir,t_ia) * residue(jdir,t_ia) / eigenvalue(t_ia)
    end forall
  enddo
  !
  ! Get the photoabsorption cross section
  do iomega=1,nomega
    photoabsorp_cross(iomega,:,:) = 4.0_dp * pi * REAL(omega(iomega),dp) / c_speedlight * dynamical_pol(iomega,:,:)
  enddo

  write(stdout,'(/,a)')     ' TRK sum rule: the two following numbers should compare well'
  write(stdout,'(a,f12.6)') ' Sum over oscillator strengths: ',trk_sumrule
  write(stdout,'(a,f12.6)') '   Number of valence electrons: ',SUM( occupation(ncore_W+1:,:) )

  write(stdout,'(/,a,f12.6)') ' Mean excitation energy (eV): ',EXP( mean_excitation / trk_sumrule ) * Ha_eV

  write(stdout,'(/,a)') ' Static dipole polarizability:'
  trace = 0.0_dp
  do idir=1,3
    write(stdout,'(3(4x,f12.6))') static_polarizability(idir,:)
    trace = trace + static_polarizability(idir,idir) / 3.0_dp
  enddo
  write(stdout,'(a,f12.6)') ' Static dipole polarizability trace: ',trace

  if( is_iomaster .AND. print_yaml_ ) then
    write(unit_yaml,'(8x,a,es18.8)') 'trk sum rule: ',trk_sumrule
    write(unit_yaml,'(8x,a,es18.8)') 'mean excitation energy: ',EXP( mean_excitation / trk_sumrule ) * Ha_eV
    write(unit_yaml,'(8x,a)')        'static polarizability:'
    do idir=1,3
      do jdir=1,3
        write(unit_yaml,'(12x,a,es18.8)') '- ',static_polarizability(idir,jdir)
      enddo
    enddo
  endif

  if( is_iomaster ) then

    open(newunit=dynpolfile,file='dynamical_dipole_polarizability.dat',form='formatted')
    open(newunit=photocrossfile,file='photoabsorption_cross_section.dat',form='formatted')
    write(dynpolfile,'(a)') '#  Imaginary part of dynamical dipole polarizability'
    write(dynpolfile,'(a)') '#  omega (eV)   Average     xx    yx    zx    xy    yy    zy    xz    yz    zz'
    write(photocrossfile,'(a)') '#  Imaginary part of dynamical dipole polarizability'
    write(photocrossfile,'(a)') '#  omega (eV)   Average     xx    yx    zx    xy    yy    zy    xz    yz    zz'
    do iomega=1,nomega
      write(dynpolfile,'(11(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,                                      &
                                           (dynamical_pol(iomega,1,1)+dynamical_pol(iomega,2,2)+dynamical_pol(iomega,3,3))/3.0_dp, &
                                           dynamical_pol(iomega,:,:)
      write(photocrossfile,'(11(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,                                      &
                                             ( photoabsorp_cross(iomega,1,1) &
                                              + photoabsorp_cross(iomega,2,2) &
                                              + photoabsorp_cross(iomega,3,3) ) / 3.0_dp, &
                                             photoabsorp_cross(iomega,:,:)
    enddo

    close(dynpolfile)
    close(photocrossfile)

  endif


  deallocate(residue)

  call stop_clock(timing_spectrum)

end subroutine optical_spectrum


!=========================================================================
subroutine stopping_power(basis,c_matrix,chi,xpy_matrix,eigenvalue)
  implicit none

  type(basis_set),intent(in)         :: basis
  real(dp),intent(in)                :: c_matrix(:,:,:)
  type(spectral_function),intent(in) :: chi
  real(dp),intent(in)                :: xpy_matrix(:,:)
  real(dp),intent(in)                :: eigenvalue(chi%npole_reso)
  !=====
  integer                            :: nstate,m_x,n_x
  integer,parameter                  :: nqradial = 500
  real(dp),parameter                 :: dqradial = 0.02_dp
!  integer,parameter                  :: nqradial = 1500
!  real(dp),parameter                 :: dqradial = 0.01_dp
  integer,parameter                  :: nq = nqradial
  integer                            :: gt
  integer                            :: t_ia,t_jb
  integer                            :: t_ia_global,t_jb_global
  integer                            :: nmat
  integer                            :: istate,astate,iaspin
  integer                            :: mpspin
  complex(dp),allocatable            :: gos_ao(:,:),gos_mo(:,:,:,:)
  complex(dp),allocatable            :: gos_tddft(:)
  integer                            :: iqs,iq,iiq,iqradial
  real(dp)                           :: fnq(chi%npole_reso)
  real(dp)                           :: qvec_list(3,nq),wq(nq)
  real(dp)                           :: qvec(3),qq
  real(dp)                           :: bethe_sumrule(nq)
  integer                            :: iv
  real(dp)                           :: stopping_cross_section(nvel_projectile)
  !real(dp)                           :: stopping_exc(nvel_projectile,chi%npole_reso)
  real(dp)                           :: vlist(3,nvel_projectile),vv
  integer                            :: stride
  integer                            :: nq_batch
  integer                            :: fstopping
  !=====


  call start_clock(timing_stopping)

  write(stdout,'(/,a)') ' Calculate the stopping power in a spherical system'
  gt = get_gaussian_type_tag(basis%gaussian_type)
  nstate = SIZE(c_matrix,DIM=2)
  m_x = SIZE(xpy_matrix,DIM=1)
  n_x = SIZE(xpy_matrix,DIM=2)

  if (nspin/=1) then
    msg='no nspin/=1 allowed'
    call issue_warning(msg)
    return
  endif

  stride = nprow_sd * npcol_sd
  write(stdout,*) 'Parallelize GOS calculation over ',stride


  !
  ! Setup the entire list of projectile velocities
  !
  do iv=1,nvel_projectile
    vlist(:,iv) = vel_projectile(:) * iv
  enddo

  !
  ! Setup the entire q-vector list
  !
  do iqradial=1,nqradial
    qq = iqradial * dqradial
    qvec_list(1,iqradial) = 0.0_dp
    qvec_list(2,iqradial) = 0.0_dp
    qvec_list(3,iqradial) = qq
    wq(iqradial)          = dqradial
  enddo

  if( print_yaml_ .AND. is_iomaster ) then
    write(unit_yaml,'(/,a)') 'stopping power:'
    write(unit_yaml,'(4x,a)') 'unit: a.u.'
    write(unit_yaml,'(4x,a)') 'q vectors:'
    do iq=1,nq
      write(unit_yaml,'(8x,a,es16.6,a,es16.6,a,es16.6,a)') '- [',qvec_list(1,iq),' , ',qvec_list(2,iq),' , ',qvec_list(3,iq),']'
    enddo
  endif

  nmat=chi%npole_reso
  allocate(gos_tddft(chi%npole_reso))

  write(stdout,'(/,1x,a,f8.3,a,f8.3)') 'Loop over q-vectors from ',NORM2(qvec_list(:,1)),' to ',NORM2(qvec_list(:,nq))
  write(stdout,'(5x,a,f8.3)') 'with increment:',dqradial
  bethe_sumrule(:) = 0.0_dp
  stopping_cross_section(:) = 0.0_dp
  !stopping_exc(:,:) = 0.0_dp
  do iqs=1,nq,stride

    nq_batch = MIN(nq,iqs+stride-1) - iqs + 1
    allocate(gos_mo(nstate,nstate,nspin,nq_batch))
    gos_mo(:,:,:,:) = 0.0_dp

    do iiq=1,nq_batch
      if( iiq /= iprow_sd + 1 + nprow_sd * ipcol_sd ) cycle
      iq = iqs + iiq - 1
      qvec(:) = qvec_list(:,iq)
      ! Get the gos oscillator strength on states
      call start_clock(timing_tmp1)
      call setup_gos_ao(basis,qvec,gos_ao)
      call stop_clock(timing_tmp1)

      call start_clock(timing_tmp2)
      do mpspin=1,nspin
        gos_mo(:,:,mpspin,iiq) = MATMUL( TRANSPOSE( c_matrix(:,:,mpspin) ) ,  MATMUL( gos_ao(:,:) , c_matrix(:,:,mpspin) ) )
      enddo
      deallocate(gos_ao)
      call stop_clock(timing_tmp2)
    enddo

    call world%sum(gos_mo)

    do iiq=1,nq_batch
      iq = iqs + iiq - 1
      qvec(:) = qvec_list(:,iq)
      call start_clock(timing_tmp3)
      gos_tddft(:) = (0.0_dp,0.0_dp)
      do t_ia=1,m_x
        t_ia_global = rowindex_local_to_global(iprow_sd,nprow_sd,t_ia)
        istate = chi%transition_table(1,t_ia_global)
        astate = chi%transition_table(2,t_ia_global)
        iaspin = chi%transition_table(3,t_ia_global)

        ! Use (i <-> j) symmetry to halve the loop
        do t_jb=1,n_x
          t_jb_global = colindex_local_to_global(ipcol_sd,npcol_sd,t_jb)

          gos_tddft(t_jb_global) = gos_tddft(t_jb_global) &
                       + gos_mo(istate,astate,iaspin,iiq) * xpy_matrix(t_ia,t_jb) * SQRT(spin_fact)
        enddo

      enddo
      call world%sum(gos_tddft)
      call stop_clock(timing_tmp3)


      fnq(:) = 2.0_dp * ABS( gos_tddft(:) )**2 * eigenvalue(:) / SUM( qvec(:)**2 )
      bethe_sumrule(iq) = SUM(fnq(:))
      !do t_ia=1,nmat
      !  write(1234,*) NORM2(qvec),eigenvalue(t_ia),fnq(t_ia)
      !enddo

      write(stdout,'(1x,a,f8.3,a,f12.6)') 'Bethe sumrule for q',NORM2(qvec(:)),':',bethe_sumrule(iq)

      do iv=1,nvel_projectile
        vv = NORM2(vlist(:,iv))
        do t_ia=1,nmat
          if( NORM2(qvec) > eigenvalue(t_ia) / vv )   &
               stopping_cross_section(iv) = stopping_cross_section(iv) + ( 4.0_dp * pi ) / vv**2  &
                                              * fnq(t_ia)  / NORM2(qvec) * wq(iq) !&
          !if( NORM2(qvec) > eigenvalue(t_ia) / vv )   &
          !     stopping_exc(iv,t_ia) = stopping_exc(iv,t_ia) + ( 4.0_dp * pi ) / vv**2  &
          !                                    * fnq(t_ia)  / NORM2(qvec) * wq(iq) !&
        enddo

      enddo
    enddo

    deallocate(gos_mo)

  enddo

  deallocate(gos_tddft)

  if( print_yaml_ .AND. is_iomaster )  then
    write(unit_yaml,'(4x,a)') 'bethe sum rule:'
    do iq=1,nq
      write(unit_yaml,'(8x,a,es16.6,a,es16.6,a)') '- [',NORM2(qvec_list(:,iq)),' , ',bethe_sumrule(iq),']'
    enddo
  endif

  write(stdout,*) 'Electronic stopping cross section: v, S0 (a.u.)'
  open(newunit=fstopping,file='stopping.dat')
  do iv=1,nvel_projectile
    write(stdout,'(1x,a,3(1x,f6.3),a,f12.6)') 'velocity ',vlist(:,iv),' : ',stopping_cross_section(iv)
    write(fstopping,'(4(2x,es18.8))') vlist(:,iv),stopping_cross_section(iv)
  enddo
  write(stdout,*)
  close(fstopping)

  !do t_ia=1,nmat
  !  write(2000+t_ia,*) '#', &
  !      chi%transition_table(1,MAXLOC(ABS(xpy_matrix(:,t_ia)),DIM=1)),&
  !      chi%transition_table(2,MAXLOC(ABS(xpy_matrix(:,t_ia)),DIM=1))
  !enddo
  !do iv=1,nvel_projectile
  !  vv = NORM2(vlist(:,iv))
  !  do t_ia=1,12 ! nmat
  !    write(stdout,'(i6,1x,2(2x,f12.6))') t_ia,vv,stopping_exc(iv,t_ia)
  !  enddo
  !enddo

  if( print_yaml_ .AND. is_iomaster )  then
    write(unit_yaml,'(4x,a)') 'stopping cross section:'
    do iv=1,nvel_projectile
      write(unit_yaml,'(8x,a,es16.6,a,es16.6,a,es16.6,a,es16.6,a)') '- [',vlist(1,iv),' , ', &
                                                                          vlist(2,iv),' , ',vlist(3,iv),' , ', &
                                                                          stopping_cross_section(iv),']'
    enddo
  endif

  call stop_clock(timing_stopping)


end subroutine stopping_power


!=========================================================================
subroutine stopping_power_3d(basis,c_matrix,chi,xpy_matrix,desc_x,eigenvalue)
  implicit none

  integer,intent(in)                 :: desc_x(NDEL)
  type(basis_set),intent(in)         :: basis
  real(dp),intent(in)                :: c_matrix(:,:,:)
  type(spectral_function),intent(in) :: chi
  real(dp),intent(in)                :: xpy_matrix(:,:)
  real(dp),intent(in)                :: eigenvalue(:)
  !=====
  real(dp),parameter                 :: QMAX = 10.0_dp ! cutoff on maximum q norm to save time
  real(dp)                           :: vlist(3,nvel_projectile),vv,v1(3),v2(3),v3(3)
  integer                            :: gt,nstate
  integer                            :: t_ia,t_jb
  integer                            :: nmat
  integer                            :: istate,astate,iaspin
  integer                            :: mpspin
  complex(dp),allocatable            :: gos_ao(:,:),gos_mo(:,:,:)
  complex(dp)                        :: gos_tddft
  real(dp)                           :: fnq
  real(dp)                           :: qvec(3),qq
  integer                            :: iv
  real(dp)                           :: stopping_cross_section(nvel_projectile)
  real(dp)                           :: stopping_exc(nvel_projectile,chi%npole_reso)
  integer                            :: stride
  integer                            :: nq_batch
  integer                            :: fstopping
  integer,parameter                  :: ncostheta=200  ! from 0 to 1
  integer,parameter                  :: nphi=6         ! from 0 to 2pi
  integer                            :: iphi,icostheta
  real(dp)                           :: phi,costheta,dphi,dcostheta
  real(dp),allocatable               :: xpy_matrix_global(:,:)
  !=====


  call start_clock(timing_stopping)

  write(stdout,'(/,a)') ' Calculate the stopping power for a 3d system'
  gt = get_gaussian_type_tag(basis%gaussian_type)
  nstate = SIZE(c_matrix,DIM=2)

  if (nspin/=1) then
    msg='no nspin/=1 allowed'
    call issue_warning(msg)
    return
  endif

  call clean_allocate('temporary non-distributed X+Y matrix',xpy_matrix_global,chi%npole_reso,chi%npole_reso)
  call gather_distributed_copy(desc_x,xpy_matrix,xpy_matrix_global)

  !
  ! Setup the entire velocity list
  !
  do iv=1,nvel_projectile
    vlist(:,iv) = vel_projectile(:) * iv
  enddo

  !
  ! Set the 3 axis (v1, v2, v3)
  ! v3 is along v
  ! v1 and v2 are orthogonal to v
  !
  v3(:) = vel_projectile(:) / NORM2(vel_projectile(:))
  ! if v3 has no component along x, then v1 will be (1 0 0)
  if( ABS(v3(1)) < 1.0e-6_dp ) then
    v1(1) = 1.0_dp
    v1(2) = 0.0_dp
    v1(3) = 0.0_dp
  else
    v1(1) = -v3(2)
    v1(2) = v3(1)
    v1(3) = 0.0_dp
    v1(:) = v1(:) / NORM2(v1(:))
  endif
  call cross_product(v1,v3,v2)
  v2(:) = -v2(:)

  stopping_exc(:,:) = 0.0_dp
  do iv=1,nvel_projectile
    write(stdout,*) 'iv/nvel_projectile',iv,' / ',nvel_projectile
    vv = NORM2(vlist(:,iv))
    do icostheta=1,ncostheta
      costheta  = REAL(icostheta,dp) / REAL(ncostheta,dp)
      dcostheta = 1.0_dp / REAL(ncostheta,dp)
      do iphi=1,nphi
        phi  = 2.0_dp * pi * REAL(iphi-1,dp) / REAL(nphi,dp)
        dphi = 2.0_dp * pi / REAL(nphi,dp)

        ! Poor man parallelization
        if( MODULO( iphi-1 + nphi*(icostheta-1), world%nproc ) /= world%rank ) cycle

        do t_jb=1,chi%npole_reso

          qq = ABS(eigenvalue(t_jb) / ( vv * costheta ))
          !qvec(1) = qq * SQRT( 1 - costheta**2 ) * COS(phi)
          !qvec(2) = qq * SQRT( 1 - costheta**2 ) * SIN(phi)
          !qvec(3) = qq * costheta
          qvec(:) = qq * SQRT( 1 - costheta**2 ) * COS(phi) * v1(:) &
                  + qq * SQRT( 1 - costheta**2 ) * SIN(phi) * v2(:) &
                  + qq * costheta * v3(:)
          if( qq > QMAX ) cycle

          call start_clock(timing_tmp1)
          call setup_gos_ao(basis,qvec,gos_ao)
          call stop_clock(timing_tmp1)

          call start_clock(timing_tmp2)
          allocate(gos_mo(nstate,nstate,nspin))
          gos_mo(:,:,:) = 0.0_dp
          do mpspin=1,nspin
            gos_mo(:,:,mpspin) = MATMUL( TRANSPOSE( c_matrix(:,:,mpspin) ) ,  MATMUL( gos_ao(:,:) , c_matrix(:,:,mpspin) ) )
          enddo
          deallocate(gos_ao)
          call stop_clock(timing_tmp2)
          ! call world%sum(gos_mo)

          call start_clock(timing_tmp3)
          gos_tddft = (0.0_dp,0.0_dp)
          do t_ia=1,chi%npole_reso
            istate = chi%transition_table(1,t_ia)
            astate = chi%transition_table(2,t_ia)
            iaspin = chi%transition_table(3,t_ia)

            gos_tddft = gos_tddft &
                           + gos_mo(istate,astate,iaspin) * xpy_matrix_global(t_ia,t_jb) * SQRT(spin_fact)
          enddo
          call stop_clock(timing_tmp3)
          deallocate(gos_mo)
          fnq = 2.0_dp * ABS( gos_tddft )**2 * eigenvalue(t_jb) / SUM( qvec(:)**2 )

          stopping_exc(iv,t_jb) = stopping_exc(iv,t_jb) + 2.0_dp / vv**2  &
                                              * fnq  * dphi * dcostheta    / ABS(costheta)

        enddo

      enddo
    enddo
  enddo ! velocity

  call world%sum(stopping_exc)
  stopping_cross_section(:) = SUM(stopping_exc(:,:),DIM=2)

  call clean_deallocate('temporary non-distributed X+Y matrix',xpy_matrix_global)

  write(stdout,*) 'Electronic stopping cross section: v, S0 (a.u.)'
  open(newunit=fstopping,file='stopping3d.dat')
  do iv=1,nvel_projectile
    write(stdout,'(1x,a,3(1x,f6.3),a,f12.6)') 'velocity ',vlist(:,iv),' : ',stopping_cross_section(iv)
    write(fstopping,'(4(2x,es18.8))') vlist(:,iv),stopping_cross_section(iv)
  enddo
  write(stdout,*)
  close(fstopping)
  do iv=1,nvel_projectile
    write(stdout,'(*(2x,es18.8))') vlist(:,iv),stopping_cross_section(iv),stopping_exc(iv,1:20)
  enddo


  if( print_yaml_ .AND. is_iomaster )  then
    write(unit_yaml,'(4x,a)') 'stopping cross section 3d:'
    do iv=1,nvel_projectile
      write(unit_yaml,'(8x,a,es16.6,a,es16.6,a,es16.6,a,es16.6,a)') '- [',vlist(1,iv),' , ', &
                                                                          vlist(2,iv),' , ',vlist(3,iv),' , ', &
                                                                          stopping_cross_section(iv),']'
    enddo
  endif

  call stop_clock(timing_stopping)

end subroutine stopping_power_3d


!=========================================================================
end module m_spectra
!=========================================================================
