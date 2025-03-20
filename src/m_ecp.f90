!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods and data for Effective Core Potentials (ECP)
!
!=========================================================================
#include "molgw.h"
module m_ecp
  use m_definitions
  use m_string_tools, only: capitalize, append_to_list, orbital_momentum_name, orbital_momentum_number
  use m_warning, only: die, issue_warning
  use m_elements
  use ISO_FORTRAN_ENV, only: IOSTAT_END


  integer, parameter :: ECP_NWCHEM = 1
  integer, parameter :: ECP_PSP6   = 2
  integer, parameter :: ECP_PSP8   = 3
  integer, parameter :: ECP_GTH    = 4

  integer, protected                :: nelement_ecp
  integer, protected, allocatable    :: element_ecp(:)

  type effective_core_potential
    integer              :: ecp_format
    integer              :: ncore          ! number of core electrons
    integer              :: necp           ! number of projectors
    integer, allocatable  :: lk(:)          ! angular momentum of the projector (-1 stands for local component)
    ! ECP
    integer, allocatable  :: nk(:)          ! r**(nk-2)
    real(dp), allocatable :: dk(:)          ! dk coefficient
    real(dp), allocatable :: zetak(:)       ! zetak coefficient (gaussian exponent)
    ! Analytic GTH
    real(dp)              :: gth_rpploc      ! r_loc^pp   following Krack's notation [Theor Chem Acc 114, 145 (2005)]
    integer               :: gth_nloc        ! number of Ci^pp coefficients in the local potential
    real(dp), allocatable  :: gth_cipp(:)     ! Ci^pp      following Krack's notation
    integer               :: gth_nl          ! number of angular momenta
    integer, allocatable   :: gth_npl(:)      ! number of projectors for angular momentum l
    real(dp), allocatable  :: gth_rl(:)       ! r^l        following Krack's notation
    real(dp), allocatable  :: gth_hijl(:, :)   ! h_ij^l     following Krack's notation
    ! KB numerical pseudo on a grid
    integer              :: mmax = 0
    real(dp), allocatable :: rad(:)
    real(dp), allocatable :: wfll(:, :)
    real(dp), allocatable :: vpspll(:, :)
    real(dp), allocatable :: ekb(:)
    real(dp), allocatable :: rhocore(:, :)
  end type


  type(effective_core_potential), allocatable :: ecp(:)

  !
  ! Grid quality description
  integer, protected      :: nradial_ecp, nangular_ecp


contains


!=========================================================================
subroutine init_ecp(ecp_elements, ecp_path, ecp_name, ecp_level_in)
  implicit none

  character(len=*), intent(in) :: ecp_elements
  character(len=*), intent(in) :: ecp_path
  character(len=*), intent(in) :: ecp_name
  integer, intent(in)          :: ecp_level_in
  !=====
  character(len=132) :: string, ecp_filename
  character(len=2)   :: element
  character(len=5)   :: amc
  integer :: ilen, inextblank, ielement_ecp, iecp
  logical :: file_exists
  !=====

  !
  ! First parse the ecp_elements line
  !
  ilen = LEN(TRIM(ecp_elements))

  ! ecp_elements is empty, no ECP needs to be setup
  if( ilen == 0 ) then
    nelement_ecp = 0
    return
  endif

  string = ecp_elements
  write(stdout, '(/,1x,a)') 'Reading ECP element list'

  !
  ! Set up the integration grid
  select case(ecp_level_in)
  case(low)       ! accuracy not guaranted, just for quick test runs
    nradial_ecp     =  12
    nangular_ecp    =   6
  case(medium)
    nradial_ecp     =  20
    nangular_ecp    =  26
  case(high)
    nradial_ecp     =  35
    nangular_ecp    =  38
  case(very_high) ! almost perfect potentials
    nradial_ecp     =  50
    nangular_ecp    = 110
  case(insane)    ! overdoing a lot
    nradial_ecp     =  200
    nangular_ecp    =  434
  case default
    call die('integration quality not recognized')
  end select


  do while( ilen > 0 )
    string = ADJUSTL(string)
    inextblank = INDEX(string, ' ')

    call append_to_list(element_number(string(1:inextblank-1)), element_ecp)

    string = string(inextblank+1:)
    ilen = LEN(TRIM(string))

  enddo

  nelement_ecp = SIZE(element_ecp)
  allocate(ecp(nelement_ecp))

  ecp(:)%ecp_format = ECP_NWCHEM
  if( INDEX(capitalize(ecp_name), 'PSP6') /= 0 ) ecp(:)%ecp_format = ECP_PSP6
  if( INDEX(capitalize(ecp_name), 'PSP8') /= 0 ) ecp(:)%ecp_format = ECP_PSP8
  if( INDEX(capitalize(ecp_name), 'GTH') /= 0 )  ecp(:)%ecp_format = ECP_GTH

  !
  ! Second, read the ECP parameters from ECP file
  !
  do ielement_ecp=1, nelement_ecp
    element = element_name(element_ecp(ielement_ecp))
    write(stdout, '(1x,a,a)') 'ECP for element: ', element

    ecp_filename = TRIM(ecp_path)//'/'//TRIM(ADJUSTL(element))//'_'//TRIM(ecp_name)
    inquire(file=TRIM(ecp_filename), exist=file_exists)
    if( .NOT. file_exists ) then
      write(stdout, '(1x,a,a)') 'Looking for file ', TRIM(ecp_filename)
      write(stdout, '(1x,a)')   'Remember the basis directory path is obtained (by priority order) from:'
      write(stdout, '(1x,a)')   '  1. the input variable basis_path'
      write(stdout, '(1x,a)')   '  2. the environment variable MOLGW_BASIS_PATH'
      write(stdout, '(1x,a)')   '  3. the location of the sources'
      call die('init_ecp: ECP file not found')
    endif

    select case(ecp(ielement_ecp)%ecp_format)
    case(ECP_NWCHEM)
      call read_ecp_file(ecp_filename, element, ecp(ielement_ecp))
    case(ECP_PSP6)
      call read_psp6_file(ecp_filename, element, ecp(ielement_ecp))
    case(ECP_PSP8)
      call read_psp8_file(ecp_filename, element, ecp(ielement_ecp))
    case(ECP_GTH)
      call read_gth_file(ecp_filename, element, ecp(ielement_ecp))
    end select

    write(stdout, '(6x,a,i3)') 'Core electrons ', ecp(ielement_ecp)%ncore
    select case(ecp(ielement_ecp)%ecp_format)
    case(ECP_PSP6, ECP_PSP8)
      write(stdout, '(6x,a)') 'l_k    KB energy (Ha)'

      do iecp=1, ecp(ielement_ecp)%necp
        if( ecp(ielement_ecp)%lk(iecp) == -1 ) then
          amc = 'local'
        else
          amc = orbital_momentum_name(ecp(ielement_ecp)%lk(iecp))
        endif
        write(stdout, '(6x,a,3x,f14.6)') amc, ecp(ielement_ecp)%ekb(iecp)
      enddo
    case(ECP_NWCHEM)
      write(stdout, '(6x,a)') 'l_k      n_k       zeta_k          d_k  '

      do iecp=1, ecp(ielement_ecp)%necp
        if( ecp(ielement_ecp)%lk(iecp) == -1 ) then
          amc = 'local'
        else
          amc = orbital_momentum_name(ecp(ielement_ecp)%lk(iecp))
        endif
        write(stdout, '(6x,a,3x,i3,2(2x,f14.6))') &
                            amc, &
                            ecp(ielement_ecp)%nk(iecp), &
                            ecp(ielement_ecp)%zetak(iecp), &
                            ecp(ielement_ecp)%dk(iecp)
      enddo
    end select

  enddo

  select case(ecp(1)%ecp_format)
  case(ECP_PSP6, ECP_PSP8)
    nradial_ecp = MINVAL(ecp(:)%mmax) - 1   ! Remove the last point for safety, usually all the projectors are zero anyway there
  end select

  select case(ecp(1)%ecp_format)
  case(ECP_GTH)
    write(stdout, '(1x,a)') 'ECP are integrated analytically'
  case default
    write(stdout, '(1x,a,i5,2x,i5)') 'ECP are integrated numerically with a grid (radial,angular): ', nradial_ecp, nangular_ecp
  end select

end subroutine init_ecp


!=========================================================================
subroutine read_ecp_file(ecp_filename, element, ecpi)
  implicit none

  character(*), intent(in)                      :: ecp_filename
  character(len=2), intent(in)                  :: element
  type(effective_core_potential), intent(inout) :: ecpi
  !=====
  integer            :: ecpunit
  integer            :: iline, i1, i2, istat
  character(len=132) :: line, amc
  integer            :: read_n
  real(dp)           :: read_zeta, read_d
  logical            :: end_of_file
  !=====

  write(stdout, *) 'read NWCHEM ECP file:', TRIM(ecp_filename)
  open(newunit=ecpunit, file=TRIM(ecp_filename), status='old', action='read')

  ! Reading an ECP file in NWCHEM format

  end_of_file = .FALSE.
  iline = 0
  line='_____'   ! Five underscores '_____' means 'advance to next line'
  do while(.NOT. end_of_file)
    iline = iline + 1
    if( line(1:5) == '_____' ) then
      read(ecpunit, '(a)', iostat=istat) line
      if( istat == IOSTAT_END ) then
        end_of_file = .TRUE.
        exit
      endif
    endif
    line = ADJUSTL(line)

    ! Remove comments if any
    if( line(1:1) == '#' ) then
      line='_____'
      cycle
    endif

    ! ECP and END should not be interpreted
    if( capitalize(line(1:3)) == 'ECP' .OR. capitalize(line(1:3)) == 'END' ) then
      line='_____'
      cycle
    endif
    i1 = INDEX(line, ' ')

    if( line(1:i1-1) /= TRIM(element) .AND. capitalize(line(1:i1-1)) /= TRIM(ADJUSTL(element)) ) then
      write(stdout, *) 'ECP file should only contain information about element '//TRIM(ADJUSTL(element))
      write(stdout, *) 'While '//line(1:i1-1)//' was found'
      call die('ECP file reading problem')
    endif

    line = ADJUSTL(line(i1+1:))

    i2 = INDEX(line, ' ')
    amc = capitalize(line(1:i2-1))
    if( amc == 'NELEC' ) then
      read(line(i2+1:), '(i10)') ecpi%ncore
      line='_____'
      cycle
    endif
    if(      amc == 'UL'  &
        .OR. amc == 'S'   &
        .OR. amc == 'P'   &
        .OR. amc == 'D'   &
        .OR. amc == 'F'   &
        .OR. amc == 'G'   &
        .OR. amc == 'H'   ) then
      istat = 0
      do while(istat == 0)
        read(ecpunit, '(a)', iostat=istat) line
        if( istat == IOSTAT_END ) then
          end_of_file = .TRUE.
          exit
        endif
        read(line, *, iostat=istat) read_n, read_zeta, read_d

        ! For the time being, only code ECP with no local potential
        if( istat == 0 ) then
          if( amc == 'UL' ) then
            call append_to_list(-1, ecpi%lk)
          else
            call append_to_list(orbital_momentum_number(amc), ecpi%lk)
          endif
          call append_to_list(read_n, ecpi%nk)
          call append_to_list(read_zeta, ecpi%zetak)
          call append_to_list(read_d, ecpi%dk)
        endif
      enddo
    else
      write(stdout, *) capitalize(line(1:i2-1)), line(1:i2-1)
      call die('problem reading ECP file')
    endif


  enddo


  ecpi%necp = SIZE(ecpi%nk)

  close(ecpunit)

end subroutine read_ecp_file


!=========================================================================
! Read psp6 ABINIT format for Haman or TM pseudos generated with FHIPP
subroutine read_psp6_file(ecp_filename, element, ecpi)
  implicit none

  character(*), intent(in)                      :: ecp_filename
  character(len=2), intent(in)                  :: element
  type(effective_core_potential), intent(inout) :: ecpi
  !=====
  integer  :: ecpunit
  integer  :: ir, il, jdum
  integer  :: pspdat, pspcod, pspxc, lmax, lloc, mmax, r2well
  real(dp) :: zatom, zion, al
  character(len=128) :: title
  !=====


  write(stdout, *) 'read PSP6 ECP file:', TRIM(ecp_filename)
  open(newunit=ecpunit, file=TRIM(ecp_filename), status='old', action='read')

  read(ecpunit, *) title
  read(ecpunit, *) zatom, zion, pspdat
  read(ecpunit, *) pspcod, pspxc, lmax, lloc, mmax, r2well
  do jdum=1, 15
    read(ecpunit, *) title
  enddo

  ecpi%ncore = NINT(zatom - zion)
  ecpi%necp  = lmax + 1
  allocate(ecpi%lk(lmax+1))
  do il=0, lmax
    ecpi%lk(il+1) = il
  enddo
  ecpi%lk(lloc+1) = -1
  ecpi%mmax       = mmax

  allocate(ecpi%rad(mmax))
  allocate(ecpi%wfll(mmax, lmax+1))
  allocate(ecpi%vpspll(mmax, lmax+1))
  allocate(ecpi%ekb(lmax+1))

  do il=1, lmax+1
    read(ecpunit, *) title
    do ir=1, mmax
      read(ecpunit, *) jdum, ecpi%rad(ir), ecpi%wfll(ir, il), ecpi%vpspll(ir, il)
    end do
  enddo

  ! Substract the local component from all the other components
  do il=1, lmax+1
    if( il == lloc + 1 ) cycle
    ecpi%vpspll(:, il) = ecpi%vpspll(:, il) - ecpi%vpspll(:, lloc+1)
  enddo

  !
  ! wfll is u_l(r)
  ! \psi_l(r) = u_l(r) / r * Y_lm(r^)
  al = LOG( ecpi%rad(2) / ecpi%rad(1) )

  !debug: Calculate the normalization
  !ecpi%ekb(:) = 0.0_dp
  !do il=1,lmax+1
  !  do ir=1,mmax-1
  !    !  ecpi%ekb(il) = ecpi%ekb(il) + (ecpi%wfll(ir+1,il)**2 + ecpi%wfll(ir,il)**2) &
  !    !                       * ( ecpi%rad(ir+1) - ecpi%rad(ir) ) / 2.0_dp
  !    ! this formula is more accurate (use the logarithmic mesh property)
  !    ecpi%ekb(il) = ecpi%ekb(il)  &
  !               + (ecpi%rad(ir+1) * ecpi%wfll(ir+1,il)**2 + ecpi%rad(ir) * ecpi%wfll(ir,il)**2) &
  !                                     * al / 2.0_dp
  !  end do
  !enddo

  ! Calculate the KB "energy"
  ecpi%ekb(:) = 0.0_dp
  do il=1, lmax+1
    if( il == lloc + 1 ) cycle
    do ir=1, mmax-1
      !  ecpi%ekb(il) = ecpi%ekb(il) &
      !         +  (  ecpi%vpspll(ir+1,il) * ecpi%wfll(ir+1,il)**2 + ecpi%vpspll(ir,il) * ecpi%wfll(ir,il)**2 ) &
      !             * ( ecpi%rad(ir+1) - ecpi%rad(ir) ) / 2.0_dp
      ecpi%ekb(il) = ecpi%ekb(il) &
            +  ( ecpi%rad(ir+1) * ecpi%vpspll(ir+1, il) * ecpi%wfll(ir+1, il)**2 &
                + ecpi%rad(ir)  * ecpi%vpspll(ir, il)   * ecpi%wfll(ir, il)**2 ) &
                        * al / 2.0_dp
    end do
    ecpi%ekb(il) = 1.0_dp / ecpi%ekb(il)
  enddo


  close(ecpunit)


end subroutine read_psp6_file


!=========================================================================
! Read psp8 ABINIT format for Don R Haman ONCVPSP pseudopotential type
subroutine read_psp8_file(ecp_filename, element, ecpi)
  implicit none

  character(*), intent(in)                      :: ecp_filename
  character(len=2), intent(in)                  :: element
  type(effective_core_potential), intent(inout) :: ecpi
  !=====
  integer, parameter :: nder=4      ! number of derivatives in the rhocore
  integer  :: ecpunit
  integer  :: ir, il, jdum, iecp, jecp, iproj
  integer  :: pspdat, pspcod, pspxc, lmax, lloc, mmax, r2well
  real(dp) :: rchrg, fchrg, qchrg
  real(dp) :: zatom, zion, al
  character(len=128) :: title
  integer  :: nproj(5), extension_switch(2)
  !=====


  write(stdout, *) 'read PSP8 ECP file:', TRIM(ecp_filename)
  open(newunit=ecpunit, file=TRIM(ecp_filename), status='old', action='read')

  ! Ne    ONCVPSP-3.3.0  r_core=   1.31204   1.70576
  ! 10.0000      8.0000      171101    zatom,zion,pspd
  ! 8   -1012   1     4   400     0    pspcod,pspxc,lmax,lloc,mmax,r2well
  ! 3.99000000  0.00000000  0.00000000    rchrg fchrg qchrg
  ! 2     2     0     0     0    nproj
  ! 1     1           extension_switch
  ! 0                         2.5609260873145D+00 -4.8346372575205D-01
  ! 1  0.0000000000000D+00  6.8997318969366D-10  3.1387163701435D-10
  ! 2  1.0000000000000D-02  9.2957245369341D-02  1.4251575380487D-02
  ! 3  2.0000000000000D-02  1.8565156365333D-01  2.8580421518622D-02
  ! 4  3.0000000000000D-02  2.7782021681486D-01  4.3063059135132D-02


  read(ecpunit, *) title
  read(ecpunit, *) zatom, zion, pspdat
  read(ecpunit, *) pspcod, pspxc, lmax, lloc, mmax, r2well
  read(ecpunit, *) rchrg, fchrg, qchrg
  read(ecpunit, *) nproj(:)
  read(ecpunit, *) extension_switch(:)
  if( ANY(extension_switch(:) /= 1) ) call die('read_psp8_file: relativistic pseudo file not implemented')

  !if( fchrg > 1.0e-8 ) call die('read_psp8_file: non linear core corrections not implemented')


  ecpi%ncore = NINT(zatom - zion)
  ecpi%necp  = SUM(nproj(:)) + 1   ! +1 corresponds to the local potential
  allocate(ecpi%lk(ecpi%necp))
  iecp = 0
  do il=0, SIZE(nproj)-1
    do iproj=1, nproj(il+1)
      iecp = iecp + 1
      ecpi%lk(iecp) = il
    enddo
  enddo
  ! Last one is the local potential
  ecpi%lk(ecpi%necp) = -1
  ecpi%mmax     = mmax

  allocate(ecpi%rad(mmax))
  allocate(ecpi%vpspll(mmax, ecpi%necp))
  allocate(ecpi%ekb(ecpi%necp))
  if( fchrg > 1.0e-8_dp ) allocate(ecpi%rhocore(mmax, nder+1))

  iecp = 0
  jecp = 0
  do il=0, lmax
    if( nproj(il+1) < 1 ) cycle
    iecp = jecp + 1
    jecp = iecp + nproj(il+1) - 1
    read(ecpunit, *) jdum, ecpi%ekb(iecp:jecp)
    do ir=1, mmax
      read(ecpunit, *) jdum, ecpi%rad(ir), ecpi%vpspll(ir, iecp:jecp)
    end do
  enddo

  ! Read local potential
  read(ecpunit, *) jdum
  do ir=1, mmax
    read(ecpunit, *) jdum, ecpi%rad(ir), ecpi%vpspll(ir, ecpi%necp)
  end do
  if( fchrg > 1.0e-8_dp ) then
    do ir=1, mmax
      read(ecpunit, *) jdum, ecpi%rad(ir), ecpi%rhocore(ir, :)
    end do
    ! in psp8 files, 4*pi*rhoc(r) is written and nobody knows why
    ecpi%rhocore(:, :) = ecpi%rhocore(:, :) / ( 4.0_dp * pi )
  endif

  close(ecpunit)

  !
  ! Checks
  ! whether the radial grid is regular else die
  if( ABS( (ecpi%rad(2)-ecpi%rad(1)) - (ecpi%rad(mmax) - ecpi%rad(mmax-1)) ) > 1.0e-5_dp ) then
    write(stdout, '(1x,a,a)') 'Non-regular grid found in ', TRIM(ecp_filename)
    call die('read_psp8_file: psp8 radial grid must be regular in this implementation')
  endif
  ! whether the first radial grid point is zero
  if( ecpi%rad(1) > 1.0e-5_dp ) then
    write(stdout, '(1x,a,a)') 'Non-zero first grid point found in ', TRIM(ecp_filename)
    call die('read_psp8_file: psp8 radial grid must with zero in this implementation')
  endif


end subroutine read_psp8_file


!=========================================================================
subroutine read_gth_file(ecp_filename, element, ecpi)
  implicit none

  character(*), intent(in)                      :: ecp_filename
  character(len=2), intent(in)                  :: element
  type(effective_core_potential), intent(inout) :: ecpi
  !=====
  integer            :: ecpunit
  integer            :: iline, i1, i2, i3, i4, istat, il
  character(len=132) :: line
  real(dp)           :: rtmp
  logical            :: end_of_file
  !=====

  write(stdout, *) 'read GTH ECP file:', TRIM(ecp_filename)
  open(newunit=ecpunit, file=TRIM(ecp_filename), status='old', action='read')

  ! Reading an GTH file in CP2K format

  ! First line is a comment
  read(ecpunit, '(a)', iostat=istat) line

  ! Second line contains 1, 2, or 3 integers that counts the number of valence electrons for s, p, d
  read(ecpunit, '(a)', iostat=istat) line

  i1=0;i2=0;i3=0
  read(line, *, iostat=istat) i1, i2, i3
  if( istat /= 0 ) then
    read(line, *, iostat=istat) i1, i2
    if( istat /= 0 ) then
      read(line, *, iostat=istat) i1
      if( istat /= 0 ) call die('read_gth_file: 2nd line is not compliant with CP2K format')
    endif
  endif

  ecpi%ncore = element_number(element) - i1 - i2 - i3

  read(ecpunit, *, iostat=istat) ecpi%gth_rpploc, ecpi%gth_nloc
  allocate(ecpi%gth_cipp(ecpi%gth_nloc))
  read(ecpunit, *, iostat=istat) i1
  ecpi%gth_nl = i1
  allocate(ecpi%gth_rl(ecpi%gth_nl))
  allocate(ecpi%gth_npl(ecpi%gth_nl))
  allocate(ecpi%gth_hijl(6, ecpi%gth_nl))
  do il=1, ecpi%gth_nl
    read(ecpunit, *, iostat=istat) rtmp, ecpi%gth_npl(il)
    if( ecpi%gth_npl(il) > 1 ) read(ecpunit, *, iostat=istat)
    if( ecpi%gth_npl(il) > 2 ) read(ecpunit, *, iostat=istat)
  enddo
  close(ecpunit)

  open(newunit=ecpunit, file=TRIM(ecp_filename), status='old', action='read')
  read(ecpunit, '(a)', iostat=istat) line
  read(ecpunit, '(a)', iostat=istat) line
  read(ecpunit, *, iostat=istat) ecpi%gth_rpploc, i1, ecpi%gth_cipp(:)
  read(ecpunit, *, iostat=istat) i1
  do il=1, ecpi%gth_nl
    select case(ecpi%gth_npl(il))
    case(1)
      read(ecpunit, *, iostat=istat) ecpi%gth_rl(il), i1, ecpi%gth_hijl(1, il)
    case(2)
      read(ecpunit, *, iostat=istat) ecpi%gth_rl(il), i1, ecpi%gth_hijl(1:2, il)
      read(ecpunit, *, iostat=istat)                    ecpi%gth_hijl(3, il)
    case(3)
      read(ecpunit, *, iostat=istat) ecpi%gth_rl(il), i1, ecpi%gth_hijl(1:3, il)
      read(ecpunit, *, iostat=istat)                    ecpi%gth_hijl(4:5, il)
      read(ecpunit, *, iostat=istat)                    ecpi%gth_hijl(  6, il)
    case default
      call die('read_gth_file: i > 3 in non-local GTH projector is not possible')
    end select
  enddo


  close(ecpunit)


end subroutine read_gth_file


!=========================================================================
end module m_ecp
!=========================================================================
