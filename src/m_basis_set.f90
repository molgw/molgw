!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to deal with the basis set and basis functions
! It works for both wavefunctions basis sets and auxiliary basis sets
!
!=========================================================================
module m_basis_set
 use m_definitions
 use m_warning
 use m_timing
 use m_mpi
 use m_elements
 use m_string_tools, only: orbital_momentum_name
 use m_atoms
 use m_ecp
 use m_gaussian
 use m_cart_to_pure


 type basis_function
   integer                      :: shell_index                ! This basis function belongs to a shell of basis functions
                                                              ! with the same exponents and angular momentum
   integer                      :: index_in_shell             !
   integer                      :: am                         ! Angular momentum number: l=0, 1, 2, 3 ...
   character(len=1)             :: amc                        ! Angular momentum letter: s, p, d, f ...
   integer                      :: nx,ny,nz                   ! Angular momentum for cartesian gaussians
   integer                      :: mm                         ! Angular momentum for pure gaussians
   integer                      :: iatom                      ! Centered on which atom
   real(dp)                     :: x0(3)                      ! Coordinates of the gaussian center
   integer                      :: ngaussian                  ! Number of primitive gausssians
   type(gaussian),allocatable   :: g(:)                       ! The primitive gaussian functions
   real(dp),allocatable         :: coeff(:)                   ! Their mixing coefficients
 end type

 type shell_type
   integer              :: am
   integer              :: ng
   real(dp),allocatable :: alpha(:)
   real(dp),allocatable :: coeff(:)
   real(dp)             :: x0(3)
   integer              :: iatom
   integer              :: istart,iend                        ! index of the shell's basis functions in the final basis set
   integer              :: istart_cart,iend_cart              ! index of the shell's basis functions in the cartesian basis set
 end type shell_type


 !
 ! A basis set is a list of basis functions
 type basis_set
   !
   ! The list
   integer                          :: ammax           ! Maximum angular momentum contained in the basis set
   integer                          :: nbf             ! Number of basis functions in the basis set
   integer                          :: nbf_cart        ! Number of underlying Cartesian functions in the basis set
   integer                          :: nshell          ! Number of shells in the basis sets
                                                       ! A shell is a group of basis functions sharing:
                                                       ! the same center,
                                                       ! the same exponents,
                                                       ! the same mixing coefficients
                                                       ! and the same angular momentum
   character(len=4)                 :: gaussian_type   ! CART or PURE
   type(basis_function),allocatable :: bfc(:)          ! Cartesian basis function
   type(basis_function),allocatable :: bff(:)          ! Final basis function (can be Cartesian or Pure)
   type(shell_type),allocatable     :: shell(:)

 end type basis_set

contains


!=========================================================================
subroutine init_basis_set(basis_path,basis_name,ecp_basis_name,gaussian_type,basis)
 implicit none

 character(len=4),intent(in) :: gaussian_type
 character(len=*),intent(in) :: basis_path
 character(len=100),intent(in) :: basis_name(natom_basis)
 character(len=100),intent(in) :: ecp_basis_name(natom_basis)
 type(basis_set),intent(out)   :: basis
!=====
 character(len=200)            :: basis_filename
 integer                       :: ibf,jbf,ng,ig
 integer                       :: ishell,ishell_file
 integer                       :: jbf_cart
 real(dp),allocatable          :: alpha(:),coeff(:)
 logical                       :: file_exists
 integer                       :: basisfile
 integer                       :: am_read,nshell_file
 logical,parameter             :: normalized=.TRUE.
 integer                       :: iatom
 integer                       :: index_in_shell
 integer                       :: nx,ny,nz,mm
 real(dp)                      :: x0(3)
!=====

 basis%nbf           = 0
 basis%nbf_cart      = 0
 basis%nshell        = 0
 basis%gaussian_type = gaussian_type

 if( TRIM(basis_name(1)) == 'none' ) return

 !
 ! LOOP OVER ATOMS
 !
 do iatom=1,natom_basis
   if( nelement_ecp > 0 ) then
     if( ANY( element_ecp(:) == zbasis(iatom) ) ) then
       basis_filename = ADJUSTL(TRIM(basis_path) // '/' // TRIM(ADJUSTL(element_name(REAL(zbasis(iatom),dp)))) &
                        // '_' // TRIM(ecp_basis_name(iatom)))
     else
       basis_filename = ADJUSTL(TRIM(basis_path) // '/' // TRIM(ADJUSTL(element_name(REAL(zbasis(iatom),dp)))) &
                       // '_' // TRIM(basis_name(iatom)))
     endif
   else
     basis_filename = ADJUSTL(TRIM(basis_path) // '/' // TRIM(ADJUSTL(element_name(REAL(zbasis(iatom),dp)))) &
                     // '_' // TRIM(basis_name(iatom)))
   endif

   inquire(file=TRIM(basis_filename),exist=file_exists)
   if(.NOT.file_exists) then
     write(stdout,'(1x,a,a)') 'Looking for file ',TRIM(basis_filename)
     write(stdout,'(1x,a)')   'Remember the basis directory path is obtained (by priority order) from:'
     write(stdout,'(1x,a)')   '  1. the input variable basis_path'
     write(stdout,'(1x,a)')   '  2. the environment variable MOLGW_BASIS_PATH'
     write(stdout,'(1x,a)')   '  3. the location of the sources'
     call die('init_basis_set: basis set file not found')
   endif

   !
   ! read first to get all the dimensions
   open(newunit=basisfile,file=TRIM(basis_filename),status='old')
   read(basisfile,*) nshell_file
   if(nshell_file<1) call die('ERROR in basis set file')
   do ishell_file=1,nshell_file
     read(basisfile,*) ng,am_read
     if(ng<1) call die('ERROR in basis set file')
     if(am_read==10) call die('Deprecated basis set file with shared exponent SP orbitals. Please split them')
     basis%nbf_cart = basis%nbf_cart + number_basis_function_am('CART'             ,am_read)
     basis%nbf      = basis%nbf      + number_basis_function_am(basis%gaussian_type,am_read)
     basis%nshell   = basis%nshell   + 1
     do ig=1,ng
       read(basisfile,*)
     enddo
   enddo
   close(basisfile)

 enddo

 allocate(basis%bfc(basis%nbf_cart))
 allocate(basis%bff(basis%nbf))
 allocate(basis%shell(basis%nshell))

 jbf         = 0
 jbf_cart    = 0
 ishell      = 0
 do iatom=1,natom_basis

   if( nelement_ecp > 0 ) then
     if( ANY( element_ecp(:) == zbasis(iatom) ) ) then
       basis_filename = ADJUSTL(TRIM(basis_path) // '/' // TRIM(ADJUSTL(element_name(REAL(zbasis(iatom),dp)))) &
                        // '_' // TRIM(ecp_basis_name(iatom)))
     else
       basis_filename = ADJUSTL(TRIM(basis_path) // '/' // TRIM(ADJUSTL(element_name(REAL(zbasis(iatom),dp)))) &
                       // '_' // TRIM(basis_name(iatom)))
     endif
   else
     basis_filename = ADJUSTL(TRIM(basis_path) // '/' // TRIM(ADJUSTL(element_name(REAL(zbasis(iatom),dp)))) &
                      // '_' //TRIM(basis_name(iatom)))
   endif

   open(newunit=basisfile,file=TRIM(basis_filename),status='old')
   read(basisfile,*) nshell_file
   do ishell_file=1,nshell_file
     read(basisfile,*) ng,am_read
     allocate(alpha(ng),coeff(ng))

     do ig=1,ng
       read(basisfile,*) alpha(ig),coeff(ig)
     enddo

     x0(:) = xbasis(:,iatom)

     !
     ! Shell setup
     !
     ishell = ishell + 1

     basis%shell(ishell)%am      = am_read
     basis%shell(ishell)%iatom   = iatom
     basis%shell(ishell)%x0(:)   = x0(:)
     basis%shell(ishell)%ng      = ng
     allocate(basis%shell(ishell)%alpha(ng))
     allocate(basis%shell(ishell)%coeff(ng))
     basis%shell(ishell)%alpha(:)    = alpha(:)
     basis%shell(ishell)%istart      = jbf + 1
     basis%shell(ishell)%iend        = jbf + number_basis_function_am(gaussian_type,am_read)
     basis%shell(ishell)%istart_cart = jbf_cart + 1
     basis%shell(ishell)%iend_cart   = jbf_cart + number_basis_function_am('CART',am_read)
     ! shell%coeff(:) is setup just after the basis functions


     !
     ! Basis function setup
     !

     !
     ! Ordering of Libint as explained in Kenny et al. J. Comput Chem. 29, 562 (2008).
     !
     nx = am_read
     ny = 0
     nz = 0
     index_in_shell = 0
     do
       ! Add the new basis function
       jbf_cart = jbf_cart + 1
       index_in_shell = index_in_shell + 1
       call init_basis_function(normalized,ng,nx,ny,nz,iatom,x0,alpha,coeff,ishell,index_in_shell,basis%bfc(jbf_cart))
       if(basis%gaussian_type == 'CART') then
         jbf = jbf + 1
         call init_basis_function(normalized,ng,nx,ny,nz,iatom,x0,alpha,coeff,ishell,index_in_shell,basis%bff(jbf))
       endif

       ! Break the loop when nz is equal to l
       if( nz == am_read ) exit

       if( nz < am_read - nx ) then
         ny = ny - 1
         nz = nz + 1
       else
         nx = nx - 1
         ny = am_read - nx
         nz = 0
       endif

     enddo

     index_in_shell = 0
     if(basis%gaussian_type == 'PURE') then
       do mm=-am_read,am_read
         jbf = jbf + 1
         index_in_shell = index_in_shell + 1
         call init_basis_function_pure(normalized,ng,am_read,mm,iatom,x0,alpha,coeff,ishell,index_in_shell,basis%bff(jbf))
       enddo
     endif

     !
     ! Include here the normalization part that does not depend on (nx,ny,nz)
     basis%shell(ishell)%coeff(:) = basis%bfc(jbf_cart-number_basis_function_am(gaussian_type,am_read)+1)%coeff(:) &
               * ( 2.0_dp / pi )**0.75_dp * 2.0_dp**am_read * alpha(:)**( 0.25_dp * ( 2.0_dp*am_read + 3.0_dp ) )

     deallocate(alpha,coeff)
   enddo
   close(basisfile)

 !
 ! END OF THE LOOP OVER ATOMS
 enddo

 ! Find the maximum angular momentum employed in the basis set
 basis%ammax = MAXVAL(basis%bfc(:)%am)

 call echo_basis_summary(basis)


end subroutine init_basis_set


!=========================================================================
subroutine echo_basis_summary(basis)
 implicit none

 type(basis_set),intent(in)    :: basis
!=====
 integer :: ibf
!=====

 write(stdout,'(/,a50,i8)') 'Total number of basis functions:',basis%nbf
 if( basis%gaussian_type == 'PURE' ) then
   write(stdout,'(a50,i8)') 'Total number of cart. functions:',basis%nbf_cart
 endif
 write(stdout,'(a50,i8)') 'Number of shells:',basis%nshell

 write(stdout,'(a50,i8)') 'Maximum angular momentum in the basis set:',basis%ammax
 write(stdout,'(a50,a8)') '                                          ',orbital_momentum_name(basis%ammax)

 if( basis%ammax > MOLGW_LMAX ) then
   write(stdout,*) 'Maximum angular momentum: ',basis%ammax
   write(stdout,*) 'while this compilation of LIBINT only supports: ',MOLGW_LMAX
   call die('echo_basis_summary: Angular momentum is too high')
 endif

 !
 ! Output all the basis functions for debugging
 if( .FALSE. ) then
   do ibf=1,basis%nbf_cart
     write(stdout,*) ' Cartesian function number',ibf
     call print_basis_function(basis%bfc(ibf))
   enddo
 endif

 write(stdout,'(a,/)') ' Basis set is fit and ready'

end subroutine echo_basis_summary


!=========================================================================
!Build a system specific auxiliary basis based the 'Auto' recipe in Gaussian
subroutine init_auxil_basis_set_auto(auxil_basis_name,basis,gaussian_type,auto_auxil_fsam,auto_auxil_lmaxinc,auxil_basis)
 implicit none

 character(len=100),intent(in) :: auxil_basis_name(natom_basis)
 character(len=4),intent(in)   :: gaussian_type
 type(basis_set),intent(in)    :: basis
 real(dp),intent(in)           :: auto_auxil_fsam
 integer,intent(in)            :: auto_auxil_lmaxinc
 type(basis_set),intent(out)   :: auxil_basis
!=====
 logical :: new_primitive,pauto
 integer :: lmax_obs,lmax_abs,am,lval
 integer :: jbf,jbf_cart
 integer :: icenter,ncenter
 integer :: iprim,jprim
 integer :: ishell,jshell,nshell_max
 integer :: icandidate,ncandidate,index_exp,nprim
 integer :: index_current,am_current,am_selected
 integer :: itrial,jtrial,ntrial
 integer :: zcenter
 real(dp) :: exponent_current,exponent_selected
 real(dp),allocatable :: exponent_candidate(:),exponent_trial(:)
 integer,allocatable  :: am_candidate(:),am_trial(:)
 logical,allocatable  :: remaining_candidate(:)
 type(shell_type),allocatable :: shell_tmp(:)
 integer :: index_in_shell,nx,ny,nz,mm
 logical,parameter :: normalized=.TRUE.
!=====

 call start_clock(timing_auto_auxil)

 pauto = TRIM( capitalize(auxil_basis_name(1)) ) == 'PAUTO'

 write(stdout,'(/,1x,a)')      'Auxiliary basis is constructed automatically'
 if( pauto ) then
   write(stdout,'(1x,a)')      '  Method: PAuto'
 else
   write(stdout,'(1x,a)')      '  Method: Auto'
 endif
 write(stdout,'(1x,a)')        '  recipe in Yang, Rendell, Frisch, J. Chem. Phys. 127, 074102 (2007)'
 write(stdout,'(1x,a25,f8.3)') 'Parameter    f_sam: ',auto_auxil_fsam
 write(stdout,'(1x,a25,i3)')   'Parameter l_MAXINC: ',auto_auxil_lmaxinc


 ncenter = MAXVAL( basis%shell(:)%iatom )
 !
 ! Evaluate the maximum number of shells possible
 nshell_max = 0
 do icenter=1,ncenter
   nprim = 0
   lmax_obs = 0
   do ishell=1,basis%nshell
     if( basis%shell(ishell)%iatom == icenter ) then
       nprim = nprim + basis%shell(ishell)%ng
       lmax_obs = MAX( lmax_obs , basis%shell(ishell)%am )
     endif
   enddo
   if( pauto ) then
     nshell_max = nshell_max + nprim**2 * get_lmax_abs(lmax_obs,zbasis(icenter))
   else
     nshell_max = nshell_max + nprim * get_lmax_abs(lmax_obs,zbasis(icenter))
   endif
 enddo


 ! Allocate a temporary storage for the shells
 allocate(shell_tmp(nshell_max))
 auxil_basis%nshell = 0

 write(stdout,'(/,1x,a)')       '========================================'
 do icenter=1,ncenter
   write(stdout,'(1x,a,i5,a,i5)') 'Center:',icenter,' / ',ncenter


   am_current = 0
   zcenter = zbasis(icenter)

   ncandidate = SUM( basis%shell(:)%ng , MASK=(basis%shell(:)%iatom == icenter) )
   if( pauto ) ncandidate = ncandidate**2
   allocate(remaining_candidate(ncandidate))
   allocate(exponent_candidate(ncandidate),exponent_trial(ncandidate))
   allocate(am_candidate(ncandidate),am_trial(ncandidate))

   !
   !  Build the candidate set
   !
   remaining_candidate(:) = .TRUE.
   index_exp = 0
   lmax_obs = 0
   do ishell=1,basis%nshell
     if( basis%shell(ishell)%iatom /= icenter ) cycle
     lmax_obs = MAX( lmax_obs , basis%shell(ishell)%am )

     if( pauto ) then
       ! PAuto recipe:
       ! All the exponents sums are considered
       ! All the angular momenta sums are considered
       do jshell=1,basis%nshell
         if( basis%shell(jshell)%iatom /= icenter ) cycle
         nprim = basis%shell(ishell)%ng * basis%shell(jshell)%ng
         do iprim=1,basis%shell(ishell)%ng
           do jprim=1,basis%shell(jshell)%ng
             index_exp = index_exp + 1
             exponent_candidate(index_exp) = basis%shell(ishell)%alpha(iprim) + basis%shell(jshell)%alpha(jprim)
             am_candidate(index_exp)       = basis%shell(ishell)%am + basis%shell(jshell)%am
           enddo
         enddo
       enddo
     else
       ! Auto recipe:
       ! All the exponents are doubled
       ! All the angular momenta are doubled
       nprim = basis%shell(ishell)%ng
       exponent_candidate(index_exp+1:index_exp+nprim) = basis%shell(ishell)%alpha(:) * 2.0_dp
       am_candidate(index_exp+1:index_exp+nprim)       = basis%shell(ishell)%am * 2
       index_exp = index_exp + nprim
     endif

   enddo

   lmax_abs = get_lmax_abs(lmax_obs,zcenter)

   write(stdout,'(1x,a,i4)') 'Maximum angular momentum for this center: ',lmax_abs


   do while( ANY(remaining_candidate) )
     ntrial = 1
     index_current    = MAXLOC( exponent_candidate(:) , MASK=remaining_candidate(:) , DIM=1 )
     exponent_current = exponent_candidate(index_current)
     exponent_trial(1) = exponent_candidate(index_current)
     am_trial(1) = am_candidate(index_current)
     remaining_candidate(index_current) = .FALSE.

     !
     ! Find the trial set
     ! Duplicate primitives will be removed later
     do icandidate=1,ncandidate
       if( remaining_candidate(icandidate) .AND. exponent_current / exponent_candidate(icandidate) < auto_auxil_fsam ) then
         ntrial = ntrial + 1
         exponent_trial(ntrial) = exponent_candidate(icandidate)
         am_trial(ntrial)       = am_candidate(icandidate)
         remaining_candidate(icandidate) = .FALSE.
       endif
     enddo

     !
     ! Find the unique primitives in the trial set
     !exponent_selected = PRODUCT( exponent_trial(1:ntrial) )**(1.0_dp/REAL(ntrial,dp))
     exponent_selected = 1.0_dp
     nprim = 0
     do itrial=1,ntrial
       new_primitive = .TRUE.
       ! Check if any of the previous trial primitive has the same angular momentum and exponent
       do jtrial=1,itrial-1
         if( am_trial(itrial) == am_trial(jtrial) &
            .AND. ABS( exponent_trial(itrial) - exponent_trial(jtrial) ) < 1.0e-6_dp ) new_primitive = .FALSE.
       enddo
       if( new_primitive ) then
         nprim = nprim + 1
         exponent_selected = exponent_selected * exponent_trial(itrial)
       endif
     enddo
     exponent_selected = (exponent_selected)**(1.0_dp/REAL(nprim,dp))

     am_selected = MAX( MAXVAL(am_trial(1:ntrial)) , am_current )

     am_selected = MIN( am_selected , lmax_abs )    ! Make sure the selected am is not larger than the allowed max
     am_current  = am_selected

     ! Save selected shell information
     do am=0,am_selected
       write(stdout,'(5x,a,es16.6,1x,i2)') '* auxiliary function (exponent, angular momentum):',exponent_selected,am
       auxil_basis%nshell = auxil_basis%nshell + 1
       shell_tmp(auxil_basis%nshell)%am    = am
       shell_tmp(auxil_basis%nshell)%ng    = 1
       allocate(shell_tmp(auxil_basis%nshell)%alpha(1))
       allocate(shell_tmp(auxil_basis%nshell)%coeff(1))
       shell_tmp(auxil_basis%nshell)%alpha(1) = exponent_selected
       shell_tmp(auxil_basis%nshell)%coeff(1) = 1.0_dp
       shell_tmp(auxil_basis%nshell)%iatom = icenter
       shell_tmp(auxil_basis%nshell)%x0(:) = xbasis(:,icenter)
     enddo

   enddo


   deallocate(am_candidate,am_trial)
   deallocate(remaining_candidate)
   deallocate(exponent_candidate,exponent_trial)
 enddo

 write(stdout,'(1x,a,/)')       '========================================'

 !
 ! Fill up the auxil_basis object
 !

 auxil_basis%gaussian_type = gaussian_type
 auxil_basis%ammax         = MAXVAL(shell_tmp(1:auxil_basis%nshell)%am)
 allocate(auxil_basis%shell(auxil_basis%nshell))
 auxil_basis%nbf           = 0
 auxil_basis%nbf_cart      = 0
 do ishell=1,auxil_basis%nshell
   am_current = shell_tmp(ishell)%am
   auxil_basis%nbf      = auxil_basis%nbf      + number_basis_function_am(gaussian_type,am_current)
   auxil_basis%nbf_cart = auxil_basis%nbf_cart + number_basis_function_am('CART'       ,am_current)
 enddo
 allocate(auxil_basis%bfc(auxil_basis%nbf_cart))
 allocate(auxil_basis%bff(auxil_basis%nbf))

 jbf = 0
 jbf_cart = 0
 do ishell=1,auxil_basis%nshell
   am_current = shell_tmp(ishell)%am

   auxil_basis%shell(ishell)%am     = shell_tmp(ishell)%am
   auxil_basis%shell(ishell)%iatom  = shell_tmp(ishell)%iatom
   auxil_basis%shell(ishell)%x0(:)  = shell_tmp(ishell)%x0(:)
   auxil_basis%shell(ishell)%ng     = shell_tmp(ishell)%ng
   allocate(auxil_basis%shell(ishell)%alpha(1))
   allocate(auxil_basis%shell(ishell)%coeff(1))
   auxil_basis%shell(ishell)%alpha(:)    = shell_tmp(ishell)%alpha(:)
   auxil_basis%shell(ishell)%istart      = jbf + 1
   auxil_basis%shell(ishell)%iend        = jbf + number_basis_function_am(gaussian_type,am_current)
   auxil_basis%shell(ishell)%istart_cart = jbf_cart + 1
   auxil_basis%shell(ishell)%iend_cart   = jbf_cart + number_basis_function_am('CART',am_current)

   !
   ! Basis function setup
   !

   !
   ! Ordering of Libint as explained in Kenny et al. J. Comput Chem. 29, 562 (2008).
   !
   nx = am_current
   ny = 0
   nz = 0
   index_in_shell = 0
   do
     ! Add the new basis function
     jbf_cart = jbf_cart + 1
     index_in_shell = index_in_shell + 1
     call init_basis_function(normalized,1,nx,ny,nz,shell_tmp(ishell)%iatom,shell_tmp(ishell)%x0,shell_tmp(ishell)%alpha,&
                              shell_tmp(ishell)%coeff,ishell,index_in_shell,auxil_basis%bfc(jbf_cart))
     if(auxil_basis%gaussian_type == 'CART') then
       jbf = jbf + 1
       call init_basis_function(normalized,1,nx,ny,nz,shell_tmp(ishell)%iatom,shell_tmp(ishell)%x0,shell_tmp(ishell)%alpha,&
                                shell_tmp(ishell)%coeff,ishell,index_in_shell,auxil_basis%bff(jbf))
     endif

     ! Break the loop when nz is equal to l
     if( nz == am_current ) exit

     if( nz < am_current - nx ) then
       ny = ny - 1
       nz = nz + 1
     else
       nx = nx - 1
       ny = am_current - nx
       nz = 0
     endif

   enddo

   index_in_shell = 0
   if(auxil_basis%gaussian_type == 'PURE') then
     do mm=-am_current,am_current
       jbf = jbf + 1
       index_in_shell = index_in_shell + 1
       call init_basis_function_pure(normalized,1,am_current,mm,shell_tmp(ishell)%iatom,shell_tmp(ishell)%x0, &
                                shell_tmp(ishell)%alpha,shell_tmp(ishell)%coeff,ishell,index_in_shell,auxil_basis%bff(jbf))
     enddo
   endif

   auxil_basis%shell(ishell)%coeff(1) = ( 2.0_dp / pi )**0.75_dp * 2.0_dp**am_current &
                                                   * shell_tmp(ishell)%alpha(1)**( 0.25_dp * ( 2.0_dp*am_current + 3.0_dp ) )

 enddo


 call echo_basis_summary(auxil_basis)

 deallocate(shell_tmp)

 call stop_clock(timing_auto_auxil)


contains

function get_lmax_abs(lmax_obs,zcenter) result(lmax_abs)
 integer,intent(in) :: zcenter,lmax_obs
 integer            :: lmax_abs
!======
 integer :: lval
!======
 if( zcenter <= 2 ) then
   lval = 0
 else if( zcenter <= 18) then
   lval = 1
 else if( zcenter <= 54) then
   lval = 2
 else
   lval = 3
 endif
 lmax_abs = MAX( lmax_obs + auto_auxil_lmaxinc , 2 * lval )

end function get_lmax_abs

end subroutine init_auxil_basis_set_auto


!=========================================================================
subroutine destroy_basis_set(basis)
 implicit none

 type(basis_set),intent(inout) :: basis
!=====
 integer :: ishell
!=====

! do ibf=1,basis%nbf_cart
!   call destroy_basis_function(basis%bfc(ibf))
! enddo
! do ibf=1,basis%nbf
!   call destroy_basis_function(basis%bff(ibf))
! enddo
 deallocate(basis%bfc)
 deallocate(basis%bff)
 do ishell=1,basis%nshell
   if(ALLOCATED(basis%shell(ishell)%alpha)) deallocate( basis%shell(ishell)%alpha )
   if(ALLOCATED(basis%shell(ishell)%coeff)) deallocate( basis%shell(ishell)%coeff )
 enddo
 deallocate(basis%shell)

end subroutine destroy_basis_set


!=========================================================================
function compare_basis_set(basis1,basis2) result(same_basis_set)
 implicit none

 logical                    :: same_basis_set
 type(basis_set),intent(in) :: basis1,basis2
!=====
 integer :: ibf
!=====

 same_basis_set = .TRUE.

 if( basis1%ammax         /= basis2%ammax         )  same_basis_set = .FALSE.
 if( basis1%nbf           /= basis2%nbf           )  same_basis_set = .FALSE.
 if( basis1%nbf_cart      /= basis2%nbf_cart      )  same_basis_set = .FALSE.
 if( basis1%nshell        /= basis2%nshell        )  same_basis_set = .FALSE.
 if( basis1%gaussian_type /= basis2%gaussian_type )  same_basis_set = .FALSE.

 ! If the basis sets already differs, then exit immediately
 if( .NOT. same_basis_set ) return

 do ibf=1,basis1%nbf
   same_basis_set = same_basis_set .AND. compare_basis_function(basis1%bfc(ibf),basis2%bfc(ibf))
 enddo


end function compare_basis_set


!=========================================================================
function compare_basis_function(bf1,bf2) result(same_basis_function)
 implicit none

 logical                         :: same_basis_function
 type(basis_function),intent(in) :: bf1,bf2
!=====
 integer                         :: ig
!=====

 same_basis_function = .TRUE.

! DO NOT compare the following commented fields. Not really necessary...
! bf1%shell_index
! bf1%amc
 if( bf1%am            /= bf2%am                        ) same_basis_function = .FALSE.
 if( bf1%nx            /= bf2%nx                        ) same_basis_function = .FALSE.
 if( bf1%ny            /= bf2%ny                        ) same_basis_function = .FALSE.
 if( bf1%nz            /= bf2%nz                        ) same_basis_function = .FALSE.
 if( bf1%iatom         /= bf2%iatom                     ) same_basis_function = .FALSE.
 if( ANY(ABS(bf1%x0(:) - bf2%x0(:)) > 1.0e-5_dp )       ) same_basis_function = .FALSE.
 if( bf1%ngaussian     /= bf2%ngaussian                 ) same_basis_function = .FALSE.

 ! If the basis functions already differs, then exit immediately
 if( .NOT. same_basis_function ) return

 do ig=1,bf1%ngaussian
   same_basis_function = same_basis_function .AND. compare_gaussian(bf1%g(ig),bf2%g(ig))
 enddo
 if( ANY(ABS(bf1%coeff(:) - bf2%coeff(:)) > 1.0e-5_dp ) ) same_basis_function = .FALSE.


end function compare_basis_function


!=========================================================================
subroutine write_basis_set(unitfile,basis)
 implicit none

 integer,intent(in)         :: unitfile
 type(basis_set),intent(in) :: basis
!=====
 integer :: ibf,ishell
!=====

 write(unitfile)  basis%ammax
 write(unitfile)  basis%nbf
 write(unitfile)  basis%nbf_cart
 write(unitfile)  basis%nshell
 write(unitfile)  basis%gaussian_type
 do ibf=1,basis%nbf_cart
   call write_basis_function(unitfile,basis%bfc(ibf))
 enddo
 do ishell=1,basis%nshell
   call write_basis_shell(unitfile,basis%shell(ishell))
 enddo


end subroutine write_basis_set


!=========================================================================
subroutine read_basis_set(unitfile,basis)
 implicit none

 integer,intent(in)          :: unitfile
 type(basis_set),intent(out) :: basis
!=====
 integer :: ibf,ishell
!=====

 read(unitfile)  basis%ammax
 read(unitfile)  basis%nbf
 read(unitfile)  basis%nbf_cart
 read(unitfile)  basis%nshell
 read(unitfile)  basis%gaussian_type
 allocate(basis%bfc(basis%nbf_cart))
 do ibf=1,basis%nbf_cart
   call read_basis_function(unitfile,basis%bfc(ibf))
 enddo
 allocate(basis%shell(basis%nshell))
 do ishell=1,basis%nshell
   call read_basis_shell(unitfile,basis%shell(ishell))
 enddo


end subroutine read_basis_set


!=========================================================================
subroutine write_basis_function(unitfile,bf)
 implicit none

 integer,intent(in)              :: unitfile
 type(basis_function),intent(in) :: bf
!=====
!=====

 write(unitfile)  bf%shell_index
 write(unitfile)  bf%am
 write(unitfile)  bf%amc
 write(unitfile)  bf%nx
 write(unitfile)  bf%ny
 write(unitfile)  bf%nz
 write(unitfile)  bf%iatom
 write(unitfile)  bf%x0(:)
 write(unitfile)  bf%ngaussian
 write(unitfile)  bf%g(:)
 write(unitfile)  bf%coeff(:)


end subroutine write_basis_function


!=========================================================================
subroutine read_basis_function(unitfile,bf)
 implicit none

 integer,intent(in)               :: unitfile
 type(basis_function),intent(out) :: bf
!=====
!=====

 read(unitfile)  bf%shell_index
 read(unitfile)  bf%am
 read(unitfile)  bf%amc
 read(unitfile)  bf%nx
 read(unitfile)  bf%ny
 read(unitfile)  bf%nz
 read(unitfile)  bf%iatom
 read(unitfile)  bf%x0(:)
 read(unitfile)  bf%ngaussian
 allocate(bf%g(bf%ngaussian))
 read(unitfile)  bf%g(:)
 allocate(bf%coeff(bf%ngaussian))
 read(unitfile)  bf%coeff(:)

end subroutine read_basis_function


!=========================================================================
subroutine write_basis_shell(unitfile,shell)
 implicit none

 integer,intent(in)          :: unitfile
 type(shell_type),intent(in) :: shell
!=====
!=====

 write(unitfile) shell%am
 write(unitfile) shell%ng
 write(unitfile) shell%alpha(:)
 write(unitfile) shell%coeff(:)
 write(unitfile) shell%x0(:)
 write(unitfile) shell%iatom
 write(unitfile) shell%istart,shell%iend
 write(unitfile) shell%istart_cart,shell%iend_cart

end subroutine write_basis_shell


!=========================================================================
subroutine read_basis_shell(unitfile,shell)
 implicit none

 integer,intent(in)           :: unitfile
 type(shell_type),intent(out) :: shell
!=====
!=====

 read(unitfile) shell%am
 read(unitfile) shell%ng
 allocate(shell%alpha(shell%ng))
 read(unitfile) shell%alpha(:)
 allocate(shell%coeff(shell%ng))
 read(unitfile) shell%coeff(:)
 read(unitfile) shell%x0(:)
 read(unitfile) shell%iatom
 read(unitfile) shell%istart,shell%iend
 read(unitfile) shell%istart_cart,shell%iend_cart

end subroutine read_basis_shell


!=========================================================================
subroutine init_basis_function(normalized,ng,nx,ny,nz,iatom,x0,alpha,coeff,shell_index,index_in_shell,bf)
 implicit none
 logical,intent(in)               :: normalized
 integer,intent(in)               :: ng,nx,ny,nz,shell_index,iatom,index_in_shell
 real(dp),intent(in)              :: x0(3),alpha(ng)
 real(dp),intent(in)              :: coeff(ng)
 type(basis_function),intent(out) :: bf
!=====
 integer                          :: ig
 real(dp)                         :: overlap
!=====

 bf%ngaussian = ng
 allocate(bf%g(bf%ngaussian))
 allocate(bf%coeff(bf%ngaussian))
 bf%nx    = nx
 bf%ny    = ny
 bf%nz    = nz
 bf%am    = nx + ny + nz
 bf%mm    = -100          ! A fake value
 bf%amc   = orbital_momentum_name(bf%am)
 bf%iatom = iatom
 bf%x0(:) = x0(:)
 bf%shell_index    = shell_index
 bf%index_in_shell = index_in_shell

 ! All the gaussians of the contraction have the same orbital momentum
 do ig=1,bf%ngaussian
   call init_gaussian_general(nx,ny,nz,alpha(ig),x0,bf%g(ig))
   bf%coeff(ig) = coeff(ig)
 enddo

 !
 ! check the normalization if requested
 if( normalized ) then
   call overlap_basis_function(bf,bf,overlap)
   if( ABS(overlap-1.0_dp) > 2.0e-5_dp ) then
     bf%coeff(:) = coeff(:) / SQRT( overlap )
!     write(stdout,*) 'normalization is different from 1.0',overlap
!     write(stdout,*) bf%nx,bf%ny,bf%nz
!     write(stdout,*) 'assuming this is a generalized contraction and rescaling coefficients'
!     write(stdout,*) bf%coeff(:)
!     write(stdout,*)
   endif
 endif


end subroutine init_basis_function


!=========================================================================
subroutine init_basis_function_pure(normalized,ng,am,mm,iatom,x0,alpha,coeff,shell_index,index_in_shell,bf)
 implicit none
 logical,intent(in)               :: normalized
 integer,intent(in)               :: ng,am,mm,shell_index,iatom,index_in_shell
 real(dp),intent(in)              :: x0(3),alpha(ng)
 real(dp),intent(in)              :: coeff(ng)
 type(basis_function),intent(out) :: bf
!=====
!=====

 bf%ngaussian = ng
 allocate(bf%g(bf%ngaussian))
 allocate(bf%coeff(bf%ngaussian))
 bf%nx    = -1
 bf%ny    = -1
 bf%nz    = -1
 bf%am    = am
 bf%mm    = mm
 bf%amc   = orbital_momentum_name(bf%am)
 bf%iatom = iatom
 bf%x0(:) = x0(:)
 bf%shell_index = shell_index
 bf%index_in_shell = index_in_shell
 bf%g(:)%alpha = alpha(:)
 bf%coeff(:)   = coeff(:)
 if( normalized ) continue

! Do not need this

!  ! All the gaussians of the contraction have the same orbital momentum
!  do ig=1,bf%ngaussian
!    call init_gaussian_general(nx,ny,nz,alpha(ig),x0,bf%g(ig))
!    bf%coeff(ig) = coeff(ig)
!  enddo
!
!  !
!  ! check the normalization if requested
!  if( normalized ) then
!    call overlap_basis_function(bf,bf,overlap)
!    if( ABS(overlap-1.0_dp) > 2.0e-5_dp ) then
! !     write(stdout,*) 'normalization is different from 1.0',overlap
! !     write(stdout,*) bf%nx,bf%ny,bf%nz
! !     write(stdout,*) 'assuming this is a generalized contraction and rescaling coefficients'
!      bf%coeff(:) = coeff(:) / SQRT( overlap )
!    endif
!  endif


end subroutine init_basis_function_pure


!=========================================================================
subroutine destroy_basis_function(bf)
 implicit none
 type(basis_function),intent(inout) :: bf
!=====

 deallocate(bf%g,bf%coeff)

end subroutine destroy_basis_function


!=========================================================================
subroutine print_basis_function(bf)
 implicit none
 type(basis_function),intent(in) :: bf
!=====
 integer :: ig
!=====

 write(stdout,*)
 write(stdout,*) '======  print out a basis function ======'
 write(stdout,'(a30,2x,1(1x,i3))')           'contraction of N gaussians',bf%ngaussian
 write(stdout,'(a30,5x,a1)')                'orbital momentum',bf%amc
 write(stdout,'(a30,1x,3(f12.6,2x))')        'centered in',bf%x0(:)
 do ig=1,bf%ngaussian
   write(stdout,'(a30,2x,1x,i3,2x,f12.6)')   'coefficient',ig,bf%coeff(ig)
 enddo
 write(stdout,*)
 do ig=1,bf%ngaussian
   call print_gaussian(bf%g(ig))
 enddo
 write(stdout,*) '====== end of basis function ======'
 write(stdout,*)

end subroutine print_basis_function


!=========================================================================
function eval_basis_function(bf,x)
 implicit none
 type(basis_function),intent(in) :: bf
 real(dp),intent(in)             :: x(3)
 real(dp)                        :: eval_basis_function
!=====
 integer                         :: ig
!=====

 eval_basis_function=0.0_dp
 do ig=1,bf%ngaussian
   eval_basis_function = eval_basis_function + eval_gaussian(bf%g(ig),x) * bf%coeff(ig)
 enddo

end function eval_basis_function


!=========================================================================
function eval_basis_function_grad(bf,x)
 implicit none
 type(basis_function),intent(in) :: bf
 real(dp),intent(in)             :: x(3)
 real(dp)                        :: eval_basis_function_grad(3)
!=====
 integer                         :: ig
!=====

 eval_basis_function_grad(:)=0.0_dp
 do ig=1,bf%ngaussian
   eval_basis_function_grad(:) = eval_basis_function_grad(:) + eval_gaussian_grad(bf%g(ig),x) * bf%coeff(ig)
 enddo

end function eval_basis_function_grad


!=========================================================================
function eval_basis_function_lapl(bf,x)
 implicit none
 type(basis_function),intent(in) :: bf
 real(dp),intent(in)             :: x(3)
 real(dp)                        :: eval_basis_function_lapl(3)
!=====
 integer                         :: ig
!=====

 eval_basis_function_lapl(:)=0.0_dp
 do ig=1,bf%ngaussian
   eval_basis_function_lapl(:) = eval_basis_function_lapl(:) + eval_gaussian_lapl(bf%g(ig),x) * bf%coeff(ig)
 enddo

end function eval_basis_function_lapl


!=========================================================================
subroutine overlap_basis_function(bf1,bf2,overlap)
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(out)            :: overlap
!=====
 integer                         :: ig,jg
 real(dp)                        :: overlap_one_gaussian
!=====

 overlap=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call overlap_recurrence(bf1%g(ig),bf2%g(jg),overlap_one_gaussian)
     overlap = overlap + overlap_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


end subroutine overlap_basis_function


!=========================================================================
subroutine overlap_three_basis_function(bf1,bf2,bf3,overlap)
 implicit none
 type(basis_function),intent(in) :: bf1,bf2,bf3
 real(dp),intent(out)            :: overlap
!=====
 type(basis_function),allocatable :: bf12(:)
 integer                          :: ibf
 real(dp)                         :: overlap_tmp
!=====

 !
 ! first multiply the two first basis functions
 call basis_function_prod(bf1,bf2,bf12)

 !
 ! then overlap the product and the third basis function
 overlap = 0.0_dp
 do ibf=1,SIZE(bf12(:))

   call overlap_basis_function(bf12(ibf),bf3,overlap_tmp)
   overlap = overlap + overlap_tmp
   call destroy_basis_function(bf12(ibf))

 enddo
 !
 ! don't forget to destroy it, else memory is leaking
 deallocate(bf12)


end subroutine overlap_three_basis_function


!=========================================================================
subroutine kinetic_basis_function(bf1,bf2,kinetic)
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(out)            :: kinetic
!=====
 integer                         :: ig,jg
 real(dp)                        :: kinetic_one_gaussian
!=====

 kinetic=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call kinetic_recurrence(bf1%g(ig),bf2%g(jg),kinetic_one_gaussian)
     kinetic = kinetic + kinetic_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


end subroutine kinetic_basis_function


!=========================================================================
subroutine nucleus_basis_function(bf1,bf2,zatom,x,nucleus_pot)
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(in)             :: zatom,x(3)
 real(dp),intent(out)            :: nucleus_pot
!=====
 integer                         :: ig,jg
 real(dp)                        :: nucleus_pot_one_gaussian
!=====

 nucleus_pot=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call nucleus_recurrence(zatom,x,bf1%g(ig),bf2%g(jg),nucleus_pot_one_gaussian)
     nucleus_pot = nucleus_pot + nucleus_pot_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


end subroutine nucleus_basis_function


!=========================================================================
subroutine basis_function_prod(bf1,bf2,bfprod)
 implicit none
 type(basis_function),intent(in)              :: bf1,bf2
 type(basis_function),allocatable,intent(out) :: bfprod(:)
!=====
 integer           :: ig,jg
 logical,parameter :: unnormalized=.FALSE.
 integer           :: fake_shell=1
 integer           :: fake_index=1
 integer           :: nbfprod,ibf
 integer           :: ix,iy,iz
 real(dp)          :: coeff_xyzg(1),c_x,c_y,c_z
 real(dp)          :: coeff
 real(dp)          :: alpha(1)
 real(dp)          :: x0(3)
 real(dp)          :: exp_fact
!=====

 nbfprod = ( 1 + bf1%nx + bf2%nx ) * ( 1 + bf1%ny + bf2%ny ) * ( 1 + bf1%nz + bf2%nz) * bf1%ngaussian * bf2%ngaussian

 allocate(bfprod(nbfprod))

 ibf = 0
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian

     alpha(1) = bf1%g(ig)%alpha + bf2%g(jg)%alpha
     coeff = bf1%coeff(ig) * bf2%coeff(jg) *  bf1%g(ig)%norm_factor * bf2%g(jg)%norm_factor
     x0(:)  = ( bf1%g(ig)%alpha * bf1%x0(:) + bf2%g(jg)%alpha * bf2%x0(:) ) / alpha(1)
     exp_fact = EXP( -bf1%g(ig)%alpha * bf2%g(jg)%alpha / alpha(1) * SUM( (bf1%x0(:)-bf2%x0(:))**2 ) )

     do ix=0,bf1%nx + bf2%nx
       c_x = c_1d(ix,bf1%nx,bf2%nx,bf1%x0(1),bf2%x0(1),x0(1))
       do iy=0,bf1%ny + bf2%ny
         c_y = c_1d(iy,bf1%ny,bf2%ny,bf1%x0(2),bf2%x0(2),x0(2))
         do iz=0,bf1%nz + bf2%nz
           c_z = c_1d(iz,bf1%nz,bf2%nz,bf1%x0(3),bf2%x0(3),x0(3))
           ibf = ibf + 1

           coeff_xyzg(1) = coeff * c_x * c_y * c_z
           call init_basis_function(unnormalized,1,ix,iy,iz,0,x0,alpha,coeff_xyzg,fake_shell,fake_index,bfprod(ibf))

           ! override the normalization
           ! the product gaussians are UNnormalized
           ! consistently with the ERI basis
           bfprod(ibf)%g(1)%norm_factor = 1.0_dp

         enddo
       enddo
     enddo

   enddo
 enddo


contains

function c_1d(ix,nx1,nx2,xi,xj,xk)
 implicit none
 integer,intent(in)  :: ix,nx1,nx2
 real(dp),intent(in) :: xi,xj,xk
 real(dp)            :: c_1d
!=====
 integer :: ii,jj
!=====

 c_1d = 0.0_dp
 do ii=0,nx1
   do jj=MAX(ix-ii,0),nx2
     c_1d = c_1d + cnk(nx1,ii) * cnk(nx2,jj) * ( xk - xi )**(nx1-ii) * ( xk - xj )**(nx2-jj)
   enddo
 enddo

end function c_1d

end subroutine basis_function_prod


!=========================================================================
subroutine basis_function_dipole(bf1,bf2,dipole)
 implicit none
 type(basis_function),intent(in)  :: bf1,bf2
 real(dp),intent(out)             :: dipole(3)
!=====
 type(basis_function)             :: bftmp
 real(dp)                         :: dipole_tmp
 logical,parameter                :: normalized=.FALSE.
 integer                          :: fake_shell=1
 integer                          :: fake_index=1
 real(dp),allocatable             :: bf2_alpha(:)
!=====

 !
 ! Calculate < phi_1 | r | phi_2 >
 ! using r = ( r - B ) + B
 !
 allocate(bf2_alpha(bf2%ngaussian))
 bf2_alpha(:) = bf2%g(:)%alpha

 ! first set up | (x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(1) = dipole_tmp
 ! first set up | Bx phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(1) = dipole(1) + dipole_tmp * bf2%x0(1)

 ! first set up | (y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(2) = dipole_tmp
 ! first set up | By phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | By phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(2) = dipole(2) + dipole_tmp * bf2%x0(2)

 ! first set up | (z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(3) = dipole_tmp
 ! first set up | Bz phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bz phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(3) = dipole(3) + dipole_tmp * bf2%x0(3)

 deallocate(bf2_alpha)

end subroutine basis_function_dipole


!=========================================================================
subroutine basis_function_quadrupole(bf1,bf2,quad)
 implicit none
 type(basis_function),intent(in)  :: bf1,bf2
 real(dp),intent(out)             :: quad(3,3)
!=====
 type(basis_function)             :: bftmp
 real(dp)                         :: quad_tmp
 logical,parameter                :: normalized=.FALSE.
 integer                          :: fake_shell=1
 integer                          :: fake_index=1
 real(dp),allocatable             :: bf2_alpha(:)
!=====

 !
 ! Calculate < phi_1 | x y | phi_2 >
 ! using x y = ( x - Bx ) * ( y - By) + Bx * ( y - By ) + By * ( x - Bx ) + Bx * By
 !

 allocate(bf2_alpha(bf2%ngaussian))
 bf2_alpha(:) = bf2%g(:)%alpha


 !
 !  terms x * y and y * x
 !

 ! first set up | (x-Bx)*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny+1,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad_tmp

 ! first set up | Bx*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad(1,2) + quad_tmp * bf2%x0(1)

 ! first set up | By*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | By*(x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad(1,2) + quad_tmp * bf2%x0(2)

 ! Overlap < phi1 | Bx*By phi2 >
 call overlap_basis_function(bf1,bf2,quad_tmp)
 quad(1,2) = quad(1,2) + quad_tmp * bf2%x0(1) * bf2%x0(2)

 quad(2,1) = quad(1,2)


 !
 !  terms x * z and z * x
 !

 ! first set up | (x-Bx)*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz+1,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad_tmp

 ! first set up | Bx*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad(1,3) + quad_tmp * bf2%x0(1)

 ! first set up | Bz*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bz*(x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad(1,3) + quad_tmp * bf2%x0(3)

 ! Overlap < phi1 | Bx*Bz phi2 >
 call overlap_basis_function(bf1,bf2,quad_tmp)
 quad(1,3) = quad(1,3) + quad_tmp * bf2%x0(1) * bf2%x0(3)

 quad(3,1) = quad(1,3)


 !
 !  terms y * z and z * y
 !

 ! first set up | (y-By)*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz+1,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By)*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad_tmp

 ! first set up | By*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | By*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad(2,3) + quad_tmp * bf2%x0(2)

 ! first set up | Bz*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bz*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad(2,3) + quad_tmp * bf2%x0(3)

 ! Overlap < phi1 | By*Bz phi2 >
 call overlap_basis_function(bf1,bf2,quad_tmp)
 quad(2,3) = quad(2,3) + quad_tmp * bf2%x0(2) * bf2%x0(3)

 quad(3,2) = quad(2,3)


 !
 !  term x * x
 !

 ! first set up | (x-Bx)*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+2,bf2%ny,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad_tmp

 ! first set up | 2Bx*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2Bx*(x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad(1,1) + quad_tmp * 2.0_dp * bf2%x0(1)

 ! first set up | Bx**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad(1,1) + quad_tmp * bf2%x0(1)**2


 !
 !  term y * y
 !

 ! first set up | (y-By)^2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+2,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad_tmp

 ! first set up | 2B*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2B*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad(2,2) + quad_tmp * 2.0_dp * bf2%x0(2)

 ! first set up | By**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | By**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad(2,2) + quad_tmp * bf2%x0(2)**2


 !
 !  term z * z
 !

 ! first set up | (z-Bz)^2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+2,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (z-Bz)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad_tmp

 ! first set up | 2B*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2B*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad(3,3) + quad_tmp * 2.0_dp * bf2%x0(3)

 ! first set up | Bz**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0, &
                          bf2%x0,bf2_alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bz**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad(3,3) + quad_tmp * bf2%x0(3)**2


 deallocate(bf2_alpha)


end subroutine basis_function_quadrupole


!=========================================================================
subroutine basis_function_gos(bf1,bf2,qvec,gos_bf1bf2)
 implicit none
 type(basis_function),intent(in)  :: bf1,bf2
 real(dp),intent(in)              :: qvec(3)
 complex(dp),intent(out)          :: gos_bf1bf2
!=====
 integer                          :: ig,jg
 complex(dp)                      :: gos_one_gaussian
!=====

 gos_bf1bf2 = 0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call evaluate_gos(bf1%g(ig),bf2%g(jg),qvec,gos_one_gaussian)
     gos_bf1bf2 = gos_bf1bf2 + gos_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo

end subroutine basis_function_gos


!=========================================================================
function basis_function_fourier(bf,qvec) RESULT(fourier_bf)
 implicit none
 type(basis_function),intent(in)  :: bf
 real(dp),intent(in)              :: qvec(3)
 complex(dp)                      :: fourier_bf
!=====
 integer                          :: ig
!=====

 fourier_bf = 0.0_dp
 do ig=1,bf%ngaussian
   fourier_bf = fourier_bf + eval_gaussian_fourier(bf%g(ig),qvec) * bf%coeff(ig)
 enddo

end function basis_function_fourier


!=========================================================================
end module m_basis_set
!=========================================================================
