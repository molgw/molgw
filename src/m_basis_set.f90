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
 use m_mpi
 use m_timing
 use m_elements
 use m_tools, only: diagonalize,invert,double_factorial,orbital_momentum_name
 use m_atoms
 use m_ecp
 use m_gaussian


 type transform
   real(dp),allocatable         :: matrix(:,:)
 end type
 integer,parameter              :: lmax_transform     = 7
 integer,parameter              :: lmax_transform_pure= 5
 type(transform)                :: cart_to_pure(0:lmax_transform)
 type(transform)                :: cart_to_pure_norm(0:lmax_transform)

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
   type(basis_function),allocatable :: bf(:)           ! Cartesian basis function
   type(basis_function),allocatable :: bff(:)          ! Final basis function (can be Cartesian or Pure)

 end type basis_set


contains


!=========================================================================
subroutine init_basis_set(basis_path,basis_name,ecp_basis_name,gaussian_type,basis)
 implicit none
 character(len=4),intent(in)   :: gaussian_type
 character(len=100),intent(in) :: basis_path
 character(len=100),intent(in) :: basis_name(natom_basis)
 character(len=100),intent(in) :: ecp_basis_name(natom_basis)
 type(basis_set),intent(out)   :: basis
!=====
 character(len=100)            :: basis_filename
 integer                       :: ibf,jbf,kbf,ng,ig,shell_index,ibf_file
 integer                       :: jbf_cart
 real(dp),allocatable          :: alpha(:),coeff(:),coeff2(:)
 logical                       :: file_exists
 integer                       :: basisfile
 integer                       :: am_read,nbf_file
 logical,parameter             :: normalized=.TRUE.
 integer                       :: iatom
 integer                       :: index_in_shell
 integer                       :: nx,ny,nz,mm
 real(dp)                      :: x0(3)
!=====

 basis%nbf           = 0
 basis%nbf_cart      = 0
 basis%gaussian_type = gaussian_type
 basis%ammax         = -1

 if(TRIM(basis_name(1))=='none') return

 !
 ! LOOP OVER ATOMS
 ! TODO could be reduced to the type of atoms in the future
 !
 do iatom=1,natom_basis

   if( nelement_ecp > 0 ) then 
     if( ANY( element_ecp(:) == basis_element(iatom) ) ) then
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(ecp_basis_name(iatom)))
     else
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
     endif
   else
     basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
   endif
  
   inquire(file=TRIM(basis_filename),exist=file_exists)
   if(.NOT.file_exists) then
     write(stdout,'(a,a)') ' Looking for file ',TRIM(basis_filename)
     call die('basis set file not found')
   endif
  
   !
   ! read first to get all the dimensions
   open(newunit=basisfile,file=TRIM(basis_filename),status='old')
   read(basisfile,*) nbf_file
   if(nbf_file<1) call die('ERROR in basis set file')
   do ibf_file=1,nbf_file
     read(basisfile,*) ng,am_read
     if(ng<1) call die('ERROR in basis set file')
     basis%nbf_cart = basis%nbf_cart + number_basis_function_am('CART'             ,am_read)
     basis%nbf      = basis%nbf      + number_basis_function_am(basis%gaussian_type,am_read)
     do ig=1,ng
       read(basisfile,*) 
     enddo
   enddo
   close(basisfile)
  
 enddo


 write(stdout,*)
 write(stdout,'(a50,i8)') 'Total number of basis functions:',basis%nbf
 if(basis%gaussian_type=='PURE') then
   write(stdout,'(a50,i8)') 'Total number of cart. functions:',basis%nbf_cart
 endif
 allocate(basis%bf(basis%nbf_cart))
 allocate(basis%bff(basis%nbf))

 jbf         = 0
 jbf_cart    = 0
 shell_index = 0
 do iatom=1,natom_basis

   if( nelement_ecp > 0 ) then 
     if( ANY( element_ecp(:) == basis_element(iatom) ) ) then
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(ecp_basis_name(iatom)))
     else
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
     endif
   else
     basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
   endif
  
   open(newunit=basisfile,file=TRIM(basis_filename),status='old')
   read(basisfile,*) nbf_file
   do ibf_file=1,nbf_file
     read(basisfile,*) ng,am_read
     allocate(alpha(ng),coeff(ng),coeff2(ng))
  
     if(am_read<10) then
       do ig=1,ng
         read(basisfile,*) alpha(ig),coeff(ig)
       enddo
     else
       do ig=1,ng
         read(basisfile,*) alpha(ig),coeff(ig),coeff2(ig)
       enddo
     endif
  
!     ! rescale the gaussian decay rate whenever zatom /= basis_element
!     if( abs( zatom(iatom) - REAL(basis_element(iatom),dp) ) > 1.d-6 ) then
!       alpha(:) = alpha(:) * ( zatom(iatom) / REAL(basis_element(iatom),dp) )**2
!       write(stdout,*) 'rescaling momentum',am_read
!       write(stdout,*) 'smallest rescaled alpha:',MINVAL(alpha(:))
!     endif
  
     x0(:) = x(:,iatom)

     shell_index = shell_index + 1

     select case(am_read)
     case(10) ! SP shared exponent
       ! Cartesian basis functions
       jbf_cart = jbf_cart + 1 ; call init_basis_function(normalized,ng,0,0,0,iatom,x0,alpha,coeff,shell_index,1,basis%bf(jbf_cart))
       shell_index = shell_index + 1
       jbf_cart = jbf_cart + 1 ; call init_basis_function(normalized,ng,1,0,0,iatom,x0,alpha,coeff2,shell_index,1,basis%bf(jbf_cart))
       jbf_cart = jbf_cart + 1 ; call init_basis_function(normalized,ng,0,1,0,iatom,x0,alpha,coeff2,shell_index,2,basis%bf(jbf_cart))
       jbf_cart = jbf_cart + 1 ; call init_basis_function(normalized,ng,0,0,1,iatom,x0,alpha,coeff2,shell_index,3,basis%bf(jbf_cart))

       ! Final basis function
       if(basis%gaussian_type == 'CART') then
         jbf     = jbf     + 1 ; call init_basis_function(normalized,ng,0,0,0,iatom,x0,alpha,coeff,shell_index,1,basis%bff(jbf))
         shell_index = shell_index + 1
         jbf     = jbf     + 1 ; call init_basis_function(normalized,ng,1,0,0,iatom,x0,alpha,coeff2,shell_index,1,basis%bff(jbf))
         jbf     = jbf     + 1 ; call init_basis_function(normalized,ng,0,1,0,iatom,x0,alpha,coeff2,shell_index,2,basis%bff(jbf))
         jbf     = jbf     + 1 ; call init_basis_function(normalized,ng,0,0,1,iatom,x0,alpha,coeff2,shell_index,3,basis%bff(jbf))
       else
         jbf     = jbf     + 1 ; call init_basis_function_pure(normalized,ng,0, 0,iatom,x0,alpha,coeff,shell_index,1,basis%bff(jbf))
         shell_index = shell_index + 1
         jbf     = jbf     + 1 ; call init_basis_function_pure(normalized,ng,1,-1,iatom,x0,alpha,coeff2,shell_index,1,basis%bff(jbf))
         jbf     = jbf     + 1 ; call init_basis_function_pure(normalized,ng,1, 0,iatom,x0,alpha,coeff2,shell_index,2,basis%bff(jbf))
         jbf     = jbf     + 1 ; call init_basis_function_pure(normalized,ng,1, 1,iatom,x0,alpha,coeff2,shell_index,3,basis%bff(jbf))
       endif

     case default

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
         call init_basis_function(normalized,ng,nx,ny,nz,iatom,x0,alpha,coeff,shell_index,index_in_shell,basis%bf(jbf_cart))
         if(basis%gaussian_type == 'CART') then
           jbf = jbf + 1
           call init_basis_function(normalized,ng,nx,ny,nz,iatom,x0,alpha,coeff,shell_index,index_in_shell,basis%bff(jbf))
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
           call init_basis_function_pure(normalized,ng,am_read,mm,iatom,x0,alpha,coeff,shell_index,index_in_shell,basis%bff(jbf))
         enddo
       endif

     end select
  
     deallocate(alpha,coeff,coeff2)
   enddo
   close(basisfile)

 !
 ! END OF THE LOOP OVER ATOMS
 enddo
 

 basis%nshell = shell_index
 write(stdout,'(a50,i8)') 'Number of shells:',basis%nshell

 ! Find the maximum angular momentum employed in the basis set
 basis%ammax=-1
 do ibf=1,basis%nbf
   basis%ammax = MAX(basis%ammax,basis%bf(ibf)%am)
 enddo
 write(stdout,'(a50,i8)') 'Maximum angular momentum in the basis set:',basis%ammax
 write(stdout,'(a50,a8)') '                                          ',orbital_momentum_name(basis%ammax)

 if(basis%ammax > lmax_transform ) then      
   write(stdout,*) 'Maximum angular momentum',basis%ammax
   call die('angular momentum too high')
 endif
 if(basis%ammax > lmax_transform_pure .AND. basis%gaussian_type == 'PURE' ) then      
   call issue_warning('Maximum angular momentum greater than the cart to pure transforms implemented')
 endif

 !
 ! finally output the basis set for debugging
 if( .FALSE. ) then
   do ibf=1,basis%nbf
     write(stdout,*) ' Cartesian function number',ibf
     call print_basis_function(basis%bf(ibf))
   enddo
 endif

 write(stdout,'(a,/)') ' Basis set is ready and fit'

end subroutine init_basis_set


!=========================================================================
subroutine destroy_basis_set(basis)
 implicit none

 type(basis_set),intent(inout) :: basis
!=====
 integer :: ibf
!=====

! do ibf=1,basis%nbf
!   call destroy_basis_function(basis%bf(ibf))
! enddo
 deallocate(basis%bf)

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
   same_basis_set = same_basis_set .AND. compare_basis_function(basis1%bf(ibf),basis2%bf(ibf))
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
 integer :: ibf
!=====

 write(unitfile)  basis%ammax         
 write(unitfile)  basis%nbf           
 write(unitfile)  basis%nbf_cart      
 write(unitfile)  basis%nshell        
 write(unitfile)  basis%gaussian_type
 do ibf=1,basis%nbf_cart
   call write_basis_function(unitfile,basis%bf(ibf))
 enddo
 

end subroutine write_basis_set


!=========================================================================
subroutine read_basis_set(unitfile,basis)
 implicit none

 integer,intent(in)          :: unitfile
 type(basis_set),intent(out) :: basis
!=====
 integer :: ibf
!=====

 read(unitfile)  basis%ammax
 read(unitfile)  basis%nbf
 read(unitfile)  basis%nbf_cart
 read(unitfile)  basis%nshell
 read(unitfile)  basis%gaussian_type
 allocate(basis%bf(basis%nbf_cart))
 do ibf=1,basis%nbf_cart
   call read_basis_function(unitfile,basis%bf(ibf))
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
!     write(stdout,*) 'normalization is different from 1.0',overlap
!     write(stdout,*) bf%nx,bf%ny,bf%nz
!     write(stdout,*) 'assuming this is a generalized contraction and rescaling coefficients'
     bf%coeff(:) = coeff(:) / SQRT( overlap )
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
 integer                          :: ig
 real(dp)                         :: overlap
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
function number_basis_function_am(gaussian_type,am)
 implicit none
 character(len=4),intent(in) :: gaussian_type
 integer,intent(in)          :: am
 integer                     :: number_basis_function_am
!=====

 select case(gaussian_type)
 case('CART')
   select case(am)
   case(0)
     number_basis_function_am = 1
   case(1)
     number_basis_function_am = 3
   case(2)
     number_basis_function_am = 6
   case(3)
     number_basis_function_am = 10
   case(4)
     number_basis_function_am = 15
   case(5)
     number_basis_function_am = 21
   case(6)
     number_basis_function_am = 28
   case(7)
     number_basis_function_am = 36
   case(10) ! stands for SP orbitals
     number_basis_function_am = 4 
   case default
     write(stdout,*) 'am=',am
     call die('number_basis_function_am: not implemented')
   end select
 case('PURE')
   if(am/=10) then
     number_basis_function_am = 2 * am + 1
   else ! stands for SP orbitals
     number_basis_function_am = 4 
   endif
 end select

end function number_basis_function_am


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
 type(basis_function)            :: bf12
 integer                         :: ig,jg
 real(dp)                        :: overlap_one_gaussian
!=====

 if(mod(bf1%nx+bf2%nx+bf3%nx,2)==1) then
   overlap=0.0_dp
   return
 endif
 if(mod(bf1%ny+bf2%ny+bf3%ny,2)==1) then
   overlap=0.0_dp
   return
 endif
 if(mod(bf1%nz+bf2%nz+bf3%nz,2)==1) then
   overlap=0.0_dp
   return
 endif
 !
 ! first multiply the two first basis functions
 call basis_function_prod(bf1,bf2,bf12)

 !
 ! then overlap the product and the third basis function
 call overlap_basis_function(bf12,bf3,overlap)

 !
 ! don't forget to destroy it, else memory is leaking
 call destroy_basis_function(bf12)


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
 type(basis_function),intent(in)  :: bf1,bf2
 type(basis_function),intent(out) :: bfprod
!=====
 integer                         :: ig,jg,kg,ng
 real(dp),allocatable            :: coeff(:),alpha(:)
 logical,parameter               :: unnormalized=.FALSE.
 real(dp)                        :: x0_dummy(3)
 integer                         :: fake_shell=1
 integer                         :: fake_index=1
!=====

 !
 ! one could save some primitive gaussians in case of bf1 * bf1
 ! however it is a very small gain
 ng = bf1%ngaussian * bf2%ngaussian
 allocate(coeff(ng),alpha(ng))
 kg=0
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     kg = kg + 1
     alpha(kg) = bf1%g(ig)%alpha + bf2%g(jg)%alpha
     coeff(kg) = bf1%coeff(ig) * bf2%coeff(jg) *  bf1%g(ig)%norm_factor * bf2%g(jg)%norm_factor 
   enddo
 enddo

 call init_basis_function(unnormalized,ng,bf1%nx+bf2%nx,bf1%ny+bf2%ny,bf1%nz+bf2%nz,0,x0_dummy,alpha,coeff,fake_shell,fake_index,bfprod)

 !
 ! override the normalization
 ! the product gaussians are UNnormalized
 ! consistently with the ERI basis
 bfprod%g(:)%norm_factor = 1.0_dp

 deallocate(coeff,alpha)

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
!=====

 ! 
 ! Calculate < phi_1 | r | phi_2 >
 ! using r = ( r - B ) + B
 !

 ! first set up | (x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(1) = dipole_tmp
 ! first set up | Bx phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bx phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(1) = dipole(1) + dipole_tmp * bf2%x0(1)

 ! first set up | (y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(2) = dipole_tmp
 ! first set up | By phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | By phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(2) = dipole(2) + dipole_tmp * bf2%x0(2)

 ! first set up | (z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(3) = dipole_tmp
 ! first set up | Bz phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bz phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(3) = dipole(3) + dipole_tmp * bf2%x0(3)


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
!=====

 ! 
 ! Calculate < phi_1 | x y | phi_2 >
 ! using x y = ( x - Bx ) * ( y - By) + Bx * ( y - By ) + By * ( x - Bx ) + Bx * By
 !



 !
 !  terms x * y and y * x
 ! 

 ! first set up | (x-Bx)*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad_tmp

 ! first set up | Bx*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad(1,2) + quad_tmp * bf2%x0(1)

 ! first set up | By*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
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
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad_tmp

 ! first set up | Bx*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad(1,3) + quad_tmp * bf2%x0(1)

 ! first set up | Bz*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
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
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By)*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad_tmp

 ! first set up | By*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | By*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad(2,3) + quad_tmp * bf2%x0(2)

 ! first set up | Bz*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
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
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+2,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad_tmp

 ! first set up | 2Bx*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2Bx*(x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad(1,1) + quad_tmp * 2.0_dp * bf2%x0(1)

 ! first set up | Bx**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bx**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad(1,1) + quad_tmp * bf2%x0(1)**2


 !
 !  term y * y
 ! 

 ! first set up | (y-By)^2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+2,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad_tmp

 ! first set up | 2B*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2B*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad(2,2) + quad_tmp * 2.0_dp * bf2%x0(2)

 ! first set up | By**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | By**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad(2,2) + quad_tmp * bf2%x0(2)**2


 !
 !  term z * z
 ! 

 ! first set up | (z-Bz)^2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+2,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (z-Bz)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad_tmp

 ! first set up | 2B*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2B*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad(3,3) + quad_tmp * 2.0_dp * bf2%x0(3)

 ! first set up | Bz**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bz**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad(3,3) + quad_tmp * bf2%x0(3)**2



end subroutine basis_function_quadrupole


!=========================================================================
subroutine gos_basis_function(bf1,bf2,qvec,gos_bf1bf2)
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


end subroutine gos_basis_function


!=========================================================================
subroutine setup_cart_to_pure_transforms(gaussian_type)
 implicit none

 character(len=4),intent(in) :: gaussian_type
!=====
 integer  :: il,ni,ii,jj,kk
 integer  :: nx,ny,nz
!=====

 write(stdout,*) 'Setting up the cartesian to pure transforms'

 if(gaussian_type == 'CART') then

   do il=0,lmax_transform
     ni = number_basis_function_am('CART',il)
     allocate(cart_to_pure     (il)%matrix(ni,ni))
     allocate(cart_to_pure_norm(il)%matrix(ni,ni))
     cart_to_pure(il)%matrix(:,:) = 0.0_dp
     do ii=1,ni
       cart_to_pure(il)%matrix(ii,ii) = 1.0_dp
     enddo
   enddo

 else
   ! Formula were read from Ref. 
   ! H.B. Schlegel and M.J. Frisch, INTERNATIONAL JOURNAL OF QUANTUM CHEMISTRY  54, 83-87 (1995).
  
   !
   ! Transform for momentum S
   allocate(cart_to_pure     (0)%matrix(1,1))
   allocate(cart_to_pure_norm(0)%matrix(1,1))
   cart_to_pure(0)%matrix(1,1) = 1.0_dp
  
   !
   ! Transform for momentum P
   allocate(cart_to_pure     (1)%matrix(3,3))
   allocate(cart_to_pure_norm(1)%matrix(3,3))
   cart_to_pure(1)%matrix(:,:) = 0.0_dp
   cart_to_pure(1)%matrix(3,2) = 1.0_dp
   cart_to_pure(1)%matrix(1,3) = 1.0_dp
   cart_to_pure(1)%matrix(2,1) = 1.0_dp
  
   !
   ! Transform for momentum D
   allocate(cart_to_pure     (2)%matrix(6,5))
   allocate(cart_to_pure_norm(2)%matrix(6,5))
   cart_to_pure(2)%matrix(:,:) =  0.0_dp
  
   cart_to_pure(2)%matrix(2,1) =  1.0_dp
  
   cart_to_pure(2)%matrix(5,2) =  1.0_dp
  
   cart_to_pure(2)%matrix(6,3) =  1.0_dp
   cart_to_pure(2)%matrix(1,3) = -0.5_dp
   cart_to_pure(2)%matrix(4,3) = -0.5_dp
  
   cart_to_pure(2)%matrix(3,4) =  1.0_dp
  
   cart_to_pure(2)%matrix(1,5) =  SQRT(3.0/4.0)
   cart_to_pure(2)%matrix(4,5) = -SQRT(3.0/4.0)
  
   !
   ! Transform for momentum F
   allocate(cart_to_pure     (3)%matrix(10,7))
   allocate(cart_to_pure_norm(3)%matrix(10,7))
   cart_to_pure(3)%matrix( :,:) =  0.0_dp
  
   cart_to_pure(3)%matrix(10,4) =  1.0_dp
   cart_to_pure(3)%matrix( 3,4) = -3.0/(2.0*SQRT(5.0))
   cart_to_pure(3)%matrix( 8,4) = -3.0/(2.0*SQRT(5.0))
  
   cart_to_pure(3)%matrix( 6,5) =  SQRT(6.0/5.0)
   cart_to_pure(3)%matrix( 1,5) = -SQRT(6.0)/4.0
   cart_to_pure(3)%matrix( 4,5) = -SQRT(6.0/5.0)/4.0
  
   cart_to_pure(3)%matrix( 9,3) =  SQRT(6.0/5.0)
   cart_to_pure(3)%matrix( 7,3) = -SQRT(6.0)/4.0
   cart_to_pure(3)%matrix( 2,3) = -SQRT(6.0/5.0)/4.0
  
   cart_to_pure(3)%matrix( 3,6) =  SQRT(3.0/4.0)
   cart_to_pure(3)%matrix( 8,6) = -SQRT(3.0/4.0)
  
   cart_to_pure(3)%matrix( 5,2) = -1.0_dp
  
   cart_to_pure(3)%matrix( 1,7) =  SQRT(10.0)/4.0
   cart_to_pure(3)%matrix( 4,7) = -SQRT(2.0)*3.0/4.0
  
   cart_to_pure(3)%matrix( 7,1) = -SQRT(10.0)/4.0
   cart_to_pure(3)%matrix( 2,1) =  SQRT(2.0)*3.0/4.0

   !
   ! Transform for momentum G
   allocate(cart_to_pure     (4)%matrix(15,9))
   allocate(cart_to_pure_norm(4)%matrix(15,9))
   cart_to_pure(4)%matrix( :,:) =  0.0_dp

   cart_to_pure(4)%matrix(15,5) =  1.0_dp
   cart_to_pure(4)%matrix( 1,5) =  3.0/8.0
   cart_to_pure(4)%matrix(11,5) =  3.0/8.0
   cart_to_pure(4)%matrix( 6,5) = -3.0*SQRT(3.0/35.0)
   cart_to_pure(4)%matrix(13,5) = -3.0*SQRT(3.0/35.0)
   cart_to_pure(4)%matrix( 4,5) =  3.0*SQRT(3.0/35.0)/4.0

   cart_to_pure(4)%matrix(10,6) =  SQRT(10.0/7.0)
   cart_to_pure(4)%matrix( 3,6) = -0.75*SQRT(10.0/7.0)
   cart_to_pure(4)%matrix( 8,6) = -0.75*SQRT( 2.0/7.0)

   cart_to_pure(4)%matrix(14,4) =  SQRT(10.0/7.0)
   cart_to_pure(4)%matrix(12,4) = -0.75*SQRT(10.0/7.0)
   cart_to_pure(4)%matrix( 5,4) = -0.75*SQRT( 2.0/7.0)

   cart_to_pure(4)%matrix( 6,7) =  1.5*SQRT( 3.0/7.0)
   cart_to_pure(4)%matrix(13,7) = -1.5*SQRT( 3.0/7.0)
   cart_to_pure(4)%matrix( 1,7) = -0.25*SQRT(5.0)
   cart_to_pure(4)%matrix(11,7) =  0.25*SQRT(5.0)

   cart_to_pure(4)%matrix( 9,3) =  3.0/SQRT(7.0)
   cart_to_pure(4)%matrix( 2,3) = -0.5*SQRT(5.0/7.0)
   cart_to_pure(4)%matrix( 7,3) = -0.5*SQRT(5.0/7.0)

   cart_to_pure(4)%matrix( 3,8) =  0.25*SQRT(10.0)
   cart_to_pure(4)%matrix( 8,8) = -0.75*SQRT(2.0)

   cart_to_pure(4)%matrix(12,2) = -0.25*SQRT(10.0)
   cart_to_pure(4)%matrix( 5,2) =  0.75*SQRT(2.0)

   cart_to_pure(4)%matrix( 1,9) =  0.125*SQRT(35.0)
   cart_to_pure(4)%matrix(11,9) =  0.125*SQRT(35.0)
   cart_to_pure(4)%matrix( 4,9) = -0.75*SQRT(3.0)

   cart_to_pure(4)%matrix( 2,1) =  SQRT(5.0/4.0)
   cart_to_pure(4)%matrix( 7,1) = -SQRT(5.0/4.0)

   !
   ! Transform for momentum H
   allocate(cart_to_pure     (5)%matrix(21,11))
   allocate(cart_to_pure_norm(5)%matrix(21,11))
   cart_to_pure(5)%matrix( :, :) =  0.0_dp

   cart_to_pure(5)%matrix(21, 6) =  1.0_dp
   cart_to_pure(5)%matrix(10, 6) = -5.0*SQRT(2.0/21.0)
   cart_to_pure(5)%matrix(19, 6) = -5.0*SQRT(2.0/21.0)
   cart_to_pure(5)%matrix( 3, 6) =  5.0/8.0*SQRT(2.0)
   cart_to_pure(5)%matrix(17, 6) = -5.0/8.0*SQRT(2.0)
   cart_to_pure(5)%matrix( 8, 6) =  0.25*SQRT(30.0/7.0)

   cart_to_pure(5)%matrix(15, 7) =  SQRT(5.0/3.0)
   cart_to_pure(5)%matrix( 6, 7) = -1.5*SQRT(5.0/7.0)
   cart_to_pure(5)%matrix(13, 7) = -1.5/SQRT(7.0)
   cart_to_pure(5)%matrix( 1, 7) = -0.125*SQRT(15.0)
   cart_to_pure(5)%matrix(11, 7) =  0.125*SQRT(5.0/3.0)
   cart_to_pure(5)%matrix( 4, 7) =  0.25*SQRT(5.0/7.0)

   cart_to_pure(5)%matrix(20, 5) =  SQRT(5.0/3.0)
   cart_to_pure(5)%matrix(18, 5) = -1.5*SQRT(5.0/7.0)
   cart_to_pure(5)%matrix( 9, 5) = -1.5/SQRT(7.0)
   cart_to_pure(5)%matrix(16, 5) = -0.125*SQRT(15.0)
   cart_to_pure(5)%matrix( 2, 5) =  0.125*SQRT(5.0/3.0)
   cart_to_pure(5)%matrix( 7, 5) =  0.25*SQRT(5.0/7.0)

   cart_to_pure(5)%matrix(10, 8) =  SQRT(5.0/4.0)
   cart_to_pure(5)%matrix(19, 8) = -SQRT(5.0/4.0)
   cart_to_pure(5)%matrix( 3, 8) = -0.25*SQRT(35.0/3.0)
   cart_to_pure(5)%matrix(17, 8) =  0.25*SQRT(35.0/3.0)

   cart_to_pure(5)%matrix(14, 4) =  SQRT(5.0/3.0)
   cart_to_pure(5)%matrix( 5, 4) =  -0.5*SQRT(5.0/3.0)
   cart_to_pure(5)%matrix(12, 4) =  -0.5*SQRT(5.0/3.0)

   cart_to_pure(5)%matrix( 6, 9) =  0.5*SQRT(10.0/3.0)
   cart_to_pure(5)%matrix(13, 9) = -0.5*SQRT(6.0)
   cart_to_pure(5)%matrix( 1, 9) = -SQRT(70.0)/16.0
   cart_to_pure(5)%matrix(11, 9) =  SQRT(70.0)/16.0
   cart_to_pure(5)%matrix( 4, 9) =  0.125*SQRT(10.0/3.0)

   cart_to_pure(5)%matrix(18, 3) = -0.5*SQRT(10.0/3.0)
   cart_to_pure(5)%matrix( 9, 3) =  0.5*SQRT(6.0)
   cart_to_pure(5)%matrix(16, 3) =  SQRT(70.0)/16.0
   cart_to_pure(5)%matrix( 2, 3) = -SQRT(70.0)/16.0
   cart_to_pure(5)%matrix( 7, 3) = -0.125*SQRT(10.0/3.0)

   cart_to_pure(5)%matrix( 3,10) =  0.125*SQRT(35.0)
   cart_to_pure(5)%matrix(17,10) =  0.125*SQRT(35.0)
   cart_to_pure(5)%matrix( 8,10) = -0.75*SQRT(3.0)

   cart_to_pure(5)%matrix( 5, 2) =  0.5*SQRT(5.0)
   cart_to_pure(5)%matrix(12, 2) = -0.5*SQRT(5.0)

   cart_to_pure(5)%matrix( 1,11) =  3.0*SQRT(14.0)/16.0
   cart_to_pure(5)%matrix(11,11) = -5.0*SQRT(14.0)/16.0
   cart_to_pure(5)%matrix( 4,11) = -5.0*SQRT(6.0)/8.0

   cart_to_pure(5)%matrix(16, 1) =  3.0*SQRT(14.0)/16.0
   cart_to_pure(5)%matrix( 2, 1) = -5.0*SQRT(14.0)/16.0
   cart_to_pure(5)%matrix( 7, 1) =  5.0*SQRT(6.0)/8.0

   !
   ! Complement with diagonal if necessary
   do il=lmax_transform_pure+1,lmax_transform
     ni = number_basis_function_am('CART',il)
     allocate(cart_to_pure     (il)%matrix(ni,ni))
     allocate(cart_to_pure_norm(il)%matrix(ni,ni))
     cart_to_pure(il)%matrix(:,:) = 0.0_dp
     do ii=1,ni
       cart_to_pure(il)%matrix(ii,ii) = 1.0_dp
     enddo
   enddo

 endif


 !
 ! Copy cart_to_pure into cart_to_pure_norm
 do il=0,lmax_transform
   cart_to_pure_norm(il)%matrix(:,:) = cart_to_pure(il)%matrix(:,:)
 enddo
 !
 ! Then introduce the normalization coefficient part that depends on (nx,ny,nz)
 do il=0,lmax_transform
   kk=0
   do ii=0,il
     nx = il - ii
     do jj=0,ii
       kk = kk + 1
       ny = ii - jj
       nz = jj
       cart_to_pure_norm(il)%matrix(kk,:) = cart_to_pure_norm(il)%matrix(kk,:) &
                 / SQRT( REAL( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) , dp ) )
         
     enddo
   enddo
 enddo


 write(stdout,*) 'Transformations set up completed'
 write(stdout,*) 

end subroutine setup_cart_to_pure_transforms


!=========================================================================
end module m_basis_set
