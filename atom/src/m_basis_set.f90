!=========================================================================
module m_basis_set
 use m_definitions
 use m_timing
 use m_warning
 use m_tools, only: element_name,diagonalize,invert
 use m_gaussian

 real(dp),parameter             :: FILTERED_EIGENVALUE=1.0d-8 ! 1.0d-6

 type basis_function
   character(len=100)           :: basis_name
   integer                      :: am
   character(len=1)             :: amc
   integer                      :: nx,ny,nz
   real(dp)                     :: x0(3)
   integer                      :: ngaussian
   type(gaussian),allocatable   :: g(:) 
   real(dp),allocatable         :: coeff(:)
 end type

 !
 ! A basis set is a list of basis functions
 ! filtering of some elements can be done by rotation 
 type basis_set
   integer                                 :: nbf
   type(basis_function),pointer            :: bf(:) 
   ! additional data needed for product basis
   integer                                 :: nbf_filtered
   integer,allocatable                     :: index_ij(:,:)
   real(dp),allocatable                    :: rotation(:,:)
 end type basis_set

contains

!=========================================================================
 subroutine init_basis_set(PRINT_VOLUME,natom,x,zatom,basis_name,basis_element,basis)
 implicit none
 integer,intent(in)            :: PRINT_VOLUME
 integer,intent(in)            :: natom
 real(dp),intent(in)           :: x(3,natom),zatom(natom)
 character(len=100),intent(in) :: basis_name
 integer,intent(in)            :: basis_element(natom)
 type(basis_set),intent(out)   :: basis
!====
 character(len=100)            :: basis_filename
 integer                       :: ibf,jbf,ng,ig
 real(dp),allocatable          :: alpha(:),coeff(:),coeff2(:)
 logical                       :: file_exists
 integer,parameter             :: basis_file=11
 integer                       :: am_tmp,nbf
 logical,parameter             :: normalized=.TRUE.
 integer                       :: iatom
 real(dp)                      :: x0(3)
!====

 basis%nbf=0
 !
 ! LOOP OVER ATOMS
 ! TODO could be reduced to the type of atoms in the future
 !
 do iatom=1,natom

   write(*,*)
   write(*,*) 'Element used for Z value:    ',TRIM(element_name(zatom(iatom)))
   write(*,*) 'Element used for the basis:  ',TRIM(element_name(REAL(basis_element(iatom),dp)))
   write(*,*) 'Basis type: ',TRIM(basis_name)
   basis_filename=TRIM(element_name(REAL(basis_element(iatom),dp)))//'_'//TRIM(basis_name)
   msg='basis file used: '//basis_filename
   call issue_warning(msg)
  
   write(*,*)
   write(*,*) 'open the basis set file ',TRIM(basis_filename)
   inquire(file=TRIM(basis_filename),exist=file_exists)
   if(.NOT.file_exists) stop'basis set file not found'
  
   !
   ! read first to get all the dimensions
   open(unit=basis_file,file=TRIM(basis_filename),status='old')
   read(basis_file,*) nbf
   if(nbf<1) stop'ERROR in basis set file'
   do ibf=1,nbf
     read(basis_file,*) ng,am_tmp
     if(ng<1) stop'ERROR in basis set file'
     basis%nbf = basis%nbf + number_basis_function_am(am_tmp)
     do ig=1,ng
       read(basis_file,*) 
     enddo
   enddo
   close(basis_file)
  
!   write(*,*) 'Number of basis functions for atom',iatom,nbf

 enddo

 write(*,*) 'Total number of basis functions',basis%nbf
 write(*,*)
 allocate(basis%bf(basis%nbf))

 jbf=0
 do iatom=1,natom

 basis_filename=TRIM(element_name(REAL(basis_element(iatom),dp)))//'_'//TRIM(basis_name)
 open(unit=basis_file,file=TRIM(basis_filename),status='old')
 read(basis_file,*) nbf
 do ibf=1,nbf
   read(basis_file,*) ng,am_tmp
   allocate(alpha(ng),coeff(ng),coeff2(ng))

   if(am_tmp<10) then
     do ig=1,ng
       read(basis_file,*) alpha(ig),coeff(ig)
     enddo
   else
     do ig=1,ng
       read(basis_file,*) alpha(ig),coeff(ig),coeff2(ig)
     enddo
   endif

   ! rescale the gaussian decay rate whenever zatom /= basis_element
   if( abs( zatom(iatom) - REAL(basis_element(iatom),dp) ) > 1.d-6 ) then
     alpha(:) = alpha(:) * ( zatom(iatom) / REAL(basis_element(iatom),dp) )**2
     write(*,*) 'rescaling momentum',am_tmp
     write(*,*) 'smallest rescaled alpha:',MINVAL(alpha(:))
   endif

   x0(:) = x(:,iatom)

   select case(am_tmp)
   case( 0)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,0,x0,alpha,coeff,basis%bf(jbf))
   case( 1)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,0,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,1,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,1,x0,alpha,coeff,basis%bf(jbf))
   case( 2)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,0,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,1,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,0,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,2,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,1,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,2,x0,alpha,coeff,basis%bf(jbf))
   case( 3)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,0,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,1,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,0,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,2,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,1,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,0,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,3,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,2,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,1,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,3,x0,alpha,coeff,basis%bf(jbf))
   case( 4)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,4,0,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,1,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,0,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,2,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,1,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,0,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,3,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,2,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,1,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,0,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,4,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,3,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,2,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,1,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,4,x0,alpha,coeff,basis%bf(jbf))
   case( 5)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,5,0,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,4,1,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,4,0,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,2,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,1,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,0,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,3,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,2,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,1,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,0,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,4,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,3,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,2,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,1,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,0,4,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,5,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,4,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,3,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,2,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,1,4,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,5,x0,alpha,coeff,basis%bf(jbf))
   case( 6)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,6,0,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,5,1,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,5,0,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,4,2,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,4,1,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,4,0,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,3,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,2,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,1,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,3,0,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,4,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,3,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,2,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,1,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,2,0,4,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,5,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,4,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,3,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,2,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,1,4,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,0,5,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,6,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,5,1,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,4,2,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,3,3,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,2,4,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,1,5,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,6,x0,alpha,coeff,basis%bf(jbf))
   case(10)
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,0,x0,alpha,coeff,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,1,0,0,x0,alpha,coeff2,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,1,0,x0,alpha,coeff2,basis%bf(jbf))
     jbf=jbf+1 ; call init_basis_function(normalized,ng,0,0,1,x0,alpha,coeff2,basis%bf(jbf))
   case default
     stop'not implemented'
   end select


   deallocate(alpha,coeff,coeff2)
 enddo
 close(basis_file)

 !
 ! END OF THE LOOP OVER ATOMS
 enddo
 

 !
 ! finally output the basis set upon request
 if(PRINT_VOLUME>5) then
   do ibf=1,basis%nbf
     write(*,*) ' basis function number',ibf
     call print_basis_function(basis%bf(ibf))
   enddo
 endif

 write(*,*) 'Basis set is ready and fit'

 end subroutine init_basis_set

!=========================================================================
 subroutine init_product_basis_set(basis,prod_basis)
 implicit none
 type(basis_set),intent(in)     :: basis
 type(basis_set),intent(out)    :: prod_basis
!====
 integer                        :: ibf,jbf,iprodbf,jprodbf
 real(dp)                       :: overlap_tmp,norm_tmp
 real(dp),allocatable           :: s_matrix(:,:),eigval(:),eigvec(:,:)
 real(dp),allocatable           :: rotation_tmp(:,:),tmp(:)
 character(len=100)             :: title
!====

 prod_basis%nbf = ( basis%nbf * ( basis%nbf + 1 ) ) / 2
 allocate(prod_basis%bf(prod_basis%nbf))
 allocate(prod_basis%index_ij(2,prod_basis%nbf))

 !
 ! Construct all products
 iprodbf = 0
 do jbf=1,basis%nbf
   do ibf=1,jbf
     iprodbf = iprodbf + 1
     prod_basis%index_ij(1,iprodbf) = ibf
     prod_basis%index_ij(2,iprodbf) = jbf
     call basis_function_prod(basis%bf(ibf),basis%bf(jbf),prod_basis%bf(iprodbf)) 
   enddo
 enddo


#ifdef AUXIL_BASIS
 allocate(s_matrix(prod_basis%nbf,prod_basis%nbf))
 allocate(eigval(prod_basis%nbf))
 allocate(eigvec(prod_basis%nbf,prod_basis%nbf))
 !
 ! first build the full overlap matrix S
 ! in order to identify the most singular eigenvectors
 do jprodbf=1,prod_basis%nbf
   do iprodbf=1,jprodbf  ! the matrix is symmetric S_ab = S_ba
     call overlap_basis_function(prod_basis%bf(iprodbf),prod_basis%bf(jprodbf),overlap_tmp)
     s_matrix(iprodbf,jprodbf) = overlap_tmp
     s_matrix(jprodbf,iprodbf) = overlap_tmp
   enddo
 enddo

 call diagonalize(prod_basis%nbf,s_matrix,eigval,eigvec)


 prod_basis%nbf_filtered = 0
 do iprodbf=1,prod_basis%nbf
   if( eigval(iprodbf) > FILTERED_EIGENVALUE ) prod_basis%nbf_filtered = prod_basis%nbf_filtered + 1
 enddo

 write(*,'(a,es12.6)') ' filtering below ',FILTERED_EIGENVALUE
 write(*,'(a,i4,a,i4)') ' Conserve ',prod_basis%nbf_filtered,' out of ',prod_basis%nbf

 allocate(prod_basis%rotation(prod_basis%nbf_filtered,prod_basis%nbf))

 jprodbf=0
 do iprodbf=1,prod_basis%nbf
   if( eigval(iprodbf) > FILTERED_EIGENVALUE ) then
     jprodbf=jprodbf+1
     prod_basis%rotation(jprodbf,:) = eigvec(:,iprodbf)    ! THIS HAS BEEN CHECKED ALREADY TWICE!
   endif
 enddo

 deallocate(s_matrix,eigval,eigvec)
#endif

 end subroutine init_product_basis_set


!=========================================================================
 subroutine destroy_basis_set(basis)
 implicit none
 type(basis_set),intent(inout) :: basis
!====

 deallocate(basis%bf)
 if(allocated(basis%index_ij)) deallocate(basis%index_ij)
 if(allocated(basis%rotation)) deallocate(basis%rotation)

 end subroutine destroy_basis_set

!=========================================================================
 subroutine init_basis_function(normalized,ng,nx,ny,nz,x0,alpha,coeff,bf)
 implicit none
 logical,intent(in)               :: normalized
 integer,intent(in)               :: ng,nx,ny,nz
 real(dp),intent(in)              :: x0(3),alpha(ng)
 real(dp),intent(in)              :: coeff(ng)
 type(basis_function),intent(out) :: bf
!====
 integer                          :: ig
 real(dp)                         :: overlap
!====

 bf%ngaussian = ng
 allocate(bf%g(bf%ngaussian))
 allocate(bf%coeff(bf%ngaussian))
 bf%nx    = nx
 bf%ny    = ny
 bf%nz    = nz
 bf%am    = nx + ny + nz
 bf%amc   = orbital_momentum_name(bf%am)
 bf%x0(:) = x0(:)

 ! All the gaussians of the contraction have the same orbital momentum
 do ig=1,bf%ngaussian
   call init_gaussian_general(nx,ny,nz,alpha(ig),x0,bf%g(ig))
   bf%coeff(ig) = coeff(ig)
 enddo

 !
 ! check the normalization if requested
 if( normalized ) then
   call overlap_basis_function(bf,bf,overlap)
   if( ABS(overlap-1.0_dp) > 2.0d-5 ) then
     write(*,*) 'normalization is different from 1.0',overlap
     write(*,*) bf%nx,bf%ny,bf%nz
     write(*,*) 'assuming this is a generalized contraction and rescaling coefficients'
     bf%coeff(:) = coeff(:) / SQRT( overlap )
   endif
 endif
 

 end subroutine init_basis_function

!=========================================================================
 subroutine destroy_basis_function(bf)
 implicit none
 type(basis_function),intent(inout) :: bf
!====
 
 deallocate(bf%g,bf%coeff)

 end subroutine destroy_basis_function

!=========================================================================
 function number_basis_function_am(am)
 integer,intent(in) :: am
 integer            :: number_basis_function_am
!=====

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
 case(10) ! stands for SP orbitals
   number_basis_function_am = 4 
 case default
   stop'number_basis_function_am: not implemented'
 end select

 end function number_basis_function_am

!=========================================================================
 subroutine print_basis_function(bf)
 implicit none
 type(basis_function),intent(in) :: bf
!====
 integer :: ig
!====

 write(*,*)
 write(*,*) '======  print out a basis function ======'
 write(*,'(a30,2x,1(x,i3))')           'contraction of N gaussians',bf%ngaussian
 write(*,'(a30,5x,a1)')                'orbital momentum',bf%amc
 write(*,'(a30,x,3(f12.6,2x))')        'centered in',bf%x0(:)
 do ig=1,bf%ngaussian
   write(*,'(a30,2x,x,i3,2x,f12.6)')   'coefficient',ig,bf%coeff(ig)
 enddo
 write(*,*)
 do ig=1,bf%ngaussian
   call print_gaussian(bf%g(ig))
 enddo
 write(*,*) '====== end of basis function ======'
 write(*,*)

 end subroutine print_basis_function

!=========================================================================
 function eval_basis_function(bf,x)
 implicit none
 type(basis_function),intent(in) :: bf
 real(dp),intent(in)             :: x(3)
 real(dp)                        :: eval_basis_function
!====
 integer                         :: ig
!====

 eval_basis_function=0.0_dp
 do ig=1,bf%ngaussian
   eval_basis_function = eval_basis_function + eval_gaussian(bf%g(ig),x) * bf%coeff(ig)
 enddo

 end function eval_basis_function

!=========================================================================
 function eval_basis_function_derivative(bf,x)
 implicit none
 type(basis_function),intent(in) :: bf
 real(dp),intent(in)             :: x(3)
 real(dp)                        :: eval_basis_function_derivative(3)
!====
 integer                         :: ig
!====

 eval_basis_function_derivative(:)=0.0_dp
 do ig=1,bf%ngaussian
   eval_basis_function_derivative(:) = eval_basis_function_derivative(:) + eval_gaussian_derivative(bf%g(ig),x) * bf%coeff(ig)
 enddo

 end function eval_basis_function_derivative


!=========================================================================
 subroutine overlap_basis_function(bf1,bf2,overlap)
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(out)            :: overlap
!====
 integer                         :: ig,jg
 real(dp)                        :: overlap_one_gaussian
!====

 overlap=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
#ifdef ATOM
     call overlap_normalized(bf1%g(ig),bf2%g(jg),overlap_one_gaussian)
#else
     call overlap_recurrence(bf1%g(ig),bf2%g(jg),overlap_one_gaussian)
#endif
     overlap = overlap + overlap_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


 end subroutine overlap_basis_function

!=========================================================================
 subroutine overlap_three_basis_function(bf1,bf2,bf3,overlap)
 implicit none
 type(basis_function),intent(in) :: bf1,bf2,bf3
 real(dp),intent(out)            :: overlap
!====
 type(basis_function)            :: bf12
 integer                         :: ig,jg
 real(dp)                        :: overlap_one_gaussian
!====

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
!====
 integer                         :: ig,jg
 real(dp)                        :: kinetic_one_gaussian
!====

 kinetic=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
#ifdef ATOM
     call kinetic_gaussian(bf1%g(ig),bf2%g(jg),kinetic_one_gaussian)
#else
     call kinetic_recurrence(bf1%g(ig),bf2%g(jg),kinetic_one_gaussian)
#endif
     kinetic = kinetic + kinetic_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


 end subroutine kinetic_basis_function

!=========================================================================
 subroutine nucleus_pot_basis_function(bf1,bf2,zatom,x,nucleus_pot)
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(in)             :: zatom,x(3)
 real(dp),intent(out)            :: nucleus_pot
!====
 integer                         :: ig,jg
 real(dp)                        :: nucleus_pot_one_gaussian
!====

 nucleus_pot=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
#ifdef ATOM
     call nucleus_pot_gaussian(bf1%g(ig),bf2%g(jg),zatom,nucleus_pot_one_gaussian)
#else
     call nucleus_recurrence(zatom,x,bf1%g(ig),bf2%g(jg),nucleus_pot_one_gaussian)
#endif
     nucleus_pot = nucleus_pot + nucleus_pot_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


 end subroutine nucleus_pot_basis_function

!=========================================================================
 subroutine basis_function_prod(bf1,bf2,bfprod)
 implicit none
 type(basis_function),intent(in)  :: bf1,bf2
 type(basis_function),intent(out) :: bfprod
!====
 integer                         :: ig,jg,kg,ng
 real(dp),allocatable            :: coeff(:),alpha(:)
 logical,parameter               :: unnormalized=.FALSE.
 real(dp)                        :: x0_dummy(3)
!====

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

 call init_basis_function(unnormalized,ng,bf1%nx+bf2%nx,bf1%ny+bf2%ny,bf1%nz+bf2%nz,x0_dummy,alpha,coeff,bfprod)

! call print_basis_function(bfprod) ; stop'DEBUG'

 !
 ! override the normalization
 ! the product gaussians are UNnormalized
 ! consistently with the ERI basis
 bfprod%g(:)%norm_factor = 1.0_dp

 deallocate(coeff,alpha)

 end subroutine basis_function_prod

!=========================================================================
end module m_basis_set
