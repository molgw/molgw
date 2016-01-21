!=========================================================================
module m_eri
 use m_definitions
 use m_mpi
 use m_memory
 use m_basis_set
 use m_timing
 use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE


 real(dp),parameter,public :: TOO_LOW_EIGENVAL=1.0e-6_dp

 integer,parameter :: BUFFER1 = 1
 integer,parameter :: BUFFER2 = 2
 !
 ! max length of a record in the ERI file
 integer,parameter,private :: line_length=1000

 real(dp),private           :: TOL_INT

 real(prec_eri),public,allocatable :: eri_4center(:)
 real(prec_eri),public,allocatable :: eri_4center_lr(:)
 real(prec_eri),public,allocatable :: eri_3center(:,:)
 real(prec_eri),public,allocatable :: eri_3center_lr(:,:)

 real(prec_eri),protected,allocatable :: eri_3center_eigen(:,:,:,:)
!FBFB LW
 real(prec_eri),protected,allocatable :: eri_3center_eigen_mixed(:,:,:,:)

 logical,protected,allocatable      :: negligible_basispair(:,:)
 logical,protected,allocatable      :: negligible_shellpair(:,:)
 integer,protected,allocatable      :: index_pair(:,:)
 integer,protected,allocatable      :: index_basis(:,:)
 integer,protected,allocatable      :: index_shellpair(:,:)
 integer,protected                  :: nshellpair

 type shell_type
   integer              :: am
   integer              :: ng
   real(dp),allocatable :: alpha(:)
   real(dp),allocatable :: coeff(:)
   real(dp)             :: x0(3)
   integer              :: istart,iend
 end type shell_type
 protected :: shell_type

 integer,protected                      :: nshell
 integer,protected                      :: nshell_auxil
 type(shell_type),protected,allocatable :: shell(:)
 type(shell_type),protected,allocatable :: shell_auxil(:)


 integer,private   :: nbf_eri         ! local copy of nbf
 integer,protected :: nsize           ! size of the eri_4center array
 integer,protected :: npair         ! number of independent pairs (i,j) with i<=j 

 integer,public    :: nauxil_3center     ! size of the 3-center matrix
                                         ! may differ from the previous one due to
                                         ! data distribution
 integer,public    :: nauxil_3center_lr  ! size of the 3-center matrix
                                         ! may differ from the previous one due to
                                         ! data distribution


! TODO write a proper interface for the call to C
! interface
!   integer(C_INT) function eval_contr_integral() bind(C)
!       info=eval_contr_integral(                &
!                               am1,am2,am3,am4, &
!                               ng1,ng2,ng3,ng4, &
!                               coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
!                               alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
!                               x01(1),x02(1),x03(1),x04(1),&
!                               rcut,
!                               int_shell(1)
!     character(kind=c_char) :: string(*)
!   end subroutine print_c
! end interface


contains


!=========================================================================
subroutine prepare_eri(basis,rcut,which_buffer)
 use m_inputparam,only: integral_level
 implicit none
!===== 
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: rcut
 integer,intent(in)         :: which_buffer
!===== 
 logical            :: file_exists
!===== 

 nbf_eri = basis%nbf

 select case(integral_level)
 case(low)       ! accuracy not guaranted, just for quick test runs
   TOL_INT = 1.0e-04_dp
 case(medium)    ! 10 meV accuracy on potentials
   TOL_INT = 1.0e-06_dp
 case(high)      !  1 meV accuracy on potentials
   TOL_INT = 1.0e-08_dp
 case(very_high) ! almost perfect potentials
   TOL_INT = 1.0e-10_dp
 case(insane)    ! No screening of any integral
   TOL_INT = 0.0_dp
 case default
   call die('integration quality not recognized')
 end select
 write(stdout,'(/,a,es9.2)') ' Tolerance on integrals set to ',TOL_INT


 if(.NOT.ALLOCATED(negligible_shellpair)) then
   call setup_shell_list(basis)
   allocate(negligible_shellpair(nshell,nshell))
   allocate(negligible_basispair(nbf_eri,nbf_eri))
   allocate(index_pair(nbf_eri,nbf_eri))
   call identify_negligible_shellpair(basis)
   call setup_shellpair()
   call setup_negligible_basispair()
 endif


 nsize = (npair*(npair+1))/2


end subroutine prepare_eri


!=========================================================================
subroutine deallocate_eri_4center()
 implicit none
!=====

 if(ALLOCATED(eri_4center)) then
   write(stdout,'(/,a)')     ' Deallocate ERI buffer'
   call clean_deallocate('4-center integrals',eri_4center)
 endif

end subroutine deallocate_eri_4center


!=========================================================================
subroutine deallocate_eri_4center_lr()
 implicit none
!=====

 if(ALLOCATED(eri_4center_lr)) then
   write(stdout,'(/,a)')     ' Deallocate LR ERI buffer'
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif

end subroutine deallocate_eri_4center_lr



!=========================================================================
subroutine deallocate_eri()
 implicit none

 integer :: ishell
!=====

 if(ALLOCATED(eri_4center)) then
   write(stdout,'(/,a)')     ' Deallocate ERI buffer'
   call clean_deallocate('4-center integrals',eri_4center)
 endif
 if(ALLOCATED(eri_4center_lr)) then
   write(stdout,'(/,a)')     ' Deallocate LR ERI buffer'
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif
 if(ALLOCATED(negligible_basispair))  deallocate(negligible_basispair)
 if(ALLOCATED(negligible_shellpair))  deallocate(negligible_shellpair)
 if(ALLOCATED(index_pair))            deallocate(index_pair)
 if(ALLOCATED(index_basis))           deallocate(index_basis)
 if(ALLOCATED(index_shellpair))       deallocate(index_shellpair)
 ! 
 ! Cleanly deallocate the shell objects
 do ishell=1,nshell
   if(ALLOCATED(shell(ishell)%alpha)) deallocate( shell(ishell)%alpha )
   if(ALLOCATED(shell(ishell)%coeff)) deallocate( shell(ishell)%coeff )
 enddo
 if(ALLOCATED(shell))                 deallocate(shell)


end subroutine deallocate_eri


!=========================================================================
function index_prod(ibf,jbf)
 implicit none
 integer,intent(in) :: ibf,jbf
 integer            :: index_prod
!=====
 integer            :: jmin,imax
!=====

 index_prod = index_pair(ibf,jbf)

end function index_prod


!=========================================================================
function index_eri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 integer            :: index_eri
!=====
! integer            :: imin,jmax,kmin,lmax
 integer            :: klmin,ijmax
 integer            :: index_ij,index_kl
!===== 

 index_ij = index_prod(ibf,jbf)
 index_kl = index_prod(kbf,lbf)

 ijmax=MAX(index_ij,index_kl)
 klmin=MIN(index_ij,index_kl)

 index_eri = (klmin-1)*npair - (klmin-1)*(klmin-2)/2 + ijmax-klmin+1

! index_eri = ibf+(jbf-1)*nbf_eri+(kbf-1)*nbf_eri**2+(lbf-1)*nbf_eri**3

end function index_eri


!=========================================================================
function eri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri = 0.0_dp
 else
   eri = eri_4center(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri


!=========================================================================
function eri_lr(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_lr
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_lr = 0.0_dp
 else
   eri_lr = eri_4center_lr(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri_lr


!=========================================================================
function eri_ri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_ri
!=====
 integer            :: index_ij,index_kl
 real(dp)           :: eri_1(1,1)
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_ri = 0.0_dp
 else
   index_ij = index_prod(ibf,jbf)
   index_kl = index_prod(kbf,lbf)
  
   eri_ri = DOT_PRODUCT( eri_3center(:,index_ij) , eri_3center(:,index_kl) )

   call xsum(eri_ri)

 endif

end function eri_ri


!=========================================================================
function eri_ri_lr(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_ri_lr
!=====
 integer            :: index_ij,index_kl
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_ri_lr = 0.0_dp
 else
   index_ij = index_prod(ibf,jbf)
   index_kl = index_prod(kbf,lbf)

   eri_ri_lr = DOT_PRODUCT( eri_3center_lr(:,index_ij) , eri_3center_lr(:,index_kl) )

   call xsum(eri_ri_lr)

 endif

end function eri_ri_lr


!=========================================================================
function eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen_ri
!=====

 eri_eigen_ri = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

 call xsum(eri_eigen_ri)

end function eri_eigen_ri


!=========================================================================
function eri_eigen_ri_paral(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen_ri_paral
!=====

 eri_eigen_ri_paral = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

end function eri_eigen_ri_paral


!=========================================================================
subroutine setup_shell_list(basis)
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer :: ibf,jbf
 integer :: ishell
!=====


 nshell = basis%nshell
 allocate(shell(nshell))

 !
 ! Set up shells information
 jbf=0
 do ishell=1,nshell
   do ibf=1,basis%nbf_cart
     if(basis%bf(ibf)%shell_index==ishell) then
       shell(ishell)%am    = basis%bf(ibf)%am
       shell(ishell)%x0(:) = basis%bf(ibf)%x0(:)
       shell(ishell)%ng    = basis%bf(ibf)%ngaussian
       allocate( shell(ishell)%alpha(shell(ishell)%ng) )
       allocate( shell(ishell)%coeff(shell(ishell)%ng) )
       shell(ishell)%alpha(:) = basis%bf(ibf)%g(:)%alpha
       !
       ! Include here the normalization part that does not depend on (nx,ny,nz)
       shell(ishell)%coeff(:) = basis%bf(ibf)%coeff(:) &
                 * ( 2.0_dp / pi )**0.75_dp * 2.0_dp**shell(ishell)%am * shell(ishell)%alpha(:)**( 0.25_dp * ( 2.0_dp*shell(ishell)%am + 3.0_dp ) )

       jbf = jbf + 1
       shell(ishell)%istart = jbf
       jbf = jbf + number_basis_function_am( basis%gaussian_type , shell(ishell)%am ) - 1
       shell(ishell)%iend   = jbf
       exit

     endif
   enddo
 enddo

end subroutine setup_shell_list


!=========================================================================
subroutine setup_shell_list_auxil(auxil_basis)
 implicit none
 
 type(basis_set),intent(in)   :: auxil_basis
!=====
 integer :: ibf,jbf
 integer :: ishell
!=====


 nshell_auxil = auxil_basis%nshell
 allocate(shell_auxil(nshell_auxil))

 !
 ! Set up shells information
 jbf=0
 do ishell=1,nshell_auxil
   do ibf=1,auxil_basis%nbf_cart
     if(auxil_basis%bf(ibf)%shell_index==ishell) then
       shell_auxil(ishell)%am    = auxil_basis%bf(ibf)%am
       shell_auxil(ishell)%x0(:) = auxil_basis%bf(ibf)%x0(:)
       shell_auxil(ishell)%ng    = auxil_basis%bf(ibf)%ngaussian
       allocate( shell_auxil(ishell)%alpha(shell_auxil(ishell)%ng) )
       allocate( shell_auxil(ishell)%coeff(shell_auxil(ishell)%ng) )
       shell_auxil(ishell)%alpha(:) = auxil_basis%bf(ibf)%g(:)%alpha
       !
       ! Include here the normalization part that does not depend on (nx,ny,nz)
       shell_auxil(ishell)%coeff(:) = auxil_basis%bf(ibf)%coeff(:) &
                 * ( 2.0_dp / pi )**0.75_dp * 2.0_dp**shell_auxil(ishell)%am * shell_auxil(ishell)%alpha(:)**( 0.25_dp * ( 2.0_dp*shell_auxil(ishell)%am + 3.0_dp ) )

       jbf = jbf + 1
       shell_auxil(ishell)%istart = jbf
       jbf = jbf + number_basis_function_am( auxil_basis%gaussian_type , shell_auxil(ishell)%am ) - 1
       shell_auxil(ishell)%iend   = jbf
       exit

     endif
   enddo
 enddo

end subroutine setup_shell_list_auxil


!=========================================================================
subroutine setup_negligible_basispair()
 implicit none
!=====
 integer :: ishell,jshell
 integer :: ibf,jbf
!=====

 negligible_basispair(:,:) = .FALSE.
 do ishell=1,nshell
   do jshell=1,nshell
     if( negligible_shellpair(ishell,jshell) ) then 
       do ibf=shell(ishell)%istart,shell(ishell)%iend
         do jbf=shell(jshell)%istart,shell(jshell)%iend
           negligible_basispair(ibf,jbf) = .TRUE.
           negligible_basispair(jbf,ibf) = .TRUE.
         enddo
       enddo
     endif
   enddo
 enddo

 npair = 0
 do jbf=1,nbf_eri
   do ibf=1,jbf
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
     endif
   enddo
 enddo
 allocate(index_basis(2,npair))

 npair = 0
 index_pair(:,:) = 0
 do jbf=1,nbf_eri
   do ibf=1,jbf
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
       index_pair(ibf,jbf) = npair
       index_pair(jbf,ibf) = npair
       index_basis(1,npair) = ibf
       index_basis(2,npair) = jbf
     endif
   enddo
 enddo


end subroutine setup_negligible_basispair


!=========================================================================
subroutine refine_negligible_basispair()
 implicit none

!=====
 integer  :: ibf,jbf,kbf,lbf
 integer  :: npair,npair_refined
 real(dp) :: max_ij
!=====

 npair         = 0
 npair_refined = 0

 do jbf=1,nbf_eri
   do ibf=1,jbf
     if( negligible_basispair(ibf,jbf) ) cycle
     npair = npair + 1

     max_ij=0.0_dp
     do lbf=1,nbf_eri
       do kbf=1,lbf
         if( negligible_basispair(kbf,lbf) ) cycle
         max_ij = MAX( max_ij , ABS(eri(ibf,jbf,kbf,lbf)) )
       enddo
     enddo
     if( max_ij < TOL_INT ) then
!       write(stdout,*) '    negl',max_ij,max_ij < TOL_INT
       negligible_basispair(ibf,jbf) = .TRUE.
     else
!       write(stdout,*) 'non negl',max_ij,max_ij < TOL_INT
       npair_refined = npair_refined + 1
     endif


   enddo
 enddo

 write(stdout,'(/,a)') ' Refining the negligible basis function pairs'
 write(stdout,'(a,x,i6)')   ' Non negligible basis function pairs stored in memory   ',npair
 write(stdout,'(a,x,i6)')   ' Non negligible basis function pairs used in calculation',npair_refined


end subroutine refine_negligible_basispair


!=========================================================================
subroutine identify_negligible_shellpair(basis)
!
! A first screening implementation
! Find negligible shell pair with
! Cauchy-Schwarz inequality
! (ij|1/r|kl)**2 <= (ij|1/r|ij) (kl|1/r|(kl) 
!
 use m_tools,only: boys_function
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer  :: info
 integer  :: iibf
 integer  :: ibf,jbf,kbf,lbf
 integer  :: n1c,n2c
 integer  :: ni,nj
 integer  :: ami,amj
 integer  :: ishell,jshell
 integer  :: ig1,ig2,ig3,ig4
 real(dp) :: zeta_12,rho,rho1,f0t(0:0),tt
 real(dp) :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
!=====
! variables used to call C
 integer(C_INT),external      :: eval_contr_integral
 integer(C_INT)               :: am1,am2
 integer(C_INT)               :: ng1,ng2
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:)
 real(C_DOUBLE)               :: x01(3),x02(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:)
 real(C_DOUBLE)               :: rcut_libint
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 call start_clock(timing_eri_screening)
 write(stdout,'(/,a)')    ' Cauchy-Schwartz screening of the 3- or 4-center integrals'

 rcut_libint = 0.0_dp

 negligible_shellpair(:,:) = .TRUE.

 do jshell=1,nshell
   ! Workload is distributed here
   if( MODULO(jshell-1,nproc) /= rank ) cycle

   do ishell=1,nshell
     ami = shell(ishell)%am
     amj = shell(jshell)%am
     !TODO: Here time could be saved by only checking ishell<= jshell
     ! But then an interexchange of indexes would have to be implemented in
     ! order to satisfy ami >= amj (required condition in libint)
     if( ami < amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     n1c = number_basis_function_am( 'CART' , ami )
     n2c = number_basis_function_am( 'CART' , amj )
     am1 = shell(ishell)%am
     am2 = shell(jshell)%am
     ng1 = shell(ishell)%ng
     ng2 = shell(jshell)%ng

     allocate(alpha1(ng1),alpha2(ng2))
     allocate(coeff1(ng1),coeff2(ng2))
     alpha1(:) = shell(ishell)%alpha(:)
     alpha2(:) = shell(jshell)%alpha(:)
     x01(:) = shell(ishell)%x0(:)
     x02(:) = shell(jshell)%x0(:)
     coeff1(:) = shell(ishell)%coeff(:)
     coeff2(:) = shell(jshell)%coeff(:)

     allocate( int_shell( n1c*n2c*n1c*n2c ) )
     allocate( integrals_cart(n1c,n2c,n1c,n2c) )
     allocate( integrals_tmp (n1c,n2c,n1c,n2c) )

     integrals_cart(:,:,:,:) = 0.0_dp

     if(ami+amj==0) then

       do ig4=1,ng2
         do ig3=1,ng1
           do ig2=1,ng2
             do ig1=1,ng1

               zeta_12 = alpha1(ig1) + alpha2(ig2)
               p(:) = ( alpha1(ig1) * x01(:) + alpha2(ig2) * x02(:) ) / zeta_12 
               q(:) = ( alpha1(ig3) * x01(:) + alpha2(ig4) * x02(:) ) / zeta_12 
               !
               ! Treat carefully the LR only integrals
               rho  = zeta_12 * zeta_12 / ( zeta_12 + zeta_12 + zeta_12*zeta_12*rcut_libint**2 )
               rho1 = zeta_12 * zeta_12 / ( zeta_12 + zeta_12 )
               tt = rho * SUM( (p(:)-q(:))**2 )
               call boys_function(f0t(0),0,tt)

               integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                     2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_12 ) * f0t(0) &
                     / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                     / zeta_12 * EXP( -alpha1(ig3)*alpha2(ig4)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) &
                     * SQRT( rho / rho1 ) &
                     * coeff1(ig1) &
                     * coeff2(ig2) &
                     * coeff1(ig3) &
                     * coeff2(ig4) * cart_to_pure_norm(0)%matrix(1,1)**4

             enddo
           enddo
         enddo
       enddo

     else

       info=eval_contr_integral(                &
                               am1,am2,am1,am2, &
                               ng1,ng2,ng1,ng2, &
                               coeff1(1),coeff2(1),coeff1(1),coeff2(1),&
                               alpha1(1),alpha2(1),alpha1(1),alpha2(1),&
                               x01(1),x02(1),x01(1),x02(1),&
                               rcut_libint, &
                               int_shell(1))


       if(info/=0) then
         write(stdout,*) am1,am2,am1,am2
         call die('ERI calculated by libint failed')
       endif

       iibf=0
       do ibf=1,n1c
         do jbf=1,n2c
           do kbf=1,n1c
             do lbf=1,n2c
               iibf=iibf+1
               integrals_cart(ibf,jbf,kbf,lbf) = int_shell(iibf)
             enddo
           enddo
         enddo
       enddo


       do lbf=1,n2c
         do kbf=1,n1c
           do jbf=1,n2c
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1c,jbf,kbf,lbf) * cart_to_pure_norm(shell(ishell)%am)%matrix(1:n1c,ibf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n2c
         do kbf=1,n1c
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2c,kbf,lbf) * cart_to_pure_norm(shell(jshell)%am)%matrix(1:n2c,jbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n2c
         do kbf=1,ni
           do jbf=1,nj
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n1c,lbf) * cart_to_pure_norm(shell(ishell)%am)%matrix(1:n1c,kbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,nj
         do kbf=1,ni
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n2c) * cart_to_pure_norm(shell(jshell)%am)%matrix(1:n2c,lbf) )
             enddo
           enddo
         enddo
       enddo
            
     endif

     do ibf=1,ni
       do jbf=1,nj
         if( ABS( integrals_cart(ibf,jbf,ibf,jbf) ) > TOL_INT**2 ) negligible_shellpair(ishell,jshell) = .FALSE.
       enddo
     enddo

     !
     ! Symmetrize
     negligible_shellpair(jshell,ishell) = negligible_shellpair(ishell,jshell)

     deallocate(integrals_cart)
     deallocate(integrals_tmp)
     deallocate(int_shell)
     deallocate(alpha1,alpha2)
     deallocate(coeff1,coeff2)

   enddo
 enddo

 call xand(negligible_shellpair)

 call stop_clock(timing_eri_screening)


end subroutine identify_negligible_shellpair


!=========================================================================
subroutine setup_shellpair()
 implicit none

 integer :: ishell,jshell
 integer :: ami,amj
 integer :: ishellpair,jshellpair
!=====

 ishellpair = 0
 jshellpair = 0
 do jshell=1,nshell
   do ishell=1,jshell ! nshell
     jshellpair = jshellpair + 1
     ! skip the identified negligible shell pairs
     if( negligible_shellpair(ishell,jshell) ) cycle
     ami = shell(ishell)%am
     amj = shell(jshell)%am
     ishellpair = ishellpair + 1

   enddo
 enddo
 nshellpair = ishellpair
 write(stdout,'(/,a,i8,a,i8)') ' Non negligible shellpairs to be computed',nshellpair,'  over a total of',jshellpair
 allocate(index_shellpair(2,nshellpair))

 ishellpair = 0
 do jshell=1,nshell
   do ishell=1,jshell ! nshell
     ! skip the identified negligible shell pairs
     if( negligible_shellpair(ishell,jshell) ) cycle
     ami = shell(ishell)%am
     amj = shell(jshell)%am
     ishellpair = ishellpair + 1
     ! Reverse if needed the order of the shell so to maximize the angular
     ! momentum of the first shell
     if( ami >= amj ) then
       index_shellpair(1,ishellpair) = ishell
       index_shellpair(2,ishellpair) = jshell
     else
       index_shellpair(1,ishellpair) = jshell
       index_shellpair(2,ishellpair) = ishell
     endif

   enddo
 enddo


end subroutine setup_shellpair


!=========================================================================
subroutine test_eri(basis)
 implicit none
 type(basis_set),intent(in)   :: basis
!=====
 integer                      :: ibf,jbf,kbf,lbf
!=====

 do jbf=1,nbf_eri
   do ibf=1,nbf_eri
     do lbf=1,nbf_eri
       do kbf=1,nbf_eri
         if( ABS(eri(ibf,jbf,kbf,lbf) - eri(kbf,lbf,ibf,jbf)) > 1.d-6 ) then
           write(stdout,*) ibf,jbf,kbf,lbf,eri(ibf,jbf,kbf,lbf)
           write(stdout,*) kbf,lbf,ibf,jbf,eri(kbf,lbf,ibf,jbf)
           write(stdout,*) ibf,basis%bf(ibf)%amc
           write(stdout,*) jbf,basis%bf(jbf)%amc
           write(stdout,*) kbf,basis%bf(kbf)%amc
           write(stdout,*) lbf,basis%bf(lbf)%amc
           call die('ERI array not symmetric')
         endif
       enddo
     enddo
   enddo
 enddo

 call die('TESTING OK')

end subroutine test_eri


!=================================================================
subroutine transform_eri_basis(nspin,c_matrix,istate,ijspin,eri_eigenstate_i)
 implicit none

 integer,intent(in)     :: nspin,istate,ijspin
 real(dp),intent(in)    :: c_matrix(nbf_eri,nbf_eri,nspin)
 real(dp),intent(inout) :: eri_eigenstate_i(nbf_eri,nbf_eri,nbf_eri,nspin)
!=====
 integer,save         :: istate_previous=0
 integer,save         :: ijspin_previous=0
 integer              :: klspin
 integer              :: ibf,jbf,kbf,lbf
 integer              :: jstate,kstate,lstate
 real(dp)             :: eri_tmp3(nbf_eri,nbf_eri,nbf_eri)
 real(dp)             :: wtime
!=====

 ! Check if the calculation can be skipped
 if( istate_previous == istate .AND. ijspin_previous == ijspin .AND. ANY(ABS(eri_eigenstate_i(:,:,:,:))>1.0e-6_dp) ) then
   return
 else
   istate_previous = istate
   ijspin_previous = ijspin
 endif


 call start_clock(timing_basis_transform)

 eri_eigenstate_i(:,:,:,:)=0.0_dp
 eri_tmp3(:,:,:)=0.0_dp

 do lbf=1,nbf_eri
   do kbf=1,nbf_eri
     do jbf=1,nbf_eri

       do ibf=1,nbf_eri
         eri_tmp3(jbf,kbf,lbf) = eri_tmp3(jbf,kbf,lbf) + eri(ibf,jbf,kbf,lbf) * c_matrix(ibf,istate,ijspin) 
       enddo


     enddo
   enddo
 enddo

 do lbf=1,nbf_eri
   do kbf=1,nbf_eri

     do jstate=1,nbf_eri
       eri_eigenstate_i(jstate,kbf,lbf,nspin) = DOT_PRODUCT( eri_tmp3(:,kbf,lbf) , c_matrix(:,jstate,ijspin) )
     enddo

   enddo
 enddo


  
 do klspin=1,nspin

   do lbf=1,nbf_eri
     do kstate=1,nbf_eri
       do jstate=1,nbf_eri
         eri_tmp3(jstate,kstate,lbf) = DOT_PRODUCT( eri_eigenstate_i(jstate,:,lbf,nspin) , c_matrix(:,kstate,klspin) )
       enddo
     enddo
   enddo

   do lstate=1,nbf_eri
     do kstate=1,nbf_eri
       do jstate=1,nbf_eri

         eri_eigenstate_i(jstate,kstate,lstate,klspin) = DOT_PRODUCT( eri_tmp3(jstate,kstate,:) , c_matrix(:,lstate,klspin) )

       enddo
     enddo
   enddo

 enddo !klspin

 call stop_clock(timing_basis_transform)

end subroutine transform_eri_basis


!=================================================================
subroutine prepare_eri_3center_eigen(c_matrix)
 use m_inputparam,only: nspin
 implicit none
 real(dp),intent(in)  :: c_matrix(nbf_eri,nbf_eri,nspin)
!=====
 integer              :: kbf,lbf
 integer              :: kstate,lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp(:,:,:)
 integer              :: ipair
!=====

 call start_clock(timing_eri_3center_eigen)

 write(stdout,'(/,a)') ' Calculate 3-center integrals on eigenstates'


 !TODO merge the 2 last indexes for prod_basis save a factor 2! (i<->j symmetry)
 call clean_allocate('3-center MO integrals',eri_3center_eigen,nauxil_3center,nbf_eri,nbf_eri,nspin)

 allocate(eri_3center_tmp(nauxil_3center,nbf_eri,nbf_eri)) 
 eri_3center_eigen(:,:,:,:) = 0.0_dp
 do klspin=1,nspin
   ! Transformation of the first index
   eri_3center_tmp(:,:,:) = 0.0_dp
   do lstate=1,nbf_eri
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       eri_3center_tmp(:,kbf,lstate) = eri_3center_tmp(:,kbf,lstate) &
                                       + c_matrix(lbf,lstate,klspin) * eri_3center(:,ipair)
       if( kbf /= lbf )  &
         eri_3center_tmp(:,lbf,lstate) = eri_3center_tmp(:,lbf,lstate) &
                                         + c_matrix(kbf,lstate,klspin) * eri_3center(:,ipair)

     enddo
   enddo

   ! Transformation of the second index
   do lstate=1,nbf_eri
     eri_3center_eigen(:,:,lstate,klspin) = MATMUL( eri_3center_tmp(:,:,lstate) , c_matrix(:,:,klspin) )
   enddo

 enddo ! klspin
 deallocate(eri_3center_tmp)

 call stop_clock(timing_eri_3center_eigen)

end subroutine prepare_eri_3center_eigen


!=================================================================
subroutine prepare_eri_3center_eigen_mixed(c_matrix)
 use m_inputparam,only: nspin
 implicit none
 real(dp),intent(in)  :: c_matrix(nbf_eri,nbf_eri,nspin)
!=====
 integer              :: kbf,lbf
 integer              :: kstate,lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp(:,:,:)
 real(dp),allocatable :: c_matrix_exx(:,:,:)
 logical              :: file_exists
!=====

 call start_clock(timing_eri_3center_eigen)

 inquire(file='fort.1000',exist=file_exists)
 if( .NOT. file_exists ) call die('fort.1000 not found')

 allocate(c_matrix_exx(nbf_eri,nbf_eri,nspin))
 open(1000,form='unformatted')
 do klspin=1,nspin
   do lstate=1,nbf_eri
     read(1000) c_matrix_exx(:,lstate,klspin)
   enddo
 enddo
 close(1000,status='delete')


 write(stdout,'(/,a)') ' Calculate 3-center integrals on MIXED eigenstates'


 !TODO merge the 2 last indexes for prod_basis save a factor 2! (i<->j symmetry)
 call clean_allocate('3-center MO integrals',eri_3center_eigen_mixed,nauxil_3center,nbf_eri,nbf_eri,nspin)

 allocate(eri_3center_tmp(nauxil_3center,nbf_eri,nbf_eri)) 
 eri_3center_eigen_mixed(:,:,:,:) = 0.0_dp
 do klspin=1,nspin
   ! Transformation of the first index
   eri_3center_tmp(:,:,:) = 0.0_dp
   do kbf=1,nbf_eri
     do lbf=1,nbf_eri
       if( negligible_basispair(kbf,lbf) ) cycle

         do lstate=1,nbf_eri
           eri_3center_tmp(:,kbf,lstate) = eri_3center_tmp(:,kbf,lstate) &
                                      + c_matrix_exx(lbf,lstate,klspin) * eri_3center(:,index_prod(kbf,lbf))
         enddo

     enddo
   enddo
   ! Transformation of the second index
   do lstate=1,nbf_eri
     eri_3center_eigen_mixed(:,:,lstate,klspin) = MATMUL( eri_3center_tmp(:,:,lstate) , c_matrix(:,:,klspin) )
   enddo

 enddo ! klspin
 deallocate(eri_3center_tmp)

 call stop_clock(timing_eri_3center_eigen)

end subroutine prepare_eri_3center_eigen_mixed


!=================================================================
subroutine destroy_eri_3center_eigen()
 implicit none
!=====

 write(stdout,'(/,a)') ' Destroy 3-center integrals on eigenstates'
 call clean_deallocate('3-center MO integrals',eri_3center_eigen)

end subroutine destroy_eri_3center_eigen


!=================================================================
subroutine destroy_eri_3center()
 implicit none
!=====

 if(ALLOCATED(eri_3center)) then
   call clean_deallocate('3-center integrals',eri_3center)
 endif

end subroutine destroy_eri_3center


!=================================================================
subroutine destroy_eri_3center_lr()
 implicit none
!=====

 if(ALLOCATED(eri_3center_lr)) then
   call clean_deallocate('3-center LR integrals',eri_3center_lr)
 endif

end subroutine destroy_eri_3center_lr


!=========================================================================
subroutine negligible_eri(tol)
 implicit none
 real(dp),intent(in) :: tol
!=====
 integer             :: icount,ibf,jbf,kbf,lbf,jcount
 integer             :: ibuffer
 real(dp)            :: integral_ij(nbf_eri,nbf_eri)
!=====

 icount=0
 do ibuffer=1,nsize
   if( ABS( eri_4center(ibuffer) ) < tol ) icount=icount+1
 enddo

 write(stdout,*) ' number of negligible integrals <',tol
 write(stdout,*) icount, ' / ',nsize,REAL(icount,dp)/REAL(nsize,dp)*100.0_dp,' [%]'


 do ibf=1,nbf_eri
   do jbf=1,nbf_eri
     integral_ij(ibf,jbf) = eri(ibf,jbf,ibf,jbf)
   enddo
 enddo

 write(stdout,*) 'testing Cauchy-Schwarz condition'
 icount=0
 jcount=0
 do ibf=1,nbf_eri
   do jbf=1,nbf_eri
     do kbf=1,nbf_eri
       do lbf=1,nbf_eri
         if( SQRT( integral_ij(ibf,jbf) * integral_ij(kbf,lbf) ) < tol ) icount = icount + 1
         if( ABS( eri(ibf,jbf,kbf,lbf) ) < tol ) jcount = jcount + 1
       enddo
     enddo
   enddo
 enddo
 write(stdout,*) ' number of negligible integrals <',tol
 write(stdout,*) icount, ' / ',nbf_eri**4,REAL(icount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'
 write(stdout,*) jcount, ' / ',nbf_eri**4,REAL(jcount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'


end subroutine negligible_eri


!=========================================================================
subroutine dump_out_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!=====
 character(len=50) :: filename
 integer           :: nline,iline,icurrent
 integer           :: erifile
!=====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 write(stdout,*) 'Dump out the ERI into file'
 write(stdout,*) 'Size of file [bytes]',REAL(nsize,dp)*prec_eri

 if( is_iomaster ) then
   open(newunit=erifile,file=TRIM(filename),form='unformatted')
   write(erifile) nsize
   write(erifile) rcut

   nline = nsize / line_length + 1
   icurrent=0
   do iline=1,nline
     write(erifile) eri_4center(icurrent+1:MIN(nsize,icurrent+line_length+1))
     icurrent = icurrent + line_length + 1
   enddo

   close(erifile)
 endif

 write(stdout,'(a,/)') ' file written'

end subroutine dump_out_eri


!=========================================================================
logical function read_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!=====
 character(len=50) :: filename
 integer           :: nline,iline,icurrent
 integer           :: integer_read
 real(dp)          :: real_read
 integer           :: erifile
!=====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 
 inquire(file=TRIM(filename),exist=read_eri)

 if(read_eri) then

   write(stdout,*) 'Try to read ERI file'
   open(newunit=erifile,file=TRIM(filename),form='unformatted',status='old')
   read(erifile) integer_read
   if(integer_read /= nsize) read_eri=.FALSE.
   read(erifile) real_read
   if(ABS(real_read-rcut) > 1.0d-6) read_eri=.FALSE.

   if(read_eri) then

     nline = nsize / line_length + 1
     icurrent=0
     do iline=1,nline
       read(erifile) eri_4center(icurrent+1:MIN(nsize,icurrent+line_length+1))
       icurrent = icurrent + line_length + 1
     enddo
     write(stdout,'(a,/)') ' ERI file read'

   else
     write(stdout,'(a,/)') ' reading aborted'
   endif

   close(erifile)

 endif


end function read_eri


!=========================================================================
end module m_eri
