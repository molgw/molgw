!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the method to prepare and store the 2-, 3-, and 4-center Coulomb integrals
!
!=========================================================================
module m_eri
 use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
 use m_definitions
 use m_mpi
 use m_memory
 use m_basis_set
 use m_timing


 real(dp),parameter,public :: TOO_LOW_EIGENVAL=1.0e-6_dp

 !
 ! max length of a record in the ERI file
 integer,parameter,private :: line_length=1000

 real(dp),private           :: TOL_INT

 real(prec_eri),public,allocatable :: eri_4center(:)
 real(prec_eri),public,allocatable :: eri_4center_lr(:)
 real(prec_eri),public,allocatable :: eri_3center(:,:)
 real(prec_eri),public,allocatable :: eri_3center_lr(:,:)

 logical,protected,allocatable      :: negligible_shellpair(:,:)
 integer,protected,allocatable      :: index_pair_1d(:)
 integer,protected,allocatable      :: index_basis(:,:)
 integer,protected,allocatable      :: index_shellpair(:,:)
 integer,protected                  :: nshellpair

 type shell_type
   integer              :: am
   integer              :: ng
   real(dp),allocatable :: alpha(:)
   real(dp),allocatable :: coeff(:)
   real(dp)             :: x0(3)
   integer              :: istart,iend   ! index of the shell's basis functions in the basis set
 end type shell_type

 integer,protected                      :: nshell
 integer,protected                      :: nshell_auxil
 type(shell_type),protected,allocatable :: shell(:)
 type(shell_type),protected,allocatable :: shell_auxil(:)
 integer,private,allocatable            :: shell_bf(:)
! integer,private,allocatable            :: index_in_shell_bf(:)


 integer,private   :: nbf_eri         ! local copy of nbf
 integer,protected :: nsize           ! size of the eri_4center array
 integer,protected :: npair         ! number of independent pairs (i,j) with i<=j 

 integer,protected :: nauxil_3center     ! size of the 3-center matrix
                                         ! may differ from the total number of 3-center integrals due to
                                         ! data distribution
 integer,protected :: nauxil_3center_lr  ! size of the 3-center matrix
                                         ! may differ from the total number of 3-center integrals due to
                                         ! data distribution

 real(dp),allocatable,public  :: eri_3center_sca(:,:)

! Parallelization information for the auxiliary basis
 integer,allocatable,protected :: iproc_ibf_auxil(:)
 integer,allocatable,protected :: ibf_auxil_g(:)       ! auxil bf index from local to global
 integer,allocatable,protected :: ibf_auxil_l(:)       ! auxil bf index from global to local
 integer,allocatable,protected :: nbf_local_iproc(:)

! Parallelization information for the auxiliary basis (LR part)
 integer,allocatable,protected :: iproc_ibf_auxil_lr(:)
 integer,allocatable,protected :: ibf_auxil_g_lr(:)
 integer,allocatable,protected :: ibf_auxil_l_lr(:)
 integer,allocatable,protected :: nbf_local_iproc_lr(:)


 interface
   function eval_contr_integral(am1,am2,am3,am4, &
                                ng1,ng2,ng3,ng4, &
                                coeff1,coeff2,coeff3,coeff4,&
                                alpha1,alpha2,alpha3,alpha4,&
                                x01,x02,x03,x04,&
                                rcut, &
                                int_shell) bind(C,name='eval_contr_integral')
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT) :: eval_contr_integral
     integer(C_INT) :: am1,am2,am3,am4
     integer(C_INT) :: ng1,ng2,ng3,ng4
     real(C_DOUBLE) :: coeff1(1),coeff2(1),coeff3(1),coeff4(1)
     real(C_DOUBLE) :: alpha1(1),alpha2(1),alpha3(1),alpha4(1)
     real(C_DOUBLE) :: x01(1),x02(1),x03(1),x04(1)
     real(C_DOUBLE) :: rcut
     real(C_DOUBLE) :: int_shell(1)
   end function eval_contr_integral
 end interface


contains


!=========================================================================
subroutine prepare_eri(basis)
 use m_inputparam,only: integral_level
 implicit none
!===== 
 type(basis_set),intent(in) :: basis
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
   call identify_negligible_shellpair(basis)
   call setup_shellpair()
   call setup_basispair()
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
 if(ALLOCATED(negligible_shellpair))   deallocate(negligible_shellpair)
 if(ALLOCATED(index_shellpair))        deallocate(index_shellpair)
 if(ALLOCATED(index_pair_1d)) call clean_deallocate('index pair',index_pair_1d)
 if(ALLOCATED(index_basis))   call clean_deallocate('index basis',index_basis)

 ! 
 ! Cleanly deallocate the shell objects
 do ishell=1,nshell
   if(ALLOCATED(shell(ishell)%alpha)) deallocate( shell(ishell)%alpha )
   if(ALLOCATED(shell(ishell)%coeff)) deallocate( shell(ishell)%coeff )
 enddo
 if(ALLOCATED(shell))                 deallocate(shell)
 if(ALLOCATED(shell_bf))              deallocate(shell_bf)
! if(ALLOCATED(index_in_shell_bf))     deallocate(index_in_shell_bf)


end subroutine deallocate_eri


!=========================================================================
subroutine deallocate_index_pair()
 implicit none

!=====
 call clean_deallocate('index pair',index_pair_1d)

end subroutine deallocate_index_pair


!=========================================================================
function index_eri(ibf,jbf,kbf,lbf)
 implicit none

 integer,intent(in) :: ibf,jbf,kbf,lbf
 integer            :: index_eri
!=====
 integer            :: klmin,ijmax
 integer            :: index_ij,index_kl
!===== 

 index_ij = index_pair(ibf,jbf)
 index_kl = index_pair(kbf,lbf)

 ijmax=MAX(index_ij,index_kl)
 klmin=MIN(index_ij,index_kl)

 index_eri = (klmin-1)*npair - (klmin-1)*(klmin-2)/2 + ijmax-klmin+1


end function index_eri


!=========================================================================
function index_pair(ibf,jbf)
 implicit none

 integer,intent(in) :: ibf,jbf
 integer            :: index_pair
!=====
 integer            :: ijmin,ijmax
!=====

 if( ibf == jbf ) then
   index_pair = ibf
 else
   ijmax=MAX(ibf,jbf)
   ijmin=MIN(ibf,jbf)

   index_pair = (ijmin-1) * (nbf_eri-1) - (ijmin-1) * (ijmin-2)/2     + ijmax - ijmin + nbf_eri

   index_pair = index_pair_1d(index_pair)
 endif


end function index_pair


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
   index_ij = index_pair(ibf,jbf)
   index_kl = index_pair(kbf,lbf)
  
   eri_ri = DOT_PRODUCT( eri_3center(:,index_ij) , eri_3center(:,index_kl) )

   call xsum_auxil(eri_ri)

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
   index_ij = index_pair(ibf,jbf)
   index_kl = index_pair(kbf,lbf)

   eri_ri_lr = DOT_PRODUCT( eri_3center_lr(:,index_ij) , eri_3center_lr(:,index_kl) )

   call xsum_auxil(eri_ri_lr)

 endif

end function eri_ri_lr


!=========================================================================
subroutine setup_shell_list(basis)
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer :: ibf,jbf,kbf
 integer :: ishell
!=====


 nshell = basis%nshell
 allocate(shell(nshell))

 !
 ! Set up shells information
 jbf=0
 do ishell=1,nshell
   do ibf=1,basis%nbf_cart
     if( basis%bf(ibf)%shell_index == ishell ) then
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

 !
 ! Set up the correspondence between basis function and shells (the inverse of
 ! the previous table more or less)
 !
 allocate(shell_bf(basis%nbf))
 ibf=1
 jbf=1
 do while (ibf<=basis%nbf_cart)
   kbf = jbf + number_basis_function_am( basis%gaussian_type , basis%bf(ibf)%am ) - 1
   shell_bf(jbf:kbf) = basis%bf(ibf)%shell_index
   jbf = kbf + 1
   ibf = ibf + number_basis_function_am( 'CART' , basis%bf(ibf)%am )
 enddo
! allocate(index_in_shell_bf(basis%nbf))
! do ibf=1,basis%nbf
!   index_in_shell_bf(ibf) = basis%bff(ibf)%index_in_shell
! enddo


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
subroutine setup_basispair()
 implicit none
!=====
 integer :: ishell,jshell
 integer :: ibf,jbf,ijbf
!=====

 npair = 0
 do jbf=1,nbf_eri
   do ibf=1,jbf
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
     endif
   enddo
 enddo

 call clean_allocate('index pair',index_pair_1d,(nbf_eri*(nbf_eri+1))/2)
 call clean_allocate('index basis',index_basis,2,npair)

 !
 ! Specific ordering where the first nbf pairs contain the diagonal terms ibf==jbf
 !
 npair = 0
 index_pair_1d(:) = 0
 do jbf=1,nbf_eri
   if( negligible_basispair(jbf,jbf) ) then
     call die('setup_negligible_basispair: this should not happen')
   endif
   npair = npair + 1
   index_pair_1d(jbf) = npair
   index_pair_1d(jbf) = npair
   index_basis(1,npair) = jbf
   index_basis(2,npair) = jbf
 enddo

 ijbf = nbf_eri

 do ibf=1,nbf_eri 
   do jbf=ibf+1,nbf_eri  ! Skip the diagonal terms since it is already included 
     ijbf = ijbf + 1
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
       index_pair_1d(ijbf) = npair
       index_basis(1,npair) = ibf
       index_basis(2,npair) = jbf
     endif
   enddo
 enddo


end subroutine setup_basispair


!=========================================================================
function negligible_basispair(ibf,jbf)
 implicit none

 integer,intent(in) :: ibf,jbf
 logical :: negligible_basispair
!=====
 integer  :: ishell,jshell
!=====


 ishell = shell_bf(ibf)
 jshell = shell_bf(jbf)

 negligible_basispair = negligible_shellpair(ishell,jshell)


end function negligible_basispair


!=========================================================================
!
! Find negligible shell pairs with
! Cauchy-Schwarz inequality: (ij|1/r|kl)**2 <= (ij|1/r|ij) (kl|1/r|(kl) 
!
!=========================================================================
subroutine identify_negligible_shellpair(basis)
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer                      :: info,ip
 integer                      :: iibf
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: n1c,n2c
 integer                      :: ni,nj
 integer                      :: ami,amj
 integer                      :: ishell,jshell
 integer                      :: ig1,ig2,ig3,ig4
 real(dp)                     :: zeta_12,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
 real(dp)                     :: workload(nproc_world)
 integer                      :: shell_proc(nshell)
!=====
! variables used to call C
 integer(C_INT)               :: am1,am2
 integer(C_INT)               :: ng1,ng2
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:)
 real(C_DOUBLE)               :: x01(3),x02(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 call start_clock(timing_eri_screening)
 write(stdout,'(/,a)')    ' Cauchy-Schwartz screening of the 3- or 4-center integrals'

 !
 ! Load balancing
 workload(:) = 0.0_dp
 do jshell=1,nshell
   amj = shell(jshell)%am
   ip = MINLOC(workload(:),DIM=1)
   !
   ! Cost function was evaluated from a few runs
   workload(ip) = workload(ip) + cost_function_eri(amj)
   shell_proc(jshell) = ip - 1
 enddo


 negligible_shellpair(:,:) = .TRUE.

 do jshell=1,nshell
   !
   ! Workload is distributed here
   if( shell_proc(jshell) /= rank_world ) cycle

   amj = shell(jshell)%am
   nj  = number_basis_function_am( basis%gaussian_type , amj )
   n2c = number_basis_function_am( 'CART' , amj )
   am2 = shell(jshell)%am
   ng2 = shell(jshell)%ng

   do ishell=1,nshell
     ami = shell(ishell)%am
     !TODO: Here time could be saved by only checking ishell<= jshell
     ! But then an interexchange of indexes would have to be implemented in
     ! order to satisfy ami >= amj (required condition in libint)
     if( ami < amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     n1c = number_basis_function_am( 'CART' , ami )
     am1 = shell(ishell)%am
     ng1 = shell(ishell)%ng

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
               rho1 = zeta_12 * zeta_12 / ( zeta_12 + zeta_12 )
               rho  = rho1
               tt = rho * SUM( (p(:)-q(:))**2 )
               call boys_function_c(f0t(0),0,tt)

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
                               0.0_C_DOUBLE, &
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

 call xand_world(negligible_shellpair)

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
   do ishell=1,jshell 
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
   do ishell=1,jshell 
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
subroutine destroy_eri_3center()
 implicit none
!=====

 if(ALLOCATED(eri_3center)) then
   call clean_deallocate('3-center integrals',eri_3center)
 endif
#ifdef SCASCA
 if(ALLOCATED(eri_3center_sca)) then
   call clean_deallocate('3-center integrals SCALAPACK',eri_3center_sca)
 endif
#endif

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
   if(ABS(real_read-rcut) > 1.0e-6_dp) read_eri=.FALSE.

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
! Rough evaluation of the CPU time to get an ERI as a function of the 
! angular momentum
! 
!=========================================================================
function cost_function_eri(am)
 implicit none
 integer,intent(in)  :: am
 real(dp)            :: cost_function_eri
!=====

 cost_function_eri = am**2 + 4.6_dp 

end function cost_function_eri


!=========================================================================
subroutine distribute_auxil_basis(nbf_auxil_basis)
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
!=====

 if( parallel_buffer ) then
   ! Use nproc instead nproc_local
  
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
   allocate(nbf_local_iproc(0:nproc_auxil-1))
  
   iproc              = nproc_auxil - 1
   nbf_local_iproc(:) = 0
   do ibf=1,nbf_auxil_basis
  
     iproc = MODULO(iproc+1,nproc_auxil)
  
     iproc_ibf_auxil(ibf) = iproc
  
     nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1
  
   enddo
  
   nauxil_3center = nbf_local_iproc(rank_auxil)
  
   allocate(ibf_auxil_g(nauxil_3center))
   allocate(ibf_auxil_l(nbf_auxil_basis))
   ibf_auxil_l(:) = 0
   ibf_local = 0
   do ibf=1,nbf_auxil_basis
     if( rank_auxil == iproc_ibf_auxil(ibf) ) then
       ibf_local = ibf_local + 1
       ibf_auxil_g(ibf_local) = ibf
       ibf_auxil_l(ibf)       = ibf_local
     endif
   enddo
  
 else

   allocate(iproc_ibf_auxil(nbf_auxil_basis))
   allocate(nbf_local_iproc(0:nproc_local-1))
  
   iproc              = nproc_local-1
   nbf_local_iproc(:) = 0
   do ibf=1,nbf_auxil_basis
  
     iproc = MODULO(iproc+1,nproc_local)
  
     iproc_ibf_auxil(ibf) = iproc
  
     nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1
  
   enddo
  
   nauxil_3center = nbf_local_iproc(rank_local)
  
   allocate(ibf_auxil_g(nauxil_3center))
   allocate(ibf_auxil_l(nbf_auxil_basis))
   ibf_auxil_l(:) = 0
   ibf_local = 0
   do ibf=1,nbf_auxil_basis
     if( rank_local == iproc_ibf_auxil(ibf) ) then
       ibf_local = ibf_local + 1
       ibf_auxil_g(ibf_local) = ibf
       ibf_auxil_l(ibf)       = ibf_local
     endif
   enddo
  
 endif

 write(stdout,'(/,a)') ' Distribute auxiliary basis functions among processors'
 do iproc=0,0
   write(stdout,'(a,i4,a,i6,a)')   ' Proc: ',iproc,' treats ',nbf_local_iproc(iproc),' auxiliary basis functions'
 enddo



end subroutine distribute_auxil_basis


!=========================================================================
subroutine distribute_auxil_basis_lr(nbf_auxil_basis)
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
!=====

 allocate(iproc_ibf_auxil_lr(nbf_auxil_basis))
 allocate(nbf_local_iproc_lr(0:nproc_auxil-1))

 iproc = nproc_auxil - 1
 nbf_local_iproc_lr(:) = 0
 do ibf=1,nbf_auxil_basis

   iproc = MODULO(iproc+1,nproc_auxil)

   iproc_ibf_auxil_lr(ibf) = iproc

   nbf_local_iproc_lr(iproc) = nbf_local_iproc_lr(iproc) + 1

 enddo

 nauxil_3center_lr = nbf_local_iproc_lr(rank_auxil)

 allocate(ibf_auxil_g_lr(nauxil_3center_lr))
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 ibf_auxil_l_lr(:) = 0
 ibf_local = 0
 do ibf=1,nbf_auxil_basis
   if( rank_auxil == iproc_ibf_auxil_lr(ibf) ) then
     ibf_local = ibf_local + 1
     ibf_auxil_g_lr(ibf_local) = ibf
     ibf_auxil_l_lr(ibf)       = ibf_local
   endif
 enddo

 write(stdout,'(/,a)') ' Distribute LR auxiliary basis functions among processors'
 do iproc=0,0
   write(stdout,'(a,i4,a,i6,a)')   ' Proc: ',iproc,' treats ',nbf_local_iproc_lr(iproc),' auxiliary basis functions'
 enddo


end subroutine distribute_auxil_basis_lr


!=========================================================================
end module m_eri
