!=========================================================================
#include "macros.h"
!=========================================================================
module m_eri
 use m_definitions
 use m_mpi
 use m_basis_set
 use m_timing

 integer,parameter :: BUFFER1 = 1
 integer,parameter :: BUFFER2 = 2
 !
 ! max length of a record in the ERI file
 integer,parameter :: line_length=1000

 real(dp),protected                 :: TOL_INT=1.0e-10_dp

 real(prec_eri),private,allocatable :: eri_buffer(:)
 real(prec_eri),private,allocatable :: eri_buffer_lr(:)
 real(prec_eri),private,allocatable :: eri_2center_m1(:,:)
 ! eri_3center_****  are only "protected" since you may need it outside for computational tricks
 real(prec_eri),protected,allocatable :: eri_3center(:,:)
 real(prec_eri),protected,allocatable :: eri_3center_eigen(:,:,:,:)

 logical,protected,allocatable      :: negligible_basispair(:,:)
 logical,private,allocatable        :: negligible_shellpair(:,:)
 integer,private,allocatable        :: index_pair(:,:)
 integer,private,allocatable        :: index_shellpair(:,:)
 integer,private                    :: nshellpair

 type shell_type
   integer              :: am
   integer              :: ng
   real(dp),allocatable :: alpha(:)
   real(dp),allocatable :: coeff(:)
   real(dp)             :: x0(3)
   integer              :: istart,iend
 end type shell_type
 integer,private                      :: nshell
 integer,private                      :: nshell_auxil
 type(shell_type),private,allocatable :: shell(:)
 type(shell_type),private,allocatable :: shell_auxil(:)


 integer,private              :: nbf_eri                ! local copy of nbf
 integer,private              :: nsize                  ! size of the eri_buffer array
 integer,private              :: nsize1                 ! number of independent pairs (i,j) with i<=j

 integer,private              :: nbf_eri_auxil          ! local copy of nbf for auxiliary basis
 integer,private              :: nsize_auxil            ! size of the eri_buffer array
 integer,private              :: nsize1_auxil           ! number of independent pairs (i,j) with i<=j


! interface
!   integet(C_INT) function eval_contr_integral() bind(C)
!       info=eval_contr_integral(                &
!                               am1,am2,am3,am4, &
!                               ng1,ng2,ng3,ng4, &
!                               coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
!                               alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
!                               x01(1),x02(1),x03(1),x04(1),&
!                               rcut,
!                               int_shell(1)
!     use ISO_C_BINDING
!     character(kind=c_char) :: string(*)
!   end subroutine print_c
! end interface


contains

!=========================================================================
subroutine prepare_eri(basis,rcut,which_buffer)
 implicit none
!===== 
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: rcut
 integer,intent(in)         :: which_buffer
!===== 
 integer            :: info
 logical            :: file_exists
!===== 

 nbf_eri = basis%nbf
 inquire(file='manual_tol_int',exist=file_exists)
 if( file_exists ) then
   open(unit=22,file='manual_tol_int',status='old')
   read(22,*) TOL_INT
   close(22)
   TOL_INT = MAX(TOL_INT,0.0_dp)
   WRITE_MASTER(msg,'(a,x,es14.4)') 'TOL_INT manually set to',TOL_INT
   call issue_warning(msg)
 else
   TOL_INT = 1.0e-10_dp
 endif

 if(.NOT.allocated(negligible_shellpair)) then
   call setup_shell_list(basis)
   allocate(negligible_shellpair(nshell,nshell))
   allocate(negligible_basispair(nbf_eri,nbf_eri))
   allocate(index_pair(nbf_eri,nbf_eri))
   call identify_negligible_shellpair(basis)
   call setup_shellpair()
   call setup_negligible_basispair()
 endif


 nsize = (nsize1*(nsize1+1))/2


end subroutine prepare_eri


!=========================================================================
subroutine allocate_eri_auxil(auxil_basis)
 implicit none
!===== 
 type(basis_set),intent(in) :: auxil_basis
!===== 
 integer            :: info
 logical            :: file_exists
!===== 

 nbf_eri_auxil = auxil_basis%nbf
 inquire(file='manual_tol_int',exist=file_exists)
 if( file_exists ) then
   open(unit=22,file='manual_tol_int',status='old')
   read(22,*) TOL_INT
   close(22)
   TOL_INT = MAX(TOL_INT,0.0_dp)
   WRITE_MASTER(msg,'(a,x,es14.4)') 'TOL_INT manually set to',TOL_INT
   call issue_warning(msg)
 else
   TOL_INT = 1.0e-10_dp
 endif

 call setup_shell_list_auxil(auxil_basis)

 nsize1_auxil = nbf_eri_auxil 
 nsize_auxil  = nsize1_auxil**2

 if(nsize_auxil<1) stop'too many or too few integrals to be stored'

 !
 ! 2-CENTER INTEGRALS 
 !
 WRITE_MASTER(*,*) 'Allocate 2-center integrals'
 call memory_statement(REAL(nsize_auxil,dp)*REAL(prec_eri/dp,dp))

 allocate(eri_2center_m1(nsize1_auxil,nsize1_auxil),stat=info)

 if(info==0) then
   WRITE_MASTER(*,*) 'success'
 else
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 !
 ! 3-CENTER INTEGRALS 
 !
 WRITE_MASTER(*,*) 'Allocate 3-center integrals'
 call memory_statement(REAL(nsize1_auxil,dp)*REAL(nsize1)*REAL(prec_eri/dp,dp))

 allocate(eri_3center(nsize1_auxil,nsize1),stat=info)
 eri_3center(:,:) = 0.0_dp

 if(info==0) then
   WRITE_MASTER(*,*) 'success'
 else
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


end subroutine allocate_eri_auxil


!=========================================================================
subroutine deallocate_eri_buffer()
 implicit none
!=====

 if(allocated(eri_buffer)) then
   WRITE_MASTER(*,'(/,a)')     ' Deallocate ERI buffer'
   call memory_statement(-REAL(nsize,dp)*REAL(prec_eri,dp)/REAL(dp,dp))
   deallocate(eri_buffer)
 endif
 if(allocated(eri_buffer_lr)) then
   WRITE_MASTER(*,'(/,a)')     ' Deallocate LR ERI buffer'
   call memory_statement(-REAL(nsize,dp)*REAL(prec_eri,dp)/REAL(dp,dp))
   deallocate(eri_buffer_lr)
 endif
 WRITE_MASTER(*,*)

end subroutine deallocate_eri_buffer


!=========================================================================
subroutine deallocate_eri()
 implicit none

 integer :: ishell
!=====

 if(allocated(eri_buffer))            deallocate(eri_buffer)
 if(allocated(eri_buffer_lr))         deallocate(eri_buffer_lr)
 if(allocated(eri_2center_m1))        deallocate(eri_2center_m1)
 if(allocated(eri_3center))           deallocate(eri_3center)
 if(allocated(negligible_basispair))  deallocate(negligible_basispair)
 if(allocated(negligible_shellpair))  deallocate(negligible_shellpair)
 if(allocated(index_pair))            deallocate(index_pair)
 if(allocated(index_shellpair))       deallocate(index_shellpair)
 ! 
 ! Cleanly deallocate the shell objects
 do ishell=1,nshell
   if(allocated(shell(ishell)%alpha)) deallocate( shell(ishell)%alpha )
   if(allocated(shell(ishell)%coeff)) deallocate( shell(ishell)%coeff )
 enddo
 if(allocated(shell))                 deallocate(shell)


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

 index_eri = (klmin-1)*nsize1 - (klmin-1)*(klmin-2)/2 + ijmax-klmin+1

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
   eri = eri_buffer(index_eri(ibf,jbf,kbf,lbf))
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
   eri_lr = eri_buffer_lr(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri_lr


!=========================================================================
function eri_ri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_ri
!=====
 integer            :: index_ij,index_kl
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_ri = 0.0_dp
 else
   index_ij = index_prod(ibf,jbf)
   index_kl = index_prod(kbf,lbf)

!     eri_ri = DOT_PRODUCT( eri_3center(:,index_ij) , MATMUL( eri_2center_m1(:,:) , eri_3center(:,index_kl) ) ) 
   eri_ri = DOT_PRODUCT( eri_3center(:,index_ij) , eri_3center(:,index_kl) )

 endif

end function eri_ri


!=========================================================================
function eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen_ri
!=====

 eri_eigen_ri = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

end function eri_eigen_ri


!=========================================================================
subroutine calculate_eri(print_eri,basis,rcut,which_buffer)
 implicit none
 logical,intent(in)           :: print_eri
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: which_buffer
!=====

 call start_clock(timing_eri_4center)

 if( .NOT. read_eri(rcut) ) call do_calculate_eri(basis,rcut,which_buffer)


 if( print_eri ) then
   call dump_out_eri(rcut)
 endif

 call stop_clock(timing_eri_4center)

end subroutine calculate_eri


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
subroutine do_calculate_eri(basis,rcut,which_buffer)
 use ISO_C_BINDING
 use m_tools,only: boys_function
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: which_buffer
!=====
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: ijshellpair,klshellpair
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ig1,ig2,ig3,ig4
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ii,i,j,k,l
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: iibf,jjbf,kkbf,llbf
 integer                      :: info
 integer                      :: ordering
 real(dp)                     :: zeta_12,zeta_34,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
!=====
! variables used to call C
 integer(C_INT),external      :: libint_init,calculate_integral
 integer(C_INT),external      :: eval_contr_integral
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: rcut_libint
!=====

 WRITE_MASTER(*,'(/,a)') ' Calculate and store all the Electron Repulsion Integrals (ERI)'
 WRITE_MASTER(*,'(a)')      ' Libint library initialized'
 WRITE_MASTER(*,'(a,i5,/)') ' Max angular momentum handled by your Libint compilation: ',libint_init()

 WRITE_MASTER(*,*) 
 WRITE_MASTER(*,*) 'Number of integrals to be stored:',nsize
 WRITE_MASTER(*,*) 'Max index size',HUGE(nsize)
 if(nsize<1) stop'too many integrals to be stored'

 WRITE_MASTER(*,*) 'Allocate 4-center integrals'
 call memory_statement(REAL(nsize,dp)*REAL(prec_eri,dp)/REAL(dp,dp))

 select case(which_buffer)
 case(BUFFER1)
   allocate(eri_buffer(nsize),stat=info)
   eri_buffer(:) = 0.0_dp
 case(BUFFER2)
   allocate(eri_buffer_lr(nsize),stat=info)
   eri_buffer_lr(:) = 0.0_dp
 end select

 if(info==0) then
   WRITE_MASTER(*,*) 'success'
 else
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif



 rcut_libint = rcut

 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   !
   ! Order the angular momenta so that libint is pleased
   ! 1) am3+am4 >= am1+am2
   ! 2) am3>=am4
   ! 3) am1>=am2
   amk = shell(kshell)%am
   aml = shell(lshell)%am


   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)

     ami = shell(ishell)%am
     amj = shell(jshell)%am
     if( amk+aml < ami+amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     nk = number_basis_function_am( basis%gaussian_type , amk )
     nl = number_basis_function_am( basis%gaussian_type , aml )


     am1 = shell(ishell)%am
     am2 = shell(jshell)%am
     am3 = shell(kshell)%am
     am4 = shell(lshell)%am
     n1c = number_basis_function_am( CARTESIAN , ami )
     n2c = number_basis_function_am( CARTESIAN , amj )
     n3c = number_basis_function_am( CARTESIAN , amk )
     n4c = number_basis_function_am( CARTESIAN , aml )
     ng1 = shell(ishell)%ng
     ng2 = shell(jshell)%ng
     ng3 = shell(kshell)%ng
     ng4 = shell(lshell)%ng
     allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
     alpha1(:) = shell(ishell)%alpha(:) 
     alpha2(:) = shell(jshell)%alpha(:)
     alpha3(:) = shell(kshell)%alpha(:)
     alpha4(:) = shell(lshell)%alpha(:)
     x01(:) = shell(ishell)%x0(:)
     x02(:) = shell(jshell)%x0(:)
     x03(:) = shell(kshell)%x0(:)
     x04(:) = shell(lshell)%x0(:)
     allocate(coeff1(shell(ishell)%ng))
     allocate(coeff2(shell(jshell)%ng))
     allocate(coeff3(shell(kshell)%ng))
     allocate(coeff4(shell(lshell)%ng))
     coeff1(:)=shell(ishell)%coeff(:)
     coeff2(:)=shell(jshell)%coeff(:)
     coeff3(:)=shell(kshell)%coeff(:)
     coeff4(:)=shell(lshell)%coeff(:)

     allocate( int_shell( n1c*n2c*n3c*n4c ) )
     allocate( integrals_cart(n1c,n2c,n3c,n4c) )
     allocate( integrals_tmp(n1c,n2c,n3c,n4c) )
     integrals_cart(:,:,:,:) = 0.0_dp


     if(am1+am2+am3+am4==0) then

       do ig4=1,ng4
         do ig3=1,ng3
           do ig2=1,ng2
             do ig1=1,ng1

               zeta_12 = alpha1(ig1) + alpha2(ig2)
               zeta_34 = alpha3(ig3) + alpha4(ig4)
               p(:) = ( alpha1(ig1) * x01(:) + alpha2(ig2) * x02(:) ) / zeta_12 
               q(:) = ( alpha3(ig3) * x03(:) + alpha4(ig4) * x04(:) ) / zeta_34 
               !
               ! Treat carefully the LR only integrals
               rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 + zeta_12*zeta_34*rcut**2 )
               rho1 = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
               tt = rho * SUM( (p(:)-q(:))**2 )
               call boys_function(f0t(0),0,tt)

               integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                     2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                     / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                     / zeta_34 * EXP( -alpha3(ig3)*alpha4(ig4)/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) &
                     * SQRT( rho / rho1 ) &
                     * coeff1(ig1) &
                     * coeff2(ig2) &
                     * coeff3(ig3) &
                     * coeff4(ig4) * cart_to_pure_norm(0)%matrix(1,1)**4

             enddo
           enddo
         enddo
       enddo

     else

        
       info=eval_contr_integral(                &
                               am1,am2,am3,am4, &
                               ng1,ng2,ng3,ng4, &
                               coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
                               alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
                               x01(1),x02(1),x03(1),x04(1),&
                               rcut_libint, &
                               int_shell(1))


       if(info/=0) then
         WRITE_MASTER(*,*) am1,am2,am3,am4
         stop 'ERI calculated by libint failed'
       endif

       iibf=0
       do ibf=1,n1c
         do jbf=1,n2c
           do kbf=1,n3c
             do lbf=1,n4c
               iibf=iibf+1
               integrals_cart(ibf,jbf,kbf,lbf) = int_shell(iibf)
             enddo
           enddo
         enddo
       enddo


       do lbf=1,n4c
         do kbf=1,n3c
           do jbf=1,n2c
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1c,jbf,kbf,lbf) * cart_to_pure_norm(shell(ishell)%am)%matrix(1:n1c,ibf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n4c
         do kbf=1,n3c
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2c,kbf,lbf) * cart_to_pure_norm(shell(jshell)%am)%matrix(1:n2c,jbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n4c
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3c,lbf) * cart_to_pure_norm(shell(kshell)%am)%matrix(1:n3c,kbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4c) * cart_to_pure_norm(shell(lshell)%am)%matrix(1:n4c,lbf) )
             enddo
           enddo
         enddo
       enddo

     endif ! is (ss|ss)
     
     do lbf=1,nl
       do kbf=1,nk
         do jbf=1,nj
           do ibf=1,ni
             if( which_buffer == BUFFER1 ) then
               eri_buffer( index_eri(shell(ishell)%istart+ibf-1, &
                                     shell(jshell)%istart+jbf-1, &
                                     shell(kshell)%istart+kbf-1, &
                                     shell(lshell)%istart+lbf-1) ) = integrals_cart(ibf,jbf,kbf,lbf)
             else
               eri_buffer_lr( index_eri(shell(ishell)%istart+ibf-1, &
                                        shell(jshell)%istart+jbf-1, &
                                        shell(kshell)%istart+kbf-1, &
                                        shell(lshell)%istart+lbf-1) ) = integrals_cart(ibf,jbf,kbf,lbf)
             endif
           enddo
         enddo
       enddo
     enddo


     deallocate(integrals_cart)
     deallocate(integrals_tmp)
     deallocate(int_shell)
     deallocate(alpha1,alpha2,alpha3,alpha4)
     deallocate(coeff1)
     deallocate(coeff2)
     deallocate(coeff3)
     deallocate(coeff4)


   enddo
 enddo


 WRITE_MASTER(*,'(a,/)') ' All ERI have been calculated'


end subroutine do_calculate_eri


!=========================================================================
subroutine calculate_eri_2center(print_eri,auxil_basis)
 use ISO_C_BINDING
 use m_tools,only: boys_function, invert
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none
 logical,intent(in)           :: print_eri
 type(basis_set),intent(in)   :: auxil_basis
!=====
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ng1,ng2,ng3,ng4
 integer                      :: ig1,ig2,ig3,ig4
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ii,i,j,k,l
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: iibf,jjbf,kkbf,llbf
 integer                      :: info
 integer                      :: ordering
 real(dp)                     :: zeta_12,zeta_34,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
 real(dp),allocatable         :: eigval(:)
!=====
! variables used to call C
 integer(C_INT),external      :: libint_init,calculate_integral
 integer(C_INT),external      :: eval_contr_integral
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: rcut_libint
!=====

 call start_clock(timing_eri_2center)

 WRITE_MASTER(*,'(/,a)')    ' Calculate, invert and store the 2-center Electron Repulsion Integrals'
 WRITE_MASTER(*,'(a)')      ' Libint library initialized'
 WRITE_MASTER(*,'(a,i5,/)') ' Max angular momentum handled by your Libint compilation: ',libint_init()

 rcut_libint = 0.0_dp

 do lshell=1,1  ! FAKE loop
   do kshell=1,nshell_auxil
     !
     ! Order the angular momenta so that libint is pleased
     ! 1) am3+am4 >= am1+am2
     ! 2) am3>=am4
     ! 3) am1>=am2
     amk = shell_auxil(kshell)%am
     aml = 0
     if( amk < aml ) cycle

     do jshell=1,1  ! FAKE loop
       do ishell=1,nshell_auxil
         ami = shell_auxil(ishell)%am
         amj = 0
         if( ami < amj ) cycle
         if( amk+aml < ami+amj ) cycle

         ni = number_basis_function_am( auxil_basis%gaussian_type , ami )
         nj = 1
         nk = number_basis_function_am( auxil_basis%gaussian_type , amk )
         nl = 1


         am1 = shell_auxil(ishell)%am
         am2 = 0
         am3 = shell_auxil(kshell)%am
         am4 = 0
         n1c = number_basis_function_am( CARTESIAN , ami )
         n2c = 1
         n3c = number_basis_function_am( CARTESIAN , amk )
         n4c = 1
         ng1 = shell_auxil(ishell)%ng
         ng2 = 1
         ng3 = shell_auxil(kshell)%ng
         ng4 = 1
         allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
         alpha1(:) = shell_auxil(ishell)%alpha(:) 
         alpha2(:) = 0.0_dp ! shell_auxil(jshell)%alpha(:)
         alpha3(:) = shell_auxil(kshell)%alpha(:)
         alpha4(:) = 0.0_dp ! shell_auxil(lshell)%alpha(:)
         x01(:) = shell_auxil(ishell)%x0(:)
         x02(:) = shell_auxil(ishell)%x0(:)
         x03(:) = shell_auxil(kshell)%x0(:)
         x04(:) = shell_auxil(kshell)%x0(:)
         allocate(coeff1(shell_auxil(ishell)%ng))
         allocate(coeff2(1))
         allocate(coeff3(shell_auxil(kshell)%ng))
         allocate(coeff4(1))
         coeff1(:)=shell_auxil(ishell)%coeff(:)
         coeff2(:)=1.0_dp
         coeff3(:)=shell_auxil(kshell)%coeff(:)
         coeff4(:)=1.0_dp

         allocate( int_shell( n1c*n2c*n3c*n4c ) )
         allocate( integrals_cart(n1c,n2c,n3c,n4c) )
         allocate( integrals_tmp(n1c,n2c,n3c,n4c) )
         integrals_cart(:,:,:,:) = 0.0_dp


         if(am1+am2+am3+am4==0) then

           do ig4=1,ng4
             do ig3=1,ng3
               do ig2=1,ng2
                 do ig1=1,ng1

                   zeta_12 = alpha1(ig1) + alpha2(ig2)
                   zeta_34 = alpha3(ig3) + alpha4(ig4)
                   p(:) = ( alpha1(ig1) * x01(:) + alpha2(ig2) * x02(:) ) / zeta_12 
                   q(:) = ( alpha3(ig3) * x03(:) + alpha4(ig4) * x04(:) ) / zeta_34 
                   !
                   ! Full range or long-range only integrals
                   rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
                   rho1 = rho
                   
                   tt = rho * SUM( (p(:)-q(:))**2 )
                   call boys_function(f0t(0),0,tt)

                   integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                         2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                         / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                         / zeta_34 * EXP( -alpha3(ig3)*alpha4(ig4)/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) &
                         * SQRT( rho / rho1 ) &
                         * coeff1(ig1)* coeff3(ig3) &
                         * cart_to_pure_norm(0)%matrix(1,1)**4

                 enddo
               enddo
             enddo
           enddo

         else


           info=eval_contr_integral(                &
                                   am1,am2,am3,am4, &
                                   ng1,ng2,ng3,ng4, &
                                   coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
                                   alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
                                   x01(1),x02(1),x03(1),x04(1),&
                                   rcut_libint, &
                                   int_shell(1))


           if(info/=0) then
             WRITE_MASTER(*,*) am1,am2,am3,am4
             stop 'ERI calculated by libint failed'
           endif

           iibf=0
           do ibf=1,n1c
             do jbf=1,n2c
               do kbf=1,n3c
                 do lbf=1,n4c
                   iibf=iibf+1
                   integrals_cart(ibf,jbf,kbf,lbf) = int_shell(iibf)
                 enddo
               enddo
             enddo
           enddo


           do lbf=1,n4c
             do kbf=1,n3c
               do jbf=1,n2c
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1c,jbf,kbf,lbf) * cart_to_pure_norm(shell_auxil(ishell)%am)%matrix(1:n1c,ibf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4c
             do kbf=1,n3c
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2c,kbf,lbf) * cart_to_pure_norm(shell_auxil(jshell)%am)%matrix(1:n2c,jbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4c
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3c,lbf) * cart_to_pure_norm(shell_auxil(kshell)%am)%matrix(1:n3c,kbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,nl
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4c) * cart_to_pure_norm(shell_auxil(lshell)%am)%matrix(1:n4c,lbf) )
                 enddo
               enddo
             enddo
           enddo


         endif
         

         do lbf=1,nl
           do kbf=1,nk
             do jbf=1,nj
               do ibf=1,ni
                 eri_2center_m1( shell_auxil(ishell)%istart+ibf-1,    &
                                 shell_auxil(kshell)%istart+kbf-1 )    = integrals_cart(ibf,jbf,kbf,lbf)
                 ! And the symmetric too
                 eri_2center_m1( shell_auxil(kshell)%istart+kbf-1,    &
                                 shell_auxil(ishell)%istart+ibf-1 )    = integrals_cart(ibf,jbf,kbf,lbf)
               enddo
             enddo
           enddo
         enddo

         deallocate(integrals_cart)
         deallocate(integrals_tmp)
         deallocate(int_shell)
         deallocate(alpha1,alpha2,alpha3,alpha4)
         deallocate(coeff1,coeff2,coeff3,coeff4)

       enddo
     enddo
   enddo
 enddo

 if( .FALSE. ) then
   !
   ! Perform in-place inversion here
   call invert(nsize1_auxil,eri_2center_m1)

   WRITE_MASTER(*,'(a)') ' All 2-center integrals have been calculated, inverted and stored'

 else

   allocate(eigval(nsize1_auxil))
   !
   ! Perform in-place diagonalization here
   call diagonalize(nsize1_auxil,eri_2center_m1,eigval)
   do jbf=1,nbf_eri_auxil
     eri_2center_m1(:,jbf) = eri_2center_m1(:,jbf) / SQRT( eigval(jbf) )
   enddo
   deallocate(eigval)

!   ! transform it back to the inverse
!   eri_2center_m1 = MATMUL( eri_2center_m1 , TRANSPOSE(eri_2center_m1) )

   WRITE_MASTER(*,'(a)') ' All 2-center integrals have been calculated, diagonalized and stored'

 endif

 call stop_clock(timing_eri_2center)

end subroutine calculate_eri_2center


!=========================================================================
subroutine calculate_eri_3center(print_eri,basis,auxil_basis)
 use ISO_C_BINDING
 use m_tools,only: boys_function
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none
 logical,intent(in)           :: print_eri
 type(basis_set),intent(in)   :: basis
 type(basis_set),intent(in)   :: auxil_basis
!=====
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: klshellpair
 integer                      :: n1,n2,n3,n4
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ng1,ng2,ng3,ng4
 integer                      :: ig1,ig2,ig3,ig4
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ii,i,j,k,l
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: iibf,jjbf,kkbf,llbf
 integer                      :: info
 integer                      :: ordering
 real(dp)                     :: zeta_12,zeta_34,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
!=====
! variables used to call C
 integer(C_INT),external      :: libint_init,calculate_integral
 integer(C_INT),external      :: eval_contr_integral
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: rcut_libint
!=====

 call start_clock(timing_eri_3center)

 WRITE_MASTER(*,'(/,a)') ' Calculate and store all the 3-center Electron Repulsion Integrals'
 WRITE_MASTER(*,'(a)')      ' Libint library initialized'
 WRITE_MASTER(*,'(a,i5,/)') ' Max angular momentum handled by your Libint compilation: ',libint_init()

 rcut_libint = 0.0_dp

 do klshellpair=1,nshellpair
     kshell = index_shellpair(1,klshellpair)
     lshell = index_shellpair(2,klshellpair)
     !
     ! Order the angular momenta so that libint is pleased
     ! 1) am3+am4 >= am1+am2
     ! 2) am3>=am4
     ! 3) am1>=am2
     amk = shell(kshell)%am
     aml = shell(lshell)%am
!     if( amk < aml ) cycle
!     if( amk < aml ) stop'SSHOULD NOT HAPPEN'

     do jshell=1,1  ! FAKE LOOP
       do ishell=1,nshell_auxil
         ami = shell_auxil(ishell)%am
         amj = 0
         if( ami < amj ) stop'PROBLEM'

         ni = number_basis_function_am( auxil_basis%gaussian_type , ami )
         nj = 1
         nk = number_basis_function_am( basis%gaussian_type , amk )
         nl = number_basis_function_am( basis%gaussian_type , aml )


         if( amk+aml >= ami+amj ) then

           am1 = shell_auxil(ishell)%am
           am2 = 0
           am3 = shell(kshell)%am
           am4 = shell(lshell)%am
           n1c = number_basis_function_am( CARTESIAN , ami )
           n2c = 1
           n3c = number_basis_function_am( CARTESIAN , amk )
           n4c = number_basis_function_am( CARTESIAN , aml )
           n1 = ni
           n2 = nj
           n3 = nk
           n4 = nl
           ng1 = shell_auxil(ishell)%ng
           ng2 = 1
           ng3 = shell(kshell)%ng
           ng4 = shell(lshell)%ng
           allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
           allocate(coeff1(ng1),coeff2(ng2),coeff3(ng3),coeff4(ng4))
           alpha1(:) = shell_auxil(ishell)%alpha(:) 
           alpha2(:) = 0.0_dp ! shell_auxil(jshell)%alpha(:)
           alpha3(:) = shell(kshell)%alpha(:)
           alpha4(:) = shell(lshell)%alpha(:)
           coeff1(:) = shell_auxil(ishell)%coeff(:)
           coeff2(:) = 1.0_dp
           coeff3(:) = shell(kshell)%coeff(:)
           coeff4(:) = shell(lshell)%coeff(:)
           x01(:) = shell_auxil(ishell)%x0(:)
           x02(:) = shell_auxil(ishell)%x0(:)
           x03(:) = shell(kshell)%x0(:)
           x04(:) = shell(lshell)%x0(:)

         else ! interexchange indexes

           am3 = shell_auxil(ishell)%am
           am4 = 0
           am1 = shell(kshell)%am
           am2 = shell(lshell)%am
           n3c = number_basis_function_am( CARTESIAN , ami )
           n4c = 1
           n1c = number_basis_function_am( CARTESIAN , amk )
           n2c = number_basis_function_am( CARTESIAN , aml )
           n3 = ni
           n4 = nj
           n1 = nk
           n2 = nl
           ng3 = shell_auxil(ishell)%ng
           ng4 = 1
           ng1 = shell(kshell)%ng
           ng2 = shell(lshell)%ng
           allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
           allocate(coeff1(ng1),coeff2(ng2),coeff3(ng3),coeff4(ng4))
           alpha3(:) = shell_auxil(ishell)%alpha(:) 
           alpha4(:) = 0.0_dp 
           alpha1(:) = shell(kshell)%alpha(:)
           alpha2(:) = shell(lshell)%alpha(:)
           coeff3(:) = shell_auxil(ishell)%coeff(:)
           coeff4(:) = 1.0_dp
           coeff1(:) = shell(kshell)%coeff(:)
           coeff2(:) = shell(lshell)%coeff(:)
           x03(:) = shell_auxil(ishell)%x0(:)
           x04(:) = shell_auxil(ishell)%x0(:)
           x01(:) = shell(kshell)%x0(:)
           x02(:) = shell(lshell)%x0(:)

         endif

         allocate( int_shell(n1c*n2c*n3c*n4c) )
         allocate( integrals_cart(n1c,n2c,n3c,n4c) )
         allocate( integrals_tmp (n1c,n2c,n3c,n4c) )
         integrals_cart(:,:,:,:) = 0.0_dp


         if(am1+am2+am3+am4==0) then

           do ig4=1,ng4
             do ig3=1,ng3
               do ig2=1,ng2
                 do ig1=1,ng1

                   zeta_12 = alpha1(ig1) + alpha2(ig2)
                   zeta_34 = alpha3(ig3) + alpha4(ig4)
                   p(:) = ( alpha1(ig1) * x01(:) + alpha2(ig2) * x02(:) ) / zeta_12 
                   q(:) = ( alpha3(ig3) * x03(:) + alpha4(ig4) * x04(:) ) / zeta_34 
                   !
                   ! Full range or long-range only integrals
                   rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
                   rho1 = rho
                   
                   tt = rho * SUM( (p(:)-q(:))**2 )
                   call boys_function(f0t(0),0,tt)

                   integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                         2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                         / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                         / zeta_34 * EXP( -alpha3(ig3)*alpha4(ig4)/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) &
                         * SQRT( rho / rho1 ) &
                         * coeff1(ig1) * coeff2(ig2) &
                         * coeff3(ig3) * coeff4(ig4)&
                         * cart_to_pure_norm(0)%matrix(1,1)**4

                 enddo
               enddo
             enddo
           enddo

         else

           info=eval_contr_integral(                &
                                   am1,am2,am3,am4, &
                                   ng1,ng2,ng3,ng4, &
                                   coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
                                   alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
                                   x01(1),x02(1),x03(1),x04(1),&
                                   rcut_libint, &
                                   int_shell(1))
    
    
           if(info/=0) then
             WRITE_MASTER(*,*) am1,am2,am3,am4
             stop 'ERI calculated by libint failed'
           endif
    
           iibf=0
           do ibf=1,n1c
             do jbf=1,n2c
               do kbf=1,n3c
                 do lbf=1,n4c
                   iibf=iibf+1
                   integrals_cart(ibf,jbf,kbf,lbf) = int_shell(iibf)
                 enddo
               enddo
             enddo
           enddo


           do lbf=1,n4c
             do kbf=1,n3c
               do jbf=1,n2c
                 do ibf=1,n1
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1c,jbf,kbf,lbf) * cart_to_pure_norm(am1)%matrix(1:n1c,ibf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4c
             do kbf=1,n3c
               do jbf=1,n2
                 do ibf=1,n1
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2c,kbf,lbf) * cart_to_pure_norm(am2)%matrix(1:n2c,jbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4c
             do kbf=1,n3
               do jbf=1,n2
                 do ibf=1,n1
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3c,lbf) * cart_to_pure_norm(am3)%matrix(1:n3c,kbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4
             do kbf=1,n3
               do jbf=1,n2
                 do ibf=1,n1
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4c) * cart_to_pure_norm(am4)%matrix(1:n4c,lbf) )
                 enddo
               enddo
             enddo
           enddo


         endif ! is (ss|ss)
         
         if(amk+aml>=ami+amj) then
           
           do lbf=1,nl
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   eri_3center( shell_auxil(ishell)%istart+ibf-1,    &
                          index_prod(shell(kshell)%istart+kbf-1,shell(lshell)%istart+lbf-1) ) = integrals_cart(ibf,jbf,kbf,lbf)
                 enddo
               enddo
             enddo
           enddo

         else

           do lbf=1,nl
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   eri_3center( shell_auxil(ishell)%istart+ibf-1,    &
                          index_prod(shell(kshell)%istart+kbf-1,shell(lshell)%istart+lbf-1) ) = integrals_cart(kbf,lbf,ibf,jbf)
                 enddo
               enddo
             enddo
           enddo

         endif


         deallocate(integrals_cart)
         deallocate(integrals_tmp)
         deallocate(int_shell)
         deallocate(alpha1,alpha2,alpha3,alpha4)
         deallocate(coeff1,coeff2,coeff3,coeff4)

       enddo
     enddo
!   enddo
! enddo
 enddo

 WRITE_MASTER(*,'(a)') ' All 3-center integrals have been calculated and stored'

 !
 ! Combine the 2-center integral into the 3-center and then get rid of them
 ! definitively
 eri_3center(:,:) = MATMUL( TRANSPOSE(eri_2center_m1) , eri_3center(:,:) )

 WRITE_MASTER(*,*) 'Now deallocate the 2-center integrals: not needed anymore'
 call memory_statement(-REAL(nsize_auxil,dp)*REAL(prec_eri/dp,dp))
 if(allocated(eri_2center_m1)) deallocate(eri_2center_m1)
 

 call stop_clock(timing_eri_3center)

end subroutine calculate_eri_3center


!=========================================================================
subroutine setup_negligible_basispair()
 implicit none
!====
 integer :: ishell,jshell
 integer :: ibf,jbf
!====

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

 nsize1 = 0
 do jbf=1,nbf_eri
   do ibf=1,jbf
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       nsize1 = nsize1 + 1
       index_pair(ibf,jbf) = nsize1
       index_pair(jbf,ibf) = nsize1
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
!       WRITE_MASTER(*,*) '    negl',max_ij,max_ij < TOL_INT
       negligible_basispair(ibf,jbf) = .TRUE.
     else
!       WRITE_MASTER(*,*) 'non negl',max_ij,max_ij < TOL_INT
       npair_refined = npair_refined + 1
     endif


   enddo
 enddo

 WRITE_MASTER(*,'(/,a)') ' Refining the negligible basis function pairs'
 WRITE_MASTER(*,'(a,x,i6)')   ' Non negligible basis function pairs stored in memory   ',npair
 WRITE_MASTER(*,'(a,x,i6)')   ' Non negligible basis function pairs used in calculation',npair_refined


end subroutine refine_negligible_basispair


!=========================================================================
subroutine identify_negligible_shellpair(basis)
!
! A first screening implementation
! Find negligible shell pair with
! Cauchy-Schwarz inequality
! (ij|1/r|kl)**2 <= (ij|1/r|ij) (kl|1/r|(kl) 
!
 use ISO_C_BINDING
 use m_tools,only: boys_function
 implicit none

 type(basis_set),intent(in)   :: basis
!====
 integer  :: info
 integer  :: iibf
 integer  :: ibf,jbf,kbf,lbf
! integer  :: n1,n2
 integer  :: n1c,n2c
 integer  :: ni,nj
 integer  :: ng1,ng2
 integer  :: ami,amj
 integer  :: ishell,jshell
 integer  :: ig1,ig2,ig3,ig4
 integer  :: nneglect
 real(dp) :: zeta_12,rho,rho1,f0t(0:0),tt
 real(dp) :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
!====
! variables used to call C
 integer(C_INT),external      :: libint_init,calculate_integral
 integer(C_INT),external      :: eval_contr_integral
 integer(C_INT)               :: am1,am2
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:)
 real(C_DOUBLE)               :: x01(3),x02(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:)
 real(C_DOUBLE)               :: rcut_libint
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====

 WRITE_MASTER(*,'(/,a)')    ' Cauchy-Schwartz screening of the 4-center integrals'
 WRITE_MASTER(*,'(a)')      ' Libint library initialized'
 WRITE_MASTER(*,'(a,i5,/)') ' Max angular momentum handled by your Libint compilation: ',libint_init()

 rcut_libint = 0.0_dp

 nneglect = 0

 do jshell=1,nshell
   do ishell=1,nshell
     ami = shell(ishell)%am
     amj = shell(jshell)%am
     !TODO: Here time could be saved by only checking ishell<= jshell
     ! But then an interexchange of indexes would have to be implemented in
     ! order to satisfy ami >= amj (required condition in libint)
     if( ami < amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     n1c = number_basis_function_am( CARTESIAN , ami )
     n2c = number_basis_function_am( CARTESIAN , amj )
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
         WRITE_MASTER(*,*) am1,am2,am1,am2
         stop 'ERI calculated by libint failed'
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

     negligible_shellpair(ishell,jshell)=.TRUE.
     do ibf=1,ni
       do jbf=1,nj
         if( ABS( integrals_cart(ibf,jbf,ibf,jbf) ) > TOL_INT ) negligible_shellpair(ishell,jshell)=.FALSE.
       enddo
     enddo

     !
     ! Symmetrize
     negligible_shellpair(jshell,ishell)=negligible_shellpair(ishell,jshell)

     if( negligible_shellpair(ishell,jshell) ) nneglect = nneglect + 1

     deallocate(integrals_cart)
     deallocate(integrals_tmp)
     deallocate(int_shell)
     deallocate(alpha1,alpha2)
     deallocate(coeff1,coeff2)

   enddo
 enddo

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
 WRITE_MASTER(*,'(/,a,i8,a,i8)') ' Non negligible shellpairs to be computed',nshellpair,'  over a total of',jshellpair
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
function libint_ordering(nx,ny,nz)
 implicit none
 integer,intent(in) :: nx,ny,nz
 integer            :: libint_ordering
!=====

 select case(nx+ny+nz)
 case(0)
                                 libint_ordering=1
 case(1)
   if(nx==1)                     libint_ordering=1
   if(ny==1)                     libint_ordering=2
   if(nz==1)                     libint_ordering=3
 case(2)
   if(nx==2          )           libint_ordering=1
   if(nx==1.AND.ny==1)           libint_ordering=2
   if(nx==1.AND.nz==1)           libint_ordering=3
   if(ny==2          )           libint_ordering=4
   if(ny==1.AND.nz==1)           libint_ordering=5
   if(nz==2          )           libint_ordering=6
 case(3)
   if(nx==3                    ) libint_ordering=1
   if(nx==2.AND.ny==1          ) libint_ordering=2
   if(nx==2.AND.nz==1          ) libint_ordering=3
   if(nx==1.AND.ny==2          ) libint_ordering=4
   if(nx==1.AND.ny==1.AND.nz==1) libint_ordering=5
   if(nx==1.AND.nz==2          ) libint_ordering=6
   if(ny==3                    ) libint_ordering=7
   if(ny==2.AND.nz==1          ) libint_ordering=8
   if(ny==1.AND.nz==2          ) libint_ordering=9
   if(nz==3                    ) libint_ordering=10
 case(4)
   if(nx==4                    ) libint_ordering=1
   if(nx==3.AND.ny==1          ) libint_ordering=2
   if(nx==3          .AND.nz==1) libint_ordering=3
   if(nx==2.AND.ny==2          ) libint_ordering=4
   if(nx==2.AND.ny==1.AND.nz==1) libint_ordering=5
   if(nx==2          .AND.nz==2) libint_ordering=6
   if(nx==1.AND.ny==3          ) libint_ordering=7
   if(nx==1.AND.ny==2          ) libint_ordering=8
   if(nx==1.AND.ny==1          ) libint_ordering=9
   if(nx==1.AND.ny==0.AND.nz==3) libint_ordering=10
   if(nx==0.AND.ny==4.AND.nz==0) libint_ordering=11
   if(nx==0.AND.ny==3.AND.nz==1) libint_ordering=12
   if(nx==0.AND.ny==2.AND.nz==2) libint_ordering=13
   if(          ny==1.AND.nz==3) libint_ordering=14
   if(                    nz==4) libint_ordering=15
 case(5)
   if(nx==5.AND.ny==0.AND.nz==0) libint_ordering=1
   if(nx==4.AND.ny==1.AND.nz==0) libint_ordering=2
   if(nx==4.AND.ny==0.AND.nz==1) libint_ordering=3
   if(nx==3.AND.ny==2.AND.nz==0) libint_ordering=4
   if(nx==3.AND.ny==1.AND.nz==1) libint_ordering=5
   if(nx==3.AND.ny==0.AND.nz==2) libint_ordering=6
   if(nx==2.AND.ny==3.AND.nz==0) libint_ordering=7
   if(nx==2.AND.ny==2.AND.nz==1) libint_ordering=8
   if(nx==2.AND.ny==1.AND.nz==2) libint_ordering=9
   if(nx==2.AND.ny==0.AND.nz==3) libint_ordering=10
   if(nx==1.AND.ny==4.AND.nz==0) libint_ordering=11
   if(nx==1.AND.ny==3.AND.nz==1) libint_ordering=12
   if(nx==1.AND.ny==2.AND.nz==2) libint_ordering=13
   if(nx==1.AND.ny==1.AND.nz==3) libint_ordering=14
   if(nx==1.AND.ny==0.AND.nz==4) libint_ordering=15
   if(nx==0.AND.ny==5.AND.nz==0) libint_ordering=16
   if(nx==0.AND.ny==4.AND.nz==1) libint_ordering=17
   if(nx==0.AND.ny==3.AND.nz==2) libint_ordering=18
   if(nx==0.AND.ny==2.AND.nz==3) libint_ordering=19
   if(nx==0.AND.ny==1.AND.nz==4) libint_ordering=20
   if(nx==0.AND.ny==0.AND.nz==5) libint_ordering=21
 case(6)
   if(nx==6.AND.ny==0.AND.nz==0) libint_ordering=1
   if(nx==5.AND.ny==1.AND.nz==0) libint_ordering=2
   if(nx==5.AND.ny==0.AND.nz==1) libint_ordering=3
   if(nx==4.AND.ny==2.AND.nz==0) libint_ordering=4
   if(nx==4.AND.ny==1.AND.nz==1) libint_ordering=5
   if(nx==4.AND.ny==0.AND.nz==2) libint_ordering=6
   if(nx==3.AND.ny==3.AND.nz==0) libint_ordering=7
   if(nx==3.AND.ny==2.AND.nz==1) libint_ordering=8
   if(nx==3.AND.ny==1.AND.nz==2) libint_ordering=9
   if(nx==3.AND.ny==0.AND.nz==3) libint_ordering=10
   if(nx==2.AND.ny==4.AND.nz==0) libint_ordering=11
   if(nx==2.AND.ny==3.AND.nz==1) libint_ordering=12
   if(nx==2.AND.ny==2.AND.nz==2) libint_ordering=13
   if(nx==2.AND.ny==1.AND.nz==3) libint_ordering=14
   if(nx==2.AND.ny==0.AND.nz==4) libint_ordering=15
   if(nx==1.AND.ny==5.AND.nz==0) libint_ordering=16
   if(nx==1.AND.ny==4.AND.nz==1) libint_ordering=17
   if(nx==1.AND.ny==3.AND.nz==2) libint_ordering=18
   if(nx==1.AND.ny==2.AND.nz==3) libint_ordering=19
   if(nx==1.AND.ny==1.AND.nz==4) libint_ordering=20
   if(nx==1.AND.ny==0.AND.nz==5) libint_ordering=21
   if(nx==0.AND.ny==6.AND.nz==0) libint_ordering=22
   if(nx==0.AND.ny==5.AND.nz==1) libint_ordering=23
   if(nx==0.AND.ny==4.AND.nz==2) libint_ordering=24
   if(nx==0.AND.ny==3.AND.nz==3) libint_ordering=25
   if(nx==0.AND.ny==2.AND.nz==4) libint_ordering=26
   if(nx==0.AND.ny==1.AND.nz==5) libint_ordering=27
   if(nx==0.AND.ny==0.AND.nz==6) libint_ordering=28
 case(7)
   if(nx==7.AND.ny==0.AND.nz==0) libint_ordering=1
   if(nx==6.AND.ny==1.AND.nz==0) libint_ordering=2
   if(nx==6.AND.ny==0.AND.nz==1) libint_ordering=3
   if(nx==5.AND.ny==2.AND.nz==0) libint_ordering=4
   if(nx==5.AND.ny==1.AND.nz==1) libint_ordering=5
   if(nx==5.AND.ny==0.AND.nz==2) libint_ordering=6
   if(nx==4.AND.ny==3.AND.nz==0) libint_ordering=7
   if(nx==4.AND.ny==2.AND.nz==1) libint_ordering=8
   if(nx==4.AND.ny==1.AND.nz==2) libint_ordering=9
   if(nx==4.AND.ny==0.AND.nz==3) libint_ordering=10
   if(nx==3.AND.ny==4.AND.nz==0) libint_ordering=11
   if(nx==3.AND.ny==3.AND.nz==1) libint_ordering=12
   if(nx==3.AND.ny==2.AND.nz==2) libint_ordering=13
   if(nx==3.AND.ny==1.AND.nz==3) libint_ordering=14
   if(nx==3.AND.ny==0.AND.nz==4) libint_ordering=15
   if(nx==2.AND.ny==5.AND.nz==0) libint_ordering=16
   if(nx==2.AND.ny==4.AND.nz==1) libint_ordering=17
   if(nx==2.AND.ny==3.AND.nz==2) libint_ordering=18
   if(nx==2.AND.ny==2.AND.nz==3) libint_ordering=19
   if(nx==2.AND.ny==1.AND.nz==4) libint_ordering=20
   if(nx==2.AND.ny==0.AND.nz==5) libint_ordering=21
   if(nx==1.AND.ny==6.AND.nz==0) libint_ordering=22
   if(nx==1.AND.ny==5.AND.nz==1) libint_ordering=23
   if(nx==1.AND.ny==4.AND.nz==2) libint_ordering=24
   if(nx==1.AND.ny==3.AND.nz==3) libint_ordering=25
   if(nx==1.AND.ny==2.AND.nz==4) libint_ordering=26
   if(nx==1.AND.ny==1.AND.nz==5) libint_ordering=27
   if(nx==1.AND.ny==0.AND.nz==6) libint_ordering=28
   if(nx==0.AND.ny==7.AND.nz==0) libint_ordering=29
   if(nx==0.AND.ny==6.AND.nz==1) libint_ordering=30
   if(nx==0.AND.ny==5.AND.nz==2) libint_ordering=31
   if(nx==0.AND.ny==4.AND.nz==3) libint_ordering=32
   if(nx==0.AND.ny==3.AND.nz==4) libint_ordering=33
   if(nx==0.AND.ny==2.AND.nz==5) libint_ordering=34
   if(nx==0.AND.ny==1.AND.nz==6) libint_ordering=35
   if(nx==0.AND.ny==0.AND.nz==7) libint_ordering=36

 case default
   stop'libint_ordering not coded for this orbital momentum'
 end select

end function libint_ordering



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
           WRITE_MASTER(*,*) ibf,jbf,kbf,lbf,eri(ibf,jbf,kbf,lbf)
           WRITE_MASTER(*,*) kbf,lbf,ibf,jbf,eri(kbf,lbf,ibf,jbf)
           WRITE_MASTER(*,*) ibf,basis%bf(ibf)%amc
           WRITE_MASTER(*,*) jbf,basis%bf(jbf)%amc
           WRITE_MASTER(*,*) kbf,basis%bf(kbf)%amc
           WRITE_MASTER(*,*) lbf,basis%bf(lbf)%amc
           stop'ERI array not symmetric'
         endif
       enddo
     enddo
   enddo
 enddo

 stop'TESTING OK'

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

!$OMP PARALLEL DEFAULT(SHARED)

!$OMP DO SCHEDULE(STATIC)
 do lbf=1,nbf_eri
   do kbf=1,nbf_eri
     do jbf=1,nbf_eri

       do ibf=1,nbf_eri
         eri_tmp3(jbf,kbf,lbf) = eri_tmp3(jbf,kbf,lbf) + eri(ibf,jbf,kbf,lbf) * c_matrix(ibf,istate,ijspin) 
       enddo


     enddo
   enddo
 enddo
!$OMP END DO


!$OMP DO SCHEDULE(STATIC)
 do lbf=1,nbf_eri
   do kbf=1,nbf_eri

     do jstate=1,nbf_eri
       eri_eigenstate_i(jstate,kbf,lbf,nspin) = DOT_PRODUCT( eri_tmp3(:,kbf,lbf) , c_matrix(:,jstate,ijspin) )
     enddo

   enddo
 enddo
!$OMP END DO

!$OMP END PARALLEL

  
 do klspin=1,nspin

!$OMP PARALLEL DEFAULT(SHARED)

!$OMP DO SCHEDULE(STATIC)
   do lbf=1,nbf_eri
     do kstate=1,nbf_eri
       do jstate=1,nbf_eri
         eri_tmp3(jstate,kstate,lbf) = DOT_PRODUCT( eri_eigenstate_i(jstate,:,lbf,nspin) , c_matrix(:,kstate,klspin) )
       enddo
     enddo
   enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
   do lstate=1,nbf_eri
     do kstate=1,nbf_eri
       do jstate=1,nbf_eri

         eri_eigenstate_i(jstate,kstate,lstate,klspin) = DOT_PRODUCT( eri_tmp3(jstate,kstate,:) , c_matrix(:,lstate,klspin) )

       enddo
     enddo
   enddo
!$OMP END DO

!$OMP END PARALLEL

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
!=====

 call start_clock(timing_eri_3center_eigen)

 WRITE_MASTER(*,'(/,a)') ' Calculate 3-center integrals on eigenstates'


 !TODO merge the 2 last indexes for prod_basis save a factor 2! (i<->j symmetry)
 allocate(eri_3center_eigen(nsize1_auxil,nbf_eri,nbf_eri,nspin))
 WRITE_MASTER(*,*) 'Allocate 3-center integrals on states'
 call memory_statement(REAL(nsize1_auxil,dp)*REAL(nbf_eri,dp)**2*REAL(prec_eri,dp)/REAL(dp,dp))

 allocate(eri_3center_tmp(nsize1_auxil,nbf_eri,nbf_eri)) 
 eri_3center_eigen(:,:,:,:) = 0.0_dp
 do klspin=1,nspin
   ! Transformation of the first index
   eri_3center_tmp(:,:,:) = 0.0_dp
   do kbf=1,nbf_eri
     do lbf=1,nbf_eri
       if( negligible_basispair(kbf,lbf) ) cycle

         do lstate=1,nbf_eri
           eri_3center_tmp(:,kbf,lstate) = eri_3center_tmp(:,kbf,lstate) &
                                      + c_matrix(lbf,lstate,klspin) * eri_3center(:,index_prod(kbf,lbf))
         enddo

     enddo
   enddo
   ! Transformation of the second index
   do lstate=1,nbf_eri
     eri_3center_eigen(:,:,lstate,klspin) = MATMUL( eri_3center_tmp(:,:,lstate) , c_matrix(:,:,klspin) )
   enddo

 enddo ! klspin
 deallocate(eri_3center_tmp)

 WRITE_MASTER(*,'(a,/)') ' Done'

 call stop_clock(timing_eri_3center_eigen)

end subroutine prepare_eri_3center_eigen


!=================================================================
subroutine destroy_eri_3center_eigen()
 implicit none
!=====

 WRITE_MASTER(*,'(/,a,/)') ' Destroy 3-center integrals on eigenstates'

 if(allocated(eri_3center_eigen)) deallocate(eri_3center_eigen)

end subroutine destroy_eri_3center_eigen


!=================================================================
subroutine destroy_eri_3center()
 implicit none
!=====

 WRITE_MASTER(*,'(/,a,/)') ' Destroy 3-center integrals'
 call memory_statement(-REAL(nsize1_auxil,dp)*REAL(nsize1)*REAL(prec_eri/dp,dp))
 if(allocated(eri_3center)) deallocate(eri_3center)

end subroutine destroy_eri_3center


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
   if( ABS( eri_buffer(ibuffer) ) < tol ) icount=icount+1
 enddo

 WRITE_MASTER(*,*) ' number of negligible integrals <',tol
 WRITE_MASTER(*,*) icount, ' / ',nsize,REAL(icount,dp)/REAL(nsize,dp)*100.0_dp,' [%]'


 do ibf=1,nbf_eri
   do jbf=1,nbf_eri
     integral_ij(ibf,jbf) = eri(ibf,jbf,ibf,jbf)
   enddo
 enddo

 WRITE_MASTER(*,*) 'testing Cauchy-Schwarz condition'
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
 WRITE_MASTER(*,*) ' number of negligible integrals <',tol
 WRITE_MASTER(*,*) icount, ' / ',nbf_eri**4,REAL(icount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'
 WRITE_MASTER(*,*) jcount, ' / ',nbf_eri**4,REAL(jcount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'


end subroutine negligible_eri


!=========================================================================
subroutine dump_out_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!====
 character(len=50) :: filename
 integer           :: nline,iline,icurrent
!====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 WRITE_MASTER(*,*) 'Dump out the ERI into file'
 WRITE_MASTER(*,*) 'Size of file [bytes]',REAL(nsize,dp)*prec_eri

 open(unit=111,file=TRIM(filename),form='unformatted')
 WRITE_MASTER(111) nsize
 WRITE_MASTER(111) rcut

 nline = nsize / line_length + 1
 icurrent=0
 do iline=1,nline
   WRITE_MASTER(111) eri_buffer(icurrent+1:MIN(nsize,icurrent+line_length+1))
   icurrent = icurrent + line_length + 1
 enddo

 close(111)

 WRITE_MASTER(*,'(a,/)') ' file written'

end subroutine dump_out_eri


!=========================================================================
logical function read_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!====
 character(len=50) :: filename
 integer           :: nline,iline,icurrent
 integer           :: integer_read
 real(dp)          :: real_read
!====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 
 inquire(file=TRIM(filename),exist=read_eri)

 if(read_eri) then

   WRITE_MASTER(*,*) 'Try to read ERI file'
   open(unit=111,file=TRIM(filename),form='unformatted',status='old')
   read(111) integer_read
   if(integer_read /= nsize) read_eri=.FALSE.
   read(111) real_read
   if(ABS(real_read-rcut) > 1.0d-6) read_eri=.FALSE.

   if(read_eri) then

     nline = nsize / line_length + 1
     icurrent=0
     do iline=1,nline
       read(111) eri_buffer(icurrent+1:MIN(nsize,icurrent+line_length+1))
       icurrent = icurrent + line_length + 1
     enddo
     WRITE_MASTER(*,'(a,/)') ' ERI file read'

   else
     WRITE_MASTER(*,'(a,/)') ' reading aborted'
   endif

   close(111)

 endif


end function read_eri


!=========================================================================
end module m_eri
