!=========================================================================
#include "macros.h"
!=========================================================================
module m_eri
 use m_definitions
 use m_mpi
 use m_basis_set

 private
 public :: eri,allocate_eri,allocate_eri_eigen,deallocate_eri,calculate_eri,transform_eri_basis_fast,transform_eri_basis_lowmem, &
           eri_lr,negligible_eri,&
           allocate_eri_auxil,calculate_eri_2center,calculate_eri_3center,&
           BUFFER1,BUFFER2
 public :: index_prod
 public :: negligible_basispair,refine_negligible_basispair

 integer,parameter :: BUFFER1 = 1
 integer,parameter :: BUFFER2 = 2
 !
 ! max length of a record in the ERI file
 integer,parameter :: line_length=1000

 real(dp)                   :: TOL_INT=1.0e-10_dp

 real(prec_eri),allocatable :: eri_buffer(:)
 real(prec_eri),allocatable :: eri_buffer_lr(:)
 real(prec_eri),allocatable :: eri_2center_m1(:,:)
 real(prec_eri),allocatable :: eri_3center(:,:)

 logical,allocatable        :: negligible_shellpair(:,:)
 logical,allocatable        :: negligible_basispair(:,:)
 integer,allocatable        :: index_pair(:,:)

 type shell_type
   integer              :: am
   integer              :: ng
   real(dp),allocatable :: alpha(:)
   real(dp),allocatable :: coeff(:)
   real(dp)             :: x0(3)
   integer              :: istart,iend
 end type shell_type
 integer                      :: nshell
 integer                      :: nshell_auxil
 type(shell_type),allocatable :: shell(:)
 type(shell_type),allocatable :: shell_auxil(:)


 integer                    :: nbf_eri                ! local copy of nbf
 integer                    :: nsize                  ! size of the eri_buffer array
 integer                    :: nsize1                 ! number of independent pairs (i,j) with i<=j

 integer                    :: nbf_eri_auxil          ! local copy of nbf for auxiliary basis
 integer                    :: nsize_auxil            ! size of the eri_buffer array
 integer                    :: nsize1_auxil           ! number of independent pairs (i,j) with i<=j

contains

!=========================================================================
subroutine allocate_eri(basis,rcut,which_buffer)
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
   call identify_negligible_shellpair(basis,rcut)
   call setup_negligible_basispair()
 endif


#if defined LOW_MEMORY2
 WRITE_MASTER(*,'(/,a)') ' Symmetrized ERI stored (8 symmetries)' 
 nsize1  = index_prod(nbf_eri,nbf_eri) 
 nsize   = index_eri(nbf_eri,nbf_eri,nbf_eri,nbf_eri)
#elif defined LOW_MEMORY3
 nsize = (nsize1*(nsize1+1))/2
#else
 WRITE_MASTER(*,'(/,a)') ' All ERI are stored'
 nsize1  = nbf_eri**2
 nsize   = nsize1**2
 stop'Not supported anymore'
#endif

 WRITE_MASTER(*,*) 'number of integrals to be stored:',nsize
 WRITE_MASTER(*,*) 'max index size',HUGE(nsize)
 if(nsize<1) stop'too many integrals to be stored'

 if(REAL(nsize,dp)*prec_eri > 1024**3 ) then
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the ERI array: ',REAL(nsize,dp)*prec_eri/1024**3,' [Gb] / proc'
 else
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the ERI array: ',REAL(nsize,dp)*prec_eri/1024**2,' [Mb] / proc'
 endif

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



end subroutine allocate_eri


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

 WRITE_MASTER(*,*) 'number of integrals to be stored:',nsize_auxil
 WRITE_MASTER(*,*) 'max index size',HUGE(nsize_auxil)
 if(nsize_auxil<1) stop'too many integrals to be stored'

 if(REAL(nsize_auxil,dp)*prec_eri > 1024**3 ) then
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the ERI array: ',REAL(nsize_auxil,dp)*prec_eri/1024**3,' [Gb] / proc'
 else
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the ERI array: ',REAL(nsize_auxil,dp)*prec_eri/1024**2,' [Mb] / proc'
 endif

 allocate(eri_2center_m1(nsize1_auxil,nsize1_auxil),stat=info)
 eri_2center_m1(:,:) = 0.0_dp

 if(info==0) then
   WRITE_MASTER(*,*) 'success'
 else
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif

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

#ifdef LOW_MEMORY3
 index_prod = index_pair(ibf,jbf)
#else
 imax=MAX(ibf,jbf)
 jmin=MIN(ibf,jbf)
 index_prod = (jmin-1)*nbf_eri - (jmin-1)*(jmin-2)/2 + imax-jmin+1
#endif

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

#if defined LOW_MEMORY2 || defined LOW_MEMORY3
 index_ij = index_prod(ibf,jbf)
 index_kl = index_prod(kbf,lbf)

 ijmax=MAX(index_ij,index_kl)
 klmin=MIN(index_ij,index_kl)

 index_eri = (klmin-1)*nsize1 - (klmin-1)*(klmin-2)/2 + ijmax-klmin+1
#else
 index_eri = ibf+(jbf-1)*nbf_eri+(kbf-1)*nbf_eri**2+(lbf-1)*nbf_eri**3
#endif


end function index_eri


!=========================================================================
function eri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri
!=====
 integer            :: index_ijkl
 integer            :: i1,i2,i3,i4
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
 integer            :: index_ijkl
 integer            :: i1,i2,i3,i4
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_lr = 0.0_dp
 else
   eri_lr = eri_buffer_lr(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri_lr


!=========================================================================
subroutine calculate_eri(print_eri,basis,rcut,which_buffer)
 implicit none
 logical,intent(in)           :: print_eri
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: which_buffer
!=====


 if( .NOT. read_eri(rcut) ) call do_calculate_eri(basis,rcut,which_buffer)


 if( print_eri ) then
   call dump_out_eri(rcut)
 endif

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
 use m_timing
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: which_buffer
!=====
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: n1,n2,n3,n4
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
 real(dp),allocatable         :: integrals_libint(:)
!=====
! variables used to call C++ 
 integer(C_INT),external      :: calculate_integral
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: omega_range
!=====

#ifndef LOW_MEMORY2 
#ifndef LOW_MEMORY3
 stop'This implementation is only compatible with preprocessing option LOW_MEMORY2 or LOW_MEMORY3'
#endif
#endif

 WRITE_MASTER(*,'(/,a)') ' Calculate and store all the Electron Repulsion Integrals (ERI)'

 if( rcut > 1.0e-6_dp ) then
   omega_range = 1.0_dp / rcut
   WRITE_MASTER(*,'(a40,x,f9.4)') ' Long-Range only integrals with rcut=',rcut
   WRITE_MASTER(*,'(a40,x,f9.4)') ' or omega=',omega_range
 else 
   omega_range = 1.0e6_dp
 endif


!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(ami,amj,amk,aml,ni,nj,nk,nl,am1,am2,am3,am4,n1,n2,n3,n4,ng1,ng2,ng3,ng4,&
!$OMP&   alpha1,alpha2,alpha3,alpha4,x01,x02,x03,x04,&
!$OMP&   integrals_libint,integrals_cart,integrals_tmp,&
!$OMP&   zeta_12,zeta_34,p,q,rho,rho1,tt,f0t,&
!$OMP&   info,iibf)

!$OMP DO SCHEDULE(DYNAMIC) 
 do lshell=1,nshell
   do kshell=1,nshell
     !
     ! Order the angular momenta so that libint is pleased
     ! 1) am3+am4 >= am1+am2
     ! 2) am3>=am4
     ! 3) am1>=am2
     amk = shell(kshell)%am
     aml = shell(lshell)%am
     if( amk < aml ) cycle
     if( negligible_shellpair(kshell,lshell) ) cycle

     do jshell=1,nshell
       do ishell=1,nshell
         ami = shell(ishell)%am
         amj = shell(jshell)%am
         if( ami < amj ) cycle
         if( amk+aml < ami+amj ) cycle
         if( negligible_shellpair(ishell,jshell) ) cycle


         ni = number_basis_function_am( basis%gaussian_type , ami )
         nj = number_basis_function_am( basis%gaussian_type , amj )
         nk = number_basis_function_am( basis%gaussian_type , amk )
         nl = number_basis_function_am( basis%gaussian_type , aml )


         am1 = shell(ishell)%am
         am2 = shell(jshell)%am
         am3 = shell(kshell)%am
         am4 = shell(lshell)%am
         n1 = number_basis_function_am( CARTESIAN , ami )
         n2 = number_basis_function_am( CARTESIAN , amj )
         n3 = number_basis_function_am( CARTESIAN , amk )
         n4 = number_basis_function_am( CARTESIAN , aml )
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

         allocate( integrals_libint( n1*n2*n3*n4 ) )
         allocate( integrals_cart(n1,n2,n3,n4) )
         allocate( integrals_tmp(n1,n2,n3,n4) )
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
                   if( rcut < 2.0e-6_dp ) then
                     rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
                     rho1 = rho
                   else
                     rho  = zeta_12 * zeta_34 * omega_range**2 / ( zeta_12*omega_range**2 + zeta_34*omega_range**2 + zeta_12*zeta_34 )
                     rho1 = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
                   endif
                   tt = rho * SUM( (p(:)-q(:))**2 )
                   call boys_function(f0t(0),0,tt)

                   integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                         2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                         / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                         / zeta_34 * EXP( -alpha3(ig3)*alpha4(ig4)/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) &
                         * SQRT( rho / rho1 ) &
                         * shell(ishell)%coeff(ig1) &
                         * shell(jshell)%coeff(ig2) &
                         * shell(kshell)%coeff(ig3) &
                         * shell(lshell)%coeff(ig4) * cart_to_pure_norm(0)%matrix(1,1)**4

                 enddo
               enddo
             enddo
           enddo

         else

           do ig4=1,ng4
             do ig3=1,ng3
               do ig2=1,ng2
                 do ig1=1,ng1

! call start_clock(timing_tmp2)
                   info=calculate_integral(omega_range,&
                                           am1,am2,am3,am4,alpha1(ig1),alpha2(ig2),alpha3(ig3),alpha4(ig4),&
                                           x01(1),x01(2),x01(3),&
                                           x02(1),x02(2),x02(3),&
                                           x03(1),x03(2),x03(3),&
                                           x04(1),x04(2),x04(3),&
                                           integrals_libint(1))

                   if(info/=0) then
                     WRITE_MASTER(*,*) am1,am2,am3,am4
                     stop 'ERI calculated by libint failed'
                   endif
! call stop_clock(timing_tmp2)
! call start_clock(timing_tmp1)

                   iibf=0
                   do ibf=1,n1
                     do jbf=1,n2
                       do kbf=1,n3
                         do lbf=1,n4
                           iibf=iibf+1
                           integrals_cart(ibf,jbf,kbf,lbf) = integrals_cart(ibf,jbf,kbf,lbf) &
                                                            + integrals_libint(iibf) * shell(ishell)%coeff(ig1) * shell(jshell)%coeff(ig2) &
                                                                                     * shell(kshell)%coeff(ig3) * shell(lshell)%coeff(ig4)
                         enddo
                       enddo
                     enddo
                   enddo
! call stop_clock(timing_tmp1)

                 enddo
               enddo
             enddo
           enddo

! call start_clock(timing_tmp3)

           do lbf=1,n4
             do kbf=1,n3
               do jbf=1,n2
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1,jbf,kbf,lbf) * cart_to_pure_norm(shell(ishell)%am)%matrix(1:n1,ibf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4
             do kbf=1,n3
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2,kbf,lbf) * cart_to_pure_norm(shell(jshell)%am)%matrix(1:n2,jbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3,lbf) * cart_to_pure_norm(shell(kshell)%am)%matrix(1:n3,kbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,nl
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4) * cart_to_pure_norm(shell(lshell)%am)%matrix(1:n4,lbf) )
                 enddo
               enddo
             enddo
           enddo

! call stop_clock(timing_tmp3)



         endif
         
! call start_clock(timing_tmp4)

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
! call stop_clock(timing_tmp4)



         deallocate(integrals_cart)
         deallocate(integrals_tmp)
         deallocate(integrals_libint)
         deallocate(alpha1,alpha2,alpha3,alpha4)

       enddo
     enddo
   enddo
 enddo
!$OMP END DO
!$OMP END PARALLEL


 WRITE_MASTER(*,'(a,/)') ' All ERI have been calculated'


end subroutine do_calculate_eri


!=========================================================================
subroutine calculate_eri_2center(print_eri,auxil_basis)
 use ISO_C_BINDING
 use m_tools,only: boys_function, invert
 use m_timing
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none
 logical,intent(in)           :: print_eri
 type(basis_set),intent(in)   :: auxil_basis
!=====
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: n1,n2,n3,n4
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
 real(dp),allocatable         :: integrals_libint(:)
!=====
! variables used to call C++ 
 integer(C_INT),external      :: calculate_integral
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: omega_range
!=====


 WRITE_MASTER(*,'(/,a)') ' Calculate, invert and store the 2-center Electron Repulsion Integrals'

 omega_range = 1.0e6_dp

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
         n1 = number_basis_function_am( CARTESIAN , ami )
         n2 = 1
         n3 = number_basis_function_am( CARTESIAN , amk )
         n4 = 1
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

         allocate( integrals_libint( n1*n2*n3*n4 ) )
         allocate( integrals_cart(n1,n2,n3,n4) )
         allocate( integrals_tmp(n1,n2,n3,n4) )
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
                         * shell_auxil(ishell)%coeff(ig1) &
                         * shell_auxil(kshell)%coeff(ig3) &
                         * cart_to_pure_norm(0)%matrix(1,1)**4

                 enddo
               enddo
             enddo
           enddo

         else

           do ig4=1,ng4
             do ig3=1,ng3
               do ig2=1,ng2
                 do ig1=1,ng1

! call start_clock(timing_tmp2)
                   info=calculate_integral(omega_range,&
                                           am1,am2,am3,am4,alpha1(ig1),alpha2(ig2),alpha3(ig3),alpha4(ig4),&
                                           x01(1),x01(2),x01(3),&
                                           x02(1),x02(2),x02(3),&
                                           x03(1),x03(2),x03(3),&
                                           x04(1),x04(2),x04(3),&
                                           integrals_libint(1))

                   if(info/=0) then
                     WRITE_MASTER(*,*) am1,am2,am3,am4
                     stop 'ERI calculated by libint failed'
                   endif
! call stop_clock(timing_tmp2)
! call start_clock(timing_tmp1)

                   iibf=0
                   do ibf=1,n1
                     do jbf=1,n2
                       do kbf=1,n3
                         do lbf=1,n4
                           iibf=iibf+1
                           integrals_cart(ibf,jbf,kbf,lbf) = integrals_cart(ibf,jbf,kbf,lbf) &
                                                            + integrals_libint(iibf) * shell_auxil(ishell)%coeff(ig1)  &
                                                                                     * shell_auxil(kshell)%coeff(ig3) 
                         enddo
                       enddo
                     enddo
                   enddo
! call stop_clock(timing_tmp1)

                 enddo
               enddo
             enddo
           enddo

! call start_clock(timing_tmp3)

           do lbf=1,n4
             do kbf=1,n3
               do jbf=1,n2
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1,jbf,kbf,lbf) * cart_to_pure_norm(shell_auxil(ishell)%am)%matrix(1:n1,ibf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4
             do kbf=1,n3
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2,kbf,lbf) * cart_to_pure_norm(shell_auxil(jshell)%am)%matrix(1:n2,jbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3,lbf) * cart_to_pure_norm(shell_auxil(kshell)%am)%matrix(1:n3,kbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,nl
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4) * cart_to_pure_norm(shell_auxil(lshell)%am)%matrix(1:n4,lbf) )
                 enddo
               enddo
             enddo
           enddo

! call stop_clock(timing_tmp3)



         endif
         
! call start_clock(timing_tmp4)

         do lbf=1,nl
           do kbf=1,nk
             do jbf=1,nj
               do ibf=1,ni
                 eri_2center_m1( shell_auxil(ishell)%istart+ibf-1,    &
                                 shell_auxil(kshell)%istart+kbf-1 )    = integrals_cart(ibf,jbf,kbf,lbf)
               enddo
             enddo
           enddo
         enddo
! call stop_clock(timing_tmp4)

         deallocate(integrals_cart)
         deallocate(integrals_tmp)
         deallocate(integrals_libint)
         deallocate(alpha1,alpha2,alpha3,alpha4)

       enddo
     enddo
   enddo
 enddo

 !
 ! Perform in-place inversion here
 call invert(nsize1,eri_2center_m1)

 WRITE_MASTER(*,'(a,/)') ' All 2-center integrals have been calculated, inverted and stored'


end subroutine calculate_eri_2center


!=========================================================================
subroutine calculate_eri_3center(print_eri,basis,auxil_basis)
 use ISO_C_BINDING
 use m_tools,only: boys_function
 use m_timing
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none
 logical,intent(in)           :: print_eri
 type(basis_set),intent(in)   :: basis
 type(basis_set),intent(in)   :: auxil_basis
!=====
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: n1,n2,n3,n4
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
 real(dp),allocatable         :: integrals_libint(:)
!=====
! variables used to call C++ 
 integer(C_INT),external      :: calculate_integral
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: omega_range
!=====


 WRITE_MASTER(*,'(/,a)') ' Calculate and store all the 3-center Electron Repulsion Integrals'

 omega_range = 1.0e6_dp


 do lshell=1,nshell
   do kshell=1,nshell
     !
     ! Order the angular momenta so that libint is pleased
     ! 1) am3+am4 >= am1+am2
     ! 2) am3>=am4
     ! 3) am1>=am2
     amk = shell(kshell)%am
     aml = shell(lshell)%am
     if( amk < aml ) cycle

     do jshell=1,1  ! FAKE LOOP
       do ishell=1,nshell_auxil
         ami = shell_auxil(ishell)%am
         amj = 0
         if( ami < amj ) cycle
         if( amk+aml < ami+amj ) cycle

         ni = number_basis_function_am( auxil_basis%gaussian_type , ami )
         nj = 1
         nk = number_basis_function_am( basis%gaussian_type , amk )
         nl = number_basis_function_am( basis%gaussian_type , aml )


         am1 = shell_auxil(ishell)%am
         am2 = 0
         am3 = shell(kshell)%am
         am4 = shell(kshell)%am
         n1 = number_basis_function_am( CARTESIAN , ami )
         n2 = 1
         n3 = number_basis_function_am( CARTESIAN , amk )
         n4 = number_basis_function_am( CARTESIAN , aml )
         ng1 = shell_auxil(ishell)%ng
         ng2 = 1
         ng3 = shell(kshell)%ng
         ng4 = shell(lshell)%ng
         allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
         alpha1(:) = shell_auxil(ishell)%alpha(:) 
         alpha2(:) = 0.0_dp ! shell_auxil(jshell)%alpha(:)
         alpha3(:) = shell(kshell)%alpha(:)
         alpha4(:) = shell(lshell)%alpha(:)
         x01(:) = shell_auxil(ishell)%x0(:)
         x02(:) = shell_auxil(ishell)%x0(:)
         x03(:) = shell(kshell)%x0(:)
         x04(:) = shell(lshell)%x0(:)

         allocate( integrals_libint( n1*n2*n3*n4 ) )
         allocate( integrals_cart(n1,n2,n3,n4) )
         allocate( integrals_tmp(n1,n2,n3,n4) )
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
                         * shell_auxil(ishell)%coeff(ig1) &
                         * shell_auxil(kshell)%coeff(ig3) &
                         * cart_to_pure_norm(0)%matrix(1,1)**4

                 enddo
               enddo
             enddo
           enddo

         else

           do ig4=1,ng4
             do ig3=1,ng3
               do ig2=1,ng2
                 do ig1=1,ng1

                   info=calculate_integral(omega_range,&
                                           am1,am2,am3,am4,alpha1(ig1),alpha2(ig2),alpha3(ig3),alpha4(ig4),&
                                           x01(1),x01(2),x01(3),&
                                           x02(1),x02(2),x02(3),&
                                           x03(1),x03(2),x03(3),&
                                           x04(1),x04(2),x04(3),&
                                           integrals_libint(1))

                   if(info/=0) then
                     WRITE_MASTER(*,*) am1,am2,am3,am4
                     stop 'ERI calculated by libint failed'
                   endif

                   iibf=0
                   do ibf=1,n1
                     do jbf=1,n2
                       do kbf=1,n3
                         do lbf=1,n4
                           iibf=iibf+1
                           integrals_cart(ibf,jbf,kbf,lbf) = integrals_cart(ibf,jbf,kbf,lbf) &
                                                            + integrals_libint(iibf) * shell_auxil(ishell)%coeff(ig1)  &
                                                                         * shell(kshell)%coeff(ig3) * shell(lshell)%coeff(ig4)
                         enddo
                       enddo
                     enddo
                   enddo

                 enddo
               enddo
             enddo
           enddo


           do lbf=1,n4
             do kbf=1,n3
               do jbf=1,n2
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1,jbf,kbf,lbf) * cart_to_pure_norm(shell_auxil(ishell)%am)%matrix(1:n1,ibf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4
             do kbf=1,n3
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2,kbf,lbf) * cart_to_pure_norm(shell_auxil(jshell)%am)%matrix(1:n2,jbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,n4
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3,lbf) * cart_to_pure_norm(shell(kshell)%am)%matrix(1:n3,kbf) )
                 enddo
               enddo
             enddo
           enddo

           do lbf=1,nl
             do kbf=1,nk
               do jbf=1,nj
                 do ibf=1,ni
                   integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4) * cart_to_pure_norm(shell(lshell)%am)%matrix(1:n4,lbf) )
                 enddo
               enddo
             enddo
           enddo


         endif
         

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


         deallocate(integrals_cart)
         deallocate(integrals_tmp)
         deallocate(integrals_libint)
         deallocate(alpha1,alpha2,alpha3,alpha4)

       enddo
     enddo
   enddo
 enddo

 WRITE_MASTER(*,'(a,/)') ' All 3-center integrals have been calculated and stored'


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

 WRITE_MASTER(*,*) 'Refining the negligible pairs'
 WRITE_MASTER(*,*) 'Non negligible pairs stored in memory   ',npair
 WRITE_MASTER(*,*) 'Non negligible pairs used in calculation',npair_refined


end subroutine refine_negligible_basispair


!=========================================================================
subroutine identify_negligible_shellpair(basis,rcut)
!
! A first screening implementation
! Find negligible shell pair with
! Cauchy-Schwarz inequality
! (ij|1/r|kl)**2 <= (ij|1/r|ij) (kl|1/r|(kl) 
!
 use ISO_C_BINDING
 use m_tools,only: boys_function
 use m_timing
 implicit none

 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
!====
 integer :: info
 integer :: iibf
 integer :: ibf,jbf,kbf,lbf
 integer :: n1,n2
 integer :: ni,nj
 integer :: ng1,ng2
 integer :: ami,amj
 integer :: ishell,jshell
 integer :: ig1,ig2,ig3,ig4
 integer :: neval,nneglect
 real(dp) :: zeta_12,rho,rho1,f0t(0:0),tt
 real(dp) :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
 real(dp),allocatable         :: integrals_libint(:)
!====
! variables used to call C++ 
 integer(C_INT),external      :: calculate_integral
 integer(C_INT)               :: am1,am2
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:)
 real(C_DOUBLE)               :: x01(3),x02(3)
 real(C_DOUBLE)               :: omega_range
!=====

 if( rcut > 1.0e-6_dp ) then
   omega_range = 1.0_dp / rcut
 else 
   omega_range = 1.0e6_dp
 endif

 neval    = 0
 nneglect = 0


 do jshell=1,nshell
   do ishell=1,nshell
     ami = shell(ishell)%am
     amj = shell(jshell)%am
     if( ami < amj ) cycle
     neval = neval + 1

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     n1 = number_basis_function_am( CARTESIAN , ami )
     n2 = number_basis_function_am( CARTESIAN , amj )
     am1 = shell(ishell)%am
     am2 = shell(jshell)%am
     ng1 = shell(ishell)%ng
     ng2 = shell(jshell)%ng

     allocate(alpha1(ng1),alpha2(ng2))
     alpha1(:) = shell(ishell)%alpha(:)
     alpha2(:) = shell(jshell)%alpha(:)
     x01(:) = shell(ishell)%x0(:)
     x02(:) = shell(jshell)%x0(:)

     allocate( integrals_libint( n1*n2*n1*n2 ) )
     allocate( integrals_cart(n1,n2,n1,n2) )
     allocate( integrals_tmp (n1,n2,n1,n2) )

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
               ! Full range or long-range only integrals
               if( rcut < 1.0e-6_dp ) then
                 rho  = zeta_12 * zeta_12 / ( zeta_12 + zeta_12 )
                 rho1 = rho
               else
                 rho  = zeta_12 * zeta_12 * omega_range**2 / ( zeta_12*omega_range**2 + zeta_12*omega_range**2 + zeta_12*zeta_12 )
                 rho1 = zeta_12 * zeta_12 / ( zeta_12 + zeta_12 )
               endif
               tt = rho * SUM( (p(:)-q(:))**2 )
               call boys_function(f0t(0),0,tt)

               integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                     2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_12 ) * f0t(0) &
                     / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                     / zeta_12 * EXP( -alpha1(ig3)*alpha2(ig4)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) &
                     * SQRT( rho / rho1 ) &
                     * shell(ishell)%coeff(ig1) &
                     * shell(jshell)%coeff(ig2) &
                     * shell(ishell)%coeff(ig3) &
                     * shell(jshell)%coeff(ig4) * cart_to_pure_norm(0)%matrix(1,1)**4

             enddo
           enddo
         enddo
       enddo

     else
       do ig4=1,ng2
         do ig3=1,ng1
           do ig2=1,ng2
             do ig1=1,ng1
               info=calculate_integral(omega_range,&
                                       am1,am2,am1,am2,alpha1(ig1),alpha2(ig2),alpha1(ig3),alpha2(ig4),&
                                       x01(1),x01(2),x01(3),&
                                       x02(1),x02(2),x02(3),&
                                       x01(1),x01(2),x01(3),&
                                       x02(1),x02(2),x02(3),&
                                       integrals_libint(1))
               if(info/=0) then
                 WRITE_MASTER(*,*) am1,am2,am1,am2
                 WRITE_MASTER(*,*) ig1,ig2,ig3,ig4
                 stop 'ERI calculated by libint failed'
               endif
               iibf=0
               do ibf=1,n1
                 do jbf=1,n2
                   do kbf=1,n1
                     do lbf=1,n2
                       iibf=iibf+1
                       integrals_cart(ibf,jbf,kbf,lbf) = integrals_cart(ibf,jbf,kbf,lbf) &
                                                        + integrals_libint(iibf) * shell(ishell)%coeff(ig1) * shell(jshell)%coeff(ig2) &
                                                                                 * shell(ishell)%coeff(ig3) * shell(jshell)%coeff(ig4)
                     enddo
                   enddo
                 enddo
               enddo

             enddo
           enddo
         enddo
       enddo

       do lbf=1,n2
         do kbf=1,n1
           do jbf=1,n2
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1,jbf,kbf,lbf) * cart_to_pure_norm(shell(ishell)%am)%matrix(1:n1,ibf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n2
         do kbf=1,n1
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2,kbf,lbf) * cart_to_pure_norm(shell(jshell)%am)%matrix(1:n2,jbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n2
         do kbf=1,ni
           do jbf=1,nj
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n1,lbf) * cart_to_pure_norm(shell(ishell)%am)%matrix(1:n1,kbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,nj
         do kbf=1,ni
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n2) * cart_to_pure_norm(shell(jshell)%am)%matrix(1:n2,lbf) )
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
     deallocate(integrals_libint)
     deallocate(alpha1,alpha2)

   enddo
 enddo

 WRITE_MASTER(*,*) 'Neglible shell pairs',nneglect,'/',neval

end subroutine identify_negligible_shellpair


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
subroutine transform_eri_basis_robust(nbf,nspin,c_matrix,eri_eigenstate)
 implicit none

 integer,intent(in) :: nspin,nbf
 real(dp),intent(in) :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: eri_eigenstate(nbf,nbf,nbf,nbf,nspin,nspin)
!=====
 integer :: ijspin,klspin
 integer :: ibf,jbf,kbf,lbf
 integer :: istate,jstate,kstate,lstate
 real(dp) :: eri_tmp(nbf,nbf,nbf,nbf)
!=====

 WRITE_MASTER(*,*) 'obtain the ERI in the eigenvector basis'

 eri_eigenstate(:,:,:,:,:,:)=0.0_dp
 do ijspin=1,nspin
   WRITE_MASTER(*,*) '== ijspin ',ijspin
   eri_tmp(:,:,:,:)=0.0_dp
   do lbf=1,nbf
     do kbf=1,nbf
       do jstate=1,nbf
         do istate=1,nbf

           do jbf=1,nbf
             do ibf=1,nbf
               eri_tmp(istate,jstate,kbf,lbf) = eri_tmp(istate,jstate,kbf,lbf) &
                  + eri(ibf,jbf,kbf,lbf) * c_matrix(ibf,istate,ijspin) * c_matrix(jbf,jstate,ijspin)
             enddo
           enddo

         enddo
       enddo
     enddo
   enddo
  
   do klspin=1,nspin
     do lstate=1,nbf
       do kstate=1,nbf
         do jstate=1,nbf
           do istate=1,nbf

             do lbf=1,nbf
               do kbf=1,nbf
                 eri_eigenstate(istate,jstate,kstate,lstate,ijspin,klspin) = eri_eigenstate(istate,jstate,kstate,lstate,ijspin,klspin) &
                                  + eri_tmp(istate,jstate,kbf,lbf) * c_matrix(kbf,kstate,klspin) * c_matrix(lbf,lstate,klspin)
               enddo
             enddo

           enddo
         enddo
       enddo
     enddo
   enddo

 enddo

end subroutine transform_eri_basis_robust


!=================================================================
subroutine transform_eri_basis_fast(nbf,nspin,c_matrix,eri_eigenstate)
 use m_timing
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none

 integer,intent(in) :: nspin,nbf
 real(dp),intent(in) :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: eri_eigenstate(nbf,nbf,nbf,nbf,nspin,nspin)
!=====
 integer :: ijspin,klspin
 integer :: ibf,jbf,kbf,lbf
 integer :: istate,jstate,kstate,lstate
 real(dp) :: eri_tmp1(nbf,nbf,nbf,nbf)
 real(dp) :: wtime
!=====

 WRITE_MASTER(*,*) 'obtain the ERI in the eigenvector basis'
 WRITE_MASTER(*,*) 'subroutine is order N^5'
! call start_clock(timing_basis_transform)

#ifdef _OPENMP
 wtime=OMP_get_wtime()
 WRITE_MASTER(*,*) 'The basis transform is using OPENMP'
#endif

 do ijspin=1,nspin
   WRITE_MASTER(*,*) '== ijspin ',ijspin
   eri_eigenstate(:,:,:,:,ijspin,:)=0.0_dp
   eri_tmp1(:,:,:,:)=0.0_dp

!$OMP PARALLEL DEFAULT(SHARED)

!$OMP DO SCHEDULE(STATIC)
   do lbf=1,nbf
     do kbf=1,nbf
       do jbf=1,nbf
         do istate=1,nbf

           do ibf=1,nbf
             eri_tmp1(istate,jbf,kbf,lbf) = eri_tmp1(istate,jbf,kbf,lbf) &
                                           +  eri(ibf,jbf,kbf,lbf) * c_matrix(ibf,istate,ijspin)
           enddo

         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP BARRIER


!$OMP DO SCHEDULE(STATIC)
   do lbf=1,nbf
     do kbf=1,nbf
       do jstate=1,nbf
         do istate=1,nbf

           eri_eigenstate(istate,jstate,kbf,lbf,ijspin,nspin) = DOT_PRODUCT( eri_tmp1(istate,:,kbf,lbf) , c_matrix(:,jstate,ijspin) )

         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP BARRIER

!$OMP END PARALLEL

  
   do klspin=1,nspin

!$OMP PARALLEL DEFAULT(SHARED)

!$OMP DO SCHEDULE(STATIC)
     do lbf=1,nbf
       do kstate=1,nbf
         do jstate=1,nbf
           do istate=1,nbf

             eri_tmp1(istate,jstate,kstate,lbf) = DOT_PRODUCT( eri_eigenstate(istate,jstate,:,lbf,ijspin,nspin) , c_matrix(:,kstate,klspin) ) 

           enddo
         enddo
       enddo
     enddo
!$OMP END DO
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
     do lstate=1,nbf
       do kstate=1,nbf
         do jstate=1,nbf
           do istate=1,nbf

             eri_eigenstate(istate,jstate,kstate,lstate,ijspin,klspin) = DOT_PRODUCT( eri_tmp1(istate,jstate,kstate,:) , c_matrix(:,lstate,klspin) )

           enddo
         enddo
       enddo
     enddo
!$OMP END DO
!$OMP BARRIER

!$OMP END PARALLEL

   enddo !klspin

 enddo !ijspin

#ifdef _OPENMP
  WRITE_MASTER(*,*) 'time (s)', OMP_get_wtime()-wtime
#endif

! call stop_clock(timing_basis_transform)
 WRITE_MASTER(*,*) 'ERI in the eigenvector basis obtained'
 WRITE_MASTER(*,*)

end subroutine transform_eri_basis_fast


!=================================================================
subroutine transform_eri_basis_lowmem(nspin,c_matrix,istate,ijspin,eri_eigenstate_i)
 use m_timing
 implicit none

 integer,intent(in)   :: nspin,istate,ijspin
 real(dp),intent(in)  :: c_matrix(nbf_eri,nbf_eri,nspin)
 real(dp),intent(out) :: eri_eigenstate_i(nbf_eri,nbf_eri,nbf_eri,nspin)
!=====
 integer              :: klspin
 integer              :: ibf,jbf,kbf,lbf
 integer              :: jstate,kstate,lstate
 real(dp)             :: eri_tmp3(nbf_eri,nbf_eri,nbf_eri)
 real(dp)             :: wtime
!=====

! call start_clock(timing_basis_transform)

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

! call stop_clock(timing_basis_transform)

end subroutine transform_eri_basis_lowmem


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
