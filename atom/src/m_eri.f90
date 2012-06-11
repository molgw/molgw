!=========================================================================
#include "macros.h"
!=========================================================================
module m_eri
 use m_definitions
 use m_mpi
 use m_basis_set

 private
 public :: eri,allocate_eri,allocate_eri_eigen,deallocate_eri,calculate_eri,transform_eri_basis_fast,transform_eri_basis_lowmem, &
           eri_lr,allocate_eri_lr,deallocate_eri_lr,negligible_eri,&
           BUFFER1,BUFFER2
 public :: index_prod

 integer,parameter :: BUFFER1 = 1
 integer,parameter :: BUFFER2 = 2
 !
 ! max length of a record in the ERI file
 integer,parameter :: line_length=1000

 real(prec_eri),allocatable :: eri_buffer(:)
 real(prec_eri),allocatable :: eri_buffer_lr(:)

#ifdef LOW_MEMORY3
 integer                    :: nsize_sparse
 real(prec_eri),allocatable :: eri_buffer_sparse(:)
 integer*2,allocatable      :: index_sparse(:,:)
#endif

 integer                    :: nbf_eri                ! local copy of nbf
 integer                    :: nsize                  ! size of the eri_buffer array
 integer                    :: nsize1                 ! number of independent pairs (i,j) with i<=j

 real(prec_eri),allocatable :: eri_eigen_buffer(:,:,:)

contains

!=========================================================================
subroutine allocate_eri(nbf)
 implicit none
!===== 
 integer,intent(in) :: nbf
!===== 
 integer            :: info
!===== 

 nbf_eri = nbf
#if LOW_MEMORY3 || LOW_MEMORY2
 WRITE_MASTER(*,'(/,a)') ' Symmetrized ERI stored'
 nsize1  = index_prod(nbf_eri,nbf_eri) 
 nsize   = index_eri(nbf_eri,nbf_eri,nbf_eri,nbf_eri)
#elif LOW_MEMORY1
 WRITE_MASTER(*,'(/,a)') ' Semi-symmetrized ERI stored (4 symmetries)'
 nsize1  = index_prod(nbf_eri,nbf_eri) 
 nsize   = nsize1*get_ntask()
#else
 WRITE_MASTER(*,'(/,a)') ' All ERI are stored'
 nsize1  = nbf_eri**2
 nsize   = nsize1**2
#endif

 WRITE_MASTER(*,*) 'number of integrals to be stored:',nsize
 WRITE_MASTER(*,*) 'max index size',HUGE(nsize)
 if(nsize<1) stop'too many integrals to be stored'

 allocate(eri_buffer(nsize),stat=info)
 if(REAL(nsize,dp)*prec_eri > 1024**3 ) then
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the ERI array: ',REAL(nsize,dp)*prec_eri/1024**3,' [Gb] / proc'
 else
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the ERI array: ',REAL(nsize,dp)*prec_eri/1024**2,' [Mb] / proc'
 endif
 if(info==0) then
   WRITE_MASTER(*,*) 'success'
 else
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif

 eri_buffer(:) = 0.0_dp

end subroutine allocate_eri

!=========================================================================
subroutine deallocate_eri()
 implicit none
!=====

 if(allocated(eri_buffer))        deallocate(eri_buffer)
#ifdef LOW_MEMORY3
 if(allocated(eri_buffer_sparse)) deallocate(eri_buffer_sparse)
 if(allocated(index_sparse))      deallocate(index_sparse)
#endif

end subroutine deallocate_eri

!=========================================================================
subroutine allocate_eri_lr(nbf)
 implicit none
!===== 
 integer,intent(in) :: nbf
!===== 
 integer            :: info
!===== 

 nbf_eri = nbf
#if LOW_MEMORY3 || LOW_MEMORY2
 WRITE_MASTER(*,'(/,a)') ' Symmetrized ERI stored'
 nsize1  = index_prod(nbf_eri,nbf_eri) 
 nsize   = index_eri(nbf_eri,nbf_eri,nbf_eri,nbf_eri)
#else
 WRITE_MASTER(*,'(/,a)') ' All ERI are stored'
 nsize1  = nbf_eri**2
 nsize   = nsize1**2
#endif

 allocate(eri_buffer_lr(nsize),stat=info)
 if(REAL(nsize,dp)*prec_eri > 1024**3 ) then
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the Long-Range ERI array: ',REAL(nsize,dp)*prec_eri/1024**3,' [Gb]'
 else
   WRITE_MASTER(*,'(a,f10.3,a)') ' Allocating the Long-Range ERI array: ',REAL(nsize,dp)*prec_eri/1024**2,' [Mb]'
 endif
 if(info==0) then
   WRITE_MASTER(*,*) 'success'
 else
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif

 eri_buffer_lr(:) = 0.0_dp

end subroutine allocate_eri_lr

!=========================================================================
subroutine deallocate_eri_lr()
 implicit none
!=====

 if(allocated(eri_buffer_lr))    deallocate(eri_buffer_lr)
#ifdef LOW_MEMORY3
 if(allocated(eri_buffer_sparse)) deallocate(eri_buffer_sparse)
 if(allocated(index_sparse))      deallocate(index_sparse)
#endif

end subroutine deallocate_eri_lr


!=========================================================================
function index_prod(ibf,jbf)
 implicit none
 integer,intent(in) :: ibf,jbf
 integer            :: index_prod
!=====
 integer            :: jmin,imax
!=====

 imax=MAX(ibf,jbf)
 jmin=MIN(ibf,jbf)
 index_prod = (jmin-1)*nbf_eri - (jmin-1)*(jmin-2)/2 + imax-jmin+1

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

#ifdef LOW_MEMORY1
 index_kl  = get_task_number(kbf,lbf)
 index_eri = index_ij + ( index_kl - 1 ) * nsize1
#else
 index_kl = index_prod(kbf,lbf)

 ijmax=MAX(index_ij,index_kl)
 klmin=MIN(index_ij,index_kl)

 index_eri = (klmin-1)*nsize1 - (klmin-1)*(klmin-2)/2 + ijmax-klmin+1
#endif

end function index_eri

!=========================================================================
function eri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri
!=====
 integer            :: ibuffer_sparse,index_ijkl
 integer            :: i1,i2,i3,i4
!=====

#if LOW_MEMORY3

 eri = 0.0_dp

 if( index_prod(ibf,jbf) >=  index_prod(kbf,lbf) ) then
   i1=MAX(ibf,jbf)
   i2=MIN(ibf,jbf)
   i3=MAX(kbf,lbf)
   i4=MIN(kbf,lbf)
 else
   i3=MAX(ibf,jbf)
   i4=MIN(ibf,jbf)
   i1=MAX(kbf,lbf)
   i2=MIN(kbf,lbf)
 endif

 do ibuffer_sparse=1,nsize_sparse

   if( index_sparse(1,ibuffer_sparse) == i1 .AND. &
       index_sparse(2,ibuffer_sparse) == i2 .AND. &
       index_sparse(3,ibuffer_sparse) == i3 .AND. &
       index_sparse(4,ibuffer_sparse) == i4 ) then
     eri = eri_buffer_sparse(ibuffer_sparse)
     exit
   endif

 enddo

#elif LOW_MEMORY2 || LOW_MEMORY1
 eri = eri_buffer(index_eri(ibf,jbf,kbf,lbf))
#else
 eri = eri_buffer(ibf+(jbf-1)*nbf_eri+(kbf-1)*nbf_eri**2+(lbf-1)*nbf_eri**3)
#endif

end function eri

!=========================================================================
function eri_lr(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_lr
!=====
 integer            :: ibuffer_sparse,index_ijkl
 integer            :: i1,i2,i3,i4
!=====

#if LOW_MEMORY2 || LOW_MEMORY1
 eri_lr = eri_buffer_lr(index_eri(ibf,jbf,kbf,lbf))
#else
 eri_lr = eri_buffer_lr(ibf+(jbf-1)*nbf_eri+(kbf-1)*nbf_eri**2+(lbf-1)*nbf_eri**3)
#endif

end function eri_lr

!=========================================================================
subroutine calculate_eri(print_volume,basis,rcut,which_buffer)
 implicit none
 integer,intent(in)           :: print_volume
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: which_buffer
!=====

 if( .NOT. read_eri(rcut) ) call do_calculate_eri(basis,rcut,which_buffer)

 if(print_volume>100) then
   call dump_out_eri(rcut)
 endif

end subroutine calculate_eri

!=========================================================================
subroutine do_calculate_eri(basis,rcut,which_buffer)
 use ISO_C_BINDING
 use m_tools,only: boys_function
 use m_timing
#ifdef OPENMP
 use omp_lib
#endif
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: which_buffer
!=====
 integer,parameter            :: NSHELLMAX=400
 integer,parameter            :: NMEMBER=28
 logical,parameter            :: DEBUG=.TRUE.
 integer                      :: info
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: ig,jg,kg,lg
 integer                      :: n1,n2,n3,n4
 integer                      :: ni,nj,nk,nl
 integer                      :: ii,i,j,k,l
 integer                      :: nshell,nint_shell
 integer                      :: ishell,jshell,kshell,lshell
 type(basis_function),pointer :: bf_current_i,bf_current_j,bf_current_k,bf_current_l
 type(gaussian),pointer       :: g_current
 integer                      :: nint_gaussian,igaussian
 real(dp),allocatable         :: int_gaussian(:)
 logical                      :: shell_already_exists
 logical                      :: need_calculation
 integer                      :: iint,nint_tot
 integer                      :: imember,jmember,kmember,lmember
 integer                      :: iindex_in_the_shell,jindex_in_the_shell,kindex_in_the_shell,lindex_in_the_shell
 integer                      :: index_integral
 integer                      :: index_tmp
 integer                      :: ami,amj,amk,aml
 real(dp),allocatable         :: int_tmp(:,:,:,:)
 real(dp)                     :: zeta_12,zeta_34,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 integer                      :: am1f,am2f,am3f,am4f
!=====
 type shell_type
   integer  :: am
   real(dp) :: alpha
   real(dp) :: x0(3)
   integer  :: nmember              !
   integer  :: index_bf(NMEMBER)    ! correspondance from shell to basis functions and primitive gaussians
   integer  :: index_g(NMEMBER)     !
 end type shell_type
 type(shell_type)             :: shell(NSHELLMAX)
!=====
! variables used to call C++ 
 integer(C_INT),external      :: calculate_integral
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: alpha1,alpha2,alpha3,alpha4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: omega_range
!=====


 WRITE_MASTER(*,'(/,a)') ' Calculate and store all the Electron Repulsion Integrals (ERI)'

 if( rcut > 1.0e-6_dp ) then
   omega_range = 1.0_dp / rcut
   WRITE_MASTER(*,'(a40,x,f9.4)') ' Long-Range only integrals with rcut=',rcut
   WRITE_MASTER(*,'(a40,x,f9.4)') ' or omega=',omega_range
 else 
   omega_range = 1.0e6_dp
 endif

 nint_gaussian=0
 ishell=0
 do ibf=1,basis%nbf
   bf_current_i => basis%bf(ibf)

   do ig=1,bf_current_i%ngaussian
     g_current => bf_current_i%g(ig)

     !
     ! check if the shell has already been created
     shell_already_exists=.FALSE.
     do jshell=1,ishell
       if( g_current%am == shell(jshell)%am .AND. ABS( g_current%alpha - shell(jshell)%alpha ) < 1.d-8 &
           .AND.  ALL( ABS( g_current%x0(:) - shell(jshell)%x0(:) ) < 1.d-8 ) ) then

         shell_already_exists=.TRUE.
         shell(jshell)%nmember = shell(jshell)%nmember + 1
         if( shell(jshell)%nmember > NMEMBER ) stop'NMEMBER hard coded parameter is too small!'
         shell(jshell)%index_bf(shell(jshell)%nmember) = ibf
         shell(jshell)%index_g (shell(jshell)%nmember) = ig

         exit
       endif
     enddo

     if(.NOT.shell_already_exists) then
       ishell=ishell+1
       if(ishell > NSHELLMAX) stop'NSHELLMAX internal parameter is too small'

       shell(ishell)%am      = g_current%am
       shell(ishell)%alpha   = g_current%alpha
       shell(ishell)%x0(:)   = g_current%x0(:)
       shell(ishell)%nmember = 1
       shell(ishell)%index_bf(shell(ishell)%nmember) = ibf
       shell(ishell)%index_g (shell(ishell)%nmember) = ig

     endif

   enddo
 enddo
 nshell=ishell
 
 call start_clock(timing_tmp2)
 WRITE_MASTER(*,*) 'number of shells',nshell

 !
 ! (ij||kl)

!!!!    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iindex_in_the_shell,jindex_in_the_shell,kindex_in_the_shell,lindex_in_the_shell,&
!!!!    !$OMP&                                 ig,jg,kg,lg,ibf,jbf,kbf,lbf,index_integral,index_tmp,&
!!!!    !$OMP&                                 ami,amj,amk,aml,am1,am2,am3,am4,alpha1,alpha2,alpha3,alpha4,&
!!!!    !$OMP&                                 n1,n2,n3,n4,ni,nj,nk,nl,nint_shell,int_shell,int_tmp,ii,info,i,j,k,l,ishell,jshell,kshell,lshell)
!!!!    
!!!!    
!!!!    !$OMP DO SCHEDULE(STATIC)

 do lshell=1,nshell
   do kshell=1,nshell

! work section HERE

     need_calculation = .FALSE.
     do lmember=1,shell(lshell)%nmember
       lbf = shell(lshell)%index_bf(lmember)
       do kmember=1,shell(kshell)%nmember
         kbf = shell(kshell)%index_bf(kmember)
         if( is_my_task(kbf,lbf) ) need_calculation=.TRUE.
       enddo
     enddo

     if( .NOT. need_calculation ) cycle


! end of work section

     do jshell=1,nshell
       do ishell=1,nshell

         call start_clock(timing_tmp3)

         ami = shell(ishell)%am
         amj = shell(jshell)%am
         amk = shell(kshell)%am
         aml = shell(lshell)%am

         !
         ! Order the angular momenta so that libint is pleased
         ! 1) am3+am4 >= am1+am2
         ! 2) am3>=am4
         ! 3) am1>=am2
         if(amk+aml>=ami+amj) then
           if(amk>=aml) then
             am3 = shell(kshell)%am
             am4 = shell(lshell)%am
             alpha3 = shell(kshell)%alpha
             alpha4 = shell(lshell)%alpha
             x03(:) = shell(kshell)%x0(:)
             x04(:) = shell(lshell)%x0(:)
           else
             am3 = shell(lshell)%am
             am4 = shell(kshell)%am
             alpha3 = shell(lshell)%alpha
             alpha4 = shell(kshell)%alpha
             x03(:) = shell(lshell)%x0(:)
             x04(:) = shell(kshell)%x0(:)
           endif
           if(ami>=amj) then
             am1 = shell(ishell)%am
             am2 = shell(jshell)%am
             alpha1 = shell(ishell)%alpha
             alpha2 = shell(jshell)%alpha
             x01(:) = shell(ishell)%x0(:)
             x02(:) = shell(jshell)%x0(:)
           else
             am1 = shell(jshell)%am
             am2 = shell(ishell)%am
             alpha1 = shell(jshell)%alpha
             alpha2 = shell(ishell)%alpha
             x01(:) = shell(jshell)%x0(:)
             x02(:) = shell(ishell)%x0(:)
           endif
         else
           if(amk>=aml) then
             am1 = shell(kshell)%am
             am2 = shell(lshell)%am
             alpha1 = shell(kshell)%alpha
             alpha2 = shell(lshell)%alpha
             x01(:) = shell(kshell)%x0(:)
             x02(:) = shell(lshell)%x0(:)
           else
             am1 = shell(lshell)%am
             am2 = shell(kshell)%am
             alpha1 = shell(lshell)%alpha
             alpha2 = shell(kshell)%alpha
             x01(:) = shell(lshell)%x0(:)
             x02(:) = shell(kshell)%x0(:)
           endif
           if(ami>=amj) then
             am3 = shell(ishell)%am
             am4 = shell(jshell)%am
             alpha3 = shell(ishell)%alpha
             alpha4 = shell(jshell)%alpha
             x03(:) = shell(ishell)%x0(:)
             x04(:) = shell(jshell)%x0(:)
           else
             am3 = shell(jshell)%am
             am4 = shell(ishell)%am
             alpha3 = shell(jshell)%alpha
             alpha4 = shell(ishell)%alpha
             x03(:) = shell(jshell)%x0(:)
             x04(:) = shell(ishell)%x0(:)
           endif
         endif

         ni = number_basis_function_am( ami )
         nj = number_basis_function_am( amj )
         nk = number_basis_function_am( amk )
         nl = number_basis_function_am( aml )

         am1f=am1
         am2f=am2
         am3f=am3
         am4f=am4
         n1 = number_basis_function_am( am1f )
         n2 = number_basis_function_am( am2f )
         n3 = number_basis_function_am( am3f )
         n4 = number_basis_function_am( am4f )
         
         nint_shell = n1 * n2 * n3 *n4
         if( allocated(int_shell) ) deallocate( int_shell )
         allocate( int_shell( nint_shell ) )

         !
         ! If all the angular momentum are S type, do the calculation immediately
         ! else call the libint library
         !
         if(am1+am2+am3+am4==0) then

!           int_shell(1) = 2.0_dp*pi**(2.5_dp) &
!                         / ( (alpha1+alpha2)*(alpha3+alpha4)*SQRT(alpha1+alpha2+alpha3+alpha4) )

           zeta_12 = alpha1 + alpha2
           zeta_34 = alpha3 + alpha4
           p(:) = ( alpha1 * x01(:) + alpha2 * x02(:) ) / zeta_12 
           q(:) = ( alpha3 * x03(:) + alpha4 * x04(:) ) / zeta_34 
           !
           ! Full range or long-range only integrals
           if( rcut < 1.0e-6_dp ) then
             rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
             rho1 = rho
           else
             rho  = zeta_12 * zeta_34 * omega_range**2 / ( zeta_12*omega_range**2 + zeta_34*omega_range**2 + zeta_12*zeta_34 )
             rho1 = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
           endif
           tt = rho * SUM( (p(:)-q(:))**2 )
           call boys_function(f0t(0),0,tt)
           int_shell(1) = 2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                 / zeta_12 * exp( -alpha1*alpha2/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                 / zeta_34 * exp( -alpha3*alpha4/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) &
                 * SQRT( rho / rho1 )

         else

           info=calculate_integral(omega_range,&
                                   am1,am2,am3,am4,alpha1,alpha2,alpha3,alpha4,&
                                   x01(1),x01(2),x01(3),&
                                   x02(1),x02(2),x02(3),&
                                   x03(1),x03(2),x03(3),&
                                   x04(1),x04(2),x04(3),&
                                   int_shell(1))
           if(info/=0) then
             WRITE_MASTER(*,*) am1,am2,am3,am4
             stop 'ERI calculated by libint failed'
           endif

         endif

         !
         ! Reorder the integrals in the original order
         ! if the ordering has been changed to please libint
         allocate(int_tmp(n1,n2,n3,n4))
         ii=0
         do i=1,n1
           do j=1,n2
             do k=1,n3
               do l=1,n4
                 ii=ii+1
                 int_tmp(i,j,k,l) = int_shell(ii)
               enddo
             enddo
           enddo
         enddo

         ii=0
         if(amk+aml>=ami+amj) then
           if(amk>=aml .AND. ami>=amj) then
             do i=1,n1
               do j=1,n2
                 do k=1,n3
                   do l=1,n4
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           else if(amk<aml .AND. ami>=amj) then
             do i=1,n1
               do j=1,n2
                 do l=1,n4
                   do k=1,n3
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           else if(amk>=aml .AND. ami<amj) then
             do j=1,n2
               do i=1,n1
                 do k=1,n3
                   do l=1,n4
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           else 
             do j=1,n2
               do i=1,n1
                 do l=1,n4
                   do k=1,n3
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           endif
         else ! amk+aml<ami+amj
           if(amk>=aml .AND. ami>=amj) then
             do k=1,n3
               do l=1,n4
                 do i=1,n1
                   do j=1,n2
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           else if(amk<aml .AND. ami>=amj) then
             do k=1,n3
               do l=1,n4
                 do j=1,n2
                   do i=1,n1
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           else if(amk>=aml .AND. ami<amj) then
             do l=1,n4
               do k=1,n3
                 do i=1,n1
                   do j=1,n2
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           else
             do l=1,n4
               do k=1,n3
                 do j=1,n2
                   do i=1,n1
                     ii=ii+1
                     int_shell(ii) = int_tmp(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
           endif
         endif
         deallocate(int_tmp)


         call stop_clock(timing_tmp3)
         call start_clock(timing_tmp4)

         !
         ! different storage according to which_buffer
         !
 
         if( which_buffer == BUFFER1 ) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iindex_in_the_shell,jindex_in_the_shell,kindex_in_the_shell,lindex_in_the_shell,&
!$OMP&     ig,jg,kg,lg,ibf,jbf,kbf,lbf,index_integral,index_tmp )

!$OMP DO SCHEDULE(STATIC)
           do lmember=1,shell(lshell)%nmember
             lbf = shell(lshell)%index_bf(lmember)
             lg  = shell(lshell)%index_g (lmember)
             lindex_in_the_shell = libint_ordering(basis%bf(lbf)%nx,basis%bf(lbf)%ny,basis%bf(lbf)%nz)
             do kmember=1,shell(kshell)%nmember
               kbf = shell(kshell)%index_bf(kmember)
               kg  = shell(kshell)%index_g (kmember)
               kindex_in_the_shell = libint_ordering(basis%bf(kbf)%nx,basis%bf(kbf)%ny,basis%bf(kbf)%nz)
               do jmember=1,shell(jshell)%nmember
                 jbf = shell(jshell)%index_bf(jmember)
                 jg  = shell(jshell)%index_g (jmember)
                 jindex_in_the_shell = libint_ordering(basis%bf(jbf)%nx,basis%bf(jbf)%ny,basis%bf(jbf)%nz)
                 do imember=1,shell(ishell)%nmember
                   ibf = shell(ishell)%index_bf(imember)
                   ig  = shell(ishell)%index_g (imember)
                   iindex_in_the_shell = libint_ordering(basis%bf(ibf)%nx,basis%bf(ibf)%ny,basis%bf(ibf)%nz)

                   index_integral = lindex_in_the_shell + (kindex_in_the_shell-1)*nl &
                                   +(jindex_in_the_shell-1)*nl*nk + (iindex_in_the_shell-1)*nl*nk*nj

#if LOW_MEMORY3 || LOW_MEMORY2 || LOW_MEMORY1
                   if(ibf<jbf) cycle
                   if(kbf<lbf) cycle
#ifndef LOW_MEMORY1
                   if(index_prod(ibf,jbf)<index_prod(kbf,lbf)) cycle
#endif
#ifdef MPI
                   if( .NOT. is_my_task(kbf,lbf) ) cycle
#endif
                   index_tmp=index_eri(ibf,jbf,kbf,lbf)
#else
                   index_tmp=ibf+(jbf-1)*nbf_eri+(kbf-1)*nbf_eri**2+(lbf-1)*nbf_eri**3
#endif
                   eri_buffer(index_tmp) = eri_buffer(index_tmp) &
                             + basis%bf(ibf)%coeff(ig) *  basis%bf(ibf)%g(ig)%norm_factor &
                             * basis%bf(jbf)%coeff(jg) *  basis%bf(jbf)%g(jg)%norm_factor &
                             * basis%bf(kbf)%coeff(kg) *  basis%bf(kbf)%g(kg)%norm_factor &
                             * basis%bf(lbf)%coeff(lg) *  basis%bf(lbf)%g(lg)%norm_factor &
                               * int_shell(index_integral) 

                 enddo
               enddo
             enddo
           enddo
!$OMP END DO
!$OMP END PARALLEL


         else    ! Store the result in the arry eri_buffer_lr


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iindex_in_the_shell,jindex_in_the_shell,kindex_in_the_shell,lindex_in_the_shell,&
!$OMP&     ig,jg,kg,lg,ibf,jbf,kbf,lbf,index_integral,index_tmp )

!$OMP DO SCHEDULE(STATIC)
           do lmember=1,shell(lshell)%nmember
             lbf = shell(lshell)%index_bf(lmember)
             lg  = shell(lshell)%index_g (lmember)
             lindex_in_the_shell = libint_ordering(basis%bf(lbf)%nx,basis%bf(lbf)%ny,basis%bf(lbf)%nz)
             do kmember=1,shell(kshell)%nmember
               kbf = shell(kshell)%index_bf(kmember)
               kg  = shell(kshell)%index_g (kmember)
               kindex_in_the_shell = libint_ordering(basis%bf(kbf)%nx,basis%bf(kbf)%ny,basis%bf(kbf)%nz)
               do jmember=1,shell(jshell)%nmember
                 jbf = shell(jshell)%index_bf(jmember)
                 jg  = shell(jshell)%index_g (jmember)
                 jindex_in_the_shell = libint_ordering(basis%bf(jbf)%nx,basis%bf(jbf)%ny,basis%bf(jbf)%nz)
                 do imember=1,shell(ishell)%nmember
                   ibf = shell(ishell)%index_bf(imember)
                   ig  = shell(ishell)%index_g (imember)
                   iindex_in_the_shell = libint_ordering(basis%bf(ibf)%nx,basis%bf(ibf)%ny,basis%bf(ibf)%nz)

                   index_integral = lindex_in_the_shell + (kindex_in_the_shell-1)*nl &
                                   +(jindex_in_the_shell-1)*nl*nk + (iindex_in_the_shell-1)*nl*nk*nj

#if LOW_MEMORY3 || LOW_MEMORY2 || LOW_MEMORY1
                   if(ibf<jbf) cycle
                   if(kbf<lbf) cycle
#ifndef LOW_MEMORY1
                   if(index_prod(ibf,jbf)<index_prod(kbf,lbf)) cycle
#endif
#ifdef MPI
                   if( .NOT. is_my_task(kbf,lbf) ) cycle
#endif

                   index_tmp=index_eri(ibf,jbf,kbf,lbf)
#else
                   index_tmp=ibf+(jbf-1)*nbf_eri+(kbf-1)*nbf_eri**2+(lbf-1)*nbf_eri**3
#endif
                   eri_buffer_lr(index_tmp) = eri_buffer_lr(index_tmp) &
                             + basis%bf(ibf)%coeff(ig) *  basis%bf(ibf)%g(ig)%norm_factor &
                             * basis%bf(jbf)%coeff(jg) *  basis%bf(jbf)%g(jg)%norm_factor &
                             * basis%bf(kbf)%coeff(kg) *  basis%bf(kbf)%g(kg)%norm_factor &
                             * basis%bf(lbf)%coeff(lg) *  basis%bf(lbf)%g(lg)%norm_factor &
                               * int_shell(index_integral) 

                 enddo
               enddo
             enddo
           enddo
!$OMP END DO
!$OMP END PARALLEL

         endif      ! if full-range or long-range


         call stop_clock(timing_tmp4)

       enddo
     enddo
   enddo
 enddo
!!!!       !$OMP END DO
 if( allocated(int_shell) ) deallocate( int_shell )
!!!!       !$OMP END PARALLEL

 call stop_clock(timing_tmp2)
 WRITE_MASTER(*,*) 'Done!'
 WRITE_MASTER(*,*)

end subroutine do_calculate_eri


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
   if(nx==1) libint_ordering=1
   if(ny==1) libint_ordering=2
   if(nz==1) libint_ordering=3
 case(2)
   if(nx==2          ) libint_ordering=1
   if(nx==1.AND.ny==1) libint_ordering=2
   if(nx==1.AND.nz==1) libint_ordering=3
   if(ny==2          ) libint_ordering=4
   if(ny==1.AND.nz==1) libint_ordering=5
   if(nz==2          ) libint_ordering=6
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
#ifdef OPENMP
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

#ifdef OPENMP
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

#ifdef OPENMP
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
 integer             :: ibuffer,ibuffer_sparse
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

#ifdef LOW_MEMORY3
 nsize_sparse = nsize - icount
 allocate(eri_buffer_sparse(nsize_sparse))
 allocate(index_sparse(4,nsize_sparse))

 if(REAL(nsize_sparse,dp)*prec_eri > 1024**3 ) then
   WRITE_MASTER(*,'(a,f8.3,a)') ' Allocating the sparse ERI array: ',REAL(nsize_sparse,dp)*prec_eri/1024**3,' [Gb]'
 else
   WRITE_MASTER(*,'(a,f8.3,a)') ' Allocating the sparse ERI array: ',REAL(nsize_sparse,dp)*prec_eri/1024**2,' [Mb]'
 endif

 ibuffer_sparse=0
 ibuffer=0
 do lbf=1,nbf_eri
   do kbf=lbf,nbf_eri
     do jbf=1,nbf_eri
       do ibf=jbf,nbf_eri

         if( index_prod(ibf,jbf) >= index_prod(kbf,lbf) ) then
           ibuffer = ibuffer + 1
           if( ABS( eri_buffer(ibuffer) ) > tol ) then
             ibuffer_sparse = ibuffer_sparse + 1
             eri_buffer_sparse(ibuffer_sparse) = eri_buffer(ibuffer)
             index_sparse(1,ibuffer_sparse) = ibf
             index_sparse(2,ibuffer_sparse) = jbf
             index_sparse(3,ibuffer_sparse) = kbf
             index_sparse(4,ibuffer_sparse) = lbf

           endif

         endif

       enddo
     enddo
   enddo
 enddo

 deallocate(eri_buffer)
#endif


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
