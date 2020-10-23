!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the method to prepare and store the 2-, 3-, and 4-center Coulomb integrals
!
!=========================================================================
module m_eri
 use m_definitions
 use m_mpi
 use m_memory
 use m_warning
 use m_scalapack
 use m_basis_set
 use m_timing
 use m_cart_to_pure
 use m_libint_tools
 use m_scalapack
 use m_inputparam,only: integral_level

 logical,parameter :: eri_pair_major = .TRUE.

 real(dp),parameter,public :: TOO_LOW_EIGENVAL=1.0e-6_dp

 !
 ! max length of a record in the ERI file
 integer,parameter,private :: line_length = 1000

 real(dp),private           :: TOL_INT


 real(dp),allocatable,public :: eri_4center(:)
 real(dp),allocatable,public :: eri_4center_lr(:)
 real(dp),allocatable,public :: eri_3center(:,:)         ! global size: npair x nauxil_2center
 real(dp),allocatable,public :: eri_3center_lr(:,:)      ! global size: npair x nauxil_2center_lr


 logical,protected,allocatable :: negligible_shellpair(:,:)
 integer,allocatable,private   :: index_pair_1d(:)
 integer,protected,allocatable :: index_basis(:,:)
 integer,protected,allocatable :: index_shellpair(:,:)
 integer,protected             :: nshellpair

 integer,private,allocatable   :: shell_bf(:)


 integer,private   :: nbf_eri               ! local copy of nbf
 ! 4-byte integers are too small to index the 4-center Coulomb integrals using no-RI/AE/5Z for 3d transition metals.
 ! Use 8-byte integers instead
 integer(kind=int8),protected :: nint_4center    ! size of the eri_4center array
 integer,protected            :: npair           ! number of independent pairs (i,j) with i<=j

 integer,public    :: nauxil_2center     ! size of the 2-center matrix
 integer,public    :: nauxil_2center_lr  ! size of the 2-center LR matrix

 integer,protected :: nauxil_3center     ! size of the 3-center matrix
                                         ! may differ from the total number of 3-center integrals due to
                                         ! data distribution
 integer,protected :: nauxil_3center_lr  ! size of the 3-center matrix
                                         ! may differ from the total number of 3-center integrals due to
                                         ! data distribution

 integer,public    :: desc_eri3(NDEL)
 integer,public    :: desc_eri3_lr(NDEL)

! Parallelization information for the auxiliary basis
 integer,allocatable,protected :: iproc_ibf_auxil(:)
 integer,allocatable,protected :: ibf_auxil_g(:)       ! auxil bf index from local to global
 integer,allocatable,protected :: ibf_auxil_l(:)       ! auxil bf index from global to local

! Parallelization information for the auxiliary basis (LR part)
 integer,allocatable,protected :: iproc_ibf_auxil_lr(:)
 integer,allocatable,protected :: ibf_auxil_g_lr(:)
 integer,allocatable,protected :: ibf_auxil_l_lr(:)


contains


!=========================================================================
subroutine prepare_eri(basis)
 implicit none
!=====
 type(basis_set),intent(in) :: basis
!=====
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
   TOL_INT = -1.0_dp
 case default
   call die('integration quality not recognized')
 end select
 write(stdout,'(/,a,es9.2)') ' Tolerance on integrals set to ',TOL_INT


 if(.NOT.ALLOCATED(negligible_shellpair)) then
   call setup_shell_index(basis)
   allocate(negligible_shellpair(basis%nshell,basis%nshell))
   call identify_negligible_shellpair(basis)
   call setup_shellpair(basis)
   call setup_basispair(basis)
 endif


 ! Carefully perform this calculation with 8-byte integers since the result can be very very large
 nint_4center = ( INT(npair,KIND=int8) * ( INT(npair,KIND=int8) + 1_int8 ) ) / 2_int8


end subroutine prepare_eri


!=========================================================================
subroutine deallocate_eri_4center()
 implicit none
!=====

 if(ALLOCATED(eri_4center)) then
   call clean_deallocate('4-center integrals',eri_4center)
 endif

end subroutine deallocate_eri_4center


!=========================================================================
subroutine deallocate_eri_4center_lr()
 implicit none
!=====

 if(ALLOCATED(eri_4center_lr)) then
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif

end subroutine deallocate_eri_4center_lr


!=========================================================================
subroutine deallocate_eri()
 implicit none

!=====

 if(ALLOCATED(eri_4center)) then
   call clean_deallocate('4-center integrals',eri_4center)
 endif
 if(ALLOCATED(eri_4center_lr)) then
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif
 if(ALLOCATED(negligible_shellpair))   deallocate(negligible_shellpair)
 if(ALLOCATED(index_shellpair))        deallocate(index_shellpair)
 if(ALLOCATED(index_pair_1d)) call clean_deallocate('index pair',index_pair_1d)
 if(ALLOCATED(index_basis))   call clean_deallocate('index basis',index_basis)

 if(ALLOCATED(shell_bf))              deallocate(shell_bf)


end subroutine deallocate_eri


!=========================================================================
pure function index_eri(ibf,jbf,kbf,lbf)
 implicit none

 integer,intent(in) :: ibf,jbf,kbf,lbf
 integer(kind=int8) :: index_eri
!=====
 integer(kind=int8) :: klmin,ijmax
 integer(kind=int8) :: index_ij,index_kl
!=====

 index_ij = index_pair(ibf,jbf)
 index_kl = index_pair(kbf,lbf)

 ijmax = MAX(index_ij,index_kl)
 klmin = MIN(index_ij,index_kl)

 index_eri = (klmin-1) * npair - ( (klmin-1) * (klmin-2) ) / 2 + ijmax - klmin + 1


end function index_eri


!=========================================================================
pure function index_pair(ibf,jbf)
 implicit none

 integer,intent(in) :: ibf,jbf
 integer            :: index_pair
!=====
 integer            :: ijmin,ijmax
!=====

 ijmax = MAX(ibf,jbf)
 ijmin = MIN(ibf,jbf)

 index_pair = (ijmin-1) * nbf_eri - ( (ijmin-2) * (ijmin-1) ) / 2  + ijmax - ijmin + 1
 index_pair = index_pair_1d(index_pair)

end function index_pair


!=========================================================================
elemental function eri(ibf,jbf,kbf,lbf)
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
subroutine setup_shell_index(basis)
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
!=====

 allocate(shell_bf(basis%nbf))
 shell_bf(:) = basis%bff(:)%shell_index

end subroutine setup_shell_index


!=========================================================================
subroutine setup_basispair(basis)
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer :: ipair
 integer :: ibf,jbf,ijbf,ijmax,ijmin
 integer :: ishell,jshell,ijshellpair
 integer :: ami,amj,ni,nj
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
 ! Specific ordering where the basis functions in a shell pair have contiguous indexes
 !
 index_pair_1d(:) = 0
 ipair = 0
 do ijshellpair=1,nshellpair
   ishell = index_shellpair(1,ijshellpair)
   jshell = index_shellpair(2,ijshellpair)
   ami = basis%shell(ishell)%am
   amj = basis%shell(jshell)%am
   ni = number_basis_function_am( basis%gaussian_type , ami )
   nj = number_basis_function_am( basis%gaussian_type , amj )

   if( ishell /= jshell ) then
     do jbf=1,nj
       do ibf=1,ni
         ipair = ipair + 1

         ijmax=MAX(basis%shell(ishell)%istart + ibf - 1,basis%shell(jshell)%istart + jbf - 1)
         ijmin=MIN(basis%shell(ishell)%istart + ibf - 1,basis%shell(jshell)%istart + jbf - 1)
         ijbf = (ijmin-1) * nbf_eri - ( (ijmin-2) * (ijmin-1) ) / 2  + ijmax - ijmin + 1

         index_pair_1d(ijbf)  = ipair
         index_basis(1,ipair) = basis%shell(ishell)%istart + ibf - 1
         index_basis(2,ipair) = basis%shell(jshell)%istart + jbf - 1
       enddo
     enddo
   else
     do jbf=1,nj
       do ibf=1,jbf
         ipair = ipair + 1

         ijmax=MAX(basis%shell(ishell)%istart + ibf - 1,basis%shell(jshell)%istart + jbf - 1)
         ijmin=MIN(basis%shell(ishell)%istart + ibf - 1,basis%shell(jshell)%istart + jbf - 1)
         ijbf = (ijmin-1) * nbf_eri - ( (ijmin-2) * (ijmin-1) ) / 2  + ijmax - ijmin + 1

         index_pair_1d(ijbf)  = ipair
         index_basis(1,ipair) = basis%shell(ishell)%istart + ibf - 1
         index_basis(2,ipair) = basis%shell(jshell)%istart + jbf - 1
       enddo
     enddo
   endif

 enddo


end subroutine setup_basispair


!=========================================================================
pure function negligible_basispair(ibf,jbf)
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
 integer                      :: ip
 integer                      :: ibf,jbf
 integer                      :: n1c,n2c
 integer                      :: ni,nj
 integer                      :: ami,amj
 integer                      :: ishell,jshell
 real(dp),allocatable         :: integrals(:,:,:,:)
 real(dp)                     :: workload(nproc_world)
 integer                      :: shell_proc(basis%nshell)
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
 do jshell=1,basis%nshell
   amj = basis%shell(jshell)%am
   ip = MINLOC(workload(:),DIM=1)
   !
   ! Cost function was evaluated from a few runs
   workload(ip) = workload(ip) + cost_function_eri(amj)
   shell_proc(jshell) = ip - 1
 enddo


 negligible_shellpair(:,:) = .TRUE.

 do jshell=1,basis%nshell
   !
   ! Workload is distributed here
   if( shell_proc(jshell) /= rank_world ) cycle

   amj = basis%shell(jshell)%am
   nj  = number_basis_function_am( basis%gaussian_type , amj )
   n2c = number_basis_function_am( 'CART' , amj )
   am2 = basis%shell(jshell)%am
   ng2 = basis%shell(jshell)%ng
   allocate(alpha2(ng2),coeff2(ng2))
   alpha2(:) = basis%shell(jshell)%alpha(:)
   coeff2(:) = basis%shell(jshell)%coeff(:)
   x02(:) = basis%shell(jshell)%x0(:)

   !$OMP PARALLEL PRIVATE(ami,ni,am1,n1c,ng1,alpha1,coeff1,x01, &
   !$OMP                  int_shell,integrals)
   !$OMP DO
   do ishell=1,basis%nshell
     ami = basis%shell(ishell)%am
     if( ami < amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     n1c = number_basis_function_am( 'CART' , ami )
     am1 = basis%shell(ishell)%am
     ng1 = basis%shell(ishell)%ng

     allocate(alpha1(ng1),coeff1(ng1))
     alpha1(:) = basis%shell(ishell)%alpha(:)
     x01(:) = basis%shell(ishell)%x0(:)
     coeff1(:) = basis%shell(ishell)%coeff(:)

     allocate( int_shell( n1c*n2c*n1c*n2c ) )

     call libint_4center(am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         0.0_C_DOUBLE,int_shell)

     call transform_libint_to_molgw(basis%gaussian_type,ami,amj,ami,amj,int_shell,integrals)

     do ibf=1,ni
       do jbf=1,nj
         if( SQRT( ABS( integrals(ibf,jbf,ibf,jbf) ) ) > TOL_INT ) negligible_shellpair(ishell,jshell) = .FALSE.
       enddo
     enddo

     !
     ! Symmetrize
     negligible_shellpair(jshell,ishell) = negligible_shellpair(ishell,jshell)

     deallocate(integrals)
     deallocate(int_shell)
     deallocate(alpha1,coeff1)

   enddo
   !$OMP END DO
   !$OMP END PARALLEL
   deallocate(alpha2,coeff2)
 enddo

 call xand_world(negligible_shellpair)

 call stop_clock(timing_eri_screening)


end subroutine identify_negligible_shellpair


!=========================================================================
subroutine setup_shellpair(basis)
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer :: ishell,jshell
 integer :: ami,amj
 integer :: ishellpair,jshellpair
 integer :: ammax,amij
!=====

 !
 ! First count the number of non-negligible shell pairs
 ishellpair = 0
 jshellpair = 0
 do jshell=1,basis%nshell
   do ishell=1,jshell
     jshellpair = jshellpair + 1
     ! skip the identified negligible shell pairs
     if( negligible_shellpair(ishell,jshell) ) cycle
     ishellpair = ishellpair + 1

   enddo
 enddo
 nshellpair = ishellpair
 write(stdout,'(/,1x,a,i8,a,i8)') 'Non negligible shellpairs to be computed',nshellpair,'  over a total of',jshellpair
 write(stdout,'(1x,a,f12.4)')     'Saving (%): ', REAL(jshellpair-nshellpair,dp)/REAL(jshellpair,dp)*100.0_dp

 allocate(index_shellpair(2,nshellpair))

 ammax = MAXVAL(basis%shell(:)%am)

 !
 ! Then order the shell pairs so to enforce the LIBINT condition:
 ! If ijshellpair < klshellpair, then ami+amj < amk+aml
 !
 ishellpair = 0
 do amij=0,2*ammax

   do jshell=1,basis%nshell
     do ishell=1,jshell
       ! skip the identified negligible shell pairs
       if( negligible_shellpair(ishell,jshell) ) cycle
       ami = basis%shell(ishell)%am
       amj = basis%shell(jshell)%am

       if( ami + amj /= amij ) cycle

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

 enddo


end subroutine setup_shellpair


!=================================================================
subroutine destroy_eri_3center_lowerlevel()
 implicit none
!=====

 if(ALLOCATED(iproc_ibf_auxil)) then
   deallocate(iproc_ibf_auxil)
 endif
 if(ALLOCATED(ibf_auxil_g)) then
   deallocate(ibf_auxil_g)
 endif
 if(ALLOCATED(ibf_auxil_l)) then
   deallocate(ibf_auxil_l)
 endif
 if(ALLOCATED(eri_3center)) then
   call clean_deallocate('3-center integrals',eri_3center)
 endif

end subroutine destroy_eri_3center_lowerlevel


!=================================================================
subroutine destroy_eri_3center_lr()
 implicit none
!=====

 if(ALLOCATED(iproc_ibf_auxil_lr)) then
   deallocate(iproc_ibf_auxil_lr)
 endif
 if(ALLOCATED(ibf_auxil_g_lr)) then
   deallocate(ibf_auxil_g_lr)
 endif
 if(ALLOCATED(ibf_auxil_l_lr)) then
   deallocate(ibf_auxil_l_lr)
 endif
 if(ALLOCATED(eri_3center_lr)) then
   call clean_deallocate('LR 3-center integrals',eri_3center_lr)
 endif

end subroutine destroy_eri_3center_lr


!=========================================================================
subroutine dump_out_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!=====
 character(len=50)  :: filename
 integer(kind=int8) :: nline,iline,icurrent
 integer            :: erifile
!=====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 write(stdout,*) 'Dump out the ERI into file'
 write(stdout,*) 'Size of file (Gbytes)',REAL(nint_4center,dp) * dp / 1024.0_dp**3

 if( is_iomaster ) then
   open(newunit=erifile,file=TRIM(filename),form='unformatted')
   write(erifile) nint_4center
   write(erifile) rcut

   nline = nint_4center / line_length + 1
   icurrent=0
   do iline=1,nline
     write(erifile) eri_4center(icurrent+1:MIN(nint_4center,icurrent+line_length+1))
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
 character(len=50)  :: filename
 integer(kind=int8) :: nline,iline,icurrent
 integer(kind=int8) :: integer_read
 real(dp)           :: real_read
 integer            :: erifile
!=====

 if(rcut < 1.0e-6_dp) then
   filename = 'molgw_eri.data'
 else
   filename = 'molgw_eri_lr.data'
 endif

 inquire(file=TRIM(filename),exist=read_eri)

 if(read_eri) then

   write(stdout,*) 'Try to read ERI file'
   open(newunit=erifile,file=TRIM(filename),form='unformatted',status='old')
   read(erifile) integer_read
   if(integer_read /= nint_4center) read_eri=.FALSE.
   read(erifile) real_read
   if(ABS(real_read-rcut) > 1.0e-6_dp) read_eri=.FALSE.

   if(read_eri) then

     nline = nint_4center / line_length + 1
     icurrent=0
     do iline=1,nline
       read(erifile) eri_4center(icurrent+1:MIN(nint_4center,icurrent+line_length+1))
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
 integer :: iproc
 integer :: ilocal,iglobal
 integer :: nbf_local_iproc(0:nproc_auxil-1)
!=====

 ! Use SCALAPACK routines to distribute the auxiliary basis
 ! Assume a processor grid: 1 x nproc_auxil

 do iproc=0,npcol_eri3_ao-1
   nbf_local_iproc(iproc) = NUMROC(nbf_auxil_basis,NB_eri3_ao,iproc,first_col,npcol_eri3_ao)
 enddo

 nauxil_3center = nbf_local_iproc(ipcol_eri3_ao)

 allocate(ibf_auxil_g(nauxil_3center))
 do ilocal=1,nauxil_3center
   ibf_auxil_g(ilocal) = INDXL2G(ilocal,NB_eri3_ao,ipcol_eri3_ao,first_col,npcol_eri3_ao)
 enddo
 allocate(ibf_auxil_l(nbf_auxil_basis))
 allocate(iproc_ibf_auxil(nbf_auxil_basis))
 do iglobal=1,nbf_auxil_basis
   ibf_auxil_l(iglobal)     = INDXG2L(iglobal,NB_eri3_ao,0,first_col,npcol_eri3_ao)
   iproc_ibf_auxil(iglobal) = INDXG2P(iglobal,NB_eri3_ao,0,first_col,npcol_eri3_ao)
 enddo


 write(stdout,'(/,a)') ' Distribute auxiliary basis functions among processors'
 write(stdout,'(1x,a,i4,a,i6,a)') 'Max auxiliary basis functions ', &
                                  MAXVAL(nbf_local_iproc(:)),' for processor ',MAXLOC(nbf_local_iproc,DIM=1)
 write(stdout,'(1x,a,i4,a,i6,a)') 'Min auxiliary basis functions ', &
                                  MINVAL(nbf_local_iproc(:)),' for processor ',MINLOC(nbf_local_iproc,DIM=1)


end subroutine distribute_auxil_basis


!=========================================================================
subroutine distribute_auxil_basis_lr(nbf_auxil_basis)
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
!=====
 integer :: iproc
 integer :: ilocal,iglobal
 integer :: nbf_local_iproc_lr(0:nproc_auxil-1)
!=====
 integer :: ibf,ibf_local
!=====


 ! Use SCALAPACK routines to distribute the auxiliary basis
 ! Assume a processor grid: 1 x nproc_auxil

 do iproc=0,npcol_eri3_ao-1
   nbf_local_iproc_lr(iproc) = NUMROC(nbf_auxil_basis,NB_eri3_ao,iproc,first_col,npcol_eri3_ao)
 enddo

 nauxil_3center_lr = nbf_local_iproc_lr(ipcol_eri3_ao)

 allocate(ibf_auxil_g_lr(nauxil_3center_lr))
 do ilocal=1,nauxil_3center_lr
   ibf_auxil_g_lr(ilocal) = INDXL2G(ilocal,NB_eri3_ao,ipcol_eri3_ao,first_col,npcol_eri3_ao)
 enddo
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 allocate(iproc_ibf_auxil_lr(nbf_auxil_basis))
 do iglobal=1,nbf_auxil_basis
   ibf_auxil_l_lr(iglobal)     = INDXG2L(iglobal,NB_eri3_ao,0,first_col,npcol_eri3_ao)
   iproc_ibf_auxil_lr(iglobal) = INDXG2P(iglobal,NB_eri3_ao,0,first_col,npcol_eri3_ao)
 enddo

 write(stdout,'(/,a)') ' Distribute LR auxiliary basis functions among processors'
 write(stdout,'(1x,a,i4,a,i6,a)') 'Max auxiliary basis functions ', &
                                  MAXVAL(nbf_local_iproc_lr(:)),' for processor ',MAXLOC(nbf_local_iproc_lr,DIM=1)
 write(stdout,'(1x,a,i4,a,i6,a)') 'Min auxiliary basis functions ', &
                                  MINVAL(nbf_local_iproc_lr(:)),' for processor ',MINLOC(nbf_local_iproc_lr,DIM=1)


end subroutine distribute_auxil_basis_lr


!=========================================================================
subroutine reshuffle_distribution_3center()
 implicit none

!=====
 integer :: info
 integer :: mlocal,nlocal
 integer :: desc3final(NDEL)
 real(dp),allocatable :: eri_3center_tmp(:,:)
!=====

#if defined(HAVE_SCALAPACK)
 write(stdout,'(/,a,i8,a,i4)') ' Final 3-center integrals distributed using a SCALAPACK grid: ',nprow_eri3_ao,' x ',npcol_eri3_ao

 if( nprow_eri3_ao == nprow_3center .AND. npcol_eri3_ao == npcol_3center .AND. MB_eri3_ao == MB_3center ) then
   write(stdout,*) 'Reshuffling not needed'
   return
 endif

 if( cntxt_eri3_ao > 0 ) then
   mlocal = NUMROC(nauxil_2center,MB_eri3_ao,iprow_eri3_ao,first_row,nprow_eri3_ao)
   nlocal = NUMROC(npair         ,NB_eri3_ao,ipcol_eri3_ao,first_col,npcol_eri3_ao)
 else
   mlocal = -1
   nlocal = -1
 endif
 call xmax_ortho(mlocal)
 call xmax_ortho(nlocal)

 if( cntxt_3center > 0 ) then
   call move_alloc(eri_3center,eri_3center_tmp)

   call DESCINIT(desc3final,nauxil_2center,npair,MB_eri3_ao,NB_eri3_ao,first_row,first_col,cntxt_eri3_ao,MAX(1,mlocal),info)

   call clean_allocate('TMP 3-center integrals',eri_3center,mlocal,nlocal)

   call PDGEMR2D(nauxil_2center,npair,eri_3center_tmp,1,1,desc_eri3,eri_3center,1,1,desc3final,cntxt_3center)
   call clean_deallocate('TMP 3-center integrals',eri_3center_tmp)
 else
   call clean_deallocate('3-center integrals',eri_3center)
   call clean_allocate('3-center integrals',eri_3center,mlocal,nlocal)
 endif

 !
 ! Propagate to the ortho MPI direction
 if( cntxt_eri3_ao <= 0 ) then
   eri_3center(:,:) = 0.0_dp
 endif
 call xsum_ortho(eri_3center)
#endif


end subroutine reshuffle_distribution_3center


!=========================================================================
end module m_eri
