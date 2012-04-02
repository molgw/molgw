module m_eri
 use m_definitions

 real(dp),allocatable :: eri_buffer(:)
 integer              :: nbf_eri                ! local copy of nbf
 integer              :: nsize                  ! size of the eri_buffer array
 integer              :: nsize1                 ! number of independent pairs (i,j) with i<=j

 real(dp),pointer :: eri_eigen_buffer(:,:,:)

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
#if SYMMETRIZED
 nsize1  = index_prod(nbf_eri,nbf_eri) !  nbf_eri*(nbf_eri+1)/2
 nsize   = index_eri(nbf_eri,nbf_eri,nbf_eri,nbf_eri)
#else
 nsize1  = nbf_eri**2
 nsize   = nsize1**2
#endif

 allocate(eri_buffer(nsize),stat=info)
 write(*,'(/,a,f8.3,a)') ' Allocating the ERI array: ',REAL(nsize,dp)*dp/1024**3,' [Gb]'
 if(info==0) then
   write(*,*) 'success'
 else
   write(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif

 eri_buffer(:) = 0.0_dp

end subroutine allocate_eri

!=========================================================================
subroutine allocate_eri_eigen(nspin)
 implicit none
!===== 
 integer,intent(in) :: nspin
!===== 
 integer            :: info
!===== 

 allocate(eri_eigen_buffer(nsize,nspin,nspin),stat=info)
 write(*,'(/,a,f8.3,a)') ' Allocating the ERI array: ',nspin**2*REAL(nsize,dp)*dp/1024**2,' [Mb]'
 if(info==0) then
   write(*,*) 'success'
 else
   write(*,*) 'failure'
   stop'Not enough memory. Ask Santa Claus for a bigger computer'
 endif

 eri_eigen_buffer(:,:,:) = 0.0_dp

end subroutine allocate_eri_eigen

!=========================================================================
subroutine deallocate_eri()
 implicit none
!=====

 deallocate(eri_buffer)

end subroutine deallocate_eri

!=========================================================================
subroutine deallocate_eri_eigen()
 implicit none
!=====

 deallocate(eri_eigen_buffer)

end subroutine deallocate_eri_eigen

!=========================================================================
function index_prod(ibf,jbf)
 implicit none
 integer,intent(in) :: ibf,jbf
 integer            :: index_prod
!=====
 integer            :: imin,jmax
!=====

 imin=MIN(ibf,jbf)
 jmax=MAX(ibf,jbf)
 index_prod = (imin-1)*nbf_eri - (imin-1)*(imin-2)/2 + jmax-imin+1

end function index_prod

!=========================================================================
function index_eri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 integer            :: index_eri
!=====
! integer            :: imin,jmax,kmin,lmax
 integer            :: ijmin,klmax
 integer            :: index_ij,index_kl
!===== 

 index_ij = index_prod(ibf,jbf)
 index_kl = index_prod(kbf,lbf)

 ijmin=MIN(index_ij,index_kl)
 klmax=MAX(index_ij,index_kl)
 index_eri = (ijmin-1)*nsize1 - (ijmin-1)*(ijmin-2)/2 + klmax-ijmin+1

end function index_eri

!=========================================================================
function eri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri
!=====

#if SYMMETRIZED
 eri = eri_buffer(index_eri(ibf,jbf,kbf,lbf))
#else
 eri = eri_buffer(ibf+(jbf-1)*nbf_eri+(kbf-1)*nbf_eri**2+(lbf-1)*nbf_eri**3)
#endif

end function eri

!=========================================================================
function eri_eigen(ibf,jbf,kbf,lbf,ijspin,klspin)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 integer,intent(in) :: ijspin,klspin
 real(dp)           :: eri_eigen
!=====

 eri_eigen = eri_eigen_buffer(index_eri(ibf,jbf,kbf,lbf),ijspin,klspin)

end function eri_eigen

#if 0
!=========================================================================
subroutine calculate_eri(basis,eri)
 use ISO_C_BINDING
 use m_definitions
 use m_timing
 use m_basis_set
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(out)         :: eri(basis%nbf,basis%nbf,basis%nbf,basis%nbf)
!=====
 integer,parameter            :: NSHELLMAX=100
 integer,parameter            :: NPRIMITIVE_MAX=20
 logical,parameter            :: DEBUG=.TRUE.
 integer                      :: info
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: ig,jg,kg,lg
 integer                      :: ni,nj,nk,nl
 integer                      :: i_current,i,j,k,l
 integer                      :: nshell,nintegrals
 integer                      :: ishell,jshell,kshell,lshell
 type(basis_function),pointer :: bf_current_i,bf_current_j,bf_current_k,bf_current_l
 type(gaussian),pointer       :: g_current
 integer                      :: nint_gaussian,igaussian
 integer                      :: index_global(basis%nbf,NPRIMITIVE_MAX)  ! no basis set will contain a function made of more than 10 gaussian
 real(dp),allocatable         :: int_gaussian(:,:,:,:)
! real*4,allocatable         :: int_gaussian(:,:,:,:)
 logical                      :: shell_already_exists
 integer                      :: ami,amj,amk,aml
 integer                      :: ii,jj,kk,ll
!=====
 type shell_type
   integer  :: am
   real(dp) :: alpha
   integer  :: index_global(NPRIMITIVE_MAX)=0
 end type shell_type
 type(shell_type)             :: shell(NSHELLMAX)
!=====
! variables used to call C++ 
 integer(C_INT),external      :: calculate_integral
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: alpha1,alpha2,alpha3,alpha4
 real(C_DOUBLE),allocatable   :: integrals(:)
!=====

 write(*,'(/,a)') ' Calculate all the Electron Repulsion Integrals (ERI) at once'

 call start_clock(timing_tmp1)
 !
 ! libint works with shells of same angular momentum and same alpha
 ! establish now a list of all gaussian shells

 nint_gaussian=0
 igaussian=0
 ishell=0
 do ibf=1,basis%nbf
   bf_current_i => basis%bf(ibf)

   do ig=1,bf_current_i%ngaussian
     g_current => bf_current_i%g(ig)

     igaussian=igaussian+1
     !
     ! check if the shell has already been created
     shell_already_exists=.FALSE.
     do jshell=1,ishell
       if( g_current%am == shell(jshell)%am .AND. ABS( g_current%alpha - shell(jshell)%alpha ) < 1.d-8 ) then
         !TODO
         shell(jshell)%index_global(libint_ordering(g_current%nx,g_current%ny,g_current%nz)) = igaussian
         shell_already_exists=.TRUE.
         index_global(ibf,ig) = igaussian
         exit
       endif
     enddo

     if(.NOT.shell_already_exists) then
       ishell=ishell+1
       if(ishell > NSHELLMAX) stop'NSHELLMAX internal parameter is too small'
       nint_gaussian = nint_gaussian + number_basis_function_am( g_current%am )

       shell(ishell)%am                                                                    = g_current%am
       shell(ishell)%alpha                                                                 = g_current%alpha
       shell(ishell)%index_global(libint_ordering(g_current%nx,g_current%ny,g_current%nz)) = igaussian
       index_global(ibf,ig)                                                                = igaussian

     endif

   enddo
 enddo
 nshell=ishell
 if(igaussian/=nint_gaussian) then
   write(*,*) igaussian,nint_gaussian
   stop'one shell is not complete'
 endif
 if(DEBUG) write(*,*) 
 if(DEBUG) write(*,*) 'found nshell=',nshell,nint_gaussian
 if(DEBUG) write(*,*) 'alpha=',shell(1:nshell)%alpha
 if(DEBUG) write(*,*) 'am   =',shell(1:nshell)%am
 if(DEBUG) write(*,*)

 allocate(int_gaussian(nint_gaussian,nint_gaussian,nint_gaussian,nint_gaussian))
 
 call stop_clock(timing_tmp1)
 call start_clock(timing_tmp2)
 ig=0; jg=0; kg=0; lg=0
 
 do lshell=1,nshell
   write(*,*) lshell,nshell
   do kshell=1,nshell
     do jshell=1,nshell
       do ishell=1,nshell

         ni = number_basis_function_am(  shell(ishell)%am )
         nj = number_basis_function_am(  shell(jshell)%am )
         nk = number_basis_function_am(  shell(kshell)%am )
         nl = number_basis_function_am(  shell(lshell)%am )
         
         nintegrals = ni * nj * nk *nl
         if( allocated(integrals) ) deallocate( integrals )
         allocate( integrals( nintegrals ) )

         !
         ! If all the angular momentum are S type, do the calculation immediately
         ! else call the libint library
         !
         ami = shell(ishell)%am
         amj = shell(jshell)%am
         amk = shell(kshell)%am
         aml = shell(lshell)%am
         if(ami+amj+amk+aml==0) then

           alpha1 = shell(ishell)%alpha
           alpha2 = shell(jshell)%alpha
           alpha3 = shell(kshell)%alpha
           alpha4 = shell(lshell)%alpha
           integrals(1) = 2.0_dp*pi**(2.5_dp) &
&                        / ( (alpha1+alpha2)*(alpha3+alpha4)*SQRT(alpha1+alpha2+alpha3+alpha4) )

           int_gaussian( shell(ishell)%index_global(1),&
                         shell(jshell)%index_global(1),&
                         shell(kshell)%index_global(1),&
                         shell(lshell)%index_global(1) ) = integrals(1)

         else
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
             else
               cycle  ! SAVE
!               am3 = shell(lshell)%am
!               am4 = shell(kshell)%am
!               alpha3 = shell(lshell)%alpha
!               alpha4 = shell(kshell)%alpha
             endif
             if(ami>=amj) then
               am1 = shell(ishell)%am
               am2 = shell(jshell)%am
               alpha1 = shell(ishell)%alpha
               alpha2 = shell(jshell)%alpha
             else
               cycle  ! SAVE
!               am1 = shell(jshell)%am
!               am2 = shell(ishell)%am
!               alpha1 = shell(jshell)%alpha
!               alpha2 = shell(ishell)%alpha
             endif
           else
               cycle  ! SAVE
!             if(amk>=aml) then
!               am1 = shell(kshell)%am
!               am2 = shell(lshell)%am
!               alpha1 = shell(kshell)%alpha
!               alpha2 = shell(lshell)%alpha
!             else
!               am1 = shell(lshell)%am
!               am2 = shell(kshell)%am
!               alpha1 = shell(lshell)%alpha
!               alpha2 = shell(kshell)%alpha
!             endif
!             if(ami>=amj) then
!               am3 = shell(ishell)%am
!               am4 = shell(jshell)%am
!               alpha3 = shell(ishell)%alpha
!               alpha4 = shell(jshell)%alpha
!             else
!               am3 = shell(jshell)%am
!               am4 = shell(ishell)%am
!               alpha3 = shell(jshell)%alpha
!               alpha4 = shell(ishell)%alpha
!             endif
           endif
!           if(DEBUG) write(*,*) 'just before the call to libint'
!           if(DEBUG) write(*,*) nintegrals
!           if(DEBUG) write(*,*) shell(ishell)%am,shell(jshell)%am,shell(kshell)%am,shell(lshell)%am
!           if(DEBUG) write(*,*) am1,am2,am3,am4
!           if(DEBUG) write(*,*) alpha1,alpha2,alpha3,alpha4
!           write(*,*) 'call to C library'

           info=calculate_integral(am1,am2,am3,am4,alpha1,alpha2,alpha3,alpha4,&
                                 integrals(1))
           if(info/=0) then
              write(*,*) 'failed'
!           else
!             write(*,*) 'success'
           endif
!          write(*,*) ishell,jshell,kshell,lshell
!         if(DEBUG) write(*,*) 'shell i j k l,    result'
!         if(DEBUG) write(*,*) ishell,jshell,kshell,lshell
!         if(DEBUG) write(*,*) integrals(1)
!         if(DEBUG) write(*,*) integrals(:)
!         if(DEBUG) write(*,*) '------------------------'

           i_current=0
           do i=1,ni
           do j=1,nj
           do k=1,nk
           do l=1,nl
             i_current=i_current+1
             ii = shell(ishell)%index_global(i)
             jj = shell(jshell)%index_global(j)
             kk = shell(kshell)%index_global(k)
             ll = shell(lshell)%index_global(l)
             int_gaussian( ii , jj , kk , ll ) = integrals(i_current)
             int_gaussian( ii , jj , ll , kk ) = integrals(i_current)
             int_gaussian( jj , ii , kk , ll ) = integrals(i_current)
             int_gaussian( jj , ii , ll , kk ) = integrals(i_current)
             int_gaussian( kk , ll , ii , jj ) = integrals(i_current)
             int_gaussian( kk , ll , jj , ii ) = integrals(i_current)
             int_gaussian( ll , kk , ii , jj ) = integrals(i_current)
             int_gaussian( ll , kk , jj , ii ) = integrals(i_current)

           enddo
           enddo
           enddo
           enddo

         endif



       enddo
     enddo
   enddo
 enddo
 if( allocated(integrals) ) deallocate( integrals )

 call stop_clock(timing_tmp2)

! if(DEBUG) write(*,*) '=========================================='
! if(DEBUG) write(*,*) int_gaussian(:,:,:,:)
! if(DEBUG) write(*,*) '=========================================='

 !
 ! calculate (ij||kl) over contractions
 
 call start_clock(timing_tmp3)
 eri(:,:,:,:) = 0.0_dp
 do lbf=1,basis%nbf
   bf_current_l => basis%bf(lbf)
   do kbf=1,basis%nbf
     bf_current_k => basis%bf(kbf)
     do jbf=1,basis%nbf
       bf_current_j => basis%bf(jbf)
       do ibf=1,basis%nbf
         bf_current_i => basis%bf(ibf)

         do lg=1,bf_current_l%ngaussian
           do kg=1,bf_current_k%ngaussian
             do jg=1,bf_current_j%ngaussian
               do ig=1,bf_current_i%ngaussian
                 eri(ibf,jbf,kbf,lbf) = eri(ibf,jbf,kbf,lbf) &
                            + bf_current_i%coeff(ig) * bf_current_j%coeff(jg) * bf_current_k%coeff(kg) * bf_current_l%coeff(lg) &
                              * int_gaussian( index_global(ibf,ig) , index_global(jbf,jg) , index_global(kbf,kg) , index_global(lbf,lg) ) &
                               * bf_current_i%g(ig)%norm_factor * bf_current_j%g(jg)%norm_factor &
                               * bf_current_k%g(kg)%norm_factor * bf_current_l%g(lg)%norm_factor
               enddo
             enddo
           enddo
         enddo


       enddo
     enddo
   enddo
 enddo
 call stop_clock(timing_tmp3)


 deallocate(int_gaussian)
          
 write(*,*) 'Done.'
 write(*,*)

end subroutine calculate_eri
#endif

!=========================================================================
subroutine calculate_eri2(basis)
 use ISO_C_BINDING
 use m_definitions
 use m_tools,only: boys_function
 use m_timing
 use m_basis_set
#ifdef OPENMP
 use omp_lib
#endif
 implicit none
 type(basis_set),intent(in)   :: basis
!=====
 integer,parameter            :: NSHELLMAX=100
 integer,parameter            :: NMEMBER=28
 integer,parameter            :: NPRIMITIVE_MAX=20
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
 integer                      :: index_global(basis%nbf,NPRIMITIVE_MAX)  
 real(dp),allocatable         :: int_gaussian(:)
 logical                      :: shell_already_exists
 integer                      :: iint,nint_tot
 integer                      :: imember,jmember,kmember,lmember
 integer                      :: iindex_in_the_shell,jindex_in_the_shell,kindex_in_the_shell,lindex_in_the_shell
 integer                      :: index_integral
 integer                      :: index_tmp
 integer                      :: ami,amj,amk,aml
 real(dp),allocatable         :: int_tmp(:,:,:,:)
 real(dp)                     :: zeta_12,zeta_34,rho,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
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
!=====


 write(*,'(/,a)') ' Calculate and store all the Electron Repulsion Integrals (ERI)'

 call start_clock(timing_tmp1)
 
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
 
 call stop_clock(timing_tmp1)
 call start_clock(timing_tmp2)
 write(*,*) 'number of shells',nshell

 !
 ! (ij||kl)

!!!!    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iindex_in_the_shell,jindex_in_the_shell,kindex_in_the_shell,lindex_in_the_shell,&
!!!!    !$OMP&                                 ig,jg,kg,lg,ibf,jbf,kbf,lbf,index_integral,index_tmp,&
!!!!    !$OMP&                                 ami,amj,amk,aml,am1,am2,am3,am4,alpha1,alpha2,alpha3,alpha4,&
!!!!    !$OMP&                                 n1,n2,n3,n4,ni,nj,nk,nl,nint_shell,int_shell,int_tmp,ii,info,i,j,k,l,ishell,jshell,kshell,lshell)
!!!!    
!!!!    
!!!!    !$OMP DO SCHEDULE(STATIC)
 do jshell=1,nshell
   do ishell=1,nshell

     do lshell=1,nshell
       do kshell=1,nshell

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

         n1 = number_basis_function_am( am1 )
         n2 = number_basis_function_am( am2 )
         n3 = number_basis_function_am( am3 )
         n4 = number_basis_function_am( am4 )
         
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
           rho = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
           tt = rho * SUM( (p(:)-q(:))**2 )
           call boys_function(f0t(0),0,tt)
           int_shell(1) = 2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                 / zeta_12 * exp( -alpha1*alpha2/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                 / zeta_34 * exp( -alpha3*alpha4/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) 

         else

           info=calculate_integral(am1,am2,am3,am4,alpha1,alpha2,alpha3,alpha4,&
                                 x01(1),x01(2),x01(3),&
                                 x02(1),x02(2),x02(3),&
                                 x03(1),x03(2),x03(3),&
                                 x04(1),x04(2),x04(3),&
                                 int_shell(1))
           if(info/=0) then
             write(*,*) am1,am2,am3,am4
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
!!!!!           !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iindex_in_the_shell,jindex_in_the_shell,kindex_in_the_shell,lindex_in_the_shell,&
!!!!!           !$OMP&     ig,jg,kg,lg,ibf,jbf,kbf,lbf,index_integral,index_tmp )
!!!!!           
!!!!!           !$OMP DO SCHEDULE(STATIC)
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

#if SYMMETRIZED
                 if(ibf>jbf) cycle
                 if(kbf>lbf) cycle
                 if(index_prod(ibf,jbf)>index_prod(kbf,lbf)) cycle

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
!!!!          !$OMP END DO
!!!!          !$OMP END PARALLEL

         call stop_clock(timing_tmp4)

       enddo
     enddo
   enddo
 enddo
!!!!       !$OMP END DO
 if( allocated(int_shell) ) deallocate( int_shell )
!!!!       !$OMP END PARALLEL

 call stop_clock(timing_tmp2)
 write(*,*) 'Done!'
 write(*,*)

end subroutine calculate_eri2


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

 case default
   stop'libint_ordering not coded for this orbital momentum'
 end select

end function libint_ordering



!=========================================================================
subroutine test_eri(basis,eri)
 use m_definitions
 use m_basis_set
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: eri(basis%nbf,basis%nbf,basis%nbf,basis%nbf)
!=====
 integer                      :: ibf,jbf,kbf,lbf
!=====

 do ibf=1,basis%nbf
   do jbf=1,basis%nbf
     do kbf=1,basis%nbf
       do lbf=1,basis%nbf
         if( ABS(eri(ibf,jbf,kbf,lbf) - eri(kbf,lbf,ibf,jbf)) > 1.d-6 ) then
           write(*,*) ibf,jbf,kbf,lbf,eri(ibf,jbf,kbf,lbf)
           write(*,*) kbf,lbf,ibf,jbf,eri(kbf,lbf,ibf,jbf)
           write(*,*) ibf,basis%bf(ibf)%amc
           write(*,*) jbf,basis%bf(jbf)%amc
           write(*,*) kbf,basis%bf(kbf)%amc
           write(*,*) lbf,basis%bf(lbf)%amc
           stop'ERI array not symmetric'
         endif
       enddo
     enddo
   enddo
 enddo


end subroutine test_eri

!=================================================================
subroutine transform_eri_basis_robust(nbf,nspin,c_matrix,eri_eigenstate)
 use m_definitions
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

 write(*,*) 'obtain the ERI in the eigenvector basis'

 eri_eigenstate(:,:,:,:,:,:)=0.0_dp
 do ijspin=1,nspin
   write(*,*) '== ijspin ',ijspin
   eri_tmp(:,:,:,:)=0.0_dp
   do lbf=1,nbf
     do kbf=1,nbf
       do jstate=1,nbf
         do istate=1,nbf

           do jbf=1,nbf
             do ibf=1,nbf
               eri_tmp(istate,jstate,kbf,lbf) = eri_tmp(istate,jstate,kbf,lbf) + eri(ibf,jbf,kbf,lbf) * c_matrix(ibf,istate,ijspin) * c_matrix(jbf,jstate,ijspin)
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
 use m_definitions
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

 write(*,*) 'obtain the ERI in the eigenvector basis'
 write(*,*) 'subroutine is order N^5'
 call start_clock(timing_basis_transform)

#ifdef OPENMP
 wtime=OMP_get_wtime()
 write(*,*) 'The basis transform is using OPENMP'
#endif

 do ijspin=1,nspin
   write(*,*) '== ijspin ',ijspin
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
  write(*,*) 'time (s)', OMP_get_wtime()-wtime
#endif

 call stop_clock(timing_basis_transform)
 write(*,*) 'ERI in the eigenvector basis obtained'
 write(*,*)

end subroutine transform_eri_basis_fast


end module m_eri
