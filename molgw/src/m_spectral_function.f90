!=========================================================================
#include "macros.h"
!=========================================================================
module m_spectral_function
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning

 !
 ! General form of any spectral function
 ! z complex number
 ! i, j running on the basis set
 ! sf_ij(z) = \sum_n L_n(i) R_n(j) / ( z - w_n )
 !
 type spectral_function 
   integer :: npole
   integer :: nprodbasis
   real(dp),allocatable :: pole(:)
   real(dp),allocatable :: residu_left(:,:)       ! first index runs on n, second index on i
   real(dp),allocatable :: residu_right(:,:)      ! first index runs on n, second index on j
 end type spectral_function

 !
 ! frozen core approximation parameters
 integer :: ncore_G=0
 integer :: ncore_W=0

 !
 ! frozen virtual approximation parameters
 integer :: nvirtual_G=HUGE(1)
 integer :: nvirtual_W=HUGE(1)

 !
 ! the boring small complex number eta: (0.0_dp,0.0001_dp) is typically over converged
 complex(dp),parameter :: ieta=(0.0_dp,0.001_dp) ! (0.0_dp,0.0001_dp)

#ifdef CRPA
 integer,parameter :: band1=1
 integer,parameter :: band2=2
#endif
contains

!=========================================================================
subroutine init_spectral_function(nbf,prod_nbf,nspin,occupation,sf)
 implicit none
 integer,intent(in)                    :: nbf,prod_nbf,nspin
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 type(spectral_function),intent(inout) :: sf
!=====
 integer                               :: ispin,ibf,jbf
 logical                               :: file_exists
!====

 !
 ! Deal with frozen core initialization
 inquire(file='manual_frozencore',exist=file_exists)
 if(file_exists) then
   !
   ! ncore_G and ncore_W contain the highest core state to be discarded
   open(13,file='manual_frozencore')
   read(13,*) ncore_G
   read(13,*) ncore_W
   close(13)
   ncore_G = MAX(ncore_G,0)
   ncore_W = MAX(ncore_W,0)
   WRITE_MASTER(msg,'(a,i4,2x,i4)') 'frozen core approximation switched on up to state (G,W) = ',ncore_G,ncore_W
   call issue_warning(msg)
 endif
 !
 ! Deal with frozen virtual initialization
 inquire(file='manual_frozenvirtual',exist=file_exists)
 if(file_exists) then
   !
   ! nvirtual_G and nvirtual_W contain the first virtual state to be discarded
   open(13,file='manual_frozenvirtual')
   read(13,*) nvirtual_G
   read(13,*) nvirtual_W
   close(13)
   WRITE_MASTER(msg,'(a,i4,2x,i4)') 'frozen virtual approximation switched on starting with state (G,W) = ',nvirtual_G,nvirtual_W
   call issue_warning(msg)
 endif




 sf%npole=0
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=1,nbf
       if( skip_transition(nspin,ibf,jbf,occupation(ibf,ispin),occupation(jbf,ispin)) ) cycle
       sf%npole = sf%npole+1
     enddo
   enddo
 enddo

 WRITE_MASTER(*,*) 
 WRITE_MASTER(*,*) 'spectral function initialized with',sf%npole,'poles'

 allocate(sf%pole(sf%npole))

#ifdef AUXIL_BASIS
 sf%nprodbasis = prod_nbf
#else
 sf%nprodbasis = prod_nbf * nspin
#endif

#ifndef CASIDA
 allocate(sf%residu_left (sf%npole,sf%nprodbasis))
 allocate(sf%residu_right(sf%npole,sf%nprodbasis))
 WRITE_MASTER(*,*) '           second index size',sf%nprodbasis

 WRITE_MASTER(*,'(a,f14.0)') ' Memory [Mb] ',REAL(SIZE(sf%residu_left(:,:)),dp)*2.0_dp/1024.0_dp**2*dp
 WRITE_MASTER(*,*)
#endif


end subroutine init_spectral_function

!=========================================================================
function skip_transition(nspin,ib1,ib2,occ1,occ2)
 implicit none
 logical             :: skip_transition
 integer,intent(in)  :: nspin,ib1,ib2
 real(dp),intent(in) :: occ1,occ2
!=====
 real(dp)            :: spin_fact
!=====
 spin_fact = REAL(-nspin+3,dp)

 skip_transition = .FALSE.

 !
 ! skip the core states if asked for a frozen-core calculation
 if( ib1 <= ncore_W .OR. ib2 <= ncore_W ) skip_transition=.TRUE.
 !
 ! skip the virtual states if asked for a frozen-virtual calculation
 if( ib1 >= nvirtual_W .OR. ib2 >= nvirtual_W ) skip_transition=.TRUE.

#ifndef CASIDA
 if( occ1 < completely_empty             .AND. occ2 < completely_empty )             skip_transition=.TRUE.
 if( occ1 > spin_fact - completely_empty .AND. occ2 > spin_fact - completely_empty ) skip_transition=.TRUE.
#else
 ! Casida case only positive transition are retained
 ! state 1 should be empty AND state 2 should be occupied
 if( occ1 > spin_fact - completely_empty .OR.  occ2 < completely_empty )             skip_transition=.TRUE.
#endif

#ifdef CRPA
 if( ib1==band1 .AND. ib2==band1 )  skip_transition=.TRUE.
 if( ib1==band1 .AND. ib2==band2 )  skip_transition=.TRUE.
 if( ib1==band2 .AND. ib2==band1 )  skip_transition=.TRUE.
 if( ib1==band2 .AND. ib2==band2 )  skip_transition=.TRUE.
#endif

end function skip_transition

!=========================================================================
subroutine destroy_spectral_function(sf)
 implicit none
 type(spectral_function),intent(inout) :: sf
!=====

 deallocate(sf%pole)
#ifndef CASIDA
 deallocate(sf%residu_left)
 deallocate(sf%residu_right)
#endif

 WRITE_MASTER(*,*) 
 WRITE_MASTER(*,*) 'spectral function destroyed'

end subroutine destroy_spectral_function


!=========================================================================
subroutine write_spectral_function(sf)
 implicit none
 type(spectral_function),intent(in) :: sf
!=====
 integer,parameter :: spectralfile=50
 integer :: iprodbasis
!=====

 WRITE_MASTER(*,'(/,a,/)') ' dumping spectral function on file' 

 open(spectralfile,file='spectral_file',form='unformatted')

 WRITE_MASTER(spectralfile) sf%npole
 WRITE_MASTER(spectralfile) sf%nprodbasis
 WRITE_MASTER(spectralfile) sf%pole(:)
 do iprodbasis=1,sf%nprodbasis
   WRITE_MASTER(spectralfile) sf%residu_left(:,iprodbasis)
 enddo
 do iprodbasis=1,sf%nprodbasis
   WRITE_MASTER(spectralfile) sf%residu_right(:,iprodbasis)
 enddo

 close(spectralfile)

end subroutine write_spectral_function


!=========================================================================
subroutine read_spectral_function(sf,reading_status)
 implicit none
 type(spectral_function),intent(inout) :: sf
 integer,intent(out)                   :: reading_status
!=====
 integer,parameter :: spectralfile=50
 integer           :: iprodbasis
 logical           :: file_exists
 integer           :: npole_read,nprodbasis_read
!=====

 WRITE_MASTER(*,'(/,a)') ' Try to read spectral function from file spectral_file' 

 inquire(file='spectral_file',exist=file_exists)
 if( .NOT. file_exists ) then
   WRITE_MASTER(*,'(a,/)') ' File does not exist'
   reading_status=1

 else

   open(spectralfile,file='spectral_file',status='old',form='unformatted')

   read(spectralfile) npole_read
   read(spectralfile) nprodbasis_read

   if( npole_read /= sf%npole .OR. nprodbasis_read /= sf%nprodbasis ) then
     WRITE_MASTER(*,'(a,/)')     ' File does not have the correct size'
     WRITE_MASTER(*,'(i5,a,i5)') npole_read,' vs ',sf%npole
     WRITE_MASTER(*,'(i5,a,i5)') nprodbasis_read,' vs ',sf%nprodbasis
     reading_status=2
   else

     read(spectralfile) sf%pole(:)
     do iprodbasis=1,sf%nprodbasis
       read(spectralfile) sf%residu_left(:,iprodbasis)
     enddo
     do iprodbasis=1,sf%nprodbasis
       read(spectralfile) sf%residu_right(:,iprodbasis)   
     enddo

     reading_status=0
     msg='reading spectral function from spectral_file'
     call issue_warning(msg)

   endif
   close(spectralfile)

 endif

end subroutine read_spectral_function



end module m_spectral_function
!=========================================================================





