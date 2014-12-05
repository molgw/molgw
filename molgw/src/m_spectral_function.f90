!=========================================================================
#include "macros.h"
!=========================================================================
module m_spectral_function
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_inputparam,only: nspin,spin_fact

 !
 ! General form of any spectral function
 ! z complex number
 ! i, j running on the basis set
 ! sf_ij(z) = \sum_n L_n(i) R_n(j) / ( z - w_n )
 !

 real(dp),parameter     :: TOL_DEG_POLE=1.0e-5_dp

 type spectral_function 
   integer              :: npole
   integer              :: nprodbasis
   real(dp),allocatable :: pole(:)
   real(dp),allocatable :: residu_left(:,:)       ! first index runs on n, second index on i
   real(dp),allocatable :: residu_right(:,:)      ! first index runs on n, second index on j
 end type spectral_function

 !
 ! frozen core approximation parameters
 integer,protected :: ncore_G=0
 integer,protected :: ncore_W=0

 !
 ! frozen virtual approximation parameters
 integer,protected :: nvirtual_G=HUGE(1)
 integer,protected :: nvirtual_W=HUGE(1)

 !
 ! the boring small complex number eta: (0.0_dp,0.0001_dp) is typically over converged
 complex(dp),parameter :: ieta=(0.0_dp,0.001_dp) ! (0.0_dp,0.0001_dp)

#ifdef CRPA
 integer,parameter :: band1=1
 integer,parameter :: band2=2
#endif

contains

!=========================================================================
subroutine init_spectral_function(nbf,occupation,ntransition)
 implicit none
 integer,intent(in)                    :: nbf
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 integer,intent(out)                   :: ntransition
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


 ntransition=0
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=1,nbf
       if( skip_transition(nspin,ibf,jbf,occupation(ibf,ispin),occupation(jbf,ispin)) ) cycle
       ntransition = ntransition + 1
     enddo
   enddo
 enddo


end subroutine init_spectral_function


!=========================================================================
subroutine allocate_spectral_function(npole,nprodbasis,sf)
 implicit none
 integer,intent(in)                    :: npole,nprodbasis
 type(spectral_function),intent(inout) :: sf
!=====

 sf%npole      = npole
 sf%nprodbasis = nprodbasis * nspin

 WRITE_MASTER(*,'(/,a,i8)') ' Spectral function initialized with npoles              : ',sf%npole
 WRITE_MASTER(*,'(a,i8)')   ' Spectral function initialized with prod basis functions: ',sf%nprodbasis

 allocate(sf%pole(sf%npole))
 allocate(sf%residu_left (sf%npole,sf%nprodbasis))
 allocate(sf%residu_right(sf%npole,sf%nprodbasis))

 WRITE_MASTER(*,'(a,f14.3,/)') ' Memory [Gb] ',REAL(sf%npole*sf%nprodbasis,dp)*2/1024.0_dp**3*dp


end subroutine allocate_spectral_function


!=========================================================================
function skip_transition(nspin,ib1,ib2,occ1,occ2)
 implicit none
 logical             :: skip_transition
 integer,intent(in)  :: nspin,ib1,ib2
 real(dp),intent(in) :: occ1,occ2
!=====

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

 if(allocated(sf%pole))         deallocate(sf%pole)
 if(allocated(sf%residu_left))  deallocate(sf%residu_left)
 if(allocated(sf%residu_right)) deallocate(sf%residu_right)

 WRITE_MASTER(*,'(/,a)') ' Spectral function destroyed'

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

   sf%npole      = npole_read
   sf%nprodbasis = nprodbasis_read

!   if( npole_read /= sf%npole .OR. nprodbasis_read /= sf%nprodbasis ) then
!     WRITE_MASTER(*,'(a,/)')     ' File does not have the correct size'
!     WRITE_MASTER(*,'(i5,a,i5)') npole_read,' vs ',sf%npole
!     WRITE_MASTER(*,'(i5,a,i5)') nprodbasis_read,' vs ',sf%nprodbasis
!     reading_status=2
!   else

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

!   endif
   close(spectralfile)

 endif

end subroutine read_spectral_function



end module m_spectral_function
!=========================================================================





