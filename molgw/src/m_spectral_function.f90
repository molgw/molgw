!=========================================================================
#include "macros.h"
!=========================================================================
module m_spectral_function
 use m_definitions
 use m_mpi
 use m_timing 
 use m_warning
 use m_inputparam

 !
 ! General form of any spectral function
 ! z complex number
 ! i, j running on the basis set
 ! sf_ij(z) = \sum_n L_n(i) R_n(j) / ( z - w_n )
 !

 type spectral_function 
   integer              :: npole
   integer              :: nprodbasis
   real(dp),allocatable :: pole(:)
   real(dp),allocatable :: residu_left(:,:)       ! first index runs on n, second index on i
   real(dp),allocatable :: residu_right(:,:)      ! first index runs on n, second index on j
 end type spectral_function

 !
 ! frozen core approximation parameters
 integer,protected :: ncore_G
 integer,protected :: ncore_W

 !
 ! frozen virtual approximation parameters
 integer,protected :: nvirtual_G
 integer,protected :: nvirtual_W

 !
 ! the boring small complex number eta: (0.0_dp,0.001_dp) is typically over converged
 ! Having a larger ieta value smoothen the oscillation far from the HOMO-LUMO gap
 complex(dp),parameter :: ieta=(0.0_dp,0.01_dp) ! (0.0_dp,0.0001_dp)

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

 ncore_G    = 0
 ncore_W    = 0
 nvirtual_G = nbf+1
 nvirtual_W = nbf+1

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
   nvirtual_G = MIN(nvirtual_G,nbf+1)
   nvirtual_W = MIN(nvirtual_W,nbf+1)
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
 sf%nprodbasis = nprodbasis

 WRITE_MASTER(*,'(/,a,i8)') ' Spectral function initialized with npoles                 : ',sf%npole
 WRITE_MASTER(*,'(a,i8)')   ' Spectral function initialized with Coulomb basis functions: ',sf%nprodbasis

 allocate(sf%pole(sf%npole))
 allocate(sf%residu_left (sf%npole,sf%nprodbasis))
 allocate(sf%residu_right(sf%npole,sf%nprodbasis))

 call memory_statement(REAL(2*sf%npole,dp)*REAL(sf%nprodbasis,dp))
 WRITE_MASTER(*,*)


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
 integer,parameter   :: spectralfile=50
 integer             :: iprodbasis,ipole
 integer             :: npole_write,ipole_write
 logical             :: file_exists
 real(dp)            :: ecut_pole
 integer,allocatable :: index_pole(:)
!=====

 WRITE_MASTER(*,'(/,a,/)') ' Writing the spectral function on file' 

 inquire(file='manual_poles',exist=file_exists)
 if(file_exists) then
   open(spectralfile,file='manual_poles',status='old')
   read(spectralfile,*) ecut_pole
   if( ecut_pole<0.0_dp ) stop'error when reading manual_poles'
   close(spectralfile)
   WRITE_ME(msg,'(a,f10.4)') 'Ouput of the spectral function with an energy cutoff [eV] ',ecut_pole*Ha_eV
   call issue_warning(msg)
 else
  ecut_pole = HUGE(1.0_dp)
 endif

 npole_write = 0
 do ipole=1,sf%npole
   if( ABS(sf%pole(ipole)) < ecut_pole ) npole_write = npole_write + 1
 enddo
 WRITE_MASTER(*,*) 'Writing',npole_write,'poles over a total of',sf%npole
 allocate(index_pole(npole_write))
 ipole_write = 0
 do ipole=1,sf%npole
   if( ABS(sf%pole(ipole)) < ecut_pole ) then
     ipole_write = ipole_write + 1
     index_pole(ipole_write) = ipole
   endif
 enddo

 open(spectralfile,file='spectral_file',form='unformatted')

 if(.NOT. file_exists) then
   WRITE_MASTER(spectralfile) calc_type%postscf_name
 else
   WRITE_ME(msg,'(a,a,f10.4)') TRIM(calc_type%postscf_name),' with cutoff above energy [eV] ',ecut_pole*Ha_eV
   WRITE_MASTER(spectralfile) msg
 endif
 WRITE_MASTER(spectralfile) npole_write
 WRITE_MASTER(spectralfile) sf%nprodbasis
 WRITE_MASTER(spectralfile) sf%pole(index_pole(:))
 do iprodbasis=1,sf%nprodbasis
   WRITE_MASTER(spectralfile) sf%residu_left(index_pole(:),iprodbasis)
 enddo
 do iprodbasis=1,sf%nprodbasis
   WRITE_MASTER(spectralfile) sf%residu_right(index_pole(:),iprodbasis)
 enddo

 close(spectralfile)
 deallocate(index_pole)

end subroutine write_spectral_function


!=========================================================================
subroutine read_spectral_function(sf,reading_status)
 implicit none
 type(spectral_function),intent(inout) :: sf
 integer,intent(out)                   :: reading_status
!=====
 integer,parameter  :: spectralfile=50
 character(len=100) :: postscf_name_read
 integer            :: iprodbasis
 logical            :: file_exists
 integer            :: npole_read,nprodbasis_read
!=====

 WRITE_MASTER(*,'(/,a)') ' Try to read spectral function from file spectral_file' 

 inquire(file='spectral_file',exist=file_exists)
 if( .NOT. file_exists ) then
   WRITE_MASTER(*,'(a,/)') ' File does not exist'
   reading_status=1
   return
 endif

 open(spectralfile,file='spectral_file',status='old',form='unformatted')

 read(spectralfile) postscf_name_read
 read(spectralfile) npole_read
 read(spectralfile) nprodbasis_read

 call allocate_spectral_function(npole_read,nprodbasis_read,sf)

! if( npole_read /= sf%npole .OR. nprodbasis_read /= sf%nprodbasis ) then
!   WRITE_MASTER(*,'(a,/)')     ' File does not have the correct size'
!   WRITE_MASTER(*,'(i5,a,i5)') npole_read,' vs ',sf%npole
!   WRITE_MASTER(*,'(i5,a,i5)') nprodbasis_read,' vs ',sf%nprodbasis
!   reading_status=2
! else

   read(spectralfile) sf%pole(:)
   do iprodbasis=1,sf%nprodbasis
     read(spectralfile) sf%residu_left(:,iprodbasis)
   enddo
   do iprodbasis=1,sf%nprodbasis
     read(spectralfile) sf%residu_right(:,iprodbasis)   
   enddo

   reading_status=0
   msg='reading spectral function from spectral_file obtained from '//TRIM(postscf_name_read)
   call issue_warning(msg)

! endif
 close(spectralfile)


end subroutine read_spectral_function



end module m_spectral_function
!=========================================================================





