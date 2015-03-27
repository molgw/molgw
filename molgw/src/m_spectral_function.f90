!=========================================================================
module m_spectral_function
 use m_definitions
 use m_mpi
 use m_timing 
 use m_warning
 use m_memory
 use m_inputparam
 use m_atoms

 !
 ! General form of any spectral function
 ! z complex number
 ! i, j running on the basis set
 ! sf_ij(z) = \sum_n L_n(i) R_n(j) / ( z - w_n )
 !

 type spectral_function 
   integer              :: npole
   integer              :: npole_reso
   integer              :: nprodbasis
   real(dp),allocatable :: transition_table(:,:)  ! correspondance table from
                                                  ! transition index to state pair indexes
   real(dp),allocatable :: pole(:)
   real(dp),allocatable :: residu_left(:,:)       ! first index runs on n, second index on i
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
 complex(dp),protected :: ieta

#ifdef CRPA
 integer,parameter :: band1=1
 integer,parameter :: band2=2
#endif

contains

!=========================================================================
subroutine init_spectral_function(nbf,occupation,sf)
 implicit none
 integer,intent(in)                    :: nbf
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 type(spectral_function),intent(out)   :: sf
!=====
 integer                               :: ijspin,ibf,jbf,itrans
 logical                               :: file_exists
!====

 ieta = (0.0_dp,1.0_dp) * pole_eta 

 ncore_G    = ncoreg
 ncore_W    = ncorew
 nvirtual_G = MIN(nvirtualg,nbf+1)
 nvirtual_W = MIN(nvirtualw,nbf+1)

 if(is_frozencore) then
   if( ncore_G == 0) ncore_G = atoms_core_states()
   if( ncore_W == 0) ncore_W = atoms_core_states()
 endif
 if( ncore_G > 0 .OR. ncore_W > 0 ) then
   write(msg,'(a,i4,2x,i4)') 'frozen core approximation switched on up to state (G,W) = ',ncore_G,ncore_W
   call issue_warning(msg)
 endif

 if( nvirtual_G <= nbf .OR. nvirtual_W <= nbf ) then
   write(msg,'(a,i4,2x,i4)') 'frozen virtual approximation switched on starting with state (G,W) = ',nvirtual_G,nvirtual_W
   call issue_warning(msg)
 endif

 !
 ! First, count the number of resonant transitions
 itrans=0
 do ijspin=1,nspin
   do ibf=1,nbf
     do jbf=1,nbf
       if( skip_transition(nspin,jbf,ibf,occupation(jbf,ijspin),occupation(ibf,ijspin)) ) cycle
       if( occupation(jbf,ijspin) - occupation(ibf,ijspin) > 0.0_dp ) cycle
       itrans = itrans + 1
     enddo
   enddo
 enddo

 ! Set the number of poles as twice the number of resonant transtions
 sf%npole_reso = itrans  
 sf%npole      = 2*itrans  
 allocate(sf%transition_table(3,sf%npole))
 ! Set the transition_table 
 itrans=0
 do ijspin=1,nspin
   do ibf=1,nbf
     do jbf=1,nbf
       if( skip_transition(nspin,jbf,ibf,occupation(jbf,ijspin),occupation(ibf,ijspin)) ) cycle
       if( occupation(jbf,ijspin) - occupation(ibf,ijspin) > 0.0_dp ) cycle
       itrans = itrans + 1
       ! Set the resonant transition table
       sf%transition_table(1,itrans) = ibf
       sf%transition_table(2,itrans) = jbf
       sf%transition_table(3,itrans) = ijspin
       ! Set the anti-resonant transition table too
       sf%transition_table(1,itrans+sf%npole/2) = jbf
       sf%transition_table(2,itrans+sf%npole/2) = ibf
       sf%transition_table(3,itrans+sf%npole/2) = ijspin
     enddo
   enddo
 enddo


end subroutine init_spectral_function


!=========================================================================
subroutine allocate_spectral_function(nprodbasis,sf)
 implicit none
 integer,intent(in)                    :: nprodbasis
 type(spectral_function),intent(inout) :: sf
!=====

 sf%nprodbasis = nprodbasis

 write(stdout,'(/,a,i8)') ' Spectral function initialized with Coulomb basis functions: ',sf%nprodbasis
 write(stdout,'(a,i8)')   ' Spectral function initialized with resonant poles         : ',sf%npole_reso

 allocate(sf%pole(sf%npole_reso))
 call clean_allocate(' left residu',sf%residu_left,sf%nprodbasis,sf%npole_reso)
 

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

 if( occ1 < completely_empty             .AND. occ2 < completely_empty )             skip_transition=.TRUE.
 if( occ1 > spin_fact - completely_empty .AND. occ2 > spin_fact - completely_empty ) skip_transition=.TRUE.

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

 if(allocated(sf%transition_table)) deallocate(sf%transition_table)
 if(allocated(sf%pole))             deallocate(sf%pole)
 if(allocated(sf%residu_left)) then
   call clean_deallocate(' left residu',sf%residu_left)
 endif

 write(stdout,'(/,a)') ' Spectral function destroyed'

end subroutine destroy_spectral_function


!=========================================================================
subroutine write_spectral_function(sf)
 implicit none
 type(spectral_function),intent(in) :: sf
!=====
 integer             :: wfile
 integer,parameter   :: tmpfile=51
 integer             :: iprodbasis,ipole
 integer             :: npole_write,ipole_write
 logical             :: file_exists
 real(dp)            :: ecut_pole
 integer,allocatable :: index_pole(:)
!=====

 write(stdout,'(/,a,/)') ' Writing the spectral function on file' 

 inquire(file='manual_poles',exist=file_exists)
 if(file_exists) then
   open(tmpfile,file='manual_poles',status='old')
   read(tmpfile,*) ecut_pole
   if( ecut_pole<0.0_dp ) stop'error when reading manual_poles'
   close(tmpfile)
   write(msg,'(a,f10.4)') 'Ouput of the spectral function with an energy cutoff [eV] ',ecut_pole*Ha_eV
   call issue_warning(msg)
 else
  ecut_pole = HUGE(1.0_dp)
 endif

 npole_write = 0
 do ipole=1,sf%npole_reso
   if( ABS(sf%pole(ipole)) < ecut_pole ) npole_write = npole_write + 1
 enddo
 write(stdout,*) 'Writing',npole_write,'poles over a total of',sf%npole_reso
 allocate(index_pole(npole_write))
 ipole_write = 0
 do ipole=1,sf%npole_reso
   if( ABS(sf%pole(ipole)) < ecut_pole ) then
     ipole_write = ipole_write + 1
     index_pole(ipole_write) = ipole
   endif
 enddo

 if( is_iomaster() ) then
   open(newunit=wfile,file='SCREENED_COULOMB',form='unformatted')

   if(.NOT. file_exists) then
     write(wfile) calc_type%postscf_name
   else
     write(msg,'(a,a,f10.4)') TRIM(calc_type%postscf_name),' with cutoff above energy [eV] ',ecut_pole*Ha_eV
     write(wfile) msg
   endif
   write(wfile) sf%nprodbasis
   write(wfile) npole_write
   write(wfile) sf%pole(index_pole(:))
   do ipole_write=1,npole_write
     write(wfile) sf%residu_left(:,index_pole(ipole_write))
   enddo

   close(wfile)
 endif

 deallocate(index_pole)

end subroutine write_spectral_function


!=========================================================================
subroutine read_spectral_function(sf,reading_status)
 implicit none
 type(spectral_function),intent(inout) :: sf
 integer,intent(out)                   :: reading_status
!=====
 integer,parameter  :: wfile=50
 character(len=100) :: postscf_name_read
 integer            :: ipole_read
 logical            :: file_exists
 integer            :: npole_read,nprodbasis_read
!=====

 write(stdout,'(/,a)') ' Try to read spectral function from file SCREENED_COULOMB' 

 inquire(file='SCREENED_COULOMB',exist=file_exists)
 if( .NOT. file_exists ) then
   write(stdout,'(a,/)') ' File does not exist'
   reading_status=1
   return
 endif

 open(wfile,file='SCREENED_COULOMB',status='old',form='unformatted')

 read(wfile) postscf_name_read
 read(wfile) nprodbasis_read
 read(wfile) npole_read

 sf%npole_reso = npole_read
 sf%npole      = npole_read * 2
 call allocate_spectral_function(nprodbasis_read,sf)

! if( npole_read /= sf%npole .OR. nprodbasis_read /= sf%nprodbasis ) then
!   write(stdout,'(a,/)')     ' File does not have the correct size'
!   write(stdout,'(i5,a,i5)') npole_read,' vs ',sf%npole
!   write(stdout,'(i5,a,i5)') nprodbasis_read,' vs ',sf%nprodbasis
!   reading_status=2
! else

   read(wfile) sf%pole(:)
   do ipole_read=1,npole_read
     read(wfile) sf%residu_left(:,ipole_read)
   enddo

   reading_status=0
   msg='reading spectral function from SCREENED_COULOMB obtained from '//TRIM(postscf_name_read)
   call issue_warning(msg)

! endif
 close(wfile)


end subroutine read_spectral_function


end module m_spectral_function
!=========================================================================
