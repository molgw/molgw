!=========================================================================
module m_spectral_function
 use m_definitions
 use m_mpi
 use m_timing 
 use m_warning
 use m_memory
 use m_inputparam
 use m_atoms
 use m_eri,only: nauxil_2center

 !
 ! General form of any spectral function
 ! z complex number
 ! i, j running on the basis set
 ! sf_ij(z) = \sum_n L_n(i) R_n(j) / ( z - w_n )
 !

 type spectral_function 
   integer              :: nprodbasis

   !
   ! Information about all the poles (the calculated ones + the approximated ones)
   !
   integer              :: npole_reso
   integer,allocatable  :: transition_table(:,:)  ! correspondance table from
                                                  ! transition index to state pair indexes

   !
   ! Information for poles actually calculated
   !
   integer              :: npole_reso_apb
   integer,allocatable  :: transition_table_apb(:,:)  ! correspondance table from
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

 integer,parameter :: nvirtual_SPA=60

 !
 ! the boring small complex number eta: (0.0_dp,0.001_dp) is typically over converged
 ! Having a larger ieta value smoothen the oscillation far from the HOMO-LUMO gap
 complex(dp),protected :: ieta


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
 allocate(sf%transition_table(3,sf%npole_reso))
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
     enddo
   enddo
 enddo

 !
 ! The same but only for the poles that will be actually calculated by
 ! diagonalization 
 ! 

 !
 ! First, count the number of resonant transitions
 itrans=0
 do ijspin=1,nspin
   do ibf=1,nbf
     do jbf=1,nbf
       if( skip_transition_apb(nspin,jbf,ibf,occupation(jbf,ijspin),occupation(ibf,ijspin)) ) cycle
       if( occupation(jbf,ijspin) - occupation(ibf,ijspin) > 0.0_dp ) cycle
       itrans = itrans + 1
     enddo
   enddo
 enddo

 ! Set the number of poles as twice the number of resonant transtions
 sf%npole_reso_apb = itrans
 allocate(sf%transition_table_apb(3,sf%npole_reso_apb))
 ! Set the transition_table 
 itrans=0
 do ijspin=1,nspin
   do ibf=1,nbf
     do jbf=1,nbf
       if( skip_transition_apb(nspin,jbf,ibf,occupation(jbf,ijspin),occupation(ibf,ijspin)) ) cycle
       if( occupation(jbf,ijspin) - occupation(ibf,ijspin) > 0.0_dp ) cycle
       itrans = itrans + 1
       ! Set the resonant transition table
       sf%transition_table_apb(1,itrans) = ibf
       sf%transition_table_apb(2,itrans) = jbf
       sf%transition_table_apb(3,itrans) = ijspin
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
 if( ib1 <= ncore_W .OR. ib2 <= ncore_W ) skip_transition = .TRUE.
 !
 ! skip the virtual states if asked for a frozen-virtual calculation
 if( ib1 >= nvirtual_W .OR. ib2 >= nvirtual_W ) skip_transition = .TRUE.

 if( occ1 < completely_empty             .AND. occ2 < completely_empty )             skip_transition = .TRUE.
 if( occ1 > spin_fact - completely_empty .AND. occ2 > spin_fact - completely_empty ) skip_transition = .TRUE.


end function skip_transition


!=========================================================================
function skip_transition_apb(nspin,ib1,ib2,occ1,occ2)
 implicit none
 logical             :: skip_transition_apb
 integer,intent(in)  :: nspin,ib1,ib2
 real(dp),intent(in) :: occ1,occ2
!=====

 skip_transition_apb = .FALSE.

 !
 ! skip the core states if asked for a frozen-core calculation
 if( ib1 <= ncore_W .OR. ib2 <= ncore_W ) skip_transition_apb = .TRUE.
 !
 ! skip the virtual states if asked for a frozen-virtual calculation
 if( ib1 >= nvirtual_SPA .OR. ib2 >= nvirtual_SPA ) skip_transition_apb = .TRUE.

 if( occ1 < completely_empty             .AND. occ2 < completely_empty )             skip_transition_apb = .TRUE.
 if( occ1 > spin_fact - completely_empty .AND. occ2 > spin_fact - completely_empty ) skip_transition_apb = .TRUE.


end function skip_transition_apb


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
 integer             :: tmpfile
 integer             :: ipole,ibf_auxil
 integer             :: npole_write,ipole_write
 integer             :: ierr,bufsize,iproc
 logical             :: file_exists
 real(dp)            :: ecut_pole
 integer,allocatable :: index_pole(:)
 real(dp),allocatable :: buffer(:)
#ifdef HAVE_MPI
 integer(kind=MPI_OFFSET_KIND) :: disp
#endif
!=====

 write(stdout,'(/,a,/)') ' Writing the spectral function on file' 

 inquire(file='manual_poles',exist=file_exists)
 if(file_exists) then
   open(newunit=tmpfile,file='manual_poles',status='old')
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

#ifndef HAVE_MPI
 if( is_iomaster ) then
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
#else
 call MPI_FILE_OPEN(MPI_COMM_WORLD,'SCREENED_COULOMB', & 
                    MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    MPI_INFO_NULL,wfile,ierr) 

 ! Only one proc has to write the poles
 disp  = 0
 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,calc_type%postscf_name,LEN(calc_type%postscf_name),MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(calc_type%postscf_name)

 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,sf%npole_reso,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(sf%npole_reso)

 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,sf%pole,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
 disp = disp + sf%npole_reso * SIZEOF(sf%pole(1))

 !
 ! Write the residu in "the" universal ordering that does not depend on the
 ! data distribution
 allocate(buffer(sf%npole_reso))
 do ibf_auxil=1,nauxil_2center
   if( rank == iproc_ibf_auxil(ibf_auxil) ) then
     buffer(:) = sf%residu_left(ibf_auxil_l(ibf_auxil),:)
     call MPI_FILE_WRITE_AT(wfile,disp,buffer,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)


   endif
   disp = disp + sf%npole_reso * SIZEOF(sf%residu_left(1,1))
 enddo
 deallocate(buffer)

 call MPI_FILE_CLOSE(wfile, ierr)

#endif

 deallocate(index_pole)

end subroutine write_spectral_function


!=========================================================================
subroutine read_spectral_function(sf,reading_status)
 implicit none
 type(spectral_function),intent(inout) :: sf
 integer,intent(out)                   :: reading_status
!=====
 integer            :: wfile
 character(len=100) :: postscf_name_read
 integer            :: ipole_read,ibf_auxil
 logical            :: file_exists
 integer            :: npole_read,nprodbasis_read
 integer            :: ierr,bufsize,iproc
 real(dp),allocatable :: buffer(:)
#ifdef HAVE_MPI
 integer(kind=MPI_OFFSET_KIND) :: disp
#endif
!=====

 write(stdout,'(/,a)') ' Try to read spectral function from file SCREENED_COULOMB' 

 inquire(file='SCREENED_COULOMB',exist=file_exists)
 if( .NOT. file_exists ) then
   write(stdout,'(a,/)') ' File does not exist'
   reading_status=1
   return
 endif

#ifndef HAVE_MPI
 open(newunit=wfile,file='SCREENED_COULOMB',status='old',form='unformatted')

 read(wfile) postscf_name_read
 read(wfile) nprodbasis_read
 read(wfile) npole_read

 sf%npole_reso = npole_read
 call allocate_spectral_function(nprodbasis_read,sf)


 read(wfile) sf%pole(:)
 do ipole_read=1,npole_read
   read(wfile) sf%residu_left(:,ipole_read)
 enddo

 reading_status=0
 msg='reading spectral function from SCREENED_COULOMB obtained from '//TRIM(postscf_name_read)
 call issue_warning(msg)

 close(wfile)
#else
 write(stdout,*) 'Reading file'
 call MPI_FILE_OPEN(MPI_COMM_WORLD,'SCREENED_COULOMB', & 
                    MPI_MODE_RDONLY, & 
                    MPI_INFO_NULL,wfile,ierr) 
 disp = 0
 call MPI_FILE_READ_AT(wfile,disp,postscf_name_read,LEN(postscf_name_read),MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(postscf_name_read)

 call MPI_FILE_READ_AT(wfile,disp,npole_read,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(npole_read)

 sf%npole_reso = npole_read
 call allocate_spectral_function(nbf_local_iproc(rank),sf)

 call MPI_FILE_READ_AT(wfile,disp,sf%pole,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
 disp = disp + sf%npole_reso * SIZEOF(sf%pole(1))

 !
 ! Read the residu from "the" universal ordering that does not depend on the
 ! data distribution
 allocate(buffer(sf%npole_reso))
 do ibf_auxil=1,nauxil_2center
   if( rank == iproc_ibf_auxil(ibf_auxil) ) then
     call MPI_FILE_READ_AT(wfile,disp,buffer,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
     sf%residu_left(ibf_auxil_l(ibf_auxil),:) = buffer(:)
   endif
   disp = disp + sf%npole_reso * SIZEOF(sf%residu_left(1,1))
 enddo
 deallocate(buffer)


 call MPI_FILE_CLOSE(wfile, ierr)

 msg='reading spectral function from SCREENED_COULOMB obtained from '//TRIM(postscf_name_read)
 call issue_warning(msg)
 reading_status=0

#endif


end subroutine read_spectral_function


end module m_spectral_function
!=========================================================================
