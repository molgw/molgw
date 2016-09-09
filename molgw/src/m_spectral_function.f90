!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the spectral decomposition of the screened Coulomb interaction W
!
!=========================================================================
module m_spectral_function
 use m_definitions
 use m_mpi
 use m_timing 
 use m_warning
 use m_memory
 use m_inputparam
 use m_atoms
 use m_eri_calculate,only: nauxil_2center

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

   !
   ! Information for poles that are not actually calculated
   !
   integer              :: npole_reso_spa
   integer,allocatable  :: transition_table_spa(:,:)  ! correspondance table from
                                                      ! transition index to state pair indexes

   real(dp),allocatable :: pole(:)
   real(dp),allocatable :: residu_left(:,:)       ! first index runs on n, second index on i

   !
   ! Static W might be stored directly in the auxiliary basis
   real(dp),allocatable :: w0(:,:)

 end type spectral_function

 !
 ! Highest occupied state
 integer,protected :: nhomo_W

 !
 ! Highest occupied state
 integer,protected :: nlumo_W

 !
 ! frozen core approximation parameter
 integer,protected :: ncore_W

 !
 ! frozen virtual approximation parameter
 integer,protected :: nvirtual_W

 !
 ! Single Pole Approximation parameter
 integer,protected :: nvirtual_SPA


contains


!=========================================================================
function index_prodstate(istate,jstate)
 implicit none
 integer,intent(in)  :: istate,jstate
 integer             :: index_prodstate
!=====
 integer             :: imin,imax
!=====

 ! Index (i,j) transformed into (I)
 ! with this ordering:
 !   (  1  2  4  7  ... )
 !   (  2  3  5  8  ... )
 !   (  4  5  6  9  ... )
 !   (  7  8  9 10  ... )
 !   ( ................ )

 imin = MIN(istate,jstate)
 imax = MAX(istate,jstate)

 index_prodstate = ( imax * (imax - 1) ) / 2 + imin


end function index_prodstate


!=========================================================================
subroutine init_spectral_function(nstate,occupation,sf)
 implicit none
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(:,:)
 type(spectral_function),intent(out)   :: sf
!=====
 integer                               :: ijspin,istate,jstate,itrans,jtrans
!=====

 if( nstate > SIZE( occupation(:,:) , DIM=1 ) ) then
   call die('init_spectral_function: nstate is too large')
 endif

 ncore_W      = ncorew
 nvirtual_W   = MIN(nvirtualw,nstate+1)
 nvirtual_SPA = MIN(nvirtualspa,nstate+1)

 if(is_frozencore) then
   if( ncore_W == 0) ncore_W = atoms_core_states()
 endif
 if( ncore_W > 0 ) then
   write(msg,'(a,i4,2x,i4)') 'frozen core approximation in W switched on up to state : ',ncore_W
   call issue_warning(msg)
 endif

 if( nvirtual_W <= nstate ) then
   write(msg,'(a,i4,2x,i4)') 'frozen virtual approximation in W switched on starting with state : ',nvirtual_W
   call issue_warning(msg)
 endif

 !
 ! Find the highest occupied state
 nhomo_W = 0
 do ijspin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ijspin) / spin_fact < completely_empty ) cycle
     nhomo_W = MAX(nhomo_W,istate)
   enddo
 enddo

 !
 ! Find the lowest occupied state
 nlumo_W = 100000
 do ijspin=1,nspin
   do istate=1,nstate
     if( (spin_fact - occupation(istate,ijspin)) / spin_fact < completely_empty) cycle
     nlumo_W = MIN(nlumo_W,istate)
   enddo
 enddo

 write(stdout,'(/,1x,a)') 'Prepare a polarizability spectral function with'
 write(stdout,'(30x,a,i6)') ' Occupied states: ',nhomo_W-ncore_W
 write(stdout,'(30x,a,i6)') '  Virtual states: ',nvirtual_W-nlumo_W
 write(stdout,'(30x,a,i6)') 'Transition space: ',(nvirtual_W-nlumo_W)*(nhomo_W-ncore_W)

 !
 ! First, count the number of resonant transitions
 itrans=0
 do ijspin=1,nspin
   do istate=1,nstate
     do jstate=1,nstate
       if( skip_transition(jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       if( occupation(jstate,ijspin) - occupation(istate,ijspin) > 0.0_dp ) cycle
       itrans = itrans + 1
     enddo
   enddo
 enddo

 sf%npole_reso = itrans  
 allocate(sf%transition_table(3,sf%npole_reso))
 ! Set the transition_table 
 itrans=0
 do ijspin=1,nspin
   do istate=1,nstate
     do jstate=1,nstate
       if( skip_transition(jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       if( occupation(jstate,ijspin) - occupation(istate,ijspin) > 0.0_dp ) cycle
       itrans = itrans + 1
       ! Set the resonant transition table
       sf%transition_table(1,itrans) = istate
       sf%transition_table(2,itrans) = jstate
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
   do istate=1,nstate
     do jstate=1,nstate
       if( skip_transition(jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       if( skip_transition_apb(jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       if( occupation(jstate,ijspin) - occupation(istate,ijspin) > 0.0_dp ) cycle
       itrans = itrans + 1
     enddo
   enddo
 enddo

 sf%npole_reso_apb = itrans
 sf%npole_reso_spa = sf%npole_reso - sf%npole_reso_apb
 allocate(sf%transition_table_apb(3,sf%npole_reso_apb))
 allocate(sf%transition_table_spa(3,sf%npole_reso_spa))
 ! Set the transition_table 
 itrans=0
 jtrans=0
 do ijspin=1,nspin
   do istate=1,nstate
     do jstate=1,nstate
       if( skip_transition(jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       if( occupation(jstate,ijspin) - occupation(istate,ijspin) > 0.0_dp ) cycle

       if( .NOT. skip_transition_apb(jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) then
         itrans = itrans + 1
         ! Set the resonant transition table
         sf%transition_table_apb(1,itrans) = istate
         sf%transition_table_apb(2,itrans) = jstate
         sf%transition_table_apb(3,itrans) = ijspin
       else
         jtrans = jtrans + 1
         ! Set the resonant transition table
         sf%transition_table_spa(1,jtrans) = istate
         sf%transition_table_spa(2,jtrans) = jstate
         sf%transition_table_spa(3,jtrans) = ijspin
       endif

     enddo
   enddo
 enddo

 if( sf%npole_reso_apb /= sf%npole_reso ) then
   write(msg,'(a,i6,2x,i6)') 'using single pole approximation with # poles instead of # ',sf%npole_reso_apb,sf%npole_reso
   call issue_warning(msg)
 endif


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
 call clean_allocate('left residu',sf%residu_left,sf%nprodbasis,sf%npole_reso)
 

end subroutine allocate_spectral_function


!=========================================================================
function skip_transition(ib1,ib2,occ1,occ2)
 implicit none
 logical             :: skip_transition
 integer,intent(in)  :: ib1,ib2
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
function skip_transition_apb(ib1,ib2,occ1,occ2)
 implicit none
 logical             :: skip_transition_apb
 integer,intent(in)  :: ib1,ib2
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

 if(ALLOCATED(sf%transition_table)) deallocate(sf%transition_table)
 if(ALLOCATED(sf%pole))             deallocate(sf%pole)
 if(ALLOCATED(sf%residu_left)) then
   call clean_deallocate('left residu',sf%residu_left)
 endif
 if(ALLOCATED(sf%w0)) then
   call clean_deallocate('static W',sf%w0)
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
 integer             :: ierr
 logical             :: file_exists
 real(dp)            :: ecut_pole
 integer,allocatable :: index_pole(:)
 real(dp),allocatable :: buffer(:)
 integer             :: iprodbasis
#ifdef HAVE_MPI
 integer(kind=MPI_OFFSET_KIND) :: disp
#endif
!=====

 write(stdout,'(/,a,/)') ' Writing the spectral function on file' 

 inquire(file='manual_poles',exist=file_exists)
 if(file_exists) then
   open(newunit=tmpfile,file='manual_poles',status='old')
   read(tmpfile,*) ecut_pole
   if( ecut_pole < 0.0_dp ) then
     call die('error when reading manual_poles')
   endif
   close(tmpfile)
   write(msg,'(a,f10.4)') 'Ouput of the spectral function with an energy cutoff (eV) ',ecut_pole*Ha_eV
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
     write(msg,'(a,a,f10.4)') TRIM(calc_type%postscf_name),' with cutoff above energy (eV) ',ecut_pole*Ha_eV
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
   call MPI_FILE_WRITE_AT(wfile,disp,sf%nprodbasis,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(sf%nprodbasis)

 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,sf%npole_reso,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(sf%npole_reso)

 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,sf%pole,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
 disp = disp + sf%npole_reso * SIZEOF(sf%pole(1))

 if( has_auxil_basis ) then
   !
   ! Write the residu in "the" universal ordering that does not depend on the
   ! data distribution
   allocate(buffer(sf%npole_reso))
   do ibf_auxil=1,nauxil_2center
     if( rank_auxil == iproc_ibf_auxil(ibf_auxil) ) then
       buffer(:) = sf%residu_left(ibf_auxil_l(ibf_auxil),:)
       call MPI_FILE_WRITE_AT(wfile,disp,buffer,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)


     endif
     disp = disp + sf%npole_reso * SIZEOF(sf%residu_left(1,1))
   enddo
   deallocate(buffer)
 else
   if(is_iomaster) then
     do iprodbasis=1,sf%nprodbasis
       call MPI_FILE_WRITE_AT(wfile,disp,sf%residu_left(iprodbasis,:),sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
       disp = disp + sf%npole_reso * SIZEOF(sf%residu_left(1,1))
     enddo
   endif
 endif

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
 integer            :: ibf_auxil
 logical            :: file_exists
 integer            :: npole_read,nprodbasis_read
 integer            :: ierr
 real(dp),allocatable :: buffer(:)
 integer            :: iprodbasis
#ifdef HAVE_MPI
 integer(kind=MPI_OFFSET_KIND) :: disp
#else
 integer :: ipole_read
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

 call MPI_FILE_READ_AT(wfile,disp,nprodbasis_read,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(nprodbasis_read)

 sf%nprodbasis = nprodbasis_read

 call MPI_FILE_READ_AT(wfile,disp,npole_read,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + SIZEOF(npole_read)

 sf%npole_reso = npole_read

 if( has_auxil_basis ) then
   call allocate_spectral_function(nbf_local_iproc(rank_auxil),sf)
 else
   call allocate_spectral_function(nprodbasis_read,sf)
 endif

 call MPI_FILE_READ_AT(wfile,disp,sf%pole,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
 disp = disp + sf%npole_reso * SIZEOF(sf%pole(1))

 if( has_auxil_basis ) then
   !
   ! Read the residu from "the" universal ordering that does not depend on the
   ! data distribution
   allocate(buffer(sf%npole_reso))
   do ibf_auxil=1,nauxil_2center
     if( rank_auxil == iproc_ibf_auxil(ibf_auxil) ) then
       call MPI_FILE_READ_AT(wfile,disp,buffer,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
       sf%residu_left(ibf_auxil_l(ibf_auxil),:) = buffer(:)
     endif
     disp = disp + sf%npole_reso * SIZEOF(sf%residu_left(1,1))
   enddo
   deallocate(buffer)
 else
   do iprodbasis=1,sf%nprodbasis
     call MPI_FILE_READ_AT(wfile,disp,sf%residu_left(iprodbasis,:),sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
     disp = disp + sf%npole_reso * SIZEOF(sf%residu_left(1,1))
   enddo
 endif



 call MPI_FILE_CLOSE(wfile, ierr)

 msg='reading spectral function from SCREENED_COULOMB obtained from '//TRIM(postscf_name_read)
 call issue_warning(msg)
 reading_status=0

#endif


end subroutine read_spectral_function


end module m_spectral_function
!=========================================================================
