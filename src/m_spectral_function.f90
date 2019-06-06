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
 use m_scalapack
 use m_inputparam
 use m_eri,only: iproc_ibf_auxil,ibf_auxil_l
 use m_eri_calculate,only: nauxil_2center,nauxil_3center
 use m_numerical_tools,only: coeffs_gausslegint

 !
 ! General form of any spectral function
 ! z complex number
 ! i, j running on the basis set
 ! sf_ij(z) = \sum_n L_n(i) R_n(j) / ( z - w_n )
 !

 type spectral_function

   integer              :: nprodbasis_total       ! total over all procs
   integer              :: nprodbasis             ! for this proc

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
   real(dp),allocatable :: residue_left(:,:)       ! first index runs on n, second index on i

   !
   ! Static or Dynamic W might be stored directly in the auxiliary basis
   real(dp),allocatable :: chi(:,:,:)
   integer              :: mchi,nchi
   integer              :: desc_chi(NDEL)
   integer              :: nomega_quad           ! Number of quadrature points
   real(dp),allocatable :: omega_quad(:)         ! quadrature points for numerical integration
   real(dp),allocatable :: weight_quad(:)        ! quadrature weight for numerical integration

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
pure function index_prodstate(istate,jstate)
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
subroutine init_spectral_function(nstate,occupation,nomega_in,sf)
 implicit none
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(:,:)
 integer,intent(in)                    :: nomega_in
 type(spectral_function),intent(out)   :: sf
!=====
 integer                               :: ijspin,istate,jstate,itrans,jtrans
 integer                               :: iomega
 integer                               :: nlumo_W_spin(nspin)
 integer                               :: nhomo_W_spin(nspin)
 real(dp),parameter                    :: alpha=1.0_dp ! 0.50_dp
 real(dp),parameter                    :: beta=1.0_dp ! 6.0_dp
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
   write(msg,'(a,i4,2x,i4)') 'frozen core approximation in W switched on up to state = ',ncore_W
   call issue_warning(msg)
 endif

 if( nvirtual_W <= nstate ) then
   write(msg,'(a,i4,2x,i4)') 'frozen virtual approximation in W switched on starting with state = ',nvirtual_W
   call issue_warning(msg)
 endif

 !
 ! Find the highest occupied state
 nhomo_W         = 0
 nhomo_W_spin(:) = 0
 do ijspin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ijspin) / spin_fact < completely_empty ) cycle
     nhomo_W              = MAX(nhomo_W,istate)
     nhomo_W_spin(ijspin) = MAX(nhomo_W_spin(ijspin),istate)
   enddo
 enddo

 !
 ! Find the lowest occupied state
 nlumo_W         = 100000
 nlumo_W_spin(:) = 100000
 do ijspin=1,nspin
   do istate=1,nstate
     if( (spin_fact - occupation(istate,ijspin)) / spin_fact < completely_empty) cycle
     nlumo_W              = MIN(nlumo_W,istate)
     nlumo_W_spin(ijspin) = MIN(nlumo_W_spin(ijspin),istate)
   enddo
 enddo

 write(stdout,'(/,1x,a)') 'Prepare a polarizability spectral function with'
 if( nspin == 1 ) then
   write(stdout,'(30x,a,i8)') ' Occupied states: ',nhomo_W-ncore_W
   write(stdout,'(30x,a,i8)') '  Virtual states: ',nvirtual_W-nlumo_W
   write(stdout,'(30x,a,i8)') 'Transition space: ',(nvirtual_W-nlumo_W)*(nhomo_W-ncore_W)
 else
   write(stdout,'(30x,a,i8,2x,i8)') ' Occupied states: ',nhomo_W_spin(:)-ncore_W
   write(stdout,'(30x,a,i8,2x,i8)') '  Virtual states: ',nvirtual_W-nlumo_W_spin(:)
   write(stdout,'(30x,a,i8)')       'Transition space: ',(nvirtual_W-nlumo_W_spin(1))*(nhomo_W_spin(1)-ncore_W) &
                                                         + (nvirtual_W-nlumo_W_spin(nspin))*(nhomo_W_spin(nspin)-ncore_W)
 endif

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

 if( has_auxil_basis ) then
   sf%nprodbasis_total = nauxil_2center
 else
   sf%nprodbasis_total = index_prodstate(nvirtual_W-1,nvirtual_W-1) * nspin
 endif

 !
 ! Set the sampling points for Chi (quadrature)
 sf%nomega_quad    = nomega_in

 if( sf%nomega_quad > 0 ) then
   allocate(sf%weight_quad(sf%nomega_quad))
   allocate(sf%omega_quad(sf%nomega_quad))

   if( sf%nomega_quad > 1 ) then
     call coeffs_gausslegint(0.0_dp,1.0_dp,sf%omega_quad,sf%weight_quad,sf%nomega_quad)

     write(stdout,'(/,1x,a)') 'Numerical integration on a grid along the imaginary axis'
     ! Variable change [0,1] -> [0,+\inf[
     write(stdout,'(a)') '    #    Frequencies (eV)    Quadrature weights'
     do iomega=1,sf%nomega_quad
       sf%weight_quad(iomega) = sf%weight_quad(iomega) / ( 2.0_dp**alpha - 1.0_dp ) * alpha * (1.0_dp -  sf%omega_quad(iomega))**(-alpha-1.0_dp) * beta
       sf%omega_quad(iomega)  =   1.0_dp / ( 2.0_dp**alpha - 1.0_dp ) * ( 1.0_dp / (1.0_dp-sf%omega_quad(iomega))**alpha - 1.0_dp ) * beta
       write(stdout,'(i5,2(2x,f14.6))') iomega,sf%omega_quad(iomega)*Ha_eV,sf%weight_quad(iomega)
     enddo
   else
     ! Only one omega means static case
     sf%weight_quad(1) = 1.0_dp
     sf%omega_quad(1)  = 0.0_dp
   endif
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
 call clean_allocate('Left residue',sf%residue_left,sf%nprodbasis,sf%npole_reso)


end subroutine allocate_spectral_function


!=========================================================================
pure function skip_transition(ib1,ib2,occ1,occ2)
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
pure function skip_transition_apb(ib1,ib2,occ1,occ2)
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
 if(ALLOCATED(sf%residue_left)) then
   call clean_deallocate('Left residue',sf%residue_left)
 endif
 if(ALLOCATED(sf%chi)) then
   call clean_deallocate('Chi',sf%chi)
 endif
 if(ALLOCATED(sf%weight_quad)) deallocate(sf%weight_quad)
 if(ALLOCATED(sf%omega_quad))  deallocate(sf%omega_quad)

 write(stdout,'(/,a)') ' Spectral function destroyed'

end subroutine destroy_spectral_function


!=========================================================================
subroutine write_spectral_function(sf)
 implicit none
 type(spectral_function),intent(in) :: sf
!=====
 integer              :: wfile
 integer              :: ibf_auxil
 integer              :: ierr
 real(dp),allocatable :: buffer(:)
 integer              :: iprodbasis
#ifdef HAVE_MPI
 integer(kind=MPI_OFFSET_KIND) :: disp
#endif
!=====
 integer :: ipole
!=====

 write(stdout,'(/,a,/)') ' Writing the spectral function on file: SCREENED_COULOMB'

 write(stdout,*) 'Number of poles written down:',sf%npole_reso


#ifndef HAVE_MPI
 if( is_iomaster ) then
   open(newunit=wfile,file='SCREENED_COULOMB',form='unformatted')

   write(wfile) calc_type%postscf_name
   write(wfile) sf%nprodbasis_total
   write(wfile) sf%npole_reso
   write(wfile) sf%pole(:)
   do ipole=1,sf%npole_reso
     write(wfile) sf%residue_left(:,ipole)
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
 disp = disp + STORAGE_SIZE(calc_type%postscf_name)

 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,sf%nprodbasis_total,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + STORAGE_SIZE(sf%nprodbasis_total)

 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,sf%npole_reso,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + STORAGE_SIZE(sf%npole_reso)

 if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile,disp,sf%pole,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
 disp = disp + sf%npole_reso * STORAGE_SIZE(sf%pole(1))


 if( has_auxil_basis ) then
   !
   ! Write the residue in "the" universal ordering that does not depend on the
   ! data distribution
   allocate(buffer(sf%npole_reso))
   do ibf_auxil=1,sf%nprodbasis_total
     if( rank_auxil == iproc_ibf_auxil(ibf_auxil) ) then

       buffer(:) = sf%residue_left(ibf_auxil_l(ibf_auxil),:)
       call MPI_FILE_WRITE_AT(wfile,disp,buffer,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

     endif
     disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1,1))
   enddo
   deallocate(buffer)
 else
   if(is_iomaster) then
     do iprodbasis=1,sf%nprodbasis_total
       call MPI_FILE_WRITE_AT(wfile,disp,sf%residue_left(iprodbasis,:),sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
       disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1,1))
     enddo
   endif
 endif

 call MPI_FILE_CLOSE(wfile, ierr)

#endif


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

#if !defined(HAVE_MPI)
 open(newunit=wfile,file='SCREENED_COULOMB',status='old',form='unformatted')

 read(wfile) postscf_name_read
 read(wfile) nprodbasis_read
 read(wfile) npole_read

 sf%npole_reso = npole_read
 call allocate_spectral_function(nprodbasis_read,sf)


 read(wfile) sf%pole(:)
 do ipole_read=1,npole_read
   read(wfile) sf%residue_left(:,ipole_read)
 enddo

 reading_status=0
 msg='reading spectral function from SCREENED_COULOMB obtained from '//TRIM(postscf_name_read)
 call issue_warning(msg)

 close(wfile)
#else

 write(stdout,*) 'Reading file SCREENED_COULOMB'
 call MPI_FILE_OPEN(MPI_COMM_WORLD,'SCREENED_COULOMB', &
                    MPI_MODE_RDONLY, &
                    MPI_INFO_NULL,wfile,ierr)

 disp = 0
 call MPI_FILE_READ_AT(wfile,disp,postscf_name_read,LEN(postscf_name_read),MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
 disp = disp + STORAGE_SIZE(postscf_name_read)

 call MPI_FILE_READ_AT(wfile,disp,nprodbasis_read,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + STORAGE_SIZE(nprodbasis_read)

 sf%nprodbasis_total = nprodbasis_read

 call MPI_FILE_READ_AT(wfile,disp,npole_read,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
 disp = disp + STORAGE_SIZE(npole_read)

 sf%npole_reso = npole_read

 if( has_auxil_basis ) then
   call allocate_spectral_function(nauxil_3center,sf)
 else
   call allocate_spectral_function(nprodbasis_read,sf)
 endif

 call MPI_FILE_READ_AT(wfile,disp,sf%pole,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
 disp = disp + sf%npole_reso * STORAGE_SIZE(sf%pole(1))


 if( has_auxil_basis ) then
   !
   ! Read the residue from "the" universal ordering that does not depend on the
   ! data distribution
   allocate(buffer(sf%npole_reso))
   do ibf_auxil=1,nauxil_2center
     if( rank_auxil == iproc_ibf_auxil(ibf_auxil) ) then
       call MPI_FILE_READ_AT(wfile,disp,buffer,sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
       sf%residue_left(ibf_auxil_l(ibf_auxil),:) = buffer(:)
     endif
     disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1,1))
   enddo
   deallocate(buffer)
 else
   do iprodbasis=1,sf%nprodbasis
     call MPI_FILE_READ_AT(wfile,disp,sf%residue_left(iprodbasis,:),sf%npole_reso,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
     disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1,1))
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
