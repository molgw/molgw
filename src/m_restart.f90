!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the reading / writing of the RESTART files
!
!=========================================================================


!=========================================================================
subroutine write_restart(restart_type,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
 use m_definitions
 use m_timing
 use m_mpi
 use m_inputparam
 use m_atoms
 use m_basis_set
 implicit none

 integer,intent(in)         :: restart_type
 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: hamiltonian_fock(basis%nbf,basis%nbf,nspin)
!=====
 integer,parameter          :: restart_version=201609
 integer                    :: restartfile
 integer                    :: ispin,istate,ibf,nstate_local
!=====

 !
 ! Only the proc iomaster writes down the RESTART file
 if( .NOT. is_iomaster) return

 call start_clock(timing_restart_file)

 select case(restart_type)
 case(SMALL_RESTART)
   write(stdout,'(/,a)') ' Writing a small RESTART file'
 case(BIG_RESTART)
   write(stdout,'(/,a)') ' Writing a big RESTART file'
 end select

 open(newunit=restartfile,file='RESTART',form='unformatted',action='write')

 ! An integer to ensure backward compatibility in the future
 write(restartfile) restart_version
 ! RESTART file type SMALL_RESTART=1 or BIG_RESTART=2
 write(restartfile) restart_type
 ! Atomic structure
 write(restartfile) natom
 write(restartfile) zatom(1:natom)
 write(restartfile) xatom(:,1:natom)
 ! Calculation type
 write(restartfile) calc_type%scf_name
 ! Basis set
 call write_basis_set(restartfile,basis)
 ! Spin channels
 write(restartfile) nspin
 ! Nstate
 write(restartfile) nstate
 ! Occupations
 write(restartfile) occupation(:,:)
 ! Eigen energies
 write(restartfile) energy(:,:)

 ! Number of states written down in the RESTART file
 if( restart_type == SMALL_RESTART ) then
   ! Identify the highest occupied state in order to
   ! save I-O in SMALL_RESTART
   nstate_local = 0
   do ispin=1,nspin
     do istate=1,nstate
       if( ANY( occupation(istate,:) > completely_empty ) ) nstate_local = istate
     enddo
   enddo
 else
   ! Or write down all the states in BIG_RESTART
   nstate_local = nstate
 endif

 write(restartfile) nstate_local

 ! Wavefunction coefficients C
 do ispin=1,nspin
   do istate=1,nstate_local
     write(restartfile) c_matrix(:,istate,ispin)
   enddo
 enddo

 if(restart_type == BIG_RESTART) then

   do ispin=1,nspin
     do ibf=1,basis%nbf
       write(restartfile) hamiltonian_fock(:,ibf,ispin)
     enddo
   enddo

 endif

 close(restartfile)
 call stop_clock(timing_restart_file)

end subroutine write_restart


!=========================================================================
subroutine read_restart(restart_type,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock,restart_filename)
 use m_definitions
 use m_mpi
 use m_inputparam
 use m_atoms
 use m_basis_set
 use m_hamiltonian_onebody
 use m_linear_algebra,only: invert
 implicit none

 integer,intent(out)         :: restart_type
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(inout)      :: occupation(nstate,nspin)
 real(dp),intent(out)        :: c_matrix(basis%nbf,nstate,nspin),energy(nstate,nspin)
 real(dp),intent(out)        :: hamiltonian_fock(basis%nbf,basis%nbf,nspin)
 character(len=*),intent(in) :: restart_filename
!=====
 integer                    :: restartfile
 integer                    :: ispin,istate,ibf,nstate_local
 logical                    :: file_exists,same_scf,same_basis,same_geometry

 integer                    :: restart_version_read
 integer                    :: restart_type_read
 character(len=100)         :: scf_name_read
 integer                    :: natom_read
 real(dp),allocatable       :: zatom_read(:),x_read(:,:)
 type(basis_set)            :: basis_read
 integer                    :: nspin_read
 integer                    :: nstate_read
 real(dp),allocatable       :: occupation_read(:,:)
 real(dp),allocatable       :: energy_read(:,:)
 real(dp),allocatable       :: c_matrix_read(:,:,:)
 real(dp),allocatable       :: hfock_read(:,:,:)
 real(dp),allocatable       :: overlap_mixedbasis(:,:)
 real(dp),allocatable       :: overlap_smallbasis(:,:)
 real(dp),allocatable       :: overlap_bigbasis_approx(:,:)
!=====

 inquire(file=restart_filename,exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,1x,a,1x,a)') TRIM(restart_filename),'file not found'
   restart_type = NO_RESTART
   return
 endif


 open(newunit=restartfile,file=restart_filename,form='unformatted',status='old',action='read')


 ! An integer to ensure backward compatibility in the future
 read(restartfile) restart_version_read

 if( restart_version_read == 201602 ) then
   write(stdout,'(1x,a,i8)') 'Old RESTART file found. Version: ',restart_version_read
   call issue_warning('RESTART file: Old version. Backward compatibility is not ensured. Skipping the reading')
   restart_type = NO_RESTART
   close(restartfile)
   return
 endif

 if( restart_version_read /= 201609 ) then
   call issue_warning('RESTART file: Version not readable. Skipping the reading')
   restart_type = NO_RESTART
   close(restartfile)
   return
 endif


 ! RESTART file type SMALL_RESTART=1 or BIG_RESTART=2
 read(restartfile) restart_type_read
 if( restart_type_read /= SMALL_RESTART .AND. restart_type_read /= BIG_RESTART ) then
   call issue_warning('RESTART file: Type not readable. Skipping the reading')
   restart_type = NO_RESTART
   close(restartfile)
   return
 endif
 restart_type = restart_type_read

 !
 ! Input keyword ignore_bigrestart enforces a small_restart
 if( ignore_bigrestart_ ) then
   restart_type = SMALL_RESTART
 endif


 ! Atomic structure
 read(restartfile) natom_read
 allocate(zatom_read(natom_read),x_read(3,natom_read))
 read(restartfile) zatom_read(1:natom_read)
 read(restartfile) x_read(:,1:natom_read)
 if( natom_read /= natom  &
  .OR. ANY( ABS( zatom_read(1:MIN(natom_read,natom)) - zatom(1:MIN(natom_read,natom)) ) > 1.0e-5_dp ) &
  .OR. ANY( ABS(   x_read(:,1:MIN(natom_read,natom)) - xatom(:,1:MIN(natom_read,natom))   ) > 1.0e-5_dp ) ) then
   same_geometry = .FALSE.
   call issue_warning('RESTART file: Geometry has changed')
 else
   same_geometry = .TRUE.
 endif
 deallocate(zatom_read,x_read)


 ! Calculation type
 read(restartfile) scf_name_read
 same_scf = ( TRIM(scf_name_read) == TRIM(calc_type%scf_name) )
 if( .NOT. same_scf) then
   call issue_warning('RESTART file: SCF type has changed')
 endif


 ! Basis set
 call read_basis_set(restartfile,basis_read)
 same_basis = compare_basis_set(basis,basis_read)
 if( .NOT. same_basis) then
   call issue_warning('RESTART file: Basis set has changed')
 endif
 if( basis%gaussian_type /= basis_read%gaussian_type ) then
   write(stdout,*) 'The basis type (cartesian or pure) cannot be changed when restarting from a previous calculation'
   call die('Erase the RESTART file or change the keyword gaussian_type and start the calculation over')
 endif


 ! Spin channels
 read(restartfile) nspin_read
 if( nspin /= nspin_read ) then
   call issue_warning('RESTART file: Number of spin channels has changed')
   call die('not implemented yet')
 endif


 ! Nstate
 read(restartfile) nstate_read
 if( nstate /= nstate_read ) then
   call issue_warning('RESTART file: Number of states has changed')
 endif


 ! Occupations
 allocate(occupation_read(nstate_read,nspin_read))
 read(restartfile) occupation_read(:,:)
 if( ANY( ABS( occupation_read(1:MIN(nstate_read,nstate),:) - occupation(1:MIN(nstate_read,nstate),:) )  > 1.0e-5_dp ) ) then
   if( temperature > 1.0e-8_dp) then
     occupation(1:MIN(nstate_read,nstate),:)=occupation_read(1:MIN(nstate_read,nstate),:)
     write(stdout,'(1xa)') "Reading occupations from a RESTART file"
     call dump_out_occupation('=== Occupations ===',nstate,nspin,occupation)
   else
     call issue_warning('RESTART file: Occupations have changed')
   endif
 endif
 deallocate(occupation_read)


 ! Eigen energies
 allocate(energy_read(nstate_read,nspin_read))
 read(restartfile) energy_read(:,:)
 energy(:,:) = 1000.0_dp
 energy(1:MIN(nstate,nstate_read),:) = energy_read(1:MIN(nstate,nstate_read),:)
 deallocate(energy_read)


 ! Number of states written down in the RESTART file
 read(restartfile) nstate_local


 ! Wavefunction coefficients C
 allocate(c_matrix_read(basis_read%nbf,nstate_local,nspin_read))
 do ispin=1,nspin_read
   do istate=1,nstate_local
     read(restartfile) c_matrix_read(:,istate,ispin)
   enddo
 enddo
 c_matrix(:,:,:) = 0.0_dp
 do ispin=1,nspin_read
   do istate=1,MIN(nstate_local,nstate)
     c_matrix(1:MIN(basis_read%nbf,basis%nbf),istate,ispin) &
          = c_matrix_read(1:MIN(basis_read%nbf,basis%nbf),istate,ispin)
   enddo
 enddo
 ! Fill the rest of the array with identity
 if( nstate_local < nstate ) then
   do ispin=1,nspin_read
     do istate=nstate_local+1,nstate
       c_matrix(istate,istate,ispin) = 1.0_dp
     enddo
   enddo
 endif


 if(restart_type == SMALL_RESTART) then

   close(restartfile)
   return

 else

   if( same_basis ) then

     do ispin=1,nspin_read
       do ibf=1,basis_read%nbf
         read(restartfile) hamiltonian_fock(:,ibf,ispin)
       enddo
     enddo

     if( same_scf .AND. same_geometry ) then
       restart_type = BIG_RESTART
     else
       restart_type = EMPTY_STATES_RESTART
     endif
     close(restartfile)
     return

   else

     allocate(hfock_read(basis_read%nbf,basis_read%nbf,nspin_read))

     do ispin=1,nspin_read
       do ibf=1,basis_read%nbf
         read(restartfile) hfock_read(:,ibf,ispin)
       enddo
     enddo

     !
     ! If the basis sets differ, then one needs to transpose the parts of the hamiltonian
     ! for a clever starting point.

     allocate(overlap_smallbasis(basis_read%nbf,basis_read%nbf))
     allocate(overlap_bigbasis_approx(basis%nbf,basis%nbf))
     allocate(overlap_mixedbasis(basis_read%nbf,basis%nbf))

     ! Calculate the overlap matrix of the read basis set, named "smallbasis"
     call setup_overlap_mixedbasis(basis_read,basis_read,overlap_smallbasis)

     ! Invert the overlap of the read basis set
     call invert(overlap_smallbasis)

     ! Get the scalar products between the old and the new basis sets
     ! Be aware: this is a rectangular matrix
     call setup_overlap_mixedbasis(basis_read,basis,overlap_mixedbasis)

     !
     ! Evaluate the estimated overlap matrix with the small basis set that is read
     ! in the RESTART file
     ! Since the diagonal of the overlap matrix should be 1 by definition, this
     ! allows us to evaluate how wrong we are and disregard the terms in the
     ! hamiltonian which are poorly evaluated in the small basis set.
     overlap_bigbasis_approx(:,:) = MATMUL( TRANSPOSE( overlap_mixedbasis), MATMUL( overlap_smallbasis , overlap_mixedbasis ) )

     overlap_mixedbasis(:,:) = MATMUL( overlap_smallbasis(:,:)  , overlap_mixedbasis(:,:) )

     do ispin=1,nspin_read
       hamiltonian_fock(:,:,ispin) = MATMUL( TRANSPOSE(overlap_mixedbasis) , MATMUL( hfock_read(:,:,ispin), overlap_mixedbasis ) )
     enddo

     !
     ! Disregard the term that are too wrong = evaluated overlap below 0.99
     ! by giving a huge value to the Hartree
     ! part of the Hamiltonian
     do ibf=1,basis%nbf
       if( ABS(overlap_bigbasis_approx(ibf,ibf) - 1.0_dp ) > 0.010_dp ) then
         hamiltonian_fock(ibf,ibf,:) = 0.0_dp
       endif
     enddo
     deallocate(overlap_smallbasis)
     deallocate(overlap_mixedbasis)
     deallocate(overlap_bigbasis_approx)

     deallocate(hfock_read)

     restart_type = BASIS_RESTART
     close(restartfile)
     return

   endif


 endif

 ! the code should never reach that point.
 call die('Internal error in the read_restart subroutine')


end subroutine read_restart


!=========================================================================
