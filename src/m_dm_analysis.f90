!=========================================================================
! This file is part of MOLGW.
!
! Copyright (C) 2010-2016  Fabien Bruneval
! MOLGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOLGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MOLGW.  If not, see <http://www.gnu.org/licenses/>.
!=========================================================================
!
! This is the main of MOLGW
!
! It consists of 3 parts:
!   1. Initialization
!   2. SCF cycles
!   3. Post-processings: GW, MP2, CI, TDDFT etc.
!
!=========================================================================
module m_dm_analysis
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_scalapack
 use m_inputparam
 use m_tools
 use m_scf
 use m_atoms
 use m_ecp
 use m_gaussian
 use m_basis_set
 use m_lbfgs
 use m_eri
 use m_eri_calculate
 use m_eri_ao_mo
 use m_dft_grid
 use m_spectral_function
 use m_hamiltonian
 use m_hamiltonian_sca
 use m_hamiltonian_onebody
 use m_hamiltonian_buffer
 use m_selfenergy_tools
 use m_virtual_orbital_space


 integer,parameter,private :: BATCH_SIZE = 128


contains


!=========================================================================
subroutine dm_dump(basis)
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer                 :: nstate
 integer                 :: file_density_matrix,file_rho_grid
 integer                 :: ispin
 integer                 :: igrid,igrid_start,igrid_end,nr,ir
 logical                 :: density_matrix_found
 real(dp),allocatable    :: p_matrix_test(:,:,:)
 real(dp),allocatable    :: c_matrix_test(:,:,:)
 real(dp),allocatable    :: occupation_test(:,:)
 real(dp)                :: normalization_test(nspin)
 real(dp)                :: chi1(nspin),chi2(nspin)
 real(dp),allocatable    :: rhor_batch_test(:,:)
 real(dp),allocatable    :: weight_batch(:)
 real(dp),allocatable    :: basis_function_r_batch(:,:)
!=====

 write(stdout,'(/,1x,a)') 'Dump the electronic density into a file'

 if( nproc_grid > 1 ) call die('dm_dump: not coded in parallel. Run with 1 core only')

 nstate = basis%nbf


 ! density matrix is tagged with suffix _test
 call clean_allocate('Density matrix',p_matrix_test,basis%nbf,basis%nbf,nspin)
 if( read_fchk /= 'NO') then
   call read_gaussian_fchk(read_fchk,'gaussian.fchk',basis,p_matrix_test)
 else
   inquire(file='DENSITY_MATRIX',exist=density_matrix_found)
   if( density_matrix_found) then
     write(stdout,'(/,1x,a)') 'Reading a MOLGW density matrix file: DENSITY_MATRIX'
     open(newunit=file_density_matrix,file='DENSITY_MATRIX',form='unformatted',action='read')
     do ispin=1,nspin
       read(file_density_matrix) p_matrix_test(:,:,ispin)
     enddo
     close(file_density_matrix)
   else
     call die('dm_dump: no correlated density matrix read or calculated though input file suggests you really want one')
   endif
 endif


 call clean_allocate('C matrix',c_matrix_test,basis%nbf,nstate,nspin)
 allocate(occupation_test(nstate,nspin))

 call get_c_matrix_from_p_matrix(p_matrix_test,c_matrix_test,occupation_test)

 write(*,*) 'Sum(occupation): ',SUM(occupation_test(:,:))

 call init_dft_grid(basis,grid_level,.FALSE.,.TRUE.,BATCH_SIZE)

 normalization_test(:) = 0.0_dp

 !
 ! Loop over batches of grid points
 !
 open(newunit=file_rho_grid,file='rho_grid.dat',form='formatted',action='write')
 do igrid_start=1,ngrid,BATCH_SIZE
   igrid_end = MIN(ngrid,igrid_start+BATCH_SIZE-1)
   nr = igrid_end - igrid_start + 1


   allocate(weight_batch(nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(rhor_batch_test(nspin,nr))

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   call get_basis_functions_r_batch(basis,igrid_start,basis_function_r_batch)
   call calc_density_r_batch(occupation_test,c_matrix_test,basis_function_r_batch,rhor_batch_test)

   ! Normalization
   normalization_test(:) = normalization_test(:) + MATMUL( rhor_batch_test(:,:) , weight_batch(:) )


   if( nspin == 1 ) then
     write(file_rho_grid,'(2(1x,es16.6))') (rhor_batch_test(:,ir),weight_batch(ir),ir=1,nr)
   else
     write(file_rho_grid,'(3(1x,es16.6))') (rhor_batch_test(:,ir),weight_batch(ir),ir=1,nr)
   endif

   deallocate(weight_batch,rhor_batch_test,basis_function_r_batch)

 enddo
 close(file_rho_grid)

 ! not coded in parallel
 !call xsum_grid(normalization_test)

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization_test(:)


 if( print_multipole_ ) then
   !
   ! Evaluate the static dipole
   call static_dipole(nstate,basis,occupation_test,c_matrix_test)
   !
   ! Evaluate the static quadrupole
   call static_quadrupole(nstate,basis,occupation_test,c_matrix_test)
 endif

 if( print_cube_ ) call plot_cube_wfn('MBPT',nstate,basis,occupation_test,c_matrix_test)


 !
 call clean_deallocate('Density matrix',p_matrix_test)
 call clean_deallocate('C matrix',c_matrix_test)
 deallocate(occupation_test)


 ! Cleanly exit the code
 call stop_clock(timing_prescf)
 call stop_clock(timing_total)

 call total_memory_statement()

 call output_timing()

 call output_all_warnings()

 write(stdout,'(/,1x,a,/)') 'This is the end'

 call finish_mpi()

 stop

end subroutine dm_dump


end module m_dm_analysis


!=========================================================================
