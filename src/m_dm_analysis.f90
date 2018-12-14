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
subroutine dm_diff(basis)
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer                 :: nstate
 integer                 :: file_density_matrix
 integer                 :: ispin
 integer                 :: igrid,igrid_start,igrid_end,nr
 logical                 :: density_matrix_found
 real(dp),allocatable    :: p_matrix_pt(:,:,:)
 real(dp),allocatable    :: c_matrix_pt(:,:,:)
 real(dp),allocatable    :: occupation_pt(:,:)
 real(dp),allocatable    :: p_matrix_ga(:,:,:)
 real(dp),allocatable    :: c_matrix_ga(:,:,:)
 real(dp),allocatable    :: occupation_ga(:,:)
 real(dp)                :: normalization_pt(nspin),normalization_ga(nspin)
 real(dp)                :: chi1(nspin),chi2(nspin)
 real(dp),allocatable    :: rhor_batch_pt(:,:),rhor_batch_ga(:,:)
 real(dp),allocatable    :: weight_batch(:)
 real(dp),allocatable    :: basis_function_r_batch(:,:)
!=====

 nstate = basis%nbf

 ! First density matrix is tagged with suffix _ga
 if( read_fchk /= 'NO') then
   call clean_allocate('Density matrix',p_matrix_ga,basis%nbf,basis%nbf,nspin)
   call read_gaussian_fchk(basis,p_matrix_ga)
 else
   inquire(file='DENSITY_MATRIX_1',exist=density_matrix_found)
   if( density_matrix_found) then
     call clean_allocate('Density matrix',p_matrix_ga,basis%nbf,basis%nbf,nspin)
     write(stdout,'(/,1x,a)') 'Reading a MOLGW density matrix file: DENSITY_MATRIX_1'
     open(newunit=file_density_matrix,file='DENSITY_MATRIX_1',form='unformatted',action='read')
     do ispin=1,nspin
       read(file_density_matrix) p_matrix_ga(:,:,ispin)
     enddo
     close(file_density_matrix)
   else
     call die('dm_diff: no correlated density matrix read or calculated though input file suggests you really want one')
   endif
 endif

 ! Second density matrix is tagged with suffix _pt
 inquire(file='DENSITY_MATRIX',exist=density_matrix_found)
 if( density_matrix_found) then
   call clean_allocate('Density matrix',p_matrix_pt,basis%nbf,basis%nbf,nspin)
   write(stdout,'(/,1x,a)') 'Reading a MOLGW density matrix file: DENSITY_MATRIX'
   open(newunit=file_density_matrix,file='DENSITY_MATRIX',form='unformatted',action='read')
   do ispin=1,nspin
     read(file_density_matrix) p_matrix_pt(:,:,ispin)
   enddo
   close(file_density_matrix)
 else
   call die('dm_diff: no correlated density matrix read or calculated though input file suggests you really want one')
 endif


 call clean_allocate('C matrix',c_matrix_pt,basis%nbf,nstate,nspin)
 call clean_allocate('C matrix',c_matrix_ga,basis%nbf,nstate,nspin)
 allocate(occupation_pt(nstate,nspin),occupation_ga(nstate,nspin))

 call get_c_matrix_from_p_matrix(p_matrix_pt,c_matrix_pt,occupation_pt)
 call get_c_matrix_from_p_matrix(p_matrix_ga,c_matrix_ga,occupation_ga)


 call init_dft_grid(basis,grid_level,.FALSE.,.TRUE.,BATCH_SIZE)

 normalization_ga(:) = 0.0_dp
 normalization_pt(:) = 0.0_dp
 chi1(:) = 0.0_dp
 chi2(:) = 0.0_dp

 !
 ! Loop over batches of grid points
 !
 do igrid_start=1,ngrid,BATCH_SIZE
   igrid_end = MIN(ngrid,igrid_start+BATCH_SIZE-1)
   nr = igrid_end - igrid_start + 1


   allocate(weight_batch(nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(rhor_batch_pt(nspin,nr))
   allocate(rhor_batch_ga(nspin,nr))

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   call get_basis_functions_r_batch(basis,igrid_start,basis_function_r_batch)
   call calc_density_r_batch(occupation_pt,c_matrix_pt,basis_function_r_batch,rhor_batch_pt)
   call calc_density_r_batch(occupation_ga,c_matrix_ga,basis_function_r_batch,rhor_batch_ga)

   ! Normalization
   normalization_pt(:) = normalization_pt(:) + MATMUL( rhor_batch_pt(:,:) , weight_batch(:) )
   normalization_ga(:) = normalization_ga(:) + MATMUL( rhor_batch_ga(:,:) , weight_batch(:) )

   chi1(:) = chi1(:) + MATMUL( ABS( rhor_batch_pt(:,:) - rhor_batch_ga(:,:) ) , weight_batch(:) )
   chi2(:) = chi2(:) + MATMUL( ABS( rhor_batch_pt(:,:) - rhor_batch_ga(:,:) )**2 , weight_batch(:) )




   deallocate(weight_batch,rhor_batch_pt,rhor_batch_ga,basis_function_r_batch)

 enddo

 call xsum_grid(normalization_ga)
 call xsum_grid(normalization_pt)
 call xsum_grid(chi1)
 call xsum_grid(chi2)

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization_pt(:)
 write(stdout,'(a,2(2x,f12.6))')   ' Number of electrons:',normalization_ga(:)

 write(stdout,'(/,a,2(2x,f12.6))') ' Chi1:',chi1(:)
 write(stdout,'(a,2(2x,f12.6))')   ' Chi2:',SQRT( chi2(:) )



 if( print_multipole_ ) then
   !
   ! Evaluate the static dipole
   call static_dipole(nstate,basis,occupation_ga,c_matrix_ga)
   !
   ! Evaluate the static quadrupole
   call static_quadrupole(nstate,basis,occupation_ga,c_matrix_ga)

   !
   ! Evaluate the static dipole
   call static_dipole(nstate,basis,occupation_pt,c_matrix_pt)
   !
   ! Evaluate the static quadrupole
   call static_quadrupole(nstate,basis,occupation_pt,c_matrix_pt)
 endif

 if( print_cube_ ) call plot_cube_wfn('MBPT',nstate,basis,occupation_pt,c_matrix_pt)


 !
 call clean_deallocate('Density matrix',p_matrix_pt)
 call clean_deallocate('Density matrix',p_matrix_ga)
 call clean_deallocate('C matrix',c_matrix_pt)
 call clean_deallocate('C matrix',c_matrix_ga)
 deallocate(occupation_pt,occupation_ga)

 write(stdout,'(/,1x,a,/)') 'This is the end'

 call die('dm_diff: exited normally')


end subroutine dm_diff


end module m_dm_analysis


!=========================================================================
