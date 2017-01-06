!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! propagation of the wavefunction in time
!
!=========================================================================
!module m_tddft_propagator
! use m_definitions
! use m_basis_set 
! use m_scf_loop
! use m_memory
! use m_hamiltonian
!contains


!=========================================================================
subroutine calculate_propagation(nstate, basis, occupation, energy, s_matrix, c_matrix,hamiltonian_kinetic,hamiltonian_nucleus)
 use m_definitions
 use m_basis_set
 use m_scf_loop
 use m_memory
 use m_hamiltonian
 use m_inputparam
 use m_dft_grid
 implicit none

!=====
 type(basis_set),intent(in)      :: basis
 integer                         :: nstate 
 real(dp),intent(in)             :: energy(nstate,nspin)
 real(dp),intent(in)             :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: c_matrix(basis%nbf,nstate,nspin) 
 real(dp),intent(in)             :: occupation(nstate,nspin)
 real(dp),intent(in)             :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer                    :: itau, jtau, ntau, idir, info
 real(dp)                   :: t_min, t_max, t_step, t_cur
 real(dp),allocatable       :: s_matrix_inv(:,:)
 real(dp),allocatable       :: p_matrix(:,:,:)
 complex(dpc), allocatable  :: l_matrix(:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al, 
 complex(dpc),allocatable   :: b_matrix(:,:) ! TDDFT Book, Springer (2006), p205
 complex(dpc),allocatable   :: unity(:,:)
 complex(dpc),allocatable   :: hamiltonian_hxc_cmplx(:,:,:)  
 complex(dpc),allocatable   :: hamiltonian_fock_cmplx(:,:,:)
 complex(dpc),allocatable   :: c_matrix_cmplx(:,:,:)
 complex(dpc),allocatable   :: p_matrix_cmplx(:,:,:)
 call init_dft_grid(grid_level)
 
 t_min=0.0_dp
 t_max=30.0_dp
 t_step=0.5_dp

allocate(s_matrix_inv(basis%nbf,basis%nbf))
allocate(p_matrix(basis%nbf,basis%nbf,nspin))
allocate(hamiltonian_hxc_cmplx(basis%nbf,basis%nbf,nspin))
allocate(hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin))
allocate(c_matrix_cmplx(basis%nbf,nstate,nspin))
allocate(p_matrix_cmplx(basis%nbf,basis%nbf,nspin))
allocate(l_matrix(basis%nbf,basis%nbf,nspin))
allocate(b_matrix(basis%nbf,basis%nbf,nspin))

 hamiltonian_fock_cmplx(:,:,1) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
   if(nspin==2) hamiltonian_fock_cmplx(:,:,nspin) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)


  c_matrix_cmplx(:,:,:)=c_matrix(:,:,:) 
  p_matrix_cmplx(:,:,:)=p_matrix(:,:,:)
  call calculate_hamiltonian_hxc_ri_cmplx(basis,nstate,basis%nbf,basis%nbf,basis%nbf,nstate,occupation, &
     c_matrix_cmplx,hamiltonian_hxc_cmplx)

  hamiltonian_fock_cmplx(:,:,:)=hamiltonian_fock_cmplx(:,:,:)+hamiltonian_hxc_cmplx(:,:,:)


 call invert(basis%nbf,s_matrix,s_matrix_inv)  
 call fill_unity(unity,basis%nbf)

!********start time loop*************
ntau=int((t_max-t_min)/t_step)
do itau=1,ntau
  t_cur=t_min+itau*t_step

  l_matrix = unity + im * t_step / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,1)  )
  b_matrix = unity - im * t_step / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,1)  )
  call invert_inplace(basis%nbf , l_matrix)
  c_matrix_cmplx = matmul(l_matrix,matmul(b_matrix,c_matrix_cmplx))

  call calculate_hamiltonian_hxc_ri_cmplx(basis,nstate,basis%nbf,basis%nbf,basis%nbf,nstate,occupation, &
             c_matrix_cmplx,hamiltonian_hxc_cmplx)
  
  
end do
!********end time loop*******************





 call destroy_dft_grid()
end subroutine calculate_propagation

subroutine fill_unity(unity,M)
        use m_definitions
        implicit none
        integer, intent(in) :: M
        complex(dpc) :: unity(M,M)
        integer :: i,j
        do i=1,M
         do j=1,M
             if (i == j) then
                     unity(i,j)=1.0_dp
             else
                     unity(i,j)=0.0_dp
             end if
          end do
        end do
end subroutine

!=========================================================================
!end module m_tddft_propagator
!=========================================================================

