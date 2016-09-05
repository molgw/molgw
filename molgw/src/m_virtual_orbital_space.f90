!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! subroutine to generate clever representations of the virtual orbital space
!
!=========================================================================
module m_virtual_orbital_space
 use m_definitions
 use m_timing
 use m_mpi
 use m_warning
 use m_memory

 real(dp),allocatable,private :: energy_ref(:,:)
 real(dp),allocatable,private :: c_matrix_ref(:,:,:)


contains

!=========================================================================
subroutine setup_virtual_smallbasis(basis,nstate,occupation,nsemax,energy,c_matrix,nstate_small)
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_hamiltonian
 use m_hamiltonian_sca
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 integer,intent(in)                    :: nsemax
 real(dp),intent(inout)                :: energy(nstate,nspin)
 real(dp),intent(inout)                :: c_matrix(basis%nbf,nstate,nspin)
 integer,intent(out)                   :: nstate_small
!=====
 integer                               :: ispin
 integer                               :: ibf
 integer                               :: istate,jstate
 integer                               :: nstate_new
 type(basis_set)                       :: basis_small
 real(dp),allocatable                  :: s_bigsmall(:,:)
 real(dp),allocatable                  :: s_small(:,:)
 real(dp),allocatable                  :: s_small_sqrt_inv(:,:)
 real(dp),allocatable                  :: h_small(:,:,:)
 real(dp),allocatable                  :: energy_small(:,:)
 real(dp),allocatable                  :: c_small(:,:,:)
 real(dp),allocatable                  :: c_big(:,:,:)
 real(dp),allocatable                  :: s_matrix(:,:)
 real(dp),allocatable                  :: h_big(:,:,:)
 real(dp),allocatable                  :: s_matrix_sqrt_inv(:,:)
 real(dp),allocatable                  :: matrix_tmp(:,:)
 real(dp),allocatable                  :: s_bar(:,:),h_bar(:,:,:),s_bar_sqrt_inv(:,:)
 real(dp),allocatable                  :: energy_bar(:,:),c_bar(:,:,:)
 integer                               :: nstate_bar
 integer                               :: nocc,nfrozen
!=====

 call start_clock(timing_fno)

 write(stdout,'(/,1x,a)') 'Prepare optimized empty states using a smaller basis set'

 call assert_experimental()

 ! Remember how to go from the small basis set to the big one
 ! 
 ! | \phi^small_a > = \sum_{bc} | \phi^big_b > * S^-1_ac * Sbs_cb
 !
 ! This is key to everything else!


 ! Initialize the small basis set
 write(stdout,*) 'Set up a smaller basis set to optimize the virtual orbital space'
 call init_basis_set(basis_path,small_basis_name,gaussian_type,basis_small)
 call issue_warning('Reduce the virtual orbital subspace by using a smaller basis set: '//TRIM(small_basis_name))

 ! Get the overlap matrix of the wavefunction basis set S: s_matrix
 allocate(s_matrix(basis%nbf,basis%nbf))
 call setup_overlap(.FALSE.,basis,s_matrix)

 ! Calculate the mixed overlap matrix Sbs: s_bigsmall
 allocate(s_bigsmall(basis%nbf,basis_small%nbf))
 call setup_overlap_mixedbasis(.FALSE.,basis,basis_small,s_bigsmall)

 call setup_sqrt_overlap(min_overlap,basis%nbf,s_matrix,nstate_new,s_matrix_sqrt_inv)
 if( nstate_new /= nstate) call die('setup_virtual_smallbasis: bug this should not happen')

 ! Calculate the overlap matrix in the small basis:
 !  tilde S = Sbs^T S^-1 Sbs
 allocate(matrix_tmp(nstate_new,basis_small%nbf))
 allocate(s_small(basis_small%nbf,basis_small%nbf))
 matrix_tmp(:,:) = MATMUL( TRANSPOSE( s_matrix_sqrt_inv ) , s_bigsmall )
 s_small(:,:) = MATMUL( TRANSPOSE( matrix_tmp ) , matrix_tmp )
 deallocate(matrix_tmp)

 ! Calculate ( tilde S )^{-1/2}
 call setup_sqrt_overlap(min_overlap,basis_small%nbf,s_small,nstate_small,s_small_sqrt_inv)


 ! Obtain the Hamiltonian in the big basis (it will be useful later)
 !
 ! H = S * C * E * C^T * S
 !
 allocate(h_big(basis%nbf,basis%nbf,nspin))
 allocate(matrix_tmp(basis%nbf,basis%nbf))
 do ispin=1,nspin
   
   ! M = E * C^T
   do istate=1,nstate
     matrix_tmp(istate,:) = energy(istate,ispin) * c_matrix(:,istate,ispin)
   enddo
   
   ! M = E * C^T * S
   matrix_tmp(1:nstate,:) = MATMUL( matrix_tmp(1:nstate,:) , s_matrix(:,:) )

   ! M = C * E * C^T * S
   matrix_tmp(:,:) = MATMUL( c_matrix(:,:,ispin) , matrix_tmp(1:nstate,:) )
 
   h_big(:,:,ispin) = MATMUL( s_matrix(:,:) , matrix_tmp(:,:) )
 enddo
 deallocate(matrix_tmp)


 ! Obtain the Hamiltonian in the small basis (it is useful now!)
 !
 ! tilde H = Sbs^T * S^-1 * H * S^-1 * Sbs
 !         = Sbs^T *    C * E * C^T  * Sbs
 !
 allocate(h_small(basis_small%nbf,basis_small%nbf,nspin))
 allocate(matrix_tmp(basis%nbf,basis%nbf))
 do ispin=1,nspin
   
   ! M = E * C^T
   do istate=1,nstate
     matrix_tmp(istate,:) = energy(istate,ispin) * c_matrix(:,istate,ispin)
   enddo
   
   ! M = E * C^T * S
   matrix_tmp(1:nstate,1:basis_small%nbf) = MATMUL( matrix_tmp(1:nstate,:) , s_bigsmall(:,:) )

   ! M = C * E * C^T * S
   matrix_tmp(:,1:basis_small%nbf) = MATMUL( c_matrix(:,:,ispin) , matrix_tmp(1:nstate,1:basis_small%nbf) )
 
   h_small(:,:,ispin) = MATMUL( TRANSPOSE(s_bigsmall(:,:)) , matrix_tmp(:,1:basis_small%nbf) )
 enddo
 deallocate(matrix_tmp)

 ! Diagonalize the small Hamiltonian in the small basis
 !
 ! tilde H * tilde C = tilde S * tilde C * tilde E
 !
 allocate(energy_small(nstate_small,nspin))
 allocate(c_small(basis_small%nbf,nstate_small,nspin))
 call diagonalize_hamiltonian_scalapack(nspin,basis_small%nbf,nstate_small,h_small,s_small_sqrt_inv,energy_small,c_small)
 call dump_out_energy('=== Energies in the initial small basis ===',&
              nstate_small,nspin,occupation(1:nstate_small,:),energy_small)

 deallocate(h_small)
 deallocate(energy_small)
 deallocate(s_small_sqrt_inv)


 !
 ! Transform the wavefunction coefficients from the small basis to the big basis
 ! tilde C -> Cbig
 allocate(c_big(basis%nbf,nstate_small,nspin))
 allocate(matrix_tmp(basis%nbf,basis%nbf))
 ! M = S^-1
 matrix_tmp(:,:) = MATMUL( s_matrix_sqrt_inv(:,:) , TRANSPOSE( s_matrix_sqrt_inv(:,:) ) )

 ! Cbig = S^-1 * Sbs * tilde C
 do ispin=1,nspin
   c_big(:,:,ispin) = MATMUL( matrix_tmp(:,:) , MATMUL( s_bigsmall(:,:) , c_small(:,:,ispin) ) )
 enddo
 deallocate(c_small)
 deallocate(matrix_tmp)
 deallocate(s_matrix_sqrt_inv)

 ! 
 ! Frozen orbitals for occupied state plus the selfenergy braket
 !
 ! Find the highest occupied state
 do istate=1,nstate
   if( ALL(occupation(istate,:) < completely_empty) ) cycle
   nocc = istate
 enddo

 ! Override the Cbig coefficients with the original C coefficients up to max(nocc,nsemax)
 nfrozen = MAX(nocc,nsemax)
 write(stdout,'(1x,a,i6)') 'Leave the first states frozen up to: ',nfrozen
 c_big(:,1:nfrozen,:) = c_matrix(:,1:nfrozen,:)

 !
 ! Final diagonalization of in the composite basis
 ! with frozen orbitals complemented with the small basis
 !
 ! Calculate the corresponding overlap matrix Sbar and hamiltonian Hbar
 allocate(s_bar(basis_small%nbf,basis_small%nbf))
 allocate(h_bar(basis_small%nbf,basis_small%nbf,nspin))
 s_bar(:,:) = MATMUL( TRANSPOSE(c_big(:,:,1)) , MATMUL( s_matrix , c_big(:,:,1) ) )
 call setup_sqrt_overlap(min_overlap,basis_small%nbf,s_bar,nstate_bar,s_bar_sqrt_inv)
 do ispin=1,nspin
   h_bar(:,:,ispin) = MATMUL( TRANSPOSE(c_big(:,:,ispin)) , MATMUL( h_big(:,:,ispin) , c_big(:,:,ispin) ) )
 enddo
 allocate(energy_bar(nstate_small,nspin))
 allocate(c_bar(basis_small%nbf,nstate_bar,nspin))
 call diagonalize_hamiltonian_scalapack(nspin,basis_small%nbf,nstate_bar,h_bar,s_bar_sqrt_inv,energy_bar,c_bar)

 do ispin=1,nspin
   c_big(:,:,ispin) = MATMUL( c_big(:,:,ispin) , c_bar(:,:,ispin) )
 enddo

 call dump_out_energy('=== Energies in the final small basis ===',&
              nstate_bar,nspin,occupation(1:nstate_bar,:),energy_bar)

 deallocate(c_bar)
 deallocate(s_bar)
 deallocate(h_bar)
 deallocate(s_bar_sqrt_inv)
 deallocate(h_big)
 deallocate(s_matrix)

 
 !
 ! Save the original c_matrix and energies
 if( .NOT. ALLOCATED(energy_ref) ) then
   allocate(energy_ref(nstate_bar,nspin))
   allocate(c_matrix_ref(basis%nbf,nstate_bar,nspin))   ! TODO clean_allocate of this

   energy_ref(:,:)     = energy(1:nstate_bar,:)
   c_matrix_ref(:,:,:) = c_matrix(:,1:nstate_bar,:)
 endif

 !
 ! And then override the c_matrix and the energy with the fresh new ones
 energy(1:nstate_bar,:)     = energy_bar(:,:)
 c_matrix(:,1:nstate_bar,:) = c_big(:,:,:)

 nstate_small = nstate_bar

 deallocate(energy_bar)
 deallocate(s_small)
 deallocate(c_big)
 deallocate(s_bigsmall)

 write(stdout,'(1x,a)') 'Optimized empty states with a smaller basis set'

 call stop_clock(timing_fno)

end subroutine setup_virtual_smallbasis


!=========================================================================
subroutine virtual_fno(basis,nstate,occupation,energy,c_matrix)
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(inout)                :: energy(nstate,nspin)
 real(dp),intent(inout)                :: c_matrix(basis%nbf,nstate,nspin)
!=====
 real(dp),parameter                    :: alpha_mp2=1.0_dp
 integer                               :: istate,jstate,astate,bstate,cstate
 integer                               :: ispin
 integer                               :: nocc,ncore,nvirtual,nvirtual_kept
 real(dp),allocatable                  :: p_matrix_mp2(:,:),ham_virtual(:,:),ham_virtual_kept(:,:)
 real(dp),allocatable                  :: occupation_mp2(:),energy_virtual_kept(:)
 real(dp)                              :: eri_ci_aj,eri_ci_bj
 real(dp)                              :: den_ca_ij,den_cb_ij
 real(dp)                              :: en_mp2
 integer                               :: nvirtualmin
 real(dp),allocatable                  :: eri_ci_i(:)
!=====

 call start_clock(timing_fno)

 write(stdout,'(/,1x,a)') 'Prepare optimized empty states with Frozen Natural Orbitals'

 call assert_experimental()

 nvirtualmin = MIN(nvirtualw,nvirtualg)

 if( nvirtualmin > nstate ) then 
   call issue_warning('virtual_fno is on, however nvirtualw and nvirtualg are not set. Skipping the Frozen Natural Orbitals generation.')
   return
 endif

 if( .NOT. has_auxil_basis ) then
   call issue_warning('virtual_fno not implemented when no auxiliary basis is provided')
   return
 endif

 if( nspin > 1 ) then
   call issue_warning('virtual_fno not implemented yet for spin polarized calculations')
   return
 endif

 ncore = ncoreg
 if(is_frozencore) then
   if( ncore == 0) ncore = atoms_core_states()
 endif


 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,1,nstate,1,nstate)


 do ispin=1,nspin
   !
   ! First, set up the list of occupied states
   nocc = ncore
   do istate=ncore+1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     nocc = istate
   enddo



   nvirtual = nstate - nocc

   allocate(p_matrix_mp2(nvirtual,nvirtual))
   p_matrix_mp2(:,:) = 0.0_dp

   allocate(eri_ci_i(nocc+1:nstate))

#if 1
   ! Approximation by Aquilante et al.
   do istate=ncore+1,nocc
     do cstate=nocc+1,nstate

       do astate=nocc+1,nstate
         eri_ci_i(astate) = eri_eigen_ri_paral(cstate,istate,ispin,astate,istate,ispin) &
                              / ( energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(cstate,ispin) )
       enddo
       call xsum_auxil(eri_ci_i)

       do bstate=nocc+1,nstate
         do astate=nocc+1,nstate

           p_matrix_mp2(astate-nocc,bstate-nocc) = &
                p_matrix_mp2(astate-nocc,bstate-nocc)  &
                   + 0.50_dp * eri_ci_i(astate) * eri_ci_i(bstate)

         enddo
       enddo

     enddo
   enddo
   deallocate(eri_ci_i)
#else

   do bstate=nocc+1,nstate
     do astate=nocc+1,nstate

#if 0
       ! Full calculation of the MP2 density matrix on virtual orbitals (See Taube and Bartlett)
       do cstate=nocc+1,nstate
         do istate=ncore+1,nocc
           do jstate=ncore+1,nocc

             den_cb_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(bstate,ispin) - energy(cstate,ispin)
             den_ca_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(astate,ispin) - energy(cstate,ispin)

             eri_ci_aj = eri_eigen_ri(cstate,istate,ispin,astate,jstate,ispin) &
                          - eri_eigen_ri(cstate,jstate,ispin,astate,istate,ispin)  * alpha_mp2
             eri_ci_bj = eri_eigen_ri(cstate,istate,ispin,bstate,jstate,ispin) &
                         - eri_eigen_ri(cstate,jstate,ispin,bstate,istate,ispin)   * alpha_mp2

             p_matrix_mp2(astate-nocc,bstate-nocc) = &
                  p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                     + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

           enddo
         enddo
       enddo
#elif 1
       ! Approximation by Aquilante et al.
       do cstate=nocc+1,nstate
         do istate=ncore+1,nocc
           den_cb_ij = energy(istate,ispin) + energy(istate,ispin) - energy(bstate,ispin) - energy(cstate,ispin)
           den_ca_ij = energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(cstate,ispin)

           eri_ci_aj = eri_eigen_ri(cstate,istate,ispin,astate,istate,ispin) 
           eri_ci_bj = eri_eigen_ri(cstate,istate,ispin,bstate,istate,ispin) 

           p_matrix_mp2(astate-nocc,bstate-nocc) = &
                p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                   + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

         enddo
       enddo
#elif 0
       do istate=ncore+1,nocc
         ! Approximation by Aquilante et al.
         den_cb_ij = energy(istate,ispin) + energy(istate,ispin) - energy(bstate,ispin) - energy(bstate,ispin)
         den_ca_ij = energy(istate,ispin) + energy(istate,ispin) - energy(astate,ispin) - energy(astate,ispin)

         eri_ci_aj = eri_eigen_ri(astate,istate,ispin,astate,istate,ispin) 
         eri_ci_bj = eri_eigen_ri(bstate,istate,ispin,bstate,istate,ispin) 

         p_matrix_mp2(astate-nocc,bstate-nocc) = &
              p_matrix_mp2(astate-nocc,bstate-nocc)  & 
                 + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )
       enddo
#else
         do istate=ncore+1,nocc
           do jstate=ncore+1,nocc

             den_cb_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(bstate,ispin) - energy(bstate,ispin)
             den_ca_ij = energy(istate,ispin) + energy(jstate,ispin) - energy(astate,ispin) - energy(astate,ispin)

             eri_ci_aj = eri_eigen_ri(astate,istate,ispin,astate,jstate,ispin) 
             eri_ci_bj = eri_eigen_ri(bstate,istate,ispin,bstate,jstate,ispin) 

             p_matrix_mp2(astate-nocc,bstate-nocc) = &
                  p_matrix_mp2(astate-nocc,bstate-nocc)  &
                     + 0.50_dp * eri_ci_aj * eri_ci_bj / ( den_cb_ij * den_ca_ij )

           enddo
         enddo

#endif

     enddo
   enddo


#endif

   allocate(occupation_mp2(nvirtual))
   call diagonalize(nvirtual,p_matrix_mp2,occupation_mp2)

!   write(stdout,*) 
!   do astate=1,nvirtual
!     write(stdout,'(i4,2x,es16.4)') astate,occupation_mp2(astate)
!   enddo
!   write(stdout,*) 

   allocate(ham_virtual(nvirtual,nvirtual))

   ham_virtual(:,:) = 0.0_dp
   do astate=1,nvirtual
     ham_virtual(astate,astate) = energy(astate+nocc,ispin)
   enddo

   nvirtual_kept = nvirtualmin - 1 - nocc

   write(stdout,'(/,1x,a,i5)')    'Max state index included ',nvirtual_kept + nocc
   write(stdout,'(1x,a,i5,a,i5)') 'Retain ',nvirtual_kept,' virtual orbitals out of ',nvirtual
   write(stdout,'(1x,a,es14.6)')  'Occupation number of the last virtual state',occupation_mp2(nvirtual - (nocc + nvirtual_kept))
   write(stdout,'(1x,a,es14.6,4x,es14.6)') '  to be compared to the first and last virtual states',occupation_mp2(nvirtual),occupation_mp2(1)
   write(stdout,*)

   allocate(ham_virtual_kept(nvirtual_kept,nvirtual_kept))
   allocate(energy_virtual_kept(nvirtual_kept))

   ham_virtual_kept(:,:)  = MATMUL( TRANSPOSE( p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) ) ,   &
                                        MATMUL( ham_virtual , p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) ) )

   deallocate(ham_virtual)

   call diagonalize(nvirtual_kept,ham_virtual_kept,energy_virtual_kept)

!   write(stdout,'(/,1x,a)') ' virtual state    FNO energy (eV)   reference energy (eV)'
!   do astate=1,nvirtual_kept
!     write(stdout,'(i4,4x,f16.5,2x,f16.5)') astate,energy_virtual_kept(astate)*Ha_eV,energy(astate+nocc,ispin)*Ha_eV
!   enddo
!   write(stdout,*)
   

   !
   ! Save the original c_matrix and energies
   if( .NOT. ALLOCATED(energy_ref) ) then
     allocate(energy_ref(nocc+1:nocc+nvirtual_kept,nspin))
     allocate(c_matrix_ref(basis%nbf,nocc+1:nocc+nvirtual_kept,nspin))   ! TODO clean_allocate of this

     energy_ref(nocc+1:nocc+nvirtual_kept,:)     = energy(nocc+1:nocc+nvirtual_kept,:)
     c_matrix_ref(:,nocc+1:nocc+nvirtual_kept,:) = c_matrix(:,nocc+1:nocc+nvirtual_kept,:)
   endif

   !
   ! And then override the c_matrix and the energy with the fresh new ones
   energy(nocc+1:nocc+nvirtual_kept,ispin)     = energy_virtual_kept(:)
   c_matrix(:,nocc+1:nocc+nvirtual_kept,ispin) = MATMUL( c_matrix(:,nocc+1:,ispin) ,  &
                                                    MATMUL( p_matrix_mp2(:,nvirtual-nvirtual_kept+1:) , &
                                                       ham_virtual_kept(:,:) ) )


   deallocate(energy_virtual_kept)
   deallocate(ham_virtual_kept)
   deallocate(p_matrix_mp2,occupation_mp2)
 enddo



 call destroy_eri_3center_eigen()

 write(stdout,'(1x,a)') 'Optimized empty states with Frozen Natural Orbitals'

 call stop_clock(timing_fno)

end subroutine virtual_fno


!=========================================================================
subroutine destroy_fno(basis,nstate,energy,c_matrix)
 use m_inputparam
 use m_tools,only: diagonalize
 use m_basis_set
 use m_eri_ao_mo
 implicit none

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(inout)                :: energy(nstate,nspin)
 real(dp),intent(inout)                :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                               :: lowerb,upperb
!=====

 write(stdout,'(/,1x,a)') 'Deallocate the Frozen Natural Orbitals'

 lowerb = LBOUND(energy_ref,DIM=1)
 upperb = UBOUND(energy_ref,DIM=1)


 energy    (lowerb:upperb,:) = energy_ref    (lowerb:upperb,:)
 c_matrix(:,lowerb:upperb,:) = c_matrix_ref(:,lowerb:upperb,:)


 deallocate(energy_ref)
 deallocate(c_matrix_ref)   ! TODO clean_allocate of this


end subroutine destroy_fno

!=========================================================================
end module m_virtual_orbital_space
