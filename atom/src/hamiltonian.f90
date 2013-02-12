#include "macros.h"
!=========================================================================
subroutine setup_overlap(print_volume,basis,s_matrix)
 use m_definitions
 use m_basis_set
 implicit none
 integer,intent(in)         :: print_volume
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(basis%nbf,basis%nbf)
!====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: overlap_tmp
!====

 WRITE_MASTER(*,*) 'Setup overlap matrix S'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am(CARTESIAN,li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am(CARTESIAN,lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call overlap_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo
     s_matrix(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                                                   MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) )


     deallocate(matrix_cart)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='=== Overlap matrix S ==='
 call dump_out_matrix(print_volume,title,basis%nbf,1,s_matrix)


end subroutine setup_overlap

!=========================================================================
subroutine setup_kinetic(print_volume,basis,hamiltonian_kinetic)
 use m_definitions
 use m_basis_set
 implicit none
 integer,intent(in)         :: print_volume
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(basis%nbf,basis%nbf)
!====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: kinetic_tmp
!====

 WRITE_MASTER(*,*) 'Setup kinetic part of the Hamiltonian'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am(CARTESIAN,li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am(CARTESIAN,lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call kinetic_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo
     hamiltonian_kinetic(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) )


     deallocate(matrix_cart)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='===  Kinetic energy contribution ==='
 call dump_out_matrix(print_volume,title,basis%nbf,1,hamiltonian_kinetic)


end subroutine setup_kinetic

!=========================================================================
subroutine setup_nucleus(print_volume,basis,hamiltonian_nucleus)
 use m_definitions
 use m_basis_set
 use m_atoms
 implicit none
 integer,intent(in)         :: print_volume
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
!====

 WRITE_MASTER(*,*) 'Setup nucleus-electron part of the Hamiltonian'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am(CARTESIAN,li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am(CARTESIAN,lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     matrix_cart(:,:) = 0.0_dp
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         do iatom=1,natom
           call nucleus_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),zatom(iatom),x(:,iatom),vnucleus_ij)
           matrix_cart(i_cart,j_cart) = matrix_cart(i_cart,j_cart) + vnucleus_ij
         enddo
       enddo
     enddo
     hamiltonian_nucleus(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) ) 


     deallocate(matrix_cart)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo


 title='===  Nucleus potential contribution ==='
 call dump_out_matrix(print_volume,title,basis%nbf,1,hamiltonian_nucleus)


end subroutine setup_nucleus

!=========================================================================
subroutine setup_hartree(print_volume,nbf,nspin,p_matrix,pot_hartree,ehartree)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_hartree(nbf,nbf,nspin)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 WRITE_MASTER(*,*) 'Calculate Hartree term'
 call start_clock(timing_hartree)

 pot_hartree(:,:,:)=0.0_dp

 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
   do jbf=1,nbf
     do ibf=1,nbf
       if( negligible_basispair(ibf,jbf) ) cycle
       if( .NOT. is_my_task(ibf,jbf) ) cycle
       do lbf=1,nbf
         !
         ! symmetry k <-> l
         do kbf=1,lbf-1 ! nbf
           if( negligible_basispair(kbf,lbf) ) cycle
           !
           ! symmetry (ij|kl) = (kl|ij) has been used to loop in the fast order
           pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                      + eri(kbf,lbf,ibf,jbf) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
         enddo
         pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                    + eri(lbf,lbf,ibf,jbf) * SUM( p_matrix(lbf,lbf,:) )
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo
 call xsum(pot_hartree)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_hartree)

 ehartree = 0.5_dp*SUM(pot_hartree(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_hartree)

end subroutine setup_hartree

!=========================================================================
subroutine setup_exchange(print_volume,nbf,nspin,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 real(dp)             :: spin_fact
 character(len=100)   :: title
!=====

 WRITE_MASTER(*,*) 'Calculate Exchange term'
 call start_clock(timing_exchange)

 spin_fact = REAL(-nspin+3,dp)

 pot_exchange(:,:,:)=0.0_dp

#if 0
 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO 
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         !
         ! symmetry k <-> l
         do kbf=1,lbf-1 ! nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           pot_exchange(ibf,jbf,ispin) = pot_exchange(ibf,jbf,ispin) &
                      - ( eri(kbf,ibf,lbf,jbf) + eri(lbf,ibf,kbf,jbf) ) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
         pot_exchange(ibf,jbf,ispin) = pot_exchange(ibf,jbf,ispin) &
                    - eri(lbf,ibf,lbf,jbf) * p_matrix(lbf,lbf,ispin) / spin_fact
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo
#else
 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO 
   do jbf=1,nbf
     do lbf=1,nbf
       if( negligible_basispair(lbf,jbf) ) cycle
       do ibf=1,nbf
         do kbf=1,nbf
           if( negligible_basispair(kbf,ibf) ) cycle
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           pot_exchange(ibf,jbf,ispin) = pot_exchange(ibf,jbf,ispin) &
                      - eri(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo
#endif

 title='=== Exchange contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange

!=========================================================================
subroutine setup_exchange_longrange(print_volume,nbf,nspin,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 real(dp)             :: spin_fact
 character(len=100)   :: title
!=====

 WRITE_MASTER(*,*) 'Calculate Long-Range Exchange term'
 call start_clock(timing_exchange)

 spin_fact = REAL(-nspin+3,dp)

 pot_exchange(:,:,:)=0.0_dp

 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           pot_exchange(ibf,jbf,ispin) = pot_exchange(ibf,jbf,ispin) &
                      - eri_lr(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo

 title='=== Exchange contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange

!=========================================================================
subroutine read_potential(print_volume,nbf,nspin,p_matrix,pot_read,eread)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_read(nbf,nbf,nspin)
 real(dp),intent(out) :: eread
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
 logical              :: file_exists
!=====

 pot_read(:,:,:)=0.0_dp

 inquire(file='manual_potential',exist=file_exists)
 if(file_exists) then
   open(unit=12,file='manual_potential',status='old')
   do ispin=1,nspin
     do jbf=1,nbf
       do ibf=1,nbf
         read(12,*) pot_read(ibf,jbf,ispin)
       enddo
     enddo
   enddo
   close(12)

 else
   stop'file not found: manual_potential'
 endif


 title='=== Read potential contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_read)

 eread = 0.5_dp*SUM(pot_read(:,:,:)*p_matrix(:,:,:))

end subroutine read_potential

!=========================================================================
subroutine setup_density_matrix(nbf,nspin,c_matrix,occupation,p_matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(in)  :: occupation(nbf,nspin)
 real(dp),intent(out) :: p_matrix(nbf,nbf,nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       p_matrix(ibf,jbf,ispin) = SUM( occupation(:,ispin) * c_matrix(ibf,:,ispin) * c_matrix(jbf,:,ispin) )
     enddo
   enddo
 enddo


end subroutine setup_density_matrix

!=========================================================================
subroutine test_density_matrix(nbf,nspin,p_matrix,s_matrix)
 use m_definitions
 use m_warning
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin),s_matrix(nbf,nbf)
!=====
 integer              :: ispin,ibf,jbf
 real(dp)             :: matrix(nbf,nbf)
 character(len=100)   :: title
!=====

 WRITE_MASTER(*,*) 'Check equality PSP = P'
 WRITE_MASTER(*,*) ' valid only for integer occupation numbers'
 do ispin=1,nspin

   !
   ! Calculate PSP
   matrix(:,:) = MATMUL( p_matrix(:,:,ispin), MATMUL( s_matrix(:,:) , p_matrix(:,:,ispin) ) )


   title='=== PSP ==='
   call dump_out_matrix(1,title,nbf,1,matrix)
   title='===  P  ==='
   call dump_out_matrix(1,title,nbf,1,p_matrix(:,:,ispin))

 enddo


end subroutine test_density_matrix

!=========================================================================
subroutine read_density_matrix(nbf,nspin,p_matrix)
 use m_definitions
 use m_warning
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(out) :: p_matrix(nbf,nbf,nspin)
!=====
 logical              :: file_exists
 integer              :: ispin,ibf,jbf
!=====

 inquire(file='manual_densitymatrix',exist=file_exists)

 if(file_exists) then
   msg='reading input density matrix'
   call issue_warning(msg)

   open(11,file='manual_densitymatrix',status='old')
   do ispin=1,nspin
     do jbf=1,nbf
       do ibf=1,nbf
         read(11,*) p_matrix(ibf,jbf,ispin) 
       enddo
     enddo
   enddo
   close(11)

 endif


end subroutine read_density_matrix

!=========================================================================
subroutine write_density_matrix(nspin,nbf,p_matrix)
 use m_definitions
 use m_warning
 use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in) :: p_matrix(nbf,nbf,nspin)
!=====
 integer              :: ispin,ibf,jbf
!=====


 WRITE_MASTER(*,*) 'output final density matrix on file'
 WRITE_MASTER(*,'(a,i5,a,i5,a,i2,/)') ' dimensions',nbf,' x ',nbf,' x ',nspin
 open(11,file='output_densitymatrix')
 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       WRITE_MASTER(11,*) p_matrix(ibf,jbf,ispin) 
     enddo
   enddo
 enddo
 close(11)


end subroutine write_density_matrix


!=========================================================================
subroutine set_occupation(electrons,magnetization,nbf,nspin,occupation)
 use m_definitions
 use m_mpi
 use m_warning
 implicit none
 real(dp),intent(in)  :: electrons,magnetization
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(out) :: occupation(nbf,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin),spin_fact
 integer              :: ibf,nlines,ilines
 logical              :: file_exists
!=====


  occupation(:,:)=0.0_dp
  spin_fact = REAL(-nspin+3,dp)

  inquire(file='manual_occupations',exist=file_exists)

  if(.NOT. file_exists) then
    remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
    if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)

    do ibf=1,nbf
      occupation(ibf,:) = MIN(remaining_electrons(:), spin_fact)
      remaining_electrons(:)  = remaining_electrons(:) - occupation(ibf,:)
    end do
  else
    WRITE_MASTER(*,*)
    WRITE_MASTER(*,*) 'occupations are read from file: manual_occupations'
    msg='reading manual occupations from file'
    call issue_warning(msg)
    open(unit=12,file='manual_occupations',status='old')
    !
    ! read nlines, all other occupations are set to zero
    read(12,*) nlines
    do ilines=1,nlines
      read(12,*) occupation(ilines,:)
    enddo
    close(12)
    WRITE_MASTER(*,*) 'occupations set, closing file'
  endif
 
  !
  ! final check
  if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0d-7 ) then
    WRITE_MASTER(*,*) 'occupation set up failed to give the right number of electrons'
    WRITE_MASTER(*,*) 'sum of occupations',SUM(occupation(:,:))
    WRITE_MASTER(*,*) 'electrons',electrons
    do ibf=1,nbf
      WRITE_MASTER(*,*) ibf,occupation(ibf,:)
    enddo
    stop'FAILURE in set_occupation'
  endif 

end subroutine set_occupation

!=========================================================================
subroutine guess_starting_c_matrix(nbf,nspin,c_matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(out) :: c_matrix(nbf,nbf,nspin)
!=====
 integer :: ibf
!=====

 !
 ! fill the c_matrix with the identity
 c_matrix(:,:,:)=0.0_dp
 do ibf=1,nbf
   c_matrix(ibf,ibf,:) = 1.0_dp
!   c_matrix(ibf,modulo(ibf,nbf)+1,:) = 1.0_dp
 enddo

end subroutine guess_starting_c_matrix

!=========================================================================
subroutine guess_starting_c_matrix_new(basis,nspin,c_matrix)
 use m_definitions
 use m_mpi
 use m_gaussian
 use m_basis_set
 implicit none
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nspin
 real(dp),intent(out)       :: c_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer  :: ibf,jbf,kbf,ig
 real(dp) :: alpha_max_bf(basis%nbf),alpha_max_remaining
!=====

 !
 ! find the sharpest gaussians
 alpha_max_bf(:)=0.0_dp
 do ibf=1,basis%nbf
   do ig=1,basis%bf(ibf)%ngaussian
     alpha_max_bf(ibf)=MAX(basis%bf(ibf)%g(ig)%alpha,alpha_max_bf(ibf))
   enddo
!   WRITE_MASTER(*,*) ibf,alpha_max_bf(ibf)
 enddo

 !
 ! fill the c_matrix 
 c_matrix(:,:,:)=0.0_dp
 do ibf=1,basis%nbf
   alpha_max_remaining=0.0_dp
   do jbf=1,basis%nbf
     if( alpha_max_bf(jbf) > alpha_max_remaining ) then
       alpha_max_remaining = alpha_max_bf(jbf)
       kbf = jbf
     endif
   enddo
   c_matrix(kbf,ibf,:) = 1.0_dp
!   WRITE_MASTER(*,*) 'chosen',ibf,kbf,alpha_max_bf(kbf)
   alpha_max_bf(kbf)   = -1.0_dp
   
 enddo

end subroutine guess_starting_c_matrix_new

!=========================================================================
subroutine setup_initial_c_matrix(print_volume,nbf,nspin,hamiltonian_nucleus,s_matrix,occupation,c_matrix)
 use m_definitions
 use m_mpi
 use m_tools
 implicit none
 integer,intent(in)         :: print_volume,nspin,nbf
 real(dp),intent(in)        :: hamiltonian_nucleus(nbf,nbf),s_matrix(nbf,nbf)
 real(dp),intent(in)        :: occupation(nbf,nspin)
 real(dp),intent(out)       :: c_matrix(nbf,nbf,nspin)
!=====
 integer                    :: ibf,jbf,kbf,lbf
 real(dp)                   :: hamiltonian(nbf,nbf),matrix(nbf,nbf),energy(nbf)
 real(dp)                   :: coeff_max,bonding
 real(dp)                   :: line(nbf)
 character(len=100)         :: title
!=====


 !
 ! Diagonalize a spin independant hamiltonian
 ! to obtain a starting point for matrix C
! hamiltonian(:,:) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
 hamiltonian(:,:) = hamiltonian_nucleus(:,:)

 title='=== bare hamiltonian ==='
 call dump_out_matrix(print_volume,title,nbf,1,hamiltonian)

 WRITE_MASTER(*,*) 'Diagonalization of an initial hamiltonian'
 call diagonalize_generalized_sym(nbf,hamiltonian,s_matrix,energy,matrix)

 title='=== Energies ==='
 call dump_out_eigenenergy(title,nbf,1,occupation,energy)


 ibf=1
 do while( ibf < nbf )
!   if( ALL( occupation(ibf,:) < completely_empty ) cycle

   jbf = ibf + 1
   !
   ! Find degenerate energies
   if( ABS( energy(ibf) - energy(jbf) ) < 1.0e-3_dp ) then
     !
     ! Find the bonding and anti bonding states
     coeff_max = MAXVAL( ABS( matrix(:,ibf) ) )
     bonding = 1.0_dp
     do kbf=1,nbf
       if( ABS( matrix(kbf,ibf) ) > coeff_max - 1.0e-3_dp ) then
         bonding = SIGN( 1.0_dp, matrix(kbf,ibf) ) * bonding
       endif
     enddo
     if( bonding < 0.0_dp ) then
       line(:)       = matrix(:,ibf)
       matrix(:,ibf) = matrix(:,jbf)
       matrix(:,jbf) = line(:)
     endif

     ibf = ibf + 2
   else
     ibf = ibf + 1
   endif

 enddo

 c_matrix(:,:,1)     = matrix(:,:)
 c_matrix(:,:,nspin) = matrix(:,:)

 matrix(:,:) = transpose( matrix(:,:) )
 title='=== C matrix ==='
 call dump_out_matrix(print_volume,title,nbf,1,matrix)


end subroutine setup_initial_c_matrix


!=========================================================================
subroutine matrix_basis_to_eigen(nspin,nbf,c_matrix,matrix_inout)
 use m_definitions
 implicit none
 integer,intent(in)      :: nspin,nbf
 real(dp),intent(in)     :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(inout)  :: matrix_inout(nbf,nbf,nspin)
!====
 integer                 :: ispin
!====


 do ispin=1,nspin
   matrix_inout(:,:,ispin) = MATMUL( TRANSPOSE( c_matrix(:,:,ispin) ) , MATMUL( matrix_inout(:,:,ispin) , c_matrix(:,:,ispin) ) )
 enddo


end subroutine matrix_basis_to_eigen
!=========================================================================
