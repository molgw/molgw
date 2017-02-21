!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian_libint
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min
 use m_basis_set
 use m_libint_tools





contains


!=========================================================================
subroutine setup_overlap_libint(print_matrix_,basis,s_matrix)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(basis%nbf,basis%nbf)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup overlap matrix S (LIBINT)'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(array_cart(ni_cart*nj_cart))

     amA = li
     amB = lj
     contrdepthA = basis%bf(ibf_cart)%ngaussian
     contrdepthB = basis%bf(jbf_cart)%ngaussian
     A(:) = basis%bf(ibf_cart)%x0(:)
     B(:) = basis%bf(jbf_cart)%x0(:)
     allocate(alphaA(contrdepthA),alphaB(contrdepthB))
     alphaA(:) = basis%bf(ibf_cart)%g(:)%alpha
     alphaB(:) = basis%bf(jbf_cart)%g(:)%alpha
     allocate(cA(contrdepthA),cB(contrdepthB))
     cA(:) = basis%bf(ibf_cart)%coeff(:) * basis%bf(ibf_cart)%g(:)%common_norm_factor
     cB(:) = basis%bf(jbf_cart)%coeff(:) * basis%bf(jbf_cart)%g(:)%common_norm_factor
     
#ifdef HAVE_LIBINT_ONEBODY
     call libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB,array_cart)
#endif

     deallocate(alphaA,alphaB,cA,cB)

     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)

     s_matrix(ibf:ibf+ni-1,jbf:jbf+nj-1) = matrix(:,:)

     deallocate(array_cart,matrix)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='=== Overlap matrix S (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix)

 call stop_clock(timing_overlap)


end subroutine setup_overlap_libint


!=========================================================================
subroutine setup_overlap_grad_libint(print_matrix_,basis,s_matrix_grad)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix_grad(basis%nbf,basis%nbf,3)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradx(:)
 real(C_DOUBLE),allocatable        :: array_cart_grady(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup gradient of the overlap matrix S (LIBINT)'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(array_cart_gradx(ni_cart*nj_cart))
     allocate(array_cart_grady(ni_cart*nj_cart))
     allocate(array_cart_gradz(ni_cart*nj_cart))

     amA = li
     amB = lj
     contrdepthA = basis%bf(ibf_cart)%ngaussian
     contrdepthB = basis%bf(jbf_cart)%ngaussian
     A(:) = basis%bf(ibf_cart)%x0(:)
     B(:) = basis%bf(jbf_cart)%x0(:)
     allocate(alphaA(contrdepthA),alphaB(contrdepthB))
     alphaA(:) = basis%bf(ibf_cart)%g(:)%alpha
     alphaB(:) = basis%bf(jbf_cart)%g(:)%alpha
     allocate(cA(contrdepthA),cB(contrdepthB))
     cA(:) = basis%bf(ibf_cart)%coeff(:) * basis%bf(ibf_cart)%g(:)%common_norm_factor
     cB(:) = basis%bf(jbf_cart)%coeff(:) * basis%bf(jbf_cart)%g(:)%common_norm_factor
     
#ifdef HAVE_LIBINT_ONEBODY
     call libint_overlap_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              array_cart_gradx,array_cart_grady,array_cart_gradz)
#else
     call die('overlap gradient not implemented without LIBINT one-body terms')
#endif

     deallocate(alphaA,alphaB,cA,cB)

     ! X
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradx,matrix)
     s_matrix_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,1) = matrix(:,:)

     ! Y
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_grady,matrix)
     s_matrix_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,2) = matrix(:,:)

     ! Z
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradz,matrix)
     s_matrix_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,3) = matrix(:,:)


     deallocate(array_cart_gradx)
     deallocate(array_cart_grady)
     deallocate(array_cart_gradz)
     deallocate(matrix)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='=== Overlap matrix S (LIBINT) X ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix_grad(:,:,1))
 title='=== Overlap matrix S (LIBINT) Y ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix_grad(:,:,2))
 title='=== Overlap matrix S (LIBINT) Z ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix_grad(:,:,3))

 call stop_clock(timing_overlap)


end subroutine setup_overlap_grad_libint


!=========================================================================
subroutine setup_kinetic_libint(print_matrix_,basis,hamiltonian_kinetic)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(basis%nbf,basis%nbf)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian (LIBINT)'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(array_cart(ni_cart*nj_cart))


     amA = li
     amB = lj
     contrdepthA = basis%bf(ibf_cart)%ngaussian
     contrdepthB = basis%bf(jbf_cart)%ngaussian
     A(:) = basis%bf(ibf_cart)%x0(:)
     B(:) = basis%bf(jbf_cart)%x0(:)
     allocate(alphaA(contrdepthA),alphaB(contrdepthB))
     alphaA(:) = basis%bf(ibf_cart)%g(:)%alpha
     alphaB(:) = basis%bf(jbf_cart)%g(:)%alpha
     allocate(cA(contrdepthA),cB(contrdepthB))
     cA(:) = basis%bf(ibf_cart)%coeff(:) * basis%bf(ibf_cart)%g(:)%common_norm_factor
     cB(:) = basis%bf(jbf_cart)%coeff(:) * basis%bf(jbf_cart)%g(:)%common_norm_factor

#ifdef HAVE_LIBINT_ONEBODY
     call libint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB,array_cart)
#endif

     deallocate(alphaA,alphaB,cA,cB)

     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)

     hamiltonian_kinetic(ibf:ibf+ni-1,jbf:jbf+nj-1) = matrix(:,:)


     deallocate(array_cart,matrix)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='===  Kinetic energy contribution (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic)

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic_libint


!=========================================================================
subroutine setup_kinetic_grad_libint(print_matrix_,basis,hamiltonian_kinetic_grad)
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic_grad(basis%nbf,basis%nbf,3)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradx(:)
 real(C_DOUBLE),allocatable        :: array_cart_grady(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup gradient of the kinetic part of the Hamiltonian (LIBINT)'

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(array_cart_gradx(ni_cart*nj_cart))
     allocate(array_cart_grady(ni_cart*nj_cart))
     allocate(array_cart_gradz(ni_cart*nj_cart))


     amA = li
     amB = lj
     contrdepthA = basis%bf(ibf_cart)%ngaussian
     contrdepthB = basis%bf(jbf_cart)%ngaussian
     A(:) = basis%bf(ibf_cart)%x0(:)
     B(:) = basis%bf(jbf_cart)%x0(:)
     allocate(alphaA(contrdepthA),alphaB(contrdepthB))
     alphaA(:) = basis%bf(ibf_cart)%g(:)%alpha
     alphaB(:) = basis%bf(jbf_cart)%g(:)%alpha
     allocate(cA(contrdepthA),cB(contrdepthB))
     cA(:) = basis%bf(ibf_cart)%coeff(:) * basis%bf(ibf_cart)%g(:)%common_norm_factor
     cB(:) = basis%bf(jbf_cart)%coeff(:) * basis%bf(jbf_cart)%g(:)%common_norm_factor

#ifdef HAVE_LIBINT_ONEBODY
     call libint_kinetic_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              array_cart_gradx,array_cart_grady,array_cart_gradz)
#else
     call die('kinetic operator gradient not implemented without LIBINT one-body terms')
#endif

     deallocate(alphaA,alphaB,cA,cB)

     ! X
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradx,matrix)
     hamiltonian_kinetic_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,1) = matrix(:,:)

     ! Y
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_grady,matrix)
     hamiltonian_kinetic_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,2) = matrix(:,:)

     ! Z
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradz,matrix)
     hamiltonian_kinetic_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,3) = matrix(:,:)


     deallocate(array_cart_gradx)
     deallocate(array_cart_grady)
     deallocate(array_cart_gradz)
     deallocate(matrix)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='===  Kinetic energy contribution (LIBINT) X ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,1))
 title='===  Kinetic energy contribution (LIBINT) Y ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,2))
 title='===  Kinetic energy contribution (LIBINT) Z ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,3))

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic_grad_libint


!=========================================================================
subroutine setup_nucleus_libint(print_matrix_,basis,hamiltonian_nucleus)
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: vnucleus_ij

 real(C_DOUBLE),allocatable        :: array_cart(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
 real(C_DOUBLE)                    :: C(3)
!=====

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian (LIBINT)'
 if( nproc_world > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(array_cart(ni_cart*nj_cart))
     array_cart(:) = 0.0_dp

     amA = li
     amB = lj
     contrdepthA = basis%bf(ibf_cart)%ngaussian
     contrdepthB = basis%bf(jbf_cart)%ngaussian
     A(:) = basis%bf(ibf_cart)%x0(:)
     B(:) = basis%bf(jbf_cart)%x0(:)
     allocate(alphaA(contrdepthA),alphaB(contrdepthB))
     alphaA(:) = basis%bf(ibf_cart)%g(:)%alpha
     alphaB(:) = basis%bf(jbf_cart)%g(:)%alpha
     allocate(cA(contrdepthA),cB(contrdepthB))
     cA(:) = basis%bf(ibf_cart)%coeff(:) * basis%bf(ibf_cart)%g(:)%common_norm_factor

     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle

       cB(:) = basis%bf(jbf_cart)%coeff(:) * basis%bf(jbf_cart)%g(:)%common_norm_factor * (-zatom(iatom))

       C(:) = x(:,iatom)
#ifdef HAVE_LIBINT_ONEBODY
       call libint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                           amB,contrdepthB,B,alphaB,cB, &
                           C,array_cart)
#endif

     enddo
     deallocate(alphaA,alphaB,cA,cB)

     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)

     hamiltonian_nucleus(ibf:ibf+ni-1,jbf:jbf+nj-1) = matrix(:,:)


     deallocate(array_cart,matrix)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 !
 ! Reduce operation
 call xsum_world(hamiltonian_nucleus)

 title='===  Nucleus potential contribution (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus_libint


!=========================================================================
subroutine setup_nucleus_grad_libint(print_matrix_,basis,hamiltonian_nucleus_grad)
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus_grad(basis%nbf,basis%nbf,natom+1,3)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp)             :: vnucleus_ij
 real(dp),allocatable :: matrixA(:,:)
 real(dp),allocatable :: matrixB(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradAx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAz(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
 real(C_DOUBLE)                    :: C(3)
!=====

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian gradient (LIBINT)'
 if( nproc_world > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif

 hamiltonian_nucleus_grad(:,:,:,:) = 0.0_dp

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(array_cart_gradAx(ni_cart*nj_cart))
     allocate(array_cart_gradAy(ni_cart*nj_cart))
     allocate(array_cart_gradAz(ni_cart*nj_cart))
     allocate(array_cart_gradBx(ni_cart*nj_cart))
     allocate(array_cart_gradBy(ni_cart*nj_cart))
     allocate(array_cart_gradBz(ni_cart*nj_cart))

     amA = li
     amB = lj
     contrdepthA = basis%bf(ibf_cart)%ngaussian
     contrdepthB = basis%bf(jbf_cart)%ngaussian
     A(:) = basis%bf(ibf_cart)%x0(:)
     B(:) = basis%bf(jbf_cart)%x0(:)
     allocate(alphaA(contrdepthA),alphaB(contrdepthB))
     alphaA(:) = basis%bf(ibf_cart)%g(:)%alpha
     alphaB(:) = basis%bf(jbf_cart)%g(:)%alpha
     allocate(cA(contrdepthA),cB(contrdepthB))
     cA(:) = basis%bf(ibf_cart)%coeff(:) * basis%bf(ibf_cart)%g(:)%common_norm_factor

     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle

       cB(:) = basis%bf(jbf_cart)%coeff(:) * basis%bf(jbf_cart)%g(:)%common_norm_factor * (-zatom(iatom))

       array_cart_gradAx(:) = 0.0_dp
       array_cart_gradAy(:) = 0.0_dp
       array_cart_gradAz(:) = 0.0_dp
       array_cart_gradBx(:) = 0.0_dp
       array_cart_gradBy(:) = 0.0_dp
       array_cart_gradBz(:) = 0.0_dp

       C(:) = x(:,iatom)
#ifdef HAVE_LIBINT_ONEBODY
       call libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                amB,contrdepthB,B,alphaB,cB, &
                                C,                           &
                                array_cart_gradAx,array_cart_gradAy,array_cart_gradAz, &
                                array_cart_gradBx,array_cart_gradBy,array_cart_gradBz)
#else
       call die('nuclear potential gradient not implemented without LIBINT one-body terms')
#endif

       ! X
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAx,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBx,matrixB)
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,iatom  ,1) = matrixA(:,:) + matrixB(:,:)
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,natom+1,1) = hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,natom+1,1) + matrixA(:,:)

       ! Y
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAy,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBy,matrixB)
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,iatom  ,2) = matrixA(:,:) + matrixB(:,:)
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,natom+1,2) = hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,natom+1,2) + matrixA(:,:)

       ! Z
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAz,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBz,matrixB)
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,iatom  ,3) = matrixA(:,:) + matrixB(:,:)
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,natom+1,3) = hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,natom+1,3) + matrixA(:,:)

     enddo
     deallocate(alphaA,alphaB,cA,cB)


     deallocate(array_cart_gradAx)
     deallocate(array_cart_gradAy)
     deallocate(array_cart_gradAz)
     deallocate(array_cart_gradBx)
     deallocate(array_cart_gradBy)
     deallocate(array_cart_gradBz)
     deallocate(matrixA)
     deallocate(matrixB)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 !
 ! Reduce operation
 call xsum_world(hamiltonian_nucleus_grad)

 title='===  Nucleus potential contribution (LIBINT) C1X==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,1))
 title='===  Nucleus potential contribution (LIBINT) C1Y==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,2))
 title='===  Nucleus potential contribution (LIBINT) C1Z==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,3))

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus_grad_libint


end module m_hamiltonian_libint
!=========================================================================
