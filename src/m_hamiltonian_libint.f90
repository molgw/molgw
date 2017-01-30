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
 use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min



 interface

   subroutine libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             overlapAB) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value  :: amA,contrdepthA
     real(C_DOUBLE),intent(in) :: A(*)
     real(C_DOUBLE),intent(in) :: alphaA(*)
     real(C_DOUBLE),intent(in) :: cA(*)
     integer(C_INT),value  :: amB,contrdepthB
     real(C_DOUBLE),intent(in) :: B(*)
     real(C_DOUBLE),intent(in) :: alphaB(*)
     real(C_DOUBLE),intent(in) :: cB(*)
     real(C_DOUBLE),intent(out) :: overlapAB(*)
     
   end subroutine libint_overlap

   subroutine libint_overlap_grad(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             overlapABx,overlapABy,overlapABz) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value  :: amA,contrdepthA
     real(C_DOUBLE),intent(in) :: A(*)
     real(C_DOUBLE),intent(in) :: alphaA(*)
     real(C_DOUBLE),intent(in) :: cA(*)
     integer(C_INT),value  :: amB,contrdepthB
     real(C_DOUBLE),intent(in) :: B(*)
     real(C_DOUBLE),intent(in) :: alphaB(*)
     real(C_DOUBLE),intent(in) :: cB(*)
     real(C_DOUBLE),intent(out) :: overlapABx(*)
     real(C_DOUBLE),intent(out) :: overlapABy(*)
     real(C_DOUBLE),intent(out) :: overlapABz(*)
     
   end subroutine libint_overlap_grad

   subroutine libint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             kineticAB) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value       :: amA,contrdepthA
     real(C_DOUBLE),intent(in)  :: A(*)
     real(C_DOUBLE),intent(in)  :: alphaA(*)
     real(C_DOUBLE),intent(in)  :: cA(*)
     integer(C_INT),value       :: amB,contrdepthB
     real(C_DOUBLE),intent(in)  :: B(*)
     real(C_DOUBLE),intent(in)  :: alphaB(*)
     real(C_DOUBLE),intent(in)  :: cB(*)
     real(C_DOUBLE),intent(out) :: kineticAB(*)

   end subroutine libint_kinetic

   subroutine libint_kinetic_grad(amA,contrdepthA,A,alphaA,cA, &
                                  amB,contrdepthB,B,alphaB,cB, &
                                  kineticABx,kineticABy,kineticABz) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value       :: amA,contrdepthA
     real(C_DOUBLE),intent(in)  :: A(*)
     real(C_DOUBLE),intent(in)  :: alphaA(*)
     real(C_DOUBLE),intent(in)  :: cA(*)
     integer(C_INT),value       :: amB,contrdepthB
     real(C_DOUBLE),intent(in)  :: B(*)
     real(C_DOUBLE),intent(in)  :: alphaB(*)
     real(C_DOUBLE),intent(in)  :: cB(*)
     real(C_DOUBLE),intent(out) :: kineticABx(*)
     real(C_DOUBLE),intent(out) :: kineticABy(*)
     real(C_DOUBLE),intent(out) :: kineticABz(*)

   end subroutine libint_kinetic_grad

   subroutine libint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             C,elecpotAB) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value         :: amA,contrdepthA
     real(C_DOUBLE),intent(in)    :: A(*)
     real(C_DOUBLE),intent(in)    :: alphaA(*)
     real(C_DOUBLE),intent(in)    :: cA(*)
     integer(C_INT),value         :: amB,contrdepthB
     real(C_DOUBLE),intent(in)    :: B(*)
     real(C_DOUBLE),intent(in)    :: alphaB(*)
     real(C_DOUBLE),intent(in)    :: cB(*)
     real(C_DOUBLE),intent(in)    :: C(*)
     real(C_DOUBLE),intent(inout) :: elecpotAB(*)
     
   end subroutine libint_elecpot

   subroutine libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                  amB,contrdepthB,B,alphaB,cB, &
                                  C,elecpotAx,elecpotAy,elecpotAz, &
                                    elecpotBx,elecpotBy,elecpotBz) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value         :: amA,contrdepthA
     real(C_DOUBLE),intent(in)    :: A(*)
     real(C_DOUBLE),intent(in)    :: alphaA(*)
     real(C_DOUBLE),intent(in)    :: cA(*)
     integer(C_INT),value         :: amB,contrdepthB
     real(C_DOUBLE),intent(in)    :: B(*)
     real(C_DOUBLE),intent(in)    :: alphaB(*)
     real(C_DOUBLE),intent(in)    :: cB(*)
     real(C_DOUBLE),intent(in)    :: C(*)
     real(C_DOUBLE),intent(inout) :: elecpotAx(*)
     real(C_DOUBLE),intent(inout) :: elecpotAy(*)
     real(C_DOUBLE),intent(inout) :: elecpotAz(*)
     real(C_DOUBLE),intent(inout) :: elecpotBx(*)
     real(C_DOUBLE),intent(inout) :: elecpotBy(*)
     real(C_DOUBLE),intent(inout) :: elecpotBz(*)
     
   end subroutine libint_elecpot_grad

 end interface



contains



!=========================================================================
subroutine setup_overlap_libint(print_matrix_,basis,s_matrix)
 use m_basis_set
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

 real(C_DOUBLE),allocatable        :: matrix_cart(:,:)
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
 write(stdout,'(/,a)') ' Setup overlap matrix S'

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

     allocate(matrix_cart(ni_cart,nj_cart))

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
                         amB,contrdepthB,B,alphaB,cB,matrix_cart)
#endif

     deallocate(alphaA,alphaB,cA,cB)


     s_matrix(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                   MATMUL( matrix_cart(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )


     deallocate(matrix_cart)
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
 use m_basis_set
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

 real(C_DOUBLE),allocatable        :: matrix_cart_gradx(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_grady(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_gradz(:,:)
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
 write(stdout,'(/,a)') ' Setup overlap matrix S'

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

     allocate(matrix_cart_gradx(ni_cart,nj_cart))
     allocate(matrix_cart_grady(ni_cart,nj_cart))
     allocate(matrix_cart_gradz(ni_cart,nj_cart))

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
                              matrix_cart_gradx,matrix_cart_grady,matrix_cart_gradz)
#endif

     deallocate(alphaA,alphaB,cA,cB)


     s_matrix_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,1) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                        MATMUL( matrix_cart_gradx(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )
     s_matrix_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,2) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                        MATMUL( matrix_cart_grady(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )
     s_matrix_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,3) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                   MATMUL( matrix_cart_gradz(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )


     deallocate(matrix_cart_gradx)
     deallocate(matrix_cart_grady)
     deallocate(matrix_cart_gradz)

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
 use m_basis_set
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

 real(C_DOUBLE),allocatable        :: matrix_cart(:,:)
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
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian'

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

     allocate(matrix_cart(ni_cart,nj_cart))


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
                         amB,contrdepthB,B,alphaB,cB,matrix_cart)
#endif

     deallocate(alphaA,alphaB,cA,cB)


     hamiltonian_kinetic(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )


     deallocate(matrix_cart)
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
 use m_basis_set
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

 real(C_DOUBLE),allocatable        :: matrix_cart_gradx(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_grady(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_gradz(:,:)
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
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian'

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

     allocate(matrix_cart_gradx(ni_cart,nj_cart))
     allocate(matrix_cart_grady(ni_cart,nj_cart))
     allocate(matrix_cart_gradz(ni_cart,nj_cart))


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
                              matrix_cart_gradx,matrix_cart_grady,matrix_cart_gradz)
#endif

     deallocate(alphaA,alphaB,cA,cB)


     hamiltonian_kinetic_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,1) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart_gradx(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )
     hamiltonian_kinetic_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,2) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart_grady(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )
     hamiltonian_kinetic_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,3) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart_gradz(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )


     deallocate(matrix_cart_gradx)
     deallocate(matrix_cart_grady)
     deallocate(matrix_cart_gradz)

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
 use m_basis_set
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
 real(dp)             :: vnucleus_ij

 real(C_DOUBLE),allocatable        :: matrix_cart(:,:)
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

     allocate(matrix_cart(ni_cart,nj_cart))
     matrix_cart(:,:) = 0.0_dp

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
                           C,matrix_cart)
#endif

     enddo
     deallocate(alphaA,alphaB,cA,cB)


     hamiltonian_nucleus(ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                              MATMUL( matrix_cart(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) ) 


     deallocate(matrix_cart)
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
 use m_basis_set
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus_grad(basis%nbf,basis%nbf,natom,3)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp)             :: vnucleus_ij

 real(C_DOUBLE),allocatable        :: matrix_cart_gradAx(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_gradAy(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_gradAz(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_gradBx(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_gradBy(:,:)
 real(C_DOUBLE),allocatable        :: matrix_cart_gradBz(:,:)
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

     allocate(matrix_cart_gradAx(ni_cart,nj_cart))
     allocate(matrix_cart_gradAy(ni_cart,nj_cart))
     allocate(matrix_cart_gradAz(ni_cart,nj_cart))
     allocate(matrix_cart_gradBx(ni_cart,nj_cart))
     allocate(matrix_cart_gradBy(ni_cart,nj_cart))
     allocate(matrix_cart_gradBz(ni_cart,nj_cart))

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

       matrix_cart_gradAx(:,:) = 0.0_dp
       matrix_cart_gradAy(:,:) = 0.0_dp
       matrix_cart_gradAz(:,:) = 0.0_dp
       matrix_cart_gradBx(:,:) = 0.0_dp
       matrix_cart_gradBy(:,:) = 0.0_dp
       matrix_cart_gradBz(:,:) = 0.0_dp

       C(:) = x(:,iatom)
#ifdef HAVE_LIBINT_ONEBODY
       call libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                amB,contrdepthB,B,alphaB,cB, &
                                C,matrix_cart_gradAx,matrix_cart_gradAy,matrix_cart_gradAz, &
                                  matrix_cart_gradBx,matrix_cart_gradBy,matrix_cart_gradBz)
#endif

       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,iatom,1) = -MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                  MATMUL( matrix_cart_gradAx(:,:) , cart_to_pure_norm(lj)%matrix(:,:)) )   &
                                                - MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) ,                   &
                                                  MATMUL( matrix_cart_gradBx(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) ) 
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,iatom,2) = -MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                  MATMUL( matrix_cart_gradAy(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )  &
                                                - MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) ,                   &  
                                                  MATMUL( matrix_cart_gradBy(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) ) 
       hamiltonian_nucleus_grad(ibf:ibf+ni-1,jbf:jbf+nj-1,iatom,3) = -MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) , &
                                                  MATMUL( matrix_cart_gradAz(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) )  &
                                                - MATMUL( TRANSPOSE(cart_to_pure_norm(li)%matrix(:,:)) ,                   &
                                                  MATMUL( matrix_cart_gradBz(:,:) , cart_to_pure_norm(lj)%matrix(:,:) ) ) 
     enddo
     deallocate(alphaA,alphaB,cA,cB)


     deallocate(matrix_cart_gradAx)
     deallocate(matrix_cart_gradAy)
     deallocate(matrix_cart_gradAz)
     deallocate(matrix_cart_gradBx)
     deallocate(matrix_cart_gradBy)
     deallocate(matrix_cart_gradBz)

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
