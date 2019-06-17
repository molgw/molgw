!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian with full distribution
! of the memory with SCALAPACK
!
! !!! EXPERIMENTAL!!!
!
!=========================================================================
module m_hamiltonian_sca
 use m_definitions
 use m_mpi
 use m_scalapack
 use m_timing
 use m_tddft_variables
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min
 use m_hamiltonian_tools
 use m_density_tools
 use m_linear_algebra,only: diagonalize


contains


!=========================================================================
subroutine matrix_cart_to_local(gaussian_type,ibf,jbf,li,lj,ni_cart,nj_cart, &
                                matrix_cart,ni,nj,matrix_local)
 use m_basis_set
 implicit none

 character(len=4),intent(in) :: gaussian_type
 integer,intent(in)          :: ibf,jbf
 integer,intent(in)          :: li,lj
 integer,intent(in)          :: ni_cart,nj_cart,ni,nj
 real(dp),intent(in)         :: matrix_cart(ni_cart,nj_cart)
 real(dp),intent(inout)      :: matrix_local(:,:)
!=====
 integer  :: gt
 integer  :: iglobal,jglobal
 integer  :: ilocal,jlocal
 real(dp) :: matrix_final(ibf:ibf+ni-1,jbf:jbf+nj-1)
!=====

#if defined(HAVE_SCALAPACK)

 gt = get_gaussian_type_tag(gaussian_type)

 matrix_final(:,:) = MATMUL( TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)) , &
                             MATMUL( matrix_cart(:,:) , cart_to_pure(lj,gt)%matrix(:,:) ) )

 do jglobal=jbf,jbf+nj-1
   jlocal = colindex_global_to_local('H',jglobal)
   if( jlocal == 0 ) cycle

   do iglobal=ibf,ibf+ni-1
     ilocal = rowindex_global_to_local('H',iglobal)
     if( ilocal == 0 ) cycle

     matrix_local(ilocal,jlocal) = matrix_final(iglobal,jglobal)

   enddo
 enddo


#endif

end subroutine matrix_cart_to_local


!=========================================================================
subroutine setup_overlap_sca(basis,s_matrix)
 use m_basis_set
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(:,:)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 real(dp),allocatable :: matrix_cart(:,:)
!=====

#if defined(HAVE_SCALAPACK)

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup overlap matrix S: SCALAPACK'

 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj        = number_basis_function_am(basis%gaussian_type,lj)
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni        = number_basis_function_am(basis%gaussian_type,li)
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call overlap_basis_function(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo

     call matrix_cart_to_local(basis%gaussian_type,ibf1,jbf1,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,s_matrix)


     deallocate(matrix_cart)
   enddo
 enddo


 call stop_clock(timing_overlap)


#endif

end subroutine setup_overlap_sca


!=========================================================================
subroutine setup_kinetic_sca(basis,hamiltonian_kinetic)
 use m_basis_set
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(:,:)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 real(dp),allocatable :: matrix_cart(:,:)
!=====

#if defined(HAVE_SCALAPACK)

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian: SCALAPACK'

 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj        = number_basis_function_am(basis%gaussian_type,lj)
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni        = number_basis_function_am(basis%gaussian_type,li)
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call kinetic_basis_function(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo

     call matrix_cart_to_local(basis%gaussian_type,ibf1,jbf1,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,hamiltonian_kinetic)

     deallocate(matrix_cart)
   enddo
 enddo


 call stop_clock(timing_hamiltonian_kin)

#endif

end subroutine setup_kinetic_sca


!=========================================================================
subroutine setup_nucleus_sca(basis,hamiltonian_nucleus)
 use m_basis_set
 use m_atoms
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus(:,:)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: natom_local
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
!=====

#if defined(HAVE_SCALAPACK)

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian: SCALAPACK'
 if( nproc_local > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_local /= MODULO(iatom-1,nproc_local) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif


 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj        = number_basis_function_am(basis%gaussian_type,lj)
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni        = number_basis_function_am(basis%gaussian_type,li)
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(matrix_cart(ni_cart,nj_cart))
     matrix_cart(:,:) = 0.0_dp
     do iatom=1,natom
       if( rank_local /= MODULO(iatom-1,nproc_local) ) cycle
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           call nucleus_basis_function(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),zvalence(iatom),xatom(:,iatom),vnucleus_ij)
           matrix_cart(i_cart,j_cart) = matrix_cart(i_cart,j_cart) + vnucleus_ij
         enddo
       enddo
     enddo

     call matrix_cart_to_local(basis%gaussian_type,ibf1,jbf1,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,hamiltonian_nucleus)


     deallocate(matrix_cart)
   enddo
 enddo

 !
 ! Reduce operation
 call xsum_local(hamiltonian_nucleus)


 call stop_clock(timing_hamiltonian_nuc)

#endif

end subroutine setup_nucleus_sca


!=========================================================================
subroutine setup_density_matrix_sca(c_matrix,occupation,p_matrix)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(out) :: p_matrix(:,:,:)
!=====
 integer  :: nbf,nstate
 integer  :: ispin,istate
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

#if defined(HAVE_SCALAPACK)
 call start_clock(timing_density_matrix)
 write(stdout,'(1x,a)') 'Build density matrix: SCALAPACK'

 nbf    = desc_c(M_)
 nstate = desc_c(N_)

 p_matrix(:,:,:) = 0.0_dp

 if( cntxt_ham > 0 ) then
   do ispin=1,nspin
     do istate=1,nstate
       if( occupation(istate,ispin) < completely_empty ) cycle
       call PDSYR('L',nbf,occupation(istate,ispin),c_matrix(:,:,ispin),1,istate,desc_c,1,p_matrix(:,:,ispin),1,1,desc_ham)
     enddo
     allocate(matrix_tmp,MOLD=p_matrix(:,:,1))
     call symmetrize_matrix_sca('L',nbf,desc_ham,p_matrix(:,:,ispin),desc_ham,matrix_tmp)
     deallocate(matrix_tmp)
   enddo
 endif

 ! Poor man distribution
 call xsum_local(p_matrix)


 call stop_clock(timing_density_matrix)
#else
 call die('setup_density_matrix_sca: is being called without SCALAPACK')
#endif

end subroutine setup_density_matrix_sca


!=========================================================================
subroutine diagonalize_hamiltonian_scalapack(hamiltonian,x_matrix,energy,c_matrix)
 implicit none

 real(dp),intent(in)  :: hamiltonian(:,:,:)
 real(dp),intent(in)  :: x_matrix(:,:)
 real(dp),intent(out) :: c_matrix(:,:,:)
 real(dp),intent(out) :: energy(:,:)
!=====
 integer :: nspin_local,nbf,nstate
 integer :: mh,nh,mc,nc,ms,ns
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
#if defined(HAVE_SCALAPACK)
 integer :: cntxt
 integer :: rank_sca,nprocs_sca
 integer :: desch(NDEL),descc(NDEL),descs(NDEL)
#endif
 integer  :: ispin
 integer  :: ilocal,jlocal,iglobal,jglobal
 integer  :: m_small,n_small
 real(dp),allocatable :: h_small(:,:),h_small2(:,:)
 real(dp),allocatable :: ham_local(:,:),c_matrix_local(:,:),s_matrix_local(:,:)
!=====

 nbf         = SIZE(c_matrix,DIM=1)
 nstate      = SIZE(c_matrix,DIM=2)
 nspin_local = SIZE(c_matrix,DIM=3)

#if defined(HAVE_SCALAPACK)

 nprow = MIN(nprow_sd,nbf/scalapack_block_min)
 npcol = MIN(npcol_sd,nbf/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 if( nprow * npcol > 1 ) then
   write(stdout,'(1x,a,i4,a,i4)') 'Generalized diagonalization using SCALAPACK with a grid',nprow,' x ',npcol
   call BLACS_PINFO( rank_sca, nprocs_sca )

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

   c_matrix(:,:,:) = 0.0_dp

   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if(cntxt > 0 ) then

     mh = NUMROC(nbf   ,block_row,iprow,first_row,nprow)
     nh = NUMROC(nbf   ,block_col,ipcol,first_col,npcol)
     mc = NUMROC(nbf   ,block_row,iprow,first_row,nprow)
     nc = NUMROC(nstate,block_col,ipcol,first_col,npcol)
     ms = NUMROC(nstate,block_row,iprow,first_row,nprow)
     ns = NUMROC(nstate,block_col,ipcol,first_col,npcol)


     call DESCINIT(desch,nbf   ,nbf   ,block_row,block_col,first_row,first_col,cntxt,MAX(1,mh),info)
     call DESCINIT(descc,nbf   ,nstate,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)
     call DESCINIT(descs,nstate,nstate,block_row,block_col,first_row,first_col,cntxt,MAX(1,ms),info)


     allocate(ham_local(mh,nh))
     allocate(c_matrix_local(mc,nc))
     allocate(s_matrix_local(mc,nc))
     allocate(h_small(ms,ns))

     !
     ! Set up the local copy of x_matrix
     do jlocal=1,nc
       jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
       do ilocal=1,mc
         iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
         s_matrix_local(ilocal,jlocal) = x_matrix(iglobal,jglobal)
       enddo
     enddo



     do ispin=1,nspin_local
       write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
       call start_clock(timing_diago_hamiltonian)

       !
       ! Set up the local copy of hamiltonian
       do jlocal=1,nh
         jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
         do ilocal=1,mh
           iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
           ham_local(ilocal,jlocal) = hamiltonian(iglobal,jglobal,ispin)
         enddo
       enddo

!       h_small(:,:) = MATMUL( TRANSPOSE(x_matrix(:,:)) , &
!                                MATMUL( hamiltonian(:,:,ispin) , x_matrix(:,:) ) )

       !
       ! H_small = ^tS^{-1/2} H S^{-1/2}
       call PDGEMM('N','N',nbf,nstate,nbf,                &
                    1.0_dp,ham_local,1,1,desch,           &
                    s_matrix_local,1,1,descc,             &
                    0.0_dp,c_matrix_local,1,1,descc)

       call PDGEMM('T','N',nstate,nstate,nbf,             &
                    1.0_dp,s_matrix_local,1,1,descc,      &
                    c_matrix_local,1,1,descc,             &
                    0.0_dp,h_small,1,1,descs)



       call diagonalize_sca(scf_diago_flavor,h_small,descs,energy(:,ispin))


!       c_matrix(:,:,ispin) = MATMUL( x_matrix(:,:) , h_small(:,:) )

       !
       ! C = S^{-1/2} C_small
       call PDGEMM('N','N',nbf,nstate,nstate,             &
                    1.0_dp,s_matrix_local,1,1,descc,      &
                    h_small,1,1,descs,                    &
                    0.0_dp,c_matrix_local,1,1,descc)


       do jlocal=1,nc
         jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
         do ilocal=1,mc
           iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
           c_matrix(iglobal,jglobal,ispin) = c_matrix_local(ilocal,jlocal)
         enddo
       enddo


      ! Nullify the eigval array for all CPUs but one, so that the all_reduce
      ! operation in the end yields the correct value
      ! Of course, using a broadcast instead would be a better solution, but I'm so lazy
       if( rank_sca /= 0 ) energy(:,ispin) = 0.0_dp


       call stop_clock(timing_diago_hamiltonian)
     enddo

     deallocate(ham_local,c_matrix_local,s_matrix_local,h_small)

     call BLACS_GRIDEXIT( cntxt )

   else
     energy(:,:) = 0.0_dp
   endif


   ! Poor man distribution TODO replace by a broadcast
   call xsum_world(energy)
   call xsum_world(c_matrix)

 else ! only one proc selected
#endif

   allocate(h_small2(nstate,nstate))

   do ispin=1,nspin_local
     write(stdout,'(1x,a,i3)') 'Generalized diagonalization for spin: ',ispin
     call start_clock(timing_diago_hamiltonian)

     allocate(h_small(nbf,nstate))
     ! h_small(:,:) = MATMUL( TRANSPOSE(x_matrix(:,:)) , &
     !                          MATMUL( hamiltonian(:,:,ispin) , x_matrix(:,:) ) )

     ! H * U
     call DGEMM('N','N',nbf,nstate,nbf,1.0d0,hamiltonian(:,:,ispin),nbf, &
                                             x_matrix,nbf,      &
                                       0.0d0,h_small,nbf)
     ! U**T * (H * U)
     call DGEMM('T','N',nstate,nstate,nbf,1.0d0,x_matrix,nbf,  &
                                                h_small,nbf,            &
                                          0.0d0,h_small2,nstate)
     deallocate(h_small)

     ! H * C' = C' * E
     call diagonalize(scf_diago_flavor,h_small2,energy(:,ispin))

     !c_matrix(:,1:nstate,ispin) = MATMUL( x_matrix(:,:) , h_small2(:,:) )
     ! C = U * C'
     call DGEMM('N','N',nbf,nstate,nstate,1.0d0,x_matrix,nbf, &
                                                h_small2,nstate,       &
                                          0.0d0,c_matrix(:,:,ispin),nbf)


     call stop_clock(timing_diago_hamiltonian)
   enddo

   deallocate(h_small2)

#if defined(HAVE_SCALAPACK)
 endif
#endif

end subroutine diagonalize_hamiltonian_scalapack


!=========================================================================
subroutine diagonalize_hamiltonian_sca(desc_h,hamiltonian,desc_sqrt,x_matrix,energy,c_matrix)
 implicit none

 integer,intent(in)   :: desc_h(NDEL),desc_sqrt(NDEL)
 real(dp),intent(in)  :: hamiltonian(:,:,:)
 real(dp),intent(in)  :: x_matrix(:,:)
 real(dp),intent(out) :: c_matrix(:,:,:)
 real(dp),intent(out) :: energy(:,:)
!=====
 integer  :: ispin_min,ispin_max
 integer  :: cntxt
 integer  :: iprow,ipcol,nprow,npcol
 integer  :: nstate,nbf
 integer  :: m_ham,n_ham
 integer  :: m_c,n_c
 integer  :: ispin
 integer  :: ilocal,jlocal,jglobal
 integer  :: m_small,n_small
 integer  :: info
 real(dp),allocatable :: h_small(:,:)
!=====

#if defined(HAVE_SCALAPACK)

 cntxt  = desc_h(CTXT_)
 nbf    = desc_h(M_)
 nstate = desc_sqrt(N_)
 m_ham  = SIZE(hamiltonian, DIM=1 )
 n_ham  = SIZE(hamiltonian, DIM=2 )
 ispin_min = LBOUND(hamiltonian,DIM=3)
 ispin_max = UBOUND(hamiltonian,DIM=3)
 m_c    = SIZE(c_matrix, DIM=1 )
 n_c    = SIZE(c_matrix, DIM=2 )

 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 if(cntxt_ham > 0 ) then
   m_small = NUMROC(nstate,desc_h(MB_),iprow,desc_h(RSRC_),nprow)
   n_small = NUMROC(nstate,desc_h(NB_),ipcol,desc_h(CSRC_),npcol)
   call DESCINIT(desc_small,nstate,nstate,desc_h(MB_),desc_h(NB_),desc_h(RSRC_),desc_h(CSRC_),cntxt,m_small,info)
   allocate(h_small(m_small,n_small))

   do ispin=ispin_min,ispin_max
     write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
     call start_clock(timing_diago_hamiltonian)

     !
     ! H_small = S^{-1/2}**T H S^{-1/2}
     ! H*S**{-1/2} --> c_matrix
     call PDGEMM('N','N',nbf,nstate,nbf,                         &
                  1.0_dp,hamiltonian(1,1,ispin),1,1,desc_h,      &
                              x_matrix,1,1,desc_sqrt,   &
                  0.0_dp,   c_matrix(1,1,ispin),1,1,desc_sqrt)

     ! S**{-1/2}**T*c_matrix --> h_small
     call PDGEMM('T','N',nstate,nstate,nbf,                      &
                  1.0_dp,x_matrix,1,1,desc_sqrt,        &
                       c_matrix(1,1,ispin),1,1,desc_sqrt,        &
                  0.0_dp,          h_small,1,1,desc_small)



     call diagonalize_sca(scf_diago_flavor,h_small,desc_small,energy(:,ispin))


     !
     ! C = S^{-1/2} C_small
     call PDGEMM('N','N',nbf,nstate,nstate,                   &
                  1.0_dp,x_matrix,1,1,desc_sqrt,     &
                                 h_small,1,1,desc_small,      &
                  0.0_dp,c_matrix(1,1,ispin),1,1,desc_sqrt)


     call stop_clock(timing_diago_hamiltonian)
   enddo

   deallocate(h_small)

 else
   c_matrix(:,:,:) = 0.0_dp
 endif

 ! Ensure that all procs know the energies
 if( rank_world /= 0 ) then
   energy(:,ispin_min:ispin_max) = 0.0_dp
 endif
 call xsum_world(energy(:,ispin_min:ispin_max))

 ! Poor man distribution
 call xsum_local(c_matrix)



#endif

end subroutine diagonalize_hamiltonian_sca


!=========================================================================
!subroutine  diagonalize_hamiltonian_sca_cmplx(ispin_min,ispin_max,desc_h, &
!               hamiltonian_cmplx,desc_sqrt,x_matrix,energy,c_matrix_cmplx)
! implicit none
!
! integer,intent(in)      :: ispin_min,ispin_max
! integer,intent(in)      :: desc_h(NDEL),desc_sqrt(NDEL)
! complex(dp),intent(in)  :: hamiltonian_cmplx(:,:,:)
! real(dp),intent(in)     :: x_matrix(:,:)
! complex(dp),intent(out)    :: c_matrix(:,:,:)
! real(dp),intent(out)    :: energy(:,:)
!!=====
! integer  :: cntxt
! integer  :: iprow,ipcol,nprow,npcol
! integer  :: nstate,nbf
! integer  :: m_ham,n_ham
! integer  :: m_c,n_c
! integer  :: ispin
! integer  :: ilocal,jlocal,jglobal
! integer  :: m_small,n_small
! integer  :: info
! complex(dp),allocatable :: h_small_cmplx(:,:)
!!=====
!
!#if defined(HAVE_SCALAPACK)
!
! cntxt  = desc_h(CTXT_A)
! nbf    = desc_h(M_A)
! nstate = desc_sqrt(N_A)
! m_ham  = SIZE(hamiltonian_cmplx, DIM=1 )
! n_ham  = SIZE(hamiltonian_cmplx, DIM=2 )
! m_c    = SIZE(c_matrix_cmplx, DIM=1 )
! n_c    = SIZE(c_matrix_cmplx, DIM=2 )
!
! call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
!
! if(cntxt_ham > 0 ) then
!   m_small = NUMROC(nstate,desc_h(MB_A),iprow,desc_h(RSRC_A),nprow)
!   n_small = NUMROC(nstate,desc_h(NB_A),ipcol,desc_h(CSRC_A),npcol)
!   call DESCINIT(desc_small,nstate,nstate,desc_h(MB_A),desc_h(NB_A),desc_h(RSRC_A),desc_h(CSRC_A),cntxt,m_small,info)
!   allocate(h_small_cmplx(m_small,n_small))
!
!   do ispin=ispin_min,ispin_max
!     write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
!     call start_clock(timing_diago_hamiltonian)
!
!     !
!     ! H_small = S^{-1/2}**T H S^{-1/2}
!     ! H*S**{-1/2} --> c_matrix
!     call PZGEMM('N','N',nbf,nstate,nbf,                         &
!                  1.0_dp,hamiltonian(1,1,ispin),1,1,desc_h,      &
!                              x_matrix,1,1,desc_sqrt,   &
!                  0.0_dp,   c_matrix(1,1,ispin),1,1,desc_sqrt)
!
!     ! S**{-1/2}**T*c_matrix --> h_small
!     call PZGEMM('T','N',nstate,nstate,nbf,                      &
!                  1.0_dp,x_matrix,1,1,desc_sqrt,        &
!                       c_matrix(1,1,ispin),1,1,desc_sqrt,        &
!                  0.0_dp,          h_small,1,1,desc_small)
!
!
!
!     call diagonalize_sca(nstate,desc_small,h_small,energy(:,ispin))
!
!
!     !
!     ! C = S^{-1/2} C_small
!     call PDGEMM('N','N',nbf,nstate,nstate,                   &
!                  1.0_dp,x_matrix,1,1,desc_sqrt,     &
!                                 h_small,1,1,desc_small,      &
!                  0.0_dp,c_matrix(1,1,ispin),1,1,desc_sqrt)
!
!
!     call stop_clock(timing_diago_hamiltonian)
!   enddo
!
!   deallocate(h_small)
!
! else
!   c_matrix(:,:,:) = 0.0_dp
! endif
!
! ! Ensure that all procs know the energies
! if( rank_world /= 0 ) then
!   energy(:,ispin_min:ispin_max) = 0.0_dp
! endif
! call xsum_world(energy(:,ispin_min:ispin_max))
!
! ! Poor man distribution
! call xsum_local(c_matrix)
!
!
!
!#endif
!
!end subroutine diagonalize_hamiltonian_sca
!



!=========================================================================
subroutine setup_sqrt_overlap_sca(TOL_OVERLAP,desc_s,s_matrix, &
                                  desc_sqrt,nstate,x_matrix)
 implicit none

 real(dp),intent(in)                :: TOL_OVERLAP
 integer,intent(in)                 :: desc_s(NDEL)
 real(dp),intent(in)                :: s_matrix(:,:)
 integer,intent(out)                :: desc_sqrt(NDEL)
 integer,intent(out)                :: nstate
 real(dp),allocatable,intent(inout) :: x_matrix(:,:)
!=====
 integer              :: cntxt
 integer              :: ms,ns,msqrt,nsqrt
 integer              :: iprow,ipcol,nprow,npcol
 integer              :: nbf
 integer              :: ibf,jbf
 integer              :: ilocal,jlocal
 integer              :: iglobal,jglobal
 integer              :: info
 real(dp),allocatable :: matrix_tmp(:,:)
 real(dp),allocatable :: s_eigval(:)
 real(dp),allocatable :: diag(:,:)
!=====

#if defined(HAVE_SCALAPACK)

 write(stdout,'(/,a)') ' Calculate overlap matrix square-root S^{1/2}: SCALAPACK'

 if( desc_s(M_) /= desc_s(N_) ) call die('setup_sqrt_overlap_sca: S matrix should be square')
 cntxt = desc_s(CTXT_)
 nbf   = desc_s(M_)
 ms    = SIZE(s_matrix, DIM=1 )
 ns    = SIZE(s_matrix, DIM=2 )
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 allocate(s_eigval(nbf))

 if( cntxt > 0 ) then

   allocate(matrix_tmp(ms,ns))
   matrix_tmp(:,:) = s_matrix(:,:)

   call diagonalize_sca(' ',matrix_tmp,desc_s,s_eigval)

   nstate = COUNT( s_eigval(:) > TOL_OVERLAP )

   !
   ! Initialize the descriptor of the rectangular matric S^{-1/2}
   msqrt = NUMROC(nbf   ,desc_s(MB_),iprow,desc_s(RSRC_),nprow)
   nsqrt = NUMROC(nstate,desc_s(NB_),ipcol,desc_s(CSRC_),npcol)
   call DESCINIT(desc_sqrt,nbf,nstate,desc_s(MB_),desc_s(NB_),desc_s(RSRC_),desc_s(CSRC_),cntxt,msqrt,info)

 else
   nstate = 0
   msqrt    = 0
   nsqrt    = 0
 endif
 ! Propagate nstate
 call xmax_world(nstate)
 call xmax_local(msqrt)
 call xmax_local(nsqrt)


 call clean_allocate('Overlap X * X**H = S**-1',x_matrix,msqrt,nsqrt)

 write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
 write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
 write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',TOL_OVERLAP
 write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',nbf

 !
 ! Filter now
 if( cntxt > 0 ) then

   ! Create the diagonal matrix than transforms a square matrix into a
   ! rectangular one
   allocate(diag(msqrt,nsqrt))
   diag(:,:) = 0.0_dp

   jglobal = 0
   do iglobal=1,nbf
     if( s_eigval(iglobal) > TOL_OVERLAP ) then
       jglobal = jglobal + 1
       if( iprow /= INDXG2P(iglobal,desc_s(MB_),0,desc_s(RSRC_),nprow) ) cycle
       if( ipcol /= INDXG2P(jglobal,desc_s(NB_),0,desc_s(CSRC_),npcol) ) cycle
       ilocal = INDXG2L(iglobal,desc_s(MB_),0,desc_s(RSRC_),nprow)
       jlocal = INDXG2L(jglobal,desc_s(NB_),0,desc_s(CSRC_),npcol)

       if( ilocal * jlocal /= 0 ) then
         diag(ilocal,jlocal) = 1.0_dp / SQRT( s_eigval(iglobal) )
       endif
     endif
   enddo


   call PDGEMM('N','N',nbf,nstate,nbf,                   &
                1.0_dp,       matrix_tmp,1,1,desc_s,     &
                                    diag,1,1,desc_sqrt,  &
                0.0_dp,x_matrix,1,1,desc_sqrt)

   deallocate(diag)
   deallocate(matrix_tmp)

 else
   x_matrix(:,:) = 0.0_dp
 endif

 call xsum_local(x_matrix)

 deallocate(s_eigval)

#endif

end subroutine setup_sqrt_overlap_sca


!=========================================================================
subroutine setup_sqrt_density_matrix_sca(nbf,m_ham,n_ham,p_matrix,p_matrix_sqrt,p_matrix_occ)
 implicit none

 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: p_matrix_sqrt(m_ham,n_ham,nspin)
 real(dp),intent(out) :: p_matrix_occ(nbf,nspin)
!=====
 integer              :: ispin,ibf
!=====

#if defined(HAVE_SCALAPACK)

 write(stdout,*) 'Calculate the square root of the density matrix: SCALAPACK'
 call start_clock(timing_sqrt_density_matrix)

 if( cntxt_ham > 0 ) then
   do ispin=1,nspin
     p_matrix_sqrt(:,:,ispin) = p_matrix(:,:,ispin)
     call diagonalize_sca(' ',p_matrix_sqrt(:,:,ispin),desc_ham,p_matrix_occ(:,ispin))

     ! Cheat on the negative eigenvalues
     ! to avoid the pathological case of non-positive definite P
     do ibf=1,nbf
       p_matrix_occ(ibf,ispin) = MAX( p_matrix_occ(ibf,ispin) , TINY(1.0_dp) )
     enddo

     call matmul_diag_sca('R',SQRT(p_matrix_occ(:,ispin)),desc_ham,p_matrix_sqrt(:,:,ispin))

   enddo
 else
   p_matrix_sqrt(:,:,:) = 0.0_dp
   p_matrix_occ(:,:)    = 0.0_dp
 endif

 call xsum_local(p_matrix_sqrt)
 call xsum_local(p_matrix_occ)

 call stop_clock(timing_sqrt_density_matrix)


#endif

end subroutine setup_sqrt_density_matrix_sca


!=========================================================================
subroutine dft_approximate_vhxc_sca(basis,m_ham,n_ham,vhxc_ij)
 use m_basis_set
 use m_dft_grid
 use m_eri_calculate
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: vhxc_ij(m_ham,n_ham)
!=====
 integer              :: idft_xc
 integer              :: igrid,ibf,jbf,ispin
 real(dp)             :: rr(3)
 real(dp)             :: normalization
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf,1)
 real(dp)             :: rhor(1)
 real(dp)             :: vxc(1),exc,excr(1)
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhgau(m_ham,n_ham)
 integer              :: iatom,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 call start_clock(timing_approx_ham)


 vhxc_ij(:,:) = 0.0_dp


 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities: SCALAPACK'

 do iatom=1,natom-nprojectile
   if( rank_local /= MODULO(iatom,nproc_local) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zvalence(iatom),zatom(iatom),coeff,alpha)

   call calculate_eri_approximate_hartree(basis,m_ham,n_ham,xatom(:,iatom),ngau,coeff,alpha,vhgau)
   vhxc_ij(:,:) = vhxc_ij(:,:) + vhgau(:,:)

   deallocate(alpha,coeff)
 enddo

 write(stdout,*) 'Simple LDA functional on a coarse grid'

 !
 ! Create a temporary grid with low quality
 ! This grid is to be destroyed at the end of the present subroutine
 call init_dft_grid(basis,low,.FALSE.,.FALSE.,1)


 normalization = 0.0_dp
 exc           = 0.0_dp
 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)
   weight = w_grid(igrid)

   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r_batch(basis,igrid,basis_function_r)

   !
   ! calculate the density at point r for spin up and spin down
   call setup_atomic_density(rr,rhor(1))

   !
   ! Normalization
   normalization = normalization + rhor(1) * weight

   call teter_lda_vxc_exc(1,rhor(1:1),vxc(1:1),excr(1:1))

   !
   ! XC energy
   exc = exc + excr(1) * weight * rhor(1)

   !
   ! HXC
   do jlocal=1,n_ham
     jglobal = colindex_local_to_global('H',jlocal)
     do ilocal=1,m_ham
       iglobal = rowindex_local_to_global('H',ilocal)

       vhxc_ij(ilocal,jlocal) =  vhxc_ij(ilocal,jlocal) &
            + weight * vxc(1) * basis_function_r(iglobal,1)  &
                           * basis_function_r(jglobal,1)

     enddo
   enddo

 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum_world(normalization)
 call xsum_world(exc)
 call xsum_local(vhxc_ij)

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization
 write(stdout,  '(a,2(2x,f12.6))') '      XC energy (Ha):',exc

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()


 call stop_clock(timing_approx_ham)

end subroutine dft_approximate_vhxc_sca


!=========================================================================
end module m_hamiltonian_sca
!=========================================================================
