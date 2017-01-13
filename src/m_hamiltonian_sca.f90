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
 use m_warning
 use m_memory
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min


contains


!=========================================================================
subroutine matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,matrix_local)
 use m_basis_set
 implicit none

 integer,intent(in)     :: ibf,jbf
 integer,intent(in)     :: li,lj
 integer,intent(in)     :: ni_cart,nj_cart,ni,nj
 integer,intent(in)     :: m_ham,n_ham
 real(dp),intent(in)    :: matrix_cart(ni_cart,nj_cart)
 real(dp),intent(inout) :: matrix_local(m_ham,n_ham)
!=====
 integer  :: iglobal,jglobal
 integer  :: ilocal,jlocal
 real(dp) :: matrix_final(ibf:ibf+ni-1,jbf:jbf+nj-1)
!=====

#ifdef HAVE_SCALAPACK

 matrix_final(:,:) = MATMUL( TRANSPOSE(cart_to_pure(li)%matrix(:,:)) , &
                             MATMUL( matrix_cart(:,:) , cart_to_pure(lj)%matrix(:,:) ) )

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
subroutine setup_overlap_sca(print_matrix_,basis,m_ham,n_ham,s_matrix)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: s_matrix(m_ham,n_ham)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
!=====

#ifdef HAVE_SCALAPACK

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup overlap matrix S: SCALAPACK'

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
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call overlap_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo

     call matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,s_matrix)


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
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix)

 call stop_clock(timing_overlap)


#endif

end subroutine setup_overlap_sca


!=========================================================================
subroutine setup_kinetic_sca(print_matrix_,basis,m_ham,n_ham,hamiltonian_kinetic)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: hamiltonian_kinetic(m_ham,n_ham)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart,iglobal,jglobal,ilocal,jlocal
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp),allocatable :: matrix_final(:,:)
!=====

#ifdef HAVE_SCALAPACK

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian: SCALAPACK'

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
     allocate(matrix_final(ni,nj))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call kinetic_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
       enddo
     enddo

     call matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,hamiltonian_kinetic)

     deallocate(matrix_cart,matrix_final)
     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo
   jbf      = 1
   jbf_cart = 1

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart

 enddo

 title='===  Kinetic energy contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic)

 call stop_clock(timing_hamiltonian_kin)

#endif

end subroutine setup_kinetic_sca


!=========================================================================
subroutine setup_nucleus_sca(print_matrix_,basis,m_ham,n_ham,hamiltonian_nucleus)
 use m_basis_set
 use m_atoms
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 real(dp),intent(out)       :: hamiltonian_nucleus(m_ham,n_ham)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
!=====

#ifdef HAVE_SCALAPACK

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
     do iatom=1,natom
       if( rank_local /= MODULO(iatom-1,nproc_local) ) cycle
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           call nucleus_basis_function(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),zatom(iatom),x(:,iatom),vnucleus_ij)
           matrix_cart(i_cart,j_cart) = matrix_cart(i_cart,j_cart) + vnucleus_ij
         enddo
       enddo
     enddo

     call matrix_cart_to_local(ibf,jbf,li,lj,ni_cart,nj_cart,matrix_cart,ni,nj,m_ham,n_ham,hamiltonian_nucleus)


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
 call xsum_local(hamiltonian_nucleus)

 title='===  Nucleus potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

#endif

end subroutine setup_nucleus_sca


!=========================================================================
subroutine setup_hartree_ri_sca(print_matrix_,nbf,m_ham,n_ham,p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: hartree_ij(m_ham,n_ham)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ilocal,jlocal
 integer              :: iglobal,jglobal
 integer              :: ibf,jbf,kbf,lbf,ispin
 integer              :: ibf_auxil,ipair
 integer              :: index_ij,index_kl
 real(dp),allocatable :: partial_sum(:)
 real(dp)             :: rtmp
 character(len=100)   :: title
!=====

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity: SCALAPACK'
 call start_clock(timing_hartree)

 allocate(partial_sum(nauxil_3center))
 partial_sum(:) = 0.0_dp

 do jlocal=1,n_ham
   jglobal = colindex_local_to_global('H',jlocal)

   do ilocal=1,m_ham
     iglobal = rowindex_local_to_global('H',ilocal)
     if( negligible_basispair(iglobal,jglobal) ) cycle

     partial_sum(:) = partial_sum(:) + eri_3center(:,index_pair(iglobal,jglobal)) * SUM( p_matrix(ilocal,jlocal,:) )  

   enddo
 enddo

 call xsum_trans(partial_sum)


 ! Hartree potential is not sensitive to spin
 hartree_ij(:,:) = 0.0_dp
 do jlocal=1,n_ham
   jglobal = colindex_local_to_global('H',jlocal)

   do ilocal=1,m_ham
     iglobal = rowindex_local_to_global('H',ilocal)
     if( negligible_basispair(iglobal,jglobal) ) cycle

     hartree_ij(ilocal,jlocal) = SUM( partial_sum(:) * eri_3center(:,index_pair(iglobal,jglobal)) )

   enddo
 enddo

 
 deallocate(partial_sum)

 !
 ! Sum up the different contribution from different procs only if needed
 call xsum_local(hartree_ij)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,hartree_ij)

 !
 ! Calculate the Hartree energy
 if( cntxt_ham > 0 ) then
   ehartree = 0.5_dp*SUM(hartree_ij(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
 else
   ehartree = 0.0_dp
 endif
 call xsum_local(ehartree)


 call stop_clock(timing_hartree)


#endif

end subroutine setup_hartree_ri_sca


!=========================================================================
subroutine setup_exchange_ri_sca(nbf,nstate,m_c,n_c,m_ham,n_ham,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 use m_eri_calculate,only: nauxil_2center
 implicit none
 integer,intent(in)   :: nbf,m_ham,n_ham
 integer,intent(in)   :: nstate,m_c,n_c
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(in)  :: c_matrix(m_c,n_c,nspin)
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: exchange_ij(m_ham,n_ham,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,ispin,istate
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 real(dp)             :: c_matrix_i(nbf)
 integer              :: iglobal,ilocal,jlocal
 integer :: mlocal,nlocal
 integer :: desc3final(NDEL)
 integer :: desc3work(NDEL)
 integer :: desctmp(NDEL)
 integer :: descx(NDEL)
 integer :: mwork,nwork
 real(dp),allocatable :: tmp_local(:,:)
 real(dp),allocatable :: matrix_tmp(:,:)
 real(dp),allocatable :: sigx(:,:)
 integer :: info
 integer :: ipair_local,ipair_global,ipair
 integer :: jbf_local,jbf_global
 integer :: mbf_local,nbf_local
!=====


 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity: SCALAPACK no buffer'
 call start_clock(timing_exchange)

 exchange_ij(:,:,:) = 0.0_dp

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'SCALAPACK local grid',nprow_3center,npcol_3center
 write(stdout,*) 'This is process:',iprow_3center,ipcol_3center
 write(stdout,*) 'Size of the non-SCALAPACK 3-center integrals:',nauxil_3center,npair

 ! nauxil x npair
 mwork = NUMROC(nauxil_2center,block_row,iprow_3center,first_row,nprow_3center)
 nwork = NUMROC(npair         ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc3work,nauxil_2center,npair,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mwork),info)

! call start_clock(timing_tmp8)
! call PDGEMR2D(nauxil_2center,npair,eri_3center,1,1,desc3final,eri_3work,1,1,desc3work,cntxt)
! call stop_clock(timing_tmp8)

 write(stdout,*) 'Size of the SCALAPACK 3-center integrals:',mwork,nwork

 allocate(tmp(mwork,nbf))

 ! nauxil x nbf
 nbf_local = NUMROC(nbf,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desctmp,nauxil_2center,nbf,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mwork),info)
 allocate(tmp_local(mwork,nbf_local))

 ! nbf x nbf
 mbf_local = NUMROC(nbf,block_row,iprow_3center,first_row,nprow_3center)
 call DESCINIT(descx,nbf,nbf,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mbf_local),info)
 allocate(sigx(mbf_local,nbf_local))

 do ispin=1,nspin

   sigx(:,:) = 0.0_dp
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle

     call start_clock(timing_tmp9)
     !
     ! First all processors must have the c_matrix for (istate, ispin)
     c_matrix_i(:) = 0.0_dp
     if( cntxt_ham > 0 ) then
       jlocal = colindex_global_to_local('H',istate)
       if( jlocal /= 0 ) then
         do ilocal=1,m_c
           iglobal = rowindex_local_to_global('H',ilocal)
           c_matrix_i(iglobal) = c_matrix(ilocal,jlocal,ispin) 
         enddo
       endif
     endif
     call xsum_world(c_matrix_i)
     call stop_clock(timing_tmp9)


     call start_clock(timing_tmp1)
     tmp(:,:) = 0.0_dp
     do ipair_local=1,nwork
       ipair_global = INDXL2G(ipair_local,block_col,ipcol_3center,first_col,npcol_3center)
       ibf = index_basis(1,ipair_global)
       jbf = index_basis(2,ipair_global)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix_i(jbf) * eri_3center_sca(:,ipair_local)
       if( ibf /= jbf )  &
         tmp(:,jbf) = tmp(:,jbf) + c_matrix_i(ibf) * eri_3center_sca(:,ipair_local)
     enddo
     call stop_clock(timing_tmp1)


     call start_clock(timing_tmp2)
     call dgsum2d(cntxt_3center,'R',' ',mwork,nbf,tmp,mwork,-1,-1)


     do jbf_local=1,nbf_local
       jbf_global = INDXL2G(jbf_local,block_col,ipcol_3center,first_col,npcol_3center)
       tmp_local(:,jbf_local) = tmp(:,jbf_global)
     enddo
     call stop_clock(timing_tmp2)


     call start_clock(timing_tmp3)
     ! Sigx(:,:) = Sigx(:,:) - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact * occ(i)
     ! C = A^T * A + C
     call PDSYRK('L','T',nbf,nauxil_2center,-occupation(istate,ispin)/spin_fact,tmp_local,1,1,desctmp,1.0_dp,sigx,1,1,descx)
     call stop_clock(timing_tmp3)


   enddo

   call start_clock(timing_tmp7)
   call PDGEMR2D(nbf,nbf,sigx,1,1,descx,exchange_ij(:,:,ispin),1,1,desc_ham,cntxt_3center)
   call stop_clock(timing_tmp7)

   if( cntxt_ham > 0 ) then
     allocate(matrix_tmp(m_ham,n_ham))
     call symmetrize_matrix_sca('L',nbf,desc_ham,exchange_ij(:,:,ispin),desc_ham,matrix_tmp)
     deallocate(matrix_tmp)
   endif

 enddo
 deallocate(tmp)
 deallocate(tmp_local)

 !
 ! Calculate the exchange energy
 if( cntxt_ham > 0 ) then
   eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )
 else
   eexchange = 0.0_dp
 endif
 call xsum_world(eexchange)

#endif
 call stop_clock(timing_exchange)



end subroutine setup_exchange_ri_sca


!=========================================================================
subroutine setup_exchange_longrange_ri_sca(print_matrix_,nbf,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: occupation(nbf,nspin)
 real(dp),intent(in)  :: c_matrix(nbf,nbf,nspin)
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
!=====

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity: SCALAPACK'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 allocate(tmp(nauxil_3center_lr,nbf))

 do ispin=1,nspin

   ! Denombrate the strictly positive eigenvalues
   nocc = COUNT( occupation(:,ispin) > completely_empty )

   do istate=1,nocc
     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix(jbf,istate,ispin) * eri_3center_lr(:,ipair) * SQRT( occupation(istate,ispin) )
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center_lr(:,ipair) * SQRT( occupation(istate,ispin) )
     enddo

     exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
                        - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
   enddo

 enddo
 deallocate(tmp)

 call xsum_local(exchange_ij)

 call dump_out_matrix(print_matrix_,'=== LR Exchange contribution ===',nbf,nspin,exchange_ij)

 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

#endif

end subroutine setup_exchange_longrange_ri_sca


!=========================================================================
subroutine setup_density_matrix_sca(nbf,nstate,m_c,n_c,c_matrix,occupation,m_ham,n_ham,p_matrix)
 implicit none
 integer,intent(in)   :: nbf,nstate
 integer,intent(in)   :: m_ham,n_ham,m_c,n_c
 real(dp),intent(in)  :: c_matrix(m_c,n_c,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(out) :: p_matrix(m_ham,n_ham,nspin)
!=====
 integer  :: ispin,istate
 real(dp) :: matrix_tmp(m_ham,n_ham)
!=====

#ifdef HAVE_SCALAPACK
 call start_clock(timing_density_matrix)
 write(stdout,'(1x,a)') 'Build density matrix: SCALAPACK'

 p_matrix(:,:,:) = 0.0_dp

 if( cntxt_ham > 0 ) then
   do ispin=1,nspin
     do istate=1,nstate
       if( occupation(istate,ispin) < completely_empty ) cycle
       call PDSYR('L',nbf,occupation(istate,ispin),c_matrix(:,:,ispin),1,istate,desc_c,1,p_matrix(:,:,ispin),1,1,desc_ham)
     enddo
     call symmetrize_matrix_sca('L',nbf,desc_ham,p_matrix(:,:,ispin),desc_ham,matrix_tmp)
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
subroutine diagonalize_hamiltonian_scalapack(nspin_local,nbf,nstate,  &
                                             hamiltonian,s_matrix_sqrt_inv,energy,c_matrix)
 implicit none

 integer,intent(in)   :: nspin_local,nbf,nstate
 real(dp),intent(in)  :: hamiltonian(nbf,nbf,nspin_local)
 real(dp),intent(in)  :: s_matrix_sqrt_inv(nbf,nstate)
 real(dp),intent(out) :: c_matrix(nbf,nstate,nspin_local)
 real(dp),intent(out) :: energy(nstate,nspin_local)
!=====
 integer :: mh,nh,mc,nc,ms,ns
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 integer :: desch(NDEL),descc(NDEL),descs(NDEL)
 real(dp),allocatable :: matrix_local(:,:)
#ifdef HAVE_SCALAPACK
 integer :: cntxt
 integer :: rank_sca,nprocs_sca
#endif
 integer  :: ispin
 integer  :: ilocal,jlocal,iglobal,jglobal
 integer  :: m_small,n_small
 real(dp),allocatable :: h_small(:,:),ham_local(:,:),c_matrix_local(:,:),s_matrix_local(:,:)
!=====

#ifdef HAVE_SCALAPACK

 nprow = MIN(nprow_sd,nbf/scalapack_block_min)
 npcol = MIN(npcol_sd,nbf/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 call BLACS_PINFO( rank_sca, nprocs_sca )

 call BLACS_GET( -1, 0, cntxt )
 call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 write(stdout,'(a,i4,a,i4)') ' Generalized diagonalization using SCALAPACK with a grid',nprow,' x ',npcol

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
   ! Set up the local copy of s_matrix_sqrt_inv
   do jlocal=1,nc
     jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
     do ilocal=1,mc
       iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
       s_matrix_local(ilocal,jlocal) = s_matrix_sqrt_inv(iglobal,jglobal)
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

!     h_small(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
!                              MATMUL( hamiltonian(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )

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



     call diagonalize_sca(nstate,descs,h_small,energy(:,ispin))


!     c_matrix(:,1:nstate,ispin) = MATMUL( s_matrix_sqrt_inv(:,:) , h_small(:,:) )

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


#else

 allocate(h_small(nstate,nstate))

 do ispin=1,nspin_local
   write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
   call start_clock(timing_diago_hamiltonian)


   h_small(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
                            MATMUL( hamiltonian(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) ) 

   call diagonalize(nstate,h_small,energy(1:nstate,ispin))

   c_matrix(:,1:nstate,ispin) = MATMUL( s_matrix_sqrt_inv(:,:) , h_small(:,:) )


   call stop_clock(timing_diago_hamiltonian)
 enddo


 deallocate(h_small)

#endif

end subroutine diagonalize_hamiltonian_scalapack


!=========================================================================
subroutine diagonalize_hamiltonian_sca(ispin_min,ispin_max,desc_h,hamiltonian,desc_sqrt,s_matrix_sqrt_inv,energy,c_matrix)
 implicit none

 integer,intent(in)   :: ispin_min,ispin_max
 integer,intent(in)   :: desc_h(NDEL),desc_sqrt(NDEL)
 real(dp),intent(in)  :: hamiltonian(:,:,:)
 real(dp),intent(in)  :: s_matrix_sqrt_inv(:,:)
 real(dp),intent(out) :: c_matrix(:,:,:)
 real(dp),intent(out) :: energy(:,:)
!=====
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

#ifdef HAVE_SCALAPACK

 cntxt  = desc_h(CTXT_A)
 nbf    = desc_h(M_A)
 nstate = desc_sqrt(N_A)
 m_ham  = SIZE(hamiltonian, DIM=1 )
 n_ham  = SIZE(hamiltonian, DIM=2 )
 m_c    = SIZE(c_matrix, DIM=1 )
 n_c    = SIZE(c_matrix, DIM=2 )

 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 if(cntxt_ham > 0 ) then
   m_small = NUMROC(nstate,desc_h(MB_A),iprow,desc_h(RSRC_A),nprow)
   n_small = NUMROC(nstate,desc_h(NB_A),ipcol,desc_h(CSRC_A),npcol)
   call DESCINIT(desc_small,nstate,nstate,desc_h(MB_A),desc_h(NB_A),desc_h(RSRC_A),desc_h(CSRC_A),cntxt,m_small,info)
   allocate(h_small(m_small,n_small))

   do ispin=ispin_min,ispin_max
     write(stdout,'(a,i3)') ' Diagonalization for spin: ',ispin
     call start_clock(timing_diago_hamiltonian)

     !
     ! H_small = S^{-1/2}**T H S^{-1/2}
     call PDGEMM('N','N',nbf,nstate,nbf,                         &
                  1.0_dp,hamiltonian(1,1,ispin),1,1,desc_h,      &
                              s_matrix_sqrt_inv,1,1,desc_sqrt,   &
                  0.0_dp,   c_matrix(1,1,ispin),1,1,desc_sqrt)

     call PDGEMM('T','N',nstate,nstate,nbf,                      &
                  1.0_dp,s_matrix_sqrt_inv,1,1,desc_sqrt,        &
                       c_matrix(1,1,ispin),1,1,desc_sqrt,        &
                  0.0_dp,          h_small,1,1,desc_small)



     call diagonalize_sca(nstate,desc_small,h_small,energy(:,ispin))


     !
     ! C = S^{-1/2} C_small 
     call PDGEMM('N','N',nbf,nstate,nstate,                   &
                  1.0_dp,s_matrix_sqrt_inv,1,1,desc_sqrt,     &
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
subroutine setup_sqrt_overlap_sca(TOL_OVERLAP,desc_s,s_matrix, &
                                  desc_sqrt,nstate,s_matrix_sqrt_inv)
 use m_tools
 implicit none

 real(dp),intent(in)                :: TOL_OVERLAP
 integer,intent(in)                 :: desc_s(NDEL)
 real(dp),intent(in)                :: s_matrix(:,:)
 integer,intent(out)                :: desc_sqrt(NDEL)
 integer,intent(out)                :: nstate
 real(dp),allocatable,intent(inout) :: s_matrix_sqrt_inv(:,:)
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

#ifdef HAVE_SCALAPACK

 write(stdout,'(/,a)') ' Calculate overlap matrix square-root S^{1/2}: SCALAPACK'

 if( desc_s(M_A) /= desc_s(N_A) ) call die('setup_sqrt_overlap_sca: S matrix should be square')
 cntxt = desc_s(CTXT_A)
 nbf   = desc_s(M_A)
 ms    = SIZE(s_matrix, DIM=1 )
 ns    = SIZE(s_matrix, DIM=2 )
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 allocate(s_eigval(nbf))

 if( cntxt > 0 ) then

   allocate(matrix_tmp(ms,ns))
   matrix_tmp(:,:) = s_matrix(:,:)

   call diagonalize_sca(nbf,desc_s,matrix_tmp,s_eigval)

   nstate = COUNT( s_eigval(:) > TOL_OVERLAP )

   ! 
   ! Initialize the descriptor of the rectangular matric S^{-1/2}
   msqrt = NUMROC(nbf   ,desc_s(MB_A),iprow,desc_s(RSRC_A),nprow)
   nsqrt = NUMROC(nstate,desc_s(NB_A),ipcol,desc_s(CSRC_A),npcol)
   call DESCINIT(desc_sqrt,nbf,nstate,desc_s(MB_A),desc_s(NB_A),desc_s(RSRC_A),desc_s(CSRC_A),cntxt,msqrt,info)
   
 else
   nstate = 0
   msqrt    = 0
   nsqrt    = 0
 endif
 ! Propagate nstate
 call xmax_world(nstate)
 call xmax_local(msqrt)
 call xmax_local(nsqrt)


 call clean_allocate('Overlap sqrt S^{-1/2}',s_matrix_sqrt_inv,msqrt,nsqrt)

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
       if( iprow /= INDXG2P(iglobal,desc_s(MB_A),0,desc_s(RSRC_A),nprow) ) cycle
       if( ipcol /= INDXG2P(jglobal,desc_s(NB_A),0,desc_s(CSRC_A),npcol) ) cycle
       ilocal = INDXG2L(iglobal,desc_s(MB_A),0,desc_s(RSRC_A),nprow)
       jlocal = INDXG2L(jglobal,desc_s(NB_A),0,desc_s(CSRC_A),npcol)

       if( ilocal * jlocal /= 0 ) then
         diag(ilocal,jlocal) = 1.0_dp / SQRT( s_eigval(iglobal) )
       endif
     endif
   enddo


   call PDGEMM('N','N',nbf,nstate,nbf,                   &
                1.0_dp,       matrix_tmp,1,1,desc_s,     &
                                    diag,1,1,desc_sqrt,  &
                0.0_dp,s_matrix_sqrt_inv,1,1,desc_sqrt)

   deallocate(diag)
   deallocate(matrix_tmp)

 else
   s_matrix_sqrt_inv(:,:) = 0.0_dp
 endif

 call xsum_local(s_matrix_sqrt_inv)

 deallocate(s_eigval)

#endif

end subroutine setup_sqrt_overlap_sca


!=========================================================================
subroutine setup_sqrt_density_matrix_sca(nbf,m_ham,n_ham,p_matrix,p_matrix_sqrt,p_matrix_occ)
 use m_tools
 implicit none

 integer,intent(in)   :: nbf,m_ham,n_ham
 real(dp),intent(in)  :: p_matrix(m_ham,n_ham,nspin)
 real(dp),intent(out) :: p_matrix_sqrt(m_ham,n_ham,nspin)
 real(dp),intent(out) :: p_matrix_occ(nbf,nspin)
!=====
 integer              :: ispin,ibf
!=====

#ifdef HAVE_SCALAPACK

 write(stdout,*) 'Calculate the square root of the density matrix: SCALAPACK'
 call start_clock(timing_sqrt_density_matrix)

 if( cntxt_ham > 0 ) then
   do ispin=1,nspin
     p_matrix_sqrt(:,:,ispin) = p_matrix(:,:,ispin)
     call diagonalize_sca(nbf,desc_ham,p_matrix_sqrt(:,:,ispin),p_matrix_occ(:,ispin))

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
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: rhor
 real(dp)             :: vxc,exc,excr
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhartree
 real(dp)             :: vhgau(m_ham,n_ham)
 integer              :: iatom,igau,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 call start_clock(timing_approx_ham)


 vhxc_ij(:,:) = 0.0_dp


 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities: SCALAPACK'

 do iatom=1,natom
   if( rank_local /= MODULO(iatom,nproc_local) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),basis_element(iatom),coeff,alpha)

   call calculate_eri_approximate_hartree(basis,m_ham,n_ham,x(:,iatom),ngau,coeff,alpha,vhgau)
   vhxc_ij(:,:) = vhxc_ij(:,:) + vhgau(:,:) 

   deallocate(alpha,coeff)
 enddo

 write(stdout,*) 'Simple LDA functional on a coarse grid'

 !
 ! Create a temporary grid with low quality
 ! This grid is to be destroyed at the end of the present subroutine
 call init_dft_grid(low)

 !
 ! If it is the first time, set up the stored arrays
 !
 if( .NOT. ALLOCATED(bfr) ) call prepare_basis_functions_r(basis)

 normalization = 0.0_dp
 exc           = 0.0_dp
 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)
   weight = w_grid(igrid)

   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r(basis,igrid,basis_function_r)

   !
   ! calculate the density at point r for spin up and spin down
   call setup_atomic_density(rr,rhor,vhartree)

   !
   ! Normalization
   normalization = normalization + rhor * weight

   call teter_lda_vxc_exc(rhor,vxc,excr)

   !
   ! XC energy
   exc = exc + excr * weight * rhor

   !
   ! HXC
   do jlocal=1,n_ham
     jglobal = colindex_local_to_global('H',jlocal)
     do ilocal=1,m_ham
       iglobal = rowindex_local_to_global('H',ilocal)

       vhxc_ij(ilocal,jlocal) =  vhxc_ij(ilocal,jlocal) &
            + weight * vxc * basis_function_r(iglobal)  &
                           * basis_function_r(jglobal)

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
