!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min


contains


!=========================================================================
subroutine setup_overlap(print_matrix_,basis,s_matrix)
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
 real(dp),allocatable :: matrix_cart(:,:)
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
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,s_matrix)

 call stop_clock(timing_overlap)


end subroutine setup_overlap


!=========================================================================
subroutine setup_overlap_mixedbasis(print_matrix_,basis1,basis2,s_matrix)
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis1,basis2
 real(dp),intent(out)       :: s_matrix(basis1%nbf,basis2%nbf)
!=====
 integer              :: ibf,jbf
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 real(dp),allocatable :: matrix_cart(:,:)
!=====

 ibf_cart = 1
 jbf_cart = 1
 ibf      = 1
 jbf      = 1
 do while(ibf_cart<=basis1%nbf_cart)
   li      = basis1%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis1%gaussian_type,li)

   do while(jbf_cart<=basis2%nbf_cart)
     lj      = basis2%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis2%gaussian_type,lj)

     allocate(matrix_cart(ni_cart,nj_cart))
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call overlap_basis_function(basis1%bf(ibf_cart+i_cart-1),basis2%bf(jbf_cart+j_cart-1),matrix_cart(i_cart,j_cart))
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


end subroutine setup_overlap_mixedbasis


!=========================================================================
subroutine setup_kinetic(print_matrix_,basis,hamiltonian_kinetic)
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
 real(dp),allocatable :: matrix_cart(:,:)
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
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic)

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic


!=========================================================================
subroutine setup_nucleus(print_matrix_,basis,hamiltonian_nucleus)
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
 real(dp),allocatable :: matrix_cart(:,:)
 real(dp)             :: vnucleus_ij
!=====

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian'
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
     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
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

 !
 ! Reduce operation
 call xsum_world(hamiltonian_nucleus)

 title='===  Nucleus potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus


!=========================================================================
subroutine setup_nucleus_ecp(basis,hamiltonian_nucleus)
 use m_basis_set
 use m_atoms
 use m_dft_grid
 use m_ecp
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: iatom
 real(dp)             :: vnucleus_ij

 integer              :: iecp
 integer              :: iproj,nproj
 integer              :: mm
 integer              :: igrid
 real(dp)             :: rr(3)
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 integer              :: iradial
 integer              :: i1,n1
 real(dp)             :: xtmp,phi,cos_theta
 real(dp)             :: wxa(nradial_ecp),xa(nradial_ecp)
 real(dp)             :: w1(nangular_ecp),x1(nangular_ecp),y1(nangular_ecp),z1(nangular_ecp)
 real(dp),allocatable :: int_fixed_r(:,:)
 real(dp),external    :: real_spherical_harmonics
 integer              :: necp,ie
 logical              :: element_has_ecp
!=====

 ! Check if there are some ECP
 if( nelement_ecp == 0 ) return


 call start_clock(timing_ecp)

 !
 ! Since there will be an allreduce operation in the end, 
 ! anticipate by dividing the input value of Hnucl by the number of procs
 if( nproc_world > 1 ) then
   hamiltonian_nucleus(:,:) = hamiltonian_nucleus(:,:) / nproc_world
 endif

 n1 = nangular_ecp
 select case(nangular_ecp)
 case(6)
   call ld0006(x1,y1,z1,w1,n1)
 case(14)
   call ld0014(x1,y1,z1,w1,n1)
 case(26)
   call ld0026(x1,y1,z1,w1,n1)
 case(38)
   call ld0038(x1,y1,z1,w1,n1)
 case(50)
   call ld0050(x1,y1,z1,w1,n1)
 case(74)
   call ld0074(x1,y1,z1,w1,n1)
 case(86)
   call ld0086(x1,y1,z1,w1,n1)
 case(110)
   call ld0110(x1,y1,z1,w1,n1)
 case(146)
   call ld0146(x1,y1,z1,w1,n1)
 case(170)
   call ld0170(x1,y1,z1,w1,n1)
 case(230)
   call ld0230(x1,y1,z1,w1,n1)
 case(302)
   call ld0302(x1,y1,z1,w1,n1)
 case(434)
   call ld0434(x1,y1,z1,w1,n1)
 case default
   write(stdout,*) 'grid points: ',nangular_ecp
   call die('setup_nucleus_ecp: Lebedev grid is not available')
 end select


 do iradial=1,nradial_ecp
   xtmp = ( iradial - 0.5_dp ) / REAL(nradial_ecp,dp)
   xa(iradial)   = -5.0_dp * log( 1.0_dp - xtmp**3)
   wxa(iradial)  = 3.0_dp * 5.0_dp * xtmp**2 / ( 1.0_dp - xtmp**3 ) / REAL(nradial_ecp,dp)
 enddo


 do iatom=1,natom
   element_has_ecp = .FALSE.
   do ie=1,nelement_ecp
     if( element_ecp(ie) == basis_element(iatom) ) then
       element_has_ecp = .TRUE.
       exit
     endif
   enddo 

   if( .NOT. element_has_ecp ) cycle

   necp = ecp(ie)%necp
     

   nproj = 0
   do iecp=1,necp
     nproj = nproj + number_basis_function_am('PURE',ecp(ie)%lk(iecp))
   enddo
  
  
   allocate(int_fixed_r(basis%nbf,nproj))
   do iradial=1,nradial_ecp
     if( MODULO(iradial-1,nproc_world) /= rank_world ) cycle

     int_fixed_r(:,:) = 0.0_dp
     do i1=1,nangular_ecp
       rr(1) = xa(iradial) * x1(i1) + x(1,iatom)
       rr(2) = xa(iradial) * y1(i1) + x(2,iatom)
       rr(3) = xa(iradial) * z1(i1) + x(3,iatom)
       call calculate_basis_functions_r(basis,rr,basis_function_r)
  
       cos_theta = z1(i1)
       phi       = ATAN2(y1(i1),x1(i1))
  
       iproj = 0
       do iecp=1,necp
         do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
           iproj = iproj + 1
           int_fixed_r(:,iproj) = int_fixed_r(:,iproj) + basis_function_r(:) &
                                     * real_spherical_harmonics(ecp(ie)%lk(iecp),mm,cos_theta,phi) &
                                        * w1(i1) * 4.0_dp * pi  
         enddo
       enddo
     enddo
  
  
     iproj = 0
     do iecp=1,necp
       do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
         iproj = iproj + 1
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             hamiltonian_nucleus(ibf,jbf) = hamiltonian_nucleus(ibf,jbf)  &
                 + int_fixed_r(ibf,iproj) * int_fixed_r(jbf,iproj) * wxa(iradial) * xa(iradial)**2  &
                    * ecp(ie)%dk(iecp) * EXP( -ecp(ie)%zetak(iecp) * xa(iradial)**2 ) * xa(iradial)**(ecp(ie)%nk(iecp)-2)
           enddo
         enddo
       enddo
     enddo
  
   enddo
  
   deallocate(int_fixed_r)

 enddo 

 call xsum_world(hamiltonian_nucleus)


 call stop_clock(timing_ecp)

end subroutine setup_nucleus_ecp


!=========================================================================
subroutine setup_hartree(print_matrix_,nbf,p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: hartree_ij(nbf,nbf)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Hartree term'
 call start_clock(timing_hartree)

 hartree_ij(:,:)=0.0_dp

 do jbf=1,nbf
   do ibf=1,nbf
     if( negligible_basispair(ibf,jbf) ) cycle
     do lbf=1,nbf
       !
       ! symmetry k <-> l
       do kbf=1,lbf-1 ! nbf
         if( negligible_basispair(kbf,lbf) ) cycle
         !
         ! symmetry (ij|kl) = (kl|ij) has been used to loop in the fast order
         hartree_ij(ibf,jbf) = hartree_ij(ibf,jbf) &
                    + eri(kbf,lbf,ibf,jbf) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
       enddo
       hartree_ij(ibf,jbf) = hartree_ij(ibf,jbf) &
                  + eri(lbf,lbf,ibf,jbf) * SUM( p_matrix(lbf,lbf,:) )
     enddo
   enddo
 enddo


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,hartree_ij)

 ehartree = 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree


!=========================================================================
subroutine setup_hartree_ri(print_matrix_,nbf,p_matrix,hartree_ij,ehartree)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: hartree_ij(nbf,nbf)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 integer              :: ibf_auxil,ipair
 integer              :: index_ij,index_kl
 real(dp),allocatable :: partial_sum(:)
 real(dp)             :: rtmp
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity'
 call start_clock(timing_hartree)


 allocate(partial_sum(nauxil_3center))
 partial_sum(:) = 0.0_dp
 do ipair=1,npair
   kbf = index_basis(1,ipair)
   lbf = index_basis(2,ipair)
   ! Factor 2 comes from the symmetry of p_matrix
   partial_sum(:) = partial_sum(:) + eri_3center(:,ipair) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
   ! Then diagonal terms have been counted twice and should be removed once.
   if( kbf == lbf ) &
     partial_sum(:) = partial_sum(:) - eri_3center(:,ipair) * SUM( p_matrix(kbf,kbf,:) )
 enddo

 ! Hartree potential is not sensitive to spin
 hartree_ij(:,:) = 0.0_dp
 do ipair=1,npair
   ibf = index_basis(1,ipair)
   jbf = index_basis(2,ipair)
   rtmp = DOT_PRODUCT( eri_3center(:,ipair) , partial_sum(:) )
   hartree_ij(ibf,jbf) = rtmp
   hartree_ij(jbf,ibf) = rtmp
 enddo

 deallocate(partial_sum)

 !
 ! Sum up the different contribution from different procs only if needed
 call xsum_auxil(hartree_ij)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,1,hartree_ij)

 ehartree = 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree_ri


!=========================================================================
subroutine setup_exchange(print_matrix_,nbf,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 write(stdout,*) 'Calculate Exchange term'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 do ispin=1,nspin
   do jbf=1,nbf
     do lbf=1,nbf
       if( negligible_basispair(lbf,jbf) ) cycle
       do kbf=1,nbf
!         if( ABS(p_matrix(kbf,lbf,ispin)) <  1.0e-12_dp ) cycle 
         do ibf=1,nbf
           if( negligible_basispair(ibf,kbf) ) cycle
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           exchange_ij(ibf,jbf,ispin) = exchange_ij(ibf,jbf,ispin) &
                      - eri(ibf,kbf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
 enddo


 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange


!=========================================================================
subroutine setup_exchange_ri(nbf,nstate,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
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

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)

 exchange_ij(:,:,:) = 0.0_dp

 allocate(tmp(nauxil_3center,nbf))

 do ispin=1,nspin

   ! Find highest occupied state
   nocc = 0
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty)  cycle
     nocc = istate
   enddo


   do istate=1,nocc
     if( MODULO( istate-1 , nproc_ortho ) /= rank_ortho ) cycle

     tmp(:,:) = 0.0_dp
     !$OMP PARALLEL PRIVATE(ibf,jbf) 
     !$OMP DO REDUCTION(+:tmp)
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix(jbf,istate,ispin) * eri_3center(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center(:,ipair)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL

     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C ; C - exchange_ij(:,:,ispin); A - tmp
     call DSYRK('L','T',nbf,nauxil_3center,-occupation(istate,ispin)/spin_fact,tmp,nauxil_3center,1.0_dp,exchange_ij(:,:,ispin),nbf)
   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo

 enddo
 deallocate(tmp)

 call xsum_world(exchange_ij)

 eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri


!=========================================================================
subroutine setup_exchange_ri_cmplx(nbf,nstate,occupation,c_matrix_cmplx,p_matrix_cmplx, &
                                   exchange_ij_cmplx,eexchange)
 use m_eri
 implicit none
 integer,intent(in)         :: nbf,nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(out)       :: eexchange
 complex(dpc),intent(in)    :: c_matrix_cmplx(nbf,nstate,nspin)
 complex(dpc),intent(in)    :: p_matrix_cmplx(nbf,nbf,nspin)
 complex(dpc),intent(out)   :: exchange_ij_cmplx(nbf,nbf,nspin)
!=====
 integer                    :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer                    :: index_ij
 integer                    :: nocc
 real(dp),allocatable       :: tmp(:,:)
 real(dp)                   :: eigval(nbf)
 integer                    :: ipair
 complex(dpc),allocatable   :: tmp_cmplx(:,:)
!=====

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)

 exchange_ij_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )
! allocate(tmp(nauxil_3center,nbf))
 allocate(tmp_cmplx(nauxil_3center,nbf))
 do ispin=1,nspin

   ! Find highest occupied state
   nocc = 0
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty)  cycle
     nocc = istate
   enddo


   do istate=1,nocc
     if( MODULO( istate-1 , nproc_ortho ) /= rank_ortho ) cycle

     tmp_cmplx(:,:) = 0.0_dp
     !$OMP PARALLEL PRIVATE(ibf,jbf) 
     !$OMP DO REDUCTION(+:tmp_cmplx)
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp_cmplx(:,ibf) = tmp_cmplx(:,ibf) + c_matrix_cmplx(jbf,istate,ispin) * eri_3center(:,ipair)
       if( ibf /= jbf ) &
            tmp_cmplx(:,jbf) = tmp_cmplx(:,jbf) + c_matrix_cmplx(ibf,istate,ispin) * eri_3center(:,ipair)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL

     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C ; C - exchange_ij(:,:,ispin); A - tmp
     call ZHERK('L','C',nbf,nauxil_3center,-occupation(istate,ispin)/spin_fact,tmp_cmplx,nauxil_3center,1.0_dp,exchange_ij_cmplx(:,:,ispin),nbf)
   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij_cmplx(ibf,jbf,ispin) = exchange_ij_cmplx(jbf,ibf,ispin)
     enddo
   enddo

 enddo ! end of loop do ispin=1,nspin
! deallocate(tmp)
! call xsum_world(exchange_ij)
 ! This interface should work also for complex exchange_ij_cmplx 
 deallocate(tmp_cmplx)
 call xsum_world(exchange_ij_cmplx)
 !!!!! CHECK THAT IM(EEXCHANGE)=0.0
 eexchange = real( 0.5_dp * SUM( exchange_ij_cmplx(:,:,:) * p_matrix_cmplx(:,:,:) ),dp)

 write(stdout,*) "This is subroutine setup_exchange_ri_cmplx and the eexchange value is: "
 write(stdout,*) eexchange, 0.5_dp * SUM(exchange_ij_cmplx(:,:,:) * p_matrix_cmplx(:,:,:))
 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri_cmplx

!=========================================================================
subroutine setup_exchange_longrange_ri(nbf,nstate,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
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

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:) = 0.0_dp

 allocate(tmp(nauxil_3center_lr,nbf))

 do ispin=1,nspin

   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty)  cycle
     if( MODULO( istate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix(jbf,istate,ispin) * eri_3center_lr(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center_lr(:,ipair)
     enddo

     !exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                   - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nauxil_3center_lr,-occupation(istate,ispin)/spin_fact,tmp,nauxil_3center_lr,1.0_dp,exchange_ij(:,:,ispin),nbf)
   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo

 enddo
 deallocate(tmp)

 call xsum_world(exchange_ij)

 eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange_ri


!=========================================================================
subroutine setup_exchange_longrange(print_matrix_,nbf,p_matrix,exchange_ij,eexchange)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
!=====

 write(stdout,*) 'Calculate Long-Range Exchange term'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           exchange_ij(ibf,jbf,ispin) = exchange_ij(ibf,jbf,ispin) &
                      - eri_lr(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
 enddo


 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange


!=========================================================================
subroutine read_external_potential(print_matrix_,nbf,nspin,p_matrix,pot_read,eread)
 use m_eri
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_read(nbf,nbf,nspin)
 real(dp),intent(out) :: eread
!=====
 integer              :: ibf,jbf,ispin
 character(len=100)   :: title
 logical              :: file_exists
 integer              :: potfile
!=====

 pot_read(:,:,:)=0.0_dp

 inquire(file='manual_potential',exist=file_exists)
 if(file_exists) then
   open(newunit=potfile,file='manual_potential',status='old')
   do ispin=1,nspin
     do jbf=1,nbf
       do ibf=1,nbf
         read(potfile,*) pot_read(ibf,jbf,ispin)
       enddo
     enddo
   enddo
   close(potfile)

 else
   call die('file not found: manual_potential')
 endif


 title='=== Read potential contribution ==='
 call dump_out_matrix(print_matrix_,title,nbf,nspin,pot_read)

 eread = 0.5_dp*SUM(pot_read(:,:,:)*p_matrix(:,:,:))

end subroutine read_external_potential


!=========================================================================
subroutine setup_density_matrix(nbf,nstate,c_matrix,occupation,p_matrix)
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(out) :: p_matrix(nbf,nbf,nspin)
!=====
 integer :: ispin,ibf,jbf
 integer :: istate
!=====

 call start_clock(timing_density_matrix)
 write(stdout,'(1x,a)') 'Build density matrix'

 p_matrix(:,:,:) = 0.0_dp
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     call DSYR('L',nbf,occupation(istate,ispin),c_matrix(:,istate,ispin),1,p_matrix(:,:,ispin),nbf)
   enddo

   ! Symmetrize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix(jbf,ibf,ispin) = p_matrix(ibf,jbf,ispin)
     enddo
   enddo
 enddo
 call stop_clock(timing_density_matrix)


end subroutine setup_density_matrix

!=========================================================================
subroutine setup_density_matrix_cmplx(nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)
 implicit none
 integer,intent(in)   :: nbf,nstate
 complex(dpc),intent(in)  :: c_matrix_cmplx(nbf,nstate,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
 complex(dpc),intent(out) :: p_matrix_cmplx(nbf,nbf,nspin)
!=====
 integer :: ispin,ibf,jbf
 integer :: istate

!=====

 call start_clock(timing_density_matrix)
 write(stdout,'(1x,a)') 'Build density matrix'

 p_matrix_cmplx(:,:,:) = 0.0_dp
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     call ZHER('L',nbf,occupation(istate,ispin),c_matrix_cmplx(:,istate,ispin),1,p_matrix_cmplx(:,:,ispin),nbf)
   enddo

   ! Symmetrize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix_cmplx(jbf,ibf,ispin) = p_matrix_cmplx(ibf,jbf,ispin)
     enddo
   enddo
 enddo
 call stop_clock(timing_density_matrix)



end subroutine setup_density_matrix_cmplx

!=========================================================================
subroutine setup_density_matrix_cmplx_slow(nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)
 implicit none
 integer,intent(in)   :: nbf,nstate
 complex(dpc),intent(in)  :: c_matrix_cmplx(nbf,nstate,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
 complex(dpc),intent(out) :: p_matrix_cmplx(nbf,nbf,nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       p_matrix_cmplx(ibf,jbf,ispin) = SUM( occupation(:,ispin) * c_matrix_cmplx(ibf,:,ispin) * conjg(c_matrix_cmplx(jbf,:,ispin)) )
     enddo
   enddo
 enddo


end subroutine setup_density_matrix_cmplx_slow

!=========================================================================
subroutine test_density_matrix(nbf,nspin,p_matrix,s_matrix)
 implicit none
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin),s_matrix(nbf,nbf)
!=====
 integer              :: ispin
 real(dp)             :: matrix(nbf,nbf)
 character(len=100)   :: title
!=====

 write(stdout,*) 'Check equality PSP = P'
 write(stdout,*) ' valid only for integer occupation numbers'
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
subroutine set_occupation(nstate,temperature,electrons,magnetization,energy,occupation)
 use m_inputparam,only: print_matrix_
 implicit none
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: electrons,magnetization,temperature
 real(dp),intent(in)  :: energy(nstate,nspin)
 real(dp),intent(out) :: occupation(nstate,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin)
 integer              :: istate,nlines,ilines
 logical              :: file_exists
 integer              :: occfile
!=====

 if( temperature < 1.0e-8_dp ) then

   occupation(:,:)=0.0_dp
  
   inquire(file='manual_occupations',exist=file_exists)
  
   if(.NOT. file_exists) then
     remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
     if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)
  
     do istate=1,nstate
       occupation(istate,:) = MIN(remaining_electrons(:), spin_fact)
       remaining_electrons(:)  = remaining_electrons(:) - occupation(istate,:)
     end do
   else
     write(stdout,*)
     write(stdout,*) 'occupations are read from file: manual_occupations'
     msg='reading manual occupations from file'
     call issue_warning(msg)
     open(newunit=occfile,file='manual_occupations',status='old')
     !
     ! read nlines, all other occupations are set to zero
     read(occfile,*) nlines
     do ilines=1,nlines
       read(occfile,*) occupation(ilines,:)
     enddo
     close(occfile)
     write(stdout,*) 'occupations set, closing file'
   endif

 else
   !
   ! Finite temperature case
   !
   call die('Finite temperature not implemented yet')

 endif

 !
 ! final check
 if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0e-7_dp ) then
   write(stdout,*) 'occupation set up failed to give the right number of electrons'
   write(stdout,*) 'sum of occupations',SUM(occupation(:,:))
   write(stdout,*) 'electrons',electrons
   do istate=1,nstate
     write(stdout,*) istate,occupation(istate,:)
   enddo
   call die('FAILURE in set_occupation')
 endif 

 call dump_out_occupation('=== Occupations ===',nstate,nspin,occupation)

end subroutine set_occupation


!=========================================================================
subroutine matrix_basis_to_eigen(nbf,nstate,c_matrix,matrix_inout)
 implicit none
 integer,intent(in)      :: nbf,nstate
 real(dp),intent(in)     :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(inout)  :: matrix_inout(nbf,nbf,nspin)
!=====
 integer                 :: ispin
!=====


 do ispin=1,nspin
   matrix_inout(1:nstate,1:nstate,ispin) = MATMUL( TRANSPOSE( c_matrix(:,:,ispin) ) , MATMUL( matrix_inout(:,:,ispin) , c_matrix(:,:,ispin) ) )
 enddo


end subroutine matrix_basis_to_eigen


!=========================================================================
subroutine evaluate_s2_operator(nbf,nstate,occupation,c_matrix,s_matrix)
 implicit none
 integer,intent(in)      :: nbf,nstate
 real(dp),intent(in)     :: occupation(nstate,nspin)
 real(dp),intent(in)     :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)     :: s_matrix(nbf,nbf)
!=====
 integer                 :: ispin,istate,jstate
 real(dp)                :: s2,s2_exact
 real(dp)                :: n1,n2,nmax,nmin
!=====

 if(nspin /= 2) return

 n1 = SUM( occupation(:,1) )
 n2 = SUM( occupation(:,2) )
 nmax = MAX(n1,n2)
 nmin = MIN(n1,n2)

 s2_exact = (nmax-nmin)/2.0_dp * ( (nmax-nmin)/2.0_dp + 1.0_dp )
 s2       = s2_exact + nmin
 do istate=1,nstate
   if( occupation(istate,1) < completely_empty ) cycle
   do jstate=1,nstate
     if( occupation(jstate,2) < completely_empty ) cycle

     s2 = s2 - ABS( DOT_PRODUCT( c_matrix(:,istate,1) , MATMUL( s_matrix(:,:) , c_matrix(:,jstate,2) ) )  &
                      * occupation(istate,1) * occupation(jstate,2) )**2

   enddo
 enddo


 write(stdout,'(/,a,f8.4)') ' Total Spin S**2: ',s2
 write(stdout,'(a,f8.4)')   ' Instead of:      ',s2_exact


end subroutine evaluate_s2_operator


!=========================================================================
subroutine level_shifting_up(nbf,nstate,s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
 implicit none
 integer,intent(in)     :: nbf,nstate
 real(dp),intent(in)    :: s_matrix(nbf,nbf)
 real(dp),intent(in)    :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)    :: occupation(nstate,nspin)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(nbf,nbf,nspin)
!=====
 integer  :: ispin
 integer  :: ibf,istate
 real(dp) :: sqrt_level_shifting(nstate)
 real(dp) :: matrix_tmp(nbf,nbf)
!=====

 write(stdout,'(/,a)')     ' Level shifting switched on'
 write(stdout,'(a,f12.6)') '   energy shift (eV):',level_shifting_energy * Ha_eV

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif


 do ispin=1,nspin
   !
   ! Shift up empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   ! 
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + matrix_tmp(:,:)

 enddo
 

end subroutine level_shifting_up


!=========================================================================
subroutine level_shifting_down(nbf,nstate,s_matrix,c_matrix,occupation,level_shifting_energy,energy,hamiltonian)
 implicit none
 integer,intent(in)     :: nbf,nstate
 real(dp),intent(in)    :: s_matrix(nbf,nbf)
 real(dp),intent(in)    :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)    :: occupation(nstate,nspin)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: energy(nstate,nspin)
 real(dp),intent(inout) :: hamiltonian(nbf,nbf,nspin)
!=====
 integer  :: ispin
 integer  :: ibf,istate
 real(dp) :: sqrt_level_shifting(nstate)
 real(dp) :: matrix_tmp(nbf,nbf)
!=====

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif

 !
 ! Shift down the energies of the virtual orbitals
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       energy(istate,ispin) = energy(istate,ispin) - level_shifting_energy
     endif
   enddo
 enddo


 do ispin=1,nspin
   !
   ! Shift down empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   ! 
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) - matrix_tmp(:,:)

 enddo
 

end subroutine level_shifting_down


!=========================================================================
subroutine setup_sqrt_overlap(TOL_OVERLAP,nbf,s_matrix,nstate,s_matrix_sqrt_inv)
 use m_tools
 implicit none

 real(dp),intent(in)                :: TOL_OVERLAP
 integer,intent(in)                 :: nbf
 real(dp),intent(in)                :: s_matrix(nbf,nbf)
 integer,intent(out)                :: nstate
 real(dp),allocatable,intent(inout) :: s_matrix_sqrt_inv(:,:)
!=====
 integer  :: ibf,jbf
 real(dp) :: s_eigval(nbf)
 real(dp) :: matrix_tmp(nbf,nbf)
!=====

 write(stdout,'(/,a)') ' Calculate overlap matrix square-root S^{1/2}'

 matrix_tmp(:,:) = s_matrix(:,:)
 ! Diagonalization with or without SCALAPACK
 call diagonalize_scalapack(scalapack_block_min,nbf,matrix_tmp,s_eigval)

 nstate = COUNT( s_eigval(:) > TOL_OVERLAP )

 call clean_allocate('Overlap sqrt S^{-1/2}',s_matrix_sqrt_inv,nbf,nstate)

 write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
 write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
 write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',TOL_OVERLAP
 write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',nbf

 ibf=0
 do jbf=1,nbf
   if( s_eigval(jbf) > TOL_OVERLAP ) then
     ibf = ibf + 1
     s_matrix_sqrt_inv(:,ibf) = matrix_tmp(:,jbf) / SQRT( s_eigval(jbf) )
   endif
 enddo


end subroutine setup_sqrt_overlap


!=========================================================================
subroutine setup_sqrt_density_matrix(nbf,p_matrix,p_matrix_sqrt,p_matrix_occ)
 use m_tools
 implicit none

 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: p_matrix_sqrt(nbf,nbf,nspin)
 real(dp),intent(out) :: p_matrix_occ(nbf,nspin)
!=====
 integer              :: ispin,ibf
!=====

 write(stdout,*) 'Calculate the square root of the density matrix'
 call start_clock(timing_sqrt_density_matrix)

 do ispin=1,nspin
   p_matrix_sqrt(:,:,ispin) = p_matrix(:,:,ispin)
   ! Diagonalization with or without SCALAPACK
   call diagonalize_scalapack(scalapack_block_min,nbf,p_matrix_sqrt(:,:,ispin),p_matrix_occ(:,ispin))
   do ibf=1,nbf
     ! this is to avoid instabilities
     if( p_matrix_occ(ibf,ispin) < 1.0e-8_dp ) then
       p_matrix_occ(ibf,ispin)    = 0.0_dp
       p_matrix_sqrt(:,ibf,ispin) = 0.0_dp
     else
       p_matrix_sqrt(:,ibf,ispin) = p_matrix_sqrt(:,ibf,ispin) * SQRT( p_matrix_occ(ibf,ispin) )
     endif
   enddo
 enddo

 call stop_clock(timing_sqrt_density_matrix)

end subroutine setup_sqrt_density_matrix


!=========================================================================
subroutine dft_exc_vxc(basis,nstate,occupation,c_matrix,p_matrix,vxc_ij,exc_xc)
 use m_definitions
 use m_mpi
 use m_timing
 use m_inputparam
 use m_basis_set
 use m_dft_grid
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: vxc_ij(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: exc_xc
!=====

 real(dp),parameter :: TOL_RHO=1.0e-10_dp
 integer  :: idft_xc
 logical  :: require_gradient,require_laplacian
 integer  :: igrid,ibf,jbf,ispin
 real(dp) :: normalization(nspin)
 real(dp) :: weight

 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: basis_function_gradr(3,basis%nbf)
 real(dp)             :: basis_function_laplr(3,basis%nbf)

 real(dp)             :: rhor(nspin)
 real(dp)             :: grad_rhor(3,nspin)
 real(dp)             :: sigma(2*nspin-1)
 real(dp)             :: tau(nspin),lapl_rhor(nspin)
 real(dp)             :: vxc_libxc(nspin)
 real(dp)             :: exc_libxc(1)
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vlapl_rho(nspin),vtau(nspin)
 real(dp)             :: dedd_r(nspin)
 real(dp)             :: dedgd_r(3,nspin)
 real(dp)             :: gradtmp(basis%nbf)
!=====

 exc_xc = 0.0_dp
 vxc_ij(:,:,:) = 0.0_dp
 if( ndft_xc == 0 ) return

 call start_clock(timing_dft)


#ifdef HAVE_LIBXC

 write(stdout,*) 'Calculate DFT XC potential'
 
 require_gradient =.FALSE.
 require_laplacian=.FALSE.
 do idft_xc=1,ndft_xc
   if( ABS(dft_xc_coef(idft_xc)) < 1.0e-6_dp ) cycle

   if(xc_f90_info_family(calc_type%xc_info(idft_xc)) == XC_FAMILY_GGA     ) require_gradient  =.TRUE.
   if(xc_f90_info_family(calc_type%xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) require_gradient  =.TRUE.
   if(xc_f90_info_family(calc_type%xc_info(idft_xc)) == XC_FAMILY_MGGA    ) require_laplacian =.TRUE.

 enddo


 !
 ! If it is the first time, then set up the stored arrays
 !
 if( .NOT. ALLOCATED(bfr) )                          call prepare_basis_functions_r(basis)
 if( require_gradient  .AND. .NOT. ALLOCATED(bfgr) ) call prepare_basis_functions_gradr(basis)
 if( require_laplacian .AND. .NOT. ALLOCATED(bfgr) ) call prepare_basis_functions_laplr(basis)

 normalization(:)=0.0_dp

 do igrid=1,ngrid

   weight = w_grid(igrid)

   !
   ! Get the functions at point r
   call get_basis_functions_r(basis,igrid,basis_function_r)
   !
   ! calculate the density at point r for spin up and spin down
   call calc_density_r(nspin,basis%nbf,nstate,occupation,c_matrix,basis_function_r,rhor)

   ! Skip all the rest if the density is too small
   if( ALL( rhor(:) < TOL_RHO )  ) cycle

   !
   ! Get the gradient and laplacian at point r
   if( require_gradient ) then
     call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
   endif
   if( require_laplacian ) then
     call get_basis_functions_laplr(basis,igrid,basis_function_gradr,basis_function_laplr)
   endif


   !
   ! Normalization
   normalization(:) = normalization(:) + rhor(:) * weight


   if( require_gradient ) then 
     call calc_density_gradr(nspin,basis%nbf,nstate,occupation,c_matrix,basis_function_r,basis_function_gradr,grad_rhor)
   endif

   if( require_laplacian ) then
     call calc_density_gradr_laplr(nspin,basis%nbf,p_matrix,basis_function_r,basis_function_gradr,basis_function_laplr,grad_rhor,tau,lapl_rhor)
   endif

   if( require_gradient .OR. require_laplacian ) then
     sigma(1) = SUM( grad_rhor(:,1)**2 )
     if(nspin==2) then
       sigma(2) = SUM( grad_rhor(:,1) * grad_rhor(:,2) )
       sigma(3) = SUM( grad_rhor(:,2)**2 )
     endif
   endif

   !
   ! LIBXC calls
   !
   dedd_r(:)    = 0.0_dp
   dedgd_r(:,:) = 0.0_dp

   do idft_xc=1,ndft_xc
     if( ABS(dft_xc_coef(idft_xc)) < 1.0e-6_dp ) cycle

     select case(xc_f90_info_family(calc_type%xc_info(idft_xc)))

     case(XC_FAMILY_LDA)
       if( dft_xc_type(idft_xc) < 1000 ) then 
         call xc_f90_lda_exc_vxc(calc_type%xc_func(idft_xc),1,rhor(1),exc_libxc(1),vxc_libxc(1))
       else
         call my_lda_exc_vxc(nspin,dft_xc_type(idft_xc),rhor,exc_libxc(1),vxc_libxc)
       endif

     case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
       if( dft_xc_type(idft_xc) < 2000 ) then 
         !
         ! Remove too small densities to stabilize the computation
         ! especially useful for Becke88
         if( ANY( rhor(:) > 1.0e-9_dp ) ) then
           call xc_f90_gga_exc_vxc(calc_type%xc_func(idft_xc),1,rhor(1),sigma(1),exc_libxc(1),vxc_libxc(1),vsigma(1))
         else
           exc_libxc(:)     = 0.0_dp
           vxc_libxc(:)     = 0.0_dp
           vsigma(:)        = 0.0_dp
         endif
       else
         call my_gga_exc_vxc_hjs(gamma_hybrid,rhor(1),sigma(1),exc_libxc(1),vxc_libxc(1),vsigma(1))
       endif

     case(XC_FAMILY_MGGA)
       call xc_f90_mgga_vxc(calc_type%xc_func(idft_xc),1,rhor(1),sigma(1),lapl_rhor(1),tau(1),vxc_libxc(1),vsigma(1),vlapl_rho(1),vtau(1))
       ! no expression for the energy
       exc_libxc(1) = 0.0_dp

     case default
       call die('functional is not LDA nor GGA nor hybrid nor meta-GGA')
     end select

     exc_xc = exc_xc + weight * exc_libxc(1) * SUM( rhor(:) ) * dft_xc_coef(idft_xc)

     dedd_r(:) = dedd_r(:) + vxc_libxc(:) * dft_xc_coef(idft_xc)

     !
     ! Set up divergence term if needed (GGA case)
     !
     if( xc_f90_info_family(calc_type%xc_info(idft_xc)) == XC_FAMILY_GGA &
        .OR. xc_f90_info_family(calc_type%xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) then
       if(nspin==1) then

         dedgd_r(:,1) = dedgd_r(:,1) + 2.0_dp * vsigma(1) * grad_rhor(:,1) * dft_xc_coef(idft_xc)

       else

         dedgd_r(:,1) = dedgd_r(:,1) + 2.0_dp * vsigma(1) * grad_rhor(:,1) * dft_xc_coef(idft_xc) &
                               + vsigma(2) * grad_rhor(:,2)

         dedgd_r(:,2) = dedgd_r(:,2) + 2.0_dp * vsigma(3) * grad_rhor(:,2) * dft_xc_coef(idft_xc) &
                               + vsigma(2) * grad_rhor(:,1)
       endif

     endif


   enddo ! loop on the XC functional


   !
   ! Eventually set up the vxc term
   !
   if( .NOT. require_gradient ) then 
     ! LDA
     do ispin=1,nspin
       call DSYR('L',basis%nbf,weight*dedd_r(ispin),basis_function_r,1,vxc_ij(:,:,ispin),basis%nbf)
     enddo

   else 
     ! GGA
     do ispin=1,nspin

       gradtmp(:) = MATMUL( dedgd_r(:,ispin) , basis_function_gradr(:,:) )
       call DSYR('L',basis%nbf,weight*dedd_r(ispin),basis_function_r,1,vxc_ij(:,:,ispin),basis%nbf)
       call DSYR2('L',basis%nbf,weight,basis_function_r,1,gradtmp,1,vxc_ij(:,:,ispin),basis%nbf)

     enddo
   endif

 enddo ! loop on the grid point




 ! Symmetrize now
 do ispin=1,nspin
   do jbf=1,basis%nbf
     do ibf=jbf+1,basis%nbf 
       vxc_ij(jbf,ibf,ispin) = vxc_ij(ibf,jbf,ispin)
     enddo
   enddo
 enddo

 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(vxc_ij)
 call xsum_grid(exc_xc)

! !
! ! Destroy operations
! do idft_xc=1,ndft_xc
!   call xc_f90_func_end(calc_type%xc_func(idft_xc))
! enddo

#else
!TODO write a call to teter to have MOLGW working without LIBXC
!   call teter_lda_vxc_exc(rhor,vxc,exc)
 write(stdout,*) 'XC energy and potential set to zero'
 write(stdout,*) 'LIBXC is not present'
#endif

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
 write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

 call stop_clock(timing_dft)

end subroutine dft_exc_vxc


!=========================================================================
subroutine dft_approximate_vhxc(basis,vhxc_ij)
 use m_basis_set
 use m_dft_grid
 use m_eri_calculate
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: vhxc_ij(basis%nbf,basis%nbf)
!=====
 integer              :: idft_xc
 integer              :: igrid,ibf,jbf,ispin
 real(dp)             :: rr(3)
 real(dp)             :: normalization
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: rhor
 real(dp)             :: vxc,excr,exc
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhartree
 real(dp)             :: vhgau(basis%nbf,basis%nbf)
 integer              :: iatom,igau,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
!=====

 call start_clock(timing_approx_ham)

 vhxc_ij(:,:) = 0.0_dp

 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities'

 do iatom=1,natom
   if( rank_grid /= MODULO(iatom,nproc_grid) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),basis_element(iatom),coeff,alpha)


   call calculate_eri_approximate_hartree(basis,basis%nbf,basis%nbf,x(:,iatom),ngau,coeff,alpha,vhgau)
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
   do jbf=1,basis%nbf
     do ibf=1,basis%nbf 
       vhxc_ij(ibf,jbf) =  vhxc_ij(ibf,jbf) + weight &
           *  vxc * basis_function_r(ibf) * basis_function_r(jbf)
     enddo
   enddo

 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(exc)
 call xsum_grid(vhxc_ij)

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization
 write(stdout,  '(a,2(2x,f12.6))') '      XC energy (Ha):',exc

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()

 call stop_clock(timing_approx_ham)

end subroutine dft_approximate_vhxc


!=========================================================================
subroutine static_dipole(nstate,basis,occupation,c_matrix)
 use m_basis_set
 use m_atoms
 implicit none

 integer,intent(in)                 :: nstate
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                            :: istate,astate,iaspin
 integer                            :: mstate,pstate,mpspin
 integer                            :: ibf,jbf
 integer                            :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart,ibf_cart,jbf_cart
 integer                            :: iatom,idir
 real(dp)                           :: dipole(3)
 real(dp),allocatable               :: dipole_basis(:,:,:)
 real(dp),allocatable               :: dipole_cart(:,:,:)
 real(dp)                           :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer :: unitfile,var_i

! call start_clock(timing_spectrum)

 write(stdout,'(/,a)') ' Calculate the static dipole'


 !
 ! First precalculate all the needed dipole in the basis set
 !
 allocate(dipole_basis(3,basis%nbf,basis%nbf))
 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   jbf_cart = 1
   jbf      = 1
   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(dipole_cart(3,ni_cart,nj_cart))


     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_dipole(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),dipole_cart(:,i_cart,j_cart))
       enddo
     enddo

     do idir=1,3
       dipole_basis(idir,ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE( cart_to_pure(li)%matrix(:,:) ) , &
             MATMUL(  dipole_cart(idir,:,:) , cart_to_pure(lj)%matrix(:,:) ) )
     enddo

     deallocate(dipole_cart)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo


 call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)

 ! Minus sign for electrons
 do idir=1,3
   dipole(idir) = -SUM( dipole_basis(idir,:,:) * SUM( p_matrix(:,:,:) , DIM=3 ) )
 enddo

!****INTERVENTIONS****
print *, "dipoe_basis WRITING"
open(newunit=unitfile, file='dipole_basis.dat')
do idir=1,3
   do var_i=1, basis%nbf
      write(unitfile,*) dipole_basis(idir,var_i,:)
   enddo
enddo
close(unitfile)
!****INTERVENTIONS****

 deallocate(dipole_basis)

 do iatom=1,natom
   dipole(:) = dipole(:) + zatom(iatom) * x(:,iatom)
 enddo

 write(stdout,'(1x,a,3(2x,f14.6))') 'Dipole (a.u.):  ',dipole(:)
 write(stdout,'(1x,a,3(2x,f14.6))') 'Dipole (Debye): ',dipole(:) * au_debye




end subroutine static_dipole


!=========================================================================
subroutine static_quadrupole(nstate,basis,occupation,c_matrix)
 use m_basis_set
 use m_atoms
 implicit none

 integer,intent(in)                 :: nstate
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                            :: istate,astate,iaspin
 integer                            :: mstate,pstate,mpspin
 integer                            :: ibf,jbf
 integer                            :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart,ibf_cart,jbf_cart
 integer                            :: iatom,idir,jdir
 real(dp)                           :: trace
 real(dp)                           :: quad(3,3)
 real(dp),allocatable               :: quad_basis(:,:,:,:)
 real(dp),allocatable               :: quad_cart(:,:,:,:)
 real(dp)                           :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====


! call start_clock(timing_spectrum)

 write(stdout,'(/,a)') ' Calculate the static quadrupole'


 !
 ! First precalculate all the needed quadrupoles in the basis set
 !
 allocate(quad_basis(3,3,basis%nbf,basis%nbf))
 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   jbf_cart = 1
   jbf      = 1
   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(quad_cart(3,3,ni_cart,nj_cart))


     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_quadrupole(basis%bf(ibf_cart+i_cart-1),basis%bf(jbf_cart+j_cart-1),quad_cart(:,:,i_cart,j_cart))
       enddo
     enddo

     do jdir=1,3
       do idir=1,3
         quad_basis(idir,jdir,ibf:ibf+ni-1,jbf:jbf+nj-1) = MATMUL( TRANSPOSE( cart_to_pure(li)%matrix(:,:) ) , &
               MATMUL(  quad_cart(idir,jdir,:,:) , cart_to_pure(lj)%matrix(:,:) ) )
       enddo
     enddo

     deallocate(quad_cart)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo


 call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)

 ! Minus sign for electrons
 do jdir=1,3
   do idir=1,3
     quad(idir,jdir) = -SUM( quad_basis(idir,jdir,:,:) * SUM( p_matrix(:,:,:) , DIM=3 ) )
   enddo
 enddo

 deallocate(quad_basis)

 do iatom=1,natom
   do jdir=1,3
     quad(:,jdir) = quad(:,jdir) + zatom(iatom) * x(:,iatom) * x(jdir,iatom)
   enddo
 enddo

 write(stdout,'(1x,a,3(2x,f14.6))') 'Quadrupole (a.u.):    ',quad(1,:)
 write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(2,:)
 write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(3,:)
 write(stdout,*)
 write(stdout,'(1x,a,3(2x,f14.6))') 'Quadrupole (Debye.A): ',quad(1,:) * au_debye * bohr_A
 write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(2,:) * au_debye * bohr_A
 write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(3,:) * au_debye * bohr_A

 trace = quad(1,1) + quad(2,2) + quad(3,3)
 quad(1,1) = quad(1,1) - trace / 3.0_dp
 quad(2,2) = quad(2,2) - trace / 3.0_dp
 quad(3,3) = quad(3,3) - trace / 3.0_dp

 write(stdout,*)
 write(stdout,'(1x,a,3(2x,f14.6))') 'Traceless quadrupole (a.u.):    ',quad(1,:)
 write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(2,:)
 write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(3,:)
 write(stdout,*)
 write(stdout,'(1x,a,3(2x,f14.6))') 'Traceless quadrupole (Debye.A): ',quad(1,:) * au_debye * bohr_A
 write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(2,:) * au_debye * bohr_A
 write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(3,:) * au_debye * bohr_A



end subroutine static_quadrupole


end module m_hamiltonian
!=========================================================================
