!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian_cmplx
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min


contains


!=========================================================================
subroutine setup_exchange_ri_cmplx(nbf,nstate,occupation,c_matrix_cmplx,p_matrix_cmplx, &
                                   exchange_ij_cmplx,eexchange)
 use m_eri
 implicit none
 integer,intent(in)         :: nbf,nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(out)       :: eexchange
 complex(dp),intent(in)    :: c_matrix_cmplx(nbf,nstate,nspin)
 complex(dp),intent(in)    :: p_matrix_cmplx(nbf,nbf,nspin)
 complex(dp),intent(out)   :: exchange_ij_cmplx(nbf,nbf,nspin)
!=====
 integer                    :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer                    :: index_ij
 integer                    :: nocc
 real(dp)                   :: eigval(nbf)
 integer                    :: ipair
 complex(dp),allocatable   :: tmp_cmplx(:,:)
!=====

! write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)

 exchange_ij_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )
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

     tmp_cmplx(:,:) = ( 0.0_dp, 0.0_dp )
     !$OMP PARALLEL PRIVATE(ibf,jbf) 
     !$OMP DO REDUCTION(+:tmp_cmplx)
     do ipair=1,npair
       ibf=index_basis(1,ipair)
       jbf=index_basis(2,ipair)
       tmp_cmplx(:,ibf) = tmp_cmplx(:,ibf) + c_matrix_cmplx(jbf,istate,ispin) * eri_3center(:,ipair)
       if( ibf /= jbf ) then
          tmp_cmplx(:,jbf) = tmp_cmplx(:,jbf) + c_matrix_cmplx(ibf,istate,ispin) * eri_3center(:,ipair)
       end if
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^H * A + C ; C - exchange_ij(:,:,ispin); A - tmp
     !    exchange_ij_cmplx(:,:,ispin) = exchange_ij_cmplx(:,:,ispin) - & 
     !            MATMUL( TRANSPOSE(tmp_cmplx(:,:)) , CONJG(tmp_cmplx(:,:)) ) * occupation(istate,ispin)/ spin_fact
     call ZHERK('L','C',nbf,nauxil_3center,-occupation(istate,ispin)/spin_fact,tmp_cmplx,nauxil_3center,1.0_dp,exchange_ij_cmplx(:,:,ispin),nbf)
       
   enddo
   exchange_ij_cmplx(:,:,ispin)=CONJG(exchange_ij_cmplx(:,:,ispin))
   !
   ! Need to make exchange_ij Hermitian (not just symmetric)
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij_cmplx(ibf,jbf,ispin) = conjg( exchange_ij_cmplx(jbf,ibf,ispin) )
     enddo
   enddo

 enddo ! end of loop do ispin=1,nspin
 deallocate(tmp_cmplx)
 ! This interface should work also for complex exchange_ij_cmplx 
 call xsum_world(exchange_ij_cmplx)
 eexchange = real( 0.5_dp * SUM( exchange_ij_cmplx(:,:,:) * conjg(p_matrix_cmplx(:,:,:)) ),dp)

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri_cmplx

!=========================================================================
subroutine setup_density_matrix_cmplx(nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)
 implicit none
 integer,intent(in)   :: nbf,nstate
 complex(dp),intent(in)  :: c_matrix_cmplx(nbf,nstate,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
 complex(dp),intent(out) :: p_matrix_cmplx(nbf,nbf,nspin)
!=====
 integer :: ispin,ibf,jbf
 integer :: istate

!=====

 call start_clock(timing_density_matrix)
! write(stdout,'(1x,a)') 'Build density matrix'

 p_matrix_cmplx(:,:,:) = ( 0.0_dp , 0.0_dp )
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     call ZHER('L',nbf,occupation(istate,ispin),c_matrix_cmplx(:,istate,ispin),1,p_matrix_cmplx(:,:,ispin),nbf)
   enddo

   ! Hermitianalize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix_cmplx(jbf,ibf,ispin) = conjg( p_matrix_cmplx(ibf,jbf,ispin) )
     enddo
   enddo
 enddo
 call stop_clock(timing_density_matrix)



end subroutine setup_density_matrix_cmplx

!=========================================================================
subroutine setup_density_matrix_cmplx_slow(nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)
 implicit none
 integer,intent(in)   :: nbf,nstate
 complex(dp),intent(in)  :: c_matrix_cmplx(nbf,nstate,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
 complex(dp),intent(out) :: p_matrix_cmplx(nbf,nbf,nspin)
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
subroutine dft_exc_vxc_cmplx(basis,nstate,occupation,c_matrix_cmplx,p_matrix,vxc_ij,exc_xc)
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
 real(dp),intent(in)        :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: vxc_ij(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: exc_xc
 complex(dp),intent(in)    :: c_matrix_cmplx(basis%nbf,nstate,nspin)
!=====

 real(dp),parameter :: TOL_RHO=1.0e-10_dp
 integer  :: idft_xc
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

! write(stdout,*) 'Calculate DFT XC potential'
 
 !
 ! If it is the first time, then set up the stored arrays
 !

 normalization(:)=0.0_dp

 do igrid=1,ngrid

   weight = w_grid(igrid)

   !
   ! Get the functions at point r
   call get_basis_functions_r(basis,igrid,basis_function_r)
   !
   ! calculate the density at point r for spin up and spin down
   call calc_density_r_cmplx(nspin,basis%nbf,nstate,occupation,c_matrix_cmplx,basis_function_r,rhor)

   ! Skip all the rest if the density is too small
   if( ALL( rhor(:) < TOL_RHO )  ) cycle

   !
   ! Get the gradient at point r
   if( dft_xc_needs_gradient ) then
     call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
   endif


   !
   ! Normalization
   normalization(:) = normalization(:) + rhor(:) * weight


   if( dft_xc_needs_gradient ) then 
     call calc_density_gradr_cmplx(nspin,basis%nbf,nstate,occupation,c_matrix_cmplx,basis_function_r,basis_function_gradr,grad_rhor)
   endif

   if( dft_xc_needs_gradient ) then
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
   if( .NOT. dft_xc_needs_gradient ) then 
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

! write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
! write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

 call stop_clock(timing_dft)

end subroutine dft_exc_vxc_cmplx

!=========================================================================
subroutine static_dipole_cmplx(nstate,basis,occupation,c_matrix_cmplx,dipole)
 use m_basis_set
 use m_atoms
 use m_timing
 implicit none

 integer,intent(in)                 :: nstate
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin)
 complex(dp),intent(in)             :: c_matrix_cmplx(basis%nbf,nstate,nspin)
 real(dp),intent(out)               :: dipole(3)
!=====
 integer                            :: gt
 integer                            :: istate,astate,iaspin
 integer                            :: mstate,pstate,mpspin
 integer                            :: ibf,jbf
 integer                            :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart,ibf_cart,jbf_cart
 integer                            :: iatom,idir
 real(dp),allocatable               :: dipole_basis(:,:,:)
 real(dp),allocatable               :: dipole_cart(:,:,:)
 complex(dp)                        :: p_matrix_cmplx(basis%nbf,basis%nbf,nspin)
!=====

! call start_clock(timing_spectrum)

! write(stdout,'(/,a)') ' Calculate the static dipole'

 gt = get_gaussian_type_tag(basis%gaussian_type)

 !
 ! First precalculate all the needed dipole in the basis set
 !
 allocate(dipole_basis(basis%nbf,basis%nbf,3))
 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bfc(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   jbf_cart = 1
   jbf      = 1
   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bfc(jbf_cart)%am
     nj_cart = number_basis_function_am('CART',lj)
     nj      = number_basis_function_am(basis%gaussian_type,lj)

     allocate(dipole_cart(ni_cart,nj_cart,3))


     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_dipole(basis%bfc(ibf_cart+i_cart-1),basis%bfc(jbf_cart+j_cart-1),dipole_cart(i_cart,j_cart,:))
       enddo
     enddo

     do idir=1,3
       dipole_basis(ibf:ibf+ni-1,jbf:jbf+nj-1,idir) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
             MATMUL(  dipole_cart(:,:,idir) , cart_to_pure(lj,gt)%matrix(:,:) ) )
     enddo

     deallocate(dipole_cart)

     jbf      = jbf      + nj
     jbf_cart = jbf_cart + nj_cart
   enddo

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo


 call setup_density_matrix_cmplx(basis%nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)

 ! Minus sign for electrons
 do idir=1,3
   dipole(idir) = real( -SUM( dipole_basis(:,:,idir) * SUM( p_matrix_cmplx(:,:,:) , DIM=3 ) ),dp)
 enddo

 deallocate(dipole_basis)

 do iatom=1,natom
   dipole(:) = dipole(:) + zatom(iatom) * x(:,iatom)
 enddo

end subroutine static_dipole_cmplx

!=========================================================================
subroutine static_dipole_fast_cmplx(basis,p_matrix_cmplx,dipole_basis,dipole)
 use m_basis_set
 use m_atoms
 use m_timing
 implicit none

!=====
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: dipole_basis(basis%nbf,basis%nbf,3)
 complex(dp),intent(in)      :: p_matrix_cmplx(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)        :: dipole(3)
!=====
 integer                     :: iatom,idir

 ! Minus sign for electrons
 do idir=1,3
   dipole(idir) = real( -SUM( dipole_basis(:,:,idir) * SUM( p_matrix_cmplx(:,:,:) , DIM=3 ) ),dp)
 enddo

 do iatom=1,natom
   dipole(:) = dipole(:) + zatom(iatom) * x(:,iatom)
 enddo

end subroutine static_dipole_fast_cmplx

end module m_hamiltonian_cmplx
!=========================================================================


