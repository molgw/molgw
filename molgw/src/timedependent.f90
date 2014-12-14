!=========================================================================
#include "macros.h"
!=========================================================================


!=========================================================================
subroutine polarizability_td(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,wpol)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_eri
 use m_dft_grid
 use m_spectral_function
 use m_inputparam
 use m_gw
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 use iso_c_binding,only: C_INT
 implicit none

 type(basis_set)                       :: basis,prod_basis,auxil_basis
 real(dp),intent(in)                   :: occupation(basis%nbf,nspin)
 real(dp),intent(in)                   :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 type(spectral_function),intent(inout) :: wpol
!=====
 integer              :: pbf,qbf,ibf,jbf,kbf,lbf,ijbf,klbf,ijbf_current,ijspin,klspin,ispin
 integer              :: istate,jstate,kstate,lstate
 integer              :: ipole
 integer              :: ntrans
 integer              :: t_ij,t_kl
 integer              :: idft_xc,igrid
 integer              :: reading_status

 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp)             :: eri_eigen_ijkl,eri_eigen_ikjl
 real(dp)             :: alpha_local,docc_ij,docc_kl
 real(dp)             :: scissor_energy(nspin)
 real(dp),allocatable :: h_2p(:,:)
 real(dp),allocatable :: eigenvalue(:),eigenvector(:,:),eigenvector_inv(:,:)
 real(dp)             :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp)             :: oscillator_strength,trk_sumrule
 real(dp)             :: energy_qp(basis%nbf,nspin)
 real(dp)             :: rr(3)
 real(dp)             :: basis_function_r       (basis%nbf)
! real(dp)             :: basis_function_r_shiftx(basis%nbf)
! real(dp)             :: basis_function_r_shifty(basis%nbf)
! real(dp)             :: basis_function_r_shiftz(basis%nbf)
 real(dp)             :: basis_function_gradr       (3,basis%nbf)
! real(dp)             :: basis_function_gradr_shiftx(3,basis%nbf)
! real(dp)             :: basis_function_gradr_shifty(3,basis%nbf)
! real(dp)             :: basis_function_gradr_shiftz(3,basis%nbf)
 real(dp)             :: rhor_r       (nspin)
! real(dp)             :: rhor_r_shiftx(nspin)
! real(dp)             :: rhor_r_shifty(nspin)
! real(dp)             :: rhor_r_shiftz(nspin)

 real(dp)             :: fxc_libxc(nspin)
 real(dp),allocatable :: fxc(:,:)
 real(dp),allocatable :: wf_r(:,:,:)
 real(dp),allocatable :: bra(:),ket(:)
 real(dp),allocatable :: bra_auxil(:,:,:,:),ket_auxil(:,:,:,:)
 real(dp),allocatable :: dipole_basis(:,:,:),dipole_state(:,:,:,:)
 real(dp),allocatable :: dipole_cart(:,:,:)
 real(dp),allocatable :: residu_left(:,:),residu_right(:,:)


 integer              :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart,ibf_cart,jbf_cart
 logical              :: require_gradient
 logical              :: is_tddft
 character(len=256)   :: string
 integer              :: iomega,idir,jdir
 integer,parameter    :: nomega=600
 complex(dp)          :: omega(nomega)
 real(dp)             :: absorp(nomega,3,3)
 real(dp)             :: static_polarizability(3,3)

#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func(ndft_xc),xc_functest
 type(xc_f90_pointer_t) :: xc_info(ndft_xc),xc_infotest
#endif
!=====

 call start_clock(timing_pola)

 if( .NOT. calc_type%is_td .AND. .NOT. calc_type%is_bse) then
   stop'BUG: this should not happend in timedependent'
 endif

 is_tddft = calc_type%is_td .AND. calc_type%is_dft

 WRITE_MASTER(*,'(/,a)') ' Calculating the polarizability for neutral excitation energies'

 ! Obtain the number of transition = the size of the matrix
 call init_spectral_function(basis%nbf,occupation,ntrans)

 allocate(h_2p(ntrans,ntrans))
 allocate(eigenvector(ntrans,ntrans))
 allocate(eigenvalue(ntrans))

 if(calc_type%is_tda) then
   msg='Tamm-Dancoff approximation is switched on'
   call issue_warning(msg)
 endif

 if(calc_type%is_td) then
   alpha_local = alpha_hybrid
 else
   alpha_local = 1.0_dp
 endif

 if(is_tddft) then
   !
   ! Prepare DFT kernel calculation with Libxc
   !
   do idft_xc=1,ndft_xc
     if(nspin==1) then
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), INT(dft_xc_type(idft_xc),C_INT), XC_UNPOLARIZED)
     else
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), INT(dft_xc_type(idft_xc),C_INT), XC_POLARIZED)
     endif
     if( MODULO(xc_f90_info_flags( xc_info(idft_xc)),XC_FLAGS_HAVE_FXC*2) < XC_FLAGS_HAVE_FXC ) then
       stop'This functional does not have the kernel implemented in Libxc'
     endif
     call xc_f90_info_name(xc_info(idft_xc),string)
     WRITE_MASTER(*,'(a,i4,a,i6,5x,a)') '   XC functional ',idft_xc,' :  ',xc_f90_info_number(xc_info(idft_xc)),&
           TRIM(string)
     if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_GGA     ) require_gradient  =.TRUE.
     if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) require_gradient  =.TRUE.
   enddo

   !
   ! calculate rho, grad rho and the kernel
   ! 
   ! Get the density matrix P from C
   call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)

   allocate(fxc(ngrid,nspin),wf_r(ngrid,basis%nbf,nspin))

   fxc(:,:) = 0.0_dp

   do igrid=1,ngrid

     rr(:) = rr_grid(:,igrid)
  
     !
     ! Get all the functions and gradients at point rr
     call get_basis_functions_r(basis,igrid,basis_function_r)
     !
     ! store the wavefunction in r
     do ispin=1,nspin
       wf_r(igrid,:,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,:,ispin) )
     enddo

!     if( require_gradient ) then
!       call get_basis_functions_gradr(basis,igrid,&
!                                      basis_function_gradr,basis_function_r_shiftx,basis_function_r_shifty,basis_function_r_shiftz,&
!                                      basis_function_gradr_shiftx,basis_function_gradr_shifty,basis_function_gradr_shiftz)
!     endif

     if( .NOT. allocated(bfr) ) stop'bfr not allocated -> weird'

     call calc_density_r(nspin,basis%nbf,p_matrix,basis_function_r,rhor_r)

     !
     ! Calculate the kernel
     ! 
     do idft_xc=1,ndft_xc
       select case(xc_f90_info_family(xc_info(idft_xc)))

       case(XC_FAMILY_LDA)
         call xc_f90_lda_fxc(xc_func(idft_xc),1_C_INT,rhor_r(1),fxc_libxc(1))
       case default
         stop'GGA kernel not yet implemented'
       end select

       ! Store the result with the weight
       ! Remove too large values for stability
       fxc(igrid,:) = fxc(igrid,:) + MIN(fxc_libxc(:),1.0E8_dp) * w_grid(igrid) * dft_xc_coef(idft_xc)

     enddo



   enddo
 endif



 if(.NOT. is_auxil_basis) then
   allocate(eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin))
 endif

 !
 ! Prepare the BSE calculation
 if( calc_type%is_bse ) then

   allocate(bra(wpol%npole),ket(wpol%npole))
   call read_energy_qp(nspin,basis%nbf,energy_qp,reading_status)
   select case(reading_status)
   case(-1)
     scissor_energy(:) = energy_qp(1,:)
     WRITE_MASTER(*,'(a,2(x,f12.6))') ' Scissor operator with value [eV]:',scissor_energy(:)*Ha_eV
     do ispin=1,nspin
       do istate=1,basis%nbf
         if( occupation(istate,ispin) > completely_empty/spin_fact ) then
           energy_qp(istate,ispin) = energy(istate,ispin)
         else
           energy_qp(istate,ispin) = energy(istate,ispin) + scissor_energy(ispin)
         endif
       enddo
     enddo
     WRITE_MASTER(*,'(/,a)') ' Scissor updated energies'
     do istate=1,basis%nbf
       WRITE_MASTER(*,'(i5,4(2x,f16.6))') istate,energy(istate,:)*Ha_eV,energy_qp(istate,:)*Ha_eV
     enddo
     WRITE_MASTER(*,*)
   case(0)
     WRITE_MASTER(*,'(a)') ' Reading OK'
   case(1,2)
     WRITE_MASTER(*,'(a,/,a)') ' Something happened during the reading of energy_qp file',' Fill up the QP energies with KS energies'
     energy_qp(:,:) = energy(:,:)
   case default
     stop'reading_status BUG'
   end select

 else
   !
   ! For any other type of calculation, just fill energy_qp array with energy
   energy_qp(:,:) = energy(:,:)
 endif

 call start_clock(timing_build_h2p)
 ! Prepare the bra and ket for BSE
 if(is_auxil_basis .AND. calc_type%is_bse) then
   allocate(bra_auxil(wpol%npole,ncore_W+1:nvirtual_W-1,ncore_W+1:nvirtual_W-1,nspin))
   allocate(ket_auxil(wpol%npole,ncore_W+1:nvirtual_W-1,ncore_W+1:nvirtual_W-1,nspin))
   do ijspin=1,nspin
     do istate=ncore_W+1,nvirtual_W-1 
       do jstate=ncore_W+1,nvirtual_W-1

         ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
         do ipole=1,wpol%npole
           bra_auxil(ipole,istate,jstate,ijspin) = DOT_PRODUCT( wpol%residu_left (ipole,:) , eri_3center_eigen(:,istate,jstate,ijspin) )
           ket_auxil(ipole,istate,jstate,ijspin) = DOT_PRODUCT( wpol%residu_right(ipole,:) , eri_3center_eigen(:,istate,jstate,ijspin) )
         enddo

       enddo
     enddo
   enddo
 endif

 h_2p(:,:)=0.0_dp
 t_ij=0
 do ijspin=1,nspin
   do istate=1,basis%nbf ! istate stands for occupied or partially occupied

     if( .NOT. is_auxil_basis) call transform_eri_basis(nspin,c_matrix,istate,ijspin,eri_eigenstate_i)

     do jstate=1,basis%nbf ! jstate stands for empty or partially empty
       if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       t_ij=t_ij+1
       docc_ij = occupation(istate,ijspin)-occupation(jstate,ijspin)


       t_kl=0
       do klspin=1,nspin
         do kstate=1,basis%nbf 

           do lstate=1,basis%nbf 
             if( skip_transition(nspin,lstate,kstate,occupation(lstate,klspin),occupation(kstate,klspin)) ) cycle
             t_kl=t_kl+1

             docc_kl = occupation(kstate,klspin)-occupation(lstate,klspin)

!             if( calc_type%is_tda .AND. ( istate > jstate ) .AND. ( kstate < lstate )) cycle
!             if( calc_type%is_tda .AND. ( istate < jstate ) .AND. ( lstate > kstate )) cycle

             if( calc_type%is_tda .AND. docc_ij * docc_kl < 0.0_dp ) cycle

             if(is_auxil_basis) then
               eri_eigen_ijkl = eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
             else
               eri_eigen_ijkl = eri_eigenstate_i(jstate,kstate,lstate,klspin)
             endif

             !
             ! Hartree part
             h_2p(t_ij,t_kl) = h_2p(t_ij,t_kl) + eri_eigen_ijkl * docc_ij

             !
             ! Add the kernel for TDDFT
             if(is_tddft) then
               h_2p(t_ij,t_kl) = h_2p(t_ij,t_kl)   &
                          + SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                                * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                                * fxc(:,ijspin) )                               & ! FIXME spin is not correct
                                * docc_ij
             endif

             !
             ! Exchange part
             if(calc_type%need_exchange .OR. calc_type%is_bse) then
               if(ijspin==klspin) then

                 if(is_auxil_basis) then
                   eri_eigen_ikjl = eri_eigen_ri(istate,kstate,ijspin,jstate,lstate,klspin)
                 else
                   eri_eigen_ikjl = eri_eigenstate_i(kstate,jstate,lstate,klspin)
                 endif

                 h_2p(t_ij,t_kl) = h_2p(t_ij,t_kl) -  eri_eigen_ikjl  &
                                   * docc_ij / spin_fact * alpha_local
               endif
             endif

             !
             ! Screened exchange part
             if(calc_type%is_bse) then
               if(ijspin==klspin) then

                 if(.NOT. is_auxil_basis) then
                   ! Here just grab the precalculated value
                   kbf = prod_basis%index_prodbasis(istate,kstate)+prod_basis%nbf*(ijspin-1)
                   bra(:) = wpol%residu_left (:,kbf)
                   kbf = prod_basis%index_prodbasis(jstate,lstate)+prod_basis%nbf*(klspin-1)
                   ket(:) = wpol%residu_right(:,kbf)
                 else
                   ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
                   bra(:) = bra_auxil(:,istate,kstate,ijspin) 
                   ket(:) = ket_auxil(:,jstate,lstate,klspin)
                 endif

                 h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  SUM( bra(:)*ket(:)/(-wpol%pole(:)) ) &
                                   * docc_ij / spin_fact   
               endif
             endif

           enddo
         enddo
       enddo !klspin


       h_2p(t_ij,t_ij) =  h_2p(t_ij,t_ij) + ( energy_qp(jstate,ijspin) - energy_qp(istate,ijspin) )


     enddo !jstate
   enddo !istate
 enddo ! ijspin

 call destroy_spectral_function(wpol)

 if( .NOT. is_auxil_basis) deallocate(eri_eigenstate_i)
 if(is_tddft)    deallocate(fxc,wf_r)
 if(calc_type%is_bse) deallocate(bra,ket)
 if(allocated(bra_auxil)) deallocate(bra_auxil)
 if(allocated(ket_auxil)) deallocate(ket_auxil)

 call stop_clock(timing_build_h2p)

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'diago 2-particle hamiltonian'
 WRITE_MASTER(*,*) 'matrix',ntrans,'x',ntrans

 call start_clock(timing_diago_h2p)
 call diagonalize_general(ntrans,h_2p,eigenvalue,eigenvector)
 call stop_clock(timing_diago_h2p)
 WRITE_MASTER(*,*) 'diago finished'
 WRITE_MASTER(*,*)
 deallocate(h_2p)
 allocate(eigenvector_inv(ntrans,ntrans))


 WRITE_MASTER(*,'(/,a,f14.8)') ' Lowest neutral excitation energy [eV]',MINVAL(ABS(eigenvalue(:)))*Ha_eV
   
 call start_clock(timing_inversion_s2p)
 call invert(ntrans,eigenvector,eigenvector_inv)
 call stop_clock(timing_inversion_s2p)

 !
 ! Calculate the spectrum now
 !

 WRITE_MASTER(*,'(/,a)') ' Calculate the optical spectrum'

 if (nspin/=1) stop'no nspin/=1 allowed'


 !
 ! First precalculate all the needed dipole in the basis set
 !
 allocate(dipole_basis(3,basis%nbf,basis%nbf))
 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am(CARTESIAN,li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   jbf_cart = 1
   jbf      = 1
   do while(jbf_cart<=basis%nbf_cart)
     lj      = basis%bf(jbf_cart)%am
     nj_cart = number_basis_function_am(CARTESIAN,lj)
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

 !
 ! Get the dipole oscillator strength with states
 allocate(dipole_state(3,basis%nbf,basis%nbf,nspin))
 dipole_state(:,:,:,:) = 0.0_dp
 do ijspin=1,nspin
   do jbf=1,basis%nbf
     do ibf=1,basis%nbf
       do jstate=1,basis%nbf
         do istate=1,basis%nbf
           dipole_state(:,istate,jstate,ijspin) = dipole_state(:,istate,jstate,ijspin) &
                                         + c_matrix(ibf,istate,ijspin) * dipole_basis(:,ibf,jbf) * c_matrix(jbf,jstate,ijspin)
         enddo
       enddo
     enddo
   enddo
 enddo

 !
 ! set the frequency mesh
 omega(1)     =MAX( 0.0_dp      ,MINVAL(ABS(eigenvalue(:)))-3.00/Ha_eV)
 omega(nomega)=MIN(20.0_dp/Ha_eV,MAXVAL(ABS(eigenvalue(:)))+3.00/Ha_eV)
 do iomega=2,nomega-1
   omega(iomega) = omega(1) + ( omega(nomega)-omega(1) ) /REAL(nomega-1,dp) * (iomega-1) 
 enddo
 ! add the broadening
 omega(:) = omega(:) + im * 0.10/Ha_eV

 allocate(residu_left (3,ntrans))
 allocate(residu_right(3,ntrans))
 residu_left (:,:) = 0.0_dp
 residu_right(:,:) = 0.0_dp
 t_ij=0
 do ijspin=1,nspin
   do istate=1,basis%nbf ! istate stands for occupied or partially occupied

     do jstate=1,basis%nbf ! jstate stands for empty or partially empty
       if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin))) cycle
       t_ij=t_ij+1
       docc_ij = occupation(istate,ijspin)-occupation(jstate,ijspin)

       do t_kl=1,ntrans
         residu_left (:,t_kl)  = residu_left (:,t_kl) &
                      + dipole_state(:,istate,jstate,ijspin) * eigenvector(t_ij,t_kl)
         residu_right(:,t_kl)  = residu_right(:,t_kl) &
                      + dipole_state(:,istate,Jstate,ijspin) * eigenvector_inv(t_kl,t_ij) &
                                       * docc_ij
       enddo

     enddo
   enddo
 enddo


 t_ij = 0
 absorp(:,:,:) = 0.0_dp
 static_polarizability(:,:) = 0.0_dp
 do ijspin=1,nspin
   do istate=1,basis%nbf ! istate stands for occupied or partially occupied

     do jstate=1,basis%nbf ! jstate stands for empty or partially empty
       if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin))) cycle
       t_ij = t_ij + 1

       do idir=1,3
         do jdir=1,3
           absorp(:,idir,jdir) = absorp(:,idir,jdir) &
                                       + residu_left(idir,t_ij) * residu_right(jdir,t_ij) &
                                         *AIMAG( -1.0_dp  / ( omega(:) - eigenvalue(t_ij) ) )
           static_polarizability(idir,jdir) = static_polarizability(idir,jdir) + residu_left(idir,t_ij) * residu_right(jdir,t_ij) / eigenvalue(t_ij)
         enddo
       enddo

     enddo
   enddo
 enddo

 WRITE_MASTER(*,'(/,a)') ' Neutral excitation energies [eV] and strengths'
 trk_sumrule=0.0_dp
 t_kl=0
 do t_ij=1,ntrans
   if(eigenvalue(t_ij) > 0.0_dp) then
     t_kl = t_kl + 1
     oscillator_strength = 2.0_dp/3.0_dp * DOT_PRODUCT(residu_left(:,t_ij),residu_right(:,t_ij)) * eigenvalue(t_ij)
     trk_sumrule = trk_sumrule + oscillator_strength
     WRITE_MASTER(*,'(i4,10(f18.8,2x))') t_kl,eigenvalue(t_ij)*Ha_eV,oscillator_strength
   endif
 enddo
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'TRK SUM RULE: the two following numbers should compare well'
 WRITE_MASTER(*,*) 'Sum over oscillator strengths',trk_sumrule
 WRITE_MASTER(*,*) 'Number of valence electrons  ',SUM( occupation(ncore_W+1:,:) )

 WRITE_MASTER(*,'(/,a)') ' Static dipole polarizability'
 do idir=1,3
   WRITE_MASTER(*,'(3(4x,f12.6))') static_polarizability(idir,:)
 enddo

 do iomega=1,nomega
   WRITE_MASTER(101,'(10(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,absorp(iomega,:,:)
 enddo 


 deallocate(residu_left,residu_right)
 deallocate(dipole_state)


 !
 ! Calculate Wp= v * chi * v 
 ! and then write it down on file
 !
 if( print_specfunc ) then
  if( is_auxil_basis) then
    call chi_to_sqrtvchisqrtv_auxil(ntrans,basis%nbf,auxil_basis%nbf,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
  else
    call chi_to_vchiv(ntrans,basis%nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
  endif
  call write_spectral_function(wpol)
 endif

 if(allocated(eigenvector))     deallocate(eigenvector)
 if(allocated(eigenvector_inv)) deallocate(eigenvector_inv)
 if(allocated(eigenvalue))      deallocate(eigenvalue)

 call stop_clock(timing_pola)

end subroutine polarizability_td

