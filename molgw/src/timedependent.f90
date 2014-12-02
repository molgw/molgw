!=========================================================================
#include "macros.h"
!=========================================================================


!=========================================================================
subroutine polarizability_td(basis,prod_basis,occupation,energy,c_matrix,wpol)
 use m_definitions
 use m_mpi
 use m_calculation_type,only: calculation_type,alpha_hybrid,alpha_hybrid_lr,dft_xc_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_eri
 use m_dft_grid
 use m_spectral_function
 use m_inputparam,only: calc_type,nspin,print_specfunc
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 use iso_c_binding,only: C_INT
 implicit none

 type(basis_set)                       :: basis,prod_basis
 real(dp),intent(in)                   :: occupation(basis%nbf,nspin)
 real(dp),intent(in)                   :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 type(spectral_function),intent(in)    :: wpol
!=====
 integer :: pbf,qbf,ibf,jbf,kbf,lbf,ijbf,klbf,ijbf_current,ijspin,klspin,ispin
 integer :: iorbital,jorbital,korbital,lorbital
 integer :: ipole
 integer :: t_ij,t_kl
 integer :: idft_xc,igrid
 integer :: reading_status

 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp),allocatable :: eri_eigenstate_k(:,:,:,:)
 real(dp)             :: spin_fact,alpha_local
 real(dp)             :: scissor_energy(nspin)
 real(dp)             :: h_2p(wpol%npole,wpol%npole)
 real(dp)             :: eigenvalue(wpol%npole),eigenvector(wpol%npole,wpol%npole),eigenvector_inv(wpol%npole,wpol%npole)
 real(dp)             :: matrix(wpol%npole,wpol%npole)
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
 real(dp),allocatable :: dipole_basis(:,:,:),dipole_state(:,:,:,:)
 real(dp),allocatable :: dipole_cart(:,:,:)
 real(dp)             :: dipole(3)
 real(dp),allocatable :: residu_left(:,:),residu_right(:,:)


 integer              :: ni,nj,li,lj,ni_cart,nj_cart,i_cart,j_cart,ibf_cart,jbf_cart
 logical              :: require_gradient
 logical              :: is_tddft
 character(len=256)   :: string
 logical,parameter    :: TDA=.FALSE.
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

 if( .NOT. calc_type%is_td .AND. .NOT. calc_type%is_bse) then
   stop'BUG: this should not happend in timedependent'
 endif

 spin_fact = REAL(-nspin+3,dp)
 is_tddft = calc_type%is_td .AND. calc_type%is_dft

 WRITE_MASTER(*,'(/,a)') ' calculating the polarizability for neutral excitation energies'

 if(TDA) then
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



 allocate(eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin))

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
       do iorbital=1,basis%nbf
         if( occupation(iorbital,ispin) > completely_empty/spin_fact ) then
           energy_qp(iorbital,ispin) = energy(iorbital,ispin)
         else
           energy_qp(iorbital,ispin) = energy(iorbital,ispin) + scissor_energy(ispin)
         endif
       enddo
     enddo
     WRITE_MASTER(*,'(/,a)') ' Scissor updated energies'
     do iorbital=1,basis%nbf
       WRITE_MASTER(*,'(i5,4(2x,f16.6))') iorbital,energy(iorbital,:)*Ha_eV,energy_qp(iorbital,:)*Ha_eV
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

 h_2p(:,:)=0.0_dp
 t_ij=0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

     call transform_eri_basis_lowmem(nspin,c_matrix,iorbital,ijspin,eri_eigenstate_i)

     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
       t_ij=t_ij+1


       t_kl=0
       do klspin=1,nspin
         do korbital=1,basis%nbf 

           do lorbital=1,basis%nbf 
             if( skip_transition(nspin,lorbital,korbital,occupation(lorbital,klspin),occupation(korbital,klspin)) ) cycle
             t_kl=t_kl+1

             if( TDA .AND. ( iorbital > jorbital ) .AND. ( korbital < lorbital )) cycle
             if( TDA .AND. ( iorbital < jorbital ) .AND. ( lorbital > korbital )) cycle

             !
             ! Hartree part
             h_2p(t_ij,t_kl) = eri_eigenstate_i(jorbital,korbital,lorbital,klspin) &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )

             !
             ! Add the kernel for TDDFT
             if(is_tddft) then
               h_2p(t_ij,t_kl) = h_2p(t_ij,t_kl)   &
                          + SUM(  wf_r(:,iorbital,ijspin) * wf_r(:,jorbital,ijspin) &
                                * wf_r(:,korbital,klspin) * wf_r(:,lorbital,klspin) &
                                * fxc(:,ijspin) )                                   & ! FIXME spin is not correct
                          * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )
             endif

             !
             ! Exchange part
             if(calc_type%need_exchange .OR. calc_type%is_bse) then
               if(ijspin==klspin) then
                 h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigenstate_i(korbital,jorbital,lorbital,klspin)  &
                          * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact * alpha_local
               endif
             endif

             !
             ! Screened exchange part
             if(calc_type%is_bse) then
               if(ijspin==klspin) then
                 kbf = prod_basis%index_prodbasis(iorbital,korbital)
                 bra(:) = wpol%residu_left (:,kbf+prod_basis%nbf*(ijspin-1))
                 kbf = prod_basis%index_prodbasis(jorbital,lorbital)
                 ket(:) = wpol%residu_right(:,kbf+prod_basis%nbf*(klspin-1))
                 h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  SUM( bra(:)*ket(:)/(-wpol%pole(:)) ) &
                          * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact   
               endif
             endif

           enddo
         enddo
       enddo !klspin


       h_2p(t_ij,t_ij) =  h_2p(t_ij,t_ij) + ( energy_qp(jorbital,ijspin) - energy_qp(iorbital,ijspin) )


     enddo !jorbital
   enddo !iorbital
 enddo ! ijspin

 if(is_tddft)    deallocate(fxc,wf_r)
 if(calc_type%is_bse) deallocate(bra,ket)

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'diago 2-particle hamiltonian'
 WRITE_MASTER(*,*) 'matrix',wpol%npole,'x',wpol%npole

 call start_clock(timing_diago_h2p)
 call diagonalize_general(wpol%npole,h_2p,eigenvalue,eigenvector)
 call stop_clock(timing_diago_h2p)
 WRITE_MASTER(*,*) 'diago finished'
 WRITE_MASTER(*,*)

! WRITE_MASTER(*,*) 'Neutral excitation energies [eV]'
! do t_ij=1,wpol%npole
!   WRITE_MASTER(*,'(1(i4,2x),20(2x,f12.6))') t_ij,eigenvalue(t_ij)*Ha_eV
! enddo

 WRITE_MASTER(*,'(/,a,f14.8)') ' Lowest neutral excitation energy [eV]',MINVAL(ABS(eigenvalue(:)))*Ha_eV
   
 call start_clock(timing_inversion_s2p)
 call invert(wpol%npole,eigenvector,eigenvector_inv)
 call stop_clock(timing_inversion_s2p)

 deallocate(eri_eigenstate_i)

 !
 ! Calculate Wp= v * chi * v 
 ! and then write it down on file
 !
 if( print_specfunc ) then
  call chi_to_vchiv(nspin,basis%nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
  call write_spectral_function(wpol)
 endif

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
 ! Get the dipole oscillator strength with orbitals
 allocate(dipole_state(3,basis%nbf,basis%nbf,nspin))
 dipole_state(:,:,:,:) = 0.0_dp
 do ijspin=1,nspin
   do jbf=1,basis%nbf
     do ibf=1,basis%nbf
       do jorbital=1,basis%nbf
         do iorbital=1,basis%nbf
           dipole_state(:,iorbital,jorbital,ijspin) = dipole_state(:,iorbital,jorbital,ijspin) &
                                         + c_matrix(ibf,iorbital,ijspin) * dipole_basis(:,ibf,jbf) * c_matrix(jbf,jorbital,ijspin)
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

 allocate(residu_left (3,wpol%npole))
 allocate(residu_right(3,wpol%npole))
 residu_left (:,:) = 0.0_dp
 residu_right(:,:) = 0.0_dp
 t_ij=0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin))) cycle
       t_ij=t_ij+1

       do t_kl=1,wpol%npole
         residu_left (:,t_kl)  = residu_left (:,t_kl) &
                      + dipole_state(:,iorbital,jorbital,ijspin) * eigenvector(t_ij,t_kl)
         residu_right(:,t_kl)  = residu_right(:,t_kl) &
                      + dipole_state(:,iorbital,Jorbital,ijspin) * eigenvector_inv(t_kl,t_ij) &
                                       * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )
       enddo

     enddo
   enddo
 enddo


 t_ij = 0
 absorp(:,:,:) = 0.0_dp
 static_polarizability(:,:) = 0.0_dp
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin))) cycle
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
 do t_ij=1,wpol%npole
   if(eigenvalue(t_ij) > 0.0_dp) then
     oscillator_strength = 2.0_dp/3.0_dp * DOT_PRODUCT(residu_left(:,t_ij),residu_right(:,t_ij)) * eigenvalue(t_ij)
     trk_sumrule = trk_sumrule + oscillator_strength
     WRITE_MASTER(*,'(i4,10(f18.8,2x))') t_ij,eigenvalue(t_ij)*Ha_eV,oscillator_strength
   endif
 enddo
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'TRK SUM RULE: the two following numbers should compare well'
 WRITE_MASTER(*,*) 'Sum over oscillator strengths',trk_sumrule
 WRITE_MASTER(*,*) 'Number of electrons          ',SUM( occupation(:,:) )

 WRITE_MASTER(*,'(/,a)') ' Static dipole polarizability'
 do idir=1,3
   WRITE_MASTER(*,'(3(4x,f12.6))') static_polarizability(idir,:)
 enddo

 do iomega=1,nomega
   WRITE_MASTER(101,'(10(e18.8,2x))') REAL(omega(iomega),dp)*Ha_eV,absorp(iomega,:,:)
 enddo 


 deallocate(residu_left,residu_right)
 deallocate(dipole_state)



end subroutine polarizability_td

