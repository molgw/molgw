!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! basic tools to play with the calculated self-energies
!
!=========================================================================
module m_selfenergy_tools
 use m_definitions
 use m_warning
 use m_mpi
 use m_inputparam
 use m_numerical_tools
 use m_hamiltonian_tools
 use m_atoms
 use m_inputparam
 use m_basis_set
 use m_dft_grid
 use m_hamiltonian_onebody
 use m_hamiltonian_wrapper
 use m_lbfgs

 !
 ! frozen core approximation parameters
 integer,protected :: ncore_G
 integer,protected :: nvirtual_G

 !
 ! Range of states to evaluate the self-energy
 integer,protected :: nsemin
 integer,protected :: nsemax

 !
 ! Highest occupied state
 integer,protected :: nhomo_G

 !
 ! Selfenergy evaluated on a frequency grid
 type selfenergy_grid
   integer                 :: nomega
   integer                 :: nomega_calc
   complex(dp),allocatable :: omega(:)
   complex(dp),allocatable :: omega_calc(:)
   real(dp),allocatable    :: weight_calc(:)
   real(dp),allocatable    :: energy0(:,:)
   complex(dp),allocatable :: sigma(:,:,:)
   complex(dp),allocatable :: sigma_calc(:,:,:)
 end type


contains


!=========================================================================
subroutine selfenergy_set_state_range(nstate_in,occupation)
 implicit none
 integer             :: nstate_in
 real(dp),intent(in) :: occupation(:,:)
!=====
 integer :: pstate
!=====

 if( nstate_in > SIZE( occupation(:,:) , DIM=1 ) ) then
   call die('selfenergy_set_state_range: nstate is too large')
 endif

 ncore_G      = ncoreg
 nvirtual_G   = MIN(nvirtualg,nstate_in+1)

 if(is_frozencore) then
   if( ncore_G == 0) ncore_G = atoms_core_states()
 endif

 if( ncore_G > 0 ) then
   write(msg,'(a,i4,2x,i4)') 'frozen core approximation in G switched on up to state = ',ncore_G
   call issue_warning(msg)
 endif

 if( nvirtual_G <= nstate_in ) then
   write(msg,'(a,i4,2x,i4)') 'frozen virtual approximation in G switched on starting with state = ',nvirtual_G
   call issue_warning(msg)
 endif

 ! Find the HOMO index
 nhomo_G = 1
 do pstate=1,nstate_in
   if( .NOT. ANY( occupation(pstate,:) < completely_empty ) ) then
     nhomo_G = MAX(nhomo_G,pstate)
   endif
 enddo

 nsemin = MAX(ncore_G+1,selfenergy_state_min,1,nhomo_G-selfenergy_state_range)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate_in,nhomo_G+selfenergy_state_range)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax

end subroutine selfenergy_set_state_range


!=========================================================================
subroutine write_selfenergy_omega(filename_root,exchange_m_vxc,occupation,energy0,se)
 implicit none

 character(len=*)    :: filename_root
 real(dp),intent(in) :: exchange_m_vxc(:,:)
 real(dp),intent(in) :: occupation(:,:),energy0(:,:)
 type(selfenergy_grid),intent(in) :: se
!=====
 integer            :: nstate
 character(len=3)   :: ctmp
 character(len=256) :: filename
 integer            :: selfenergyfile,selfenergyfile_cmplx
 integer            :: pstate,pspin
 integer            :: iomega
 real(dp)           :: spectral_function_w(nspin),sign_occ(nspin)
!=====

 ! Just the master writes
 if( .NOT. is_iomaster ) return

 nstate = SIZE(exchange_m_vxc,DIM=1)

 write(stdout,'(/,1x,a)') 'Write Sigma(omega) on file'

 !
 ! omega is defined with respect to energy0_a
 ! Absolute omega is omega + energy0_a
 !
 do pstate=nsemin,nsemax
   write(ctmp,'(i3.3)') pstate
   filename = TRIM(filename_root) // '_state' // TRIM(ctmp) // '.dat'
   write(stdout,'(1x,a,a)') 'Writing selfenergy in file: ', TRIM(filename)
   open(newunit=selfenergyfile,file=filename)

   write(selfenergyfile,'(a)') &
    '#       omega (eV)          Re SigmaC (eV)     Im SigmaC (eV)    omega - e_gKS - Vxc + SigmaX (eV)     A (eV^-1)'

   do pspin=1,nspin
     sign_occ(:) = SIGN( 1.0_dp , occupation(pstate,pspin) - spin_fact * 0.50_dp )
   enddo

   do iomega=-se%nomega,se%nomega
     spectral_function_w(:) = 1.0_dp / pi * ABS(   &
                                       AIMAG( 1.0_dp   &
                                         / ( se%energy0(pstate,:) + se%omega(iomega) - energy0(pstate,:) + ieta * sign_occ(:) &
                                               - exchange_m_vxc(pstate,:) - se%sigma(iomega,pstate,:) ) ) )

     write(selfenergyfile,'(20(f16.8,2x))') ( se%omega(iomega) + se%energy0(pstate,:) )*Ha_eV,           &
                                            REAL(se%sigma(iomega,pstate,:),dp) * Ha_eV,                  &
                                            AIMAG(se%sigma(iomega,pstate,:)) * Ha_eV,                    &
                                            ( REAL(se%omega(iomega),dp) + se%energy0(pstate,:)           &
                                             - energy0(pstate,:) - exchange_m_vxc(pstate,:) ) * Ha_eV,   &
                                            spectral_function_w(:) / Ha_eV
   enddo
   write(selfenergyfile,*)
   close(selfenergyfile)

   if( se%nomega_calc > 0 ) then
     filename = TRIM(filename_root) // '_cmplx_state' // TRIM(ctmp) // '.dat'
     write(stdout,'(1x,a,a)') 'Writing selfenergy for complex frequencies in file: ', TRIM(filename)
     open(newunit=selfenergyfile_cmplx,file=filename)
     write(selfenergyfile_cmplx,'(a)') &
      '#       omega (eV)          Re SigmaC (eV)     Im SigmaC (eV)    omega - e_gKS - Vxc + SigmaX (eV)     A (eV^-1)'
     do iomega=1,se%nomega_calc
       write(selfenergyfile_cmplx,'(20(f16.8,2x))') ( se%omega_calc(iomega) + se%energy0(pstate,:) )*Ha_eV,   &
                                                    se%sigma_calc(iomega,pstate,:) * Ha_eV,                   &
                                                    0.0_dp,0.0_dp
     enddo
     close(selfenergyfile_cmplx)
   endif

 enddo


end subroutine write_selfenergy_omega


!=========================================================================
subroutine find_qp_energy_linearization(se,exchange_m_vxc,energy0,energy_qp_z,zz)
 implicit none

 type(selfenergy_grid),intent(in) :: se
 real(dp),intent(in)              :: exchange_m_vxc(:,:),energy0(:,:)
 real(dp),intent(out)             :: energy_qp_z(:,:)
 real(dp),intent(out),optional    :: zz(:,:)
!=====
 integer  :: nstate
 integer  :: pstate,pspin
 real(dp) :: zz_p(nspin)
!=====

 nstate = SIZE(exchange_m_vxc,DIM=1)

 ! First, a dummy initialization
 energy_qp_z(:,:) = energy0(:,:)

 ! Then overwrite the interesting energy with the calculated GW one
 !$OMP PARALLEL
 !$OMP DO PRIVATE(pspin,zz_p)
 do pstate=nsemin,nsemax

   if( se%nomega > 0 .AND. PRESENT(zz) ) then
     zz_p(:) = REAL( se%sigma(1,pstate,:) - se%sigma(-1,pstate,:) ,dp) / REAL( se%omega(1) - se%omega(-1) ,dp)
     zz_p(:) = 1.0_dp / ( 1.0_dp - zz_p(:) )
     ! Constrain Z to be in [0:1] to avoid crazy values
     ! Z falls out of [0:1] when a weak self-energy pole is very close to Eks.
     ! Z out of [0:1] is an indicator for whether it happened or not.
     do pspin=1,nspin
       zz_p(pspin) = MIN( MAX(zz_p(pspin),0.0_dp) , 1.0_dp )
     enddo

     zz(pstate,:)          = zz_p(:)
     energy_qp_z(pstate,:) = se%energy0(pstate,:)  &
             + zz_p(:) * ( energy0(pstate,:) - se%energy0(pstate,:)  &
                            + REAL(se%sigma(0,pstate,:),dp) + exchange_m_vxc(pstate,:) )

   else

     energy_qp_z(pstate,:) = energy0(pstate,:) + REAL(se%sigma(0,pstate,:),dp) + exchange_m_vxc(pstate,:)

   endif

 enddo
 !$OMP END DO
 !$OMP END PARALLEL

end subroutine find_qp_energy_linearization


!=========================================================================
subroutine find_qp_energy_graphical(se,exchange_m_vxc,energy0,energy_qp_g)
 implicit none

 type(selfenergy_grid),intent(in) :: se
 real(dp),intent(in)  :: exchange_m_vxc(:,:),energy0(:,:)
 real(dp),intent(out) :: energy_qp_g(:,:)
!=====
 integer,parameter :: NFIXED_POINTS_MAX = 4
 integer  :: nfixed,ifixed
 integer  :: nstate,pstate,pspin
 real(dp) :: energy_fixed_point(NFIXED_POINTS_MAX,nsemin:nsemax,nspin),z_weight(NFIXED_POINTS_MAX,nsemin:nsemax,nspin)
 real(dp) :: equation_lhs(-se%nomega:se%nomega),equation_rhs(-se%nomega:se%nomega)
!=====

 nstate = SIZE(exchange_m_vxc,DIM=1)


 if( se%nomega > 0 ) then

   write(stdout,'(/,1x,a)') 'Graphical solution to the QP equation'
   write(stdout,'(1x,a,sp,f8.3,a,f8.3)') 'Scanning range around the input energy (eV): ', &
                                REAL(se%omega(-se%nomega),dp)*Ha_eV,'  --',REAL(se%omega(se%nomega),dp)*Ha_eV
   write(stdout,'(1x,a,i6)')             'Number of discretization points:             ',SIZE(se%omega(:))

   ! First, a dummy initialization
   energy_qp_g(:,:) = 0.0_dp
   z_weight(:,:,:)  = 0.0_dp
   energy_fixed_point(:,:,:) = 0.0_dp

   ! Then overwrite the interesting energy with the calculated GW one
   do pstate=nsemin,nsemax

     if( MODULO(pstate-nsemin,world%nproc) /= world%rank ) cycle

     do pspin=1,nspin
       !
       ! QP equation:
       ! E_GW = E0 + \omega =  E_gKS + \Sigma_c(E0+\omega) + \Sigma_x - v_xc
       !
       equation_lhs(:) = REAL(se%omega(:),dp)+se%energy0(pstate,pspin)
       equation_rhs(:) = REAL(se%sigma(:,pstate,pspin),dp) + exchange_m_vxc(pstate,pspin) + energy0(pstate,pspin)
       call find_fixed_point(equation_lhs,equation_rhs,energy_fixed_point(:,pstate,pspin),z_weight(:,pstate,pspin))

       ! If no reasonable QP energy found, then set E_qp to E_gKS for safety
       if( z_weight(1,pstate,pspin) < 1.0e-6_dp ) then
         energy_qp_g(pstate,pspin) = energy0(pstate,pspin)
       else
         energy_qp_g(pstate,pspin) = energy_fixed_point(1,pstate,pspin)
       endif

     enddo

   enddo

   call world%sum(energy_qp_g)
   call world%sum(z_weight)
   call world%sum(energy_fixed_point)

   ! Master IO node outputs the solution details
   write(stdout,'(/,1x,a)') 'state spin    QP energy (eV)  QP spectral weight'
   do pstate=nsemin,nsemax
     do pspin=1,nspin
       nfixed = COUNT( z_weight(:,pstate,pspin) > 0.0_dp )
       if( nfixed > 0 ) then
         write(stdout,'(1x,i5,2x,i3,*(3x,f12.6,3x,f12.6))') pstate,pspin, &
                   ( energy_fixed_point(ifixed,pstate,pspin)*Ha_eV,z_weight(ifixed,pstate,pspin), ifixed = 1,nfixed)
       else
         write(stdout,'(1x,i5,2x,i3,a)') pstate,pspin,'      has no graphical solution in the calculated range'
       endif
     enddo
   enddo

   if( ANY(z_weight(1,:,:) < 0.0_dp) ) then
     call issue_warning('At least one state had no graphical solution in the calculated range.'  &
                     // ' Increase nomega_sigma or step_sigma.')
   endif

 else

   do pstate=nsemin,nsemax
     energy_qp_g(pstate,:) = energy0(pstate,:) + REAL(se%sigma(0,pstate,:),dp) + exchange_m_vxc(pstate,:)
   enddo

 endif

 ! Unchanged energies for the states that were not calculated (outside range [nsemin:nsemax])
 energy_qp_g(:nsemin-1,:) = energy0(:nsemin-1,:)
 energy_qp_g(nsemax+1:,:) = energy0(nsemax+1:,:)



end subroutine find_qp_energy_graphical


!=========================================================================
subroutine find_fixed_point(xx,fx,energy_fixed_point,z_weight)
 implicit none
 real(dp),intent(in) :: xx(:)
 real(dp),intent(in) :: fx(:)
 real(dp),intent(out) :: energy_fixed_point(:),z_weight(:)
!=====
 integer  :: nx,ix
 integer  :: nfpmx,ifixed,jfixed
 real(dp) :: gx(SIZE(xx))
 real(dp) :: gpxi
 real(dp) :: z_zero
 integer  :: enumerate(SIZE(energy_fixed_point))
!=====

 nx = SIZE(xx)
 nfpmx = SIZE(energy_fixed_point)
 do ifixed=1,nfpmx
   enumerate(ifixed) = ifixed
 enddo

 ! Negative value to indicate something bad happened
 z_weight(:)           = -1.0_dp
 energy_fixed_point(:) =  0.0_dp

 !
 ! g(x) contains f(x) - x
 gx(:) = fx(:) - xx(:)
 ifixed = 0

 do ix=1,nx-1
   ! Check for sign changes => fixed point
   if( gx(ix) * gx(ix+1) < 0.0_dp ) then
     ! Evaluate the weight Z of the pole
     !  Z = ( 1 - d\Sigma / d\omega )^-1
     !  Z should be in [0,1]
     z_zero = 1.0_dp / ( 1.0_dp - ( fx(ix+1) - fx(ix) ) / ( xx(ix+1) - xx(ix) ) )

     ! Z not valid
     if( z_zero < 0.00001_dp ) cycle
     ! Compare the new Z with the largest Z's found at this stage
     if( z_zero > z_weight(nfpmx) ) then

       jfixed = MINVAL( enumerate(:), z_zero > z_weight(:) )
       z_weight(jfixed+1:nfpmx)           = z_weight(jfixed:nfpmx-1)
       energy_fixed_point(jfixed+1:nfpmx) = energy_fixed_point(jfixed:nfpmx-1)
       z_weight(jfixed) = z_zero
       gpxi = ( gx(ix+1) - gx(ix) ) / ( xx(ix+1) - xx(ix) )
       energy_fixed_point(jfixed) = xx(ix) - gx(ix) / gpxi

     endif

   endif
 enddo

end subroutine find_fixed_point


!=========================================================================
subroutine output_qp_energy(calcname,energy0,exchange_m_vxc,ncomponent,se,energy1,energy2,zz)
 implicit none

 character(len=*)             :: calcname
 integer                      :: ncomponent
 real(dp),intent(in)          :: energy0(:,:),exchange_m_vxc(:,:)
 type(selfenergy_grid),intent(in) :: se
 real(dp),intent(in)          :: energy1(:,:)
 real(dp),intent(in),optional :: energy2(:,:),zz(:,:)
!=====
 integer           :: pstate,ii
 character(len=36) :: sigc_label
!=====

 if( ncomponent > 2 ) call die('output_qp_energy: too many components. Not implemented yet')
 if( ncomponent < 1 ) call die('output_qp_energy: too few components. Something is not correct in the coding.')

 if( ncomponent == 1 ) then
   sigc_label = 'SigC'
 else
   if( nspin == 1 ) then
     sigc_label = 'SigC_1      SigC_2'
   else
     sigc_label = 'SigC_1                  SigC_2'
   endif
 endif


 write(stdout,'(/,1x,a,1x,a)') TRIM(calcname),'eigenvalues (eV)'

 if( PRESENT(zz) .AND. PRESENT(energy2) ) then

   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,10x,a,8x,a,5x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'Z','E_qp^lin','E_qp^graph'
   else
     write(stdout,'(3x,a,15x,a,22x,a,19x,a,24x,a,21x,a,17x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'Z','E_qp^lin','E_qp^graph'
     write(stdout,'(12x,14(a4,9x))') (' up ','down',ii=1,5+ncomponent)
   endif

   do pstate=nsemin,nsemax

     write(stdout,'(i4,1x,20(1x,f12.6))') pstate,energy0(pstate,:)*Ha_eV,  &
                                          exchange_m_vxc(pstate,:)*Ha_eV,  &
                                          REAL(se%sigma(0,pstate,:)*Ha_eV,dp),  &
                                          zz(pstate,:),                    &
                                          energy1(pstate,:)*Ha_eV,         &
                                          energy2(pstate,:)*Ha_eV
   enddo

 else

   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,9x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'E_qp'
   else
     write(stdout,'(3x,a,15x,a,22x,a,20x,a,22x,a)') '#','E0','SigX-Vxc',TRIM(sigc_label),'E_qp'
     write(stdout,'(12x,10(a4,9x))') (' up ','down',ii=1,3+ncomponent)
   endif

   do pstate=nsemin,nsemax

     write(stdout,'(i4,1x,20(1x,f12.6))') pstate,energy0(pstate,:)*Ha_eV,       &
                                          exchange_m_vxc(pstate,:)*Ha_eV,       &
                                          REAL(se%sigma(0,pstate,:)*Ha_eV,dp),  &
                                          energy1(pstate,:)*Ha_eV
   enddo



 endif

end subroutine output_qp_energy


!=========================================================================
subroutine output_qp_energy_yaml(calcname,energy0,exchange_m_vxc,se,energy1,energy2,zz)
 implicit none

 character(len=*)             :: calcname
 real(dp),intent(in)          :: energy0(:,:),exchange_m_vxc(:,:)
 type(selfenergy_grid),intent(in) :: se
 real(dp),intent(in)          :: energy1(:,:)
 real(dp),intent(in),optional :: energy2(:,:),zz(:,:)
!=====
 integer          :: pstate,ii,pspin
 character(len=6) :: char6
!=====

 if( .NOT. ( print_yaml_ .AND. is_iomaster ) ) return

 write(unit_yaml,'(/,a,a)') lower(calcname),' selfenergy:'
 write(unit_yaml,'(4x,a)') 'correlation:'
 write(unit_yaml,'(8x,a)') 'unit: eV'
 do pspin=1,nspin
   write(unit_yaml,'(8x,a,i2,a)') 'spin channel',pspin,':'
   do pstate=nsemin,nsemax
     write(char6,'(i6)') pstate
     write(unit_yaml,'(12x,a6,a,1x,es18.8)') ADJUSTL(char6),':',REAL(se%sigma(0,pstate,pspin),dp) * Ha_eV
   enddo
 enddo

 write(unit_yaml,'(4x,a)') 'exchange minus vxc:'
 write(unit_yaml,'(8x,a)') 'unit: eV'
 do pspin=1,nspin
   write(unit_yaml,'(8x,a,i2,a)') 'spin channel',pspin,':'
   do pstate=nsemin,nsemax
     write(char6,'(i6)') pstate
     write(unit_yaml,'(12x,a6,a,1x,es18.8)') ADJUSTL(char6),':',REAL(exchange_m_vxc(pstate,pspin),dp) * Ha_eV
   enddo
 enddo

 if( PRESENT(zz) ) then
   write(unit_yaml,'(4x,a)') 'renormalization factor:'
   do pspin=1,nspin
     write(unit_yaml,'(8x,a,i2,a)') 'spin channel',pspin,':'
     do pstate=nsemin,nsemax
       write(char6,'(i6)') pstate
       write(unit_yaml,'(12x,a6,a,1x,es18.8)') ADJUSTL(char6),':',zz(pstate,pspin)
     enddo
   enddo
 endif


end subroutine output_qp_energy_yaml


!=========================================================================
subroutine init_selfenergy_grid(selfenergy_technique,energy0,se)
 implicit none

 integer,intent(in)                  :: selfenergy_technique
 real(dp),intent(in)                 :: energy0(:,:)
 type(selfenergy_grid),intent(inout) :: se
!=====
 real(dp),parameter   :: alpha=1.0_dp
 real(dp),parameter   :: beta=1.0_dp
 integer              :: iomega,pstate
 real(dp),allocatable :: omega_gaussleg(:)
 real(dp)             :: efermi
 integer              :: iunittmp
 logical              :: manual_efermi
!=====

 se%nomega_calc = 0
 se%nomega      = 0

 inquire(file='manual_efermi',exist=manual_efermi)
 if(manual_efermi) then
   open(newunit=iunittmp,file='manual_efermi',action='read')
   read(iunittmp,*) efermi
   close(iunittmp)
   !
   ! efermi needs to be in the HOMO-LUMO gap
   if( efermi < MAXVAL(energy0(nhomo_G,:)) .OR. efermi > MINVAL(energy0(nhomo_G+1,:)) ) then
     write(stdout,*) 'efermi is out of the HOMO-LUMO gap:',efermi,&
                     MAXVAL(energy0(nhomo_G,:)),MINVAL(energy0(nhomo_G+1,:))
     call die('init_selfenergy_grid: efermi needs to be in the HOMO-LUMO gap')
   endif
 else
   efermi = 0.5_dp * ( MAXVAL(energy0(nhomo_G,:)) + MINVAL(energy0(nhomo_G+1,:)) )
 endif

 select case(selfenergy_technique)

 case(EVSC,static_selfenergy)

   allocate(se%omega(-se%nomega:se%nomega))
   se%omega(0) = 0.0_dp

 case(imaginary_axis_pade)
   !
   ! Set the final sampling points for Sigma on the real axis
   se%nomega = nomega_sigma/2
   allocate(se%omega(-se%nomega:se%nomega))
   do iomega=-se%nomega,se%nomega
     se%omega(iomega) = efermi + step_sigma * iomega
   enddo

   !
   ! Set the calculated sampling points for Sigma on the imaginary axis
   se%nomega_calc = nomega_sigma_calc
   allocate(se%omega_calc(se%nomega_calc))
   do iomega=1,se%nomega_calc
     se%omega_calc(iomega) = efermi + step_sigma_calc * (iomega-1) * im
   enddo

 case(imaginary_axis_homolumo)
   !
   ! Set the final sampling points for Sigma on the real axis
   se%nomega = nomega_sigma/2
   allocate(se%omega(-se%nomega:se%nomega))
   do iomega=-se%nomega,se%nomega
     se%omega(iomega) = efermi + step_sigma * iomega
   enddo

   !
   ! Set the initial sampling points for Sigma on the real axis inside the HOMO-LUMO gap
   se%nomega_calc = nomega_sigma_calc
   allocate(se%omega_calc(se%nomega_calc))
   do iomega=1,se%nomega_calc
     se%omega_calc(iomega) = MAXVAL(energy0(nhomo_G,:)) + 0.01_dp &
                             + (iomega-1) / (se%nomega_calc-1.0_dp) &
                                 * ( MINVAL(energy0(nhomo_G+1,:)) - MAXVAL(energy0(nhomo_G,:)) - 0.02_dp)
   enddo

 case(imaginary_axis_integral)
   !
   ! No final sampling points for Sigma on the real axis
   se%nomega = 0
   allocate(se%omega(-se%nomega:se%nomega))
   !
   ! Set the calculated sampling points for Sigma on the imaginary axis
   ! so to have a Gauss-Legendre type quadrature
   se%nomega_calc = nomega_sigma
   allocate(se%omega_calc(se%nomega_calc))
   allocate(se%weight_calc(se%nomega_calc))
   allocate(omega_gaussleg(se%nomega_calc))
   call coeffs_gausslegint(0.0_dp,1.0_dp,omega_gaussleg,se%weight_calc,se%nomega_calc)

   ! Variable change [0,1] -> [0,+\inf[
   do iomega=1,se%nomega_calc
     se%weight_calc(iomega) = se%weight_calc(iomega) / ( 2.0_dp**alpha - 1.0_dp ) &
                                 * alpha * (1.0_dp -  omega_gaussleg(iomega))**(-alpha-1.0_dp) * beta
     se%omega_calc(iomega)  = im / ( 2.0_dp**alpha - 1.0_dp ) &
                               * ( 1.0_dp / (1.0_dp - omega_gaussleg(iomega))**alpha - 1.0_dp ) * beta
   enddo
   deallocate(omega_gaussleg)

 case(one_shot)
   !
   ! Most standard case:
   se%nomega = nomega_sigma/2
   allocate(se%omega(-se%nomega:se%nomega))
   do iomega=-se%nomega,se%nomega
     se%omega(iomega) = step_sigma * iomega
   enddo

 end select


 !
 ! Set the central point of the grid
 !
 allocate(se%energy0(nsemin:nsemax,nspin))

 select case(selfenergy_technique)
 case(imaginary_axis_pade,imaginary_axis_homolumo)
   ! in this case the central point is already included in the complex frequency se%omega_calc
   se%energy0(nsemin:nsemax,:) = 0.0_dp
 case default
   se%energy0(nsemin:nsemax,:) = energy0(nsemin:nsemax,:)
 end select

 allocate(se%sigma(-se%nomega:se%nomega,nsemin:nsemax,nspin))

 select case(selfenergy_technique)
 case(imaginary_axis_pade,imaginary_axis_integral,imaginary_axis_homolumo)
   allocate(se%sigma_calc(se%nomega_calc,nsemin:nsemax,nspin))
 end select

end subroutine init_selfenergy_grid


!=========================================================================
subroutine destroy_selfenergy_grid(se)
 implicit none
   type(selfenergy_grid),intent(inout) :: se
  !=====

   se%nomega      = 0
   se%nomega_calc = 0
   if( ALLOCATED(se%omega) )  deallocate(se%omega)
   if( ALLOCATED(se%omega_calc) ) deallocate(se%omega_calc)
   deallocate(se%energy0)
   deallocate(se%sigma)
   if( ALLOCATED(se%sigma_calc) )  deallocate(se%sigma_calc)
   if( ALLOCATED(se%weight_calc) ) deallocate(se%weight_calc)

  end subroutine destroy_selfenergy_grid


  !=========================================================================
  subroutine setup_exchange_m_vxc(basis,occupation,energy,c_matrix,hamiltonian_fock,exchange_m_vxc)
   implicit none

   type(basis_set),intent(in) :: basis
   real(dp),intent(in)        :: occupation(:,:)
   real(dp),intent(in)        :: energy(:,:)
   real(dp),intent(in)        :: c_matrix(:,:,:)
   real(dp),intent(in)        :: hamiltonian_fock(:,:,:)
   real(dp),intent(out)       :: exchange_m_vxc(:,:,:)
  !=====
   integer              :: nstate
   integer              :: ispin,pstate
   real(dp)             :: exc
   real(dp),allocatable :: occupation_tmp(:,:)
   real(dp),allocatable :: p_matrix_tmp(:,:,:)
   real(dp),allocatable :: hxc_val(:,:,:),hexx_val(:,:,:),hxmxc(:,:,:)
  !=====

   call start_clock(timing_x_m_vxc)
   write(stdout,*) 'Calculate the (Sigma_x - Vxc) matrix'

   nstate = SIZE(occupation,DIM=1)

   !
   ! Testing the core/valence splitting
   !
   if(dft_core > 0) then
     if( beta_hybrid > 0.001 ) then
       call die('setup_exchange_m_vxc: RSH not implemented for DFT core-valence splitting')
     endif
     write(msg,'(a,i4,2x,i4)') 'DFT core-valence interaction switched on up to state = ',dft_core
     call issue_warning(msg)

     allocate(occupation_tmp(nstate,nspin))
     allocate(p_matrix_tmp(basis%nbf,basis%nbf,nspin))
     allocate(hxc_val(basis%nbf,basis%nbf,nspin))
     allocate(hexx_val(basis%nbf,basis%nbf,nspin))
     allocate(hxmxc(basis%nbf,basis%nbf,nspin))

     ! Override the occupation of the core electrons
     occupation_tmp(:,:) = occupation(:,:)
     occupation_tmp(1:dft_core,:) = 0.0_dp

     if( calc_type%is_dft ) then
       call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.FALSE.,BATCH_SIZE)
       call setup_density_matrix(c_matrix,occupation_tmp,p_matrix_tmp)
       call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation_tmp,c_matrix,hxc_val,exc)
       call destroy_dft_grid()
     endif

     call calculate_exchange(basis,p_matrix_tmp,hexx_val,occupation=occupation_tmp,c_matrix=c_matrix)

     hxc_val(:,:,:) = hxc_val(:,:,:) + alpha_hybrid * hexx_val(:,:,:)
     hxmxc(:,:,:) = hexx_val(:,:,:) - hxc_val(:,:,:)

     deallocate(occupation_tmp,p_matrix_tmp)

     !
     ! Calculate the matrix Sigma_x - Vxc
     ! for the forthcoming GW corrections
     !
     call matrix_ao_to_mo(c_matrix,hxmxc,exchange_m_vxc)

     deallocate(hxc_val,hexx_val,hxmxc)

   else

     !
     ! Calculate the matrix < p | Sigma_x - Vxc | q >
     ! this is equal to < p | F - H | q > and < p | H | q > = e_p \delta_{pq}

     call matrix_ao_to_mo(c_matrix,hamiltonian_fock,exchange_m_vxc)

     do ispin=1,nspin
       do pstate=1,nstate
         exchange_m_vxc(pstate,pstate,ispin) = exchange_m_vxc(pstate,pstate,ispin) - energy(pstate,ispin)
       enddo
     enddo

   endif

   call stop_clock(timing_x_m_vxc)

  end subroutine setup_exchange_m_vxc


  !=========================================================================
  subroutine apply_qs_approximation(s_matrix,c_matrix,selfenergy)
   implicit none

   real(dp),intent(in)    :: s_matrix(:,:),c_matrix(:,:,:)
   real(dp),intent(inout) :: selfenergy(:,:,:)
  !=====
   integer :: ispin
  !=====

   !
   ! Kotani's Hermitianization trick
   !
   do ispin=1,nspin
     selfenergy(:,:,ispin) = 0.5_dp * ( selfenergy(:,:,ispin) + TRANSPOSE(selfenergy(:,:,ispin)) )

     ! Transform the matrix elements back to the AO basis
     ! do not forget the overlap matrix S
     ! C^T S C = I
     ! the inverse of C is C^T S
     ! the inverse of C^T is S C
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , &
                               MATMUL( selfenergy(nsemin:nsemax,nsemin:nsemax,ispin), &
                                 MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )

   enddo

  end subroutine apply_qs_approximation


  !=========================================================================
  subroutine self_energy_fit(se)
   implicit none

   type(selfenergy_grid),intent(inout) :: se
 !=====
 integer :: pstate,pspin
 integer :: iomega
 integer :: ilbfgs,ipp,ii,info
 real(dp),parameter :: dq=1.0e-4_dp
 integer,parameter :: nlbfgs=100
 integer,parameter :: npp=6
 integer,parameter :: nparam=6
 real(dp)          :: coeff(nparam*npp)
 real(dp)          :: coeffdq(nparam*npp)
 real(dp)          :: dchi2(nparam*npp)
 real(dp)          :: chi2
 type(lbfgs_state) :: lbfgs_plan
 !=====

 do pspin=1,nspin
   do pstate=nsemin,nsemax

     !
     ! Optimization: first guess
     do ipp=1,npp
       coeff(1+(ipp-1)*nparam) = 0.5_dp  * ipp
       coeff(2+(ipp-1)*nparam) = 0.5_dp  * ipp
       coeff(3+(ipp-1)*nparam) = 0.01_dp * ipp
       coeff(4+(ipp-1)*nparam) = 0.01_dp * ipp
       coeff(5+(ipp-1)*nparam) = 0.3_dp  / ipp
       coeff(6+(ipp-1)*nparam) = 0.3_dp  / ipp
     enddo
     chi2 = eval_chi2(coeff)

     call lbfgs_init(lbfgs_plan,nparam*npp,5,diag_guess=1.0_dp)
     do ilbfgs=1,nlbfgs
       write(stdout,*) 'chi2=',chi2
       do ipp=1,npp
         do ii=1,nparam
           coeffdq(:) = coeff(:)
           coeffdq(ii+(ipp-1)*nparam) = coeffdq(ii+(ipp-1)*nparam) + dq
           dchi2(ii+(ipp-1)*nparam) = ( eval_chi2(coeffdq) - eval_chi2(coeff) ) / dq
         enddo
       enddo
       info = lbfgs_execute(lbfgs_plan,coeff,chi2,dchi2)
       chi2 = eval_chi2(coeff)
       if( chi2 < 5.0e-8_dp ) exit
     enddo
     call lbfgs_destroy(lbfgs_plan)
     write(stdout,'(/,1x,a)') '=========================='
     write(stdout,'(1x,a)') '   #           Coeff              Re Pole            Im Pole'
     do ipp=1,npp
       write(stdout,'(1x,i4,4(2x,f18.8))') 2*ipp-1,                   &
                                           coeff(5+(ipp-1)*nparam)**2, &
                                           coeff(1+(ipp-1)*nparam)**2, &
                                          -coeff(3+(ipp-1)*nparam)**2
       write(stdout,'(1x,i4,4(2x,f18.8))') 2*ipp,                     &
                                           coeff(6+(ipp-1)*nparam)**2, &
                                          -coeff(2+(ipp-1)*nparam)**2, &
                                           coeff(4+(ipp-1)*nparam)**2
     enddo
     write(stdout,'(1x,a)') '=========================='


     do iomega=1,se%nomega_calc
       write(500+pstate,'(6(1x,f18.8))') (se%omega_calc(iomega) + se%energy0(pstate,pspin))*Ha_eV, &
                                         se%sigma_calc(iomega,pstate,pspin)*Ha_eV, &
                                         eval_func(coeff, se%omega_calc(iomega) + se%energy0(pstate,pspin) )*Ha_eV
     enddo


     ! Extrapolate to the final desired omega's
     do iomega=-se%nomega,se%nomega
       se%sigma(iomega,pstate,pspin) = eval_func(coeff, se%omega(iomega) + se%energy0(pstate,pspin) )
     enddo

   enddo
 enddo




contains


function eval_func(coeff_in,zz)
 implicit none
 real(dp),intent(in)    :: coeff_in(nparam*npp)
 complex(dp),intent(in) :: zz
 complex(dp)            :: eval_func
!=====
 integer :: ipp
!=====

 eval_func = 0.0_dp
 do ipp=1,npp
   eval_func = eval_func + coeff_in(5+(ipp-1)*nparam)**2  &
                             / ( zz - ( coeff_in(1+(ipp-1)*nparam)**2 - im * coeff_in(3+(ipp-1)*nparam)**2 ) ) &
                         + coeff_in(6+(ipp-1)*nparam)**2  &
                             / ( zz + ( coeff_in(2+(ipp-1)*nparam)**2 - im * coeff_in(4+(ipp-1)*nparam)**2 ) )
 enddo

end function eval_func


function eval_chi2(coeff_in)
 implicit none
 real(dp),intent(in)    :: coeff_in(nparam*npp)
 real(dp)               :: eval_chi2
!=====
 integer  :: iomega_calc
 real(dp) :: weight
 real(dp) :: norm
!=====

 eval_chi2 = 0.0_dp
 norm      = 0.0_dp
 do iomega_calc=1,se%nomega_calc
   weight = 1.0_dp / ABS(1.0_dp+se%omega_calc(iomega_calc))**2
   eval_chi2 = eval_chi2         &
                + ABS( se%sigma_calc(iomega_calc,pstate,pspin) &
                      - eval_func(coeff_in, se%omega_calc(iomega_calc) + se%energy0(pstate,pspin) ) )**2 &
                      * weight
   norm = norm + weight

 enddo
 eval_chi2 = eval_chi2 / norm


end function eval_chi2


end subroutine self_energy_fit


!=========================================================================
subroutine self_energy_pade(se)
 implicit none

 type(selfenergy_grid),intent(inout) :: se
!=====
 integer  :: pstate,pspin
 integer  :: iomega,iomega_calc
 real(dp) :: sign_eta
 complex(dp) :: omega_sym(2*se%nomega_calc-1)
 complex(dp) :: sigma_sym(2*se%nomega_calc-1)
!=====

 do pspin=1,nspin
   do pstate=nsemin,nsemax

     ! First create the symmetric sigma
     ! using sigma(-iw) = sigma(iw)*
     omega_sym(se%nomega_calc) = se%omega_calc(1)
     sigma_sym(se%nomega_calc) = se%sigma_calc(1,pstate,pspin)
     do iomega_calc=2,se%nomega_calc
       omega_sym(se%nomega_calc+iomega_calc-1) = se%omega_calc(iomega_calc)
       sigma_sym(se%nomega_calc+iomega_calc-1) = se%sigma_calc(iomega_calc,pstate,pspin)
       omega_sym(se%nomega_calc-iomega_calc+1) = CONJG(se%omega_calc(iomega_calc))
       sigma_sym(se%nomega_calc-iomega_calc+1) = CONJG(se%sigma_calc(iomega_calc,pstate,pspin))
     enddo


     do iomega=-se%nomega,se%nomega
       sign_eta = -SIGN( 1.0_dp , se%omega(iomega)%re - se%omega_calc(1)%re )
       se%sigma(iomega,pstate,pspin) = pade(omega_sym,sigma_sym, se%omega(iomega) + ieta * sign_eta )
     enddo
   enddo
 enddo


end subroutine self_energy_pade


!=========================================================================
! Fit a low-order polynomial on the available real frequencies for sigma
! and make an extrapolation/interpolation to the requested frequencies
subroutine self_energy_polynomial(se)
 implicit none

 type(selfenergy_grid),intent(inout) :: se
!=====
 integer :: pspin,pstate,iomega
 real(dp) :: a0,a1,a2     ! Polynomial coefficients a0 + a1*x + a2*x**2 + ...
!=====

 do pspin=1,nspin
   do pstate=nsemin,nsemax

     a0 = se%sigma_calc(se%nomega_calc/2+1,pstate,pspin)%re
     a1 = ( se%sigma_calc(se%nomega_calc,pstate,pspin)%re - se%sigma_calc(1,pstate,pspin)%re ) &
            / ( se%omega_calc(se%nomega_calc)%re - se%omega_calc(1)%re )
     a2 = (            se%sigma_calc(se%nomega_calc,pstate,pspin)%re   &
            +          se%sigma_calc(1,pstate,pspin)%re                 &
            - 2.0_dp * se%sigma_calc(se%nomega_calc/2+1,pstate,pspin)%re ) &
            / ( se%omega_calc(se%nomega_calc)%re - se%omega_calc(1)%re )**2 * 2.0_dp

     do iomega=-se%nomega,se%nomega
       se%sigma(iomega,pstate,pspin) = a0 + a1 * (se%omega(iomega)-se%omega_calc(se%nomega_calc/2+1)) &
                                          + a2 * (se%omega(iomega)-se%omega_calc(se%nomega_calc/2+1))**2
     enddo

   enddo
 enddo


end subroutine self_energy_polynomial


!=========================================================================
! Prediction the CBS limit for GW
!
! No this is not black magic...
!
! Using the simplest model from Bruneval et al. JCTC 2020
! convergence error
! Delta E_i = A_basis + B_basis * ln( t_i )
! where t_i = < \phi_i | -\nabla^2 / 2 | \phi_i >
subroutine selfenergy_convergence_prediction(basis,c_matrix,eqp)
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(in)        :: eqp(:,:)
 !=====
 integer  :: pspin,pstate,iatom
 logical  :: basis_recognized
 real(dp) :: abasis,bbasis
 real(dp) :: hkin(basis%nbf,basis%nbf)
 real(dp) :: t_i(nsemin:nsemax,nspin)
 real(dp) :: deltae
 !=====

 !
 ! Be careful this routine is in eV !
 !
 write(stdout,'(/,1x,a)') 'Estimate the Complete Basis Set limit for free'
 write(stdout,*)          '  see Bruneval, Maliyov, Lapointe, Marinica, JCTC 16, 4399 (2020)'
 write(stdout,*)          '      https://doi.org/10.1021/acs.jctc.0c00433'

 !
 ! Retrieve the linear regression parameters trained on a benchmark of organic molecules (GW@BHLYP level)
 !
 !
 basis_recognized = .TRUE.
 do iatom=2,SIZE(basis_name(:))
   if( TRIM(basis_name(iatom)) /= TRIM(basis_name(1)) ) basis_recognized = .FALSE.
 enddo

 select case(TRIM(basis_name(1)))
 case('cc-pVDZ')
   abasis =  0.9383
   bbasis = -0.4500
 case('cc-pVTZ')
   abasis =  0.4863
   bbasis = -0.2156
 case('cc-pVQZ')
   abasis =  0.2703
   bbasis = -0.1066
 case('cc-pV5Z')
   abasis =  0.1271
   bbasis = -0.0483
 case('cc-pV6Z')
   abasis =  0.1271 * 0.5    ! evaluated
   bbasis = -0.0483 * 0.5    ! evaluated
 case('aug-cc-pVDZ')
   abasis =  0.5351
   bbasis = -0.3061
 case('aug-cc-pVTZ')
   abasis =   0.5063
   bbasis =  -0.2086
 case('aug-cc-pVQZ')
   abasis =  0.5063  * 0.5    ! evaluated
   bbasis = -0.2086  * 0.5    ! evaluated
 case('aug-cc-pV5Z')
   abasis =  0.5063  * 0.25    ! evaluated
   bbasis = -0.2086  * 0.25    ! evaluated
 case('aug-cc-pV6Z')
   abasis =  0.5063  * 0.125    ! evaluated
   bbasis = -0.2086  * 0.125    ! evaluated
 case default
   basis_recognized = .FALSE.
 end select

 if( .NOT. basis_recognized ) then
    write(stdout,*) 'Automatic extrapolation to CBS not possible because the basis set is not recognized'
    write(stdout,*) 'only fitted for a Dunning basis cc-pVXZ or aug-cc-pVXZ'
    write(stdout,*) 'only fitted for the same basis on all atoms'
    return
 endif

 call setup_kinetic(basis,hkin)

 do pspin=1,nspin
   do pstate=nsemin,nsemax
     t_i(pstate,pspin) = DOT_PRODUCT( c_matrix(:,pstate,pspin) ,  MATMUL( hkin(:,:) ,  c_matrix(:,pstate,pspin) ) ) * Ha_eV
   enddo
 enddo


 write(stdout,'(/,1x,a,a)')           'For basis: ',basis_name(1)
 write(stdout,'(1x,a,f7.4,a,f7.4,a)') 'Magical formula reads Delta E_i = ',abasis,' + ', &
                                      bbasis,' x LOG( <i|-\nabla^2/2|i>  (eV) ) (eV)'
 write(stdout,'(5x,a)')               'Formula was trained for organic molecules and GW@BHLYP'
 write(stdout,'(5x,a)')               'Accuracy is correct for cc-pVDZ and excellent above'
 write(stdout,'(5x,a)')               'Accuracy is excellent for aug-cc-pVDZ and above'

 write(stdout,'(/,1x,a)') 'Extrapolation to CBS (eV)'
 write(stdout,'(25x,a,a,a)') '<i|-\nabla^2/2|i>    Delta E_i     E_i(',TRIM(basis_name(1)),')      E_i(CBS)'
 do pspin=1,nspin
   do pstate=nsemin,nsemax
      deltae = abasis + bbasis * LOG( t_i(pstate,pspin) )
      write(stdout,'(1x,a,i4,a,i2,a,*(4x,f12.6))') 'state ',pstate,' spin ',pspin,' : ', &
                                            t_i(pstate,pspin),    &
                                            deltae, &
                                            eqp(pstate,pspin)*Ha_eV, &
                                            eqp(pstate,pspin)*Ha_eV + deltae
   enddo
 enddo


end subroutine selfenergy_convergence_prediction


!=========================================================================
end module m_selfenergy_tools
!=========================================================================
