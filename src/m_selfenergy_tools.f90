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
   integer                 :: nomegai
   complex(dp),allocatable :: omega(:)
   complex(dp),allocatable :: omegai(:)
   real(dp),allocatable    :: energy0(:,:)
   complex(dp),allocatable :: sigma(:,:,:)
   complex(dp),allocatable :: sigmai(:,:,:)
 end type


contains


!=========================================================================
subroutine selfenergy_set_state_range(nstate,occupation)
 implicit none
 integer             :: nstate
 real(dp),intent(in) :: occupation(:,:)
!=====
 integer :: istate
!=====

 if( nstate > SIZE( occupation(:,:) , DIM=1 ) ) then
   call die('selfenergy_set_state_range: nstate is too large')
 endif

 ncore_G      = ncoreg
 nvirtual_G   = MIN(nvirtualg,nstate+1)

 if(is_frozencore) then
   if( ncore_G == 0) ncore_G = atoms_core_states()
 endif

 if( ncore_G > 0 ) then
   write(msg,'(a,i4,2x,i4)') 'frozen core approximation in G switched on up to state = ',ncore_G
   call issue_warning(msg)
 endif

 if( nvirtual_G <= nstate ) then
   write(msg,'(a,i4,2x,i4)') 'frozen virtual approximation in G switched on starting with state = ',nvirtual_G
   call issue_warning(msg)
 endif

 ! Find the HOMO index
 nhomo_G = 1
 do istate=1,nstate
   if( .NOT. ANY( occupation(istate,:) < completely_empty ) ) then
     nhomo_G = MAX(nhomo_G,istate)
   endif
 enddo

 nsemin = MAX(ncore_G+1,selfenergy_state_min,1,nhomo_G-selfenergy_state_range)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate,nhomo_G+selfenergy_state_range)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax

end subroutine selfenergy_set_state_range


!=========================================================================
subroutine write_selfenergy_omega(filename_root,nstate,exchange_m_vxc,energy0,se)
 implicit none

 character(len=*)    :: filename_root
 integer,intent(in)  :: nstate
 real(dp),intent(in) :: exchange_m_vxc(nstate,nspin)
 real(dp),intent(in) :: energy0(nstate,nspin)
 type(selfenergy_grid),intent(in) :: se
!=====
 character(len=3)   :: ctmp
 character(len=256) :: filename
 integer :: selfenergyfile
 integer :: pstate
 integer :: iomega
 real(dp) :: spectral_function_w(nspin)
!=====

 ! Just the master writes
 if( .NOT. is_iomaster ) return

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

   write(selfenergyfile,'(a)') '#       omega (eV)          Re SigmaC (eV)     Im SigmaC (eV)    omega - e_gKS - Vxc + SigmaX (eV)     A (eV^-1)'

   do iomega=-se%nomega,se%nomega
     spectral_function_w(:) = 1.0_dp / pi * ABS(   &
                                       AIMAG( 1.0_dp   &
                                         / ( se%energy0(pstate,:)+se%omega(iomega) - energy0(pstate,:)  &
                                               - exchange_m_vxc(pstate,:) - se%sigma(iomega,pstate,:) ) ) )

     write(selfenergyfile,'(20(f16.8,2x))') ( se%omega(iomega) + se%energy0(pstate,:) )*Ha_eV,             &
                                            REAL(se%sigma(iomega,pstate,:),dp) * Ha_eV,                    &
                                            AIMAG(se%sigma(iomega,pstate,:)) * Ha_eV,                      &
                                            (REAL(se%omega(iomega),dp)+se%energy0(pstate,:) - energy0(pstate,:) - exchange_m_vxc(pstate,:) )*Ha_eV, &
                                            spectral_function_w(:) / Ha_eV
   enddo
   if( se%nomegai > 0 ) then
     do iomega=-se%nomegai,se%nomegai
       write(selfenergyfile,'(20(f16.8,2x))') ( se%omegai(iomega) + se%energy0(pstate,:) )*Ha_eV,     &
                                              REAL(se%sigmai(iomega,pstate,:),dp) * Ha_eV,            &
                                              AIMAG(se%sigmai(iomega,pstate,:)) * Ha_eV,              &
                                              0.0_dp,0.0_dp
     enddo
   endif
   write(selfenergyfile,*)
   close(selfenergyfile)

 enddo


end subroutine write_selfenergy_omega


!=========================================================================
subroutine find_qp_energy_linearization(se,nstate,exchange_m_vxc,energy0,energy_qp_z,zz)
 implicit none

 type(selfenergy_grid),intent(in) :: se
 integer,intent(in)               :: nstate
 real(dp),intent(in)              :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out)             :: energy_qp_z(nstate,nspin)
 real(dp),intent(out),optional    :: zz(nsemin:nsemax,nspin)
!=====
 integer  :: pstate,pspin
 real(dp) :: zz_p(nspin)
!=====

 ! First, a dummy initialization
 energy_qp_z(:,:) = energy0(:,:)

 ! Then overwrite the interesting energy with the calculated GW one
 do pstate=nsemin,nsemax

   if( se%nomega > 0 .AND. PRESENT(zz) ) then
     zz_p(:) = ( REAL(se%sigma(1,pstate,:),dp) - REAL(se%sigma(-1,pstate,:),dp) ) / ( se%omega(1) - se%omega(-1) )
     zz_p(:) = 1.0_dp / ( 1.0_dp - zz_p(:) )
     ! Contrain Z to be in [0:1] to avoid crazy values
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


end subroutine find_qp_energy_linearization


!=========================================================================
subroutine find_qp_energy_graphical(se,nstate,exchange_m_vxc,energy0,energy_qp_g)
 implicit none

 type(selfenergy_grid),intent(in) :: se
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out) :: energy_qp_g(nstate,nspin)
!=====
 integer  :: pstate,pspin
 integer  :: info
 real(dp) :: sigma_xc_m_vxc(-se%nomega:se%nomega)
!=====

 ! First, a dummy initialization
 energy_qp_g(:,:) = 0.0_dp

 ! Then overwrite the interesting energy with the calculated GW one
 do pstate=nsemin,nsemax

   if( MODULO(pstate-nsemin,nproc_world) /= rank_world ) cycle

   do pspin=1,nspin
     sigma_xc_m_vxc(:) = REAL(se%sigma(:,pstate,pspin),dp) + exchange_m_vxc(pstate,pspin)
     !
     ! QP equation:
     ! E_GW = E0 + \omega =  E_gKS + \Sigma_c(E0+\omega) + \Sigma_x - v_xc
     !
     energy_qp_g(pstate,pspin) = find_fixed_point(se%nomega,REAL(se%omega,dp)+se%energy0(pstate,pspin),sigma_xc_m_vxc(:)+energy0(pstate,pspin),info)
     if( info /= 0 ) then
       write(stdout,'(1x,a,i4,2x,i4)') 'Unreliable graphical solution of the QP equation for state,spin: ',pstate,pspin
       call issue_warning('No fixed point found for the QP equation. Try to increase nomega_sigma or step_sigma.')
     endif
   enddo

 enddo

 call xsum_world(energy_qp_g)

 energy_qp_g(:nsemin-1,:) = energy0(:nsemin-1,:)
 energy_qp_g(nsemax+1:,:) = energy0(nsemax+1:,:)



end subroutine find_qp_energy_graphical


!=========================================================================
function find_fixed_point(nx,xx,fx,info) result(fixed_point)
 implicit none
 integer,intent(in)  :: nx
 real(dp),intent(in) :: xx(-nx:nx)
 real(dp),intent(in) :: fx(-nx:nx)
 integer,intent(out) :: info
 real(dp)            :: fixed_point
!=====
 integer             :: ix,imin1,imin2
 real(dp)            :: rmin
 real(dp)            :: gx(-nx:nx)
 real(dp)            :: gpx
!=====

 !
 ! g(x) contains f(x) - x
 gx(:) = fx(:) - xx(:)

 ! Find the sign change in g(x) which is the closest to ix=0
 ! Search positive x
 imin1 = 1000000
 do ix=0,nx-1
   if( gx(ix) * gx(ix+1) < 0.0_dp ) then
     imin1 = ix
     exit
   endif
 enddo
 ! Search negative x
 imin2 = 1000000
 do ix=0,-nx+1,-1
   if( gx(ix) * gx(ix-1) < 0.0_dp ) then
     imin2 = ix
     exit
   endif
 enddo

 if( imin1 == 1000000 .AND. imin2 == 1000000 ) then

   if( gx(0) > 0.0_dp ) then 
     info =  1
     fixed_point = xx(nx)
   else
     info = -1
     fixed_point = xx(-nx)
   endif

 else
   info = 0
   ! 
   ! If sign change found, interpolate the abscissa between 2 grid points
   if( ABS(imin1) <= ABS(imin2) )  then
     gpx = ( gx(imin1+1) - gx(imin1) ) / ( xx(imin1+1) - xx(imin1) )
     fixed_point = xx(imin1) - gx(imin1) / gpx 
   else
     gpx = ( gx(imin2) - gx(imin2-1) ) / ( xx(imin2) - xx(imin2-1) )
     fixed_point = xx(imin2-1) - gx(imin2-1) / gpx 
   endif
 endif


end function find_fixed_point


!=========================================================================
subroutine output_qp_energy(calcname,nstate,energy0,exchange_m_vxc,ncomponent,se,energy1,energy2,zz)
 implicit none
 
 character(len=*)             :: calcname
 integer                      :: nstate,ncomponent
 real(dp),intent(in)          :: energy0(nstate,nspin),exchange_m_vxc(nstate,nspin)
 type(selfenergy_grid),intent(in) :: se
 real(dp),intent(in)          :: energy1(nstate,nspin)
 real(dp),intent(in),optional :: energy2(nstate,nspin),zz(nsemin:nsemax,nspin)
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
subroutine init_selfenergy_grid(selfenergy_technique,nstate,energy0,se)
 use m_atoms
 implicit none

 integer,intent(in)                  :: selfenergy_technique,nstate
 real(dp),intent(in)                 :: energy0(nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer            :: iomega,pstate
!=====

 se%nomegai = 0
 se%nomega  = 0

 select case(selfenergy_technique)
 
 case(EVSC,static)

   allocate(se%omega(-se%nomega:se%nomega))
   se%omega(0) = 0.0_dp

 case(imaginary_axis)
   !
   ! Set the final sampling points for Sigma
   se%nomega = nomega_sigma/2
   allocate(se%omega(-se%nomega:se%nomega))
   do iomega=-se%nomega,se%nomega
     se%omega(iomega) = step_sigma * iomega
   enddo

   !
   ! Set the calculated sampling points for Sigma
   se%nomegai = nomega_sigma / 2
   allocate(se%omegai(-se%nomegai:se%nomegai))
   do iomega=-se%nomegai,se%nomegai
     se%omegai(iomega) = step_sigma * iomega * im * 3.0_dp
   enddo


 case(one_shot)
   select case(calc_type%selfenergy_approx)
   case(GSIGMA)
     se%nomega = 1

   case default
     se%nomega = nomega_sigma/2
     allocate(se%omega(-se%nomega:se%nomega))
     do iomega=-se%nomega,se%nomega
       se%omega(iomega) = step_sigma * iomega
     enddo

   end select

 end select


 !
 ! Set the central point of the grid
 allocate(se%energy0(nsemin:nsemax,nspin))

 select case(selfenergy_technique)
 case(imaginary_axis)
   ! Find the center of the HOMO-LUMO gap
   forall(pstate=nsemin:nsemax)
     se%energy0(pstate,:) = 0.5_dp * ( energy0(nhomo_G,:) + energy0(nhomo_G+1,:) )
   end forall
   forall(pstate=nsemin:nhomo_G)
     se%energy0(pstate,:) = energy0(nhomo_G,:) + 0.05_dp
   end forall
   forall(pstate=nhomo_G+1:nsemax)
     se%energy0(pstate,:) = energy0(nhomo_G+1,:) - 0.05_dp
   end forall
 case default
   se%energy0(nsemin:nsemax,:) = energy0(nsemin:nsemax,:)
 end select

 !
 ! Set the central point of the grid
 allocate(se%sigma(-se%nomega:se%nomega,nsemin:nsemax,nspin))
 select case(selfenergy_technique)
 case(imaginary_axis)
   allocate(se%sigmai(-se%nomegai:se%nomegai,nsemin:nsemax,nspin))
 end select


end subroutine init_selfenergy_grid


!=========================================================================
subroutine destroy_selfenergy_grid(se)
 implicit none
 type(selfenergy_grid),intent(inout) :: se
!=====

 se%nomega  = 0
 se%nomegai = 0
 if( ALLOCATED(se%omega) )  deallocate(se%omega)
 if( ALLOCATED(se%omegai) ) deallocate(se%omegai)
 deallocate(se%energy0)
 deallocate(se%sigma)
 if( ALLOCATED(se%sigmai) ) deallocate(se%sigmai)

end subroutine destroy_selfenergy_grid


!=========================================================================
subroutine setup_exchange_m_vxc_diag(basis,nstate,occupation,energy,c_matrix,hamiltonian_fock,exchange_m_vxc_diag)
 use m_inputparam
 use m_basis_set
 use m_dft_grid
 use m_hamiltonian
 use m_hamiltonian_sca
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: hamiltonian_fock(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: exchange_m_vxc_diag(nstate,nspin)
!=====
 integer              :: ispin,istate
 real(dp)             :: exc,eexx
 real(dp),allocatable :: occupation_tmp(:,:)
 real(dp),allocatable :: p_matrix_tmp(:,:,:)
 real(dp),allocatable :: hxc_val(:,:,:),hexx_val(:,:,:),hxmxc(:,:,:)
!=====

 !
 ! Testing the core/valence splitting
 !
 if(dft_core > 0) then
   if( alpha_hybrid_lr > 0.001 ) then
     call die('RSH not implemented yet')
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
     call init_dft_grid(grid_level)
     call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation_tmp,p_matrix_tmp)
     call dft_exc_vxc(basis,nstate,occupation_tmp,c_matrix,p_matrix_tmp,hxc_val,exc)
     call destroy_dft_grid()
   endif

   if( .NOT. has_auxil_basis ) then
     call setup_exchange(print_matrix_,basis%nbf,p_matrix_tmp,hexx_val,eexx)
   else
     if( parallel_ham ) then
       call die('setup_exchange_m_vxc_diag: case not implemented')
     else
       call setup_exchange_ri(basis%nbf,nstate,occupation_tmp,c_matrix,p_matrix_tmp,hexx_val,eexx)
     endif
   endif

   hxc_val(:,:,:) = hxc_val(:,:,:) + alpha_hybrid * hexx_val(:,:,:)
   hxmxc(:,:,:) = hexx_val(:,:,:) - hxc_val(:,:,:) 

   deallocate(occupation_tmp,p_matrix_tmp)

   !
   ! Calculate the diagonal of the matrix Sigma_x - Vxc
   ! for the forthcoming GW corrections
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,nstate

        exchange_m_vxc_diag(istate,ispin) =  DOT_PRODUCT(  c_matrix(:,istate,ispin) , &
                                                MATMUL( hxmxc(:,:,ispin) , c_matrix(:,istate,ispin) ) )
     enddo
   enddo
   deallocate(hxc_val,hexx_val,hxmxc)

 else

   !
   ! Calculate the diagonal of the matrix Sigma_x - Vxc
   ! for the forthcoming GW corrections
   exchange_m_vxc_diag(:,:) = 0.0_dp
   do ispin=1,nspin
     do istate=1,nstate

        exchange_m_vxc_diag(istate,ispin) =  DOT_PRODUCT(  c_matrix(:,istate,ispin) , &
                                                MATMUL( hamiltonian_fock(:,:,ispin) , c_matrix(:,istate,ispin) ) ) &
                                              - energy(istate,ispin)
     enddo
   enddo

 endif


end subroutine setup_exchange_m_vxc_diag


!=========================================================================
subroutine apply_qs_approximation(nbf,nstate,s_matrix,c_matrix,selfenergy)
 implicit none

 integer,intent(in)     :: nbf,nstate
 real(dp),intent(in)    :: s_matrix(nbf,nbf),c_matrix(nbf,nstate,nspin)
 real(dp),intent(inout) :: selfenergy(nbf,nbf,nspin)
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
subroutine self_energy_fit(nstate,energy,se)
 use m_lbfgs
 implicit none

 integer,intent(in)                  :: nstate
 real(dp),intent(in)                 :: energy(nstate,nspin)
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer :: pstate,pspin
 integer :: iomega
 integer :: ilbfgs,ipp,ii,info
 real(dp),parameter :: dq=1.0e-5_dp
 integer,parameter :: nlbfgs=100
 integer,parameter :: npp=4
 real(dp)          :: coeff(6*npp)
 real(dp)          :: coeffdq(6*npp)
 real(dp)          :: dchi2(6*npp)
 real(dp)          :: chi2
 real(dp)          :: shift(6*npp)
 type(lbfgs_state) :: lbfgs_plan
!=====

 do pspin=1,nspin
   do pstate=nsemin,nsemax

     !
     ! Optimization: first guess
     do ipp=1,npp
       coeff(1+(ipp-1)*6) = 0.5_dp  * ipp
       coeff(2+(ipp-1)*6) = 0.5_dp  * ipp
       coeff(3+(ipp-1)*6) = 0.01_dp * ipp
       coeff(4+(ipp-1)*6) = 0.01_dp * ipp
       coeff(5+(ipp-1)*6) = 0.3_dp  / ipp
       coeff(6+(ipp-1)*6) = 0.3_dp  / ipp
     enddo
     chi2 = eval_chi2(coeff)

     call lbfgs_init(lbfgs_plan,6*npp,5,diag_guess=1.0_dp)
     do ilbfgs=1,nlbfgs
       write(stdout,*) 'chi2=',chi2
       do ipp=1,npp
         do ii=1,6
           coeffdq(:) = coeff(:)
           coeffdq(ii+(ipp-1)*6) = coeffdq(ii+(ipp-1)*6) + dq 
           dchi2(ii+(ipp-1)*6) = ( eval_chi2(coeffdq) - eval_chi2(coeff) ) / dq
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
       write(stdout,'(1x,i4,3(2x,f18.8))') 2*ipp-1,                   &
                                           coeff(5+(ipp-1)*6)**2, &
                                           coeff(1+(ipp-1)*6)**2, &
                                          -coeff(3+(ipp-1)*6)**2
       write(stdout,'(1x,i4,3(2x,f18.8))') 2*ipp,                     &
                                           coeff(6+(ipp-1)*6)**2, &
                                          -coeff(2+(ipp-1)*6)**2, &
                                           coeff(4+(ipp-1)*6)**2
     enddo
     write(stdout,'(1x,a)') '=========================='


     do iomega=-se%nomegai,se%nomegai
       write(500+pstate,'(6(1x,f18.8))') (se%omegai(iomega) + se%energy0(pstate,pspin))*Ha_eV, &
                                         se%sigmai(iomega,pstate,pspin)*Ha_eV, &
                                         eval_func(coeff, se%omegai(iomega) + se%energy0(pstate,pspin) )*Ha_eV
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
 real(dp),intent(in)    :: coeff_in(6*npp)
 complex(dp),intent(in) :: zz
 complex(dp)            :: eval_func
!=====
 integer :: ipp
!=====

 eval_func = 0.0_dp
 do ipp=1,npp
   eval_func = eval_func + coeff_in(5+(ipp-1)*6)**2 &
                             / ( zz - ( coeff_in(1+(ipp-1)*6)**2 - im * coeff_in(3+(ipp-1)*6)**2 ) ) &
                         + coeff_in(6+(ipp-1)*6)**2 &
                             / ( zz + ( coeff_in(2+(ipp-1)*6)**2 - im * coeff_in(4+(ipp-1)*6)**2 ) )
 enddo

end function eval_func


function eval_chi2(coeff_in)
 implicit none
 real(dp),intent(in)    :: coeff_in(6*npp)
 real(dp)               :: eval_chi2
!=====
 integer  :: iomega
 real(dp) :: weight
 real(dp) :: norm
!=====

 eval_chi2 = 0.0_dp
 norm      = 0.0_dp
 do iomega=-se%nomegai,se%nomegai
   weight = 1.0_dp / ABS(1.0_dp+se%omegai(iomega))**2
   eval_chi2 = eval_chi2         &
                + ABS( se%sigmai(iomega,pstate,pspin) &
                      - eval_func(coeff_in, se%omegai(iomega) + se%energy0(pstate,pspin) ) )**2 &
                      * weight
   norm = norm + weight
                              
 enddo 
 eval_chi2 = eval_chi2 / norm


end function eval_chi2


end subroutine self_energy_fit

!=========================================================================
end module m_selfenergy_tools
!=========================================================================
