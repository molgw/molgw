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
 use m_inputparam,only: nspin,ncoreg,nvirtualg,is_frozencore, &
                        selfenergy_state_min,selfenergy_state_max,selfenergy_state_range

 !
 ! frozen core approximation parameters
 integer,protected :: ncore_G
 integer,protected :: nvirtual_G

 integer,protected :: nsemin
 integer,protected :: nsemax

 !
 ! Highest occupied state
 integer,protected :: nhomo_G



contains


!=========================================================================
subroutine write_selfenergy_omega(filename_root,nstate,energy0,exchange_m_vxc,nomegai,omegai,nsemin,nsemax,selfenergy_omega)
 implicit none

 character(len=*)    :: filename_root
 integer,intent(in)  :: nstate,nsemin,nsemax,nomegai
 real(dp),intent(in) :: energy0(nstate,nspin),exchange_m_vxc(nstate,nspin)
 real(dp),intent(in) :: omegai(nomegai)
 real(dp),intent(in) :: selfenergy_omega(nomegai,nsemin:nsemax,nspin)
!=====
 character(len=3)   :: ctmp
 character(len=256) :: filename
 integer :: selfenergyfile
 integer :: astate
 integer :: iomegai
!=====

 ! Just the master writes
 if( .NOT. is_iomaster ) return

 write(stdout,'(/,1x,a)') 'Write Sigma(omega) on file'

 !
 ! omega is defined with respect to energy0_a
 ! Absolute omega is omega + energy0_a
 !
 do astate=nsemin,nsemax
   write(ctmp,'(i3.3)') astate
   filename = TRIM(filename_root) // '_state' // TRIM(ctmp) // '.dat'
   write(stdout,'(1x,a,a)') 'Writing selfenergy in file: ', TRIM(filename)
   open(newunit=selfenergyfile,file=filename)

   write(selfenergyfile,'(a)') '# omega (eV)             SigmaC (eV)    omega - e_gKS - Vxc + SigmaX (eV)     A (eV^-1) '

   do iomegai=1,nomegai
     write(selfenergyfile,'(20(f12.6,2x))') ( omegai(iomegai) + energy0(astate,:) )*Ha_eV,               &
                                        selfenergy_omega(iomegai,astate,:)*Ha_eV,                    &
                                        ( omegai(iomegai) - exchange_m_vxc(astate,:) )*Ha_eV,     &
                                        1.0_dp/pi/ABS( omegai(iomegai) - exchange_m_vxc(astate,:) &
                                                - selfenergy_omega(iomegai,astate,:) ) / Ha_eV
   enddo
   write(selfenergyfile,*)
   close(selfenergyfile)

 enddo


end subroutine write_selfenergy_omega


!=========================================================================
subroutine find_qp_energy_linearization(nomegai,omegai,nsemin,nsemax,selfenergy_omega,nstate,exchange_m_vxc,energy0,energy_qp_z,zz)
 implicit none

 integer,intent(in)            :: nomegai
 integer,intent(in)            :: nsemin,nsemax
 integer,intent(in)            :: nstate
 real(dp),intent(in)           :: omegai(-nomegai:nomegai)
 real(dp),intent(in)           :: selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nspin)
 real(dp),intent(in)           :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out)          :: energy_qp_z(nstate,nspin)
 real(dp),intent(out),optional :: zz(nsemin:nsemax,nspin)
!=====
 integer  :: astate,aspin
 real(dp) :: zz_a(nspin)
!=====

 ! First, a dummy initialization
 energy_qp_z(:,:) = energy0(:,:)

 ! Then overwrite the interesting energy with the calculated GW one
 do astate=nsemin,nsemax

   if( nomegai > 0 .AND. PRESENT(zz) ) then
     zz_a(:) = ( selfenergy_omega(1,astate,:) - selfenergy_omega(-1,astate,:) ) / ( omegai(1) - omegai(-1) )
     zz_a(:) = 1.0_dp / ( 1.0_dp - zz_a(:) )
     ! Contrain Z to be in [0:1] to avoid crazy values
     do aspin=1,nspin
       zz_a(aspin) = MIN( MAX(zz_a(aspin),0.0_dp) , 1.0_dp )
     enddo

     zz(astate,:)          = zz_a(:)
     energy_qp_z(astate,:) = energy0(astate,:) + zz_a(:) * ( selfenergy_omega(0,astate,:) + exchange_m_vxc(astate,:) )

   else

     energy_qp_z(astate,:) = energy0(astate,:) + selfenergy_omega(0,astate,:) + exchange_m_vxc(astate,:)

   endif

 enddo


end subroutine find_qp_energy_linearization


!=========================================================================
subroutine find_qp_energy_graphical(nomegai,omegai,nsemin,nsemax,selfenergy_omega,nstate,exchange_m_vxc,energy0,energy_qp_g)
 implicit none

 integer,intent(in)   :: nomegai
 integer,intent(in)   :: nsemin,nsemax
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: omegai(-nomegai:nomegai)
 real(dp),intent(in)  :: selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nspin)
 real(dp),intent(in)  :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out) :: energy_qp_g(nstate,nspin)
!=====
 integer  :: astate,aspin
 integer  :: info
 real(dp) :: sigma_xc_m_vxc(-nomegai:nomegai)
!=====

 ! First, a dummy initialization
 energy_qp_g(:,:) = 0.0_dp

 ! Then overwrite the interesting energy with the calculated GW one
 do astate=nsemin,nsemax

   if( MODULO(astate-nsemin,nproc_world) /= rank_world ) cycle

   do aspin=1,nspin
     sigma_xc_m_vxc(:) = selfenergy_omega(:,astate,aspin) + exchange_m_vxc(astate,aspin)
     !
     ! QP equation:
     ! E_GW - E_gKS = \Sigma_c(E_GW) + \Sigma_x - v_xc
     !
     ! Remember: omegai = E_GW - E_gKS
     energy_qp_g(astate,aspin) = find_fixed_point(nomegai,omegai,sigma_xc_m_vxc,info) + energy0(astate,aspin)
     if( info /= 0 ) then
       write(stdout,'(1x,a,i4,2x,i4)') 'Unreliable graphical solution of the QP equation for state,spin: ',astate,aspin
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
subroutine output_qp_energy(calcname,nstate,nsemin,nsemax,energy0,exchange_m_vxc,selfenergy,energy1,energy2,zz)
 implicit none
 
 character(len=*)             :: calcname
 integer                      :: nstate,nsemin,nsemax
 real(dp),intent(in)          :: energy0(nstate,nspin),exchange_m_vxc(nstate,nspin)
 real(dp),intent(in)          :: selfenergy(nsemin:nsemax,nspin),energy1(nstate,nspin)
 real(dp),intent(in),optional :: energy2(nstate,nspin),zz(nsemin:nsemax,nspin)
!=====
 integer  :: astate
!=====


 write(stdout,'(/,1x,a,1x,a)') TRIM(calcname),'eigenvalues (eV)'

 if( PRESENT(zz) .AND. PRESENT(energy2) ) then

   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,10x,a,11x,a,10x,a)') '#','E0','SigX-Vxc','SigC','Z','E_Z','E_qp'
   else
     write(stdout,'(3x,a,15x,a,22x,a,19x,a,24x,a,23x,a,23x,a)') '#','E0','SigX-Vxc','SigC','Z','E_Z','E_qp'
     write(stdout,'(12x,12(a4,9x))') ' up ','down',' up ','down',' up ','down',' up ','down',' up ','down',' up ','down'
   endif

   do astate=nsemin,nsemax

     write(stdout,'(i4,1x,20(1x,f12.6))') astate,energy0(astate,:)*Ha_eV,   &
                                          exchange_m_vxc(astate,:)*Ha_eV,   &
                                          selfenergy(astate,:)*Ha_eV,       &
                                          zz(astate,:),                     &
                                          energy1(astate,:)*Ha_eV,          &
                                          energy2(astate,:)*Ha_eV
   enddo

 else

   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,9x,a)') '#','E0','SigX-Vxc','SigC','E_qp'
   else
     write(stdout,'(3x,a,15x,a,22x,a,20x,a,22x,a)') '#','E0','SigX-Vxc','SigC','E_qp'
     write(stdout,'(12x,8(a4,9x))') ' up ','down',' up ','down',' up ','down',' up ','down'
   endif

   do astate=nsemin,nsemax

     write(stdout,'(i4,1x,20(1x,f12.6))') astate,energy0(astate,:)*Ha_eV,       &
                                          exchange_m_vxc(astate,:)*Ha_eV,       &
                                          selfenergy(astate,:)*Ha_eV,           &
                                          energy1(astate,:)*Ha_eV
   enddo



 endif

end subroutine output_qp_energy


!=========================================================================
subroutine selfenergy_set_state_ranges(nstate,occupation)
 use m_atoms
 implicit none
!=====
 integer,intent(in)  :: nstate
 real(dp),intent(in) :: occupation(nstate,nspin)
!=====
 integer :: istate
 integer :: nsemax_tmp
!=====

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


end subroutine selfenergy_set_state_ranges


!=========================================================================
subroutine setup_exchange_m_vxc_diag(basis,nstate,m_ham,n_ham, &
     occupation,c_matrix,hamiltonian_exx,hamiltonian_xc,exchange_m_vxc_diag)
 use m_inputparam
 use m_basis_set
 use m_dft_grid
 use m_hamiltonian
 use m_hamiltonian_sca
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate,m_ham,n_ham
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(inout)     :: hamiltonian_exx(m_ham,n_ham,nspin),hamiltonian_xc(m_ham,n_ham,nspin)
 real(dp),intent(out)       :: exchange_m_vxc_diag(nstate,nspin)
!=====
 integer  :: ispin,istate
 real(dp) :: exc,eexx
 real(dp),allocatable :: occupation_tmp(:,:)
 real(dp),allocatable :: p_matrix_tmp(:,:,:),p_matrix_occ(:,:),p_matrix_sqrt(:,:,:)
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
   allocate(p_matrix_sqrt(basis%nbf,basis%nbf,nspin))
   allocate(p_matrix_occ(basis%nbf,nspin))
   ! Override the occupation of the core electrons
   occupation_tmp(:,:) = occupation(:,:)
   do istate=1,dft_core
     occupation_tmp(istate,:) = 0.0_dp
   enddo

   if( calc_type%is_dft ) then
     call init_dft_grid(grid_level)
     call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation_tmp,p_matrix_tmp)
     call setup_sqrt_density_matrix(basis%nbf,p_matrix_tmp,p_matrix_sqrt,p_matrix_occ)
     call dft_exc_vxc(basis,p_matrix_occ,p_matrix_sqrt,p_matrix_tmp,hamiltonian_xc,exc)
     call destroy_dft_grid()
   endif

   if( .NOT. is_full_auxil ) then
     call setup_exchange(print_matrix_,basis%nbf,p_matrix_tmp,hamiltonian_exx,eexx)
   else
     if( parallel_ham ) then
       call setup_exchange_ri_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix_occ,p_matrix_sqrt,p_matrix_tmp,hamiltonian_exx,eexx)
     else
       call setup_exchange_ri(print_matrix_,basis%nbf,p_matrix_occ,p_matrix_sqrt,p_matrix_tmp,hamiltonian_exx,eexx)
     endif
   endif

   deallocate(occupation_tmp,p_matrix_tmp,p_matrix_sqrt,p_matrix_occ)

   hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + alpha_hybrid * hamiltonian_exx(:,:,:)
 endif



 !
 ! Calculate the diagonal of the matrix Sigma_x - Vxc
 ! for the forthcoming GW corrections
 exchange_m_vxc_diag(:,:) = 0.0_dp
 do ispin=1,nspin
   do istate=1,nstate

      exchange_m_vxc_diag(istate,ispin) =  DOT_PRODUCT(  c_matrix(:,istate,ispin) , &
                                              MATMUL( ( hamiltonian_exx(:,:,ispin) - hamiltonian_xc(:,:,ispin) ) , &
                                                        c_matrix(:,istate,ispin) ) )
   enddo
 enddo

end subroutine setup_exchange_m_vxc_diag


!=========================================================================
end module m_selfenergy_tools
!=========================================================================
