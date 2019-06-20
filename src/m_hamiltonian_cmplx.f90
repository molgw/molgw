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
 use m_basis_set
 use m_eri_calculate
 use m_density_tools

contains


!=========================================================================
subroutine setup_exchange_ri_cmplx(occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 implicit none
 real(dp),intent(in)     :: occupation(:,:)
 complex(dp),intent(in)  :: c_matrix(:,:,:)
 complex(dp),intent(in)  :: p_matrix(:,:,:)
 complex(dp),intent(out) :: exchange_ij(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer                 :: nauxil_local
 integer                 :: nbf,nstate
 integer                 :: nocc
 integer                 :: ibf,jbf,ispin,istate
 complex(dp),allocatable :: tmp_cmplx(:,:),c_t_cmplx(:,:)
 integer                 :: ipair,ipair_local,iauxil
 integer                 :: ibf_auxil_first,nbf_auxil_core
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Complex Exchange term with Resolution-of-Identity'

 exchange_ij(:,:,:) = (0.0_dp, 0.0_dp)

 ! Find highest occupied state
 nocc = get_number_occupied_states(occupation)

 nbf    = SIZE(exchange_ij,DIM=1)
 nstate = SIZE(occupation(:,:),DIM=1)

 nauxil_local = nauxil_3center

 allocate(tmp_cmplx(nocc,nbf))
 allocate(c_t_cmplx(nocc,nbf))

 do ispin=1,nspin

   !$OMP PARALLEL DO
   do ibf=1,nbf
     c_t_cmplx(:,ibf) = c_matrix(ibf,1:nocc,ispin) * SQRT( occupation(1:nocc,ispin) / spin_fact )
   enddo
   !$OMP END PARALLEL DO

   do iauxil=1,nauxil_local
     if( MODULO( iauxil - 1 , nproc_ortho ) /= rank_ortho ) cycle
     tmp_cmplx(:,:) = (0.0_dp, 0.0_dp)
     !$OMP PARALLEL PRIVATE(ibf,jbf,ipair)
     !$OMP DO REDUCTION(+:tmp_cmplx)
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp_cmplx(:,ibf) = tmp_cmplx(:,ibf) + c_t_cmplx(:,jbf) * eri_3center(ipair,iauxil)
       tmp_cmplx(:,jbf) = tmp_cmplx(:,jbf) + c_t_cmplx(:,ibf) * eri_3center(ipair,iauxil)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( CONJG(TRANSPOSE(tmp(:,:))) , tmp(:,:) ) / spin_fact
     ! C = A^H * A + C
     call ZHERK('L','C',nbf,nocc,-1.0_dp,tmp_cmplx,nocc,1.0_dp,exchange_ij(1,1,ispin),nbf)

   enddo
 enddo

 deallocate(c_t_cmplx)
 deallocate(tmp_cmplx)


 !
 ! Need to symmetrize exchange_ij
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo
 enddo
 call xsum_world(exchange_ij)

 eexchange = 0.5_dp * REAL( SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) ) , dp)

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri_cmplx


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
   dipole(:) = dipole(:) + zatom(iatom) * xatom(:,iatom)
 enddo

end subroutine static_dipole_fast_cmplx


!=========================================================================
subroutine calc_density_in_disc_cmplx_dft_grid(basis,occupation,c_matrix_cmplx,num,time_cur)
 use m_inputparam
 use m_dft_grid
 use m_atoms
 use m_density_tools
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(:,:)
 complex(dp),intent(in)     :: c_matrix_cmplx(:,:,:)
 integer,intent(in)         :: num
 real(dp),intent(in)        :: time_cur
!=====
 real(dp)             :: length
 integer              :: nstate,ndisc
 integer              :: ibf,jbf,ispin
 integer              :: idft_xc
 integer              :: igrid_start,igrid_end,ir,nr
 character(len=200)   :: file_name(2)
 real(dp)             :: z_min,z_max
 integer              :: file_out(2),igrid
 !vectors in the plane
 real(dp)             :: vec_r(3)
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: tmp_batch(:,:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: rhor_batch(:,:)
 real(dp)             :: dz_disc
 real(dp),allocatable :: charge_disc(:,:)
 real(dp),allocatable :: charge_out(:,:)
 integer(dp),allocatable :: count_z_section(:,:)
 integer              :: idisc,i_max_atom,nocc
 logical              :: file_exists
 integer              :: imanual
!=====

 call start_clock(timing_calc_dens_disc)

 nocc = SIZE(c_matrix_cmplx(:,:,:),DIM=2)

 if( excit_type%form==EXCIT_PROJECTILE ) then
   i_max_atom=natom-nprojectile
 else
   i_max_atom=natom
 endif

 inquire(file='manual_disc_dft_grid',exist=file_exists)
 if(file_exists) then
   open(newunit=imanual,file='manual_disc_dft_grid',status='old')
   read(imanual,*) ndisc
   read(imanual,*) length
   close(imanual)
 else
   ndisc=100
   length=10.0_dp
   call issue_warning('calc_density_in_disc_cmplx_dft_grid: manual file was not found')
 endif

 z_min =MIN(MINVAL( xatom(3,1:i_max_atom) ),MINVAL( xbasis(3,:) )) - length
 z_max =MAX(MAXVAL( xatom(3,1:i_max_atom) ),MAXVAL( xbasis(3,:) )) + length

 nstate = SIZE(occupation,DIM=1)

 write(stdout,*) 'Calculate electronic density in discs'

 !
 ! Loop over batches of grid points
 !

 allocate(charge_disc(ndisc,nspin))
 allocate(charge_out(2,nspin))
 allocate(count_z_section(ndisc,nspin))
 charge_disc(:,:) = 0.0_dp
 charge_out(:,:) = 0.0_dp
 count_z_section(:,:) = 0

 dz_disc=(z_max-z_min)/ndisc

 do igrid_start=1,ngrid,BATCH_SIZE
   igrid_end = MIN(ngrid,igrid_start+BATCH_SIZE-1)
   nr = igrid_end - igrid_start + 1

   allocate(weight_batch(nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(rhor_batch(nspin,nr))

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   call get_basis_functions_r_batch(basis,igrid_start,basis_function_r_batch)

   call calc_density_r_batch(occupation,c_matrix_cmplx,basis_function_r_batch,rhor_batch)

   do ir=1,nr
     igrid = igrid_start + ir - 1
     vec_r=rr_grid(1:3,igrid)
     do ispin=1,nspin
       idisc = INT((vec_r(3)-z_min)/dz_disc) + 1
       if( idisc > 0 .AND. idisc <= ndisc .AND. (vec_r(1)**2+vec_r(2)**2)**0.5_dp <= r_disc ) then
         charge_disc(idisc,ispin)=charge_disc(idisc,ispin)+rhor_batch(ispin,ir) * weight_batch(ir)
         count_z_section(idisc,ispin)=count_z_section(idisc,ispin)+1
       end if
       if( idisc <= 0 ) then
         charge_out(1,ispin)=charge_out(1,ispin)+rhor_batch(ispin,ir) * weight_batch(ir)
       end if
       if( idisc > ndisc ) then
         charge_out(2,ispin)=charge_out(2,ispin)+rhor_batch(ispin,ir) * weight_batch(ir)
       end if
     enddo ! loop on the ispin
   enddo ! loop on ir

   deallocate(weight_batch)
   deallocate(basis_function_r_batch)
   deallocate(rhor_batch)

 enddo ! loop on the batches
 call xsum_grid(charge_disc(:,:))

 if( is_iomaster ) then

   do ispin=1,nspin
     write(file_name(ispin),'(a,i4.4,a,i1,a,i3.3,f0.3,a)') 'disc_dens_',num, "_s_",ispin,"_r_",INT(r_disc),r_disc-INT(r_disc),".dat"
     open(newunit=file_out(ispin),file=file_name(ispin))
   enddo

   do ispin=1,nspin
     write(file_out(ispin),'(a,F12.6,a,3F12.6)') '# Time: ',time_cur, '  Projectile position (A): ',xatom(:,natom+nghost)*bohr_A
     do idisc=1,ndisc
       write(file_out(ispin),'(F16.4,F19.10,i6)') (z_min+idisc*dz_disc)*bohr_A,charge_disc(idisc,ispin),count_z_section(idisc,ispin)
     end do
     close(file_out(ispin))
   end do
   write(stdout,'(A,2F12.6)') "Charge out of region, left and right:",charge_out(:,1)

 end if
 !
 ! Sum up the contributions from all procs only if needed

 deallocate(charge_disc)

 call stop_clock(timing_calc_dens_disc)

end subroutine calc_density_in_disc_cmplx_dft_grid


end module m_hamiltonian_cmplx
!=========================================================================
