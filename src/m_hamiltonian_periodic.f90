!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval
!
! This module contains
!
!=========================================================================
#include "molgw.h"
#if !defined(NO_LIBINT)
#include<libint2/libint2_params.h>
#endif

module m_hamiltonian_periodic
  use m_definitions
  use m_timing
  use m_mpi
  use m_scalapack
  use m_warning
  use m_memory
  use m_cart_to_pure
  use m_inputparam
  use m_basis_set
  use m_libint_tools
  use m_libcint_tools
  use m_io
  use m_atoms
  use m_dft_grid
  use m_ecp
  use m_eri
  use m_hamiltonian_tools
#if defined(HAVE_FFTW3)
  use m_fftw3
#endif



contains


!=========================================================================
subroutine setup_hartree_periodic(basis, p_matrix, hartree_ao, ehartree)
  implicit none

  type(basis_set), intent(in) :: basis
  class(*), intent(in)  :: p_matrix(:, :, :)
  real(dp), intent(out) :: hartree_ao(:, :)
  real(dp), intent(out) :: ehartree
  !=====
  integer              :: nbf
  integer              :: ibf, kbf, lbf
  integer              :: ipair
  real(dp), allocatable :: x_vector(:)
  integer              :: timing_xxdft_hartree
  real(dp), allocatable :: pmat(:)

  integer(C_INT), parameter :: ngx=60
  integer(C_INT), parameter :: ngy=60
  integer(C_INT), parameter :: ngz=60
  integer, parameter :: nr = ngx * ngy * ngz
  type(C_PTR) :: plan
  type(C_PTR) :: pr, pg, pvr, pvg
  real(C_DOUBLE), pointer :: rhor(:, :, :), vr(:, :, :)
  complex(C_DOUBLE_COMPLEX), pointer ::  rhog(:, :, :), vg(:, :, :)
  integer :: igx, igy, igz, ix, iy, iz, ir
  real(dp) :: rr(3), rgrid(3, nr), rhogrid(nr, nspin)
  real(dp) :: dr(3, 3)
  real(dp), allocatable :: c_matrix(:, :, :), occupation(:, :)
  !=====

#if !defined(HAVE_FFTW3)
  call die('setup_hartree_periodic: requires FFTW at compilation time')
#else

  nbf = SIZE(hartree_ao(:, :), DIM=1)

  select type(p_matrix)
  type is(real(dp))
    timing_xxdft_hartree   = timing_hartree
  type is(complex(dp))
    timing_xxdft_hartree   = timing_tddft_hartree
  class default
    call die("setup_hartree_periodic: p_matrix is neither real nor complex")
  end select

  call start_clock(timing_xxdft_hartree)

  write(stdout, *) 'Calculate periodic Hartree term'

  ir = 0
  do iz=0, ngz - 1
    do iy=0, ngy - 1
      do ix=0, ngx - 1
        ir = ir + 1
        rgrid(:, ir) = MATMUL( aprim(:, :), [DBLE(ix)/ngx, DBLE(iy)/ngy, DBLE(iz)/ngz] )
      enddo
    enddo
  enddo

  select type(p_matrix)
  type is(real(dp))
    call get_c_matrix_from_p_matrix(p_matrix, c_matrix, occupation)
  end select
  call start_clock(timing_tmp1)
  call calculate_density_periodic(basis, c_matrix, occupation, rgrid, rhogrid)
  call stop_clock(timing_tmp1)
  call output_timing_line('tmp1' , timing_tmp1 , 0)

  
  dr(:, 1) = aprim(:, 1) / ngx
  dr(:, 2) = aprim(:, 2) / ngy
  dr(:, 3) = aprim(:, 3) / ngz
  call write_cube_file('toto.cube', ngx, ngy, ngz, dr, rhogrid(:, 1), comment='test')

  
  do ir=1, nr 
    write(1001, '(1x,i8,1x,3(f12.6,1x),4x,f16.6)') ir, rgrid(:, ir), rhogrid(ir, 1)
  enddo


  ! Preparing calls to FFT
  pr = fftw_alloc_real(int(ngx * ngy * ngz, C_SIZE_T))
  call C_F_POINTER(pr, rhor, [ngx, ngy, ngz])
  pg = fftw_alloc_complex(int( (ngx/2+1) * ngy * ngz, C_SIZE_T))
  call C_F_POINTER(pg, rhog, [(ngx/2+1), ngy, ngz])

  rhor(:, :, :) = RESHAPE(rhogrid(:, 1), [ngx, ngy, ngz])
  ! reversed order ngz, ngy, ngx
  plan = fftw_plan_dft_r2c_3d(ngz, ngy, ngx, rhor, rhog, FFTW_ESTIMATE)
  call fftw_execute_dft_r2c(plan, rhor, rhog)
  call fftw_destroy_plan(plan)

  rhog(:, :, :) = rhog(:, :, :) / REAL(nr, KIND=dp)

  write(stdout, *) 'number of electrons', SUM(rhor) * volume / REAL(nr, KIND=dp)
  write(stdout, *) 'number of electrons', rhog(1, 1, 1)%re * volume


  pvg = fftw_alloc_complex(int( (ngx/2+1) * ngy * ngz, C_SIZE_T))
  call C_F_POINTER(pvg, vg, [(ngx/2+1), ngy, ngz])

  do igz=1, ngz
    do igy=1, ngy
      do igx=1, ngx/2+1
        write(1000,'(3(i3,1x),2(1x,es16.6))') igx, igy, igz, rhog(igx, igy, igz)
        if( igx * igy * igz > 1 ) then
          vg(igx, igy, igz) = vcoul(igx, igy, igz) * rhog(igx, igy, igz)
        else
          vg(igx, igy, igz) = 0.0_dp
        endif
      enddo
    enddo
  enddo

  pvr = fftw_alloc_real(int(ngx * ngy * ngz, C_SIZE_T))
  call C_F_POINTER(pvr, vr, [ngx, ngy, ngz])

  plan = fftw_plan_dft_c2r_3d(ngz, ngy, ngx, vg, vr, FFTW_ESTIMATE)
  call fftw_execute_dft_c2r(plan, vg, vr)
  call fftw_destroy_plan(plan)

  write(stdout,*) 'FBFB Eh', 0.5d0 * SUM( vr * rhor ) * volume / REAL(nr, KIND=dp)

  rhogrid(:, :) = RESHAPE(vr(:, :, :), [nr, 1])
  call calculate_hao_periodic(basis, rgrid, rhogrid, hartree_ao)

  call fftw_free(pr)
  call fftw_free(pg)
  call fftw_free(pvr)
  call fftw_free(pvg)

  call dump_out_matrix(.TRUE., '=== Hartree contribution (FFTW) ===', hartree_ao)


  hartree_ao(:, :) = 0.0_dp

  allocate(pmat(npair))
  allocate(x_vector(nauxil_local))

  select type(p_matrix)
  type is(real(dp))
    do ipair=1, npair
      kbf = index_basis(1, ipair)
      lbf = index_basis(2, ipair)
      pmat(ipair) = SUM(p_matrix(kbf, lbf, :)) * 2.0_dp
    enddo
  type is(complex(dp))
    do ipair=1, npair
      kbf = index_basis(1, ipair)
      lbf = index_basis(2, ipair)
      ! As all pairs contribute twice for (k, l) and (l, k) and as P is hermitian,
      ! only the real part survives
      pmat(ipair) = SUM(p_matrix(kbf, lbf, :)%re) * 2.0_dp
    enddo
  end select

  ! X_P = \sum_{\alpha \beta} P_{\alpha \beta} * ( \alpha \beta | P )
  call DGEMV('T', npair, nauxil_local, 1.0d0, eri_3center, npair, pmat, 1, 0.0d0, x_vector, 1)
  ! v_H_{alpha beta} = \sum_P ( alpha beta | P ) * X_P
  call DGEMV('N', npair, nauxil_local, 1.0d0, eri_3center, npair, x_vector, 1, 0.0d0, pmat, 1)

  !$OMP PARALLEL PRIVATE(kbf, lbf)
  !$OMP DO
  do ipair=1, npair
    kbf = index_basis(1, ipair)
    lbf = index_basis(2, ipair)
    hartree_ao(kbf, lbf) = pmat(ipair)
    hartree_ao(lbf, kbf) = pmat(ipair)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Do not forget that the eri_3center(ibf, ibf | P ) included a factor 0.50
  do ibf=1, nbf
    hartree_ao(ibf, ibf) = hartree_ao(ibf, ibf) * 2.0_dp
  enddo
  deallocate(x_vector, pmat)

  !
  ! Sum up the different contribution from different procs
  call world%sum(hartree_ao)
  hartree_ao(:, :) = hartree_ao(:, :) / REAL(poorman%nproc, dp)

  call dump_out_matrix(.TRUE., '=== Hartree contribution (ref) ===', hartree_ao)

  select type(p_matrix)
  type is(real(dp))
    ehartree = 0.5_dp * SUM( hartree_ao(:, :) * SUM(p_matrix(:, :, :), DIM=3) )
  type is(complex(dp))
    ehartree = 0.5_dp * SUM( hartree_ao(:, :) * SUM(REAL(p_matrix(:, :, :), dp), DIM=3) )
  end select

  write(stdout,*) 'Eh 2 gaussian', ehartree

  call stop_clock(timing_xxdft_hartree)
#endif

  call die("STOP: ENOUGH")

contains

pure function vcoul(igx, igy, igz) 
  implicit none
  integer, intent(in) :: igx, igy, igz
  real(dp) :: vcoul
  !=====
  integer  :: igxm, igym, igzm
  real(dp) :: g2
  !=====
  
  igxm = MERGE( igx - 1, igx - 1 - ngx, igx <= ngx/2)
  igym = MERGE( igy - 1, igy - 1 - ngy, igy <= ngy/2)
  igzm = MERGE( igz - 1, igz - 1 - ngz, igz <= ngz/2)
  g2 = SUM( MATMUL(bprim(:, :) , [igxm, igym ,igzm])**2 )

  vcoul = 4.0_dp * pi / g2  !* (1.0 -  COS( SQRT(g2) * 37.0 ))

end function vcoul

end subroutine setup_hartree_periodic


!=========================================================================
subroutine calculate_density_periodic(basis, c_matrix, occupation, rr, rhor)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in)  :: c_matrix(:, :, :)
  real(dp), intent(in)  :: occupation(:, :)
  real(dp), intent(in)  :: rr(:, :)
  real(dp), intent(out) :: rhor(:, :)
  !=====
  integer :: nx=1
  integer :: ispin, ir, nr, nstate, ix, iy, iz
  real(dp) :: shift(3)
  real(dp), allocatable :: rr_shifted(:, :)
  real(dp), allocatable :: bfr(:, :)
  real(dp), allocatable :: phir(:, :)
  !=====

  nr = SIZE(rr, DIM=2)
  nstate = SIZE(c_matrix, DIM=2)

  allocate(rr_shifted, MOLD=rr)
  allocate(bfr(basis%nbf,nr))
  allocate(phir(nr, nstate))

  rhor(:, :) = 0.0_dp
  do iz=-nx,nx
    do iy=-nx,nx
      do ix=-nx,nx
        shift(:) = MATMUL( aprim(:, :), [ix, iy, iz] )

        do ir=1, nr
          rr_shifted(:, ir) = rr(:, ir) - shift(:)
          call calculate_basis_functions_r(basis, rr_shifted(:, ir), bfr(:, ir))
        enddo

        do ispin=1, nspin
          !phir(:, :) = MATMUL( TRANPOSE(bfr(:, :)) , c_matrix(:, :, ispin) )
          call DGEMM('T', 'N', nr, nstate, basis%nbf, 1.0d0, bfr, basis%nbf, c_matrix(:,:,ispin), basis%nbf, 0.0d0, phir, nr)

          rhor(:, :) = rhor(:, :) + MATMUL( phir(:, :)**2, occupation(:, :) )
        enddo

      enddo
    enddo
  enddo

  deallocate(rr_shifted)
  deallocate(bfr, phir)

end subroutine calculate_density_periodic


!=========================================================================
subroutine calculate_hao_periodic(basis, rr, vr, h_ao)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in)  :: rr(:, :)
  real(dp), intent(in) :: vr(:, :)
  real(dp), intent(out) :: h_ao(:, :)
  !=====
  integer :: nx=1
  integer :: ispin, ir, nr, nstate, ix, iy, iz, ibf, jbf
  real(dp) :: shift(3)
  real(dp), allocatable :: rr_shifted(:, :)
  real(dp), allocatable :: bfr(:, :)
  !=====

  nr = SIZE(rr, DIM=2)

  allocate(rr_shifted, MOLD=rr)
  allocate(bfr(basis%nbf,nr))

  h_ao(:, :) = 0.0_dp
  do iz=-nx,nx
    do iy=-nx,nx
      do ix=-nx,nx
        shift(:) = MATMUL( aprim(:, :), [ix, iy, iz] )

        do ir=1, nr
          rr_shifted(:, ir) = rr(:, ir) - shift(:)
          call calculate_basis_functions_r(basis, rr_shifted(:, ir), bfr(:, ir))
        enddo

        do ispin=1, nspin
          do jbf=1, basis%nbf
            do ibf=1, basis%nbf
              h_ao(ibf, jbf) = h_ao(ibf, jbf) + SUM( bfr(ibf, :) * vr(:, ispin) * bfr(jbf, :) )
            enddo
          enddo
        enddo

      enddo
    enddo
  enddo
  h_ao(:, :) = h_ao(:, :) * volume / nr

  deallocate(rr_shifted)
  deallocate(bfr)

end subroutine calculate_hao_periodic


end module m_hamiltonian_periodic


!=========================================================================
