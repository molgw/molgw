!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains the calculation of the vertex functions 
! within different flavors: TDHF
!
!=========================================================================
#include "molgw.h"
module m_tdhf_selfenergy
  use m_definitions
  use m_mpi
  use m_timing
  use m_inputparam
  use m_warning
  use m_basis_set
  use m_spectral_function
  use m_eri_ao_mo
  use m_selfenergy_tools
  use m_tddft_fxc
  use m_linear_response


contains


!              call DGEMM('N','N',nauxil_global,nvirtual_G-nhomo_G-1,nauxil_global, &
!                         1.0_dp,chi_wwp(:,:),nauxil_global, &
!                                eri3_a_m,nauxil_global, &
!                         0.0_dp,tmp,nauxil_global)
!              call DGEMM('T','N',nvirtual_G-ncore_G-1,nvirtual_G-nhomo_G-1,nauxil_global, &
!                         1.0_dp,eri3_r_i(:,:),nauxil_global, &
!                                tmp(:,:),nauxil_global, &
!                         0.0_dp,braket2_ra(:,:),nvirtual_G-ncore_G-1)

!=========================================================================
!
!
subroutine tdhf_selfenergy(basis,occupation,energy,c_matrix,se)
  implicit none

  type(basis_set)                    :: basis
  real(dp),intent(in)                :: occupation(:,:),energy(:,:)
  real(dp),intent(in)                :: c_matrix(:,:,:)
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer                 :: nstate,nmat,imat
  integer                 :: istate,astate,jstate,bstate,iaspin,spole
  integer                 :: pstate,iomega_sigma
  real(dp),allocatable    :: x_matrix(:,:),y_matrix(:,:)
  real(dp),allocatable    :: a_matrix(:,:),b_matrix(:,:)
  real(dp)                :: erpa_tmp,egw_tmp
  type(spectral_function) :: wpol
  complex(dp),allocatable :: sigma_tdhf(:,:,:)
  real(dp) :: eri_iajp,eri_ijap,eri_ajip
  real(dp) :: eri_iabp,eri_ibap,eri_abip
  real(dp),allocatable :: eri_tmp1o(:,:), eri_tmp2o(:,:), eri_tmp3o(:,:)
  real(dp),allocatable :: eri_tmp1v(:,:), eri_tmp2v(:,:), eri_tmp3v(:,:)
  real(dp),allocatable :: num_tmp1o(:,:), num_tmp2o(:,:)
  real(dp),allocatable :: num_tmp1v(:,:), num_tmp2v(:,:)
  real(dp),allocatable :: xpy_matrix(:,:)
  real(dp) :: num1,num2
  real(dp),allocatable    :: axpby_matrix(:,:),aypbx_matrix(:,:)
  !=====

  call start_clock(timing_gw_self)

  if( .NOT. has_auxil_basis ) call die('tdhf_selfenergy: not implemented without an auxiliary basis')
  if( nspin > 1 ) call die('tdhf_selfenergy: not implemented for spin unrestricted')

  nstate = SIZE(energy,DIM=1)

  write(stdout,'(/,1x,a)') 'Calculate Sigma_TDHF (Bruneval-Foerster formula'



  call wpol%init(nstate,occupation,0)
  nmat = wpol%npole_reso

  call clean_allocate('X matrix',x_matrix,nmat,nmat)
  call clean_allocate('Y matrix',y_matrix,nmat,nmat)
  call clean_allocate('A matrix',a_matrix,nmat,nmat)
  call clean_allocate('B matrix',b_matrix,nmat,nmat)

  ! Get A and B, X and Y
  call polarizability(.FALSE.,.TRUE.,basis,occupation,energy,c_matrix,erpa_tmp,egw_tmp,wpol, &
                       a_matrix=a_matrix,b_matrix=b_matrix,x_matrix=x_matrix,y_matrix=y_matrix)

  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)

  ! Enforce PT2 manually
  if( .FALSE. ) then
    call issue_warning("Enforce PT2 for debug")
    x_matrix(:,:) = 0.0d0
    y_matrix(:,:) = 0.0d0
    do imat=1,nmat
      spole = imat
      istate = wpol%transition_table(1,imat)
      astate = wpol%transition_table(2,imat)
      wpol%pole(spole) = energy(astate,1) - energy(istate,1)
      x_matrix(imat,spole) = 1.0d0
    enddo
  endif
  if( .FALSE. ) then
    call issue_warning("Print X for debug")
    
    write(101,*) "# i, a, spole, x_matrix(ia,spole),y_matrix(ia,spole)"
    do spole=1,nmat
      do imat=1,nmat
        istate = wpol%transition_table(1,imat)
        astate = wpol%transition_table(2,imat)
        write(101,*) istate,astate,spole,x_matrix(imat,spole),y_matrix(imat,spole)
      enddo

      write(*,*) '1?',spole,SUM( x_matrix(:,spole) * x_matrix(:,spole) - y_matrix(:,spole) * y_matrix(:,spole) )
      write(*,*) '0?',spole,SUM( -y_matrix(:,spole) * x_matrix(:,spole) + x_matrix(:,spole) * y_matrix(:,spole) )
    enddo
    allocate(xpy_matrix(nmat,nmat))
    xpy_matrix = MATMUL( TRANSPOSE(x_matrix), x_matrix ) -  MATMUL( TRANSPOSE(y_matrix), y_matrix )
    call dump_out_matrix(.TRUE.,'===  X**T * X - Y**T * Y ===',xpy_matrix)
    xpy_matrix = MATMUL(-TRANSPOSE(y_matrix), x_matrix ) +  MATMUL( TRANSPOSE(x_matrix), y_matrix )
    call dump_out_matrix(.TRUE.,'=== -Y**T * X + x**T * Y ===',xpy_matrix)
    deallocate(xpy_matrix)
  endif

  allocate(sigma_tdhf(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(xpy_matrix(nmat,nmat))
  xpy_matrix(:,:) = x_matrix(:,:) + y_matrix(:,:)
  allocate(eri_tmp1o(nmat,ncore_G+1:nhomo_G))
  allocate(eri_tmp2o(nmat,ncore_G+1:nhomo_G))
  allocate(eri_tmp3o(nmat,ncore_G+1:nhomo_G))
  allocate(eri_tmp1v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(eri_tmp2v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(eri_tmp3v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(num_tmp1o(nmat,ncore_G+1:nhomo_G))
  allocate(num_tmp2o(nmat,ncore_G+1:nhomo_G))
  allocate(num_tmp1v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(num_tmp2v(nmat,nhomo_G+1:nvirtual_G-1))
  sigma_tdhf(:,:,:) = 0.0_dp

  do pstate=nsemin,nsemax

    !
    ! Ocuppied states
    !
    do jstate=ncore_G+1,nhomo_G
      do imat=1,nmat
        istate = wpol%transition_table(1,imat)
        astate = wpol%transition_table(2,imat)

        eri_tmp1o(imat,jstate) = eri_eigen(istate,astate,1,jstate,pstate,1)
        eri_tmp2o(imat,jstate) = eri_eigen(istate,jstate,1,astate,pstate,1)  ! set to zero to recover GW
        eri_tmp3o(imat,jstate) = eri_eigen(astate,jstate,1,istate,pstate,1)  ! set to zero to recover GW
      enddo
    enddo

    ! Fortran version
    !num_tmp2o(:,:) = MATMUL( TRANSPOSE(xpy_matrix(:,:)) , eri_tmp1o(:,:) )
    !num_tmp1o(:,:) = 2.0 * num_tmp2o(:,:) &
    !                - MATMUL( TRANSPOSE(x_matrix(:,:)) , eri_tmp3o(:,:) ) & ! Here the role of X and Y is swapped
    !                - MATMUL( TRANSPOSE(y_matrix(:,:)) , eri_tmp2o(:,:) )   ! as compared to Vacondio

    ! BLAS version
    call DGEMM('T','N',nmat,nhomo_G-ncore_G,nmat, &
                   1.0d0,xpy_matrix(:,:),nmat,eri_tmp1o(:,:),nmat,&
                   0.0d0,num_tmp2o(:,:),nmat)

    num_tmp1o(:,:) = 2.0 * num_tmp2o(:,:)
    call DGEMM('T','N',nmat,nhomo_G-ncore_G,nmat, &
                  -1.0d0,x_matrix(:,:),nmat,eri_tmp3o(:,:),nmat,&
                   1.0d0,num_tmp1o(:,:),nmat)
    call DGEMM('T','N',nmat,nhomo_G-ncore_G,nmat, &
                  -1.0d0,y_matrix(:,:),nmat,eri_tmp2o(:,:),nmat,&
                   1.0d0,num_tmp1o(:,:),nmat)


    do jstate=ncore_G+1,nhomo_G
      do spole=1,wpol%npole_reso
        sigma_tdhf(:,pstate,1) = sigma_tdhf(:,pstate,1) &
                       +  num_tmp1o(spole,jstate) * num_tmp2o(spole,jstate) &
                          / ( se%omega(:) + se%energy0(pstate,1) - energy(jstate,1) + wpol%pole(spole) - ieta )

      enddo ! loop over spole
    enddo ! loop over jstate

    !
    ! Virtual states
    !
    do bstate=nhomo_G+1,nvirtual_G-1
      do imat=1,nmat
        istate = wpol%transition_table(1,imat)
        astate = wpol%transition_table(2,imat)

        eri_tmp1v(imat,bstate) = eri_eigen(istate,astate,1,bstate,pstate,1)
        eri_tmp2v(imat,bstate) = eri_eigen(istate,bstate,1,astate,pstate,1)  ! set to zero to recover GW
        eri_tmp3v(imat,bstate) = eri_eigen(astate,bstate,1,istate,pstate,1)  ! set to zero to recover GW
      enddo
    enddo

    ! Fortran version
    !num_tmp2v(:,:) = MATMUL( TRANSPOSE(xpy_matrix(:,:)) , eri_tmp1v(:,:) )
    !num_tmp1v(:,:) = 2.0 * num_tmp2v(:,:) &
    !                - MATMUL( TRANSPOSE(x_matrix(:,:)) , eri_tmp2v(:,:) ) & ! Here the role of X and Y is *conserved*
    !                - MATMUL( TRANSPOSE(y_matrix(:,:)) , eri_tmp3v(:,:) )   ! as compared to Vacondio

    ! BLAS version
    call DGEMM('T','N',nmat,nvirtual_G-nhomo_G-1,nmat, &
                   1.0d0,xpy_matrix(:,:),nmat,eri_tmp1v(:,:),nmat,&
                   0.0d0,num_tmp2v(:,:),nmat)
    num_tmp1v(:,:) = 2.0 * num_tmp2v(:,:)
    call DGEMM('T','N',nmat,nvirtual_G-nhomo_G-1,nmat, &
                  -1.0d0,x_matrix(:,:),nmat,eri_tmp2v(:,:),nmat,&
                   1.0d0,num_tmp1v(:,:),nmat)
    call DGEMM('T','N',nmat,nvirtual_G-nhomo_G-1,nmat, &
                  -1.0d0,y_matrix(:,:),nmat,eri_tmp3v(:,:),nmat,&
                   1.0d0,num_tmp1v(:,:),nmat)

    do bstate=nhomo_G+1,nvirtual_G-1
      do spole=1,wpol%npole_reso
        sigma_tdhf(:,pstate,1) = sigma_tdhf(:,pstate,1) &
                       +  num_tmp1v(spole,bstate) * num_tmp2v(spole,bstate) &
                          / ( se%omega(:) + se%energy0(pstate,1) - energy(bstate,1) - wpol%pole(spole) + ieta )
      enddo ! loop over spole
    enddo ! loop over bstate


  enddo ! loop over pstate

  se%sigma(:,:,:) = sigma_tdhf(:,:,:)

  call clean_deallocate('X matrix',x_matrix)
  call clean_deallocate('Y matrix',y_matrix)

  call clean_deallocate('A matrix',a_matrix)
  call clean_deallocate('B matrix',b_matrix)
  call clean_deallocate('AX+BY matrix',axpby_matrix)
  call clean_deallocate('AY+BX matrix',aypbx_matrix)

  deallocate(sigma_tdhf)
  call destroy_eri_3center_eigen()
  call wpol%destroy()

  call stop_clock(timing_gw_self)

end subroutine tdhf_selfenergy


!=========================================================================
! Vacondio et al.'s expression as I interprete it
!
subroutine tdhf_vacondio_selfenergy(basis,occupation,energy,c_matrix,se)
  implicit none

  type(basis_set)                    :: basis
  real(dp),intent(in)                :: occupation(:,:),energy(:,:)
  real(dp),intent(in)                :: c_matrix(:,:,:)
  type(selfenergy_grid),intent(inout) :: se
  !=====
  integer                 :: nstate,nmat,imat
  integer                 :: istate,astate,jstate,bstate,iaspin,spole
  integer                 :: pstate,iomega_sigma
  real(dp),allocatable    :: x_matrix(:,:),y_matrix(:,:)
  real(dp),allocatable    :: a_matrix(:,:),b_matrix(:,:)
  real(dp)                :: erpa_tmp,egw_tmp
  type(spectral_function) :: wpol
  complex(dp),allocatable :: sigma_tdhf(:,:,:)
  real(dp) :: eri_iajp,eri_ijap,eri_ajip
  real(dp) :: eri_iabp,eri_ibap,eri_abip
  real(dp),allocatable :: eri_tmp1(:,:), eri_tmp2(:,:), eri_tmp3(:,:)
  real(dp),allocatable :: eri_tmp1v(:,:), eri_tmp2v(:,:), eri_tmp3v(:,:)
  real(dp),allocatable :: num_tmp1(:,:), num_tmp2(:,:)
  real(dp),allocatable :: num_tmp1v(:,:), num_tmp2v(:,:)
  real(dp),allocatable :: xpy_matrix(:,:)
  real(dp) :: num1,num2
  real(dp),allocatable    :: axpby_matrix(:,:),aypbx_matrix(:,:)
  !=====

  call start_clock(timing_gw_self)

  if( .NOT. has_auxil_basis ) call die('tdhf_selfenergy: not implemented without an auxiliary basis')
  if( nspin > 1 ) call die('tdhf_selfenergy: not implemented for spin unrestricted')

  nstate = SIZE(energy,DIM=1)

  write(stdout,'(/,1x,a)') 'Calculate Sigma_TDHF (Bruneval-Foerster formula'

  call wpol%init(nstate,occupation,0)
  nmat = wpol%npole_reso

  call clean_allocate('X matrix',x_matrix,nmat,nmat)
  call clean_allocate('Y matrix',y_matrix,nmat,nmat)
  call clean_allocate('A matrix',a_matrix,nmat,nmat)
  call clean_allocate('B matrix',b_matrix,nmat,nmat)

  ! Get A and B, X and Y
  call polarizability(.FALSE.,.TRUE.,basis,occupation,energy,c_matrix,erpa_tmp,egw_tmp,wpol, &
                       a_matrix=a_matrix,b_matrix=b_matrix,x_matrix=x_matrix,y_matrix=y_matrix)

  !!
  !! A -> A' in Kresse's notation
  !! Remove the energy difference on the diagonal
  !call remove_a_energy_diag(energy,wpol,a_matrix)


  !call clean_allocate('AX+BY matrix',axpby_matrix,nmat,nmat)
  !call clean_allocate('AY+BX matrix',aypbx_matrix,nmat,nmat)
  !! AX + BY
  !axpby_matrix(:,:) = MATMUL( a_matrix, x_matrix ) + MATMUL( b_matrix, y_matrix )
  !! AY + BX
  !aypbx_matrix(:,:) = MATMUL( a_matrix, y_matrix ) + MATMUL( b_matrix, x_matrix )

  call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)

  allocate(sigma_tdhf(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(xpy_matrix(nmat,nmat))
  xpy_matrix(:,:) = x_matrix(:,:) + y_matrix(:,:)
  allocate(eri_tmp1(nmat,ncore_G+1:nhomo_G))
  allocate(eri_tmp2(nmat,ncore_G+1:nhomo_G))
  allocate(eri_tmp3(nmat,ncore_G+1:nhomo_G))
  allocate(eri_tmp1v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(eri_tmp2v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(eri_tmp3v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(num_tmp1(nmat,ncore_G+1:nhomo_G))
  allocate(num_tmp2(nmat,ncore_G+1:nhomo_G))
  allocate(num_tmp1v(nmat,nhomo_G+1:nvirtual_G-1))
  allocate(num_tmp2v(nmat,nhomo_G+1:nvirtual_G-1))
  sigma_tdhf(:,:,:) = 0.0_dp

  do pstate=nsemin,nsemax

    do jstate=ncore_G+1,nhomo_G
      do imat=1,nmat
        istate = wpol%transition_table(1,imat)
        astate = wpol%transition_table(2,imat)

        eri_tmp1(imat,jstate) = eri_eigen(istate,astate,1,jstate,pstate,1)
        eri_tmp2(imat,jstate) = eri_eigen(istate,jstate,1,astate,pstate,1)  ! set to zero to recover GW
        eri_tmp3(imat,jstate) = eri_eigen(astate,jstate,1,istate,pstate,1)  ! set to zero to recover GW
     enddo
   enddo

   num_tmp2(:,:) = MATMUL( TRANSPOSE(xpy_matrix(:,:)) , eri_tmp1(:,:) )

   num_tmp1(:,:) = 2.0 * num_tmp2(:,:) &
                   - MATMUL( TRANSPOSE(x_matrix(:,:)) , eri_tmp2(:,:) ) &
                   - MATMUL( TRANSPOSE(y_matrix(:,:)) , eri_tmp3(:,:) )


  !    do spole=1,wpol%npole_reso

  !      num1 = 0.0
  !      num2 = 0.0
  !      do imat=1,nmat
  !        istate = wpol%transition_table(1,imat)
  !        astate = wpol%transition_table(2,imat)

  !        eri_iajp = eri_eigen(istate,astate,1,jstate,pstate,1)
  !        eri_ijap = eri_eigen(istate,jstate,1,astate,pstate,1)  ! set to zero to recover GW
  !        eri_ajip = eri_eigen(astate,jstate,1,istate,pstate,1)  ! set to zero to recover GW

  !        num1 = num1 + 2.0 * x_matrix(imat,spole) * eri_iajp - x_matrix(imat,spole) * eri_ijap  &
  !                    + 2.0 * y_matrix(imat,spole) * eri_iajp - y_matrix(imat,spole) * eri_ajip 
  !        num2 = num2 + ( x_matrix(imat,spole) + y_matrix(imat,spole) ) * eri_iajp

  !      enddo

    do jstate=ncore_G+1,nhomo_G
      do spole=1,wpol%npole_reso
        sigma_tdhf(:,pstate,1) = sigma_tdhf(:,pstate,1) &
                       +  num_tmp1(spole,jstate) * num_tmp2(spole,jstate) &
                          / ( se%omega(:) + se%energy0(pstate,1) - energy(jstate,1) + wpol%pole(spole) -ieta )

      enddo ! loop over spole
    enddo ! loop over jstate

    do bstate=nhomo_G+1,nvirtual_G-1
      do imat=1,nmat
        istate = wpol%transition_table(1,imat)
        astate = wpol%transition_table(2,imat)

        eri_tmp1v(imat,bstate) = eri_eigen(istate,astate,1,bstate,pstate,1)
        eri_tmp2v(imat,bstate) = eri_eigen(istate,bstate,1,astate,pstate,1)  ! set to zero to recover GW
        eri_tmp3v(imat,bstate) = eri_eigen(astate,bstate,1,istate,pstate,1)  ! set to zero to recover GW
     enddo
   enddo
   num_tmp2v(:,:) = MATMUL( TRANSPOSE(xpy_matrix(:,:)) , eri_tmp1v(:,:) )
   num_tmp1v(:,:) = 2.0 * num_tmp2v(:,:) &
                   - MATMUL( TRANSPOSE(x_matrix(:,:)) , eri_tmp2v(:,:) ) &
                   - MATMUL( TRANSPOSE(y_matrix(:,:)) , eri_tmp3v(:,:) )
    do bstate=nhomo_G+1,nvirtual_G-1
      do spole=1,wpol%npole_reso
        sigma_tdhf(:,pstate,1) = sigma_tdhf(:,pstate,1) &
                       +  num_tmp1v(spole,bstate) * num_tmp2v(spole,bstate) &
                          / ( se%omega(:) + se%energy0(pstate,1) - energy(bstate,1) - wpol%pole(spole) +ieta )
      enddo ! loop over spole
    enddo ! loop over bstate

    !do bstate=nhomo_G+1,nvirtual_G-1
    !  do spole=1,wpol%npole_reso

    !    num1 = 0.0
    !    num2 = 0.0
    !    do imat=1,nmat
    !      istate = wpol%transition_table(1,imat)
    !      astate = wpol%transition_table(2,imat)

    !      eri_iabp = eri_eigen(istate,astate,1,bstate,pstate,1)
    !      eri_ibap = eri_eigen(istate,bstate,1,astate,pstate,1)  ! set to zero to recover GW
    !      eri_abip = eri_eigen(astate,bstate,1,istate,pstate,1)  ! set to zero to recover GW

    !      num1 = num1 + 2.0 * x_matrix(imat,spole) * eri_iabp - x_matrix(imat,spole) * eri_ibap  &
    !                  + 2.0 * y_matrix(imat,spole) * eri_iabp - y_matrix(imat,spole) * eri_abip
    !      num2 = num2 + ( x_matrix(imat,spole) + y_matrix(imat,spole) ) * eri_iabp

    !    enddo

    !    sigma_tdhf(:,pstate,1) = sigma_tdhf(:,pstate,1) &
    !                   +  num1 * num2 &
    !                      / ( se%omega(:) + se%energy0(pstate,1) - energy(bstate,1) - wpol%pole(spole) +ieta )


    !  enddo ! loop over spole
    !enddo ! loop over bstate

  enddo ! loop over pstate

  se%sigma(:,:,:) = sigma_tdhf(:,:,:)

  call clean_deallocate('X matrix',x_matrix)
  call clean_deallocate('Y matrix',y_matrix)

  call clean_deallocate('A matrix',a_matrix)
  call clean_deallocate('B matrix',b_matrix)
  call clean_deallocate('AX+BY matrix',axpby_matrix)
  call clean_deallocate('AY+BX matrix',aypbx_matrix)

  deallocate(sigma_tdhf)
  call destroy_eri_3center_eigen()
  call wpol%destroy()

  call stop_clock(timing_gw_self)

end subroutine tdhf_vacondio_selfenergy


end module m_tdhf_selfenergy
!=========================================================================
