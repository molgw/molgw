!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval, Mauricio Rodriguez-Mayorga
!
! This module contains
! the routines to calculate the ACFD total energy
!
!=========================================================================
module m_acfd
  use m_definitions
  use m_timing
  use m_warning
  use m_memory
  use m_inputparam
  use m_mpi
  use m_scalapack
  use m_cart_to_pure
  use m_basis_set
  use m_spectral_function
  use m_spectra
  use m_scf
  use m_linear_response
  use m_numerical_tools,only: coeffs_gausslegint




contains


!=========================================================================
subroutine acfd_total_energy(basis,nstate,occupation,energy,c_matrix,en_mbpt)
  type(basis_set),intent(in)               :: basis
  integer,intent(in)                       :: nstate
  real(dp),intent(in)                      :: occupation(:,:),energy(:,:)
  real(dp),intent(in)                      :: c_matrix(:,:,:)
  type(energy_contributions),intent(inout) :: en_mbpt
  !=====
  type(spectral_function)    :: wpol
  real(dp)                   :: erpa_singlet,erpa_triplet,egw_tmp
  integer,parameter          :: nlambda=20
  real(dp)                   :: wlambda(nlambda),lambda(nlambda)
  real(dp),allocatable       :: x_matrix(:,:),y_matrix(:,:)
  real(dp),allocatable       :: a_matrix(:,:),b_matrix(:,:)
  integer                    :: nmat,desc_x(NDEL)
  !=====

  select case(TRIM(postscf))
  case('RPA')
    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol,enforce_spin_multiplicity=1)
    call destroy_spectral_function(wpol)
    en_mbpt%rpa = erpa_singlet
    write(stdout,'(a,2x,f19.10)') ' RPA Energy      (Ha):',en_mbpt%rpa

    en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPA Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  case('RPAX','RPAX-II')
    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol,enforce_spin_multiplicity=1)
    call destroy_spectral_function(wpol)
    en_mbpt%rpa = 0.50_dp * erpa_singlet

    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_triplet,egw_tmp,wpol,enforce_spin_multiplicity=3)
    call destroy_spectral_function(wpol)
    en_mbpt%rpa = en_mbpt%rpa + 1.50_dp * erpa_triplet
    write(stdout,'(a,2x,f19.10)') ' RPAx Energy      (Ha):',en_mbpt%rpa

    en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPAx Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  case('RPA-I')
    write(stdout,'(/,1x,a,i4)') 'RPA with integration over lambda with nlambda: ',nlambda 
    call coeffs_gausslegint(0.0_dp,1.0_dp,lambda,wlambda,nlambda)

    call init_spectral_function(nstate,occupation,0,wpol)
    nmat = wpol%npole_reso
    m_x = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
    n_x = NUMROC(nmat,block_col,ipcol_sd,first_col,npcol_sd)
    call DESCINIT(desc_x,nmat,nmat,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,m_x),info)
  
    call clean_allocate('X matrix',x_matrix,m_x,n_x)
    call clean_allocate('Y matrix',y_matrix,m_x,n_x)

    call clean_allocate('A matrix',a_matrix,m_x,n_x)
    call clean_allocate('B matrix',b_matrix,m_x,n_x)
    ! Get A and B
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                        enforce_spin_multiplicity=1,lambda=1.0_dp,a_matrix=a_matrix,b_matrix=b_matrix)
    call destroy_spectral_function(wpol)

    en_mbpt%rpa = 0.0_dp
    do ilambda=1,nlambda
      write(stdout,'(1x,a,i4,a,i4)') '=== Lambda',ilambda,' / ',nlambda
      call init_spectral_function(nstate,occupation,0,wpol)
      call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                          enforce_spin_multiplicity=1,lambda=lambda(ilambda),x_matrix=x_matrix,y_matrix=y_matrix)
      call destroy_spectral_function(wpol)

      call calculate_ec_acft(desc_x,a_matrix,b_matrix,x_matrix,y_matrix,erpa_singlet)

      en_mbpt%rpa = en_mbpt%rpa + erpa_singlet * wlambda(ilambda)
    enddo
    if(has_auxil_basis) call destroy_eri_3center_eigen()

    call clean_deallocate('X matrix',x_matrix)
    call clean_deallocate('Y matrix',y_matrix)
    call clean_deallocate('A matrix',a_matrix)
    call clean_deallocate('B matrix',b_matrix)

    write(stdout,'(a,2x,f19.10)') ' RPA Energy      (Ha):',en_mbpt%rpa

    en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPA Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  case('RPAX-I')
    write(stdout,'(/,1x,a,i4)') 'RPAx with integration over lambda with nlambda: ',nlambda 
    call coeffs_gausslegint(0.0_dp,1.0_dp,lambda,wlambda,nlambda)

    call init_spectral_function(nstate,occupation,0,wpol)
    nmat = wpol%npole_reso
    m_x = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
    n_x = NUMROC(nmat,block_col,ipcol_sd,first_col,npcol_sd)
    call DESCINIT(desc_x,nmat,nmat,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,m_x),info)
  
    call clean_allocate('X matrix',x_matrix,m_x,n_x)
    call clean_allocate('Y matrix',y_matrix,m_x,n_x)

    call clean_allocate('A matrix',a_matrix,m_x,n_x)
    call clean_allocate('B matrix',b_matrix,m_x,n_x)
    !
    ! Singlet
    !
    ! Get A and B
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                        enforce_spin_multiplicity=1,lambda=1.0_dp,a_matrix=a_matrix,b_matrix=b_matrix)
    call destroy_spectral_function(wpol)

    en_mbpt%rpa = 0.0_dp
    do ilambda=1,nlambda
      write(stdout,'(1x,a,i4,a,i4)') '=== Lambda',ilambda,' / ',nlambda
      call init_spectral_function(nstate,occupation,0,wpol)
      call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_singlet,egw_tmp,wpol, &
                          enforce_spin_multiplicity=1,lambda=lambda(ilambda),x_matrix=x_matrix,y_matrix=y_matrix)
      call destroy_spectral_function(wpol)

      call calculate_ec_acft(desc_x,a_matrix,b_matrix,x_matrix,y_matrix,erpa_singlet)

      en_mbpt%rpa = en_mbpt%rpa + erpa_singlet * wlambda(ilambda)
    enddo
    call destroy_spectral_function(wpol)
    !
    ! Triplet
    !
    ! Get A and B
    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_triplet,egw_tmp,wpol, &
                        enforce_spin_multiplicity=3,lambda=1.0_dp,a_matrix=a_matrix,b_matrix=b_matrix)
    call destroy_spectral_function(wpol)

    do ilambda=1,nlambda
      write(stdout,'(1x,a,i4,a,i4)') '=== Lambda',ilambda,' / ',nlambda
      call init_spectral_function(nstate,occupation,0,wpol)
      call polarizability(.FALSE.,.FALSE.,basis,occupation,energy,c_matrix,erpa_triplet,egw_tmp,wpol, &
                          enforce_spin_multiplicity=3,lambda=lambda(ilambda),x_matrix=x_matrix,y_matrix=y_matrix)
      call destroy_spectral_function(wpol)

      call calculate_ec_acft(desc_x,a_matrix,b_matrix,x_matrix,y_matrix,erpa_triplet)

      en_mbpt%rpa = en_mbpt%rpa + 3.0_dp * erpa_triplet * wlambda(ilambda)
    enddo


    if(has_auxil_basis) call destroy_eri_3center_eigen()

    call clean_deallocate('X matrix',x_matrix)
    call clean_deallocate('Y matrix',y_matrix)
    call clean_deallocate('A matrix',a_matrix)
    call clean_deallocate('B matrix',b_matrix)

    write(stdout,'(a,2x,f19.10)') ' RPA Energy      (Ha):',en_mbpt%rpa

    en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    write(stdout,*)
    write(stdout,'(a,2x,f19.10)') ' RPA Total Energy (Ha):',en_mbpt%total
    write(stdout,*)

  case default
    call die('acfd_total_energy: postscf option not recognized')
  end select
  call print_energy_yaml('mbpt energy',en_mbpt)

end subroutine acfd_total_energy


!=========================================================================
subroutine calculate_ec_acft(desc_x,a_matrix,b_matrix,x_matrix,y_matrix,erpa)
  implicit none

  integer,intent(in)   :: desc_x(NDEL)
  real(dp),intent(in)  :: x_matrix(:,:),y_matrix(:,:),a_matrix(:,:),b_matrix(:,:)
  real(dp),intent(out) :: erpa
  !=====
  integer :: imat,nmat
  real(dp),allocatable :: m_matrix(:,:)
  !=====

  allocate(m_matrix,MOLD=a_matrix)
  nmat=SIZE(a_matrix,DIM=1)
  m_matrix(:,:) = -a_matrix(:,:)
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( x_matrix , MATMUL(a_matrix,TRANSPOSE(x_matrix)) )
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( y_matrix , MATMUL(a_matrix,TRANSPOSE(y_matrix)) )
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( x_matrix , MATMUL(b_matrix,TRANSPOSE(y_matrix)) )
  m_matrix(:,:) = m_matrix(:,:) + MATMUL( y_matrix , MATMUL(b_matrix,TRANSPOSE(x_matrix)) )

  erpa = 0.0_dp
  do imat=1,nmat
    erpa = erpa + m_matrix(imat,imat)
  enddo



end subroutine calculate_ec_acft


end module m_acfd
!=========================================================================
