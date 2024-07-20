!==================================================================
! This file is part of MOLGW.
! Author: M. Rodriguez-Mayorga
!
! This file contains
! - The construction of the X2C Hamiltonian 
!=========================================================================
#include "molgw.h"
subroutine x2c_init(basis)
  use m_definitions
  use m_warning
  use m_timing
  use m_mpi
  use m_elements
  use m_string_tools, only: orbital_momentum_name,append_to_list
  use m_atoms
  use m_ecp
  use m_gaussian
  use m_cart_to_pure
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_hamiltonian_onebody
  use m_libcint_tools
  implicit none

  type(basis_set),intent(inout)  :: basis
  !====
  type(basis_set)                :: basis_nrel
  character(len=100)             :: basis_name_1
  character(len=100),allocatable :: basis_name_nrel(:)
  logical                        :: this_is_large
  logical,allocatable            :: is_large(:),is_large_tmp(:)
  integer                        :: istring,ibf,jbf,ishell,jshell,igaus,ngaus,ngaus_nrl
  integer                        :: nshell,nshell_nrel,nbasis,ntyp,shell_typ,shell_typ_nrl
  real(dp)                       :: eext
  real(dp),allocatable           :: scalar_s_matrix(:,:)
  real(dp),allocatable           :: scalar_nucleus(:,:)
  real(dp),allocatable           :: scalar_nabla_ao(:,:,:) ! < alpha | nabla_r | beta >
  !=====

  call start_clock(timing_x2c)

#if defined(HAVE_LIBCINT)

  allocate(basis_name_nrel(ncenter_basis))

  write(stdout,'(/,a)') ' X2C Hamiltonian construction'
  write(stdout,'(a)')   ' ============================'

  write(stdout,'(/,a,/)') ' Computing large component reference basis'
  !! Find the reference non-rel. basis
  write(basis_name_1,'(a)') trim(basis_name(1))
  do istring=1,len(basis_name_1)
    if( (basis_name_1(istring:istring)=='_' .and. basis_name_1(istring+1:istring+1)=='r') .and.   &
    &  (  basis_name_1(istring+2:istring+2)=='e' .and. basis_name_1(istring+3:istring+3)=='l') ) then
      basis_name_1=basis_name_1(1:istring-1)
    endif
  enddo
  basis_name_nrel(1:ncenter_basis)=basis_name_1
  
  !! Initialize the non-rel. basis (used to find the large component AOs)
  call init_basis_set(basis_path,basis_name_nrel,ecp_basis_name,gaussian_type, &
                        even_tempered_alpha,even_tempered_beta,even_tempered_n_list,basis_nrel)
  nshell=SIZE(basis%shell(:)%icenter)
  nshell_nrel=SIZE(basis_nrel%shell(:)%icenter)
  allocate(is_large_tmp(nshell))
  is_large_tmp(:)=.false.
  do ishell=1,nshell
   shell_typ=basis%shell(ishell)%am   ! 0 for s, 1 for p, 2 for d, 3 for f,...,
   ngaus=basis%shell(ishell)%ng
   do jshell=1,nshell_nrel
    shell_typ_nrl=basis_nrel%shell(jshell)%am   ! 0 for s, 1 for p, 2 for d, 3 for f,...,
    ngaus_nrl=basis_nrel%shell(jshell)%ng
    if((shell_typ==shell_typ_nrl .and. ngaus==ngaus_nrl).and.(.not.is_large_tmp(ishell))) then
     this_is_large=.true.
     do igaus=1,ngaus
       if( basis%shell(ishell)%coeff(igaus) /= basis_nrel%shell(jshell)%coeff(igaus) .or. &
        &  basis%shell(ishell)%alpha(igaus) /= basis_nrel%shell(jshell)%alpha(igaus) ) then
        this_is_large=.false.
       endif
     enddo
     is_large_tmp(ishell)=this_is_large
    endif 
   enddo
  enddo
  
  !! We already found the large component shells. No longer need the large comp. basis
  call destroy_basis_set(basis_nrel)
 
  !! Define atomic basis as large or small
  allocate(is_large(basis%nbf))
  nbasis=0
  is_large=.false.
  do ishell=1,nshell
   shell_typ=basis%shell(ishell)%am
   if(shell_typ==0) then
    ntyp=1
   else if(shell_typ==1) then
    ntyp=3
   else if(shell_typ==2) then
    ntyp=6
   else if(shell_typ==3) then
    ntyp=10
   else if(shell_typ==4) then
    ntyp=15
   else if(shell_typ==5) then
    ntyp=21
   else if(shell_typ==6) then
    ntyp=28
   else 
    call issue_warning('Shell type >6 in X2C is not implemented!') 
   endif
   do ibf=1,ntyp
     is_large(ibf+nbasis)=is_large_tmp(ishell)
   enddo
   nbasis=nbasis+ntyp
  enddo
  deallocate(is_large_tmp)

  !! Calculate all integrals for H^X2C
  call init_libcint(basis)
  !! Calculate overlap matrix S
  call clean_allocate('Scalar overlap matrix S',scalar_s_matrix,basis%nbf,basis%nbf)
  call clean_allocate('Scalar nabla_ao operator M',scalar_nabla_ao,basis%nbf,basis%nbf,3)
  call clean_allocate('Scalar nucleus operator V',scalar_nucleus,basis%nbf,basis%nbf)
  !! S only depends onto the basis set
  call setup_overlap(basis,scalar_s_matrix)
  !! Nucleus-electron interaction
  call setup_nucleus(basis,scalar_nucleus)
  !! External electric field
  call setup_electric_field(basis,scalar_nucleus,eext)
  !! <alpha | nabla_r | beta >
  call setup_nabla_ao(basis,scalar_nabla_ao)
  !! destroy libcint info
  call destroy_libcint(basis)



  call clean_deallocate('Scalar overlap matrix S',scalar_s_matrix)
  call clean_deallocate('Scalar nabla_ao operator M',scalar_nabla_ao)
  call clean_deallocate('Scalar nucleus operator V',scalar_nucleus)
  deallocate(is_large)
  write(stdout,'(/,a)') ' Completed X2C Hamiltonian construction'
  write(stdout,'(a,/)') ' ======================================'

#endif

  deallocate(basis_name_nrel)
  call stop_clock(timing_x2c)

end subroutine x2c_init

!==================================================================
