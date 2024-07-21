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
  integer                        :: istring,ibf,jbf,ishell,jshell,igaus,ngaus,ngaus_nrl
  integer                        :: nshell,nshell_nrel,nbasis,ntyp,shell_typ,shell_typ_nrl
  real(dp)                       :: eext
  real(dp)                       :: Vext_pq(4),S_pq(4),Dz_pq(4),Dy_pq(4),Dx_pq(4)
  complex(dp)                    :: H4c_me
  logical,allocatable            :: is_large(:),is_large_tmp(:)
  real(dp),allocatable           :: scalar_s_matrix(:,:)
  real(dp),allocatable           :: scalar_nucleus(:,:)
  real(dp),allocatable           :: scalar_nabla_ao(:,:,:) ! < alpha | nabla_r | beta >
  real(dp),allocatable           :: s_matrix_4c(:,:)
  complex(dp),allocatable        :: H_4c_ukb_mat(:,:)
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
  basis_name_nrel(:)=basis_name_1
  
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

  !! Calculate all integrals for UKB H_4C
  call init_libcint(basis)
  !! Calculate overlap matrix S
  call clean_allocate('Scalar overlap matrix S',scalar_s_matrix,basis%nbf,basis%nbf)
  call clean_allocate('Scalar nabla operator D',scalar_nabla_ao,basis%nbf,basis%nbf,3)
  call clean_allocate('Scalar nucleus operator V',scalar_nucleus,basis%nbf,basis%nbf)
  !! S only depends onto the basis set
  call setup_overlap(basis,scalar_s_matrix)
  !! Nucleus-electron interaction
  call setup_nucleus(basis,scalar_nucleus)
  !! External electric field
  call setup_electric_field(basis,scalar_nucleus,eext)
  !! <alpha | nabla_r | beta >
  call setup_nabla_ao(basis,scalar_nabla_ao)
  write(stdout,'(a)') ' '
  !! destroy libcint info
  call destroy_libcint(basis)

  !! Initialize the H_4C UKB
  write(stdout,'(/,a)') 'Building the 4C Hamiltonian'
  nbasis=2*basis%nbf
  call clean_allocate('H_4C in UKB',H_4c_ukb_mat,nbasis,nbasis)
  call clean_allocate('4C overlap matrix S',s_matrix_4c,nbasis,nbasis)
  allocate(is_large_tmp(nbasis))
  H_4c_ukb_mat(:,:)=complex_zero
  s_matrix_4c(:,:)=zero
  do ibf=1,basis%nbf

   is_large_tmp(2*ibf-1)=is_large(ibf)
   is_large_tmp(2*ibf)=is_large(ibf)

   do jbf=1,basis%nbf
    ! Set H_4c_ukb_mat and s_matrix_4c
    if(is_large(ibf)) then  ! L
     if(is_large(jbf)) then  ! LL
      ! s_matrix_4c
      s_matrix_4c(2*ibf-1,2*jbf-1)=scalar_s_matrix(ibf,jbf)
      s_matrix_4c(2*ibf,2*jbf)=scalar_s_matrix(ibf,jbf)
      ! p1 and q1
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Vext_pq(1)=scalar_nucleus(ibf,jbf)
      S_pq(1)=scalar_s_matrix(ibf,jbf)
      H_4c_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p1 and q2
      H_4c_ukb_mat(2*ibf-1,2*jbf)=complex_zero
      ! p2 and q1
      H_4c_ukb_mat(2*ibf,2*jbf-1)=complex_zero
      ! p2 and q2
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Vext_pq(2)=scalar_nucleus(ibf,jbf)
      S_pq(2)=scalar_s_matrix(ibf,jbf)
      H_4c_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     else                    ! LS
      ! p1 and q3
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(1)=scalar_nabla_ao(ibf,jbf,3)
      H_4c_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p1 and q4
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(1)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(1)=scalar_nabla_ao(ibf,jbf,2)
      H_4c_ukb_mat(2*ibf-1,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p2 and q3
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(2)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(2)=scalar_nabla_ao(ibf,jbf,2)
      H_4c_ukb_mat(2*ibf,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p2 and q4
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(2)=scalar_nabla_ao(ibf,jbf,3)
      H_4c_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     endif
    else                    ! S 
     if(is_large(jbf)) then  ! SL
      ! p3 and q1
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(3)=scalar_nabla_ao(ibf,jbf,3)
      H_4c_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p3 and q2
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(3)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(3)=scalar_nabla_ao(ibf,jbf,2)
      H_4c_ukb_mat(2*ibf-1,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p4 and q1
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(4)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(4)=scalar_nabla_ao(ibf,jbf,2)
      H_4c_ukb_mat(2*ibf,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p4 and q2
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(4)=scalar_nabla_ao(ibf,jbf,3)
      H_4c_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     else                    ! SS
      ! s_matrix_4c
      s_matrix_4c(2*ibf-1,2*jbf-1)=scalar_s_matrix(ibf,jbf)
      s_matrix_4c(2*ibf,2*jbf)=scalar_s_matrix(ibf,jbf)
      ! p3 and q3
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Vext_pq(3)=scalar_nucleus(ibf,jbf)
      S_pq(3)=scalar_s_matrix(ibf,jbf)
      H_4c_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p3 and q4
      H_4c_ukb_mat(2*ibf-1,2*jbf)=complex_zero
      ! p4 and q3
      H_4c_ukb_mat(2*ibf,2*jbf-1)=complex_zero
      ! p4 and q4
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Vext_pq(4)=scalar_nucleus(ibf,jbf)
      S_pq(4)=scalar_s_matrix(ibf,jbf)
      H_4c_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     endif
    endif  
 
   enddo

  enddo
  write(stdout,'(a,/)') 'Completed the 4C Hamiltonian'

  !! No longer need the scalar matrices (integrals) TODO: I will probably need < alpha | nabla_r | beta > in restoring the RKB...
  deallocate(is_large)
  call clean_deallocate('Scalar overlap matrix S',scalar_s_matrix)
  call clean_deallocate('Scalar nabla operator D',scalar_nabla_ao)
  call clean_deallocate('Scalar nucleus operator V',scalar_nucleus)

  !! Shuffle the matrices to the
  !! ( LL  LS )
  !! ( SL  SS )
  !! shape
  call shuffle_complex(nbasis,is_large_tmp,H_4c_ukb_mat)
  call shuffle_real(nbasis,is_large_tmp,s_matrix_4c)
  deallocate(is_large_tmp)

  !! TODO
  !! - Lowdin orthogonalize
  !! - Set RKB (including orthogonalize)
  !! - Diag. to find the eigenvectors (Coef)
  !! - Use Coef to find R decoupling matrix
  !! - Transform to the X2C matrices

  call clean_deallocate('H_4C in UKB',H_4c_ukb_mat)
  call clean_deallocate('4C overlap matrix S',s_matrix_4c)

  write(stdout,'(/,a)') ' Completed X2C Hamiltonian construction'
  write(stdout,'(a,/)') ' ======================================'

#endif

  deallocate(basis_name_nrel)
  call stop_clock(timing_x2c)

end subroutine x2c_init

!==================================================================
function H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq) result(H_val)
  use m_definitions
  implicit none
  real(dp),intent(in)  :: Vext_pq(4)
  real(dp),intent(in)  :: S_pq(4)
  real(dp),intent(in)  :: Dz_pq(4)
  real(dp),intent(in)  :: Dy_pq(4)
  real(dp),intent(in)  :: Dx_pq(4)
  complex(dp)          :: H_val
  !=====
  !=====

  H_val =                           Vext_pq(1) + Vext_pq(2) + Vext_pq(3) + Vext_pq(4)  & 
    &   +c_speedlight*c_speedlight*(   S_pq(1) +    S_pq(2) -    S_pq(3) -    S_pq(4)) &
    &   -c_speedlight*im          *(  Dz_pq(1) -   Dz_pq(2) +   Dz_pq(3) -   Dz_pq(4)) &
    &   -c_speedlight             *(  Dy_pq(1) -   Dy_pq(2) +   Dy_pq(3) -   Dy_pq(4)) &
    &   -c_speedlight*im          *(  Dx_pq(1) +   Dx_pq(2) +   Dx_pq(3) +   Dx_pq(4))

end function H4c_me

!==================================================================

subroutine shuffle_complex(nbasis,is_large,matrix)
  use m_definitions
  implicit none
  integer,intent(in)        :: nbasis
  logical,intent(in)        :: is_large(nbasis)
  complex(dp),intent(inout) :: matrix(nbasis,nbasis)
  !=====
  integer                   :: ibf,jbf
  complex(dp),allocatable   :: tmp_matrix(:,:)
  !=====

  allocate(tmp_matrix(nbasis,nbasis))
  jbf=1
  do ibf=1,nbasis
   if(is_large(ibf)) then
    tmp_matrix(jbf,:)=matrix(ibf,:)
    jbf=jbf+1  
   endif
  enddo
  do ibf=1,nbasis
   if(.not.is_large(ibf)) then
    tmp_matrix(jbf,:)=matrix(ibf,:)
    jbf=jbf+1  
   endif
  enddo
  jbf=1
  do ibf=1,nbasis
   if(is_large(ibf)) then
    matrix(:,jbf)=tmp_matrix(:,ibf)
    jbf=jbf+1  
   endif
  enddo
  do ibf=1,nbasis
   if(.not.is_large(ibf)) then
    matrix(:,jbf)=tmp_matrix(:,ibf)
    jbf=jbf+1  
   endif
  enddo
  deallocate(tmp_matrix)

end subroutine shuffle_complex

!==================================================================

subroutine shuffle_real(nbasis,is_large,matrix)
  use m_definitions
  implicit none
  integer,intent(in)     :: nbasis
  logical,intent(in)     :: is_large(nbasis)
  real(dp),intent(inout) :: matrix(nbasis,nbasis)
  !=====
  integer                :: ibf,jbf
  real(dp),allocatable   :: tmp_matrix(:,:)
  !=====

  allocate(tmp_matrix(nbasis,nbasis))
  jbf=1
  do ibf=1,nbasis
   if(is_large(ibf)) then
    tmp_matrix(jbf,:)=matrix(ibf,:)
    jbf=jbf+1  
   endif
  enddo
  do ibf=1,nbasis
   if(.not.is_large(ibf)) then
    tmp_matrix(jbf,:)=matrix(ibf,:)
    jbf=jbf+1  
   endif
  enddo
  jbf=1
  do ibf=1,nbasis
   if(is_large(ibf)) then
    matrix(:,jbf)=tmp_matrix(:,ibf)
    jbf=jbf+1  
   endif
  enddo
  do ibf=1,nbasis
   if(.not.is_large(ibf)) then
    matrix(:,jbf)=tmp_matrix(:,ibf)
    jbf=jbf+1  
   endif
  enddo
  deallocate(tmp_matrix)

end subroutine shuffle_real

!==================================================================
