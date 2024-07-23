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
  use m_hamiltonian_tools
  use m_libcint_tools
  implicit none

  type(basis_set),intent(inout)  :: basis
  !====
  type(basis_set)                :: basis_nrel
  character(len=100)             :: basis_name_1
  character(len=100),allocatable :: basis_name_nrel(:)
  logical                        :: found_basis_name
  logical                        :: this_is_large
  integer                        :: istring,ibf,jbf,iibf,jjbf,ishell,jshell,igaus,ngaus,ngaus_nrl,nstate_large
  integer                        :: nshell,nshell_nrel,nbasis,nbasis_L,nbasis_S,ntyp,shell_typ,shell_typ_nrl
  integer                        :: nstate_rkb
  integer                        :: info,lwork
  real(dp)                       :: eext
  real(dp)                       :: Vext_pq(4),S_pq(4),Dz_pq(4),Dy_pq(4),Dx_pq(4)
  complex(dp)                    :: MpSqL_me,H4c_me ! Functions used to build matrix elements
  logical,allocatable            :: is_large(:),is_large_4c(:)
  integer,allocatable            :: ipiv(:)
  real(dp),allocatable           :: W(:)
  real(dp),allocatable           :: scalar_s_matrix(:,:)
  real(dp),allocatable           :: scalar_nucleus(:,:)
  real(dp),allocatable           :: scalar_nabla_ao(:,:,:) ! < alpha | nabla_r | beta >
  real(dp),allocatable           :: s_matrix_large(:,:)
  real(dp),allocatable           :: x_matrix_large(:,:)
  real(dp),allocatable           :: s_matrix_4c(:,:)
  complex(dp),allocatable        :: Work(:)
  complex(dp),allocatable        :: U_mat(:,:),Tmp_matrix(:,:)
  complex(dp),allocatable        :: c_matrix_ukb2rkb(:,:) ! NOTE: This is like AO^sph to AO^cart in non-rel. calcs.
  complex(dp),allocatable        :: H_4c_ukb_mat(:,:)
  complex(dp),allocatable        :: MpSqL_matrix(:,:)     ! < AO^S | sigma . p | AO^L > matrix (the projector)
  complex(dp),allocatable        :: c_matrix_small(:,:)
  complex(dp),allocatable        :: s_matrix_small(:,:)
  complex(dp),allocatable        :: x_matrix_small(:,:)
  complex(dp),allocatable        :: c_matrix(:,:)
  complex(dp),allocatable        :: s_matrix(:,:)
  complex(dp),allocatable        :: x_matrix(:,:)
  complex(dp),allocatable        :: H_4c_rkb_mat(:,:)
  !=====

  call start_clock(timing_x2c)

#if defined(HAVE_LIBCINT)

  allocate(basis_name_nrel(ncenter_basis))

  write(stdout,'(/,a)') ' X2C Hamiltonian construction'
  write(stdout,'(a)')   ' ============================'

  write(stdout,'(/,a,/)') ' Computing large component reference basis'
  !! Find the reference non-rel. basis
  found_basis_name=.false.
  write(basis_name_1,'(a)') trim(basis_name(1))
  do istring=1,len(basis_name_1)
   if(.not.found_basis_name) then
     if( (basis_name_1(istring:istring)=='_' .and. basis_name_1(istring+1:istring+1)=='r') .and.   &
     &  (  basis_name_1(istring+2:istring+2)=='e' .and. basis_name_1(istring+3:istring+3)=='l') ) then
       basis_name_1=basis_name_1(1:istring-1)
       found_basis_name=.true.
     endif
   endif
  enddo
  basis_name_nrel(:)=basis_name_1
  
  !! Initialize the non-rel. basis (used to find the large component AOs)
  call init_basis_set(basis_path,basis_name_nrel,ecp_basis_name,gaussian_type, &
                        even_tempered_alpha,even_tempered_beta,even_tempered_n_list,basis_nrel)
  nshell=SIZE(basis%shell(:)%icenter)
  nshell_nrel=SIZE(basis_nrel%shell(:)%icenter)
  allocate(is_large_4c(nshell))
  is_large_4c(:)=.false.
  do ishell=1,nshell
   shell_typ=basis%shell(ishell)%am   ! 0 for s, 1 for p, 2 for d, 3 for f,...,
   ngaus=basis%shell(ishell)%ng
   do jshell=1,nshell_nrel
    shell_typ_nrl=basis_nrel%shell(jshell)%am   ! 0 for s, 1 for p, 2 for d, 3 for f,...,
    ngaus_nrl=basis_nrel%shell(jshell)%ng
    if((shell_typ==shell_typ_nrl .and. ngaus==ngaus_nrl).and.(.not.is_large_4c(ishell))) then
     this_is_large=.true.
     do igaus=1,ngaus
       if( basis%shell(ishell)%coeff(igaus) /= basis_nrel%shell(jshell)%coeff(igaus) .or. &
        &  basis%shell(ishell)%alpha(igaus) /= basis_nrel%shell(jshell)%alpha(igaus) ) then
        this_is_large=.false.
       endif
     enddo
     is_large_4c(ishell)=this_is_large
    endif 
   enddo
  enddo

  !! Find if nstate_large=basis_nrel%nbf before proceeding
  call clean_allocate('Large overlap matrix S',s_matrix_large,basis_nrel%nbf,basis_nrel%nbf)
  call init_libcint(basis_nrel)
  call setup_overlap(basis_nrel,s_matrix_large)
  call setup_x_matrix(min_overlap,s_matrix_large,nstate_large,x_matrix_large)
  call clean_deallocate('Large overlap matrix S',s_matrix_large)
  call destroy_libcint(basis_nrel)
  
  !! We already found the large component shells. No longer need the large comp. basis
  call destroy_basis_set(basis_nrel)

  !! Do not proceed if nstate/=basis_nrel%nbf because we need to build for each AO^L an AO^S in restricted-KB (see below).
  if(nstate_large/=basis_nrel%nbf) then
    write(stdout,'(/,a,i10,i10,/)') 'Nstate and basis_large ',nstate_large,basis_nrel%nbf
    call die("X2C requires basis sets that are linearly independent.")
  endif
 
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
     is_large(ibf+nbasis)=is_large_4c(ishell)
   enddo
   nbasis=nbasis+ntyp
  enddo
  deallocate(is_large_4c)

  !! Calculate all integrals for unrestricted-KB H_4C
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

  !! Initialize the H_4C unrestricted-KB
  write(stdout,'(/,a)') ' Building the H_4C in unrestricted-KB (UKB)'
  nbasis=2*basis%nbf
  call clean_allocate('H_4C in UKB',H_4c_ukb_mat,nbasis,nbasis)
  call clean_allocate('4C UKB overlap matrix S',s_matrix_4c,nbasis,nbasis)
  allocate(is_large_4c(nbasis))
  H_4c_ukb_mat(:,:)=complex_zero
  s_matrix_4c(:,:)=zero
  nbasis_L=0; nbasis_S=0;
  do ibf=1,basis%nbf
 
   if(is_large(ibf)) then
    nbasis_L=nbasis_L+1
   else
    nbasis_S=nbasis_S+1
   endif

   is_large_4c(2*ibf-1)=is_large(ibf)
   is_large_4c(2*ibf)=is_large(ibf)

   do jbf=1,basis%nbf
    ! Set H_4c_ukb_mat and s_matrix_4c
    if(is_large(ibf)) then  ! L
     if(is_large(jbf)) then  ! LL
      ! s_matrix_4c
      s_matrix_4c(2*ibf-1,2*jbf-1)=scalar_s_matrix(ibf,jbf)
      s_matrix_4c(2*ibf,2*jbf)=scalar_s_matrix(ibf,jbf)
      ! x_matrix_4c (S^-1/2 large)
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
  nbasis_L=2*nbasis_L; nbasis_S=2*nbasis_S; nstate_rkb=2*nbasis_L;
  write(stdout,'(a,i10)') ' UKB Nbasis Large',nbasis_L
  write(stdout,'(a,i10)') ' UKB Nbasis Small',nbasis_S
  write(stdout,'(a)') ' Doing Lowdin orthonormalization of the Large component'
  write(stdout,'(a)') ' Filling (S^-1/2)^Large in restricted-KB (RKB) X'
  call clean_allocate('Full RKB X matrix',x_matrix,nstate_rkb,nstate_rkb)
  call clean_allocate('Full RKB S matrix',s_matrix,nstate_rkb,nstate_rkb)
  x_matrix=complex_zero; s_matrix=complex_zero;
  do ibf=1,nbasis_L/2
   do jbf=1,nbasis_L/2
    x_matrix(2*ibf-1,2*jbf-1)=x_matrix_large(ibf,jbf)
    x_matrix(2*ibf,2*jbf)=x_matrix_large(ibf,jbf)
    s_matrix(2*ibf-1,2*jbf-1)=scalar_s_matrix(ibf,jbf)
    s_matrix(2*ibf,2*jbf)=scalar_s_matrix(ibf,jbf)
   enddo
  enddo
  write(stdout,'(a,/)') ' Completed the H_4C in UKB'

  !! No longer need these scalar and pure large component matrices (integrals)
  call clean_deallocate('Large X * X^H = S^-1',x_matrix_large)
  call clean_deallocate('Scalar overlap matrix S',scalar_s_matrix)
  call clean_deallocate('Scalar nucleus operator V',scalar_nucleus)

  !! Shuffle the matrices to the
  !! ( LL  LS )
  !! ( SL  SS )
  !! shape
  call shuffle_complex(nbasis,is_large_4c,H_4c_ukb_mat)
  call shuffle_real(nbasis,is_large_4c,s_matrix_4c)
  call shuffle_real(basis%nbf,is_large,scalar_nabla_ao(:,:,1))
  call shuffle_real(basis%nbf,is_large,scalar_nabla_ao(:,:,2))
  call shuffle_real(basis%nbf,is_large,scalar_nabla_ao(:,:,3))
  deallocate(is_large_4c,is_large)

  !! Set RKB for each AO^L -> build an AO^S
  !! Use the M_pS,qL matrix < pS | sigma . p | qL > = \sum_tS S_pS,tS C_tS,qL then C=S^-1 M_pSqL
  write(stdout,'(/,a)') ' Imposing RKB'
  call clean_allocate('M_pSqL matrix ',MpSqL_matrix,nbasis_S,nbasis_L)
  call clean_allocate('Small overlap matrix ',s_matrix_small,nbasis_S,nbasis_S)
  call clean_allocate('Coefficients matrix small ',c_matrix_small,nbasis_S,nbasis_L)
  s_matrix_small(1:nbasis_S,1:nbasis_S)=s_matrix_4c(nbasis_L+1:nbasis,nbasis_L+1:nbasis)
  lwork=-1
  allocate(ipiv(nbasis_S),Work(1))
  call zgetrf(nbasis_S,nbasis_S,s_matrix_small,nbasis_S,ipiv,info)
  if(info/=0) call die("Error in XC2 computing S_small^-1 in zgetrf")
  call zgetri(nbasis_S,s_matrix_small,nbasis_S,ipiv,Work,lwork,info)
  lwork=nint(real(Work(1)))
  deallocate(Work)
  allocate(Work(lwork))
  call zgetri(nbasis_S,s_matrix_small,nbasis_S,ipiv,Work,lwork,info)
  if(info/=0) call die("Error in XC2 computing S_small^-1 in zgetri")
  deallocate(ipiv,Work)
  do ibf=1,nbasis_S
   do jbf=1,nbasis_L
    Dx_pq=zero;Dy_pq=zero;Dz_pq=zero;
    iibf=ibf+nbasis_L
    jjbf=jbf
    ! S1
    if(mod(iibf,2)/=0) then
     iibf=(iibf+1)/2
     ! L1 
     if(mod(jjbf,2)/=0) then
      jjbf=(jjbf+1)/2
      Dz_pq(1)=scalar_nabla_ao(iibf,jjbf,3)
     ! L2 
     else
      jjbf=jjbf/2
      Dx_pq(1)=scalar_nabla_ao(iibf,jjbf,1)
      Dy_pq(1)=scalar_nabla_ao(iibf,jjbf,2)
     endif
    ! S2 
    else
     iibf=iibf/2
     ! L1 
     if(mod(jjbf,2)/=0) then
      jjbf=(jjbf+1)/2
      Dx_pq(2)=scalar_nabla_ao(iibf,jjbf,1)
      Dy_pq(2)=scalar_nabla_ao(iibf,jjbf,2)
     ! L2 
     else
      jjbf=jjbf/2
      Dz_pq(2)=scalar_nabla_ao(iibf,jjbf,3)
     endif
    endif
    MpSqL_matrix(ibf,jbf)=MpSqL_me(Dx_pq,Dy_pq,Dz_pq)
   enddo
  enddo
  c_matrix_small=matmul(s_matrix_small,MpSqL_matrix) !! C = S^-1 M_pSqL
  call clean_deallocate('M_pSqL matrix ',MpSqL_matrix)
  call clean_deallocate('Scalar nabla operator D',scalar_nabla_ao)
  write(stdout,'(a,/)') ' Completed imposing RKB'

  !! Build RKB Hamiltonian
  write(stdout,'(/,a)') ' Building the H_4C in RKB'
  write(stdout,'(a,i10)') ' RKB Nbasis ',nstate_rkb
  call clean_allocate('H_4C in RKB',H_4c_rkb_mat,nstate_rkb,nstate_rkb)
  call clean_allocate('UKB to RKB coefficients',c_matrix_ukb2rkb,nbasis_L+nbasis_S,nstate_rkb)
   ! c_matrix_ukb2rk
   ! ( 1  0              )
   ! ( 0  c_matrix_small )
  c_matrix_ukb2rkb=complex_zero
  do ibf=1,nbasis_L
   c_matrix_ukb2rkb(ibf,ibf)=1.0e0
  enddo
  c_matrix_ukb2rkb(nbasis_L+1:,nbasis_L+1:)=c_matrix_small(:,:)
  H_4c_rkb_mat=matmul(conjg(transpose(c_matrix_ukb2rkb)),matmul(H_4c_ukb_mat,c_matrix_ukb2rkb))
  call clean_deallocate('H_4C in UKB',H_4c_ukb_mat)
  write(stdout,'(a,/)') ' Built the H_4C in RKB'

  !! - Lowdin orthogonalize restricted-KB Small component basis
  write(stdout,'(/,a)') ' Orthonomalizing the RKB small component'
   ! C^dagger s_matrix_small C = S_SS with S_SS being used to define x_matrix_small=(S^SS)^-1/2
  call clean_allocate('(RKB S_SS)^-1/2 matrix',x_matrix_small,nbasis_L,nbasis_L) 
  s_matrix_small(1:nbasis_S,1:nbasis_S)=s_matrix_4c(nbasis_L+1:nbasis,nbasis_L+1:nbasis)
  x_matrix_small=matmul(transpose(conjg(c_matrix_small)),matmul(s_matrix_small,c_matrix_small))
  s_matrix(nbasis_L+1:,nbasis_L+1:)=x_matrix_small(:,:) ! NOTE: save in s_matrix
  allocate(W(nbasis_L),U_mat(nbasis_L,nbasis_L),Tmp_matrix(nbasis_L,nbasis_L)) 
  call diagonalize(' ',x_matrix_small,W,U_mat)
  Tmp_matrix=complex_zero
  do ibf=1,nbasis_L
   Tmp_matrix(ibf,ibf)=1.0e0/(sqrt(W(ibf))+1.0e-10) 
  enddo
  x_matrix_small=matmul(matmul(U_mat,Tmp_matrix),transpose(conjg(U_mat)))
   ! C = C (RKB S_SS)^-1/2 
  c_matrix_small=matmul(c_matrix_small,x_matrix_small)
  x_matrix(nbasis_L+1:,nbasis_L+1:)=x_matrix_small(:,:) ! NOTE: save in x_matrix
   ! Check that it is orthonormal
  x_matrix_small=matmul(transpose(conjg(c_matrix_small)),matmul(s_matrix_small,c_matrix_small))
  do ibf=1,nbasis_L
   do jbf=1,nbasis_L
     if(ibf/=jbf) then
      if(abs(x_matrix_small(ibf,jbf))>1e-8) then
       write(stdout,'(a,i10,i10,2f20.5)') 'Error S comp. orthonorm. ',ibf,jbf,x_matrix_small(ibf,jbf)
      endif
     else
      if(abs(x_matrix_small(ibf,jbf)-1.0e0)>1e-8) then
       write(stdout,'(a,i10,i10,2f20.5)') 'Error S comp. orthonorm. ',ibf,jbf,x_matrix_small(ibf,jbf)
      endif
     endif
   enddo
  enddo
  deallocate(W,U_mat,Tmp_matrix) 
  call clean_deallocate('Small overlap matrix ',s_matrix_small)
  call clean_deallocate('(RBK S_SS)^-1/2 matrix ',x_matrix_small)
  call clean_deallocate('Coefficients matrix small ',c_matrix_small)
  write(stdout,'(a,/)') ' Completed the orthonomalization of the RKB small component'

  !! Diagonalize the H_4C in RKB
   ! H^RKB C = S C e -> H_4c_rkb_mat c_matrix = s_matrix c_matrix e
   ! Building S^-1/2 H^RKB S^-1/2  with S^-1/2 = x_matrix
  write(stdout,'(/,a)') ' Diagonalizing the H_4C in RKB'
  call clean_allocate('Full RKB wavefunctions C',c_matrix,nstate_rkb,nstate_rkb)
  H_4c_rkb_mat=matmul(conjg(transpose(x_matrix)),matmul(H_4c_rkb_mat,x_matrix))
  allocate(W(nstate_rkb),U_mat(nstate_rkb,nstate_rkb)) 
  call diagonalize(' ',H_4c_rkb_mat,W,U_mat)
  W(:)=W(:)-c_speedlight*c_speedlight ! NOTE: People add -c^2 I_4 to each < 4c_AO_basis_p | H | 4c_AO_basis_q > term
  write(stdout,'(a,i10)') ' Note: -c^2 was added the energy of each state'
  write(stdout,'(1x,a)') '=== Energies ==='
  write(stdout,'(a)') '   #                               (Ha)                       &
  &                         (eV)      '
  write(stdout,'(a)') '                         bar                       ubar       &
  &               bar                       ubar '
  do ibf=1,nstate_rkb/2
   write(stdout,'(1x,i5,2(2(1x,f25.5)))') ibf,W(2*ibf-1),W(2*ibf),W(2*ibf-1)*Ha_eV,W(2*ibf)*Ha_eV 
   if(ibf==nbasis_L/2) write(stdout,'(a)') '  --------------------------------------------------------------'
  enddo
  c_matrix=matmul(x_matrix,U_mat) ! NOTE: As we do in m_scf_loop.f90 we multiply S^-1/2 U = c_matrix
  deallocate(W,U_mat)
  write(stdout,'(a,/)') ' Diagonalized the H_4C in RKB'

  !! TODO
  !! - Use c_matrix to find the R decoupling matrix
  !! - Transform to the X2C matrices ( get  ---> H^X2C to be used in SCF calcs.) 

  call clean_deallocate('Full RKB wavefunctions C',c_matrix)
  call clean_deallocate('Full RKB X matrix',x_matrix)
  call clean_deallocate('Full RKB S matrix',s_matrix)
  call clean_deallocate('H_4C in RKB',H_4c_rkb_mat)
  call clean_deallocate('4C UKB overlap matrix S',s_matrix_4c) !! Can we remove it before? do we need it for X2C decoupling?
  call clean_deallocate('UKB to RKB coefficients',c_matrix_ukb2rkb)

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
function MpSqL_me(Dx_pq,Dy_pq,Dz_pq) result(M_pSqL)
  use m_definitions
  implicit none
  real(dp),intent(in)  :: Dx_pq(4)
  real(dp),intent(in)  :: Dy_pq(4)
  real(dp),intent(in)  :: Dz_pq(4)
  complex(dp)          :: M_pSqL
  !=====
  !=====

  M_pSqL = -im * ( Dx_pq(1) + Dx_pq(2) ) & 
    &      -     ( Dy_pq(1) - Dy_pq(2) ) &
    &      -im * ( Dz_pq(1) - Dz_pq(2) ) 

end function MpSqL_me

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
