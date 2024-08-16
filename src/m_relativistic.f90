!==================================================================
! This file is part of MOLGW.
! Author: M. Rodriguez-Mayorga
!
! This file contains
! - The construction of the Relativistic Hamiltonians (4C and X2C) 
!=========================================================================
#include "molgw.h"
module m_relativistic
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

 public::relativistic_init,check_CdaggerSC_I

 private::shuffle_real,shuffle_complex,H4c_me,MpSqL_me

contains

!==================================================================
function H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq) result(H_val)
  real(dp),intent(in)  :: Vext_pq(4)
  real(dp),intent(in)  :: S_pq(4)
  real(dp),intent(in)  :: Dz_pq(4)
  real(dp),intent(in)  :: Dy_pq(4)
  real(dp),intent(in)  :: Dx_pq(4)
  complex(dp)          :: H_val
  !=====
  !=====

  H_val =                      Vext_pq(1)  + Vext_pq(2) + Vext_pq(3) + Vext_pq(4)  & 
    &   +2.0d0*c_speedlight*c_speedlight*(              -    S_pq(3) -    S_pq(4)) & ! -c^2 I_4x4 added
    &   -c_speedlight*im      *(  Dz_pq(1) -   Dz_pq(2) +   Dz_pq(3) -   Dz_pq(4)) &
    &   -c_speedlight         *(  Dy_pq(1) -   Dy_pq(2) +   Dy_pq(3) -   Dy_pq(4)) &
    &   -c_speedlight*im      *(  Dx_pq(1) +   Dx_pq(2) +   Dx_pq(3) +   Dx_pq(4))

  ! Symmetric (+) and (-) energy states for the Dirac Hamitonian (i.e. Vext_pq = 0)
  !  &  +c_speedlight*c_speedlight*(  S_pq(1) +   S_pq(2) -  S_pq(3) -    S_pq(4)) &

end function H4c_me

!==================================================================
function MpSqL_me(Dx_pq,Dy_pq,Dz_pq) result(M_pSqL)
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

!=========================================================================
subroutine relativistic_init(basis,is_x2c,electrons_in,nstate,c_matrix,s_matrix,x_matrix,H_rel_rkb_mat)

  type(basis_set),intent(inout)  :: basis
  logical,intent(in)             :: is_x2c
  real(dp),intent(in)            :: electrons_in
  integer,intent(out)            :: nstate
  complex(dp),allocatable,dimension(:,:),intent(inout)::c_matrix
  complex(dp),allocatable,dimension(:,:),intent(inout)::s_matrix
  complex(dp),allocatable,dimension(:,:),intent(inout)::x_matrix
  complex(dp),allocatable,dimension(:,:),intent(inout)::H_rel_rkb_mat(:,:)
  !====
  type(basis_set)                :: basis_nrel
  character(len=100)             :: basis_name_1
  character(len=100),allocatable :: basis_name_nrel(:)
  logical                        :: found_basis_name
  logical                        :: this_is_large
  integer                        :: nbf_large
  integer                        :: istring,ibf,jbf,iibf,jjbf,ishell,jshell,kshell,igaus,ngaus,ngaus_nrl,nstate_large
  integer                        :: nshell,nshell_nrel,nbasis,nbasis_L,nbasis_S,ntyp,shell_typ,shell_typ_nrl
  integer                        :: nstate_rkb,ielectrons
  integer                        :: info,lwork
  real(dp)                       :: eext
  real(dp)                       :: Vext_pq(4),S_pq(4),Dz_pq(4),Dy_pq(4),Dx_pq(4)
  logical,allocatable            :: is_large(:),is_large_4c(:)
  integer,allocatable            :: ipiv(:)
  real(dp),allocatable           :: W(:),E_state(:)
  real(dp),allocatable           :: scalar_s_matrix(:,:)
  real(dp),allocatable           :: scalar_nucleus(:,:)
  real(dp),allocatable           :: scalar_nabla_ao(:,:,:) ! < alpha | nabla_r | beta >
  real(dp),allocatable           :: s_matrix_large(:,:)
  real(dp),allocatable           :: x_matrix_large(:,:)
  real(dp),allocatable           :: s_matrix_4c(:,:)
  complex(dp),allocatable        :: Work(:)
  complex(dp),allocatable        :: A_mat(:,:),B_mat(:,:),R_mat(:,:)
  complex(dp),allocatable        :: U_mat(:,:),V_mat(:,:),Tmp_matrix(:,:)
  complex(dp),allocatable        :: c_matrix_ukb2rkb(:,:) ! NOTE: This is like AO^sph to AO^cart in non-rel. calcs.
  complex(dp),allocatable        :: H_rel_ukb_mat(:,:)
  complex(dp),allocatable        :: MpSqL_matrix(:,:)     ! < AO^S | sigma . p | AO^L > matrix (the projector)
  complex(dp),allocatable        :: c_matrix_small(:,:)
  complex(dp),allocatable        :: s_matrix_small(:,:)
  complex(dp),allocatable        :: x_matrix_small(:,:)
  complex(dp),allocatable        :: H_rel_rkb_ortho_mat(:,:)
  !=====

  call start_clock(timing_relativistic)


#if defined(HAVE_LIBCINT)


  allocate(basis_name_nrel(ncenter_basis))

  write(stdout,'(/,a)') ' Relativitistic Hamiltonian construction'
  write(stdout,'(a,/)') ' ======================================='

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
  nbf_large=basis_nrel%nbf
  deallocate(basis_name_nrel)
  nshell=SIZE(basis%shell(:)%icenter)
  nshell_nrel=SIZE(basis_nrel%shell(:)%icenter)
  allocate(is_large_4c(nshell))
  is_large_4c(:)=.false.
  kshell=0
  do ishell=1,nshell
   shell_typ=basis%shell(ishell)%am   ! 0 for s, 1 for p, 2 for d, 3 for f,...,
   ngaus=basis%shell(ishell)%ng
   if(kshell<nshell_nrel) then
    do jshell=1,nshell_nrel
     shell_typ_nrl=basis_nrel%shell(jshell)%am   ! 0 for s, 1 for p, 2 for d, 3 for f,...,
     ngaus_nrl=basis_nrel%shell(jshell)%ng
     if((shell_typ==shell_typ_nrl .and. ngaus==ngaus_nrl).and.(.not.is_large_4c(ishell))) then
      this_is_large=.true.
      do igaus=1,ngaus
        if( abs(basis%shell(ishell)%coeff(igaus) - basis_nrel%shell(jshell)%coeff(igaus))>1e-4 .or. &
         &  abs(basis%shell(ishell)%alpha(igaus) - basis_nrel%shell(jshell)%alpha(igaus))>1e-4 ) then
         this_is_large=.false.
        endif
      enddo
      is_large_4c(ishell)=this_is_large
      if(is_large_4c(ishell)) kshell=kshell+1
     endif
    enddo
   endif
  enddo

  !! Find if nstate_large=basis_nrel%nbf before proceeding
  call clean_allocate('Large overlap matrix S',s_matrix_large,basis_nrel%nbf,basis_nrel%nbf)
  call init_libcint(basis_nrel)
  call setup_overlap(basis_nrel,s_matrix_large)
  call setup_x_matrix(min_overlap,s_matrix_large,nstate_large,x_matrix_large)
  call destroy_libcint(basis_nrel)
  
  !! We already found the large component shells. No longer need the large comp. basis
  call destroy_basis_set(basis_nrel)

  !! Do not proceed if nstate/=basis_nrel%nbf because we need to build for each AO^L an AO^S in restricted-KB (see below).
  if(nstate_large/=nbf_large) then
    write(stdout,'(/,a,i10,i10,/)') 'Nstate and basis_large ',nstate_large,basis_nrel%nbf
    call die("Relativistic requires basis sets that are linearly independent.")
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
    call die('Shell type >6 in Relativistic is not implemented!') 
   endif
   do ibf=1,ntyp
     is_large(ibf+nbasis)=is_large_4c(ishell)
   enddo
   nbasis=nbasis+ntyp
  enddo
  deallocate(is_large_4c)

  !! Calculate all integrals for unrestricted-KB H_rel
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

  !! Initialize the H_rel unrestricted-KB
  write(stdout,'(/,a)') ' Building the H_rel in unrestricted-KB (UKB)'
  nbasis=2*basis%nbf
  call clean_allocate('H_rel in UKB',H_rel_ukb_mat,nbasis,nbasis)
  call clean_allocate('4C UKB overlap matrix S',s_matrix_4c,nbasis,nbasis)
  allocate(is_large_4c(nbasis))
  H_rel_ukb_mat(:,:)=complex_zero
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
    ! Set H_rel_ukb_mat and s_matrix_4c
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
      H_rel_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p1 and q2
      H_rel_ukb_mat(2*ibf-1,2*jbf)=complex_zero
      ! p2 and q1
      H_rel_ukb_mat(2*ibf,2*jbf-1)=complex_zero
      ! p2 and q2
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Vext_pq(2)=scalar_nucleus(ibf,jbf)
      S_pq(2)=scalar_s_matrix(ibf,jbf)
      H_rel_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     else                    ! LS
      ! p1 and q3
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(1)=scalar_nabla_ao(ibf,jbf,3)
      H_rel_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p1 and q4
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(1)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(1)=scalar_nabla_ao(ibf,jbf,2)
      H_rel_ukb_mat(2*ibf-1,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p2 and q3
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(2)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(2)=scalar_nabla_ao(ibf,jbf,2)
      H_rel_ukb_mat(2*ibf,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p2 and q4
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(2)=scalar_nabla_ao(ibf,jbf,3)
      H_rel_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     endif
    else                    ! S 
     if(is_large(jbf)) then  ! SL
      ! p3 and q1
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(3)=scalar_nabla_ao(ibf,jbf,3)
      H_rel_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p3 and q2
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(3)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(3)=scalar_nabla_ao(ibf,jbf,2)
      H_rel_ukb_mat(2*ibf-1,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p4 and q1
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dx_pq(4)=scalar_nabla_ao(ibf,jbf,1)
      Dy_pq(4)=scalar_nabla_ao(ibf,jbf,2)
      H_rel_ukb_mat(2*ibf,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p4 and q2
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Dz_pq(4)=scalar_nabla_ao(ibf,jbf,3)
      H_rel_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     else                    ! SS
      ! s_matrix_4c
      s_matrix_4c(2*ibf-1,2*jbf-1)=scalar_s_matrix(ibf,jbf)
      s_matrix_4c(2*ibf,2*jbf)=scalar_s_matrix(ibf,jbf)
      ! p3 and q3
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Vext_pq(3)=scalar_nucleus(ibf,jbf)
      S_pq(3)=scalar_s_matrix(ibf,jbf)
      H_rel_ukb_mat(2*ibf-1,2*jbf-1)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
      ! p3 and q4
      H_rel_ukb_mat(2*ibf-1,2*jbf)=complex_zero
      ! p4 and q3
      H_rel_ukb_mat(2*ibf,2*jbf-1)=complex_zero
      ! p4 and q4
      Vext_pq=zero;S_pq=zero;Dz_pq=zero;Dy_pq=zero;Dx_pq=zero;
      Vext_pq(4)=scalar_nucleus(ibf,jbf)
      S_pq(4)=scalar_s_matrix(ibf,jbf)
      H_rel_ukb_mat(2*ibf,2*jbf)=H4c_me(Vext_pq,S_pq,Dz_pq,Dy_pq,Dx_pq)
     endif
    endif  
 
   enddo

  enddo
  nbasis_L=2*nbasis_L; nbasis_S=2*nbasis_S; nstate_rkb=2*nbasis_L;
  write(stdout,'(a,i10)') ' UKB Nbasis Large',nbasis_L
  write(stdout,'(a,i10)') ' UKB Nbasis Small',nbasis_S
  write(stdout,'(a)') ' Doing Lowdin orthonormalization for the Large component'
  write(stdout,'(a)') ' Filling (S^-1/2)^Large in restricted-KB (RKB) X'
  call clean_allocate('Full RKB X matrix',x_matrix,nstate_rkb,nstate_rkb)
  call clean_allocate('Full RKB S matrix',s_matrix,nstate_rkb,nstate_rkb)
  x_matrix=complex_zero; s_matrix=complex_zero;
  do ibf=1,nbasis_L/2
   do jbf=1,nbasis_L/2
    x_matrix(2*ibf-1,2*jbf-1)=x_matrix_large(ibf,jbf)
    x_matrix(2*ibf,2*jbf)=x_matrix_large(ibf,jbf)
    s_matrix(2*ibf-1,2*jbf-1)=s_matrix_large(ibf,jbf)
    s_matrix(2*ibf,2*jbf)=s_matrix_large(ibf,jbf)
   enddo
  enddo
  write(stdout,'(a,/)') ' Completed the H_rel in UKB'
  !! Do not proceed if nbasis_L/2 /= nbf_large because we did not find all the large component basis
  if(nbasis_L/2 /= nbf_large) then
    write(stdout,'(/,a,i10,i10,/)') 'Nbasis_L/2 and basis%nbf large ',nbasis_L/2,nbf_large
    call die("Relativistic requires large comp. basis sets in files foo_rel to be also present in files foo.")
  endif

  !! No longer need these scalar and pure large component matrices (integrals)
  call clean_deallocate('Large overlap matrix S',s_matrix_large)
  call clean_deallocate('Large X * X^H = S^-1',x_matrix_large)
  call clean_deallocate('Scalar overlap matrix S',scalar_s_matrix)
  call clean_deallocate('Scalar nucleus operator V',scalar_nucleus)

  !! Shuffle the matrices to the
  !! ( LL  LS )
  !! ( SL  SS )
  !! shape
  call shuffle_complex(nbasis,is_large_4c,H_rel_ukb_mat)
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
  if(info/=0) then
   call issue_warning("Error in Relativistic computing S_small^-1 in zgetrf")
   s_matrix_small(1:nbasis_S,1:nbasis_S)=s_matrix_4c(nbasis_L+1:nbasis,nbasis_L+1:nbasis)
   allocate(W(nbasis_S),U_mat(nbasis_S,nbasis_S))
   call diagonalize(' ',s_matrix_small,W,U_mat)
   s_matrix_small=COMPLEX_ZERO
   do ibf=1,nbasis_S
    if(abs(W(ibf))<1.0e-8) then
     write(stdout,'(a,i5,f20.8)') ' Eigenvalue lower than 1e-8 in S_mall^-1',ibf,W(ibf)
    endif
    s_matrix_small(ibf,ibf)=1.0d0/(W(ibf)+1.0e-10)
   enddo
   s_matrix_small=matmul(matmul(U_mat,s_matrix_small),transpose(conjg(U_mat)))
   deallocate(W,U_mat)
  else
   call zgetri(nbasis_S,s_matrix_small,nbasis_S,ipiv,Work,lwork,info)
   lwork=nint(real(Work(1)))
   deallocate(Work)
   allocate(Work(lwork))
   call zgetri(nbasis_S,s_matrix_small,nbasis_S,ipiv,Work,lwork,info)
   if(info/=0) then
    call issue_warning("Error in Relativistic computing S_small^-1 in zgetri")
    s_matrix_small(1:nbasis_S,1:nbasis_S)=s_matrix_4c(nbasis_L+1:nbasis,nbasis_L+1:nbasis)
    allocate(W(nbasis_S),U_mat(nbasis_S,nbasis_S))
    call diagonalize(' ',s_matrix_small,W,U_mat)
    s_matrix_small=COMPLEX_ZERO
    do ibf=1,nbasis_S
     if(abs(W(ibf))<1.0e-8) then
      write(stdout,'(a,i5,f20.8)') ' Eigenvalue lower than 1e-8 in S_mall^-1',ibf,W(ibf)
     endif
     s_matrix_small(ibf,ibf)=1.0d0/(W(ibf)+1.0e-10)
    enddo
    s_matrix_small=matmul(matmul(U_mat,s_matrix_small),transpose(conjg(U_mat)))
    deallocate(W,U_mat)
   endif
  endif
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
  write(stdout,'(/,a)') ' Building the H_rel in RKB'
  write(stdout,'(a,i10)') ' RKB Nbasis ',nstate_rkb
  call clean_allocate('H_rel in RKB',H_rel_rkb_mat,nstate_rkb,nstate_rkb)
  call clean_allocate('UKB to RKB coefficients',c_matrix_ukb2rkb,nbasis_L+nbasis_S,nstate_rkb) ! ML+MS x 2ML
   ! c_matrix_ukb2rk
   ! ( 1  0              )
   ! ( 0  c_matrix_small )
  c_matrix_ukb2rkb=complex_zero
  do ibf=1,nbasis_L
   c_matrix_ukb2rkb(ibf,ibf)=1.0d0
  enddo
  c_matrix_ukb2rkb(nbasis_L+1:,nbasis_L+1:)=c_matrix_small(:,:)
  H_rel_rkb_mat=matmul(conjg(transpose(c_matrix_ukb2rkb)),matmul(H_rel_ukb_mat,c_matrix_ukb2rkb))
  call clean_deallocate('H_rel in UKB',H_rel_ukb_mat)
  write(stdout,'(a,/)') ' Completed the H_rel in RKB'

  !! - Lowdin orthogonalize restricted-KB Small component basis
  write(stdout,'(/,a)') ' Orthonomalizing the RKB small component'
   ! C^dagger s_matrix_small C = S_SS with S_SS being used to define x_matrix_small=(S^SS)^-1/2
  call clean_allocate('(RKB S_SS)^-1/2 matrix',x_matrix_small,nbasis_L,nbasis_L) 
  s_matrix_small(1:nbasis_S,1:nbasis_S)=s_matrix_4c(nbasis_L+1:nbasis,nbasis_L+1:nbasis)
  x_matrix_small=matmul(transpose(conjg(c_matrix_small)),matmul(s_matrix_small,c_matrix_small)) ! (C_S)^dagger S_SS C_S ? I
  s_matrix(nbasis_L+1:,nbasis_L+1:)=x_matrix_small(:,:) ! NOTE: save in s_matrix
  allocate(W(nbasis_L),U_mat(nbasis_L,nbasis_L),Tmp_matrix(nbasis_L,nbasis_L)) 
  call diagonalize(' ',x_matrix_small,W,U_mat)
  Tmp_matrix=complex_zero
  do ibf=1,nbasis_L
   if(abs(W(ibf))<1.0e-8) then
    write(stdout,'(a,i5,f20.8)') ' Eigenvalue lower than 1e-8 in (RKB S_SS)^-1/2',ibf,W(ibf)
   endif
   Tmp_matrix(ibf,ibf)=1.0d0/(sqrt(W(ibf))+1.0e-10) 
  enddo
  x_matrix_small=matmul(matmul(U_mat,Tmp_matrix),transpose(conjg(U_mat)))
  x_matrix(nbasis_L+1:,nbasis_L+1:)=x_matrix_small(:,:) ! NOTE: save in x_matrix
   ! C = C (RKB S_SS)^-1/2 
  c_matrix_small=matmul(c_matrix_small,x_matrix_small)
   ! Check that it is orthonormal -> X^dagger (C_S_UKB)^dagger S_SS_UKB C_S_UKB X = X^dagger S_SS_RKB X
  x_matrix_small=matmul(transpose(conjg(c_matrix_small)),matmul(s_matrix_small,c_matrix_small)) ! (C_S)^dagger S_SS C_S = I
  do ibf=1,nbasis_L
   do jbf=1,nbasis_L
     if(ibf/=jbf) then
      if(abs(x_matrix_small(ibf,jbf))>1d-8) then
       write(stdout,'(a,i10,i10,2f20.5)') 'Error S comp. orthonorm. ',ibf,jbf,x_matrix_small(ibf,jbf)
      endif
     else
      if(abs(x_matrix_small(ibf,jbf)-1.0d0)>1d-8) then
       write(stdout,'(a,i10,i10,2f20.5)') 'Error S comp. orthonorm. ',ibf,jbf,x_matrix_small(ibf,jbf)
      endif
     endif
   enddo
  enddo
  deallocate(W,U_mat,Tmp_matrix) 
  call clean_deallocate('4C UKB overlap matrix S',s_matrix_4c)
  call clean_deallocate('Small overlap matrix ',s_matrix_small)
  call clean_deallocate('(RBK S_SS)^-1/2 matrix ',x_matrix_small)
  call clean_deallocate('Coefficients matrix small ',c_matrix_small)
  write(stdout,'(a,/)') ' Completed the orthonomalization of the RKB small component'

  !! Diagonalize the H_rel in RKB
   ! H^RKB C = S C e -> H_rel_rkb_mat c_matrix = s_matrix c_matrix e
   ! Building
   !          (S^-1/2)^dagger H^RKB S^-1/2 
   ! with S^-1/2 = x_matrix and ((S^-1/2)^dagger S = S^1/2 because ((S^-1/2)^dagger S S^-1/2 = I
  write(stdout,'(/,a)') ' Diagonalizing the H_rel in RKB'
  call clean_allocate('H_rel in RKB ortho',H_rel_rkb_ortho_mat,nstate_rkb,nstate_rkb)
  call clean_allocate('Full RKB wavefunctions C',c_matrix,nstate_rkb,nstate_rkb)
  H_rel_rkb_ortho_mat=matmul(conjg(transpose(x_matrix)),matmul(H_rel_rkb_mat,x_matrix))
  allocate(W(nstate_rkb),U_mat(nstate_rkb,nstate_rkb),E_state(nstate_rkb)) 
  call diagonalize(' ',H_rel_rkb_ortho_mat,W,U_mat)
  E_state(1:nbasis_L)=W(nbasis_L+1:2*nbasis_L)
  E_state(nbasis_L+1:2*nbasis_L)=W(1:nbasis_L) 
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
  call clean_deallocate('H_rel in RKB ortho',H_rel_rkb_ortho_mat)
  write(stdout,'(a,/)') ' Diagonalized the H_rel in RKB'
  ielectrons=nint(electrons_in)
  write(stdout,'(/,/,a25,1x,f19.10,/)') 'Rel Hcore Energy (Ha):',sum(E_state(1:ielectrons)) 

  !! NOTE: At this state we have all the fixed one-body contributions to H_rel (i.e. kinetic+external potential)
  !!       to do Dirac-HF/DFT-SCF 4c-calculations. We have a guess (i.e. c_matrix ) that can be used to build
  !!       the electronic density n(r) = ( MOs(r) 1-RDM MOs(s)^dagger ) with
  !!                              (  MOs(r) ) = ( AO^RKB (r) )*c_matrix
  !!                                          = ( AO^L+S (r) )*c_matrix_ukb2rkb*c_matrix
  !!       (  MOs(r) )_{1 x nbasis_L} is a row vector
  !!       (  AO^L+S (R) )_{1 x nbasis_L+nbasis_S} is a row vector
  !!       c_matrix_ukb2rkb _{nbasis_L+nbasis_S x 2*nbasis_L}
  !!       c_matrix _{2*nbasis_L x  2*nbasis_L}
  !!       Notice that in SCF calcs. we would only occupy the lowest (+) energy states [i.e. starting with nbasis_L+1]
  !!

  if(.not.is_x2c) then ! 4c-calculations

   deallocate(E_state)

   ! Set the crucial nstate
   nstate=nstate_rkb

   ! TODO this array should be passed to molgw.f90 for 4C calcs.
   call clean_deallocate('UKB to RKB coefficients',c_matrix_ukb2rkb)
   
   write(stdout,'(/,a)') ' Completed Relativistic Hamiltonian construction'
   write(stdout,'(a,/)') ' ==============================================='

  else ! X2C

   !! Build the decoupling matrix R from ->  C_L(+) (C_L(+))^dagger R = C_L(+) (C_S(+))^dagger ==> A R = B
    ! Rearrange the C matrix to have first the (+) energy states
    ! NOTE: Our decoupling is not fully cancelling C_L(-) because we apply it not to the the diag. H_rel_rkb_ortho_mat
    !       but to H C = S C e. Our
   write(stdout,'(/,a)') ' Computing X2C decoupling matrix R'
   allocate(Tmp_matrix(nstate_rkb,nstate_rkb))
   Tmp_matrix(:,1:nbasis_L)=c_matrix(:,nbasis_L+1:nstate_rkb)
   Tmp_matrix(:,nbasis_L+1:nstate_rkb)=c_matrix(:,1:nbasis_L)
   c_matrix=Tmp_matrix
   deallocate(Tmp_matrix)
   allocate(Tmp_matrix(nbasis_L,nbasis_L))
   call clean_allocate('A matrix ',A_mat,nbasis_L,nbasis_L)
   call clean_allocate('B matrix ',B_mat,nbasis_L,nbasis_L)
   call clean_allocate('R matrix ',R_mat,nbasis_L,nbasis_L)
   Tmp_matrix(:,:)=c_matrix(1:nbasis_L,1:nbasis_L) ! C^(+,L)
   A_mat=matmul(Tmp_matrix,transpose(conjg(Tmp_matrix)))
   R_mat(:,:)=c_matrix(nbasis_L+1:,1:nbasis_L)      ! C^(+,S)
   B_mat=matmul(Tmp_matrix,transpose(conjg(R_mat)))
   allocate(W(nbasis_L),U_mat(nbasis_L,nbasis_L))
   call diagonalize(' ',A_mat,W,U_mat)
   R_mat=COMPLEX_ZERO
   do ibf=1,nbasis_L
    if(abs(W(ibf))<1.0e-8) then
     write(stdout,'(a,i5,f20.8)') ' Eigenvalue lower than 1e-8 in A^-1',ibf,W(ibf)
    endif
    R_mat(ibf,ibf)=1.0d0/(W(ibf)+1.0e-10)
   enddo
   A_mat=matmul(matmul(U_mat,R_mat),transpose(conjg(U_mat)))
   R_mat=matmul(A_mat,B_mat)
   deallocate(W,U_mat)
!   lwork=-1
!   allocate(ipiv(nbasis_L),Work(1))
!   call zgetrf(nbasis_L,nbasis_L,A_mat,nbasis_L,ipiv,info)
!   if(info/=0) call die("Error in Relativistic computing A^-1 in zgetrf")
!   call zgetri(nbasis_L,A_mat,nbasis_L,ipiv,Work,lwork,info)
!   lwork=nint(real(Work(1)))
!   deallocate(Work)
!   allocate(Work(lwork))
!   call zgetri(nbasis_L,A_mat,nbasis_L,ipiv,Work,lwork,info)
!   R_mat=matmul(A_mat,B_mat)
    ! Check that R solves C_L(+) (C_L(+))^dagger R = C_L(+) (C_S(+))^dagger => A R = B
   Tmp_matrix(:,:)=c_matrix(1:nbasis_L,1:nbasis_L) ! C^(+,L)
   A_mat=matmul(Tmp_matrix,transpose(conjg(Tmp_matrix)))
   A_mat=matmul(A_mat,R_mat)
   A_mat=A_mat-B_mat ! A contains A R - B
   do ibf=1,nbasis_L
    do jbf=1,nbasis_L
      if(abs(A_mat(ibf,jbf))>1d-5) then
       write(stdout,'(a,i10,i10,2f20.5)') 'Error computing R decoup. matrix. ',ibf,jbf,A_mat(ibf,jbf)
      endif
    enddo
   enddo
   R_mat=transpose(conjg(R_mat))
   deallocate(Tmp_matrix)!,Work,ipiv)
   write(stdout,'(a)') ' Checked that the R matrix solves the system of equations'
   write(stdout,'(a,/)') ' Completed R matrix construction'
   
    ! Compute normalization matrices
    ! in A we put 1/ sqrt[ I + R^dagger R ]
    ! in B we put 1/ sqrt[ I + R R^dagger ]
   write(stdout,'(/,a)') ' Computing normalization factors for the transformation matrix'
   write(stdout,'(a)') ' 1/ sqrt[ I + R^dagger R ] '
   write(stdout,'(a)') ' 1/ sqrt[ I + R R^dagger ] '
   A_mat=matmul(transpose(conjg(R_mat)),R_mat)
   B_mat=matmul(R_mat,transpose(conjg(R_mat)))
   do ibf=1,nbasis_L
    A_mat(ibf,ibf)=A_mat(ibf,ibf)+1.0d0
    B_mat(ibf,ibf)=B_mat(ibf,ibf)+1.0d0
   enddo
   allocate(W(nbasis_L),U_mat(nbasis_L,nbasis_L),Tmp_matrix(nbasis_L,nbasis_L))
   Tmp_matrix=complex_zero; U_mat=complex_zero; W=complex_zero;
   call diagonalize(' ',A_mat,W,U_mat)
   do ibf=1,nbasis_L
    if(abs(W(ibf))<1.0e-8) then
     write(stdout,'(a,i5,f20.8)') ' Eigenvalue lower than 1e-8 in 1/ sqrt[ I + R^dagger R ]',ibf,W(ibf)
    endif
    Tmp_matrix(ibf,ibf)=1.0d0/(sqrt(W(ibf))+1.0e-10) 
   enddo  
   A_mat=matmul(matmul(U_mat,Tmp_matrix),transpose(conjg(U_mat)))
   Tmp_matrix=complex_zero; U_mat=complex_zero; W=complex_zero;
   call diagonalize(' ',B_mat,W,U_mat)
   do ibf=1,nbasis_L
    if(abs(W(ibf))<1.0e-8) then
     write(stdout,'(a,i5,f20.8)') ' Eigenvalue lower than 1e-8 in 1/ sqrt[ I + R R^dagger ]',ibf,W(ibf)
    endif
    Tmp_matrix(ibf,ibf)=1.0d0/(sqrt(W(ibf))+1.0e-10) 
   enddo  
   B_mat=matmul(matmul(U_mat,Tmp_matrix),transpose(conjg(U_mat)))
   deallocate(W,U_mat,Tmp_matrix)
   write(stdout,'(a,/)') ' Completed normalization factors for the transformation matrix'
   
   write(stdout,'(/,a)') ' Computing decoupling unitary matrix'
   call clean_allocate('U transformation matrix ',U_mat,nstate_rkb,nstate_rkb)
   allocate(Tmp_matrix(nstate_rkb,nstate_rkb))
   U_mat=complex_zero
   Tmp_matrix=complex_zero
    ! Normalization matrix
   Tmp_matrix(1:nbasis_L,1:nbasis_L)=A_mat(:,:)  
   Tmp_matrix(nbasis_L+1:,nbasis_L+1:)=B_mat(:,:)
    ! Decoupling matrix
   U_mat(nbasis_L+1:,1:nbasis_L)=-R_mat(:,:)
   R_mat=transpose(conjg(R_mat))
   U_mat(1:nbasis_L,nbasis_L+1:)=R_mat
   do ibf=1,nstate_rkb
    U_mat(ibf,ibf)=1.0d0
   enddo
   U_mat=matmul(Tmp_matrix,U_mat)
   Tmp_matrix=matmul(U_mat,transpose(conjg(U_mat)))
   do ibf=1,nbasis_rkb
    do jbf=1,nbasis_rkb
     if(ibf==jbf) then
      if(abs(Tmp_matrix(ibf,jbf)-1.0d0)>1d-8) then
       write(stdout,'(a,i5,i5,2f10.5)') ' Error in U matrix',ibf,jbf,s_matrix(ibf,jbf)
      endif
     else
      if(abs(Tmp_matrix(ibf,jbf))>1d-8) then
       write(stdout,'(a,i5,i5,2f10.5)') ' Error in U matrix',ibf,jbf,s_matrix(ibf,jbf)
      endif
     endif
    enddo
   enddo
   deallocate(Tmp_matrix)
   write(stdout,'(a,/)') ' Checked that the U decoupling matrix is unitary'
   call clean_deallocate('A matrix ',A_mat)
   call clean_deallocate('B matrix ',B_mat)
   call clean_deallocate('R matrix ',R_mat)
   write(stdout,'(a,/)') ' Completed decoupling unitary matrix construction'
   
   !! Transform to the X2C matrices (get ---> H^X2C, S^X2C, C^X2C, X^X2C, etc. to be used in SCF calcs.) 
   write(stdout,'(/,a)') ' Decoupling 4C -> X2C'
    ! H^RKB C = S C e 
    ! H^RKB  U^dagger  U C = S U^dagger   U C   e
    ! U H^RKB U^dagger U C = U S U^dagger U C   e
    !   H^x2c          C^x2c = S^x2c      C^x2c e    
   allocate(Tmp_matrix(nstate_rkb,nstate_rkb))
    ! C^x2c = U C
   Tmp_matrix=matmul(U_mat,c_matrix)
   deallocate(c_matrix)
   allocate(c_matrix(nbasis_L,nbasis_L))
   c_matrix(:,:)=Tmp_matrix(1:nbasis_L,1:nbasis_L)
    ! S^x2c = U S U^dagger
   Tmp_matrix=matmul(matmul(U_mat,s_matrix),transpose(conjg(U_mat)))
   deallocate(s_matrix)
   allocate(s_matrix(nbasis_L,nbasis_L))
   s_matrix(:,:)=Tmp_matrix(1:nbasis_L,1:nbasis_L)
    ! X = (S^-1/2)^x2c 
   deallocate(x_matrix)
   allocate(x_matrix(nbasis_L,nbasis_L))
   x_matrix=s_matrix
   allocate(W(nbasis_L),V_mat(nbasis_L,nbasis_L))
   call diagonalize(' ',x_matrix,W,V_mat)
   x_matrix=complex_zero
   do ibf=1,nbasis_L
    if(abs(W(ibf))<1.0e-8) then
     write(stdout,'(a,i5,f20.8)') ' Eigenvalue lower than 1e-8 in X matrix calc.',ibf,W(ibf)
    endif
    x_matrix(ibf,ibf)=1.0d0/(sqrt(W(ibf))+1.0e-10) 
   enddo
   x_matrix=matmul(matmul(V_mat,x_matrix),transpose(conjg(V_mat)))
   deallocate(W,V_mat)
    ! H^x2c = U H^RKB U^dagger
   Tmp_matrix=matmul(matmul(U_mat,H_rel_rkb_mat),transpose(conjg(U_mat)))
   deallocate(H_rel_rkb_mat)
   allocate(H_rel_rkb_mat(nbasis_L,nbasis_L))
   H_rel_rkb_mat(:,:)=Tmp_matrix(1:nbasis_L,1:nbasis_L)
   deallocate(Tmp_matrix)
   write(stdout,'(a,/)') ' Completed decoupling 4C -> X2C'

!!  Use this to show that H^X2C reproduces the (+) energy states written in c_matrix
!!  and test other properties of transformed operators
!!  block
!!  complex(dp),allocatable :: tmp0(:,:)
!!  complex(dp),allocatable :: tmp1(:,:)
!!  complex(dp),allocatable :: tmp2(:,:)
!!  allocate(tmp0(nbasis_L,nbasis_L))
!!  allocate(tmp1(nbasis_L,nbasis_L))
!!  allocate(tmp2(nbasis_L,nbasis_L))
!!  write(stdout,'(a)') ' Checking that H^x2c C^x2c = S^x2c C^x2c e(+) after decoupling'
!!  tmp2=complex_zero
!!  do ibf=1,nbasis_L
!!   tmp2(ibf,ibf)=E_state(ibf)
!!  enddo
!!  tmp0=matmul(H_rel_rkb_mat,c_matrix)
!!  tmp1=matmul(matmul(s_matrix,c_matrix),tmp2)
!!  tmp2=tmp0-tmp1
!!  do ibf=1,nbasis_L
!!   do jbf=1,nbasis_L
!!    if(abs(tmp2(ibf,jbf))>1d-8) write(stdout,*) ibf,jbf,tmp2(ibf,jbf)
!!   enddo
!!  enddo
!!  write(stdout,'(a)') ' Checking Hermiticity of H^x2c'
!!  do ibf=1,nbasis_L
!!   do jbf=1,nbasis_L
!!    if(abs(H_rel_rkb_mat(ibf,jbf)-conjg(H_rel_rkb_mat(jbf,ibf)))>1d-8) then
!!     write(stdout,*) ibf,jbf,H_rel_rkb_mat(ibf,jbf),H_rel_rkb_mat(jbf,ibf)
!!    endif
!!   enddo
!!  enddo
!!  write(stdout,'(a)') ' Checking Hermiticity of S^x2c'
!!  do ibf=1,nbasis_L
!!   do jbf=1,nbasis_L
!!    if(abs(s_matrix(ibf,jbf)-conjg(s_matrix(jbf,ibf)))>1d-8) then
!!     write(stdout,*) ibf,jbf,s_matrix(ibf,jbf),s_matrix(jbf,ibf)
!!    endif
!!   enddo
!!  enddo
!!  write(stdout,'(a)') ' Checking (X^x2c)^dagger S^x2c X^x2c = I'
!!  tmp0=matmul(conjg(transpose(x_matrix)),matmul(s_matrix,x_matrix))
!!  do ibf=1,nbasis_L
!!   do jbf=1,nbasis_L
!!    if(ibf/=jbf) then
!!     if(abs(tmp0(ibf,jbf))>1d-8) write(stdout,*) ibf,jbf,tmp0(ibf,jbf)
!!    else
!!     if(abs(tmp0(ibf,jbf)-1.0d0)>1d-8) write(stdout,*) ibf,jbf,tmp0(ibf,jbf)
!!    endif
!!   enddo
!!  enddo
!!  deallocate(tmp0)
!!  deallocate(tmp1)
!!  deallocate(tmp2)
!!  endblock    
  
   deallocate(E_state)
   
   ! TODO these arrays should be passed to molgw.f90 for full picture change
   call clean_deallocate('U decoupling matrix ',U_mat)
   call clean_deallocate('UKB to RKB coefficients',c_matrix_ukb2rkb)

   ! Set the crucial nstate
   nstate=nbasis_L
   write(stdout,'(/,a)') ' Diagonalizing the H^X2C (for testing)'
   call clean_allocate('H_X2C ortho',H_rel_rkb_ortho_mat,nstate,nstate)
   H_rel_rkb_ortho_mat=matmul(conjg(transpose(x_matrix)),matmul(H_rel_rkb_mat,x_matrix))
   allocate(U_mat(nstate,nstate),W(nstate)) 
   call diagonalize(' ',H_rel_rkb_ortho_mat,W,U_mat)
   write(stdout,'(1x,a)') '=== Energies ==='
   write(stdout,'(a)') '   #                               (Ha)                       &
   &                         (eV)      '
   write(stdout,'(a)') '                         bar                       ubar       &
   &               bar                       ubar '
   do ibf=1,nstate/2
    write(stdout,'(1x,i5,2(2(1x,f25.5)))') ibf,W(2*ibf-1),W(2*ibf),W(2*ibf-1)*Ha_eV,W(2*ibf)*Ha_eV 
    if(ibf==ielectrons/2) write(stdout,'(a)') '  --------------------------------------------------------------'
   enddo
   deallocate(W,U_mat)
   call clean_deallocate('H_X2C ortho',H_rel_rkb_ortho_mat)

   write(stdout,'(/,a)') ' Completed X2C Hamiltonian construction'
   write(stdout,'(a,/)') ' ======================================'

  endif


#endif

  call stop_clock(timing_relativistic)

end subroutine relativistic_init

!==================================================================
subroutine shuffle_complex(nbasis,is_large,matrix)
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

!
! Check deviation from the identity of ( C^x2c )^dagger S C^x2c with the actual overlap matrix S
! and overwrite S and X matrix if the MAE > 1e-6 to preserve orthonormality
!==================================================================
subroutine check_CdaggerSC_I(basis,c_matrix_rel,s_matrix_rel,x_matrix_rel,s_matrix,x_matrix)

  type(basis_set),intent(inout)  :: basis
  real(dp),intent(in)            :: s_matrix(:,:),x_matrix(:,:)
  complex(dp),intent(inout)      :: c_matrix_rel(:,:),s_matrix_rel(:,:),x_matrix_rel(:,:)
  !====
  integer                 :: istate,jstate,nstate
  real(dp)                :: err_x2c_coef
  complex(dp),allocatable :: tmp_matrix(:,:)
  !====
  
  nstate=2*basis%nbf

  write(stdout,'(/,a)') ' Checking (C^x2c)^dagger S C^x2c = I'
  allocate(tmp_matrix(nstate,nstate))
  tmp_matrix=COMPLEX_ZERO
  do istate=1,nstate/2
    do jstate=1,nstate/2
       tmp_matrix(2*istate-1,2*jstate-1)=s_matrix(istate,jstate)
       tmp_matrix(2*istate  ,2*jstate  )=s_matrix(istate,jstate)
    enddo
  enddo
  tmp_matrix=matmul(transpose(conjg(c_matrix_rel)),matmul(tmp_matrix,c_matrix_rel))
  err_x2c_coef=0.0_dp
  do istate=1,nstate
    do jstate=1,nstate
       if(istate==jstate) then
         if(abs(tmp_matrix(istate,jstate)-1.0_dp)>1e-6) then
           err_x2c_coef=err_x2c_coef+abs(tmp_matrix(istate,jstate)-1.0_dp)
         endif
       else
         if(abs(tmp_matrix(istate,jstate))>1e-6) then
           err_x2c_coef=err_x2c_coef+abs(tmp_matrix(istate,jstate))
         endif
       endif
    enddo
  enddo
  deallocate(tmp_matrix)
  err_x2c_coef=err_x2c_coef/(nstate*nstate)
  write(stdout,'(a,f10.6)') ' MAE in (C^x2c)^dagger S C^x2c = I',err_x2c_coef
  if(err_x2c_coef>1e-6) then
    ! We prefer to enforce orthonormality for the C^x2c states
    write(stdout,'(a)') ' The MAE > 1e-6, overwriting S and X matrices before doing the SCF procedure'
    s_matrix_rel=COMPLEX_ZERO
    x_matrix_rel=COMPLEX_ZERO
    do istate=1,nstate/2
      do jstate=1,nstate/2
         s_matrix_rel(2*istate-1,2*jstate-1)=s_matrix(istate,jstate)
         s_matrix_rel(2*istate  ,2*jstate  )=s_matrix(istate,jstate)
         x_matrix_rel(2*istate-1,2*jstate-1)=x_matrix(istate,jstate)
         x_matrix_rel(2*istate  ,2*jstate  )=x_matrix(istate,jstate)
      enddo
    enddo
  endif
  

end subroutine check_CdaggerSC_I
!==================================================================

end module m_relativistic
!==================================================================
