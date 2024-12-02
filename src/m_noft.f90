!==================================================================
! This file is part of MOLGW.
! Author: Mauricio Rodriguez-Mayorga
!
! This file contains
! - NOFT energy opt.
!=========================================================================
#include "molgw.h"
module m_noft
  use m_definitions
  use m_mpi
  use m_cart_to_pure
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_hamiltonian_tools
  use m_hamiltonian_onebody
  use m_hamiltonian_twobodies
  use m_hamiltonian_wrapper
  use m_noft_driver


  logical,parameter,private       :: noft_verbose=.FALSE.,noft_1_spin=.TRUE.
  logical                         :: noft_edft=.FALSE.,noft_fcidump_in=.FALSE.
  integer,private                 :: nstate_noft,nstate_frozen,irs_noft
  real(dp)                        :: ExcDFT,E_t_vext
  real(dp),allocatable,private    :: AhCORE(:,:),T_Vext(:)            ! hCORE matrix (T+Ven) in AO basis
  complex(dp),allocatable,private :: AhCORE_cmplx(:,:)
  type(basis_set),pointer         :: basis_pointer


contains

!=========================================================================
subroutine noft_energy(basis,occupation,Enoft,Vnn,Aoverlap,c_matrix,c_matrix_rel,hkin,hnuc,hkin_nuc_rel)
  implicit none

  type(basis_set),intent(in),target :: basis
  real(dp),intent(inout)    :: occupation(:,:)
  real(dp),intent(in)       :: Vnn
  real(dp),intent(out)      :: Enoft
  real(dp),intent(in),optional        :: Aoverlap(:,:)
  real(dp),intent(in),optional        :: hkin(:,:),hnuc(:,:)
  real(dp),intent(inout),optional     :: c_matrix(:,:,:)
  complex(dp),intent(in),optional     :: hkin_nuc_rel(:,:)
  complex(dp),intent(inout),optional  :: c_matrix_rel(:,:)
  !====
  logical                   :: file_exists=.false.
  integer                   :: istate,lwork,info,iao,jao
  integer                   :: iorb,iorb1,istat,nelectrons,iunit=366
  integer                   :: imethorb,imethocc,nstate_occ,nstate_beta,nstate_alpha,nstate_coupled
  integer                   :: iNOTupdateOCC,iNOTupdateORB,iprintdmn,iprintswdmn,iprintints,ireadOCC,ireadCOEF
  integer                   :: ireadFdiag,ireadGAMMAs,ista,inof
  real(dp)                  :: ran_num,coeff_old
  real(dp),allocatable      :: occ(:,:),energy(:,:),occ_print(:,:)
  real(dp),allocatable      :: NO_COEF(:,:)
  real(dp),allocatable      :: tmp_mat0(:,:),tmp_mat(:,:),Work(:)
  complex(dp),allocatable   :: tmp_mat0_cmplx(:,:),tmp_mat_cmplx(:,:)
  complex(dp),allocatable   :: NO_COEF_cmplx(:,:)
  character(len=100)        :: msgw
  character(len=200)        :: ofile_name
  !=====

  basis_pointer => basis

  ! Init clock
  call start_clock(timing_noft_energy)

  ! Write hearder and set name for $name.noft file
  write(ofile_name,'(2a)') trim(output_name),'noft'
  write(stdout,'(/,a)') ' =================================================='
  write(stdout,'(a)')   ' NOFT calculation'
  write(stdout,'(3a)')  ' writting NOFT results to ',trim(ofile_name),' file.'
  write(stdout,'(a)')   ' =================================================='
  write(stdout,'(/,a)') ' '

  Enoft = zero; occupation = zero; ExcDFT = zero;
  nstate_noft = SIZE(c_matrix,DIM=2) ! Number of lin. indep. molecular orbitals

  ! These varibles will remain fixed for a while
  imethorb=1;imethocc=1;iNOTupdateOCC=0;iNOTupdateORB=0;
  iprintdmn=0;iprintswdmn=0;iprintints=0;ireadOCC=0;ireadCOEF=0;
  ireadFdiag=0;ireadGAMMAs=0;ista=0;inof=8;
  if(noft_NOTupdateOCC=='yes') iNOTupdateOCC=1
  if(noft_NOTupdateORB=='yes') iNOTupdateORB=1
  if(noft_printdmn=='yes') iprintdmn=1
  if(noft_printswdmn=='yes') iprintswdmn=1
  if(noft_printints=='yes') iprintints=1
  if(noft_readOCC=='yes') ireadOCC=1
  if(noft_readCOEF=='yes') ireadCOEF=1
  if(noft_readFdiag=='yes') ireadFdiag=1
  if(noft_readGAMMAS=='yes') ireadGAMMAs=1
  if(noft_NR_OCC=='yes') imethocc=2
  if(noft_QC_ORB=='yes') imethorb=2

  select case(capitalize(noft_functional))
  case('GNOFS')
    inof=8
    ista=3
  case('GNOF')
    inof=8
  case('PNOF7S')
    inof=7
    ista=1
  case('PNOF7')
    inof=7
  case('PNOF5')
    inof=5
  case('MULLER')
    inof=100
  case('POWER')
    inof=101
  case('CA')
    inof=102
  case('CGA')
    inof=103
  case('GU')
    inof=104
  case('HF')
    inof=0
  case('PCCD')
    inof=-1
  case default
    call die('noft_energy: NOFT functional not recognized. Check your input')
  end select
  if(TRIM(postscf)=='PCCD') inof=-1  ! Last check that it is not pCCD

  noft_fcidump_in=(noft_fcidump=='yes')
  if( noft_fcidump_in .and. has_auxil_basis ) then
    noft_fcidump_in=.false.
    write(msgw,'(a)') 'The FCIDUMP file is not available with RI.'
    call issue_warning(msgw)
  endif
  if( noft_fcidump_in .and. noft_complex=='yes' ) then
    noft_fcidump_in=.false.
    write(msgw,'(a)') 'The FCIDUMP file is not available with complex orbitals.'
    call issue_warning(msgw)
  endif

  !
  ! Perform a relativistic or a non-relativistic NOFT calculations
  !
  if(TRIM(x2c) == 'yes') then ! relativistic

   nstate_noft=2*basis%nbf

   call clean_allocate('AhCORE_cmplx',AhCORE_cmplx,nstate_noft,nstate_noft)
   call clean_allocate('NO_COEF_cmplx',NO_COEF_cmplx,nstate_noft,nstate_noft)
   call clean_allocate('NO_energies',energy,nstate_noft,1)
   AhCORE_cmplx=hkin_nuc_rel; NO_COEF_cmplx=c_matrix_rel;

   ! Recover GUESS=core if required
   if(TRIM(init_hamiltonian)=='CORE') then
     call clean_allocate('tmp_mat0_cmplx',tmp_mat0_cmplx,nstate_noft,nstate_noft,noft_verbose)
     call clean_allocate('tmp_mat_cmplx',tmp_mat_cmplx,nstate_noft,nstate_noft,noft_verbose)
     tmp_mat_cmplx=matmul(conjg(transpose(NO_COEF_cmplx)),matmul(AhCORE_cmplx,NO_COEF_cmplx)) ! Hcore^MO basis
     call diagonalize(' ',tmp_mat_cmplx,energy(:,1),tmp_mat0_cmplx)
     NO_COEF_cmplx=matmul(NO_COEF_cmplx,tmp_mat0_cmplx) ! rotate to the basis where Hcore (=Hcore^X2C) is diagonal
     call clean_deallocate('tmp_mat0_cmplx',tmp_mat0_cmplx,noft_verbose)
     call clean_deallocate('tmp_mat_cmplx',tmp_mat_cmplx,noft_verbose)
   endif

   ! TODO: call noft suborutine from module standalone, do SCF, etc, etc
   c_matrix_rel=NO_COEF_cmplx
    
   ! Clean arrays
   call clean_deallocate('NO_energies',energy)
   call clean_deallocate('AhCORE_cmplx',AhCORE_cmplx)
   call clean_deallocate('NO_COEF_cmplx',NO_COEF_cmplx)

  else ! non-relativistic

   ! Allocate arrays and initialize them
   if(noft_complex=='yes') then
     call clean_allocate('AhCORE_cmplx',AhCORE_cmplx,basis%nbf,basis%nbf)
   else
     call clean_allocate('AhCORE',AhCORE,basis%nbf,basis%nbf)
   endif
   call clean_allocate('NO_occ',occ,basis%nbf,1)
   call clean_allocate('NO_energies',energy,basis%nbf,1)
   if(noft_complex=='yes') then
     call clean_allocate('NO_COEF_cmplx',NO_COEF_cmplx,basis%nbf,basis%nbf)
   else
     call clean_allocate('NO_COEF',NO_COEF,basis%nbf,basis%nbf)
   endif
   occ(:,:)    = zero
   energy(:,:) = zero
   ! Save Atomic Orbital hCORE integrals
   if(noft_complex=='yes') then
     AhCORE_cmplx(:,:) = hkin(:,:) + hnuc(:,:)
   else
     AhCORE(:,:) = hkin(:,:) + hnuc(:,:)
   endif
   
   ! Initially copy c_matrix (HF orbs) to NO_COEF
   if(noft_complex=='yes') then
     NO_COEF_cmplx(:,:)=complex_zero
     write(stdout,'(/,a)') ' Adding Random Imaginary Phases '
     write(stdout,'(a,/)') ' ------------------------------ '
     do istate=1,nstate_noft
       call random_number(ran_num) ! For complex orbs, each one has its own random phase (to have real and imaginary orbs)
       if(noft_nophases=='yes') ran_num=0.0e0
       write(stdout,'(a,I10,a,f8.5,a,f8.5,a)') ' MO',istate,': (',real(exp(im*ran_num)),',',aimag(exp(im*ran_num)),')'
       NO_COEF_cmplx(:,istate)=exp(im*ran_num)*c_matrix(:,istate,1)
     enddo
     write(stdout,*) ' '
   else
     NO_COEF(:,:)=zero
     do istate=1,nstate_noft
       NO_COEF(:,istate)=c_matrix(:,istate,1)
     enddo
   endif
   ! Replace NO_COEF by Hcore orbs for the initial GUESS?
   if(TRIM(init_hamiltonian)=='CORE') then
     call clean_allocate('tmp_mat0',tmp_mat0,basis%nbf,basis%nbf,noft_verbose)
     call clean_allocate('tmp_mat',tmp_mat,basis%nbf,basis%nbf,noft_verbose)
     tmp_mat0(:,:)=zero
     do istate=1,nstate_noft
       tmp_mat0(:,istate)=c_matrix(:,istate,1)
     enddo
     allocate(Work(1))
     tmp_mat=matmul(hkin,tmp_mat0)+matmul(hnuc,tmp_mat0)
     tmp_mat=matmul(transpose(tmp_mat0),tmp_mat)
     lwork=-1
     call DSYEV('V','L',basis%nbf,tmp_mat,basis%nbf,energy(:,1),Work,lwork,info)
     lwork=nint(Work(1))
     if(info==0) then
       deallocate(Work)
       allocate(Work(lwork))
       energy=zero
       call DSYEV('V','L',basis%nbf,tmp_mat,basis%nbf,energy(:,1),Work,lwork,info)
     endif
     if(noft_complex=='yes') then
       NO_COEF_cmplx=matmul(NO_COEF_cmplx,tmp_mat)
       write(stdout,'(/,a)') ' Adding New Random Imaginary Phases '
       write(stdout,'(a,/)') ' ---------------------------------- '
       do istate=1,nstate_noft
         call random_number(ran_num) ! For complex orbs, each one has its own random phase (to have real and imaginary orbs)
         if(noft_nophases=='yes') ran_num=0.0e0
         write(stdout,'(a,I10,a,f8.5,a,f8.5,a)') ' MO',istate,': (',real(exp(im*ran_num)),',',aimag(exp(im*ran_num)),')'
         NO_COEF_cmplx(:,istate)=exp(im*ran_num)*NO_COEF_cmplx(:,istate)
       enddo
       write(stdout,*) ' '
     else
       NO_COEF=matmul(NO_COEF,tmp_mat)
     endif
     write(stdout,'(/,a,/)') ' Approximate Hamiltonian Hcore used as GUESS in NOFT calc.'
     call clean_deallocate('tmp_mat0',tmp_mat0,noft_verbose)
     call clean_deallocate('tmp_mat',tmp_mat,noft_verbose)
     deallocate(Work)
   endif
   
   ! Build NO_COEF_cmplx with phases times NO_COEF (real)
   if(TRIM(init_hamiltonian)=='NOFT' .and. noft_complex=='yes') then
     inquire(file='NO_COEF_BIN',exist=file_exists)
     if( file_exists ) then
       open(unit=iunit,form='unformatted',file='NO_COEF_BIN',iostat=istat,status='old')
       if(istat==0) then
        do
         read(iunit,iostat=istat) iorb,iorb1,coeff_old
         if(istat/=0) then
          exit
         endif
         if(((iorb/=0).and.(iorb1/=0)).and.iorb*iorb1<=basis%nbf*basis%nbf) then
          NO_COEF_cmplx(iorb,iorb1)=coeff_old
         else
          exit
         endif
        enddo
       endif
       write(stdout,'(/,a)') ' Adding New Random Imaginary Phases '
       write(stdout,'(a,/)') ' ---------------------------------- '
       do istate=1,nstate_noft
         call random_number(ran_num) ! For complex orbs, each one has its own random phase (to have real and imaginary orbs)
         if(noft_nophases=='yes') ran_num=0.0e0
         write(stdout,'(a,I10,a,f7.5,a,f7.5,a)') ' MO',istate,': (',real(exp(im*ran_num)),',',aimag(exp(im*ran_num)),')'
         NO_COEF_cmplx(:,istate)=exp(im*ran_num)*NO_COEF_cmplx(:,istate)
       enddo
       write(stdout,*) ' '
       write(stdout,'(/,a,/)') ' Reading the NO_COEF_BIN file to set the initial complex NO_COEF.'
     else
       write(stdout,'(/,a,/)') ' Did not find NO_COEF_BIN file to set the initial complex NO_COEF.'
     endif
   endif
   
   ! Not ready for open-shell calcs. (TODO)
   nelectrons = NINT(electrons)
   nstate_coupled = noft_ncoupled - 1
   nstate_frozen = (nelectrons-2*noft_npairs)/2
   nstate_beta = nstate_frozen+noft_npairs
   nstate_alpha = nstate_beta
   do
     nstate_occ=nstate_frozen+noft_npairs*(nstate_coupled+1)
     if(nstate_occ<=nstate_noft) then
       exit
     else
       nstate_coupled=nstate_coupled-1
     endif
   enddo
   
   !
   ! Setup the grids for the quadrature of DFT potential/energy
   irs_noft=0
   if( calc_type%is_dft .and. noft_dft=='yes' ) then
     if( nspin /= 2 ) call die('molgw: RS-NOFT calculations need nspin=2')
     if(noft_rsinter=='yes') then
       irs_noft=1
     else
       irs_noft=2
     endif
     if( .not.calc_type%need_exchange_lr ) then
       write(msgw,'(a)') 'LR exchange is needed for RS-NOFT.'
       call die(msgw)
     endif
     write(stdout,'(a)') ' '
     write(stdout,'(a,f10.5)') ' RS-NOFT amount of exact exchange (alpha_hybrid)           ',alpha_hybrid
     write(stdout,'(a,f10.5)') ' RS-NOFT amount of long-range exact exchange (beta_hybrid) ',beta_hybrid
     write(stdout,'(a,f10.5)') ' RS-NOFT error function parameter (gamma_hybrid)           ',gamma_hybrid
     call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
   endif
   
   ! Call module initialization and run NOFT calc.
   if(noft_complex=='yes') then
   
     if(noft_restart=='yes') then
       call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
        imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
        Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF_cmplx=NO_COEF_cmplx,lowmemERI=(noft_lowmemERI=='yes'),&
        restart=(noft_restart=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,&
        ireadFdiag=ireadFdiag,iNOTupdateOCC=iNOTupdateOCC,iNOTupdateORB=iNOTupdateORB,Lpower=noft_Lpower,&
        fcidump=noft_fcidump_in,irange_sep=irs_noft,hessian=(noft_hessian=='yes'))
     else
       call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
        imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
        Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF_cmplx=NO_COEF_cmplx,lowmemERI=(noft_lowmemERI=='yes'),&
        Lpower=noft_Lpower,fcidump=noft_fcidump_in,irange_sep=irs_noft,hessian=(noft_hessian=='yes'))
     endif
   
   else
   
     if(noft_restart=='yes') then
       call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
        imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
        Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),&
        restart=(noft_restart=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,&
        ireadFdiag=ireadFdiag,iNOTupdateOCC=iNOTupdateOCC,iNOTupdateORB=iNOTupdateORB,Lpower=noft_Lpower,&
        fcidump=noft_fcidump_in,irange_sep=irs_noft,hessian=(noft_hessian=='yes'))
     else
       call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
        imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
        Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),&
        Lpower=noft_Lpower,fcidump=noft_fcidump_in,irange_sep=irs_noft,hessian=(noft_hessian=='yes'))
     endif
   
   endif
   
   if( irs_noft/=0 ) then ! Compute total Energy for range-sep NOFT switching off the hamiltonian_xc (we only need hCORE=T+Vext).
     noft_edft=.true.     ! So, we restart but we will not update orbs nor occs.
     call clean_allocate('T_Vext',T_Vext,basis%nbf,noft_verbose)
     write(ofile_name,'(a)') 'tmp_dft_noft'
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
      imethocc,imethorb,noft_nscf,0,0,0,noft_ithresh_lambda,noft_ndiis,Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),&
      mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),restart=.true.,ireadGAMMAS=1,ireadOCC=0,&
      ireadCOEF=1,ireadFdiag=1,iNOTupdateOCC=1,iNOTupdateORB=1,Lpower=noft_Lpower,fcidump=(noft_fcidump=='yes'),irange_sep=irs_noft)
     Enoft=Enoft+ExcDFT
     write(stdout,'(/,a,2x,f19.10)')   ' Nucleus-Nucleus (Ha):',Vnn
     write(stdout,'(a,2x,f19.10)')     ' Hcore Energy (Ha)   :',sum(T_Vext(:)*occ(:,1))
     write(stdout,'(a,2x,f19.10,/)')   ' XC Energy (Ha)      :',ExcDFT
     call system("rm tmp_dft_noft")
     call clean_deallocate('T_Vext',T_Vext,noft_verbose)
     call destroy_dft_grid()
   endif
   
   ! If required print post-procesing files
   occupation(1:nstate_occ,1)=occ(1:nstate_occ,1)
   if( irs_noft==0 ) then
     if(noft_complex=='yes') then
       if(print_wfn_files_ ) then
         call clean_allocate('Occ_print',occ_print,nstate_noft,1,noft_verbose)
         occ_print(1:nstate_noft,1)=occ(1:nstate_noft,1)
         ! Update c_matrix with real part of optimized NO_COEF
         do istate=1,nstate_noft
           c_matrix(:,istate,1)=real(NO_COEF_cmplx(:,istate))
         enddo
         call print_wfn_file('NOFT_RE',basis,occ_print,c_matrix,Enoft,energy)
         ! Update c_matrix with imaginary part of optimized NO_COEF
         do istate=1,nstate_noft
           c_matrix(:,istate,1)=aimag(NO_COEF_cmplx(:,istate))
         enddo
         call print_wfn_file('NOFT_IM',basis,occ_print,c_matrix,Enoft,energy)
         call clean_deallocate('Occ_print',occ_print,noft_verbose)
       endif
     else
       ! Update c_matrix with optimized NO_COEF
       do istate=1,nstate_noft
         c_matrix(:,istate,1)=NO_COEF(:,istate)
       enddo
       ! Select the post-procesing files
       if(print_wfn_ .or. print_cube_ .or. print_wfn_files_ ) then
         call clean_allocate('Occ_print',occ_print,nstate_noft,1,noft_verbose)
         occ_print(1:nstate_noft,1)=occ(1:nstate_noft,1)
         if( print_wfn_ )  call plot_wfn(basis,c_matrix)
         if( print_wfn_ )  call plot_rho('NOFT',basis,occ_print,c_matrix)
         if( print_cube_ ) call plot_cube_wfn('NOFT',basis,occ_print,c_matrix)
         if( print_wfn_files_ ) call print_wfn_file('NOFT',basis,occ_print,c_matrix,Enoft,energy)
         call clean_deallocate('Occ_print',occ_print,noft_verbose)
       endif
     endif
   endif
   
   ! Deallocate arrays and print the normal termination
   if(noft_complex=='yes') then
     call clean_deallocate('AhCORE_cmplx',AhCORE_cmplx)
   else
     call clean_deallocate('AhCORE',AhCORE)
   endif
   call clean_deallocate('NO_occ',occ)
   call clean_deallocate('NO_energies',energy)
   if(noft_complex=='yes') then
     call clean_deallocate('NO_COEF_cmplx',NO_COEF_cmplx)
   else
     call clean_deallocate('NO_COEF',NO_COEF)
   endif

  endif

  write(stdout,'(/,a)') ' =================================================='
  write(stdout,'(a)')   ' NOFT SCF loop ends here'
  write(stdout,'(a)')   ' =================================================='
  write(stdout,'(/,a)') ' '

  ! Stop clock
  call stop_clock(timing_noft_energy)

end subroutine noft_energy


!==================================================================
subroutine mo_ints(nbf,nstate_occ,nstate_kji,Occ,DM2_JK,NO_COEF,hCORE,ERImol,ERImolJsr,ERImolLsr,&
     &             NO_COEF_cmplx,hCORE_cmplx,ERImol_cmplx,ERImolJsr_cmplx,ERImolLsr_cmplx,all_ERIs,&
     &             Edft_xc,do_xc_dft)
  implicit none

  logical,optional,intent(in)     :: all_ERIs,do_xc_dft
  integer,intent(in)              :: nbf,nstate_occ,nstate_kji
  real(dp),intent(in)             :: Occ(nstate_occ)
  real(dp),optional,intent(inout) :: Edft_xc
  real(dp),optional,intent(in)    :: DM2_JK(2,nstate_occ,nstate_occ)
  real(dp),optional,intent(in)    :: NO_COEF(nbf,nbf)
  real(dp),optional,intent(inout) :: hCORE(nbf,nbf)
  real(dp),optional,intent(inout) :: ERImol(nbf,nstate_kji,nstate_kji,nstate_kji)
  real(dp),optional,intent(inout) :: ERImolJsr(nbf,nstate_kji,nstate_kji)
  real(dp),optional,intent(inout) :: ERImolLsr(nbf,nstate_kji,nstate_kji)
  complex(dp),optional,intent(in)    :: NO_COEF_cmplx(nbf,nbf)
  complex(dp),optional,intent(inout) :: hCORE_cmplx(nbf,nbf)
  complex(dp),optional,intent(inout) :: ERImol_cmplx(nbf,nstate_kji,nstate_kji,nstate_kji)
  complex(dp),optional,intent(inout) :: ERImolJsr_cmplx(nbf,nstate_kji,nstate_kji)
  complex(dp),optional,intent(inout) :: ERImolLsr_cmplx(nbf,nstate_kji,nstate_kji)
  !====
  logical                    :: all_ERIs_in=.false.,long_range=.true.,do_xc_dft_tmp=.true.
  integer                    :: istate,jstate,pstate,qstate,ispin
  character(len=100)         :: msgw
  real(dp)                   :: ERI_lkji
  complex(dp)                :: ERI_lkji_cmplx
  real(dp),allocatable       :: occupation(:,:)
  real(dp),allocatable       :: tmp_c_matrix(:,:,:),hamiltonian_xc(:,:,:)
  complex(dp),allocatable    :: tmp_c_matrix_cmplx(:,:,:)
  !=====

  if(present(all_ERIs)) all_ERIs_in=all_ERIs

  ! Comment: Despite the arrays are of size nbf x nbf, we use nstate_noft = num. lin. indep. states in the ERI  transformation.
  ! Doing this, we save some time in the loops because nstate_noft <= nbf

  ! Molecular hCORE including all one-body interactions
  if(noft_complex=='yes') then

    ! Build 3D array for complex c_matrix and init hCORE_cmplx
    call clean_allocate('tmp_c_matrix',tmp_c_matrix_cmplx,nbf,nstate_noft,nspin,noft_verbose)
    hCORE_cmplx(:,:)=complex_zero;tmp_c_matrix_cmplx(:,:,:)=complex_zero;
    do istate=1,nstate_noft
      do ispin=1,nspin
        tmp_c_matrix_cmplx(:,istate,ispin)=NO_COEF_cmplx(:,istate)
      enddo
    enddo

    ! Add the sr-NOFT term
    if( noft_NOTvxc=='yes ' ) do_xc_dft_tmp=.false.
    if( present(do_xc_dft) )  do_xc_dft_tmp=.true.
    if( (irs_noft/=0) .and. (.not.noft_edft .and. do_xc_dft_tmp) ) then
      ! Prepare the DFT contribution (takes part only during orb. optimization and is switched off for final energy calculation)
      call clean_allocate('occupation',occupation,nbf,nspin,noft_verbose)
      call clean_allocate('hamiltonian_xc',hamiltonian_xc,nbf,nbf,nspin,noft_verbose)
      ! MRM: The first call of mo_ints contains occ(1:Nfrozen+Npairs)=2.0 
      occupation(:,:)=zero; hamiltonian_xc(:,:,:)=zero;
      if( ANY(Occ(:nstate_occ)>completely_empty) ) then
        if ( nspin==1 ) then ! In principle, this option should not be used because we need nspin=2 to use Pi(r)
          occupation(:nstate_occ,1)=2.0e0*Occ(:nstate_occ)
        else
          do istate=1,nstate_occ
            do ispin=1,nspin
              occupation(istate,ispin)=Occ(istate)
            enddo
          enddo
        endif
        call dft_exc_vxc_batch(BATCH_SIZE,basis_pointer,occupation,tmp_c_matrix_cmplx,hamiltonian_xc,ExcDFT,dm2_JK=DM2_JK)
        if(present(Edft_xc)) Edft_xc=ExcDFT
      endif   
      hamiltonian_xc(:,:,1)=SUM(hamiltonian_xc(:,:,:),DIM=3)
      if ( nspin==2 ) hamiltonian_xc(:,:,1)=0.5e0*hamiltonian_xc(:,:,1)
      hCORE_cmplx=matmul(conjg(transpose(NO_COEF_cmplx(:,:))),matmul(hamiltonian_xc(:,:,1),NO_COEF_cmplx(:,:)))
      call clean_deallocate('hamiltonian_xc',hamiltonian_xc,noft_verbose)
      call clean_deallocate('occupation',occupation,noft_verbose)
    endif      

    ! T+Vext+V_xc(?) part
    hCORE_cmplx=hCORE_cmplx+matmul(conjg(transpose(NO_COEF_cmplx)),matmul(AhCORE_cmplx,NO_COEF_cmplx))
    if( noft_edft ) then
      do istate=1,nstate_noft
        T_Vext(istate)=real(hCORE_cmplx(istate,istate))
      enddo
    endif

  else

    ! Build 3D array for c_matrix and init hCORE
    call clean_allocate('tmp_c_matrix',tmp_c_matrix,nbf,nstate_noft,nspin,noft_verbose)
    hCORE(:,:)=zero;tmp_c_matrix(:,:,:)=zero;
    do istate=1,nstate_noft
      do ispin=1,nspin
        tmp_c_matrix(:,istate,ispin)=NO_COEF(:,istate)
      enddo
    enddo

    ! Add the sr-NOFT term
    if( noft_NOTvxc=='yes ' ) do_xc_dft_tmp=.false.
    if( present(do_xc_dft) )  do_xc_dft_tmp=.true.
    if( (irs_noft/=0) .and. (.not.noft_edft .and. do_xc_dft_tmp) ) then
      ! Prepare the DFT contribution (takes part only during orb. optimization and is switched off for final energy calculation)
      call clean_allocate('occupation',occupation,nbf,nspin,noft_verbose)
      call clean_allocate('hamiltonian_xc',hamiltonian_xc,nbf,nbf,nspin,noft_verbose)
      ! MRM: The first call of mo_ints contains occ(1:Nfrozen+Npairs)=2.0
      occupation(:,:)=zero; hamiltonian_xc(:,:,:)=zero;
      if( ANY(Occ(:nstate_occ)>completely_empty) ) then
        if ( nspin==1 ) then ! In principle, this option should not be used because we need nspin=2 to use Pi(r)
          occupation(:nstate_occ,1)=2.0e0*Occ(:nstate_occ)
        else
          do istate=1,nstate_occ
            do ispin=1,nspin
              occupation(istate,ispin)=Occ(istate)
            enddo
          enddo
        endif
        call dft_exc_vxc_batch(BATCH_SIZE,basis_pointer,occupation,tmp_c_matrix,hamiltonian_xc,ExcDFT,dm2_JK=DM2_JK)
        if(present(Edft_xc)) Edft_xc=ExcDFT
      endif
      hamiltonian_xc(:,:,1)=SUM(hamiltonian_xc(:,:,:),DIM=3)
      if ( nspin==2 ) hamiltonian_xc(:,:,1)=0.5e0*hamiltonian_xc(:,:,1)
      hCORE=matmul(transpose(NO_COEF(:,:)),matmul(hamiltonian_xc(:,:,1),NO_COEF(:,:)))
      call clean_deallocate('hamiltonian_xc',hamiltonian_xc,noft_verbose)
      call clean_deallocate('occupation',occupation,noft_verbose)
    endif

    ! T+Vext+V_xc(?) part
    hCORE=hCORE+matmul(transpose(NO_COEF),matmul(AhCORE,NO_COEF))
    if( noft_edft ) then
      do istate=1,nstate_noft
        T_Vext(istate)=hCORE(istate,istate)
      enddo
    endif

  endif

  ! Molecular ERImol including all two-body interactions (maybe also including sr-ERImol)
  if( irs_noft/=0 ) then

    if(present(ERImol) .and. present(ERImolJsr) .and. present(ERImolLsr)) then
      ERImol(:,:,:,:)=zero; ERImolJsr(:,:,:)=zero; ERImolLsr(:,:,:)=zero
      if(has_auxil_basis) then ! RI case
        call calculate_eri_3center_eigen(tmp_c_matrix,1,nstate_noft,1,nstate_kji,verbose=noft_verbose,long_range=long_range, &
         &   only_one_spin=noft_1_spin)
        ! <lk| [alpha+beta*erf(gamma r12)]/r12 |ji> format used for ERImol
        ! Hartree : <li|ji>^sr = <li| 1/r12 |ji> - <li| [alpha+beta*erf(gamma r12)]/r12 |ji>
        ! Time-rev: <lk|ii>^sr = <lk| 1/r12 |ii> - <lk| [alpha+beta*erf(gamma r12)]/r12 |ii>
        ! Exchange: Not needed a <li|ij>^sr term to be passed to the NOFT module.
        do istate=1,nstate_occ
          do jstate=1,nstate_occ
            do pstate=1,nstate_noft
              ! Hartree: <pi|ji> format used for ERImol
              ERI_lkji=eri_eigen_ri(pstate,jstate,1,istate,istate,1)
              ERImol(pstate,istate,jstate,istate)=alpha_hybrid*ERI_lkji &
               +beta_hybrid*eri_eigen_ri_lr(pstate,jstate,1,istate,istate,1)
              ERImolJsr(pstate,istate,jstate)=ERI_lkji-ERImol(pstate,istate,jstate,istate)
              ! Exchange: <pj|ji> format used for ERImol
              ERI_lkji=eri_eigen_ri(pstate,jstate,1,jstate,istate,1)
              ERImol(pstate,jstate,jstate,istate)=alpha_hybrid*ERI_lkji &
               +beta_hybrid*eri_eigen_ri_lr(pstate,jstate,1,jstate,istate,1)
              ! Time-rev: <pi|jj> format used for ERImol
              ERI_lkji=eri_eigen_ri(pstate,jstate,1,istate,jstate,1)
              ERImol(pstate,istate,jstate,jstate)=alpha_hybrid*ERI_lkji &
               +beta_hybrid*eri_eigen_ri_lr(pstate,jstate,1,istate,jstate,1)
              ERImolLsr(pstate,istate,jstate)=ERI_lkji-ERImol(pstate,istate,jstate,jstate)
            enddo
          enddo
        enddo
        call destroy_eri_3center_eigen(verbose=noft_verbose,long_range=long_range)
      else            ! Normal case (not using RI)
        !TODO
        write(msgw,'(a)') 'LR exchange requires RI for RS-NOFT (hint: include the RI basis).'
        call die(msgw)
      endif
    endif

    if(present(ERImol_cmplx) .and. present(ERImolJsr_cmplx) .and. present(ERImolLsr_cmplx)) then
      ERImol_cmplx(:,:,:,:)=complex_zero; ERImolJsr_cmplx(:,:,:)=complex_zero; ERImolLsr_cmplx(:,:,:)=complex_zero
      if(has_auxil_basis) then ! RI case
        call calculate_eri_3center_eigen_cmplx(tmp_c_matrix_cmplx,1,nstate_noft,1,nstate_kji,verbose=noft_verbose,&
         &   long_range=long_range,only_one_spin=noft_1_spin)
        ! <lk| [alpha+beta*erf(gamma r12)]/r12 |ji> format used for ERImol
        ! Hartree : <li|ji>^sr = <li| 1/r12 |ji> - <li| [alpha+beta*erf(gamma r12)]/r12 |ji>
        ! Time-rev: <lk|ii>^sr = <lk| 1/r12 |ii> - <lk| [alpha+beta*erf(gamma r12)]/r12 |ii>
        ! Exchange: Not needed a <li|ij>^sr term to be passed to the NOFT module.
        do istate=1,nstate_occ
          do jstate=1,nstate_occ
            do pstate=1,nstate_noft
              ! Hartree: <pi|ji> format used for ERImol
              ERI_lkji_cmplx=eri_eigen_ri_cmplx(pstate,jstate,1,istate,istate,1)
              ERImol_cmplx(pstate,istate,jstate,istate)=alpha_hybrid*ERI_lkji_cmplx &
               +beta_hybrid*eri_eigen_ri_lr_cmplx(pstate,jstate,1,istate,istate,1)
              ERImolJsr_cmplx(pstate,istate,jstate)=ERI_lkji_cmplx-ERImol_cmplx(pstate,istate,jstate,istate)
              ! Exchange: <pj|ji> format used for ERImol
              ERI_lkji_cmplx=eri_eigen_ri_cmplx(pstate,jstate,1,jstate,istate,1)
              ERImol_cmplx(pstate,jstate,jstate,istate)=alpha_hybrid*ERI_lkji_cmplx &
               +beta_hybrid*eri_eigen_ri_lr_cmplx(pstate,jstate,1,jstate,istate,1)
              ! Time-rev: <pi|jj> format used for ERImol 
              ERI_lkji_cmplx=eri_eigen_ri_cmplx(pstate,jstate,1,istate,jstate,1)
              ERImol_cmplx(pstate,istate,jstate,jstate)=alpha_hybrid*ERI_lkji_cmplx &
               +beta_hybrid*eri_eigen_ri_lr_cmplx(pstate,jstate,1,istate,jstate,1)
              ERImolLsr_cmplx(pstate,istate,jstate)=ERI_lkji-ERImol_cmplx(pstate,istate,jstate,jstate)
            enddo
          enddo
        enddo
        call destroy_eri_3center_eigen_cmplx(verbose=noft_verbose,long_range=long_range)
      else            ! Normal case (not using RI)
        !TODO
        write(msgw,'(a)') 'LR exchange requires RI for RS-NOFT (hint: include the RI basis).'
        call die(msgw)
      endif
    endif

  else

    if(noft_complex=='yes') then

      if(present(ERImol_cmplx)) then
        ERImol_cmplx(:,:,:,:)=complex_zero
        if(has_auxil_basis) then ! RI case
          call calculate_eri_3center_eigen_cmplx(tmp_c_matrix_cmplx,1,nstate_noft,1,nstate_kji,verbose=noft_verbose)
          if(all_ERIs_in .and. nstate_noft==nstate_kji) then
           do istate=1,nstate_noft
            do jstate=1,nstate_noft
             do pstate=1,nstate_noft
              do qstate=1,nstate_noft
               ERImol_cmplx(qstate,pstate,jstate,istate)=eri_eigen_ri_cmplx(qstate,jstate,1,pstate,istate,1)
              enddo
             enddo
            enddo
           enddo
          else
           do istate=1,nstate_occ
             do jstate=1,nstate_occ
               do pstate=1,nstate_noft
                 ERImol_cmplx(pstate,istate,jstate,istate)=eri_eigen_ri_cmplx(pstate,jstate,1,istate,istate,1) ! <li|ji> format used for ERImol
                 ERImol_cmplx(pstate,jstate,jstate,istate)=eri_eigen_ri_cmplx(pstate,jstate,1,jstate,istate,1) ! <lj|ji> format used for ERImol
                 ERImol_cmplx(pstate,istate,jstate,jstate)=eri_eigen_ri_cmplx(pstate,jstate,1,istate,jstate,1) ! <li|jj> format used for ERImol
               enddo
             enddo
           enddo
          endif
          call destroy_eri_3center_eigen_cmplx(noft_verbose)
        else            ! Normal case (not using RI)
          call form_erimol(nbf,nstate_noft,nstate_kji,c_matrix_cmplx=tmp_c_matrix_cmplx,ERImol_cmplx=ERImol_cmplx)
        endif
      endif

    else

      if(present(ERImol)) then
        ERImol(:,:,:,:)=zero
        if(has_auxil_basis) then ! RI case
          call calculate_eri_3center_eigen(tmp_c_matrix,1,nstate_noft,1,nstate_kji,verbose=noft_verbose)
          if(all_ERIs_in .and. nstate_noft==nstate_kji) then
           do istate=1,nstate_noft
            do jstate=1,nstate_noft
             do pstate=1,nstate_noft
              do qstate=1,nstate_noft
               ERImol(qstate,pstate,jstate,istate)=eri_eigen_ri(qstate,jstate,1,pstate,istate,1)
              enddo
             enddo
            enddo
           enddo
          else
           do istate=1,nstate_occ
             do jstate=1,nstate_occ
               do pstate=1,nstate_noft
                 ERImol(pstate,istate,jstate,istate)=eri_eigen_ri(pstate,jstate,1,istate,istate,1) ! <li|ji> format used for ERImol
                 ERImol(pstate,jstate,jstate,istate)=eri_eigen_ri(pstate,jstate,1,jstate,istate,1) ! <lj|ji> format used for ERImol
                 ERImol(pstate,istate,jstate,jstate)=eri_eigen_ri(pstate,jstate,1,istate,jstate,1) ! <li|jj> format used for ERImol
               enddo
             enddo
           enddo
          endif
          call destroy_eri_3center_eigen(noft_verbose)
        else            ! Normal case (not using RI)
          call form_erimol(nbf,nstate_noft,nstate_kji,c_matrix=tmp_c_matrix,ERImol=ERImol)
        endif
      endif

    endif

  endif

  ! Deallocate tmp_c_matrix
  if(noft_complex=='yes') then
    call clean_deallocate('tmp_c_matrix',tmp_c_matrix_cmplx,noft_verbose)
  else
    call clean_deallocate('tmp_c_matrix',tmp_c_matrix,noft_verbose)
  endif

end subroutine mo_ints

!==================================================================
subroutine mo_ints_x2c(nstate_occ,nstate_kji,NO_COEF_x2c,hCORE_x2c,ERImol_x2c)
  implicit none

  integer,intent(in)        :: nstate_occ,nstate_kji
  complex(dp),intent(in)    :: NO_COEF_x2c(nstate_noft,nstate_noft)
  complex(dp),intent(inout) :: hCORE_x2c(nstate_noft,nstate_noft)
  complex(dp),intent(inout) :: ERImol_x2c(nstate_noft,nstate_kji,nstate_kji,nstate_kji)
  !====
  integer                    :: istate,jstate,pstate,qstate
  !=====

  ! Molecular hCORE including all one-body interactions
  ! T+Vext part
  hCORE_x2c=matmul(conjg(transpose(NO_COEF_x2c)),matmul(AhCORE_cmplx,NO_COEF_x2c))

  ERImol_x2c(:,:,:,:)=complex_zero
  if(has_auxil_basis) then ! RI case
    call calculate_eri_x2c(NO_COEF_x2c,nstate_noft,nstate_min=1,nstate_max=nstate_noft,mstate_min=1,mstate_max=nstate_kji)
    if(nstate_noft==nstate_kji) then
     do istate=1,nstate_noft
      do jstate=1,nstate_noft
       do pstate=1,nstate_noft
        do qstate=1,nstate_noft
         ERImol_x2c(qstate,pstate,jstate,istate)=eri_eigen_ri_x2c(qstate,jstate,pstate,istate)
        enddo
       enddo
      enddo
     enddo
    else
     do istate=1,nstate_occ
       do jstate=1,nstate_occ    ! Delta_ij contributions (Hartee and exchange terms incl. ubar and bar)
         do pstate=1,nstate_noft
           ! For 2nd, 4th, and 5th terms in Eq. 64 Rel-RDMFT paper SciPost Chem. 1, 004 (2022)
           ERImol_x2c(pstate,istate,jstate,istate)=eri_eigen_ri_x2c(pstate,jstate,istate,istate) ! <li|ji> format used for ERImol
           ERImol_x2c(pstate,jstate,jstate,istate)=eri_eigen_ri_x2c(pstate,jstate,jstate,istate) ! <lj|ji> format used for ERImol
         enddo
       enddo
       do jstate=1,nstate_occ/2  ! Pi_ij contributions (i.e. <\ubar{i} \bar{i} | \ubar{j} \bar{j}>, <\bar{i} \ubar{i} | \bar{j} \ubar{j}>,
         do pstate=1,nstate_noft ! <\ubar{i} \bar{i} | \bar{j} \ubar{j}>, and <\bar{i} \ubar{i} | \ubar{j} \bar{j}>)
           ! For 6th and 7th terms in Eq. 64 Rel-RDMFT paper SciPost Chem. 1, 004 (2022) NOTE: Focus on the types of 2-RDM elements NOT on the integrals of Eq. 64
           ERImol_x2c(pstate,istate,2*jstate-1,2*jstate  )=eri_eigen_ri_x2c(pstate,2*jstate-1,istate,2*jstate  ) ! <li|\ubar{j} \bar{j} > format used for ERImol
           ERImol_x2c(pstate,istate,2*jstate  ,2*jstate-1)=eri_eigen_ri_x2c(pstate,2*jstate  ,istate,2*jstate-1) ! <li|\bar{j} \ubar{j} > format used for ERImol
         enddo
       enddo
     enddo
    endif
    call destroy_eri_3center_eigen_x2c()
  else            ! TODO: Normal case (not using RI)
    call die("NOFT X2C needs an auxiliary basis")
  endif

end subroutine mo_ints_x2c

!==================================================================
end module m_noft
!==================================================================
