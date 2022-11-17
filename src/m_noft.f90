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


 logical,parameter,private    :: noft_verbose=.FALSE.
 logical                      :: rs_noft=.FALSE.,noft_edft=.FALSE.
 integer,private              :: nstate_noft,nstate_frozen 
 real(dp)                     :: ExcDFT,E_t_vext
 real(dp),allocatable,private :: AhCORE(:,:),T_Vext(:)            ! hCORE matrix (T+Ven) in AO basis
 type(basis_set),pointer      :: basis_pointer


contains

!=========================================================================
subroutine noft_energy(basis,c_matrix,occupation,hkin,hnuc,Aoverlap,Enoft,Vnn)
 implicit none

 type(basis_set),intent(in),target :: basis
 real(dp),intent(inout)    :: c_matrix(:,:,:)
 real(dp),intent(in)       :: Aoverlap(:,:)
 real(dp),intent(in)       :: hkin(:,:),hnuc(:,:)
 real(dp),intent(inout)    :: occupation(:,:)
 real(dp),intent(in)       :: Vnn
 real(dp),intent(out)      :: Enoft
 character(len=100)        :: msgw
!====
 integer                   :: istate,lwork,info
 integer                   :: nelectrons
 integer                   :: imethorb,imethocc,iERItyp,nstate_occ,nstate_beta,nstate_alpha,nstate_coupled
 integer                   :: iNOTupdateOCC,iNOTupdateORB,iprintdmn,iprintswdmn,iprintints,ireadOCC,ireadCOEF 
 integer                   :: ireadFdiag,ireadGAMMAs,ista,inof 
 real(dp)                  :: ran_num
 real(dp),allocatable      :: occ(:,:),energy(:,:),occ_print(:,:)
 real(dp),allocatable      :: NO_COEF(:,:)
 real(dp),allocatable      :: tmp_mat0(:,:),tmp_mat(:,:),Work(:) 
 complex(dp),allocatable   :: NO_COEF_cmplx(:,:)
 complex(dp),allocatable   :: tmp_mat_cmplx(:,:)
 character(len=200)        :: ofile_name
!=====

 basis_pointer=>basis

 ! Init clock
 call start_clock(timing_noft_energy)

 ! Write hearder and set name for $name.noft file
 write(ofile_name,'(2a)') trim(output_name),'noft'
 write(stdout,'(/,a)') ' =================================================='
 write(stdout,'(a)')   ' NOFT calculation'
 write(stdout,'(3a)')  ' writting NOFT results to ',trim(ofile_name),' file.'
 write(stdout,'(a)')   ' =================================================='
 write(stdout,'(/,a)') ' '

 !
 ! Setup the grids for the quadrature of DFT potential/energy
 if( calc_type%is_dft .and. noft_dft=='yes' ) then
   rs_noft=.true.
   if( .not.calc_type%need_exchange_lr ) then
     write(msgw,'(a)') 'LR exchange is needed for rs-NOFT. The calculation will proceed without range-sep.'
     call issue_warning(msgw)
     rs_noft=.false.
   endif
   call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
 endif

 Enoft = zero; occupation = zero; ExcDFT = zero;
 nstate_noft = SIZE(c_matrix,DIM=2) ! Number of lin. indep. molecular orbitals

 ! These varibles will remain fixed for a while
 ! iERItyp=1 -> use notation <ij|kl>
 imethorb=1;iERItyp=1;imethocc=1;iNOTupdateOCC=0;iNOTupdateORB=0;
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
 if(noft_sta=='yes') ista=1

 select case(capitalize(noft_functional))
 case('GNOF')
   inof=8
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
 case default
   call die('noft_energy: NOFT functional not recognized. Check your input')
 end select

 ! Allocate arrays and initialize them 
 call clean_allocate('AhCORE',AhCORE,basis%nbf,basis%nbf)
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
 AhCORE(:,:) = hkin(:,:) + hnuc(:,:)

 ! Initially copy c_matrix (HF orbs) to NO_COEF
 if(noft_complex=='yes') then
   NO_COEF_cmplx(:,:)=complex_zero
!   call random_number(ran_num) ! For complex orbs with same random phase (i.e. to have real and imaginary orbs)
   ran_num=pi/four
   do istate=1,nstate_noft
!     call random_number(ran_num) ! For complex orbs each having its own random phase (i.e. to have real and imaginary orbs)
     NO_COEF_cmplx(:,istate)=exp(im*ran_num)*c_matrix(:,istate,1)
!     NO_COEF_cmplx(:,istate)=im*c_matrix(:,istate,1)
!     NO_COEF_cmplx(:,istate)=c_matrix(:,istate,1)
   enddo
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
   allocate(Work(1))
   tmp_mat0=matmul(AhCORE,NO_COEF)
   tmp_mat=matmul(transpose(NO_COEF),tmp_mat0)
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
    call clean_allocate('tmp_mat_cmplx',tmp_mat_cmplx,basis%nbf,basis%nbf,noft_verbose)
    tmp_mat_cmplx=matmul(NO_COEF_cmplx,tmp_mat)
    NO_COEF_cmplx=tmp_mat_cmplx
    call clean_deallocate('tmp_mat_cmplx',tmp_mat_cmplx,noft_verbose)
   else
     tmp_mat0=matmul(NO_COEF,tmp_mat)
     NO_COEF=tmp_mat0
   endif
   write(stdout,'(/,a,/)') ' Approximate Hamiltonian Hcore used as GUESS in NOFT calc.'
   call clean_deallocate('tmp_mat0',tmp_mat0,noft_verbose)
   call clean_deallocate('tmp_mat',tmp_mat,noft_verbose)
   deallocate(Work)
 endif 

 ! Not ready for open-shell calcs. (TODO)
 nelectrons=NINT(electrons)
 nstate_coupled=noft_ncoupled-1
 nstate_frozen=(nelectrons-2*noft_npairs)/2
 nstate_beta=nstate_frozen+noft_npairs
 nstate_alpha=nstate_beta
 do
  nstate_occ=nstate_frozen+noft_npairs*(nstate_coupled+1)
  if(nstate_occ<=nstate_noft) then
    exit
  else
    nstate_coupled=nstate_coupled-1
  endif
 enddo 

 ! Call module initialization and run NOFT calc.
 if(noft_complex=='yes') then
   if(noft_restart=='yes') then
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
     & iERItyp,imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
     & Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF_cmplx=NO_COEF_cmplx,lowmemERI=(noft_lowmemERI=='yes'),&
     & restart=(noft_restart=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,&
     & ireadFdiag=ireadFdiag,iNOTupdateOCC=iNOTupdateOCC,iNOTupdateORB=iNOTupdateORB,Lpower=noft_Lpower,&
     & fcidump=(noft_fcidump=='yes'),range_sep=rs_noft)
   else
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
     & iERItyp,imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
     & Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF_cmplx=NO_COEF_cmplx,lowmemERI=(noft_lowmemERI=='yes'),&
     & Lpower=noft_Lpower,fcidump=(noft_fcidump=='yes'),range_sep=rs_noft)
   endif
 else
   if(noft_restart=='yes') then
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
     & iERItyp,imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
     & Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),&
     & restart=(noft_restart=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,&
     & ireadFdiag=ireadFdiag,iNOTupdateOCC=iNOTupdateOCC,iNOTupdateORB=iNOTupdateORB,Lpower=noft_Lpower,&
     & fcidump=(noft_fcidump=='yes'),range_sep=rs_noft)
   else
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
     & iERItyp,imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
     & Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),&
     & Lpower=noft_Lpower,fcidump=(noft_fcidump=='yes'),range_sep=rs_noft)
   endif
 endif
 
 if( rs_noft ) then ! Compute total Energy for range-sep NOFT switching off the hamiltonian_xc (we only need hCORE=T+Vext).
   noft_edft=.true. ! this is why we restart but will not update orbs or occs.
   call clean_allocate('T_Vext',T_Vext,basis%nbf,noft_verbose)
   write(ofile_name,'(a)') 'tmp_dft_noft'
   call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
   & iERItyp,imethocc,imethorb,noft_nscf,0,0,0,noft_ithresh_lambda,noft_ndiis,Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),&
   & mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),restart=.true.,ireadGAMMAS=1,ireadOCC=1,&
   & ireadCOEF=1,ireadFdiag=1,iNOTupdateOCC=1,iNOTupdateORB=1,Lpower=noft_Lpower,fcidump=(noft_fcidump=='yes'),range_sep=rs_noft)
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

 ! Deallocate arrays and print the normal termination 
 call clean_deallocate('AhCORE',AhCORE)
 call clean_deallocate('NO_occ',occ)
 call clean_deallocate('NO_energies',energy)
 if(noft_complex=='yes') then
   call clean_deallocate('NO_COEF_cmplx',NO_COEF_cmplx)
 else
   call clean_deallocate('NO_COEF',NO_COEF)
 endif
 write(stdout,'(/,a)') ' =================================================='
 write(stdout,'(a)')   ' NOFT SCF loop ends here'
 write(stdout,'(a)')   ' =================================================='
 write(stdout,'(/,a)') ' '

 ! Stop clock
 call stop_clock(timing_noft_energy)

end subroutine noft_energy


!==================================================================
subroutine mo_ints(nbf,nstate_occ,nstate_kji,Occ,NO_COEF,hCORE,ERImol,ERImolv,ERImolH,NO_COEF_cmplx,hCORE_cmplx,&
& ERImol_cmplx,ERImolv_cmplx)
 implicit none

 integer,intent(in)              :: nbf,nstate_occ,nstate_kji
 real(dp),intent(in)             :: Occ(nstate_occ)
 real(dp),optional,intent(in)    :: NO_COEF(nbf,nbf)
 real(dp),optional,intent(inout) :: hCORE(nbf,nbf)
 real(dp),optional,intent(inout) :: ERImol(nbf,nstate_kji,nstate_kji,nstate_kji)
 real(dp),optional,intent(inout) :: ERImolH(nbf,nstate_kji,nstate_kji)
 real(dp),optional,intent(inout) :: ERImolv(nbf*nstate_kji*nstate_kji*nstate_kji)
 complex(dp),optional,intent(in)    :: NO_COEF_cmplx(nbf,nbf)
 complex(dp),optional,intent(inout) :: hCORE_cmplx(nbf,nbf)
 complex(dp),optional,intent(inout) :: ERImol_cmplx(nbf,nstate_kji,nstate_kji,nstate_kji)
 complex(dp),optional,intent(inout) :: ERImolv_cmplx(nbf*nstate_kji*nstate_kji*nstate_kji)
!====
 logical                    :: long_range=.true.
 integer                    :: istate,jstate,kstate,lstate
 character(len=100)         :: msgw
 real(dp)                   :: full_eri
 real(dp),allocatable       :: occupation(:,:)!,hamiltonian_hartree(:,:)
 real(dp),allocatable       :: tmp_c_matrix(:,:,:),hamiltonian_xc(:,:,:)!,p_matrix(:,:,:)
 complex(dp),allocatable    :: tmp_c_matrix_cmplex(:,:,:)
!=====

! Comment: Despite the arrays are of size nbf x nbf, we use nstate_noft = num. lin. indep. states in the ERI  transformation. 
! Doing this, we save some time in the loops because nstate_noft <= nbf

 if( rs_noft ) then

   if(noft_complex=='yes') then

     !TODO
     write(msgw,'(a)') 'LR complex exchange is needed for rs-NOFT, but not coded.'
     call issue_warning(msgw)

   else

     ! Build 3D array for c_matrix and init hCORE
     call clean_allocate('tmp_c_matrix',tmp_c_matrix,nbf,nstate_noft,1,noft_verbose)
     call clean_allocate('occupation',occupation,nbf,1,noft_verbose)
     occupation(:,:)=zero; occupation(:nstate_occ,1)=two*Occ(:nstate_occ);
     tmp_c_matrix(:,:,:)=zero; hCORE(:,:)=zero; 
     do istate=1,nstate_noft
      tmp_c_matrix(:,istate,1)=NO_COEF(:,istate)
     enddo

     if(.not.noft_edft) then
       ! Prepare the DFT contribution (takes part only during orb. optimization and is switched off for final energy calculation)
       call clean_allocate('hamiltonian_xc',hamiltonian_xc,nbf,nbf,1,noft_verbose)
       hamiltonian_xc(:,:,:)=zero;
       call dft_exc_vxc_batch(BATCH_SIZE,basis_pointer,occupation,tmp_c_matrix,hamiltonian_xc,ExcDFT)
       hCORE=matmul(transpose(NO_COEF(:,:)),matmul(hamiltonian_xc(:,:,1),NO_COEF(:,:)))
       call clean_deallocate('hamiltonian_xc',hamiltonian_xc,noft_verbose)
       
       ! MRM: Actually, we don't need the Vhartree in AO basis... But, we could use it in the future.
       ! Prepare the Vhartree contribution
         !call clean_allocate('density matrix P',p_matrix,nbf,nbf,1,noft_verbose)
         !call clean_allocate('hamiltonian_hartree',hamiltonian_hartree,nbf,nbf,noft_verbose)
         !p_matrix(:,:,:)=zero; hamiltonian_hartree(:,:)=zero;
         !call setup_density_matrix(tmp_c_matrix,occupation,p_matrix)
         !call calculate_hartree(basis_pointer,p_matrix,hamiltonian_hartree,Ehartree)
         !hCORE=hCORE+matmul(transpose(NO_COEF(:,:)),matmul(hamiltonian_hartree(:,:),NO_COEF(:,:)))
         !call clean_deallocate('hamiltonian_hartree',hamiltonian_hartree,noft_verbose)
         !call clean_deallocate('density matrix P',p_matrix,noft_verbose)
     endif    
 
     call clean_deallocate('occupation',occupation,noft_verbose)

     ! hCORE = T + Vext (and add the DFT if needed) 
     hCORE=hCORE+matmul(transpose(NO_COEF(:,:)),matmul(AhCORE(:,:),NO_COEF(:,:)))
     if(noft_edft) then
       do istate=1,nstate_noft
         T_Vext(istate)=hCORE(istate,istate)
       enddo
     endif

     ! ERI terms
     if(present(ERImol) .and. present(ERImolH)) then
       ERImol(:,:,:,:)=zero; ERImolH(:,:,:)=zero;
       if(has_auxil_basis) then ! RI case
         call calculate_eri_3center_eigen(tmp_c_matrix,1,nstate_noft,1,nstate_kji,verbose=noft_verbose,long_range=long_range)
         do istate=1,nstate_occ
           do jstate=1,nstate_occ
             do kstate=1,nstate_occ
               do lstate=1,nstate_noft
                 full_eri=eri_eigen_ri(lstate,jstate,1,kstate,istate,1)
                 ERImol(lstate,kstate,jstate,istate)=alpha_hybrid*full_eri & ! <lk| [alpha+beta*erf(gamma r12)]/r12 |ji> format used for ERImol
                 & +beta_hybrid*eri_eigen_ri_lr(lstate,jstate,1,kstate,istate,1) 
                 if(kstate==istate) then ! Hartree: <li|ji>^Hartree = <li| 1/r12 |ji> - <li| [alpha+beta*erf(gamma r12)]/r12 |ji>
                   ERImolH(lstate,istate,jstate)=full_eri-ERImol(lstate,istate,jstate,istate) ! <li|ji> format used for ERImol
                 endif
               enddo
             enddo
           enddo
         enddo
         call destroy_eri_3center_eigen(verbose=noft_verbose,long_range=long_range)
       else            ! Normal case (not using RI)
         !TODO
         write(msgw,'(a)') 'LR exchange without RI is needed for rs-NOFT, but not coded.'
         call issue_warning(msgw)
       endif
     endif

     call clean_deallocate('tmp_c_matrix',tmp_c_matrix,noft_verbose)

   endif

 else

   if(noft_complex=='yes') then

     ! hCORE part
     hCORE_cmplx(:,:)=complex_zero
     hCORE_cmplx=matmul(conjg(transpose(NO_COEF_cmplx)),matmul(AhCORE,NO_COEF_cmplx))

     ! ERI terms
     if(present(ERImol_cmplx)) then
       ERImol_cmplx(:,:,:,:)=complex_zero
       call clean_allocate('tmp_c_matrix',tmp_c_matrix_cmplex,nbf,nstate_noft,1,noft_verbose)
       do istate=1,nstate_noft
        tmp_c_matrix_cmplex(:,istate,1)=NO_COEF_cmplx(:,istate)
       enddo
       if(has_auxil_basis) then ! RI case
         call calculate_eri_3center_eigen_cmplx(tmp_c_matrix_cmplex,1,nstate_noft,1,nstate_kji,verbose=noft_verbose)
         do istate=1,nstate_occ
           do jstate=1,nstate_occ
             do kstate=1,nstate_occ
               do lstate=1,nstate_noft
                 ERImol_cmplx(lstate,kstate,jstate,istate)=eri_eigen_ri_cmplx(lstate,jstate,1,kstate,istate,1) ! <lk|ji> format used for ERImol
               enddo
             enddo
           enddo
         enddo
         call destroy_eri_3center_eigen_cmplx(noft_verbose)
       else            ! Normal case (not using RI) 
         call form_erimol(nbf,nstate_noft,nstate_kji,c_matrix_cmplx=tmp_c_matrix_cmplex,ERImol_cmplx=ERImol_cmplx)
       endif
       call clean_deallocate('tmp_c_matrix',tmp_c_matrix_cmplex,noft_verbose)
     endif

   else

     ! hCORE part
     hCORE(:,:)=zero
     hCORE=matmul(transpose(NO_COEF),matmul(AhCORE,NO_COEF))

     ! ERI terms
     if(present(ERImol)) then
       ERImol(:,:,:,:)=zero
       call clean_allocate('tmp_c_matrix',tmp_c_matrix,nbf,nstate_noft,1,noft_verbose)
       do istate=1,nstate_noft
        tmp_c_matrix(:,istate,1)=NO_COEF(:,istate)
       enddo
       if(has_auxil_basis) then ! RI case
         call calculate_eri_3center_eigen(tmp_c_matrix,1,nstate_noft,1,nstate_kji,verbose=noft_verbose)
         do istate=1,nstate_occ
           do jstate=1,nstate_occ
             do kstate=1,nstate_occ
               do lstate=1,nstate_noft
                 ERImol(lstate,kstate,jstate,istate)=eri_eigen_ri(lstate,jstate,1,kstate,istate,1) ! <lk|ji> format used for ERImol
               enddo
             enddo
           enddo
         enddo
         call destroy_eri_3center_eigen(noft_verbose)
       else            ! Normal case (not using RI)
         call form_erimol(nbf,nstate_noft,nstate_kji,c_matrix=tmp_c_matrix,ERImol=ERImol)
       endif
       call clean_deallocate('tmp_c_matrix',tmp_c_matrix,noft_verbose)
     endif

   endif

 endif

end subroutine mo_ints


!==================================================================
end module m_noft
!==================================================================
