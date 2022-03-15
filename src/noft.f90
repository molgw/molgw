!==================================================================
! This file is part of MOLGW.
! Author: Mauricio Rodriguez-Mayorga
!
! This file contains
! - NOFT energy opt. with Resolution-of-Identity
!=========================================================================
#include "molgw.h"
subroutine noft_energy(nelectrons,nstate,basis,c_matrix,AhCORE_in,AOverlap_in,Enoft,Vnn)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_hamiltonian_onebody
 use m_noft_driver
 implicit none

 integer,intent(in)         :: nelectrons,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(inout)     :: AhCORE_in(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: AOverlap_in(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: Vnn
 real(dp),intent(out)       :: Enoft
!====
 integer                    :: istate,lwork,info,verbose=-1
 integer                    :: imethorb,imethocc,iERItyp,nstate_occ,nstate_frozen,nstate_beta,nstate_alpha,nstate_coupled
 integer                    :: iNOTupdateOCC,iNOTupdateORB,iprintdmn,iprintswdmn,iprintints,ireadOCC,ireadCOEF 
 integer                    :: ireadFdiag,ireadGAMMAs,ista,inof 
 real(dp)                   :: ran_num
 real(dp),allocatable       :: occ(:,:),energy(:,:),occ_print(:,:)
 real(dp),allocatable       :: NO_COEF(:,:)
 real(dp),allocatable       :: tmp_mat0(:,:),tmp_mat(:,:),Work(:) 
 complex(dp)                :: ran_num_cmplx
 complex(dp),allocatable    :: NO_COEF_cmplx(:,:)
 complex(dp),allocatable    :: tmp_mat_cmplx(:,:)
 character(len=200)         :: ofile_name
 external::mo_ints
!=====

 ! Init clock
 call start_clock(timing_noft_energy)

 ! Write hearder and set name for $name.noft file
 write(ofile_name,'(2a)') trim(output_name),'noft'
 write(stdout,'(/,a)') ' =================================================='
 write(stdout,'(a)')   ' NOFT calculation'
 write(stdout,'(3a)')  ' writting NOFT results to ',trim(ofile_name),' file.'
 write(stdout,'(a)')   ' =================================================='
 write(stdout,'(/,a)') ' '

 Enoft = zero
 nstate_nof=nstate  ! Number of lin. indep. molecular orbitals
 ! These varibles will remain fixed for a while
 ! iERItyp=1 -> use notation <ij|kl>
 imethorb=1;iERItyp=1;imethocc=1;iNOTupdateOCC=0;iNOTupdateORB=0;
 iprintdmn=0;iprintswdmn=0;iprintints=0;ireadOCC=0;ireadCOEF=0;
 ireadFdiag=0;ireadGAMMAs=0;ista=1;inof=7;
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
 if(noft_functional=='PNOF7') then
   inof=7
 elseif(noft_functional=='PNOF5') then
   inof=5
 elseif(noft_functional=='MULLER') then
   inof=-1
 elseif(noft_functional=='POWER') then
   inof=-2
 else ! HF
   inof=0 
 endif
 ! Allocate arrays and initialize them 
 call clean_allocate('AhCORE',AhCORE,basis%nbf,basis%nbf)
 call clean_allocate('Aoverlap',Aoverlap,basis%nbf,basis%nbf)
 call clean_allocate('NO_occ',occ,basis%nbf,1)
 call clean_allocate('NO_energies',energy,basis%nbf,1)
 if(noft_complex=='yes') then
   call clean_allocate('NO_COEF_cmplx',NO_COEF_cmplx,basis%nbf,basis%nbf)
 else
   call clean_allocate('NO_COEF',NO_COEF,basis%nbf,basis%nbf)
 endif
 occ=zero; energy=zero;
  ! Save Atomic hCORE integrals and atomic overlaps
 AhCORE(:,:)=AhCORE_in(:,:)
 Aoverlap(:,:)=AOverlap_in(:,:) 
  ! Initially copy c_matrix (HF orbs) to NO_COEF
 if(noft_complex=='yes') then
   NO_COEF_cmplx(:,:)=complex_zero
   call random_number(ran_num)
   ran_num_cmplx=exp(im*ran_num)
   do istate=1,nstate_nof
!     call random_number(ran_num)
!     NO_COEF_cmplx(:,istate)=exp(im*ran_num)*c_matrix(:,istate,1)
!     NO_COEF_cmplx(:,istate)=ran_num_cmplx*c_matrix(:,istate,1)
!     NO_COEF_cmplx(:,istate)=im*c_matrix(:,istate,1)
     NO_COEF_cmplx(:,istate)=c_matrix(:,istate,1)     ! This is a faked complex orbs usage because orbs are still real. TODO
   enddo
 else
   NO_COEF(:,:)=zero
   do istate=1,nstate_nof
     NO_COEF(:,istate)=c_matrix(:,istate,1)
   enddo
 endif
  ! Replace NO_COEF by Hcore orbs for the initial GUESS?
 if(TRIM(init_hamiltonian)=='CORE') then
   call clean_allocate('tmp_mat0',tmp_mat0,basis%nbf,basis%nbf,verbose)
   call clean_allocate('tmp_mat',tmp_mat,basis%nbf,basis%nbf,verbose)
   allocate(Work(1))
   tmp_mat0=matmul(AhCORE_in,NO_COEF)
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
    call clean_allocate('tmp_mat_cmplx',tmp_mat_cmplx,basis%nbf,basis%nbf,verbose)
    tmp_mat_cmplx=matmul(NO_COEF_cmplx,tmp_mat)
    NO_COEF_cmplx=tmp_mat_cmplx
    call clean_deallocate('tmp_mat_cmplx',tmp_mat_cmplx,verbose)
   else
     tmp_mat0=matmul(NO_COEF,tmp_mat)
     NO_COEF=tmp_mat0
   endif
   write(stdout,'(/,a,/)') ' Approximate Hamiltonian Hcore used as GUESS in NOFT calc.'
   call clean_deallocate('tmp_mat0',tmp_mat0,verbose)
   call clean_deallocate('tmp_mat',tmp_mat,verbose)
   deallocate(Work)
 endif 

 ! Not ready for open-shell calcs. (TODO)
 nstate_coupled=noft_ncoupled-1
 nstate_frozen=(nelectrons-2*noft_npairs)/2
 nstate_beta=nstate_frozen+noft_npairs
 nstate_alpha=nstate_beta
 do
  nstate_occ=nstate_frozen+noft_npairs*(nstate_coupled+1)
  if(nstate_occ<=nstate_nof) then
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
     & fcidump=(noft_fcidump=='yes'))
   else
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
     & iERItyp,imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
     & Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF_cmplx=NO_COEF_cmplx,lowmemERI=(noft_lowmemERI=='yes'),&
     & Lpower=noft_Lpower,fcidump=(noft_fcidump=='yes'))
   endif
 else
   if(noft_restart=='yes') then
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
     & iERItyp,imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
     & Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),&
     & restart=(noft_restart=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,&
     & ireadFdiag=ireadFdiag,iNOTupdateOCC=iNOTupdateOCC,iNOTupdateORB=iNOTupdateORB,Lpower=noft_Lpower,&
     & fcidump=(noft_fcidump=='yes'))
   else
     call run_noft(inof,ista,basis%nbf,nstate_occ,nstate_frozen,noft_npairs,nstate_coupled,nstate_beta,nstate_alpha,&
     & iERItyp,imethocc,imethorb,noft_nscf,iprintdmn,iprintswdmn,iprintints,noft_ithresh_lambda,noft_ndiis,&
     & Enoft,noft_tolE,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(noft_lowmemERI=='yes'),&
     & Lpower=noft_Lpower,fcidump=(noft_fcidump=='yes'))
   endif
 endif
 
 ! If required print post-procesing files 
 if(noft_complex=='yes') then
   if(print_wfn_files_ ) then
     call clean_allocate('Occ_print',occ_print,nstate_nof,1,verbose)
     occ_print(1:nstate_nof,1)=occ(1:nstate_nof,1)
     ! Update c_matrix with real part of optimized NO_COEF
     do istate=1,nstate_nof 
       c_matrix(:,istate,1)=real(NO_COEF_cmplx(:,istate))
     enddo
     call print_wfn_file('NOFT_RE',basis,occ_print,c_matrix,Enoft,energy)
     ! Update c_matrix with imaginary part of optimized NO_COEF
     do istate=1,nstate_nof 
       c_matrix(:,istate,1)=aimag(NO_COEF_cmplx(:,istate))
     enddo
     call print_wfn_file('NOFT_IM',basis,occ_print,c_matrix,Enoft,energy)
     call clean_deallocate('Occ_print',occ_print,verbose)
   endif
 else
   ! Update c_matrix with optimized NO_COEF
   do istate=1,nstate_nof 
     c_matrix(:,istate,1)=NO_COEF(:,istate)
   enddo
   ! Select the post-procesing files 
   if(print_wfn_ .or. print_cube_ .or. print_wfn_files_ ) then
     call clean_allocate('Occ_print',occ_print,nstate_nof,1,verbose)
     occ_print(1:nstate_nof,1)=occ(1:nstate_nof,1)
     if( print_wfn_ )  call plot_wfn(basis,c_matrix)
     if( print_wfn_ )  call plot_rho('NOFT',basis,occ_print,c_matrix)
     if( print_cube_ ) call plot_cube_wfn('NOFT',basis,occ_print,c_matrix)
     if( print_wfn_files_ ) call print_wfn_file('NOFT',basis,occ_print,c_matrix,Enoft,energy)
     call clean_deallocate('Occ_print',occ_print,verbose)
   endif
 endif

 ! Deallocate arrays and print the normal termination 
 call clean_deallocate('AhCORE',AhCORE)
 call clean_deallocate('Aoverlap',Aoverlap)
 call clean_deallocate('NO_occ',occ)
 call clean_deallocate('NO_energies',energy)
 if(noft_complex=='yes') then
   call clean_deallocate('NO_COEF_cmplx',NO_COEF_cmplx)
 else
   call clean_deallocate('NO_COEF',NO_COEF)
 endif
 write(stdout,'(/,a)') ' =================================================='
 write(stdout,'(a)')   ' NOFT SCF done'
 write(stdout,'(a)')   ' =================================================='
 write(stdout,'(/,a)') ' '

 ! Stop clock
 call stop_clock(timing_noft_energy)

end subroutine noft_energy

subroutine mo_ints(nbf,nstate_occ,nstate_kji,NO_COEF,hCORE,ERImol,ERImolv,NO_COEF_cmplx,hCORE_cmplx,ERImol_cmplx,ERImolv_cmplx)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_hamiltonian_onebody
 implicit none

 integer,intent(in)         :: nbf,nstate_occ,nstate_kji
 real(dp),optional,intent(in)    :: NO_COEF(nbf,nbf)
 real(dp),optional,intent(inout) :: hCORE(nbf,nbf)
 real(dp),optional,intent(inout) :: ERImol(nbf,nstate_kji,nstate_kji,nstate_kji)
 real(dp),optional,intent(inout) :: ERImolv(nbf*nstate_kji*nstate_kji*nstate_kji)
 complex(dp),optional,intent(in)    :: NO_COEF_cmplx(nbf,nbf)
 complex(dp),optional,intent(inout) :: hCORE_cmplx(nbf,nbf)
 complex(dp),optional,intent(inout) :: ERImol_cmplx(nbf,nstate_kji,nstate_kji,nstate_kji)
 complex(dp),optional,intent(inout) :: ERImolv_cmplx(nbf*nstate_kji*nstate_kji*nstate_kji)
!====
 integer                    :: istate,jstate,kstate,lstate,verbose=-1
 real(dp),allocatable       :: tmp_hcore(:,:)
 real(dp),allocatable       :: tmp_c_matrix(:,:,:)
 complex(dp),allocatable    :: tmp_hcore_cmplx(:,:)
 complex(dp),allocatable    :: tmp_c_matrix_cmplex(:,:,:)
!=====

 if(noft_complex=='yes') then

   ! hCORE part
   call clean_allocate('tmp_hcore',tmp_hcore_cmplx,nbf,nbf,verbose)
   hCORE_cmplx(:,:)=complex_zero; tmp_hcore_cmplx(:,:)=complex_zero;
   tmp_hcore_cmplx=matmul(AhCORE,NO_COEF_cmplx)
   hCORE_cmplx=matmul(conjg(transpose(NO_COEF_cmplx)),tmp_hcore_cmplx)
   call clean_deallocate('tmp_hcore',tmp_hcore_cmplx,verbose)

   ! ERI terms
   if(present(ERImol_cmplx)) then
     ERImol_cmplx(:,:,:,:)=complex_zero
     call clean_allocate('tmp_c_matrix',tmp_c_matrix_cmplex,nbf,nstate_nof,1,verbose)
     do istate=1,nstate_nof
      tmp_c_matrix_cmplex(:,istate,1)=NO_COEF_cmplx(:,istate)
     enddo
     if(noft_ri) then ! RI case
       call calculate_eri_3center_eigen_cmplx(tmp_c_matrix_cmplex,1,nstate_nof,1,nstate_kji,verbose=verbose)
       do istate=1,nstate_occ
         do jstate=1,nstate_occ
           do kstate=1,nstate_occ
             do lstate=1,nstate_nof
               ERImol_cmplx(lstate,kstate,jstate,istate)=eri_eigen_ri_cmplx(lstate,jstate,1,kstate,istate,1) ! <lk|ji> format used for ERImol
             enddo
           enddo
         enddo
       enddo
       call destroy_eri_3center_eigen_cmplx(verbose)
     else            ! Normal case (not using RI) 
       call form_erimol(nbf,nstate_nof,nstate_kji,c_matrix_cmplx=tmp_c_matrix_cmplex,ERImol_cmplx=ERImol_cmplx)
     endif
     call clean_deallocate('tmp_c_matrix',tmp_c_matrix_cmplex,verbose)
   endif

 else

   ! hCORE part
   call clean_allocate('tmp_hcore',tmp_hcore,nbf,nbf,verbose)
   hCORE(:,:)=zero; tmp_hcore(:,:)=zero;
   tmp_hcore=matmul(AhCORE,NO_COEF)
   hCORE=matmul(transpose(NO_COEF),tmp_hcore)
   call clean_deallocate('tmp_hcore',tmp_hcore,verbose)

   ! ERI terms
   if(present(ERImol)) then
     ERImol(:,:,:,:)=zero
     call clean_allocate('tmp_c_matrix',tmp_c_matrix,nbf,nstate_nof,1,verbose)
     do istate=1,nstate_nof
      tmp_c_matrix(:,istate,1)=NO_COEF(:,istate)
     enddo
     if(noft_ri) then ! RI case
       call calculate_eri_3center_eigen(tmp_c_matrix,1,nstate_nof,1,nstate_kji,verbose=verbose)
       do istate=1,nstate_occ
         do jstate=1,nstate_occ
           do kstate=1,nstate_occ
             do lstate=1,nstate_nof
               ERImol(lstate,kstate,jstate,istate)=eri_eigen_ri(lstate,jstate,1,kstate,istate,1) ! <lk|ji> format used for ERImol
             enddo
           enddo
         enddo
       enddo
       call destroy_eri_3center_eigen(verbose)
     else            ! Normal case (not using RI)
       call form_erimol(nbf,nstate_nof,nstate_kji,c_matrix=tmp_c_matrix,ERImol=ERImol)
     endif
     call clean_deallocate('tmp_c_matrix',tmp_c_matrix,verbose)
   endif

 endif

end subroutine mo_ints

!==================================================================
