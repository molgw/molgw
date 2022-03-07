!==================================================================
! This file is part of MOLGW.
! Author: Mauricio Rodriguez-Mayorga
!
! This file contains
! - NOFT energy opt. with Resolution-of-Identity
!=========================================================================
#include "molgw.h"
subroutine noft_energy(Nelect,nstate,basis,c_matrix,AhCORE_in,AOverlap_in,Enoft,Vnn)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_hamiltonian_onebody
 use m_noft_driver
 implicit none

 integer,intent(in)         :: Nelect,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(inout)     :: AhCORE_in(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: AOverlap_in(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: Vnn
 real(dp),intent(out)       :: Enoft
!====
 integer                    :: istate,lwork,info,verbose=-1
 integer                    :: imethorb,iERItyp,NBF_occ,Nfrozen,Nbeta,Nalpha,Nvcoupled
 real(dp)                   :: ran_num
 real(dp),allocatable       :: occ(:,:),energy(:,:),occ_print(:,:)
 real(dp),allocatable       :: NO_COEF(:,:)
 real(dp),allocatable       :: tmp_mat0(:,:),tmp_mat(:,:),Work(:) 
! complex(dp)                :: ran_numC
 complex(dp),allocatable    :: NO_COEFc(:,:)
 complex(dp),allocatable    :: tmp_mat0C(:,:)
 character(len=200)         :: ofile_name
 external::mo_ints
!=====

 ! Init clock
 call start_clock(timing_noft_energy)

 ! Write hearder and set name for $name.noft file
 write(ofile_name,'(2a)') trim(output_name),'noft'
 write(stdout,'(/,a)') ' =================================================='
 write(stdout,'(a)')   ' RI-NOFT calculation'
 write(stdout,'(3a)')  ' writting NOFT results to ',trim(ofile_name),' file.'
 write(stdout,'(a)')   ' =================================================='
 write(stdout,'(/,a)') ' '

 Enoft = zero
 nbf_noft=nstate  ! Number of lin. indep. molecular orbitals
 ! These varibles will remain fixed for a while
 ! iERItyp=1 -> use notation <ij|kl>
 imethorb=1;iERItyp=1;

 ! Allocate arrays and initialize them 
 call clean_allocate('AhCORE',AhCORE,basis%nbf,basis%nbf)
 call clean_allocate('Aoverlap',Aoverlap,basis%nbf,basis%nbf)
 call clean_allocate('NO_occ',occ,basis%nbf,1)
 call clean_allocate('NO_energies',energy,basis%nbf,1)
 if(complexnoft=='yes') then
   call clean_allocate('NO_COEFc',NO_COEFc,basis%nbf,basis%nbf)
 else
   call clean_allocate('NO_COEF',NO_COEF,basis%nbf,basis%nbf)
 endif
 occ=zero; energy=zero;
  ! Save Atomic hCORE integrals and atomic overlaps
 AhCORE(:,:)=AhCORE_in(:,:)
 Aoverlap(:,:)=AOverlap_in(:,:) 
  ! Initially copy c_matrix (HF orbs) to NO_COEF
 if(complexnoft=='yes') then
   NO_COEFc(:,:)=complex_zero
!   call random_number(ran_num)
!   ran_numC=exp(im*ran_num)
   do istate=1,nbf_noft
!     call random_number(ran_num)
     !NO_COEFc(:,istate)=exp(im*ran_num)*c_matrix(:,istate,1)
     !NO_COEFc(:,istate)=ran_numC*c_matrix(:,istate,1)
     !NO_COEFc(:,istate)=im*c_matrix(:,istate,1)  ! OK
     NO_COEFc(:,istate)=c_matrix(:,istate,1)      ! OK
   enddo
 else
   NO_COEF(:,:)=zero
   do istate=1,nbf_noft
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
   if(complexnoft=='yes') then
    call clean_allocate('tmp_mat0C',tmp_mat0C,basis%nbf,basis%nbf,verbose)
    tmp_mat0C=matmul(NO_COEFc,tmp_mat)
    NO_COEFc=tmp_mat0C
    call clean_deallocate('tmp_mat0C',tmp_mat0C,verbose)
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
 Nvcoupled=ncoupled-1
 Nfrozen=(Nelect-2*Npairs)/2
 Nbeta=Nfrozen+Npairs
 Nalpha=Nbeta
 do
  NBF_occ=Nfrozen+Npairs*(Nvcoupled+1)
  if(NBF_occ<=nbf_noft) then
    exit
  else
    Nvcoupled=Nvcoupled-1
  endif
 enddo 

 ! Call module initialization and run NOFT calc.
 if(complexnoft=='yes') then
   if(restartnoft=='yes') then
     call run_noft(INOF,Ista,basis%nbf,NBF_occ,Nfrozen,Npairs,Nvcoupled,Nbeta,Nalpha,iERItyp,&
     & imethocc,imethorb,nscf_nof,iprintdmn,iprintswdmn,iprintints,ithresh_lambda,ndiis_nof,&
     & Enoft,tolE_nof,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEFc=NO_COEFc,lowmemERI=(lowmemERI=='yes'),&
     & restart=(restartnoft=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,&
     & ireadFdiag=ireadFdiag,iNOTupdateOCC=iNOTupdateOCC,iNOTupdateORB=iNOTupdateORB,Lpower=Lpower,&
     & fcidump=(fcidump=='yes'))
   else
     call run_noft(INOF,Ista,basis%nbf,NBF_occ,Nfrozen,Npairs,Nvcoupled,Nbeta,Nalpha,iERItyp,&
     & imethocc,imethorb,nscf_nof,iprintdmn,iprintswdmn,iprintints,ithresh_lambda,ndiis_nof,&
     & Enoft,tolE_nof,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEFc=NO_COEFc,lowmemERI=(lowmemERI=='yes'),&
     & Lpower=Lpower,fcidump=(fcidump=='yes'))
   endif
 else
   if(restartnoft=='yes') then
     call run_noft(INOF,Ista,basis%nbf,NBF_occ,Nfrozen,Npairs,Nvcoupled,Nbeta,Nalpha,iERItyp,&
     & imethocc,imethorb,nscf_nof,iprintdmn,iprintswdmn,iprintints,ithresh_lambda,ndiis_nof,&
     & Enoft,tolE_nof,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(lowmemERI=='yes'),&
     & restart=(restartnoft=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,&
     & ireadFdiag=ireadFdiag,iNOTupdateOCC=iNOTupdateOCC,iNOTupdateORB=iNOTupdateORB,Lpower=Lpower,&
     & fcidump=(fcidump=='yes'))
   else
     call run_noft(INOF,Ista,basis%nbf,NBF_occ,Nfrozen,Npairs,Nvcoupled,Nbeta,Nalpha,iERItyp,&
     & imethocc,imethorb,nscf_nof,iprintdmn,iprintswdmn,iprintints,ithresh_lambda,ndiis_nof,&
     & Enoft,tolE_nof,Vnn,Aoverlap,occ(:,1),mo_ints,ofile_name,NO_COEF=NO_COEF,lowmemERI=(lowmemERI=='yes'),&
     & Lpower=Lpower,fcidump=(fcidump=='yes'))
   endif
 endif
 
 ! If required print post-procesing files 
 if(complexnoft=='yes') then
   if(print_wfn_files_ ) then
     call clean_allocate('Occ_print',occ_print,nbf_noft,1,verbose)
     occ_print(1:nbf_noft,1)=occ(1:nbf_noft,1)
     ! Update c_matrix with real part of optimized NO_COEF
     do istate=1,nbf_noft 
       c_matrix(:,istate,1)=real(NO_COEFc(:,istate))
     enddo
     call print_wfn_file('NOFT_RE',basis,occ_print,c_matrix,Enoft,energy)
     ! Update c_matrix with imaginary part of optimized NO_COEF
     do istate=1,nbf_noft 
       c_matrix(:,istate,1)=aimag(NO_COEFc(:,istate))
     enddo
     call print_wfn_file('NOFT_IM',basis,occ_print,c_matrix,Enoft,energy)
     call clean_deallocate('Occ_print',occ_print,verbose)
   endif
 else
   ! Update c_matrix with optimized NO_COEF
   do istate=1,nbf_noft 
     c_matrix(:,istate,1)=NO_COEF(:,istate)
   enddo
   ! Select the post-procesing files 
   if(print_wfn_ .or. print_cube_ .or. print_wfn_files_ ) then
     call clean_allocate('Occ_print',occ_print,nbf_noft,1,verbose)
     occ_print(1:nbf_noft,1)=occ(1:nbf_noft,1)
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
 if(complexnoft=='yes') then
   call clean_deallocate('NO_COEFc',NO_COEFc)
 else
   call clean_deallocate('NO_COEF',NO_COEF)
 endif
 write(stdout,'(/,a)') ' =================================================='
 write(stdout,'(a)')   ' RI-NOFT SCF done'
 write(stdout,'(a)')   ' =================================================='
 write(stdout,'(/,a)') ' '

 ! Stop clock
 call stop_clock(timing_noft_energy)

end subroutine noft_energy

subroutine mo_ints(nbf,nbf_occ,nbf_kji,NO_COEF,hCORE,ERImol,ERImolv,NO_COEFc,hCOREc,ERIcmol,ERIcmolv)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_hamiltonian_onebody
 implicit none

 integer,intent(in)         :: nbf,nbf_occ,nbf_kji
 real(dp),optional,intent(in)    :: NO_COEF(nbf,nbf)
 real(dp),optional,intent(inout) :: hCORE(nbf,nbf)
 real(dp),optional,intent(inout) :: ERImol(nbf,nbf_kji,nbf_kji,nbf_kji)
 real(dp),optional,intent(inout) :: ERImolv(nbf*nbf_kji*nbf_kji*nbf_kji)
 complex(dp),optional,intent(in)    :: NO_COEFc(nbf,nbf)
 complex(dp),optional,intent(inout) :: hCOREc(nbf,nbf)
 complex(dp),optional,intent(inout) :: ERIcmol(nbf,nbf_kji,nbf_kji,nbf_kji)
 complex(dp),optional,intent(inout) :: ERIcmolv(nbf*nbf_kji*nbf_kji*nbf_kji)
!====
 integer                    :: istate,jstate,kstate,lstate,verbose=-1
 real(dp),allocatable       :: tmp_hcore(:,:)
 real(dp),allocatable       :: tmp_c_matrix(:,:,:)
 complex(dp),allocatable    :: tmp_hcoreC(:,:)
 complex(dp),allocatable    :: tmp_c_matrixC(:,:,:)
!=====

 if(complexnoft=='yes') then

   ! hCORE part
   call clean_allocate('tmp_hcore',tmp_hcoreC,nbf,nbf,verbose)
   hCOREC(:,:)=complex_zero; tmp_hcoreC(:,:)=complex_zero;
   tmp_hcoreC=matmul(AhCORE,NO_COEFc)
   hCOREC=matmul(conjg(transpose(NO_COEFc)),tmp_hcoreC)
   call clean_deallocate('tmp_hcore',tmp_hcoreC,verbose)

   ! ERI terms
   if(present(ERIcmol)) then
     ERIcmol(:,:,:,:)=complex_zero
     call clean_allocate('tmp_c_matrix',tmp_c_matrixC,nbf,nbf_noft,1,verbose)
     do istate=1,nbf_noft
      tmp_c_matrixC(:,istate,1)=NO_COEFc(:,istate)
     enddo
     if(noft_ri) then ! RI case
       call calculate_eri_3center_eigenC(tmp_c_matrixC,1,nbf_noft,1,nbf_kji,verbose=verbose)
       do istate=1,nbf_occ
         do jstate=1,nbf_occ
           do kstate=1,nbf_occ
             do lstate=1,nbf_noft
               ERIcmol(lstate,kstate,jstate,istate)=eri_eigen_riC(lstate,jstate,1,kstate,istate,1) ! <lk|ji> format used for ERImol
             enddo
           enddo
         enddo
       enddo
       call destroy_eri_3center_eigenC(verbose)
     else            ! Normal case 
      ! TODO
     endif
     call clean_deallocate('tmp_c_matrix',tmp_c_matrixC,verbose)
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
     call clean_allocate('tmp_c_matrix',tmp_c_matrix,nbf,nbf_noft,1,verbose)
     do istate=1,nbf_noft
      tmp_c_matrix(:,istate,1)=NO_COEF(:,istate)
     enddo
     if(noft_ri) then ! RI case
       call calculate_eri_3center_eigen(tmp_c_matrix,1,nbf_noft,1,nbf_kji,verbose=verbose)
       do istate=1,nbf_occ
         do jstate=1,nbf_occ
           do kstate=1,nbf_occ
             do lstate=1,nbf_noft
               ERImol(lstate,kstate,jstate,istate)=eri_eigen_ri(lstate,jstate,1,kstate,istate,1) ! <lk|ji> format used for ERImol
             enddo
           enddo
         enddo
       enddo
       call destroy_eri_3center_eigen(verbose)
     else            ! Normal case 
      ! TODO
     endif
     call clean_deallocate('tmp_c_matrix',tmp_c_matrix,verbose)
   endif

 endif

end subroutine mo_ints

!==================================================================
