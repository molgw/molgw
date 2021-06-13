!==================================================================
! This file is part of MOLGW.
! Author: Mauricio Rodriguez-Mayorga
!
! This file contains
! - NOFT energy opt. with Resolution-of-Identity
!=========================================================================
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
 integer                    :: istate
 real(dp),allocatable       :: NO_COEF(:,:),occ(:,:),energy(:,:),occ_print(:,:) 
 integer::imethorb,iERItyp,NBF_occ,Nfrozen,Nbeta,Nalpha,Nvcoupled
 character(len=200)         :: ofile_name
 external::mo_ints
!=====

 ! Init clock
 call start_clock(timing_noft_energy)

 write(stdout,'(/,a)') ' RI-NOFT calculation'

 ofile_name='molgw.noft'
 Enoft = 0.0_dp
 nbf_noft=nstate  ! Number of lin. indep. molecular orbitals
 ! These can be fixed for a while... 
 !  iERItyp=1 -> use notation <ij|kl>
 imethorb=1;iERItyp=1;

 ! Allocate arrays and initialize them 
 call clean_allocate('AhCORE',AhCORE,basis%nbf,basis%nbf)
 call clean_allocate('Aoverlap',Aoverlap,basis%nbf,basis%nbf)
 call clean_allocate('NO_COEF',NO_COEF,basis%nbf,basis%nbf)
 call clean_allocate('NO_occ',occ,basis%nbf,1)
 call clean_allocate('NO_energies',energy,basis%nbf,1)
 occ=0.0_dp; energy=0.0_dp;
  ! Save Atomic hCORE integrals and atomic overlaps
 AhCORE(:,:)=AhCORE_in(:,:)
 Aoverlap(:,:)=AOverlap_in(:,:) 
  ! Copy current guess NOs to NO_COEF array
 NO_COEF(:,:)=0.0_dp
 do istate=1,nbf_noft
   NO_COEF(:,istate)=c_matrix(:,istate,1)
 enddo

 ! Not ready for open-shell calcs. (TODO)
 Nvcoupled=ncoupled-1
 Nfrozen=(Nelect-2*Npairs)/2
 Nbeta=Nfrozen+Npairs
 Nalpha=Nbeta
 do
  NBF_occ=Nfrozen+Npairs*(Nvcoupled+1)
  if(NBF_occ<nbf_noft) then
   exit
  else
   Nvcoupled=Nvcoupled-1
  endif
 enddo 

 ! Call module initialization and run NOFT calc.
 if(restartnoft=='yes') then
   call run_noft(INOF,Ista,basis%nbf,NBF_occ,Nfrozen,Npairs,Nvcoupled,Nbeta,Nalpha,iERItyp,&
   & imethocc,imethorb,nscf_nof,iprintdmn,iprintints,ithresh_lambda,ndiis_nof,Enoft,tolE_nof,Vnn,NO_COEF,&
   & Aoverlap,occ(:,1),mo_ints,ofile_name,lowmemERI=(lowmemERI=='yes'),&
   & restart=(restartnoft=='yes'),ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,ireadFdiag=ireadFdiag)
 else
   call run_noft(INOF,Ista,basis%nbf,NBF_occ,Nfrozen,Npairs,Nvcoupled,Nbeta,Nalpha,iERItyp,&
   & imethocc,imethorb,nscf_nof,iprintdmn,iprintints,ithresh_lambda,ndiis_nof,Enoft,tolE_nof,Vnn,NO_COEF,&
   & Aoverlap,occ(:,1),mo_ints,ofile_name,lowmemERI=(lowmemERI=='yes'))
 endif
 
 ! Update c_matrix with optimized NO_COEF
 do istate=1,nbf_noft 
   c_matrix(:,istate,1)=NO_COEF(:,istate)
 enddo

 ! If required print post-procesing files 
 if(print_wfn_ .or. print_cube_ .or. print_wfn_files_ ) then
   call clean_allocate('Occ_print',occ_print,nbf_noft,1)
   occ_print(1:nbf_noft,1)=occ(1:nbf_noft,1)
   if( print_wfn_ )  call plot_wfn(basis,c_matrix)
   if( print_wfn_ )  call plot_rho(basis,occ_print,c_matrix)
   if( print_cube_ ) call plot_cube_wfn('NOFT',basis,occ_print,c_matrix)
   if( print_wfn_files_ ) call print_wfn_file('NOFT',basis,occ_print,c_matrix,Enoft,energy)
   call clean_deallocate('Occ_print',occ_print)
 endif

 ! Deallocate arrays 
 call clean_deallocate('AhCORE',AhCORE)
 call clean_deallocate('Aoverlap',Aoverlap)
 call clean_deallocate('NO_COEF',NO_COEF)
 call clean_deallocate('NO_occ',occ)
 call clean_deallocate('NO_energies',energy)

 ! Stop clock
 call stop_clock(timing_noft_energy)

end subroutine noft_energy

subroutine mo_ints(nbf,nbf_occ,nbf_kji,NO_COEF,hCORE,ERImol,ERImolv)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_hamiltonian_onebody
 implicit none

 integer,intent(in)         :: nbf,nbf_occ,nbf_kji
 real(dp),intent(in)        :: NO_COEF(nbf,nbf)
 real(dp),intent(inout)     :: hCORE(nbf,nbf)
 real(dp),optional,intent(inout) :: ERImol(nbf,nbf_kji,nbf_kji,nbf_kji)
 real(dp),optional,intent(inout) :: ERImolv(nbf*nbf_kji*nbf_kji*nbf_kji)
!====
 integer                    :: istate,jstate,kstate,lstate
 real(dp),allocatable       :: tmp_hcore(:,:)
 real(dp),allocatable       :: tmp_c_matrix(:,:,:)
!=====

 ! hCORE part
 call clean_allocate('tmp_hcore',tmp_hcore,nbf,nbf)
 hCORE(:,:)=0.0_dp; tmp_hcore(:,:)=0.0_dp;
 tmp_hcore=matmul(AhCORE,NO_COEF)
 hCORE=matmul(transpose(NO_COEF),tmp_hcore)
 call clean_deallocate('tmp_hcore',tmp_hcore)

 ! ERI terms
 if(present(ERImol)) then
   ERImol(:,:,:,:)=0.0_dp
   call clean_allocate('tmp_c_matrix',tmp_c_matrix,nbf,nbf_noft,1)
   do istate=1,nbf_noft
    tmp_c_matrix(:,istate,1)=NO_COEF(:,istate)
   enddo
   if(noft_ri) then ! RI case
     call calculate_eri_3center_eigen(tmp_c_matrix,1,nbf_noft,1,nbf_noft)
     do istate=1,nbf_occ
       do jstate=1,nbf_occ
         do kstate=1,nbf_occ
           do lstate=1,nbf_noft
             ERImol(lstate,kstate,jstate,istate)=eri_eigen_ri(lstate,jstate,1,kstate,istate,1) ! <lk|ji> format used for ERImol
           enddo
         enddo
       enddo
     enddo
     call destroy_eri_3center_eigen()
   else            ! Normal case 
    ! TODO
   endif
   call clean_deallocate('tmp_c_matrix',tmp_c_matrix)
 endif

end subroutine mo_ints

!==================================================================
