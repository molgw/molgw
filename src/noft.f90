!==================================================================
! This file is part of MOLGW.
! Author: Mauricio Rodriguez-Mayorga
!
! This file contains
! - NOFT energy opt. with Resolution-of-Identity
!=========================================================================
subroutine noft_energy(nstate,basis,c_matrix,AhCORE_in,AOverlap_in,enoft,Vnn)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore
 use m_hamiltonian_onebody
 use m_noft_driver
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(inout)     :: AhCORE_in(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: AOverlap_in(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: Vnn
 real(dp),intent(out)       :: enoft
!====
 integer                    :: istate
 real(dp),allocatable       :: NO_COEF(:,:),occ(:,:),energy(:,:) 
 logical::lrestart=.true.
 integer::INOF,Ista,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta,Nalpha,itermax,NTHRESHL,NDIIS
 integer::imethorb,imethocc,iprintdmn,iprintints
 real(dp)::tolE
 external::mo_ints
!=====

 ! Init clock
 call start_clock(timing_noft_energy)

 write(stdout,'(/,a)') ' RI-NOFT calculation'

 nbf_noft=nstate
 enoft = 0.0_dp
 ! These can be fixed for a while...
 itermax=1000;NTHRESHL=4;NDIIS=6;tolE=1.0d-8;imethorb=1;

 ! Allocate arrays and initialize them 
 call clean_allocate('AhCORE',AhCORE,basis%nbf,basis%nbf)
 call clean_allocate('Aoverlap',Aoverlap,basis%nbf,basis%nbf)
 call clean_allocate('NO_COEF',NO_COEF,basis%nbf,basis%nbf)
 call clean_allocate('NO_occ',occ,nbf_noft,1)
 call clean_allocate('NO_energies',energy,nbf_noft,1)
 occ=0.0_dp; energy=0.0_dp;
  ! Save Atomic hCORE integrals and atomic overlaps
 AhCORE(:,:)=AhCORE_in(:,:)
 Aoverlap(:,:)=AOverlap_in(:,:) 
 NO_COEF(:,:)=0.0_dp
  ! Copy current guess NOs to NO_COEF array
 do istate=1,nbf_noft
   NO_COEF(:,istate)=c_matrix(:,istate,1)
 enddo

 ! To be input parameters! TODO 
 INOF=7;Ista=1;NBF_occ=10;Nfrozen=0;Npairs=1;Ncoupled=9;Nbeta=1;Nalpha=Nbeta;imethocc=1;iprintdmn=1;iprintints=0 

 ! Call module initialization and run NOFT calc.
 call run_noft(INOF,Ista,nbf_noft,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta,Nalpha,1,imethocc,imethorb,itermax,&
 & iprintdmn,iprintints,NTHRESHL,NDIIS,enoft,tolE,Vnn,NO_COEF,Aoverlap,occ(:,1),mo_ints,&
 & restart=LRESTART,ireadGAMMAS=1,ireadOCC=1,ireadCOEF=1,ireadFdiag=1)

 ! Update c_matrix with optimized NO_COEF
 do istate=1,nbf_noft
   c_matrix(:,istate,1)=NO_COEF(:,istate)
 enddo

 ! If required print post-procesing files 
 if( print_wfn_ )  call plot_wfn(basis,c_matrix)
 if( print_wfn_ )  call plot_rho(basis,occ,c_matrix)
 if( print_cube_ ) call plot_cube_wfn('NOFT',basis,occ,c_matrix)
 if( print_wfn_files_ ) call print_wfn_file('NOFT',basis,occ,c_matrix,enoft,energy)

 ! Deallocate arrays 
 call clean_deallocate('AhCORE',AhCORE)
 call clean_deallocate('Aoverlap',Aoverlap)
 call clean_deallocate('NO_COEF',NO_COEF)
 call clean_deallocate('NO_occ',occ)
 call clean_deallocate('NO_energies',energy)

 ! Stop clock
 call stop_clock(timing_noft_energy)

end subroutine noft_energy

subroutine mo_ints(nbf,NO_COEF,hCORE,ERImol)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_hamiltonian_onebody
 implicit none

 integer,intent(in)         :: nbf
 real(dp),intent(in)        :: NO_COEF(nbf,nbf)
 real(dp),intent(inout)     :: hCORE(nbf,nbf)
 real(dp),intent(inout)     :: ERImol(nbf,nbf,nbf,nbf)
!====
 integer                    :: istate,jstate,kstate,lstate
 real(dp),allocatable       :: tmp_hcore(:,:)
 real(dp),allocatable       :: tmp_c_matrix(:,:,:)
!=====

 ! hCORE part
 call clean_allocate('tmp_hcore',tmp_hcore,nbf,nbf)
 hCORE(:,:)=0.0d0; tmp_hcore(:,:)=0.0d0;
 tmp_hcore=matmul(AhCORE,NO_COEF)
 hCORE=matmul(transpose(NO_COEF),tmp_hcore)
 call clean_deallocate('tmp_hcore',tmp_hcore)

 ! ERI terms
 ERImol(:,:,:,:)=0.0d0
 call clean_allocate('tmp_c_matrix',tmp_c_matrix,nbf,nbf_noft,1)
 do istate=1,nbf_noft
  tmp_c_matrix(:,istate,1)=NO_COEF(:,istate)
 enddo
 if(noft_ri) then ! RI case
   call calculate_eri_3center_eigen(tmp_c_matrix,1,nbf_noft,1,nbf_noft)
   do istate=1,nbf_noft
     do jstate=1,nbf_noft
       do kstate=1,nbf_noft
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

end subroutine mo_ints

!==================================================================
