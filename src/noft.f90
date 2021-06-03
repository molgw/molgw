!==================================================================
! This file is part of MOLGW.
! Author: Mauricio Rodriguez-Mayorga
!
! This file contains
! - NOFT energy opt. with Resolution-of-Identity
!=========================================================================
subroutine noft_energy_ri(nstate,basis,c_matrix,hCORE,Overlap,enoft)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore
 use m_hamiltonian_onebody
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(inout)     :: hCORE(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: Overlap(basis%nbf,basis%nbf)
 real(dp),intent(out)       :: enoft
!====
 integer                    :: istate,jstate,kstate,lstate
 real(dp)                   :: tmp_iajb,energy_denom
 real(dp),allocatable       :: NO_COEF(:,:)
 real(dp),allocatable       :: ERImol(:,:,:,:)
 integer                    :: nocc(nspin)
 integer                    :: ncore
 external                   :: mo_ints
!=====

 call start_clock(timing_noft_energy)

 write(stdout,'(/,a)') ' RI-NOFT calculation'

 allocate(AhCORE(basis%nbf,basis%nbf),Aoverlap(basis%nbf,basis%nbf),NO_COEF(basis%nbf,basis%nbf))
 allocate(ERImol(basis%nbf,basis%nbf,basis%nbf,basis%nbf))

 AhCORE(:,:)=hCORE(:,:)
 Aoverlap(:,:)=overlap(:,:) 
 NO_COEF(:,:)=0.0_dp
 do istate=1,basis%nbf
  NO_COEF(istate,1:nstate)=c_matrix(istate,1:nstate,1)
 enddo

 ncore = ncoreg
 if(is_frozencore) then
   if( ncore == 0) ncore = atoms_core_states()
 endif

 call calculate_eri_3center_eigen(c_matrix,ncore+1,nstate,ncore+1,nstate)

 enoft = 0.0_dp
 
 call mo_ints(basis%nbf,NO_COEF,hCORE,ERImol)
 !tmp_iajb = eri_eigen_ri(istate,astate,iaspin,jstate,bstate,jbspin)
 tmp_iajb = eri_eigen_ri(1,2,1,1,2,1)

 enoft = enoft + 0.5_dp * energy_denom * tmp_iajb**2


 write(stdout,'(a,f16.10,/)') ' NOFT :',enoft

 call destroy_eri_3center_eigen()
 deallocate(AhCORE,Aoverlap,NO_COEF)
 deallocate(ERImol)
 call stop_clock(timing_noft_energy)

end subroutine noft_energy_ri

subroutine mo_ints(nbf,NO_COEF,hcoreMOL,ERImol)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_hamiltonian_onebody
 implicit none

 integer,intent(in)         :: nbf
 real(dp),intent(in)        :: NO_COEF(nbf,nbf)
 real(dp),intent(inout)     :: hcoreMOL(nbf,nbf)
 real(dp),intent(inout)     :: ERImol(nbf,nbf,nbf,nbf)
!====
 integer                    :: istate

 write(*,*) 'IN HERE'
 do istate=1,nbf
  write(*,*) AhCORE(istate,:)
 enddo

end subroutine mo_ints

!==================================================================
