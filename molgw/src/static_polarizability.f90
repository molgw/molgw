!=========================================================================
! This file contains
! the routine to calculate the static polarizability within RPA
!=========================================================================


!=========================================================================
subroutine static_polarizability(basis,auxil_basis,occupation,energy,wpol_out)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_tools
 use m_basis_set
 use m_spectral_function
 use m_eri
 implicit none

 type(basis_set),intent(in)            :: basis,auxil_basis
 real(dp),intent(in)                   :: occupation(basis%nbf,nspin)
 real(dp),intent(in)                   :: energy(basis%nbf,nspin)
 type(spectral_function),intent(inout) :: wpol_out
!=====
 integer                   :: t_ij,t_kl
 integer                   :: istate,jstate,ijspin
 integer                   :: jbf_auxil,ibf_auxil,ibf_auxil_local
 real(dp),allocatable      :: vsqchi0vsq(:,:),eri_3center_ij(:)
 real(dp)                  :: docc,denom
!=====

 call start_clock(timing_pola_static)

 write(stdout,'(/,a)') ' Calculate the static polarizability within RPA'

 if( .NOT.has_auxil_basis) then
   call die('static_polarizability requires an auxiliary basis')
 endif

 call clean_allocate('static W',wpol_out%w0,nauxil_2center,nauxil_2center)
 
 call clean_allocate('temp chi0 matrix',vsqchi0vsq,nauxil_2center,nauxil_2center)

 allocate(eri_3center_ij(auxil_basis%nbf))

 !
 ! First evaluate v^{1/2} \chi_0 v^{1/2}
 !
 ! Loop over resonant transitions
 vsqchi0vsq(:,:) = 0.0_dp
 do t_ij=1,wpol_out%npole_reso
   istate = wpol_out%transition_table_apb(1,t_ij)
   jstate = wpol_out%transition_table_apb(2,t_ij)
   ijspin = wpol_out%transition_table_apb(3,t_ij)

   docc = occupation(jstate,ijspin) - occupation(istate,ijspin)
   ! Factor 2.0 comes from resonant+antiresonant
   denom = -2.0_dp * docc / ( energy(istate,ijspin) - energy(jstate,ijspin) )

   !
   ! Communicate the needed 3-center integrals
   eri_3center_ij(:) = 0.0_dp
   do ibf_auxil_local=1,auxil_basis%nbf_local
      ibf_auxil = ibf_auxil_g(ibf_auxil_local)
      eri_3center_ij(ibf_auxil) = eri_3center_eigen(ibf_auxil_local,istate,jstate,ijspin)
   enddo
   call xsum(eri_3center_ij)


   do jbf_auxil=1,nauxil_2center
     vsqchi0vsq(:,jbf_auxil) = vsqchi0vsq(:,jbf_auxil) &
          + eri_3center_ij(:) * eri_3center_ij(jbf_auxil) * denom
   enddo

 enddo

 !
 ! Second calculate v^{1/2} \chi v^{1/2} = ( 1 -  v^{1/2} \chi_0 v^{1/2} )^{-1} 
 !                                             * v^{1/2} \chi_0 v^{1/2}
 !
 wpol_out%w0(:,:) = -vsqchi0vsq(:,:)
 forall(jbf_auxil=1:nauxil_2center)
   wpol_out%w0(jbf_auxil,jbf_auxil) = 1.0_dp + wpol_out%w0(jbf_auxil,jbf_auxil)
 end forall

 call invert(nauxil_2center,wpol_out%w0)

 wpol_out%w0(:,:) = MATMUL( wpol_out%w0(:,:) , vsqchi0vsq(:,:) )


 call clean_deallocate('temp chi0 matrix',vsqchi0vsq)
 deallocate(eri_3center_ij)


 call stop_clock(timing_pola_static)

end subroutine static_polarizability


!=========================================================================
