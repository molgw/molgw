!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the routine to calculate the static polarizability within RPA
!
!=========================================================================


!=========================================================================
subroutine static_polarizability(nstate,basis,occupation,energy,wpol_out)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_tools
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                    :: nstate
 type(basis_set),intent(in)            :: basis
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin)
 type(spectral_function),intent(inout) :: wpol_out
!=====
 integer                   :: t_ij
 integer                   :: istate,jstate,ijspin
 integer                   :: jbf_auxil,ibf_auxil,ibf_auxil_local
 real(dp),allocatable      :: vsqchi0vsq(:,:)
 real(dp)                  :: eri_3center_ij(nauxil_2center)
 real(dp)                  :: docc,denom
!=====

 call start_clock(timing_pola_static)

 write(stdout,'(/,a)') ' Calculate the static polarizability within RPA'

 if( .NOT. has_auxil_basis ) then
   call die('static_polarizability requires an auxiliary basis')
 endif

 call clean_allocate('static W',wpol_out%w0,nauxil_2center,nauxil_2center)
 
 call clean_allocate('temp chi0 matrix',vsqchi0vsq,nauxil_2center,nauxil_2center)


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
   do ibf_auxil_local=1,nauxil_3center
     ibf_auxil = ibf_auxil_g(ibf_auxil_local)
     eri_3center_ij(ibf_auxil) = eri_3center_eigen(ibf_auxil_local,istate,jstate,ijspin)
   enddo
   call xsum(eri_3center_ij)


   do jbf_auxil=1,nauxil_2center
     if( MODULO( jbf_auxil , nproc ) /= rank ) cycle 
     vsqchi0vsq(:,jbf_auxil) = vsqchi0vsq(:,jbf_auxil) &
          + eri_3center_ij(:) * eri_3center_ij(jbf_auxil) * denom
   enddo

 enddo

 call xsum(vsqchi0vsq)


 !
 ! Second calculate v^{1/2} \chi v^{1/2} = ( 1 -  v^{1/2} \chi_0 v^{1/2} )^{-1} 
 !                                             * v^{1/2} \chi_0 v^{1/2}
 !
 wpol_out%w0(:,:) = -vsqchi0vsq(:,:)
 forall(jbf_auxil=1:nauxil_2center)
   wpol_out%w0(jbf_auxil,jbf_auxil) = 1.0_dp + wpol_out%w0(jbf_auxil,jbf_auxil)
 end forall


 ! TODO I should use SCALAPACK for the next two operations
 call invert(nauxil_2center,wpol_out%w0)
 wpol_out%w0(:,:) = MATMUL( wpol_out%w0(:,:) , vsqchi0vsq(:,:) )


 call clean_deallocate('temp chi0 matrix',vsqchi0vsq)


 call stop_clock(timing_pola_static)

end subroutine static_polarizability


!=========================================================================
subroutine dynamical_polarizability(nstate,basis,occupation,energy,omega,wpol_in,vsqchi0vsqm1)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_tools
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                    :: nstate
 type(basis_set),intent(in)            :: basis
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin)
 real(dp),intent(in)                   :: omega
 type(spectral_function),intent(in)    :: wpol_in
 real(dp),intent(out)                  :: vsqchi0vsqm1(nauxil_2center,nauxil_2center)
!=====
 integer                   :: t_ij
 integer                   :: istate,jstate,ijspin
 integer                   :: jbf_auxil,ibf_auxil,ibf_auxil_local
 real(dp)                  :: eri_3center_ij(nauxil_2center)
 real(dp)                  :: docc,denom
 real(dp) :: tmp(nauxil_2center,nauxil_2center)
!=====

 call start_clock(timing_pola_static)

! write(stdout,'(/,a)') ' Calculate the dynamical polarizability within RPA'

 if( .NOT. has_auxil_basis ) then
   call die('dynamical_polarizability requires an auxiliary basis')
 endif


 !
 ! First evaluate v^{1/2} \chi_0 v^{1/2}
 !
 ! Loop over resonant transitions
 vsqchi0vsqm1(:,:) = 0.0_dp
 do t_ij=1,wpol_in%npole_reso
   istate = wpol_in%transition_table_apb(1,t_ij)
   jstate = wpol_in%transition_table_apb(2,t_ij)
   ijspin = wpol_in%transition_table_apb(3,t_ij)

   docc = occupation(istate,ijspin) - occupation(jstate,ijspin)
   denom = REAL( docc / ( omega - energy(jstate,ijspin) + energy(istate,ijspin) + ieta  ) &
                -docc / ( omega + energy(jstate,ijspin) - energy(istate,ijspin) - ieta  ) , dp )

   !
   ! Communicate the needed 3-center integrals
   eri_3center_ij(:) = 0.0_dp
   do ibf_auxil_local=1,nauxil_3center
     ibf_auxil = ibf_auxil_g(ibf_auxil_local)
     eri_3center_ij(ibf_auxil) = eri_3center_eigen(ibf_auxil_local,istate,jstate,ijspin)
   enddo
   call xsum(eri_3center_ij)


   do jbf_auxil=1,nauxil_2center
     if( MODULO( jbf_auxil , nproc ) /= rank ) cycle 
     vsqchi0vsqm1(:,jbf_auxil) = vsqchi0vsqm1(:,jbf_auxil) &
          + eri_3center_ij(:) * eri_3center_ij(jbf_auxil) * denom
   enddo

 enddo

 call xsum(vsqchi0vsqm1)

! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(1,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(2,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(3,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(4,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(5,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(6,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(7,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(8,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(9,:)

! call invert(nauxil_2center,vsqchi0vsqm1)

 call diagonalize(nauxil_2center,vsqchi0vsqm1,eri_3center_ij)
 forall(jbf_auxil=1:nauxil_2center)
   tmp(:,jbf_auxil) = vsqchi0vsqm1(:,jbf_auxil) * eri_3center_ij(jbf_auxil)
 endforall
 vsqchi0vsqm1(:,:) = MATMUL( vsqchi0vsqm1(:,:) , TRANSPOSE(tmp) )
! write(stdout,*) 'eigenvalues'
! write(stdout,*) eri_3center_ij(:)
! write(stdout,*)
! write(stdout,*) 1.0_dp/eri_3center_ij(:)
!
! write(stdout,*)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(1,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(2,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(3,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(4,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(5,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(6,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(7,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(8,:)
! write(stdout,'(20(es14.6,2x))') vsqchi0vsqm1(9,:)


 call stop_clock(timing_pola_static)

end subroutine dynamical_polarizability


!=========================================================================
