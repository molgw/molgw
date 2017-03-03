!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the routine to calculate the static polarizability within RPA
!
!=========================================================================
subroutine static_polarizability(nstate,occupation,energy,wpol_out)
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
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin)
 type(spectral_function),intent(inout) :: wpol_out
!=====
 integer                   :: t_ia
 integer                   :: istate,astate,iaspin
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

 call clean_allocate('Static W',wpol_out%w0,nauxil_2center,nauxil_2center)
 
 call clean_allocate('temp chi0 matrix',vsqchi0vsq,nauxil_2center,nauxil_2center)


 !
 ! First evaluate v^{1/2} \chi_0 v^{1/2}
 !
 ! Loop over resonant transitions
 vsqchi0vsq(:,:) = 0.0_dp
 do t_ia=1,wpol_out%npole_reso_apb
   istate = wpol_out%transition_table_apb(1,t_ia)
   astate = wpol_out%transition_table_apb(2,t_ia)
   iaspin = wpol_out%transition_table_apb(3,t_ia)

   docc = occupation(astate,iaspin) - occupation(istate,iaspin)
   ! Factor 2.0 comes from resonant+antiresonant
   denom = -2.0_dp * docc / ( energy(istate,iaspin) - energy(astate,iaspin) )

   !
   ! Communicate the needed 3-center integrals
   eri_3center_ij(:) = 0.0_dp
   do ibf_auxil_local=1,nauxil_3center
     ibf_auxil = ibf_auxil_g(ibf_auxil_local)
     eri_3center_ij(ibf_auxil) = eri_3center_eigen(ibf_auxil_local,istate,astate,iaspin)
   enddo
   call xsum_auxil(eri_3center_ij)


   do jbf_auxil=1,nauxil_2center
     if( MODULO( jbf_auxil , nproc_auxil ) /= rank_auxil ) cycle 
     vsqchi0vsq(:,jbf_auxil) = vsqchi0vsq(:,jbf_auxil) &
          + eri_3center_ij(:) * eri_3center_ij(jbf_auxil) * denom
   enddo

 enddo

 call xsum_auxil(vsqchi0vsq)


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
subroutine dynamical_polarizability(nstate,occupation,energy,omega,wpol_in,vsqchi0vsqm1)
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
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin)
 real(dp),intent(in)                   :: omega
 type(spectral_function),intent(in)    :: wpol_in
 real(dp),intent(out)                  :: vsqchi0vsqm1(nauxil_2center,nauxil_2center)
!=====
 integer                   :: t_ia
 integer                   :: istate,astate,iaspin
 integer                   :: jbf_auxil,ibf_auxil,ibf_auxil_local
 real(dp)                  :: eri_3center_ij(nauxil_2center)
 real(dp)                  :: docc,denom
 real(dp)                  :: tmp(nauxil_2center,nauxil_2center)
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
 do t_ia=1,wpol_in%npole_reso_apb
   istate = wpol_in%transition_table_apb(1,t_ia)
   astate = wpol_in%transition_table_apb(2,t_ia)
   iaspin = wpol_in%transition_table_apb(3,t_ia)

   docc = occupation(istate,iaspin) - occupation(astate,iaspin)
   denom = REAL( docc / ( omega - energy(astate,iaspin) + energy(istate,iaspin) + ieta  ) &
                -docc / ( omega + energy(astate,iaspin) - energy(istate,iaspin) - ieta  ) , dp )

   !
   ! Communicate the needed 3-center integrals
   eri_3center_ij(:) = 0.0_dp
   do ibf_auxil_local=1,nauxil_3center
     ibf_auxil = ibf_auxil_g(ibf_auxil_local)
     eri_3center_ij(ibf_auxil) = eri_3center_eigen(ibf_auxil_local,istate,astate,iaspin)
   enddo
   call xsum_auxil(eri_3center_ij)


   do jbf_auxil=1,nauxil_2center
     if( MODULO( jbf_auxil , nproc_auxil ) /= rank_auxil ) cycle 
     vsqchi0vsqm1(:,jbf_auxil) = vsqchi0vsqm1(:,jbf_auxil) &
          + eri_3center_ij(:) * eri_3center_ij(jbf_auxil) * denom
   enddo

 enddo

 call xsum_auxil(vsqchi0vsqm1)

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
subroutine dynamical_polarizability_sca(nstate,occupation,energy,nomega,omega,wpol_in,desc_chi,mlocal,nlocal,chi)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_scalapack
 use m_inputparam
 use m_mpi
 use m_tools
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin)
 integer,intent(in)                    :: nomega
 real(dp),intent(in)                   :: omega(nomega)
 type(spectral_function),intent(in)    :: wpol_in
 integer,intent(in)                    :: desc_chi(NDEL)
 integer,intent(in)                    :: mlocal,nlocal
 real(dp),intent(out)                  :: chi(mlocal,nlocal)
!=====
 integer              :: iomega
 integer              :: ilocal,jlocal
 integer              :: iglobal,jglobal
 integer              :: t_ia
 integer              :: istate,astate,iaspin
 integer              :: info
 real(dp)             :: docc,de,factor_sqrt
 real(dp),allocatable :: eri_3center_t(:,:)
 real(dp),allocatable :: chi0(:,:)
 real(dp),allocatable :: one_m_chi0(:,:)
 real(dp),allocatable :: one_m_chi0m1(:,:)
 integer              :: desc_eri3_t(NDEL)
!=====

 call start_clock(timing_pola_static)

#ifdef HAVE_SCALAPACK
 write(stdout,'(/,a)') ' Calculate the dynamical polarizability within RPA: SCALAPACK'

 if( .NOT. has_auxil_basis ) then
   call die('dynamical_polarizability_sca requires an auxiliary basis')
 endif

 if( nspin == 2 ) then
   call die('dynamical_polarizability_sca: nspin==2 not implemented')
 endif



 do iomega=1,nomega

   call clean_allocate('TMP 3-center MO integrals',eri_3center_t,nauxil_3center,wpol_in%npole_reso_apb)

   call DESCINIT(desc_eri3_t,nauxil_2center,wpol_in%npole_reso_apb,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,nauxil_3center),info)

   !
   ! First evaluate v^{1/2} \chi_0 v^{1/2}
   !
   ! Loop over resonant transitions
   do t_ia=1,wpol_in%npole_reso_apb
     istate = wpol_in%transition_table_apb(1,t_ia)
     astate = wpol_in%transition_table_apb(2,t_ia)
     iaspin = wpol_in%transition_table_apb(3,t_ia)

     docc = occupation(istate,iaspin) - occupation(astate,iaspin)
     de   = energy(astate,iaspin)     - energy(istate,iaspin)
     factor_sqrt = SQRT( 2.0_dp * docc * de / ( omega(iomega)**2 + de**2 ) )

     eri_3center_t(:,t_ia) = eri_3center_eigen(:,istate,astate,iaspin) * factor_sqrt

   enddo

   call clean_allocate('Chi0',chi0,mlocal,nlocal)

   call PDSYRK('L','N',nauxil_2center,wpol_in%npole_reso_apb,1.0_dp,eri_3center_t,1,1,desc_eri3_t,0.0_dp,chi0,1,1,desc_chi)
   chi0(:,:) = -chi0(:,:)


   call clean_deallocate('TMP 3-center MO integrals',eri_3center_t)

   call clean_allocate('1-Chi0',one_m_chi0,mlocal,nlocal)

   ! Symmetrize chi0
   call symmetrize_matrix_sca('L',nauxil_2center,desc_chi,chi0,desc_chi,one_m_chi0)


   one_m_chi0(:,:) = -chi0(:,:)
   do jlocal=1,nlocal
     jglobal = INDXL2G(jlocal,NBLOCK_AUXIL,ipcol_auxil,first_col,npcol_auxil)
     do ilocal=1,mlocal
       iglobal = INDXL2G(ilocal,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
       if( iglobal == jglobal ) one_m_chi0(ilocal,jlocal) = one_m_chi0(ilocal,jlocal) + 1.0_dp
     enddo
   enddo

   call clean_allocate('(1-Chi0)**-1',one_m_chi0m1,mlocal,nlocal)


   call invert_sca(desc_chi,one_m_chi0,one_m_chi0m1)

   call PDGEMM('N','N',nauxil_2center,nauxil_2center,nauxil_2center, &
               1.0_dp,one_m_chi0m1,1,1,desc_chi,    &
                      chi0        ,1,1,desc_chi,    &
               0.0_dp,chi         ,1,1,desc_chi)



   call clean_deallocate('1-Chi0',one_m_chi0)
   call clean_deallocate('(1-Chi0)**-1',one_m_chi0m1)
   call clean_deallocate('Chi0',chi0)


 enddo

#endif

 call stop_clock(timing_pola_static)

end subroutine dynamical_polarizability_sca


!=========================================================================
