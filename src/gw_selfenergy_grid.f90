!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the calculation of the GW self-energy and RPA polarizability
! on a grid of imaginary frequencies
!
!=========================================================================
subroutine polarizability_grid_scalapack(basis,nstate,occupation,energy,c_matrix,erpa,wpol)
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

 type(basis_set),intent(in)            :: basis
 integer,intent(in)                    :: nstate
 real(dp),intent(in)                   :: occupation(nstate,nspin)
 real(dp),intent(in)                   :: energy(nstate,nspin)
 real(dp),intent(in)                   :: c_matrix(nstate,basis%nbf,nspin)
 real(dp),intent(out)                  :: erpa
 type(spectral_function),intent(inout) :: wpol
!=====
 integer              :: iomega
 integer              :: ilocal,jlocal
 integer              :: iglobal,jglobal
 integer              :: t_ia
 integer              :: istate,astate,iaspin
 integer              :: info
 real(dp)             :: docc,de,factor_sqrt
 real(dp),allocatable :: eri3_t(:,:)
 real(dp),allocatable :: eri3_sca(:,:)
 real(dp),allocatable :: chi0(:,:)
 real(dp),allocatable :: one_m_chi0(:,:)
 real(dp),allocatable :: one_m_chi0m1(:,:)
 real(dp)             :: eigval(nauxil_2center)
 integer              :: desc_eri3_t(NDEL)
 integer              :: desc_eri3_final(NDEL)
 integer              :: meri3,neri3
!=====

 call start_clock(timing_pola_dynamic)

 write(stdout,'(/,1x,a)') 'Calculation of RPA polarizability on imaginary axis grid'
#ifdef HAVE_SCALAPACK
 write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow_sd,' x ',npcol_sd
#endif


 if( wpol%nomega_quad < 1 ) call die('polarizability_grid_sca: manual_imag_axis file should provide a positive integral number of frequencies')

 if( .NOT. has_auxil_basis ) then
   call die('dynamical_polarizability_sca requires an auxiliary basis')
 endif



 wpol%nprodbasis = nauxil_3center
 wpol%mchi = NUMROC(nauxil_2center,block_row,iprow_sd,first_row,nprow_sd)
 wpol%nchi = NUMROC(nauxil_2center,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(wpol%desc_chi,nauxil_2center,nauxil_2center,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,wpol%mchi),info)
 call clean_allocate('Chi',wpol%chi,wpol%mchi,wpol%nchi,wpol%nomega_quad)

 write(stdout,'(1x,a,i7,a,i7)') 'Matrix sizes   ',nauxil_2center,' x ',nauxil_2center
 write(stdout,'(1x,a,i7,a,i7)') 'Distributed in ',wpol%mchi,' x ',wpol%nchi

 if( has_auxil_basis ) call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1)



 !
 ! Get the processor grid included in the input wpol%desc_chi
 meri3 = NUMROC(nauxil_2center     ,wpol%desc_chi(MB_A),iprow_sd,wpol%desc_chi(RSRC_A),nprow_sd)
 neri3 = NUMROC(wpol%npole_reso_apb,wpol%desc_chi(NB_A),ipcol_sd,wpol%desc_chi(CSRC_A),npcol_sd)
 call DESCINIT(desc_eri3_final,nauxil_2center,wpol%npole_reso_apb,wpol%desc_chi(MB_A),wpol%desc_chi(NB_A), &
               wpol%desc_chi(RSRC_A),wpol%desc_chi(CSRC_A),wpol%desc_chi(CTXT_A),MAX(1,meri3),info)

#ifdef HAVE_SCALAPACK
 call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)
#endif
 call clean_allocate('TMP 3-center MO integrals',eri3_t,nauxil_3center,wpol%npole_reso_apb)
 call clean_allocate('Chi0',chi0,wpol%mchi,wpol%nchi)
 call clean_allocate('1-Chi0',one_m_chi0,wpol%mchi,wpol%nchi)
 call clean_allocate('(1-Chi0)**-1',one_m_chi0m1,wpol%mchi,wpol%nchi)

 call DESCINIT(desc_eri3_t,nauxil_2center,wpol%npole_reso_apb,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,nauxil_3center),info)


 erpa = 0.0_dp
 do iomega=1,wpol%nomega_quad

   write(stdout,'(1x,a,i4,a,i4)') 'Loop on frequencies: ',iomega,' / ',wpol%nomega_quad

   !
   ! First evaluate v^{1/2} \chi_0 v^{1/2}
   !
   ! Loop over resonant transitions
   do t_ia=1,wpol%npole_reso_apb
     istate = wpol%transition_table_apb(1,t_ia)
     astate = wpol%transition_table_apb(2,t_ia)
     iaspin = wpol%transition_table_apb(3,t_ia)

     docc = occupation(istate,iaspin) - occupation(astate,iaspin)
     de   = energy(astate,iaspin)     - energy(istate,iaspin)
     factor_sqrt = SQRT( 2.0_dp * docc * de / ( wpol%omega_quad(iomega)**2 + de**2 ) )

     eri3_t(:,t_ia) = eri_3center_eigen(:,istate,astate,iaspin) * factor_sqrt

   enddo

#ifdef HAVE_SCALAPACK
   call PDGEMR2D(nauxil_2center,wpol%npole_reso_apb,eri3_t,1,1,desc_eri3_t, &
                                                  eri3_sca,1,1,desc_eri3_final,wpol%desc_chi(CTXT_A))
#endif

#ifdef HAVE_SCALAPACK
   call PDSYRK('L','N',nauxil_2center,wpol%npole_reso_apb,1.0_dp,eri3_sca,1,1,desc_eri3_final,0.0_dp,chi0,1,1,wpol%desc_chi)
#else
   call DSYRK('L','N',nauxil_2center,wpol%npole_reso_apb,1.0_dp,eri3_t,nauxil_2center,0.0_dp,chi0,nauxil_2center)
#endif
   chi0(:,:) = -chi0(:,:)



   ! Symmetrize chi0
   call symmetrize_matrix_sca('L',nauxil_2center,wpol%desc_chi,chi0,wpol%desc_chi,one_m_chi0)


   one_m_chi0(:,:) = -chi0(:,:)
   do jlocal=1,wpol%nchi
     jglobal = colindex_local_to_global_descriptor(wpol%desc_chi,jlocal)
     do ilocal=1,wpol%mchi
       iglobal = rowindex_local_to_global_descriptor(wpol%desc_chi,ilocal)
       if( iglobal == jglobal ) one_m_chi0(ilocal,jlocal) = one_m_chi0(ilocal,jlocal) + 1.0_dp
     enddo
   enddo


   one_m_chi0m1(:,:) = one_m_chi0(:,:)

   ! Diagonalize (1-chi0) in order to have RPA total energy.
   ! might be time-consuming and not necessary for most applications
   ! TODO: add an input variable to trigger or not its calculation
   call diagonalize_sca(nauxil_2center,wpol%desc_chi,one_m_chi0m1,eigval)
   erpa = erpa + SUM( LOG(eigval(:)) + 1.0_dp - eigval(:) ) / (2.0_dp * pi) * wpol%weight_quad(iomega)


   call invert_sca(wpol%desc_chi,one_m_chi0,one_m_chi0m1)

#ifdef HAVE_SCALAPACK
   call PDGEMM('N','N',nauxil_2center,nauxil_2center,nauxil_2center, &
               1.0_dp,one_m_chi0m1        ,1,1,wpol%desc_chi,    &
                      chi0                ,1,1,wpol%desc_chi,    &
               0.0_dp,wpol%chi(:,:,iomega),1,1,wpol%desc_chi)
#else
   call DGEMM('N','N',nauxil_2center,nauxil_2center,nauxil_2center, &
              1.0_dp,one_m_chi0m1,nauxil_2center, &
                     chi0        ,nauxil_2center, &
              0.0_dp,wpol%chi(:,:,iomega),nauxil_2center)
#endif


 enddo

#ifdef HAVE_SCALAPACK
 call clean_deallocate('TMP 3-center MO integrals',eri3_sca)
#endif
 call clean_deallocate('TMP 3-center MO integrals',eri3_t)
 call clean_deallocate('1-Chi0',one_m_chi0)
 call clean_deallocate('(1-Chi0)**-1',one_m_chi0m1)
 call clean_deallocate('Chi0',chi0)

 call destroy_eri_3center_eigen()

 write(stdout,'(/,1x,a,f16.10)') 'RPA correlation energy (Ha): ',erpa

 call stop_clock(timing_pola_dynamic)


end subroutine polarizability_grid_scalapack


!=========================================================================
subroutine gw_selfenergy_imag_scalapack(basis,nstate,occupation,energy,c_matrix,wpol,se)
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
 use m_selfenergy_tools
 implicit none

 type(basis_set),intent(in)          :: basis
 integer,intent(in)                  :: nstate
 real(dp),intent(in)                 :: occupation(nstate,nspin)
 real(dp),intent(in)                 :: energy(nstate,nspin)
 real(dp),intent(in)                 :: c_matrix(nstate,basis%nbf,nspin)
 type(spectral_function),intent(in)  :: wpol
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer              :: iomegas
 integer              :: iomega
 integer              :: ilocal,jlocal
 integer              :: iglobal,jglobal
 integer              :: t_ia
 integer              :: istate,astate,iaspin
 integer              :: info
 real(dp)             :: docc,de,factor_sqrt
 real(dp),allocatable :: eri3_t(:,:)
 real(dp),allocatable :: eri3_sca(:,:)
 real(dp),allocatable :: chi_eri3_sca(:,:)
 integer              :: desc_eri3_t(NDEL)
 integer              :: iprow,ipcol,nprow,npcol
 integer              :: desc_eri3_final(NDEL)
 integer              :: meri3,neri3
 integer              :: mstate,pstate,mpspin
 integer              :: prange,plocal
!=====


 if( .NOT. has_auxil_basis ) then
   call die('gw_selfenergy_imag_sca requires an auxiliary basis')
 endif

 call start_clock(timing_self)

 write(stdout,'(/,1x,a)') 'GW self-energy on a grid of imaginary frequencies'
 call BLACS_GRIDINFO(wpol%desc_chi(CTXT_A),nprow,npcol,iprow,ipcol)

#ifdef HAVE_SCALAPACK
 write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow,' x ',npcol
#endif


 if( has_auxil_basis ) call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_G+1,nvirtual_G-1,nsemin,nsemax)


 prange = nvirtual_G - ncore_G - 1

 ! Get the processor grid included in the input wpol%desc_chi
 meri3 = NUMROC(nauxil_2center,wpol%desc_chi(MB_A),iprow,wpol%desc_chi(RSRC_A),nprow)
 neri3 = NUMROC(prange        ,wpol%desc_chi(NB_A),ipcol,wpol%desc_chi(CSRC_A),npcol)
 call DESCINIT(desc_eri3_final,nauxil_2center,prange,wpol%desc_chi(MB_A),wpol%desc_chi(NB_A), &
               wpol%desc_chi(RSRC_A),wpol%desc_chi(CSRC_A),wpol%desc_chi(CTXT_A),MAX(1,meri3),info)

 call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)
 call clean_allocate('TMP 3-center MO integrals',chi_eri3_sca,meri3,neri3)

 call DESCINIT(desc_eri3_t,nauxil_2center,prange,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,nauxil_3center),info)

 se%sigmai(:,:,:) = 0.0_dp

 do mpspin=1,nspin
   do mstate=nsemin,nsemax

#ifdef HAVE_SCALAPACK
     call PDGEMR2D(nauxil_2center,prange,eri_3center_eigen(:,:,mstate,mpspin),1,1,desc_eri3_t, &
                                                                     eri3_sca,1,1,desc_eri3_final,wpol%desc_chi(CTXT_A))
#else
     eri3_sca(:,1:prange) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,mstate,mpspin)
#endif


     do iomega=1,wpol%nomega_quad
#ifdef HAVE_SCALAPACK
       call PDGEMM('N','N',nauxil_2center,prange,nauxil_2center,     &
                   1.0_dp,wpol%chi(:,:,iomega),1,1,wpol%desc_chi,    &
                          eri3_sca            ,1,1,desc_eri3_final,  &
                   0.0_dp,chi_eri3_sca        ,1,1,desc_eri3_final)
#else
       call DGEMM('N','N',nauxil_2center,prange,nauxil_2center,  &
                  1.0_dp,wpol%chi(:,:,iomega),nauxil_2center,    &
                         eri3_sca            ,nauxil_2center,    &
                  0.0_dp,chi_eri3_sca        ,nauxil_2center)
#endif

       do iomegas=0,se%nomegai
         do plocal=1,neri3
           pstate = INDXL2G(plocal,wpol%desc_chi(NB_A),ipcol,wpol%desc_chi(CSRC_A),npcol) + ncore_G
           se%sigmai(iomegas,mstate,mpspin) = se%sigmai(iomegas,mstate,mpspin) &
                         - wpol%weight_quad(iomega) * (  1.0_dp / ( ( se%energy0(mstate,mpspin) + se%omegai(iomegas) - energy(pstate,mpspin) ) &
                                                                  + im * wpol%omega_quad(iomega) )   &
                                                       + 1.0_dp / ( ( se%energy0(mstate,mpspin) + se%omegai(iomegas) - energy(pstate,mpspin) )  &
                                                                  - im * wpol%omega_quad(iomega) )  ) &
                            * DOT_PRODUCT( eri3_sca(:,plocal) , chi_eri3_sca(:,plocal) )  &
                                  /  (2.0_dp * pi)
         enddo
       enddo

     enddo

   enddo
 enddo
 call xsum_world(se%sigmai)

 forall(iomegas=1:se%nomegai)
   se%sigmai(-iomegas,:,:) = CONJG( se%sigmai(iomegas,:,:) )
 end forall

 call clean_deallocate('TMP 3-center MO integrals',eri3_sca)
 call clean_deallocate('TMP 3-center MO integrals',chi_eri3_sca)

 call destroy_eri_3center_eigen()

 call stop_clock(timing_self)

end subroutine gw_selfenergy_imag_scalapack


!=========================================================================
