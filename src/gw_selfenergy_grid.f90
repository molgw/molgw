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
 real(dp),intent(in)                   :: c_matrix(basis%nbf,nstate,nspin)
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
 real(dp),allocatable :: chi0(:,:)
 real(dp),allocatable :: one_m_chi0(:,:)
 real(dp),allocatable :: one_m_chi0m1(:,:)
 real(dp)             :: eigval(nauxil_2center)
 real(dp)             :: eint
 integer              :: desc_eri3_t(NDEL)
 integer              :: desc_eri3_final(NDEL)
 integer              :: meri3,neri3
#ifdef HAVE_SCALAPACK
 real(dp),allocatable :: eri3_sca(:,:)
#endif
!=====

 call start_clock(timing_rpa_dynamic)

 write(stdout,'(/,1x,a)') 'Calculation of RPA polarizability on imaginary axis grid'
#ifdef HAVE_SCALAPACK
 write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow_sd,' x ',npcol_sd
#endif


 if( wpol%nomega_quad < 1 ) call die('polarizability_grid_sca: nomega_imag input variable should be greater than 1')

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

 if( has_auxil_basis ) call calculate_eri_3center_eigen(c_matrix,ncore_W+1,nhomo_W,nlumo_W,nvirtual_W-1)



 !
 ! Get the processor grid included in the input wpol%desc_chi
 meri3 = NUMROC(nauxil_2center     ,wpol%desc_chi(MB_),iprow_sd,wpol%desc_chi(RSRC_),nprow_sd)
 neri3 = NUMROC(wpol%npole_reso_apb,wpol%desc_chi(NB_),ipcol_sd,wpol%desc_chi(CSRC_),npcol_sd)
 call DESCINIT(desc_eri3_final,nauxil_2center,wpol%npole_reso_apb,wpol%desc_chi(MB_),wpol%desc_chi(NB_), &
               wpol%desc_chi(RSRC_),wpol%desc_chi(CSRC_),wpol%desc_chi(CTXT_),MAX(1,meri3),info)

#ifdef HAVE_SCALAPACK
 call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)
#endif
 call clean_allocate('TMP 3-center MO integrals',eri3_t,nauxil_3center,wpol%npole_reso_apb)
 call clean_allocate('Chi0',chi0,wpol%mchi,wpol%nchi)
 call clean_allocate('1-Chi0',one_m_chi0,wpol%mchi,wpol%nchi)
 call clean_allocate('(1-Chi0)**-1',one_m_chi0m1,wpol%mchi,wpol%nchi)

 call DESCINIT(desc_eri3_t,nauxil_2center,wpol%npole_reso_apb,MB_auxil,NB_auxil,first_row,first_col,cntxt_auxil,MAX(1,nauxil_3center),info)


 erpa = 0.0_dp
 eint = 0.0_dp
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
                                                    eri3_sca,1,1,desc_eri3_final,wpol%desc_chi(CTXT_))
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
   ! might be a bit time-consuming but we only calculate the eigenvalues
   call diagonalize_eigval_sca(postscf_diago_flavor,one_m_chi0m1,wpol%desc_chi,eigval)
   erpa = erpa + SUM( LOG(eigval(:)) + 1.0_dp - eigval(:) ) / (2.0_dp * pi) * wpol%weight_quad(iomega)
   eint = eint + SUM( -( 1.0_dp - eigval(:) ) / eigval(:) + 1.0_dp - eigval(:) ) / (2.0_dp * pi) * wpol%weight_quad(iomega)

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
 write(stdout,'(1x,a,f16.10)')   'GW  correlation energy (Ha): ',eint

 call stop_clock(timing_rpa_dynamic)


end subroutine polarizability_grid_scalapack


!=========================================================================
subroutine gw_selfenergy_imag_scalapack(basis,nstate,energy,c_matrix,wpol,se)
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
 real(dp),intent(in)                 :: energy(nstate,nspin)
 real(dp),intent(in)                 :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in)  :: wpol
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer              :: iomegas
 integer              :: iomega
 integer              :: info
 real(dp),allocatable :: eri3_sca(:,:)
 real(dp),allocatable :: chi_eri3_sca(:,:)
 real(dp)             :: v_chi_v_p
 integer              :: desc_eri3_t(NDEL)
 integer              :: iprow,ipcol,nprow,npcol
 integer              :: desc_eri3_final(NDEL)
 integer              :: meri3,neri3
 integer              :: mstate,pstate,mpspin
 integer              :: prange,plocal
 complex(dp),allocatable :: sigmaigw(:,:,:)
!=====


 if( .NOT. has_auxil_basis ) then
   call die('gw_selfenergy_imag_sca requires an auxiliary basis')
 endif

 call start_clock(timing_gw_self)

 write(stdout,'(/,1x,a)') 'GW self-energy on a grid of imaginary frequencies'

 nprow = 1
 npcol = 1
#ifdef HAVE_SCALAPACK
 ! Get the processor grid included in the input wpol%desc_chi
 call BLACS_GRIDINFO(wpol%desc_chi(CTXT_),nprow,npcol,iprow,ipcol)
 write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow,' x ',npcol
#endif


 if( has_auxil_basis ) call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,nsemin,nsemax)


 prange = nvirtual_G - ncore_G - 1

 meri3 = NUMROC(nauxil_2center,wpol%desc_chi(MB_),iprow,wpol%desc_chi(RSRC_),nprow)
 neri3 = NUMROC(prange        ,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol)
 call DESCINIT(desc_eri3_final,nauxil_2center,prange,wpol%desc_chi(MB_),wpol%desc_chi(NB_), &
               wpol%desc_chi(RSRC_),wpol%desc_chi(CSRC_),wpol%desc_chi(CTXT_),MAX(1,meri3),info)

 call clean_allocate('TMP 3-center MO integrals',eri3_sca,meri3,neri3)
 call clean_allocate('TMP 3-center MO integrals',chi_eri3_sca,meri3,neri3)

 call DESCINIT(desc_eri3_t,nauxil_2center,prange,MB_auxil,NB_auxil,first_row,first_col,cntxt_auxil,MAX(1,nauxil_3center),info)


 ! OPENMP does not want to reduce se%sigmai then work with a temporary array sigmaigw
 allocate(sigmaigw(0:se%nomegai,nsemin:nsemax,nspin))
 sigmaigw(:,:,:) = 0.0_dp

 do mpspin=1,nspin
   do mstate=nsemin,nsemax

#ifdef HAVE_SCALAPACK
     call PDGEMR2D(nauxil_2center,prange,eri_3center_eigen(:,:,mstate,mpspin),1,1,desc_eri3_t, &
                                                                     eri3_sca,1,1,desc_eri3_final,wpol%desc_chi(CTXT_))
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

       !$OMP PARALLEL PRIVATE(pstate, v_chi_v_p)
       !$OMP DO REDUCTION(+:sigmaigw)
       do plocal=1,neri3
         pstate = INDXL2G(plocal,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol) + ncore_G

         v_chi_v_p = DOT_PRODUCT( eri3_sca(:,plocal) , chi_eri3_sca(:,plocal) )

         sigmaigw(:,mstate,mpspin) = sigmaigw(:,mstate,mpspin) &
                       - wpol%weight_quad(iomega) * (  1.0_dp / ( ( se%energy0(mstate,mpspin) + se%omegai(0:) - energy(pstate,mpspin) ) &
                                                                + im * wpol%omega_quad(iomega) )   &
                                                     + 1.0_dp / ( ( se%energy0(mstate,mpspin) + se%omegai(0:) - energy(pstate,mpspin) )  &
                                                                - im * wpol%omega_quad(iomega) )  ) &
                          * v_chi_v_p /  (2.0_dp * pi)
       enddo
       !$OMP END DO
       !$OMP END PARALLEL

     enddo

   enddo
 enddo
 call xsum_world(sigmaigw)

 se%sigmai(0:,:,:) = sigmaigw(0:,:,:)
 forall(iomegas=1:se%nomegai)
   se%sigmai(-iomegas,:,:) = CONJG( sigmaigw(iomegas,:,:) )
 end forall

 call clean_deallocate('TMP 3-center MO integrals',eri3_sca)
 call clean_deallocate('TMP 3-center MO integrals',chi_eri3_sca)

 call destroy_eri_3center_eigen()

 call stop_clock(timing_gw_self)

end subroutine gw_selfenergy_imag_scalapack


!=========================================================================
subroutine gw_selfenergy_imag_scalapack_full(basis,nstate,energy,c_matrix,wpol,se)
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
 real(dp),intent(in)                 :: energy(nstate,nspin)
 real(dp),intent(in)                 :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in)  :: wpol
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer              :: iomegas
 integer              :: iomega
 integer              :: info
 real(dp),allocatable :: eri3_sca_p(:,:),eri3_sca_q(:,:)
 real(dp),allocatable :: chi_eri3_sca_q(:,:)
 real(dp)             :: v_chi_v_pq
 integer              :: desc_eri3_t(NDEL)
 integer              :: iprow,ipcol,nprow,npcol
 integer              :: desc_eri3_final(NDEL)
 integer              :: meri3,neri3
 integer              :: mstate,pstate,qstate,mpspin
 integer              :: prange,mlocal
 complex(dp),allocatable :: sigmaigw(:,:,:,:)
!=====


 if( .NOT. has_auxil_basis ) then
   call die('gw_selfenergy_imag_sca requires an auxiliary basis')
 endif

 call start_clock(timing_gw_self)

 write(stdout,'(/,1x,a)') 'GW self-energy on a grid of imaginary frequencies'

 nprow = 1
 npcol = 1
#ifdef HAVE_SCALAPACK
 ! Get the processor grid included in the input wpol%desc_chi
 call BLACS_GRIDINFO(wpol%desc_chi(CTXT_),nprow,npcol,iprow,ipcol)
 write(stdout,'(1x,a,i4,a,i4)') 'SCALAPACK grid',nprow,' x ',npcol
#endif


 if( has_auxil_basis ) call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,nsemin,nsemax)


 prange = nvirtual_G - ncore_G - 1

 meri3 = NUMROC(nauxil_2center,wpol%desc_chi(MB_),iprow,wpol%desc_chi(RSRC_),nprow)
 neri3 = NUMROC(prange        ,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol)
 call DESCINIT(desc_eri3_final,nauxil_2center,prange,wpol%desc_chi(MB_),wpol%desc_chi(NB_), &
               wpol%desc_chi(RSRC_),wpol%desc_chi(CSRC_),wpol%desc_chi(CTXT_),MAX(1,meri3),info)

 call clean_allocate('TMP 3-center MO integrals',eri3_sca_p,meri3,neri3)
 call clean_allocate('TMP 3-center MO integrals',eri3_sca_q,meri3,neri3)
 call clean_allocate('TMP 3-center MO integrals',chi_eri3_sca_q,meri3,neri3)

 call DESCINIT(desc_eri3_t,nauxil_2center,prange,MB_auxil,NB_auxil,first_row,first_col,cntxt_auxil,MAX(1,nauxil_3center),info)


 ! OPENMP does not want to reduce se%sigmai then work with a temporary array sigmaigw
 allocate(sigmaigw(-se%nomegai:se%nomegai,nsemin:nsemax,nsemin:nsemax,nspin))
 sigmaigw(:,:,:,:) = 0.0_dp

 do mpspin=1,nspin
   do qstate=nsemin,nsemax

     eri3_sca_q(:,1:prange) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,qstate,mpspin)


     do iomega=1,wpol%nomega_quad

       call DGEMM('N','N',nauxil_2center,prange,nauxil_2center,  &
                  1.0_dp,wpol%chi(:,:,iomega),nauxil_2center,    &
                         eri3_sca_q          ,nauxil_2center,    &
                  0.0_dp,chi_eri3_sca_q      ,nauxil_2center)

   do pstate=nsemin,nsemax
     eri3_sca_p(:,1:prange) = eri_3center_eigen(:,ncore_G+1:nvirtual_G-1,pstate,mpspin)

       !$OMP PARALLEL PRIVATE(mstate, v_chi_v_pq)
       !$OMP DO REDUCTION(+:sigmaigw)
       do mlocal=1,neri3
         mstate = INDXL2G(mlocal,wpol%desc_chi(NB_),ipcol,wpol%desc_chi(CSRC_),npcol) + ncore_G

         v_chi_v_pq = DOT_PRODUCT( eri3_sca_p(:,mlocal) , chi_eri3_sca_q(:,mlocal) )

         sigmaigw(:,pstate,qstate,mpspin) = sigmaigw(:,pstate,qstate,mpspin) &
                       - wpol%weight_quad(iomega) * (  1.0_dp / ( ( se%omegai(:) - energy(mstate,mpspin) ) &
                                                                + im * wpol%omega_quad(iomega) )   &
                                                     + 1.0_dp / ( ( se%omegai(:) - energy(mstate,mpspin) )  &
                                                                - im * wpol%omega_quad(iomega) )  ) &
                          * v_chi_v_pq /  (2.0_dp * pi)
       enddo
       !$OMP END DO
       !$OMP END PARALLEL

   enddo

     enddo

   enddo
 enddo
 call xsum_world(sigmaigw)

! do iomegas=1,se%nomegai
!   sigmaigw(-iomegas,:,:,:) = CONJG( sigmaigw(iomegas,:,:,:) )
! enddo


 write(stdout,*) 'Calculate GW density matrix'

 do iomegas=-2,2
   write(stdout,*) iomegas,AIMAG(se%omegai(iomegas))
   write(stdout,*) sigmaigw(iomegas,4:7,4:7,1)
   write(stdout,*)
 enddo

 pstate=7
 qstate=7
 do iomegas=-se%nomegai,se%nomegai
   write(stdout,*) iomegas,AIMAG(se%omegai(iomegas)),sigmaigw(iomegas,pstate,qstate,1) / (2.0_dp * pi)  &
                       / ( se%omegai(iomegas) - energy(pstate,1) ) / ( se%omegai(iomegas) - energy(qstate,1) )
 enddo
 do pstate=4,7
! do qstate=4,7
 qstate = pstate
 write(stdout,*) pstate,qstate,SUM( sigmaigw(:,pstate,qstate,1) / pi  &
                       / ( se%omegai(:) - energy(pstate,1) ) / ( se%omegai(:) - energy(qstate,1) ) ) * ABS(se%omegai(1) - se%omegai(0)) 
! enddo
 enddo

 stop 'ENOUGH'
 
 call clean_deallocate('TMP 3-center MO integrals',eri3_sca_p)
 call clean_deallocate('TMP 3-center MO integrals',eri3_sca_q)
 call clean_deallocate('TMP 3-center MO integrals',chi_eri3_sca_q)

 call destroy_eri_3center_eigen()


 call stop_clock(timing_gw_self)

end subroutine gw_selfenergy_imag_scalapack_full


!=========================================================================
