!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the GW self-energy
! within different flavors: G0W0, GnW0, GnWn, COHSEX, QSGW
!
!=========================================================================
subroutine gw_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,wpol,se,energy_gw)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_timing
 use m_inputparam
 use m_warning,only: issue_warning
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_tools,only: coeffs_gausslegint,diagonalize
 use m_selfenergy_tools
 implicit none

 integer,intent(in)                  :: nstate,selfenergy_approx
 type(basis_set)                     :: basis
 real(dp),intent(in)                 :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)                 :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in)  :: wpol
 type(selfenergy_grid),intent(inout) :: se
 real(dp),intent(out)                :: energy_gw
!=====
 integer               :: iomega
 integer               :: ipstate
 integer               :: pstate,bstate
 integer               :: istate,ispin,ipole
 real(dp),allocatable  :: bra(:,:)
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: fact_full_a,fact_empty_a
!=====

 call start_clock(timing_gw_self)

 write(stdout,*)
 select case(selfenergy_approx)
 case(GW)
   write(stdout,*) 'Perform a one-shot G0W0 calculation'
 case(ONE_RING)
   write(stdout,*) 'Perform a one-shot one-ring calculation'
 case(COHSEX)
   write(stdout,*) 'Perform a COHSEX calculation'
 case(GnW0)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnW0 calculation'
 case(GnWn)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnWn calculation'
 case default
   write(stdout,*) 'type:',selfenergy_approx
   call die('gw_selfenergy: calculation type unknown')
 end select


 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,nsemin,nsemax,ncore_G+1,nvirtual_G-1)
 endif


 call clean_allocate('Temporary array',bra,1,wpol%npole_reso,nsemin,nsemax)


 energy_gw = 0.0_dp


 se%sigma(:,:,:) = 0.0_dp

 do ispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !INNER LOOP of G

     if( MODULO( istate - (ncore_G+1) , nproc_ortho) /= rank_ortho ) cycle

     !
     ! Prepare the bra and ket with the knowledge of index istate and pstate
     if( .NOT. has_auxil_basis) then
       !$OMP PARALLEL
       !$OMP DO PRIVATE(ipstate)
       ! Here just grab the precalculated value
       do pstate=nsemin,nsemax
         ipstate = index_prodstate(istate,pstate) + (ispin-1) * index_prodstate(nvirtual_W-1,nvirtual_W-1)
         bra(:,pstate) = wpol%residue_left(ipstate,:)
       enddo
       !$OMP END DO
       !$OMP END PARALLEL
     else
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,nsemin:nsemax)     = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,nsemin:nsemax,istate,ispin) )
       call xsum_auxil(bra)
     endif



     ! The application of residue theorem only retains the pole in given
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ispin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ispin)) / spin_fact


     do ipole=1,wpol%npole_reso

       if( MODULO( ipole - 1 , nproc_auxil ) /= rank_auxil ) cycle


       select case(selfenergy_approx)

       case(GW,GnW0,GnWn,ONE_RING)
         ! ymbyun 2018/07/11
         ! For now, only G0W0/GnW0/GnWn are parallelized.
         ! COLLAPSE(2) is bad for bra(:,:) in terms of memory affinity.
         ! However, it is good for G0W0 with # of threads > |nsemax - nsemin| (e.g. when only HOMO and LUMO energies are needed).
         !$OMP PARALLEL
         !$OMP DO COLLAPSE(2)
         !
         ! calculate only the diagonal !
         do pstate=nsemin,nsemax
           do iomega=-se%nomega,se%nomega
             se%sigma(iomega,pstate,ispin) = se%sigma(iomega,pstate,ispin) &
                      + bra(ipole,pstate) * bra(ipole,pstate)                                          &
                        * ( fact_full_i  / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(istate,ispin) + wpol%pole(ipole) - ieta )  &
                          + fact_empty_i / ( se%energy0(pstate,ispin) + se%omega(iomega) - energy(istate,ispin) - wpol%pole(ipole) + ieta ) )
           enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
       case(COHSEX)

         do pstate=nsemin,nsemax
           !
           ! SEX
           !
           se%sigma(0,pstate,ispin) = se%sigma(0,pstate,ispin) &
                      + bra(ipole,pstate) * bra(ipole,pstate) &
                            * fact_full_i / wpol%pole(ipole) * 2.0_dp

           !
           ! COH
           !
           se%sigma(0,pstate,ispin) = se%sigma(0,pstate,ispin) &
                      - bra(ipole,pstate) * bra(ipole,pstate) &
                            / wpol%pole(ipole)

         enddo

       case default
         call die('BUG')
       end select

     enddo !ipole

   enddo !istate
 enddo !ispin

 ! Sum up the contribution from different poles (= different procs)
 call xsum_world(se%sigma)


 write(stdout,'(a)') ' Sigma_c(omega) is calculated'



 call clean_deallocate('Temporary array',bra)

 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 endif


 call stop_clock(timing_gw_self)


end subroutine gw_selfenergy


!=========================================================================
subroutine gw_selfenergy_scalapack(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,wpol,se)
 use m_definitions
 use m_timing
 use m_warning,only: issue_warning
 use m_mpi
 use m_scalapack
 use m_inputparam
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_tools,only: coeffs_gausslegint
 use m_selfenergy_tools
 implicit none

 integer,intent(in)                  :: nstate,selfenergy_approx
 type(basis_set)                     :: basis
 real(dp),intent(in)                 :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)                 :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in)  :: wpol
 type(selfenergy_grid),intent(inout) :: se
!=====
 integer                 :: pstate,pspin
 integer                 :: iomega
 integer                 :: istate,ipole
 real(dp)                :: fact_full_i,fact_empty_i
 integer                 :: desc_wauxil(NDEL),desc_wsd(NDEL)
 integer                 :: desc_3auxil(NDEL),desc_3sd(NDEL)
 integer                 :: desc_bra(NDEL)
 integer                 :: mlocal,nlocal
 integer                 :: ilocal,jlocal,jglobal
 integer                 :: info
 real(dp),allocatable    :: eri_3tmp_auxil(:,:),eri_3tmp_sd(:,:)
 real(dp),allocatable    :: wresidue_sd(:,:)
 real(dp),allocatable    :: bra(:,:)
 complex(dp),allocatable :: sigmagw(:,:,:)
!=====

 if(.NOT. has_auxil_basis) return

#ifdef HAVE_SCALAPACK
 call start_clock(timing_gw_self)

 write(stdout,*)
 select case(selfenergy_approx)
 case(ONE_RING)
   write(stdout,*) 'Perform a one-shot one_ring calculation: SCALAPACK'
 case(GW)
   write(stdout,*) 'Perform a one-shot G0W0 calculation: SCALAPACK'
 case(GnW0)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnW0 calculation: SCALAPACK'
 case(GnWn)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnWn calculation: SCALAPACK'
 case default
   write(stdout,*) 'type:',selfenergy_approx
   call die('gw_selfenergy_scalapack: calculation type unknown')
 end select


 call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,nsemin,nsemax)



 !
 ! SCALAPACK preparation for W
 !  wpol%residue_left
 mlocal = NUMROC(nauxil_2center ,MB_auxil,iprow_auxil,first_row,nprow_auxil)
 nlocal = NUMROC(wpol%npole_reso,NB_auxil,ipcol_auxil,first_col,npcol_auxil)
 call DESCINIT(desc_wauxil,nauxil_2center,wpol%npole_reso,MB_auxil,NB_auxil,first_row,first_col,cntxt_auxil,MAX(1,mlocal),info)
 !
 ! Change data distribution
 ! from cntxt_auxil to cntxt_sd
 mlocal = NUMROC(nauxil_2center ,block_row,iprow_sd,first_row,nprow_sd)
 nlocal = NUMROC(wpol%npole_reso,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desc_wsd,nauxil_2center,wpol%npole_reso,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mlocal),info)
 call clean_allocate('TMP distributed W',wresidue_sd,mlocal,nlocal)
 call PDGEMR2D(nauxil_2center,wpol%npole_reso,wpol%residue_left,1,1,desc_wauxil, &
                                                    wresidue_sd,1,1,desc_wsd,cntxt_sd)

 ! TODO maybe I should deallocate here wpol%residue_left

 ! Temporary array sigmagw is created because OPENMP does not want to work directly with se%sigma
 allocate(sigmagw,SOURCE=se%sigma)
 sigmagw(:,:,:)  = 0.0_dp

 do pspin=1,nspin
   do pstate=nsemin,nsemax


     !
     ! SCALAPACK preparation for the 3-center integrals
     !
     mlocal = NUMROC(nauxil_2center      ,MB_auxil,iprow_auxil,first_row,nprow_auxil)
     nlocal = NUMROC(nvirtual_G-ncore_G-1,NB_auxil,ipcol_auxil,first_col,npcol_auxil)
     call DESCINIT(desc_3auxil,nauxil_2center,nvirtual_G-ncore_G-1,MB_auxil,NB_auxil,first_row,first_col,cntxt_auxil,MAX(1,mlocal),info)

     if( cntxt_auxil > 0 ) then
       call clean_allocate('TMP 3center eigen',eri_3tmp_auxil,mlocal,nlocal)
       do jlocal=1,nlocal
         jglobal = INDXL2G(jlocal,NB_auxil,ipcol_auxil,first_col,npcol_auxil) + ncore_G
         do ilocal=1,mlocal
           eri_3tmp_auxil(ilocal,jlocal) = eri_3center_eigen(ilocal,jglobal,pstate,pspin)
         enddo
       enddo
     else
       call clean_allocate('TMP 3center eigen',eri_3tmp_auxil,1,1)
     endif
     !
     ! Change data distribution
     ! from cntxt_auxil to cntxt_sd
     mlocal = NUMROC(nauxil_2center      ,block_row,iprow_sd,first_row,nprow_sd)
     nlocal = NUMROC(nvirtual_G-ncore_G-1,block_col,ipcol_sd,first_col,npcol_sd)
     call DESCINIT(desc_3sd,nauxil_2center,nvirtual_G-ncore_G-1,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mlocal),info)
     call clean_allocate('TMP 3center eigen',eri_3tmp_sd,mlocal,nlocal)
     call PDGEMR2D(nauxil_2center,nvirtual_G-ncore_G-1,eri_3tmp_auxil,1,1,desc_3auxil, &
                                                          eri_3tmp_sd,1,1,desc_3sd,cntxt_sd)
     call clean_deallocate('TMP 3center eigen',eri_3tmp_auxil)


     !
     ! Prepare a SCALAPACKed bra that is to contain  wresidue**T * v**1/2
     mlocal = NUMROC(wpol%npole_reso     ,block_row,iprow_sd,first_row,nprow_sd)
     nlocal = NUMROC(nvirtual_G-ncore_G-1,block_col,ipcol_sd,first_col,npcol_sd)
     call DESCINIT(desc_bra,wpol%npole_reso,nvirtual_G-ncore_G-1,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mlocal),info)
     call clean_allocate('Temporary array',bra,mlocal,nlocal)

     ! And calculate it
     call PDGEMM('T','N',wpol%npole_reso,nvirtual_G-ncore_G-1,nauxil_2center, &
                             1.0_dp,wresidue_sd,1,1,desc_wsd,    &
                                    eri_3tmp_sd,1,1,desc_3sd,    &
                             0.0_dp,bra        ,1,1,desc_bra)
     call clean_deallocate('TMP 3center eigen',eri_3tmp_sd)



     !$OMP PARALLEL PRIVATE(istate,ipole,fact_full_i,fact_empty_i)
     !$OMP DO REDUCTION(+:sigmagw)
     do jlocal=1,nlocal
       istate = INDXL2G(jlocal,block_col,ipcol_sd,first_col,npcol_sd) + ncore_G
       do ilocal=1,mlocal
         ipole = INDXL2G(ilocal,block_row,iprow_sd,first_row,nprow_sd)


         ! The application of residue theorem only retains the pole in given
         ! quadrants.
         ! The positive poles of W go with the poles of occupied states in G
         ! The negative poles of W go with the poles of empty states in G
         fact_full_i   = occupation(istate,pspin) / spin_fact
         fact_empty_i = (spin_fact - occupation(istate,pspin)) / spin_fact

         sigmagw(:,pstate,pspin) = sigmagw(:,pstate,pspin) &
                  + bra(ilocal,jlocal) * bra(ilocal,jlocal)                    &
                    * ( fact_full_i  / ( se%energy0(pstate,pspin) + se%omega(:) - energy(istate,pspin) + wpol%pole(ipole) - ieta )   &
                      + fact_empty_i / ( se%energy0(pstate,pspin) + se%omega(:) - energy(istate,pspin) - wpol%pole(ipole) + ieta )  )
       enddo  !ilocal -> ipole
     enddo !jlocal -> istate
     !$OMP END DO
     !$OMP END PARALLEL

     call clean_deallocate('Temporary array',bra)

   enddo !pstate
 enddo !pspin

 ! Sum up the contribution from different poles (= different procs)
 call xsum_world(sigmagw)

 se%sigma(:,:,:) = sigmagw(:,:,:)
 deallocate(sigmagw)

 write(stdout,'(a)') ' Sigma_c(omega) is calculated'

 call clean_deallocate('TMP distributed W',wresidue_sd)
 call destroy_eri_3center_eigen()

 call stop_clock(timing_gw_self)

#endif

end subroutine gw_selfenergy_scalapack


!=========================================================================
subroutine gw_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,wpol,selfenergy)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_timing
 use m_inputparam
 use m_warning,only: issue_warning
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_tools,only: coeffs_gausslegint
 use m_selfenergy_tools
 implicit none

 integer,intent(in)                 :: nstate
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
!=====
 integer               :: ipstate,pstate,qstate,istate
 integer               :: ispin,ipole
 real(dp),allocatable  :: bra(:,:)
 real(dp)              :: fact_full_i,fact_empty_i
!=====

 call start_clock(timing_gw_self)

 write(stdout,*)
 select case(calc_type%selfenergy_approx)
 case(ONE_RING)
   write(stdout,*) 'Perform a QP self-consistent one-ring calculation (QS1-ring)'
 case(GW)
   write(stdout,*) 'Perform a QP self-consistent GW calculation (QSGW)'
 case(COHSEX)
   write(stdout,*) 'Perform a self-consistent COHSEX calculation'
 case default
   call die('gw_selfenergy_qs: calculation type unknown')
 end select


 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,nsemin,nsemax,ncore_G+1,nvirtual_G-1)
 endif

 call clean_allocate('Temporary array',bra,1,wpol%npole_reso,nsemin,nsemax)


 selfenergy(:,:,:)  = 0.0_dp
 do ispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !INNER LOOP of G

     if( MODULO( istate - (ncore_G+1) , nproc_ortho) /= rank_ortho ) cycle

     !
     ! Prepare the bra and ket with the knowledge of index istate and pstate
     if( .NOT. has_auxil_basis) then
       ! Here just grab the precalculated value
       do pstate=nsemin,nsemax
         ipstate = index_prodstate(istate,pstate) + (ispin-1) * index_prodstate(nvirtual_W-1,nvirtual_W-1)
         bra(:,pstate) = wpol%residue_left(ipstate,:)
       enddo
     else
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,nsemin:nsemax)     = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,nsemin:nsemax,istate,ispin) )
       call xsum_auxil(bra)
     endif


     ! The application of residue theorem only retains the pole in given
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ispin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ispin)) / spin_fact


     do ipole=1,wpol%npole_reso

       if( MODULO( ipole - 1 , nproc_auxil ) /= rank_auxil ) cycle

       select case(calc_type%selfenergy_approx)

       case(GW)

         do qstate=nsemin,nsemax
           do pstate=nsemin,nsemax

             selfenergy(pstate,qstate,ispin) = selfenergy(pstate,qstate,ispin) &
                        + bra(ipole,pstate) * bra(ipole,qstate)                            &
                          * ( REAL(  fact_full_i  / ( energy(qstate,ispin) - ieta  - energy(istate,ispin) + wpol%pole(ipole) ) , dp ) &
                            + REAL(  fact_empty_i / ( energy(qstate,ispin) + ieta  - energy(istate,ispin) - wpol%pole(ipole) ) , dp ) )

           enddo
         enddo

       case(COHSEX)
         do qstate=nsemin,nsemax
           do pstate=nsemin,nsemax
             !
             ! SEX
             !
             selfenergy(pstate,qstate,ispin) = selfenergy(pstate,qstate,ispin) &
                        + bra(ipole,pstate) * bra(ipole,qstate)                                &
                              * fact_full_i / wpol%pole(ipole) * 2.0_dp

             !
             ! COH
             !
             selfenergy(pstate,qstate,ispin) = selfenergy(pstate,qstate,ispin) &
                        - bra(ipole,pstate) * bra(ipole,qstate) &
                              / wpol%pole(ipole)
           enddo
         enddo

       case default
         call die('BUG')
       end select

     enddo !ipole

   enddo !istate
 enddo !ispin

 ! Sum up the contribution from different poles (= different procs)
 call xsum_world(selfenergy)


 ! Kotani's hermitianization trick
 call apply_qs_approximation(s_matrix,c_matrix,selfenergy)


 call clean_deallocate('Temporary array',bra)

 if(has_auxil_basis) call destroy_eri_3center_eigen()


 call stop_clock(timing_gw_self)


end subroutine gw_selfenergy_qs


!=========================================================================
