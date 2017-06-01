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

 integer,intent(in)                 :: nstate,selfenergy_approx
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 type(selfenergy_grid),intent(inout) :: se
 real(dp),intent(out)               :: energy_gw
!=====
 integer               :: iomega
 integer               :: ipstate
 integer               :: pstate,bstate
 integer               :: istate,ispin,ipole
 real(dp),allocatable  :: bra(:,:)
 real(dp),allocatable  :: bra_exx(:,:)
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: fact_full_a,fact_empty_a
 real(dp)              :: energy_lw(nstate,nspin)
 character(len=3)      :: ctmp
 integer               :: selfenergyfile
! GW tilde
 real(dp),allocatable  :: vsqchi0vsqm1(:,:)
 real(dp)              :: omega_m_ei,bra2
! LW devel
 complex(dp),allocatable :: omegac(:)
 complex(dp),allocatable :: selfenergy_omegac(:,:,:,:)
 complex(dp),allocatable :: matrix(:,:),eigvec(:,:)
 real(dp),allocatable    :: eigval(:),x1(:),weights(:)
 real(dp)                :: tr_log_gsigma,tr_gsigma,rdiag,mu
 real(dp),allocatable    :: c_matrix_exx(:,:,:)
!=====

 call start_clock(timing_self)

 write(stdout,*)
 select case(selfenergy_approx)
 case(LW,LW2)
   write(stdout,*) 'Perform the calculation of Tr[ ln ( 1 - GSigma ) ]'
 case(GSIGMA)
   write(stdout,*) 'Perform the calculation of Tr[ SigmaG ]'
 case(GV)
   write(stdout,*) 'Perform a perturbative HF calculation'
 case(GW)
   write(stdout,*) 'Perform a one-shot G0W0 calculation'
 case(ONE_RING)
   write(stdout,*) 'Perform a one-shot one-ring calculation'
 case(COHSEX)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case(GnW0)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnW0 calculation'
 case(GnWn)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnWn calculation'
 case default
   write(stdout,*) 'type:',selfenergy_approx
   call die('gw_selfenergy: calculation type unknown')
 end select


 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,nsemin,nsemax,ncore_G+1,nvirtual_G-1)
   if( calc_type%selfenergy_approx == LW .OR. calc_type%selfenergy_approx == LW2 .OR. calc_type%selfenergy_approx == GSIGMA ) &
       call calculate_eri_3center_eigen_mixed(basis%nbf,nstate,c_matrix)
 endif


 call clean_allocate('Temporary array',bra,1,wpol%npole_reso,nsemin,nsemax)

 if(selfenergy_approx==LW .OR. selfenergy_approx==LW2 .OR. selfenergy_approx==GSIGMA) then
   call clean_allocate('Temporary array for LW',bra_exx,1,wpol%npole_reso,nsemin,nsemax)
 endif


 energy_gw = 0.0_dp


 !
 ! The ones with imaginary frequencies
 select case(selfenergy_approx)
 case(LW,LW2)
   allocate(omegac(1:se%nomega))
   
   allocate(x1(se%nomega),weights(se%nomega))  
   call coeffs_gausslegint(0.0_dp,1.0_dp,x1,weights,se%nomega)

   mu =-0.10_dp
   write(stdout,*) 'mu is set manually here to',mu
   do iomega=1,se%nomega
     omegac(iomega) = x1(iomega) / ( 1.0_dp - x1(iomega) ) * im + mu
     weights(iomega) = weights(iomega) / (1.0_dp - x1(iomega))**2
     write(stdout,'(i4,3(2x,f12.4))') iomega,REAL(omegac(iomega)),AIMAG(omegac(iomega)),weights(iomega)
   enddo
   deallocate(x1)
 end select


 if( selfenergy_approx == LW .OR. selfenergy_approx == LW2 .OR. selfenergy_approx == GSIGMA ) then
   call issue_warning('reading G\tilde')
   open(1001,form='unformatted')
   read(1001) energy_lw(:,:)
   close(1001,status='delete')
 endif


 !
 ! Which calculation type needs a complex sigma?
 !
 select case(selfenergy_approx)
 case(LW,LW2)                     ! matrix complex
   allocate(selfenergy_omegac(1:se%nomega,nsemin:nsemax,nsemin:nsemax,nspin))
   selfenergy_omegac(:,:,:,:) = 0.0_dp
 end select


 se%sigma(:,:,:)  = 0.0_dp

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
       if( selfenergy_approx==LW .OR. selfenergy_approx==LW2 .OR. selfenergy_approx==GSIGMA) then
         bra_exx(:,nsemin:nsemax) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen_mixed(:,istate,nsemin:nsemax,ispin) )
         call xsum_auxil(bra_exx)
       endif
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

       case(LW,LW2)

         do bstate=nsemin,nsemax
           do pstate=nsemin,nsemax

             do iomega=1,se%nomega
               selfenergy_omegac(iomega,pstate,bstate,ispin) = selfenergy_omegac(iomega,pstate,bstate,ispin) &
                        + bra_exx(ipole,pstate) * bra_exx(ipole,bstate) &
                          * (  fact_full_i  / ( omegac(iomega) - energy(istate,ispin) + wpol%pole(ipole)  )    &
                             + fact_empty_i / ( omegac(iomega) - energy(istate,ispin) - wpol%pole(ipole) ) )   &
                          / ( omegac(iomega) - energy_lw(pstate,ispin) )
             enddo

           enddo
         enddo

       case(GSIGMA)

         !
         ! calculate only the diagonal !
         do pstate=nsemin,nsemax
           fact_full_a   = occupation(pstate,ispin) / spin_fact
           fact_empty_a  = (spin_fact - occupation(pstate,ispin)) / spin_fact
           se%sigma(0,pstate,ispin) = se%sigma(0,pstate,ispin) &
                    - bra(ipole,pstate) * bra(ipole,pstate) &
                      * ( fact_full_i  * fact_empty_a / ( energy(pstate,ispin)  - energy(istate,ispin) + wpol%pole(ipole) - ieta )   &
                        - fact_empty_i * fact_full_a  / ( energy(pstate,ispin)  - energy(istate,ispin) - wpol%pole(ipole) + ieta )   )
           se%sigma(1,pstate,ispin) = se%sigma(1,pstate,ispin) &
                    - bra_exx(ipole,pstate) * bra_exx(ipole,pstate) &
                      * ( fact_full_i  * fact_empty_a / ( energy_lw(pstate,ispin)  - energy(istate,ispin) + wpol%pole(ipole) - ieta )   &
                        - fact_empty_i * fact_full_a  / ( energy_lw(pstate,ispin)  - energy(istate,ispin) - wpol%pole(ipole) + ieta )  )
         enddo


       case(GW,GnW0,GnWn,ONE_RING)

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

       case(COHSEX)

         do pstate=nsemin,nsemax
           !
           ! SEX
           !
           se%sigma(0,pstate,ispin) = se%sigma(0,pstate,ispin) &
                      + bra(ipole,pstate) * bra(ipole,pstate) &
                            * fact_full_i / wpol%pole(ipole) * 2.0_dp  &
                            * beta_cohsex

           !
           ! COH
           !
           se%sigma(0,pstate,ispin) = se%sigma(0,pstate,ispin) &
                      - bra(ipole,pstate) * bra(ipole,pstate) &
                            / wpol%pole(ipole)                &
                            * alpha_cohsex

         enddo

       case(GV)
         !
         ! Do nothing: no correlation in this case
         !
       case default 
         call die('BUG')
       end select

     enddo !ipole

   enddo !istate
 enddo !ispin

 ! Sum up the contribution from different poles (= different procs)
 call xsum_world(se%sigma)

 if( ALLOCATED(selfenergy_omegac) ) then
   call xsum_world(selfenergy_omegac)
 endif


 write(stdout,'(a)') ' Sigma_c(omega) is calculated'


 !
 ! Only the EXPERIMENTAL features need to do some post-processing here!
 select case(selfenergy_approx)
 case(GSIGMA) !==========================================================

   energy_gw = 0.5_dp * SUM(REAL(se%sigma(1,:,:),dp)) * spin_fact
   write(stdout,*) 'Tr[Sigma tilde G]:',2.0_dp*energy_gw

   energy_gw = 0.5_dp * SUM(REAL(se%sigma(0,:,:))) * spin_fact
   write(stdout,*) '       Tr[SigmaG]:',2.0_dp*energy_gw

 case(LW)

   allocate(matrix(nsemin:nsemax,nsemin:nsemax))
   allocate(eigvec(nsemin:nsemax,nsemin:nsemax))
   allocate(eigval(nsemax-nsemin+1))

   tr_log_gsigma = 0.0_dp
   tr_gsigma     = 0.0_dp

   do ispin=1,nspin
     do iomega=1,se%nomega

       rdiag = 0.d0
       do istate=nsemin,nsemax
         rdiag = rdiag + REAL(selfenergy_omegac(iomega,istate,istate,ispin),dp) * 2.0_dp
       enddo

       matrix(:,:) = selfenergy_omegac(iomega,:,:,ispin) + CONJG(TRANSPOSE( selfenergy_omegac(iomega,:,:,ispin) )) &
                    - MATMUL( selfenergy_omegac(iomega,:,:,ispin) , CONJG(TRANSPOSE( selfenergy_omegac(iomega,:,:,ispin) )) )

       call diagonalize(nsemax-nsemin+1,matrix,eigval,eigvec)

       tr_gsigma     = tr_gsigma     + rdiag                         * spin_fact / (2.0 * pi) * weights(iomega)
       tr_log_gsigma = tr_log_gsigma + SUM(LOG( 1.0_dp - eigval(:))) * spin_fact / (2.0 * pi) * weights(iomega)

     enddo
   enddo


   write(stdout,*) 'Tr[Log(1-GSigma)] ',tr_log_gsigma 
   write(stdout,*) 'Tr[GSigma]        ',tr_gsigma 
   write(stdout,*) 'Sum               ',(tr_log_gsigma+tr_gsigma) 
   deallocate(matrix,eigvec,eigval)

 case(LW2)

   allocate(matrix(basis%nbf,basis%nbf))
   allocate(eigvec(basis%nbf,basis%nbf))
   allocate(eigval(basis%nbf))

   tr_log_gsigma = 0.0_dp
   tr_gsigma     = 0.0_dp

   do iomega=1,se%nomega

     matrix(:,:) = MATMUL( selfenergy_omegac(iomega,:,:,1) , selfenergy_omegac(iomega,:,:,1) )   &
                + MATMUL( TRANSPOSE(CONJG(selfenergy_omegac(iomega,:,:,1))) , TRANSPOSE(CONJG(selfenergy_omegac(iomega,:,:,1))) )

     rdiag=0.d0
     do istate=1,basis%nbf
       rdiag = rdiag - REAL(matrix(istate,istate),dp) * spin_fact / (2.0 * pi)   * 0.5_dp  ! -1/2 comes from the log expansion
     enddo


     tr_gsigma = tr_gsigma + rdiag * weights(iomega)

   enddo


   write(stdout,*) 'Tr[tGSigma*tGSigma]        ',tr_gsigma 
   deallocate(matrix,eigvec,eigval)

 end select



 call clean_deallocate('Temporary array',bra)
 if(ALLOCATED(bra_exx)) call clean_deallocate('Temporary array for LW',bra_exx)

 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
   if( calc_type%selfenergy_approx == LW .OR. calc_type%selfenergy_approx == LW2 .OR. calc_type%selfenergy_approx == GSIGMA ) &
       call calculate_eri_3center_eigen_mixed(basis%nbf,nstate,c_matrix)
 endif

 if(ALLOCATED(omegac)) deallocate(omegac)
 if(ALLOCATED(weights)) deallocate(weights)
 if(ALLOCATED(selfenergy_omegac)) deallocate(selfenergy_omegac)


 call stop_clock(timing_self)


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
 integer               :: pstate,pspin
 integer               :: iomega
 integer               :: istate,ipole
 real(dp)              :: fact_full_i,fact_empty_i
 integer               :: selfenergyfile
 integer               :: desc_wauxil(NDEL),desc_wsd(NDEL)
 integer               :: desc_3auxil(NDEL),desc_3sd(NDEL)
 integer               :: desc_bra(NDEL)
 integer               :: mlocal,nlocal
 integer               :: ilocal,jlocal,iglobal,jglobal
 integer               :: info
 real(dp),allocatable  :: eri_3tmp_auxil(:,:),eri_3tmp_sd(:,:)
 real(dp),allocatable  :: wresidue_sd(:,:)
 real(dp),allocatable  :: bra(:,:)
!=====

 if(.NOT. has_auxil_basis) return

#ifdef HAVE_SCALAPACK
 call start_clock(timing_self)

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


 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore_G+1,nvirtual_G-1,nsemin,nsemax)



 !
 ! SCALAPACK preparation for W
 !  wpol%residue_left
 mlocal = NUMROC(nauxil_2center ,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
 nlocal = NUMROC(wpol%npole_reso,NBLOCK_AUXIL,ipcol_auxil,first_col,npcol_auxil)
 call DESCINIT(desc_wauxil,nauxil_2center,wpol%npole_reso,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,mlocal),info)
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
 

 se%sigma(:,:,:)  = 0.0_dp

 do pspin=1,nspin
   do pstate=nsemin,nsemax


     !
     ! SCALAPACK preparation for the 3-center integrals
     !
     mlocal = NUMROC(nauxil_2center      ,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
     nlocal = NUMROC(nvirtual_G-ncore_G-1,NBLOCK_AUXIL,ipcol_auxil,first_col,npcol_auxil)
     call DESCINIT(desc_3auxil,nauxil_2center,nvirtual_G-ncore_G-1,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,mlocal),info)

     if( cntxt_auxil > 0 ) then
       call clean_allocate('TMP 3center eigen',eri_3tmp_auxil,mlocal,nlocal)
       do jlocal=1,nlocal
         jglobal = INDXL2G(jlocal,NBLOCK_AUXIL,ipcol_auxil,first_col,npcol_auxil) + ncore_G
         do ilocal=1,mlocal
           eri_3tmp_auxil(ilocal,jlocal) = eri_3center_eigen(ilocal,jglobal,pstate,pspin)
         enddo
       enddo
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

         do iomega=-se%nomega,se%nomega
           se%sigma(iomega,pstate,pspin) = se%sigma(iomega,pstate,pspin) &
                    + bra(ilocal,jlocal) * bra(ilocal,jlocal)                    & 
                      * ( fact_full_i  / ( se%energy0(pstate,pspin) + se%omega(iomega) - energy(istate,pspin) + wpol%pole(ipole) - ieta )   &
                        + fact_empty_i / ( se%energy0(pstate,pspin) + se%omega(iomega) - energy(istate,pspin) - wpol%pole(ipole) + ieta )  )
         enddo !iomega
       enddo  !ilocal -> ipole
     enddo !jlocal -> istate

     call clean_deallocate('Temporary array',bra)

   enddo !pstate
 enddo !pspin

 ! Sum up the contribution from different poles (= different procs)
 call xsum_world(se%sigma)


 write(stdout,'(a)') ' Sigma_c(omega) is calculated'

 call clean_deallocate('TMP distributed W',wresidue_sd)
 call destroy_eri_3center_eigen()

 call stop_clock(timing_self)

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

 call start_clock(timing_self)

 write(stdout,*)
 select case(calc_type%selfenergy_approx)
 case(ONE_RING)
   write(stdout,*) 'Perform a QP self-consistent one-ring calculation (QS1-ring)'
 case(GW)
   write(stdout,*) 'Perform a QP self-consistent GW calculation (QSGW)'
 case(COHSEX)
   write(stdout,*) 'Perform a self-consistent COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case default
   call die('gw_selfenergy_qs: calculation type unknown')
 end select


 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,nsemin,nsemax,ncore_G+1,nvirtual_G-1)
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
                              * fact_full_i / wpol%pole(ipole) * 2.0_dp                        &
                              * beta_cohsex

             !
             ! COH
             !
             selfenergy(pstate,qstate,ispin) = selfenergy(pstate,qstate,ispin) &
                        - bra(ipole,pstate) * bra(ipole,qstate) & 
                              / wpol%pole(ipole)                &
                              * alpha_cohsex
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
 call apply_qs_approximation(basis%nbf,nstate,s_matrix,c_matrix,selfenergy)


 call clean_deallocate('Temporary array',bra)

 if(has_auxil_basis) call destroy_eri_3center_eigen()


 call stop_clock(timing_self)


end subroutine gw_selfenergy_qs


!=========================================================================
