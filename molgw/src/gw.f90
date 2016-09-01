!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the GW self-energy
! within different flavors: G0W0, GnW0, GnWn, COHSEX
!
!=========================================================================
subroutine gw_selfenergy(nstate,gwmethod,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy,energy_gw)
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam
 use m_warning,only: issue_warning,msg
 use m_basis_set
 use m_spectral_function
 use m_eri_ao_mo
 use m_tools,only: coeffs_gausslegint
 use m_selfenergy_tools
 implicit none

 integer,intent(in)                 :: nstate,gwmethod
 type(basis_set)                    :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)               :: energy_gw
!=====
 integer               :: nprodbasis
 integer               :: iomegai
 complex(dp),allocatable :: omegac(:)
 complex(dp),allocatable :: selfenergy_omegac(:,:,:,:)
 real(dp),allocatable  :: selfenergy_omega(:,:,:,:)
 real(dp),allocatable  :: zz(:,:)
 integer               :: iastate
 integer               :: astate,bstate
 integer               :: istate,ispin,ipole
 real(dp),allocatable  :: bra(:,:)
 real(dp),allocatable  :: bra_exx(:,:)
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: fact_full_a,fact_empty_a
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin),energy_qp_z(nstate,nspin)
 character(len=3)      :: ctmp
 integer               :: reading_status
 integer               :: selfenergyfile
! GW tilde
 real(dp),allocatable  :: vsqchi0vsqm1(:,:)
 real(dp)              :: omega_m_ei,bra2
! LW devel
 complex(dp),allocatable :: matrix(:,:),eigvec(:,:)
 real(dp),allocatable    :: eigval(:),x1(:),weights(:)
 real(dp)                :: tr_log_gsigma,tr_gsigma,rdiag,mu
 real(dp),allocatable    :: c_matrix_exx(:,:,:)
!=====

 call start_clock(timing_self)

 write(stdout,*)
 select case(gwmethod)
 case(LW,LW2)
   write(stdout,*) 'Perform the calculation of Tr[ ln ( 1 - GSigma ) ]'
 case(GSIGMA)
   write(stdout,*) 'Perform the calculation of Tr[ SigmaG ]'
 case(GSIGMA3)
   write(stdout,*) 'Perform the calculation of Tr[ SigmaG ] numerical for test only'
 case(GV)
   write(stdout,*) 'Perform a perturbative HF calculation'
 case(QS)
   write(stdout,*) 'Perform a QP self-consistent GW calculation (QSGW)'
 case(GWTILDE)
   write(stdout,*) 'Perform a one-shot GWtilde calculation'
 case(G0W0)
   write(stdout,*) 'Perform a one-shot G0W0 calculation'
 case(G0W0_IOMEGA)
   write(stdout,*) 'Perform a one-shot G0W0 calculation'
 case(COHSEX)
   write(stdout,*) 'Perform a COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case(QSCOHSEX)
   write(stdout,*) 'Perform a self-consistent COHSEX calculation'
   if( ABS(alpha_cohsex - 1.0_dp) > 1.0e-4_dp .OR. ABS(beta_cohsex - 1.0_dp) > 1.0e-4_dp ) then
     write(stdout,'(a,2(2x,f12.6))') ' Tuned COHSEX with parameters alpha, beta: ',alpha_cohsex,beta_cohsex
   endif
 case(GnW0)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnW0 calculation'
 case(GnWn)
   write(stdout,*) 'Perform an eigenvalue self-consistent GnWn calculation'
 end select


 nprodbasis = index_prodstate(nstate,nstate)

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,nsemin,nsemax,ncore_G+1,nvirtual_G-1)
   if( calc_type%gwmethod == LW .OR. calc_type%gwmethod == LW2 .OR. calc_type%gwmethod == GSIGMA ) &
       call calculate_eri_3center_eigen_mixed(basis%nbf,nstate,c_matrix)
 endif

 call clean_allocate('Temporary array',bra,1,wpol%npole_reso,nsemin,nsemax)

 if(gwmethod==LW .OR. gwmethod==LW2 .OR. gwmethod==GSIGMA) then
   call clean_allocate('Temporary array for LW',bra_exx,1,wpol%npole_reso,nsemin,nsemax)
 endif


 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 call issue_warning('small complex number is '//msg)


 !
 ! The ones with imaginary frequencies
 select case(gwmethod)
 case(LW,LW2,G0W0_IOMEGA)
   allocate(omegac(1:nomegai))
   
   allocate(x1(nomegai),weights(nomegai))  
   call coeffs_gausslegint(0.0_dp,1.0_dp,x1,weights,nomegai)

   mu =-0.10_dp
   write(stdout,*) 'mu is set manually here to',mu
   do iomegai=1,nomegai
     omegac(iomegai) = x1(iomegai) / ( 1.0_dp - x1(iomegai) ) * im + mu
     weights(iomegai) = weights(iomegai) / (1.0_dp - x1(iomegai))**2
     write(stdout,'(i4,3(2x,f12.4))') iomegai,REAL(omegac(iomegai)),AIMAG(omegac(iomegai)),weights(iomegai)
   enddo
   deallocate(x1)
 end select

 !
 ! Which calculation type needs to update energy_qp
 !
 if( calc_type%read_energy_qp ) then
   call read_energy_qp(nstate,energy_qp,reading_status)
   if(reading_status/=0) then
     call issue_warning('File energy_qp not found: assuming 1st iteration')
     energy_qp(:,:) = energy(:,:)
   endif

 else
   energy_qp(:,:) = energy(:,:)

 endif

 if( gwmethod == LW .OR. gwmethod == LW2 .OR. gwmethod == GSIGMA ) then
   call issue_warning('reading G\tilde')
   open(1001,form='unformatted')
   read(1001) energy_qp(:,:)
   close(1001,status='delete')
 endif


 !
 ! Which calculation type needs the diagonal only?
 ! Which calculation type needs a complex sigma?
 !
 select case(gwmethod)
 case(QS,QSCOHSEX,GSIGMA3) ! matrix real
   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nsemin:nsemax,nspin))

 case(LW,LW2,G0W0_IOMEGA)                     ! matrix complex
   allocate(selfenergy_omegac(1:nomegai,nsemin:nsemax,nsemin:nsemax,nspin))

 case default
   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,1,nspin))

 end select

 if( ALLOCATED(selfenergy_omega) )  selfenergy_omega(:,:,:,:)  = 0.0_dp
 if( ALLOCATED(selfenergy_omegac) ) selfenergy_omegac(:,:,:,:) = 0.0_dp



 do ispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !INNER LOOP of G

     if( MODULO( istate -(ncore_G+1) , nproc_ortho) /= rank_ortho ) cycle

     !
     ! Prepare the bra and ket with the knowledge of index istate and astate
     if( .NOT. has_auxil_basis) then
       ! Here just grab the precalculated value
       do astate=nsemin,nsemax
         iastate = index_prodstate(istate,astate) + nprodbasis * (ispin-1)
         bra(:,astate) = wpol%residu_left(iastate,:)
       end do
     else
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,nsemin:nsemax)     = MATMUL( TRANSPOSE(wpol%residu_left(:,:)) , eri_3center_eigen(:,nsemin:nsemax,istate,ispin) )
       call xsum_auxil(bra)
       if( gwmethod==LW .OR. gwmethod==LW2 .OR. gwmethod==GSIGMA) then
         bra_exx(:,nsemin:nsemax) = MATMUL( TRANSPOSE(wpol%residu_left(:,:)) , eri_3center_eigen_mixed(:,istate,nsemin:nsemax,ispin) )
         call xsum_auxil(bra_exx)
       endif
     endif


     ! The application of residu theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ispin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ispin)) / spin_fact


     do ipole=1,wpol%npole_reso


       select case(gwmethod)

       case(QS)

         do bstate=nsemin,nsemax
           do astate=nsemin,nsemax

             selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
                        + bra(ipole,astate) * bra(ipole,bstate)                            &  
                          * ( REAL(  fact_full_i  / ( energy_qp(bstate,ispin) - ieta  - energy_qp(istate,ispin) + wpol%pole(ipole) ) , dp ) &
                            + REAL(  fact_empty_i / ( energy_qp(bstate,ispin) + ieta  - energy_qp(istate,ispin) - wpol%pole(ipole) ) , dp ) )

           enddo
         enddo

       case(LW,LW2)

         if(ipole==1) write(stdout,*) 'FBFB LW',istate

         do bstate=nsemin,nsemax
           do astate=nsemin,nsemax

             do iomegai=1,nomegai
               selfenergy_omegac(iomegai,astate,bstate,ispin) = selfenergy_omegac(iomegai,astate,bstate,ispin) &
                        + bra_exx(ipole,astate) * bra_exx(ipole,bstate) &
                          * (  fact_full_i  / ( omegac(iomegai) - energy(istate,ispin) + wpol%pole(ipole)  )    &
                             + fact_empty_i / ( omegac(iomegai) - energy(istate,ispin) - wpol%pole(ipole) ) )   &
                          / ( omegac(iomegai) - energy_qp(astate,ispin) )
             enddo

           enddo
         enddo

       case(GSIGMA3)

         do bstate=nsemin,nsemax
           do astate=nsemin,nsemax

             fact_full_a   = occupation(astate,ispin) / spin_fact
             fact_empty_a  = (spin_fact - occupation(astate,ispin)) / spin_fact
             selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
                      - bra(ipole,astate) * bra(ipole,bstate) &
                        * ( REAL(  fact_full_i  * fact_empty_a / ( energy_qp(astate,ispin)  - energy(istate,ispin) + wpol%pole(ipole) - ieta )  , dp )  &
                          - REAL(  fact_empty_i * fact_full_a  / ( energy_qp(astate,ispin)  - energy(istate,ispin) - wpol%pole(ipole) + ieta )  , dp ) )

           enddo
         enddo


       case(GSIGMA)

         do astate=nsemin,nsemax
           fact_full_a   = occupation(astate,ispin) / spin_fact
           fact_empty_a  = (spin_fact - occupation(astate,ispin)) / spin_fact
           !
           ! calculate only the diagonal !
           selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                    - bra(ipole,astate) * bra(ipole,astate) &
                      * ( REAL(  fact_full_i  * fact_empty_a / ( energy(astate,ispin)  - energy(istate,ispin) + wpol%pole(ipole) - ieta )  , dp )  &
                        - REAL(  fact_empty_i * fact_full_a  / ( energy(astate,ispin)  - energy(istate,ispin) - wpol%pole(ipole) + ieta )  , dp ) )
           selfenergy_omega(1,astate,1,ispin) = selfenergy_omega(1,astate,1,ispin) &
                    - bra_exx(ipole,astate) * bra_exx(ipole,astate) &
                      * ( REAL(  fact_full_i  * fact_empty_a / ( energy_qp(astate,ispin)  - energy(istate,ispin) + wpol%pole(ipole) - ieta )  , dp )  &
                        - REAL(  fact_empty_i * fact_full_a  / ( energy_qp(astate,ispin)  - energy(istate,ispin) - wpol%pole(ipole) + ieta )  , dp ) )
         enddo


       case(G0W0,GnW0,GnWn)

         do astate=nsemin,nsemax
           !
           ! calculate only the diagonal !
           do iomegai=-nomegai,nomegai
             selfenergy_omega(iomegai,astate,1,ispin) = selfenergy_omega(iomegai,astate,1,ispin) &
                      + bra(ipole,astate) * bra(ipole,astate)                                          & 
                        * ( REAL(  fact_full_i  / ( energy_qp(astate,ispin) + omegai(iomegai) - energy_qp(istate,ispin) + wpol%pole(ipole) - ieta )  , dp )  &
                          + REAL(  fact_empty_i / ( energy_qp(astate,ispin) + omegai(iomegai) - energy_qp(istate,ispin) - wpol%pole(ipole) + ieta )  , dp ) )
           enddo
         enddo

       case(GWTILDE)

         allocate(vsqchi0vsqm1(nauxil_2center,nauxil_2center))

         do astate=nsemin,nsemax
           do iomegai=-nomegai,nomegai
             omega_m_ei = energy_qp(astate,ispin) + omegai(iomegai) - energy_qp(istate,ispin)
             call dynamical_polarizability(nstate,occupation,energy_qp,omega_m_ei,wpol,vsqchi0vsqm1)

             bra2 = DOT_PRODUCT( eri_3center_eigen(:,astate,istate,ispin) , MATMUL( vsqchi0vsqm1(:,:) , wpol%residu_left(:,ipole)) )

             selfenergy_omega(iomegai,astate,1,ispin) = selfenergy_omega(iomegai,astate,1,ispin) &
                      + bra2 * bra(ipole,astate)                                          & 
                        * ( REAL(  fact_full_i  / ( omega_m_ei + wpol%pole(ipole) - ieta )  , dp )  &
                          + REAL(  fact_empty_i / ( omega_m_ei - wpol%pole(ipole) + ieta )  , dp ) )
           enddo
         enddo

         deallocate(vsqchi0vsqm1)

       case(G0W0_IOMEGA)

         do astate=nsemin,nsemax
           !
           ! calculate only the diagonal !
             do iomegai=1,nomegai
               selfenergy_omegac(iomegai,astate,1,ispin) = selfenergy_omegac(iomegai,astate,1,ispin) &
                      + bra(ipole,astate) * bra(ipole,astate)                                          & 
                        * ( fact_full_i  / ( omegac(iomegai) - energy_qp(istate,ispin) + wpol%pole(ipole) ) &  !  - ieta )    &
                          + fact_empty_i / ( omegac(iomegai) - energy_qp(istate,ispin) - wpol%pole(ipole) ) )  ! + ieta )  )
           enddo
         enddo


       case(COHSEX)

         do astate=nsemin,nsemax
           !
           ! SEX
           !
           selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                      + bra(ipole,astate) * bra(ipole,astate) &
                            * fact_full_i / wpol%pole(ipole) * 2.0_dp  &
                            * beta_cohsex

           !
           ! COH
           !
           selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                      - bra(ipole,astate) * bra(ipole,astate) &
                            / wpol%pole(ipole)                &
                            * alpha_cohsex

         enddo

       case(QSCOHSEX) 
         do bstate=nsemin,nsemax
           do astate=nsemin,nsemax
             !
             ! SEX
             !
             selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
                        + bra(ipole,astate) * bra(ipole,bstate)                                & 
                              * fact_full_i / wpol%pole(ipole) * 2.0_dp                        &
                              * beta_cohsex

             !
             ! COH
             !
             selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
                        - bra(ipole,astate) * bra(ipole,bstate) & 
                              / wpol%pole(ipole)                &
                              * alpha_cohsex
           enddo
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
 if( ALLOCATED(selfenergy_omega) ) then
   call xsum_ortho(selfenergy_omega)
 endif
 if( ALLOCATED(selfenergy_omegac) ) then
   call xsum_ortho(selfenergy_omegac)
 endif


 write(stdout,'(a)') ' Sigma_c(omega) is calculated'

 !
 ! Kotani's Hermitianization trick
 !
 if(gwmethod==QS) then
   do ispin=1,nspin
     selfenergy_omega(0,:,:,ispin) = 0.5_dp * ( selfenergy_omega(0,:,:,ispin) + TRANSPOSE(selfenergy_omega(0,:,:,ispin)) )
   enddo
 endif


 select case(gwmethod)
 case(QS,QSCOHSEX)
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )
   enddo

 case(GV) !==========================================================

   selfenergy(:,:,:) = selfenergy_omega(0,:,:,:)
   
   call find_qp_energy_linearization(selfenergy_omega(:,:,1,:),nstate,exchange_m_vxc_diag,energy_qp,energy_qp_new)
   call output_qp_energy('GV',nstate,energy_qp,exchange_m_vxc_diag,1,selfenergy_omega(0,:,1,:),energy_qp_new)

   call write_energy_qp(nstate,energy_qp_new)


 case(COHSEX) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall
   
   call find_qp_energy_linearization(selfenergy_omega(:,:,1,:),nstate,exchange_m_vxc_diag,energy_qp,energy_qp_new)
   call output_qp_energy('COHSEX',nstate,energy_qp,exchange_m_vxc_diag,1,selfenergy_omega(0,:,1,:),energy_qp_new)

   call write_energy_qp(nstate,energy_qp_new)


 case(G0W0_IOMEGA) !==========================================================

   if( is_iomaster ) then

     do astate=nsemin,nsemax
       write(ctmp,'(i3.3)') astate
       write(stdout,'(1x,a,a)') 'Printing file: ','selfenergy_omegac_state'//TRIM(ctmp)
       open(newunit=selfenergyfile,file='selfenergy_omegac_state'//TRIM(ctmp))
       write(selfenergyfile,'(a)') '# Re omega (eV)  Im omega (eV)  SigmaC (eV) '
       do iomegai=1,nomegai
         write(selfenergyfile,'(20(f12.6,2x))') omegac(iomegai)*Ha_eV,                                     &
                                            selfenergy_omegac(iomegai,astate,1,:)*Ha_eV
       enddo
       write(selfenergyfile,*)
       close(selfenergyfile)
     enddo

   endif



 case(G0W0,GWTILDE) !==========================================================

   if( calc_type%gwmethod == G0W0GAMMA0 .OR. calc_type%gwmethod == G0W0SOX0 ) then
     if( is_iomaster ) then
       open(newunit=selfenergyfile,file='g0w0.dat',form='unformatted')
       do ispin=1,nspin
         do astate=nsemin,nsemax
           write(selfenergyfile) selfenergy_omega(:,astate,1,ispin)
         enddo
       enddo
       close(selfenergyfile)
     endif
   endif

   if( print_sigma_) then
     call write_selfenergy_omega('selfenergy_gw',nstate,energy_qp,exchange_m_vxc_diag,selfenergy_omega(:,:,1,:))
   endif


   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall

   allocate(zz(nsemin:nsemax,nspin))
   call find_qp_energy_linearization(selfenergy_omega(:,:,1,:),nstate,exchange_m_vxc_diag,energy,energy_qp_z,zz)
   call find_qp_energy_graphical(selfenergy_omega(:,:,1,:),nstate,exchange_m_vxc_diag,energy,energy_qp_new)

   call output_qp_energy('G0W0',nstate,energy_qp,exchange_m_vxc_diag,1,selfenergy_omega(0,:,1,:),energy_qp_z,energy_qp_new,zz)
   deallocate(zz)

   if( .NOT. (calc_type%gwmethod == G0W0GAMMA0 .OR. calc_type%gwmethod == G0W0SOX0) )  then
     call write_energy_qp(nstate,energy_qp_new)
   endif



 case(GnW0,GnWn) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall

   call find_qp_energy_linearization(selfenergy_omega(:,:,1,:),nstate,exchange_m_vxc_diag,energy,energy_qp_new)

   if( gwmethod == GnW0 ) then
     call output_qp_energy('GnW0',nstate,energy,exchange_m_vxc_diag,1,selfenergy_omega(0,:,1,:),energy_qp_new)
   else
     call output_qp_energy('GnWn',nstate,energy,exchange_m_vxc_diag,1,selfenergy_omega(0,:,1,:),energy_qp_new)
   endif

   call write_energy_qp(nstate,energy_qp_new)

 case(GSIGMA) !==========================================================

   energy_gw = 0.5_dp * SUM(selfenergy_omega(1,:,1,:)) * spin_fact
   write(stdout,*) 'Tr[Sigma tilde G]:',2.0_dp*energy_gw

   energy_gw = 0.5_dp * SUM(selfenergy_omega(0,:,1,:)) * spin_fact
   write(stdout,*) '       Tr[SigmaG]:',2.0_dp*energy_gw

 case(GSIGMA3) !==========================================================

   energy_gw = 0.0_dp
   do astate=nsemin,nsemax
     energy_gw = energy_gw + 0.5_dp * SUM(selfenergy_omega(0,astate,astate,:)) * spin_fact
   enddo
   write(stdout,*) 'Tr[SigmaG]:',2.0_dp*energy_gw


   allocate(c_matrix_exx(basis%nbf,nsemin:nsemax,nspin))
   open(1000,form='unformatted')
   do ispin=1,nspin
     do astate=nsemin,nsemax
       read(1000) c_matrix_exx(:,astate,ispin)
     enddo
   enddo
   close(1000)
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( TRANSPOSE(c_matrix_exx(:,:,ispin)) , MATMUL( selfenergy(:,:,ispin), c_matrix_exx(:,:,ispin) ) )
!     selfenergy(:,:,ispin) = MATMUL( TRANSPOSE(c_matrix(:,:,ispin)) , MATMUL( selfenergy(:,:,ispin), c_matrix(:,:,ispin) ) )
   enddo
   deallocate(c_matrix_exx)

   energy_gw = 0.0_dp
   do astate=nsemin,nsemax
     energy_gw = energy_gw + 0.5_dp * SUM(selfenergy(astate,astate,:)) * spin_fact
   enddo
   write(stdout,*) 'Tr[SigmaG]:',2.0_dp*energy_gw


 case(LW)

   allocate(matrix(nsemin:nsemax,nsemin:nsemax))
   allocate(eigvec(nsemin:nsemax,nsemin:nsemax))
   allocate(eigval(nsemax-nsemin+1))

   tr_log_gsigma = 0.0_dp
   tr_gsigma     = 0.0_dp

   do ispin=1,nspin
     do iomegai=1,nomegai

       rdiag = 0.d0
       do istate=nsemin,nsemax
         rdiag = rdiag + REAL(selfenergy_omegac(iomegai,istate,istate,ispin),dp) * 2.0_dp
       enddo

       matrix(:,:) = selfenergy_omegac(iomegai,:,:,ispin) + CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,ispin) )) &
                    - MATMUL( selfenergy_omegac(iomegai,:,:,ispin) , CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,ispin) )) )

       call diagonalize(nsemax-nsemin+1,matrix,eigval,eigvec)

       tr_gsigma     = tr_gsigma     + rdiag                         * spin_fact / (2.0 * pi) * weights(iomegai)
       tr_log_gsigma = tr_log_gsigma + SUM(LOG( 1.0_dp - eigval(:))) * spin_fact / (2.0 * pi) * weights(iomegai)

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

   do iomegai=1,nomegai

     matrix(:,:) = MATMUL( selfenergy_omegac(iomegai,:,:,1) , selfenergy_omegac(iomegai,:,:,1) )   &
                + MATMUL( TRANSPOSE(CONJG(selfenergy_omegac(iomegai,:,:,1))) , TRANSPOSE(CONJG(selfenergy_omegac(iomegai,:,:,1))) )

     rdiag=0.d0
     do istate=1,basis%nbf
       rdiag = rdiag - REAL(matrix(istate,istate),dp) * spin_fact / (2.0 * pi)   * 0.5_dp  ! -1/2 comes from the log expansion
     enddo


     tr_gsigma = tr_gsigma + rdiag * weights(iomegai)

   enddo


   write(stdout,*) 'Tr[tGSigma*tGSigma]        ',tr_gsigma 
   deallocate(matrix,eigvec,eigval)

 end select


 !
 ! Output the new HOMO and LUMO energies
 !
 select case(gwmethod)
 case(G0W0,GV,GnW0,GnWn)
   call output_new_homolumo('GW',nstate,occupation,energy_qp_new,nsemin,nsemax)
 case(COHSEX)
   call output_new_homolumo('COHSEX',nstate,occupation,energy_qp_new,nsemin,nsemax)
 end select




 call clean_deallocate('Temporary array',bra)
 if(ALLOCATED(bra_exx)) call clean_deallocate('Temporary array for LW',bra_exx)

 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
   if( calc_type%gwmethod == LW .OR. calc_type%gwmethod == LW2 .OR. calc_type%gwmethod == GSIGMA ) &
       call calculate_eri_3center_eigen_mixed(basis%nbf,nstate,c_matrix)
 endif

 if(ALLOCATED(omegac)) deallocate(omegac)
 if(ALLOCATED(weights)) deallocate(weights)
 if(ALLOCATED(selfenergy_omega)) deallocate(selfenergy_omega)
 if(ALLOCATED(selfenergy_omegac)) deallocate(selfenergy_omegac)


 call stop_clock(timing_self)


end subroutine gw_selfenergy


!=========================================================================
