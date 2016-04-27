!=========================================================================
! This file is part of MOLGW.
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
 logical               :: file_exists=.FALSE.
 integer               :: nprodbasis
 integer               :: homo
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 complex(dp),allocatable  :: omegac(:)
 real(dp),allocatable     :: selfenergy_omega(:,:,:,:)
 complex(dp),allocatable  :: selfenergy_omegac(:,:,:,:)
 real(dp),allocatable  :: sigma_xc_m_vxc_diag(:)
 integer               :: ndim2
 integer               :: bbf,ibf,iastate
 integer               :: astate,bstate
 integer               :: istate,ispin,ipole
 real(dp),allocatable  :: bra(:,:)
 real(dp),allocatable  :: bra_exx(:,:)
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: fact_full_a,fact_empty_a
 real(dp)              :: zz_a(nspin)
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: zz(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin),energy_qp_z(nstate,nspin)
 real(dp)              :: energy_qp_z_a(nspin),energy_qp_omega(nspin)
 character(len=3)      :: ctmp
 integer               :: reading_status
 integer               :: selfenergyfile
 integer               :: nsemin,nsemax
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
 case(G0W0)
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
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix)
   if( calc_type%gwmethod == LW .OR. calc_type%gwmethod == LW2 .OR. calc_type%gwmethod == GSIGMA ) &
       call calculate_eri_3center_eigen_mixed(basis%nbf,nstate,c_matrix)
 endif

 !
 ! Set the range of states on which to evaluate the self-energy
 nsemin = MAX(ncore_G+1   ,selfenergy_state_min,1)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax
 call clean_allocate('Temporary array',bra,1,wpol%npole_reso,nsemin,nsemax)

 if(gwmethod==LW .OR. gwmethod==LW2 .OR. gwmethod==GSIGMA) then
   call clean_allocate('Temporary array for LW',bra_exx,1,wpol%npole_reso,nsemin,nsemax)
 endif


 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)


 select case(gwmethod)
 case(GV,QS,COHSEX,GSIGMA3,QSCOHSEX,GnW0,GnWn)
   nomegai = 0
   allocate(omegai(-nomegai:nomegai))
   omegai(0) = 0.0_dp

 case(GSIGMA)
   nomegai = 1
   allocate(omegai(-nomegai:nomegai))
   omegai(:) = 0.0_dp

 case(LW,LW2)
   nomegai =  32

   allocate(omegac(1:nomegai))
   
   allocate(x1(nomegai),weights(nomegai))  
   call coeffs_gausslegint(0.0_dp,1.0_dp,x1,weights,nomegai)

   mu =-0.10_dp
   write(stdout,*) 'mu is set manually here to',mu
   do iomegai=1,nomegai
     omegac(iomegai) = x1(iomegai) / ( 1.0_dp - x1(iomegai) ) * im + mu
     weights(iomegai) = weights(iomegai) / (1.0_dp - x1(iomegai))**2
     write(stdout,*) iomegai,AIMAG(omegac(iomegai)),weights(iomegai)
   enddo
   deallocate(x1)

 case default
   nomegai = nomega_sigma/2
   allocate(omegai(-nomegai:nomegai))
   do iomegai=-nomegai,nomegai
     omegai(iomegai) = step_sigma * iomegai
   enddo

 end select

 !
 ! Which calculation type needs to update energy_qp
 !
 select case(gwmethod)
 case(GnW0,GnWn,GSIGMA3)
   call read_energy_qp(nstate,energy_qp,reading_status)
   if(reading_status/=0) then
     call issue_warning('File energy_qp not found: assuming 1st iteration')
     energy_qp(:,:) = energy(:,:)
   endif

 case(LW,LW2,GSIGMA)
   call issue_warning('reading G\tilde')
   open(1001,form='unformatted')
   read(1001) energy_qp(:,:)
   close(1001,status='delete')

 case default
   energy_qp(:,:) = energy(:,:)

 end select

 !
 ! Which calculation type needs the diagonal only?
 ! Which calculation type needs a complex sigma?
 !
 select case(gwmethod)
 case(COHSEX,GV,GSIGMA,G0W0,GnW0,GnWn)   ! diagonal real
   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,1,nspin))

 case(QS,QSCOHSEX,GSIGMA3) ! matrix real
   allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nsemin:nsemax,nspin))

 case(LW,LW2)                     ! matrix complex
   allocate(selfenergy_omegac(1:nomegai,nsemin:nsemax,nsemin:nsemax,nspin))

 case default
   call die('GW case does not exist. Should not happen')
 end select
 if( ALLOCATED(selfenergy_omega) )  selfenergy_omega(:,:,:,:)  = 0.0_dp
 if( ALLOCATED(selfenergy_omegac) ) selfenergy_omegac(:,:,:,:) = 0.0_dp



 do ispin=1,nspin
   do istate=1,nstate !INNER LOOP of G
     if(gwmethod==LW .or. gwmethod==LW2) write(stdout,*) 'FBFB LW',istate
     !
     ! Apply the frozen core and frozen virtual approximation to G
     if(istate <= ncore_G)    cycle
     if(istate >= nvirtual_G) cycle

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
       call xsum(bra)
       if( gwmethod==LW .OR. gwmethod==LW2 .OR. gwmethod==GSIGMA) then
         bra_exx(:,nsemin:nsemax) = MATMUL( TRANSPOSE(wpol%residu_left(:,:)) , eri_3center_eigen_mixed(:,istate,nsemin:nsemax,ispin) )
         call xsum(bra_exx)
       endif
     endif


     ! The application of residu theorem only retains the pole in certain
     ! quadrants.
     ! The positive poles of W go with the poles of occupied states in G
     ! The negative poles of W go with the poles of empty states in G
     fact_full_i   = occupation(istate,ispin) / spin_fact
     fact_empty_i = (spin_fact - occupation(istate,ispin)) / spin_fact


     do ipole=1,wpol%npole_reso

       if( rank /= MODULO(ipole-1,nproc) ) cycle

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
   call xsum(selfenergy_omega)
 endif
 if( ALLOCATED(selfenergy_omegac) ) then
   call xsum(selfenergy_omegac)
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
 case(QS)
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,nsemin:nsemax,ispin) ) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,nsemin:nsemax,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(QSCOHSEX)
   !FBFB
   write(stdout,*)  'FBFB    #    sigc  '
   do ispin=1,nspin
     do astate=nsemin,nsemax
       write(stdout,*) astate,selfenergy_omega(0,astate,astate,ispin)
     enddo
   enddo
   write(stdout,*) 
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
   
   write(stdout,'(/,a)') ' Gv     Eigenvalues (eV)'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do astate=nsemin,nsemax
     zz_a(:) = 1.0_dp 
     energy_qp_new(astate,:) = energy_qp(astate,:) + selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,     &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV, &
                                           zz_a(:),energy_qp_new(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nstate,energy_qp_new)


 case(COHSEX) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall
   
   write(stdout,'(/,a)') ' COHSEX Eigenvalues (eV)'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do astate=nsemin,nsemax
     zz_a(:) = 1.0_dp 
     energy_qp_new(astate,:) = energy_qp(astate,:) + selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,     &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV, &
                                           zz_a(:),energy_qp_new(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nstate,energy_qp_new)


 case(G0W0) !==========================================================

   if(print_sigma_ .AND. is_iomaster ) then

     do astate=nsemin,nsemax
       write(ctmp,'(i3.3)') astate
       open(newunit=selfenergyfile,file='selfenergy_omega_state'//TRIM(ctmp))
       write(selfenergyfile,'(a)') '# omega (eV)             SigmaC (eV)    omega - e_KS - Vxc + SigmaX (eV)     A (eV^-1) '
       do iomegai=-nomegai,nomegai
         write(selfenergyfile,'(20(f12.6,2x))') ( omegai(iomegai) + energy_qp(astate,:) )*Ha_eV,               &
                                            selfenergy_omega(iomegai,astate,1,:)*Ha_eV,                    &
                                            ( omegai(iomegai) - exchange_m_vxc_diag(astate,:) )*Ha_eV,     &
                                            1.0_dp/pi/ABS( omegai(iomegai) - exchange_m_vxc_diag(astate,:) &
                                                    - selfenergy_omega(iomegai,astate,1,:) ) / Ha_eV
       enddo
       write(selfenergyfile,*)
       close(selfenergyfile)
     enddo

   endif


   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall

   zz(:,:) = 0.0_dp
   energy_qp_z(:,:) = 0.0_dp
   energy_qp_new(:,:) = 0.0_dp

   ! Then overwrite the interesting energy with the calculated GW one
   do astate=nsemin,nsemax

     if( MODULO(astate-nsemin,nproc) /= rank ) cycle

     zz_a(:) = ( selfenergy_omega(1,astate,1,:) - selfenergy_omega(-1,astate,1,:) ) / ( omegai(1) - omegai(-1) )
     zz_a(:) = 1.0_dp / ( 1.0_dp - zz_a(:) )
     ! Contrain Z to be in [0:1] to avoid crazy values
     do ispin=1,nspin
       zz_a(ispin) = MIN( MAX(zz_a(ispin),0.0_dp) , 1.0_dp )
     enddo

     energy_qp_z_a(:) = energy_qp(astate,:) + zz_a(:) * ( selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:) )

     allocate(sigma_xc_m_vxc_diag(-nomegai:nomegai))
     do ispin=1,nspin
       sigma_xc_m_vxc_diag(:) = selfenergy_omega(:,astate,1,ispin) + exchange_m_vxc_diag(astate,ispin)
       energy_qp_omega(ispin) = find_fixed_point(nomegai,omegai,sigma_xc_m_vxc_diag) + energy_qp(astate,ispin) 
     enddo
     deallocate(sigma_xc_m_vxc_diag)

     zz(astate,:)            = zz_a(:)
     energy_qp_z(astate,:)   = energy_qp_z_a(:)
     energy_qp_new(astate,:) = energy_qp_omega(:) 
   enddo

   call xsum(zz)
   call xsum(energy_qp_z)
   call xsum(energy_qp_new)

   energy_qp_new(:nsemin-1,:) = energy(:nsemin-1,:)
   energy_qp_new(nsemax+1:,:) = energy(nsemax+1:,:)

   write(stdout,'(/,a)') ' G0W0 Eigenvalues (eV)'
   if(nspin==1) then
     write(stdout,'(a)') '   #          E0        SigX-Vxc      SigC          Z         G0W0_Z         G0W0_qp'
   else
     write(stdout,'(a)') &
       '   #                E0                      SigX-Vxc                    SigC                       Z                       G0W0_Z                      G0W0_qp'
   endif

   do astate=nsemin,nsemax
     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,          & 
                                        exchange_m_vxc_diag(astate,:)*Ha_eV,       &
                                        selfenergy_omega(0,astate,1,:)*Ha_eV,      &
                                        zz(astate,:),                              & 
                                        energy_qp_z(astate,:)*Ha_eV,               &
                                        energy_qp_new(astate,:)*Ha_eV
   enddo


   call write_energy_qp(nstate,energy_qp_new)



 case(GnW0,GnWn) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=nsemin:nsemax)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall


   if( gwmethod==GnW0) then
     write(stdout,'(/,a)') ' GnW0 Eigenvalues (eV)'
   else
     write(stdout,'(/,a)') ' GnWn Eigenvalues (eV)'
   endif
   if(nspin==1) then
     write(stdout,'(a)') '  #          E0        SigX-Vxc      SigC          Z          GW(n-1)       GW(n)'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       GW'
   endif

   ! First give energy_qp_new a meaningful default value
   energy_qp_new(:,:) = energy(:,:)
   ! Then overwrite the interesting energy with the calculated GW one
   do astate=nsemin,nsemax
     zz_a(:) = 1.0_dp 

     energy_qp_new(astate,:) = energy(astate,:) + selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy(astate,:)*Ha_eV,                    &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,       &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV,      &
                                                 zz_a(:),                                     &
                                                 energy_qp(astate,:)*Ha_eV,                 &
                                                 energy_qp_new(astate,:)*Ha_eV
   enddo

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
 case(G0W0,GV,COHSEX,GnW0,GnWn)
   do istate=1,nstate
     if( ANY(occupation(istate,:) > completely_empty) ) homo = istate
   enddo
   write(stdout,*)
   if( homo >= nsemin .AND. homo <= nsemax ) then
     write(stdout,'(a,2(2x,f12.6))') ' GW HOMO (eV):',energy_qp_new(homo,:)*Ha_eV
   endif
   if( homo+1 >= nsemin .AND. homo+1 <= nsemax ) then
     write(stdout,'(a,2(2x,f12.6))') ' GW LUMO (eV):',energy_qp_new(homo+1,:)*Ha_eV
   endif
 end select

 call clean_deallocate('Temporary array',bra)
 if(ALLOCATED(bra_exx)) call clean_deallocate('Temporary array for LW',bra_exx)

 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
   if( calc_type%gwmethod == LW .OR. calc_type%gwmethod == LW2 .OR. calc_type%gwmethod == GSIGMA ) &
       call calculate_eri_3center_eigen_mixed(basis%nbf,nstate,c_matrix)
 endif

 if(ALLOCATED(omegai)) deallocate(omegai)
 if(ALLOCATED(omegac)) deallocate(omegac)
 if(ALLOCATED(weights)) deallocate(weights)
 if(ALLOCATED(selfenergy_omega)) deallocate(selfenergy_omega)
 if(ALLOCATED(selfenergy_omegac)) deallocate(selfenergy_omegac)

 call stop_clock(timing_self)


contains


function find_fixed_point(nx,xx,fx) result(fixed_point)
 implicit none
 integer,intent(in)  :: nx
 real(dp),intent(in) :: xx(-nx:nx)
 real(dp),intent(in) :: fx(-nx:nx)
 real(dp)            :: fixed_point
!=====
 integer             :: ix,imin
 real(dp)            :: rmin
 real(dp)            :: gx(-nx:nx)
 real(dp)            :: gpx
!=====


 gx(:) = fx(:) - xx(:)

 rmin = HUGE(1.0_dp)
 do ix=-nx,nx
   if( ABS(gx(ix)) < rmin ) then
     rmin = ABS(gx(ix))
     imin = ix
   endif
 enddo 


 if( imin == -nx .OR. imin == nx) then
   fixed_point = xx(imin)
 else 
   if( gx(imin)*gx(imin+1) < 0.0_dp )  then
     gpx = ( gx(imin+1) - gx(imin) ) / ( xx(imin+1) - xx(imin) )
     fixed_point = xx(imin) - gx(imin) / gpx 
   else if( gx(imin)*gx(imin-1) < 0.0_dp )  then
     gpx = ( gx(imin) - gx(imin-1) ) / ( xx(imin) - xx(imin-1) )
     fixed_point = xx(imin-1) - gx(imin-1) / gpx 
   else
     fixed_point = xx(imin)
   endif
 endif



end function find_fixed_point


end subroutine gw_selfenergy


!=========================================================================
