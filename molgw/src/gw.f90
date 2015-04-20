!=========================================================================
subroutine gw_selfenergy(gwmethod,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam
 use m_warning,only: issue_warning,msg
 use m_basis_set
 use m_spectral_function
 use m_eri,only: eri_3center_eigen
 use m_tools,only: coeffs_gausslegint
 implicit none

 integer,intent(in)                 :: gwmethod
 type(basis_set)                    :: basis,prod_basis
 real(dp),intent(in)                :: occupation(basis%nbf,nspin),energy(basis%nbf,nspin),exchange_m_vxc_diag(basis%nbf,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
!=====
 logical               :: file_exists=.FALSE.
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 complex(dp),allocatable  :: omegac(:)
 real(dp),allocatable  :: selfenergy_omega(:,:,:,:)
 complex(dp),allocatable  :: selfenergy_omegac(:,:,:,:)
 real(dp),allocatable  :: sigma_xc_m_vxc_diag(:)
 integer               :: ndim2
 integer               :: bbf,ibf,kbf
 integer               :: astate,bstate
 integer               :: istate,ispin,ipole
 real(dp)              :: bra(wpol%npole_reso,basis%nbf)
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: fact_full_a,fact_empty_a
 real(dp)              :: zz(nspin)
 real(dp)              :: energy_qp(basis%nbf,nspin)
 real(dp)              :: energy_qp_new(basis%nbf,nspin)
 real(dp)              :: energy_qp_z(nspin),energy_qp_omega(nspin)
 character(len=3)      :: ctmp
 integer               :: reading_status
 integer               :: selfenergyfile
! LW devel
 complex(dp),allocatable :: matrix(:,:),eigvec(:,:)
 real(dp),allocatable :: eigval(:),x1(:),weights(:)
 real(dp) :: tr_log_gsigma,tr_gsigma,rdiag,mu
!=====

 call start_clock(timing_self)


 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)

 write(stdout,*)
 select case(gwmethod)
 case(LW)
   write(stdout,*) 'perform the calculation of Tr[ ln ( 1 - GSigma ) ]'
 case(GSIGMA)
   write(stdout,*) 'perform the calculation of Tr[ SigmaG ]'
 case(GSIGMA2)
   write(stdout,*) 'perform the calculation of Tr[ SigmaG ] numerical for test only'
 case(GV)
   write(stdout,*) 'perform a perturbative HF calculation'
 case(QS)
   write(stdout,*) 'perform a QP self-consistent GW calculation'
 case(G0W0)
   write(stdout,*) 'perform a one-shot G0W0 calculation'
 case(COHSEX)
   write(stdout,*) 'perform a COHSEX calculation'
 case(QSCOHSEX)
   write(stdout,*) 'perform a self-consistent COHSEX calculation'
 case(GnW0)
   write(stdout,*) 'perform an eigenvalue self-consistent GnW0 calculation'
 case(GnWn)
   write(stdout,*) 'perform an eigenvalue self-consistent GnWn calculation'
 end select

 if(gwmethod==GV .OR. gwmethod==GSIGMA .OR. gwmethod==QS .OR. gwmethod==COHSEX  &
   .OR. gwmethod==QSCOHSEX .OR.  gwmethod==GnW0 .OR. gwmethod==GnWn) then
   nomegai = 0
   allocate(omegai(-nomegai:nomegai))
   omegai(0) = 0.0_dp
 else if(gwmethod==LW .OR. gwmethod==GSIGMA2) then
   nomegai =  16  !16000

   allocate(omegac(1:nomegai))
   
   allocate(x1(nomegai),weights(nomegai))  !FBFB deallocate weights
   call coeffs_gausslegint(0.0_dp,1.0_dp,x1,weights,nomegai)

   mu = 0.00_dp
   write(stdout,*) 'mu is set manually here to',mu
   do iomegai=1,nomegai
     omegac(iomegai) = x1(iomegai) / ( 1.0_dp - x1(iomegai) ) * im + mu
     weights(iomegai) = weights(iomegai) / (1.0_dp - x1(iomegai))**2
     write(stdout,*) iomegai,AIMAG(omegac(iomegai)),weights(iomegai)
   enddo
   deallocate(x1)


 else
   nomegai = nomega_sigma/2
   allocate(omegai(-nomegai:nomegai))
   do iomegai=-nomegai,nomegai
     omegai(iomegai) = step_sigma * iomegai
   enddo
 endif

 if( gwmethod==GnW0 .OR. gwmethod==GnWn .OR. gwmethod==GSIGMA .OR. gwmethod==GSIGMA2 .OR. gwmethod==LW ) then
   call read_energy_qp(nspin,basis%nbf,energy_qp,reading_status)
   if(reading_status/=0) then
     call issue_warning('File energy_qp not found: assuming 1st iteration')
     energy_qp(:,:) = energy(:,:)
   endif
 else
   energy_qp(:,:) = energy(:,:)
 endif

 if( gwmethod==GV .OR. gwmethod==GSIGMA .OR. gwmethod==G0W0 .OR.  gwmethod==GSIGMA2 &
    .OR. gwmethod==GnW0 .OR. gwmethod==GnWn ) then
   ndim2=1
 else
   ndim2=basis%nbf
 endif

 if( gwmethod /= LW .AND. gwmethod /=GSIGMA2 ) then
   allocate(selfenergy_omega(-nomegai:nomegai,basis%nbf,ndim2,nspin))
   selfenergy_omega(:,:,:,:) = 0.0_dp
 else
!   allocate(selfenergy_omegac(-nomegai:nomegai,basis%nbf,ndim2,nspin))
   allocate(selfenergy_omegac(1:nomegai,basis%nbf,ndim2,nspin))
   selfenergy_omegac(:,:,:,:) = 0.0_dp
 endif

 do ispin=1,nspin
   do istate=1,basis%nbf !INNER LOOP of G
     if(gwmethod==LW) write(stdout,*) 'FBFB',istate
     !
     ! Apply the frozen core and frozen virtual approximation to G
     if(istate <= ncore_G)    cycle
     if(istate >= nvirtual_G) cycle

     !
     ! Prepare the bra and ket with the knowledge of index istate and astate
     if( .NOT. has_auxil_basis) then
       ! Here just grab the precalculated value
       do astate=1,basis%nbf
         kbf = prod_basis%index_prodbasis(istate,astate) + prod_basis%nbf*(ispin-1)
         bra(:,astate) = wpol%residu_left(kbf,:)
       enddo
     else
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,:) = MATMUL( TRANSPOSE(wpol%residu_left(:,:)) , eri_3center_eigen(:,:,istate,ispin) )
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

         do bstate=1,basis%nbf
           do astate=1,basis%nbf

             selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
                        + bra(ipole,astate) * bra(ipole,bstate)                            &  
                          * ( REAL(  fact_full_i   / ( energy_qp(bstate,ispin) - ieta  - energy_qp(istate,ispin) + wpol%pole(ipole) ) , dp ) &
                            + REAL(  fact_empty_i / ( energy_qp(bstate,ispin) + ieta  - energy_qp(istate,ispin) - wpol%pole(ipole) ) , dp ) )

           enddo
         enddo

       case(LW)

         do bstate=1,basis%nbf
           do astate=1,basis%nbf

!             do iomegai=-nomegai,nomegai
             do iomegai=       1,nomegai
               selfenergy_omegac(iomegai,astate,bstate,ispin) = selfenergy_omegac(iomegai,astate,bstate,ispin) &
                        + bra(ipole,astate) * bra(ipole,bstate) &
                          * (  fact_full_i  / ( omegac(iomegai) - energy(istate,ispin) + wpol%pole(ipole)  )    &
                             + fact_empty_i / ( omegac(iomegai) - energy(istate,ispin) - wpol%pole(ipole) ) )   &
                          / ( omegac(iomegai) - energy_qp(astate,ispin) )
             enddo

           enddo
         enddo

       case(GSIGMA2)

         do astate=1,basis%nbf

!           do iomegai=-nomegai,nomegai
           do iomegai=       1,nomegai
             selfenergy_omegac(iomegai,astate,1,ispin) = selfenergy_omegac(iomegai,astate,1,ispin) &
                      + bra(ipole,astate) * bra(ipole,astate) &
                        * (  fact_full_i  / ( omegac(iomegai) - energy(istate,ispin) + wpol%pole(ipole)  )    &
                           + fact_empty_i / ( omegac(iomegai) - energy(istate,ispin) - wpol%pole(ipole) ) )   &
                        / ( omegac(iomegai) - energy_qp(astate,ispin) )
           enddo

         enddo



       case(GSIGMA)

         do astate=1,basis%nbf
           fact_full_a   = occupation(astate,ispin) / spin_fact
           fact_empty_a  = (spin_fact - occupation(astate,ispin)) / spin_fact
           !
           ! calculate only the diagonal !
           selfenergy_omega(0,astate,1,ispin) = selfenergy_omega(0,astate,1,ispin) &
                    + bra(ipole,astate) * bra(ipole,astate) &
                      * ( REAL(  fact_full_i * fact_empty_a   / ( energy_qp(astate,ispin)  - energy(istate,ispin) + wpol%pole(ipole) - ieta )  , dp )  &
                        - REAL(  fact_empty_i * fact_full_a  / ( energy_qp(astate,ispin)  - energy(istate,ispin) - wpol%pole(ipole) + ieta )  , dp ) )
         enddo


       case(G0W0,GnW0,GnWn)

         do astate=1,basis%nbf
           !
           ! calculate only the diagonal !
           do iomegai=-nomegai,nomegai
             selfenergy_omega(iomegai,astate,1,ispin) = selfenergy_omega(iomegai,astate,1,ispin) &
                      + bra(ipole,astate) * bra(ipole,astate)                                          & 
                        * ( REAL(  fact_full_i   / ( energy_qp(astate,ispin) + omegai(iomegai) - energy_qp(istate,ispin) + wpol%pole(ipole) - ieta )  , dp )  &
                          + REAL(  fact_empty_i / ( energy_qp(astate,ispin) + omegai(iomegai) - energy_qp(istate,ispin) - wpol%pole(ipole) + ieta )  , dp ) )
           enddo
         enddo

       case(COHSEX,QSCOHSEX) 

         !
         ! SEX
         !
         do bstate=1,basis%nbf
           do astate=1,basis%nbf
             selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
                        - bra(ipole,astate) * bra(ipole,bstate)                            & 
                              * fact_full_i / (-wpol%pole(ipole))
           enddo
         enddo
         !
         ! COH
         !
         do bstate=1,basis%nbf
           do astate=1,basis%nbf
             selfenergy_omega(0,astate,bstate,ispin) = selfenergy_omega(0,astate,bstate,ispin) &
                        - bra(ipole,astate) * bra(ipole,bstate) & 
                              * fact_empty_i / wpol%pole(ipole)
           enddo
         enddo

       case(GV)
         !
         ! Do nothing: no correlation in this case
         !
       case default 
         stop'BUG'
       end select

     enddo !ipole

   enddo !istate
 enddo !ispin

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
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,ispin) ) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(QSCOHSEX)
   ! Transform the matrix elements back to the AO basis
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,ispin)) , MATMUL( selfenergy_omega(0,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(GV) !==========================================================

   selfenergy(:,:,:) = selfenergy_omega(0,:,:,:)
   
   write(stdout,'(/,a)') ' COHSEX Eigenvalues [eV]'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do astate=1,basis%nbf
     zz(:) = 1.0_dp 
     energy_qp_new(astate,:) = energy_qp(astate,:) + selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,     &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV, &
                                           zz(:),energy_qp_new(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp_new)


 case(COHSEX) !==========================================================

   selfenergy(:,:,:) = selfenergy_omega(0,:,:,:)
   
   write(stdout,'(/,a)') ' COHSEX Eigenvalues [eV]'
   if(nspin==1) then
     write(stdout,*) '  #          E0        SigX-Vxc      SigC          Z         COHSEX'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       COHSEX'
   endif
   do astate=1,basis%nbf
     zz(:) = 1.0_dp 
     energy_qp_new(astate,:) = energy_qp(astate,:) + selfenergy_omega(0,astate,astate,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,               &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,     &
                                                 selfenergy_omega(0,astate,astate,:)*Ha_eV, &
                                           zz(:),energy_qp_new(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp_new)


 case(G0W0) !==========================================================

   if(print_sigma_ .AND. is_iomaster ) then

     do astate=1,basis%nbf
       write(ctmp,'(i3.3)') astate
       open(newunit=selfenergyfile,file='selfenergy_omega_state'//TRIM(ctmp))
       write(selfenergyfile,'(a)') '# omega + e_KS (eV)     SigmaC (eV)    omega + e_KS - Vxc + SigmaX (eV)     A (eV^-1) '
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
   forall(astate=1:basis%nbf)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall

   write(stdout,'(/,a)') ' G0W0 Eigenvalues [eV]'
   if(nspin==1) then
     write(stdout,'(a)') '   #          E0        SigX-Vxc      SigC          Z         G0W0_Z         G0W0_qp'
   else
     write(stdout,'(a)') &
       '   #                E0                      SigX-Vxc                    SigC                       Z                       G0W0_Z                      G0W0_qp'
   endif
   do astate=1,basis%nbf
     zz(:) = ( selfenergy_omega(1,astate,1,:) - selfenergy_omega(-1,astate,1,:) ) / ( omegai(1) - omegai(-1) )
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )
     ! Contrain Z to be in [0:1] to avoid crazy values
     do ispin=1,nspin
       zz(ispin) = MIN( MAX(zz(ispin),0.0_dp) , 1.0_dp )
     enddo

     energy_qp_z(:) = energy_qp(astate,:) + zz(:) * ( selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:) )

     allocate(sigma_xc_m_vxc_diag(-nomegai:nomegai))
     do ispin=1,nspin
       sigma_xc_m_vxc_diag(:) = selfenergy_omega(:,astate,1,ispin) + exchange_m_vxc_diag(astate,ispin)
       energy_qp_omega(ispin) = find_fixed_point(nomegai,omegai,sigma_xc_m_vxc_diag) + energy_qp(astate,ispin) 
     enddo
     deallocate(sigma_xc_m_vxc_diag)



     write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,                 & 
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,       &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV,      &
                                                 zz(:),                                     & 
                                                 energy_qp_z(:)*Ha_eV,                      &
                                                 energy_qp_omega(:)*Ha_eV

     energy_qp_new(astate,:) = energy_qp_omega(:) 
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp_new)



 case(GnW0,GnWn) !==========================================================

   ! Only had the diagonal calculated...
   selfenergy(:,:,:) = 0.0_dp
   forall(astate=1:basis%nbf)
     selfenergy(astate,astate,:) = selfenergy_omega(0,astate,1,:)
   end forall


   if( gwmethod==GnW0) then
     write(stdout,'(/,a)') ' GnW0 Eigenvalues [eV]'
   else
     write(stdout,'(/,a)') ' GnWn Eigenvalues [eV]'
   endif
   if(nspin==1) then
     write(stdout,'(a)') '  #          E0        SigX-Vxc      SigC          Z          GW(n-1)       GW(n)'
   else
     write(stdout,'(a)') '  #                E0                      SigX-Vxc                    SigC                       Z                       GW'
   endif
   do astate=1,basis%nbf
     zz(:) = 1.0_dp 

     energy_qp_new(astate,:) = energy(astate,:) + selfenergy_omega(0,astate,1,:) + exchange_m_vxc_diag(astate,:)

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy(astate,:)*Ha_eV,                    &
                                                 exchange_m_vxc_diag(astate,:)*Ha_eV,       &
                                                 selfenergy_omega(0,astate,1,:)*Ha_eV,      &
                                                 zz(:),                                     &
                                                 energy_qp(astate,:)*Ha_eV,                 &
                                                 energy_qp_new(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp_new)

 case(GSIGMA) !==========================================================

   write(stdout,*) 'Tr[SigmaG]:',SUM(selfenergy_omega(0,:,1,:)) * spin_fact

 case(GSIGMA2) !==========================================================

   write(stdout,*) 'Tr[SigmaG]: numerical eval'

 rdiag = 0.0_dp
! do iomegai=-nomegai,nomegai
 do iomegai=       1,nomegai
   rdiag = rdiag + SUM(selfenergy_omegac(iomegai,:,1,:)) * spin_fact / (2.0 * pi) * weights(iomegai)
   write(124,'(3(f18.8,2x))') AIMAG(omegac(iomegai)),SUM(selfenergy_omegac(iomegai,:,1,:)) * spin_fact / (2.0 * pi)  * 2.0_dp
 enddo
 write(stdout,*) 'Result:', rdiag   * 2.0_dp

 case(LW)

   write(stdout,*) 'FBFB LW'
   allocate(matrix(basis%nbf,basis%nbf))
   allocate(eigvec(basis%nbf,basis%nbf))
   allocate(eigval(basis%nbf))

   tr_log_gsigma = 0.0_dp
   tr_gsigma     = 0.0_dp

!   do iomegai=-nomegai,nomegai
!     matrix(:,:) = 0.5_dp *(  selfenergy_omegac(iomegai,:,:,1) + CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,1) ))  )
   do iomegai=       1,nomegai
     matrix(:,:) =  selfenergy_omegac(iomegai,:,:,1) + CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,1) ))

     call diagonalize(basis%nbf,matrix,eigval,eigvec)
     write(125,'(3(f18.8,2x))') AIMAG(omegac(iomegai)),SUM(eigval(:)) * spin_fact / (2.0 * pi)
   enddo

!   do iomegai=-nomegai,nomegai
   do iomegai=       1,nomegai

     rdiag=0.d0
     do istate=1,basis%nbf
       rdiag = rdiag + REAL(selfenergy_omegac(iomegai,istate,istate,1),dp) * spin_fact / (2.0 * pi)   * 2.0_dp
     enddo

!     matrix(:,:) = 0.5_dp *(  selfenergy_omegac(iomegai,:,:,1) + CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,1) )) &
!                  - MATMUL( selfenergy_omegac(iomegai,:,:,1) , CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,1) )) )  )
  
     matrix(:,:) = selfenergy_omegac(iomegai,:,:,1) + CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,1) )) &
                  - MATMUL( selfenergy_omegac(iomegai,:,:,1) , CONJG(TRANSPOSE( selfenergy_omegac(iomegai,:,:,1) )) )

     call diagonalize(basis%nbf,matrix,eigval,eigvec)


     write(126,'(3(f18.8,2x))') AIMAG(omegac(iomegai)),-SUM(LOG( 1.0_dp - eigval(:))) * spin_fact / (2.0 * pi)
     write(127,'(3(f18.8,2x))') AIMAG(omegac(iomegai)),rdiag+SUM(LOG( 1.0_dp - eigval(:))) * spin_fact / (2.0 * pi) 
     tr_gsigma = tr_gsigma + rdiag * weights(iomegai)
     tr_log_gsigma = tr_log_gsigma - SUM(LOG( 1.0_dp - eigval(:))) * spin_fact / (2.0 * pi) * weights(iomegai)

   enddo


   write(stdout,*) 'Log     ',tr_log_gsigma 
   write(stdout,*) 'GSigma  ',tr_gsigma 
   write(stdout,*) 'Diff    ',(tr_log_gsigma-tr_gsigma) 
   deallocate(matrix,eigvec,eigval)

 end select


 if(ALLOCATED(omegai)) deallocate(omegai)
 if(ALLOCATED(omegac)) deallocate(omegac)
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
