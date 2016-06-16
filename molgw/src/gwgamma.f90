!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the GW self-energy with vertex function
! within different flavors: G0W0GAMMA0
!
!=========================================================================
subroutine gwgamma_selfenergy(nstate,gwmethod,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy,energy_gw)
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
 real(dp),allocatable     :: selfenergy_omega(:,:,:,:)
 real(dp),allocatable     :: selfenergy_omega_gw(:,:,:,:)
 real(dp),allocatable     :: selfenergy_omega_gamma(:,:,:,:)
 real(dp),allocatable     :: selfenergy_omega_sox(:,:,:,:)
 real(dp),allocatable  :: sigma_xc_m_vxc_diag(:)
 integer               :: ndim2
 integer               :: astate,bstate,cstate
 integer               :: istate,jstate,kstate,ispin,stpole
 integer               :: iastate,jbstate,kcstate
 integer               :: kcmstate
 integer               :: mstate
 real(dp),allocatable  :: bra(:,:)
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: fact_full_a,fact_empty_a
 real(dp)              :: vcoul,vcoul1,vcoul2
 real(dp)              :: zz_a(nspin)
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: zz(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin),energy_qp_z(nstate,nspin)
 real(dp)              :: energy_qp_z_a(nspin),energy_qp_omega(nspin)
 character(len=3)      :: ctmp
 integer               :: reading_status
 integer               :: selfenergyfile
 integer               :: nsemin,nsemax
!=====

 call start_clock(timing_self)

 write(stdout,*)
 select case(gwmethod)
 case(G0W0GAMMA0)
   write(stdout,*) 'Perform a one-shot G0W0GAMMA0 calculation'
 end select

 nprodbasis = index_prodstate(nstate,nstate)

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix)
 else
   stop'NOT implemented'
 endif

 !
 ! Set the range of states on which to evaluate the self-energy
 nsemin = MAX(ncore_G+1   ,selfenergy_state_min,1)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax
 call clean_allocate('Temporary array',bra,nstate,nstate)


 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)


 nomegai = nomega_sigma/2
 allocate(omegai(-nomegai:nomegai))
 do iomegai=-nomegai,nomegai
   omegai(iomegai) = step_sigma * iomegai
 enddo


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
 case default
   energy_qp(:,:) = energy(:,:)

 end select

 !
 !
 allocate(selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,1,nspin))
 allocate(selfenergy_omega_gamma(-nomegai:nomegai,nsemin:nsemax,1,nspin))
 allocate(selfenergy_omega_sox(-nomegai:nomegai,nsemin:nsemax,1,nspin))

 if( ALLOCATED(selfenergy_omega_gamma) ) selfenergy_omega_gamma(:,:,:,:)  = 0.0_dp
 if( ALLOCATED(selfenergy_omega_sox) )   selfenergy_omega_sox(:,:,:,:)  = 0.0_dp

 write(stdout,*) 'Calculate SOX'

 do ispin=1,nspin

   !==========================
   do kstate=ncore_G+1,nvirtual_G-1
     if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle
     do istate=ncore_G+1,nvirtual_G-1
       if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
       do bstate=ncore_G+1,nvirtual_G-1
         if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty) cycle

         do mstate=nsemin,nsemax

           vcoul1 = DOT_PRODUCT( eri_3center_eigen(:,mstate,istate,ispin) , eri_3center_eigen(:,bstate,kstate,ispin) )
           vcoul2 = DOT_PRODUCT( eri_3center_eigen(:,istate,bstate,ispin) , eri_3center_eigen(:,kstate,mstate,ispin) )
           !
           ! calculate only the diagonal !
           do iomegai=-nomegai,nomegai
             selfenergy_omega_sox(iomegai,mstate,1,ispin) = selfenergy_omega_sox(iomegai,mstate,1,ispin) &
                 - vcoul1 * vcoul2            &
                   *  REAL(  1.0_dp / ( energy_qp(mstate,ispin) + omegai(iomegai) - energy_qp(istate,ispin) - energy_qp(kstate,ispin) + energy_qp(bstate,ispin) - ieta )  , dp ) 
           enddo
         enddo

       enddo
     enddo
   enddo


   !==========================
   do cstate=ncore_G+1,nvirtual_G-1
     if( (spin_fact - occupation(cstate,ispin)) / spin_fact < completely_empty) cycle
     do jstate=ncore_G+1,nvirtual_G-1
       if( occupation(jstate,ispin) / spin_fact < completely_empty ) cycle
       do astate=ncore_G+1,nvirtual_G-1
         if( (spin_fact - occupation(astate,ispin)) / spin_fact < completely_empty) cycle

         do mstate=nsemin,nsemax

           vcoul1 = DOT_PRODUCT( eri_3center_eigen(:,mstate,astate,ispin) , eri_3center_eigen(:,jstate,cstate,ispin) )
           vcoul2 = DOT_PRODUCT( eri_3center_eigen(:,astate,jstate,ispin) , eri_3center_eigen(:,cstate,mstate,ispin) )
           !
           ! calculate only the diagonal !
           do iomegai=-nomegai,nomegai
             selfenergy_omega_sox(iomegai,mstate,1,ispin) = selfenergy_omega_sox(iomegai,mstate,1,ispin) &
                 - vcoul1 * vcoul2            &
                   *  REAL(  1.0_dp / ( energy_qp(mstate,ispin) + omegai(iomegai) - energy_qp(astate,ispin) - energy_qp(cstate,ispin) + energy_qp(jstate,ispin) + ieta )  , dp ) 
           enddo
         enddo

       enddo
     enddo
   enddo


 enddo


#if 0
 write(stdout,*) 'Calculate SOSEX'

 do ispin=1,nspin
   do stpole=1,wpol%npole_reso

     write(stdout,*) 'FBFB poles:',stpole,' / ',wpol%npole_reso

     do kcstate=1,nstate
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,kcstate)     = MATMUL( wpol%residu_left(:,stpole) , eri_3center_eigen(:,:,kcstate,ispin) )
     enddo
     call xsum(bra)


     !==========================
     do kstate=ncore_G+1,nvirtual_G-1
       if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle
       do istate=ncore_G+1,nvirtual_G-1
         if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
         do bstate=ncore_G+1,nvirtual_G-1
           if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty) cycle

           !
           ! calculate only the diagonal !
           do mstate=nsemin,nsemax

             vcoul = DOT_PRODUCT( eri_3center_eigen(:,istate,kstate,ispin) , eri_3center_eigen(:,bstate,mstate,ispin) )
             do iomegai=-nomegai,nomegai
               selfenergy_omega_gamma(iomegai,mstate,1,ispin) = selfenergy_omega_gamma(iomegai,mstate,1,ispin) &
                        - bra(mstate,kstate) * bra(istate,bstate) * vcoul                          &  
                          *  REAL(  1.0_dp / ( energy_qp(mstate,ispin) + omegai(iomegai) - energy_qp(kstate,ispin) + wpol%pole(stpole) - ieta )  , dp )  &
                          *  REAL(  1.0_dp / ( -wpol%pole(stpole) + energy_qp(istate,ispin) - energy_qp(bstate,ispin) + ieta )  , dp ) 
             enddo
           enddo

         enddo
       enddo
     enddo

     !==========================
     do cstate=ncore_G+1,nvirtual_G-1
       if( (spin_fact - occupation(cstate,ispin)) / spin_fact < completely_empty ) cycle
       do istate=ncore_G+1,nvirtual_G-1
         if( occupation(istate,ispin) / spin_fact < completely_empty ) cycle
         do bstate=ncore_G+1,nvirtual_G-1
           if( (spin_fact - occupation(bstate,ispin)) / spin_fact < completely_empty ) cycle

           !
           ! calculate only the diagonal !
           do mstate=nsemin,nsemax

             vcoul = DOT_PRODUCT( eri_3center_eigen(:,istate,cstate,ispin) , eri_3center_eigen(:,bstate,mstate,ispin) ) 
             do iomegai=-nomegai,nomegai
               selfenergy_omega_gamma(iomegai,mstate,1,ispin) = selfenergy_omega_gamma(iomegai,mstate,1,ispin) &
                        - bra(mstate,cstate) * bra(istate,bstate) * vcoul                          &  
                          *  REAL(  1.0_dp / ( energy_qp(mstate,ispin) + omegai(iomegai) - energy_qp(cstate,ispin) - wpol%pole(stpole) + ieta )  , dp )  &
                          *  REAL(  1.0_dp / ( energy_qp(mstate,ispin) + omegai(iomegai) - energy_qp(cstate,ispin) + energy_qp(istate,ispin) - energy_qp(bstate,ispin) + ieta )  , dp ) 
             enddo
           enddo

         enddo
       enddo
     enddo

     !==========================
     do kstate=ncore_G+1,nvirtual_G-1
       if( occupation(kstate,ispin) / spin_fact < completely_empty ) cycle
       do astate=ncore_G+1,nvirtual_G-1
         if( (spin_fact - occupation(astate,ispin)) / spin_fact < completely_empty  ) cycle
         do jstate=ncore_G+1,nvirtual_G-1
           if( occupation(jstate,ispin) / spin_fact < completely_empty ) cycle

           !
           ! calculate only the diagonal !
           do mstate=nsemin,nsemax

             vcoul = DOT_PRODUCT( eri_3center_eigen(:,astate,kstate,ispin) , eri_3center_eigen(:,jstate,mstate,ispin) )
             do iomegai=-nomegai,nomegai
               selfenergy_omega_gamma(iomegai,mstate,1,ispin) = selfenergy_omega_gamma(iomegai,mstate,1,ispin) &
                        - bra(mstate,kstate) * bra(astate,jstate) * vcoul                          &  
                          *  REAL(  1.0_dp / ( energy_qp(mstate,ispin) + omegai(iomegai) - energy_qp(kstate,ispin) + energy_qp(astate,ispin) - energy_qp(jstate,ispin)  - ieta )  , dp )  &
                          *  REAL(  1.0_dp / ( energy_qp(jstate,ispin) - energy_qp(astate,ispin) - wpol%pole(stpole) + ieta )  , dp ) 
             enddo
           enddo

         enddo
       enddo
     enddo




   enddo !stpole
 enddo !ispin
#endif

 ! Sum up the contribution from different poles (= different procs)
 if( ALLOCATED(selfenergy_omega_gamma) ) then
   call xsum(selfenergy_omega_gamma)
 endif
 if( ALLOCATED(selfenergy_omega_sox) ) then
   call xsum(selfenergy_omega_sox)
 endif


 write(stdout,'(a)') ' Sigma_c(omega) is calculated'

 allocate(selfenergy_omega_gw(-nomegai:nomegai,nsemin:nsemax,1,nspin))

 open(newunit=selfenergyfile,file='g0w0.dat',status='old',form='unformatted')
 do ispin=1,nspin
   do astate=nsemin,nsemax
     read(selfenergyfile) selfenergy_omega_gw(:,astate,1,ispin)
   enddo
 enddo
 close(selfenergyfile,status='delete')


 forall(astate=nsemin:nsemax)
   selfenergy_omega(:,astate,1,:) = selfenergy_omega_gw(:,astate,1,:) + selfenergy_omega_sox(:,astate,1,:) + selfenergy_omega_gamma(:,astate,1,:)
 end forall


 if(print_sigma_ .AND. is_iomaster ) then

   do astate=nsemin,nsemax
     write(ctmp,'(i3.3)') astate
     write(stdout,'(x,a,a)') 'Printing file: ','selfenergygamma_omega_state'//TRIM(ctmp)
     open(newunit=selfenergyfile,file='selfenergygamma_omega_state'//TRIM(ctmp))
     write(selfenergyfile,'(a)') '# omega (eV)             SigmaC (eV)    omega - e_KS - Vxc + SigmaX (eV)     A (eV^-1) '
     do iomegai=-nomegai,nomegai
       write(selfenergyfile,'(20(f12.6,2x))') ( omegai(iomegai) + energy_qp(astate,:) )*Ha_eV,            &
                                          selfenergy_omega(iomegai,astate,1,:) * Ha_eV,   &
                                          selfenergy_omega_gw(iomegai,astate,1,:) * Ha_eV,   &
                                          selfenergy_omega_sox(iomegai,astate,1,:) * Ha_eV,   &
                                          selfenergy_omega_gamma(0,astate,1,:) * Ha_eV,   &
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

 write(stdout,'(/,a)') ' G0W0Gamma0 Eigenvalues (eV)'
 if(nspin==1) then
   write(stdout,'(a)') '   #          E0        SigX-Vxc    SigC_G0W0    SigC_SOX   SigC_Gamma0   SigC_TOT      Z         G0W0_Z         G0W0_qp'
 else
   write(stdout,'(a)') &
     '   #                E0                      SigX-Vxc                    SigC                       Z                       G0W0_Z                      G0W0_qp'
 endif

 do astate=nsemin,nsemax
   write(stdout,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,          & 
                                      exchange_m_vxc_diag(astate,:)*Ha_eV,       &
                                      selfenergy_omega_gw(0,astate,1,:)*Ha_eV,           &
                                      selfenergy_omega_sox(0,astate,1,:)*Ha_eV,  &
                                      selfenergy_omega_gamma(0,astate,1,:)*Ha_eV, &
                                      selfenergy_omega(0,astate,1,:)*Ha_eV, &
                                      zz(astate,:),                              & 
                                      energy_qp_z(astate,:)*Ha_eV,               &
                                      energy_qp_new(astate,:)*Ha_eV
 enddo


 call write_energy_qp(nstate,energy_qp_new)




 !
 ! Output the new HOMO and LUMO energies
 !
 select case(gwmethod)
 case(G0W0GAMMA0)
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

 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
   if( calc_type%gwmethod == LW .OR. calc_type%gwmethod == LW2 .OR. calc_type%gwmethod == GSIGMA ) &
       call calculate_eri_3center_eigen_mixed(basis%nbf,nstate,c_matrix)
 endif

 if(ALLOCATED(omegai)) deallocate(omegai)
 if(ALLOCATED(selfenergy_omega)) deallocate(selfenergy_omega)

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


end subroutine gwgamma_selfenergy






!=========================================================================
subroutine gwgamma2_selfenergy(nstate,gwmethod,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy,energy_gw)
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
 complex(dp),allocatable  :: omegai(:)
 complex(dp),allocatable     :: selfenergy_omega(:,:,:,:)
 real(dp),allocatable     :: selfenergy_omega_gw(:,:,:,:)
 complex(dp),allocatable     :: selfenergy_omega_sox(:,:,:,:)
 real(dp),allocatable  :: sigma_xc_m_vxc_diag(:)
 integer               :: ndim2
 integer               :: astate,bstate,cstate
 integer               :: istate,jstate,kstate,ispin,stpole
 integer               :: iastate,jbstate,kcstate
 integer               :: kcmstate
 integer               :: mstate
 real(dp),allocatable  :: bra(:,:)
 real(dp)              :: fact_full_i,fact_empty_i
 real(dp)              :: fact_full_a,fact_empty_a
 real(dp)              :: vcoul,vcoul1,vcoul2
 real(dp)              :: zz_a(nspin)
 real(dp)              :: energy_qp(nstate,nspin)
 real(dp)              :: zz(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin),energy_qp_z(nstate,nspin)
 real(dp)              :: energy_qp_z_a(nspin),energy_qp_omega(nspin)
 character(len=3)      :: ctmp
 integer               :: reading_status
 integer               :: selfenergyfile
 integer               :: nsemin,nsemax
 real(dp)              :: mu
 integer               :: qstate,rstate,sstate,pstate
 real(dp)              :: eq,er,es
 real(dp)              :: fq,fr
 real(dp),allocatable  :: x1(:),weights(:)
 real(dp)              :: omega
!=====

 call start_clock(timing_self)

 write(stdout,*)
 select case(gwmethod)
 case(G0W0GAMMA0)
   write(stdout,*) 'Perform a one-shot G0W0GAMMA0 imaginary calculation'
 end select

 nprodbasis = index_prodstate(nstate,nstate)

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix)
 else
   stop'NOT implemented'
 endif

 !
 ! Set the range of states on which to evaluate the self-energy
 nsemin = MAX(ncore_G+1   ,selfenergy_state_min,1)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate)

 write(stdout,'(a,i4,a,i4)') ' Calculate state range from ',nsemin,' to ',nsemax
 call clean_allocate('Temporary array',bra,nstate,nstate)


 energy_gw = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)


 nomegai =  24

 allocate(omegai(1:nomegai))

 allocate(x1(nomegai),weights(nomegai))
 call coeffs_gausslegint(0.0_dp,1.0_dp,x1,weights,nomegai)

 mu =  0.00_dp / Ha_eV
 write(stdout,*) 'mu is set manually here to (Ha)',mu
 write(stdout,*) 'mu is set manually here to (eV)',mu*Ha_eV
 do iomegai=1,nomegai
   omegai(iomegai) = x1(iomegai) / ( 1.0_dp - x1(iomegai) ) * im + mu
   weights(iomegai) = weights(iomegai) / (1.0_dp - x1(iomegai))**2
   write(stdout,'(i4,3(2x,f12.4))') iomegai,REAL(omegai(iomegai)),AIMAG(omegai(iomegai)),weights(iomegai)
 enddo
 deallocate(x1)



 energy_qp(:,:) = energy(:,:)


 !
 !
 allocate(selfenergy_omega_sox(0:0,nsemin:nsemax,1,nspin))

 if( ALLOCATED(selfenergy_omega_sox) )   selfenergy_omega_sox(:,:,:,:)  = 0.0_dp

 write(stdout,*) 'Calculate SOX'

 omega = -4.00 / Ha_eV

 do iomegai=1,nomegai

   do ispin=1,nspin

       !==========================
       do qstate=ncore_G+1,nvirtual_G-1
         fq = occupation(qstate,ispin) / spin_fact
         eq = energy_qp(qstate,ispin)
         do rstate=ncore_G+1,nvirtual_G-1
           er = energy_qp(rstate,ispin)
           do sstate=ncore_G+1,nvirtual_G-1
             es = energy_qp(sstate,ispin)

             do pstate=nsemin,nsemax

               vcoul1 = DOT_PRODUCT( eri_3center_eigen(:,pstate,qstate,ispin) , eri_3center_eigen(:,rstate,sstate,ispin) ) &
                       -DOT_PRODUCT( eri_3center_eigen(:,pstate,sstate,ispin) , eri_3center_eigen(:,rstate,qstate,ispin) )

               vcoul2 = DOT_PRODUCT( eri_3center_eigen(:,pstate,qstate,ispin) , eri_3center_eigen(:,rstate,sstate,ispin) ) 
               vcoul2 = vcoul1

               ! Positive omegai
               selfenergy_omega_sox(0,pstate,1,ispin) = selfenergy_omega_sox(0,pstate,1,ispin) &
                     - vcoul1 * vcoul2 * weights(iomegai) / (2.0_dp * pi ) * (fq-fr)             &
                       / ( omega + omegai(iomegai) - es )                   &
                       / ( omegai(iomegai) - er + eq )

               ! Negative omegai
               selfenergy_omega_sox(0,pstate,1,ispin) = selfenergy_omega_sox(0,pstate,1,ispin) &
                     - vcoul1 * vcoul2 * weights(iomegai) / (2.0_dp * pi ) * (fq-fr)             &
                       / ( omega - omegai(iomegai) - es )                   &
                       / ( -omegai(iomegai) - er + eq )

             enddo

           enddo
         enddo
       enddo

   enddo !ispin
 enddo !iomegai

 write(stdout,'(a)') ' Sigma_c(omega) is calculated'

 do ispin=1,nspin
 do pstate=nsemin,nsemax
   write(stdout,'(x,a,10(2x,f16.6))') 'Result:',omega*Ha_eV,&
                      selfenergy_omega_sox(0,pstate,1,ispin)*Ha_eV
 enddo
 enddo

 stop'ENOUGH'




 do iomegai=1,nomegai

   do ispin=1,nspin

     do stpole=1,wpol%npole_reso

       write(stdout,*) 'FBFB poles:',stpole,' / ',wpol%npole_reso

       do kcstate=1,nstate
         ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
         bra(:,kcstate)     = MATMUL( wpol%residu_left(:,stpole) , eri_3center_eigen(:,:,kcstate,ispin) )
       enddo
       call xsum(bra)

       !==========================
       do qstate=ncore_G+1,nvirtual_G-1
         fq = occupation(qstate,ispin) / spin_fact
         eq = energy_qp(qstate,ispin)
         do rstate=ncore_G+1,nvirtual_G-1
           er = energy_qp(rstate,ispin)
           do sstate=ncore_G+1,nvirtual_G-1
             es = energy_qp(sstate,ispin)

             do pstate=nsemin,nsemax

               vcoul1 = DOT_PRODUCT( eri_3center_eigen(:,pstate,qstate,ispin) , eri_3center_eigen(:,rstate,sstate,ispin) ) &
                       -DOT_PRODUCT( eri_3center_eigen(:,pstate,sstate,ispin) , eri_3center_eigen(:,rstate,qstate,ispin) )

               ! Positive omegai
               vcoul2 = bra(pstate,qstate) * bra(rstate,sstate) &
                       * ( 1.0_dp / ( omegai(iomegai) - wpol%pole(stpole) ) - 1.0_dp / ( omegai(iomegai) + wpol%pole(stpole) ) )   &
                                   + DOT_PRODUCT( eri_3center_eigen(:,pstate,qstate,ispin) , eri_3center_eigen(:,rstate,sstate,ispin) ) 

               selfenergy_omega_sox(0,pstate,1,ispin) = selfenergy_omega_sox(0,pstate,1,ispin) &
                     - vcoul1 * vcoul2 * weights(iomegai) / (2.0_dp * pi ) * (fq-fr)             &
                       / ( energy_qp(pstate,ispin)+0.10_dp + omegai(iomegai) - es )                   &
                       / ( omegai(iomegai) - er + eq )

               ! Negative omegai
               vcoul2 = bra(pstate,qstate) * bra(rstate,sstate) &
                       * ( 1.0_dp / ( -omegai(iomegai) - wpol%pole(stpole) ) - 1.0_dp / ( -omegai(iomegai) + wpol%pole(stpole) ) )   &
                                   + DOT_PRODUCT( eri_3center_eigen(:,pstate,qstate,ispin) , eri_3center_eigen(:,rstate,sstate,ispin) ) 

               selfenergy_omega_sox(0,pstate,1,ispin) = selfenergy_omega_sox(0,pstate,1,ispin) &
                     - vcoul1 * vcoul2 * weights(iomegai) / (2.0_dp * pi ) * (fq-fr)             &
                       / ( energy_qp(pstate,ispin)+0.10_dp - omegai(iomegai) - es )                   &
                       / ( -omegai(iomegai) - er + eq )

             enddo

           enddo
         enddo
       enddo

     enddo !stpole

   enddo !ispin
 enddo !iomegai


 write(stdout,'(a)') ' Sigma_c(omega) is calculated'

 
 do ispin=1,nspin
 do pstate=nsemin,nsemax
   write(stdout,'(x,a,10(2x,f16.6))') 'Result:',(energy_qp(pstate,ispin)+0.10_dp)*Ha_eV,&
                      selfenergy_omega_sox(0,pstate,1,ispin)*Ha_eV
 enddo
 enddo

 stop'ENOUGH'


 call stop_clock(timing_self)




end subroutine gwgamma2_selfenergy



!=========================================================================
