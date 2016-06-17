!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the MP2 self-energy evaluation
!=========================================================================
subroutine mp2_selfenergy(method,nstate,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,selfenergy,emp2)
 use m_definitions
 use m_mpi
 use m_warning
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_spectral_function
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: method,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(out)       :: selfenergy(basis%nbf,basis%nbf,nspin),emp2
!=====
#ifdef CHI0
 logical,parameter     :: ring_only=.true.
#else
 logical,parameter     :: ring_only=.false.
#endif
 integer               :: astate,bstate,bstate2
 integer               :: homo
 real(dp),allocatable  :: selfenergy_ring(:,:,:,:)
 real(dp),allocatable  :: selfenergy_sox(:,:,:,:)
 real(dp)              :: selfenergy_final(nstate,nstate,nspin)
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 real(dp),allocatable  :: sigma_xc_m_vxc_diag(:)
 integer               :: istate,jstate,kstate
 integer               :: abispin,jkspin,ispin
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek
 real(dp)              :: omega
 real(dp)              :: zz(nspin)
 real(dp)              :: fact_real,fact_energy
 real(dp)              :: emp2_sox,emp2_ring
 logical               :: file_exists
 character(len=3)      :: ctmp
 integer               :: iomegafile
 real(dp),allocatable  :: eri_eigenstate_i(:,:,:,:)
 integer               :: nket1,nket2
 integer               :: nsemin,nsemax
 real(dp)              :: coul_ibjk,coul_ijkb,coul_iakj
 real(dp)              :: energy_qp_z(nstate,nspin)
 real(dp)              :: energy_qp_new(nstate,nspin)
 real(dp)              :: energy_qp_omega(nspin)
 real(dp)              :: ehomo,elumo
 type(spectral_function) :: sf_dummy
!=====

 call start_clock(timing_mp2_self)

 ! This call does not allocate anything fortunately
 call init_spectral_function(nstate,occupation,sf_dummy)

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)

 !
 ! Set the range of states on which to evaluate the self-energy
 nsemin = MAX(ncore_G+1   ,selfenergy_state_min,1)
 nsemax = MIN(nvirtual_G-1,selfenergy_state_max,nstate)

 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 select case(method)
 case(QS)
   write(stdout,*) 'with the QP self-consistent approach'
   nket1 = nsemin
   nket2 = nsemax
 case(perturbative)
   write(stdout,*) 'with the perturbative approach'
   nket1 = 1
   nket2 = 1
 end select


 if(method==QS) then
   nomegai=0
   allocate(omegai(-nomegai:nomegai))
   omegai(0) =  0.00_dp
 else
   nomegai = nomega_sigma/2
   allocate(omegai(-nomegai:nomegai))
   do iomegai=-nomegai,nomegai
     omegai(iomegai) = step_sigma * iomegai
   enddo
 endif


 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix)
 else
   allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
 endif



 allocate(selfenergy_ring(-nomegai:nomegai,nsemin:nsemax,nket1:nket2,nspin))
 allocate(selfenergy_sox (-nomegai:nomegai,nsemin:nsemax,nket1:nket2,nspin))


 selfenergy_ring(:,:,:,:) = 0.0_dp
 selfenergy_sox(:,:,:,:)  = 0.0_dp

 do abispin=1,nspin
   do istate=ncore_G+1,nvirtual_G-1 !LOOP of the first Green's function

     if( .NOT. has_auxil_basis ) then
       call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,istate,abispin,eri_eigenstate_i)
     endif

     fi = occupation(istate,abispin)
     ei = energy(istate,abispin)

     do astate=nsemin,nsemax ! external loop ( bra )
       do bstate=nsemin,nsemax   ! external loop ( ket )
         if( astate /= bstate .AND. method == perturbative ) cycle
         if( method == perturbative ) then
           bstate2 = 1
         else
           bstate2 = bstate
         endif

         do jkspin=1,nspin
           do jstate=ncore_G+1,nvirtual_G-1  !LOOP of the second Green's function
             fj=occupation(jstate,jkspin)
             ej=energy(jstate,jkspin)
  
             do kstate=ncore_G+1,nvirtual_G-1 !LOOP of the third Green's function
               fk=occupation(kstate,jkspin)
               ek=energy(kstate,jkspin)
  
               fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**3
               fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**3
  
               if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle

               if( has_auxil_basis ) then
                 coul_iakj = eri_eigen_ri(istate,astate,abispin,kstate,jstate,jkspin)
                 coul_ibjk = eri_eigen_ri(istate,bstate,abispin,jstate,kstate,jkspin)
                 if( abispin == jkspin ) then
                   coul_ijkb = eri_eigen_ri(istate,jstate,abispin,kstate,bstate,abispin) 
                 endif
               else
                 coul_iakj = eri_eigenstate_i(astate,kstate,jstate,jkspin)
                 coul_ibjk = eri_eigenstate_i(bstate,jstate,kstate,jkspin)
                 if( abispin == jkspin ) then
                   coul_ijkb = eri_eigenstate_i(jstate,kstate,bstate,abispin) 
                 endif
               endif

               do iomegai=-nomegai,nomegai
                 omega = energy(bstate,abispin) + omegai(iomegai)
  
                 fact_real = REAL( fact_occ1 / (omega-ei+ej-ek+ieta) + fact_occ2 / (omega-ei+ej-ek-ieta) , dp)
                 fact_energy = REAL( fact_occ1 / (energy(astate,abispin)-ei+ej-ek+ieta) , dp )
  
                 selfenergy_ring(iomegai,astate,bstate2,abispin) = selfenergy_ring(iomegai,astate,bstate2,abispin) &
                          + fact_real * coul_iakj * coul_ibjk * spin_fact

                 if(iomegai==0 .AND. astate==bstate .AND. occupation(astate,abispin)>completely_empty) then
                   emp2_ring = emp2_ring + occupation(astate,abispin) &
                                         * fact_energy * coul_iakj * coul_ibjk * spin_fact
                 endif
  
                 if( abispin == jkspin ) then

                   selfenergy_sox(iomegai,astate,bstate2,abispin) = selfenergy_sox(iomegai,astate,bstate2,abispin) &
                            - fact_real * coul_iakj * coul_ijkb

                   if(iomegai==0 .AND. astate==bstate .AND. occupation(astate,abispin)>completely_empty) then
                     emp2_sox = emp2_sox - occupation(astate,abispin) &
                               * fact_energy * coul_iakj * coul_ijkb
                   endif

                 endif
    
    
               enddo ! iomega
    
             enddo
           enddo
         enddo
       enddo 
     enddo
   enddo 
 enddo ! abispin

 emp2_ring = 0.5_dp * emp2_ring
 emp2_sox  = 0.5_dp * emp2_sox
 emp2 = emp2_ring + emp2_sox
 write(stdout,'(/,a)')       ' MP2 Energy'
 write(stdout,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
 write(stdout,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
 write(stdout,'(a,f14.8,/)') ' MP2 correlation :',emp2


 if(method==perturbative) then

   if( print_sigma_) then
     call write_selfenergy_omega('selfenergy_mp2',nstate,energy,exchange_m_vxc_diag,SIZE(omegai),omegai, &
                                 nsemin,nsemax,selfenergy_ring(:,:,1,:)+selfenergy_sox(:,:,1,:))
   endif


   write(stdout,*) '=============================='
   write(stdout,*) ' selfenergy RING + SOX'
   write(stdout,'(a)') ' #         Energies           Sigx-Vxc       one-ring              SOX             MP2                Z       eMP2_Z       eMP2_qp'

   do astate=nsemin,nsemax
     zz(:) = ( selfenergy_ring(1,astate,1,:)+selfenergy_sox(1,astate,1,:) &
                 - selfenergy_ring(-1,astate,1,:)-selfenergy_sox(-1,astate,1,:) ) / REAL( omegai(1)-omegai(-1) ,dp)
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )

     energy_qp_z(astate,:) = energy(astate,:) + zz(:) * ( selfenergy_ring(0,astate,1,:) + selfenergy_sox(0,astate,1,:) + exchange_m_vxc_diag(astate,:) )

     allocate(sigma_xc_m_vxc_diag(-nomegai:nomegai))
     do ispin=1,nspin
       sigma_xc_m_vxc_diag(:) = selfenergy_ring(:,astate,1,ispin) + selfenergy_sox(:,astate,1,ispin) + exchange_m_vxc_diag(astate,ispin)
       energy_qp_omega(ispin) = find_fixed_point(nomegai,omegai,sigma_xc_m_vxc_diag) + energy(astate,ispin)
     enddo
     deallocate(sigma_xc_m_vxc_diag)
     energy_qp_new(astate,:) = energy_qp_omega(:) 


     write(stdout,'(i4,20(3x,f14.5))') astate,energy(astate,:)*Ha_eV,&
           exchange_m_vxc_diag(astate,:)*Ha_eV,&
           selfenergy_ring(0,astate,1,:) * Ha_eV,&
           selfenergy_sox (0,astate,1,:) * Ha_eV,&
         ( selfenergy_ring(0,astate,1,:)+selfenergy_sox(0,astate,1,:) ) * Ha_eV,&
         zz(:),&
         energy_qp_z(astate,:) * Ha_eV, &
         energy_qp_new(astate,:) * Ha_eV
   enddo
   write(stdout,*) '=============================='

   selfenergy(:,:,:) = 0.0_dp
   if(.NOT.ring_only) then
     do astate=nsemin,nsemax
       selfenergy(astate,1,:) = selfenergy_ring(0,astate,1,:) + selfenergy_sox(0,astate,1,:)
     enddo
   else
     do astate=nsemin,nsemax
       selfenergy(astate,1,:) = selfenergy_ring(0,astate,1,:)
     enddo
   endif

   call output_new_homolumo('MP2',nstate,occupation,energy_qp_z,nsemin,nsemax,ehomo,elumo)

 endif


 if(method==QS) then
 ! QS trick of Faleev-Kotani-vanSchilfgaarde
   do abispin=1,nspin
     if(.NOT.ring_only) then
       selfenergy_final(:,:,abispin) = 0.5_dp * ( selfenergy_ring(0,:,:,abispin) + selfenergy_sox(0,:,:,abispin) &
                                               +  TRANSPOSE( selfenergy_ring(0,:,:,abispin) + selfenergy_sox(0,:,:,abispin) ) )
     else
       write(stdout,*) 'ring_only'
       selfenergy_final(:,:,abispin) = 0.5_dp * ( selfenergy_ring(0,:,:,abispin)  &
                                               +  TRANSPOSE( selfenergy_ring(0,:,:,abispin) ) )
     endif
   enddo

   ! Transform the matrix elements back to the non interacting states
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do abispin=1,nspin
     selfenergy(:,:,abispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,abispin) ) , MATMUL( selfenergy_final(:,:,abispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,abispin)), s_matrix(:,:) ) ) )
   enddo
 endif

 if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
 deallocate(omegai)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)

 call stop_clock(timing_mp2_self)

end subroutine mp2_selfenergy


!=========================================================================
