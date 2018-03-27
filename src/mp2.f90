!==================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! - MP2 total energy with or without Resolution-of-Identity
! - Single excitation contribution to total energy
! - Full CI for 2 electrons
!=========================================================================
subroutine mp2_energy_ri(nstate,basis,occupation,energy,c_matrix,emp2)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)       :: emp2
!====
 integer                    :: astate,bstate,istate,jstate
 integer                    :: iaspin,jbspin
 real(dp)                   :: energy_denom
 real(dp)                   :: tmp_iajb,tmp_ibja
 real(dp)                   :: contrib1,contrib2
 real(dp)                   :: fact
 integer                    :: nocc(nspin)
 integer                    :: ncore,nstate_mp2
!=====

 call start_clock(timing_mp2_energy)

 write(stdout,'(/,a)') ' RI-MP2 correlation calculation'


 ncore = ncoreg
 if(is_frozencore) then
   if( ncore == 0) ncore = atoms_core_states()
 endif

 call calculate_eri_3center_eigen(c_matrix,ncore+1,nstate,ncore+1,nstate)

 nstate_mp2 = MIN( nvirtualg-1, nstate )



 emp2 = 0.0_dp
 contrib1 = 0.0_dp
 contrib2 = 0.0_dp


 do iaspin=1,nspin
   !
   ! First, set up the list of occupied states
   nocc(iaspin) = ncore
   do istate=ncore+1,nstate
     if( occupation(istate,iaspin) < completely_empty ) cycle
     nocc(iaspin) = istate
   enddo
 enddo


 do iaspin=1,nspin

   do istate=ncore+1,nocc(iaspin)


     write(stdout,'(i4,2x,i4,a,i4)') iaspin,istate-ncore,' / ',nocc-ncore

     do jbspin=1,nspin

       do jstate=ncore+1,nocc(jbspin)
       
         do astate=ncore+1,nstate_mp2
           if( occupation(astate,iaspin) > spin_fact - completely_empty ) cycle

           do bstate=ncore+1,nstate_mp2
             if( occupation(bstate,jbspin) > spin_fact - completely_empty ) cycle

             fact =  occupation(istate,iaspin) * ( spin_fact - occupation(astate,iaspin) ) &
&                   *occupation(jstate,jbspin) * ( spin_fact - occupation(bstate,jbspin) ) / spin_fact**2

             energy_denom = energy(istate,iaspin) + energy(jstate,jbspin) &
&                                    - energy(astate,iaspin) - energy(bstate,jbspin)
             ! Avoid the zero denominators
             if( ABS(energy_denom) < 1.d-18) then
               write(stdout,*) 'you skipped something'
               cycle
             endif

             energy_denom =  fact / energy_denom 

             tmp_iajb = eri_eigen_ri(istate,astate,iaspin,jstate,bstate,jbspin)

             contrib1 = contrib1 + 0.5_dp * energy_denom * tmp_iajb**2

             if(iaspin==jbspin) then
               tmp_ibja = eri_eigen_ri(istate,bstate,iaspin,jstate,astate,jbspin)
               contrib2 = contrib2 - 0.5_dp * energy_denom * tmp_iajb*tmp_ibja / spin_fact
             endif

           enddo
         enddo
       enddo
     enddo !jbspin

   enddo ! istate

 enddo !iaspin


 emp2 = contrib1 + contrib2
 write(stdout,'(/,a)')       ' MP2 contributions'
 write(stdout,'(a,f16.10)')   ' 2-ring diagram  :',contrib1
 write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib2
 write(stdout,'(a,f16.10,/)') ' MP2 correlation :',emp2

 call destroy_eri_3center_eigen()
 call stop_clock(timing_mp2_energy)

end subroutine mp2_energy_ri


!=========================================================================
subroutine mp3_energy_ri(nstate,basis,occupation,energy,c_matrix,emp3)
 use m_definitions
 use m_mpi
 use m_cart_to_pure
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)       :: emp3
!====
 integer                    :: astate,bstate,cstate,dstate,istate,jstate,kstate,lstate
 integer                    :: iaspin,jbspin,kcspin
 integer                    :: cspin,dspin,kspin,lspin
 real(dp)                   :: energy_denom
 real(dp)                   :: tmp_iajb,tmp_ibja
 real(dp)                   :: contrib1,contrib2,contrib3
 real(dp)                   :: fact
 real(dp)                   :: denom1,denom2,numer1,numer2
 integer                    :: nocc(nspin)
 integer                    :: ncore,nstate_mp3
 real(dp)                   :: t_ijab_tilde,x_ijab,t_ijcd,t_klab,t_kjac,t_kiac,t_ikac
!=====

 call start_clock(timing_mp2_energy)

 write(stdout,'(/,a)') ' RI-MP3 correlation calculation'

 if( nspin > 1 ) call die('MP3 not implemented for unrestricted calculations')

 ncore = ncoreg
 if(is_frozencore) then
   if( ncore == 0) ncore = atoms_core_states()
 endif

 call calculate_eri_3center_eigen(c_matrix,ncore+1,nstate,ncore+1,nstate)

 nstate_mp3 = MIN( nvirtualg-1, nstate )



 emp3 = 0.0_dp
 contrib1 = 0.0_dp
 contrib2 = 0.0_dp
 contrib3 = 0.0_dp


 do iaspin=1,nspin
   !
   ! First, set up the list of occupied states
   nocc(iaspin) = ncore
   do istate=ncore+1,nstate
     if( occupation(istate,iaspin) < completely_empty ) cycle
     nocc(iaspin) = istate
   enddo
 enddo


#if 0
 do iaspin=1,nspin
   do istate=ncore+1,nocc(iaspin)
     do astate=nocc(iaspin)+1,nstate_mp3

       do jbspin=1,nspin
         do jstate=ncore+1,nocc(iaspin)
           do bstate=nocc(iaspin)+1,nstate_mp3

             denom1 = energy(istate,iaspin) + energy(jstate,jbspin) &
                                   - energy(astate,iaspin) - energy(bstate,jbspin)
             numer1 = eri_eigen_ri(astate,istate,iaspin,bstate,jstate,jbspin)
             if( iaspin == jbspin ) numer1 = numer1 - eri_eigen_ri(astate,jstate,iaspin,bstate,istate,jbspin) / spin_fact

             !
             ! Contrib1 
             cspin = iaspin
             dspin = jbspin
             do cstate=nocc(iaspin)+1,nstate_mp3
               do dstate=nocc(iaspin)+1,nstate_mp3

                 denom2 = energy(istate,iaspin) + energy(jstate,jbspin) &
                                       - energy(cstate,iaspin) - energy(dstate,jbspin)
                 numer2 = eri_eigen_ri(cstate,istate,iaspin,dstate,jstate,jbspin)
                 if( cspin == dspin ) numer2 = numer2 - eri_eigen_ri(cstate,jstate,iaspin,dstate,istate,jbspin) / spin_fact

                 contrib1 = contrib1 + 0.125_dp * numer1 / denom1 * numer2 / denom2 &
                                      * eri_eigen_ri(astate,cstate,iaspin,bstate,dstate,jbspin)
                 if( cspin == dspin ) &
                   contrib1 = contrib1 - 0.125_dp * numer1 / denom1 * numer2 / denom2 &
                                      * eri_eigen_ri(astate,dstate,iaspin,bstate,cstate,jbspin) / spin_fact 

               enddo
             enddo

             !
             ! Contrib2 
             kspin = iaspin
             lspin = jbspin
             do kstate=ncore+1,nocc(iaspin)
               do lstate=ncore+1,nocc(iaspin)

                 denom2 = energy(kstate,iaspin) + energy(lstate,jbspin) &
                                       - energy(astate,iaspin) - energy(bstate,jbspin)
                 numer2 = eri_eigen_ri(astate,kstate,iaspin,bstate,lstate,jbspin)
                 if( kspin == lspin ) numer2 = numer2 - eri_eigen_ri(astate,lstate,iaspin,bstate,kstate,jbspin) / spin_fact

                 contrib2 = contrib2 + 0.125_dp * numer1 / denom1 * numer2 / denom2 &
                                      * eri_eigen_ri(kstate,istate,iaspin,lstate,jstate,jbspin) 
                 if( kspin == lspin ) &
                   contrib2 = contrib2 - 0.125_dp * numer1 / denom1 * numer2 / denom2 &
                                      * eri_eigen_ri(kstate,jstate,iaspin,lstate,istate,jbspin) / spin_fact 

               enddo
             enddo

             !
             ! Contrib3
             do kcspin=1,nspin
               do kstate=ncore+1,nocc(iaspin)
                 do cstate=nocc(iaspin)+1,nstate_mp3

                   denom2 = energy(kstate,iaspin) + energy(jstate,jbspin) &
                                         - energy(cstate,iaspin) - energy(bstate,jbspin)
                   numer2 = eri_eigen_ri(cstate,kstate,kcspin,bstate,jstate,jbspin)
                   if( kcspin == jbspin ) numer2 = numer2 - eri_eigen_ri(cstate,jstate,kcspin,kstate,bstate,jbspin) / spin_fact

                   contrib3 = contrib3 + numer1 / denom1 * numer2 / denom2 &
                                        * eri_eigen_ri(astate,istate,iaspin,kstate,cstate,kcspin)
                   if( iaspin == kcspin ) &
                     contrib3 = contrib3 - numer1 / denom1 * numer2 / denom2 &
                                            * eri_eigen_ri(astate,cstate,iaspin,kstate,istate,kcspin) / spin_fact


                 enddo
               enddo
             enddo



       enddo
     enddo
   enddo
 enddo
#else

 ! From Helgaker's book
 write(stdout,*) 'From Helgaker book'
 iaspin=1
 jbspin=1

 do istate=ncore+1,nocc(iaspin)
   do astate=nocc(iaspin)+1,nstate_mp3
     do jstate=ncore+1,nocc(iaspin)
       do bstate=nocc(iaspin)+1,nstate_mp3

         t_ijab_tilde = - 2.0_dp * ( 2.0_dp * eri_eigen_ri(istate,astate,iaspin,jstate,bstate,jbspin)  &
                                     - eri_eigen_ri(istate,bstate,iaspin,jstate,astate,jbspin) ) &
                                       / ( energy(astate,iaspin) + energy(bstate,jbspin) - energy(istate,iaspin) - energy(jstate,jbspin) )

         !FIXME: Commenting this is inconsistent with Helgaker however yields the correct numerical results
         !if( istate == jstate .AND. astate == bstate ) t_ijab_tilde = t_ijab_tilde * 2.0_dp

         x_ijab = 0.0_dp
         !
         ! Contrib1 
         do cstate=nocc(iaspin)+1,nstate_mp3
           do dstate=nocc(iaspin)+1,nstate_mp3

             t_ijcd = - eri_eigen_ri(istate,cstate,iaspin,jstate,dstate,jbspin)  &
                          / ( energy(cstate,iaspin) + energy(dstate,jbspin) - energy(istate,iaspin) - energy(jstate,jbspin) )

             x_ijab = x_ijab + 0.5_dp * eri_eigen_ri(astate,cstate,iaspin,bstate,dstate,jbspin) * t_ijcd

           enddo
         enddo

         !
         ! Contrib2 
         do kstate=ncore+1,nocc(iaspin)
           do lstate=ncore+1,nocc(iaspin)

             t_klab = - eri_eigen_ri(kstate,astate,iaspin,lstate,bstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(bstate,jbspin) - energy(kstate,iaspin) - energy(lstate,jbspin) )

             x_ijab = x_ijab + 0.5_dp * eri_eigen_ri(kstate,istate,iaspin,lstate,jstate,jbspin) * t_klab

           enddo
         enddo

         !
         ! Contrib3
         do kstate=ncore+1,nocc(iaspin)
           do cstate=nocc(iaspin)+1,nstate_mp3

             t_ikac = - eri_eigen_ri(istate,astate,iaspin,kstate,cstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(cstate,jbspin) - energy(istate,iaspin) - energy(kstate,jbspin) )

             x_ijab = x_ijab + ( 2.0_dp * eri_eigen_ri(bstate,jstate,iaspin,kstate,cstate,jbspin) &
                                  - eri_eigen_ri(bstate,cstate,iaspin,kstate,jstate,jbspin) ) * t_ikac


             t_kjac = - eri_eigen_ri(kstate,astate,iaspin,jstate,cstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(cstate,jbspin) - energy(kstate,iaspin) - energy(jstate,jbspin) )
             x_ijab = x_ijab - eri_eigen_ri(bstate,cstate,iaspin,kstate,istate,jbspin) * t_kjac

             t_kiac = - eri_eigen_ri(kstate,astate,iaspin,istate,cstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(cstate,jbspin) - energy(kstate,iaspin) - energy(istate,jbspin) )
             x_ijab = x_ijab - eri_eigen_ri(bstate,jstate,iaspin,kstate,cstate,jbspin) * t_kiac


           enddo
         enddo

         contrib1 = contrib1 + t_ijab_tilde * x_ijab

       enddo
     enddo
   enddo
 enddo
#endif



 emp3 = contrib1 + contrib2 + contrib3
! write(stdout,'(/,a)')       ' MP3 contributions'
! write(stdout,'(a,f16.10)')   ' 2-ring diagram  :',contrib1
! write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib2
! write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib3
! write(stdout,'(a,f16.10,/)') ' MP3 correlation :',emp3

 call destroy_eri_3center_eigen()
 call stop_clock(timing_mp2_energy)

end subroutine mp3_energy_ri


!==================================================================
subroutine mp2_energy(nstate,basis,occupation,c_matrix,energy,emp2)
 use m_definitions
 use m_mpi
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam,only: nspin,spin_fact,ncoreg
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)       :: emp2
!=====
 integer                    :: astate,bstate,istate,jstate
 integer                    :: ibf,jbf,abf,bbf,iaspin,jbspin
 real(dp)                   :: energy_denom
 real(dp)                   :: tmp_ixjx(basis%nbf,basis%nbf)
 real(dp)                   :: tmp_iajx(basis%nbf),tmp_ixja(basis%nbf)
 real(dp)                   :: tmp_iajb,tmp_ibja
 real(dp)                   :: contrib1,contrib2
 real(dp)                   :: fact
 real(dp),allocatable       :: tmp_ixxx(:,:,:)
 integer                    :: nocc
!=====

 call start_clock(timing_mp2_energy)

 write(stdout,*) 'starting the MP2 calculation'


 emp2 = 0.0_dp
 contrib1 = 0.0_dp
 contrib2 = 0.0_dp


 allocate(tmp_ixxx(basis%nbf,basis%nbf,basis%nbf))

 do iaspin=1,nspin
   !
   ! First, set up the list of occupied states
   nocc = ncoreg
   do istate=ncoreg+1,nstate
     if( occupation(istate,iaspin) < completely_empty ) cycle
     nocc = istate
   enddo


   do istate=ncoreg+1,nocc


     write(stdout,'(i4,2x,i4,a,i4)') iaspin,istate-ncoreg,' / ',nocc

     tmp_ixxx(:,:,:) = 0.0_dp
     do bbf=1,basis%nbf
       do jbf=1,basis%nbf
         if( negligible_basispair(jbf,bbf) ) cycle
         do abf=1,basis%nbf
           do ibf=1,basis%nbf
             if( negligible_basispair(ibf,abf) ) cycle
             tmp_ixxx(abf,jbf,bbf) = tmp_ixxx(abf,jbf,bbf) &
&               + c_matrix(ibf,istate,iaspin) * eri(ibf,abf,jbf,bbf)
           enddo
         enddo
       enddo
     enddo

     do jbspin=1,nspin
       do jstate=ncoreg+1,nstate
         if( occupation(jstate,jbspin) < completely_empty ) cycle
       
         tmp_ixjx(:,:) = 0.0_dp
         do bbf=1,basis%nbf
           do jbf=1,basis%nbf
             do abf=1,basis%nbf
               tmp_ixjx(abf,bbf) = tmp_ixjx(abf,bbf) + c_matrix(jbf,jstate,jbspin) * tmp_ixxx(abf,jbf,bbf)
             enddo
           enddo
         enddo

         do astate=1,nstate
           if( occupation(astate,iaspin) > spin_fact - completely_empty ) cycle

           tmp_iajx(:) = 0.0_dp
           do bbf=1,basis%nbf
             do abf=1,basis%nbf
               tmp_iajx(bbf) = tmp_iajx(bbf) + c_matrix(abf,astate,iaspin) * tmp_ixjx(abf,bbf)
             enddo
           enddo

           if(iaspin==jbspin) then 
             tmp_ixja(:) = 0.0_dp
             do abf=1,basis%nbf
               do bbf=1,basis%nbf
                 tmp_ixja(bbf) = tmp_ixja(bbf) + c_matrix(abf,astate,iaspin) * tmp_ixjx(bbf,abf)
               enddo
             enddo
           endif

           do bstate=1,nstate
             if( occupation(bstate,jbspin) > spin_fact - completely_empty ) cycle

             fact =  occupation(istate,iaspin) * ( spin_fact - occupation(astate,iaspin) ) &
&                   *occupation(jstate,jbspin) * ( spin_fact - occupation(bstate,jbspin) ) / spin_fact**2

             energy_denom = energy(istate,iaspin) + energy(jstate,jbspin) &
&                                    - energy(astate,iaspin) - energy(bstate,jbspin)
             ! Avoid the zero denominators
             if( ABS(energy_denom) < 1.d-18) then
               write(stdout,*) 'you skipped something'
               cycle
             endif

             energy_denom =  fact / energy_denom 

             tmp_iajb = SUM( tmp_iajx(:) * c_matrix(:,bstate,jbspin) )

             contrib1 = contrib1 + 0.5_dp * energy_denom * tmp_iajb**2 

             if(iaspin==jbspin) then
               tmp_ibja = SUM( tmp_ixja(:) * c_matrix(:,bstate,jbspin) )
               contrib2 = contrib2 - 0.5_dp * energy_denom * tmp_iajb*tmp_ibja / spin_fact
             endif

           enddo
         enddo
       enddo
     enddo !jbspin

   enddo ! istate

 enddo !iaspin

 deallocate(tmp_ixxx)


 emp2 = contrib1 + contrib2
 write(stdout,'(/,a)')       ' MP2 contributions'
 write(stdout,'(a,f16.10)')   ' 2-ring diagram  :',contrib1
 write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib2
 write(stdout,'(a,f16.10,/)') ' MP2 correlation :',emp2

 call stop_clock(timing_mp2_energy)

end subroutine mp2_energy


!==================================================================
subroutine single_excitations(nstate,nbf,energy,occupation,c_matrix,fock_matrix,energy_se)
 use m_definitions
 use m_timing
 use m_inputparam,only: nspin,spin_fact
 use m_hamiltonian
 implicit none

 integer,intent(in)         :: nstate,nbf
 real(dp),intent(in)        :: energy(nstate,nspin),occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)        :: fock_matrix(nbf,nbf,nspin)
 real(dp),intent(out)       :: energy_se
!=====
 integer                    :: ier
 integer                    :: ispin
 integer                    :: istate,astate
 real(dp),allocatable       :: fock_matrix_eigen(:,:,:)
!=====

 call start_clock(timing_single_excitation)

 energy_se = 0.0_dp

 allocate(fock_matrix_eigen(nbf,nbf,nspin),stat=ier)
 ier = ABS(ier)
 call xmax_world(ier)
 if( ier /= 0 ) then
   write(stdout,*) 'Skip this step that is too memory demanding'
   return
 endif
 
 fock_matrix_eigen(:,:,:) = fock_matrix(:,:,:)
 !
 ! Rotate the Fock matrix to the eigenstate basis
 call matrix_basis_to_eigen(c_matrix,fock_matrix_eigen)


 do ispin=1,nspin
   ! loop on occupied states
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     ! loop on virtual states
     do astate=1,nstate
       if( occupation(astate,ispin) > spin_fact - completely_empty ) cycle
       energy_se = energy_se + fock_matrix_eigen(istate,astate,ispin)**2 / ( energy(istate,ispin) - energy(astate,ispin) ) * spin_fact
     enddo
   enddo
 enddo

 deallocate(fock_matrix_eigen)

 call stop_clock(timing_single_excitation)

end subroutine single_excitations


!==================================================================
subroutine full_ci_2electrons_spin(print_wfn_,nstate,spinstate,basis,h_1e,c_matrix,nuc_nuc)
 use m_definitions
 use m_mpi
 use m_warning
 use m_tools
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam,only: nspin,has_auxil_basis
 use m_dft_grid
 implicit none
!
 integer,parameter :: cip=dp
 integer,parameter :: nx=4000
!
 logical,intent(in)         :: print_wfn_
 integer,intent(in)         :: spinstate,nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: h_1e(basis%nbf,basis%nbf),c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: nuc_nuc
!=====
 integer,parameter     :: neig=6
 integer,parameter     :: nblock=1
 integer,parameter     :: ncycle=12
 integer               :: gt
 integer               :: bigm,bigm_max
 integer               :: ieig,jeig,keig,neigc,icycle,iblock,jblock
 real(dp),allocatable  :: bb(:,:),qq(:,:),atilde(:,:),ab(:,:)
 real(dp),allocatable  :: bb_s(:,:),atilde_s(:,:),ab_s(:,:)
 real(dp),allocatable  :: lambda(:),alphavec(:,:)
 real(dp)              :: h_1e_hf(nstate,nstate)
 integer               :: nconf,iconf,jconf,kconf
 integer               :: ibf,jbf
 integer               :: istate,jstate,kstate,lstate
 integer               :: istate1,istate2,jstate1,jstate2,ispin1,ispin2,jspin1,jspin2
 real(cip),allocatable :: hamiltonian(:,:)
 real(cip),allocatable :: energy(:),eigenvector(:,:)
 real(cip),allocatable :: test1(:),test2(:),test3(:),hphi(:),gr1(:),gr2(:)
 real(cip)             :: delta,eigen,hessian,norm_gr1
 integer               :: iline,ix
 real(dp)              :: rhor(nx),rhor_hf(nx),rr(3)
 real(dp)              :: rhor_t(nx)
 integer,parameter     :: ny=nx,nz=nx
 integer               :: iy,iz
 real(dp)              :: wx(nx),wy(ny),wz(nz)
 real(dp)              :: norm
 real(dp)              :: basis_function_r(basis%nbf)
 real(dp)              :: eval_wfn(nstate)
 real(dp),allocatable  :: eri_hf_i(:,:,:,:)
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,1x,a,/)') 'Enter full CI subroutine'

 if( nspin /= 1) call die('CI is only implemented starting from spin-restricted SCF')

 if( .NOT. has_auxil_basis ) then
   allocate(eri_hf_i(nstate,nstate,nstate,nspin))
 else
   call calculate_eri_3center_eigen(c_matrix)
 endif

 gt = get_gaussian_type_tag(basis%gaussian_type)

 write(stdout,*) 'Obtain the one-electron Hamiltonian in the HF basis'

 h_1e_hf(:,:) = MATMUL( TRANSPOSE(c_matrix(:,:,1)) , MATMUL( h_1e(:,:) , c_matrix(:,:,1) ) )

 select case(spinstate)
 case(0)
   write(stdout,*) 'calculate spin singlet state'
 case(1)
   write(stdout,*) 'calculate spin triplet state'
 case default
   call die('BUG: spin state not possible')
 end select

 nconf = ( 2 * nstate * (2 * nstate - 1) ) / 2
 write(stdout,*)
 write(stdout,*) 'CI matrix lower than',nconf,' x ',nconf
 allocate(hamiltonian(nconf,nconf))
 hamiltonian(:,:) = 0.0_dp
 do iconf=1,nconf
   hamiltonian(iconf,iconf) = nuc_nuc
 enddo

 iconf=0
 do istate1=1,nstate

   if( .NOT. has_auxil_basis ) then
     call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,istate1,1,eri_hf_i)
   endif

   do ispin1=1,2

     do istate2=istate1,nstate
       do ispin2=1,2
         if(istate1==istate2 .AND. (ispin2==1 .OR. ispin1==2) ) cycle
         !
         ! S^2 selection
!         if(ABS(ispin1-ispin2)==spinstate) cycle
         !
         ! for two electrons, the two electron wavefunction is even in space
         ! the two one-particle wavefunctions have then to have the parity 
!TODO         if(.NOT.symmetry(istate1,istate2)) cycle
         iconf=iconf+1

         jconf=0
         do jstate1=1,nstate
           do jspin1=1,2
             do jstate2=jstate1,nstate
               do jspin2=1,2
                 if(jstate1==jstate2 .AND. (jspin2==1 .OR. jspin1==2) ) cycle
                 !
                 ! S^2 selection
!                 if(ABS(jspin1-jspin2)==spinstate) cycle
                 !
                 ! for two electrons, the two electron wavefunction is even in space
                 ! the two one-particle wavefunctions have then to have the parity 
!TODO                 if(.NOT.symmetry(jstate1,jstate2)) cycle
                 jconf=jconf+1

!         write(stdout,'(10(i4,x))') jconf,jstate1,jspin1,jstate2,jspin2

                 if( istate2==jstate2 .AND. ispin2==jspin2 .AND. ispin1==jspin1 ) hamiltonian(iconf,jconf) = hamiltonian(iconf,jconf) + h_1e_hf(istate1,jstate1)
                 if( istate1==jstate1 .AND. ispin1==jspin1 .AND. ispin2==jspin2 ) hamiltonian(iconf,jconf) = hamiltonian(iconf,jconf) + h_1e_hf(istate2,jstate2)
                 if( istate2==jstate1 .AND. ispin2==jspin1 .AND. ispin1==jspin2 ) hamiltonian(iconf,jconf) = hamiltonian(iconf,jconf) - h_1e_hf(istate1,jstate2)
                 if( istate1==jstate2 .AND. ispin1==jspin2 .AND. ispin2==jspin1 ) hamiltonian(iconf,jconf) = hamiltonian(iconf,jconf) - h_1e_hf(istate2,jstate1)


                 !
                 ! Not so painful implementation of the determinant rules as shown in
                 ! p. 70 of "Modern Quantum Chemistry" by A. Szabo and N. S. Ostlung

                 if( has_auxil_basis ) then
                   if( ispin1==jspin1 .AND. ispin2==jspin2 ) &
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_eigen_ri(istate1,jstate1,1,istate2,jstate2,1)
                   if( ispin1==jspin2 .AND. ispin2==jspin1 ) &
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) - eri_eigen_ri(istate1,jstate2,1,istate2,jstate1,1)
                 else
                   if( ispin1==jspin1 .AND. ispin2==jspin2 ) &
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate1,istate2,jstate2,1)
                   if( ispin1==jspin2 .AND. ispin2==jspin1 ) &
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) - eri_hf_i(jstate2,istate2,jstate1,1)
                 endif

               enddo
             enddo
           enddo
         enddo

       enddo
     enddo
   enddo
 enddo

 ! I dont need the integrals anymore
 if( .NOT. has_auxil_basis ) then
   deallocate(eri_hf_i)
 else
   call destroy_eri_3center_eigen()
 endif

 ! Adjust the real size of the CI hamiltonian
 nconf=iconf
 write(stdout,*)
 write(stdout,*) 'CI matrix finally is',nconf,' x ',nconf

 write(stdout,'(a,f16.10)') ' Single determinant ground state energy (Ha): ',hamiltonian(1,1)
! write(stdout,*) '=========== H_1e ============== '
! do istate=1,nstate
!   write(stdout,'(i4,2x,20(1x,f12.6))') iconf,h_1e_hf(istate,1:nstate)
! enddo
! write(stdout,*) '=========== full H ============== '
! do iconf=1,nconf
!   write(stdout,'(i4,2x,20(1x,f12.6))') iconf,hamiltonian(iconf,1:nconf)
! enddo

 allocate(energy(nconf),eigenvector(nconf,nconf))
 !
 ! resize the matrices
 eigenvector(:,:) = hamiltonian(1:nconf,1:nconf)
 deallocate(hamiltonian)
 allocate(hamiltonian(nconf,nconf))
 hamiltonian(:,:) = eigenvector(:,:)


 ! full LAPACK diago
 write(stdout,*) 
 write(stdout,*) 'starting the diago'
 call diagonalize(nconf,hamiltonian,energy,eigenvector)
 write(stdout,*) 'Full diago DONE'
 write(stdout,*) energy(1:MIN(neig,nconf))
 write(stdout,'(a,i4,2x,20(1x,f7.4))') ' Full diago ',1,eigenvector(1:MIN(20,nconf),1)
 write(stdout,'(a,i4,2x,20(1x,f7.4))') ' Full diago ',2,eigenvector(1:MIN(20,nconf),2)
 write(stdout,'(a,i4,2x,20(1x,f7.4))') ' Full diago ',3,eigenvector(1:MIN(20,nconf),3)
 write(stdout,'(a,i4,2x,20(1x,f7.4))') ' Full diago ',4,eigenvector(1:MIN(20,nconf),4)
 write(stdout,'(a,i4,2x,20(1x,f7.4))') ' Full diago ',5,eigenvector(1:MIN(20,nconf),5)


 write(stdout,*)
 write(stdout,*) 'normalization',SUM(eigenvector(:,1)**2)
 write(stdout,*)
 write(stdout,'(a30,2x,f16.10)') 'CI ground-state energy (Ha):',energy(1)
 write(stdout,'(a30,2x,f16.10)') 'correlation energy (Ha):',energy(1)-hamiltonian(1,1)
 write(stdout,*)
  
 deallocate(hamiltonian)


 !
 ! Plot the ground state density if requested
 !
 if( print_wfn_ ) then
   write(stdout,*)
   write(stdout,*) 'calculate the density'
  
  
   rhor(:)=0.0_dp
   rhor_hf(:)=0.0_dp
   do ix=1,nx
     rr(1)= ( DBLE(ix-1)/DBLE(nx-1) - 0.5 ) * 10.0
     rr(2)= 0.0
     rr(3)= 0.0
  

     call calculate_basis_functions_r(basis,rr,basis_function_r)
  
     eval_wfn(:) = MATMUL( basis_function_r(:) ,  c_matrix(:,:,1) )

    
    do kconf=1,1
    
     iconf=0
     do istate1=1,nstate
       do ispin1=1,2
         do istate2=istate1,nstate
           do ispin2=1,2
             if(istate1==istate2 .AND. (ispin2==1 .OR. ispin1==2) ) cycle
  !           !
  !           ! S^2 selection
  !           if(ABS(ispin1-ispin2)==spinstate) cycle
             !
             ! for two electrons, the two electron wavefunction is even in space
             ! the two one-particle wavefunctions have then to have the parity 
    !TODO         if(.NOT.symmetry(istate1,istate2)) cycle
             iconf=iconf+1
    
             jconf=0
             do jstate1=1,nstate
               do jspin1=1,2
                 do jstate2=jstate1,nstate
                   do jspin2=1,2
                     if(jstate1==jstate2 .AND. (jspin2==1 .OR. jspin1==2) ) cycle
    !                 !
    !                 ! S^2 selection
    !                 if(ABS(jspin1-jspin2)==spinstate) cycle
                     !
                     ! for two electrons, the two electron wavefunction is even in space
                     ! the two one-particle wavefunctions have then to have the parity 
    !TODO                 if(.NOT.symmetry(jstate1,jstate2)) cycle
                     jconf=jconf+1
    
    
                       if( istate2==jstate2 .AND. ispin2==jspin2 .AND. ispin1==jspin1 ) then
                         rhor(ix) = rhor(ix)  &
                          + eigenvector(iconf,kconf) * eigenvector(jconf,kconf) &
                           * eval_wfn(istate1) * eval_wfn(jstate1) 
                       endif
                       if( istate1==jstate1 .AND. ispin1==jspin1 .AND. ispin2==jspin2 ) then 
                         rhor(ix) = rhor(ix)  &
                          + eigenvector(iconf,kconf) * eigenvector(jconf,kconf) &
                           * eval_wfn(istate2) * eval_wfn(jstate2) 
                       endif
                       if( istate2==jstate1 .AND. ispin2==jspin1 .AND. ispin1==jspin2 ) then
                         rhor(ix) = rhor(ix)  &
                          - eigenvector(iconf,kconf) * eigenvector(jconf,kconf) &
                           * eval_wfn(istate1) * eval_wfn(jstate2)
                       endif
                       if( istate1==jstate2 .AND. ispin1==jspin2 .AND. ispin2==jspin1 ) then
                         rhor(ix) = rhor(ix)  &
                          - eigenvector(iconf,kconf) * eigenvector(jconf,kconf) &
                           * eval_wfn(istate2) * eval_wfn(jstate1)
                       endif
  
                       !
                       ! HARTREE-FOCK PART
                       if( iconf==kconf .AND. jconf==kconf ) then
                         if( istate2==jstate2 .AND. ispin2==jspin2 .AND. ispin1==jspin1 ) then
                           rhor_hf(ix) = rhor_hf(ix)  &
                            + eval_wfn(istate1) * eval_wfn(jstate1) 
                         endif
                         if( istate1==jstate1 .AND. ispin1==jspin1 .AND. ispin2==jspin2 ) then 
                           rhor_hf(ix) = rhor_hf(ix)  &
                            + eval_wfn(istate2) * eval_wfn(jstate2) 
                         endif
                         if( istate2==jstate1 .AND. ispin2==jspin1 .AND. ispin1==jspin2 ) then
                           rhor_hf(ix) = rhor_hf(ix)  &
                            - eval_wfn(istate1) * eval_wfn(jstate2)
                         endif
                         if( istate1==jstate2 .AND. ispin1==jspin2 .AND. ispin2==jspin1 ) then
                           rhor_hf(ix) = rhor_hf(ix)  &
                            - eval_wfn(istate2) * eval_wfn(jstate1)
                         endif
                       endif
    
    
    
                   enddo
                 enddo
               enddo
             enddo
    
           enddo
         enddo
       enddo
     enddo
  
     enddo  !kconf
  
   enddo !ix
   rhor(:)    = rhor   (:) * 0.5_dp  ! divide by two to compare with phi^star * phi which is normalized to unity
   rhor_hf(:) = rhor_hf(:) * 0.5_dp  ! divide by two to compare with phi^star * phi which is normalized to unity
  
  ! write(stdout,*) 'NORM',norm / DBLE(nx*ny*nz)
  ! write(stdout,*) 'norm',SUM(rhor(:)*wx(:))
  
  ! rhor_t(:)=0.0_dp
  ! do ix=1,nx
  ! do iy=1,ny
  ! do iz=1,nz
  !   rr(1)= x(ix)
  !   rr(2)= y(iy)
  !   rr(3)= z(iz)
  !   rhor_t(ix)=rhor_t(ix)+ eval_basis_function(basis%bfc(1),rr)**2 * wy(iy) * wz(iz)
  ! enddo !iz
  ! enddo !iy
  ! enddo !ix
  ! write(stdout,*) 'norm',SUM(rhor_t(:)*wx(:))
  ! do ix=1,nx
  !   rr(1)= x(ix)
  !   write(11,'(5(e12.6,2x))') rr(1),rhor_t(ix)
  ! enddo
  
   do ix=1,nx
  !   rr(1)= x(ix)
     rr(1)= (DBLE(ix-1)/DBLE(nx-1)-0.5)*10.00
     rr(2)= 0.0
     rr(3)= 0.0
     write(10,'(5(e14.6,2x))') rr(1),rhor(ix),rhor_hf(ix)
   enddo

 endif

 ! finalize the CI calculation
 deallocate(energy,eigenvector)

 call stop_clock(timing_full_ci)

end subroutine full_ci_2electrons_spin


!==================================================================
