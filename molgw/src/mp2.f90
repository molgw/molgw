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

 call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,ncore+1,nstate,ncore+1,nstate)

 ncore = ncoreg
 if(is_frozencore) then
   if( ncore == 0) ncore = atoms_core_states()
 endif

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
subroutine single_excitations(nstate,nbf,energy,occupation,c_matrix,fock_matrix)
 use m_definitions
 use m_timing
 use m_inputparam,only: nspin,spin_fact
 use m_scf,only: en
 use m_hamiltonian
 implicit none

 integer,intent(in)         :: nstate,nbf
 real(dp),intent(in)        :: energy(nstate,nspin),occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(inout)     :: fock_matrix(nbf,nbf,nspin)
!=====
 integer                    :: ispin
 integer                    :: istate,astate
!=====

 call start_clock(timing_single_excitation)

 !
 ! Rotate the Fock matrix to the eigenstate basis
 call matrix_basis_to_eigen(nbf,nstate,c_matrix,fock_matrix)


 en%se = 0.0_dp
 do ispin=1,nspin
   ! loop on occupied states
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     ! loop on virtual states
     do astate=1,nstate
       if( occupation(astate,ispin) > spin_fact - completely_empty ) cycle
       en%se = en%se + fock_matrix(istate,astate,ispin)**2 / ( energy(istate,ispin) - energy(astate,ispin) ) * spin_fact
     enddo
   enddo
 enddo

 call stop_clock(timing_single_excitation)

end subroutine single_excitations


!==================================================================
subroutine full_ci_2electrons_spin(print_wfn_,nstate,spinstate,basis,h_1e,c_matrix,nuc_nuc)
 use m_definitions
 use m_mpi
 use m_tools
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam,only: nspin,has_auxil_basis
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
 real(dp)              :: xxx(nx),y(ny),z(nz)
 real(dp)              :: wx(nx),wy(ny),wz(nz)
 real(dp)              :: norm
 integer               :: ibf_cart,li,ni,ni_cart,i_cart
 real(dp),allocatable  :: basis_function_r_cart(:)
 real(dp)              :: basis_function_r(basis%nbf)
 real(dp)              :: eval_wfn(nstate)
 real(dp),allocatable  :: eri_hf_i(:,:,:,:)
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,x,a,/)') 'Enter full CI subroutine'

 if( nspin /= 1) call die('CI is only implemented starting from spin-restricted SCF')

 if( .NOT. has_auxil_basis ) then
   allocate(eri_hf_i(nstate,nstate,nstate,nspin))
 else
   call calculate_eri_3center_eigen(basis%nbf,nstate,c_matrix,1,nstate,1,nstate)
 endif


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
!   write(stdout,'(i4,2x,20(x,f12.6))') iconf,h_1e_hf(istate,1:nstate)
! enddo
! write(stdout,*) '=========== full H ============== '
! do iconf=1,nconf
!   write(stdout,'(i4,2x,20(x,f12.6))') iconf,hamiltonian(iconf,1:nconf)
! enddo

 allocate(energy(nconf),eigenvector(nconf,nconf))
 !
 ! resize the matrices
 eigenvector(:,:) = hamiltonian(1:nconf,1:nconf)
 deallocate(hamiltonian)
 allocate(hamiltonian(nconf,nconf))
 hamiltonian(:,:) = eigenvector(:,:)

 if(.FALSE.) then 
   ! home made steepest descent algo
   write(stdout,*)
   write(stdout,*) 'hamiltonian too big to be diagonalized'
   write(stdout,*) 'Slow steepest descent algorithm'
   write(stdout,*) 'home cooked'
   allocate(test1(nconf),test2(nconf),hphi(nconf),gr1(nconf),gr2(nconf))
   test1(:) = 0.001
   test1(1) = 1. 
   test1(:) = test1(:) / SQRT( SUM( test1(:)**2 ) )
   delta = 0.10
  
   do iline=1,500
  
     hphi = MATMUL(eigenvector,test1)
     eigen = DOT_PRODUCT(test1,hphi)
     gr1(:) = 2. * ( hphi(:) - eigen * test1(:) )
     norm_gr1 = SQRT(DOT_PRODUCT(gr1,gr1))
!     write(stdout,*) 'x1',norm_gr1
  
!     write(stdout,*) 'guessed delta',delta
     test2(:) = test1(:) - delta * gr1(:) / norm_gr1
     test2(:) = test2(:) / SQRT( SUM( test2(:)**2 ) )
  
     hphi = MATMUL(eigenvector,test2)
     eigen = DOT_PRODUCT(test2,hphi)
     gr2(:) = 2. * ( hphi(:) - eigen * test2(:) )
!     write(stdout,*) 'x2',DOT_PRODUCT(gr2,gr1)/norm_gr1
     hessian = DOT_PRODUCT(gr1,gr2-gr1) / ( delta * norm_gr1 )
     delta = -norm_gr1 / hessian
!     write(stdout,*) 'optimal delta',delta
     test2(:) = test1(:) - delta * gr1(:) / norm_gr1
     test2(:) = test2(:) / SQRT( SUM( test2(:)**2 ) )
  
     if(MODULO(iline,20)==0) then
       write(stdout,*) 'diff',iline,eigen,SUM(ABS(test2(:)-test1(:)))/DBLE(nconf)
     endif
!     write(stdout,*) '==================='
     if( SUM(ABS(test2(:)-test1(:)))/DBLE(nconf) < 1.e-12_dp ) exit
     test1(:) = test2(:)
  
   enddo

   eigenvector(:,1) =test1(:)
   energy(1) = eigen
   
   deallocate(test1,test2,hphi,gr1,gr2)

 endif ! OLD stepest descent

 if( nconf>100 ) then
   write(stdout,*) 
   write(stdout,*) 'Davidson diago'
   write(stdout,*) 'trial vectors'

   bigm_max = neig + (ncycle-1) * neig
   bigm     = neig

   allocate(bb(nconf,bigm_max))
   allocate(atilde(bigm_max,bigm_max))
   allocate(ab(nconf,bigm_max))

   allocate(qq(nconf,neig))

   !
   ! Initialize with stupid coefficients
   bb(:,1:bigm)=0.001_dp
   do ieig=1,bigm
     bb(ieig,ieig) = 1.0_dp
   enddo
   do ieig=1,bigm
     ! orthogonalize to previous vectors
     do keig=1,ieig-1
       bb(:,ieig) = bb(:,ieig) - bb(:,keig) * DOT_PRODUCT( bb(:,ieig) , bb(:,keig) ) 
     enddo
     ! normalize
     bb(:,ieig) = bb(:,ieig) / NORM2( bb(:,ieig) )
   enddo


   ab(:,1:bigm) = MATMUL( hamiltonian(:,:) , bb(:,1:bigm) )


   do icycle=1,ncycle

     bigm = icycle * neig
     write(*,*) 'icycle bigm',icycle,bigm


     atilde(1:bigm,1:bigm) = MATMUL( TRANSPOSE(bb(:,1:bigm)) , ab(:,1:bigm) )
!     write(*,*) 
!     write(*,*) '==============================='
!     do ieig=1,bigm
!       write(*,'(i4,2x,20(2x,f14.6))') ieig,atilde(ieig,1:bigm)
!     enddo
!     write(*,*) '==============================='


     allocate(lambda(bigm),alphavec(bigm,bigm))
     call diagonalize(bigm,atilde(1:bigm,1:bigm),lambda,alphavec)

     write(stdout,*) 'icycle',icycle,lambda(1:neig)

     if( icycle == ncycle ) then
       energy(1:neig) = lambda(1:neig)
       eigenvector(:,1:neig) = MATMUL( bb(:,1:bigm) , alphavec(1:bigm,1:neig) )
       deallocate(lambda,alphavec)
       exit
     endif

     do ieig=1,neig

       qq(:,ieig) = MATMUL( ab(:,1:bigm) ,  alphavec(1:bigm,ieig) ) &
                     - lambda(ieig) * MATMUL( bb(:,1:bigm) , alphavec(1:bigm,ieig) )

       write(stdout,'(a,i4,x,i4,x,e12.4)') ' Residual norm for eigenvalue,cycle',ieig,icycle,NORM2(qq(:,ieig))

       do iconf=1,nconf
         bb(iconf,bigm+ieig) = qq(iconf,ieig) / ( lambda(ieig) - hamiltonian(iconf,iconf) )
       enddo
     enddo

     !
     ! orthogonalize to all previous
     do ieig=bigm+1,bigm+neig
       do keig=1,ieig-1
         bb(:,ieig) = bb(:,ieig) - bb(:,keig) * DOT_PRODUCT( bb(:,ieig) , bb(:,keig) )
       enddo
       ! normalize
       bb(:,ieig) = bb(:,ieig) / NORM2( bb(:,ieig) )
     enddo


     ab(:,bigm+1:bigm+neig) = MATMUL( hamiltonian(:,:) , bb(:,bigm+1:bigm+neig) )



     deallocate(lambda,alphavec)

   enddo ! icycle

   deallocate(ab,atilde)

   

   write(stdout,*) 
   write(stdout,*) 'Davidson diago DONE'
   write(stdout,*) energy(1:min(neig,nconf))

   deallocate(bb,qq)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Davidson   ',1,eigenvector(1:min(20,nconf),1)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Davidson   ',2,eigenvector(1:min(20,nconf),2)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Davidson   ',3,eigenvector(1:min(20,nconf),3)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Davidson   ',4,eigenvector(1:min(20,nconf),4)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Davidson   ',5,eigenvector(1:min(20,nconf),5)

 endif


 if( nconf>800 ) then
   write(stdout,*) 
   write(stdout,*) 'Davidson diago'
   write(stdout,*) 'trial vectors'

   allocate(bb_s(nconf,neig+ncycle*nblock))
   allocate(qq(nconf,nblock))
   !
   ! Initialize with stupid coefficients
   bb_s(:,1:neig)=0.001_dp
   do ieig=1,neig
     bb_s(ieig,ieig) = 1.0_dp
     ! orthogonalize to previous vectors
     do keig=1,ieig-1
       bb_s(:,ieig) = bb_s(:,ieig) - bb_s(:,keig) * DOT_PRODUCT( bb_s(:,ieig) , bb_s(:,keig) )
     enddo
     ! normalize
     bb_s(:,ieig) = bb_s(:,ieig) / NORM2( bb_s(:,ieig) )
   enddo

   do jeig=1,neig,nblock

     do icycle=1,ncycle

       neigc = neig + ( icycle - 1 ) * nblock

       allocate(bb(nconf,neigc))
       bb(:,:) = bb_s(:,1:neigc)

       allocate(ab(nconf,neigc),atilde(neigc,neigc))

       ab(:,:) = MATMUL( hamiltonian(:,:) , bb(:,:) )

       atilde(:,:) = MATMUL( TRANSPOSE(bb(:,:)) , ab(:,:) )

       allocate(lambda(neigc),alphavec(neigc,neigc))
       call diagonalize(neigc,atilde,lambda,alphavec)
       write(stdout,*) 'icycle',icycle,lambda(1:neig)

       do iblock=1,nblock
         qq(:,iblock) = MATMUL( ab ,  alphavec(:,jeig+iblock-1) ) &
                 - lambda(jeig+iblock-1) * MATMUL ( bb , alphavec(:,jeig+iblock-1) )

         write(stdout,'(a,i4,x,i4,x,e12.4)') ' Residual norm for eigenvalue,cycle',&
                       icycle,jeig+iblock-1,NORM2(qq(:,iblock))

         do iconf=1,nconf
           qq(iconf,iblock) = qq(iconf,iblock) / ( lambda(jeig+iblock-1) - hamiltonian(iconf,iconf) )
         enddo
         ! orthogonalize
         do ieig=1,neigc
           qq(:,iblock) = qq(:,iblock) - bb(:,ieig) * DOT_PRODUCT( bb(:,ieig) , qq(:,iblock) )
         enddo
         do jblock=1,iblock-1
           qq(:,iblock) = qq(:,iblock) - qq(:,jblock) * DOT_PRODUCT( qq(:,jblock) , qq(:,iblock) )
         enddo
         qq(:,iblock) = qq(:,iblock) / NORM2( qq(:,iblock) )

       enddo

       if(icycle<ncycle) then
         do iblock=1,nblock
           bb_s(:,neigc+iblock) = qq(:,iblock)
         enddo
       else
         do iblock=1,nblock
           do ieig=jeig,neig
             bb_s(:,ieig+iblock-1) = MATMUL( bb , alphavec(:,ieig+iblock-1) )
           enddo
           energy(jeig+iblock-1) = lambda(jeig+iblock-1)
         enddo
       endif
       deallocate(ab,atilde,lambda,alphavec,bb)

     enddo ! icycle
     write(stdout,*) 
   enddo ! jeig

   write(stdout,*) 'diago DONE'
   write(stdout,*) energy(1:min(neig,nconf))
   eigenvector(:,1:neig) = bb_s(:,1:neig)

   deallocate(bb_s,qq)

  else
   ! full LAPACK diago
   write(stdout,*) 
   write(stdout,*) 'starting the diago'
   call diagonalize(nconf,hamiltonian,energy,eigenvector)
   write(stdout,*) 'Full diago DONE'
   write(stdout,*) energy(1:min(neig,nconf))
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Full diago ',1,eigenvector(1:min(20,nconf),1)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Full diago ',2,eigenvector(1:min(20,nconf),2)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Full diago ',3,eigenvector(1:min(20,nconf),3)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Full diago ',4,eigenvector(1:min(20,nconf),4)
   write(stdout,'(a,i4,2x,20(x,f7.4))') ' Full diago ',5,eigenvector(1:min(20,nconf),5)

 endif

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
  

     !
     ! First precalculate all the needed basis function evaluations at point rr
     !
     ibf_cart = 1
     ibf      = 1
     do while(ibf_cart<=basis%nbf_cart)
       li      = basis%bf(ibf_cart)%am
       ni_cart = number_basis_function_am('CART',li)
       ni      = number_basis_function_am(basis%gaussian_type,li)
  
       allocate(basis_function_r_cart(ni_cart))
  
       do i_cart=1,ni_cart
         basis_function_r_cart(i_cart) = eval_basis_function(basis%bf(ibf_cart+i_cart-1),rr)
       enddo
       basis_function_r(ibf:ibf+ni-1) = MATMUL(  basis_function_r_cart(:) , cart_to_pure(li)%matrix(:,:) )
       deallocate(basis_function_r_cart)
  
       ibf      = ibf      + ni
       ibf_cart = ibf_cart + ni_cart
     enddo
     !
     ! Precalculation done!
     !

  
     eval_wfn(:) = 0.0_dp
     do istate=1,nstate
       do ibf=1,basis%nbf
         eval_wfn(istate) = eval_wfn(istate) + c_matrix(ibf,istate,1) * basis_function_r(ibf)
       enddo
     enddo
    
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
  !   rhor_t(ix)=rhor_t(ix)+ eval_basis_function(basis%bf(1),rr)**2 * wy(iy) * wz(iz)
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

