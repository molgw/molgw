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
 implicit none

 integer,intent(in)  :: method,nstate
 type(basis_set)     :: basis
 real(dp),intent(in) :: occupation(nstate,nspin),energy(nstate,nspin),exchange_m_vxc_diag(nstate,nspin)
 real(dp),intent(in) :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in) :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(out) :: selfenergy(basis%nbf,basis%nbf,nspin),emp2
!=====

 integer     :: bbf,ibf,kbf,lbf
 integer     :: astate,bstate
 real(dp)    :: homo,lumo,fact_full,fact_empty

#ifdef CHI0
 logical,parameter    :: ring_only=.true.
#else
 logical,parameter    :: ring_only=.false.
#endif

 real(dp),allocatable :: selfenergy_ring(:,:,:,:)
 real(dp),allocatable :: selfenergy_sox(:,:,:,:)
 real(dp)             :: selfenergy_final(basis%nbf,basis%nbf,nspin)

 integer :: nomegai
 integer :: iomegai
 real(dp),allocatable :: omegai(:)

 complex(dp)           :: ieta
 integer               :: istate,jstate,kstate
 integer               :: abispin,jkspin
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek
 real(dp)              :: eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin)
 real(dp)              :: omega
 real(dp)              :: zz(nspin)
 real(dp)              :: fact_real,fact_nega
 real(dp)              :: emp2_sox,emp2_ring
 logical               :: file_exists
 character(len=3)      :: ctmp
 integer               :: selfenergyfile,iomegafile
!=====

 call start_clock(timing_mp2_self)

 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp

 ieta = im * pole_eta
 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)

 write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
 select case(method)
 case(QS)
   write(stdout,*) 'with the QP self-consistent approach'
 case(perturbative)
   write(stdout,*) 'with the perturbative approach'
 end select


 if(method==QS) then
   nomegai=1
   allocate(omegai(nomegai))
   omegai(1) =  0.00_dp
 else
   ! look for manual omegas
   inquire(file='manual_omega',exist=file_exists)
   if(file_exists) then
     msg='reading frequency file for self-energy evaluation'
     call issue_warning(msg)
     open(newunit=iomegafile,file='manual_omega',status='old')
     read(iomegafile,*) nomegai
     allocate(omegai(nomegai))
     read(iomegafile,*) omegai(1)
     read(iomegafile,*) omegai(nomegai)
     close(iomegafile)
     do iomegai=2,nomegai-1
       omegai(iomegai) = omegai(1) + (omegai(nomegai)-omegai(1)) * (iomegai-1) / DBLE(nomegai-1)
     enddo
   else
     nomegai=3
     allocate(omegai(nomegai))
     omegai(1) = -0.01_dp
     omegai(2) =  0.00_dp
     omegai(3) =  0.01_dp
   endif
 endif

 allocate(selfenergy_ring(nomegai,nstate,nstate,nspin))
 allocate(selfenergy_sox(nomegai,nstate,nstate,nspin))

 selfenergy_ring(:,:,:,:) = 0.0_dp
 selfenergy_sox(:,:,:,:)  = 0.0_dp

#ifdef _OPENMP
 write(stdout,*) 'OPENMP is used for the MP2 self-energy'
#endif
 do abispin=1,nspin
   do istate=1,nstate !LOOP of the first Green's function

     call calculate_eri_4center_eigen(basis%nbf,nstate,c_matrix,istate,abispin,eri_eigenstate_i)

     fi = occupation(istate,abispin)
     ei = energy(istate,abispin)

     do astate=1,nstate ! external loop ( bra )
       do bstate=1,nstate ! external loop ( ket )


         do iomegai=1,nomegai
           omega = energy(bstate,abispin) + omegai(iomegai)
    
    
    
           do jkspin=1,nspin
             do jstate=1,nstate !LOOP of the second Green's function
               fj=occupation(jstate,jkspin)
               ej=energy(jstate,jkspin)
    
               do kstate=1,nstate !LOOP of the third Green's function
                 fk=occupation(kstate,jkspin)
                 ek=energy(kstate,jkspin)
    
                 fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**2
                 fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**2
    
                 if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle
    
                 fact_real = REAL( fact_occ1 / (omega-ei+ej-ek+ieta) + fact_occ2 / (omega-ei+ej-ek-ieta) , dp)
                 fact_nega = REAL( fact_occ1 / (energy(astate,abispin)-ei+ej-ek+ieta) , dp )
    
                 selfenergy_ring(iomegai,astate,bstate,abispin) = selfenergy_ring(iomegai,astate,bstate,abispin) &
                          + fact_real * eri_eigenstate_i(astate,kstate,jstate,jkspin) &
                                         * eri_eigenstate_i(bstate,jstate,kstate,jkspin)

                 if(iomegai==1 .AND. astate==bstate .AND. occupation(astate,abispin)>completely_empty) then
                   emp2_ring = emp2_ring + occupation(astate,abispin) &
                                         * fact_nega * eri_eigenstate_i(astate,kstate,jstate,jkspin) &
                                                     * eri_eigenstate_i(bstate,jstate,kstate,jkspin)
                 endif
    
                 if( abispin == jkspin ) then
                   selfenergy_sox(iomegai,astate,bstate,abispin) = selfenergy_sox(iomegai,astate,bstate,abispin) &
                            - fact_real * eri_eigenstate_i(astate,jstate,kstate,abispin) &
                                           * eri_eigenstate_i(jstate,kstate,bstate,abispin) / spin_fact

                   if(iomegai==1 .AND. astate==bstate .AND. occupation(astate,abispin)>completely_empty) then
                     emp2_sox = emp2_sox - occupation(astate,abispin) &
                               * fact_nega * eri_eigenstate_i(astate,jstate,kstate,jkspin) &
                                           * eri_eigenstate_i(jstate,kstate,bstate,abispin) / spin_fact
                   endif
                 endif
    
    
               enddo
    
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

   if(file_exists .AND. is_iomaster ) then
     do astate=1,MIN(2,nstate)
       write(ctmp,'(i3.3)') astate
       open(newunit=selfenergyfile,file='selfenergy_omega_state'//TRIM(ctmp))
       do iomegai=1,nomegai
         write(selfenergyfile,'(20(f12.6,2x))') DBLE(omegai(iomegai))+energy(astate,:),&
&                         DBLE(selfenergy_ring(iomegai,astate,astate,:)),&
&                         DBLE(selfenergy_sox(iomegai,astate,astate,:)),&
&                         DBLE(selfenergy_ring(iomegai,astate,astate,:))+&
&                         DBLE(selfenergy_sox(iomegai,astate,astate,:)),&
&                         DBLE(omegai(iomegai))-exchange_m_vxc_diag(astate,:)
       enddo
       write(selfenergyfile,*)
     enddo
     close(selfenergyfile)
     call die('output the self energy in a file')
   endif


   write(stdout,*) '=============================='
   write(stdout,*) ' selfenergy RING + SOX'
   write(stdout,'(a)') ' #         Energies           Sigx-Vxc       one-ring              SOX             MP2              Z            QP-eigenvalues'
   do bstate=1,nstate 
     zz(:) = REAL( selfenergy_ring(3,bstate,bstate,:)+selfenergy_sox(3,bstate,bstate,:) &
                 - selfenergy_ring(1,bstate,bstate,:)-selfenergy_sox(1,bstate,bstate,:) ) / REAL( omegai(3)-omegai(1) )
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )
     write(stdout,'(i4,18(3x,f14.5))') bstate,energy(bstate,:)*Ha_eV,&
           exchange_m_vxc_diag(bstate,:)*Ha_eV,&
           selfenergy_ring(2,bstate,bstate,:) * Ha_eV,&
           selfenergy_sox (2,bstate,bstate,:) * Ha_eV,&
         ( selfenergy_ring(2,bstate,bstate,:)+selfenergy_sox(2,bstate,bstate,:) ) * Ha_eV,&
         zz(:),&
         ( energy(bstate,:) + exchange_m_vxc_diag(bstate,:) + selfenergy_ring(2,bstate,bstate,:) + selfenergy_sox(2,bstate,bstate,:) ) * Ha_eV
   enddo
   write(stdout,*) '=============================='

   selfenergy(:,:,:) = 0.0_dp
   if(.NOT.ring_only) then
     do bstate=1,nstate
       selfenergy(bstate,bstate,:) = selfenergy_ring(2,bstate,bstate,:) + selfenergy_sox(2,bstate,bstate,:)
     enddo
   else
     do bstate=1,nstate
       selfenergy(bstate,bstate,:) = selfenergy_ring(2,bstate,bstate,:)
     enddo
   endif

 endif


 if(method==QS) then
 ! QS trick of Faleev-Kotani-vanSchilfgaarde
   do abispin=1,nspin
     if(.NOT.ring_only) then
       selfenergy_final(:,:,abispin) = 0.5_dp * ( selfenergy_ring(1,:,:,abispin) + selfenergy_sox(1,:,:,abispin) &
                                               +  TRANSPOSE( selfenergy_ring(1,:,:,abispin) + selfenergy_sox(1,:,:,abispin) ) )
     else
       write(stdout,*) 'ring_only'
       selfenergy_final(:,:,abispin) = 0.5_dp * ( selfenergy_ring(1,:,:,abispin)  &
                                               +  TRANSPOSE( selfenergy_ring(1,:,:,abispin) ) )
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

 deallocate(omegai)
 deallocate(selfenergy_ring)
 deallocate(selfenergy_sox)

 call stop_clock(timing_mp2_self)

end subroutine mp2_selfenergy
!=========================================================================
