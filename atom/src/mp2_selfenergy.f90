!=========================================================================
#include "macros.h"
!=========================================================================
subroutine mp2_selfenergy(method,nspin,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,selfenergy,emp2)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_warning
 use m_basis_set
 use m_eri
 implicit none

 integer,intent(in)  :: method,nspin
 type(basis_set)     :: basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin),energy(basis%nbf,nspin),exchange_m_vxc_diag(basis%nbf,nspin)
 real(dp),intent(in) :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in) :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(out) :: selfenergy(basis%nbf,basis%nbf,nspin),emp2
!=====

 integer     :: bbf,ibf,kbf,lbf
 integer     :: aorbital,borbital
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

 complex(dp),parameter :: ieta=(0.0_dp,0.0001_dp)
 integer               :: iorbital,jorbital,korbital
 integer               :: abispin,jkspin
 real(dp)              :: spin_fact
 real(dp)              :: fact_occ1,fact_occ2
 real(dp)              :: fi,fj,fk,ei,ej,ek
#if defined LOW_MEMORY2 || defined LOW_MEMORY3
 real(dp)              :: eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin)
#else
 real(dp)              :: eri_eigenstate(basis%nbf,basis%nbf,basis%nbf,basis%nbf,nspin,nspin)
#endif
 real(dp)              :: omega
 real(dp)              :: zz(nspin)
 real(dp)              :: fact_real,fact_nega
 real(dp)              :: emp2_sox,emp2_ring
 logical               :: file_exists
 character(len=3)      :: ctmp
!=====

 spin_fact = REAL(-nspin+3,dp)
 emp2_ring = 0.0_dp
 emp2_sox  = 0.0_dp

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'Perform the second-order self-energy calculation'
 select case(method)
 case(QS)
   WRITE_MASTER(*,*) 'with the QP self-consistent approach'
 case(perturbative)
   WRITE_MASTER(*,*) 'with the perturbative approach'
 end select


#ifndef LOW_MEMORY2
#ifndef LOW_MEMORY3
 call transform_eri_basis_fast(basis%nbf,nspin,c_matrix,eri_eigenstate)
#endif
#endif


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
     open(12,file='manual_omega',status='old')
     read(12,*) nomegai
     allocate(omegai(nomegai))
     read(12,*) omegai(1)
     read(12,*) omegai(nomegai)
     close(12)
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

 allocate(selfenergy_ring(nomegai,basis%nbf,basis%nbf,nspin))
 allocate(selfenergy_sox(nomegai,basis%nbf,basis%nbf,nspin))

 selfenergy_ring(:,:,:,:) = 0.0_dp
 selfenergy_sox(:,:,:,:)  = 0.0_dp

#ifdef _OPENMP
 WRITE_MASTER(*,*) 'OPENMP is used for the MP2 self-energy'
#endif
 do abispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED) &
#if defined LOW_MEMORY2 || defined LOW_MEMORY3
!$OMP PRIVATE(omega,fi,ei,fj,ej,fk,ek,fact_occ1,fact_occ2,fact_real,fact_nega,eri_eigenstate_i) 
#else
!$OMP PRIVATE(omega,fi,ei,fj,ej,fk,ek,fact_occ1,fact_occ2,fact_real,fact_nega) 
#endif
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:emp2_ring,emp2_sox,selfenergy_ring,selfenergy_sox)
   do iorbital=1,basis%nbf !LOOP of the first Green's function
#if defined LOW_MEMORY2 || defined LOW_MEMORY3
     call transform_eri_basis_lowmem(nspin,c_matrix,iorbital,abispin,eri_eigenstate_i)
#endif
     do aorbital=1,basis%nbf ! external loop ( bra )


       do borbital=1,basis%nbf ! external loop ( ket )

         do iomegai=1,nomegai
           omega = energy(borbital,abispin) + omegai(iomegai)
    
    
           fi=occupation(iorbital,abispin)
           ei=energy(iorbital,abispin)
    
           do jkspin=1,nspin
             do jorbital=1,basis%nbf !LOOP of the second Green's function
               fj=occupation(jorbital,jkspin)
               ej=energy(jorbital,jkspin)
    
               do korbital=1,basis%nbf !LOOP of the third Green's function
                 fk=occupation(korbital,jkspin)
                 ek=energy(korbital,jkspin)
    
                 fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**2
                 fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**2
    
                 if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle
    
                 fact_real = REAL( fact_occ1 / (omega-ei+ej-ek+ieta) + fact_occ2 / (omega-ei+ej-ek-ieta) , dp)
                 fact_nega = REAL( fact_occ1 / (energy(aorbital,abispin)-ei+ej-ek+ieta) , dp )
    
#if defined LOW_MEMORY2 || defined LOW_MEMORY3
                 selfenergy_ring(iomegai,aorbital,borbital,abispin) = selfenergy_ring(iomegai,aorbital,borbital,abispin) &
                          + fact_real * eri_eigenstate_i(aorbital,korbital,jorbital,jkspin) &
                                         * eri_eigenstate_i(borbital,jorbital,korbital,jkspin)

                 if(iomegai==1 .AND. aorbital==borbital .AND. occupation(aorbital,abispin)>completely_empty) then
                   emp2_ring = emp2_ring + occupation(aorbital,abispin) &
                                         * fact_nega * eri_eigenstate_i(aorbital,korbital,jorbital,jkspin) &
                                                     * eri_eigenstate_i(borbital,jorbital,korbital,jkspin)
                 endif
    
                 if( abispin == jkspin ) then
                   selfenergy_sox(iomegai,aorbital,borbital,abispin) = selfenergy_sox(iomegai,aorbital,borbital,abispin) &
                            - fact_real * eri_eigenstate_i(aorbital,jorbital,korbital,abispin) &
                                           * eri_eigenstate_i(jorbital,korbital,borbital,abispin) / spin_fact

                   if(iomegai==1 .AND. aorbital==borbital .AND. occupation(aorbital,abispin)>completely_empty) then
                     emp2_sox = emp2_sox - occupation(aorbital,abispin) &
                               * fact_nega * eri_eigenstate_i(aorbital,jorbital,korbital,jkspin) &
                                           * eri_eigenstate_i(jorbital,korbital,borbital,abispin) / spin_fact
                   endif
                 endif
#else
                 selfenergy_ring(iomegai,aorbital,borbital,abispin) = selfenergy_ring(iomegai,aorbital,borbital,abispin) &
                          + fact_real * eri_eigenstate(aorbital,iorbital,korbital,jorbital,abispin,jkspin) &
                                         * eri_eigenstate(jorbital,korbital,iorbital,borbital,jkspin,abispin)

                 if(iomegai==1 .AND. aorbital==borbital .AND. occupation(aorbital,abispin)>completely_empty) then
                   emp2_ring = emp2_ring + occupation(aorbital,abispin) &
                                         * fact_nega * eri_eigenstate(aorbital,iorbital,korbital,jorbital,abispin,jkspin) &
                                                     * eri_eigenstate(jorbital,korbital,iorbital,borbital,jkspin,abispin)
                 endif
    
                 if( abispin == jkspin ) then
                   selfenergy_sox(iomegai,aorbital,borbital,abispin) = selfenergy_sox(iomegai,aorbital,borbital,abispin) &
                            - fact_real * eri_eigenstate(aorbital,iorbital,jorbital,korbital,abispin,abispin) &
                                           * eri_eigenstate(iorbital,jorbital,korbital,borbital,abispin,abispin) / spin_fact

                   if(iomegai==1 .AND. aorbital==borbital .AND. occupation(aorbital,abispin)>completely_empty) then
                     emp2_sox = emp2_sox - occupation(aorbital,abispin) &
                               * fact_nega * eri_eigenstate(aorbital,iorbital,jorbital,korbital,abispin,jkspin) &
                                           * eri_eigenstate(iorbital,jorbital,korbital,borbital,abispin,abispin) / spin_fact
                   endif
                 endif
#endif

    
    
               enddo
    
             enddo
           enddo
         enddo
       enddo 
     enddo
   enddo 
!$OMP END DO
!$OMP END PARALLEL
 enddo ! abispin

 emp2_ring = 0.5_dp * emp2_ring
 emp2_sox  = 0.5_dp * emp2_sox
 emp2 = emp2_ring + emp2_sox
 WRITE_MASTER(*,'(/,a)')       ' MP2 Energy'
 WRITE_MASTER(*,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
 WRITE_MASTER(*,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
 WRITE_MASTER(*,'(a,f14.8,/)') ' MP2 correlation :',emp2

 if(method==perturbative) then

   if(file_exists) then
     do aorbital=1,MIN(2,basis%nbf)
       write(ctmp,'(i3.3)') aorbital
       open(20+aorbital,file='selfenergy_omega_state'//TRIM(ctmp))
       do iomegai=1,nomegai
         WRITE_MASTER(20+aorbital,'(20(f12.6,2x))') DBLE(omegai(iomegai))+energy(aorbital,:),&
&                         DBLE(selfenergy_ring(iomegai,aorbital,aorbital,:)),&
&                         DBLE(selfenergy_sox(iomegai,aorbital,aorbital,:)),&
&                         DBLE(selfenergy_ring(iomegai,aorbital,aorbital,:))+&
&                         DBLE(selfenergy_sox(iomegai,aorbital,aorbital,:)),&
&                         DBLE(omegai(iomegai))-exchange_m_vxc_diag(aorbital,:)
       enddo
       WRITE_MASTER(20+aorbital,*)
     enddo
     close(20+aorbital)
     stop'output the self energy in a file'
   endif


   WRITE_MASTER(*,*) '=============================='
   WRITE_MASTER(*,*) ' selfenergy RING + SOX'
   WRITE_MASTER(*,'(a)') ' #         Energies           Sigx-Vxc       one-ring              SOX             MP2              Z            QP-eigenvalues'
   do borbital=1,basis%nbf 
     zz(:) = REAL( selfenergy_ring(3,borbital,borbital,:)+selfenergy_sox(3,borbital,borbital,:) &
                 - selfenergy_ring(1,borbital,borbital,:)-selfenergy_sox(1,borbital,borbital,:) ) / REAL( omegai(3)-omegai(1) )
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )
     WRITE_MASTER(*,'(i4,18(3x,f14.5))') borbital,energy(borbital,:)*Ha_eV,&
           exchange_m_vxc_diag(borbital,:)*Ha_eV,&
           selfenergy_ring(2,borbital,borbital,:) * Ha_eV,&
           selfenergy_sox (2,borbital,borbital,:) * Ha_eV,&
         ( selfenergy_ring(2,borbital,borbital,:)+selfenergy_sox(2,borbital,borbital,:) ) * Ha_eV,&
         zz(:),&
         ( energy(borbital,:) + exchange_m_vxc_diag(borbital,:) + selfenergy_ring(2,borbital,borbital,:) + selfenergy_sox(2,borbital,borbital,:) ) * Ha_eV
   enddo
   WRITE_MASTER(*,*) '=============================='

   selfenergy(:,:,:) = 0.0_dp
   if(.NOT.ring_only) then
     do borbital=1,basis%nbf
       selfenergy(borbital,borbital,:) = selfenergy_ring(2,borbital,borbital,:) + selfenergy_sox(2,borbital,borbital,:)
     enddo
   else
     do borbital=1,basis%nbf
       selfenergy(borbital,borbital,:) = selfenergy_ring(2,borbital,borbital,:)
     enddo
   endif

 endif


 if(method==QS) then
 ! QS trick of Faleev-Kotani-vanSchilfgaarde
   do abispin=1,nspin
     if(.NOT.ring_only) then
       selfenergy_final(:,:,abispin) = 0.5_dp * ( selfenergy_ring(1,:,:,abispin) + selfenergy_sox(1,:,:,abispin) &
                                               +  transpose( selfenergy_ring(1,:,:,abispin) + selfenergy_sox(1,:,:,abispin) ) )
     else
       WRITE_MASTER(*,*) 'ring_only'
       selfenergy_final(:,:,abispin) = 0.5_dp * ( selfenergy_ring(1,:,:,abispin)  &
                                               +  transpose( selfenergy_ring(1,:,:,abispin) ) )
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


end subroutine mp2_selfenergy
!=========================================================================
