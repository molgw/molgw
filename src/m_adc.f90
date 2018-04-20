!==================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
!
!==================================================================
module m_adc
 use m_definitions
 use m_tools
 use m_basis_set
 use m_inputparam
 use m_eri_ao_mo
 use m_selfenergy_tools



contains

!==================================================================
subroutine adc2(basis,nstate,occupation,energy,c_matrix)
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)  :: nstate
 real(dp),intent(in) :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in) :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer,parameter :: ncycle=1
 logical,parameter :: TWOPH_TDA=.FALSE.
 integer :: icycle
 integer :: ispin,pspin
 integer :: istate,jstate,kstate,lstate ! occupied
 integer :: astate,bstate,cstate,dstate ! empty
 integer :: pstate,qstate,rstate,sstate ! occupied OR empty
 integer :: nhomo,nlumo,nhole,npart
 integer :: bmat,imat,jmat
 real(dp),allocatable :: b_matrix(:,:)
 real(dp),allocatable :: x_matrix(:,:)
 real(dp),allocatable :: eigenvalue(:)
 real(dp),allocatable :: g_pq(:,:),siginf_pq(:,:)
!=====

 !
 ! Routine closely follows the notations of Schirmer, Cederbaum, Walter, Phys. Rev. A 28, 1242 (1983).
 !
 write(stdout,'(/,1x,a)') 'Entering ADC(2) routine'

 call selfenergy_set_state_range(nstate,occupation)

 if( nspin /= 1 ) call die('adc2: not coded for spin-unrestricted calculations')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix)
   call die('adc2: not coded with RI')
 else
   call calculate_eri_4center_eigen_uks(c_matrix,1,nvirtual_G-1)  ! TODO set the nstate_min to a more finely tuned value
 endif

 ispin = 1

 ! Find highest occupied state
 nhomo = 0
 do istate=1,nstate
   if( occupation(istate,1) < completely_empty)  cycle
   nhomo = istate
 enddo
 nlumo = nhomo + 1

 nhole = nhomo - ncore_G
 npart = nvirtual_G - 1 - nhomo

 write(stdout,'(1x,a,i4,5x,a,i4,5x,a,i4)') 'Total states: ',(nhole+npart)*nspin,'Hole states: ',nhole*nspin,'Particle states: ',npart*nspin

 bmat = ( nhole + npart + nhole * (npart * (npart+1)) / 2 + (nhole * (nhole+1)) / 2 * npart ) * nspin
 write(stdout,*) 'ADC B matrix',bmat,' x ',bmat

 call clean_allocate('B matrix',b_matrix,bmat,bmat)
 call clean_allocate('Eigenvectors',x_matrix,bmat,bmat)
 allocate(eigenvalue(bmat))
 allocate(g_pq(ncore_G+1:nvirtual_G-1,ncore_G+1:nvirtual_G-1))
 allocate(siginf_pq(ncore_G+1:nvirtual_G-1,ncore_G+1:nvirtual_G-1))


 siginf_pq(:,:) = 0.0_dp

 do icycle=1,ncycle
   write(stdout,'(/,/,1x,a,i4)') 'Cycle: ',icycle
   !
   ! A) \epsilon part
   ! 
   b_matrix(:,:) = 0.0_dp
   b_matrix(ncore_G+1:nvirtual_G-1,ncore_G+1:nvirtual_G-1) = siginf_pq(:,:)
   jmat = 0
   do pstate=ncore_G+1,nvirtual_G-1
     do pspin=1,nspin
       jmat = jmat + 1
       b_matrix(jmat,jmat) = energy(pstate,pspin)
     enddo
   enddo
  
   !
   ! B) U part
   ! 
   do pstate=ncore_G+1,nvirtual_G-1
     imat = pstate - ncore_G
     jmat = nhole + npart
     !
     ! 2 particles - 1 hole
     do istate=ncore_G+1,nhomo
       do astate=nlumo,nvirtual_G-1
         do bstate=astate,nvirtual_G-1
           jmat = jmat + 1
           b_matrix(imat,jmat) = eri_physics_as(pstate,istate,astate,bstate)
         enddo
       enddo
     enddo
     !
     ! 2 holes - 1 particle
     do istate=ncore_G+1,nhomo
       do jstate=istate,nhomo
         do astate=nlumo,nvirtual_G-1
           jmat = jmat + 1
           b_matrix(imat,jmat) = eri_physics_as(pstate,astate,istate,jstate)
         enddo
       enddo
     enddo
  
   enddo
  
   !
   ! C) K part
   ! 
   !
   ! Part I: 2 particles - 1 hole
   jmat = nhole + npart
   do istate=ncore_G+1,nhomo
     do astate=nlumo,nvirtual_G-1
       do bstate=astate,nvirtual_G-1
         jmat = jmat + 1
         b_matrix(jmat,jmat) = energy(astate,ispin) + energy(bstate,ispin) - energy(istate,ispin)
  !       write(*,*) jmat,b_matrix(jmat,jmat)*Ha_eV
       enddo
     enddo
   enddo
   !
   ! Part II: 2 holes - 1 particle
   do istate=ncore_G+1,nhomo
     do jstate=istate,nhomo
       do astate=nlumo,nvirtual_G-1
         jmat = jmat + 1
         b_matrix(jmat,jmat) = energy(istate,ispin) + energy(jstate,ispin) - energy(astate,ispin)
  !       write(*,*) jmat,b_matrix(jmat,jmat)*Ha_eV
       enddo
     enddo
   enddo
  
  
   !
   ! D) C part
   ! 
   if( TWOPH_TDA ) then
     ! Part I: 2ph-TDA
     imat = nhole + npart
     do istate=ncore_G+1,nhomo
       do astate=nlumo,nvirtual_G-1
         do bstate=astate,nvirtual_G-1
           imat = imat + 1
           jmat = nhole + npart
    
           do jstate=ncore_G+1,nhomo
             do cstate=nlumo,nvirtual_G-1
               do dstate=cstate,nvirtual_G-1
                 jmat = jmat + 1
                 b_matrix(imat,jmat) = b_matrix(imat,jmat) + eri_physics_as(astate,bstate,cstate,dstate) * delta(istate,jstate) &
                                                           - eri_physics_as(jstate,bstate,istate,dstate) * delta(astate,cstate) &
                                                           - eri_physics_as(jstate,astate,istate,cstate) * delta(bstate,dstate) &
                                                           + eri_physics_as(bstate,astate,cstate,dstate) * delta(istate,jstate) &
                                                           - eri_physics_as(jstate,astate,istate,dstate) * delta(bstate,cstate) &
                                                           - eri_physics_as(jstate,bstate,istate,cstate) * delta(astate,dstate)
               enddo
             enddo
           enddo
    
         enddo
       enddo
     enddo
    
     ! Part II: 2ph-TDA
     imat = nhole + npart + nhole * (npart * (npart+1)) / 2
     do istate=ncore_G+1,nhomo
       do jstate=istate,nhomo
         do astate=nlumo,nvirtual_G-1
           imat = imat + 1
           jmat = nhole + npart + nhole * (npart * (npart+1)) / 2
    
           do kstate=ncore_G+1,nhomo
             do lstate=kstate,nhomo
               do bstate=nlumo,nvirtual_G-1
                 jmat = jmat + 1
                 b_matrix(imat,jmat) = b_matrix(imat,jmat) - eri_physics_as(istate,jstate,kstate,lstate) * delta(astate,bstate) &
                                                           + eri_physics_as(bstate,jstate,astate,lstate) * delta(istate,kstate) &
                                                           + eri_physics_as(bstate,istate,astate,kstate) * delta(jstate,lstate) &
                                                           - eri_physics_as(jstate,istate,kstate,lstate) * delta(astate,bstate) &
                                                           + eri_physics_as(bstate,istate,astate,lstate) * delta(jstate,kstate) &
                                                           + eri_physics_as(bstate,jstate,astate,kstate) * delta(istate,lstate)
               enddo
             enddo
           enddo
    
         enddo
       enddo
     enddo
  
   endif
  
  
  
   ! Symmetrize the matrix B
   do imat=1,bmat
     do jmat=imat+1,bmat
       b_matrix(jmat,imat) = b_matrix(imat,jmat)
     enddo
   enddo
  
  
  
   call diagonalize(b_matrix,eigenvalue,x_matrix)
  
   write(stdout,'(/,1x,a)') '================================='
   write(stdout,'(1x,a)') 'G Poles (eV)'
   do imat=1,bmat
     if( NORM2(x_matrix(1:nhole+npart,imat)) > 0.10_dp ) & 
       write(stdout,'(2(2x,f12.6))') eigenvalue(imat) * Ha_eV,NORM2(x_matrix(1:nhole+npart,imat))
   enddo
   write(stdout,'(1x,a,/)') '================================='
  
  
   siginf_pq(:,:) = 0.0_dp
   do istate=ncore_G+1,nhomo
     do pstate=ncore_G+1,nvirtual_G-1
       do qstate=ncore_G+1,nvirtual_G-1
         siginf_pq(pstate,qstate) = siginf_pq(pstate,qstate) - eri_physics_as(pstate,istate,qstate,istate)
       enddo
     enddo
   enddo
   g_pq(:,:) = 0.0_dp
   do imat=1,bmat
     if( eigenvalue(imat) > 0.0_dp ) cycle
     forall(pstate=ncore_G+1:nvirtual_G-1,qstate=ncore_G+1:nvirtual_G-1)
       g_pq(pstate,qstate) = g_pq(pstate,qstate) + x_matrix(pstate-ncore_G,imat) * x_matrix(qstate-ncore_G,imat)
     end forall
   enddo

   do rstate=ncore_G+1,nvirtual_G-1
   do sstate=ncore_G+1,nvirtual_G-1
     do pstate=ncore_G+1,nvirtual_G-1
       do qstate=ncore_G+1,nvirtual_G-1
         siginf_pq(pstate,qstate) = siginf_pq(pstate,qstate) + eri_physics_as(pstate,rstate,qstate,sstate) * g_pq(rstate,sstate)
       enddo
     enddo
   enddo
   enddo


 enddo ! big cycle loop

 call clean_deallocate('B matrix',b_matrix)
 call clean_deallocate('Eigenvectors',x_matrix)
 deallocate(g_pq)
 deallocate(eigenvalue)

end subroutine adc2


!==================================================================
pure function eri_physics_as(istate,jstate,kstate,lstate)
 implicit none
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp) :: eri_physics_as
!=====
!=====

 eri_physics_as = 2.0_dp * eri_4center_eigen_uks(istate,kstate,jstate,lstate)  &
                  - eri_4center_eigen_uks(istate,lstate,jstate,kstate)

end function eri_physics_as


!==================================================================
pure function delta(istate,jstate)
 implicit none
 integer,intent(in) :: istate,jstate
 real(dp) :: delta
!=====
!=====

 delta = 0.0_dp
 if( istate == jstate ) delta = 1.0_dp

end function delta

!==================================================================
end module m_adc
!==================================================================
