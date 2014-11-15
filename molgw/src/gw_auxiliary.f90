!=========================================================================
#include "macros.h"
!=========================================================================

#ifdef AUXIL_BASIS
!=========================================================================
subroutine polarizability_rpa_aux(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,rpa_correlation,wpol)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_spectral_function
 implicit none

 integer,intent(in)  :: nspin
 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin)
 real(dp),intent(in) :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in) :: sinv_v_sinv(prod_basis%nbf,prod_basis%nbf)
 real(dp),intent(out) :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol
!=====
 integer :: pbf,qbf,ibf,jbf,ijspin,klspin
 integer :: iorbital,jorbital,korbital,lorbital
 integer :: ipole
 integer :: t_ij,t_kl

 real(dp) :: spin_fact 
 real(dp) :: overlap_tmp
 real(dp) :: h_2p(wpol%npole,wpol%npole)
 real(dp) :: eigenvalue(wpol%npole),eigenvector(wpol%npole,wpol%npole),eigenvector_inv(wpol%npole,wpol%npole)
 real(dp) :: matrix(wpol%npole,wpol%npole)
 real(dp) :: three_bf_overlap(prod_basis%nbf,basis%nbf,basis%nbf,nspin)
 real(dp) :: coulomb_left(prod_basis%nbf)
 real(dp) :: exchange_left(prod_basis%nbf)
 real(dp) :: left(wpol%npole,prod_basis%nbf),right(wpol%npole,prod_basis%nbf)

 logical :: TDHF=.FALSE.
!=====
 spin_fact = REAL(-nspin+3,dp)

 WRITE_MASTER(*,'(/,a)') ' calculating CHI alla rpa'
 if(TDHF) then
   msg='calculating the TDHF polarizability'
   call issue_warning(msg)
 endif

 WRITE_MASTER(*,*) 'three wf overlaps'
 call start_clock(timing_overlap3)
 !
 ! set up overlaps (aux function * orbital * orbital)
 three_bf_overlap(:,:,:,:) = 0.0_dp
   do pbf=1,prod_basis%nbf
     do ibf=1,basis%nbf
       do jbf=1,basis%nbf
         call overlap_three_basis_function(prod_basis%bf(pbf),basis%bf(ibf),basis%bf(jbf),overlap_tmp)

 do ijspin=1,nspin


         do iorbital=1,basis%nbf
           do jorbital=1,basis%nbf
             three_bf_overlap(pbf,iorbital,jorbital,ijspin) = three_bf_overlap(pbf,iorbital,jorbital,ijspin) &
&                + c_matrix(ibf,iorbital,ijspin) * c_matrix(jbf,jorbital,ijspin) &
&                    * overlap_tmp 
           enddo
         enddo

 enddo

       enddo
     enddo
   enddo
 call stop_clock(timing_overlap3)
 WRITE_MASTER(*,*) 'three wf overlaps: DONE'

 t_ij=0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
       if(iorbital==jorbital) cycle  ! intra state transitions are not allowed!
       if( abs(occupation(jorbital,ijspin)-occupation(iorbital,ijspin))<completely_empty ) cycle
       t_ij=t_ij+1

       coulomb_left = matmul( three_bf_overlap(:,iorbital,jorbital,ijspin), sinv_v_sinv )


       t_kl=0
       do klspin=1,nspin
         do korbital=1,basis%nbf 
if(TDHF)          exchange_left = matmul( three_bf_overlap(:,iorbital,korbital,ijspin), sinv_v_sinv )
           do lorbital=1,basis%nbf 
             if(korbital==lorbital) cycle  ! intra state transitions are not allowed!
             if( abs(occupation(lorbital,klspin)-occupation(korbital,klspin))<completely_empty ) cycle
             t_kl=t_kl+1

             h_2p(t_ij,t_kl) = dot_product(coulomb_left, three_bf_overlap(:,korbital,lorbital,klspin)) &
&                       * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )

if(TDHF) then
        if(ijspin==klspin) then
             h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -dot_product(exchange_left, three_bf_overlap(:,jorbital,lorbital,ijspin)) &
&                     * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact 
        endif
endif

           enddo
         enddo
       enddo !klspin

       h_2p(t_ij,t_ij) =  h_2p(t_ij,t_ij) + ( energy(jorbital,ijspin) - energy(iorbital,ijspin) )

!       WRITE_MASTER(*,'(4(i4,2x),2(2x,f12.6))') t_ij,iorbital,jorbital,ijspin,( energy(jorbital,ijspin) - energy(iorbital,ijspin) ),h_2p(t_ij,t_ij)
     enddo !jorbital
   enddo !iorbital
 enddo ! ijspin

 WRITE_MASTER(*,*) 'diago 2-particle hamiltonian'
 WRITE_MASTER(*,*) 'matrix',wpol%npole,'x',wpol%npole
 call start_clock(timing_diago_h2p)
 call diagonalize_general(wpol%npole,h_2p,eigenvalue,eigenvector)
 call stop_clock(timing_diago_h2p)
 WRITE_MASTER(*,*) 'diago finished'
   
 call start_clock(timing_inversion_s2p)
 call invert(wpol%npole,eigenvector,eigenvector_inv)
 call stop_clock(timing_inversion_s2p)

 left(:,:) = 0.0_dp 
 right(:,:) = 0.0_dp 
 do pbf=1,prod_basis%nbf
   t_ij=0
   do ijspin=1,nspin
     do iorbital=1,basis%nbf 
       do jorbital=1,basis%nbf
         if(iorbital==jorbital) cycle  ! intra state transitions are not allowed!
         if( abs(occupation(jorbital,ijspin)-occupation(iorbital,ijspin))<completely_empty ) cycle
         t_ij=t_ij+1
         left(:,pbf) = left(:,pbf) + three_bf_overlap(pbf,iorbital,jorbital,ijspin) * eigenvector(t_ij,:)
         right(:,pbf) = right(:,pbf) + three_bf_overlap(pbf,iorbital,jorbital,ijspin) * eigenvector_inv(:,t_ij) &
                                              * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )
       enddo
     enddo
   enddo
 enddo

 wpol%residu_left  = transpose( matmul( sinv_v_sinv , transpose(left) ) ) 
 wpol%residu_right = matmul( right , sinv_v_sinv )  
 wpol%pole = eigenvalue


 rpa_correlation=0.0_dp

end subroutine polarizability_rpa_aux


!=========================================================================
subroutine gw_selfenergy_aux(method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_basis_set
 use m_spectral_function
 implicit none

 integer,intent(in)  :: method,nspin
 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin),energy(basis%nbf,nspin),exchange_m_vxc_diag(basis%nbf,nspin)
 real(dp),intent(in) :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in) :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out) :: selfenergy(basis%nbf,basis%nbf,nspin)

 integer :: nomegai
 integer :: iomegai
 complex(dp),allocatable :: omegai(:)
 complex(dp),allocatable :: selfenergy_tmp(:,:,:,:)

 logical     :: file_exists=.FALSE.

 integer     :: bbf,ibf,kbf,lbf
 integer     :: aorbital,borbital
 integer     :: iorbital,ispin,ipole
 real(dp)    :: spin_fact,overlap_tmp
 real(dp)    :: three_state_overlap_i(prod_basis%nbf,basis%nbf)
 real(dp)    :: bra(wpol%npole,basis%nbf),ket(wpol%npole,basis%nbf)
 real(dp)    :: fact_full,fact_empty
 real(dp)    :: zz(nspin)
 complex(dp) :: energy_complex
!=====
 spin_fact = REAL(-nspin+3,dp)

 WRITE_MASTER(*,*)
 select case(method)
 case(QS)
   WRITE_MASTER(*,*) 'perform a QP self-consistent GW calculation'
 case(perturbative)
   WRITE_MASTER(*,*) 'perform a one-shot G0W0 calculation'
 end select


 if(method==QS) then
   nomegai=1
 else
   ! look for manual omegas
   inquire(file='manual_omega',exist=file_exists)
   WRITE_MASTER(*,*) file_exists
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

 allocate(selfenergy_tmp(nomegai,basis%nbf,basis%nbf,nspin))
 

 selfenergy_tmp(:,:,:,:) = (0.0_dp,0.0_dp)

 do ispin=1,nspin
   do iorbital=1,basis%nbf !INNER LOOP of G

     three_state_overlap_i(:,:) = 0.0_dp
     do bbf=1,basis%nbf
       do kbf=1,prod_basis%nbf
         do ibf=1,basis%nbf

           call overlap_three_basis_function(prod_basis%bf(kbf),basis%bf(bbf),basis%bf(ibf),overlap_tmp)

           do borbital=1,basis%nbf
             three_state_overlap_i(kbf,borbital) = three_state_overlap_i(kbf,borbital) &
&                  + c_matrix(ibf,iorbital,ispin) * c_matrix(bbf,borbital,ispin)  &
&                      * overlap_tmp 
           enddo


         enddo
       enddo 
     enddo 

     bra = matmul(wpol%residu_left,three_state_overlap_i)
     ket = matmul(wpol%residu_right,three_state_overlap_i)

     do ipole=1,wpol%npole
       if( wpol%pole(ipole) < 0.0_dp ) then
         fact_empty = (spin_fact - occupation(iorbital,ispin)) / spin_fact
         fact_full = 0.0_dp
       else
         fact_empty = 0.0_dp
         fact_full = occupation(iorbital,ispin) / spin_fact
       endif
       if( ABS(fact_empty - fact_full) < 0.0001 ) cycle

       if(method==QS) then

         do borbital=1,basis%nbf
           !
           ! Take care about the energies that may lie in the vicinity of the
           ! poles of Sigma
           if( ABS( energy(borbital,ispin) - energy(iorbital,ispin) + wpol%pole(ipole) ) < aimag(ieta) ) then
             energy_complex = energy(borbital,ispin)  + (0.0_dp,1.0_dp) &
&                * ( aimag(ieta) - ABS( energy(borbital,ispin) - energy(iorbital,ispin) + wpol%pole(ipole) ) )
           else
             energy_complex = energy(borbital,ispin) 
           endif
!           WRITE_MASTER(*,'(i4,x,i4,x,i4,x,10(2x,f12.6))') borbital,iorbital,ipole,&
!&                energy(borbital,ispin),energy(iorbital,ispin),wpol%pole(ipole),aimag(energy_complex)
           do aorbital=1,basis%nbf

             selfenergy_tmp(1,aorbital,borbital,ispin) = selfenergy_tmp(1,aorbital,borbital,ispin) &
&                       - bra(ipole,aorbital) * ket(ipole,borbital) &
&                         * ( fact_empty / (       energy_complex  - energy(iorbital,ispin) + wpol%pole(ipole) ) &
&                            -fact_full  / ( CONJG(energy_complex) - energy(iorbital,ispin) + wpol%pole(ipole) ) )

           enddo
         enddo

       else if(method==perturbative) then

         do aorbital=1,basis%nbf
           !
           ! calculate only the diagonal !
!           do borbital=1,basis%nbf
           borbital=aorbital
             do iomegai=1,nomegai
               !
               ! Take care about the energies that may lie in the vicinity of the
               ! poles of Sigma
               if( ABS( energy(borbital,ispin) + omegai(iomegai) - energy(iorbital,ispin) + wpol%pole(ipole) ) < aimag(ieta) ) then
!                 WRITE_MASTER(*,*) 'avoid pole'
                 energy_complex = energy(borbital,ispin)  + (0.0_dp,1.0_dp) &
&                    * ( aimag(ieta) - ABS( energy(borbital,ispin) + omegai(iomegai) - energy(iorbital,ispin) + wpol%pole(ipole) ) )
               else
                 energy_complex = energy(borbital,ispin) 
               endif
               selfenergy_tmp(iomegai,aorbital,borbital,ispin) = selfenergy_tmp(iomegai,aorbital,borbital,ispin) &
&                       - bra(ipole,aorbital) * ket(ipole,borbital) &
&                         * ( fact_empty / (       energy_complex         + omegai(iomegai) &
&                                          - energy(iorbital,ispin) + wpol%pole(ipole)     ) &
&                            -fact_full  / ( CONJG(energy_complex)        + omegai(iomegai) &
&                                          - energy(iorbital,ispin) + wpol%pole(ipole)     ) )
             enddo
!           enddo
         enddo

       else
         stop'BUG'
       endif

     enddo !ipole

   enddo !iorbital
 enddo !ispin


 if(method==QS) then
   !
   ! Kotani's Hermitianization trick
   do ispin=1,nspin
     selfenergy_tmp(1,:,:,ispin) = 0.5_dp * REAL( ( selfenergy_tmp(1,:,:,ispin) + transpose(selfenergy_tmp(1,:,:,ispin)) ) )
   enddo
!   WRITE_MASTER(*,*) '  diagonal on the eigenvector basis'
!   WRITE_MASTER(*,*) '  #        e_GW         Sigc                   '
!   do aorbital=1,basis%nbf
!     WRITE_MASTER(*,'(i4,x,12(x,f12.6))') aorbital,energy(aorbital,:),REAL(selfenergy_tmp(1,aorbital,aorbital,:),dp)
!   enddo

   ! Transform the matrix elements back to the non interacting states
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,ispin) ) , MATMUL( selfenergy_tmp(1,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,ispin)), s_matrix(:,:) ) ) )
   enddo

 else if(method==perturbative) then

   if(file_exists) then
     open(13,file='selfenergy_omega')
     do iomegai=1,nomegai
       WRITE_MASTER(13,'(20(f12.6,2x))') DBLE(omegai(iomegai)),( DBLE(selfenergy_tmp(iomegai,aorbital,aorbital,:)), aorbital=1,5 )
     enddo
     close(13)
     stop'output the self energy in a file'
   endif
   selfenergy(:,:,:) = REAL( selfenergy_tmp(2,:,:,:) )
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,*) 'G0W0 Eigenvalues [Ha]'
   if(nspin==1) then
     WRITE_MASTER(*,*) '  #          E0         Sigc          Z         G0W0'
   else
     WRITE_MASTER(*,'(a)') '  #                E0                       Sigc                       Z                        G0W0'
   endif
   do aorbital=1,basis%nbf
     zz(:) = REAL( selfenergy_tmp(3,aorbital,aorbital,:) - selfenergy_tmp(1,aorbital,aorbital,:) ) / REAL( omegai(3)-omegai(1) )
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )

     WRITE_MASTER(*,'(i4,x,12(x,f12.6))') aorbital,energy(aorbital,:),REAL(selfenergy_tmp(2,aorbital,aorbital,:),dp),&
           zz(:),energy(aorbital,:)+zz(:)*REAL(selfenergy_tmp(2,aorbital,aorbital,:) + exchange_m_vxc_diag(aorbital,:) ,dp)
   enddo

   WRITE_MASTER(*,*)
   WRITE_MASTER(*,*) 'G0W0 Eigenvalues [eV]'
   if(nspin==1) then
     WRITE_MASTER(*,*) '  #          E0         Sigc          Z         G0W0'
   else
     WRITE_MASTER(*,'(a)') '  #                E0                       Sigc                       Z                        G0W0'
   endif
   do aorbital=1,basis%nbf
     zz(:) = REAL( selfenergy_tmp(3,aorbital,aorbital,:) - selfenergy_tmp(1,aorbital,aorbital,:) ) / REAL( omegai(3)-omegai(1) )
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )

     WRITE_MASTER(*,'(i4,x,12(x,f12.6))') aorbital,energy(aorbital,:)*Ha_eV,REAL(selfenergy_tmp(2,aorbital,aorbital,:),dp)*Ha_eV,&
           zz(:),( energy(aorbital,:)+zz(:)*REAL(selfenergy_tmp(2,aorbital,aorbital,:) + exchange_m_vxc_diag(aorbital,:) ,dp) )*Ha_eV
   enddo

 endif

 deallocate(omegai)
 deallocate(selfenergy_tmp)

end subroutine gw_selfenergy_aux


#endif
!=========================================================================
