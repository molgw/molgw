!=========================================================================
#include "macros.h"
!=========================================================================


!=========================================================================
subroutine polarizability_rpa(nspin,basis,prod_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none

 integer,intent(in)  :: nspin
 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin)
 real(dp),intent(in) :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out) :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol
!=====
 integer :: pbf,qbf,ibf,jbf,kbf,lbf,ijbf,klbf,ijbf_current,ijspin,klspin
 integer :: iorbital,jorbital,korbital,lorbital
 integer :: ipole
 integer :: t_ij,t_kl
 integer :: reading_status

#if defined LOW_MEMORY2 || defined LOW_MEMORY3
 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp),allocatable :: eri_eigenstate_k(:,:,:,:)
#else
 real(dp) :: eri_eigenstate(basis%nbf,basis%nbf,basis%nbf,basis%nbf,nspin,nspin)
#endif
 real(dp) :: spin_fact 
 real(dp) :: h_2p(wpol%npole,wpol%npole)
 real(dp) :: eigenvalue(wpol%npole),eigenvector(wpol%npole,wpol%npole),eigenvector_inv(wpol%npole,wpol%npole)
 real(dp) :: matrix(wpol%npole,wpol%npole)
 real(dp) :: rtmp
 real(dp) :: alpha1,alpha2

 logical :: TDHF=.FALSE.
!=====
 spin_fact = REAL(-nspin+3,dp)
 rpa_correlation = 0.0_dp

 WRITE_MASTER(*,'(/,a)') ' calculating CHI alla rpa'

 call read_spectral_function(wpol,reading_status)
 if( reading_status == 0 ) then
   WRITE_MASTER(*,'(a,/)') ' no need to calculate W: already done'
   return
 endif

 inquire(file='manual_tdhf',exist=TDHF)
 if(TDHF) then
   open(unit=18,file='manual_tdhf',status='old')
   read(18,*) alpha1
   read(18,*) alpha2
   close(18)
   WRITE_ME(msg,'(a,f12.6,3x,f12.6)') 'calculating the TDHF polarizability with alphas  ',alpha1,alpha2
   call issue_warning(msg)
 else
   alpha1=0.0_dp
   alpha2=0.0_dp
 endif


#if defined LOW_MEMORY2 || defined LOW_MEMORY3
 allocate(eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin))
#else
 call transform_eri_basis_fast(basis%nbf,nspin,c_matrix,eri_eigenstate)
#endif


 t_ij=0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

#if defined LOW_MEMORY2 || LOW_MEMORY3
 call start_clock(timing_tmp1)
       call transform_eri_basis_lowmem(nspin,c_matrix,iorbital,ijspin,eri_eigenstate_i)
 call stop_clock(timing_tmp1)
#ifdef CRPA
       if( iorbital==band1 .OR. iorbital==band2) then
         do jorbital=band1,band2
           do korbital=band1,band2
             do lorbital=band1,band2
               write(1000,'(4(i6,x),2x,f16.8)') iorbital,jorbital,korbital,lorbital,eri_eigenstate_i(jorbital,korbital,lorbital,1)
             enddo
           enddo
         enddo
       endif
#endif
#endif

     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
       t_ij=t_ij+1



       t_kl=0
       do klspin=1,nspin
         do korbital=1,basis%nbf 

           do lorbital=1,basis%nbf 
             if( skip_transition(nspin,lorbital,korbital,occupation(lorbital,klspin),occupation(korbital,klspin)) ) cycle
             t_kl=t_kl+1

#ifndef CHI0

#if defined LOW_MEMORY2 || defined LOW_MEMORY3

!!FBFB
!             if( t_kl /= t_ij .AND. & 
!                   ( iorbital >= 30 .OR. jorbital >= 30 .OR. korbital >= 30 .OR. lorbital >= 30 ) ) then
!               h_2p(t_ij,t_kl) = 0.0_dp
!             else
             h_2p(t_ij,t_kl) = eri_eigenstate_i(jorbital,korbital,lorbital,klspin) &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )
!             endif


#else

             h_2p(t_ij,t_kl) = eri_eigenstate(iorbital,jorbital,korbital,lorbital,ijspin,klspin) &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )
#endif

#else
             h_2p(t_ij,t_kl) = 0.0_dp
#endif


           if(TDHF) then
             if(ijspin==klspin) then
#if defined LOW_MEMORY2 || defined LOW_MEMORY3
               h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigenstate_i(korbital,jorbital,lorbital,klspin)  &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact * alpha1
#else
               h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigenstate(iorbital,korbital,jorbital,lorbital,ijspin,klspin)  &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact * alpha1
#endif
             endif
           endif

           enddo
         enddo
       enddo !klspin

       h_2p(t_ij,t_ij) =  h_2p(t_ij,t_ij) + ( energy(jorbital,ijspin) - energy(iorbital,ijspin) )

       rpa_correlation = rpa_correlation - 0.25_dp * ABS( h_2p(t_ij,t_ij) )

!       WRITE_MASTER(*,'(4(i4,2x),2(2x,f12.6))') t_ij,iorbital,jorbital,ijspin,( energy(jorbital,ijspin) - energy(iorbital,ijspin) ),h_2p(t_ij,t_ij)
     enddo !jorbital
   enddo !iorbital
 enddo ! ijspin

 WRITE_MASTER(*,*) 'diago 2-particle hamiltonian'
 WRITE_MASTER(*,*) 'matrix',wpol%npole,'x',wpol%npole
! WRITE_MASTER(*,*) '=============='
! t_ij=0
! do ijspin=1,nspin
!   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied
!     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
!!INTRA       if(iorbital==jorbital) cycle  ! intra state transitions are not allowed!
!!SKIP       if( abs(occupation(jorbital,ijspin)-occupation(iorbital,ijspin))<completely_empty ) cycle
!       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
!
!       t_ij=t_ij+1
!       WRITE_MASTER(*,'(3(i4,2x),20(2x,f12.6))') ijspin,iorbital,jorbital,h_2p(t_ij,:)
!     enddo
!   enddo
! enddo
! WRITE_MASTER(*,*) '=============='
! do t_ij=1,wpol%npole
!   WRITE_MASTER(*,'(1(i4,2x),20(2x,f12.6))') t_ij,h_2p(t_ij,:)
! enddo
 call start_clock(timing_diago_h2p)
 call diagonalize_general(wpol%npole,h_2p,eigenvalue,eigenvector)
 call stop_clock(timing_diago_h2p)
 WRITE_MASTER(*,*) 'diago finished'
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'calculate the RPA energy using the Tamm-Dancoff decomposition'
 WRITE_MASTER(*,*) 'formula (23) from F. Furche J. Chem. Phys. 129, 114105 (2008)'
 rpa_correlation = rpa_correlation + 0.25_dp * SUM( ABS(eigenvalue(:)) )
 WRITE_MASTER(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',rpa_correlation

 WRITE_MASTER(*,'(/,a,f14.8)') ' Lowest neutral excitation energy [eV]',MINVAL(ABS(eigenvalue(:)))*Ha_eV

! do t_ij=1,wpol%npole
!   WRITE_MASTER(*,'(1(i4,2x),20(2x,f12.6))') t_ij,eigenvalue(t_ij)
! enddo
   
 call start_clock(timing_inversion_s2p)
 call invert(wpol%npole,eigenvector,eigenvector_inv)
 call stop_clock(timing_inversion_s2p)

#if defined LOW_MEMORY2 || defined LOW_MEMORY3
 deallocate(eri_eigenstate_i)
 allocate(eri_eigenstate_k(basis%nbf,basis%nbf,basis%nbf,nspin))
#endif

 wpol%pole = eigenvalue
 wpol%residu_left (:,:) = 0.0_dp
 wpol%residu_right(:,:) = 0.0_dp
 t_kl=0
 do klspin=1,nspin
   do kbf=1,basis%nbf 
#if defined LOW_MEMORY2 || LOW_MEMORY3
 call start_clock(timing_tmp1)
     call transform_eri_basis_lowmem(nspin,c_matrix,kbf,klspin,eri_eigenstate_k)
 call stop_clock(timing_tmp1)
#endif


     do lbf=1,basis%nbf
       if( skip_transition(nspin,lbf,kbf,occupation(lbf,klspin),occupation(kbf,klspin)) ) cycle
       t_kl=t_kl+1


       do ijspin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(ibf,jbf,ijbf_current)
         do ijbf=1,prod_basis%nbf
           ibf = prod_basis%index_ij(1,ijbf)
           jbf = prod_basis%index_ij(2,ijbf)

           ijbf_current = ijbf+prod_basis%nbf*(ijspin-1)


#if defined LOW_MEMORY2 || defined LOW_MEMORY3
           wpol%residu_left (:,ijbf_current)  = wpol%residu_left (:,ijbf_current) &
                        + eri_eigenstate_k(lbf,ibf,jbf,ijspin) *  eigenvector(t_kl,:)
           wpol%residu_right(:,ijbf_current)  = wpol%residu_right(:,ijbf_current) &
                        + eri_eigenstate_k(lbf,ibf,jbf,ijspin) * eigenvector_inv(:,t_kl) &
                                         * ( occupation(kbf,klspin)-occupation(lbf,klspin) )
#else
           wpol%residu_left (:,ijbf_current)  = wpol%residu_left (:,ijbf_current) &
                        + eri_eigenstate(ibf,jbf,kbf,lbf,ijspin,klspin) *  eigenvector(t_kl,:)
           wpol%residu_right(:,ijbf_current)  = wpol%residu_right(:,ijbf_current) &
                        + eri_eigenstate(kbf,lbf,ibf,jbf,klspin,ijspin) * eigenvector_inv(:,t_kl) &
                                         * ( occupation(kbf,klspin)-occupation(lbf,klspin) )
#endif


         enddo
!$OMP END DO
!$OMP END PARALLEL
       enddo
     enddo
   enddo
 enddo


#ifdef CRPA
 do ijbf=1,prod_basis%nbf
   iorbital = prod_basis%index_ij(1,ijbf)
   jorbital = prod_basis%index_ij(2,ijbf)
   if(iorbital /=band1 .AND. iorbital /=band2) cycle
   if(jorbital /=band1 .AND. jorbital /=band2) cycle
   do klbf=1,prod_basis%nbf
     korbital = prod_basis%index_ij(1,klbf)
     lorbital = prod_basis%index_ij(2,klbf)
     if(korbital /=band1 .AND. korbital /=band2) cycle
     if(lorbital /=band1 .AND. lorbital /=band2) cycle
     rtmp=0.0_dp
     do ipole=1,wpol%npole
       rtmp = rtmp + wpol%residu_left(ipole,ijbf) * wpol%residu_right(ipole,klbf) / ( -wpol%pole(ipole) )
     enddo
     write(1001,'(4(i6,x),2x,f16.8)') iorbital,jorbital,korbital,lorbital,rtmp
   enddo
 enddo
#endif
 


#if defined LOW_MEMORY2 || defined LOW_MEMORY3
 deallocate(eri_eigenstate_k)
#endif

end subroutine polarizability_rpa


!=========================================================================
subroutine polarizability_casida(nspin,basis,prod_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none

 integer,intent(in)  :: nspin
 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin)
 real(dp),intent(in) :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out) :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol
!=====
 integer :: pbf,qbf,ibf,jbf,kbf,lbf,ijbf,klbf,ijbf_current,ijspin,klspin
 integer :: iorbital,jorbital,korbital,lorbital
 integer :: ipole
 integer :: t_ij,t_kl
 integer :: t_ij_local,t_kl_local
 integer :: transition_index(basis%nbf,basis%nbf,nspin)

 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp),allocatable :: eri_eigenstate_ij(:,:,:)

 real(dp) :: spin_fact 
! real(dp) :: apb(wpol%npole,wpol%npole)
! real(dp) :: amb_diag(wpol%npole)
 real(dp),allocatable :: apb(:,:)
 real(dp),allocatable :: amb_diag(:)
 real(dp) :: eigenvalue(wpol%npole)
! not used real(dp) :: eigenvector_inv(wpol%npole,wpol%npole)
! real(dp) :: eigenvector(wpol%npole,wpol%npole)

 logical :: TDHF=.FALSE.

 integer          :: mlocal,nlocal
#ifdef HAVE_SCALAPACK
 integer          :: desc(ndel)
#endif
!FIXME
 integer :: task(basis%nbf,nspin)
 integer :: itask
!=====

 spin_fact = REAL(-nspin+3,dp)
 rpa_correlation = 0.0_dp

#ifndef LOW_MEMORY2
#ifndef LOW_MEMORY3
 stop 'NOT implemented'
#endif
#endif

#ifdef HAVE_SCALAPACK
 call init_desc(wpol%npole,desc,mlocal,nlocal)
#else
 mlocal = wpol%npole
 nlocal = wpol%npole
#endif


 !
 ! The implementation closely follows the notation of F. Furche in JCP 132, 234114 (2010).
 !
 WRITE_MASTER(*,'(/,a)') ' calculating CHI a la Casida'
 if(TDHF) then
   msg='calculating the TDHF polarizability'
   call issue_warning(msg)
 endif

 allocate(eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin))

 !
 ! First set up the diagonal for the full matrix
 ! and transition indexing

 allocate(amb_diag(wpol%npole))
 amb_diag(:)=0.0_dp
 transition_index(:,:,:)=0
 t_ij=0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied
     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
       t_ij=t_ij+1
       transition_index(iorbital,jorbital,ijspin) = t_ij

       amb_diag(t_ij) = energy(jorbital,ijspin) - energy(iorbital,ijspin)

     enddo
   enddo
 enddo


 !
 ! Second set up the body of the local sub matrix (in scalapack sense)
 ! 
 allocate(apb(mlocal,nlocal))
 apb(:,:) = 0.0_dp

 !
 ! Prepare the SCALAPACK distribution of the ROWS
 !
 task(:,:) = 0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied
     if( occupation(iorbital,ijspin) < completely_empty .OR. iorbital <= ncore_W ) cycle
     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty

       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle

       t_ij=transition_index(iorbital,jorbital,ijspin)

       t_ij_local = rowindex_global_to_local(t_ij)
       if(t_ij_local==0) cycle

       task(iorbital,ijspin) = task(iorbital,ijspin) + 1

     enddo
   enddo
 enddo

#if 0
 WRITE_ME(111+rank,*) 'rank',rank,nproc
 do ijspin=1,nspin
   do iorbital=1,basis%nbf
     WRITE_ME(111+rank,*) task(iorbital,ijspin)
   enddo
 enddo
#endif

 !
 ! transitions ROWS
 !
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

     if( occupation(iorbital,ijspin) < completely_empty .OR. iorbital <= ncore_W ) cycle
     !
     ! SCALAPACK parallelization
     !
     if( task(iorbital,ijspin) == 0 ) cycle

     call transform_eri_basis_lowmem(nspin,c_matrix,iorbital,ijspin,eri_eigenstate_i)


     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty

       if( skip_transition(nspin,jorbital,iorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle

       t_ij=transition_index(iorbital,jorbital,ijspin)

       t_ij_local = rowindex_global_to_local(t_ij)
       if(t_ij_local==0) cycle

       !
       ! transitions COLS
       !
       do klspin=1,nspin
         do korbital=1,basis%nbf 

           do lorbital=1,basis%nbf 
             if( skip_transition(nspin,lorbital,korbital,occupation(lorbital,klspin),occupation(korbital,klspin)) ) cycle

             t_kl=transition_index(korbital,lorbital,klspin)

             t_kl_local = colindex_global_to_local(t_kl)
             if(t_kl_local==0) cycle

             apb(t_ij_local,t_kl_local) =  2.0_dp * eri_eigenstate_i(jorbital,korbital,lorbital,klspin) &
                                                      * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )


!TODO check TDHF implementation
!           if(TDHF) then
!             if(ijspin==klspin) then
!               h_2p(t_ij_local,t_kl_local) =  h_2p(t_ij_local,t_kl_local) -  eri_eigenstate_i(korbital,jorbital,lorbital,klspin)  &
!                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact 
!             endif
!           endif

             if(t_ij==t_kl) then
               if( task(iorbital,ijspin) == 0 ) then
                 WRITE_ME(*,'(a,10(2x,i5))') ' === should have skipped',rank,iorbital,ijspin,t_ij,t_ij_local
               endif
               apb(t_ij_local,t_kl_local) = apb(t_ij_local,t_kl_local) + amb_diag(t_ij)
               rpa_correlation = rpa_correlation - 0.25_dp * ( amb_diag(t_ij) + apb(t_ij_local,t_kl_local) )
!               WRITE_ME(*,'(a,5(2x,i5),2(2x,e16.6))') ' === should have skipped',rank,iorbital,t_ij,wpol%npole,t_ij_local,amb_diag(t_ij),apb(t_ij_local,t_kl_local)
             endif

           enddo
         enddo
       enddo ! klspin



     enddo !jorbital
   enddo !iorbital
 enddo ! ijspin


 !
 ! Reduce the rpa_correlation from all cpus
 !
 call sum_sca(rpa_correlation)

 !
 ! apb is transformed here into (a-b)**1/2 * (a+b) * (a-b)**1/2
 !
 do t_kl_local=1,nlocal
   t_kl = colindex_local_to_global(t_kl_local)

   do t_ij_local=1,mlocal
     t_ij = rowindex_local_to_global(t_ij_local)
     apb(t_ij_local,t_kl_local) = SQRT( amb_diag(t_ij) ) * apb(t_ij_local,t_kl_local) * SQRT( amb_diag(t_kl) )
   enddo

 enddo

 WRITE_MASTER(*,*) 'Diago Casida matrix'
 WRITE_MASTER(*,*) 'Matrix size:',wpol%npole,'x',wpol%npole

 !
 ! Symmetric in-place diagonalization
 call start_clock(timing_diago_h2p)

#ifdef HAVE_SCALAPACK

 call diagonalize_sca(desc,wpol%npole,mlocal,nlocal,apb,eigenvalue)

#else
! call diagonalize(wpol%npole,matrix,eigenvalue,eigenvector)
 call diagonalize_wo_vectors(wpol%npole,apb,eigenvalue)
#endif

 call stop_clock(timing_diago_h2p)

!!TEST sort eigenvalues
! WRITE_MASTER(*,*) 'eigenvalue'
! WRITE_MASTER(*,*) eigenvalue(:)
! t_kl=wpol%npole
! do t_ij=1,wpol%npole-1
!   if( ABS( eigenvalue(t_ij) - eigenvalue(t_ij+1) ) < 1.0e-5_dp ) then
!      write(*,*) 'remove',t_ij,t_ij+1
!      t_kl=t_kl-1
!   endif
! enddo
! write(*,*) 't_kl / npole',t_kl ,wpol%npole
!!END


 WRITE_MASTER(*,*) 'diago finished'
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'calculate the RPA energy using the Tamm-Dancoff decomposition'
 WRITE_MASTER(*,*) 'formula (9) from J. Chem. Phys. 132, 234114 (2010)'
 rpa_correlation = rpa_correlation + 0.50_dp * SUM( SQRT(eigenvalue(:)) )
 WRITE_MASTER(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',rpa_correlation




 deallocate(eri_eigenstate_i)
 deallocate(apb,amb_diag)

end subroutine polarizability_casida


!=========================================================================
subroutine gw_selfenergy(gwmethod,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_basis_set
 use m_spectral_function
 implicit none

 integer,intent(in)  :: gwmethod,nspin
 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin),energy(basis%nbf,nspin),exchange_m_vxc_diag(basis%nbf,nspin)
 real(dp),intent(in) :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in) :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out) :: selfenergy(basis%nbf,basis%nbf,nspin)

 logical               :: file_exists=.FALSE.
 logical               :: write_sigma_omega=.FALSE.
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 real(dp),allocatable  :: selfenergy_tmp(:,:,:,:)

 integer     :: bbf,ibf,kbf,lbf
 integer     :: aorbital,borbital
 integer     :: iorbital,ispin,ipole
 real(dp)    :: spin_fact,overlap_tmp
 real(dp)    :: bra(wpol%npole,basis%nbf),ket(wpol%npole,basis%nbf)
 real(dp)    :: fact_full,fact_empty
 real(dp)    :: zz(nspin)
 real(dp)    :: energy_re
 real(dp)    :: energy_qp(basis%nbf,nspin)
 character(len=3) :: ctmp
!=====
 spin_fact = REAL(-nspin+3,dp)

 WRITE_ME(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)

 WRITE_MASTER(*,*)
 select case(gwmethod)
 case(QS)
   WRITE_MASTER(*,*) 'perform a QP self-consistent GW calculation'
 case(perturbative)
   WRITE_MASTER(*,*) 'perform a one-shot G0W0 calculation'
 case(COHSEX)
   WRITE_MASTER(*,*) 'perform a COHSEX calculation'
 case(QSCOHSEX)
   WRITE_MASTER(*,*) 'perform a self-consistent COHSEX calculation'
 end select

 if(gwmethod==QS .OR. gwmethod==COHSEX .OR. gwmethod==QSCOHSEX) then
   nomegai=1
   allocate(omegai(nomegai))
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
     write_sigma_omega=.TRUE.
   else
     nomegai=3
     allocate(omegai(nomegai))
     omegai(1) = -0.01_dp
     omegai(2) =  0.00_dp
     omegai(3) =  0.01_dp
   endif
 endif

 allocate(selfenergy_tmp(nomegai,basis%nbf,basis%nbf,nspin))
 

 selfenergy_tmp(:,:,:,:) = 0.0_dp

 do ispin=1,nspin
   do iorbital=1,basis%nbf !INNER LOOP of G

     !
     ! Apply the frozen core and frozen virtual approximation to G
     if(iorbital <= ncore_G)    cycle
     if(iorbital >= nvirtual_G) cycle

     !
     ! Prepare the bra and ket with the knowledge of index iorbital and aorbital
     do aorbital=1,basis%nbf
       kbf = prod_basis%index_prodbasis(iorbital,aorbital)
       bra(:,aorbital) = wpol%residu_left (:,kbf+prod_basis%nbf*(ispin-1))
       ket(:,aorbital) = wpol%residu_right(:,kbf+prod_basis%nbf*(ispin-1))
     enddo

     do ipole=1,wpol%npole

       if( wpol%pole(ipole) < 0.0_dp ) then
         fact_empty = (spin_fact - occupation(iorbital,ispin)) / spin_fact
         fact_full = 0.0_dp
       else
         fact_empty = 0.0_dp
         fact_full = occupation(iorbital,ispin) / spin_fact
       endif
       if( ABS(fact_empty - fact_full) < 0.0001 ) cycle

       select case(gwmethod)
       case(QS)

         do borbital=1,basis%nbf
           do aorbital=1,basis%nbf

             selfenergy_tmp(1,aorbital,borbital,ispin) = selfenergy_tmp(1,aorbital,borbital,ispin) &
                        - bra(ipole,aorbital) * ket(ipole,borbital) &
                          * REAL(   fact_empty / ( energy(borbital,ispin) + ieta  - energy(iorbital,ispin) + wpol%pole(ipole) ) &
                                   -fact_full  / ( energy(borbital,ispin) - ieta  - energy(iorbital,ispin) + wpol%pole(ipole) ) , dp )

           enddo
         enddo

       case(perturbative)

         do aorbital=1,basis%nbf
           !
           ! calculate only the diagonal !
           borbital=aorbital
!           do borbital=1,basis%nbf
             do iomegai=1,nomegai
               selfenergy_tmp(iomegai,aorbital,borbital,ispin) = selfenergy_tmp(iomegai,aorbital,borbital,ispin) &
                        - bra(ipole,aorbital) * ket(ipole,borbital) &
                          * REAL(  fact_empty / ( energy(borbital,ispin) + ieta + omegai(iomegai) - energy(iorbital,ispin) + wpol%pole(ipole)     ) &
                                  -fact_full  / ( energy(borbital,ispin) - ieta + omegai(iomegai) - energy(iorbital,ispin) + wpol%pole(ipole)     )  , dp )
             enddo
!           enddo
         enddo

       case(COHSEX,QSCOHSEX) 

         !
         ! SEX
         !
         do borbital=1,basis%nbf
           do aorbital=1,basis%nbf
             selfenergy_tmp(1,aorbital,borbital,ispin) = selfenergy_tmp(1,aorbital,borbital,ispin) &
                        - bra(ipole,aorbital) * ket(ipole,borbital) &
                          * fact_full  / (-wpol%pole(ipole))    
           enddo
         enddo
         !
         ! COH
         !
         do borbital=1,basis%nbf
           do aorbital=1,basis%nbf
             selfenergy_tmp(1,aorbital,borbital,ispin) = selfenergy_tmp(1,aorbital,borbital,ispin) &
                        + bra(ipole,aorbital) * ket(ipole,borbital) &
                          * fact_empty  / (-wpol%pole(ipole))
           enddo
         enddo

       case default 
         stop'BUG'
       end select

     enddo !ipole

   enddo !iorbital
 enddo !ispin

 !
 ! Kotani's Hermitianization trick
 !
 if(gwmethod==QS) then
   do ispin=1,nspin
     selfenergy_tmp(1,:,:,ispin) = 0.5_dp * ( selfenergy_tmp(1,:,:,ispin) + transpose(selfenergy_tmp(1,:,:,ispin)) )
   enddo
 endif

 select case(gwmethod)
 case(QS)
   ! Transform the matrix elements back to the non interacting states
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,ispin) ) , MATMUL( selfenergy_tmp(1,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,ispin)), s_matrix(:,:) ) ) )
   enddo

 case(QSCOHSEX)
   ! Transform the matrix elements back to the non interacting states
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,ispin)) , MATMUL( selfenergy_tmp(1,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(COHSEX) !==========================================================

   selfenergy(:,:,:) = REAL( selfenergy_tmp(1,:,:,:) )
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,*) 'COHSEX Eigenvalues [eV]'
   if(nspin==1) then
     WRITE_MASTER(*,*) '  #          E0        Sigx-Vxc      Sigc          Z         G0W0'
   else
     WRITE_MASTER(*,'(a)') '  #                E0                      Sigx-Vxc                    Sigc                       Z                       G0W0'
   endif
   do aorbital=1,basis%nbf
     zz(:) = 1.0_dp 
     energy_qp(aorbital,:) = energy(aorbital,:)+zz(:)*REAL(selfenergy_tmp(1,aorbital,aorbital,:) + exchange_m_vxc_diag(aorbital,:))

     WRITE_MASTER(*,'(i4,x,20(x,f12.6))') aorbital,energy(aorbital,:)*Ha_eV,exchange_m_vxc_diag(aorbital,:)*Ha_eV,REAL(selfenergy_tmp(1,aorbital,aorbital,:),dp)*Ha_eV,&
           zz(:),energy_qp(aorbital,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp)

 case(perturbative) !==========================================================

   if(write_sigma_omega) then
     do aorbital=1,basis%nbf ! MIN(2,basis%nbf)
       WRITE_ME(ctmp,'(i3.3)') aorbital
       open(200+aorbital,file='selfenergy_omega_state'//TRIM(ctmp))
       do iomegai=1,nomegai
         WRITE_MASTER(200+aorbital,'(20(f12.6,2x))') ( DBLE(omegai(iomegai))+energy(aorbital,:) )*Ha_eV,&
                                                     ( DBLE(selfenergy_tmp(iomegai,aorbital,aorbital,:)) )*Ha_eV,&
                                                     ( DBLE(omegai(iomegai))-exchange_m_vxc_diag(aorbital,:) )*Ha_eV,&
                                                     ( 1.0_dp/pi/ABS( DBLE(omegai(iomegai))-exchange_m_vxc_diag(aorbital,:) - DBLE(selfenergy_tmp(iomegai,aorbital,aorbital,:)) ) ) / Ha_eV
       enddo
       WRITE_MASTER(200+aorbital,*)
     enddo
     close(200+aorbital)
     stop'output the self energy in a file'
   endif

   selfenergy(:,:,:) = REAL( selfenergy_tmp(2,:,:,:) )
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,*) 'G0W0 Eigenvalues [eV]'
   if(nspin==1) then
     WRITE_MASTER(*,*) '  #          E0        Sigx-Vxc      Sigc          Z         G0W0'
   else
     WRITE_MASTER(*,'(a)') '  #                E0                      Sigx-Vxc                    Sigc                       Z                       G0W0'
   endif
   do aorbital=1,basis%nbf
     zz(:) = REAL( selfenergy_tmp(3,aorbital,aorbital,:) - selfenergy_tmp(1,aorbital,aorbital,:) ) / REAL( omegai(3)-omegai(1) )
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )
     energy_qp(aorbital,:) = energy(aorbital,:)+zz(:)*REAL(selfenergy_tmp(2,aorbital,aorbital,:) + exchange_m_vxc_diag(aorbital,:))

     WRITE_MASTER(*,'(i4,x,20(x,f12.6))') aorbital,energy(aorbital,:)*Ha_eV,exchange_m_vxc_diag(aorbital,:)*Ha_eV,REAL(selfenergy_tmp(2,aorbital,aorbital,:),dp)*Ha_eV,&
           zz(:),energy_qp(aorbital,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp)

 end select

 deallocate(omegai)
 deallocate(selfenergy_tmp)


end subroutine gw_selfenergy


!=========================================================================
