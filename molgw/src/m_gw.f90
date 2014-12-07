!=========================================================================
#include "macros.h"
!=========================================================================
module m_gw

contains

!=========================================================================
subroutine polarizability_rpa(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_eri
 use m_spectral_function
 use m_inputparam
 implicit none

 type(basis_set)                       :: basis,prod_basis,auxil_basis
 real(dp),intent(in)                   :: occupation(basis%nbf,nspin)
 real(dp),intent(in)                   :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)                  :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol
!=====
 integer :: ibf,jbf,ijbf,klbf,ijspin,klspin
 integer :: istate,jstate,kstate,lstate
 integer :: ipole
 integer :: ntrans
 integer :: t_ij,t_kl
 integer :: reading_status

 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp)             :: eri_eigen_ijkl
 real(dp)             :: eri_eigen_ikjl
 real(dp),allocatable :: h_2p(:,:),eigenvector(:,:),eigenvector_inv(:,:)
 real(dp),allocatable :: eigenvalue(:)
 real(dp)             :: rtmp
 real(dp)             :: alpha1,alpha2

 logical :: TDHF=.FALSE.
!=====

 call start_clock(timing_pola)

 rpa_correlation = 0.0_dp

 WRITE_MASTER(*,'(/,a)') ' calculating CHI alla rpa'

 ! Obtain the number of transition = the size of the matrix
 call init_spectral_function(basis%nbf,occupation,ntrans)

 call read_spectral_function(wpol,reading_status)
 if( reading_status == 0 ) then
   WRITE_MASTER(*,'(a,/)') ' no need to calculate W: already done'
   return
 endif

 allocate(h_2p(ntrans,ntrans))
 allocate(eigenvector(ntrans,ntrans))
 allocate(eigenvalue(ntrans))

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


 if(is_auxil_basis) then
   call prepare_eri_3center_eigen(c_matrix)
 else
   allocate(eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin))
 endif


 call start_clock(timing_build_h2p)
 t_ij=0
 do ijspin=1,nspin
   do istate=1,basis%nbf ! istate stands for occupied or partially occupied

      if( .NOT. is_auxil_basis ) call transform_eri_basis(nspin,c_matrix,istate,ijspin,eri_eigenstate_i)


#ifdef CRPA
       ! attempt of coding CRPA
       if( istate==band1 .OR. istate==band2) then
         do jstate=band1,band2
           do kstate=band1,band2
             do lstate=band1,band2
               write(1000,'(4(i6,x),2x,f16.8)') istate,jstate,kstate,lstate,eri_eigenstate_i(jstate,kstate,lstate,1)
             enddo
           enddo
         enddo
       endif
#endif

     do jstate=1,basis%nbf ! jstate stands for empty or partially empty
       if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       t_ij=t_ij+1



       t_kl=0
       do klspin=1,nspin
         do kstate=1,basis%nbf 

           do lstate=1,basis%nbf 
             if( skip_transition(nspin,lstate,kstate,occupation(lstate,klspin),occupation(kstate,klspin)) ) cycle
             t_kl=t_kl+1

             if(is_auxil_basis) then
               eri_eigen_ijkl = eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
             else
               eri_eigen_ijkl = eri_eigenstate_i(jstate,kstate,lstate,klspin)
             endif
             
#ifndef CHI0

!!FBFB
!             if( t_kl /= t_ij .AND. & 
!                   ( istate >= 30 .OR. jstate >= 30 .OR. kstate >= 30 .OR. lstate >= 30 ) ) then
!               h_2p(t_ij,t_kl) = 0.0_dp
!             else
             h_2p(t_ij,t_kl) = eri_eigen_ijkl * ( occupation(istate,ijspin)-occupation(jstate,ijspin) )
!             endif


#else
             h_2p(t_ij,t_kl) = 0.0_dp
#endif


           if(TDHF) then
             if(ijspin==klspin) then
               if(is_auxil_basis) then
                 eri_eigen_ikjl = eri_eigen_ri(istate,kstate,ijspin,jstate,lstate,klspin)
               else
                 eri_eigen_ikjl = eri_eigenstate_i(kstate,jstate,lstate,klspin)
               endif
               h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigen_ikjl       &
                        * ( occupation(istate,ijspin)-occupation(jstate,ijspin) ) / spin_fact * alpha1
             endif
           endif

           enddo
         enddo
       enddo !klspin

       h_2p(t_ij,t_ij) =  h_2p(t_ij,t_ij) + ( energy(jstate,ijspin) - energy(istate,ijspin) )

       rpa_correlation = rpa_correlation - 0.25_dp * ABS( h_2p(t_ij,t_ij) )

     enddo !jstate
   enddo !istate
 enddo ! ijspin
 if(allocated(eri_eigenstate_i)) deallocate(eri_eigenstate_i)

 call stop_clock(timing_build_h2p)


 WRITE_MASTER(*,*) 'diago 2-particle hamiltonian'
 WRITE_MASTER(*,*) 'matrix',ntrans,'x',ntrans

 call start_clock(timing_diago_h2p)
 call diagonalize_general(ntrans,h_2p,eigenvalue,eigenvector)
 call stop_clock(timing_diago_h2p)
 WRITE_MASTER(*,*) 'diago finished'

 deallocate(h_2p)
 allocate(eigenvector_inv(ntrans,ntrans))

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'calculate the RPA energy using the Tamm-Dancoff decomposition'
 WRITE_MASTER(*,*) 'formula (23) from F. Furche J. Chem. Phys. 129, 114105 (2008)'
 rpa_correlation = rpa_correlation + 0.25_dp * SUM( ABS(eigenvalue(:)) )
 WRITE_MASTER(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',rpa_correlation

 WRITE_MASTER(*,'(/,a,f14.8)') ' Lowest neutral excitation energy [eV]',MINVAL(ABS(eigenvalue(:)))*Ha_eV

! do t_ij=1,ntrans
!   WRITE_MASTER(*,'(1(i4,2x),20(2x,f12.6))') t_ij,eigenvalue(t_ij)
! enddo
   
 call start_clock(timing_inversion_s2p)
 call invert(ntrans,eigenvector,eigenvector_inv)
 call stop_clock(timing_inversion_s2p)


 !
 ! Finally calculate v * \chi * v and store it in object wpol
 ! Deallocation is made inside chi_to_vchiv
 if(is_auxil_basis) then
   call chi_to_vchiv_auxil(ntrans,basis%nbf,auxil_basis%nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 else
   call chi_to_vchiv(ntrans,basis%nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 endif

 if(is_auxil_basis) call destroy_eri_3center_eigen()

 ! If requested write the spectral function on file
 if( print_specfunc ) call write_spectral_function(wpol)


#ifdef CRPA
 ! Constrained RPA attempt
 do ijbf=1,prod_basis%nbf
   istate = prod_basis%index_ij(1,ijbf)
   jstate = prod_basis%index_ij(2,ijbf)
   if(istate /=band1 .AND. istate /=band2) cycle
   if(jstate /=band1 .AND. jstate /=band2) cycle
   do klbf=1,prod_basis%nbf
     kstate = prod_basis%index_ij(1,klbf)
     lstate = prod_basis%index_ij(2,klbf)
     if(kstate /=band1 .AND. kstate /=band2) cycle
     if(lstate /=band1 .AND. lstate /=band2) cycle
     rtmp=0.0_dp
     do ipole=1,wpol%npole
       rtmp = rtmp + wpol%residu_left(ipole,ijbf) * wpol%residu_right(ipole,klbf) / ( -wpol%pole(ipole) )
     enddo
     write(1001,'(4(i6,x),2x,f16.8)') istate,jstate,kstate,lstate,rtmp
   enddo
 enddo
#endif
 
 call stop_clock(timing_pola)

end subroutine polarizability_rpa


!=========================================================================
subroutine chi_to_vchiv(ntrans,nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_inputparam,only: nspin,is_auxil_basis
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,ntrans
 type(basis_set),intent(in)            :: prod_basis
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 real(dp),intent(in)                   :: c_matrix(nbf,nbf,nspin)
 real(dp),allocatable,intent(inout)    :: eigenvector    (:,:)
 real(dp),allocatable,intent(inout)    :: eigenvector_inv(:,:)
 real(dp),allocatable,intent(inout)    :: eigenvalue     (:)
 type(spectral_function),intent(inout) :: wpol
!=====
 integer                    :: t_kl,klspin,ijspin
 integer                    :: istate,jstate,kstate,lstate,ijstate,ijstate_spin
 integer                    :: ipole
 real(dp)                   :: docc_kl
 real(dp)                   :: eri_eigen_klij
 real(dp),allocatable       :: eri_eigenstate_k(:,:,:,:)
!=====

 call start_clock(timing_buildw)

 if( .NOT. is_auxil_basis ) allocate(eri_eigenstate_k(nbf,nbf,nbf,nspin))

 call allocate_spectral_function(ntrans,prod_basis%nbf,wpol)

 wpol%pole(:) = eigenvalue(:)

 wpol%residu_left (:,:) = 0.0_dp
 wpol%residu_right(:,:) = 0.0_dp
 t_kl=0
 do klspin=1,nspin
   do kstate=1,nbf 

     if( .NOT. is_auxil_basis ) call transform_eri_basis(nspin,c_matrix,kstate,klspin,eri_eigenstate_k)

     do lstate=1,nbf
       if( skip_transition(nspin,lstate,kstate,occupation(lstate,klspin),occupation(kstate,klspin)) ) cycle
       t_kl=t_kl+1

       docc_kl = occupation(kstate,klspin)-occupation(lstate,klspin)

       do ijspin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(istate,jstate,ijstate_spin)
         do ijstate=1,prod_basis%nbf
           istate = prod_basis%index_ij(1,ijstate)
           jstate = prod_basis%index_ij(2,ijstate)

           ijstate_spin = ijstate+prod_basis%nbf*(ijspin-1)

           if(is_auxil_basis) then
             eri_eigen_klij = eri_eigen_ri(kstate,lstate,klspin,istate,jstate,ijspin)
           else
             eri_eigen_klij = eri_eigenstate_k(lstate,istate,jstate,ijspin)
           endif

           wpol%residu_left (:,ijstate_spin)  = wpol%residu_left (:,ijstate_spin) &
                        + eri_eigen_klij *  eigenvector(t_kl,:)
           wpol%residu_right(:,ijstate_spin)  = wpol%residu_right(:,ijstate_spin) &
                        + eri_eigen_klij * eigenvector_inv(:,t_kl) * docc_kl

         enddo
!$OMP END DO
!$OMP END PARALLEL
       enddo
     enddo
   enddo
 enddo

 deallocate(eigenvector)
 deallocate(eigenvector_inv)
 deallocate(eigenvalue)
 if(allocated(eri_eigenstate_k)) deallocate(eri_eigenstate_k)

 call stop_clock(timing_buildw)

end subroutine chi_to_vchiv


!=========================================================================
subroutine chi_to_vchiv_auxil(ntrans,nbf,nbf_auxil,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_inputparam,only: nspin,is_auxil_basis
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none
 
 integer,intent(in)         :: nbf,nbf_auxil,ntrans
 type(basis_set),intent(in) :: prod_basis
 real(dp),intent(in)        :: occupation(nbf,nspin)
 real(dp),intent(in)        :: c_matrix(nbf,nbf,nspin)
 real(dp),allocatable,intent(inout)    :: eigenvector    (:,:)
 real(dp),allocatable,intent(inout)    :: eigenvector_inv(:,:)
 real(dp),allocatable,intent(inout)    :: eigenvalue     (:)
 type(spectral_function),intent(inout) :: wpol
!=====
 integer                    :: t_ij,t_kl,klspin,ijspin
 integer                    :: istate,jstate,kstate,lstate,ijstate,ijstate_spin
 real(dp)                   :: docc_kl
 real(dp),allocatable       :: res_left(:,:),res_right(:,:)
!=====

 call start_clock(timing_buildw)

 allocate(res_left(nbf_auxil,ntrans))
 allocate(res_right(nbf_auxil,ntrans))

 res_left (:,:) = 0.0_dp
 res_right(:,:) = 0.0_dp
 t_kl=0
 do klspin=1,nspin
   do kstate=1,nbf 

     do lstate=1,nbf
       if( skip_transition(nspin,lstate,kstate,occupation(lstate,klspin),occupation(kstate,klspin)) ) cycle
       t_kl=t_kl+1

       docc_kl = occupation(kstate,klspin)-occupation(lstate,klspin)

       t_ij=0
       do ijspin=1,nspin
         do istate=1,nbf
      
           do jstate=1,nbf
             if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin))) cycle
             t_ij=t_ij+1


             res_left (:,t_ij) = res_left (:,t_ij) + eri_3center_eigen(:,kstate,lstate,klspin) * eigenvector(t_kl,t_ij)
             res_right(:,t_ij) = res_right(:,t_ij) + eri_3center_eigen(:,kstate,lstate,klspin) * eigenvector_inv(t_ij,t_kl) *docc_kl


           enddo
         enddo
       enddo

     enddo
   enddo
 enddo

 deallocate(eigenvector)
 deallocate(eigenvector_inv)

 call allocate_spectral_function(ntrans,prod_basis%nbf,wpol)

 wpol%pole(:) = eigenvalue(:)
 deallocate(eigenvalue)

 wpol%residu_left (:,:) = 0.0_dp
 wpol%residu_right(:,:) = 0.0_dp

 do ijspin=1,nspin
   do ijstate=1,prod_basis%nbf
     istate = prod_basis%index_ij(1,ijstate)
     jstate = prod_basis%index_ij(2,ijstate)
     ijstate_spin = ijstate+prod_basis%nbf*(ijspin-1)

     t_kl=0
     do klspin=1,nspin
       do kstate=1,nbf
    
         do lstate=1,nbf
           if( skip_transition(nspin,lstate,kstate,occupation(lstate,klspin),occupation(kstate,klspin))) cycle
           t_kl=t_kl+1

           wpol%residu_left (t_kl,ijstate_spin) = DOT_PRODUCT( res_left (:,t_kl) , eri_3center_eigen(:,istate,jstate,ijspin) )
           wpol%residu_right(t_kl,ijstate_spin) = DOT_PRODUCT( res_right(:,t_kl) , eri_3center_eigen(:,istate,jstate,ijspin) )

         enddo
       enddo

     enddo

   enddo
 enddo

 deallocate(res_left)
 deallocate(res_right)

 call stop_clock(timing_buildw)

end subroutine chi_to_vchiv_auxil


!=========================================================================
subroutine polarizability_casida(basis,prod_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_eri
 use m_spectral_function
 use m_inputparam,only: nspin,spin_fact
 implicit none

 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin)
 real(dp),intent(in) :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out) :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol
!=====
 integer :: pbf,qbf,kbf,lbf,klbf,ijspin,klspin
 integer :: istate,jstate,kstate,lstate
 integer :: ipole
 integer :: t_ij,t_kl
 integer :: t_ij_local,t_kl_local
 integer :: transition_index(basis%nbf,basis%nbf,nspin)

 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp),allocatable :: eri_eigenstate_ij(:,:,:)

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

 call start_clock(timing_pola)

 rpa_correlation = 0.0_dp

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

 call start_clock(timing_build_h2p)
 !
 ! First set up the diagonal for the full matrix
 ! and transition indexing

 allocate(amb_diag(wpol%npole))
 amb_diag(:)=0.0_dp
 transition_index(:,:,:)=0
 t_ij=0
 do ijspin=1,nspin
   do istate=1,basis%nbf ! istate stands for occupied or partially occupied
     do jstate=1,basis%nbf ! jstate stands for empty or partially empty
       if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle
       t_ij=t_ij+1
       transition_index(istate,jstate,ijspin) = t_ij

       amb_diag(t_ij) = energy(jstate,ijspin) - energy(istate,ijspin)

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
   do istate=1,basis%nbf ! istate stands for occupied or partially occupied
     if( occupation(istate,ijspin) < completely_empty .OR. istate <= ncore_W ) cycle
     do jstate=1,basis%nbf ! jstate stands for empty or partially empty

       if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle

       t_ij=transition_index(istate,jstate,ijspin)

       t_ij_local = rowindex_global_to_local(t_ij)
       if(t_ij_local==0) cycle

       task(istate,ijspin) = task(istate,ijspin) + 1

     enddo
   enddo
 enddo


 !
 ! transitions ROWS
 !
 do ijspin=1,nspin
   do istate=1,basis%nbf ! istate stands for occupied or partially occupied

     if( occupation(istate,ijspin) < completely_empty .OR. istate <= ncore_W ) cycle
     !
     ! SCALAPACK parallelization
     !
     if( task(istate,ijspin) == 0 ) cycle

     call transform_eri_basis(nspin,c_matrix,istate,ijspin,eri_eigenstate_i)


     do jstate=1,basis%nbf ! jstate stands for empty or partially empty

       if( skip_transition(nspin,jstate,istate,occupation(jstate,ijspin),occupation(istate,ijspin)) ) cycle

       t_ij=transition_index(istate,jstate,ijspin)

       t_ij_local = rowindex_global_to_local(t_ij)
       if(t_ij_local==0) cycle

       !
       ! transitions COLS
       !
       do klspin=1,nspin
         do kstate=1,basis%nbf 

           do lstate=1,basis%nbf 
             if( skip_transition(nspin,lstate,kstate,occupation(lstate,klspin),occupation(kstate,klspin)) ) cycle

             t_kl=transition_index(kstate,lstate,klspin)

             t_kl_local = colindex_global_to_local(t_kl)
             if(t_kl_local==0) cycle

             apb(t_ij_local,t_kl_local) =  2.0_dp * eri_eigenstate_i(jstate,kstate,lstate,klspin) &
                                                      * ( occupation(istate,ijspin)-occupation(jstate,ijspin) )


!TODO check TDHF implementation
!           if(TDHF) then
!             if(ijspin==klspin) then
!               h_2p(t_ij_local,t_kl_local) =  h_2p(t_ij_local,t_kl_local) -  eri_eigenstate_i(kstate,jstate,lstate,klspin)  &
!                        * ( occupation(istate,ijspin)-occupation(jstate,ijspin) ) / spin_fact 
!             endif
!           endif

             if(t_ij==t_kl) then
               if( task(istate,ijspin) == 0 ) then
!                 WRITE_ME(*,'(a,10(2x,i5))') ' === should have skipped',rank,istate,ijspin,t_ij,t_ij_local
               endif
               apb(t_ij_local,t_kl_local) = apb(t_ij_local,t_kl_local) + amb_diag(t_ij)
               rpa_correlation = rpa_correlation - 0.25_dp * ( amb_diag(t_ij) + apb(t_ij_local,t_kl_local) )
!               WRITE_ME(*,'(a,5(2x,i5),2(2x,e16.6))') ' === should have skipped',rank,istate,t_ij,wpol%npole,t_ij_local,amb_diag(t_ij),apb(t_ij_local,t_kl_local)
             endif

           enddo
         enddo
       enddo ! klspin



     enddo !jstate
   enddo !istate
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

 call stop_clock(timing_build_h2p)

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


 WRITE_MASTER(*,*) 'diago finished'
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'calculate the RPA energy using the Tamm-Dancoff decomposition'
 WRITE_MASTER(*,*) 'formula (9) from J. Chem. Phys. 132, 234114 (2010)'
 rpa_correlation = rpa_correlation + 0.50_dp * SUM( SQRT(eigenvalue(:)) )
 WRITE_MASTER(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',rpa_correlation




 deallocate(eri_eigenstate_i)
 deallocate(apb,amb_diag)

 call stop_clock(timing_pola)

end subroutine polarizability_casida


!=========================================================================
subroutine gw_selfenergy(gwmethod,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning,only: issue_warning
 use m_basis_set
 use m_spectral_function
 use m_inputparam,only: nspin,spin_fact
 implicit none

 integer,intent(in)  :: gwmethod
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

 integer     :: bbf,ibf,kbf
 integer     :: astate,bstate
 integer     :: istate,ispin,ipole
 real(dp)    :: overlap_tmp
 real(dp)    :: bra(wpol%npole,basis%nbf),ket(wpol%npole,basis%nbf)
 real(dp)    :: fact_full,fact_empty
 real(dp)    :: zz(nspin)
 real(dp)    :: energy_re
 real(dp)    :: energy_qp(basis%nbf,nspin)
 character(len=3) :: ctmp
!=====

 call start_clock(timing_self)


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
   do istate=1,basis%nbf !INNER LOOP of G

     !
     ! Apply the frozen core and frozen virtual approximation to G
     if(istate <= ncore_G)    cycle
     if(istate >= nvirtual_G) cycle

     !
     ! Prepare the bra and ket with the knowledge of index istate and astate
     do astate=1,basis%nbf
       kbf = prod_basis%index_prodbasis(istate,astate)
       bra(:,astate) = wpol%residu_left (:,kbf+prod_basis%nbf*(ispin-1))
       ket(:,astate) = wpol%residu_right(:,kbf+prod_basis%nbf*(ispin-1))
     enddo

     do ipole=1,wpol%npole

       if( wpol%pole(ipole) < 0.0_dp ) then
         fact_empty = (spin_fact - occupation(istate,ispin)) / spin_fact
         fact_full = 0.0_dp
       else
         fact_empty = 0.0_dp
         fact_full = occupation(istate,ispin) / spin_fact
       endif
       if( ABS(fact_empty - fact_full) < 0.0001 ) cycle

       select case(gwmethod)
       case(QS)

         do bstate=1,basis%nbf
           do astate=1,basis%nbf

             selfenergy_tmp(1,astate,bstate,ispin) = selfenergy_tmp(1,astate,bstate,ispin) &
                        - bra(ipole,astate) * ket(ipole,bstate) &
                          * REAL(   fact_empty / ( energy(bstate,ispin) + ieta  - energy(istate,ispin) + wpol%pole(ipole) ) &
                                   -fact_full  / ( energy(bstate,ispin) - ieta  - energy(istate,ispin) + wpol%pole(ipole) ) , dp )

           enddo
         enddo

       case(perturbative)

         do astate=1,basis%nbf
           !
           ! calculate only the diagonal !
           bstate=astate
!           do bstate=1,basis%nbf
             do iomegai=1,nomegai
               selfenergy_tmp(iomegai,astate,bstate,ispin) = selfenergy_tmp(iomegai,astate,bstate,ispin) &
                        - bra(ipole,astate) * ket(ipole,bstate) &
                          * REAL(  fact_empty / ( energy(bstate,ispin) + ieta + omegai(iomegai) - energy(istate,ispin) + wpol%pole(ipole)     ) &
                                  -fact_full  / ( energy(bstate,ispin) - ieta + omegai(iomegai) - energy(istate,ispin) + wpol%pole(ipole)     )  , dp )
             enddo
!           enddo
         enddo

       case(COHSEX,QSCOHSEX) 

         !
         ! SEX
         !
         do bstate=1,basis%nbf
           do astate=1,basis%nbf
             selfenergy_tmp(1,astate,bstate,ispin) = selfenergy_tmp(1,astate,bstate,ispin) &
                        - bra(ipole,astate) * ket(ipole,bstate) &
                          * fact_full  / (-wpol%pole(ipole))    
           enddo
         enddo
         !
         ! COH
         !
         do bstate=1,basis%nbf
           do astate=1,basis%nbf
             selfenergy_tmp(1,astate,bstate,ispin) = selfenergy_tmp(1,astate,bstate,ispin) &
                        + bra(ipole,astate) * ket(ipole,bstate) &
                          * fact_empty  / (-wpol%pole(ipole))
           enddo
         enddo

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
   do astate=1,basis%nbf
     zz(:) = 1.0_dp 
     energy_qp(astate,:) = energy(astate,:)+zz(:)*REAL(selfenergy_tmp(1,astate,astate,:) + exchange_m_vxc_diag(astate,:))

     WRITE_MASTER(*,'(i4,x,20(x,f12.6))') astate,energy(astate,:)*Ha_eV,exchange_m_vxc_diag(astate,:)*Ha_eV,REAL(selfenergy_tmp(1,astate,astate,:),dp)*Ha_eV,&
           zz(:),energy_qp(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp)

 case(perturbative) !==========================================================

   if(write_sigma_omega) then

     do astate=1,basis%nbf ! MIN(2,basis%nbf)
       WRITE_ME(ctmp,'(i3.3)') astate
       open(200+astate,file='selfenergy_omega_state'//TRIM(ctmp))
       do iomegai=1,nomegai
         WRITE_MASTER(200+astate,'(20(f12.6,2x))') ( DBLE(omegai(iomegai))+energy(astate,:) )*Ha_eV,&
                                                     ( DBLE(selfenergy_tmp(iomegai,astate,astate,:)) )*Ha_eV,&
                                                     ( DBLE(omegai(iomegai))-exchange_m_vxc_diag(astate,:) )*Ha_eV,&
                                                     ( 1.0_dp/pi/ABS( DBLE(omegai(iomegai))-exchange_m_vxc_diag(astate,:) &
                                                                   - DBLE(selfenergy_tmp(iomegai,astate,astate,:)) ) ) / Ha_eV
       enddo
       WRITE_MASTER(200+astate,*)
     enddo
     close(200+astate)

   else

     selfenergy(:,:,:) = REAL( selfenergy_tmp(2,:,:,:) )
     WRITE_MASTER(*,*)
     WRITE_MASTER(*,*) 'G0W0 Eigenvalues [eV]'
     if(nspin==1) then
       WRITE_MASTER(*,*) '  #          E0        Sigx-Vxc      Sigc          Z         G0W0'
     else
       WRITE_MASTER(*,'(a)') '  #                E0                      Sigx-Vxc                    Sigc                       Z                       G0W0'
     endif
     do astate=1,basis%nbf
       zz(:) = REAL( selfenergy_tmp(3,astate,astate,:) - selfenergy_tmp(1,astate,astate,:) ) / REAL( omegai(3)-omegai(1) )
       zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )
       energy_qp(astate,:) = energy(astate,:)+zz(:)*REAL(selfenergy_tmp(2,astate,astate,:) + exchange_m_vxc_diag(astate,:))

       WRITE_MASTER(*,'(i4,x,20(x,f12.6))') astate,energy(astate,:)*Ha_eV,exchange_m_vxc_diag(astate,:)*Ha_eV,REAL(selfenergy_tmp(2,astate,astate,:),dp)*Ha_eV,&
             zz(:),energy_qp(astate,:)*Ha_eV
     enddo

     call write_energy_qp(nspin,basis%nbf,energy_qp)
   endif

 end select

 deallocate(omegai)
 deallocate(selfenergy_tmp)

 call stop_clock(timing_self)

end subroutine gw_selfenergy


end module m_gw
!=========================================================================
