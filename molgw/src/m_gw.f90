!=========================================================================
#include "macros.h"
!=========================================================================
module m_gw
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam

contains

!=========================================================================
subroutine polarizability_rpa(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_warning,only: issue_warning
 use m_basis_set
 use m_spectral_function
 use m_tools
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
 integer :: t_ij
 integer :: reading_status

 real(dp),allocatable :: eigenvector(:,:),eigenvector_inv(:,:)
 real(dp),allocatable :: eigenvalue(:)
 real(dp)             :: rtmp
!=====

 call start_clock(timing_pola)

 WRITE_MASTER(*,'(/,a)') ' calculating CHI alla rpa'

 ! Obtain the number of transition = the size of the matrix
 call init_spectral_function(basis%nbf,occupation,wpol)

 call read_spectral_function(wpol,reading_status)
 if( reading_status == 0 ) then
   WRITE_MASTER(*,'(a,/)') ' no need to calculate W: already done'
   return
 endif

 !
 ! Build the 2-particle hamiltonian in transition space
 !
 allocate(eigenvector(wpol%npole,wpol%npole))
 allocate(eigenvalue(wpol%npole))

 select case(transition_ordering)
 case(SAVE_CPU)
   call build_h2p    (basis%nbf,c_matrix,occupation,energy,wpol,eigenvalue,eigenvector,rpa_correlation)
 case(CASIDA)
   call build_h2p_sym(basis%nbf,c_matrix,occupation,energy,wpol,eigenvalue,eigenvector,rpa_correlation)
 end select

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'Calculate RPA correlation energy'
 WRITE_MASTER(*,*) 'with Eq. (23) of F. Furche J. Chem. Phys. 129, 114105 (2008)'
 WRITE_MASTER(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',rpa_correlation


 allocate(eigenvector_inv(wpol%npole,wpol%npole))


 WRITE_MASTER(*,'(/,a,f14.8)') ' Lowest neutral excitation energy [eV]',MINVAL(ABS(eigenvalue(:)))*Ha_eV

! do t_ij=1,wpol%npole
!   WRITE_MASTER(124,'(1(i4,2x),20(2x,f12.6))') t_ij,eigenvalue(t_ij)*Ha_eV
! enddo
   
 call start_clock(timing_inversion_s2p)
 call invert(wpol%npole,eigenvector,eigenvector_inv)
 call stop_clock(timing_inversion_s2p)


 !
 ! Finally calculate v * \chi * v and store it in object wpol
 ! Deallocation is made inside chi_to_vchiv
 if(is_auxil_basis) then
   call chi_to_sqrtvchisqrtv_auxil(basis%nbf,auxil_basis%nbf,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 else
   call chi_to_vchiv(basis%nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 endif

 ! If requested write the spectral function on file
 if( print_specfunc ) call write_spectral_function(wpol)

 call stop_clock(timing_pola)

end subroutine polarizability_rpa


!=========================================================================
subroutine chi_to_vchiv(nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf
 type(basis_set),intent(in)            :: prod_basis
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 real(dp),intent(in)                   :: c_matrix(nbf,nbf,nspin)
 real(dp),allocatable,intent(inout)    :: eigenvector    (:,:)
 real(dp),allocatable,intent(inout)    :: eigenvector_inv(:,:)
 real(dp),allocatable,intent(inout)    :: eigenvalue     (:)
 type(spectral_function),intent(inout) :: wpol
!=====
 integer                               :: t_kl,klspin,ijspin
 integer                               :: istate,jstate,kstate,lstate,ijstate,ijstate_spin
 integer                               :: klstate_min
 integer                               :: klstate_max
 integer                               :: ipole
 real(dp)                              :: docc_kl
 real(dp)                              :: eri_eigen_klij
 real(dp),allocatable                  :: eri_eigenstate_klmin(:,:,:,:)
!=====

 call start_clock(timing_buildw)

 WRITE_MASTER(*,'(/,a)') ' Build W = v * chi * v'

 if( .NOT. is_auxil_basis ) then
   allocate(eri_eigenstate_klmin(nbf,nbf,nbf,nspin))
   ! Set this to zero and then enforce the calculation of the first array of Coulomb integrals
   eri_eigenstate_klmin(:,:,:,:) = 0.0_dp
 endif

 call allocate_spectral_function(prod_basis%nbf*nspin,wpol)

 wpol%pole(:) = eigenvalue(:)

 wpol%residu_left (:,:) = 0.0_dp
 wpol%residu_right(:,:) = 0.0_dp
 do t_kl=1,wpol%npole
   kstate = wpol%transition_table(1,t_kl)
   lstate = wpol%transition_table(2,t_kl)
   klspin = wpol%transition_table(3,t_kl)
   docc_kl= occupation(kstate,klspin)-occupation(lstate,klspin)

   if( .NOT. is_auxil_basis ) then
     klstate_min = MIN(kstate,lstate)
     klstate_max = MAX(kstate,lstate)
     call transform_eri_basis(nspin,c_matrix,klstate_min,klspin,eri_eigenstate_klmin)
   endif

   do ijspin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(istate,jstate,ijstate_spin,eri_eigen_klij)
     do ijstate=1,prod_basis%nbf
       istate = prod_basis%index_ij(1,ijstate)
       jstate = prod_basis%index_ij(2,ijstate)

       ijstate_spin = ijstate+prod_basis%nbf*(ijspin-1)

       if(is_auxil_basis) then
         eri_eigen_klij = eri_eigen_ri(kstate,lstate,klspin,istate,jstate,ijspin)
       else
         eri_eigen_klij = eri_eigenstate_klmin(klstate_max,istate,jstate,ijspin)
       endif

       wpol%residu_left (:,ijstate_spin)  = wpol%residu_left (:,ijstate_spin) &
                    + eri_eigen_klij * eigenvector(t_kl,:)
       wpol%residu_right(:,ijstate_spin)  = wpol%residu_right(:,ijstate_spin) &
                    + eri_eigen_klij * eigenvector_inv(:,t_kl) * docc_kl

     enddo
!$OMP END DO
!$OMP END PARALLEL
   enddo
 enddo

 deallocate(eigenvector)
 deallocate(eigenvector_inv)
 deallocate(eigenvalue)
 if(allocated(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)

 call stop_clock(timing_buildw)

end subroutine chi_to_vchiv


!=========================================================================
subroutine chi_to_sqrtvchisqrtv_auxil(nbf,nbf_auxil,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,nbf_auxil
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 real(dp),intent(in)                   :: c_matrix(nbf,nbf,nspin)
 real(dp),allocatable,intent(inout)    :: eigenvector    (:,:)
 real(dp),allocatable,intent(inout)    :: eigenvector_inv(:,:)
 real(dp),allocatable,intent(inout)    :: eigenvalue     (:)
 type(spectral_function),intent(inout) :: wpol
!=====
 integer                               :: t_ij,t_kl,klspin,ijspin
 integer                               :: kstate,lstate
 real(dp)                              :: docc_kl
 real(dp),allocatable                  :: eri_3center_2index(:,:)
!=====

 call start_clock(timing_buildw)

 WRITE_MASTER(*,'(/,a)') ' Build v^{1/2} * chi * v^{1/2}'

 call allocate_spectral_function(nbf_auxil,wpol)
 wpol%pole(:) = eigenvalue(:)
 deallocate(eigenvalue)

 allocate(eri_3center_2index(nbf_auxil,wpol%npole))
 do t_kl=1,wpol%npole
   kstate = wpol%transition_table(1,t_kl)
   lstate = wpol%transition_table(2,t_kl)
   klspin = wpol%transition_table(3,t_kl)
   eri_3center_2index(:,t_kl) = eri_3center_eigen(:,kstate,lstate,klspin)
 enddo
 wpol%residu_left(:,:)  = TRANSPOSE( MATMUL( eri_3center_2index , eigenvector ) )

 do t_kl=1,wpol%npole
   kstate = wpol%transition_table(1,t_kl)
   lstate = wpol%transition_table(2,t_kl)
   klspin = wpol%transition_table(3,t_kl)
   docc_kl = occupation(kstate,klspin)-occupation(lstate,klspin)
   eri_3center_2index(:,t_kl) = eri_3center_2index(:,t_kl) * docc_kl
 enddo
 wpol%residu_right(:,:) = MATMUL( eigenvector_inv , TRANSPOSE( eri_3center_2index ) )

 deallocate(eri_3center_2index)


 deallocate(eigenvector)
 deallocate(eigenvector_inv)

 call stop_clock(timing_buildw)

end subroutine chi_to_sqrtvchisqrtv_auxil


!=========================================================================
subroutine polarizability_casida(basis,prod_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_definitions
 use m_mpi
 use m_timing 
 use m_warning,only: issue_warning
 use m_tools
 use m_basis_set
 use m_eri
 use m_spectral_function
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

 allocate(eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin))

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
 use m_timing 
 use m_warning,only: issue_warning,msg
 use m_basis_set
 use m_spectral_function
 use m_eri,only: eri_3center_eigen
 implicit none

 integer,intent(in)                 :: gwmethod
 type(basis_set)                    :: basis,prod_basis
 real(dp),intent(in)                :: occupation(basis%nbf,nspin),energy(basis%nbf,nspin),exchange_m_vxc_diag(basis%nbf,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
!=====
 logical               :: file_exists=.FALSE.
 logical               :: write_sigma_omega=.FALSE.
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 real(dp),allocatable  :: selfenergy_tmp(:,:,:,:)
 integer               :: bbf,ibf,kbf
 integer               :: astate,bstate
 integer               :: istate,ispin,ipole
 real(dp)              :: overlap_tmp
 real(dp)              :: bra(wpol%npole,basis%nbf),ket(wpol%npole,basis%nbf)
 real(dp)              :: fact_full,fact_empty
 real(dp)              :: zz(nspin)
 real(dp)              :: energy_re
 real(dp)              :: energy_qp(basis%nbf,nspin)
 character(len=3)      :: ctmp
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
     if( .NOT. is_auxil_basis) then
       ! Here just grab the precalculated value
       do astate=1,basis%nbf
         kbf = prod_basis%index_prodbasis(istate,astate) + prod_basis%nbf*(ispin-1)
         bra(:,astate) = wpol%residu_left (:,kbf)
         ket(:,astate) = wpol%residu_right(:,kbf)
       enddo
     else
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,:) = MATMUL( wpol%residu_left (:,:) , eri_3center_eigen(:,:,istate,ispin) )
       ket(:,:) = MATMUL( wpol%residu_right(:,:) , eri_3center_eigen(:,:,istate,ispin) )
     endif

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

     do astate=1,basis%nbf
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
       ! Contrain Z to be in [0:1] to avoid crazy values
       do ispin=1,nspin
         zz(ispin) = MIN( MAX(zz(ispin),0.0_dp) , 1.0_dp )
       enddo
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


!=========================================================================
subroutine build_h2p(nbf,c_matrix,occupation,energy,wpol,eigenvalue,eigenvector,rpa_correlation)
 use m_tools
 use m_spectral_function
 use m_eri
 implicit none

 integer,intent(in)                 :: nbf
 real(dp),intent(in)                :: occupation(nbf,nspin),energy(nbf,nspin)
 real(dp),intent(in)                :: c_matrix(nbf,nbf,nspin)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: eigenvalue(wpol%npole),eigenvector(wpol%npole,wpol%npole)
 real(dp),intent(out)               :: rpa_correlation
!=====
 integer              :: t_ij,t_kl
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijstate_min
 integer              :: ijspin,klspin
 real(dp),allocatable :: eri_eigenstate_ijmin(:,:,:,:)
 real(dp),allocatable :: h_2p(:,:)
 real(dp)             :: eri_eigen_ijkl
 real(dp)             :: eri_eigen_ikjl
 real(dp)             :: alpha,docc_ij

 logical              :: is_ij
 logical              :: TDHF=.FALSE.
!=====

 call start_clock(timing_build_h2p)

 inquire(file='manual_tdhf',exist=TDHF)
 if(TDHF) then
   open(unit=18,file='manual_tdhf',status='old')
   read(18,*) alpha
   close(18)
   WRITE_ME(msg,'(a,f12.6,3x,f12.6)') 'calculating the TDHF polarizability with alphas  ',alpha
   call issue_warning(msg)
 else
   alpha=0.0_dp
 endif


 if( .NOT. is_auxil_basis) then
   allocate(eri_eigenstate_ijmin(nbf,nbf,nbf,nspin))
   ! Set this to zero and then enforce the calculation of the first series of Coulomb integrals
   eri_eigenstate_ijmin(:,:,:,:) = 0.0_dp
 endif

 allocate(h_2p(wpol%npole,wpol%npole))

 do t_ij=1,wpol%npole
   istate = wpol%transition_table(1,t_ij)
   jstate = wpol%transition_table(2,t_ij)
   ijspin = wpol%transition_table(3,t_ij)
   docc_ij= occupation(istate,ijspin)-occupation(jstate,ijspin)

   if( .NOT. is_auxil_basis ) then
     ijstate_min = MIN(istate,jstate)
     is_ij = (ijstate_min == istate)
     call transform_eri_basis(nspin,c_matrix,ijstate_min,ijspin,eri_eigenstate_ijmin)
   endif

   do t_kl=1,wpol%npole
     kstate = wpol%transition_table(1,t_kl)
     lstate = wpol%transition_table(2,t_kl)
     klspin = wpol%transition_table(3,t_kl)

     if(is_auxil_basis) then
       eri_eigen_ijkl = eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
     else
       if(is_ij) then ! treating (i,j)
         eri_eigen_ijkl = eri_eigenstate_ijmin(jstate,kstate,lstate,klspin)
       else           ! treating (j,i)
         eri_eigen_ijkl = eri_eigenstate_ijmin(istate,kstate,lstate,klspin)
       endif
     endif


     h_2p(t_ij,t_kl) = eri_eigen_ijkl * docc_ij

     if(TDHF) then
       if(ijspin==klspin) then
         if(is_auxil_basis) then
           eri_eigen_ikjl = eri_eigen_ri(istate,kstate,ijspin,jstate,lstate,klspin)
         else
           if(is_ij) then
             eri_eigen_ikjl = eri_eigenstate_ijmin(kstate,jstate,lstate,klspin)
           else
             eri_eigen_ikjl = eri_eigenstate_ijmin(lstate,istate,kstate,klspin)
           endif
         endif
         h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigen_ikjl       &
                                                 * docc_ij / spin_fact * alpha
       endif
     endif

   enddo ! t_kl

   h_2p(t_ij,t_ij) =  h_2p(t_ij,t_ij) + ( energy(jstate,ijspin) - energy(istate,ijspin) )

 enddo ! t_ij

 if(allocated(eri_eigenstate_ijmin)) deallocate(eri_eigenstate_ijmin)

 call stop_clock(timing_build_h2p)

 rpa_correlation = 0.0_dp
 do t_ij=1,wpol%npole
   rpa_correlation = rpa_correlation - 0.25_dp * ABS( h_2p(t_ij,t_ij) )
 enddo

 WRITE_MASTER(*,*) 'diago 2-particle hamiltonian'
 WRITE_MASTER(*,*) 'matrix',wpol%npole,'x',wpol%npole

 call start_clock(timing_diago_h2p)
 call diagonalize_general(wpol%npole,h_2p,eigenvalue,eigenvector)
 call stop_clock(timing_diago_h2p)
 WRITE_MASTER(*,*) 'diago finished'

 deallocate(h_2p)

 rpa_correlation = rpa_correlation + 0.25_dp * SUM( ABS(eigenvalue(:)) )


end subroutine build_h2p


!=========================================================================
subroutine build_h2p_sym(nbf,c_matrix,occupation,energy,wpol,eigenvalue,eigenvector,rpa_correlation)
 use m_spectral_function
 use m_eri
 use m_tools 
 implicit none

 integer,intent(in)                 :: nbf
 real(dp),intent(in)                :: occupation(nbf,nspin),energy(nbf,nspin)
 real(dp),intent(in)                :: c_matrix(nbf,nbf,nspin)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: eigenvalue(wpol%npole)
 real(dp),intent(out)               :: eigenvector(wpol%npole,wpol%npole)
 real(dp),intent(out)               :: rpa_correlation
!=====
 integer              :: t_ij,t_kl
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
 integer              :: ijstate_min
 real(dp),allocatable :: eri_eigenstate_ijmin(:,:,:,:)
 real(dp)             :: eri_eigen_ijkl
 real(dp)             :: eri_eigen_ikjl,eri_eigen_iljk
 real(dp)             :: alpha
 real(dp),allocatable :: apb_matrix(:,:),amb_matrix(:,:),amb_eigval(:),bigomega(:)
 real(dp),allocatable :: amb_matrix_sqrt(:,:),amb_matrix_sqrtm1(:,:)
 real(dp),allocatable :: cc_matrix(:,:)
 integer              :: nmat
 logical              :: is_ij
! real(dp),allocatable :: bigx(:),bigy(:)
 real(dp),allocatable :: bigx(:,:),bigy(:,:),cc_matrix_bigomega(:,:)

 logical              :: TDHF=.FALSE.
!=====

 call start_clock(timing_build_h2p)

 inquire(file='manual_tdhf',exist=TDHF)
 if(TDHF) then
   open(unit=18,file='manual_tdhf',status='old')
   read(18,*) alpha
   close(18)
   WRITE_ME(msg,'(a,f12.6,3x,f12.6)') 'calculating the TDHF polarizability with alphas  ',alpha
   call issue_warning(msg)
 else
   alpha=0.0_dp
 endif

 nmat=wpol%npole/2
 allocate(apb_matrix(nmat,nmat))
 allocate(amb_matrix(nmat,nmat))

 if( .NOT. is_auxil_basis) then
   allocate(eri_eigenstate_ijmin(nbf,nbf,nbf,nspin))
   ! Set this to zero and then enforce the calculation of the first series of
   ! Coulomb integrals
   eri_eigenstate_ijmin(:,:,:,:) = 0.0_dp
 endif

 do t_ij=1,nmat ! only resonant transition
   istate = wpol%transition_table(1,t_ij)
   jstate = wpol%transition_table(2,t_ij)
   ijspin = wpol%transition_table(3,t_ij)

   if( .NOT. is_auxil_basis ) then
     ijstate_min = MIN(istate,jstate)
     is_ij = (ijstate_min == istate)
     call transform_eri_basis(nspin,c_matrix,ijstate_min,ijspin,eri_eigenstate_ijmin)
   endif

   do t_kl=1,nmat ! only resonant transition
     kstate = wpol%transition_table(1,t_kl)
     lstate = wpol%transition_table(2,t_kl)
     klspin = wpol%transition_table(3,t_kl)


     if(is_auxil_basis) then
       eri_eigen_ijkl = eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
     else
       if(is_ij) then ! treating (i,j)
         eri_eigen_ijkl = eri_eigenstate_ijmin(jstate,kstate,lstate,klspin)
       else           ! treating (j,i)
         eri_eigen_ijkl = eri_eigenstate_ijmin(istate,kstate,lstate,klspin)
       endif
     endif

     apb_matrix(t_ij,t_kl) = 2.0_dp * eri_eigen_ijkl * spin_fact
     amb_matrix(t_ij,t_kl) = 0.0_dp

     if(TDHF) then
       if(ijspin==klspin) then
         if(is_auxil_basis) then
           eri_eigen_ikjl = eri_eigen_ri(istate,kstate,ijspin,jstate,lstate,klspin)
           eri_eigen_iljk = eri_eigen_ri(istate,lstate,ijspin,jstate,lstate,klspin)
         else
           if(is_ij) then
             eri_eigen_ikjl = eri_eigenstate_ijmin(kstate,jstate,lstate,klspin)
             eri_eigen_iljk = eri_eigenstate_ijmin(lstate,jstate,kstate,klspin)
           else
             eri_eigen_ikjl = eri_eigenstate_ijmin(lstate,istate,kstate,klspin)
             eri_eigen_iljk = eri_eigenstate_ijmin(kstate,istate,lstate,klspin)
           endif
         endif
         apb_matrix(t_ij,t_kl) = apb_matrix(t_ij,t_kl) - eri_eigen_ikjl * alpha - eri_eigen_iljk * alpha
         amb_matrix(t_ij,t_kl) = amb_matrix(t_ij,t_kl) - eri_eigen_ikjl * alpha + eri_eigen_iljk * alpha
       endif
     endif

   enddo 

   apb_matrix(t_ij,t_ij) =  apb_matrix(t_ij,t_ij) + ( energy(jstate,ijspin) - energy(istate,ijspin) )
   amb_matrix(t_ij,t_ij) =  amb_matrix(t_ij,t_ij) + ( energy(jstate,ijspin) - energy(istate,ijspin) )

 enddo 

 if(allocated(eri_eigenstate_ijmin)) deallocate(eri_eigenstate_ijmin)
 call stop_clock(timing_build_h2p)


 do t_ij=1,nmat
   rpa_correlation = rpa_correlation - 0.25_dp * ( apb_matrix(t_ij,t_ij) + amb_matrix(t_ij,t_ij) )
 enddo

 allocate(cc_matrix(nmat,nmat))
 allocate(amb_eigval(nmat))

 WRITE_MASTER(*,*) 'Diago 2-particle hamiltonian a la Casida'

 !
 ! Calculate (A-B)^{1/2}
 ! First diagonalize (A-B):
 ! (A-B) R = R D
 ! (A-B) is real symmetric, hence R is orthogonal R^{-1} = tR
 ! (A-B)       = R D tR 
 ! (A-B)^{1/2} = R D^{1/2} tR 
 call start_clock(timing_diago_h2p)
 WRITE_MASTER(*,'(a,i8,a,i8)') ' Diago to get (A - B)^{1/2}                   ',nmat,' x ',nmat
 call diagonalize(nmat,amb_matrix,amb_eigval)
 call stop_clock(timing_diago_h2p)

 allocate(amb_matrix_sqrt(nmat,nmat),amb_matrix_sqrtm1(nmat,nmat))
 do t_kl=1,nmat
   amb_matrix_sqrt  (:,t_kl) = amb_matrix(:,t_kl)*SQRT(amb_eigval(t_kl))
   amb_matrix_sqrtm1(:,t_kl) = amb_matrix(:,t_kl)/SQRT(amb_eigval(t_kl))
 enddo
 deallocate(amb_eigval)

 amb_matrix = TRANSPOSE( amb_matrix )
 amb_matrix_sqrt  (:,:) = MATMUL( amb_matrix_sqrt(:,:)   , amb_matrix(:,:) )
 amb_matrix_sqrtm1(:,:) = MATMUL( amb_matrix_sqrtm1(:,:) , amb_matrix(:,:) )
 
 cc_matrix(:,:) = MATMUL( amb_matrix_sqrt , MATMUL( apb_matrix , amb_matrix_sqrt)  )

! WRITE_MASTER(*,*) 'CC ',matrix_is_symmetric(nmat,cc_matrix)

 deallocate(apb_matrix,amb_matrix)

 allocate(bigomega(nmat))
 call start_clock(timing_diago_h2p)
 WRITE_MASTER(*,'(a,i8,a,i8)') ' Diago (A - B)^{1/2} * (A + B) * (A - B)^{1/2}',nmat,' x ',nmat
 call diagonalize(nmat,cc_matrix,bigomega)
 call stop_clock(timing_diago_h2p)

 bigomega(:) = SQRT(bigomega(:))

 rpa_correlation = rpa_correlation + 0.50_dp * SUM( bigomega(:) )

 allocate(bigx(nmat,nmat),bigy(nmat,nmat))
 allocate(cc_matrix_bigomega(nmat,nmat))

 do t_kl=1,nmat
   cc_matrix_bigomega(:,t_kl) = cc_matrix(:,t_kl) * bigomega(t_kl)
   ! Resonant
   eigenvalue(t_kl)      =  bigomega(t_kl)
   ! AntiResonant
   eigenvalue(t_kl+nmat) = -bigomega(t_kl)
 enddo
 bigx(:,:) = 0.5_dp * MATMUL( amb_matrix_sqrt(:,:)   , cc_matrix(:,:) )
 bigy(:,:) = 0.5_dp * MATMUL( amb_matrix_sqrtm1(:,:) , cc_matrix_bigomega(:,:) )

 ! Resonant
 eigenvector(1:nmat       ,1:nmat)        = bigx(:,:)+bigy(:,:)
 eigenvector(nmat+1:2*nmat,1:nmat)        = bigx(:,:)-bigy(:,:)
 ! AntiResonant
 eigenvector(1:nmat       ,nmat+1:2*nmat) = bigx(:,:)-bigy(:,:)
 eigenvector(nmat+1:2*nmat,nmat+1:2*nmat) = bigx(:,:)+bigy(:,:)

 deallocate(cc_matrix_bigomega)


 deallocate(amb_matrix_sqrt,amb_matrix_sqrtm1,bigomega)
 deallocate(bigx,bigy)
 deallocate(cc_matrix)



end subroutine build_h2p_sym


!=========================================================================
end module m_gw
!=========================================================================
