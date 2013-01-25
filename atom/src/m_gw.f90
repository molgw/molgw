!=========================================================================
#include "macros.h"
!=========================================================================
module m_gw
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_timing 
 use m_warning

 !
 ! General form of any spectral function
 ! z complex number
 ! i, j running on the basis set
 ! sf_ij(z) = \sum_n L_n(i) R_n(j) / ( z - w_n )
 !
 type spectral_function 
   integer :: npole
   integer :: nprodbasis
   real(dp),allocatable :: pole(:)
   real(dp),allocatable :: residu_left(:,:)       ! first index runs on n, second index on i
   real(dp),allocatable :: residu_right(:,:)      ! first index runs on n, second index on j
 end type spectral_function

 !
 ! frozen core approximation parameters
 integer :: ncore_G=0
 integer :: ncore_W=0

 !
 ! frozen virtual approximation parameters
 integer :: nvirtual_G=HUGE(dp)
 integer :: nvirtual_W=HUGE(dp)

#ifdef CRPA
 integer,parameter :: band1=1
 integer,parameter :: band2=2
#endif

contains

!=========================================================================
subroutine init_spectral_function(nbf,prod_nbf,nspin,occupation,sf)
 implicit none
 integer,intent(in)                    :: nbf,prod_nbf,nspin
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 type(spectral_function),intent(inout) :: sf
!=====
 integer                               :: ispin,ibf,jbf
 logical                               :: file_exists
!=====

 !
 ! Deal with frozen core initialization
 inquire(file='manual_frozencore',exist=file_exists)
 if(file_exists) then
   !
   ! ncore_G and ncore_W contain the highest core state to be discarded
   open(13,file='manual_frozencore')
   read(13,*) ncore_G
   read(13,*) ncore_W
   close(13)
   ncore_G = MAX(ncore_G,0)
   ncore_W = MAX(ncore_W,0)
   WRITE_MASTER(msg,'(a,i4,2x,i4)') 'frozen core approximation switched on up to state (G,W) = ',ncore_G,ncore_W
   call issue_warning(msg)
 endif
 !
 ! Deal with frozen virtual initialization
 inquire(file='manual_frozenvirtual',exist=file_exists)
 if(file_exists) then
   !
   ! nvirtual_G and nvirtual_W contain the first virtual state to be discarded
   open(13,file='manual_frozenvirtual')
   read(13,*) nvirtual_G
   read(13,*) nvirtual_W
   close(13)
   nvirtual_G = MAX(nvirtual_G,0)
   nvirtual_W = MAX(nvirtual_W,0)
   WRITE_MASTER(msg,'(a,i4,2x,i4)') 'frozen virtual approximation switched on up to state (G,W) = ',nvirtual_G,nvirtual_W
   call issue_warning(msg)
 endif




 sf%npole=0
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=1,nbf
       if( skip_transition(nspin,ibf,jbf,occupation(ibf,ispin),occupation(jbf,ispin)) ) cycle
       sf%npole = sf%npole+1
     enddo
   enddo
 enddo

 WRITE_MASTER(*,*) 
 WRITE_MASTER(*,*) 'spectral function initialized with',sf%npole,'poles'

 allocate(sf%pole(sf%npole))

#ifdef AUXIL_BASIS
 sf%nprodbasis = prod_nbf
#else
 sf%nprodbasis = prod_nbf * nspin
#endif

 allocate(sf%residu_left (sf%npole,sf%nprodbasis))
 allocate(sf%residu_right(sf%npole,sf%nprodbasis))
 WRITE_MASTER(*,*) '           second index size',sf%nprodbasis

 WRITE_MASTER(*,'(a,f14.0)') ' Memory [Mb] ',REAL(SIZE(sf%residu_left(:,:)),dp)*2.0_dp/1024.0_dp**2*dp
 WRITE_MASTER(*,*)


end subroutine init_spectral_function

!=========================================================================
function skip_transition(nspin,ib1,ib2,occ1,occ2)
 implicit none
 logical             :: skip_transition
 integer,intent(in)  :: nspin,ib1,ib2
 real(dp),intent(in) :: occ1,occ2
!=====
 real(dp)            :: spin_fact
!=====
 spin_fact = REAL(-nspin+3,dp)

 skip_transition = .FALSE.

 !
 ! skip the core states if asked for a frozen-core calculation
 if( ib1 <= ncore_W .OR. ib2 <= ncore_W ) skip_transition=.TRUE.
 !
 ! skip the virtual states if asked for a frozen-virtual calculation
 if( ib1 >= nvirtual_W .OR. ib2 >= nvirtual_W ) skip_transition=.TRUE.

#ifndef CASIDA
 if( occ1 < completely_empty             .AND. occ2 < completely_empty )             skip_transition=.TRUE.
 if( occ1 > spin_fact - completely_empty .AND. occ2 > spin_fact - completely_empty ) skip_transition=.TRUE.
#else
 ! Casida case only positive transition are retained
 if( occ1 > spin_fact - completely_empty .OR.  occ2 < completely_empty )             skip_transition=.TRUE.
#endif

#ifdef CRPA
 if( ib1==band1 .AND. ib2==band1 )  skip_transition=.TRUE.
 if( ib1==band1 .AND. ib2==band2 )  skip_transition=.TRUE.
 if( ib1==band2 .AND. ib2==band1 )  skip_transition=.TRUE.
 if( ib1==band2 .AND. ib2==band2 )  skip_transition=.TRUE.
#endif

end function skip_transition

!=========================================================================
subroutine destroy_spectral_function(sf)
 implicit none
 type(spectral_function),intent(inout) :: sf
!=====

 deallocate(sf%pole)
 deallocate(sf%residu_left)
 deallocate(sf%residu_right)

 WRITE_MASTER(*,*) 
 WRITE_MASTER(*,*) 'spectral function destroyed'

end subroutine destroy_spectral_function


!=========================================================================
subroutine write_spectral_function(sf)
 implicit none
 type(spectral_function),intent(in) :: sf
!=====
 integer,parameter :: spectralfile=50
 integer :: iprodbasis
!=====

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'dumping spectral function on file' 

 open(spectralfile,file='spectral_file',form='unformatted')

 WRITE_MASTER(spectralfile) sf%npole
 WRITE_MASTER(spectralfile) sf%nprodbasis
 WRITE_MASTER(spectralfile) sf%pole(:)
 do iprodbasis=1,sf%nprodbasis
   WRITE_MASTER(spectralfile) sf%residu_left(:,iprodbasis)
 enddo
 do iprodbasis=1,sf%nprodbasis
   WRITE_MASTER(spectralfile) sf%residu_right(:,iprodbasis)
 enddo

 close(spectralfile)

end subroutine write_spectral_function


!=========================================================================
subroutine read_spectral_function(sf)
 implicit none
 type(spectral_function),intent(inout) :: sf
!=====
 integer,parameter :: spectralfile=50
 integer :: iprodbasis
!=====

 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'reading spectral function from file' 

 open(spectralfile,file='spectral_file',status='old',form='unformatted')

 read(spectralfile) sf%npole
 read(spectralfile) sf%nprodbasis
 read(spectralfile) sf%pole(:)
 do iprodbasis=1,sf%nprodbasis
   read(spectralfile) sf%residu_left(:,iprodbasis)
 enddo
 do iprodbasis=1,sf%nprodbasis
   read(spectralfile) sf%residu_right(:,iprodbasis)   
 enddo

 close(spectralfile)

 WRITE_MASTER(*,*) 'Reading done'

end subroutine read_spectral_function


#ifdef AUXIL_BASIS
!=========================================================================
subroutine polarizability_rpa(nspin,basis,prod_basis,occupation,energy,c_matrix,sinv_v_sinv,rpa_correlation,wpol,w_pol)
 use m_tools
 use m_basis_set
 implicit none

 integer,intent(in)  :: nspin
 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin)
 real(dp),intent(in) :: energy(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in) :: sinv_v_sinv(prod_basis%nbf,prod_basis%nbf)
 real(dp),intent(out) :: rpa_correlation
 type(spectral_function),intent(inout) :: wpol
 complex(dp),optional :: w_pol(prod_basis%nbf_filtered,prod_basis%nbf_filtered,NOMEGA)
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

end subroutine polarizability_rpa
#endif

#ifndef AUXIL_BASIS
!=========================================================================
subroutine polarizability_rpa_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_tools
 use m_basis_set
 use m_eri
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

#ifndef LOW_MEMORY2
 real(dp) :: eri_eigenstate(basis%nbf,basis%nbf,basis%nbf,basis%nbf,nspin,nspin)
#else
 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp),allocatable :: eri_eigenstate_k(:,:,:,:)
#endif
 real(dp) :: spin_fact 
 real(dp) :: h_2p(wpol%npole,wpol%npole)
 real(dp) :: eigenvalue(wpol%npole),eigenvector(wpol%npole,wpol%npole),eigenvector_inv(wpol%npole,wpol%npole)
 real(dp) :: matrix(wpol%npole,wpol%npole)
 real(dp) :: rtmp

 logical :: TDHF=.FALSE.
!=====
 spin_fact = REAL(-nspin+3,dp)
 rpa_correlation = 0.0_dp

 WRITE_MASTER(*,'(/,a)') ' calculating CHI alla rpa'
 if(TDHF) then
   msg='calculating the TDHF polarizability'
   call issue_warning(msg)
 endif

#ifndef LOW_MEMORY2
 call transform_eri_basis_fast(basis%nbf,nspin,c_matrix,eri_eigenstate)
#else
 allocate(eri_eigenstate_i(basis%nbf,basis%nbf,basis%nbf,nspin))
#endif


 t_ij=0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

#ifdef LOW_MEMORY2
       call transform_eri_basis_lowmem(nspin,c_matrix,iorbital,ijspin,eri_eigenstate_i)
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
       if( skip_transition(nspin,iorbital,jorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
       t_ij=t_ij+1



       t_kl=0
       do klspin=1,nspin
         do korbital=1,basis%nbf 

           do lorbital=1,basis%nbf 
             if( skip_transition(nspin,korbital,lorbital,occupation(lorbital,klspin),occupation(korbital,klspin)) ) cycle
             t_kl=t_kl+1

#ifndef CHI0
#ifndef LOW_MEMORY2
             h_2p(t_ij,t_kl) = eri_eigenstate(iorbital,jorbital,korbital,lorbital,ijspin,klspin) &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )
#else
             h_2p(t_ij,t_kl) = eri_eigenstate_i(jorbital,korbital,lorbital,klspin) &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )

!!Testing the formula for partial occupations of WT Yang.
!             h_2p(t_ij,t_kl) = eri_eigenstate_i(jorbital,korbital,lorbital,klspin) &
!                        * SIGN( 1.0_dp , occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) &
!                        * SQRT(    MAX( occupation(iorbital,ijspin) , occupation(jorbital,ijspin) )   &
!                               *   MAX( occupation(korbital,klspin) , occupation(lorbital,klspin) )   &
!                               * ( spin_fact - MIN( occupation(iorbital,ijspin) , occupation(jorbital,ijspin) )  )&
!                               * ( spin_fact - MIN( occupation(korbital,klspin) , occupation(lorbital,klspin) )  )  )
!!Testing my formula
!             h_2p(t_ij,t_kl) = eri_eigenstate_i(jorbital,korbital,lorbital,klspin) &
!                        * SIGN( 1.0_dp , occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) &
!                        *         MAX( occupation(iorbital,ijspin) , occupation(jorbital,ijspin) ) &
!                              * ( spin_fact - MIN( occupation(iorbital,ijspin) , occupation(jorbital,ijspin) )  ) / spin_fact


#endif
#else
             h_2p(t_ij,t_kl) = 0.0_dp
#endif


           if(TDHF) then
             if(ijspin==klspin) then
#ifndef LOW_MEMORY2
               h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigenstate(iorbital,korbital,jorbital,lorbital,ijspin,klspin)  &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact 
#else
               h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigenstate_i(korbital,jorbital,lorbital,klspin)  &
                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact 
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
!       if( skip_transition(nspin,iorbital,jorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
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

! do t_ij=1,wpol%npole
!   WRITE_MASTER(*,'(1(i4,2x),20(2x,f12.6))') t_ij,eigenvalue(t_ij)
! enddo
   
 call start_clock(timing_inversion_s2p)
 call invert(wpol%npole,eigenvector,eigenvector_inv)
 call stop_clock(timing_inversion_s2p)

#ifdef LOW_MEMORY2
 deallocate(eri_eigenstate_i)
 allocate(eri_eigenstate_k(basis%nbf,basis%nbf,basis%nbf,nspin))
#endif

 wpol%pole = eigenvalue
 wpol%residu_left (:,:) = 0.0_dp
 wpol%residu_right(:,:) = 0.0_dp
 t_kl=0
 do klspin=1,nspin
   do kbf=1,basis%nbf 
#ifdef LOW_MEMORY2
     call transform_eri_basis_lowmem(nspin,c_matrix,kbf,klspin,eri_eigenstate_k)
#endif


     do lbf=1,basis%nbf
       if( skip_transition(nspin,kbf,lbf,occupation(lbf,klspin),occupation(kbf,klspin)) ) cycle
       t_kl=t_kl+1


       do ijspin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(ibf,jbf,ijbf_current)
         do ijbf=1,prod_basis%nbf
           ibf = prod_basis%index_ij(1,ijbf)
           jbf = prod_basis%index_ij(2,ijbf)

           ijbf_current = ijbf+prod_basis%nbf*(ijspin-1)


#ifndef LOW_MEMORY2
           wpol%residu_left (:,ijbf_current)  = wpol%residu_left (:,ijbf_current) &
                        + eri_eigenstate(ibf,jbf,kbf,lbf,ijspin,klspin) *  eigenvector(t_kl,:)
           wpol%residu_right(:,ijbf_current)  = wpol%residu_right(:,ijbf_current) &
                        + eri_eigenstate(kbf,lbf,ibf,jbf,klspin,ijspin) * eigenvector_inv(:,t_kl) &
                                         * ( occupation(kbf,klspin)-occupation(lbf,klspin) )
#else
           wpol%residu_left (:,ijbf_current)  = wpol%residu_left (:,ijbf_current) &
                        + eri_eigenstate_k(lbf,ibf,jbf,ijspin) *  eigenvector(t_kl,:)
           wpol%residu_right(:,ijbf_current)  = wpol%residu_right(:,ijbf_current) &
                        + eri_eigenstate_k(lbf,ibf,jbf,ijspin) * eigenvector_inv(:,t_kl) &
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
 


#ifdef LOW_MEMORY2
 deallocate(eri_eigenstate_k)
#endif

end subroutine polarizability_rpa_noaux

!=========================================================================
subroutine polarizability_casida_noaux(nspin,basis,prod_basis,occupation,energy,c_matrix,rpa_correlation,wpol)
 use m_tools
 use m_basis_set
 use m_eri
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

 real(dp),allocatable :: eri_eigenstate_i(:,:,:,:)
 real(dp),allocatable :: eri_eigenstate_k(:,:,:,:)

 real(dp) :: spin_fact 
 real(dp) :: apb(wpol%npole,wpol%npole)
 real(dp) :: amb_diag(wpol%npole)
 real(dp) :: eigenvalue(wpol%npole),eigenvector(wpol%npole,wpol%npole),eigenvector_inv(wpol%npole,wpol%npole)
 real(dp) :: matrix(wpol%npole,wpol%npole)

 logical :: TDHF=.FALSE.
!=====
 spin_fact = REAL(-nspin+3,dp)
 rpa_correlation = 0.0_dp

#ifndef LOW_MEMORY2
 stop 'NOT implemented'
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

 apb(:,:) = 0.0_dp

 t_ij=0
 do ijspin=1,nspin
   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied

     !
     ! TODO: Add check if it is OK like this 
     if( occupation(iorbital,ijspin) < completely_empty ) cycle
     call transform_eri_basis_lowmem(nspin,c_matrix,iorbital,ijspin,eri_eigenstate_i)


     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty

       if( skip_transition(nspin,iorbital,jorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
       t_ij=t_ij+1

       amb_diag(t_ij) = energy(jorbital,ijspin) - energy(iorbital,ijspin)
       apb(t_ij,t_ij) = amb_diag(t_ij)

       t_kl=0
       do klspin=1,nspin
         do korbital=1,basis%nbf 

           do lorbital=1,basis%nbf 
             if( skip_transition(nspin,korbital,lorbital,occupation(lorbital,klspin),occupation(korbital,klspin)) ) cycle
             t_kl=t_kl+1

             apb(t_ij,t_kl) = apb(t_ij,t_kl) + 2.0_dp * eri_eigenstate_i(jorbital,korbital,lorbital,klspin) &
                                                      * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) )


!TODO check TDHF implementation
!           if(TDHF) then
!             if(ijspin==klspin) then
!               h_2p(t_ij,t_kl) =  h_2p(t_ij,t_kl) -  eri_eigenstate_i(korbital,jorbital,lorbital,klspin)  &
!                        * ( occupation(iorbital,ijspin)-occupation(jorbital,ijspin) ) / spin_fact 
!             endif
!           endif

           enddo
         enddo
       enddo ! klspin



       rpa_correlation = rpa_correlation - 0.25_dp * ( amb_diag(t_ij) + apb(t_ij,t_ij) )

     enddo !jorbital
   enddo !iorbital
 enddo ! ijspin

 do t_kl=1,wpol%npole
   do t_ij=1,wpol%npole
     matrix(t_ij,t_kl) = SQRT( amb_diag(t_ij) ) * apb(t_ij,t_kl) * SQRT( amb_diag(t_kl) )
   enddo
 enddo

 WRITE_MASTER(*,*) 'diago Casida matrix'
 WRITE_MASTER(*,*) 'matrix',wpol%npole,'x',wpol%npole
! WRITE_MASTER(*,*) '=============='
! t_ij=0
! do ijspin=1,nspin
!   do iorbital=1,basis%nbf ! iorbital stands for occupied or partially occupied
!     do jorbital=1,basis%nbf ! jorbital stands for empty or partially empty
!!INTRA       if(iorbital==jorbital) cycle  ! intra state transitions are not allowed!
!!SKIP       if( abs(occupation(jorbital,ijspin)-occupation(iorbital,ijspin))<completely_empty ) cycle
!       if( skip_transition(nspin,iorbital,jorbital,occupation(jorbital,ijspin),occupation(iorbital,ijspin)) ) cycle
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

 !
 ! Symmetric in-place diagonalization
 call start_clock(timing_diago_h2p)
! call diagonalize_general(wpol%npole,matrix,eigenvalue,eigenvector)
 call diagonalize(wpol%npole,matrix,eigenvalue)
 call stop_clock(timing_diago_h2p)
 WRITE_MASTER(*,*) 'diago finished'
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) 'calculate the RPA energy using the Tamm-Dancoff decomposition'
 WRITE_MASTER(*,*) 'formula (9) from J. Chem. Phys. 132, 234114 (2010)'
 rpa_correlation = rpa_correlation + 0.50_dp * SUM( SQRT(eigenvalue(:)) )
 WRITE_MASTER(*,'(/,a,f14.8)') ' RPA energy [Ha]: ',rpa_correlation


 deallocate(eri_eigenstate_i)
!§  stop'ENOUGH'
!§  allocate(eri_eigenstate_k(basis%nbf,basis%nbf,basis%nbf,nspin))
!§ 
!§  wpol%pole(:)           = SQRT(eigenvalue(:))
!§  wpol%residu_left (:,:) = 0.0_dp
!§  wpol%residu_right(:,:) = 0.0_dp
!§  t_kl=0
!§  do klspin=1,nspin
!§    do kbf=1,basis%nbf 
!§ 
!§      call transform_eri_basis_lowmem(nspin,c_matrix,kbf,klspin,eri_eigenstate_k)
!§ 
!§      do lbf=1,basis%nbf
!§        if( skip_transition(nspin,kbf,lbf,occupation(lbf,klspin),occupation(kbf,klspin)) ) cycle
!§        t_kl=t_kl+1
!§ 
!§ 
!§        do ijspin=1,nspin
!§ !$OMP PARALLEL DEFAULT(SHARED)
!§ !$OMP DO PRIVATE(ibf,jbf,ijbf_current)
!§          do ijbf=1,prod_basis%nbf
!§            ibf = prod_basis%index_ij(1,ijbf)
!§            jbf = prod_basis%index_ij(2,ijbf)
!§ 
!§            ijbf_current = ijbf+prod_basis%nbf*(ijspin-1)
!§ 
!§ 
!§            wpol%residu_left (:,ijbf_current)  = wpol%residu_left (:,ijbf_current) &
!§                         + eri_eigenstate_k(lbf,ibf,jbf,ijspin) *  eigenvector(t_kl,:)
!§            wpol%residu_right(:,ijbf_current)  = wpol%residu_right(:,ijbf_current) &
!§                         + eri_eigenstate_k(lbf,ibf,jbf,ijspin) * eigenvector_inv(:,t_kl) &
!§                                          * ( occupation(kbf,klspin)-occupation(lbf,klspin) )
!§ 
!§ 
!§          enddo
!§ !$OMP END DO
!§ !$OMP END PARALLEL
!§        enddo
!§      enddo
!§    enddo
!§  enddo
!§ 
!§  deallocate(eri_eigenstate_k)

end subroutine polarizability_casida_noaux
#endif


#ifdef AUXIL_BASIS
!=========================================================================
subroutine gw_selfenergy(method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_basis_set
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
 complex(dp) :: eta=(0.0_dp,0.0001_dp)
! complex(dp) :: eta=(0.0_dp,1.50_dp)
! complex(dp) :: eta=(0.0_dp,0.0_dp)
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
           if( ABS( energy(borbital,ispin) - energy(iorbital,ispin) + wpol%pole(ipole) ) < aimag(eta) ) then
             energy_complex = energy(borbital,ispin)  + (0.0_dp,1.0_dp) &
&                * ( aimag(eta) - ABS( energy(borbital,ispin) - energy(iorbital,ispin) + wpol%pole(ipole) ) )
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
               if( ABS( energy(borbital,ispin) + omegai(iomegai) - energy(iorbital,ispin) + wpol%pole(ipole) ) < aimag(eta) ) then
!                 WRITE_MASTER(*,*) 'avoid pole'
                 energy_complex = energy(borbital,ispin)  + (0.0_dp,1.0_dp) &
&                    * ( aimag(eta) - ABS( energy(borbital,ispin) + omegai(iomegai) - energy(iorbital,ispin) + wpol%pole(ipole) ) )
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

end subroutine gw_selfenergy
#endif

#ifndef AUXIL_BASIS
!=========================================================================
subroutine gw_selfenergy_noaux(method,nspin,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_basis_set
 implicit none

 integer,intent(in)  :: method,nspin
 type(basis_set)     :: basis,prod_basis
 real(dp),intent(in) :: occupation(basis%nbf,nspin),energy(basis%nbf,nspin),exchange_m_vxc_diag(basis%nbf,nspin)
 real(dp),intent(in) :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in) :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out) :: selfenergy(basis%nbf,basis%nbf,nspin)

 logical               :: file_exists=.FALSE.
 integer               :: nomegai
 integer               :: iomegai
 complex(dp),parameter :: ieta=(0.0_dp,0.0001_dp) ! (0.0_dp,0.0001_dp)
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
 character(len=3) :: ctmp
!=====
 spin_fact = REAL(-nspin+3,dp)

 write(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)

 WRITE_MASTER(*,*)
 select case(method)
 case(QS)
   WRITE_MASTER(*,*) 'perform a QP self-consistent GW calculation'
 case(perturbative)
   WRITE_MASTER(*,*) 'perform a one-shot G0W0 calculation'
 end select

 if(method==QS) then
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

     bra(:,:)=0.0_dp
     ket(:,:)=0.0_dp
     !
     !TODO the search could be improved
     do kbf=1,prod_basis%nbf
       if( prod_basis%index_ij(1,kbf)==iorbital ) then
         aorbital = prod_basis%index_ij(2,kbf)
         bra(:,aorbital) = wpol%residu_left (:,kbf+prod_basis%nbf*(ispin-1))
         ket(:,aorbital) = wpol%residu_right(:,kbf+prod_basis%nbf*(ispin-1))
       else if( prod_basis%index_ij(2,kbf)==iorbital ) then
         aorbital = prod_basis%index_ij(1,kbf)
         bra(:,aorbital) = wpol%residu_left (:,kbf+prod_basis%nbf*(ispin-1))
         ket(:,aorbital) = wpol%residu_right(:,kbf+prod_basis%nbf*(ispin-1))
       endif
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

       if(method==QS) then

         do borbital=1,basis%nbf
           do aorbital=1,basis%nbf

             selfenergy_tmp(1,aorbital,borbital,ispin) = selfenergy_tmp(1,aorbital,borbital,ispin) &
                        - bra(ipole,aorbital) * ket(ipole,borbital) &
                          * REAL(   fact_empty / ( energy(borbital,ispin) + ieta  - energy(iorbital,ispin) + wpol%pole(ipole) ) &
                                   -fact_full  / ( energy(borbital,ispin) - ieta  - energy(iorbital,ispin) + wpol%pole(ipole) ) , dp )

           enddo
         enddo

       else if(method==perturbative) then

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
     selfenergy_tmp(1,:,:,ispin) = 0.5_dp * ( selfenergy_tmp(1,:,:,ispin) + transpose(selfenergy_tmp(1,:,:,ispin)) )
   enddo

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
     do aorbital=1,basis%nbf ! MIN(2,basis%nbf)
       write(ctmp,'(i3.3)') aorbital
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
   WRITE_MASTER(*,*) 'G0W0 Eigenvalues [Ha]'
   if(nspin==1) then
     WRITE_MASTER(*,*) '  #          E0        Sigx-Vxc      Sigc          Z         G0W0'
   else
     WRITE_MASTER(*,'(a)') '  #                E0                      Sigx-Vxc                    Sigc                       Z                       G0W0'
   endif
   do aorbital=1,basis%nbf
     zz(:) = REAL( selfenergy_tmp(3,aorbital,aorbital,:) - selfenergy_tmp(1,aorbital,aorbital,:) ) / REAL( omegai(3)-omegai(1) )
     zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )

     WRITE_MASTER(*,'(i4,x,20(x,f12.6))') aorbital,energy(aorbital,:),exchange_m_vxc_diag(aorbital,:),REAL(selfenergy_tmp(2,aorbital,aorbital,:),dp),&
           zz(:),energy(aorbital,:)+zz(:)*REAL(selfenergy_tmp(2,aorbital,aorbital,:) + exchange_m_vxc_diag(aorbital,:) ,dp)
   enddo

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

     WRITE_MASTER(*,'(i4,x,20(x,f12.6))') aorbital,energy(aorbital,:)*Ha_eV,exchange_m_vxc_diag(aorbital,:)*Ha_eV,REAL(selfenergy_tmp(2,aorbital,aorbital,:),dp)*Ha_eV,&
           zz(:),( energy(aorbital,:)+zz(:)*REAL(selfenergy_tmp(2,aorbital,aorbital,:) + exchange_m_vxc_diag(aorbital,:) ,dp) )*Ha_eV
   enddo

 endif

 deallocate(omegai)
 deallocate(selfenergy_tmp)


end subroutine gw_selfenergy_noaux
#endif


!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================

#if 0
!=========================================================================
subroutine cohsex_selfenergy(nstate,nspin,occupation,energy,c_matrix,w_pol,selfenergy)
 use m_tools,only: coeffs_gausslegint,diagonalize
 implicit none

 integer,intent(in)  :: nstate,nspin
 real(dp),intent(in) :: energy(nstate,nspin),c_matrix(nstate,nstate,nspin)
 complex(dp),intent(in) :: w_pol(nstate_pola,nstate_pola,NOMEGA)
 real(dp),intent(in) :: occupation(nstate,nspin)
 real(dp),intent(out) :: selfenergy(nstate,nstate,nspin)

 integer :: astate,bstate,istate,kstate,lstate
 integer :: iorbital,ispin
 real(dp) :: three_state_overlap_kai,three_state_overlap_ib(nstate_pola)
 complex(dp) :: ket(nstate_pola,nstate)

 real(dp) :: density_matrix(nstate,nstate,nspin)
 real(dp) :: partial(nstate_pola,nstate,nstate)

 real(dp) :: fact
 complex(dp) :: selfenergy_tmp(nstate,nstate,nspin)
!=====



 WRITE_MASTER(*,*) 'building COHSEX self-energy'
 WRITE_MASTER(*,*)


! Poor man implementation of the COHSEX:
! the delta function in the COH is represented by a sum over states...
! and also the rich man implementation

 selfenergy_tmp(:,:,:) = (0.0_dp,0.0_dp)

! SEX
 if(.TRUE.) then

   WRITE_MASTER(*,*) 'calculate first the SEX term'

   density_matrix(:,:,:) = 0.0_dp
   do kstate=1,nstate
     do lstate=1,nstate
       do iorbital=1,nstate
         density_matrix(kstate,lstate,:) = density_matrix(kstate,lstate,:) &
&            + occupation(iorbital,:) * c_matrix(kstate,iorbital,:) * c_matrix(lstate,iorbital,:)
       enddo
     enddo
   enddo

   selfenergy_tmp(:,:,:) = 0.0_dp

   do bstate=1,nstate
     do lstate=1,nstate
       partial(:,lstate,bstate) = matmul( three_state_overlap_product_pola(:,lstate,bstate), w_pol(:,:,1) ) 
     enddo
   enddo

   do bstate=1,nstate
     do astate=1,nstate
       do kstate=1,nstate
         do lstate=1,nstate
           selfenergy_tmp(astate,bstate,:) = selfenergy_tmp(astate,bstate,:) &
&             - dot_product( three_state_overlap_product_pola(:,kstate,astate) ,  partial(:,lstate,bstate) ) &
&                * density_matrix(kstate,lstate,:) 
         enddo
       enddo
     enddo
   enddo

 else
   WRITE_MASTER(*,*) 'WARNING: SEX term skipped!'
 endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! the following calculation of the SEX term is wrong and I don't know why
 if(.FALSE.) then

 do ispin=1,nspin
   do iorbital=1,nstate !INNER LOOP of G

!     fact = 0.5_dp - occupation(iorbital,ispin) / spin_fact
!SEX     WRITE_MASTER(*,*) 'Screened exchange only'

     if( occupation(iorbital,ispin)<completely_empty ) cycle
     fact =  - occupation(iorbital,ispin) / spin_fact

     ket(:,:) = 0.0_dp
     do bstate=1,nstate

       three_state_overlap_ib = matmul( three_state_overlap_product(:,:,bstate) , c_matrix(:,iorbital,ispin) )

!       ket(:,bstate) = ket(:,bstate) + fact * matmul( w_pol(:,:,1) ,
!       three_state_overlap_ib(:) )
       ket(:,bstate) = matmul( w_pol(:,:,1) , three_state_overlap_ib(:) )

     enddo !INNER LOOP of G

     do astate=1,nstate 

       do kstate=1,nstate_pola
         three_state_overlap_kai = 0.0_dp
         do istate=1,nstate
           three_state_overlap_kai = three_state_overlap_kai &
&                + c_matrix(istate,iorbital,ispin) * three_state_overlap_product_pola(kstate,astate,istate)
         enddo
         do bstate=1,nstate 
           selfenergy_tmp(astate,bstate,ispin) = selfenergy_tmp(astate,bstate,ispin) + three_state_overlap_kai*ket(kstate,bstate) * fact
         enddo
       enddo

     enddo

   enddo ! iorbital
 enddo ! ispin

 endif ! WRONG SEX calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! COH
 if(.TRUE.) then
   WRITE_MASTER(*,*) 'then build the COulomb Hole'
   do astate=1,nstate
     do bstate=1,nstate
       do kstate=1,nstate_pola
         do lstate=1,nstate_pola
           selfenergy_tmp(astate,bstate,:) = selfenergy_tmp(astate,bstate,:) &
&            + 0.5_dp * w_pol(kstate,lstate,1) &
&                  * four_wf_product(basis(astate),basis(bstate),basis_pola(kstate),basis_pola(lstate))
         enddo
       enddo
     enddo
   enddo
 else 
   WRITE_MASTER(*,*) 'WARNING: COH term skipped!'
 endif

! selfenergy_tmp(:,:,:,:) = -selfenergy_tmp(:,:,:,:) / ( 2.0_dp * pi )

 selfenergy(:,:,:) = selfenergy_tmp(:,:,:)

! WRITE_MASTER(*,'(a,6(2x,f12.6))') 'COHSEX [Ha]:',selfenergy(1,1,:)
! WRITE_MASTER(*,'(a,6(2x,f12.6))') 'COHSEX [Ha]:',selfenergy(2,2,:)
! WRITE_MASTER(*,'(a,6(2x,f12.6))') 'COHSEX [Ha]:',selfenergy(1,3,:)
! WRITE_MASTER(*,'(a,6(2x,f12.6))') 'COHSEX [Ha]:',selfenergy(3,1,:)
! WRITE_MASTER(*,'(a,6(2x,f12.6))') 'COHSEX [Ha]:',selfenergy(3,3,:)


end subroutine cohsex_selfenergy

!=========================================================================
function screened_exchange_energy(nstate,nspin,density_matrix,w_pol)
 implicit none

 integer,intent(in)  :: nstate,nspin
 real(dp),intent(in) :: density_matrix(nstate,nstate,nspin)
 complex(dp),intent(in) :: w_pol(nstate_pola,nstate_pola,NOMEGA)
 real(dp) :: screened_exchange_energy
!=====
 integer :: kstate,lstate,astate,bstate
 real(dp) :: pkl,ppq
 real(dp) :: partial(nstate_pola,nstate,nstate)
!=====

 screened_exchange_energy = 0.0_dp

 do bstate=1,nstate
   do lstate=1,nstate
     partial(:,lstate,bstate) = matmul( three_state_overlap_product_pola(:,lstate,bstate), w_pol(:,:,1) ) 
   enddo
 enddo

 do astate=1,nstate
   do bstate=1,nstate
     do lstate=1,nstate
       do kstate=1,nstate

         screened_exchange_energy = screened_exchange_energy &
&          - 0.5_dp * dot_product( three_state_overlap_product_pola(:,kstate,astate) ,  partial(:,lstate,bstate) ) &
&            * SUM( density_matrix(astate,bstate,:) * density_matrix(lstate,kstate,:) ) / spin_fact

       enddo
     enddo
   enddo
 enddo

end function screened_exchange_energy


#endif



end module m_gw
!=========================================================================





