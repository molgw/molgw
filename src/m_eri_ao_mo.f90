!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to perform the Atomic Orbital to Molecular Orbital transform
!
!=========================================================================
module m_eri_ao_mo
 use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
 use m_definitions
 use m_mpi
 use m_memory
 use m_basis_set
 use m_timing
 use m_eri


 real(dp),protected,allocatable :: eri_3center_eigen(:,:,:,:)
 real(dp),protected,allocatable :: eri_3center_eigen_lr(:,:,:,:)
 real(dp),protected,allocatable :: eri_3center_eigen_mixed(:,:,:,:)


contains


!=========================================================================
function eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen_ri
!=====

 eri_eigen_ri = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

 call xsum_auxil(eri_eigen_ri)

end function eri_eigen_ri


!=========================================================================
function eri_eigen_ri_paral(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen_ri_paral
!=====

 eri_eigen_ri_paral = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

end function eri_eigen_ri_paral


!=================================================================
subroutine calculate_eri_4center_eigen(nbf,nstate,c_matrix,istate,ijspin,eri_eigenstate_i)
 use m_inputparam,only: nspin
 implicit none

 integer,intent(in)     :: nbf,nstate
 integer,intent(in)     :: istate,ijspin
 real(dp),intent(in)    :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(inout) :: eri_eigenstate_i(nstate,nstate,nstate,nspin)
!=====
 integer,save         :: istate_previous=0
 integer,save         :: ijspin_previous=0
 integer              :: klspin
 integer              :: ibf,jbf,kbf,lbf
 integer              :: jstate,kstate,lstate
 real(dp)             :: eri_tmp3(nbf,nbf,nbf),eri_tmp2(nbf,nbf,nbf)
!=====

 ! Check if the calculation can be skipped
 if( istate_previous == istate .AND. ijspin_previous == ijspin .AND. ANY(ABS(eri_eigenstate_i(:,:,:,:))>1.0e-6_dp) ) then
   return
 else
   istate_previous = istate
   ijspin_previous = ijspin
 endif


 call start_clock(timing_basis_transform)

 eri_eigenstate_i(:,:,:,:)=0.0_dp
 eri_tmp2(:,:,:)=0.0_dp
 eri_tmp3(:,:,:)=0.0_dp

 do lbf=1,nbf
   do kbf=1,nbf
     do jbf=1,nbf

       do ibf=1,nbf
         eri_tmp3(jbf,kbf,lbf) = eri_tmp3(jbf,kbf,lbf) + eri(ibf,jbf,kbf,lbf) * c_matrix(ibf,istate,ijspin) 
       enddo


     enddo
   enddo
 enddo

 do lbf=1,nbf
   do kbf=1,nbf

     do jstate=1,nstate
       eri_tmp2(jstate,kbf,lbf) = DOT_PRODUCT( eri_tmp3(:,kbf,lbf) , c_matrix(:,jstate,ijspin) )
     enddo

   enddo
 enddo


  
 do klspin=1,nspin

   do lbf=1,nbf
     do kstate=1,nstate
       do jstate=1,nstate
         eri_tmp3(jstate,kstate,lbf) = DOT_PRODUCT( eri_tmp2(jstate,:,lbf) , c_matrix(:,kstate,klspin) )
       enddo
     enddo
   enddo

   do lstate=1,nstate
     do kstate=1,nstate
       do jstate=1,nstate

         eri_eigenstate_i(jstate,kstate,lstate,klspin) = DOT_PRODUCT( eri_tmp3(jstate,kstate,:) , c_matrix(:,lstate,klspin) )

       enddo
     enddo
   enddo

 enddo !klspin


 call stop_clock(timing_basis_transform)

end subroutine calculate_eri_4center_eigen


!=================================================================
subroutine calculate_eri_3center_eigen(nbf,nstate,c_matrix,mstate_min,mstate_max,nstate_min,nstate_max)
 use m_inputparam,only: nspin
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
 integer,intent(in)   :: mstate_min,mstate_max,nstate_min,nstate_max
!=====
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp_l(:,:)
 integer              :: ipair
!=====

 call start_clock(timing_eri_3center_eigen)

 write(stdout,'(/,a)') ' Calculate 3-center integrals on eigenstates'


 !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
 call clean_allocate('3-center MO integrals',eri_3center_eigen,1,nauxil_3center,mstate_min,mstate_max,nstate_min,nstate_max,1,nspin)
 eri_3center_eigen(:,:,:,:) = 0.0_dp

 call clean_allocate('TMP 3-center ints',eri_3center_tmp_l,nauxil_3center,nbf)

 do klspin=1,nspin

   do lstate=nstate_min,nstate_max
     if( MODULO( lstate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     eri_3center_tmp_l(:,:) = 0.0_dp

     ! Transformation of the first index
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       eri_3center_tmp_l(:,kbf) = eri_3center_tmp_l(:,kbf) &
                                       + c_matrix(lbf,lstate,klspin) * eri_3center(:,ipair)
       if( kbf /= lbf ) &
       eri_3center_tmp_l(:,lbf) = eri_3center_tmp_l(:,lbf) &
                                       + c_matrix(kbf,lstate,klspin) * eri_3center(:,ipair)
     enddo


     ! Transformation of the second index
     eri_3center_eigen(:,mstate_min:mstate_max,lstate,klspin) = MATMUL( eri_3center_tmp_l(:,:) , c_matrix(:,mstate_min:mstate_max,klspin) )

   enddo

 enddo ! klspin

 call clean_deallocate('TMP 3-center ints',eri_3center_tmp_l)

 call xsum_ortho(eri_3center_eigen)

 call stop_clock(timing_eri_3center_eigen)

end subroutine calculate_eri_3center_eigen


!=================================================================
subroutine calculate_eri_3center_eigen_lr(nbf,nstate,c_matrix)
 use m_inputparam,only: nspin
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
!=====
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp_l(:,:)
 integer              :: ipair
!=====

 call start_clock(timing_eri_3center_eigen)

 write(stdout,'(/,a)') ' Calculate LR 3-center integrals on eigenstates'


 !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
 call clean_allocate('LR 3-center MO integrals',eri_3center_eigen_lr,nauxil_3center_lr,nstate,nstate,nspin)
 eri_3center_eigen_lr(:,:,:,:) = 0.0_dp

 allocate(eri_3center_tmp_l(nauxil_3center_lr,nbf))

 do klspin=1,nspin

   do lstate=1,nstate
     if( MODULO( lstate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     eri_3center_tmp_l(:,:) = 0.0_dp

     ! Transformation of the first index
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       eri_3center_tmp_l(:,kbf) = eri_3center_tmp_l(:,kbf) &
                                       + c_matrix(lbf,lstate,klspin) * eri_3center_lr(:,ipair)
       if( kbf /= lbf )  &
         eri_3center_tmp_l(:,lbf) = eri_3center_tmp_l(:,lbf) &
                                         + c_matrix(kbf,lstate,klspin) * eri_3center_lr(:,ipair)

     enddo

   ! Transformation of the second index
     eri_3center_eigen_lr(:,:,lstate,klspin) = MATMUL( eri_3center_tmp_l(:,:) , c_matrix(:,:,klspin) )

   enddo

 enddo ! klspin
 deallocate(eri_3center_tmp_l)

 call xsum_ortho(eri_3center_eigen)

 call stop_clock(timing_eri_3center_eigen)

end subroutine calculate_eri_3center_eigen_lr


!=================================================================
subroutine calculate_eri_3center_eigen_mixed(nbf,nstate,c_matrix)
 use m_inputparam,only: nspin
 implicit none

 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
!=====
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp(:,:,:)
 real(dp),allocatable :: c_matrix_exx(:,:,:)
 logical              :: file_exists
!=====

 call start_clock(timing_eri_3center_eigen)

 inquire(file='fort.1000',exist=file_exists)
 if( .NOT. file_exists ) call die('fort.1000 not found')

 allocate(c_matrix_exx(nbf,nstate,nspin))
 open(1000,form='unformatted')
 do klspin=1,nspin
   do lstate=1,nstate
     read(1000) c_matrix_exx(:,lstate,klspin)
   enddo
 enddo
 close(1000,status='delete')


 write(stdout,'(/,a)') ' Calculate 3-center integrals on MIXED eigenstates'


 !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
 call clean_allocate('3-center mixed MO integrals',eri_3center_eigen_mixed,nauxil_3center,nstate,nstate,nspin)
 eri_3center_eigen_mixed(:,:,:,:) = 0.0_dp

 allocate(eri_3center_tmp(nauxil_3center,nbf,nstate)) 

 !TODO fix all this mess here to make it more similar to the previous subroutine
 do klspin=1,nspin
   ! Transformation of the first index
   eri_3center_tmp(:,:,:) = 0.0_dp
   do kbf=1,nbf
     do lbf=1,nbf
       if( negligible_basispair(kbf,lbf) ) cycle

         do lstate=1,nstate
           eri_3center_tmp(:,kbf,lstate) = eri_3center_tmp(:,kbf,lstate) &
                                      + c_matrix_exx(lbf,lstate,klspin) * eri_3center(:,index_pair(kbf,lbf))
         enddo

     enddo
   enddo
   ! Transformation of the second index
   do lstate=1,nstate
     eri_3center_eigen_mixed(:,:,lstate,klspin) = MATMUL( eri_3center_tmp(:,:,lstate) , c_matrix(:,:,klspin) )
   enddo

 enddo ! klspin
 deallocate(eri_3center_tmp)

 call stop_clock(timing_eri_3center_eigen)

end subroutine calculate_eri_3center_eigen_mixed


!=================================================================
subroutine destroy_eri_3center_eigen()
 implicit none
!=====

 write(stdout,'(/,a)') ' Destroy 3-center integrals on eigenstates'
 call clean_deallocate('3-center MO integrals',eri_3center_eigen)

end subroutine destroy_eri_3center_eigen


!=================================================================
subroutine destroy_eri_3center_eigen_lr()
 implicit none
!=====

 write(stdout,'(/,a)') ' Destroy LR 3-center integrals on eigenstates'
 call clean_deallocate('LR 3-center MO integrals',eri_3center_eigen_lr)

end subroutine destroy_eri_3center_eigen_lr


!=================================================================
subroutine destroy_eri_3center_eigen_mixed()
 implicit none
!=====

 write(stdout,'(/,a)') ' Destroy 3-center mixed integrals on eigenstates'
 call clean_deallocate('3-center mixed MO integrals',eri_3center_eigen_mixed)

end subroutine destroy_eri_3center_eigen_mixed


!=========================================================================
end module m_eri_ao_mo
