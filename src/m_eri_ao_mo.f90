!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to perform the Atomic Orbital to Molecular Orbital transform
!
!=========================================================================
module m_eri_ao_mo
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_scalapack
 use m_memory
 use m_warning
 use m_basis_set
 use m_timing
 use m_eri
 use m_inputparam,only: nspin,has_auxil_basis


 real(dp),protected,allocatable :: eri_3center_eigen(:,:,:,:)
 real(dp),protected,allocatable :: eri_3center_eigen_lr(:,:,:,:)
 real(dp),protected,allocatable :: eri_3center_eigen_mixed(:,:,:,:)

 real(dp),protected,allocatable :: eri_4center_eigen_uks(:,:,:,:)


contains


!=========================================================================
function eri_eigen(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen
!=====

 if(has_auxil_basis) then
   eri_eigen = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )
   call xsum_auxil(eri_eigen)
 else
   eri_eigen = eri_4center_eigen_uks(istate,jstate,kstate,lstate)
 endif


end function eri_eigen


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
pure function eri_eigen_ri_paral(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen_ri_paral
!=====

 eri_eigen_ri_paral = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

end function eri_eigen_ri_paral


!=================================================================
subroutine calculate_eri_4center_eigen(nbf,nstate,c_matrix,istate,ijspin,eri_eigenstate_i)
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


 call start_clock(timing_eri_4center_eigen)

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


 call stop_clock(timing_eri_4center_eigen)

end subroutine calculate_eri_4center_eigen


!=================================================================
subroutine calculate_eri_4center_eigen_uks(c_matrix,nstate_min,nstate_max)
 implicit none

 real(dp),intent(in)    :: c_matrix(:,:,:)
 integer,intent(in)     :: nstate_min,nstate_max
!=====
 integer              :: nbf,nstate
 integer              :: ibf,jbf,kbf,lbf
 integer              :: istate,jstate
 real(dp),allocatable :: eri_tmp3(:,:,:),eri_tmp2(:,:,:),eri_tmp1(:,:)
 integer,allocatable  :: id(:)
!=====

 write(stdout,'(/,1x,a)') 'Calculate all the 4-center MO integrals at once'

 if( nspin /= 1 ) call die('calculate_eri_4center_eigen_uks: requires spin-restricted calculation')

 nbf    = SIZE(c_matrix,DIM=1)
 nstate = SIZE(c_matrix,DIM=2)

 call start_clock(timing_eri_4center_eigen)

 call clean_allocate('4-center MO integrals',eri_4center_eigen_uks, &
                     nstate_min,nstate_max,nstate_min,nstate_max,nstate_min,nstate_max,nstate_min,nstate_max)
 eri_4center_eigen_uks(:,:,:,:) = 0.0_dp

 allocate(eri_tmp3(nbf,nbf,nbf),eri_tmp2(nstate_min:nstate_max,nbf,nbf),eri_tmp1(nstate_min:nstate_max,nbf))

 allocate(id(nbf))
 forall(ibf=1:nbf)
   id(ibf) = ibf
 end forall


 do istate=nstate_min,nstate_max
   if( MODULO( istate - nstate_min , nproc_ortho ) /= rank_ortho ) cycle

   do lbf=1,nbf
     do kbf=1,nbf
       do jbf=1,nbf
         eri_tmp3(jbf,kbf,lbf) = DOT_PRODUCT( c_matrix(1:nbf,istate,1) ,  eri(id(1:nbf),jbf,kbf,lbf) )
       enddo
     enddo
   enddo


   do lbf=1,nbf
     eri_tmp2(nstate_min:nstate_max,1:nbf,lbf) = MATMUL( TRANSPOSE(c_matrix(:,nstate_min:nstate_max,1)) , eri_tmp3(:,1:nbf,lbf) )
   enddo

   do jstate=nstate_min,nstate_max
     eri_tmp1(nstate_min:nstate_max,1:nbf) = MATMUL( TRANSPOSE(c_matrix(:,nstate_min:nstate_max,1)) ,  eri_tmp2(jstate,:,1:nbf) )

     eri_4center_eigen_uks(istate,jstate,nstate_min:nstate_max,nstate_min:nstate_max) =   &
                       MATMUL( eri_tmp1(:,:) , c_matrix(:,nstate_min:nstate_max,1) )
   enddo

 enddo

 call xsum_ortho(eri_4center_eigen_uks)


 deallocate(eri_tmp1,eri_tmp2,eri_tmp3)
 deallocate(id)

 call stop_clock(timing_eri_4center_eigen)

end subroutine calculate_eri_4center_eigen_uks


!=================================================================
subroutine destroy_eri_4center_eigen_uks()
 implicit none
!=====
!=====

 call clean_deallocate('4-center MO integrals',eri_4center_eigen_uks)

end subroutine destroy_eri_4center_eigen_uks


!=================================================================
subroutine calculate_eri_3center_eigen(c_matrix,mstate_min,mstate_max,nstate_min,nstate_max)
 implicit none
 real(dp),intent(in)         :: c_matrix(:,:,:)
 integer,optional,intent(in) :: mstate_min,mstate_max,nstate_min,nstate_max
!=====
 integer              :: nbf,nstate
 integer              :: mstate_min_,mstate_max_,nstate_min_,nstate_max_
 integer              :: mstate_count_,nstate_count_
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp_l(:,:)
 integer              :: ipair
!=====

 call start_clock(timing_eri_3center_eigen)

 write(stdout,'(/,a)') ' Calculate 3-center integrals on eigenstates'

 nbf    = SIZE(c_matrix,DIM=1)
 nstate = SIZE(c_matrix,DIM=2)

 if( PRESENT(mstate_min) ) then
   mstate_min_ = mstate_min
 else
   mstate_min_ = 1
 endif

 if( PRESENT(mstate_max) ) then
   mstate_max_ = mstate_max
 else
   mstate_max_ = nstate
 endif

 if( PRESENT(nstate_min) ) then
   nstate_min_ = nstate_min
 else
   nstate_min_ = 1
 endif

 if( PRESENT(nstate_max) ) then
   nstate_max_ = nstate_max
 else
   nstate_max_ = nstate
 endif

 mstate_count_ = mstate_max_ - mstate_min_ + 1
 nstate_count_ = nstate_max_ - nstate_min_ + 1

 !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
 call clean_allocate('3-center MO integrals',eri_3center_eigen,1,nauxil_3center,mstate_min_,mstate_max_,nstate_min_,nstate_max_,1,nspin)
 eri_3center_eigen(:,:,:,:) = 0.0_dp

 call clean_allocate('TMP 3-center ints',eri_3center_tmp_l,nauxil_3center,nbf)

 do klspin=1,nspin

   do lstate=nstate_min_,nstate_max_
     if( MODULO( lstate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     eri_3center_tmp_l(:,:) = 0.0_dp

     ! Transformation of the first index
     !$OMP PARALLEL PRIVATE(kbf,lbf)
     !$OMP DO REDUCTION(+:eri_3center_tmp_l)
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       eri_3center_tmp_l(:,kbf) = eri_3center_tmp_l(:,kbf) &
                                       + c_matrix(lbf,lstate,klspin) * eri_3center(:,ipair)
       if( kbf /= lbf ) &
       eri_3center_tmp_l(:,lbf) = eri_3center_tmp_l(:,lbf) &
                                       + c_matrix(kbf,lstate,klspin) * eri_3center(:,ipair)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL


     ! Transformation of the second index
     !eri_3center_eigen(:,mstate_min_:mstate_max_,lstate,klspin) = MATMUL( eri_3center_tmp_l(:,:) , c_matrix(:,mstate_min_:mstate_max_,klspin) )

     call DGEMM('N','N',nauxil_3center,mstate_count_,nbf,1.0d0,eri_3center_tmp_l,nauxil_3center,   &
                                                               c_matrix(1,mstate_min_,klspin),nbf, &
                                                         0.0d0,eri_3center_eigen(1,mstate_min_,lstate,klspin),nauxil_3center)

   enddo

 enddo ! klspin

 call clean_deallocate('TMP 3-center ints',eri_3center_tmp_l)

 call xsum_ortho(eri_3center_eigen)

 call stop_clock(timing_eri_3center_eigen)

end subroutine calculate_eri_3center_eigen


!=================================================================
subroutine calculate_eri_3center_eigen_lr(c_matrix)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
!=====
 integer              :: nbf,nstate
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp_l(:,:)
 integer              :: ipair
!=====

 call start_clock(timing_eri_3center_eigen)

 write(stdout,'(/,a)') ' Calculate LR 3-center integrals on eigenstates'

 if( npcol_3center > 1 ) call die('calculate_eri_3center_eigen_lr: incompatible with npcol_3center > 1')

 nbf    = SIZE(c_matrix,DIM=1)
 nstate = SIZE(c_matrix,DIM=2)

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
subroutine calculate_eri_3center_eigen_mixed(c_matrix)
 implicit none

 real(dp),intent(in)  :: c_matrix(:,:,:)
!=====
 integer              :: nbf,nstate
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp(:,:,:)
 real(dp),allocatable :: c_matrix_exx(:,:,:)
 logical              :: file_exists
!=====

 call start_clock(timing_eri_3center_eigen)

 nbf    = SIZE(c_matrix,DIM=1)
 nstate = SIZE(c_matrix,DIM=2)

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
