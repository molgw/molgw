!=========================================================================
subroutine setup_hartree(PRINT_VOLUME,nbf,nspin,p_matrix,pot_hartree,ehartree)
 use m_definitions
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: PRINT_VOLUME
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_hartree(nbf,nbf,nspin)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
 integer              :: ii,ii_tmp
!=====

 call start_clock(timing_hartree)

 pot_hartree(:,:,:)=0.0_dp

!######ifndef LOW_MEMORY
#if 1      
 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ij|kl) = (kl|ij) has been used to loop in the fast order
           pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                      + eri(kbf,lbf,ibf,jbf) * SUM( p_matrix(kbf,lbf,:) )
         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo
#else
 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO REDUCTION(+:pot_hartree) PRIVATE(ii_tmp,ibf,jbf,kbf,lbf) 
 do ii=1,nsize_sparse
   ii_tmp = index_sparse(ii)
   ii_tmp = ii-1
   lbf = ii_tmp/nbf**3 + 1
   ii_tmp = ii_tmp - (lbf-1)*nbf**3
   kbf = ii_tmp/nbf**2 + 1
   ii_tmp = ii_tmp - (kbf-1)*nbf**2
   jbf = ii_tmp/nbf    + 1
   ii_tmp = ii_tmp - (jbf-1)*nbf
   ibf = ii_tmp + 1

   pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                + eri_buffer(ii) * SUM( p_matrix(kbf,lbf,:) )
 enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo
#endif

 title='=== Hartree contribution ==='
 call dump_out_matrix(PRINT_VOLUME,title,nbf,nspin,pot_hartree)

 ehartree = 0.5_dp*SUM(pot_hartree(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_hartree)

end subroutine setup_hartree

!=========================================================================
subroutine setup_exchange(PRINT_VOLUME,nbf,nspin,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: PRINT_VOLUME
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_exchange(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 real(dp)             :: spin_fact
 character(len=100)   :: title
!=====

 call start_clock(timing_exchange)

 spin_fact = REAL(-nspin+3,dp)

 pot_exchange(:,:,:)=0.0_dp

 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           pot_exchange(ibf,jbf,ispin) = pot_exchange(ibf,jbf,ispin) &
                      - eri(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact
         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo

 title='=== Exchange contribution ==='
 call dump_out_matrix(PRINT_VOLUME,title,nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange
