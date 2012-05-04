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
 integer              :: ii
 real(dp)             :: eri_tmp
!=====

 call start_clock(timing_hartree)

 pot_hartree(:,:,:)=0.0_dp

#ifndef LOW_MEMORY3
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
!$OMP DO REDUCTION(+:pot_hartree) PRIVATE(ibf,jbf,kbf,lbf) 
 do ii=1,nsize_sparse
   ibf = index_sparse(1,ii)
   jbf = index_sparse(2,ii)
   kbf = index_sparse(3,ii)
   lbf = index_sparse(4,ii)
   eri_tmp = eri_buffer_sparse(ii)

   pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                + eri_tmp * SUM( p_matrix(kbf,lbf,:) )

   if( ibf/=kbf .AND. jbf/=lbf ) then
     pot_hartree(kbf,lbf,ispin) = pot_hartree(kbf,lbf,ispin) &
                  + eri_tmp * SUM( p_matrix(ibf,jbf,:) )
     if(ibf/=jbf) then
       pot_hartree(jbf,ibf,ispin) = pot_hartree(jbf,ibf,ispin) &
                  + eri_tmp * SUM( p_matrix(kbf,lbf,:) )
       pot_hartree(kbf,lbf,ispin) = pot_hartree(kbf,lbf,ispin) &
                    + eri_tmp * SUM( p_matrix(jbf,ibf,:) )
       if(kbf/=lbf) then
         pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                      + eri_tmp * SUM( p_matrix(lbf,kbf,:) )
         pot_hartree(jbf,ibf,ispin) = pot_hartree(jbf,ibf,ispin) &
                      + eri_tmp * SUM( p_matrix(lbf,kbf,:) )
         pot_hartree(lbf,kbf,ispin) = pot_hartree(lbf,kbf,ispin) &
                      + eri_tmp * SUM( p_matrix(ibf,jbf,:) )
         pot_hartree(lbf,kbf,ispin) = pot_hartree(lbf,kbf,ispin) &
                      + eri_tmp * SUM( p_matrix(jbf,ibf,:) )
       endif
     else
       if(kbf/=lbf) then
         pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                      + eri_tmp * SUM( p_matrix(lbf,kbf,:) )
         pot_hartree(lbf,kbf,ispin) = pot_hartree(lbf,kbf,ispin) &
                      + eri_tmp * SUM( p_matrix(ibf,jbf,:) )
       endif
     endif
   else
     if(ibf/=jbf) then
       pot_hartree(jbf,ibf,ispin) = pot_hartree(jbf,ibf,ispin) &
                  + eri_tmp * SUM( p_matrix(kbf,lbf,:) )
       if(kbf/=lbf) then
         pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                    + eri_tmp * SUM( p_matrix(lbf,kbf,:) )
         pot_hartree(jbf,ibf,ispin) = pot_hartree(jbf,ibf,ispin) &
                    + eri_tmp * SUM( p_matrix(lbf,kbf,:) )
       endif
     else
       if(kbf/=lbf) then
         pot_hartree(ibf,jbf,ispin) = pot_hartree(ibf,jbf,ispin) &
                    + eri_tmp * SUM( p_matrix(lbf,kbf,:) )
       endif
     endif
   endif

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

!=========================================================================
subroutine setup_exchange_shortrange(PRINT_VOLUME,nbf,nspin,p_matrix,pot_exchange,eexchange)
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
                      - ( eri(kbf,ibf,lbf,jbf) - eri_scr(kbf,ibf,lbf,jbf) ) &     !  short-range only 
                       * p_matrix(kbf,lbf,ispin) / spin_fact 
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

end subroutine setup_exchange_shortrange
