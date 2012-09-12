!=========================================================================
subroutine setup_hartree(print_volume,nbf,nspin,p_matrix,pot_hartree,ehartree)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_hartree(nbf,nbf,nspin)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
!=====

 call start_clock(timing_hartree)

 pot_hartree(:,:,:)=0.0_dp

 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
   do jbf=1,nbf
     do ibf=1,nbf
       if( .NOT. is_my_task(ibf,jbf) ) cycle
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
 call xsum(pot_hartree)


 title='=== Hartree contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_hartree)

 ehartree = 0.5_dp*SUM(pot_hartree(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_hartree)

end subroutine setup_hartree

!=========================================================================
subroutine setup_exchange(print_volume,nbf,nspin,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
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
     do lbf=1,nbf
       if( .NOT. is_my_task(lbf,jbf) ) cycle
       do ibf=1,nbf
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
 call xsum(pot_exchange)


 title='=== Exchange contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange

!=========================================================================
subroutine setup_exchange_longrange(print_volume,nbf,nspin,p_matrix,pot_exchange,eexchange)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
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
                      - eri_lr(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo

 title='=== Exchange contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_exchange)

 eexchange = 0.5_dp*SUM(pot_exchange(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange

!=========================================================================
subroutine read_potential(print_volume,nbf,nspin,p_matrix,pot_read,eread)
 use m_definitions
 use m_mpi
 use m_timing
 use m_eri
 implicit none
 integer,intent(in)   :: print_volume
 integer,intent(in)   :: nbf,nspin
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: pot_read(nbf,nbf,nspin)
 real(dp),intent(out) :: eread
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
 character(len=100)   :: title
 logical              :: file_exists
!=====

 pot_read(:,:,:)=0.0_dp

 inquire(file='manual_potential',exist=file_exists)
 if(file_exists) then
   open(unit=12,file='manual_potential',status='old')
   do ispin=1,nspin
     do jbf=1,nbf
       do ibf=1,nbf
         read(12,*) pot_read(ibf,jbf,ispin)
       enddo
     enddo
   enddo
   close(12)

 else
   stop'file not found: manual_potential'
 endif


 title='=== Read potential contribution ==='
 call dump_out_matrix(print_volume,title,nbf,nspin,pot_read)

 eread = 0.5_dp*SUM(pot_read(:,:,:)*p_matrix(:,:,:))

end subroutine read_potential
