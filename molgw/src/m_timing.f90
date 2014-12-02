!=========================================================================
#include "macros.h"
!=========================================================================
module m_timing
 use m_definitions
 use m_mpi


 integer,parameter :: timing_total             = 1
 integer,parameter :: timing_dft               = 2
 integer,parameter :: timing_pola              = 3
 integer,parameter :: timing_self              = 4
 integer,parameter :: timing_prodbasis         = 5
 integer,parameter :: timing_eri               = 6
 integer,parameter :: timing_exchange          = 7
 integer,parameter :: timing_hartree           = 8
 integer,parameter :: timing_overlap3          = 9
 integer,parameter :: timing_diago_h2p         = 10
 integer,parameter :: timing_inversion_s2p     = 11
 integer,parameter :: timing_mp2_energy        = 12
 integer,parameter :: timing_mp2_self          = 13
 integer,parameter :: timing_basis_transform   = 14
 integer,parameter :: timing_single_excitation = 15
 integer,parameter :: timing_eri_auxil         = 16
 
 integer,parameter :: timing_tmp1          = 41
 integer,parameter :: timing_tmp2          = 42
 integer,parameter :: timing_tmp3          = 43
 integer,parameter :: timing_tmp4          = 44
 integer,parameter :: timing_tmp5          = 45
 integer,parameter :: timing_tmp6          = 46
 integer,parameter :: timing_tmp7          = 47
 integer,parameter :: timing_tmp8          = 48
 integer,parameter :: timing_tmp9          = 49

 integer,parameter :: NTIMING=50
 integer           :: count_rate,count_max
 logical           :: time_running(NTIMING)
 real(dp)          :: time_start(NTIMING)
 real(dp)          :: timing(NTIMING)
 integer(dp)       :: calls(NTIMING)
 
contains

subroutine init_timing()
 implicit none
 
 time_running(:) = .FALSE.
 timing(:)       = 0.0_dp
 calls(:)        = 0

 call system_clock(COUNT_RATE=count_rate,COUNT_MAX=count_max)

end subroutine

subroutine start_clock(itiming)
 implicit none
 integer,intent(in) :: itiming
!=====
 real(dp)           :: time_tmp
 integer            :: count_tmp
!=====
 
 if(time_running(itiming)) then
   WRITE_MASTER(*,*) 'clock # is already started:',itiming
   stop'error in start clock'
 endif

 time_running(itiming)=.TRUE.

 call system_clock(COUNT=count_tmp)
 time_start(itiming) = count_tmp
 calls(itiming) = calls(itiming) + 1

end subroutine start_clock
 
subroutine stop_clock(itiming)
 implicit none
 integer,intent(in) :: itiming
!=====
 real(dp)           :: time_tmp
 integer            :: count_tmp
!===== 
  
 if(.NOT.time_running(itiming)) stop'error in start clock'

 time_running(itiming)=.FALSE. 

 call system_clock(COUNT=count_tmp)
 timing(itiming) = timing(itiming) + MODULO( count_tmp - NINT(time_start(itiming)) , count_max) / REAL(count_rate,dp)

end subroutine stop_clock

subroutine output_timing()
 implicit none
!=====

 WRITE_MASTER(*,'(/,a)') '                 --- Timings in (s) and # of calls ---'
 WRITE_MASTER(*,'(a30,2x,f12.2)') 'Total time',timing(timing_total)
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a30,2x,f12.2)')    '4-center integrals' ,timing(timing_eri)
 if( calls(timing_eri_auxil) > 0 ) then
   WRITE_MASTER(*,'(a30,2x,f12.2)')    '2- and 3-center integrals' ,timing(timing_eri_auxil)
 endif
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'Hartree'         ,timing(timing_hartree),calls(timing_hartree)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'Exchange'        ,timing(timing_exchange),calls(timing_exchange)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'DFT xc'          ,timing(timing_dft),calls(timing_dft)
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'ERI basis transform' ,timing(timing_basis_transform),calls(timing_basis_transform)
 WRITE_MASTER(*,*)

 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'Total RPA polarization' ,timing(timing_pola),calls(timing_pola)

 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'Diago 2 particle H' ,timing(timing_diago_h2p),calls(timing_diago_h2p)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'Invert 2 particle S' ,timing(timing_inversion_s2p),calls(timing_inversion_s2p)
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'GW self-energy'  ,timing(timing_self),calls(timing_self)
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'MP2 energy'      ,timing(timing_mp2_energy),calls(timing_mp2_energy)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'MP2 self-energy' ,timing(timing_mp2_self),calls(timing_mp2_self)
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'Single Excit.'   ,timing(timing_single_excitation),calls(timing_single_excitation)
 WRITE_MASTER(*,'(a)') '                 ----------------------'

 !
 ! developer's timings
 if( ANY(timing(timing_tmp1:timing_tmp9)>1.d-5) ) then
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp1   ' ,timing(timing_tmp1),calls(timing_tmp1)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp2   ' ,timing(timing_tmp2),calls(timing_tmp2)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp3   ' ,timing(timing_tmp3),calls(timing_tmp3)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp4   ' ,timing(timing_tmp4),calls(timing_tmp4)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp5   ' ,timing(timing_tmp5),calls(timing_tmp5)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp6   ' ,timing(timing_tmp6),calls(timing_tmp6)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp7   ' ,timing(timing_tmp7),calls(timing_tmp7)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp8   ' ,timing(timing_tmp8),calls(timing_tmp8)
   WRITE_MASTER(*,'(a30,2x,f12.2,2x,i8)') 'timing tmp9   ' ,timing(timing_tmp9),calls(timing_tmp9)
   WRITE_MASTER(*,*)
 endif

end subroutine

end module m_timing
