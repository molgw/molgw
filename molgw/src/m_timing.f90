!=========================================================================
module m_timing
 use m_definitions


 integer,parameter :: timing_total             = 90

 integer,parameter :: timing_prescf            = 81
 integer,parameter :: timing_scf               = 82
 integer,parameter :: timing_postscf           = 83

 integer,parameter :: timing_dft               = 2
 integer,parameter :: timing_pola              = 3
 integer,parameter :: timing_self              = 4
 integer,parameter :: timing_prodbasis         = 5
 integer,parameter :: timing_eri_4center       = 6
 integer,parameter :: timing_exchange          = 7
 integer,parameter :: timing_hartree           = 8
 integer,parameter :: timing_overlap3          = 9
 integer,parameter :: timing_diago_h2p         = 10
 integer,parameter :: timing_inversion_s2p     = 11
 integer,parameter :: timing_mp2_energy        = 12
 integer,parameter :: timing_mp2_self          = 13
 integer,parameter :: timing_basis_transform   = 14
 integer,parameter :: timing_single_excitation = 15
 integer,parameter :: timing_eri_2center       = 16
 integer,parameter :: timing_eri_3center       = 17
 integer,parameter :: timing_eri_3center_eigen = 18
 integer,parameter :: timing_buildw            = 19
 integer,parameter :: timing_build_h2p         = 20
 integer,parameter :: timing_restart_file      = 21
 integer,parameter :: timing_diago_hamiltonian = 22
 integer,parameter :: timing_hamiltonian_nuc   = 23
 integer,parameter :: timing_hamiltonian_kin   = 24
 integer,parameter :: timing_build_common      = 25
 integer,parameter :: timing_build_tddft       = 26
 integer,parameter :: timing_build_bse         = 27
 integer,parameter :: timing_spectrum          = 28
 integer,parameter :: timing_eri_screening     = 29
 
 integer,parameter :: timing_tmp1              = 91
 integer,parameter :: timing_tmp2              = 92
 integer,parameter :: timing_tmp3              = 93
 integer,parameter :: timing_tmp4              = 94
 integer,parameter :: timing_tmp5              = 95
 integer,parameter :: timing_tmp6              = 96
 integer,parameter :: timing_tmp7              = 97
 integer,parameter :: timing_tmp8              = 98
 integer,parameter :: timing_tmp9              = 99

 integer,parameter :: NTIMING=100
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
   write(stdout,*) 'clock # is already started:',itiming
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
  
 if(.NOT.time_running(itiming)) then
   write(stdout,*) 'clock # has not been started:',itiming
   stop'error in start clock'
 endif


 time_running(itiming)=.FALSE. 

 call system_clock(COUNT=count_tmp)
 timing(itiming) = timing(itiming) + MODULO( count_tmp - NINT(time_start(itiming)) , count_max) / REAL(count_rate,dp)

end subroutine stop_clock

subroutine output_timing()
 implicit none
!=====

 write(stdout,'(/,a,/)') '                 --- Timings in (s) and # of calls ---'

 write(stdout,'(a30,6x,f12.2)')  'Total time',timing(timing_total)
 write(stdout,'(/,a,/)') '                 ----------------------'

 write(stdout,'(a30,6x,f12.2)')  'Total pre SCF',timing(timing_prescf)
 write(stdout,'(a30,6x,f12.2)')      'Total SCF',timing(timing_scf)
 write(stdout,'(a30,6x,f12.2)') 'Total post SCF',timing(timing_postscf)
 write(stdout,'(/,a,/)') '                 ----------------------'

 write(stdout,'(a30,6x,f12.2,2x,i8)')          'integral screening' ,timing(timing_eri_screening),calls(timing_eri_screening)
 write(stdout,'(a30,6x,f12.2,2x,i8)')          '4-center integrals' ,timing(timing_eri_4center),calls(timing_eri_4center)
 if( calls(timing_eri_2center) > 0 ) then
   write(stdout,'(a30,6x,f12.2,2x,i8)')        '       2-center integrals' ,timing(timing_eri_2center),calls(timing_eri_2center)
   write(stdout,'(a30,6x,f12.2,2x,i8)')        '       3-center integrals' ,timing(timing_eri_3center),calls(timing_eri_3center)
 endif
 write(stdout,'(a30,6x,f12.2,2x,i8)')          'Kinetic Hamiltonian',timing(timing_hamiltonian_kin),calls(timing_hamiltonian_kin)
 write(stdout,'(a30,6x,f12.2,2x,i8)')  'Electron-nuclei Hamiltonian',timing(timing_hamiltonian_nuc),calls(timing_hamiltonian_nuc)

 write(stdout,*)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'Hartree'          ,timing(timing_hartree),calls(timing_hartree)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'Exchange'         ,timing(timing_exchange),calls(timing_exchange)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'DFT xc'           ,timing(timing_dft),calls(timing_dft)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'Hamiltonian diago',timing(timing_diago_hamiltonian),calls(timing_diago_hamiltonian)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'RESTART file writing',timing(timing_restart_file),calls(timing_restart_file)
 write(stdout,*)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'Single Excit.'   ,timing(timing_single_excitation),calls(timing_single_excitation)
 write(stdout,*)

 write(stdout,'(a30,6x,f12.2,2x,i8)') 'Total chi polarization' ,timing(timing_pola),calls(timing_pola)
 if( calls(timing_eri_3center_eigen) > 0 ) then
   write(stdout,'(a30,6x,f12.2,2x,i8)')'Rotation 3-center integrals' ,timing(timing_eri_3center_eigen),calls(timing_eri_3center_eigen)
 endif
 if( calls(timing_basis_transform) > 0 ) then
   write(stdout,'(a32,4x,f12.2,2x,i8)') 'ERI basis transform' ,timing(timing_basis_transform),calls(timing_basis_transform)
 endif
 write(stdout,'(a32,4x,f12.2,2x,i8)') '    Build 2 particle H' ,timing(timing_build_h2p),calls(timing_build_h2p)
 write(stdout,'(a34,2x,f12.2,2x,i8)') '           Common part' ,timing(timing_build_common),calls(timing_build_common)
 write(stdout,'(a34,2x,f12.2,2x,i8)') '            TDDFT part' ,timing(timing_build_tddft),calls(timing_build_tddft)
 write(stdout,'(a34,2x,f12.2,2x,i8)') '              BSE part' ,timing(timing_build_bse),calls(timing_build_bse)
 write(stdout,'(a32,4x,f12.2,2x,i8)') '    Diago 2 particle H' ,timing(timing_diago_h2p),calls(timing_diago_h2p)
 write(stdout,'(a32,4x,f12.2,2x,i8)') '   Invert 2 particle S' ,timing(timing_inversion_s2p),calls(timing_inversion_s2p)
 write(stdout,'(a32,4x,f12.2,2x,i8)') '               Build W' ,timing(timing_buildw),calls(timing_buildw)
 write(stdout,'(a32,4x,f12.2,2x,i8)') '      Optical spectrum' ,timing(timing_spectrum),calls(timing_spectrum)
 write(stdout,*)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'GW self-energy'  ,timing(timing_self),calls(timing_self)
 write(stdout,*)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'MP2 energy'      ,timing(timing_mp2_energy),calls(timing_mp2_energy)
 write(stdout,'(a30,6x,f12.2,2x,i8)') 'MP2 self-energy' ,timing(timing_mp2_self),calls(timing_mp2_self)
 write(stdout,'(a)') '                 ----------------------'

 !
 ! developer's timings
 if( ANY(timing(timing_tmp1:timing_tmp9)>1.d-5) ) then
   write(stdout,*)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp1   ' ,timing(timing_tmp1),calls(timing_tmp1)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp2   ' ,timing(timing_tmp2),calls(timing_tmp2)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp3   ' ,timing(timing_tmp3),calls(timing_tmp3)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp4   ' ,timing(timing_tmp4),calls(timing_tmp4)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp5   ' ,timing(timing_tmp5),calls(timing_tmp5)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp6   ' ,timing(timing_tmp6),calls(timing_tmp6)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp7   ' ,timing(timing_tmp7),calls(timing_tmp7)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp8   ' ,timing(timing_tmp8),calls(timing_tmp8)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing tmp9   ' ,timing(timing_tmp9),calls(timing_tmp9)
   write(stdout,*)
 endif

end subroutine

end module m_timing
