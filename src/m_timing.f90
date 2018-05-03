!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the timings of the code
!
!=========================================================================
module m_timing
 use m_definitions
 use m_warning,only: die

 integer,parameter :: NTIMING=140

 integer,parameter :: timing_total               =  1

 integer,parameter :: timing_prescf              = 81
 integer,parameter :: timing_scf                 = 82
 integer,parameter :: timing_postscf             = 83

 integer,parameter :: timing_dft                 =  2
 integer,parameter :: timing_pola                =  3
 integer,parameter :: timing_gw_self             =  4
 integer,parameter :: timing_overlap             =  5
 integer,parameter :: timing_eri_4center         =  6
 integer,parameter :: timing_exchange            =  7
 integer,parameter :: timing_hartree             =  8
 integer,parameter :: timing_sqrt_density_matrix =  9
 integer,parameter :: timing_diago_h2p           = 10
 integer,parameter :: timing_rpa_static          = 11
 integer,parameter :: timing_mp2_energy          = 12
 integer,parameter :: timing_pt_self             = 13
 integer,parameter :: timing_eri_4center_eigen   = 14
 integer,parameter :: timing_single_excitation   = 15
 integer,parameter :: timing_eri_2center         = 16
 integer,parameter :: timing_eri_3center         = 17
 integer,parameter :: timing_eri_3center_eigen   = 18
 integer,parameter :: timing_vchiv               = 19
 integer,parameter :: timing_build_h2p           = 20
 integer,parameter :: timing_restart_file        = 21
 integer,parameter :: timing_diago_hamiltonian   = 22
 integer,parameter :: timing_hamiltonian_nuc     = 23
 integer,parameter :: timing_hamiltonian_kin     = 24
 integer,parameter :: timing_build_common        = 25
 integer,parameter :: timing_build_tddft         = 26
 integer,parameter :: timing_build_bse           = 27
 integer,parameter :: timing_spectrum            = 28
 integer,parameter :: timing_eri_screening       = 29
 integer,parameter :: timing_hamiltonian_ecp     = 30
 integer,parameter :: timing_sca_distr1          = 31
 integer,parameter :: timing_grid_generation     = 32
 integer,parameter :: timing_diis                = 33
 integer,parameter :: timing_approx_ham          = 34
 integer,parameter :: timing_sca_distr2          = 35
 integer,parameter :: timing_fno                 = 36
 integer,parameter :: timing_full_ci             = 37
 integer,parameter :: timing_gwgamma_self        = 38
 integer,parameter :: timing_ecp                 = 39
 integer,parameter :: timing_density_matrix      = 40
 integer,parameter :: timing_force               = 41
 integer,parameter :: timing_rpa_dynamic         = 42
 integer,parameter :: timing_ci_selfenergy       = 43
 integer,parameter :: timing_ham_ci              = 44
 integer,parameter :: timing_ci_diago            = 45
 integer,parameter :: timing_ci_write            = 46
 integer,parameter :: timing_ci_config           = 47
 integer,parameter :: timing_zeroes_ci           = 48
 integer,parameter :: timing_density_matrix_cmplx= 49
 integer,parameter :: timing_aomo_pola           = 50
 integer,parameter :: timing_aomo_ci             = 51
 integer,parameter :: timing_aomo_self           = 52
 
 integer,parameter :: timing_tmp0                = 90
 integer,parameter :: timing_tmp1                = 91
 integer,parameter :: timing_tmp2                = 92
 integer,parameter :: timing_tmp3                = 93
 integer,parameter :: timing_tmp4                = 94
 integer,parameter :: timing_tmp5                = 95
 integer,parameter :: timing_tmp6                = 96
 integer,parameter :: timing_tmp7                = 97
 integer,parameter :: timing_tmp8                = 98
 integer,parameter :: timing_tmp9                = 99

 integer,parameter :: timing_tddft_loop          = 110
 integer,parameter :: timing_tddft_fourier       = 111
 integer,parameter :: timing_tddft_one_iter      = 112
 integer,parameter :: timing_tddft_propagation   = 113
 integer,parameter :: timing_tddft_hamiltonian_fock = 114
 integer,parameter :: timing_print_cube_rho_tddft= 115
 integer,parameter :: timing_restart_tddft_file  = 116
 integer,parameter :: timing_propagate_diago     = 117
 integer,parameter :: timing_propagate_matmul    = 118
 integer,parameter :: timing_print_line_rho_tddft= 119

 integer           :: count_rate,count_max
 logical           :: time_running(NTIMING)
 real(dp)          :: time_start(NTIMING)
 real(dp)          :: timing(NTIMING)
 integer(dp)       :: calls(NTIMING)
 
contains


!=========================================================================
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
 integer            :: count_tmp
!=====
 
 if(time_running(itiming)) then
   write(stdout,*) 'clock # is already started:',itiming
   call die('error in start clock')
 endif

 time_running(itiming)=.TRUE.

 call system_clock(COUNT=count_tmp)
 time_start(itiming) = count_tmp
 calls(itiming) = calls(itiming) + 1

end subroutine start_clock
 

!=========================================================================
subroutine stop_clock(itiming)
 implicit none
 integer,intent(in) :: itiming
!=====
 integer            :: count_tmp
!===== 
  
 if(.NOT.time_running(itiming)) then
   write(stdout,*) 'clock # has not been started:',itiming
   call die('error in stop clock')
 endif


 time_running(itiming)=.FALSE. 

 call system_clock(COUNT=count_tmp)
 timing(itiming) = timing(itiming) + MODULO( count_tmp - NINT(time_start(itiming)) , count_max) / REAL(count_rate,dp)

end subroutine stop_clock


!=========================================================================
subroutine output_timing()
 implicit none
!=====
 logical,parameter :: ONLY_ONE_CALL = .TRUE.
!=====

 write(stdout,'(/,a,/)') '                 --- Timings in (s) and # of calls ---'

 call output_timing_line('Total time' ,timing_total ,0,ONLY_ONE_CALL)

 write(stdout,'(/,a,/)') '                 -------------------------------------'

<<<<<<< HEAD
 write(stdout,'(a30,6x,f12.2)')  'Total pre SCF',timing(timing_prescf)
 write(stdout,'(a30,6x,f12.2)')      'Total SCF',timing(timing_scf)
 write(stdout,'(a30,6x,f12.2)') 'Total post SCF',timing(timing_postscf) 
=======
 call output_timing_line('Total pre SCF' ,timing_prescf ,0,ONLY_ONE_CALL)
 call output_timing_line('Total SCF'     ,timing_scf    ,0,ONLY_ONE_CALL)
 call output_timing_line('Total post SCF',timing_postscf,0,ONLY_ONE_CALL)

>>>>>>> origin
 write(stdout,'(/,a,/)') '                 -------------------------------------'
 write(stdout,'(a,/)')   '                             Pre SCF'

 call output_timing_line('Integral pre-screening',timing_eri_screening,1)
 call output_timing_line('4-center integrals',timing_eri_4center,1)
 call output_timing_line('2-center integrals',timing_eri_2center,1)
 call output_timing_line('3-center integrals',timing_eri_3center,1)
 call output_timing_line('Overlap matrix S',timing_overlap,1)
 call output_timing_line('Approximate guess Hamiltonian',timing_approx_ham,1)
 call output_timing_line('Kinetic Hamiltonian',timing_hamiltonian_kin,1)
 call output_timing_line('Electron-nucleus Hamiltonian',timing_hamiltonian_nuc,1)
 call output_timing_line('ECP Hamiltonian',timing_hamiltonian_ecp,1)

 write(stdout,'(/,a,/)') '                 -------------------------------------'
 write(stdout,'(a,/)')   '                                 SCF'

 call output_timing_line('DFT Grid generation',timing_grid_generation,1)
 call output_timing_line('Density matrix',timing_density_matrix,1)
 call output_timing_line('Hartree potential',timing_hartree,1)
 call output_timing_line('Exchange operator',timing_exchange,1)
 call output_timing_line('DFT xc potential',timing_dft,1)
 call output_timing_line('Hamiltonian diagonalization',timing_diago_hamiltonian,1)
 call output_timing_line('Pulay DIIS mixing',timing_diis,1)
 call output_timing_line('RESTART file writing',timing_restart_file,1)
 call output_timing_line('Singles correction',timing_single_excitation,1)
 call output_timing_line('Virtual FNO generation',timing_fno,1)
 call output_timing_line('Forces',timing_force,1)

<<<<<<< HEAD
 write(stdout,*)
 write(stdout,'(a)') '                 -------------------------------------'
 write(stdout,*)

 if( calls(timing_sca_distr1)+calls(timing_sca_distr2) > 0 ) then
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing SCA1   ' ,timing(timing_sca_distr1),calls(timing_sca_distr2)
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'timing SCA2   ' ,timing(timing_sca_distr1),calls(timing_sca_distr2)
 endif

 if( calls(timing_tddft_loop) > 0 ) then
   write(stdout,'(a30,6x,f12.2,2x,i8)') 'TD-DFT Loop'                    ,timing(timing_tddft_loop),calls(timing_tddft_loop)
   write(stdout,'(a32,4x,f12.2,2x,i8)') 'Propagation for TD-DFT'         ,timing(timing_tddft_propagation),calls(timing_tddft_propagation)
   write(stdout,'(a34,2x,f12.2,2x,i8)') 'Diago in tddft propagation'     ,timing(timing_propagate_diago),calls(timing_propagate_diago)
   write(stdout,'(a34,2x,f12.2,2x,i8)') 'Matmul in tddft propagation'    ,timing(timing_propagate_matmul),calls(timing_propagate_matmul)

   write(stdout,'(a32,4x,f12.2,2x,i8)') 'Hamiltonian_fock calculation'   ,timing(timing_tddft_hamiltonian_fock),calls(timing_tddft_hamiltonian_fock)
   write(stdout,'(a34,2x,f12.2,2x,i8)') 'Complex density matrix'         ,timing(timing_density_matrix_cmplx),calls(timing_density_matrix_cmplx)

   write(stdout,'(a30,4x,f12.2,2x,i8)') 'RESTART_TDDFT file writing'     ,timing(timing_restart_tddft_file),calls(timing_restart_tddft_file)

   if(calls(timing_print_cube_rho_tddft)>0) then
      write(stdout,'(a32,4x,f12.2,2x,i8)') 'Cube density file writing'   ,timing(timing_print_cube_rho_tddft),calls(timing_print_cube_rho_tddft)
   end if
   if(calls(timing_print_line_rho_tddft)>0) then
      write(stdout,'(a32,4x,f12.2,2x,i8)') 'Line density file writing'   ,timing(timing_print_line_rho_tddft),calls(timing_print_line_rho_tddft)
   end if
 end if

=======
 write(stdout,'(/,a,/)') '                 -------------------------------------'
 write(stdout,'(a,/)')   '                            Post SCF'

 call output_timing_line('Response function chi',timing_pola,1)
 call output_timing_line('3-center AO to MO transform',timing_eri_3center_eigen,2)
 call output_timing_line('4-center AO to MO transform',timing_eri_4center_eigen,2)
 call output_timing_line('Static polarization for BSE',timing_rpa_static,2)
 call output_timing_line('Dynamic polarization',timing_rpa_dynamic,2)
 call output_timing_line('Build 2-particle Hamiltonian',timing_build_h2p,2)
 call output_timing_line('RPA part',timing_build_common,3)
 call output_timing_line('TDDFT part',timing_build_tddft,3)
 call output_timing_line('BSE part',timing_build_bse,3)
 call output_timing_line('Diago 2 particle H',timing_diago_h2p,2)
 call output_timing_line('Build W',timing_vchiv,2)
 call output_timing_line('Optical spectrum',timing_spectrum,2)

 call output_timing_line('GW self-energy',timing_gw_self,1)
 call output_timing_line('PT self-energy',timing_pt_self,1)
 call output_timing_line('GWGamma self-energy',timing_gwgamma_self,1)
 call output_timing_line('MP2 energy',timing_mp2_energy,1)


 call output_timing_line('Full CI for few electrons',timing_full_ci,1)
 call output_timing_line('Setup CI configurations',timing_ci_config,2)
 call output_timing_line('Screen CI Hamiltonian zeroes',timing_zeroes_ci,2)
 call output_timing_line('Build CI Hamiltonian',timing_ham_ci,2)
 call output_timing_line('CI diagonalization',timing_ci_diago,2)
 call output_timing_line('CI eigenvector file writing',timing_ci_write,2)
 call output_timing_line('CI self-energy',timing_ci_selfenergy,2)

 write(stdout,'(/,a,/)') '                 -------------------------------------'



 !
 ! Developer's timings for temporary use only!
>>>>>>> origin
 !
 if( ANY( calls(timing_tmp0:timing_tmp9) > 0 ) .OR. calls(timing_sca_distr1) > 0 .OR. calls(timing_sca_distr2) > 0 ) then
   write(stdout,'(a,/)')   '                            Testing'
   call output_timing_line('timing SCALAPACK tmp1',timing_sca_distr1,1)
   call output_timing_line('timing SCALAPACK tmp2',timing_sca_distr2,1)
   call output_timing_line('Tmp timing 0',timing_tmp0,1)
   call output_timing_line('Tmp timing 1',timing_tmp1,1)
   call output_timing_line('Tmp timing 2',timing_tmp2,1)
   call output_timing_line('Tmp timing 3',timing_tmp3,1)
   call output_timing_line('Tmp timing 4',timing_tmp4,1)
   call output_timing_line('Tmp timing 5',timing_tmp5,1)
   call output_timing_line('Tmp timing 6',timing_tmp6,1)
   call output_timing_line('Tmp timing 7',timing_tmp7,1)
   call output_timing_line('Tmp timing 8',timing_tmp8,1)
   call output_timing_line('Tmp timing 9',timing_tmp9,1)
   write(stdout,'(/,a,/)') '                 -------------------------------------'
 endif

<<<<<<< HEAD
 write(stdout,'(/,a)') '                 -------------------------------------'

=======
>>>>>>> origin

end subroutine output_timing


!=========================================================================
subroutine output_timing_line(title,itiming,level,only_one_call)
 implicit none

 character(len=*),intent(in) :: title
 integer,intent(in)          :: itiming
 integer,intent(in)          :: level
 logical,intent(in),optional :: only_one_call
!=====
 integer,parameter           :: max_length = 40
 character(len=50)           :: prepend
 integer                     :: lt,lp
!=====

 ! No writing if the timing counter has never been used
 if( calls(itiming) < 1 ) return

 lt = LEN_TRIM(title)

 select case(level)
 case(0)
   prepend = ''
 case(1)
   prepend = '|'
 case(2)
   prepend = '      |'
 case(3)
   prepend = '            |'
 end select

 lp = LEN_TRIM(prepend)

 prepend = TRIM(prepend) // REPEAT('-',max_length-lt-lp)

 if( PRESENT(only_one_call) ) then
   write(stdout,'(1x,a41,4x,f12.2)') TRIM(title),timing(itiming)
 else
   write(stdout,'(1x,a,1x,a,4x,f12.2,2x,i8)') TRIM(prepend),TRIM(title),timing(itiming),calls(itiming)
 endif


end subroutine output_timing_line


!=========================================================================
end module m_timing
