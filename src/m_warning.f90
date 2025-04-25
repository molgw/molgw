!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the management of the warning and error messages.
!
!=========================================================================
#include "molgw.h"
module m_warning
  use m_definitions


  integer, parameter          :: NWARNINGMAX=100
  integer, private            :: nwarning
  character(len=128), private :: warning_list(NWARNINGMAX)

  character(len=128)         :: msg


contains


!=========================================================================
subroutine init_warning()
  implicit none

  nwarning=0
  warning_list(:)=''

end subroutine init_warning


!=========================================================================
subroutine issue_warning(msgw)
  implicit none
  character(len=*), intent(in) :: msgw
  !=====

  write(stdout, '(/,a,a)') ' WARNING: ', TRIM(msgw)

  !
  ! Eliminate the storage of identical warnings
  if(ANY(warning_list(1:nwarning)==TRIM(msgw))) return

  nwarning = nwarning + 1
  warning_list(nwarning) = TRIM(msgw)

end subroutine issue_warning


!=========================================================================
subroutine output_all_warnings(unit_yaml)
  implicit none

  integer, intent(in), optional :: unit_yaml
  !=====
  integer :: iwarning
  !=====

  if(nwarning>0) then
    write(stdout, '(/,a,/)') ' SUMMARY of all the WARNINGS issued during the run:'
    do iwarning=1, nwarning
      write(stdout, '(i2,a,5x,a)') iwarning, '.', warning_list(iwarning)
    enddo
    write(stdout, '(/)')
  endif
  if( PRESENT(unit_yaml) ) then
    write(unit_yaml, '(/,a)')  'warnings:'
    do iwarning=1, nwarning
      write(unit_yaml, '(4x,a,a)') '- ', warning_list(iwarning)
    enddo
  endif

end subroutine output_all_warnings


!=========================================================================
subroutine die(msg)
  implicit none
  character(*), intent(in) :: msg
  !=====
  !=====

  write(stdout, '(/,a)') '=============================='
  write(stdout, '(a,a)') 'STOP: ',msg
  write(stdout, '(a,/)') '=============================='
  write(stderr, '(a,a)') 'STOP: ',msg
  error stop

end subroutine die


!=========================================================================
subroutine assert_experimental()
  implicit none
  !=====
  !=====

#if defined(ACTIVATE_EXPERIMENTAL)
  call issue_warning('Activating an experimental part of the code. Hopefully you know what you are doing')
#else
  write(stdout, *) 'This part of the code is experimental. It should not be used unless you know what you are doing'
  call die('Compilation with -DACTIVATE_EXPERIMENTAL is required')
#endif

end subroutine assert_experimental


!=========================================================================
end module m_warning
!=========================================================================
