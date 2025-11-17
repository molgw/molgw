!=========================================================================
! This file is part of MOLGW.
! Author: Marc Sarraute
!
! This module contains
! the drivers to link with ELPA
!=========================================================================
#include "molgw.h"
module m_elpa
  use m_definitions
  use m_timing
  use m_warning
  use m_memory
  use m_mpi

#if defined(HAVE_ELPA)
  use elpa
#endif

  implicit none

  private

#if defined(HAVE_ELPA)

!Public procedures
!Had to choose names different from those provided by elpa
  public :: elpa_func_allocate            ! Allocate a ELPA handle and set up MPI information
  public :: elpa_func_deallocate          ! Deallocate a ELPA handle
  public :: elpa_func_error_handler       ! Manage errors (print a readable message)
  public :: elpa_func_get_communicators   ! Get rows and cols communicators (not supposed to be called directly)
  public :: elpa_func_set_matrix          ! Set matrix specifications in a ELPA handle
  public :: elpa_func_solve_evp           ! Solve the diagonalization problem (use a ELPA handle)

  interface elpa_func_solve_evp
    module procedure elpa_func_solve_evp_real
    module procedure elpa_func_solve_evp_complex
  end interface elpa_func_solve_evp

!ELPA generalized handle
  type, public :: elpa_hdl_t
    logical :: is_allocated = .FALSE.
    logical :: matrix_is_set = .FALSE.
    class(elpa_t), pointer :: elpa
    integer :: mpi_comm_parent
    integer :: elpa_comm_rows, elpa_comm_cols
    integer :: process_row, process_col
    integer :: gpu = 0
  end type elpa_hdl_t

#endif

contains  !==============================================================================
!!***

#if defined(HAVE_ELPA)


!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_allocate
!! NAME
!!  elpa_func_allocate
!!
!! FUNCTION
!!  Allocate a ELPA handle and set it up with communicators specification
!!
!! INPUTS
!!  [blacs_ctx]= -- optional -- Blacs context
!!  [gpu]= -- optional -- Flag (0 or 1): use GPU version (currently only NVidia)
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)= ELPA handle
!!
!! SOURCE

subroutine elpa_func_allocate(elpa_hdl, gpu, blacs_ctx)

 !Arguments ------------------------------------
  integer, intent(in), optional :: gpu, blacs_ctx
  type(elpa_hdl_t), intent(inout) :: elpa_hdl

 !Local variables-------------------------------
  integer :: err, l_gpu, l_blacs_ctx
  logical :: debug_mode = .FALSE.
  character(len=10) :: varname

 ! *********************************************************************

  err = ELPA_OK
  ! if optional parameter is present, use it
  ! else use default value, i.e. don't use GPU
  l_gpu = 0
  if (PRESENT(gpu)) then
    if(gpu /=0 ) l_gpu = 1
  end if

  elpa_hdl%elpa => elpa_allocate(err)
  call elpa_func_error_handler(err_code=err, err_msg='Error in initialization')

  if(l_gpu == 1) then
    varname = "nvidia-gpu"
    if (err == ELPA_OK) call elpa_hdl%elpa%set(varname, l_gpu, err)
    call elpa_func_error_handler(err_code=err, err_msg='Error when enabling GPU on ELPA')
  end if

  if(debug_mode) then
    if (err == ELPA_OK) call elpa_hdl%elpa%set("debug", 1, err) 
    call elpa_func_error_handler(err_code=err, err_msg='Error when enabling debug on ELPA')
  end if

  if(PRESENT(blacs_ctx)) then
    if (err == ELPA_OK) call elpa_hdl%elpa%set("blacs_context", INT(blacs_ctx, kind=C_INT), err)
    call elpa_func_error_handler(err_code=err, err_varname=varname)
  end if

  elpa_hdl%is_allocated = .TRUE.

end subroutine elpa_func_allocate
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_deallocate
!! NAME
!!  elpa_func_deallocate_matrix
!!
!! FUNCTION
!!  Deallocate a ELPA handle
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)= ELPA handle
!!
!! SOURCE

subroutine elpa_func_deallocate(elpa_hdl)

  type(elpa_hdl_t), intent(inout) :: elpa_hdl
  !=====
  !=====
 
 
  call elpa_deallocate(elpa_hdl%elpa)
 
  elpa_hdl%matrix_is_set = .FALSE.
  elpa_hdl%is_allocated = .FALSE.

end subroutine elpa_func_deallocate
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_error_handler
!! NAME
!!  elpa_func_error_handler
!!
!! FUNCTION
!!  Handle ELPA errors
!!
!! INPUTS
!!  [err_code]= --optional-- Error code
!!  [err_msg]= --optional-- Generic error message
!!  [err_varname]= -- optional-- Name of the ELPA variable related to the error
!!
!! OUTPUT
!!  No output, only printing
!!
!! SOURCE

subroutine elpa_func_error_handler(err_code, err_msg, err_varname)

!Arguments ------------------------------------
 integer, optional :: err_code
 character(len=*), optional :: err_msg, err_varname

!Local variables-------------------------------
 integer :: err_code_
 character(len=500) :: msg
 character(len=100) :: err_strg

! *********************************************************************

 err_code_ = -100
 if (PRESENT(err_code)) err_code_ = err_code
 if (err_code_ == ELPA_OK) return

 err_strg = ''
 if (err_code_ == ELPA_ERROR) err_strg='ELPA_ERROR'
 if (err_code_ == ELPA_ERROR_ENTRY_READONLY) err_strg='ELPA_ERROR_ENTRY_READONLY'
 if (err_code_ == ELPA_ERROR_ENTRY_NOT_FOUND) err_strg='ELPA_ERROR_ENTRY_NOT_FOUND'
 if (err_code_ == ELPA_ERROR_ENTRY_ALREADY_SET) err_strg='ELPA_ERROR_ENTRY_ALREADY_SET'
 if (err_code_ == ELPA_ERROR_ENTRY_INVALID_VALUE) err_strg='ELPA_ERROR_ENTRY_INVALID_VALUE'
 if (err_code_ == ELPA_ERROR_ENTRY_NO_STRING_REPRESENTATION) err_strg='ELPA_ERROR_NO_STRING_REPRESENTATION'
 if (err_code_ == ELPA_ERROR_SETUP) err_strg='ELPA_ERROR_SETUP'
 if (err_code_ == ELPA_ERROR_CRITICAL) err_strg='ELPA_ERROR_CRITICAL'
 if (err_code_ == ELPA_ERROR_API_VERSION) err_strg='ELPA_ERROR_API_VERSION'
 if (err_code_ == ELPA_ERROR_AUTOTUNE_API_VERSION) err_strg='ELPA_ERROR_AUTOTUNE_API_VERSION'
 if (err_code_ == ELPA_ERROR_AUTOTUNE_OBJECT_CHANGED) err_strg='ELPA_ERROR_AUTOTUNE_OBJECT_CHANGED'
 if (err_code_ == ELPA_ERROR_CANNOT_OPEN_FILE) err_strg='ELPA_ERROR_CANNOT_OPEN_FILE'
 if (err_code_ == ELPA_ERROR_DURING_COMPUTATION) err_strg='ELPA_ERROR_DURING_COMPUTATION'

 write(msg, '(a)') 'ELPA library error!'
 if (PRESENT(err_msg)) then
   if (TRIM(err_msg) /= "") write(msg, '(3a)') TRIM(msg), CHAR(10), TRIM(err_msg)
 end if
 if (PRESENT(err_varname)) then
   if (TRIM(err_varname) /= "") write(msg, '(4a)') TRIM(msg), CHAR(10), 'Variable: ', TRIM(err_varname)
 end if
 if (TRIM(err_strg) /= "") write(msg, '(4a)') TRIM(msg), CHAR(10), 'Error code: ', TRIM(err_strg)
   call die(msg)

end subroutine elpa_func_error_handler
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_get_communicators
!! NAME
!!  elpa_func_get_communicators
!!
!! FUNCTION
!!  Wrapper to elpa_get_communicators ELPA function
!!
!! INPUTS
!!  mpi_comm_parent=Global communicator for the calculations (in)
!!  process_row=Row coordinate of the calling process in the process grid (in)
!!  process_col=Column coordinate of the calling process in the process grid (in)
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)= ELPA handle
!!
!! SOURCE

subroutine elpa_func_get_communicators(elpa_hdl, mpi_comm_parent, process_row, process_col)

!Arguments ------------------------------------
 integer, intent(in)  :: mpi_comm_parent, process_row, process_col
 type(elpa_hdl_t), intent(inout) :: elpa_hdl

!Local variables-------------------------------
 integer  :: err
 character(len=20) :: varname
 character(len=200):: msg

! *********************************************************************

 err = ELPA_OK
 varname = ''

 if (.NOT. elpa_hdl%is_allocated) then
   call die('ELPA handle not allocated!')
 end if

 if (err == ELPA_OK) then
   varname = 'mpi_comm_parent'
   call elpa_hdl%elpa%set(TRIM(varname), mpi_comm_parent, err)
 end if
 if (err == ELPA_OK) then
   varname = 'process_row'
   call elpa_hdl%elpa%set(TRIM(varname), process_row, err)
 end if
 if (err==ELPA_OK) then
   varname = 'process_col'
   call elpa_hdl%elpa%set(TRIM(varname), process_col, err)
 end if
 if (err == ELPA_OK) then
   err = elpa_hdl%elpa%setup()
   write(msg, '(3a)') "Error when setting '", varname, "'during ELPA communicators setup"
   call elpa_func_error_handler(err_code=err, err_msg=msg)
 endif

 elpa_hdl%mpi_comm_parent = mpi_comm_parent
 elpa_hdl%process_row = process_row
 elpa_hdl%process_col = process_col

 call elpa_func_error_handler(err_code=err, err_msg='Error in elpa_get_communicators', err_varname=varname)

end subroutine elpa_func_get_communicators
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_set_matrix
!! NAME
!!  elpa_func_set_matrix
!!
!! FUNCTION
!!  Set parameters decribing a matrix and it's MPI distribution
!!  in a ELPA handle
!!
!! INPUTS
!!  na=Order of matrix A
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  nev=Number of eigenvalues needed.
!!  local_nrows=Leading dimension of A
!!  local_ncols=Local columns of matrixes A and Q (eigenvectors)
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_set_matrix(elpa_hdl, na, nblk, nev, local_nrows, local_ncols)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: na, nblk, nev, local_nrows, local_ncols
 type(elpa_hdl_t), intent(inout) :: elpa_hdl
!arrays

!Local variables-------------------------------
 integer :: err
 character(len=15) :: varname

! *********************************************************************

 err = ELPA_OK
 varname = ''

 if (.NOT. elpa_hdl%is_allocated) then
   call die('ELPA handle not allocated!')
 end if

 if (err == ELPA_OK) then
   varname = "na"
   call elpa_hdl%elpa%set(TRIM(varname), na, err)
 end if
 if (err == ELPA_OK) then
   varname = "nblk"
   call elpa_hdl%elpa%set(TRIM(varname), nblk, err)
 end if
 if (err == ELPA_OK) then
   varname = "nev"
   call elpa_hdl%elpa%set(TRIM(varname), nev, err)
 end if
 if (err == ELPA_OK) then
   varname = "local_nrows"
   call elpa_hdl%elpa%set(TRIM(varname), local_nrows, err)
 end if
 if (err == ELPA_OK) then
   varname = "local_ncols"
   call elpa_hdl%elpa%set(TRIM(varname), local_ncols, err)
 end if

 call elpa_func_error_handler(err_code=err, err_msg='Error during matrix initialization', err_varname=varname)

 elpa_hdl%matrix_is_set = .TRUE.

end subroutine elpa_func_set_matrix
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_solve_evp_real
!! NAME
!!  elpa_func_solve_evp_real
!!
!! FUNCTION
!!  Wrapper to elpa_solve_evp_real ELPA function
!!
!! INPUTS
!!  nev=Number of eigenvalues needed.
!!  use_two_stage=if TRUE, use ELPA 2-stage solver (better on CPU, worse on GPU)
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(local_nrows, local_ncols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na, na))
!!                     even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!!  aa(local_nrows, local_ncols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_solve_evp_real(elpa_hdl, aa, qq, ev, nev, use_two_stage)

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: nev
 type(elpa_hdl_t), intent(inout) :: elpa_hdl
 logical :: use_two_stage
!arrays
 real(dp), intent(inout) :: aa(:, :)
 real(dp), intent(out) :: ev(:), qq(:, :)

!Local variables-------------------------------
 integer :: err
 logical  :: success

! *********************************************************************

 success = .TRUE.
 err = ELPA_OK

 if (.NOT. elpa_hdl%is_allocated) then
   call die('ELPA handle not allocated!')
 end if
 if (.NOT. elpa_hdl%matrix_is_set) then
   call die('Matrix not set in ELPA handle!')
 end if
 if(SIZE(aa) /= elpa_hdl%elpa%local_nrows * elpa_hdl%elpa%local_ncols) call die('BUG: matrix A has wrong sizes!')
 if(SIZE(qq) /= elpa_hdl%elpa%local_nrows * elpa_hdl%elpa%local_ncols) call die('BUG: matrix Q has wrong sizes!')
 if(SIZE(ev) /= elpa_hdl%elpa%na) call die('BUG: matrix EV has wrong sizes!')

 if(use_two_stage) then
   if (err == ELPA_OK) call elpa_hdl%elpa%set("solver", ELPA_SOLVER_2STAGE, err)
 else
   if (err == ELPA_OK) call elpa_hdl%elpa%set("solver", ELPA_SOLVER_1STAGE, err)
 end if
 if (err == ELPA_OK) call elpa_hdl%elpa%eigenvectors(aa, ev, qq, err)
 success = (err == ELPA_OK)

 if (.NOT. success) call elpa_func_error_handler(err_msg='Error in solve_evp_real!')

end subroutine elpa_func_solve_evp_real
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_solve_evp_complex
!! NAME
!!  elpa_func_solve_evp_complex
!!
!! FUNCTION
!!  Wrapper to elpa_solve_evp_complex ELPA function
!!
!! INPUTS
!!  nev=Number of eigenvalues needed.
!!  use_two_stage=if TRUE, use ELPA 2-stage solver (better on CPU, worse on GPU)
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(local_nrows, local_ncols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na, na))
!!                      even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!!  aa(local_nrows, local_ncols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_solve_evp_complex(elpa_hdl, aa, qq, ev, nev, use_two_stage)

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: nev
 type(elpa_hdl_t), intent(inout) :: elpa_hdl
 logical :: use_two_stage
!arrays
 complex(dp), intent(inout) :: aa(:, :)
 real(dp), intent(out) :: ev(:)
 complex(dp), intent(out) :: qq(:, :)

!Local variables-------------------------------
 integer :: err
 logical  :: success

! *********************************************************************

 success = .TRUE.
 err = ELPA_OK

 if (.NOT. elpa_hdl%is_allocated) then
   call die('ELPA handle not allocated!')
 end if
 if (.NOT. elpa_hdl%matrix_is_set) then
   call die('Matrix not set in ELPA handle!')
 end if
 if(SIZE(aa) /= elpa_hdl%elpa%local_nrows * elpa_hdl%elpa%local_ncols) call die('BUG: matrix A has wrong sizes!')
 if(SIZE(qq) /= elpa_hdl%elpa%local_nrows * elpa_hdl%elpa%local_ncols) call die('BUG: matrix Q has wrong sizes!')
 if(SIZE(ev) /= elpa_hdl%elpa%na) call die('BUG: matrix EV has wrong sizes!')

 if(use_two_stage) then
   if (err == ELPA_OK) call elpa_hdl%elpa%set("solver", ELPA_SOLVER_2STAGE, err)
 else
   if (err == ELPA_OK) call elpa_hdl%elpa%set("solver", ELPA_SOLVER_1STAGE, err)
 end if
 if (err == ELPA_OK) call elpa_hdl%elpa%eigenvectors(aa, ev, qq, err)
 success = (err == ELPA_OK)

 if (.NOT. success) call elpa_func_error_handler(err_msg='Error in solve_evp_complex!')

end subroutine elpa_func_solve_evp_complex
!!***

!----------------------------------------------------------------------

#endif

end module m_elpa
!!***
