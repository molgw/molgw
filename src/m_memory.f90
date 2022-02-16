!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the memory management
!=========================================================================
module m_memory
  use m_definitions
  use m_warning,only: die

  real(dp),private :: total_memory = 0.0_dp    ! Total memory occupied
                                               ! by the big arrays in Mb
  real(dp),private :: peak_memory  = 0.0_dp    ! Max memory occupied
                                               ! by the big arrays in Mb

  interface clean_allocate
    module procedure clean_allocate_i1d
    module procedure clean_allocate_i2d
    module procedure clean_allocate_s1d
    module procedure clean_allocate_s2d
    module procedure clean_allocate_s3d
    module procedure clean_allocate_s4d
    module procedure clean_allocate_1d
    module procedure clean_allocate_1d_int8
    module procedure clean_allocate_2d
    module procedure clean_allocate_2d_range
    module procedure clean_allocate_3d
    module procedure clean_allocate_3d_range
    module procedure clean_allocate_4d
    module procedure clean_allocate_4d_range
    module procedure clean_allocate_c1d
    module procedure clean_allocate_c2d
    module procedure clean_allocate_c3d
    module procedure clean_allocate_c4d
  end interface

  interface clean_deallocate
    module procedure clean_deallocate_i1d
    module procedure clean_deallocate_i2d
    module procedure clean_deallocate_s1d
    module procedure clean_deallocate_s2d
    module procedure clean_deallocate_s3d
    module procedure clean_deallocate_s4d
    module procedure clean_deallocate_1d
    module procedure clean_deallocate_2d
    module procedure clean_deallocate_3d
    module procedure clean_deallocate_4d
    module procedure clean_deallocate_c1d
    module procedure clean_deallocate_c2d
    module procedure clean_deallocate_c3d
    module procedure clean_deallocate_c4d
  end interface


contains


!=========================================================================
subroutine total_memory_statement()
  implicit none
  !=====
  !=====
 
  write(stdout,'(/,a,/)') '                 ----------------------------'
  write(stdout,'(a)')     '                 --- Memory usage summary ---'
 
  write(stdout,'(/,a)') ' Memory that was not deallocated properly'
  write(stdout,'(a30,f9.3)') ' Memory (Gb): ',total_memory / 1024._dp
 
  write(stdout,'(/,a)') ' Maximum memory used during the run'
  write(stdout,'(a30,f9.3)') ' Peak memory (Gb): ',peak_memory / 1024._dp
 
  write(stdout,*)

end subroutine total_memory_statement


!=========================================================================
function get_peak_memory() RESULT(peak)
  implicit none
 
  real(dp) :: peak
  !=====
  !=====
 
  peak = peak_memory / 1024.0_dp


end function get_peak_memory


!=========================================================================
subroutine clean_allocate_i1d(array_name,array,n1)
  implicit none

  character(len=*),intent(in)       :: array_name
  integer,allocatable,intent(inout) :: array(:)
  integer,intent(in)                :: n1
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(4,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1),stat=info)

  if(info/=0) then
    write(stdout,'(a,a)')    ' Failure when allocating ',array_name
    write(stdout,'(a,f9.3)') ' with size (Gb) ',mem_mb / 1024.0_dp
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_i1d


!=========================================================================
subroutine clean_allocate_i2d(array_name,array,n1,n2)
  implicit none

  character(len=*),intent(in)       :: array_name
  integer,allocatable,intent(inout) :: array(:,:)
  integer,intent(in)                :: n1,n2
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(4,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2),stat=info)

  if(info/=0) then
    write(stdout,'(a,a)')    ' Failure when allocating ',array_name
    write(stdout,'(a,f9.3)') ' with size (Gb) ',mem_mb / 1024.0_dp
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_i2d


!=========================================================================
subroutine clean_allocate_1d(array_name,array,n1)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(dp),allocatable,intent(inout) :: array(:)
  integer,intent(in)                 :: n1
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_1d


!=========================================================================
subroutine clean_allocate_1d_int8(array_name,array,n1)
  implicit none

  character(len=*),intent(in)        :: array_name
  real(dp),allocatable,intent(inout) :: array(:)
  integer(kind=int8),intent(in)      :: n1
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated ->'//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif

  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_1d_int8


!=========================================================================
subroutine clean_allocate_2d(array_name,array,n1,n2,verbose)
  implicit none

  character(len=*),intent(in)        :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:)
  integer,intent(in)                 :: n1,n2
  integer,optional,intent(in)        :: verbose
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  if(PRESENT(verbose)) then
   if(verbose/=-1) call write_memory_allocate(array_name,mem_mb)
  else
   call write_memory_allocate(array_name,mem_mb)
  endif

end subroutine clean_allocate_2d


!=========================================================================
subroutine clean_allocate_2d_range(array_name,array,n1s,n1f,n2s,n2f,verbose)
  implicit none

  character(len=*),intent(in)        :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:)
  integer,intent(in)                 :: n1s,n1f,n2s,n2f
  integer,optional,intent(in)        :: verbose
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1f-n1s+1,dp) * REAL(n2f-n2s+1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1s:n1f,n2s:n2f),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  if(PRESENT(verbose)) then
   if(verbose/=-1) call write_memory_allocate(array_name,mem_mb)
  else
   call write_memory_allocate(array_name,mem_mb)
  endif

end subroutine clean_allocate_2d_range


!=========================================================================
subroutine clean_allocate_3d(array_name,array,n1,n2,n3,verbose)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:,:)
  integer,intent(in)                 :: n1,n2,n3
  integer,optional,intent(in)        :: verbose
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2,n3),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  if(PRESENT(verbose)) then
   if(verbose/=-1) call write_memory_allocate(array_name,mem_mb)
  else
   call write_memory_allocate(array_name,mem_mb)
  endif

end subroutine clean_allocate_3d


!=========================================================================
subroutine clean_allocate_3d_range(array_name,array,n1s,n1f,n2s,n2f,n3s,n3f)
  implicit none

  character(len=*),intent(in)        :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:,:)
  integer,intent(in)                 :: n1s,n1f,n2s,n2f,n3s,n3f
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1f-n1s+1,dp) * REAL(n2f-n2s+1,dp) * REAL(n3f-n3s+1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1s:n1f,n2s:n2f,n3s:n3f),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_3d_range


!=========================================================================
subroutine clean_allocate_4d(array_name,array,n1,n2,n3,n4,verbose)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:,:,:)
  integer,intent(in)                 :: n1,n2,n3,n4
  integer,optional,intent(in)        :: verbose
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2,n3,n4),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  if(PRESENT(verbose)) then
   if(verbose/=-1) call write_memory_allocate(array_name,mem_mb)
  else
   call write_memory_allocate(array_name,mem_mb)
  endif

end subroutine clean_allocate_4d


!=========================================================================
subroutine clean_allocate_4d_range(array_name,array,n1s,n1f,n2s,n2f,n3s,n3f,n4s,n4f,verbose)
  implicit none

  character(len=*),intent(in)        :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:,:,:)
  integer,intent(in)                 :: n1s,n1f,n2s,n2f,n3s,n3f,n4s,n4f
  integer,optional,intent(in)        :: verbose
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(dp,dp) * REAL(n1f-n1s+1,dp) * REAL(n2f-n2s+1,dp) * REAL(n3f-n3s+1,dp) * REAL(n4f-n4s+1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1s:n1f,n2s:n2f,n3s:n3f,n4s:n4f),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  if(PRESENT(verbose)) then
   if(verbose/=-1) call write_memory_allocate(array_name,mem_mb)
  else
   call write_memory_allocate(array_name,mem_mb)
  endif

end subroutine clean_allocate_4d_range


!=========================================================================
subroutine clean_allocate_c1d(array_name,array,n1)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:)
  integer,intent(in)                    :: n1
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_c1d


!=========================================================================
subroutine clean_allocate_c2d(array_name,array,n1,n2)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:,:)
  integer,intent(in)                    :: n1,n2
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_c2d


!=========================================================================
subroutine clean_allocate_c3d(array_name,array,n1,n2,n3)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:,:,:)
  integer,intent(in)                    :: n1,n2,n3
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2,n3),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_c3d


!=========================================================================
subroutine clean_allocate_c4d(array_name,array,n1,n2,n3,n4)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:,:,:,:)
  integer,intent(in)                    :: n1,n2,n3,n4
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2,n3,n4),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_c4d


!=========================================================================
subroutine clean_deallocate_i1d(array_name,array)
  implicit none

  character(len=*),intent(in)       :: array_name
  integer,allocatable,intent(inout) :: array(:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:),DIM=1)

  mem_mb = REAL(4,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_i1d


!=========================================================================
subroutine clean_deallocate_i2d(array_name,array)
  implicit none

  character(len=*),intent(in)       :: array_name
  integer,allocatable,intent(inout) :: array(:,:)
  !=====
  real(dp)            :: mem_mb
  integer(kind=int8)  :: n1,n2
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:),DIM=1,KIND=int8)
  n2 = SIZE(array(:,:),DIM=2,KIND=int8)

  mem_mb = REAL(4,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_i2d


!=========================================================================
subroutine clean_deallocate_1d(array_name,array)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(dp),allocatable,intent(inout) :: array(:)
  !=====
  real(dp)            :: mem_mb
  integer(kind=int8)  :: n1
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:),KIND=int8)

  mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_1d


!=========================================================================
subroutine clean_deallocate_2d(array_name,array,verbose)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:)
  integer,optional,intent(in)        :: verbose
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:),DIM=1)
  n2 = SIZE(array(:,:),DIM=2)

  mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb
  
  if(PRESENT(verbose)) then
   if(verbose/=-1)call write_memory_deallocate(array_name,mem_mb)
  else
   call write_memory_deallocate(array_name,mem_mb)
  endif

end subroutine clean_deallocate_2d


!=========================================================================
subroutine clean_deallocate_3d(array_name,array,verbose)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:,:)
  integer,optional,intent(in)        :: verbose
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2,n3
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:,:),DIM=1)
  n2 = SIZE(array(:,:,:),DIM=2)
  n3 = SIZE(array(:,:,:),DIM=3)

  mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  if(PRESENT(verbose)) then
   if(verbose/=-1)call write_memory_deallocate(array_name,mem_mb)
  else
   call write_memory_deallocate(array_name,mem_mb)
  endif

end subroutine clean_deallocate_3d


!=========================================================================
subroutine clean_deallocate_4d(array_name,array,verbose)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(dp),allocatable,intent(inout) :: array(:,:,:,:)
  integer,optional,intent(in)        :: verbose
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2,n3,n4
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:,:,:),DIM=1)
  n2 = SIZE(array(:,:,:,:),DIM=2)
  n3 = SIZE(array(:,:,:,:),DIM=3)
  n4 = SIZE(array(:,:,:,:),DIM=4)

  mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp)/ 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  if(PRESENT(verbose)) then
   if(verbose/=-1)call write_memory_deallocate(array_name,mem_mb)
  else
   call write_memory_deallocate(array_name,mem_mb)
  endif

end subroutine clean_deallocate_4d


!=========================================================================
subroutine clean_deallocate_c1d(array_name,array)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:))

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_c1d


!=========================================================================
subroutine clean_deallocate_c2d(array_name,array)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:,:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:),DIM=1)
  n2 = SIZE(array(:,:),DIM=2)

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_c2d


!=========================================================================
subroutine clean_deallocate_c3d(array_name,array)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:,:,:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2,n3
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:,:),DIM=1)
  n2 = SIZE(array(:,:,:),DIM=2)
  n3 = SIZE(array(:,:,:),DIM=3)

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_c3d


!=========================================================================
subroutine clean_deallocate_c4d(array_name,array)
  implicit none

  character(len=*),intent(in)           :: array_name
  complex(dp),allocatable,intent(inout) :: array(:,:,:,:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2,n3,n4
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:,:,:),DIM=1)
  n2 = SIZE(array(:,:,:,:),DIM=2)
  n3 = SIZE(array(:,:,:,:),DIM=3)
  n4 = SIZE(array(:,:,:,:),DIM=4)

  mem_mb = 2._dp*REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp)/ 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_c4d


!=========================================================================
subroutine clean_allocate_s1d(array_name,array,n1)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:)
  integer,intent(in)                 :: n1
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(sp,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_s1d


!=========================================================================
subroutine clean_allocate_s2d(array_name,array,n1,n2)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:,:)
  integer,intent(in)                 :: n1,n2
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_s2d


!=========================================================================
subroutine clean_allocate_s3d(array_name,array,n1,n2,n3)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:,:,:)
  integer,intent(in)                 :: n1,n2,n3
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2,n3),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_s3d


!=========================================================================
subroutine clean_allocate_s4d(array_name,array,n1,n2,n3,n4)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:,:,:,:)
  integer,intent(in)                 :: n1,n2,n3,n4
  !=====
  integer             :: info
  real(dp)            :: mem_mb
  !=====

  if( ALLOCATED(array) ) then
    call die('clean_allocate: Cannot allocate. This array is already allocated -> '//TRIM(array_name))
  endif

  mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp) / 1024._dp**2

  ! The allocation itself
  allocate(array(n1,n2,n3,n4),stat=info)

  if(info/=0) then
    write(stdout,*) 'failure'
    call die('clean_allocate: Not enough memory. Buy a bigger computer')
  endif


  total_memory = total_memory + mem_mb
  peak_memory = MAX(peak_memory,total_memory)

  call write_memory_allocate(array_name,mem_mb)

end subroutine clean_allocate_s4d


!=========================================================================
subroutine clean_deallocate_s1d(array_name,array)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:))

  mem_mb = REAL(sp,dp) * REAL(n1,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_s1d


!=========================================================================
subroutine clean_deallocate_s2d(array_name,array)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:,:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:),DIM=1)
  n2 = SIZE(array(:,:),DIM=2)

  mem_mb = REAL(sp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_s2d


!=========================================================================
subroutine clean_deallocate_s3d(array_name,array)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:,:,:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2,n3
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:,:),DIM=1)
  n2 = SIZE(array(:,:,:),DIM=2)
  n3 = SIZE(array(:,:,:),DIM=3)

  mem_mb = REAL(sp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_s3d


!=========================================================================
subroutine clean_deallocate_s4d(array_name,array)
  implicit none

  character(len=*),intent(in)       :: array_name
  real(sp),allocatable,intent(inout) :: array(:,:,:,:)
  !=====
  real(dp)            :: mem_mb
  integer             :: n1,n2,n3,n4
  !=====

  if( .NOT. ALLOCATED(array) ) return

  n1 = SIZE(array(:,:,:,:),DIM=1)
  n2 = SIZE(array(:,:,:,:),DIM=2)
  n3 = SIZE(array(:,:,:,:),DIM=3)
  n4 = SIZE(array(:,:,:,:),DIM=4)

  mem_mb = REAL(sp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp)/ 1024._dp**2

  ! The allocation itself
  deallocate(array)

  total_memory = total_memory - mem_mb

  call write_memory_deallocate(array_name,mem_mb)

end subroutine clean_deallocate_s4d


!=========================================================================
subroutine write_memory_allocate(array_name,mem_mb)
  implicit none
  !=====
  character(len=*),intent(in)       :: array_name
  real(dp),intent(in)               :: mem_mb
  !=====

  write(stdout,'(1x,a,a50,a,sp,f8.3,s,a,a,f8.3,a)')   &
        '  Allocate: ',TRIM(array_name),        &
            ',  Mem: ', mem_mb/1024.0_dp      ,' (Gb)',   &
      ',  Total Mem: ', total_memory/1024.0_dp,' (Gb)'

end subroutine write_memory_allocate


!=========================================================================
subroutine write_memory_deallocate(array_name,mem_mb)
  implicit none
  !=====
  character(len=*),intent(in)       :: array_name
  real(dp),intent(in)               :: mem_mb
  !=====

  write(stdout,'(1x,a,a50,a,sp,f8.3,s,a,a,f8.3,a)')   &
        'Deallocate: ',TRIM(array_name),        &
            ',  Mem: ',-mem_mb/1024.0_dp      ,' (Gb)',   &
      ',  Total Mem: ', total_memory/1024.0_dp,' (Gb)'

end subroutine write_memory_deallocate


end module m_memory
!=========================================================================
