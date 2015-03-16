!=========================================================================
#include "macros.h"
!=========================================================================
module m_memory
 use m_definitions
 use m_mpi

 real(dp),private :: total_memory=0.0_dp     ! Total memory occupied 
                                             ! by the big arrays in Mb
 real(dp),private :: peak_memory=0.0_dp      ! Max memory occupied 
                                             ! by the big arrays in Mb

 interface clean_allocate
  module procedure clean_allocate_s1d
  module procedure clean_allocate_s2d
  module procedure clean_allocate_s3d
  module procedure clean_allocate_s4d
  module procedure clean_allocate_1d
  module procedure clean_allocate_2d
  module procedure clean_allocate_3d
  module procedure clean_allocate_4d
 end interface

 interface clean_deallocate
  module procedure clean_deallocate_s1d
  module procedure clean_deallocate_s2d
  module procedure clean_deallocate_s3d
  module procedure clean_deallocate_s4d
  module procedure clean_deallocate_1d
  module procedure clean_deallocate_2d
  module procedure clean_deallocate_3d
  module procedure clean_deallocate_4d
 end interface


contains


!=========================================================================
subroutine total_memory_statement()
 implicit none
!=====

 WRITE(*,'(/,a)') ' Total memory that was not deallocated properly'
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE(*,'(/,a)') ' Maximum memory used during the run'
 if( peak_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') '  Peak memory [Mb]: ',peak_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') '  Peak memory [Gb]: ',peak_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine total_memory_statement


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

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_allocate_1d


!=========================================================================
subroutine clean_allocate_2d(array_name,array,n1,n2)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:)
 integer,intent(in)                 :: n1,n2
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1,n2),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_allocate_2d


!=========================================================================
subroutine clean_allocate_3d(array_name,array,n1,n2,n3)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:)
 integer,intent(in)                 :: n1,n2,n3
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1,n2,n3),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_allocate_3d


!=========================================================================
subroutine clean_allocate_4d(array_name,array,n1,n2,n3,n4)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:,:)
 integer,intent(in)                 :: n1,n2,n3,n4
!=====
 integer             :: info
 real(dp)            :: mem_mb
!=====

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1,n2,n3,n4),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_allocate_4d


!=========================================================================
subroutine clean_deallocate_1d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1
!=====

 n1 = SIZE(array(:))

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_1d


!=========================================================================
subroutine clean_deallocate_2d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2
!=====

 n1 = SIZE(array(:,:),DIM=1)
 n2 = SIZE(array(:,:),DIM=2)

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_2d


!=========================================================================
subroutine clean_deallocate_3d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3
!=====

 n1 = SIZE(array(:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:),DIM=3)

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_3d


!=========================================================================
subroutine clean_deallocate_4d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(dp),allocatable,intent(inout) :: array(:,:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3,n4
!=====

 n1 = SIZE(array(:,:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:,:),DIM=3)
 n4 = SIZE(array(:,:,:,:),DIM=4)

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp)/ 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_4d

!FBFBFBFBFBBF

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

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

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

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1,n2),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

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

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1,n2,n3),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

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

 WRITE_MASTER(*,'(a,x,a)') ' Allocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Memory [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 allocate(array(n1,n2,n3,n4),stat=info)

 if(info/=0) then
   WRITE_MASTER(*,*) 'failure'
   stop'Not enough memory. Buy a bigger computer'
 endif


 total_memory = total_memory + mem_mb
 peak_memory = MAX(peak_memory,total_memory)
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_allocate_s4d


!=========================================================================
subroutine clean_deallocate_s1d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1
!=====

 n1 = SIZE(array(:))

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_s1d


!=========================================================================
subroutine clean_deallocate_s2d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2
!=====

 n1 = SIZE(array(:,:),DIM=1)
 n2 = SIZE(array(:,:),DIM=2)

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) * REAL(n2,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_s2d


!=========================================================================
subroutine clean_deallocate_s3d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3
!=====

 n1 = SIZE(array(:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:),DIM=3)

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) / 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_s3d


!=========================================================================
subroutine clean_deallocate_s4d(array_name,array)
 implicit none

 character(len=*),intent(in)       :: array_name
 real(sp),allocatable,intent(inout) :: array(:,:,:,:)
!=====
 integer             :: info
 real(dp)            :: mem_mb
 integer             :: n1,n2,n3,n4
!=====

 n1 = SIZE(array(:,:,:,:),DIM=1)
 n2 = SIZE(array(:,:,:,:),DIM=2)
 n3 = SIZE(array(:,:,:,:),DIM=3)
 n4 = SIZE(array(:,:,:,:),DIM=4)

 WRITE_MASTER(*,'(a,x,a)') ' Deallocate',TRIM(array_name)
 mem_mb = REAL(dp,dp) * REAL(n1,dp) *REAL(n2,dp) * REAL(n3,dp) * REAL(n4,dp)/ 1024._dp**2
 if( mem_mb < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Mb]: ',mem_mb
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Free [Gb]: ',mem_mb / 1024._dp
 endif

 ! The allocation itself
 deallocate(array)

 total_memory = total_memory - mem_mb
 if( total_memory < 500._dp ) then
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Mb]: ',total_memory
 else
   WRITE_MASTER(*,'(a30,f9.3)') ' Total memory [Gb]: ',total_memory / 1024._dp
 endif

 WRITE_MASTER(*,*)

end subroutine clean_deallocate_s4d


end module m_memory
!=========================================================================
