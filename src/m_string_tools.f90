!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! some unsorted basic operation on strings and lists
!
!=========================================================================
module m_string_tools
 use m_definitions
 use m_warning,only: die

 interface append_to_list
   module procedure append_to_list_i
   module procedure append_to_list_r
 end interface


contains


!=========================================================================
pure function capitalize(str)
 implicit none
 character(*), intent(in) :: str
 character(LEN(str))      :: capitalize
!=====
 character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
 character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
 integer :: ic, ii
!=====

 capitalize = str
 do ii=1,LEN_TRIM(str)
   ic = INDEX(low,str(ii:ii))
   if (ic > 0) capitalize(ii:ii) = cap(ic:ic)
 enddo

end function capitalize


!=========================================================================
function orbital_momentum_number(amc)
 character(len=1),intent(in) :: amc
 integer :: orbital_momentum_number
!=====

 select case(capitalize(amc))
 case('S')
   orbital_momentum_number = 0
 case('P')
   orbital_momentum_number = 1
 case('D')
   orbital_momentum_number = 2
 case('F')
   orbital_momentum_number = 3
 case('G')
   orbital_momentum_number = 4
 case('H')
   orbital_momentum_number = 5
 case('I')
   orbital_momentum_number = 6
 case('K')
   orbital_momentum_number = 7
 case default
   write(stdout,*) amc,capitalize(amc)
   call die('orbital_momentum_number: keyword unknown')
 end select


end function orbital_momentum_number


!=========================================================================
pure function orbital_momentum_name(am)
 integer,intent(in) :: am
 character(len=1) :: orbital_momentum_name
!=====

 select case(am)
 case(0)
   orbital_momentum_name='s'
 case(1)
   orbital_momentum_name='p'
 case(2)
   orbital_momentum_name='d'
 case(3)
   orbital_momentum_name='f'
 case(4)
   orbital_momentum_name='g'
 case(5)
   orbital_momentum_name='h'
 case(6)
   orbital_momentum_name='i'
 case(7)
   orbital_momentum_name='k'
 case(8)
   orbital_momentum_name='l'
 case(9)
   orbital_momentum_name='m'
 case(10)
   orbital_momentum_name='n'
 case(11)
   orbital_momentum_name='o'
 case(12)
   orbital_momentum_name='q'
 case default
   orbital_momentum_name='x'
 end select

end function orbital_momentum_name


!=========================================================================
subroutine append_to_list_i(new_element,list)
 implicit none

 integer,intent(in)                :: new_element
 integer,allocatable,intent(inout) :: list(:)
!=====
 integer :: nsize
 integer,allocatable :: list_old(:)
!=====

 if( ALLOCATED(list) ) then
   nsize = SIZE(list)
 else
   nsize = 0
 endif

 ! Copy old list and free the list
 allocate(list_old(nsize))
 if( nsize > 0 ) then
   list_old(1:nsize) =list(1:nsize)
   deallocate(list)
 endif

 allocate(list(nsize+1))
 if( nsize > 0 ) then
   list(1:nsize) = list_old(1:nsize)
 endif
 list(nsize+1) = new_element


end subroutine append_to_list_i


!=========================================================================
subroutine append_to_list_r(new_element,list)
 implicit none

 real(dp),intent(in)                :: new_element
 real(dp),allocatable,intent(inout) :: list(:)
!=====
 integer :: nsize
 real(dp),allocatable :: list_old(:)
!=====

 if( ALLOCATED(list) ) then
   nsize = SIZE(list)
 else
   nsize = 0
 endif

 ! Copy old list and free the list
 allocate(list_old(nsize))
 if( nsize > 0 ) then
   list_old(1:nsize) =list(1:nsize)
   deallocate(list)
 endif

 allocate(list(nsize+1))
 if( nsize > 0 ) then
   list(1:nsize) = list_old(1:nsize)
 endif
 list(nsize+1) = new_element


end subroutine append_to_list_r


!=========================================================================
pure function get_number_of_elements(string) result(num)
 implicit none
 character(len=*),intent(in)  :: string
 integer                      :: num
!=====
 integer   :: ip,pos
!=====

 pos = 1
 num = 0

 do
   ip = VERIFY(string(pos:),' ')  !-- Find next non-blank
   if( ip == 0 ) exit             !-- No word found
   num = num + 1                 !-- Found something
   pos = pos + ip - 1             !-- Move to start of the word
   ip = SCAN(string(pos:),' ')    !-- Find next blank
   if( ip == 0 ) exit             !-- No blank found
   pos = pos + ip - 1             !-- Move to the blank
 end do

end function get_number_of_elements


!=========================================================================
subroutine string_to_integers(string_in,iarray)
 implicit none

 character(len=*),intent(in) :: string_in
 integer,intent(inout)       :: iarray(:)
!=====
 character(LEN(string_in)) :: string
 integer                   :: ilen,inextblank,ii
!=====

 string = string_in

 ilen = LEN(TRIM(string))
 ii = 0
 do while( ilen > 0 )
   string = ADJUSTL(string)
   inextblank = INDEX(string,' ')
   ii = ii + 1
   if( ii > SIZE(iarray) ) exit
   read(string(1:inextblank-1),*) iarray(ii)
   string = string(inextblank+1:)
   ilen = LEN(TRIM(string))
 enddo

end subroutine string_to_integers


!=========================================================================
subroutine string_to_reals(string_in,rarray)
 implicit none

 character(len=*),intent(in) :: string_in
 real(dp),intent(inout)      :: rarray(:)
!=====
 character(LEN(string_in)) :: string
 integer            :: ilen,inextblank,ii
!=====

 string = string_in

 ilen = LEN(TRIM(string))
 ii = 0
 do while( ilen > 0 )
   string = ADJUSTL(string)
   inextblank = INDEX(string,' ')
   ii = ii + 1
   if( ii > SIZE(rarray) ) exit
   read(string(1:inextblank-1),*) rarray(ii)
   string = string(inextblank+1:)
   ilen = LEN(TRIM(string))
 enddo

end subroutine string_to_reals


!=========================================================================
function get_number_of_lines(filename) result(nlines)
 implicit none
 character(len=*),intent(in)  :: filename
 character(len=100)           :: cur_string
 integer                      :: nlines
!=====
 integer   :: file_unit,io
!=====

 nlines=0
 open(newunit=file_unit, file = filename)
 do
   read(file_unit,'(A)',iostat=io)cur_string
   if ( io /= 0 ) exit
   nlines=nlines+1
 end do
 close(file_unit)

end function get_number_of_lines


!=========================================================================
end module m_string_tools


!=========================================================================
