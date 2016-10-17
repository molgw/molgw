!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods and data for Effective Core Potentials (ECP)
!
!=========================================================================
module m_ecp
 use m_definitions
 use m_tools,only: capitalize,append_to_list,orbital_momentum_name,orbital_momentum_number
 use m_warning,only: die,issue_warning
 use m_elements


 integer,protected                :: nelement_ecp
 integer,protected,allocatable    :: element_ecp(:)

 type effective_core_potential
   integer              :: nelec
   integer              :: necp
   integer,allocatable  :: lk(:)
   integer,allocatable  :: nk(:)
   real(dp),allocatable :: dk(:)
   real(dp),allocatable :: zetak(:)
 end type

 type(effective_core_potential),allocatable :: ecp(:)



contains


subroutine init_ecp(ecp_elements,ecp_path,ecp_name)
 implicit none

 character(len=*),intent(in) :: ecp_elements
 character(len=*),intent(in) :: ecp_path
 character(len=*),intent(in) :: ecp_name
!=====
 character(len=132) :: string,ecp_filename
 character(len=2)   :: element
 integer :: ecpfile
 integer :: ilen,inextblank,ielement_ecp,iecp
 logical :: file_exists
!=====

 !
 ! First parse the ecp_elements line
 !
 ilen = LEN(TRIM(ecp_elements))

 ! ecp_elements is empty, no ECP needs to be setup
 if( ilen == 0 ) return

 string = ecp_elements
 write(stdout,'(/,1x,a)') 'Reading ECP element list'

 do while( ilen > 0 )
   string = ADJUSTL(string)
   inextblank = INDEX(string,' ')

   call append_to_list(element_number(string(1:inextblank-1)),element_ecp)

   string = string(inextblank+1:)
   ilen = LEN(TRIM(string))

 enddo 

 nelement_ecp = SIZE(element_ecp)
 allocate(ecp(nelement_ecp))


 !
 ! Second, read the ECP parameters from ECP file
 !
 do ielement_ecp=1,nelement_ecp
   element = element_name(REAL(element_ecp(ielement_ecp),dp))
   write(stdout,'(1x,a,a)') 'ECP for element: ',element

   ecp_filename = TRIM(ecp_path)//TRIM(element)//'_'//TRIM(ecp_name)
   inquire(file=TRIM(ecp_filename),exist=file_exists)
   if( .NOT. file_exists ) then
     write(stdout,*) 'Looking for file: ',ecp_filename
     call die('File not found')
   endif
   open(newunit=ecpfile,file=TRIM(ecp_filename),status='old',action='read')

   call read_ecp_file(ecpfile,element,ecp(ielement_ecp))
   close(ecpfile)

   write(stdout,'(6x,a,i3)') 'Core electrons ',ecp(ielement_ecp)%nelec
   write(stdout,'(6x,a)') 'l_k  n_k       zeta_k          d_k  '

   do iecp=1,ecp(ielement_ecp)%necp
     write(stdout,'(6x,a,2x,i3,2(2x,f14.6))') &
                         orbital_momentum_name(ecp(ielement_ecp)%lk(iecp)), &
                         ecp(ielement_ecp)%nk(iecp), &
                         ecp(ielement_ecp)%zetak(iecp), &
                         ecp(ielement_ecp)%dk(iecp) 
   enddo

 enddo



end subroutine init_ecp


!=========================================================================
subroutine read_ecp_file(ecpunit,element,ecpi)
 use ISO_FORTRAN_ENV,only: iostat_end
 implicit none

 integer,intent(in)                           :: ecpunit
 character(len=2),intent(in)                  :: element
 type(effective_core_potential),intent(inout) :: ecpi
!=====
 integer :: iline,i1,i2,istat
 integer :: nelec
 character(len=132) :: line,amc
 character(len=2)   :: element_read
 integer  :: read_n
 real(dp) :: read_zeta,read_d
 logical  :: end_of_file
!=====

 ! Reading an ECP file in NWCHEM format

 end_of_file = .FALSE.
 iline = 0
 line='_____'   ! Five underscores '_____' means 'advance to next line'
 do while(.NOT. end_of_file)
   iline = iline + 1
   if( line(1:5) == '_____' ) then 
     read(ecpunit,'(a)',iostat=istat) line
     if( istat == iostat_end ) then
!       write(*,*) 'End of file reached'
       end_of_file = .TRUE.
       exit
     endif
   endif
   line = ADJUSTL(line)

   ! Remove comments if any
   if( line(1:1) == '#' ) then
     line='_____'
     cycle
   endif

   ! ECP and END should not be interpreted
   if( capitalize(line(1:3)) == 'ECP' .OR. capitalize(line(1:3)) == 'END' ) then
     line='_____'
     cycle
   endif
   i1 = INDEX(line,' ')
!   write(*,*) i1,'_'//line(1:i1-1)//'_'

   if( line(1:i1-1) /= TRIM(element) .AND. capitalize(line(1:i1-1)) /= TRIM(element) ) then
     write(stdout,*) 'ECP file should only contain information about element '//TRIM(element)
     write(stdout,*) 'While '//line(1:i1-1)//' was found'
     call die('ECP file reading problem')
   endif
!   if( i1-1 > 2 ) write(*,*) 'in line ',iline,'the first symbol should be an element'

   line = ADJUSTL(line(i1+1:))

   i2 = INDEX(line,' ')
!   write(*,*) 'key  _'//line(1:i2-1)//'_'
   amc = capitalize(line(1:i2-1))
   if( amc == 'NELEC' ) then
     read(line(i2+1:),'(i10)') ecpi%nelec
     line='_____'
     cycle
   endif
   if(      amc == 'UL'  &
       .OR. amc == 'S'   &
       .OR. amc == 'P'   &
       .OR. amc == 'D'   &
       .OR. amc == 'F'   &
       .OR. amc == 'G'   ) then
     istat = 0
     do while(istat == 0)
       read(ecpunit,'(a)',iostat=istat) line
       if( istat == iostat_end ) then
         write(*,*) 'End of file reached'
         end_of_file = .TRUE.
         exit
       endif
       read(line,'(i3,e16.8,e16.8)',iostat=istat) read_n,read_zeta,read_d

       ! For the time being, only code ECP with no local potential
       if( amc == 'UL' ) then
         if( ABS(read_d) > 1.0e-4_dp ) then
           call die('For the time being, only coded for ECP with no local potential')
         else 
           cycle
         endif
       endif

       if( istat == 0 ) then
         call append_to_list(orbital_momentum_number(amc),ecpi%lk)
         call append_to_list(read_n,ecpi%nk)
         call append_to_list(read_zeta,ecpi%zetak)
         call append_to_list(read_d,ecpi%dk)
       endif
     end do
   else
     write(*,*) capitalize(line(1:i2-1)),line(1:i2-1)
     call die('problem reading ECP file')
   endif
     

 enddo


 ecpi%necp = SIZE(ecpi%nk)




end subroutine read_ecp_file


end module m_ecp
!=========================================================================
