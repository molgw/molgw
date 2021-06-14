!!****m* DoNOF/m_output
!! NAME
!! m_output
!!
!! FUNCTION
!!  Initialize the output file.nof, appen lines, and close it .
!!
!! COPYRIGHT
!!
!! NOTES
!!
!! SOURCE

module m_nofoutput

 implicit none
 
 character(len=100)::output_file

 public::init_output,write_output,write_header,write_footer

contains

!!***
!!****f* DoNOF/init_output
!! NAME
!! init_output
!!
!! FUNCTION
!!  Initialize the output file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine init_output(output_file_in)
!Arguments ------------------------------------
!scalars
!arrays
 character(len=100)::output_file_in
!Local variables ------------------------------
!scalars
 integer::iunit=313
!arrays

!************************************************************************

 output_file=trim(output_file_in)
 open(unit=iunit,form='formatted',file=output_file)
 
end subroutine init_output
!!***

!!***
!!****f* DoNOF/write_output
!! NAME
!! write_output
!!
!! FUNCTION
!!  Write string to output file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine write_output(msg)
!Arguments ------------------------------------
!scalars
!arrays
 character(len=200),intent(inout)::msg
!Local variables ------------------------------
!scalars
 integer::iunit=313
!arrays

!************************************************************************

 write(iunit,'(a)') trim(msg)
 
end subroutine write_output
!!***

!!***
!!****f* DoNOF/close_output
!! NAME
!! close_output
!!
!! FUNCTION
!!  Close the output file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine close_output()
!Arguments ------------------------------------
!scalars
!arrays
!Local variables ------------------------------
!scalars
 integer::iunit=313
!arrays

!************************************************************************

 close(iunit)
  
end subroutine close_output
!!***

!!***
!!****f* DoNOF/write_header
!! NAME
!! write_header
!!
!! FUNCTION
!!  Writes the header on the output file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine write_header()
!Arguments ------------------------------------
!scalars
!arrays
!Local variables ------------------------------
!scalars
!arrays
 character(len=200)::msg
 character(8)::date
 character(10)::time
 character(5)::zone
 integer,dimension(8)::tvalues

!************************************************************************

 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' -------------------------------------------'
 call write_output(msg)
 write(msg,'(a)') ' Entering RUN-NOF module for NOFT calcs.'
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' Developed by: Dr. M. Rodriguez-Mayorga '
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') '  First version: VU Amsterdam 2021 '
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 call date_and_time(date=date,time=time,zone=zone,values=tvalues)
 write(msg,'(a)') ' Starting date and time'
 call write_output(msg)
 write(msg,'(a,a2,a,a2,a,a4,a,i2,a,i2,a,i2)') " ",date(7:8),"/",date(5:6),"/",date(1:4)," ",tvalues(5),":",&
 & tvalues(6),":",tvalues(7)
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' -------------------------------------------'
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 
end subroutine write_header
!!***

!!***
!!****f* DoNOF/write_footer
!! NAME
!! write_footer
!!
!! FUNCTION
!!  Writes the header on the output file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine write_footer()
!Arguments ------------------------------------
!scalars
!arrays
!Local variables ------------------------------
!scalars
!arrays
 character(len=200)::msg
 character(8)::date
 character(10)::time
 character(5)::zone
 integer,dimension(8)::tvalues

!************************************************************************

 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' -------------------------------------------'
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 call date_and_time(date=date,time=time,zone=zone,values=tvalues)
 write(msg,'(a)') ' Final date and time'
 call write_output(msg)
 write(msg,'(a,a2,a,a2,a,a4,a,i2,a,i2,a,i2)') " ",date(7:8),"/",date(5:6),"/",date(1:4)," ",tvalues(5),":",&
 & tvalues(6),":",tvalues(7)
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' Normal termination of RUN-NOF module.'
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') '   |        <^>        | '
 call write_output(msg)
 write(msg,'(a)') '   ||===I||(-@-)||I===|| '
 call write_output(msg)
 write(msg,'(a)') '   |        \_/        | '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') '   |                   | '
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' "Your feeble skills are no match for the '
 call write_output(msg)
 write(msg,'(a)') ' power of the dark side." Emperor Palpatine '
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a)') ' -------------------------------------------'
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)

end subroutine write_footer
!!***

end module m_nofoutput
!!***
