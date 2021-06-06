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

 public::init_output,write_output

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

end module m_nofoutput
!!***
