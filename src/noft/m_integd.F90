!!****m* DoNOF/m_integd
!! NAME
!! basic integrals variables
!!
!! FUNCTION
!! This module contains definitions of the arrays containing all integrals used by NOFT module.
!!
!! COPYRIGHT
!! This file is distributed under the terms of the
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

module m_integd

 use m_vars
 implicit none

!!****t* m_integd/integ_t
!! NAME
!! integ_t
!!
!! FUNCTION
!! Datatype storing integral arrays needed
!!
!! SOURCE

 type,public :: integ_t

 integer::iERItyp=0              ! Type of ERI notation to use DoNOF=0, Physicist = 1, Chemist = 2   
! arrays 
 real(dp),allocatable,dimension(:)::ERI_J,ERI_K
 real(dp),allocatable,dimension(:,:)::hCORE,Overlap
 real(dp),allocatable,dimension(:,:,:,:)::ERImol

 contains 
   procedure :: free => integ_free
   ! Destructor.

   procedure :: eritoeriJK => eri_to_eriJK  
   ! ERImol to ERI_J and ERI_K.

   procedure :: print_int => print_ints 
   ! Print hCORE and ERImol integrals in their current status

 end type integ_t

 public :: integ_init     ! Main creation method.
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/integ_init
!! NAME
!! integ_init
!!
!! FUNCTION
!!  Initialize the data type integ_t 
!!
!! INPUTS
!! NBF_tot=Number of total orbitals
!! NBF_occ=Number of orbitals that are occupied
!! Overlap_in=S_ao overlap matrix in atomic orbs.
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine integ_init(INTEGd,NBF_tot,NBF_occ,iERItyp_in,Overlap_in)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_tot,NBF_occ,iERItyp_in
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(NBF_tot,NBF_tot),intent(in)::Overlap_in
!Local variables ------------------------------
!scalars
 integer::NBF_ldiag
 real(dp)::totMEM
!arrays
!************************************************************************

 INTEGd%iERItyp=iERItyp_in
 NBF_ldiag=NBF_occ*(NBF_occ+1)/2
 ! Calculate memory needed
 totMEM=2*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_tot*NBF_tot*NBF_tot
 totMEM=8*totMEM       ! Bytes
 totMEM=totMEM*1.0d-6  ! Bytes to Mb  
 if(totMEM>1.0d3) then     ! Mb to Gb
  write(*,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM*1.0d-3,' Gb'
 elseif(totMEM<1.0d0) then ! Mb to Kb
  write(*,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM*1.0d3,' Kb'
 else                      ! Mb
  write(*,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM,' Mb'
 endif
 ! Allocate arrays
 allocate(INTEGd%ERI_J(NBF_ldiag),INTEGd%ERI_K(NBF_ldiag))
 allocate(INTEGd%hCORE(NBF_tot,NBF_tot),INTEGd%Overlap(NBF_tot,NBF_tot))
 allocate(INTEGd%ERImol(NBF_tot,NBF_tot,NBF_tot,NBF_tot))
 INTEGd%Overlap=Overlap_in

end subroutine integ_init
!!***

!!***
!!****f* DoNOF/integ_free
!! NAME
!! integ_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type integ_t 
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

subroutine integ_free(INTEGd)
!Arguments ------------------------------------
!scalars
 class(integ_t),intent(inout)::INTEGd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 deallocate(INTEGd%hCORE) 
 deallocate(INTEGd%Overlap) 
 deallocate(INTEGd%ERImol) 
 deallocate(INTEGd%ERI_J) 
 deallocate(INTEGd%ERI_K) 

end subroutine integ_free
!!***

!!***
!!****f* DoNOF/eri_to_eriJK
!! NAME
!! eri_to_eriJK
!!
!! FUNCTION
!!  Get ERI_J, ERI_K from ERI
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

subroutine eri_to_eriJK(INTEGd,NBF_occ)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_occ
 class(integ_t),intent(inout)::INTEGd
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2
!arrays
!************************************************************************

 iorb2=1
 do iorb=1,NBF_occ
  do iorb1=1,iorb
   if(INTEGd%iERItyp==0) then
    INTEGd%ERI_J(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb1,iorb) ! J in DoNOF {ij|ji}
    INTEGd%ERI_K(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb,iorb1) ! K in DoNOF {ij|ij}
   elseif(INTEGd%iERItyp==1) then
    INTEGd%ERI_J(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb,iorb1) ! J in <ij|ij>
    INTEGd%ERI_K(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb1,iorb) ! K in <ij|ji>
   elseif(INTEGd%iERItyp==2) then
    INTEGd%ERI_J(iorb2)=INTEGd%ERImol(iorb,iorb,iorb1,iorb1) ! J in (ii|jj)
    INTEGd%ERI_K(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb1,iorb) ! K in (ij|ji)
   else 
    ! Nth
   endif
   iorb2=iorb2+1
  enddo
 enddo

end subroutine eri_to_eriJK
!!***

!!***
!!****f* DoNOF/print_ints
!! NAME
!! print_ints
!!
!! FUNCTION
!!  Print hCORE and ERImol integrals in their current status 
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

subroutine print_ints(INTEGd,NBF_tot)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_tot
 class(integ_t),intent(in)::INTEGd
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iunit=312
 real(dp)::tol8=1.0d-8
!arrays
!************************************************************************
 
 ! Print ERImol
 open(unit=iunit,form='unformatted',file='ERImol')
 do iorb=1,NBF_tot
  do iorb1=1,NBF_tot
   do iorb2=1,NBF_tot
    do iorb3=1,NBF_tot
     if(dabs(INTEGd%ERImol(iorb,iorb1,iorb2,iorb3))>tol8) then
      if(INTEGd%iERItyp==0) then
       write(iunit) iorb,iorb1,iorb3,iorb2,INTEGd%ERImol(iorb,iorb1,iorb2,iorb3) ! DoNOF {ij|lk} 
      elseif(INTEGd%iERItyp==1) then
       write(iunit) iorb,iorb1,iorb2,iorb3,INTEGd%ERImol(iorb,iorb1,iorb2,iorb3) ! <ij|kl>
      elseif(INTEGd%iERItyp==2) then
       write(iunit) iorb,iorb2,iorb1,iorb3,INTEGd%ERImol(iorb,iorb1,iorb2,iorb3) ! (ik|jl)
      else
       ! Nth
      endif
     endif
    enddo
   enddo
  enddo
 enddo 
 write(iunit) 0,0,0,0,0.0d0
 write(iunit) 0,0,0,0,0.0d0
 close(iunit)

 ! Print hCORE
 open(unit=iunit,form='unformatted',file='hCORE')
 do iorb=1,NBF_tot
  do iorb1=1,NBF_tot
   if(dabs(INTEGd%hCORE(iorb,iorb1))>tol8) then
    write(iunit) iorb,iorb1,INTEGd%hCORE(iorb,iorb1)
   endif
  enddo
 enddo
 write(iunit) 0,0,0.0d0
 close(iunit)

end subroutine print_ints
!!***

end module m_integd
!!***
