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

 use m_nofoutput
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

 integer::iERItyp=0              ! Type of ERI notation to use DoNOF=0, Physicist=1, Chemist=2, Vectorial(Phys)=-1   
 integer::NBF_jkl=0              ! Size of the basis for the <:j|kl> terms
 integer::NBF2,NBF3,NBF4         ! Sizes used in the vectorial allocation of ERIs
! arrays 
 real(dp),allocatable,dimension(:)::ERI_J,ERI_K,ERI_L
 real(dp),allocatable,dimension(:)::ERImolv
 real(dp),allocatable,dimension(:,:)::hCORE,Overlap
 real(dp),allocatable,dimension(:,:,:,:)::ERImol

 contains 
   procedure :: free => integ_free
   ! Destructor.

   procedure :: eritoeriJKL => eri_to_eriJKL
   ! ERImol to ERI_J, ERI_K, and ERI_L.

   procedure :: print_int => print_ints 
   ! Print hCORE and ERImol integrals in their current status

   procedure :: print_dump => print_fcidump
   ! Print FCIDUMP integrals file in their current status

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

subroutine integ_init(INTEGd,NBF_tot,NBF_occ,iERItyp_in,Overlap_in,lowmemERI)
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::lowmemERI
 integer,intent(in)::NBF_tot,NBF_occ,iERItyp_in
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(NBF_tot,NBF_tot),intent(in)::Overlap_in
!Local variables ------------------------------
!scalars
 integer::NBF_ldiag
 real(dp)::totMEM
!arrays
 character(len=200)::msg
!************************************************************************

 INTEGd%iERItyp=iERItyp_in
 INTEGd%NBF_jkl=NBF_tot
 NBF_ldiag=NBF_occ*(NBF_occ+1)/2
 ! Calculate memory needed
 if(present(lowmemERI)) then
  if(lowmemERI) then
   totMEM=3*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_occ*NBF_occ*NBF_occ
   INTEGd%NBF_jkl=NBF_occ
  else
   totMEM=3*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_tot*NBF_tot*NBF_tot
  endif
 else
  totMEM=3*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_tot*NBF_tot*NBF_tot
 endif
 totMEM=8*totMEM       ! Bytes
 totMEM=totMEM*tol6    ! Bytes to Mb  
 if(totMEM>thousand) then  ! Mb to Gb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM*tol3,' Gb'
 elseif(totMEM<one) then   ! Mb to Kb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM*thousand,' Kb'
 else                      ! Mb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM,' Mb'
 endif
 call write_output(msg)
 ! Allocate arrays
 allocate(INTEGd%ERI_J(NBF_ldiag),INTEGd%ERI_K(NBF_ldiag),INTEGd%ERI_L(NBF_ldiag))
 allocate(INTEGd%hCORE(NBF_tot,NBF_tot),INTEGd%Overlap(NBF_tot,NBF_tot))
 if(INTEGd%iERItyp/=-1) then ! Allocate ERImol all types except vectorial
  allocate(INTEGd%ERImol(NBF_tot,INTEGd%NBF_jkl,INTEGd%NBF_jkl,INTEGd%NBF_jkl))
  allocate(INTEGd%ERImolv(1))
 else                        ! Allocate vectorial ERI case
  INTEGd%NBF2=NBF_tot
  INTEGd%NBF3=INTEGd%NBF2*INTEGd%NBF_jkl
  INTEGd%NBF4=INTEGd%NBF3*INTEGd%NBF_jkl
  allocate(INTEGd%ERImol(1,1,1,1))
  allocate(INTEGd%ERImolv(NBF_tot*INTEGd%NBF_jkl*INTEGd%NBF_jkl*INTEGd%NBF_jkl))
 endif
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
 deallocate(INTEGd%ERImolv) 
 deallocate(INTEGd%ERI_J) 
 deallocate(INTEGd%ERI_K) 
 deallocate(INTEGd%ERI_L) 

end subroutine integ_free
!!***

!!***
!!****f* DoNOF/eri_to_eriJKL
!! NAME
!! eri_to_eriJKL
!!
!! FUNCTION
!!  Get ERI_J, ERI_K, and ERI_L from ERI
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

subroutine eri_to_eriJKL(INTEGd,NBF_occ)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_occ
 class(integ_t),intent(inout)::INTEGd
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorbm1,iorb1m1
!arrays
!************************************************************************

 iorb2=1
 do iorb=1,NBF_occ
  do iorb1=1,iorb
   if(INTEGd%iERItyp==0) then
    INTEGd%ERI_J(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb1,iorb) ! J in DoNOF {ij|ji}
    INTEGd%ERI_K(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb,iorb1) ! K in DoNOF {ij|ij}
    INTEGd%ERI_L(iorb2)=INTEGd%ERImol(iorb,iorb,iorb1,iorb1) ! L in DoNOF {ii|jj}
   elseif(INTEGd%iERItyp==1) then
    INTEGd%ERI_J(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb,iorb1) ! J in <ij|ij>
    INTEGd%ERI_K(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb1,iorb) ! K in <ij|ji>
    INTEGd%ERI_L(iorb2)=INTEGd%ERImol(iorb,iorb,iorb1,iorb1) ! L in <ii|jj>
   elseif(INTEGd%iERItyp==2) then
    INTEGd%ERI_J(iorb2)=INTEGd%ERImol(iorb,iorb,iorb1,iorb1) ! J in (ii|jj)
    INTEGd%ERI_K(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb1,iorb) ! K in (ij|ji)
    INTEGd%ERI_L(iorb2)=INTEGd%ERImol(iorb,iorb1,iorb,iorb1) ! L in (ij|ij)
   elseif(INTEGd%iERItyp==-1) then
    iorbm1=iorb-1
    iorb1m1=iorb1-1
    INTEGd%ERI_J(iorb2)=INTEGd%ERImolv(iorbm1+iorb1m1*INTEGd%NBF2+iorbm1*INTEGd%NBF3+iorb1m1*INTEGd%NBF4+1) ! J in <i+j*NBF2+i*NBF3+j*NBF4> 
    INTEGd%ERI_K(iorb2)=INTEGd%ERImolv(iorbm1+iorb1m1*INTEGd%NBF2+iorb1m1*INTEGd%NBF3+iorbm1*INTEGd%NBF4+1) ! K in <i+j*NBF2+j*NBF3+i*NBF4> 
    INTEGd%ERI_L(iorb2)=INTEGd%ERImolv(iorbm1+iorbm1*INTEGd%NBF2+iorb1m1*INTEGd%NBF3+iorb1m1*INTEGd%NBF4+1) ! L in <i+i*NBF2+j*NBF3+j*NBF4> 
   else 
    ! Nth
   endif
   iorb2=iorb2+1
  enddo
 enddo

end subroutine eri_to_eriJKL
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

subroutine print_ints(INTEGd)
!Arguments ------------------------------------
!scalars
 class(integ_t),intent(in)::INTEGd
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iunit=312
 integer::iorbm1,iorb1m1,iorb2m1,iorb3m1
 real(dp)::ERIval
!arrays
!************************************************************************
 
 ! Print ERImol
 open(unit=iunit,form='unformatted',file='ERImol')
 do iorb=1,INTEGd%NBF_jkl
  do iorb1=1,INTEGd%NBF_jkl
   do iorb2=1,INTEGd%NBF_jkl
    do iorb3=1,INTEGd%NBF_jkl
     if(INTEGd%iERItyp/=-1) then
      ERIval=INTEGd%ERImol(iorb,iorb1,iorb2,iorb3)
      if(dabs(ERIval)>tol8) then
       if(INTEGd%iERItyp==0) then
        write(iunit) iorb,iorb1,iorb3,iorb2,ERIval ! DoNOF {ij|lk} 
       elseif(INTEGd%iERItyp==1) then
        write(iunit) iorb,iorb1,iorb2,iorb3,ERIval ! <ij|kl>
       elseif(INTEGd%iERItyp==2) then
        write(iunit) iorb,iorb2,iorb1,iorb3,ERIval ! (ik|jl)
       else
        ! Nth
       endif
      endif
     else
      iorbm1=iorb-1
      iorb1m1=iorb1-1
      iorb2m1=iorb2-1
      iorb3m1=iorb3-1
      ERIval=INTEGd%ERImolv(iorbm1+iorb1m1*INTEGd%NBF2+iorb2m1*INTEGd%NBF3+iorb3m1*INTEGd%NBF4+1) ! <i+j*NBF2+k*NBF3+l*NBF4> 
      if(dabs(ERIval)>tol8) then
       write(iunit) iorb,iorb1,iorb2,iorb3,ERIval
      endif
     endif
    enddo
   enddo
  enddo
 enddo 
 write(iunit) 0,0,0,0,zero
 write(iunit) 0,0,0,0,zero
 close(iunit)

 ! Print hCORE
 open(unit=iunit,form='unformatted',file='hCORE')
 do iorb=1,INTEGd%NBF_jkl
  do iorb1=1,INTEGd%NBF_jkl
   if(dabs(INTEGd%hCORE(iorb,iorb1))>tol8) then
    write(iunit) iorb,iorb1,INTEGd%hCORE(iorb,iorb1)
   endif
  enddo
 enddo
 write(iunit) 0,0,zero
 close(iunit)

end subroutine print_ints
!!***

!!***
!!****f* DoNOF/print_fcidump
!! NAME
!! print_fcidump
!!
!! FUNCTION
!!  Print FCIDUMP integrals in their current status 
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

subroutine print_fcidump(INTEGd,Nel,Vnn)
!Arguments ------------------------------------
 integer,intent(in)::Nel
 real(dp),intent(in)::Vnn
!scalars
 class(integ_t),intent(in)::INTEGd
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iunit=312
 integer::iorbm1,iorb1m1,iorb2m1,iorb3m1
 real(dp)::ERIval
!arrays
!************************************************************************
 
 ! Print FCIDUMP
 open(unit=iunit,file='FCIDUMP')
 write(iunit,'(a12,i4,a8,i4)') '$ FCI NORB =',INTEGd%NBF_jkl,' NELEC =',Nel
 write(iunit,'(a10,i10)') '  MEMORY =',1000000
 write(iunit,'(a1)') '$ '
 do iorb=1,INTEGd%NBF_jkl
  do iorb1=1,iorb
   do iorb2=iorb,INTEGd%NBF_jkl
    do iorb3=1,iorb2
     if(INTEGd%iERItyp/=-1) then
      if(INTEGd%iERItyp==0) then
       ERIval=INTEGd%ERImol(iorb,iorb2,iorb3,iorb1) ! DoNOF {ij|lk}
      elseif(INTEGd%iERItyp==1) then
       ERIval=INTEGd%ERImol(iorb,iorb2,iorb1,iorb3) ! <ij|kl>
      elseif(INTEGd%iERItyp==2) then
       ERIval=INTEGd%ERImol(iorb,iorb1,iorb2,iorb3) ! (ik|jl)
      else
       ERIval=zero                                  ! Nth
      endif
     else
      iorbm1=iorb-1
      iorb1m1=iorb1-1
      iorb2m1=iorb2-1
      iorb3m1=iorb3-1
      ERIval=INTEGd%ERImolv(iorbm1+iorb2m1*INTEGd%NBF2+iorb1m1*INTEGd%NBF3+iorb3m1*INTEGd%NBF4+1) ! <i+j*NBF2+k*NBF3+l*NBF4> 
     endif
     if(dabs(ERIval)>tol8) write(iunit,'(f15.8,4i4)') ERIval,iorb,iorb1,iorb2,iorb3
    enddo
   enddo
  enddo
 enddo 
 do iorb=1,INTEGd%NBF_jkl
  do iorb1=1,iorb
   if(dabs(INTEGd%hCORE(iorb,iorb1))>tol8) then
    write(iunit,'(f15.8,4i4)') INTEGd%hCORE(iorb,iorb1),iorb,iorb1,0,0
   endif
  enddo
 enddo
 write(iunit,'(f15.8,4i4)') Vnn,0,0,0,0
 close(iunit)

end subroutine print_fcidump
!!***

end module m_integd
!!***
