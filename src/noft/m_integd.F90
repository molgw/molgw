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
 use m_definitions

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

 logical::complex_ints=.false.   ! Set to true in case complex integrals are employed
 integer::irange_sep=0           ! rs-NOFT calcs. 0=no, 1=intra, 2=ex_corr
 integer::iERItyp=0              ! Type of ERI notation to use DoNOF=0, Physicist=1, Chemist=2, Vectorial(Phys)=-1   
 integer::NBF_jkl=0              ! Size of the basis for the <:j|kl> terms
 integer::NBF2,NBF3,NBF4         ! Sizes used in the vectorial allocation of ERIs
! arrays 
 real(dp),allocatable,dimension(:)::ERI_J,ERI_K,ERI_L
 real(dp),allocatable,dimension(:)::ERI_Jsr,ERI_Lsr
 real(dp),allocatable,dimension(:)::ERImolv
 real(dp),allocatable,dimension(:,:)::hCORE,Overlap
 real(dp),allocatable,dimension(:,:,:)::ERImolJsr
 real(dp),allocatable,dimension(:,:,:)::ERImolLsr
 real(dp),allocatable,dimension(:,:,:,:)::ERImol
 complex(dp),allocatable,dimension(:)::ERI_J_cmplx,ERI_K_cmplx,ERI_L_cmplx
 complex(dp),allocatable,dimension(:)::ERImolv_cmplx
 complex(dp),allocatable,dimension(:,:)::hCORE_cmplx
 complex(dp),allocatable,dimension(:,:,:,:)::ERImol_cmplx

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

subroutine integ_init(INTEGd,NBF_tot,NBF_occ,iERItyp_in,Overlap_in,complex_ints_in,irs_noft,lowmemERI)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::complex_ints_in
 logical,optional,intent(in)::lowmemERI
 integer,intent(in)::NBF_tot,NBF_occ,iERItyp_in,irs_noft
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(NBF_tot,NBF_tot),intent(in)::Overlap_in
!Local variables ------------------------------
!scalars
 integer::NBF_ldiag
 real(dp)::totMEM
!arrays
 character(len=200)::msg
!************************************************************************

 if(irs_noft/=0) INTEGd%irange_sep=irs_noft
 INTEGd%iERItyp=iERItyp_in
 INTEGd%NBF_jkl=NBF_tot
 NBF_ldiag=NBF_occ*(NBF_occ+1)/2
 INTEGd%complex_ints=complex_ints_in
 ! Calculate memory needed
 if(INTEGd%complex_ints) then
  if(present(lowmemERI)) then
   if(lowmemERI) then
    totMEM=3*2*NBF_ldiag+NBF_tot*NBF_tot+4*NBF_tot*NBF_tot+16*NBF_tot*NBF_occ*NBF_occ*NBF_occ
    INTEGd%NBF_jkl=NBF_occ
   else
    totMEM=3*2*NBF_ldiag+NBF_tot*NBF_tot+4*NBF_tot*NBF_tot+16*NBF_tot*NBF_tot*NBF_tot*NBF_tot
   endif
  else
   totMEM=3*2*NBF_ldiag+NBF_tot*NBF_tot+4*NBF_tot*NBF_tot+16*NBF_tot*NBF_tot*NBF_tot*NBF_tot
  endif
 else
  if(present(lowmemERI)) then
   if(lowmemERI) then
    totMEM=5*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_occ*NBF_occ*NBF_occ
    INTEGd%NBF_jkl=NBF_occ
    if(INTEGd%irange_sep/=0) then
     totMEM=totMEM+2*NBF_tot*NBF_occ*NBF_occ
    endif
   else
    totMEM=5*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_tot*NBF_tot*NBF_tot
    if(INTEGd%irange_sep/=0) then
     totMEM=totMEM+2*NBF_tot*NBF_tot*NBF_tot
    endif
   endif
  else
   totMEM=5*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_tot*NBF_tot*NBF_tot
   if(INTEGd%irange_sep/=0) then
    totMEM=totMEM+2*NBF_tot*NBF_tot*NBF_tot
   endif
  endif
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
 if(INTEGd%complex_ints) then
  allocate(INTEGd%ERI_J_cmplx(NBF_ldiag),INTEGd%ERI_K_cmplx(NBF_ldiag),INTEGd%ERI_L_cmplx(NBF_ldiag))
  allocate(INTEGd%hCORE_cmplx(NBF_tot,NBF_tot),INTEGd%Overlap(NBF_tot,NBF_tot))
  if(INTEGd%iERItyp/=-1) then ! Allocate ERImol all types except vectorial
   allocate(INTEGd%ERImol_cmplx(NBF_tot,INTEGd%NBF_jkl,INTEGd%NBF_jkl,INTEGd%NBF_jkl))
   allocate(INTEGd%ERImolv_cmplx(1))
  else                        ! Allocate vectorial ERI case
   INTEGd%NBF2=NBF_tot
   INTEGd%NBF3=INTEGd%NBF2*INTEGd%NBF_jkl
   INTEGd%NBF4=INTEGd%NBF3*INTEGd%NBF_jkl
   allocate(INTEGd%ERImol_cmplx(1,1,1,1))
   allocate(INTEGd%ERImolv_cmplx(NBF_tot*INTEGd%NBF_jkl*INTEGd%NBF_jkl*INTEGd%NBF_jkl))
  endif
 else
  if(INTEGd%irange_sep/=0) then
   allocate(INTEGd%ERImolJsr(NBF_tot,INTEGd%NBF_jkl,INTEGd%NBF_jkl))
   allocate(INTEGd%ERImolLsr(NBF_tot,INTEGd%NBF_jkl,INTEGd%NBF_jkl))
  else 
   allocate(INTEGd%ERImolJsr(1,1,1))
   allocate(INTEGd%ERImolLsr(1,1,1))
  endif
  allocate(INTEGd%ERI_J(NBF_ldiag),INTEGd%ERI_K(NBF_ldiag),INTEGd%ERI_L(NBF_ldiag))
  allocate(INTEGd%ERI_Jsr(NBF_ldiag),INTEGd%ERI_Lsr(NBF_ldiag))
  INTEGd%ERI_J=zero; INTEGd%ERI_K=zero; INTEGd%ERI_L=zero;
  INTEGd%ERI_Jsr=zero; INTEGd%ERI_Lsr=zero;
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

 if(INTEGd%complex_ints) then
  deallocate(INTEGd%hCORE_cmplx) 
  deallocate(INTEGd%ERImol_cmplx) 
  deallocate(INTEGd%ERImolv_cmplx) 
  deallocate(INTEGd%ERI_J_cmplx) 
  deallocate(INTEGd%ERI_K_cmplx) 
  deallocate(INTEGd%ERI_L_cmplx)
 else
  deallocate(INTEGd%hCORE) 
  deallocate(INTEGd%ERImol) 
  deallocate(INTEGd%ERImolJsr) 
  deallocate(INTEGd%ERImolLsr) 
  deallocate(INTEGd%ERImolv) 
  deallocate(INTEGd%ERI_J) 
  deallocate(INTEGd%ERI_K) 
  deallocate(INTEGd%ERI_L)
  deallocate(INTEGd%ERI_Jsr)
  deallocate(INTEGd%ERI_Lsr)
 endif 
 deallocate(INTEGd%Overlap) 

end subroutine integ_free
!!***

!!***
!!****f* DoNOF/eri_to_eriJKL
!! NAME
!! eri_to_eriJKL
!!
!!  For complex orbs. with time-reversal symmetry [i.e. states p_alpha = (p_beta)* ]: 
!!      < p_alpha q_beta | r_alpha s_beta > =  < p_alpha s_alpha | r_alpha q_beta >
!!             and      Lij integral (alpha beta) = Kij integral (alpha alpha)
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

 if(INTEGd%complex_ints) then
  iorb2=1
  do iorb=1,NBF_occ
   do iorb1=1,iorb
    if(INTEGd%iERItyp==0) then
     INTEGd%ERI_J_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb1,iorb) ! J in DoNOF {ij|ji}
     INTEGd%ERI_K_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb,iorb1) ! K in DoNOF {ij|ij}
     INTEGd%ERI_L_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb,iorb1) ! L in DoNOF {ii|jj} Lij is Kij with time-rev symm
    elseif(INTEGd%iERItyp==1) then
     INTEGd%ERI_J_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb,iorb1) ! J in <ij|ij>
     INTEGd%ERI_K_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb1,iorb) ! K in <ij|ji>
     INTEGd%ERI_L_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb1,iorb) ! L in <ii|jj> Lij is Kij with time-rev symm
    elseif(INTEGd%iERItyp==2) then
     INTEGd%ERI_J_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb,iorb1,iorb1) ! J in (ii|jj)
     INTEGd%ERI_K_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb1,iorb) ! K in (ij|ji)
     INTEGd%ERI_L_cmplx(iorb2)=INTEGd%ERImol_cmplx(iorb,iorb1,iorb1,iorb) ! L in (ij|ij) Lij is Kij with time-rev symm
    elseif(INTEGd%iERItyp==-1) then
     iorbm1=iorb-1
     iorb1m1=iorb1-1
     INTEGd%ERI_J_cmplx(iorb2)=INTEGd%ERImolv_cmplx(iorbm1+iorb1m1*INTEGd%NBF2+iorbm1*INTEGd%NBF3+iorb1m1*INTEGd%NBF4+1)!J in <i+j*NBF2+i*NBF3+j*NBF4> 
     INTEGd%ERI_K_cmplx(iorb2)=INTEGd%ERImolv_cmplx(iorbm1+iorb1m1*INTEGd%NBF2+iorb1m1*INTEGd%NBF3+iorbm1*INTEGd%NBF4+1)!K in <i+j*NBF2+j*NBF3+i*NBF4> 
     INTEGd%ERI_L_cmplx(iorb2)=INTEGd%ERImolv_cmplx(iorbm1+iorb1m1*INTEGd%NBF2+iorb1m1*INTEGd%NBF3+iorbm1*INTEGd%NBF4+1)!Lij is Kij with time-rev symm
    !INTEGd%ERI_L_cmplx(iorb2)=INTEGd%ERImolv_cmplx(iorbm1+iorbm1*INTEGd%NBF2+iorb1m1*INTEGd%NBF3+iorb1m1*INTEGd%NBF4+1)!L in <i+i*NBF2+j*NBF3+j*NBF4> 
    else 
     ! Nth
    endif
    iorb2=iorb2+1
   enddo
  enddo
 else
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
    if(INTEGd%irange_sep/=0) then
     INTEGd%ERI_Jsr(iorb2)=INTEGd%ERImolJsr(iorb,iorb1,iorb) ! Hartree <ij|ij>
     INTEGd%ERI_Lsr(iorb2)=INTEGd%ERImolLsr(iorb,iorb,iorb1) ! Hartree <ii|jj>
    endif
    iorb2=iorb2+1
   enddo
  enddo
 endif

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
 complex(dp)::ERIval_cmplx
!arrays
!************************************************************************
 
 ! Print ERImol
 open(unit=iunit,form='unformatted',file='ERImol')
 if(INTEGd%complex_ints) then
  do iorb=1,INTEGd%NBF_jkl
   do iorb1=1,INTEGd%NBF_jkl
    do iorb2=1,INTEGd%NBF_jkl
     do iorb3=1,INTEGd%NBF_jkl
      if(INTEGd%iERItyp/=-1) then
       ERIval_cmplx=INTEGd%ERImol_cmplx(iorb,iorb1,iorb2,iorb3)
       if(cdabs(ERIval_cmplx)>tol8) then
        if(INTEGd%iERItyp==0) then
         write(iunit) iorb,iorb1,iorb3,iorb2,ERIval_cmplx ! DoNOF {ij|lk} 
        elseif(INTEGd%iERItyp==1) then
         write(iunit) iorb,iorb1,iorb2,iorb3,ERIval_cmplx ! <ij|kl>
        elseif(INTEGd%iERItyp==2) then
         write(iunit) iorb,iorb2,iorb1,iorb3,ERIval_cmplx ! (ik|jl)
        else
         ! Nth
        endif
       endif
      else
       iorbm1=iorb-1
       iorb1m1=iorb1-1
       iorb2m1=iorb2-1
       iorb3m1=iorb3-1
       ERIval_cmplx=INTEGd%ERImolv_cmplx(iorbm1+iorb1m1*INTEGd%NBF2+iorb2m1*INTEGd%NBF3+iorb3m1*INTEGd%NBF4+1) ! <i+j*NBF2+k*NBF3+l*NBF4> 
       if(cdabs(ERIval_cmplx)>tol8) then
        write(iunit) iorb,iorb1,iorb2,iorb3,ERIval_cmplx
       endif
      endif
     enddo
    enddo
   enddo
  enddo 
  write(iunit) 0,0,0,0,complex_zero
  write(iunit) 0,0,0,0,complex_zero
 else
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
 endif
 close(iunit)

 ! Print hCORE
 open(unit=iunit,form='unformatted',file='hCORE')
 if(INTEGd%complex_ints) then
  do iorb=1,INTEGd%NBF_jkl
   do iorb1=1,INTEGd%NBF_jkl
    if(cdabs(INTEGd%hCORE_cmplx(iorb,iorb1))>tol8) then
     write(iunit) iorb,iorb1,INTEGd%hCORE_cmplx(iorb,iorb1)
    endif
   enddo
  enddo
  write(iunit) 0,0,complex_zero
 else
  do iorb=1,INTEGd%NBF_jkl
   do iorb1=1,INTEGd%NBF_jkl
    if(dabs(INTEGd%hCORE(iorb,iorb1))>tol8) then
     write(iunit) iorb,iorb1,INTEGd%hCORE(iorb,iorb1)
    endif
   enddo
  enddo
  write(iunit) 0,0,zero
 endif
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
 character(len=200)::msg
!************************************************************************
 
 ! Print FCIDUMP
 if(.not.INTEGd%complex_ints) then
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
 else ! TODO
  write(msg,'(a)') 'FCIDUMP file not programed/produced for complex orbitals'
  call write_output(msg)
 endif

end subroutine print_fcidump
!!***

end module m_integd
!!***
