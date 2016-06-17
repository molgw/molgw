!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! basic tools to play with the calculated self-energies
!
!=========================================================================
module m_selfenergy_tools
 use m_definitions
 use m_warning
 use m_mpi
 use m_inputparam,only: nspin
! use m_memory
! use m_timing
! use m_scalapack
! use m_tools,only: invert


contains


!=========================================================================
subroutine write_selfenergy_omega(filename_root,nstate,energy0,exchange_m_vxc,nomegai,omegai,nsemin,nsemax,selfenergy_omega)
 implicit none

 character(len=*)    :: filename_root
 integer,intent(in)  :: nstate,nsemin,nsemax,nomegai
 real(dp),intent(in) :: energy0(nstate,nspin),exchange_m_vxc(nstate,nspin)
 real(dp),intent(in) :: omegai(nomegai)
 real(dp),intent(in) :: selfenergy_omega(nomegai,nsemin:nsemax,nspin)
!=====
 character(len=3)   :: ctmp
 character(len=256) :: filename
 integer :: selfenergyfile
 integer :: astate
 integer :: iomegai
!=====

 ! Just the master writes
 if( .NOT. is_iomaster ) return

 write(stdout,'(/,x,a)') 'Write Sigma(omega) on file'

 !
 ! omega is defined with respect to energy0_a
 ! Absolute omega is omega + energy0_a
 !
 do astate=nsemin,nsemax
   write(ctmp,'(i3.3)') astate
   filename = TRIM(filename_root) // '_state' // TRIM(ctmp) // '.dat'
   write(stdout,'(x,a,a)') 'Writing selfenergy in file: ', TRIM(filename)
   open(newunit=selfenergyfile,file=filename)

   write(selfenergyfile,'(a)') '# omega (eV)             SigmaC (eV)    omega - e_gKS - Vxc + SigmaX (eV)     A (eV^-1) '

   do iomegai=1,nomegai
     write(selfenergyfile,'(20(f12.6,2x))') ( omegai(iomegai) + energy0(astate,:) )*Ha_eV,               &
                                        selfenergy_omega(iomegai,astate,:)*Ha_eV,                    &
                                        ( omegai(iomegai) - exchange_m_vxc(astate,:) )*Ha_eV,     &
                                        1.0_dp/pi/ABS( omegai(iomegai) - exchange_m_vxc(astate,:) &
                                                - selfenergy_omega(iomegai,astate,:) ) / Ha_eV
   enddo
   write(selfenergyfile,*)
   close(selfenergyfile)

 enddo


end subroutine write_selfenergy_omega


!=========================================================================
function find_fixed_point(nx,xx,fx) result(fixed_point)
 implicit none
 integer,intent(in)  :: nx
 real(dp),intent(in) :: xx(-nx:nx)
 real(dp),intent(in) :: fx(-nx:nx)
 real(dp)            :: fixed_point
!=====
 integer             :: ix,imin
 real(dp)            :: rmin
 real(dp)            :: gx(-nx:nx)
 real(dp)            :: gpx
!=====


 gx(:) = fx(:) - xx(:)

 rmin = HUGE(1.0_dp)
 do ix=-nx,nx
   if( ABS(gx(ix)) < rmin ) then
     rmin = ABS(gx(ix))
     imin = ix
   endif
 enddo 


 if( imin == -nx .OR. imin == nx) then
   fixed_point = xx(imin)
 else 
   if( gx(imin)*gx(imin+1) < 0.0_dp )  then
     gpx = ( gx(imin+1) - gx(imin) ) / ( xx(imin+1) - xx(imin) )
     fixed_point = xx(imin) - gx(imin) / gpx 
   else if( gx(imin)*gx(imin-1) < 0.0_dp )  then
     gpx = ( gx(imin) - gx(imin-1) ) / ( xx(imin) - xx(imin-1) )
     fixed_point = xx(imin-1) - gx(imin-1) / gpx 
   else
     fixed_point = xx(imin)
   endif
 endif


end function find_fixed_point




end module m_selfenergy_tools


!=========================================================================
