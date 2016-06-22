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
subroutine find_qp_energy_linearization(nomegai,omegai,nsemin,nsemax,selfenergy_omega,nstate,exchange_m_vxc,energy0,energy_qp_z,zz)
 implicit none

 integer,intent(in)            :: nomegai
 integer,intent(in)            :: nsemin,nsemax
 integer,intent(in)            :: nstate
 real(dp),intent(in)           :: omegai(-nomegai:nomegai)
 real(dp),intent(in)           :: selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nspin)
 real(dp),intent(in)           :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out)          :: energy_qp_z(nstate,nspin)
 real(dp),intent(out),optional :: zz(nsemin:nsemax,nspin)
!=====
 integer  :: astate,aspin
 real(dp) :: zz_a(nspin)
!=====

 ! First, a dummy initialization
 energy_qp_z(:,:) = energy0(:,:)

 ! Then overwrite the interesting energy with the calculated GW one
 do astate=nsemin,nsemax

   if( nomegai > 0 .AND. PRESENT(zz) ) then
     zz_a(:) = ( selfenergy_omega(1,astate,:) - selfenergy_omega(-1,astate,:) ) / ( omegai(1) - omegai(-1) )
     zz_a(:) = 1.0_dp / ( 1.0_dp - zz_a(:) )
     ! Contrain Z to be in [0:1] to avoid crazy values
     do aspin=1,nspin
       zz_a(aspin) = MIN( MAX(zz_a(aspin),0.0_dp) , 1.0_dp )
     enddo

     zz(astate,:)          = zz_a(:)
     energy_qp_z(astate,:) = energy0(astate,:) + zz_a(:) * ( selfenergy_omega(0,astate,:) + exchange_m_vxc(astate,:) )

   else

     energy_qp_z(astate,:) = energy0(astate,:) + selfenergy_omega(0,astate,:) + exchange_m_vxc(astate,:)

   endif

 enddo


end subroutine find_qp_energy_linearization


!=========================================================================
subroutine find_qp_energy_graphical(nomegai,omegai,nsemin,nsemax,selfenergy_omega,nstate,exchange_m_vxc,energy0,energy_qp_g)
 implicit none

 integer,intent(in)   :: nomegai
 integer,intent(in)   :: nsemin,nsemax
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: omegai(-nomegai:nomegai)
 real(dp),intent(in)  :: selfenergy_omega(-nomegai:nomegai,nsemin:nsemax,nspin)
 real(dp),intent(in)  :: exchange_m_vxc(nstate,nspin),energy0(nstate,nspin)
 real(dp),intent(out) :: energy_qp_g(nstate,nspin)
!=====
 integer  :: astate,aspin
 real(dp) :: sigma_xc_m_vxc(-nomegai:nomegai)
!=====

 ! First, a dummy initialization
 energy_qp_g(:,:) = 0.0_dp

 ! Then overwrite the interesting energy with the calculated GW one
 do astate=nsemin,nsemax

   if( MODULO(astate-nsemin,nproc) /= rank ) cycle

   do aspin=1,nspin
     sigma_xc_m_vxc(:) = selfenergy_omega(:,astate,aspin) + exchange_m_vxc(astate,aspin)
     energy_qp_g(astate,aspin) = find_fixed_point(nomegai,omegai,sigma_xc_m_vxc) + energy0(astate,aspin)
   enddo

 enddo

 call xsum(energy_qp_g)

 energy_qp_g(:nsemin-1,:) = energy0(:nsemin-1,:)
 energy_qp_g(nsemax+1:,:) = energy0(nsemax+1:,:)



end subroutine find_qp_energy_graphical


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


!=========================================================================
subroutine output_qp_energy(calcname,nstate,nsemin,nsemax,energy0,exchange_m_vxc,selfenergy,energy1,energy2,zz)
 implicit none
 
 character(len=*)             :: calcname
 integer                      :: nstate,nsemin,nsemax
 real(dp),intent(in)          :: energy0(nstate,nspin),exchange_m_vxc(nstate,nspin)
 real(dp),intent(in)          :: selfenergy(nsemin:nsemax,nspin),energy1(nstate,nspin)
 real(dp),intent(in),optional :: energy2(nstate,nspin),zz(nsemin:nsemax,nspin)
!=====
 integer  :: astate
!=====


 write(stdout,'(/,x,a,x,a)') TRIM(calcname),'eigenvalues (eV)'

 if( PRESENT(zz) .AND. PRESENT(energy2) ) then

   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,10x,a,11x,a,10x,a)') '#','E0','SigX-Vxc','SigC','Z','E_Z','E_qp'
   else
     write(stdout,'(3x,a,15x,a,22x,a,19x,a,24x,a,23x,a,23x,a)') '#','E0','SigX-Vxc','SigC','Z','E_Z','E_qp'
     write(stdout,'(12x,12(a4,9x))') ' up ','down',' up ','down',' up ','down',' up ','down',' up ','down',' up ','down'
   endif

   do astate=nsemin,nsemax

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy0(astate,:)*Ha_eV,   &
                                        exchange_m_vxc(astate,:)*Ha_eV,   &
                                        selfenergy(astate,:)*Ha_eV,       &
                                        zz(astate,:),                     &
                                        energy1(astate,:)*Ha_eV,          &
                                        energy2(astate,:)*Ha_eV
   enddo

 else

   if(nspin==1) then
     write(stdout,'(3x,a,8x,a,9x,a,7x,a,9x,a)') '#','E0','SigX-Vxc','SigC','E_qp'
   else
     write(stdout,'(3x,a,15x,a,22x,a,20x,a,22x,a)') '#','E0','SigX-Vxc','SigC','E_qp'
     write(stdout,'(12x,8(a4,9x))') ' up ','down',' up ','down',' up ','down',' up ','down'
   endif

   do astate=nsemin,nsemax

     write(stdout,'(i4,x,20(x,f12.6))') astate,energy0(astate,:)*Ha_eV,       &
                                        exchange_m_vxc(astate,:)*Ha_eV,       &
                                        selfenergy(astate,:)*Ha_eV,           &
                                        energy1(astate,:)*Ha_eV
   enddo



 endif

end subroutine output_qp_energy


end module m_selfenergy_tools


!=========================================================================
