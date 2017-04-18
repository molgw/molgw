program spectrum
 use m_definitions
 use m_warning
 use ISO_C_BINDING
 implicit none
!=====
 integer                    :: ntau, itau,idir,iomega,nomega
 integer                    :: file_transforms, file_dipolar_spectra 
 integer                    :: file_dipole_time, file_excit_field, io
 integer                    :: num_rows
 real(dp)                   :: time_cur, time_min, time_sim
 real(dp)                   :: damp_factor
 real(dp)                   :: omega_factor
 complex(dp),allocatable    :: m_times(:)
 complex(dp),allocatable    :: dipole_time_ref(:,:)
 complex(dp),allocatable    :: dipole_time_damped(:,:)
 complex(dp),allocatable    :: trans_m_excit_field(:,:)
 complex(dp),allocatable    :: m_excit_field(:,:)
 complex(dp),allocatable    :: trans_dipole_time(:,:)
 logical                    :: file_exists
 character(len=100)         :: name_dipole_time, name_excitation
 character(len=100)         :: name_dipolar_spectra, name_transforms
 character(len=100)         :: test_string
 type(C_PTR)                :: plan

#ifdef HAVE_FFTW
#include "fftw3.f03"

! ntau=int((time_sim-time_min)/time_step)+1

 damp_factor=250.0_dp

 name_dipole_time='dipole_time.dat'
 name_excitation='excitation.dat'
 name_dipolar_spectra='dipolar_spectra.dat'
 name_transforms='transforms.dat'

 inquire(file=name_dipole_time,exist=file_exists)
 if(.NOT. file_exists) then
   call die('No dipole_time file found')
 endif

 inquire(file=name_excitation,exist=file_exists)
 if(.NOT. file_exists) then
   call die('No excitation file found')
 endif

 inquire(file=name_dipolar_spectra,exist=file_exists)
 if(.NOT. file_exists) then
   call issue_warning('Dipolar spectra file will be rewritten')
 endif

 ntau = 0 
 open (newunit=file_dipole_time, file = name_dipole_time)
 do
   if( ntau==0  ) then
     read(file_dipole_time,'(A)',iostat=io)test_string
     call get_number_of_elements(test_string,num_rows)
     write(stdout,*) test_string,num_rows
     if(num_rows/=4) call die('Dipole file must contain 4 columns: time dx dy dz')
   else 
     read(file_dipole_time,*,iostat=io)
   end if
   if (io/=0) exit
   ntau = ntau + 1
 end do

 write(stdout,*) ntau

 do itau=1,ntau
   read(file_dipole_time,*,iostat=io)time_min
   write(stdout,*) time_min
   if (io/=0) exit
   ntau = ntau + 1
 end do
 close (file_dipole_time)


 stop

 allocate(m_times(ntau))

 time_cur=time_min
 do itau=1,ntau
   dipole_time_damped(itau,:)=dipole_time_ref(itau,:)*exp(-m_times(itau) / damp_factor )
 end do

 do idir=1,3
   plan = fftw_plan_dft_1d(ntau,dipole_time_damped(1:ntau,idir),trans_dipole_time(1:ntau,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,dipole_time_damped(1:ntau,idir),trans_dipole_time(1:ntau,idir))
   call fftw_destroy_plan(plan)

 end do

 trans_dipole_time =  trans_dipole_time / ntau

 do idir=1,3
   plan = fftw_plan_dft_1d(ntau,m_excit_field(1:ntau,idir),trans_m_excit_field(1:ntau,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,m_excit_field(1:ntau,idir),trans_m_excit_field(1:ntau,idir))
   call fftw_destroy_plan(plan)
 end do

 trans_m_excit_field(:,:)=trans_m_excit_field(:,:) / ntau

 open(newunit=file_transforms,file=name_transforms)
 open(newunit=file_dipolar_spectra,file=name_dipolar_spectra)
 write(file_dipolar_spectra,*) "# omega(eV), Average, xx, xy, xz, yx, yy, yz, zx, zy, zz"
 write(file_transforms,*) "# omega(eV), |E_x(omega)|, real(E_x(omega)), aimag(E_x(omega)), |d_x(omega)|, real(d_x(omega)), aimag(d_x(omega))"

 nomega=ntau
 do iomega=1,nomega/2 ! here iomega correspons to frequency
   omega_factor = 4.0_dp * pi * 2 * pi * (iomega-1) / time_sim / c_speedlight
   if(2*pi * iomega / time_sim * Ha_eV < 500.0_dp) then
     write(file_dipolar_spectra,"(x,es16.8E3)",advance='no')  2*pi * iomega / time_sim * Ha_eV

     write(file_dipolar_spectra,"(3x,es16.8E3)",advance='no') &
             SUM( aimag((trans_dipole_time(iomega,:)) / (trans_m_excit_field(iomega,:))) ) / 3.0_dp * omega_factor

     do idir=1,3
!       if( excit_dir(idir) > 1.0E-10_dp ) then
         write(file_dipolar_spectra,"(3(3x,es16.8E3))",advance='no') aimag((trans_dipole_time(iomega,:))/(trans_m_excit_field(iomega,idir))) * omega_factor
!       else
!         write(file_dipolar_spectra,"(3(es16.8E3))",advance='no') -1.0_dp, -1.0_dp, -1.0_dp
!       end if
       if(idir==3) write(file_dipolar_spectra,*)
     end do
     write(file_transforms,*) pi * iomega / time_sim * Ha_eV, abs(trans_m_excit_field(iomega,1)), real(trans_m_excit_field(iomega,1)), aimag(trans_m_excit_field(iomega,1)), abs(trans_dipole_time(iomega,1)), real(trans_dipole_time(iomega,1)), aimag(trans_dipole_time(iomega,1))
   end if
 end do

 close(file_transforms)
 close(file_dipolar_spectra)
#else
 call issue_warning("tddft: calculate_propagation; fftw is not present")
#endif

end program spectrum

