!=========================================================================
#include "macros.h"
!=========================================================================
module m_inputparam
 use m_definitions
 use m_mpi
 use m_atoms
 use m_calculation_type
 use m_basis_set
 use m_scf

 integer,protected                   :: nspin
 integer,protected                   :: nscf
 real(dp),protected                  :: alpha_mixing
 character(len=100),protected        :: basis_name
 integer,protected                   :: gaussian_type
 real(dp),protected                  :: electrons
 real(dp),protected                  :: magnetization
 type(calculation_type),protected    :: calc_type
 character(len=100),protected        :: quadrature_name

 logical,protected            :: print_matrix
 logical,protected            :: print_basis
 logical,protected            :: print_eri
 logical,protected            :: print_densitymatrix
 logical,protected            :: print_wfn
 logical,protected            :: print_specfunc


contains


!=========================================================================
subroutine read_inputparameter_molecule()
! use m_tools
 implicit none

 character(len=100)                 :: read_char
 character(len=100)                 :: read_line
 character(len=100)                 :: line_wocomment
 character(len=100)                 :: mixing_name
 character(len=100)                 :: dft_accuracy
 character(len=100)                 :: gaussian_name
 integer                            :: print_volume
 integer                            :: ipos,jpos
 integer                            :: istat,iatom,jatom
 real(dp)                           :: charge,length_factor
!=====

 read(*,*) read_char
 call init_calculation_type(calc_type,read_char)

 WRITE_MASTER(*,'(/,a,/)')    ' Summary of the input parameters '
 WRITE_MASTER(*,'(a20,2x,a)') ' calculation type: ',TRIM(read_char)

 read(*,*) nspin,charge,magnetization
 if(nspin/=1 .AND. nspin/=2) stop'nspin in incorrect'
 if(magnetization<-1.d-5)    stop'magnetization is negative'
 if(magnetization>1.d-5 .AND. nspin==1) stop'magnetization is non-zero and nspin is 1'

 read(*,fmt='(a100)') read_line 
 ipos=index(read_line,'#',back=.false.)
 if(ipos==0) ipos=101
 line_wocomment(:ipos-1)=read_line(:ipos-1)
 line_wocomment=ADJUSTL(line_wocomment)
 jpos=index(line_wocomment,' ',back=.false.)
 basis_name = line_wocomment(:jpos-1)
 line_wocomment(1:ipos-jpos-1)=line_wocomment(jpos+1:ipos-1)
 select case(TRIM(ADJUSTL(line_wocomment(1:ipos-jpos-1))))
 case('PURE','pure')
   gaussian_type=PURE
   gaussian_name='PURE'
 case('CART','cart')
   gaussian_type=CARTESIAN
   gaussian_name='CARTESIAN'
 case default
   stop'Error in input line 3: second keyword should either PURE or CART'
 end select

 read(*,*) nscf,alpha_mixing,mixing_name,dft_accuracy
 if(nscf<1) stop'nscf too small'
 if(alpha_mixing<0.0 .OR. alpha_mixing > 1.0 ) stop'alpha_mixing should be inside [0,1]'
 select case(TRIM(mixing_name))
 case('SIMPLE')
   mixing_scheme = simple_mixing
 case('PULAY')
   mixing_scheme = pulay_mixing
 case default
   stop'mixing scheme not recognized'
 end select

 select case(TRIM(dft_accuracy))
 case('VERYLOW','verylow','VL','vl')
   quadrature_name = 'very low'
 case('LOW','low','L','l')
   quadrature_name = 'low'
 case('MEDIUM','medium','M','m')
   quadrature_name = 'medium'
 case('HIGH','high','H','h')
   quadrature_name = 'high'
 case('VERYHIGH','veryhigh','VH','vh')
   quadrature_name = 'very high'
 case default
   stop'integration quality not recognized'
 end select

 read(*,*) print_volume

 print_matrix        = MODULO(print_volume       ,2)==1
 print_basis         = MODULO(print_volume/10    ,2)==1
 print_eri           = MODULO(print_volume/100   ,2)==1
 print_densitymatrix = MODULO(print_volume/1000  ,2)==1
 print_wfn           = MODULO(print_volume/10000 ,2)==1
 print_specfunc      = MODULO(print_volume/100000,2)==1

 
 read(*,*) natom,read_char
 if(natom<1) stop'natom<1'

 !
 ! lengths are stored internally in bohr
 !
 select case(TRIM(read_char))
 case('A','a','Angstrom','ANGSTROM','angstrom')
   length_factor=1.0_dp/bohr_A
 case('bohr','BOHR','Bohr','AU','au','a.u.','a.u','A.U','A.U.')
   length_factor=1.0_dp
 case default
   stop'units for lengths in input file not understood'
 end select

 allocate(zatom(natom),x(3,natom),basis_element(natom))
 do iatom=1,natom
   read(*,*) zatom(iatom),x(:,iatom)
 enddo
 x(:,:) = x(:,:) * length_factor

 !
 ! check for atoms too close
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( SQRT( SUM( (x(:,iatom)-x(:,jatom))**2 ) ) < 0.2 ) then
       WRITE_MASTER(*,*) 'atoms',iatom,jatom
       WRITE_MASTER(*,*) 'are closer than 0.2 bohr'
       stop'stop here'
     endif
   enddo
 enddo


 basis_element(:)=NINT(zatom(:))

 electrons = SUM(zatom(:)) - charge


 !
 ! summarize input parameters
 WRITE_MASTER(*,'(a25,i3)')   ' Natom: ',natom
 WRITE_MASTER(*,'(a25,f8.4)') ' Electrons: ',electrons
 WRITE_MASTER(*,'(a25,f8.4)') ' Charge: ',charge
 WRITE_MASTER(*,'(a25,f8.4)') ' Magnetization: ',magnetization
 WRITE_MASTER(*,'(a25,2x,a)') ' Basis set: ',basis_name
 WRITE_MASTER(*,'(a25,2x,a)') ' Gaussian type: ',gaussian_name
 WRITE_MASTER(*,'(a25,i3)')   ' Spin polarization: ',nspin
 WRITE_MASTER(*,'(a25,i3)')   ' SCF steps: ',nscf
 WRITE_MASTER(*,'(a25,f8.4)') ' Mixing: ',alpha_mixing
 WRITE_MASTER(*,'(a25,2x,a)') ' Quadrature accuracy: ',quadrature_name
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a19)')      ' Print volume:'
 WRITE_MASTER(*,'(a30,l3)')   ' - matrices details:   ',print_matrix        
 WRITE_MASTER(*,'(a30,l3)')   ' - basis set details:  ',print_basis
 WRITE_MASTER(*,'(a30,l3)')   ' - ERI file:           ',print_eri           
 WRITE_MASTER(*,'(a30,l3)')   ' - density matrix file:',print_densitymatrix 
 WRITE_MASTER(*,'(a30,l3)')   ' - plot some wfns:     ',print_wfn           
 WRITE_MASTER(*,'(a30,l3)')   ' - dump spectral functs',print_specfunc      


 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) '================================'
 WRITE_MASTER(*,*) '      atom list'
 WRITE_MASTER(*,*) '                       bohr                                        angstrom'
 do iatom=1,natom
   WRITE_MASTER(*,'(2x,a2,3(x,f12.6),6x,3(x,f12.6))') element_name(zatom(iatom)),x(:,iatom),x(:,iatom)*bohr_A
 enddo


 WRITE_MASTER(*,*) '================================'
 WRITE_MASTER(*,*)

end subroutine read_inputparameter_molecule

end module m_inputparam

