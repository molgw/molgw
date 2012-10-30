!=========================================================================
#include "macros.h"
!=========================================================================
subroutine header()
 use m_definitions
 use m_mpi
 implicit none
 integer :: values(8) 
!=====

 WRITE_MASTER(*,'(x,70("="))') 
 WRITE_MASTER(*,'(/,/,12x,a,/,/)') ' Welcome to the fascinating world of MOLGW'
 WRITE_MASTER(*,'(x,70("="))') 
 WRITE_MASTER(*,*)

 call date_and_time(VALUES=values)

 WRITE_MASTER(*,'(a,i2.2,a,i2.2,a,i4.4)') ' Today is ',values(3),'/',values(2),'/',values(1)
 WRITE_MASTER(*,'(a,i2.2,a,i2.2)')        ' It is now ',values(5),':',values(6)
 select case(values(5))
 case(03,04,05,06,07)
   WRITE_MASTER(*,*) 'And it is too early to work. Go back to sleep'
 case(22,23,00,01,02)
   WRITE_MASTER(*,*) 'And it is too late to work. Go to bed and have a sleep'
 case(12,13)
   WRITE_MASTER(*,*) 'Go and get some good food'
 case(17)
   WRITE_MASTER(*,*) 'Dont forget to go and get the kids'
 case default
   WRITE_MASTER(*,*) 'And it is perfect time to work'
 end select

end subroutine header

!=========================================================================
subroutine dump_out_array(is_energy,title,n,nspin,array)
 use m_definitions
 use m_mpi
 implicit none
 logical,intent(in)            :: is_energy
 character(len=100),intent(in) :: title
 integer,intent(in)            :: n,nspin
 real(dp),intent(in)           :: array(n,nspin)
!=====
 integer,parameter :: MAXSIZE=300
!=====
 integer :: i,ispin
!=====

 WRITE_MASTER(*,'(/,x,a)') TRIM(title)

 if(is_energy) then
   if(nspin==1) then
     WRITE_MASTER(*,'(a)') '   #       [Ha]         [eV]      '
   else
     WRITE_MASTER(*,'(a)') '   #              [Ha]                      [eV]      '
     WRITE_MASTER(*,'(a)') '           spin 1       spin 2       spin 1       spin 2'
   endif
   do i=1,MIN(n,MAXSIZE)
     WRITE_MASTER(*,'(x,i3,2(2(x,f12.5)),2x)') i,array(i,:),array(i,:)*Ha_eV
   enddo
 else
   if(nspin==2) then
     WRITE_MASTER(*,'(a)') '           spin 1       spin 2 '
   endif
   do i=1,MIN(n,MAXSIZE)
     WRITE_MASTER(*,'(x,i3,2(2(x,f12.5)),2x)') i,array(i,:)
   enddo
 endif
 WRITE_MASTER(*,*)

end subroutine dump_out_array

!=========================================================================
subroutine dump_out_matrix(print_volume,title,n,nspin,matrix)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)            :: print_volume       
 character(len=100),intent(in) :: title
 integer,intent(in)            :: n,nspin
 real(dp),intent(in)           :: matrix(n,n,nspin)
!=====
 integer,parameter :: MAXSIZE=25
!=====
 integer :: i,ispin
!=====

 if(MODULO(print_volume,2)<1) return

 WRITE_MASTER(*,'(/,x,a)') TRIM(title)

 do ispin=1,nspin
   if(nspin==2) then
     WRITE_MASTER(*,'(a,i1)') ' spin polarization # ',ispin
   endif
   do i=1,MIN(n,MAXSIZE)
     WRITE_MASTER(*,'(x,i3,100(x,f12.5))') i,matrix(i,1:MIN(n,MAXSIZE),ispin)
   enddo
   WRITE_MASTER(*,*)
 enddo
 WRITE_MASTER(*,*)

end subroutine dump_out_matrix

!=========================================================================
subroutine output_homolumo(nbf,nspin,occupation,energy,homo,lumo)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)  :: nbf,nspin
 real(dp),intent(in) :: occupation(nbf,nspin),energy(nbf,nspin)
 real(dp),intent(out) :: homo(nspin),lumo(nspin)
 real(dp) :: homo_tmp,lumo_tmp
 integer :: ispin,ibf


 do ispin=1,nspin
   homo_tmp=-1.d+5
   lumo_tmp= 1.d+5
   do ibf=1,nbf
     if(occupation(ibf,ispin)>completely_empty) then
       homo_tmp = MAX( homo_tmp , energy(ibf,ispin) )
     endif

     if(occupation(ibf,ispin)<1.0_dp - completely_empty ) then
       lumo_tmp = MIN( lumo_tmp , energy(ibf,ispin) )
     endif

   enddo
   homo(ispin) = homo_tmp
   lumo(ispin) = lumo_tmp
 enddo


 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a,2(3x,f12.6))') ' HOMO energy    [Ha]:',homo(:)
 WRITE_MASTER(*,'(a,2(3x,f12.6))') ' LUMO energy    [Ha]:',lumo(:)
 WRITE_MASTER(*,'(a,2(3x,f12.6))') ' HOMO-LUMO gap  [Ha]:',lumo(:)-homo(:)
 WRITE_MASTER(*,*)


end subroutine output_homolumo

!=========================================================================
subroutine read_inputparameter_molecule(calc_type,nspin,nscf,alpha_mixing,print_volume,&
      basis_name,electrons,magnetization)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_warning
 use m_tools
 use m_atoms
 use m_scf
 implicit none

 type(calculation_type),intent(out) :: calc_type
 integer,intent(out)                :: nspin,nscf
 real(dp),intent(out)               :: alpha_mixing
 integer,intent(out)                :: print_volume  
 character(len=100),intent(out)     :: basis_name
 real(dp),intent(out)               :: electrons 
 real(dp),intent(out)               :: magnetization
!=====                              
 character(len=100)                 :: read_char
 character(len=100)                 :: read_line
 character(len=100)                 :: line_wocomment
 character(len=100)                 :: mixing_name
 integer                            :: ipos
 integer                            :: istat,iatom,jatom
 real(dp)                           :: charge,length_factor
!=====

 read(*,*) read_char
 call init_calculation_type(calc_type,read_char)

 WRITE_MASTER(*,'(/,a,/)')    ' Summary of the input parameters '
 WRITE_MASTER(*,'(a20,2x,a)') ' calculation type: ',TRIM(read_char)

 read(*,*) nspin,charge,magnetization
 if(nspin/=1 .AND. nspin/=2) stop'nspin in incorrect'
! if(magnetization<-1.d-5)    stop'magnetization is negative'
 if(magnetization>1.d-5 .AND. nspin==1) stop'magnetization is non-zero and nspin is 1'

! basis_element=NINT(zatom)
 read(*,fmt='(a100)') read_line !,basis_element
 ipos=index(read_line,'#',back=.false.)
 if(ipos==0) ipos=101
 line_wocomment(:ipos-1)=read_line(:ipos-1)
 line_wocomment=ADJUSTL(line_wocomment)
 ipos=index(line_wocomment,' ',back=.false.)
 basis_name = line_wocomment(:ipos-1)
! line_wocomment(1:)=line_wocomment(ipos+1:)
! read(line_wocomment,*,iostat=istat) basis_element
! if(istat==0) then
!   msg='override the automatic basis set selection with the rescaled basis for element '//element_name(REAL(basis_element,dp))
!   call issue_warning(msg)
! endif

 read(*,*) nscf,alpha_mixing,mixing_name
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

 read(*,*) print_volume

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
 WRITE_MASTER(*,'(a20,i3)')   ' Natom: ',natom
 WRITE_MASTER(*,'(a20,f8.4)') ' Electrons: ',electrons
 WRITE_MASTER(*,'(a20,f8.4)') ' Charge: ',charge
 WRITE_MASTER(*,'(a20,f8.4)') ' Magnetization: ',magnetization
 WRITE_MASTER(*,'(a20,2x,a)') ' Basis set: ',basis_name
 WRITE_MASTER(*,'(a20,i3)')   ' Spin polarization: ',nspin
 WRITE_MASTER(*,'(a20,i3)')   ' SCF steps: ',nscf
 WRITE_MASTER(*,'(a20,f8.4)') ' Mixing: ',alpha_mixing
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a19)')      ' Print volume:'
 WRITE_MASTER(*,'(a30,l3)')   ' - matrices details:   ',MODULO(print_volume      ,2)==1
 WRITE_MASTER(*,'(a30,l3)')   ' - basis set details:  ',MODULO(print_volume/10   ,2)==1
 WRITE_MASTER(*,'(a30,l3)')   ' - ERI file:           ',MODULO(print_volume/100  ,2)==1
 WRITE_MASTER(*,'(a30,l3)')   ' - density matrix file:',MODULO(print_volume/1000 ,2)==1
 WRITE_MASTER(*,'(a30,l3)')   ' - plot some wfns:     ',MODULO(print_volume/10000,2)==1


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

!=========================================================================
subroutine plot_wfn(nspin,basis,c_matrix)
 use m_definitions
 use m_mpi
 use m_basis_set,only: basis_set,eval_basis_function
 implicit none
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer,parameter          :: nr=10000
 real(dp),parameter         :: length=15.0_dp
 integer                    :: ir,ibf
 real(dp)                   :: x(3),phi1,phi2,phase1,phase2
!=====

 phase1=1.0_dp
 phase2=1.0_dp
 x(:)=0.0_dp

 do ir=1,nr
   x(1)= length * 2.0_dp * ( (ir-1.0_dp) / REAL(nr-1,dp) - 0.5_dp )
   phi1=0.0_dp
   phi2=0.0_dp

   do ibf=1,basis%nbf
     phi1=phi1+eval_basis_function(basis%bf(ibf),x)*c_matrix(ibf,1,2)
     phi2=phi2+eval_basis_function(basis%bf(ibf),x)*c_matrix(ibf,2,2)
   enddo

   !
   ! turn the wfns so that they are all positive at a given point
   if(ir==1) then
     if(phi1<0.0_dp) phase1=-1.0_dp
     if(phi2<0.0_dp) phase2=-1.0_dp
   endif

   WRITE_MASTER(101,*) x(1),phi1*phase1,phi2*phase2
   WRITE_MASTER(102,*) x(1),phi1**2,phi2**2
 enddo

end subroutine plot_wfn

!=========================================================================
