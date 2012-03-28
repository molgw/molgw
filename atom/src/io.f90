module m_io

contains
!=========================================================================
subroutine dump_out_array(is_energy,title,n,nspin,array)
 use m_definitions
 implicit none
 logical,intent(in)            :: is_energy
 character(len=100),intent(in) :: title
 integer,intent(in)            :: n,nspin
 real(dp),intent(in)           :: array(n,nspin)
!=====
 integer,parameter :: MAXSIZE=200
!=====
 integer :: i,ispin
!=====

 write(*,'(/,x,a)') TRIM(title)

 if(is_energy) then
   if(nspin==1) then
     write(*,'(a)') '   #       [Ha]         [eV]      '
   else
     write(*,'(a)') '   #              [Ha]                      [eV]      '
     write(*,'(a)') '           spin 1       spin 2       spin 1       spin 2'
   endif
   do i=1,MIN(n,MAXSIZE)
     write(*,'(x,i3,2(2(x,f12.5)),2x)') i,array(i,:),array(i,:)*Ha_eV
   enddo
 else
   if(nspin==2) then
     write(*,'(a)') '           spin 1       spin 2 '
   endif
   do i=1,MIN(n,MAXSIZE)
     write(*,'(x,i3,2(2(x,f12.5)),2x)') i,array(i,:)
   enddo
 endif
 write(*,*)

end subroutine dump_out_array

!=========================================================================
subroutine dump_out_matrix(print_volume,title,n,nspin,matrix)
 use m_definitions
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

 if(print_volume<10) return

 write(*,'(/,x,a)') TRIM(title)

 do ispin=1,nspin
   if(nspin==2) write(*,'(a,i1)') ' spin polarization # ',ispin
   do i=1,MIN(n,MAXSIZE)
     write(*,'(x,i3,100(x,f12.5))') i,matrix(i,1:MIN(n,MAXSIZE),ispin)
   enddo
   write(*,*)
 enddo
 write(*,*)

end subroutine dump_out_matrix

!=========================================================================
subroutine output_homolumo(nbf,nspin,occupation,energy,homo,lumo)
 use m_definitions
 implicit none
 integer,intent(in)  :: nbf,nspin
 real(dp),intent(in) :: occupation(nbf,nspin),energy(nbf,nspin)
 real(dp),intent(out) :: homo,lumo
 real(dp) :: homo_tmp,lumo_tmp
 integer :: ispin,ibf

 homo_tmp=-1.d+5
 lumo_tmp= 1.d+5

 do ispin=1,nspin
   do ibf=1,nbf
     if(occupation(ibf,ispin)>completely_empty) then
       homo_tmp = MAX( homo_tmp , energy(ibf,ispin) )
     endif

     if(occupation(ibf,ispin)<1.0_dp - completely_empty ) then
       lumo_tmp = MIN( lumo_tmp , energy(ibf,ispin) )
     endif

   enddo
 enddo

 homo = homo_tmp
 lumo = lumo_tmp

 write(*,*)
 write(*,'(a,x,f12.6)') 'HOMO energy    [Ha]:',homo_tmp
 write(*,'(a,x,f12.6)') 'LUMO energy    [Ha]:',lumo_tmp
 write(*,'(a,x,f12.6)') 'HOMO-LUMO gap  [Ha]:',lumo_tmp-homo_tmp
 write(*,*)


end subroutine output_homolumo

!=========================================================================
subroutine read_inputparameter(calc_type,nspin,nscf,alpha_mixing,print_volume,basis_name,zatom,electrons,magnetization,basis_element)
 use m_definitions
 use m_calculation_type
 use m_warning
 use m_tools
 implicit none

 type(calculation_type),intent(out) :: calc_type
 integer,intent(out)                :: nspin,nscf
 integer,intent(out)                :: basis_element
 real(dp),intent(out)               :: alpha_mixing
 integer,intent(out)                :: print_volume  
 character(len=100),intent(out)     :: basis_name
 real(dp),intent(out)               :: zatom
 real(dp),intent(out)               :: electrons 
 real(dp),intent(out)               :: magnetization
!=====                              
 character(len=100)                 :: read_char
 character(len=100)                 :: read_line
 character(len=100)                 :: line_wocomment
 integer                            :: ipos
 integer                            :: istat
!=====


 read(*,*) read_char
 call init_calculation_type(calc_type,read_char)

 write(*,'(/,a,/)')    ' Summary of the input parameters '
 write(*,'(a20,2x,a)') ' calculation type: ',TRIM(read_char)

 read(*,*) nspin,zatom,electrons,magnetization
 if(nspin/=1 .AND. nspin/=2) stop'nspin in incorrect'
 if(zatom<1.d-5)             stop'zatom is negative'
 if(electrons<1.d-5)         stop'electrons is negative'
 if(magnetization<-1.d-5)    stop'magnetization is negative'
 if(magnetization>1.d-5 .AND. nspin==1) stop'magnetization is non-zero and nspin is 1'

 basis_element=NINT(zatom)
 read(*,fmt='(a100)') read_line !,basis_element
 ipos=index(read_line,'#',back=.false.)
 if(ipos==0) ipos=101
 line_wocomment(:ipos-1)=read_line(:ipos-1)
 line_wocomment=ADJUSTL(line_wocomment)
 ipos=index(line_wocomment,' ',back=.false.)
 basis_name = line_wocomment(:ipos-1)
 line_wocomment(1:)=line_wocomment(ipos+1:)
 read(line_wocomment,*,iostat=istat) basis_element
 if(istat==0) then
   msg='override the automatic basis set selection with the rescaled basis for element '//element_name(REAL(basis_element,dp))
   call issue_warning(msg)
 endif

 read(*,*) nscf,alpha_mixing
 if(nscf<1) stop'nscf too small'
 if(alpha_mixing<0.0 .OR. alpha_mixing > 1.0 ) stop'alpha_mixing should be inside [0,1]'

 read(*,*) print_volume

 !
 ! summarize input parameters
 write(*,'(a20,f8.4)') ' Atom Z: ',zatom
 write(*,'(a20,f8.4)') ' Electrons: ',electrons
 write(*,'(a20,f8.4)') ' Magnetization: ',magnetization
 write(*,'(a20,2x,a)') ' Basis set: ',basis_name
 write(*,'(a20,i3)')   ' Spin polarization: ',nspin
 write(*,'(a20,i3)')   ' SCF steps: ',nscf


 write(*,*)

end subroutine read_inputparameter

!=========================================================================
subroutine read_inputparameter_molecule(calc_type,nspin,nscf,alpha_mixing,print_volume,&
      basis_name,zatom,electrons,magnetization,basis_element,natom,x)
 use m_definitions
 use m_calculation_type
 use m_warning
 use m_tools
 implicit none

 type(calculation_type),intent(out) :: calc_type
 integer,intent(out)                :: nspin,nscf,natom
 integer,pointer,intent(inout)        :: basis_element(:)
 real(dp),intent(out)               :: alpha_mixing
 integer,intent(out)                :: print_volume  
 character(len=100),intent(out)     :: basis_name
 real(dp),pointer,intent(inout)   :: zatom(:)
 real(dp),pointer,intent(inout)   :: x(:,:)
 real(dp),intent(out)               :: electrons 
 real(dp),intent(out)               :: magnetization
!=====                              
 character(len=100)                 :: read_char
 character(len=100)                 :: read_line
 character(len=100)                 :: line_wocomment
 integer                            :: ipos
 integer                            :: istat,iatom
 real(dp)                           :: charge
!=====

 write(*,*) 'inside read_input'
 read(*,*) read_char
 call init_calculation_type(calc_type,read_char)

 write(*,'(/,a,/)')    ' Summary of the input parameters '
 write(*,'(a20,2x,a)') ' calculation type: ',TRIM(read_char)

 read(*,*) nspin,charge,magnetization
 if(nspin/=1 .AND. nspin/=2) stop'nspin in incorrect'
 if(magnetization<-1.d-5)    stop'magnetization is negative'
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

 read(*,*) nscf,alpha_mixing
 if(nscf<1) stop'nscf too small'
 if(alpha_mixing<0.0 .OR. alpha_mixing > 1.0 ) stop'alpha_mixing should be inside [0,1]'

 read(*,*) print_volume

 read(*,*) natom
 if(natom<1) stop'natom<1'

 allocate(zatom(natom),x(3,natom),basis_element(natom))
 do iatom=1,natom
   read(*,*) zatom(iatom),x(:,iatom)
 enddo
 basis_element(:)=NINT(zatom(:))

 electrons = SUM(zatom(:)) + charge



 !
 ! summarize input parameters
 write(*,'(a20,i3)')   ' Natom: ',natom
 write(*,'(a20,f8.4)') ' Electrons: ',electrons
 write(*,'(a20,f8.4)') ' Charge: ',charge
 write(*,'(a20,f8.4)') ' Magnetization: ',magnetization
 write(*,'(a20,2x,a)') ' Basis set: ',basis_name
 write(*,'(a20,i3)')   ' Spin polarization: ',nspin
 write(*,'(a20,i3)')   ' SCF steps: ',nscf
 write(*,'(a20,f8.4)')   ' Mixing: ',alpha_mixing


 write(*,*)

end subroutine read_inputparameter_molecule

!=========================================================================
subroutine plot_wfn(nspin,basis,c_matrix)
 use m_definitions
 use m_basis_set
 implicit none
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer,parameter          :: nr=4000
 integer                    :: ir,ibf
 real(dp)                   :: x(3),phi1,phi2,phase1,phase2
!=====

 phase1=1.0_dp
 phase2=1.0_dp
 x(:)=0.0_dp

 do ir=1,nr
   x(1)=DBLE(ir-1)*5.0d-3
   phi1=0.0_dp
   phi2=0.0_dp

   do ibf=1,basis%nbf
     phi1=phi1+eval_basis_function(basis%bf(ibf),x)*c_matrix(ibf,1,1)
     phi2=phi2+eval_basis_function(basis%bf(ibf),x)*c_matrix(ibf,1,nspin)
   enddo

   !
   ! turn the wfns so that they are all positive at r=0
   if(ir==1) then
     if(phi1<0.0_dp) phase1=-1.0_dp
     if(phi2<0.0_dp) phase2=-1.0_dp
   endif

   write(101,*) x(1),phi1*phase1,phi2*phase2
 enddo

end subroutine plot_wfn
!=========================================================================
end module m_io
