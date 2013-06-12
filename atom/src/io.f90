!=========================================================================
#include "macros.h"
!=========================================================================
subroutine header()
 use m_definitions
 use m_mpi
 use m_warning
 implicit none
 integer           :: values(8) 
 character(len=12) :: chartmp
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


 WRITE_MASTER(*,*) 'Compilation options'
#ifdef HAVE_LIBXC
 call xc_f90_version(values(1),values(2))
 WRITE_ME(chartmp,'(i2,a,i2)') values(1),'.',values(2)
 msg='LIBXC version '//TRIM(chartmp)
 call issue_warning(msg)
#endif
#ifdef CHI0
 msg='CHI0 option has been swichted on at compilation time'
 call issue_warning(msg)
#endif
#ifdef OPENMP
 write(msg,'(i6)') OMP_get_max_threads()
 msg='OPENMP option is activated with threads number'//msg
 call issue_warning(msg)
#endif
#ifdef LOW_MEMORY2
 msg='LOW_MEMORY version 2 option has been swichted on at compilation time'
 call issue_warning(msg)
#endif
#ifdef LOW_MEMORY3
 msg='LOW_MEMORY version 3 option has been swichted on at compilation time'
 call issue_warning(msg)
#endif
#ifdef CASIDA
 msg='CASIDA option has been swichted on at compilation time'
 call issue_warning(msg)
#endif
#ifdef HAVE_MPI
 msg='Running with MPI'
 call issue_warning(msg)
#endif
#ifdef SCALAPACK
 msg='Running with SCALAPACK'
 call issue_warning(msg)
#endif


end subroutine header

!=========================================================================
subroutine dump_out_occupation(title,n,nspin,occupation)
 use m_definitions
 use m_mpi
 implicit none
 character(len=100),intent(in) :: title
 integer,intent(in)            :: n,nspin
 real(dp),intent(in)           :: occupation(n,nspin)
!=====
 integer,parameter :: MAXSIZE=1000
!=====
 integer :: i,ispin
!=====

 WRITE_MASTER(*,'(/,x,a)') TRIM(title)

 if(nspin==2) then
   WRITE_MASTER(*,'(a)') '           spin 1       spin 2 '
 endif
 do i=1,MIN(n,MAXSIZE)
   WRITE_MASTER(*,'(x,i3,2(2(x,f12.5)),2x)') i,occupation(i,:)
 enddo

 WRITE_MASTER(*,*)

end subroutine dump_out_occupation

!=========================================================================
subroutine dump_out_eigenenergy(title,n,nspin,occupation,energy)
 use m_definitions
 use m_mpi
 implicit none
 character(len=100),intent(in) :: title
 integer,intent(in)            :: n,nspin
 real(dp),intent(in)           :: occupation(n,nspin),energy(n,nspin)
!=====
 integer,parameter :: MAXSIZE=300
!=====
 integer  :: i,ispin
 real(dp) :: spin_fact
!=====

 spin_fact = REAL(-nspin+3,dp)

 WRITE_MASTER(*,'(/,x,a)') TRIM(title)

 if(nspin==1) then
   WRITE_MASTER(*,'(a)') '   #       [Ha]         [eV]      '
 else
   WRITE_MASTER(*,'(a)') '   #              [Ha]                      [eV]      '
   WRITE_MASTER(*,'(a)') '           spin 1       spin 2       spin 1       spin 2'
 endif
 do i=1,MIN(n,MAXSIZE)
   WRITE_MASTER(*,'(x,i3,2(2(x,f12.5)),2x)') i,energy(i,:),energy(i,:)*Ha_eV
   if(i<n) then
     if( ANY( occupation(i+1,:) < spin_fact/2.0_dp .AND. occupation(i,:) > spin_fact/2.0 ) ) then 
        if(nspin==1) then
          WRITE_MASTER(*,'(a)') '  -----------------------------'
        else
          WRITE_MASTER(*,'(a)') '  -------------------------------------------------------'
        endif
     endif
   endif
 enddo

 WRITE_MASTER(*,*)

end subroutine dump_out_eigenenergy

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
      basis_name,gaussian_type,electrons,magnetization,nradial_grid,nangular_grid)
 use m_definitions
 use m_mpi
 use m_calculation_type
 use m_warning
 use m_tools
 use m_atoms
 use m_scf
 use m_basis_set
 implicit none

 type(calculation_type),intent(out) :: calc_type
 integer,intent(out)                :: nspin,nscf
 real(dp),intent(out)               :: alpha_mixing
 integer,intent(out)                :: print_volume  
 character(len=100),intent(out)     :: basis_name
 integer,intent(out)                :: gaussian_type
 integer,intent(out)                :: nradial_grid,nangular_grid
 real(dp),intent(out)               :: electrons 
 real(dp),intent(out)               :: magnetization
!=====                              
 character(len=100)                 :: read_char
 character(len=100)                 :: read_line
 character(len=100)                 :: line_wocomment
 character(len=100)                 :: mixing_name
 character(len=100)                 :: dft_accuracy
 character(len=100)                 :: gaussian_name
 character(len=100)                 :: quadrature_name
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
   nradial_grid  =  6
   nangular_grid =  6 
   quadrature_name = 'very low'
 case('LOW','low','L','l')
   nradial_grid  = 10
   nangular_grid = 14 
   quadrature_name = 'low'
 case('MEDIUM','medium','M','m')
   nradial_grid  = 20
   nangular_grid = 26 
   quadrature_name = 'medium'
 case('HIGH','high','H','h')
   nradial_grid  = 40
   nangular_grid = 50 
   quadrature_name = 'high'
 case('VERYHIGH','veryhigh','VH','vh')
   nradial_grid  = 80
   nangular_grid = 86 
   quadrature_name = 'very high'
 case default
   stop'integration quality not recognized'
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
 use m_atoms
 use m_basis_set
 implicit none
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer,parameter          :: nr=20000
 real(dp),parameter         :: length=10.0_dp
 integer                    :: ir,ibf
 integer                    :: istate1,istate2,istate,ispin
 real(dp)                   :: rr(3)
 real(dp),allocatable       :: phi(:,:),phase(:,:)
 real(dp)                   :: u(3),a(3)
 logical                    :: file_exists
 real(dp)                   :: xmin,xmax
 real(dp)                   :: basis_function_r(basis%nbf)
 integer                    :: ibf_cart,ni_cart,ni,li,i_cart
 real(dp),allocatable       :: basis_function_r_cart(:)
!=====

 WRITE_MASTER(*,*) 
 WRITE_MASTER(*,*) 'Plotting some selected wavefunctions'
 inquire(file='manual_plotwfn',exist=file_exists)
 if(file_exists) then
   open(100,file='manual_plotwfn',status='old')
   read(100,*) istate1,istate2
   read(100,*) u(:)
   read(100,*) a(:)
   close(100)
 else
   istate1=1
   istate2=2
   u(:)=0.0_dp
   u(1)=1.0_dp
   a(:)=0.0_dp
 endif
 u(:) = u(:) / SQRT(SUM(u(:)**2))
 allocate(phase(istate1:istate2,nspin),phi(istate1:istate2,nspin))
 WRITE_MASTER(*,'(a,2(2x,i4))')   ' states:   ',istate1,istate2
 WRITE_MASTER(*,'(a,3(2x,f8.3))') ' direction:',u(:)
 WRITE_MASTER(*,'(a,3(2x,f8.3))') ' origin:   ',a(:)

 xmin = MINVAL( u(1)*x(1,:) + u(2)*x(2,:) + u(3)*x(3,:) ) - length
 xmax = MAXVAL( u(1)*x(1,:) + u(2)*x(2,:) + u(3)*x(3,:) ) + length

 phase(:,:)=1.0_dp

 do ir=1,nr
   rr(:) = ( xmin + (ir-1)*(xmax-xmin)/REAL(nr-1,dp) ) * u(:) + a(:)

   phi(:,:) = 0.0_dp
   
   !
   ! First precalculate all the needed basis function evaluations at point rr
   !
   ibf_cart = 1
   ibf      = 1
   do while(ibf_cart<=basis%nbf_cart)
     li      = basis%bf(ibf_cart)%am
     ni_cart = number_basis_function_am(CARTESIAN,li)
     ni      = number_basis_function_am(basis%gaussian_type,li)

     allocate(basis_function_r_cart(ni_cart))

     do i_cart=1,ni_cart
       basis_function_r_cart(i_cart) = eval_basis_function(basis%bf(ibf_cart+i_cart-1),rr)
     enddo
     basis_function_r(ibf:ibf+ni-1) = MATMUL(  basis_function_r_cart(:) , cart_to_pure(li)%matrix(:,:) )
     deallocate(basis_function_r_cart)

     ibf      = ibf      + ni
     ibf_cart = ibf_cart + ni_cart
   enddo
   !
   ! Precalculation done!
   !

   do ispin=1,nspin
     phi(istate1:istate2,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,istate1:istate2,ispin) )
   enddo

   !
   ! turn the wfns so that they are all positive at a given point
   if(ir==1) then
     do ispin=1,nspin
       do istate=istate1,istate2
         if( phi(istate,ispin) < 0.0_dp ) phase(istate,ispin) = -1.0_dp
       enddo
     enddo
   endif

   WRITE_MASTER(101,'(50(e16.8,2x))') DOT_PRODUCT(rr(:),u(:)),phi(:,:)*phase(:,:)
   WRITE_MASTER(102,'(50(e16.8,2x))') DOT_PRODUCT(rr(:),u(:)),phi(:,:)**2

 enddo

 deallocate(phase,phi)

end subroutine plot_wfn

!=========================================================================
subroutine plot_cube_wfn(nspin,basis,c_matrix)
 use m_definitions
 use m_mpi
 use m_atoms
 use m_basis_set
 implicit none
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer                    :: nx
 integer                    :: ny
 integer                    :: nz
 real(dp),parameter         :: length=4.0_dp
 integer                    :: ir,ibf
 integer                    :: istate1,istate2,istate,ispin
 real(dp)                   :: rr(3)
 real(dp),allocatable       :: phi(:,:),phase(:,:)
 real(dp)                   :: u(3),a(3)
 logical                    :: file_exists
 real(dp)                   :: xmin,xmax,ymin,ymax,zmin,zmax
 real(dp)                   :: basis_function_r(basis%nbf)
 integer                    :: ix,iy,iz,iatom
 integer                    :: ibf_cart,ni_cart,ni,li,i_cart
 real(dp),allocatable       :: basis_function_r_cart(:)
 integer                    :: file_unit
 character(len=200)         :: file_name
!=====

 WRITE_MASTER(*,*) 
 WRITE_MASTER(*,*) 'Plotting some selected wavefunctions in a cube file'
 inquire(file='manual_cubewfn',exist=file_exists)
 if(file_exists) then
   open(100,file='manual_cubewfn',status='old')
   read(100,*) istate1,istate2
   read(100,*) nx,ny,nz
   close(100)
 else
   istate1=1
   istate2=2
   nx=40
   ny=40
   nz=40
 endif
 allocate(phase(istate1:istate2,nspin),phi(istate1:istate2,nspin))
 WRITE_MASTER(*,'(a,2(2x,i4))')   ' states:   ',istate1,istate2

 xmin = MINVAL( x(1,:) ) - length
 xmax = MAXVAL( x(1,:) ) + length
 ymin = MINVAL( x(2,:) ) - length
 ymax = MAXVAL( x(2,:) ) + length
 zmin = MINVAL( x(3,:) ) - length
 zmax = MAXVAL( x(3,:) ) + length

 do istate=istate1,istate2
   do ispin=1,nspin
     file_unit=1000+istate-istate1+(ispin-1)*(istate2-istate1+1)
     WRITE_ME(file_name,'(a,i3.3,a,i1,a)') 'wfn_',istate,'_',ispin,'.cube'
     open(unit=file_unit,file=file_name)
     WRITE_MASTER(file_unit,'(a)') 'cube file generated from MOLGW'
     WRITE_MASTER(file_unit,'(a,i4)') 'wavefunction ',istate1
     WRITE_MASTER(file_unit,'(i6,3(f12.6,2x))') natom,xmin,ymin,zmin
     WRITE_MASTER(file_unit,'(i6,3(f12.6,2x))') nx,(xmax-xmin)/REAL(nx,dp),0.,0.
     WRITE_MASTER(file_unit,'(i6,3(f12.6,2x))') ny,0.,(ymax-ymin)/REAL(ny,dp),0.
     WRITE_MASTER(file_unit,'(i6,3(f12.6,2x))') nz,0.,0.,(zmax-zmin)/REAL(nz,dp)
     do iatom=1,natom
       WRITE_MASTER(file_unit,'(i6,4(2x,f12.6))') NINT(zatom(iatom)),0.0,x(:,iatom)
     enddo
   enddo
 enddo

 phase(:,:)=1.0_dp

 do ix=1,nx
   rr(1) = ( xmin + (ix-1)*(xmax-xmin)/REAL(nx,dp) ) 
   do iy=1,ny
     rr(2) = ( ymin + (iy-1)*(ymax-ymin)/REAL(ny,dp) ) 
     do iz=1,nz
       rr(3) = ( zmin + (iz-1)*(zmax-zmin)/REAL(nz,dp) ) 


       phi(:,:) = 0.0_dp
       
       !
       ! First precalculate all the needed basis function evaluations at point rr
       !
       ibf_cart = 1
       ibf      = 1
       do while(ibf_cart<=basis%nbf_cart)
         li      = basis%bf(ibf_cart)%am
         ni_cart = number_basis_function_am(CARTESIAN,li)
         ni      = number_basis_function_am(basis%gaussian_type,li)
    
         allocate(basis_function_r_cart(ni_cart))
    
         do i_cart=1,ni_cart
           basis_function_r_cart(i_cart) = eval_basis_function(basis%bf(ibf_cart+i_cart-1),rr)
         enddo
         basis_function_r(ibf:ibf+ni-1) = MATMUL(  basis_function_r_cart(:) , cart_to_pure(li)%matrix(:,:) )
         deallocate(basis_function_r_cart)
    
         ibf      = ibf      + ni
         ibf_cart = ibf_cart + ni_cart
       enddo
       !
       ! Precalculation done!
       !

       do ispin=1,nspin
         phi(istate1:istate2,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,istate1:istate2,ispin) )
       enddo

       !
       ! turn the wfns so that they are all positive at a given point
       if(ir==1) then
         do ispin=1,nspin
           do istate=istate1,istate2
             if( phi(istate,ispin) < 0.0_dp ) phase(istate,ispin) = -1.0_dp
           enddo
         enddo
       endif

       do istate=istate1,istate2
         do ispin=1,nspin
           file_unit=1000+istate-istate1+(ispin-1)*(istate2-istate1+1)
           WRITE_MASTER(file_unit,'(50(e16.8,2x))') phi(istate,ispin)*phase(istate,ispin)
         enddo
       enddo

     enddo
   enddo
 enddo

 deallocate(phase,phi)

 do istate=istate1,istate2
   do ispin=1,nspin
     file_unit=1000+istate-istate1+(ispin-1)*(istate2-istate1+1)
     close(file_unit)
   enddo
 enddo

end subroutine plot_cube_wfn

!=========================================================================
