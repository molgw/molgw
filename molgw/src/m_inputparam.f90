!=========================================================================
#include "macros.h"
!=========================================================================
module m_inputparam
 use m_definitions
 use m_mpi
 use m_atoms
 use m_basis_set
 use m_scf
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif

 !
 ! Method definitions
 integer,parameter :: perturbative = 101
 integer,parameter :: QS           = 102
 integer,parameter :: COHSEX       = 103
 integer,parameter :: QSCOHSEX     = 104

 type calculation_type
   character(len=100) :: calc_name
   character(len=100) :: scf_name
   character(len=100) :: postscf_name
   logical            :: is_dft
   logical            :: need_exchange
   logical            :: need_rpa
   logical            :: is_lr_mbpt
   logical            :: is_screened_hybrid
   logical            :: is_gw
   logical            :: is_mp2
   logical            :: is_ci
   logical            :: read_potential
   logical            :: is_bse,is_td,is_tda
   integer            :: gwmethod                    ! perturbative or quasiparticle self-consistent
 end type calculation_type

 integer,protected                :: nspin
 real(dp),protected               :: spin_fact
 integer,protected                :: nscf
 real(dp),protected               :: alpha_mixing
 character(len=100),protected     :: basis_name
 character(len=100),protected     :: auxil_basis_name
 integer,protected                :: gaussian_type
 real(dp),protected               :: electrons
 real(dp),protected               :: magnetization
 type(calculation_type),protected :: calc_type
 character(len=100),protected     :: quadrature_name
 logical,protected                :: is_auxil_basis

 logical,protected                :: print_matrix
 logical,protected                :: print_basis
 logical,protected                :: print_eri
 logical,protected                :: ignore_big_restart
 logical,protected                :: print_wfn
 logical,protected                :: print_specfunc

 real(dp)                         :: alpha_hybrid    = 1.0_dp
 real(dp)                         :: alpha_hybrid_lr = 0.0_dp
 real(dp)                         :: rcut            = 0.0_dp

 integer,protected                :: ndft_xc      = 0
 integer,protected,allocatable    :: dft_xc_type(:)
 real(dp),protected,allocatable   :: dft_xc_coef(:)


contains


!=========================================================================
subroutine init_calculation_type(calc_type,input_key)
 implicit none
!=====
 type(calculation_type),intent(out)   :: calc_type
 character(len=100),intent(in)        :: input_key
!=====
 integer                              :: ipos
 character(len=100)                   :: key1,key2
!=====

 msg='calculation name: '//TRIM(input_key)
 call issue_warning(msg)
 !
 ! default values
 calc_type%calc_name           =  TRIM(input_key)
 calc_type%is_dft              = .FALSE.
 calc_type%need_exchange       = .FALSE.
 calc_type%need_rpa            = .FALSE.
 calc_type%is_lr_mbpt          = .FALSE.
 calc_type%is_screened_hybrid  = .FALSE.
 calc_type%is_gw               = .FALSE.
 calc_type%is_mp2              = .FALSE.
 calc_type%is_ci               = .FALSE.
 calc_type%is_bse              = .FALSE.
 calc_type%is_td               = .FALSE.
 calc_type%is_tda              = .FALSE.
 calc_type%gwmethod            = 0
 calc_type%read_potential      = .FALSE.
 calc_type%postscf_name        = 'None'
 

 ipos=index(input_key,'+',.TRUE.)

 key1=''
 key2=''

 !
 ! If it exists, first read the last part of the calculation specifier
 if(ipos/=0) then
   key1(:ipos-1) = input_key(:ipos-1)
   key2(1:) = input_key(ipos+1:)

   calc_type%postscf_name =  TRIM(key2)

   select case(TRIM(key2))
   case('GW')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = perturbative
   case('COHSEX')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = COHSEX
   case('LRGW')
     calc_type%is_gw      =.TRUE.
     calc_type%gwmethod   = perturbative
     calc_type%is_lr_mbpt = .TRUE.
   case('MP2')
     calc_type%is_mp2   =.TRUE.
     calc_type%gwmethod = perturbative
   case('CI')
     calc_type%is_ci =.TRUE.
   case('BSE')
     calc_type%is_bse =.TRUE.
   case('BSE-TDA')
     calc_type%is_bse =.TRUE.
     calc_type%is_tda =.TRUE.
   case('TD')
     calc_type%is_td =.TRUE.
   case('TDA')
     calc_type%is_td  =.TRUE.
     calc_type%is_tda =.TRUE.
   case default
     stop'error reading calculation type part 2'
   end select
 else
   key1 = input_key
 endif

 calc_type%scf_name =  TRIM(key1)

 !
 ! Then read the first part of the calculation specifier
 select case(TRIM(key1))
 case('CI')
   calc_type%is_ci =.TRUE.
   calc_type%need_exchange = .TRUE.
 case('HF')
   calc_type%need_exchange = .TRUE.  
 case('MP2')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_mp2        = .TRUE.
   calc_type%gwmethod      = QS
 case('GW')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_gw         = .TRUE.
   calc_type%gwmethod      = QS
 case('COHSEX')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_gw         = .TRUE.
   calc_type%gwmethod      = QSCOHSEX
 case('VIN')
   calc_type%read_potential= .TRUE.  
 case default
   !
   ! If the calculation type is none of the above, let's assume it is DFT-type
   calc_type%is_dft=.TRUE.
   call init_dft_type(key1,calc_type)
 end select

end subroutine init_calculation_type


!=========================================================================
subroutine init_dft_type(key,calc_type)
 implicit none
!=====
 character(len=100),intent(in)          :: key
 type(calculation_type),intent(inout)   :: calc_type
!=====


 select case(TRIM(key))
 case('LDAx','PBEx','PBEhx','Bx','PW91x','BJx','RPPx',&
      'BHANDH','BHANDHLYP','B3LYP','PBE0','HSE03','HSE06','HSE08','HCTH','CAM-B3LYP','TD-CAM-B3LYP')
   ndft_xc=1
 case('LDA','SPL','VWN','VWN_RPA','PBE','PBEh','BLYP','PW91')
   ndft_xc=2
 case('TESTHSE')
   ndft_xc=1
 case('TESTPBE0','TESTLDA0')
   ndft_xc=2
 case default
   WRITE_MASTER(*,*) 'error reading calculation type'
   WRITE_MASTER(*,*) TRIM(key)
   stop' is unknown'
 end select

 allocate(dft_xc_type(ndft_xc))
 allocate(dft_xc_coef(ndft_xc))
 !
 ! default is one, otherwise it is modified later
 dft_xc_coef(:) = 1.0_dp

 select case(TRIM(key))
#ifdef HAVE_LIBXC
 !
 ! LDA functionals
 case('LDAx')
   dft_xc_type(1) = XC_LDA_X
 case('SPL')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PZ
 case('LDA')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PW
 case('VWN')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_VWN
 case('VWN_RPA')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_VWN_RPA
 !
 ! GGA functionals
 case('PBEx')
   dft_xc_type(1) = XC_GGA_X_PBE
 case('PBE')
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_C_PBE
 case('PBEhx')
   dft_xc_type(1) = XC_GGA_X_WPBEH
 case('PBEh')
   dft_xc_type(1) = XC_GGA_X_WPBEH
   dft_xc_type(2) = XC_GGA_C_PBE
 case('Bx')
   dft_xc_type(1) = XC_GGA_X_B88
 case('BLYP')
   dft_xc_type(1) = XC_GGA_X_B88
   dft_xc_type(2) = XC_GGA_C_LYP
 case('PW91x')
   dft_xc_type(1) = XC_GGA_X_PW91
 case('PW91')
   dft_xc_type(1) = XC_GGA_X_PW91
   dft_xc_type(2) = XC_GGA_C_PW91
 case('HCTH')
   dft_xc_type(1) = XC_GGA_XC_HCTH_407
 case('TH')
   dft_xc_type(1) = XC_GGA_XC_TH1
 !
 ! Meta-GGA functionals
 case('BJx')
   dft_xc_type(1) = XC_MGGA_X_BJ06
 case('RPPx')
   dft_xc_type(1) = XC_MGGA_X_RPP09
 !
 ! Hybrid functionals
 case('BHANDH')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_BHANDH
   alpha_hybrid = 0.50_dp
 case('BHANDHLYP')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_BHANDHLYP
   alpha_hybrid = 0.50_dp
 case('B3LYP')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_B3LYP
   alpha_hybrid = 0.20_dp
 case('PBE0')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_PBEH
   alpha_hybrid = 0.25_dp
 case('HSE03')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_HSE03
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   rcut         = 1.0_dp / ( 0.15_dp / SQRT(2.0_dp) )
 case('HSE06')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_HSE06
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   rcut         = 1.0_dp / 0.11_dp
 case('HSE08')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.
   dft_xc_type(1) = XC_HYB_GGA_XC_HJS_PBE
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   rcut           = 1.0_dp / 0.11_dp
 case('CAM-B3LYP')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.
   dft_xc_type(1)  = XC_HYB_GGA_XC_CAM_B3LYP
   alpha_hybrid    =  0.19_dp 
   alpha_hybrid_lr =  0.46_dp 
   rcut            =  1.0_dp / 0.33_dp  
 case('TD-CAM-B3LYP')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.
   dft_xc_type(1)  = XC_HYB_GGA_XC_TUNED_CAM_B3LYP
   alpha_hybrid    =  0.0799_dp 
   alpha_hybrid_lr =  0.9201_dp
   rcut            =  1.0_dp / 0.150_dp  
 ! Testing
 case('TESTHSE')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_GGA_X_wPBEh
   dft_xc_coef(2) = -0.25_dp
   alpha_hybrid   = 0.25_dp
   rcut           = 1.0_dp / 0.11_dp
 case('TESTLDA0')
   calc_type%need_exchange       = .TRUE.
   alpha_hybrid   = 0.25_dp
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PW
   dft_xc_coef(1) =  1.00_dp - alpha_hybrid
   dft_xc_coef(2) =  1.00_dp
 case('TESTPBE0')
   calc_type%need_exchange       = .TRUE.
   alpha_hybrid   = 0.25_dp
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_C_PBE
   dft_xc_coef(1) =  1.00_dp - alpha_hybrid
   dft_xc_coef(2) =  1.00_dp
#endif
 case default
   stop'error reading calculation type part 1'
 end select


end subroutine init_dft_type


!=========================================================================
subroutine output_calculation_type(calc_type)
 implicit none
 type(calculation_type),intent(in)   :: calc_type
!=====

  WRITE_MASTER(*,*) 'need_exchange               ',calc_type%need_exchange
  WRITE_MASTER(*,*) 'need_rpa                    ',calc_type%need_rpa
  WRITE_MASTER(*,*) 'is_gw                       ',calc_type%is_gw
  WRITE_MASTER(*,*) 'is_mp2                      ',calc_type%is_mp2
  WRITE_MASTER(*,*) 'is_ci                       ',calc_type%is_ci
  WRITE_MASTER(*,*) 'is_td                       ',calc_type%is_td
  WRITE_MASTER(*,*) 'is_tda                      ',calc_type%is_tda
  WRITE_MASTER(*,*) 'is_bse                      ',calc_type%is_bse
  WRITE_MASTER(*,*) 'method                      ',calc_type%gwmethod  

end subroutine output_calculation_type


!=========================================================================
subroutine read_inputparameter_molecule()
! use m_tools
 implicit none

 character(len=100)                 :: read_char
 character(len=200)                 :: read_line
 character(len=100)                 :: line_wocomment
 character(len=100)                 :: mixing_name
 character(len=100)                 :: quadrature_accuracy
 character(len=100)                 :: field1,field2,field3
 integer                            :: print_volume
 integer                            :: ipos,jpos,kpos
 integer                            :: istat,iatom,jatom
 real(dp)                           :: charge,length_factor
!=====

 !
 ! Reading line 1
 read(*,*) read_char
 call init_calculation_type(calc_type,read_char)

 !
 ! Reading line 2
 read(*,*) nspin,charge,magnetization
 if(nspin/=1 .AND. nspin/=2) stop'nspin in incorrect'
 if(magnetization<-1.d-5)    stop'magnetization is negative'
 if(magnetization>1.d-5 .AND. nspin==1) stop'magnetization is non-zero and nspin is 1'

 spin_fact = REAL(-nspin+3,dp)

 !
 ! Reading line 3
 ! 2 mandatory fields + 1 optional field
 read(*,fmt='(a200)') read_line 
 ipos=index(read_line,'#',back=.false.)
 if(ipos==0) ipos=201
 line_wocomment(:ipos-1)=read_line(:ipos-1)
 line_wocomment=ADJUSTL(line_wocomment)
 jpos=index(line_wocomment,' ',back=.false.)
 basis_name = line_wocomment(:jpos-1)

 line_wocomment(1:ipos-jpos-1)=line_wocomment(jpos+1:ipos-1)
 line_wocomment(ipos-jpos:)=' '
 line_wocomment=ADJUSTL(line_wocomment)

 kpos=index(line_wocomment,' ',back=.false.)
 field2=line_wocomment(1:kpos-1)
 select case(TRIM(field2))
 case('PURE','pure')
   gaussian_type=PURE
 case('CART','cart')
   gaussian_type=CARTESIAN
 case default
   stop'Error in input line 3: second keyword should either PURE or CART'
 end select

 ! Optional field: auxiliary basis
 line_wocomment(1:ipos-kpos-1)=line_wocomment(kpos+1:ipos-1)
 line_wocomment(ipos-kpos:)=' '
 line_wocomment=ADJUSTL(line_wocomment)
 kpos=index(line_wocomment,' ',back=.false.)
 if(kpos==0 .OR. kpos==1) then
   auxil_basis_name = 'none'
 else
   auxil_basis_name = line_wocomment(:kpos-1)
 endif

 if( auxil_basis_name=='none' .OR. auxil_basis_name=='NONE' ) then
   auxil_basis_name='none'
   is_auxil_basis = .FALSE.
 else
   is_auxil_basis = .TRUE.
 endif
 


 !
 ! Reading line 4
 read(*,*) nscf,alpha_mixing,mixing_name,quadrature_accuracy
 if(nscf<1) stop'nscf too small'
 if(alpha_mixing<0.0 .OR. alpha_mixing > 1.0 ) stop'alpha_mixing should be inside [0,1]'
 select case(TRIM(mixing_name))
 case('SIMPLE')
   mixing_scheme = simple_mixing
 case('PULAY')
   mixing_scheme = pulay_mixing
 case('FBFB')
   mixing_scheme = FBFB_mixing
 case default
   WRITE_MASTER(*,*) TRIM(mixing_name)
   stop'mixing scheme not recognized'
 end select

 select case(TRIM(quadrature_accuracy))
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
 case('INSANE','insane','I','i')
   quadrature_name = 'insane'
 case default
   stop'integration quality not recognized'
 end select

 !
 ! Reading line 5
 read(*,*) print_volume

 print_matrix        = MODULO(print_volume       ,2)==1
 print_basis         = MODULO(print_volume/10    ,2)==1
 print_eri           = MODULO(print_volume/100   ,2)==1
 ignore_big_restart  = MODULO(print_volume/1000  ,2)==1
 print_wfn           = MODULO(print_volume/10000 ,2)==1
 print_specfunc      = MODULO(print_volume/100000,2)==1

 
 !
 ! Reading line 6 and following
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
 ! Summarize input parameters
 WRITE_MASTER(*,'(/,a,/)')    ' Summary of the input parameters '
 WRITE_MASTER(*,'(a25,2x,a)') ' Calculation type: ',calc_type%calc_name
 WRITE_MASTER(*,'(a25,2x,a)') '         SCF type: ',calc_type%scf_name
 WRITE_MASTER(*,'(a25,2x,a)') '    Post SCF type: ',calc_type%postscf_name
 WRITE_MASTER(*,'(a25,i3)')   ' Natom: ',natom
 WRITE_MASTER(*,'(a25,f8.4)') ' Electrons: ',electrons
 WRITE_MASTER(*,'(a25,f8.4)') ' Charge: ',charge
 WRITE_MASTER(*,'(a25,f8.4)') ' Magnetization: ',magnetization
 WRITE_MASTER(*,'(a25,2x,a)') ' Basis set: ',basis_name
 WRITE_MASTER(*,'(a25,2x,a)') ' Auxiliary basis set: ',auxil_basis_name
 if(gaussian_type==PURE) then
   WRITE_MASTER(*,'(a25,2x,a)') ' Gaussian type: ','pure'
 else
   WRITE_MASTER(*,'(a25,2x,a)') ' Gaussian type: ','cartesian'
 endif
 WRITE_MASTER(*,'(a25,i3)')   ' Spin polarization: ',nspin
 WRITE_MASTER(*,'(a25,i3)')   ' SCF steps: ',nscf
 WRITE_MASTER(*,'(a25,f8.4)') ' Mixing: ',alpha_mixing
 WRITE_MASTER(*,'(a25,2x,a)') ' Quadrature accuracy: ',quadrature_name
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a19)')      ' IO options:'
 WRITE_MASTER(*,'(a30,l3)')   ' - matrices details:   ',print_matrix        
 WRITE_MASTER(*,'(a30,l3)')   ' - basis set details:  ',print_basis
 WRITE_MASTER(*,'(a30,l3)')   ' - ERI file:           ',print_eri           
 WRITE_MASTER(*,'(a30,l3)')   ' - ignore big RESTART: ',ignore_big_restart
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

