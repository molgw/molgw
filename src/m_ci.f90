!==================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! full ci calculations for a few electrons  (1, 2, 3, 4, 5)
! using occupation number vectors
!
!==================================================================
module m_ci
 use m_definitions
 use m_mpi
 use m_warning
 use m_tools
 use m_memory
 use m_timing
 use m_scalapack
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam


 integer,private              :: nfrozen_ci
 integer,private              :: nstate_ci

 real(dp),allocatable,private :: h_1body(:,:)

 real(dp),allocatable,private :: energy_1e(:)
 real(dp),allocatable,private :: energy_2e(:)
 real(dp),allocatable,private :: energy_3e(:)
 real(dp),allocatable,private :: energy_4e(:)
 real(dp),allocatable,private :: energy_5e(:)

 real(dp),allocatable,private :: eigvec_1e(:,:)
 real(dp),allocatable,private :: eigvec_2e(:,:)
 real(dp),allocatable,private :: eigvec_3e(:,:)
 real(dp),allocatable,private :: eigvec_4e(:,:)
 real(dp),allocatable,private :: eigvec_5e(:,:)

 type, private :: configurations
   integer             :: nelec
   integer             :: nconf
   integer             :: sz
   integer,allocatable :: sporb_occ(:,:)        ! spinor-orbitals with occupation equal to 1
 end type

 type(configurations),private :: conf_1e
 type(configurations),private :: conf_2e
 type(configurations),private :: conf_3e
 type(configurations),private :: conf_4e
 type(configurations),private :: conf_5e


contains 


!==================================================================
! Sign introduced by creation/annihilation operators
! as defined in Hellgaker's book (chapter 1 box 1)
pure function gamma_sign(on,isporb)
 implicit none

 integer,intent(in) :: on(:)
 integer,intent(in) :: isporb
 integer            :: gamma_sign
!=====
!=====

 gamma_sign = 1 - 2 * MODULO( COUNT( on(1:isporb-1) == 1 ) , 2 ) 

end function gamma_sign


!==================================================================
! From spin-orbital index, get the spin = +1 or -1
elemental function sporb_to_spin(isporb) result(ispin)
 implicit none

 integer,intent(in) :: isporb
 integer :: ispin
!=====
!=====

 ispin = 2 * MODULO( isporb , 2 ) - 1

end function sporb_to_spin


!==================================================================
! From spin-orbital index, get the orbital index 
elemental function sporb_to_state(isporb) result(istate)
 implicit none

 integer,intent(in) :: isporb
 integer :: istate
!=====
!=====

 istate = ( isporb + 1 ) / 2

end function sporb_to_state


!==================================================================
! From occupied spin-orbitals to an occupation-number vector
pure function sporb_to_on(sporb_occ) result(on)
 implicit none

 integer,intent(in) :: sporb_occ(:)
 integer :: on(2*nstate_ci)
!=====
 integer :: ielec
!=====

 on(:) = 0
 do ielec=1,SIZE(sporb_occ)
   on(sporb_occ(ielec)) = 1
 enddo

end function sporb_to_on


!==================================================================
subroutine prepare_ci(nstate_in,nfrozen_in,h_1e,c_matrix)
 implicit none

 integer,intent(in)  :: nstate_in,nfrozen_in
 real(dp),intent(in) :: h_1e(:,:),c_matrix(:,:,:)
!=====
!=====

 nstate_ci  = nstate_in
 nfrozen_ci = nfrozen_in

 write(stdout,'(/,1x,a,i4,a,i4)') 'Prepare CI with active states ranging from ',nfrozen_ci+1,' to ',nstate_ci

 ! Calculate the one-electron hamiltonian on the eigenstate basis
 call build_1e_hamiltonian(c_matrix,h_1e)


end subroutine prepare_ci


!==================================================================
subroutine setup_configurations_ci(nelec,spinstate,conf)
 implicit none

 integer,intent(in)               :: nelec
 integer,intent(in)               :: spinstate
 type(configurations),intent(out) :: conf
!=====
 integer :: iconf
 integer :: ispin(nelec),istate(nelec)
 integer :: sporb(nelec)
 integer :: isporb,jsporb,ksporb,lsporb,msporb
!=====

 conf%nelec = nelec
 conf%sz    = spinstate

 select case(conf%sz)
 case(-100)
   write(stdout,*) 'Any spin state'
 case(0)
   write(stdout,*) 'Spin singlet'
 case(1)
   write(stdout,*) 'Spin doublet'
 case(2)
   write(stdout,*) 'Spin triplet'
 case(3)
   write(stdout,*) 'Spin quadruplet'
 case default
   call die('setup_configurations_ci: spin case not recognized')
 end select

 select case(conf%nelec)
 case(1)
   conf%nconf = 0
   do isporb=2*nfrozen_ci+1,2*nstate_ci
     sporb(1)  = isporb
     ispin(:)  = sporb_to_spin(  sporb(:) )
     istate(:) = sporb_to_state( sporb(:) )
     if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
       conf%nconf = conf%nconf + 1
     endif
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do isporb=2*nfrozen_ci+1,2*nstate_ci
     sporb(1)  = isporb
     ispin(:)  = sporb_to_spin(  sporb(:) )
     istate(:) = sporb_to_state( sporb(:) )
     if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
       iconf = iconf + 1
       conf%sporb_occ(:,iconf) = sporb(:)
     endif
   enddo


 case(2)
   conf%nconf = 0
   do jsporb=2*nfrozen_ci+1,2*nstate_ci
     do isporb=jsporb+1,2*nstate_ci
       sporb(1)  = isporb
       sporb(2)  = jsporb
       ispin(:)  = sporb_to_spin(  sporb(:) )
       istate(:) = sporb_to_state( sporb(:) )
       if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
         conf%nconf = conf%nconf + 1
       endif
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do jsporb=2*nfrozen_ci+1,2*nstate_ci
     do isporb=jsporb+1,2*nstate_ci
       sporb(1)  = isporb
       sporb(2)  = jsporb
       ispin(:)  = sporb_to_spin( sporb(:) )
       istate(:) = sporb_to_state( sporb(:) )
       if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
         iconf = iconf + 1
         conf%sporb_occ(:,iconf) = sporb(:)
       endif
     enddo
   enddo


 case(3)
   conf%nconf = 0
   do ksporb=2*nfrozen_ci+1,2*nstate_ci
     do jsporb=ksporb+1,2*nstate_ci
       do isporb=jsporb+1,2*nstate_ci
         sporb(1)  = isporb
         sporb(2)  = jsporb
         sporb(3)  = ksporb
         ispin(:)  = sporb_to_spin(  sporb(:) )
         istate(:) = sporb_to_state( sporb(:) )
         if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
           conf%nconf = conf%nconf + 1
         endif
       enddo
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do ksporb=2*nfrozen_ci+1,2*nstate_ci
     do jsporb=ksporb+1,2*nstate_ci
       do isporb=jsporb+1,2*nstate_ci
         sporb(1)  = isporb
         sporb(2)  = jsporb
         sporb(3)  = ksporb
         ispin(:)  = sporb_to_spin( sporb(:) )
         istate(:) = sporb_to_state( sporb(:) )
         if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
           iconf = iconf + 1
           conf%sporb_occ(:,iconf) = sporb(:)
         endif
       enddo
     enddo
   enddo

 case(4)
   conf%nconf = 0
   do lsporb=2*nfrozen_ci+1,2*nstate_ci
     do ksporb=lsporb+1,2*nstate_ci
       do jsporb=ksporb+1,2*nstate_ci
         do isporb=jsporb+1,2*nstate_ci
           sporb(1)  = isporb
           sporb(2)  = jsporb
           sporb(3)  = ksporb
           sporb(4)  = lsporb
           ispin(:)  = sporb_to_spin(  sporb(:) )
           istate(:) = sporb_to_state( sporb(:) )
           if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
             conf%nconf = conf%nconf + 1
           endif
         enddo
       enddo
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do lsporb=2*nfrozen_ci+1,2*nstate_ci
     do ksporb=lsporb+1,2*nstate_ci
       do jsporb=ksporb+1,2*nstate_ci
         do isporb=jsporb+1,2*nstate_ci
           sporb(1)  = isporb
           sporb(2)  = jsporb
           sporb(3)  = ksporb
           sporb(4)  = lsporb
           ispin(:)  = sporb_to_spin( sporb(:) )
           istate(:) = sporb_to_state( sporb(:) )
           if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
             iconf = iconf + 1
             conf%sporb_occ(:,iconf) = sporb(:)
           endif
         enddo
       enddo
     enddo
   enddo

 case(5)
   conf%nconf = 0
   do msporb=2*nfrozen_ci+1,2*nstate_ci
     do lsporb=msporb+1,2*nstate_ci
       do ksporb=lsporb+1,2*nstate_ci
         do jsporb=ksporb+1,2*nstate_ci
           do isporb=jsporb+1,2*nstate_ci
             sporb(1)  = isporb
             sporb(2)  = jsporb
             sporb(3)  = ksporb
             sporb(4)  = lsporb
             sporb(5)  = msporb
             ispin(:)  = sporb_to_spin(  sporb(:) )
             istate(:) = sporb_to_state( sporb(:) )
             if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
               conf%nconf = conf%nconf + 1
             endif
           enddo
         enddo
       enddo
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do msporb=2*nfrozen_ci+1,2*nstate_ci
     do lsporb=msporb+1,2*nstate_ci
       do ksporb=lsporb+1,2*nstate_ci
         do jsporb=ksporb+1,2*nstate_ci
           do isporb=jsporb+1,2*nstate_ci
             sporb(1)  = isporb
             sporb(2)  = jsporb
             sporb(3)  = ksporb
             sporb(4)  = lsporb
             sporb(5)  = msporb
             ispin(:)  = sporb_to_spin( sporb(:) )
             istate(:) = sporb_to_state( sporb(:) )
             if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
               iconf = iconf + 1
               conf%sporb_occ(:,iconf) = sporb(:)
             endif
           enddo
         enddo
       enddo
     enddo
   enddo


 case default
   call die('setup_configurations_ci: number of electrons not coded as of today')
 end select

 write(stdout,'(/,1x,a,i2,a,i6)') 'Electron count:',conf%nelec,'           Configurations: ',conf%nconf


end subroutine setup_configurations_ci


!==================================================================
subroutine build_1e_hamiltonian(c_matrix,h_1e)
 implicit none

 real(dp),intent(in) :: c_matrix(:,:,:)
 real(dp),intent(in) :: h_1e(:,:)
!=====
!=====

 if( SIZE(c_matrix,DIM=2) /= nstate_ci        ) call die('build_1e_hamiltonian: nstate_ci not initialized yet')
 if( SIZE(c_matrix,DIM=1) /= SIZE(h_1e,DIM=2) ) call die('build_1e_hamiltonian: c_matrix and h_1e have inconsistent sizes')
 if( SIZE(h_1e,DIM=1)     /= SIZE(h_1e,DIM=2) ) call die('build_1e_hamiltonian: h_1e is not square')

 if( ALLOCATED(h_1body) ) return

 write(stdout,*) 'Obtain the one-electron Hamiltonian in the HF basis'

 allocate(h_1body(nstate_ci,nstate_ci))
 h_1body(:,:) = MATMUL( TRANSPOSE(c_matrix(:,:,1)) , MATMUL( h_1e(:,:) , c_matrix(:,:,1) ) )


end subroutine build_1e_hamiltonian


!==================================================================
subroutine build_ci_hamiltonian(conf,h_ci)
 implicit none

 type(configurations),intent(in) :: conf
 real(dp)                        :: h_ci(conf%nconf,conf%nconf)
!=====
 integer :: ielec,jelec
 integer :: iconf,jconf
 integer :: iisporb(conf%nelec),jjsporb(conf%nelec)
 integer :: iistate(conf%nelec),jjstate(conf%nelec)
 integer :: iispin(conf%nelec),jjspin(conf%nelec)
 integer :: on_i(2*nstate_ci),on_j(2*nstate_ci)
 logical :: mask(2*nstate_ci)
 integer :: isporb,jsporb,istate,jstate,ispin,jspin
 integer :: ksporb,lsporb,kstate,lstate,kspin,lspin
!=====

 call start_clock(timing_tmp1)
 !
 ! Follow the second-quantization notations from Hellgaker book Chapter 1.
 ! Use occupation number vectors on_i(:) and on_j(:) filled with 0's and three 1's.
 !

 h_ci(:,:) = 0.0_dp

 do jconf=1,conf%nconf
   jjsporb(:) = conf%sporb_occ(:,jconf)
   jjspin(:)  = sporb_to_spin(  jjsporb(:) )
   jjstate(:) = sporb_to_state( jjsporb(:) )

   on_j(:) = sporb_to_on(jjsporb)


   ! Only fills the upper triangle
   do iconf=1,jconf
     iisporb(:) = conf%sporb_occ(:,iconf)
     iispin(:)  = sporb_to_spin(  iisporb(:) )
     iistate(:) = sporb_to_state( iisporb(:) )

     on_i(:) = sporb_to_on(iisporb)


     !
     ! 1-body part
     !

     !
     ! Exact same ON-vector
     if( iconf == jconf ) then 
       do ielec=1,conf%nelec
         h_ci(iconf,jconf) = h_ci(iconf,jconf) + h_1body(iistate(ielec),iistate(ielec))
       enddo
     endif

     !
     ! ON-vectors differ by one occupation number
     if( COUNT( on_j(:) - on_i(:) == 1 ) == 1 ) then
       jsporb = MAXLOC( on_j(:) - on_i(:) , DIM=1 )
       isporb = MINLOC( on_j(:) - on_i(:) , DIM=1 )
       istate = sporb_to_state( isporb )
       jstate = sporb_to_state( jsporb )

       if( sporb_to_spin( isporb ) == sporb_to_spin( jsporb ) ) &
         h_ci(iconf,jconf) = h_ci(iconf,jconf)  &
                        + h_1body(istate,jstate) * gamma_sign(on_j,jsporb) * gamma_sign(on_i,isporb)

     endif


     !
     ! 2-body part
     !

     !
     ! Exact same ON-vector
     if( iconf == jconf ) then 
       do jelec=1,conf%nelec
         jstate = jjstate(jelec)

         do ielec=1,conf%nelec
           istate = iistate(ielec)
           h_ci(iconf,jconf) = h_ci(iconf,jconf)  &
                      + 0.5_dp * eri_eigen(istate,istate,1,jstate,jstate,1)

           if( iispin(ielec) == jjspin(jelec) )  &
             h_ci(iconf,jconf) = h_ci(iconf,jconf)  &
                        - 0.5_dp * eri_eigen(istate,jstate,1,jstate,istate,1)

         enddo
       enddo
     endif

     !
     ! ON-vectors differ by one occupation number
     if( COUNT( on_j(:) - on_i(:) == 1 ) == 1 ) then
       isporb = MINLOC( on_j(:) - on_i(:) , DIM=1 )
       jsporb = MAXLOC( on_j(:) - on_i(:) , DIM=1 )
       istate = sporb_to_state( isporb )
       jstate = sporb_to_state( jsporb )
       ispin  = sporb_to_spin(  isporb )
       jspin  = sporb_to_spin(  jsporb )

       do ksporb=2*nfrozen_ci+1,2*nstate_ci
         if( on_i(ksporb) == 0 ) cycle
         if( on_j(ksporb) == 0 ) cycle
         kstate = sporb_to_state( ksporb )
         kspin  = sporb_to_spin(  ksporb )

         if( ispin == jspin ) &
           h_ci(iconf,jconf) = h_ci(iconf,jconf)  &
                      + eri_eigen(istate,jstate,1,kstate,kstate,1)   &
                           * gamma_sign(on_i,isporb) * gamma_sign(on_j,jsporb)

         if( ispin == kspin .AND. jspin == kspin ) &
           h_ci(iconf,jconf) = h_ci(iconf,jconf)  &
                      - eri_eigen(istate,kstate,1,kstate,jstate,1)  &
                           * gamma_sign(on_i,isporb) * gamma_sign(on_j,jsporb)
       enddo

     endif

     !
     ! ON-vectors differ by two occupation numbers
     if( COUNT( on_j(:) - on_i(:) == 1 ) == 2 ) then
       ! Find the two indexes k < l
       ksporb = MAXLOC( on_j(:) - on_i(:) , DIM=1 )
       mask(:)      = .TRUE.
       mask(ksporb) = .FALSE.
       lsporb = MAXLOC( on_j(:) - on_i(:) , DIM=1 , MASK=mask)

       ! Find the two indexes i < j
       isporb = MINLOC( on_j(:) - on_i(:) , DIM=1 )
       mask(:)      = .TRUE.
       mask(isporb) = .FALSE.
       jsporb = MINLOC( on_j(:) - on_i(:) , DIM=1 , MASK=mask)

       istate = sporb_to_state( isporb )
       jstate = sporb_to_state( jsporb )
       kstate = sporb_to_state( ksporb )
       lstate = sporb_to_state( lsporb )
       ispin  = sporb_to_spin(  isporb )
       jspin  = sporb_to_spin(  jsporb )
       kspin  = sporb_to_spin(  ksporb )
       lspin  = sporb_to_spin(  lsporb )

       if( ispin == kspin .AND. jspin == lspin ) &
         h_ci(iconf,jconf) = h_ci(iconf,jconf)  &
                    + eri_eigen(istate,kstate,1,jstate,lstate,1)           &
                         * gamma_sign(on_i,isporb) * gamma_sign(on_i,jsporb)  &
                         * gamma_sign(on_j,ksporb) * gamma_sign(on_j,lsporb)

       if( ispin == lspin .AND. jspin == kspin ) &
         h_ci(iconf,jconf) = h_ci(iconf,jconf)  &
                    - eri_eigen(istate,lstate,1,jstate,kstate,1)           &
                         * gamma_sign(on_i,isporb) * gamma_sign(on_i,jsporb)  &
                         * gamma_sign(on_j,ksporb) * gamma_sign(on_j,lsporb)

     endif


     !
     ! Symmetrize here
     h_ci(jconf,iconf) = h_ci(iconf,jconf)

   enddo

 enddo

 write(*,*) '============================='
 write(*,*) sporb_to_on(conf%sporb_occ(:,1))
 write(*,*) h_ci(1,1)
 write(*,*) '============================='
#if 0
 !
 ! Frozen core
 do iconf=1,conf%nconf
   ! one-body part
   do istate=1,nfrozen_ci
     h_ci(iconf,iconf) = h_ci(iconf,iconf) + h_1body(istate,istate) * 2.0_dp
   enddo
   ! two-body part
   do jstate=1,nfrozen_ci
     do istate=1,nfrozen_ci
       h_ci(iconf,iconf) = h_ci(iconf,iconf)  &
                  + 2.0_dp * eri_eigen(istate,istate,1,jstate,jstate,1)

       h_ci(iconf,iconf) = h_ci(iconf,iconf)  &
                  - 1.0_dp * eri_eigen(istate,jstate,1,jstate,istate,1)

     enddo
   enddo
 enddo
 write(*,*) h_ci(1,1)
! stop 'enough'
#endif

 call stop_clock(timing_tmp1)

end subroutine build_ci_hamiltonian


!==================================================================
subroutine full_ci_2electrons_selfenergy(occupation)
 use m_selfenergy_tools
 implicit none

 real(dp),intent(in)        :: occupation(:,:)
!=====
 integer,parameter          :: nomega=5000
 integer,parameter          :: ns=-1
 integer                    :: ielec,jelec
 integer                    :: is
 integer                    :: iconf,jconf
 integer                    :: on_i(2*nstate_ci),on_j(2*nstate_ci)
 integer                    :: on_tmp(2*nstate_ci)
 integer                    :: iisporb(3),iistate(3),iispin(3)
 integer                    :: jjsporb(2),jjstate(2),jjspin(2)
 integer                    :: isporb,istate,ispin
 integer                    :: ns_occ,ns_virt
 real(dp),allocatable       :: fs_occ(:,:),fs_virt(:,:)
 real(dp),allocatable       :: es_occ(:),es_virt(:)
 integer                    :: iomega
 real(dp)                   :: omega
 complex(dp)                :: gi_w
 character(len=3)           :: ctmp3
 character(len=1)           :: ctmp1
 integer                    :: unit_gf
!=====

 call start_clock(timing_ci_selfenergy)

 write(stdout,'(/,1x,a)') 'Full CI self-energy for 2 electrons'


 !
 ! Set the range of states on which to evaluate the self-energy
 call selfenergy_set_state_range(nstate_ci,occupation)

 if( .NOT. ALLOCATED(energy_1e) ) call die('full_ci_2electrons_selfenergy: previous calculation for 1e needed')
 if( .NOT. ALLOCATED(energy_2e) ) call die('full_ci_2electrons_selfenergy: previous calculation for 2e needed')
 if( .NOT. ALLOCATED(energy_3e) ) call die('full_ci_2electrons_selfenergy: previous calculation for 3e needed')

 write(stdout,'(1x,a,sp,i4)') 'Previous CI calculation had spin state Sz(1e): ',conf_1e%sz
 write(stdout,'(1x,a,i3)')    'Previous CI calculation had spin state Sz(2e): ',conf_2e%sz
 write(stdout,'(1x,a,sp,i4)') 'Previous CI calculation had spin state Sz(3e): ',conf_3e%sz

 !
 ! Choose how many Lehmann excitations to calculate
 ! If ns is negative, calculate all of them
 if( ns > 0 ) then
   ns_occ  = ns
   ns_virt = ns
 else
   ns_occ  = conf_1e%nconf
   ns_virt = conf_3e%nconf
 endif

 allocate(fs_occ(2*nstate_ci,ns_occ))
 allocate(fs_virt(2*nstate_ci,ns_virt))
 allocate(es_occ(ns_occ))
 allocate(es_virt(ns_virt))

 !
 ! Build the Lehman amplitude for occupied states
 !
 fs_occ(:,:) = 0.0_dp
 do jconf=1,conf_2e%nconf
   jjsporb(:) = conf_2e%sporb_occ(:,jconf)
   jjspin(:)  = sporb_to_spin(  jjsporb(:) )
   jjstate(:) = sporb_to_state( jjsporb(:) )

   on_j(:) = sporb_to_on(jjsporb)

   do iconf=1,conf_1e%nconf
     iisporb(:1) = conf_1e%sporb_occ(:,iconf)
     iispin(:1)  = sporb_to_spin(  iisporb(:1) )
     iistate(:1) = sporb_to_state( iisporb(:1) )

     on_i(:) = sporb_to_on(iisporb(:1))

     do isporb=2*nfrozen_ci+1,2*nstate_ci
       if( isporb == iisporb(1) ) cycle
       on_tmp(:) = on_i(:)
       on_tmp(isporb) = on_tmp(isporb) + 1
       if( ALL( on_j(:) - on_tmp(:) == 0 ) ) then
         fs_occ(isporb,:) = fs_occ(isporb,:) + eigvec_1e(iconf,:ns_occ) * eigvec_2e(jconf,1) * gamma_sign(on_i,isporb)
       endif
     enddo

   enddo
 enddo

 write(stdout,*) '====================='
 write(stdout,*) 'Occupied states'
 write(stdout,*)

 do is=1,ns_occ
   es_occ(is) = energy_2e(1) - energy_1e(is)
!   write(stdout,'(/,1x,a,i4,1x,f10.6)') '=== Excitation (eV): ',is,es_occ(is) * Ha_eV
  
!   write(stdout,'(1x,a,i2,a)') '=== fs_occ(',is,')'
   do isporb=2*nfrozen_ci,2*nstate_ci
     ispin  = sporb_to_spin(isporb)
     istate = sporb_to_state(isporb)
!     if( ABS(fs_occ(isporb,is)) > 1.0e-4_dp ) write(stdout,*) istate,ispin,fs_occ(isporb,is)
   enddo
  
!   write(stdout,'(1x,a,/)') '========='

 enddo



 !
 ! Build the Lehman amplitude for virtual states
 !
 fs_virt(:,:) = 0.0_dp
 do jconf=1,conf_2e%nconf
   jjsporb(:) = conf_2e%sporb_occ(:,jconf)
   jjspin(:)  = sporb_to_spin(  jjsporb(:) )
   jjstate(:) = sporb_to_state( jjsporb(:) )

   on_j(:) = sporb_to_on(jjsporb)

   do iconf=1,conf_3e%nconf
     iisporb(:3) = conf_3e%sporb_occ(:,iconf)
     iispin(:3)  = sporb_to_spin(  iisporb(:3) )
     iistate(:3) = sporb_to_state( iisporb(:3) )

     on_i(:) = sporb_to_on(iisporb(:3))


     do isporb=2*nfrozen_ci+1,2*nstate_ci
       ! Cannot create an electron in a state that is already occupied
       if( isporb == jjsporb(1) .OR. isporb == jjsporb(2) ) cycle
       on_tmp(:) = on_j(:)
       on_tmp(isporb) = on_tmp(isporb) + 1
       if( ALL( on_i(:) - on_tmp(:) == 0 ) ) then
         fs_virt(isporb,:) = fs_virt(isporb,:) + eigvec_3e(iconf,:ns_virt) * eigvec_2e(jconf,1) * gamma_sign(on_j,isporb)
       endif
     enddo


   enddo
 enddo

 write(stdout,*) '====================='
 write(stdout,*) 'Virtual states'
 write(stdout,*)

 do is=1,ns_virt
   es_virt(is) = energy_3e(is) - energy_2e(1)
!   write(stdout,'(/,1x,a,i4,1x,f10.6)') '=== Excitation (eV): ',is,es_virt(is) * Ha_eV
  
!   write(stdout,'(1x,a,i2,a)') '=== fs_virt(',is,')'
   do isporb=2*nfrozen_ci+1,2*nstate_ci
     ispin  = sporb_to_spin( isporb )
     istate = sporb_to_state( isporb )
!     if( ABS(fs_virt(isporb,is)) > 1.0e-4_dp ) write(stdout,*) istate,ispin,fs_virt(isporb,is)
   enddo
  
!   write(stdout,'(1x,a,/)') '========='

 enddo

 do isporb=2*nsemin-1,2*nsemax
   ispin = sporb_to_spin(isporb)
   istate = sporb_to_state(isporb)

   write(ctmp3,'(i3.3)') istate
   write(ctmp1,'(i1)') MODULO( isporb-1 , 2 ) + 1

   open(newunit=unit_gf,file='exact_greens_function_state'//ctmp3//'_spin'//ctmp1//'.dat',action='write')
   do iomega=-nomega,nomega
     omega = iomega / REAL(nomega,dp) * 100.0_dp / Ha_eV
     gi_w = 0.0_dp

     do is=1,ns_virt
       gi_w = gi_w + fs_virt(isporb,is)**2 / ( omega - es_virt(is) - ieta )
     enddo
     do is=1,ns_occ
       gi_w = gi_w + fs_occ(isporb,is)**2 / ( omega - es_occ(is) + ieta )
     enddo

     write(unit_gf,'(4(1x,es18.8))') omega * Ha_eV, gi_w / Ha_eV, ABS(AIMAG(gi_w)) / pi / Ha_eV
   enddo
   close(unit_gf)
 enddo



 deallocate(fs_virt,fs_occ)
 deallocate(es_virt,es_occ)

 call stop_clock(timing_ci_selfenergy)


end subroutine full_ci_2electrons_selfenergy


!==================================================================
subroutine full_ci_4electrons_selfenergy(occupation)
 use m_selfenergy_tools
 implicit none

 real(dp),intent(in)        :: occupation(:,:)
!=====
 integer,parameter          :: nomega=5000
 integer,parameter          :: ns=-1
 integer                    :: ielec,jelec
 integer                    :: is
 integer                    :: iconf,jconf
 integer                    :: on_i(2*nstate_ci),on_j(2*nstate_ci)
 integer                    :: on_tmp(2*nstate_ci)
 integer                    :: iisporb(5),iistate(5),iispin(5)
 integer                    :: jjsporb(4),jjstate(4),jjspin(4)
 integer                    :: isporb,istate,ispin
 integer                    :: ns_occ,ns_virt
 real(dp),allocatable       :: fs_occ(:,:),fs_virt(:,:)
 real(dp),allocatable       :: es_occ(:),es_virt(:)
 integer                    :: iomega
 real(dp)                   :: omega
 complex(dp)                :: gi_w
 character(len=3)           :: ctmp3
 character(len=1)           :: ctmp1
 integer                    :: unit_gf
!=====

 call start_clock(timing_ci_selfenergy)

 write(stdout,'(/,1x,a)') 'Full CI self-energy for 4 electrons'


 !
 ! Set the range of states on which to evaluate the self-energy
 call selfenergy_set_state_range(nstate_ci,occupation)

 if( .NOT. ALLOCATED(energy_3e) ) call die('full_ci_2electrons_selfenergy: previous calculation for 3e needed')
 if( .NOT. ALLOCATED(energy_4e) ) call die('full_ci_2electrons_selfenergy: previous calculation for 4e needed')
 if( .NOT. ALLOCATED(energy_5e) ) call die('full_ci_2electrons_selfenergy: previous calculation for 5e needed')

 write(stdout,'(1x,a,sp,i4)') 'Previous CI calculation had spin state Sz(3e): ',conf_3e%sz
 write(stdout,'(1x,a,i3)')    'Previous CI calculation had spin state Sz(4e): ',conf_4e%sz
 write(stdout,'(1x,a,sp,i4)') 'Previous CI calculation had spin state Sz(5e): ',conf_5e%sz

 !
 ! Choose how many Lehmann excitations to calculate
 ! If ns is negative, calculate all of them
 if( ns > 0 ) then
   ns_occ  = ns
   ns_virt = ns
 else
   ns_occ  = conf_3e%nconf
   ns_virt = conf_5e%nconf
 endif

 allocate(fs_occ(2*nstate_ci,ns_occ))
 allocate(fs_virt(2*nstate_ci,ns_virt))
 allocate(es_occ(ns_occ))
 allocate(es_virt(ns_virt))

 !
 ! Build the Lehman amplitude for occupied states
 !
 fs_occ(:,:) = 0.0_dp
 do jconf=1,conf_4e%nconf
   jjsporb(:) = conf_4e%sporb_occ(:,jconf)
   jjspin(:)  = sporb_to_spin(  jjsporb(:) )
   jjstate(:) = sporb_to_state( jjsporb(:) )

   on_j(:) = sporb_to_on(jjsporb)

   do iconf=1,conf_3e%nconf
     iisporb(:3) = conf_3e%sporb_occ(:,iconf)
     iispin(:3)  = sporb_to_spin(  iisporb(:3) )
     iistate(:3) = sporb_to_state( iisporb(:3) )

     on_i(:) = sporb_to_on(iisporb(:3))

     !
     ! Evaluate for any s, < 4 , 0 | a_i^+ | 3 , s >
     !
     do isporb=2*nfrozen_ci+1,2*nstate_ci
       if( isporb == iisporb(1) ) cycle
       if( isporb == iisporb(2) ) cycle
       if( isporb == iisporb(3) ) cycle
       on_tmp(:) = on_i(:)
       on_tmp(isporb) = on_tmp(isporb) + 1
       if( ALL( on_j(:) - on_tmp(:) == 0 ) ) then
         fs_occ(isporb,:) = fs_occ(isporb,:) + eigvec_3e(iconf,:ns_occ) * eigvec_4e(jconf,1) * gamma_sign(on_i,isporb)
       endif
     enddo

   enddo
 enddo

 write(stdout,*) '====================='
 write(stdout,*) 'Occupied states'
 write(stdout,*)

 do is=1,ns_occ
   es_occ(is) = energy_4e(1) - energy_3e(is)
   write(stdout,'(/,1x,a,i4,1x,f12.6)') '=== Excitation (eV): ',is,es_occ(is) * Ha_eV
  
!   write(stdout,'(1x,a,i2,a)') '=== fs_occ(',is,')'
   do isporb=2*nfrozen_ci,2*nstate_ci
     ispin  = sporb_to_spin(isporb)
     istate = sporb_to_state(isporb)
!     if( ABS(fs_occ(isporb,is)) > 1.0e-4_dp ) write(stdout,*) istate,ispin,fs_occ(isporb,is)
   enddo
  
!   write(stdout,'(1x,a,/)') '========='

 enddo



 !
 ! Build the Lehman amplitude for virtual states
 !
 fs_virt(:,:) = 0.0_dp
 do jconf=1,conf_4e%nconf
   jjsporb(:) = conf_4e%sporb_occ(:,jconf)
   jjspin(:)  = sporb_to_spin(  jjsporb(:) )
   jjstate(:) = sporb_to_state( jjsporb(:) )

   on_j(:) = sporb_to_on(jjsporb)

   do iconf=1,conf_5e%nconf
     iisporb(:5) = conf_5e%sporb_occ(:,iconf)
     iispin(:5)  = sporb_to_spin(  iisporb(:5) )
     iistate(:5) = sporb_to_state( iisporb(:5) )

     on_i(:) = sporb_to_on(iisporb(:5))


     !
     ! Evaluate for any s, < 5 , s | a_i^+ | 4 , 0 >
     !
     do isporb=2*nfrozen_ci+1,2*nstate_ci
       ! Cannot create an electron in a state that is already occupied
       if( isporb == jjsporb(1) .OR. isporb == jjsporb(2) ) cycle
       if( isporb == jjsporb(3) .OR. isporb == jjsporb(4) ) cycle
       on_tmp(:) = on_j(:)
       on_tmp(isporb) = on_tmp(isporb) + 1
       if( ALL( on_i(:) - on_tmp(:) == 0 ) ) then
         fs_virt(isporb,:) = fs_virt(isporb,:) + eigvec_5e(iconf,:ns_virt) * eigvec_4e(jconf,1) * gamma_sign(on_j,isporb)
       endif
     enddo


   enddo
 enddo

 write(stdout,*) '====================='
 write(stdout,*) 'Virtual states'
 write(stdout,*)

 do is=1,ns_virt
   es_virt(is) = energy_5e(is) - energy_4e(1)
   write(stdout,'(/,1x,a,i4,1x,f12.6)') '=== Excitation (eV): ',is,es_virt(is) * Ha_eV
  
!   write(stdout,'(1x,a,i2,a)') '=== fs_virt(',is,')'
   do isporb=2*nfrozen_ci+1,2*nstate_ci
     ispin  = sporb_to_spin( isporb )
     istate = sporb_to_state( isporb )
!     if( ABS(fs_virt(isporb,is)) > 1.0e-4_dp ) write(stdout,*) istate,ispin,fs_virt(isporb,is)
   enddo
  
!   write(stdout,'(1x,a,/)') '========='

 enddo

 do isporb=2*nsemin-1,2*nsemax
   ispin = sporb_to_spin(isporb)
   istate = sporb_to_state(isporb)

   write(ctmp3,'(i3.3)') istate
   write(ctmp1,'(i1)') MODULO( isporb-1 , 2 ) + 1

   open(newunit=unit_gf,file='exact_greens_function_state'//ctmp3//'_spin'//ctmp1//'.dat',action='write')
   do iomega=-nomega,nomega
     omega = iomega / REAL(nomega,dp) * 100.0_dp / Ha_eV
     gi_w = 0.0_dp

     do is=1,ns_virt
       gi_w = gi_w + fs_virt(isporb,is)**2 / ( omega - es_virt(is) - ieta )
     enddo
     do is=1,ns_occ
       gi_w = gi_w + fs_occ(isporb,is)**2 / ( omega - es_occ(is) + ieta )
     enddo

     write(unit_gf,'(4(1x,es18.8))') omega * Ha_eV, gi_w / Ha_eV, ABS(AIMAG(gi_w)) / pi / Ha_eV
   enddo
   close(unit_gf)
 enddo



 deallocate(fs_virt,fs_occ)
 deallocate(es_virt,es_occ)

 call stop_clock(timing_ci_selfenergy)

end subroutine full_ci_4electrons_selfenergy


!==================================================================
subroutine full_ci_1electron_on(save_coefficients,spinstate,nuc_nuc)
 implicit none

 logical,intent(in)         :: save_coefficients
 integer,intent(in)         :: spinstate
 real(dp),intent(in)        :: nuc_nuc
!=====
 real(dp),allocatable       :: h_ci(:,:)
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,1x,a,/)') 'Full CI for 1 electron'


 call setup_configurations_ci(1,spinstate,conf_1e)


 call clean_allocate('CI hamiltonian',h_ci,conf_1e%nconf,conf_1e%nconf)

 call build_ci_hamiltonian(conf_1e,h_ci)

 allocate(energy_1e(conf_1e%nconf))
 allocate(eigvec_1e(conf_1e%nconf,conf_1e%nconf))

 call diagonalize(conf_1e%nconf,h_ci,energy_1e,eigvec_1e)

 write(stdout,'(/,1x,a,f19.10)')   '     Uncorrelated energy (Ha): ',h_ci(1,1)
 write(stdout,'(1x,a,f19.10,/)')   '      Correlation energy (Ha): ',energy_1e(1) - h_ci(1,1)
 write(stdout,'(1x,a,f19.10)')     '     Ground-state energy (Ha): ',energy_1e(1) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '1st excited-state energy (Ha): ',energy_1e(2) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '2nd excited-state energy (Ha): ',energy_1e(3) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '3rd excited-state energy (Ha): ',energy_1e(4) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '4th excited-state energy (Ha): ',energy_1e(5) + nuc_nuc

 if( .NOT. save_coefficients ) then
   deallocate(energy_1e)
   deallocate(eigvec_1e)
 endif

 call clean_deallocate('CI hamiltonian',h_ci)

 call stop_clock(timing_full_ci)


end subroutine full_ci_1electron_on


!==================================================================
subroutine full_ci_2electrons_on(save_coefficients,spinstate,nuc_nuc)
 implicit none

 logical,intent(in)         :: save_coefficients
 integer,intent(in)         :: spinstate
 real(dp),intent(in)        :: nuc_nuc
!=====
 real(dp),allocatable       :: h_ci(:,:)
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,1x,a,/)') 'Full CI for 2 electrons'


 call setup_configurations_ci(2,spinstate,conf_2e)


 call clean_allocate('CI hamiltonian',h_ci,conf_2e%nconf,conf_2e%nconf)
 call build_ci_hamiltonian(conf_2e,h_ci)

 allocate(energy_2e(conf_2e%nconf))
 allocate(eigvec_2e(conf_2e%nconf,conf_2e%nconf))

 call diagonalize(conf_2e%nconf,h_ci,energy_2e,eigvec_2e)

 write(stdout,'(/,1x,a,f19.10)')   '     Uncorrelated energy (Ha): ',h_ci(1,1)
 write(stdout,'(1x,a,f19.10,/)')   '      Correlation energy (Ha): ',energy_2e(1) - h_ci(1,1)
 write(stdout,'(1x,a,f19.10)')     '     Ground-state energy (Ha): ',energy_2e(1) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '1st excited-state energy (Ha): ',energy_2e(2) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '2nd excited-state energy (Ha): ',energy_2e(3) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '3rd excited-state energy (Ha): ',energy_2e(4) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '4th excited-state energy (Ha): ',energy_2e(5) + nuc_nuc

 if( .NOT. save_coefficients ) then
   deallocate(energy_2e)
   deallocate(eigvec_2e)
 endif

 call clean_deallocate('CI hamiltonian',h_ci)

 call stop_clock(timing_full_ci)


end subroutine full_ci_2electrons_on


!==================================================================
subroutine full_ci_3electrons_on(save_coefficients,spinstate,nuc_nuc)
 implicit none

 logical,intent(in)         :: save_coefficients
 integer,intent(in)         :: spinstate
 real(dp),intent(in)        :: nuc_nuc
!=====
 real(dp),allocatable       :: h_ci(:,:)
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,1x,a,/)') 'Full CI for 3 electrons'


 call setup_configurations_ci(3,spinstate,conf_3e)


 call clean_allocate('CI hamiltonian',h_ci,conf_3e%nconf,conf_3e%nconf)

 call build_ci_hamiltonian(conf_3e,h_ci)

 allocate(energy_3e(conf_3e%nconf))
 allocate(eigvec_3e(conf_3e%nconf,conf_3e%nconf))

 call diagonalize(conf_3e%nconf,h_ci,energy_3e,eigvec_3e)

 write(stdout,'(/,1x,a,f19.10)')   '     Uncorrelated energy (Ha): ',h_ci(1,1)
 write(stdout,'(1x,a,f19.10,/)')   '      Correlation energy (Ha): ',energy_3e(1) - h_ci(1,1)
 write(stdout,'(1x,a,f19.10)')     '     Ground-state energy (Ha): ',energy_3e(1) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '1st excited-state energy (Ha): ',energy_3e(2) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '2nd excited-state energy (Ha): ',energy_3e(3) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '3rd excited-state energy (Ha): ',energy_3e(4) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '4th excited-state energy (Ha): ',energy_3e(5) + nuc_nuc

 if( .NOT. save_coefficients ) then
   deallocate(energy_3e)
   deallocate(eigvec_3e)
 endif

 call clean_deallocate('CI hamiltonian',h_ci)

 call stop_clock(timing_full_ci)


end subroutine full_ci_3electrons_on


!==================================================================
subroutine full_ci_4electrons_on(save_coefficients,spinstate,nuc_nuc)
 implicit none

 logical,intent(in)         :: save_coefficients
 integer,intent(in)         :: spinstate
 real(dp),intent(in)        :: nuc_nuc
!=====
 real(dp),allocatable       :: h_ci(:,:)
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,1x,a,/)') 'Full CI for 4 electrons'


 call setup_configurations_ci(4,spinstate,conf_4e)


 call clean_allocate('CI hamiltonian',h_ci,conf_4e%nconf,conf_4e%nconf)

 call build_ci_hamiltonian(conf_4e,h_ci)

 allocate(energy_4e(conf_4e%nconf))
 allocate(eigvec_4e(conf_4e%nconf,conf_4e%nconf))

#ifdef HAVE_SCALAPACK
 call diagonalize(conf_4e%nconf,h_ci,energy_4e,eigvec_4e)
#else
 call diagonalize_scalapack(scalapack_block_min,conf_4e%nconf,h_ci,energy_4e)
 eigvec_4e(:,:) = h_ci(:,:)
#endif

 write(stdout,'(/,1x,a,f19.10)')   '     Uncorrelated energy (Ha): ',h_ci(1,1)
 write(stdout,'(1x,a,f19.10,/)')   '      Correlation energy (Ha): ',energy_4e(1) - h_ci(1,1)
 write(stdout,'(1x,a,f19.10)')     '     Ground-state energy (Ha): ',energy_4e(1) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '1st excited-state energy (Ha): ',energy_4e(2) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '2nd excited-state energy (Ha): ',energy_4e(3) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '3rd excited-state energy (Ha): ',energy_4e(4) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '4th excited-state energy (Ha): ',energy_4e(5) + nuc_nuc

 if( .NOT. save_coefficients ) then
   deallocate(energy_4e)
   deallocate(eigvec_4e)
 endif

 call clean_deallocate('CI hamiltonian',h_ci)

 call stop_clock(timing_full_ci)


end subroutine full_ci_4electrons_on


!==================================================================
subroutine full_ci_5electrons_on(save_coefficients,spinstate,nuc_nuc)
 implicit none

 logical,intent(in)         :: save_coefficients
 integer,intent(in)         :: spinstate
 real(dp),intent(in)        :: nuc_nuc
!=====
 real(dp),allocatable       :: h_ci(:,:)
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,1x,a,/)') 'Full CI for 5 electrons'


 call setup_configurations_ci(5,spinstate,conf_5e)


 call clean_allocate('CI hamiltonian',h_ci,conf_5e%nconf,conf_5e%nconf)

 call build_ci_hamiltonian(conf_5e,h_ci)

 allocate(energy_5e(conf_5e%nconf))
 allocate(eigvec_5e(conf_5e%nconf,conf_5e%nconf))

 call start_clock(timing_tmp0)
#ifdef HAVE_SCALAPACK
 call diagonalize(conf_5e%nconf,h_ci,energy_5e,eigvec_5e)
#else
 call diagonalize_scalapack(scalapack_block_min,conf_5e%nconf,h_ci,energy_5e)
 eigvec_5e(:,:) = h_ci(:,:)
#endif
 call stop_clock(timing_tmp0)

 write(stdout,'(/,1x,a,f19.10)')   '     Uncorrelated energy (Ha): ',h_ci(1,1)
 write(stdout,'(1x,a,f19.10,/)')   '      Correlation energy (Ha): ',energy_5e(1) - h_ci(1,1)
 write(stdout,'(1x,a,f19.10)')     '     Ground-state energy (Ha): ',energy_5e(1) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '1st excited-state energy (Ha): ',energy_5e(2) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '2nd excited-state energy (Ha): ',energy_5e(3) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '3rd excited-state energy (Ha): ',energy_5e(4) + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '4th excited-state energy (Ha): ',energy_5e(5) + nuc_nuc

 if( .NOT. save_coefficients ) then
   deallocate(energy_5e)
   deallocate(eigvec_5e)
 endif

 call clean_deallocate('CI hamiltonian',h_ci)

 call stop_clock(timing_full_ci)


end subroutine full_ci_5electrons_on


end module m_ci
!==================================================================
