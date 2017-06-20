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
 use m_selfenergy_tools


 integer,private              :: nfrozen_ci
 integer,private              :: nstate_ci

 real(dp),allocatable,private :: h_1body(:,:)

 type, private :: configurations
   integer             :: nelec
   integer             :: nconf
   integer             :: nstate
   integer             :: sz
   integer,allocatable :: sporb_occ(:,:)        ! spinor-orbitals with occupation equal to 1
 end type

 type(configurations),target,private :: conf_0  ! Neutral    configuration: N electrons
 type(configurations),target,private :: conf_p  ! +1-charged configuration: N-1 electrons
 type(configurations),target,private :: conf_m  ! -1-charged configuration: N+1 electrons

 integer,private                     :: desc_0(NDEL)
 integer,private                     :: desc_p(NDEL)
 integer,private                     :: desc_m(NDEL)

 real(dp),allocatable,private        :: energy_0(:)
 real(dp),allocatable,private        :: energy_p(:)
 real(dp),allocatable,private        :: energy_m(:)

 real(dp),allocatable,private        :: eigvec_0(:,:)
 real(dp),allocatable,private        :: eigvec_p(:,:)
 real(dp),allocatable,private        :: eigvec_m(:,:)

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
subroutine destroy_ci()
 implicit none
!=====
!=====
 
 call clean_deallocate('Eigenstate-Fock operator',h_1body)

 if( ALLOCATED(conf_0%sporb_occ) ) deallocate(conf_0%sporb_occ)
 if( ALLOCATED(conf_p%sporb_occ) ) deallocate(conf_p%sporb_occ)
 if( ALLOCATED(conf_m%sporb_occ) ) deallocate(conf_m%sporb_occ)

 if( ALLOCATED(eigvec_0) ) call clean_deallocate('CI eigenvectors',eigvec_0)
 if( ALLOCATED(eigvec_p) ) call clean_deallocate('CI eigenvectors',eigvec_p)
 if( ALLOCATED(eigvec_m) ) call clean_deallocate('CI eigenvectors',eigvec_m)

 if( ALLOCATED(energy_0) ) deallocate(energy_0)
 if( ALLOCATED(energy_p) ) deallocate(energy_p)
 if( ALLOCATED(energy_m) ) deallocate(energy_m)


end subroutine destroy_ci


!==================================================================
subroutine setup_configurations_ci(nelec,spinstate,ci_type_in,conf)
 implicit none

 integer,intent(in)               :: nelec
 integer,intent(in)               :: spinstate
 character(len=*),intent(in)      :: ci_type_in
 type(configurations),intent(out) :: conf
!=====
 integer :: iconf
 integer :: ispin(nelec),istate(nelec)
 integer :: sporb(nelec)
 integer :: isporb,jsporb,ksporb,lsporb,msporb,nsporb,osporb,psporb
 integer :: excitation_max
 integer :: on_0(2*nstate_ci)
!=====

 conf%nelec = nelec
 conf%sz    = spinstate

 if( MODULO(conf%nelec,2) /= MODULO(conf%sz,2) ) then
   write(stdout,*) 'nelec=',conf%nelec
   write(stdout,*) 'spinz=',conf%sz
   call die('setup_configurations_ci: spin multiplicity is not compatible with the number of electrons')
 endif

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

 select case(ci_type_in)
 case('ALL')
   excitation_max = 100
 case('CISD')
   excitation_max = 2
 case('CISDT')
   excitation_max = 3
 case('CISDTQ')
   excitation_max = 4
 case default
   call die('setup_configurations_ci: ci_type not understood')
 end select

 ! Set the ground state vector
 on_0(:) = 0
 on_0(1:conf%nelec) = 1

 ! First, freeze the core electrons in the lowest available spin-orbitals
 do isporb=1,2*nfrozen_ci
   sporb(isporb) = isporb
 enddo

 ! Then populate the rest
 select case(conf%nelec-2*nfrozen_ci)
 case(1)
   conf%nconf = 0
   do isporb=2*nfrozen_ci+1,2*nstate_ci
     sporb(2*nfrozen_ci+1) = isporb
     ispin(:)  = sporb_to_spin(  sporb(:) )
     istate(:) = sporb_to_state( sporb(:) )
     if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
       conf%nconf = conf%nconf + 1
     endif
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do isporb=2*nfrozen_ci+1,2*nstate_ci
     sporb(2*nfrozen_ci+1) = isporb
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
       sporb(2*nfrozen_ci+1)  = isporb
       sporb(2*nfrozen_ci+2)  = jsporb
       ispin(:)  = sporb_to_spin(  sporb(:) )
       istate(:) = sporb_to_state( sporb(:) )
       if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
         if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
           conf%nconf = conf%nconf + 1
         endif
       endif
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do jsporb=2*nfrozen_ci+1,2*nstate_ci
     do isporb=jsporb+1,2*nstate_ci
       sporb(2*nfrozen_ci+1)  = isporb
       sporb(2*nfrozen_ci+2)  = jsporb
       ispin(:)  = sporb_to_spin( sporb(:) )
       istate(:) = sporb_to_state( sporb(:) )
       if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
         if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
           iconf = iconf + 1
           conf%sporb_occ(:,iconf) = sporb(:)
         endif
       endif
     enddo
   enddo


 case(3)
   conf%nconf = 0
   do ksporb=2*nfrozen_ci+1,2*nstate_ci
     do jsporb=ksporb+1,2*nstate_ci
       do isporb=jsporb+1,2*nstate_ci
         sporb(2*nfrozen_ci+1) = isporb
         sporb(2*nfrozen_ci+2) = jsporb
         sporb(2*nfrozen_ci+3) = ksporb
         ispin(:)  = sporb_to_spin(  sporb(:) )
         istate(:) = sporb_to_state( sporb(:) )
         if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
           if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
             conf%nconf = conf%nconf + 1
           endif
         endif
       enddo
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do ksporb=2*nfrozen_ci+1,2*nstate_ci
     do jsporb=ksporb+1,2*nstate_ci
       do isporb=jsporb+1,2*nstate_ci
         sporb(2*nfrozen_ci+1) = isporb
         sporb(2*nfrozen_ci+2) = jsporb
         sporb(2*nfrozen_ci+3) = ksporb
         ispin(:)  = sporb_to_spin( sporb(:) )
         istate(:) = sporb_to_state( sporb(:) )
         if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
           if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
             iconf = iconf + 1
             conf%sporb_occ(:,iconf) = sporb(:)
           endif
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
           sporb(2*nfrozen_ci+1) = isporb
           sporb(2*nfrozen_ci+2) = jsporb
           sporb(2*nfrozen_ci+3) = ksporb
           sporb(2*nfrozen_ci+4) = lsporb
           ispin(:)  = sporb_to_spin(  sporb(:) )
           istate(:) = sporb_to_state( sporb(:) )
           if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
             if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
               conf%nconf = conf%nconf + 1
             endif
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
           sporb(2*nfrozen_ci+1) = isporb
           sporb(2*nfrozen_ci+2) = jsporb
           sporb(2*nfrozen_ci+3) = ksporb
           sporb(2*nfrozen_ci+4) = lsporb
           ispin(:)  = sporb_to_spin( sporb(:) )
           istate(:) = sporb_to_state( sporb(:) )
           if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
             if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
               iconf = iconf + 1
               conf%sporb_occ(:,iconf) = sporb(:)
             endif
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
             sporb(2*nfrozen_ci+1) = isporb
             sporb(2*nfrozen_ci+2) = jsporb
             sporb(2*nfrozen_ci+3) = ksporb
             sporb(2*nfrozen_ci+4) = lsporb
             sporb(2*nfrozen_ci+5) = msporb
             ispin(:)  = sporb_to_spin(  sporb(:) )
             istate(:) = sporb_to_state( sporb(:) )
             if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
               if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                 conf%nconf = conf%nconf + 1
               endif
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
             sporb(2*nfrozen_ci+1) = isporb
             sporb(2*nfrozen_ci+2) = jsporb
             sporb(2*nfrozen_ci+3) = ksporb
             sporb(2*nfrozen_ci+4) = lsporb
             sporb(2*nfrozen_ci+5) = msporb
             ispin(:)  = sporb_to_spin( sporb(:) )
             istate(:) = sporb_to_state( sporb(:) )
             if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
               if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                 iconf = iconf + 1
                 conf%sporb_occ(:,iconf) = sporb(:)
               endif
             endif
           enddo
         enddo
       enddo
     enddo
   enddo

 case(6)
   conf%nconf = 0
   do nsporb=2*nfrozen_ci+1,2*nstate_ci
     do msporb=nsporb+1,2*nstate_ci
       do lsporb=msporb+1,2*nstate_ci
         do ksporb=lsporb+1,2*nstate_ci
           do jsporb=ksporb+1,2*nstate_ci
             do isporb=jsporb+1,2*nstate_ci
               sporb(2*nfrozen_ci+1) = isporb
               sporb(2*nfrozen_ci+2) = jsporb
               sporb(2*nfrozen_ci+3) = ksporb
               sporb(2*nfrozen_ci+4) = lsporb
               sporb(2*nfrozen_ci+5) = msporb
               sporb(2*nfrozen_ci+6) = nsporb
               ispin(:)  = sporb_to_spin(  sporb(:) )
               istate(:) = sporb_to_state( sporb(:) )
               if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
                 if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                   conf%nconf = conf%nconf + 1
                 endif
               endif
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do nsporb=2*nfrozen_ci+1,2*nstate_ci
     do msporb=nsporb+1,2*nstate_ci
       do lsporb=msporb+1,2*nstate_ci
         do ksporb=lsporb+1,2*nstate_ci
           do jsporb=ksporb+1,2*nstate_ci
             do isporb=jsporb+1,2*nstate_ci
               sporb(2*nfrozen_ci+1) = isporb
               sporb(2*nfrozen_ci+2) = jsporb
               sporb(2*nfrozen_ci+3) = ksporb
               sporb(2*nfrozen_ci+4) = lsporb
               sporb(2*nfrozen_ci+5) = msporb
               sporb(2*nfrozen_ci+6) = nsporb
               ispin(:)  = sporb_to_spin( sporb(:) )
               istate(:) = sporb_to_state( sporb(:) )
               if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
                 if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                   iconf = iconf + 1
                   conf%sporb_occ(:,iconf) = sporb(:)
                 endif
               endif
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo

 case(7)
   conf%nconf = 0
   do osporb=2*nfrozen_ci+1,2*nstate_ci
     do nsporb=osporb+1,2*nstate_ci
       do msporb=nsporb+1,2*nstate_ci
         do lsporb=msporb+1,2*nstate_ci
           do ksporb=lsporb+1,2*nstate_ci
             do jsporb=ksporb+1,2*nstate_ci
               do isporb=jsporb+1,2*nstate_ci
                 sporb(2*nfrozen_ci+1) = isporb
                 sporb(2*nfrozen_ci+2) = jsporb
                 sporb(2*nfrozen_ci+3) = ksporb
                 sporb(2*nfrozen_ci+4) = lsporb
                 sporb(2*nfrozen_ci+5) = msporb
                 sporb(2*nfrozen_ci+6) = nsporb
                 sporb(2*nfrozen_ci+7) = osporb
                 ispin(:)  = sporb_to_spin(  sporb(:) )
                 istate(:) = sporb_to_state( sporb(:) )
                 if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
                   if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                     conf%nconf = conf%nconf + 1
                   endif
                 endif
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do osporb=2*nfrozen_ci+1,2*nstate_ci
     do nsporb=osporb+1,2*nstate_ci
       do msporb=nsporb+1,2*nstate_ci
         do lsporb=msporb+1,2*nstate_ci
           do ksporb=lsporb+1,2*nstate_ci
             do jsporb=ksporb+1,2*nstate_ci
               do isporb=jsporb+1,2*nstate_ci
                 sporb(2*nfrozen_ci+1) = isporb
                 sporb(2*nfrozen_ci+2) = jsporb
                 sporb(2*nfrozen_ci+3) = ksporb
                 sporb(2*nfrozen_ci+4) = lsporb
                 sporb(2*nfrozen_ci+5) = msporb
                 sporb(2*nfrozen_ci+6) = nsporb
                 sporb(2*nfrozen_ci+7) = osporb
                 ispin(:)  = sporb_to_spin( sporb(:) )
                 istate(:) = sporb_to_state( sporb(:) )
                 if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
                   if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                     iconf = iconf + 1
                     conf%sporb_occ(:,iconf) = sporb(:)
                   endif
                 endif
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo

 case(8)
   conf%nconf = 0
   do psporb=2*nfrozen_ci+1,2*nstate_ci
     do osporb=psporb+1,2*nstate_ci
       do nsporb=osporb+1,2*nstate_ci
         do msporb=nsporb+1,2*nstate_ci
           do lsporb=msporb+1,2*nstate_ci
             do ksporb=lsporb+1,2*nstate_ci
               do jsporb=ksporb+1,2*nstate_ci
                 do isporb=jsporb+1,2*nstate_ci
                   sporb(2*nfrozen_ci+1) = isporb
                   sporb(2*nfrozen_ci+2) = jsporb
                   sporb(2*nfrozen_ci+3) = ksporb
                   sporb(2*nfrozen_ci+4) = lsporb
                   sporb(2*nfrozen_ci+5) = msporb
                   sporb(2*nfrozen_ci+6) = nsporb
                   sporb(2*nfrozen_ci+7) = osporb
                   sporb(2*nfrozen_ci+8) = psporb
                   ispin(:)  = sporb_to_spin(  sporb(:) )
                   istate(:) = sporb_to_state( sporb(:) )
                   if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
                     if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                       conf%nconf = conf%nconf + 1
                     endif
                   endif
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
   allocate(conf%sporb_occ(conf%nelec,conf%nconf))
   iconf = 0
   do psporb=2*nfrozen_ci+1,2*nstate_ci
     do osporb=psporb+1,2*nstate_ci
       do nsporb=osporb+1,2*nstate_ci
         do msporb=nsporb+1,2*nstate_ci
           do lsporb=msporb+1,2*nstate_ci
             do ksporb=lsporb+1,2*nstate_ci
               do jsporb=ksporb+1,2*nstate_ci
                 do isporb=jsporb+1,2*nstate_ci
                   sporb(2*nfrozen_ci+1) = isporb
                   sporb(2*nfrozen_ci+2) = jsporb
                   sporb(2*nfrozen_ci+3) = ksporb
                   sporb(2*nfrozen_ci+4) = lsporb
                   sporb(2*nfrozen_ci+5) = msporb
                   sporb(2*nfrozen_ci+6) = nsporb
                   sporb(2*nfrozen_ci+7) = osporb
                   sporb(2*nfrozen_ci+8) = psporb
                   ispin(:)  = sporb_to_spin( sporb(:) )
                   istate(:) = sporb_to_state( sporb(:) )
                   if( SUM(ispin(:)) == spinstate .OR. spinstate == -100 ) then
                     if( COUNT( sporb_to_on(sporb) - on_0(:) == 1 ) <= excitation_max ) then
                       iconf = iconf + 1
                       conf%sporb_occ(:,iconf) = sporb(:)
                     endif
                   endif
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo


 case default
   write(stdout,*) 'Active elecrons:',conf%nelec-2*nfrozen_ci
   call die('setup_configurations_ci: number of active electrons not coded as of today')
 end select

 write(stdout,'(/,1x,a,i2,a,i2,a,i8)') 'Electron count: ',conf%nelec, &
                                       '         Active electrons: ',conf%nelec-2*nfrozen_ci, &
                                       '           Configurations: ',conf%nconf


end subroutine setup_configurations_ci


!==================================================================
subroutine build_1e_hamiltonian(c_matrix,h_1e)
 implicit none

 real(dp),intent(in) :: c_matrix(:,:,:)
 real(dp),intent(in) :: h_1e(:,:)
!=====
 integer :: nstate
!=====

 if( SIZE(c_matrix,DIM=1) /= SIZE(h_1e,DIM=2) ) call die('build_1e_hamiltonian: c_matrix and h_1e have inconsistent sizes')
 if( SIZE(h_1e,DIM=1)     /= SIZE(h_1e,DIM=2) ) call die('build_1e_hamiltonian: h_1e is not square')
 nstate = SIZE(c_matrix,DIM=2)

 if( ALLOCATED(h_1body) ) return

 write(stdout,*) 'Obtain the one-electron Hamiltonian in the HF basis'

 call clean_allocate('Eigenstate-Fock operator',h_1body,nstate,nstate)
 h_1body(:,:) = MATMUL( TRANSPOSE(c_matrix(:,:,1)) , MATMUL( h_1e(:,:) , c_matrix(:,:,1) ) )


end subroutine build_1e_hamiltonian


!==================================================================
pure function hamiltonian_ci(iisporb,jjsporb) RESULT(h_ci_ij)
 implicit none

 integer,intent(in) :: iisporb(:)
 integer,intent(in) :: jjsporb(:)
 real(dp)           :: h_ci_ij
!=====
 integer :: nelec
 integer :: ielec,jelec
 integer,allocatable :: iistate(:),jjstate(:)
 integer,allocatable :: iispin(:),jjspin(:)
 integer :: on_i(2*nstate_ci),on_j(2*nstate_ci)
 logical :: mask(2*nstate_ci)
 integer :: isporb,jsporb,istate,jstate,ispin,jspin
 integer :: ksporb,lsporb,kstate,lstate,kspin,lspin
 integer :: excitation_order
!=====

 h_ci_ij = 0.0_dp

 nelec = SIZE(iisporb)
 allocate(iistate(nelec),jjstate(nelec))
 allocate(iispin(nelec),jjspin(nelec))

 iispin(:)  = sporb_to_spin(  iisporb(:) )
 iistate(:) = sporb_to_state( iisporb(:) )
 on_i(:)    = sporb_to_on( iisporb(:) )
 jjspin(:)  = sporb_to_spin(  jjsporb(:) )
 jjstate(:) = sporb_to_state( jjsporb(:) )
 on_j(:)    = sporb_to_on( jjsporb(:) )

 excitation_order = COUNT( on_j(:) - on_i(:) == 1 )

 select case(excitation_order)

 !
 ! Exact same ON-vector
 case(0)
   !
   ! 1-body part
   !
   do ielec=1,nelec
     h_ci_ij = h_ci_ij + h_1body(iistate(ielec),iistate(ielec))
   enddo
   !
   ! 2-body part
   !
   do jelec=1,nelec
     jstate = jjstate(jelec)

     do ielec=1,nelec
       istate = iistate(ielec)
       h_ci_ij = h_ci_ij  &
                  + 0.5_dp * eri_eigen(istate,istate,1,jstate,jstate,1)

       if( iispin(ielec) == jjspin(jelec) )  &
         h_ci_ij = h_ci_ij  &
                    - 0.5_dp * eri_eigen(istate,jstate,1,jstate,istate,1)

     enddo
   enddo


 !
 ! ON-vectors differ by one occupation number
 case(1)
   !
   ! 1-body part
   !
   jsporb = MAXLOC( on_j(:) - on_i(:) , DIM=1 )
   isporb = MINLOC( on_j(:) - on_i(:) , DIM=1 )
   istate = sporb_to_state( isporb )
   jstate = sporb_to_state( jsporb )

   if( sporb_to_spin( isporb ) == sporb_to_spin( jsporb ) ) &
     h_ci_ij = h_ci_ij  &
                    + h_1body(istate,jstate) * gamma_sign(on_j,jsporb) * gamma_sign(on_i,isporb)

   !
   ! 2-body part
   !
   isporb = MINLOC( on_j(:) - on_i(:) , DIM=1 )
   jsporb = MAXLOC( on_j(:) - on_i(:) , DIM=1 )
   istate = sporb_to_state( isporb )
   jstate = sporb_to_state( jsporb )
   ispin  = sporb_to_spin(  isporb )
   jspin  = sporb_to_spin(  jsporb )

   do ksporb=1,2*nstate_ci
     if( on_i(ksporb) == 0 ) cycle
     if( on_j(ksporb) == 0 ) cycle
     kstate = sporb_to_state( ksporb )
     kspin  = sporb_to_spin(  ksporb )

     if( ispin == jspin ) &
       h_ci_ij = h_ci_ij  &
                  + eri_eigen(istate,jstate,1,kstate,kstate,1)   &
                       * gamma_sign(on_i,isporb) * gamma_sign(on_j,jsporb)

     if( ispin == kspin .AND. jspin == kspin ) &
       h_ci_ij = h_ci_ij  &
                  - eri_eigen(istate,kstate,1,kstate,jstate,1)  &
                       * gamma_sign(on_i,isporb) * gamma_sign(on_j,jsporb)
   enddo


 !
 ! ON-vectors differ by two occupation numbers
 case(2)
   !
   ! 2-body part
   !
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
     h_ci_ij = h_ci_ij  &
                + eri_eigen(istate,kstate,1,jstate,lstate,1)           &
                     * gamma_sign(on_i,isporb) * gamma_sign(on_i,jsporb)  &
                     * gamma_sign(on_j,ksporb) * gamma_sign(on_j,lsporb)

   if( ispin == lspin .AND. jspin == kspin ) &
     h_ci_ij = h_ci_ij  &
                - eri_eigen(istate,lstate,1,jstate,kstate,1)           &
                     * gamma_sign(on_i,isporb) * gamma_sign(on_i,jsporb)  &
                     * gamma_sign(on_j,ksporb) * gamma_sign(on_j,lsporb)

 end select

 deallocate(iistate,jjstate)
 deallocate(iispin,jjspin)

end function hamiltonian_ci


!==================================================================
subroutine build_ci_hamiltonian(conf,desc_ham,h_ci)
 implicit none

 type(configurations),intent(in) :: conf
 integer,intent(in)              :: desc_ham(NDEL)
 real(dp),intent(out)            :: h_ci(:,:)
!=====
 integer :: mconf,nconf
 integer :: iconf,jconf
 integer :: iconf_global,jconf_global
!=====

 call start_clock(timing_ham_ci)

 write(stdout,'(1x,a)') 'Build CI hamiltonian'
 !
 ! Follow the second-quantization notations from Hellgaker book Chapter 1.
 ! Use occupation number vectors on_i(:) and on_j(:) filled with 0's and three 1's.
 !

 mconf = SIZE(h_ci,DIM=1)
 nconf = SIZE(h_ci,DIM=2)

 do jconf=1,nconf
   jconf_global = colindex_local_to_global(desc_ham,jconf)

   do iconf=1,mconf
     iconf_global = rowindex_local_to_global(desc_ham,iconf)
!     if( iconf_global > jconf_global ) cycle  !FBFB TODO: activate this

     h_ci(iconf,jconf) = hamiltonian_ci(conf%sporb_occ(:,iconf_global),conf%sporb_occ(:,jconf_global))

!     h_ci(jconf,iconf) = h_ci(iconf,jconf)   ! TODO: SCALAPACKization of this

   enddo

 enddo

 call stop_clock(timing_ham_ci)

end subroutine build_ci_hamiltonian


!==================================================================
subroutine full_ci_nelectrons_selfenergy()
 implicit none

!=====
 type(selfenergy_grid) :: se
 integer               :: is
 integer               :: iconf,jconf,kconf
 integer               :: iconf_global,jconf_global,kconf_global
 integer               :: on_i(2*nstate_ci),on_j(2*nstate_ci)
 integer               :: on_tmp(2*nstate_ci)
 integer,allocatable   :: iisporb(:),iistate(:),iispin(:)
 integer,allocatable   :: jjsporb(:),jjstate(:),jjspin(:)
 integer               :: isporb,istate,ispin
 integer               :: ns_occ,ns_virt
 real(dp),allocatable  :: fs_occ(:,:),fs_virt(:,:)
 real(dp),allocatable  :: es_occ(:),es_virt(:)
 integer               :: iomega
 complex(dp)           :: gi_w
 character(len=3)      :: ctmp3
 character(len=1)      :: ctmp1
 integer               :: unit_gf
 real(dp)              :: eigvec0(conf_0%nconf)
 real(dp)              :: energy0_dummy(nsemax,nspin)
!=====

 call start_clock(timing_ci_selfenergy)

 ns_occ  = 0
 ns_virt = 0
 write(stdout,'(/,1x,a,i4)') 'Full CI self-energy for electron count: ',conf_0%nelec
 write(stdout,'(1x,a,i3,a,i3)')    'Previous CI calculation had spin state Sz(',conf_0%nelec,'): ',conf_0%sz

 if( ALLOCATED(conf_p%sporb_occ) ) then
   ns_occ  = conf_p%nstate
   write(stdout,'(1x,a,i4)') 'Hole part will evaluated with excitations: ',ns_occ
   write(stdout,'(1x,a,i3,a,sp,i4)') 'Previous CI calculation had spin state Sz(',conf_p%nelec,'): ',conf_p%sz
 endif
 if( ALLOCATED(conf_m%sporb_occ) ) then
   ns_virt = conf_m%nstate
   write(stdout,'(1x,a,i4)') 'Electron part will evaluated with excitations: ',ns_virt
   write(stdout,'(1x,a,i3,a,sp,i4)') 'Previous CI calculation had spin state Sz(',conf_m%nelec,'): ',conf_m%sz
 endif



 !
 ! Gather the N-electron CI coefficients since only one is needed
 eigvec0(:) = 0.0_dp
 do jconf=1,SIZE(eigvec_0,DIM=2)
   jconf_global = colindex_local_to_global(desc_0,jconf)
   if( jconf_global /= 1 ) cycle
   do iconf=1,SIZE(eigvec_0,DIM=1)
     iconf_global = rowindex_local_to_global(desc_0,iconf)
     eigvec0(iconf_global) = eigvec_0(iconf,jconf)
   enddo
 enddo
 call xsum_world(eigvec0)

 allocate(jjsporb(conf_0%nelec))
 allocate(jjspin(conf_0%nelec))
 allocate(jjstate(conf_0%nelec))


 !
 ! Build the Lehman amplitude for occupied states
 !
 if( ns_occ > 0 ) then
   allocate(fs_occ(2*nstate_ci,ns_occ))
   allocate(es_occ(ns_occ))

   write(stdout,*) '====================='
   write(stdout,*) 'Occupied states'

   allocate(iisporb(conf_p%nelec))
   allocate(iispin(conf_p%nelec))
   allocate(iistate(conf_p%nelec))

   fs_occ(:,:) = 0.0_dp
   do jconf=1,conf_0%nconf
     jjsporb(:) = conf_0%sporb_occ(:,jconf)
     jjspin(:)  = sporb_to_spin(  jjsporb(:) )
     jjstate(:) = sporb_to_state( jjsporb(:) )

     on_j(:) = sporb_to_on(jjsporb)

     do iconf=1,SIZE(eigvec_p,DIM=1)
       iconf_global = rowindex_local_to_global(desc_p,iconf)
       iisporb(:) = conf_p%sporb_occ(:,iconf_global)
       iispin(:)  = sporb_to_spin(  iisporb(:) )
       iistate(:) = sporb_to_state( iisporb(:) )

       on_i(:) = sporb_to_on(iisporb(:))

       !
       ! Evaluate for any s, < N , 0 | a_i^+ | N-1 , s >
       !
       do isporb=1,2*nstate_ci
         on_tmp(:) = on_i(:)
         on_tmp(isporb) = on_tmp(isporb) + 1
         if( ALL( on_j(:) - on_tmp(:) == 0 ) ) then
           do kconf=1,SIZE(eigvec_p,DIM=2)
             kconf_global = colindex_local_to_global(desc_p,kconf)
             if( kconf_global > ns_occ ) cycle
             fs_occ(isporb,kconf_global) = fs_occ(isporb,kconf_global) &
                                     + eigvec_p(iconf,kconf) * eigvec0(jconf) * gamma_sign(on_i,isporb)
           enddo
         endif
       enddo

     enddo
   enddo

   call xsum_world(fs_occ)

   do is=1,ns_occ
     es_occ(is) = energy_0(1) - energy_p(is)
     !write(100,*) is,es_occ(is) * Ha_eV
   enddo
   write(stdout,'(1x,a,f12.6,4x,f12.6,/)')   '-IP (eV) and Weight: ',es_occ(1) * Ha_eV,fs_occ(2*nfrozen_ci+conf_0%nelec,1)**2

   deallocate(iisporb,iispin,iistate)

 endif


 !
 ! Build the Lehman amplitude for virtual states
 !
 if( ns_virt > 0 ) then
   allocate(fs_virt(2*nstate_ci,ns_virt))
   allocate(es_virt(ns_virt))

   write(stdout,*) '====================='
   write(stdout,*) 'Virtual states'

   allocate(iisporb(conf_m%nelec))
   allocate(iispin(conf_m%nelec))
   allocate(iistate(conf_m%nelec))

   fs_virt(:,:) = 0.0_dp
   do jconf=1,conf_0%nconf
     jjsporb(:) = conf_0%sporb_occ(:,jconf)
     jjspin(:)  = sporb_to_spin(  jjsporb(:) )
     jjstate(:) = sporb_to_state( jjsporb(:) )

     on_j(:) = sporb_to_on(jjsporb)

     do iconf=1,SIZE(eigvec_m,DIM=1)
       iconf_global = rowindex_local_to_global(desc_m,iconf)
       iisporb(:) = conf_m%sporb_occ(:,iconf_global)
       iispin(:)  = sporb_to_spin(  iisporb(:) )
       iistate(:) = sporb_to_state( iisporb(:) )

       on_i(:) = sporb_to_on(iisporb(:))


       !
       ! Evaluate for any s, < N+1 , s | a_i^+ | N , 0 >
       !
       do isporb=1,2*nstate_ci
         on_tmp(:) = on_j(:)
         on_tmp(isporb) = on_tmp(isporb) + 1
         if( ALL( on_i(:) - on_tmp(:) == 0 ) ) then
           do kconf=1,SIZE(eigvec_m,DIM=2)
             kconf_global = colindex_local_to_global(desc_m,kconf)
             if( kconf_global > ns_virt ) cycle
             fs_virt(isporb,kconf_global) = fs_virt(isporb,kconf_global)  &
                             + eigvec_m(iconf,kconf) * eigvec0(jconf) * gamma_sign(on_j,isporb)
           enddo
         endif
       enddo


     enddo
   enddo

   call xsum_world(fs_virt)

   do is=1,ns_virt
     es_virt(is) = energy_m(is) - energy_0(1)
   enddo

   write(stdout,'(1x,a,f12.6,4x,f12.6,/)') '-EA (eV) and Weight: ',es_virt(1) * Ha_eV,fs_virt(2*nfrozen_ci+conf_0%nelec+1,1)**2

   deallocate(iisporb,iispin,iistate)

 endif

 deallocate(jjsporb,jjspin,jjstate)

 !
 ! Setup the frequency grid
 ! 
 do istate=nsemin,nsemax
   energy0_dummy(istate,:) = 0.0_dp
 enddo
 call init_selfenergy_grid(one_shot,nsemax,energy0_dummy,se)


 !
 ! Output CI Green's function
 !
 do isporb=2*nsemin-1,2*nsemax
   ispin = sporb_to_spin(isporb)
   istate = sporb_to_state(isporb)

   write(ctmp3,'(i3.3)') istate
   write(ctmp1,'(i1)') MODULO( isporb-1 , 2 ) + 1

   open(newunit=unit_gf,file='exact_greens_function_state'//ctmp3//'_spin'//ctmp1//'.dat',action='write')
   do iomega=-se%nomega,se%nomega

     gi_w = 0.0_dp
     do is=1,ns_virt
       gi_w = gi_w + fs_virt(isporb,is)**2 / ( se%omega(iomega) - es_virt(is) - ieta )
     enddo
     do is=1,ns_occ
       gi_w = gi_w + fs_occ(isporb,is)**2 / ( se%omega(iomega) - es_occ(is) + ieta )
     enddo

     write(unit_gf,'(8(1x,es18.8))') se%omega(iomega) * Ha_eV, gi_w / Ha_eV, ABS(AIMAG(gi_w)) / pi / Ha_eV
   enddo
   close(unit_gf)
 enddo

 !
 ! Output CI Green's function summed over spin
 !
 do isporb=2*nsemin-1,2*nsemax,2
   istate = sporb_to_state(isporb)

   write(ctmp3,'(i3.3)') istate

   open(newunit=unit_gf,file='exact_greens_function_state'//ctmp3//'.dat',action='write')
   do iomega=-se%nomega,se%nomega

     gi_w = 0.0_dp
     do is=1,ns_virt
       gi_w = gi_w + fs_virt(isporb  ,is)**2 / ( se%omega(iomega) - es_virt(is) - ieta )
       gi_w = gi_w + fs_virt(isporb+1,is)**2 / ( se%omega(iomega) - es_virt(is) - ieta )
     enddo
     do is=1,ns_occ
       gi_w = gi_w + fs_occ(isporb  ,is)**2 / ( se%omega(iomega) - es_occ(is) + ieta )
       gi_w = gi_w + fs_occ(isporb+1,is)**2 / ( se%omega(iomega) - es_occ(is) + ieta )
     enddo

     write(unit_gf,'(8(1x,es18.8))') se%omega(iomega) * Ha_eV, gi_w / Ha_eV, ABS(AIMAG(gi_w)) / pi / Ha_eV
   enddo
   close(unit_gf)
 enddo



 if( ALLOCATED(fs_virt) ) deallocate(fs_virt,es_virt)
 if( ALLOCATED(fs_occ ) ) deallocate(fs_occ,es_occ)

 call destroy_selfenergy_grid(se)

 call stop_clock(timing_ci_selfenergy)

end subroutine full_ci_nelectrons_selfenergy


!==================================================================
subroutine full_ci_nelectrons_on(save_coefficients,nelectron,spinstate,nuc_nuc)
 implicit none

 integer,intent(in)         :: save_coefficients
 integer,intent(in)         :: nelectron,spinstate
 real(dp),intent(in)        :: nuc_nuc
!=====
 logical,parameter            :: incore=.FALSE.
 character(len=12)            :: filename_eigvec
 integer,parameter            :: mb=128
 integer                      :: desc_ham(NDEL)
 integer                      :: mham,nham
 integer                      :: desc_vec(NDEL)
 integer                      :: mvec,nvec
 integer                      :: info,read_status
 real(dp)                     :: ehf
 real(dp),allocatable         :: h_ci(:,:)
 type(configurations),pointer :: conf
 real(dp),allocatable         :: energy(:)
 real(dp),allocatable         :: eigvec(:,:)

 type(configurations)         :: conf_sd
 real(dp),allocatable         :: eigvec_sd(:,:)
 integer                      :: desc_vec_sd(NDEL)
 integer                      :: mvec_sd,nvec_sd
!=====

 call start_clock(timing_full_ci)

 write(stdout,'(/,1x,a,i4,/)') 'Full CI for electron count: ',nelectron

 select case(save_coefficients)
 case(0)
   conf => conf_0
   filename_eigvec = 'EIGVEC_CI_0' 
 case(-1)
   conf => conf_m
   filename_eigvec = 'EIGVEC_CI_M' 
 case(1)
   conf => conf_p
   filename_eigvec = 'EIGVEC_CI_P' 
 case default
   call die('full_ci_nelectrons_on: error')
 end select

 call setup_configurations_ci(nelectron,spinstate,ci_type,conf)
 if( save_coefficients == 0 ) then
   conf%nstate = 1                         ! MIN(ci_nstate,conf%nconf)
 else
   conf%nstate = MIN(ci_nstate,conf%nconf) ! conf%nconf
 endif

 allocate(energy(conf%nconf))
 ehf = hamiltonian_ci(conf%sporb_occ(:,1),conf%sporb_occ(:,1))

 if( incore .OR. conf%nstate == conf%nconf ) then
   mham = NUMROC(conf%nconf,block_row,iprow_sd,first_row,nprow_sd)
   nham = NUMROC(conf%nconf,block_col,ipcol_sd,first_col,npcol_sd)
   call DESCINIT(desc_ham,conf%nconf,conf%nconf,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mham),info)
  
   if( nprow_sd * npcol_sd > 1 ) then
     write(stdout,'(1x,a,i5,a,i5)') 'Use SCALAPACK proc grid: ',nprow_sd,' x ',npcol_sd
     write(stdout,'(1x,a,i5,a,i5)') '        Sub-matrix size: ',mham,' x ',nham
   endif
  
   call clean_allocate('CI hamiltonian',h_ci,mham,nham)
  
  
   call build_ci_hamiltonian(conf,desc_ham,h_ci)


   mvec = NUMROC(conf%nconf,block_row,iprow_sd,first_row,nprow_sd)
   nvec = NUMROC(conf%nstate,block_col,ipcol_sd,first_col,npcol_sd)
   call DESCINIT(desc_vec,conf%nconf,conf%nstate,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mvec),info)
   call clean_allocate('CI eigenvectors',eigvec,mvec,nvec)
 
 else
   !
   ! cntxt_auxil is a row-only distribution: ( Ncore x 1 )
   ! use block_size = mb
   mvec = NUMROC(conf%nconf,mb,iprow_auxil,first_row,nprow_auxil)
   nvec = NUMROC(conf%nstate,mb,ipcol_auxil,first_col,npcol_auxil)
   call DESCINIT(desc_vec,conf%nconf,conf%nstate,mb,mb,first_row,first_col,cntxt_auxil,MAX(1,mvec),info)
   call clean_allocate('CI eigenvectors',eigvec,mvec,nvec)
   eigvec(:,:) = 0.0_dp
 endif

 call start_clock(timing_ci_diago)
 if( conf%nstate == conf%nconf ) then
   write(stdout,'(1x,a,i8,a,i8)') 'Full diagonalization of CI hamiltonian',conf%nconf,' x ',conf%nconf
   call diagonalize_sca(conf%nconf,desc_ham,h_ci,energy,desc_ham,eigvec)
 else
   write(stdout,'(1x,a,i8,a,i8)') 'Partial diagonalization of CI hamiltonian',conf%nconf,' x ',conf%nconf
   if( incore ) then
     call diagonalize_davidson_sca(toldav,desc_ham,h_ci,conf%nstate,energy,desc_vec,eigvec)
   else

   
     call read_eigvec_ci(filename_eigvec,conf,desc_vec,eigvec,read_status)

     if( read_status /= 0 ) then
       call setup_configurations_ci(nelectron,spinstate,'CISD',conf_sd)
       conf_sd%nstate = conf%nstate
       mvec_sd = NUMROC(conf_sd%nconf,mb,iprow_auxil,first_row,nprow_auxil)
       nvec_sd = NUMROC(conf_sd%nstate,mb,ipcol_auxil,first_col,npcol_auxil)
       call DESCINIT(desc_vec_sd,conf_sd%nconf,conf_sd%nstate,mb,mb,first_row,first_col,cntxt_auxil,MAX(1,mvec),info)
       call clean_allocate('CISD eigenvectors',eigvec_sd,mvec_sd,nvec_sd)
       eigvec_sd(:,:) = 0.0_dp
       call diagonalize_davidson_ci(toldav,'',conf_sd,conf_sd%nstate,energy(1:conf_sd%nconf),desc_vec_sd,eigvec_sd)
       call translate_eigvec_ci(conf_sd,desc_vec_sd,eigvec_sd,conf,desc_vec,eigvec)
       call clean_deallocate('CISD eigenvectors',eigvec_sd)
     endif

     call diagonalize_davidson_ci(toldav,filename_eigvec,conf,conf%nstate,energy,desc_vec,eigvec)

   endif
 endif

 call stop_clock(timing_ci_diago)

 call clean_deallocate('CI hamiltonian',h_ci)

 write(stdout,'(/,1x,a,f19.10)')   '        Uncorrelated energy (Ha): ',ehf + nuc_nuc
 write(stdout,'(1x,a,f19.10,/)')   '         Correlation energy (Ha): ',energy(1) - ehf
 write(stdout,'(1x,a,f19.10)')     '     CI ground-state energy (Ha): ',energy(1) + nuc_nuc
 if( conf%nstate >= 5 ) then
   write(stdout,'(1x,a,f19.10)')     'CI 1st excited-state energy (Ha): ',energy(2) + nuc_nuc
   write(stdout,'(1x,a,f19.10)')     'CI 2nd excited-state energy (Ha): ',energy(3) + nuc_nuc
   write(stdout,'(1x,a,f19.10)')     'CI 3rd excited-state energy (Ha): ',energy(4) + nuc_nuc
   write(stdout,'(1x,a,f19.10)')     'CI 4th excited-state energy (Ha): ',energy(5) + nuc_nuc
 endif



 select case(save_coefficients)
 case(0)
   desc_0(:) = desc_vec
   call move_alloc(energy,energy_0)
   call move_alloc(eigvec,eigvec_0)
 case(-1)
   desc_m(:) = desc_vec
   call move_alloc(energy,energy_m)
   call move_alloc(eigvec,eigvec_m)
 case(1)
   desc_p(:) = desc_vec
   call move_alloc(energy,energy_p)
   call move_alloc(eigvec,eigvec_p)
 case default
   call clean_deallocate('CI Eigenvectors',eigvec)
   deallocate(energy)
 end select


 call stop_clock(timing_full_ci)


end subroutine full_ci_nelectrons_on


!=========================================================================
subroutine diagonalize_davidson_ci(tolerance,filename,conf,neig_calc,eigval,desc_vec,eigvec)
 implicit none

 real(dp),intent(in)             :: tolerance
 character(len=*),intent(in)     :: filename
 type(configurations),intent(in) :: conf
 integer,intent(in)              :: neig_calc
 integer,intent(in)              :: desc_vec(NDEL)
 real(dp),intent(out)            :: eigval(conf%nconf)
 real(dp),intent(out)            :: eigvec(:,:)
!=====
 integer,parameter    :: dim_factor=1
 integer              :: neig_dim
 integer              :: ncycle
 integer              :: iconf,jconf,kconf,iconf_global,kconf_global
 integer              :: mm,mm_max
 integer              :: ieig,icycle
 integer              :: mvec,nvec
 real(dp),allocatable :: bb(:,:),atilde(:,:),ab(:,:),qq(:,:)
 real(dp),allocatable :: lambda(:),alphavec(:,:)
 real(dp),allocatable :: ab_i(:),ab_iblock(:,:)
 real(dp),allocatable :: h_i(:),h_iblock(:,:)
 integer              :: iconf_min,iconf_max
 real(dp)             :: residual_norm,norm2_i
 integer              :: desc_bb(NDEL)
 integer              :: mbb,nbb
 integer              :: desc_at(NDEL)
 integer              :: mat,nat
 integer              :: desc_qq(NDEL)
 integer              :: mqq,nqq
 integer              :: info
 integer              :: ilocal,jlocal,iglobal,jglobal
 real(dp)             :: ham_diag(conf%nconf)
 real(dp)             :: rtmp
 integer              :: cntxt,nprow,npcol,iprow,ipcol
 integer              :: mb,nb,rdest,mb_local
!=====

 write(stdout,'(/,1x,a,i5)') 'Davidson diago for eigenvector count: ',neig_calc

 mvec = SIZE(eigvec(:,:),DIM=1)
 nvec = SIZE(eigvec(:,:),DIM=2)
 if( nvec /= neig_calc ) call die('diagonalize_davidson_ci: distribution grid not allowed')

 !
 ! Local scalapack information all obtained from the input desc_vec
 cntxt = desc_vec(CTXT_)
 mb    = desc_vec(MB_)
 nb    = desc_vec(NB_)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 write(stdout,'(1x,a,i4,a,i4)') 'Distribution:',nprow,' x ',npcol
 write(stdout,'(1x,a,i6,a,i8)') 'Local number of configurations: ',mvec,'  over a total of ',conf%nconf
 write(stdout,'(1x,a,i4)')      'Block size:',mb

 eigval(:) = 0.0_dp

 !
 ! All procs have the full diagonal
 do iconf=1,conf%nconf
   ham_diag(iconf) = hamiltonian_ci(conf%sporb_occ(:,iconf),conf%sporb_occ(:,iconf))
 enddo

 neig_dim = neig_calc * dim_factor
 ncycle = 30
 mm     = 0
 mm_max = neig_dim * ncycle
 if( mm_max > conf%nconf ) then
   ncycle = conf%nconf / neig_dim
   mm_max = neig_dim * ncycle
 endif

 mbb = NUMROC(conf%nconf,mb,iprow,first_row,nprow)
 nbb = NUMROC(mm_max    ,nb,ipcol,first_col,npcol)
 call DESCINIT(desc_bb,conf%nconf,mm_max,mb,nb,first_row,first_col,cntxt,MAX(1,mbb),info)

 call clean_allocate('Trial vectors',bb,mbb,nbb)
 call clean_allocate('Hamiltonian applications',ab,mbb,nbb)

 mqq = NUMROC(conf%nconf,mb,iprow,first_row,nprow)
 nqq = NUMROC(neig_dim  ,nb,ipcol,first_col,npcol)
 call DESCINIT(desc_qq,conf%nconf,neig_dim,mb,nb,first_row,first_col,cntxt,MAX(1,mqq),info)
 call clean_allocate('Error vectors',qq,mqq,nqq)


 !
 ! Initialize with stupid coefficients
 do jlocal=1,nbb
   jglobal = colindex_local_to_global(desc_bb,jlocal)
   do ilocal=1,mbb
     iglobal = rowindex_local_to_global(desc_bb,ilocal)
     bb(ilocal,jlocal) = MIN( EXP( -REAL(iglobal,dp) ) , 0.1_dp )
     if( iglobal == jglobal ) bb(ilocal,jlocal) = 1.0_dp
   enddo
 enddo
 !
 ! Then override them in case eigvec isn't empty
 rtmp = ABS(eigvec(1,1))
 call xsum_auxil(rtmp)
 if( rtmp > 1.0e-12_dp ) then 
   write(stdout,*) 'Found existing eigenvectors'
   bb(:,1:neig_calc) = eigvec(:,:)
 endif
 call orthogonalize_sca(desc_bb,1,neig_dim,bb)



 ! Block per block algo
 do iconf_min=1,conf%nconf,mb
   iconf_max = MIN( iconf_min + mb - 1 , conf%nconf)
   mb_local = iconf_max - iconf_min + 1

   rdest = INDXG2P(iconf_min,mb,0,first_row,nprow)

   allocate(h_iblock(iconf_min:iconf_max,mvec))
   allocate(ab_iblock(iconf_min:iconf_max,neig_dim))

   forall(iconf=iconf_min:iconf_max,kconf=1:mvec)
     h_iblock(iconf,kconf) = hamiltonian_ci(conf%sporb_occ(:,iconf), &
                                            conf%sporb_occ(:,indxl2g_pure(kconf,mb,iprow,first_row,nprow)))
   end forall

   !ab_iblock(:,:) = MATMUL( h_iblock(:,:) , bb(:,mm+1:mm+neig_dim) )
   call DGEMM('N','N',mb_local,neig_dim,mvec,     &
              1.0_dp,h_iblock,mb_local,       &
                     bb(1,mm+1),mvec,         &
              0.0_dp,ab_iblock,mb_local)

   call DGSUM2D(cntxt,'A',' ',mb_local,neig_dim,ab_iblock,1,rdest,0)
   if( iprow == rdest ) then
     iconf = rowindex_global_to_local(desc_bb,iconf_min)
     ab(iconf:iconf+mb_local-1,mm+1:mm+neig_dim) = ab_iblock(:,:)
   endif
   deallocate(h_iblock)
   deallocate(ab_iblock)
 enddo



 do icycle=1,ncycle

   mm = icycle * neig_dim

   allocate(lambda(mm))

   mat = NUMROC(mm,mb,iprow,first_row,nprow)
   nat = NUMROC(mm,nb,ipcol,first_col,npcol)
   call DESCINIT(desc_at,mm,mm,mb,nb,first_row,first_col,cntxt,MAX(1,mat),info)
   allocate(atilde(mat,nat),alphavec(mat,nat))

   !atilde(1:mm,1:mm) = MATMUL( TRANSPOSE(bb(:,1:mm)) , ab(:,1:mm) )
   call PDGEMM('T','N',mm,mm,conf%nconf,1.0_dp,bb,1,1,desc_bb,ab,1,1,desc_bb,0.0_dp,atilde,1,1,desc_at)

   call diagonalize_sca(mm,desc_at,atilde,lambda,desc_at,alphavec)

   deallocate(atilde)

   !write(stdout,*) 'icycle',icycle,lambda(1:mm)

   ! qq = bb * alphavec
   call PDGEMM('N','N',conf%nconf,neig_dim,mm,1.0_dp,bb,1,1,desc_bb,alphavec,1,1,desc_at,0.0_dp,qq,1,1,desc_qq)
   eigvec(:,1:neig_calc) = qq(:,1:neig_calc)
   eigval(1:neig_calc)   = lambda(1:neig_calc)

   ! Write down found CI eigenvectors at each cycle
   call write_eigvec_ci(filename,conf,desc_vec,eigvec)

   ! qq = qq * Lambda
   do ieig=1,neig_dim
     call PDSCAL(conf%nconf,lambda(ieig),qq,1,ieig,desc_qq,1)
   enddo


   ! qq = ab * alphavec - lambda * bb * alphavec
   call PDGEMM('N','N',conf%nconf,neig_dim,mm,1.0_dp,ab,1,1,desc_bb,alphavec,1,1,desc_at,-1.0_dp,qq,1,1,desc_qq)


   deallocate(alphavec)

   residual_norm = 0.0_dp
   do ieig=1,neig_calc
     norm2_i = 0.0_dp
     call PDNRM2(conf%nconf,norm2_i,qq,1,ieig,desc_qq,1)
     residual_norm = MAX( residual_norm , norm2_i )
   enddo
   call xmax_world(residual_norm)

   write(stdout,'(1x,a,i4,1x,i5,1x,es12.4,1x,f19.10)') 'Cycle, Subspace dim, Max residual norm, Electronic energy: ', &
                                                      icycle,mm,residual_norm,lambda(1)


   !
   ! Convergence reached... or not
   if( icycle == ncycle .OR. residual_norm < tolerance ) then
     exit
   endif


   !
   ! New trial vectors
   !
   do jlocal=1,nqq
     jglobal = colindex_local_to_global(desc_qq,jlocal)
     do ilocal=1,mqq
       iglobal = rowindex_local_to_global(desc_qq,ilocal)
       qq(ilocal,jlocal) = qq(ilocal,jlocal) / ( lambda(jglobal) - ham_diag(iglobal) )
     enddo
   enddo

   call PDGEMR2D(conf%nconf,neig_dim,qq,1,1,desc_qq,bb,1,mm+1,desc_bb,cntxt)

   call orthogonalize_sca(desc_bb,mm+1,mm+neig_dim,bb)



   !call start_clock(timing_tmp4)

   ! Block per block algo
   do iconf_min=1,conf%nconf,mb
     iconf_max = MIN( iconf_min + mb - 1 , conf%nconf)
     mb_local = iconf_max - iconf_min + 1

     rdest = INDXG2P(iconf_min,mb,0,first_row,nprow)

     allocate(h_iblock(iconf_min:iconf_max,mvec))
     allocate(ab_iblock(iconf_min:iconf_max,neig_dim))

     !call start_clock(timing_tmp5)
     forall(iconf=iconf_min:iconf_max,kconf=1:mvec)
       h_iblock(iconf,kconf) = hamiltonian_ci(conf%sporb_occ(:,iconf), &
                                              conf%sporb_occ(:,indxl2g_pure(kconf,mb,iprow,first_row,nprow)))
     end forall
     !call stop_clock(timing_tmp5)

     !call start_clock(timing_tmp6)
     !ab_iblock(:,:) = MATMUL( h_iblock(:,:) , bb(:,mm+1:mm+neig_dim) )
     call DGEMM('N','N',mb_local,neig_dim,mvec,     &
                1.0_dp,h_iblock,mb_local,       &
                       bb(1,mm+1),mvec,         &
                0.0_dp,ab_iblock,mb_local)
     !call stop_clock(timing_tmp6)

     !call start_clock(timing_tmp7)
     call DGSUM2D(cntxt,'A',' ',mb_local,neig_dim,ab_iblock,1,rdest,0)
     if( iprow == rdest ) then
       iconf = rowindex_global_to_local(desc_bb,iconf_min)
       ab(iconf:iconf+mb_local-1,mm+1:mm+neig_dim) = ab_iblock(:,:)
     endif
     !call stop_clock(timing_tmp7)
     deallocate(h_iblock)
     deallocate(ab_iblock)
   enddo

   !call stop_clock(timing_tmp4)


   deallocate(lambda)


 enddo ! icycle

 if( ALLOCATED(lambda) ) deallocate(lambda)
 call clean_deallocate('Trial vectors',bb)
 call clean_deallocate('Hamiltonian applications',ab)
 call clean_deallocate('Error vectors',qq)


end subroutine diagonalize_davidson_ci



!==================================================================
subroutine translate_eigvec_ci(conf_in,desc_vec_in,eigvec_in,conf_out,desc_vec_out,eigvec_out)
 implicit none

 type(configurations),intent(in) :: conf_in,conf_out
 integer,intent(in)              :: desc_vec_in(NDEL),desc_vec_out(NDEL)
 real(dp),intent(in)             :: eigvec_in(:,:)
 real(dp),intent(out)            :: eigvec_out(:,:)
!=====
 integer :: iconf_in,iconf_out,iconf_in_local,iconf_out_local
 integer :: neig
 logical :: found
 real(dp),allocatable :: eigvec_tmp(:)
!=====

 neig = SIZE(eigvec_in(:,:),DIM=2)
 allocate(eigvec_tmp(neig))

 write(stdout,'(/,1x,a)') 'Tranpose the CI eigvectors'
 write(stdout,'(1x,a,i8)') 'Initial space: ',conf_in%nconf
 write(stdout,'(1x,a,i8)') '  Final space: ',conf_out%nconf

 ! Preliminary checks
 if( conf_in%nelec /= conf_out%nelec ) then
   call die('translate_eigvec_ci: configuration spaces have different number of electrons')
 endif


 do iconf_in=1,conf_in%nconf
   found = .FALSE.
   iconf_in_local = rowindex_global_to_local(desc_vec_in,iconf_in)
   if( iconf_in_local /= 0 ) then
     eigvec_tmp(:) = eigvec_in(iconf_in_local,1:neig)
   else
     eigvec_tmp(:) = 0.0_dp
   endif
   call xsum_auxil(eigvec_tmp)

   do iconf_out=1,conf_out%nconf
     if( ALL( conf_in%sporb_occ(:,iconf_in) == conf_out%sporb_occ(:,iconf_out) ) ) then
       found = .TRUE.

       iconf_out_local = rowindex_global_to_local(desc_vec_out,iconf_out)
       if( iconf_out_local /= 0 ) &
         eigvec_out(iconf_out_local,1:neig) = eigvec_tmp(:)

     endif
   enddo

   if( .NOT. found) then
     call die('translate_eigvec_ci: final configuration spaces does not contain the initial configuration space')
   endif
 enddo

 deallocate(eigvec_tmp)

end subroutine translate_eigvec_ci


!==================================================================
subroutine write_eigvec_ci(filename,conf,desc_vec,eigvec)
 implicit none

 character(len=*),intent(in)     :: filename
 type(configurations),intent(in) :: conf
 integer,intent(in)              :: desc_vec(NDEL)
 real(dp),intent(in)             :: eigvec(:,:)
!=====
 integer :: iconf,iconf_local
 integer :: neig
 integer :: cifile
 real(dp),allocatable :: eigvec_tmp(:)
!=====

 if( LEN(filename) == 0 ) return

 call start_clock(timing_ci_write)

 neig = SIZE(eigvec(:,:),DIM=2)
 allocate(eigvec_tmp(neig))

 if( neig /= conf%nstate ) call die('write_eigvec_ci: distribution grid should be row-only')

 write(stdout,'(1x,a,a,5x,a,f12.3)') 'Write CI eigvectors in ',TRIM(filename), &
                                     'File size [Mb]: ',REAL(conf%nconf,dp) * REAL(neig,dp) * 8.0_dp / 1024.0_dp**2

 if( is_iomaster ) then
   open(newunit=cifile,file=TRIM(filename),form='unformatted',status='unknown',action='write')
   write(cifile) conf%nconf
   write(cifile) conf%nstate
 endif

 do iconf=1,conf%nconf
   iconf_local = rowindex_global_to_local(desc_vec,iconf)
   if( iconf_local /= 0 ) then
     eigvec_tmp(:) = eigvec(iconf_local,1:conf%nstate)
   else
     eigvec_tmp(:) = 0.0_dp
   endif
   call xsum_auxil(iomaster,eigvec_tmp)
   if( is_iomaster ) write(cifile) eigvec_tmp(:)
 enddo
 if( is_iomaster ) close(cifile)

 deallocate(eigvec_tmp)

 call stop_clock(timing_ci_write)

end subroutine write_eigvec_ci


!==================================================================
subroutine read_eigvec_ci(filename,conf,desc_vec,eigvec,read_status)
 implicit none

 character(len=*),intent(in)     :: filename
 type(configurations),intent(in) :: conf
 integer,intent(in)              :: desc_vec(NDEL)
 integer,intent(out)             :: read_status
 real(dp),intent(out)            :: eigvec(:,:)
!=====
 logical :: file_exists
 integer :: iconf,iconf_local
 integer :: neig
 integer :: cifile
 real(dp),allocatable :: eigvec_tmp(:)
 integer :: nconf_read,nstate_read
!=====

 neig = SIZE(eigvec(:,:),DIM=2)
 eigvec(:,:) = 0.0_dp
 read_status = 1

 if( neig /= conf%nstate ) call die('read_eigvec_ci: distribution grid should be row-only')

 inquire(file=TRIM(filename),exist=file_exists)
 if( .NOT. file_exists ) then
   write(stdout,'(a,a,a)') 'File ',TRIM(filename),' not found'
   return
 endif

 write(stdout,'(/,1x,a,a)')   'Read the CI eigvectors from file: ',TRIM(filename)

 open(newunit=cifile,file=TRIM(filename),form='unformatted',status='unknown',action='read')

 read(cifile) nconf_read
 read(cifile) nstate_read

 if( nconf_read /= conf%nconf .OR. nstate_read /= conf%nstate ) then
   call issue_warning(TRIM(filename)//' file not compatible: ignore it')
   return
 endif


 allocate(eigvec_tmp(conf%nstate))

 do iconf=1,conf%nconf
   iconf_local = rowindex_global_to_local(desc_vec,iconf)
   read(cifile) eigvec_tmp(:)

   if( iconf_local /= 0 ) then
     eigvec(iconf_local,1:conf%nstate) = eigvec_tmp(:)
   endif

 enddo
 close(cifile)

 deallocate(eigvec_tmp)

 ! Set the value to 0 when the reading was carried out properly
 read_status = 0

end subroutine read_eigvec_ci


!==================================================================
end module m_ci
!==================================================================
