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

 integer,parameter,private     :: key_int=8


 integer,private              :: nfrozen_ci
 integer,private              :: nstate_ci

 real(dp),private             :: h_frozen_frozen
 real(dp),allocatable,private :: h_frozen_active(:)

 real(dp),allocatable,private :: h_1body(:,:)

 type, private :: configurations
   integer             :: nelec
   integer             :: nelec_valence
   integer             :: nelec_valence_up
   integer             :: nelec_valence_down
   integer             :: nconf
   integer             :: nstate
   integer             :: sz
   integer(kind=key_int),allocatable :: keyud(:,:)  ! its bits encode the occupation: 0 or 1
                                                    ! first element is spin up
                                                    ! second element is spin down
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

 type, private ::  sparse_matrix
   real(dp)             :: nnz_total
   integer              :: nnz
   real(dp),allocatable :: val(:)
   integer,allocatable  :: row_ind(:)
   integer,allocatable  :: col_ptr(:)
 end type


contains


!==================================================================
! Sign introduced by creation/annihilation operators
! as defined in Hellgaker's book (chapter 1 box 1)
!==================================================================
pure function gamma_sign_keyud(keyud,istate,ispin)
 implicit none

 integer(key_int),intent(in) :: keyud(2)
 integer,intent(in) :: istate,ispin
 integer            :: gamma_sign_keyud
!=====
!=====

 if( ispin == 1 ) then
   gamma_sign_keyud = POPPAR( IAND(keyud(1),MASKR(istate-nfrozen_ci)) ) * 2 - 1
 else
   gamma_sign_keyud = ( MODULO( POPPAR(keyud(1)) + POPPAR( IAND(keyud(2),MASKR(istate-nfrozen_ci)) ) , 2) ) * 2 - 1
 endif

end function gamma_sign_keyud


!==================================================================
pure function get_keyud(sporbup,sporbdown) result(keyud)
 implicit none

 integer,intent(in)    :: sporbup(:),sporbdown(:)
 integer(kind=key_int) :: keyud(2)
!=====
 integer :: ielec
!=====

 keyud(:) = 0_key_int
 do ielec=1,SIZE(sporbup)
   keyud(1) = keyud(1) + SHIFTL(1_key_int,sporbup(ielec)-1-nfrozen_ci)
 enddo
 do ielec=1,SIZE(sporbdown)
   keyud(2) = keyud(2) + SHIFTL(1_key_int,sporbdown(ielec)-1-nfrozen_ci)
 enddo

end function get_keyud


!==================================================================
subroutine increment_sporb(sporb)
 implicit none

 integer,intent(inout) :: sporb(:)
!=====
 integer :: nelec,ielec,jelec
 logical :: room
!=====

 nelec = SIZE(sporb)

 sporb(nelec) = sporb(nelec) + 1

 if( sporb(nelec) > nstate_ci ) then
   room = .FALSE.
   do ielec=nelec-1,1,-1
     if( sporb(ielec) < nstate_ci - (nelec - ielec) ) then
       room = .TRUE.
       sporb(ielec) = sporb(ielec) + 1
       do jelec=ielec+1,nelec
         sporb(jelec) = sporb(jelec-1) + 1
       enddo
       exit
     endif
   enddo
   ! If there is no room to increment again, then return a conventional zero
   if( .NOT. room ) sporb(:) = 0
 endif

end subroutine increment_sporb


!==================================================================
pure function get_spinz_from_keyud(keyud) result(spinz)
 implicit none

 integer(kind=key_int),intent(in) :: keyud(2)
 integer :: spinz
!=====
!=====

 spinz = POPCNT(keyud(1)) - POPCNT(keyud(2))

end function get_spinz_from_keyud


!==================================================================
! From key, get the spins = 1 or 2
subroutine get_spins_from_keyud(keyud,ispin)
 implicit none

 integer(kind=key_int),intent(in) :: keyud(2)
 integer,intent(out)              :: ispin(:)
!=====
 integer :: ipos,ielec
!=====

 ispin(:) = 2
 do ielec=1,POPCNT(keyud(1))
   ispin(ielec) = 1
 enddo

end subroutine get_spins_from_keyud


!==================================================================
! From key, get the states
subroutine get_states_from_keyud(keyud,istate)
 implicit none

 integer(kind=key_int),intent(in) :: keyud(2)
 integer,intent(out)              :: istate(:)
!=====
 integer :: ipos,ielec
!=====

 ielec = 0
 do ipos=0,BIT_SIZE(keyud(1))-1
   if( BTEST(keyud(1),ipos) ) then
     ielec = ielec + 1
     istate(ielec) = ipos + 1 + nfrozen_ci
   endif
 enddo
 do ipos=0,BIT_SIZE(keyud(2))-1
   if( BTEST(keyud(2),ipos) ) then
     ielec = ielec + 1
     istate(ielec) = ipos + 1 + nfrozen_ci
   endif
 enddo

end subroutine get_states_from_keyud


!==================================================================
subroutine prepare_ci(basis,nstate_in,nfrozen_in,c_matrix)
 use m_hamiltonian_onebody
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate_in,nfrozen_in
 real(dp),intent(in)        :: c_matrix(:,:,:)
!=====
 real(dp) :: h_1e(basis%nbf,basis%nbf)
 real(dp) :: hnuc(basis%nbf,basis%nbf)
!=====

 nstate_ci  = nstate_in
 nfrozen_ci = nfrozen_in

 write(stdout,'(/,1x,a,i4,a,i4)') 'Prepare CI with active states ranging from ',nfrozen_ci+1,' to ',nstate_ci

 write(stdout,'(1x,a,i4)') '    Max active space size: ',BIT_SIZE(0_key_int)
 write(stdout,'(1x,a,i4)') 'Current active space size: ',nstate_ci -nfrozen_ci
 if( nstate_ci -nfrozen_ci > BIT_SIZE(0_key_int) ) call die('prepare_ci: current active space too large')

 call setup_kinetic(basis,h_1e)
 call setup_nucleus(basis,hnuc)
 if( nelement_ecp > 0 ) then
   call setup_nucleus_ecp(basis,hnuc)
 endif
 h_1e(:,:) = h_1e(:,:) + hnuc(:,:)

 ! Calculate the one-electron hamiltonian on the eigenstate basis
 call build_1e_hamiltonian(c_matrix,h_1e)


end subroutine prepare_ci


!==================================================================
subroutine destroy_ci()
 implicit none
!=====
!=====

 call clean_deallocate('Eigenstate-Fock operator',h_1body)

 if( ALLOCATED(conf_0%keyud) ) deallocate(conf_0%keyud)
 if( ALLOCATED(conf_p%keyud) ) deallocate(conf_p%keyud)
 if( ALLOCATED(conf_m%keyud) ) deallocate(conf_m%keyud)

 if( ALLOCATED(eigvec_0) ) call clean_deallocate('CI eigenvectors',eigvec_0)
 if( ALLOCATED(eigvec_p) ) call clean_deallocate('CI eigenvectors',eigvec_p)
 if( ALLOCATED(eigvec_m) ) call clean_deallocate('CI eigenvectors',eigvec_m)

 if( ALLOCATED(energy_0) ) deallocate(energy_0)
 if( ALLOCATED(energy_p) ) deallocate(energy_p)
 if( ALLOCATED(energy_m) ) deallocate(energy_m)

 if( ALLOCATED(h_frozen_active) ) deallocate(h_frozen_active)

end subroutine destroy_ci


!==================================================================
subroutine setup_configurations_ci(nelec,spinstate,ci_type_in,conf)
 implicit none

 integer,intent(in)               :: nelec
 integer,intent(in)               :: spinstate
 character(len=*),intent(in)      :: ci_type_in
 type(configurations),intent(out) :: conf
!=====
 integer :: iconf,jconf
 integer :: excitation_max
 integer :: kstate,lstate
 integer :: ielec,isporb,ispin,istate
 integer :: nref,iref
 integer(kind=key_int) :: keyud_test(2)
 integer(kind=key_int),allocatable :: keyud_ref(:,:)
 integer :: fu
 logical :: file_exists
 character(len=142) :: string
 character(len=2)   :: ctmp2
 integer :: ref_sporb(2*(nstate_ci-nfrozen_ci+1))
 integer,allocatable :: sporbup(:),sporbdown(:)
 integer :: order_up,order_down
!=====

 call start_clock(timing_ci_config)

 write(stdout,'(/,1x,a,a)') 'Setup CI space with excitations: ',ci_type_in

 conf%nelec              = nelec
 conf%sz                 = spinstate
 conf%nelec_valence      = conf%nelec - 2 * nfrozen_ci
 conf%nelec_valence_up   = ( conf%nelec_valence + spinstate ) / 2
 conf%nelec_valence_down = ( conf%nelec_valence - spinstate ) / 2

 if( MODULO(conf%nelec,2) /= MODULO(conf%sz,2) .AND. .NOT. conf%sz == -100 ) then
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
 case('CISDTQ5')
   excitation_max = 5
 case('CISDTQ56')
   excitation_max = 6
 case default
   call die('setup_configurations_ci: ci_type not understood')
 end select


 !
 ! Set the reference states
 write(ctmp2,'(i2.2)') nelec
 inquire(file='manual_ci_ref_'//ctmp2,exist=file_exists)
 if(file_exists) then
   open(newunit=fu,file='manual_ci_ref_'//ctmp2,status='old',action='read')
   read(fu,'(i8)') nref
   allocate(keyud_ref(2,nref))

   do iref=1,nref
     ref_sporb(:) = 0
     read(fu,'(a)')  string
     call string_to_integers(string,ref_sporb)
     if( SUM(ref_sporb(:)) /= conf%nelec_valence ) &
       call die('setup_configuration_ci: manual_ci_ref_xx file does not have the correct number of valence electrons')
     keyud_ref(:,iref) = 0_key_int
     do isporb=1,SIZE(ref_sporb)
       ispin = MODULO(isporb-1,2) + 1
       istate = (isporb+1)/2
       if( ref_sporb(isporb) == 1 ) &
         keyud_ref(ispin,iref) = keyud_ref(ispin,iref) + SHIFTL(1_key_int,istate-1)
     enddo
     if( get_spinz_from_keyud(keyud_ref(:,iref)) /= spinstate ) then
       call die('setup_configuration_ci: manual_ci_ref_xx file does not represent the correct spin state')
     endif
   enddo
   close(fu)

 else
   nref = 1
   allocate(keyud_ref(2,nref))
   keyud_ref(:,1) = 0_key_int
   do ielec=1,conf%nelec_valence
     ispin = 2 - MODULO(ielec,2)
     keyud_ref(ispin,1) = keyud_ref(ispin,1) + SHIFTL(1_key_int,(ielec-1)/2)
   enddo
 endif
 write(stdout,'(1x,a,i4)') 'Number of reference states: ',nref
 do iref=1,nref
   write(stdout,'(1x,a,*(1x,b64.64))') 'Reference state key: ',keyud_ref(:,iref)
 enddo


 allocate(sporbup(conf%nelec_valence_up))
 allocate(sporbdown(conf%nelec_valence_down))

 !
 ! First, determine the total number of CI configurations
 do ielec=1,conf%nelec_valence_up
   sporbup(ielec) = ielec + nfrozen_ci
 enddo

 iconf = 0
 ! Spin up loop
 do
   do ielec=1,conf%nelec_valence_down
     sporbdown(ielec) = ielec + nfrozen_ci
   enddo

   keyud_test = get_keyud(sporbup,sporbdown)
   order_up = keys_diff_order(keyud_test(1),keyud_ref(1,:))

   if( order_up <= 2*excitation_max ) then

     ! Spin down loop
     do
       keyud_test = get_keyud(sporbup,sporbdown)
       order_down = keys_diff_order(keyud_test(2),keyud_ref(2,:))

       ! Check if the excitation is not too high
       if( order_up + order_down <= 2*excitation_max ) then
         iconf = iconf + 1
       endif
       if( SIZE(sporbdown) > 0 ) then
         call increment_sporb(sporbdown)
         if( sporbdown(1) == 0 ) exit
       else
         exit
       endif
     enddo

   endif

   call increment_sporb(sporbup)
   if( sporbup(1) == 0 ) exit
 enddo
 conf%nconf = iconf

 allocate(conf%keyud(2,conf%nconf))

 !
 ! Second, generate and store all the keys
 do ielec=1,conf%nelec_valence_up
   sporbup(ielec) = ielec + nfrozen_ci
 enddo

 iconf = 0
 ! Spin up loop
 do
   do ielec=1,conf%nelec_valence_down
     sporbdown(ielec) = ielec + nfrozen_ci
   enddo

   keyud_test = get_keyud(sporbup,sporbdown)
   order_up = keys_diff_order(keyud_test(1),keyud_ref(1,:))

   if( order_up <= 2*excitation_max ) then

     ! Spin down loop
     do
       keyud_test = get_keyud(sporbup,sporbdown)
       order_down = keys_diff_order(keyud_test(2),keyud_ref(2,:))

       ! Check if the excitation is not too high
       if( order_up + order_down <= 2*excitation_max ) then
         iconf = iconf + 1
         conf%keyud(:,iconf) = keyud_test(:)
       endif

       if( SIZE(sporbdown) > 0 ) then
         call increment_sporb(sporbdown)
         if( sporbdown(1) == 0 ) exit
       else
         exit
       endif
     enddo

   endif

   call increment_sporb(sporbup)
   if( sporbup(1) == 0 ) exit
 enddo

 deallocate(sporbup,sporbdown)
 deallocate(keyud_ref)


 if( iconf /= conf%nconf ) then
   write(stdout,*) 'iconf=',iconf,'conf%nconf=',conf%nconf
   call die('setup_configurations_ci: bug')
 endif


 write(stdout,'(/,3(1x,a,i2,/),1x,a,i8)')       &
                        'Valence electrons: ',conf%nelec_valence,      &
                        '     Up electrons: ',conf%nelec_valence_up,   &
                        '   Down electrons: ',conf%nelec_valence_down, &
                        '   Configurations: ',conf%nconf

 !
 ! Calculate the frozen-frozen contribution to the Hamiltonian
 !
 h_frozen_frozen = 0.0_dp
 do kstate=1,nfrozen_ci
   h_frozen_frozen = h_frozen_frozen + 2.0_dp * h_1body(kstate,kstate)
   do lstate=1,nfrozen_ci
     h_frozen_frozen = h_frozen_frozen + 2.0_dp * eri_eigen(kstate,kstate,1,lstate,lstate,1) - eri_eigen(kstate,lstate,1,kstate,lstate,1)
   enddo
 enddo
 !
 ! Calculate the frozen-active contribution to the Hamiltonian
 !
 if( .NOT. ALLOCATED(h_frozen_active) ) then
   allocate(h_frozen_active(nfrozen_ci+1:nstate_ci))

   h_frozen_active(:) = 0.0_dp
   do lstate=nfrozen_ci+1,nstate_ci
     do kstate=1,nfrozen_ci
       h_frozen_active(lstate) = h_frozen_active(lstate) &
                     + 2.0_dp * eri_eigen(kstate,kstate,1,lstate,lstate,1) - 1.0_dp * eri_eigen(kstate,lstate,1,kstate,lstate,1)
     enddo
   enddo

 endif

 call stop_clock(timing_ci_config)

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
function keyud_diff_order(keyud1,keyud2) RESULT(order)
 implicit none

 integer(kind=key_int),intent(in) :: keyud1(2),keyud2(2)
 integer                          :: order
!=====
!=====

 order = POPCNT( IEOR(keyud1(1),keyud2(1)) ) + POPCNT( IEOR(keyud1(2),keyud2(2)) )

end function keyud_diff_order


!==================================================================
function key_diff_order(key1,key2) RESULT(order)
 implicit none

 integer(kind=key_int),intent(in) :: key1,key2
 integer                          :: order
!=====
!=====

 order = POPCNT( IEOR(key1,key2) )

end function key_diff_order


!==================================================================
function keysud_diff_order(keyud1,keysud2) RESULT(order)
 implicit none

 integer(kind=key_int),intent(in) :: keyud1(2),keysud2(:,:)
 integer                          :: order
!=====
 integer :: iref
!=====

 order = 1000
 do iref=1,SIZE(keysud2,DIM=2)
   order = MIN(order, keyud_diff_order(keyud1,keysud2(:,iref)) )
 enddo

end function keysud_diff_order


!==================================================================
function keys_diff_order(key1,keys2) RESULT(order)
 implicit none

 integer(kind=key_int),intent(in) :: key1,keys2(:)
 integer                          :: order
!=====
 integer :: iref
!=====

 order = 1000
 do iref=1,SIZE(keys2)
   order = MIN(order, key_diff_order(key1,keys2(iref)) )
 enddo

end function keys_diff_order


!==================================================================
function hamiltonian_ci(keyudi,keyudj) RESULT(h_ci_ij)
 implicit none

 integer(kind=key_int),intent(in) :: keyudi(2),keyudj(2)
 real(dp)           :: h_ci_ij
!=====
 integer :: nelec_valence
 integer :: ielec,jelec
 integer :: isporb,jsporb,istate,jstate,ispin,jspin
 integer :: ksporb,lsporb,kstate,lstate,kspin,lspin
 integer :: msporb,mstate,mspin
 integer :: sign_factor
 integer(kind=key_int) :: keyud_ieor(2),keyud_iand(2),keyudii(2),keyudjj(2)
 integer,allocatable :: iistate(:),jjstate(:)
 integer,allocatable :: iispin(:),jjspin(:)
 integer,allocatable :: iisporb(:),jjsporb(:)
!=====

 !
 ! Follow the second-quantization rules from Hellgaker book Chapter 1.
 !

 h_ci_ij = 0.0_dp

 select case(keyud_diff_order(keyudi,keyudj))
 case(0)
   !
   ! Exact same ON-vector

   nelec_valence = POPCNT(keyudi(1)) + POPCNT(keyudi(2))


   allocate(iistate(nelec_valence),jjstate(nelec_valence))
   allocate(iispin(nelec_valence),jjspin(nelec_valence))

   call get_spins_from_keyud(keyudi,iispin)
   call get_states_from_keyud(keyudi,iistate)
   call get_spins_from_keyud(keyudj,jjspin)
   call get_states_from_keyud(keyudj,jjstate)

   h_ci_ij = h_frozen_frozen
   do ielec=1,nelec_valence
     istate = iistate(ielec)
     h_ci_ij = h_ci_ij + h_frozen_active(istate)
   enddo

   !
   ! 1-body part
   !
   do ielec=1,nelec_valence
     h_ci_ij = h_ci_ij + h_1body(iistate(ielec),iistate(ielec))
   enddo
   !
   ! 2-body part
   !
   do jelec=1,nelec_valence
     jstate = jjstate(jelec)

     do ielec=1,nelec_valence
       istate = iistate(ielec)
       h_ci_ij = h_ci_ij  &
                  + 0.5_dp * eri_eigen(istate,istate,1,jstate,jstate,1)

       if( iispin(ielec) == jjspin(jelec) )  &
         h_ci_ij = h_ci_ij  &
                    - 0.5_dp * eri_eigen(istate,jstate,1,jstate,istate,1)

     enddo
   enddo

   deallocate(iistate,jjstate)
   deallocate(iispin,jjspin)


 case(2)
   !
   ! ON-vectors differ by one occupation number

   !
   ! Identify the spin-orbitals that differ between the two Slater determinants
   ! If only one differs, bra: i      |  ket: k
   keyud_ieor(:) = IEOR(keyudi(:),keyudj(:))
   keyud_iand(:) = IAND(keyudi(:),keyudj(:))
   keyudii(:)    = IAND(keyud_ieor(:),keyudi(:))
   keyudjj(:)    = IAND(keyud_ieor(:),keyudj(:))


   if( keyud_ieor(1) /= 0_key_int ) then
     ispin = 1
     kspin = 1
     istate = TRAILZ(keyudii(1)) + nfrozen_ci + 1
     kstate = TRAILZ(keyudjj(1)) + nfrozen_ci + 1
   else
     ispin = -1
     kspin = -1
     istate = TRAILZ(keyudii(2)) + nfrozen_ci + 1
     kstate = TRAILZ(keyudjj(2)) + nfrozen_ci + 1
   endif

   sign_factor = gamma_sign_keyud(keyudi,istate,ispin) * gamma_sign_keyud(keyudj,kstate,kspin)


   if( ispin == kspin ) then

     !
     ! 1-body part
     !
     h_ci_ij = h_1body(istate,kstate) * sign_factor

     !
     ! 2-body part
     !

     ! Frozen states contribution
     do mstate=1,nfrozen_ci
       h_ci_ij = h_ci_ij + 2.0_dp * eri_eigen(istate,kstate,1,mstate,mstate,1) * sign_factor   &
                                  - eri_eigen(istate,mstate,1,kstate,mstate,1) * sign_factor
     enddo


     do mstate=nfrozen_ci+1,nfrozen_ci+BIT_SIZE(keyud_iand(1))
       if( BTEST(keyud_iand(1),mstate-nfrozen_ci-1) ) then
         mspin  = 1
         h_ci_ij = h_ci_ij + eri_eigen(istate,kstate,1,mstate,mstate,1) * sign_factor

         if( ispin == mspin ) &
           h_ci_ij = h_ci_ij - eri_eigen(istate,mstate,1,mstate,kstate,1)  * sign_factor
       endif
     enddo

     do mstate=nfrozen_ci+1,nfrozen_ci+BIT_SIZE(keyud_iand(2))
       if( BTEST(keyud_iand(2),mstate-nfrozen_ci-1) ) then
         mspin  = -1
         h_ci_ij = h_ci_ij + eri_eigen(istate,kstate,1,mstate,mstate,1) * sign_factor

         if( ispin == mspin ) &
           h_ci_ij = h_ci_ij - eri_eigen(istate,mstate,1,mstate,kstate,1)  * sign_factor
       endif
     enddo

   endif


 case(4)
   !
   ! ON-vectors differ by two occupation numbers

   !
   ! Identify the spin-orbitals that differ between the two Slater determinants
   ! If only one differs, bra: i      |  ket: k
   ! If two differ,       bra: i < j  |  ket: k < l

   keyud_ieor(:) = IEOR(keyudi,keyudj)
   keyudii(:)    = IAND(keyud_ieor,keyudi)
   keyudjj(:)    = IAND(keyud_ieor,keyudj)

   select case(POPCNT(keyud_ieor(1)))
   case(2)
     ispin =  1
     jspin = -1
     kspin =  1
     lspin = -1
     istate = TRAILZ(keyudii(1)) + nfrozen_ci + 1
     kstate = TRAILZ(keyudjj(1)) + nfrozen_ci + 1
     jstate = TRAILZ(keyudii(2)) + nfrozen_ci + 1
     lstate = TRAILZ(keyudjj(2)) + nfrozen_ci + 1
   case(0)
     ispin = -1
     jspin = -1
     kspin = -1
     lspin = -1
     istate = TRAILZ(keyudii(2)) + nfrozen_ci + 1
     kstate = TRAILZ(keyudjj(2)) + nfrozen_ci + 1
     jstate = BIT_SIZE(0_key_int) - LEADZ(keyudii(2)) + nfrozen_ci
     lstate = BIT_SIZE(0_key_int) - LEADZ(keyudjj(2)) + nfrozen_ci
   case(4)
     ispin = 1
     jspin = 1
     kspin = 1
     lspin = 1
     istate = TRAILZ(keyudii(1)) + nfrozen_ci + 1
     kstate = TRAILZ(keyudjj(1)) + nfrozen_ci + 1
     jstate = BIT_SIZE(0_key_int) - LEADZ(keyudii(1)) + nfrozen_ci
     lstate = BIT_SIZE(0_key_int) - LEADZ(keyudjj(1)) + nfrozen_ci
   end select

   sign_factor = gamma_sign_keyud(keyudi,istate,ispin) * gamma_sign_keyud(keyudi,jstate,jspin)  &
               * gamma_sign_keyud(keyudj,kstate,kspin) * gamma_sign_keyud(keyudj,lstate,lspin)

   if( ispin == kspin .AND. jspin == lspin ) &
     h_ci_ij = eri_eigen(istate,kstate,1,jstate,lstate,1) * sign_factor

   if( ispin == lspin .AND. jspin == kspin ) &
     h_ci_ij = h_ci_ij - eri_eigen(istate,lstate,1,jstate,kstate,1) * sign_factor


 end select


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

 mconf = SIZE(h_ci,DIM=1)
 nconf = SIZE(h_ci,DIM=2)

 do jconf=1,nconf
   jconf_global = colindex_local_to_global(desc_ham,jconf)

   do iconf=1,mconf   !TODO use symmetry of H to reduce the calculation by 2
     iconf_global = rowindex_local_to_global(desc_ham,iconf)

     h_ci(iconf,jconf) = hamiltonian_ci(conf%keyud(:,iconf_global),conf%keyud(:,jconf_global))

   enddo
 enddo

 call stop_clock(timing_ham_ci)

end subroutine build_ci_hamiltonian


!==================================================================
subroutine build_ci_hamiltonian_sparse(conf,desc,h)
 implicit none

 type(configurations),intent(in)   :: conf
 integer,intent(in)                :: desc(NDEL)
 type(sparse_matrix),intent(inout) :: h
!=====
 integer :: ii,mvec,nnztmp
 integer :: iconf,jconf
 integer :: jconf_global
 integer :: cntxt,nprow,npcol,iprow,ipcol
 real(dp) :: h_ij
!=====


 write(stdout,'(1x,a)') 'Build CI hamiltonian with sparse storage'


 cntxt = desc(CTXT_)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 mvec = NUMROC(conf%nconf,desc(MB_),iprow,first_row,nprow)

 call start_clock(timing_zeroes_ci)
 !
 ! Find the maximum size of the sparse CI hamiltonian
 h%nnz = 0
 do jconf=1,mvec
   jconf_global = rowindex_local_to_global(desc,jconf)
   do iconf=1,jconf_global-1   ! Use symmetry of H
     if( keyud_diff_order(conf%keyud(:,iconf),conf%keyud(:,jconf_global)) <= 4 ) h%nnz = h%nnz + 1
   enddo
 enddo
 h%nnz_total = REAL( h%nnz , dp )
 call xsum_auxil(h%nnz_total)
 write(stdout,'(1x,a,i10)')  ' CI hamiltonian elements on this proc: ',h%nnz
 nnztmp = h%nnz
 call xmin_auxil(nnztmp)
 write(stdout,'(1x,a,i10)')  'Min CI hamiltonian elements on a proc: ',nnztmp
 nnztmp = h%nnz
 call xmax_auxil(nnztmp)
 write(stdout,'(1x,a,i10)')  'Max CI hamiltonian elements on a proc: ',nnztmp
 ! total number of terms in the triangular matrix is N * (N-1) / 2
 write(stdout,'(1x,a,f8.3)') 'CI hamiltonian sparsity (%): ',h%nnz_total / REAL(conf%nconf,dp) / REAL(conf%nconf-1,dp) * 200.0_dp

 call stop_clock(timing_zeroes_ci)


 call start_clock(timing_ham_ci)

 call clean_allocate('Sparce CI values',h%val,h%nnz)
 call clean_allocate('Sparce CI indexes',h%row_ind,h%nnz)
 allocate(h%col_ptr(conf%nconf+1))

 ii = 0
 h%col_ptr(1) = ii + 1
 do jconf=1,mvec
   jconf_global = rowindex_local_to_global(desc,jconf)

   do iconf=1,jconf_global-1  ! Use symmetry of H

     h_ij  = hamiltonian_ci(conf%keyud(:,iconf),conf%keyud(:,jconf_global))

     if( ABS(h_ij) < 1.e-8_dp ) cycle

     ii = ii + 1
     h%val(ii) = h_ij
     h%row_ind(ii) = iconf

   enddo
   h%col_ptr(jconf+1) = ii + 1

 enddo

 h%nnz_total = REAL( h%col_ptr(mvec+1) ,dp)
 call xsum_auxil(h%nnz_total)
 write(stdout,'(1x,a,f8.3)') 'CI hamiltonian sparsity (%): ',h%nnz_total / REAL(conf%nconf,dp) / REAL(conf%nconf-1,dp) * 200.0_dp


 call stop_clock(timing_ham_ci)

end subroutine build_ci_hamiltonian_sparse


!==================================================================
subroutine full_ci_nelectrons_selfenergy()
 implicit none

!=====
 type(selfenergy_grid) :: se
 integer               :: is
 integer               :: iconf,jconf,kconf
 integer               :: iconf_global,jconf_global,kconf_global
 integer               :: istate,ispin
 integer               :: ns_occ,ns_virt
 real(dp),allocatable  :: fs_occ(:,:,:),fs_virt(:,:,:)
 real(dp),allocatable  :: es_occ(:),es_virt(:)
 integer               :: iomega
 complex(dp)           :: gi_w
 character(len=3)      :: ctmp3
 character(len=1)      :: ctmp1
 integer               :: unit_gf
 real(dp)              :: eigvec0(conf_0%nconf)
 real(dp)              :: energy0_dummy(nsemax,nspin)
 integer(kind=key_int) :: keyudi(2),keyudj(2)
!=====

 call start_clock(timing_ci_selfenergy)

 ns_occ  = 0
 ns_virt = 0
 write(stdout,'(/,1x,a,i4)') 'Full CI self-energy for electron count: ',conf_0%nelec
 write(stdout,'(1x,a,i3,a,i3)')    'Previous CI calculation had spin state Sz(',conf_0%nelec,'): ',conf_0%sz

 if( ALLOCATED(conf_p%keyud) ) then
   ns_occ  = conf_p%nstate
   write(stdout,'(1x,a,i4)') 'Hole     part will evaluated with excitations: ',ns_occ
   write(stdout,'(1x,a,i3,a,sp,i4)') 'Previous CI calculation had spin state Sz(',conf_p%nelec,'): ',conf_p%sz
 endif
 if( ALLOCATED(conf_m%keyud) ) then
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


 !
 ! Build the Lehman amplitude for occupied states
 !
 if( ns_occ > 0 ) then
   allocate(fs_occ(nstate_ci,2,ns_occ))
   allocate(es_occ(ns_occ))

   write(stdout,*) '====================='
   write(stdout,*) 'Occupied states'


   fs_occ(:,:,:) = 0.0_dp
   do jconf=1,conf_0%nconf
     keyudj(:) = conf_0%keyud(:,jconf)

     do iconf=1,SIZE(eigvec_p,DIM=1)
       iconf_global = rowindex_local_to_global(desc_p,iconf)
       keyudi(:) = conf_p%keyud(:,iconf_global)

       if( keyud_diff_order(keyudi,keyudj) /= 1 ) cycle

       ! Find the only spin and orbital that differ
       if( POPCNT(IEOR(keyudi(1),keyudj(1))) /= 0) then
         ispin  = 1
       else
         ispin  = 2
       endif
       istate = TRAILZ(IEOR(keyudi(ispin),keyudj(ispin))) + nfrozen_ci + 1

       !
       ! Evaluate for any s, < N , 0 | a_i^+ | N-1 , s >
       !
       do kconf=1,SIZE(eigvec_p,DIM=2)
         kconf_global = colindex_local_to_global(desc_p,kconf)
         if( kconf_global > ns_occ ) cycle
         fs_occ(istate,ispin,kconf_global) = fs_occ(istate,ispin,kconf_global) &
                                 + eigvec_p(iconf,kconf) * eigvec0(jconf) * gamma_sign_keyud(keyudi,istate,ispin)
       enddo

     enddo
   enddo

   call xsum_world(fs_occ)

   do is=1,ns_occ
     es_occ(is) = energy_0(1) - energy_p(is)
   enddo
   write(stdout,'(1x,a,f12.6,4x,f12.6)') '-IP (eV) and Weight: ',es_occ(1) * Ha_eV,fs_occ(conf_0%nelec/2,2,1)**2
   write(stdout,'(1x,a,f12.6,/)')        'Lowest excit.  (eV): ',es_occ(ns_occ) * Ha_eV

 endif


 !
 ! Build the Lehman amplitude for virtual states
 !
 if( ns_virt > 0 ) then
   allocate(fs_virt(nstate_ci,2,ns_virt))
   allocate(es_virt(ns_virt))

   write(stdout,*) '====================='
   write(stdout,*) 'Virtual states'


   fs_virt(:,:,:) = 0.0_dp
   do jconf=1,conf_0%nconf
     keyudj(:) = conf_0%keyud(:,jconf)

     do iconf=1,SIZE(eigvec_m,DIM=1)
       iconf_global = rowindex_local_to_global(desc_m,iconf)
       keyudi(:) = conf_m%keyud(:,iconf_global)


       if( keyud_diff_order(keyudi,keyudj) /= 1 ) cycle

       ! Find the only spin and orbital that differ
       if( POPCNT(IEOR(keyudi(1),keyudj(1))) /= 0) then
         ispin  = 1
       else
         ispin  = 2
       endif
       istate = TRAILZ(IEOR(keyudi(ispin),keyudj(ispin))) + nfrozen_ci + 1

       !
       ! Evaluate for any s, < N+1 , s | a_i^+ | N , 0 >
       !
       do kconf=1,SIZE(eigvec_m,DIM=2)
         kconf_global = colindex_local_to_global(desc_m,kconf)
         if( kconf_global > ns_virt ) cycle
         fs_virt(istate,ispin,kconf_global) = fs_virt(istate,ispin,kconf_global)  &
                         + eigvec_m(iconf,kconf) * eigvec0(jconf) * gamma_sign_keyud(keyudj,istate,ispin)
       enddo


     enddo
   enddo

   call xsum_world(fs_virt)

   do is=1,ns_virt
     es_virt(is) = energy_m(is) - energy_0(1)
   enddo

   write(stdout,'(1x,a,f12.6,4x,f12.6,/)') '-EA (eV) and Weight: ',es_virt(1) * Ha_eV,fs_virt(conf_0%nelec/2+1,2,1)**2

 endif


 !
 ! Setup the frequency grid
 !
 do istate=nsemin,nsemax
   energy0_dummy(istate,:) = 0.0_dp
 enddo
 call init_selfenergy_grid(one_shot,energy0_dummy,se)


 !
 ! Output CI Green's function
 !
 do ispin=1,2
   do istate=nsemin,nsemax

     write(ctmp3,'(i3.3)') istate
     write(ctmp1,'(i1)') ispin

     open(newunit=unit_gf,file='exact_greens_function_state'//ctmp3//'_spin'//ctmp1//'.dat',action='write')
     do iomega=-se%nomega,se%nomega

       gi_w = 0.0_dp
       do is=1,ns_virt
         gi_w = gi_w + fs_virt(istate,ispin,is)**2 / ( se%omega(iomega) - es_virt(is) - ieta )
       enddo
       do is=1,ns_occ
         gi_w = gi_w + fs_occ(istate,ispin,is)**2 / ( se%omega(iomega) - es_occ(is) + ieta )
       enddo

       write(unit_gf,'(8(1x,es18.8))') se%omega(iomega) * Ha_eV, gi_w / Ha_eV, ABS(AIMAG(gi_w)) / pi / Ha_eV
     enddo
     close(unit_gf)
   enddo
 enddo

 !
 ! Output CI Green's function summed over spin
 !
 do istate=nsemin,nsemax

   write(ctmp3,'(i3.3)') istate

   open(newunit=unit_gf,file='exact_greens_function_state'//ctmp3//'.dat',action='write')
   do iomega=-se%nomega,se%nomega

     gi_w = 0.0_dp
     do is=1,ns_virt
       gi_w = gi_w + SUM(fs_virt(istate,:,is)**2) / ( se%omega(iomega) - es_virt(is) - ieta )
     enddo
     do is=1,ns_occ
       gi_w = gi_w + SUM(fs_occ(istate,:,is)**2) / ( se%omega(iomega) - es_occ(is) + ieta )
     enddo

     write(unit_gf,'(8(1x,es18.8))') se%omega(iomega) * Ha_eV, gi_w / Ha_eV, ABS(AIMAG(gi_w)) / pi / Ha_eV
   enddo
   close(unit_gf)

   open(newunit=unit_gf,file='exact_spectral_function_state'//ctmp3//'.dat',action='write')
   do is=1,ns_virt
     write(unit_gf,'(8(1x,es18.8))') es_virt(is) * Ha_eV, SUM(fs_virt(istate,:,is)**2)
   enddo
   do is=1,ns_occ
     write(unit_gf,'(8(1x,es18.8))') es_occ(is) * Ha_eV, SUM(fs_occ(istate,:,is)**2)
   enddo
   close(unit_gf)

 enddo



 if( ALLOCATED(fs_virt) ) deallocate(fs_virt,es_virt)
 if( ALLOCATED(fs_occ ) ) deallocate(fs_occ,es_occ)

 call destroy_selfenergy_grid(se)

 call stop_clock(timing_ci_selfenergy)

end subroutine full_ci_nelectrons_selfenergy


#ifdef HAVE_SCALAPACK
!==================================================================
subroutine full_ci_nelectrons(save_coefficients,nelectron,spinstate,nuc_nuc)
 implicit none

 integer,intent(in)         :: save_coefficients
 integer,intent(in)         :: nelectron,spinstate
 real(dp),intent(in)        :: nuc_nuc
!=====
 logical,parameter            :: incore=.TRUE.
 character(len=12)            :: filename_eigvec
 integer,parameter            :: mb_max=64
 integer                      :: mb_sd,mb
 integer                      :: desc_ham(NDEL)
 integer                      :: mham,nham
 integer                      :: desc_vec(NDEL)
 integer                      :: mvec,nvec
 integer                      :: info,read_status
 integer                      :: nstate_read
 real(dp)                     :: residual_norm
 real(dp)                     :: ehf
 real(dp),allocatable         :: h_ci(:,:)
 type(configurations),pointer :: conf
 real(dp),allocatable         :: energy(:)
 real(dp),allocatable         :: eigvec(:,:)

 type(configurations)         :: conf_sd
 real(dp),allocatable         :: eigvec_sd(:,:)
 integer                      :: desc_vec_sd(NDEL)
 integer                      :: mvec_sd,nvec_sd

 type(sparse_matrix)          :: h
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
   call die('full_ci_nelectrons: error')
 end select

 call setup_configurations_ci(nelectron,spinstate,ci_type,conf)
 if( save_coefficients == 0 ) then
   conf%nstate = MIN(ci_nstate,conf%nconf)
 else
   conf%nstate = MIN(ci_nstate_self,conf%nconf)
 endif

 allocate(energy(conf%nconf))
 ehf = hamiltonian_ci(conf%keyud(:,1),conf%keyud(:,1))

 if( conf%nstate == conf%nconf ) then
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
   ! use a calculated block_size = mb
   mb = MIN( mb_max , 2**FLOOR( LOG(REAL(conf%nconf/nprow_auxil,dp)) / LOG(2.0_dp) ) )
   mvec = NUMROC(conf%nconf,mb,iprow_auxil,first_row,nprow_auxil)
   nvec = NUMROC(conf%nstate,mb,ipcol_auxil,first_col,npcol_auxil)
   call DESCINIT(desc_vec,conf%nconf,conf%nstate,mb,mb,first_row,first_col,cntxt_auxil,MAX(1,mvec),info)
   call clean_allocate('CI eigenvectors',eigvec,mvec,nvec)
   eigvec(:,:) = 0.0_dp
 endif

 if( conf%nstate == conf%nconf ) then
   call start_clock(timing_ci_diago)
   write(stdout,'(1x,a,i8,a,i8)') 'Full diagonalization of CI hamiltonian',conf%nconf,' x ',conf%nconf
   call diagonalize_sca(conf%nconf,desc_ham,h_ci,energy,desc_ham,eigvec)
   call stop_clock(timing_ci_diago)
 else
   write(stdout,'(1x,a,i8,a,i8)') 'Partial diagonalization of CI hamiltonian',conf%nconf,' x ',conf%nconf
   if( incore ) then


     call read_eigvec_ci(filename_eigvec,conf,desc_vec,eigvec,energy,nstate_read,residual_norm,read_status)
     if( residual_norm < toldav .AND. nstate_read == conf%nstate .AND. read_status == 0 ) then
       write(stdout,'(1x,a)') 'All CI eigenvectors are obtained from a previous converged calculation. Skip diagonalization'
       write(stdout,'(1x,a,es12.4,/)') 'Residual norm: ',residual_norm
     else
       call build_ci_hamiltonian_sparse(conf,desc_vec,h)
       call start_clock(timing_ci_diago)
       call diagonalize_davidson_ci(toldav,filename_eigvec,conf,conf%nstate,energy,desc_vec,eigvec,h)
       call stop_clock(timing_ci_diago)
       call clean_deallocate('Sparce CI values',h%val)
       call clean_deallocate('Sparce CI indexes',h%row_ind)
       deallocate(h%col_ptr)
     endif

   else


     call read_eigvec_ci(filename_eigvec,conf,desc_vec,eigvec,energy,nstate_read,residual_norm,read_status)

     if( read_status /= 0 ) then
       call setup_configurations_ci(nelectron,spinstate,'CISD',conf_sd)
       conf_sd%nstate = conf%nstate
       mb_sd = MIN( mb_max , 2**FLOOR( LOG(REAL(conf_sd%nconf/nprow_auxil,dp)) / LOG(2.0_dp) ) )

       mvec_sd = NUMROC(conf_sd%nconf,mb_sd,iprow_auxil,first_row,nprow_auxil)
       nvec_sd = NUMROC(conf_sd%nstate,mb_sd,ipcol_auxil,first_col,npcol_auxil)
       call DESCINIT(desc_vec_sd,conf_sd%nconf,conf_sd%nstate,mb_sd,mb_sd,first_row,first_col,cntxt_auxil,MAX(1,mvec),info)
       call clean_allocate('CISD eigenvectors',eigvec_sd,mvec_sd,nvec_sd)
       eigvec_sd(:,:) = 0.0_dp
       call start_clock(timing_ci_diago)
       call diagonalize_davidson_ci(toldav,'',conf_sd,conf_sd%nstate,energy(1:conf_sd%nconf),desc_vec_sd,eigvec_sd)
       call stop_clock(timing_ci_diago)
       call translate_eigvec_ci(conf_sd,desc_vec_sd,eigvec_sd,conf,desc_vec,eigvec)
       call clean_deallocate('CISD eigenvectors',eigvec_sd)
     endif

     call start_clock(timing_ci_diago)
     call diagonalize_davidson_ci(toldav,filename_eigvec,conf,conf%nstate,energy,desc_vec,eigvec)
     call stop_clock(timing_ci_diago)

   endif
 endif


 call clean_deallocate('CI hamiltonian',h_ci)

 write(stdout,'(/,1x,a,f19.10)')   '        Uncorrelated energy (Ha): ',ehf + nuc_nuc
 write(stdout,'(1x,a,f19.10)')     '         Correlation energy (Ha): ',energy(1) - ehf
 write(stdout,'(1x,a,f19.10,/)')   '     CI ground-state energy (Ha): ',energy(1) + nuc_nuc
 if( conf%nstate > 1 ) &
   write(stdout,'(1x,a,f19.10,5x,f12.6)')     'CI 1st excited-state energy (Ha), diff (eV): ',energy(2) + nuc_nuc, (energy(2)-energy(1)) * Ha_eV
 if( conf%nstate > 2 ) &
   write(stdout,'(1x,a,f19.10,5x,f12.6)')     'CI 2nd excited-state energy (Ha), diff (eV): ',energy(3) + nuc_nuc, (energy(3)-energy(1)) * Ha_eV
 if( conf%nstate > 3 ) &
   write(stdout,'(1x,a,f19.10,5x,f12.6)')     'CI 3rd excited-state energy (Ha), diff (eV): ',energy(4) + nuc_nuc, (energy(4)-energy(1)) * Ha_eV
 if( conf%nstate > 4 ) &
   write(stdout,'(1x,a,f19.10,5x,f12.6)')     'CI 4th excited-state energy (Ha), diff (eV): ',energy(5) + nuc_nuc, (energy(5)-energy(1)) * Ha_eV



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


end subroutine full_ci_nelectrons



!=========================================================================
subroutine diagonalize_davidson_ci(tolerance,filename,conf,neig_calc,eigval,desc_vec,eigvec,h)
 implicit none

 real(dp),intent(in)             :: tolerance
 character(len=*),intent(in)     :: filename
 type(configurations),intent(in) :: conf
 integer,intent(in)              :: neig_calc
 integer,intent(in)              :: desc_vec(NDEL)
 real(dp),intent(out)            :: eigval(conf%nconf)
 real(dp),intent(out)            :: eigvec(:,:)
 type(sparse_matrix),intent(in),optional :: h
!=====
 integer              :: nstep
 integer              :: iconf
 integer              :: mm,mm_max
 integer              :: ieig,icycle
 integer              :: mvec,nvec
 real(dp),allocatable :: bb(:,:),atilde(:,:),ab(:,:),qq(:,:)
 real(dp),allocatable :: lambda(:),alphavec(:,:)
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
 integer              :: mb,nb
 real(dp)             :: value_min
 integer              :: index_min
 integer              :: index_trial(neig_calc)
 logical              :: still_available(conf%nconf)
!=====

#if defined(HAVE_SCALAPACK)
 write(stdout,'(/,1x,a,i5)') 'Davidson diago for eigenvector count: ',neig_calc
 if( PRESENT(h) ) then
   write(stdout,'(1x,a,f10.3)') 'Sparse matrix incore implementation with storage size (Mb): ', &
                                h%nnz * 8.0_dp / 1024.0_dp**2
 endif

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
   ham_diag(iconf) = hamiltonian_ci(conf%keyud(:,iconf),conf%keyud(:,iconf))
 enddo

 nstep = nstep_dav
 mm     = 0
 mm_max = neig_calc * nstep
 if( mm_max > conf%nconf ) then
   nstep = conf%nconf / neig_calc
   mm_max = neig_calc * nstep
   write(stdout,'(1x,a)')    'Diagonalization problem is too small'
   write(stdout,'(1x,a,i4)') 'Number of Davidson steps has been reduced to ',nstep
 endif

 !
 ! Find the smallest diagonal terms => good candidates for initial trial vectors
 still_available(:) = .TRUE.
 do ieig=1,neig_calc
   index_min = 0
   value_min = 10000.0_dp
   do iconf=1,conf%nconf
     if( .NOT. still_available(iconf) ) cycle
     if( ham_diag(iconf) < value_min ) then
       value_min = ham_diag(iconf)
       index_min = iconf
     endif
   enddo
   index_trial(ieig) = index_min
   still_available(index_min) = .FALSE.
 enddo

 mbb = NUMROC(conf%nconf,mb,iprow,first_row,nprow)
 nbb = NUMROC(mm_max    ,nb,ipcol,first_col,npcol)
 call DESCINIT(desc_bb,conf%nconf,mm_max,mb,nb,first_row,first_col,cntxt,MAX(1,mbb),info)

 call clean_allocate('Trial vectors',bb,mbb,nbb)
 call clean_allocate('Hamiltonian applications',ab,mbb,nbb)

 mqq = NUMROC(conf%nconf,mb,iprow,first_row,nprow)
 nqq = NUMROC(neig_calc ,nb,ipcol,first_col,npcol)
 call DESCINIT(desc_qq,conf%nconf,neig_calc,mb,nb,first_row,first_col,cntxt,MAX(1,mqq),info)
 call clean_allocate('Error vectors',qq,mqq,nqq)


 !
 ! Initialize with stupid coefficients
 do ieig=1,neig_calc
   do ilocal=1,mbb
     iglobal = rowindex_local_to_global(desc_bb,ilocal)
     bb(ilocal,ieig) = MIN( EXP( -REAL(iglobal,dp) ) , 0.1_dp )
     if( iglobal == index_trial(ieig) ) bb(ilocal,ieig) = 1.0_dp
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
 call orthogonalize_sca(desc_bb,1,neig_calc,bb)


 ! Obtain Ab matrix
 call get_ab()


 do icycle=1,nstep

   mm = icycle * neig_calc

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
   call PDGEMM('N','N',conf%nconf,neig_calc,mm,1.0_dp,bb,1,1,desc_bb,alphavec,1,1,desc_at,0.0_dp,qq,1,1,desc_qq)
   eigvec(:,1:neig_calc) = qq(:,1:neig_calc)
   eigval(1:neig_calc)   = lambda(1:neig_calc)

   ! qq = qq * Lambda
   do ieig=1,neig_calc
     call PDSCAL(conf%nconf,lambda(ieig),qq,1,ieig,desc_qq,1)
   enddo


   ! qq = ab * alphavec - lambda * bb * alphavec
   call PDGEMM('N','N',conf%nconf,neig_calc,mm,1.0_dp,ab,1,1,desc_bb,alphavec,1,1,desc_at,-1.0_dp,qq,1,1,desc_qq)


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

   ! Write down found CI eigenvectors at each cycle
   call write_eigvec_ci(filename,conf,desc_vec,eigvec,eigval,residual_norm)


   !
   ! Convergence reached... or not
   if( residual_norm < tolerance ) then
     write(stdout,'(1x,a,es12.4)') 'Desired accuracy has been reached: ',toldav
     exit
   endif
   if( icycle == nstep ) then
     write(stdout,'(1x,a,es12.4)') 'Maximum number of Davidson steps performed without reaching desired accuracy: ',toldav
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

   call PDGEMR2D(conf%nconf,neig_calc,qq,1,1,desc_qq,bb,1,mm+1,desc_bb,cntxt)

   call orthogonalize_sca(desc_bb,mm+1,mm+neig_calc,bb)


   ! Obtain Ab matrix
   call get_ab()

   deallocate(lambda)


 enddo ! icycle

 if( ALLOCATED(lambda) ) deallocate(lambda)
 call clean_deallocate('Trial vectors',bb)
 call clean_deallocate('Hamiltonian applications',ab)
 call clean_deallocate('Error vectors',qq)


contains

subroutine get_ab()
 implicit none

!=====
 integer              :: iconf_min,iconf_max
 integer              :: ii,jj,rdest,mb_local
 integer              :: iconf,iconf_global
 integer              :: kconf,kconf_global
 real(dp),allocatable :: ab_iblock(:,:),h_iblock(:,:)
 real(dp),allocatable :: ab_i(:)
 real(dp),allocatable :: bb_i(:)
!=====

 if( .NOT. PRESENT(h) ) then
   do iconf_min=1,conf%nconf,mb
     iconf_max = MIN( iconf_min + mb - 1 , conf%nconf)
     mb_local = iconf_max - iconf_min + 1

     rdest = INDXG2P(iconf_min,mb,0,first_row,nprow)

     allocate(h_iblock(iconf_min:iconf_max,mvec))
     allocate(ab_iblock(iconf_min:iconf_max,neig_calc))

     do kconf=1,mvec
       kconf_global = indxl2g_pure(kconf,mb,iprow,first_row,nprow)
       do iconf=iconf_min,iconf_max
         h_iblock(iconf,kconf) = hamiltonian_ci(conf%keyud(:,iconf),conf%keyud(:,kconf_global))
       enddo
     enddo

     !ab_iblock(:,:) = MATMUL( h_iblock(:,:) , bb(:,mm+1:mm+neig_calc) )
     call DGEMM('N','N',mb_local,neig_calc,mvec,     &
                1.0_dp,h_iblock,mb_local,       &
                       bb(1,mm+1),mvec,         &
                0.0_dp,ab_iblock,mb_local)

     call DGSUM2D(cntxt,'A',' ',mb_local,neig_calc,ab_iblock,1,rdest,0)

     if( iprow == rdest ) then
       iconf = rowindex_global_to_local(desc_bb,iconf_min)
       ab(iconf:iconf+mb_local-1,mm+1:mm+neig_calc) = ab_iblock(:,:)
     endif
     deallocate(h_iblock)
     deallocate(ab_iblock)
   enddo

 else

   allocate(ab_i(conf%nconf),bb_i(conf%nconf))
   do ieig=mm+1,mm+neig_calc

     bb_i(:) = 0.0_dp
     do iconf=1,mvec
       iconf_global = rowindex_local_to_global(desc_bb,iconf)
       bb_i(iconf_global) =  bb(iconf,ieig)
     enddo

     call DGSUM2D(cntxt,'A',' ',conf%nconf,1,bb_i,1,-1,-1)

     ab_i(:) = 0.0_dp
     ! Use symmetry of H here
     do kconf=1,mvec
       kconf_global = indxl2g_pure(kconf,mb,iprow,first_row,nprow)
       ab_i(kconf_global) = ham_diag(kconf_global) * bb(kconf,ieig)

       do jj=h%col_ptr(kconf),h%col_ptr(kconf+1)-1
         ab_i(kconf_global) = ab_i(kconf_global) + h%val(jj) * bb_i(h%row_ind(jj))
       enddo
     enddo

     do kconf=1,mvec
       do ii=h%col_ptr(kconf),h%col_ptr(kconf+1)-1
         ab_i(h%row_ind(ii)) = ab_i(h%row_ind(ii)) + h%val(ii) * bb(kconf,ieig)
       enddo
     enddo

     call DGSUM2D(cntxt,'A',' ',conf%nconf,1,ab_i,1,-1,-1)

     do iconf=1,mvec
       iconf_global = rowindex_local_to_global(desc_bb,iconf)
       ab(iconf,ieig) = ab_i(iconf_global)
     enddo

   enddo
   deallocate(ab_i,bb_i)

 endif

end subroutine get_ab

#else
 call die('diagonalize_davidson_ci: only works with a compilation using SCALAPACK')
#endif

end subroutine diagonalize_davidson_ci
#endif


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

 if( conf_in%nelec_valence /= conf_out%nelec_valence ) then
   call die('translate_eigvec_ci: configuration spaces have different number of active electrons')
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
     if( ALL( conf_in%keyud(:,iconf_in) == conf_out%keyud(:,iconf_out) ) ) then
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
subroutine write_eigvec_ci(filename,conf,desc_vec,eigvec,eigval,residual_norm)
 implicit none

 character(len=*),intent(in)     :: filename
 type(configurations),intent(in) :: conf
 integer,intent(in)              :: desc_vec(NDEL)
 real(dp),intent(in)             :: eigvec(:,:)
 real(dp),intent(in)             :: eigval(:)
 real(dp),intent(in)             :: residual_norm
!=====
 integer :: iconf,iconf_local
 integer :: neig,ieig
 integer :: cifile
 real(dp),allocatable :: eigvec_tmp(:)
!=====

 if( LEN(filename) == 0 ) return

 call start_clock(timing_ci_write)

 neig = SIZE(eigvec(:,:),DIM=2)
 allocate(eigvec_tmp(conf%nconf))

 if( neig /= conf%nstate ) call die('write_eigvec_ci: distribution grid should be row-only')

 write(stdout,'(1x,a,a,5x,a,f12.3)') 'Write CI eigvectors in ',TRIM(filename), &
                                     'File size (Mb): ',REAL(conf%nconf,dp) * REAL(neig,dp) * 8.0_dp / 1024.0_dp**2

 if( is_iomaster ) then
   open(newunit=cifile,file=TRIM(filename),form='unformatted',status='unknown',action='write')
   write(cifile) conf%nconf
   write(cifile) conf%nstate
   write(cifile) residual_norm
   write(cifile) eigval(1:conf%nstate)
 endif

 do ieig=1,neig
   eigvec_tmp(:) = 0.0_dp
   do iconf_local=1,SIZE(eigvec(:,:),DIM=1)
     iconf = rowindex_local_to_global(desc_vec,iconf_local)
     eigvec_tmp(iconf) = eigvec(iconf_local,ieig)
   enddo
   call xsum_auxil(iomaster,eigvec_tmp)
   if( is_iomaster ) write(cifile) eigvec_tmp(:)
 enddo
 if( is_iomaster ) close(cifile)

 deallocate(eigvec_tmp)

 call stop_clock(timing_ci_write)

end subroutine write_eigvec_ci


!==================================================================
subroutine read_eigvec_ci(filename,conf,desc_vec,eigvec,eigval,nstate_read,residual_norm,read_status)
 implicit none

 character(len=*),intent(in)     :: filename
 type(configurations),intent(in) :: conf
 integer,intent(in)              :: desc_vec(NDEL)
 integer,intent(out)             :: read_status
 real(dp),intent(out)            :: eigvec(:,:)
 real(dp),intent(out)            :: eigval(:)
 integer,intent(out)             :: nstate_read
 real(dp),intent(out)            :: residual_norm
!=====
 logical :: file_exists
 integer :: iconf,iconf_local
 integer :: neig,ieig
 integer :: cifile
 real(dp),allocatable :: eigvec_tmp(:)
 integer :: nconf_read
!=====

 neig = SIZE(eigvec(:,:),DIM=2)
 eigvec(:,:) = 0.0_dp
 eigval(:)   = 0.0_dp
 read_status = 1
 residual_norm = 1.0e8_dp
 nstate_read = 0

 if( neig /= conf%nstate ) call die('read_eigvec_ci: distribution grid should be row-only')

 inquire(file=TRIM(filename),exist=file_exists)
 if( .NOT. file_exists ) then
   write(stdout,'(1x,a,a,a)') 'File ',TRIM(filename),' does not exist'
   return
 endif

 write(stdout,'(/,1x,a,a)')   'Read the CI eigvectors from file: ',TRIM(filename)

 open(newunit=cifile,file=TRIM(filename),form='unformatted',status='unknown',action='read')

 read(cifile) nconf_read
 read(cifile) nstate_read
 read(cifile) residual_norm
 read(cifile) eigval(1:nstate_read)

 if( nconf_read /= conf%nconf ) then
   call issue_warning(TRIM(filename)//' file not compatible: ignore it')
   close(cifile)
   return
 endif

 if( nstate_read < conf%nstate ) then
   call issue_warning(TRIM(filename)//' file contains fewer states than required. Read and use them anyway')
   close(cifile)
   return
 endif


 allocate(eigvec_tmp(conf%nconf))

 nstate_read = MIN(conf%nstate,nstate_read)
 do ieig=1,nstate_read
   read(cifile) eigvec_tmp(:)
   do iconf_local=1,SIZE(eigvec(:,:),DIM=1)
     iconf = rowindex_local_to_global(desc_vec,iconf_local)

     eigvec(iconf_local,ieig) = eigvec_tmp(iconf)

   enddo
 enddo
 close(cifile)


 deallocate(eigvec_tmp)

 ! Set the value to 0 when the reading was carried out properly
 read_status = 0

end subroutine read_eigvec_ci


!==================================================================
end module m_ci
!==================================================================
