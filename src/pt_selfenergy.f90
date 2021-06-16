!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the perturbation theory to 2nd order evaluation of the self-energy
!
!=========================================================================
subroutine pt2_selfenergy(selfenergy_approx,nstate,basis,occupation,energy,c_matrix,se,emp2)
  use m_definitions
  use m_mpi
  use m_warning
  use m_timing
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_selfenergy_tools
  implicit none

  integer,intent(in)         :: selfenergy_approx,nstate
  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
  type(selfenergy_grid),intent(inout) :: se
  real(dp),intent(out)       :: emp2
  !=====
  integer                 :: pstate,qstate
  complex(dp),allocatable :: selfenergy_ring(:,:,:)
  complex(dp),allocatable :: selfenergy_sox(:,:,:)
  integer                 :: iomega
  integer                 :: istate,jstate,kstate
  integer                 :: pqispin,jkspin
  real(dp)                :: fact_occ1,fact_occ2
  real(dp)                :: fi,fj,fk,ei,ej,ek
  complex(dp)             :: omega
  complex(dp)             :: fact_comp
  real(dp)                :: fact_energy
  real(dp)                :: emp2_sox,emp2_ring
  real(dp),allocatable    :: eri_eigenstate_i(:,:,:,:)
  real(dp)                :: coul_iqjk,coul_ijkq,coul_ipkj
  !=====

  call start_clock(timing_pt_self)

  emp2_ring = 0.0_dp
  emp2_sox  = 0.0_dp


  write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
  write(stdout,*) 'with the perturbative approach'



  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
  else
    allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
  endif



  allocate(selfenergy_ring(-se%nomega:se%nomega,nsemin:nsemax,nspin))
  allocate(selfenergy_sox (-se%nomega:se%nomega,nsemin:nsemax,nspin))


  selfenergy_ring(:,:,:) = 0.0_dp
  selfenergy_sox(:,:,:)  = 0.0_dp

  do pqispin=1,nspin
    do istate=ncore_G+1,nvirtual_G-1 !LOOP of the first Green's function
      if( MODULO( istate - (ncore_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      if( .NOT. has_auxil_basis ) then
        call calculate_eri_4center_eigen(c_matrix,istate,pqispin,eri_eigenstate_i)
      endif

      fi = occupation(istate,pqispin)
      ei = energy(istate,pqispin)

      !$OMP PARALLEL
      !$OMP DO PRIVATE(qstate,fj,ej,fk,ek,fact_occ1,fact_occ2,coul_ipkj,coul_iqjk,coul_ijkq,omega,fact_comp,fact_energy) &
      !$OMP REDUCTION(+:emp2_ring,emp2_sox)
      do pstate=nsemin,nsemax ! external loop ( bra )
        qstate=pstate         ! external loop ( ket )

        do jkspin=1,nspin
          do jstate=ncore_G+1,nvirtual_G-1  !LOOP of the second Green's function
            fj = occupation(jstate,jkspin)
            ej = energy(jstate,jkspin)

            do kstate=ncore_G+1,nvirtual_G-1 !LOOP of the third Green's function
              fk = occupation(kstate,jkspin)
              ek = energy(kstate,jkspin)

              fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**3
              fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**3

              if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle

              if( has_auxil_basis ) then
                coul_ipkj = eri_eigen_ri(istate,pstate,pqispin,kstate,jstate,jkspin)
                coul_iqjk = eri_eigen_ri(istate,qstate,pqispin,jstate,kstate,jkspin)
                if( pqispin == jkspin ) then
                  coul_ijkq = eri_eigen_ri(istate,jstate,pqispin,kstate,qstate,pqispin)
                endif
              else
                coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
                coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
                if( pqispin == jkspin ) then
                  coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin)
                endif
              endif

              do iomega=-se%nomega,se%nomega
                omega = energy(qstate,pqispin) + se%omega(iomega)

                fact_comp   = fact_occ1 / ( omega - ei + ej - ek + ieta ) &
                            + fact_occ2 / ( omega - ei + ej - ek - ieta )
                fact_energy = REAL( fact_occ1 / (energy(pstate,pqispin) - ei + ej - ek + ieta) , dp )

                selfenergy_ring(iomega,pstate,pqispin) = selfenergy_ring(iomega,pstate,pqispin) &
                         + fact_comp * coul_ipkj * coul_iqjk * spin_fact

                if(iomega==0 .AND. occupation(pstate,pqispin)>completely_empty) then
                  emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                        * fact_energy * coul_ipkj * coul_iqjk * spin_fact
                endif

                if( pqispin == jkspin ) then

                  selfenergy_sox(iomega,pstate,pqispin) = selfenergy_sox(iomega,pstate,pqispin) &
                           - fact_comp * coul_ipkj * coul_ijkq

                  if(iomega==0 .AND. occupation(pstate,pqispin)>completely_empty) then
                    emp2_sox = emp2_sox - occupation(pstate,pqispin) &
                              * fact_energy * coul_ipkj * coul_ijkq
                  endif

                endif


              enddo ! iomega

            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo
  enddo ! pqispin

  call ortho%sum(selfenergy_ring)
  call ortho%sum(selfenergy_sox)
  call ortho%sum(emp2_ring)
  call ortho%sum(emp2_sox)

  emp2_ring = 0.5_dp * emp2_ring
  emp2_sox  = 0.5_dp * emp2_sox

  if( selfenergy_approx == ONE_RING ) then
    emp2_sox = 0.0_dp
    selfenergy_sox(:,:,:) = 0.0_dp
  endif
  if( selfenergy_approx == SOX ) then
    emp2_ring = 0.0_dp
    selfenergy_ring(:,:,:) = 0.0_dp
  endif

  if( nsemin <= ncore_G+1 .AND. nsemax >= nhomo_G ) then
    emp2 = emp2_ring + emp2_sox
    write(stdout,'(/,a)')       ' MP2 Energy'
    write(stdout,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
    write(stdout,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
    write(stdout,'(a,f14.8,/)') ' MP2 correlation :',emp2
  else
    emp2 = 0.0_dp
  endif

  !$OMP PARALLEL
  !$OMP WORKSHARE
  se%sigma(:,:,:) = selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:)
  !$OMP END WORKSHARE
  !$OMP END PARALLEL

  write(stdout,'(/,1x,a)') ' Spin  State      1-ring             SOX              PT2'
  do pqispin=1,nspin
    do pstate=nsemin,nsemax
      write(stdout,'(1x,i2,2x,i4,*(2x,f12.5))') pqispin,pstate,&
                                                REAL(selfenergy_ring(0,pstate,pqispin),dp)*Ha_eV,&
                                                REAL(selfenergy_sox(0,pstate,pqispin),dp)*Ha_eV,&
                                                REAL(se%sigma(0,pstate,pqispin),dp)*Ha_eV

    enddo
  enddo

  if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
  deallocate(selfenergy_ring)
  deallocate(selfenergy_sox)
  if(has_auxil_basis) call destroy_eri_3center_eigen()

  call stop_clock(timing_pt_self)

end subroutine pt2_selfenergy


!=========================================================================
subroutine onering_selfenergy(nstate,basis,occupation,energy,c_matrix,se,emp2)
  use m_definitions
  use m_mpi
  use m_warning
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_spectral_function
  use m_selfenergy_tools
  implicit none

  integer,intent(in)         :: nstate
  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
  type(selfenergy_grid),intent(inout) :: se
  real(dp),intent(out)       :: emp2
  !=====
  type(spectral_function) :: vchi0v
  !=====

  call start_clock(timing_pt_self)

  if( .NOT. has_auxil_basis ) &
    call die('onering_selfenergy: only implemented when an auxiliary basis is available')

  emp2 = 0.0_dp


  write(stdout,'(/,a)') ' Perform the one-ring self-energy calculation'
  write(stdout,*) 'with the perturbative approach'

  call init_spectral_function(nstate,occupation,0,vchi0v)

  call polarizability_onering(basis,nstate,energy,c_matrix,vchi0v)

#if defined(HAVE_SCALAPACK)
  call gw_selfenergy_scalapack(ONE_RING,nstate,basis,occupation,energy,c_matrix,vchi0v,se)
#else
  call gw_selfenergy(ONE_RING,nstate,basis,occupation,energy,c_matrix,vchi0v,se)
#endif

  call destroy_spectral_function(vchi0v)

  call stop_clock(timing_pt_self)


end subroutine onering_selfenergy


!=========================================================================
subroutine pt2_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,selfenergy,emp2)
  use m_definitions
  use m_mpi
  use m_warning
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_selfenergy_tools
  implicit none

  integer,intent(in)         :: nstate
  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
  real(dp),intent(in)        :: s_matrix(basis%nbf,basis%nbf)
  real(dp),intent(out)       :: selfenergy(basis%nbf,basis%nbf,nspin)
  real(dp),intent(out)       :: emp2
  !=====
  integer                 :: pstate,qstate
  complex(dp),allocatable :: selfenergy_ring(:,:,:)
  complex(dp),allocatable :: selfenergy_sox(:,:,:)
  integer                 :: istate,jstate,kstate
  integer                 :: pqispin,jkspin
  real(dp)                :: fact_occ1,fact_occ2
  real(dp)                :: fi,fj,fk,ei,ej,ek,ep,eq
  complex(dp)             :: fact_comp
  real(dp)                :: fact_energy
  real(dp)                :: emp2_sox,emp2_ring
  real(dp),allocatable    :: eri_eigenstate_i(:,:,:,:)
  real(dp)                :: coul_iqjk,coul_ijkq,coul_ipkj
  !=====

  call start_clock(timing_pt_self)

  emp2_ring = 0.0_dp
  emp2_sox  = 0.0_dp


  write(stdout,'(/,a)') ' Perform the second-order self-energy calculation'
  write(stdout,*) 'with the QP self-consistent approach'



  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
  else
    allocate(eri_eigenstate_i(nstate,nstate,nstate,nspin))
  endif



  allocate(selfenergy_ring (nsemin:nsemax,nsemin:nsemax,nspin))
  allocate(selfenergy_sox  (nsemin:nsemax,nsemin:nsemax,nspin))


  selfenergy_ring(:,:,:) = 0.0_dp
  selfenergy_sox(:,:,:)  = 0.0_dp

  do pqispin=1,nspin
    do istate=ncore_G+1,nvirtual_G-1 !LOOP of the first Green's function
      if( MODULO( istate - (ncore_G+1) , ortho%nproc ) /= ortho%rank ) cycle

      if( .NOT. has_auxil_basis ) then
        call calculate_eri_4center_eigen(c_matrix,istate,pqispin,eri_eigenstate_i)
      endif

      fi = occupation(istate,pqispin)
      ei = energy(istate,pqispin)

      !$OMP PARALLEL
      !$OMP DO PRIVATE(fj,ej,fk,ek,fact_occ1,fact_occ2,coul_ipkj,coul_iqjk,coul_ijkq,ep,eq,fact_comp,fact_energy)   &
      !$OMP REDUCTION(+:emp2_ring,emp2_sox) COLLAPSE(2)
      do pstate=nsemin,nsemax ! external loop ( bra )
        do qstate=nsemin,nsemax   ! external loop ( ket )

          do jkspin=1,nspin
            do jstate=ncore_G+1,nvirtual_G-1  !LOOP of the second Green's function
              fj = occupation(jstate,jkspin)
              ej = energy(jstate,jkspin)

              do kstate=ncore_G+1,nvirtual_G-1 !LOOP of the third Green's function
                fk = occupation(kstate,jkspin)
                ek = energy(kstate,jkspin)

                fact_occ1 = (spin_fact-fi) *            fj  * (spin_fact-fk) / spin_fact**3
                fact_occ2 =            fi  * (spin_fact-fj) *            fk  / spin_fact**3

                if( fact_occ1 < completely_empty .AND. fact_occ2 < completely_empty ) cycle

                if( has_auxil_basis ) then
                  coul_ipkj = eri_eigen_ri(istate,pstate,pqispin,kstate,jstate,jkspin)
                  coul_iqjk = eri_eigen_ri(istate,qstate,pqispin,jstate,kstate,jkspin)
                  if( pqispin == jkspin ) then
                    coul_ijkq = eri_eigen_ri(istate,jstate,pqispin,kstate,qstate,pqispin)
                  endif
                else
                  coul_ipkj = eri_eigenstate_i(pstate,kstate,jstate,jkspin)
                  coul_iqjk = eri_eigenstate_i(qstate,jstate,kstate,jkspin)
                  if( pqispin == jkspin ) then
                    coul_ijkq = eri_eigenstate_i(jstate,kstate,qstate,pqispin)
                  endif
                endif

                ep = energy(pstate,pqispin)
                eq = energy(qstate,pqispin)

                fact_comp   = fact_occ1 / ( eq - ei + ej - ek + ieta) &
                            + fact_occ2 / ( eq - ei + ej - ek - ieta)
                fact_energy = REAL( fact_occ1 / ( ep - ei + ej - ek + ieta) , dp )

                selfenergy_ring(pstate,qstate,pqispin) = selfenergy_ring(pstate,qstate,pqispin) &
                         + fact_comp * coul_ipkj * coul_iqjk * spin_fact

                if( pstate == qstate .AND. occupation(pstate,pqispin) > completely_empty ) then
                  emp2_ring = emp2_ring + occupation(pstate,pqispin) &
                                        * fact_energy * coul_ipkj * coul_iqjk * spin_fact
                endif

                if( pqispin == jkspin ) then

                  selfenergy_sox(pstate,qstate,pqispin) = selfenergy_sox(pstate,qstate,pqispin) &
                           - fact_comp * coul_ipkj * coul_ijkq

                  if( pstate == qstate .AND. occupation(pstate,pqispin) > completely_empty ) then
                    emp2_sox = emp2_sox - occupation(pstate,pqispin) &
                              * fact_energy * coul_ipkj * coul_ijkq
                  endif

                endif



              enddo
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo
  enddo ! pqispin

  call ortho%sum(selfenergy_ring)
  call ortho%sum(selfenergy_sox)
  call ortho%sum(emp2_ring)
  call ortho%sum(emp2_sox)

  emp2_ring = 0.5_dp * emp2_ring
  emp2_sox  = 0.5_dp * emp2_sox
  if( nsemin <= ncore_G+1 .AND. nsemax >= nhomo_G ) then
    emp2 = emp2_ring + emp2_sox
    write(stdout,'(/,a)')       ' MP2 Energy'
    write(stdout,'(a,f14.8)')   ' 2-ring diagram  :',emp2_ring
    write(stdout,'(a,f14.8)')   ' SOX diagram     :',emp2_sox
    write(stdout,'(a,f14.8,/)') ' MP2 correlation :',emp2
  else
    emp2 = 0.0_dp
  endif

  !$OMP PARALLEL
  !$OMP WORKSHARE
  selfenergy(:,:,:) = REAL( selfenergy_ring(:,:,:) + selfenergy_sox(:,:,:) ,dp)
  !$OMP END WORKSHARE
  !$OMP END PARALLEL

  call apply_qs_approximation(s_matrix,c_matrix,selfenergy)


  if( ALLOCATED(eri_eigenstate_i) ) deallocate(eri_eigenstate_i)
  deallocate(selfenergy_ring)
  deallocate(selfenergy_sox)
  if(has_auxil_basis) call destroy_eri_3center_eigen()

  call stop_clock(timing_pt_self)

end subroutine pt2_selfenergy_qs


!=========================================================================
subroutine pt3_selfenergy(selfenergy_approx,selfenergy_technique,nstate,basis,occupation,energy,c_matrix,se,emp3)
  use m_definitions
  use m_mpi
  use m_warning
  use m_timing
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_selfenergy_tools
  implicit none

  integer,intent(in)         :: selfenergy_approx,selfenergy_technique
  integer,intent(in)         :: nstate
  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
  type(selfenergy_grid),intent(inout) :: se
  real(dp),intent(out)       :: emp3
  !=====
  integer,parameter       :: ONERING=1,SOX_=2
  integer,parameter       :: Ah=3,Ax=4
  integer,parameter       :: Cd=5,Cx=6
  integer,parameter       :: DRINGS=7,Deh=8
  integer,parameter       :: TWORINGS=9
  integer                 :: pstate,qstate
  complex(dp),allocatable :: selfenergy(:,:,:,:)
  complex(dp),allocatable :: selfenergy1(:),selfenergy2(:),selfenergy0(:)
  real(dp)                :: selfenergy3,selfenergy4
  integer                 :: iomega
  integer                 :: istate,jstate,kstate,lstate
  integer                 :: astate,bstate,cstate,dstate
  integer                 :: pqspin
  complex(dp)             :: omega
  complex(dp)             :: denom1,denom2
  real(dp)                :: num1,num2,num3
  real(dp)                :: num1a,num1b,num2a,num2b,num3a,num3b
  real(dp)                :: eri_pqjk,eri_pjqk,eri_pqbc,eri_pbqc,eri_pqjc,eri_pjqc,eri_pqkb,eri_pkqb
  real(dp)                :: eri_paib,eri_pbia,eri_iajc,eri_ijac,eri_qcjb,eri_qbjc
  real(dp)                :: eri_pija,eri_pjia,eri_qija,eri_qaib
  real(dp)                :: eri_pkla,eri_plka
  real(dp)                :: eri_qjbc,eri_qjkb,eri_qkjb,eri_qbik,eri_qkib
  real(dp)                :: eri_kaib,eri_icja,eri_acib,eri_ikja
  !=====

  call start_clock(timing_pt_self)

  ! Emp3 is not calculated so far
  emp3 = 0.0_dp

  write(stdout,'(/,a)') ' Perform the third-order self-energy calculation'
  select case(TRIM(pt3_a_diagrams))
  case('YES')
    write(stdout,'(1x,a)') 'Include all the 2nd and 3rd order diagrams'
  case('NO')
    write(stdout,'(1x,a)') 'Discard the A diagrams: the PT2 density matrix effect on Hartree and Exchange'
  case('ONLY')
    write(stdout,'(1x,a)') 'Only calculate the A diagrams: the PT2 density matrix effect on Hartree and Exchange'
  case default
    call die('pt3_selfenergy: pt3_a_diagrams option not valid')
  end select

  if( nspin /= 1 ) call die('pt3_selfenergy: only implemented for spin restricted calculations')

  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
  else
    call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
  endif


  allocate(selfenergy(-se%nomega:se%nomega,ONERING:TWORINGS,nsemin:nsemax,nspin))
  allocate(selfenergy0(-se%nomega:se%nomega))
  allocate(selfenergy1(-se%nomega:se%nomega))
  allocate(selfenergy2(-se%nomega:se%nomega))

  write(stdout,'(/,2x,a5,*(a11,3x))') 'state','PT2','1-ring','SOX', &
                                              'PT3','A diagrams','A hartree','A exch', &
                                              'C diagrams','ladder','ladder X', &
                                              'D diagrams','2-rings','2-rings & X','el-hole & X'

  selfenergy(:,:,:,:) = 0.0_dp

  pqspin = 1
  do pstate=nsemin,nsemax
    qstate = pstate

    ! Check for degenerate states
    if( pstate > nsemin .AND. ABS(energy(pstate,pqspin)-energy(pstate-1,pqspin)) < 1.0e-5_dp ) then
      write(stdout,*) 'State ',pstate,'is degenerated with state',pstate-1,'=> skip calculation'
      selfenergy(:,:,pstate,pqspin) = selfenergy(:,:,pstate-1,pqspin)
      cycle
    endif

    !$OMP PARALLEL &
    !$OMP& PRIVATE(eri_pqjk,eri_pjqk,eri_pqbc,eri_pbqc,eri_pqjc,eri_pjqc,eri_pqkb,eri_pkqb, &
    !$OMP&         eri_paib,eri_pbia,eri_iajc,eri_ijac,eri_qcjb,eri_qbjc, &
    !$OMP&         eri_pija,eri_pjia,eri_qija,eri_qaib,eri_pkla,eri_plka, &
    !$OMP&         eri_qjbc,eri_qjkb,eri_qkjb,eri_qbik,eri_qkib, &
    !$OMP&         eri_kaib,eri_icja,eri_acib,eri_ikja, &
    !$OMP&         denom1,denom2,num1,num2,num3,num1a,num1b,num2a,num2b,num3a,num3b,omega)

    !
    ! A diagrams family
    !
    if( pt3_a_diagrams == 'ONLY' .OR. pt3_a_diagrams == 'YES' ) then


      do astate=nhomo_G+1,nvirtual_G-1
        if( MODULO( astate - (nhomo_G+1) , ortho%nproc ) /= ortho%rank ) cycle

        selfenergy3 = 0.0_dp
        selfenergy4 = 0.0_dp

        !$OMP DO REDUCTION(+:selfenergy3,selfenergy4)
        do bstate=nhomo_G+1,nvirtual_G-1

          ! A1   i,j,k   a,b
          do istate=ncore_G+1,nhomo_G
            do jstate=ncore_G+1,nhomo_G
              num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) &
                    - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
              denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)

              do kstate=ncore_G+1,nhomo_G
                denom2 = energy(kstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
                eri_pqjk = eri_eigen(pstate,qstate,pqspin,jstate,kstate,pqspin)
                eri_pjqk = eri_eigen(pstate,jstate,pqspin,kstate,qstate,pqspin)

                eri_kaib = eri_eigen(astate,kstate,pqspin,bstate,istate,pqspin)

                selfenergy3 = selfenergy3 - 2.0_dp * eri_pqjk * num2 * eri_kaib / ( denom1 * denom2 )
                selfenergy4 = selfenergy4 +          eri_pjqk * num2 * eri_kaib / ( denom1 * denom2 )
              enddo
            enddo
          enddo

          ! A2   i,j   a,b,c
          do cstate=nhomo_G+1,nvirtual_G-1
            eri_pqbc = eri_eigen(pstate,qstate,pqspin,cstate,bstate,pqspin)
            eri_pbqc = eri_eigen(pstate,bstate,pqspin,qstate,cstate,pqspin)

            do istate=ncore_G+1,nhomo_G
              do jstate=ncore_G+1,nhomo_G
                denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
                denom2 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(cstate,pqspin)

                num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) &
                      - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
                eri_icja = eri_eigen(istate,cstate,pqspin,jstate,astate,pqspin)

                selfenergy3 = selfenergy3 + 2.0_dp * eri_pqbc * num2 * eri_icja / ( denom1 * denom2 )
                selfenergy4 = selfenergy4 -          eri_pbqc * num2 * eri_icja / ( denom1 * denom2 )
              enddo
            enddo
          enddo

          ! A3,A4   i,j   a,b,c
          do istate=ncore_G+1,nhomo_G
            do jstate=ncore_G+1,nhomo_G
              num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) &
                    - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
              denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)

              do cstate=nhomo_G+1,nvirtual_G-1
                denom2 = energy(jstate,pqspin) - energy(cstate,pqspin)
                eri_pqjc = eri_eigen(pstate,qstate,pqspin,jstate,cstate,pqspin)
                eri_pjqc = eri_eigen(pstate,jstate,pqspin,qstate,cstate,pqspin)

                eri_acib = eri_eigen(astate,cstate,pqspin,istate,bstate,pqspin)

                selfenergy3 = selfenergy3 + 4.0_dp * eri_pqjc * num2 * eri_acib / ( denom1 * denom2 )
                selfenergy4 = selfenergy4 - 2.0_dp * eri_pjqc * num2 * eri_acib / ( denom1 * denom2 )
              enddo
            enddo
          enddo

          ! A5,A6   i,j,k   a,b
          do istate=ncore_G+1,nhomo_G
            do jstate=ncore_G+1,nhomo_G
              num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) &
                    - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
              denom1 = energy(jstate,pqspin) +  energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)

              do kstate=ncore_G+1,nhomo_G
                denom2 = energy(kstate,pqspin) - energy(bstate,pqspin)
                eri_pqkb = eri_eigen(pstate,qstate,pqspin,kstate,bstate,pqspin)
                eri_pkqb = eri_eigen(pstate,kstate,pqspin,bstate,qstate,pqspin)

                eri_ikja = eri_eigen(istate,kstate,pqspin,jstate,astate,pqspin)

                selfenergy3 = selfenergy3 - 4.0_dp * eri_pqkb  * num2 * eri_ikja / ( denom1 * denom2 )
                selfenergy4 = selfenergy4 + 2.0_dp * eri_pkqb  * num2 * eri_ikja / ( denom1 * denom2 )
              enddo
            enddo
          enddo
        enddo ! bstate
        !$OMP END DO
        !$OMP SINGLE
        selfenergy(:,Ah,pstate,pqspin) = selfenergy(:,Ah,pstate,pqspin) + selfenergy3
        selfenergy(:,Ax,pstate,pqspin) = selfenergy(:,Ax,pstate,pqspin) + selfenergy4
        !$OMP END SINGLE

      enddo ! astate

    endif


    if( pt3_a_diagrams == 'NO' .OR. pt3_a_diagrams == 'YES' ) then


      !
      ! B diagrams family (= 2nd order diagrams)
      !
      do astate=nhomo_G+1,nvirtual_G-1
        if( MODULO( astate - (nhomo_G+1) , ortho%nproc ) /= ortho%rank ) cycle

        selfenergy1(:) = (0.0_dp,0.0_dp)
        selfenergy2(:) = (0.0_dp,0.0_dp)

        ! B1 i,j    a
        !$OMP DO COLLAPSE(2) REDUCTION(+:selfenergy1,selfenergy2)
        do istate=ncore_G+1,nhomo_G
          do jstate=ncore_G+1,nhomo_G

            eri_pija = eri_eigen(pstate,istate,pqspin,jstate,astate,pqspin)
            eri_pjia = eri_eigen(pstate,jstate,pqspin,istate,astate,pqspin)
            eri_qija = eri_eigen(qstate,istate,pqspin,jstate,astate,pqspin)
            do iomega=-se%nomega,se%nomega
              omega = energy(qstate,pqspin) + se%omega(iomega)
              denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(jstate,pqspin) - ieta
              selfenergy1(iomega) = selfenergy1(iomega) + 2.0_dp * eri_pija * eri_qija / denom1
              selfenergy2(iomega) = selfenergy2(iomega) -          eri_pjia * eri_qija / denom1
            enddo

          enddo
        enddo
        !$OMP END DO

        !$OMP SINGLE
        selfenergy(:,ONERING,pstate,pqspin) = selfenergy(:,ONERING,pstate,pqspin) + selfenergy1(:)
        selfenergy(:,SOX_,pstate,pqspin)    = selfenergy(:,SOX_,pstate,pqspin)    + selfenergy2(:)
        !$OMP END SINGLE

        selfenergy1(:) = (0.0_dp,0.0_dp)
        selfenergy2(:) = (0.0_dp,0.0_dp)

        ! B2 i    a,b
        !$OMP DO COLLAPSE(2) REDUCTION(+:selfenergy1,selfenergy2)
        do istate=ncore_G+1,nhomo_G
          do bstate=nhomo_G+1,nvirtual_G-1

            eri_paib = eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin)
            eri_pbia = eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)
            eri_qaib = eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin)
            do iomega=-se%nomega,se%nomega
              omega = energy(qstate,pqspin) + se%omega(iomega)
              denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
              selfenergy1(iomega) = selfenergy1(iomega) + 2.0_dp * eri_paib * eri_qaib / denom1
              selfenergy2(iomega) = selfenergy2(iomega) -          eri_pbia * eri_qaib / denom1
            enddo

          enddo
        enddo  ! istate
        !$OMP END DO

        !$OMP SINGLE
        selfenergy(:,ONERING,pstate,pqspin) = selfenergy(:,ONERING,pstate,pqspin) + selfenergy1(:)
        selfenergy(:,SOX_,pstate,pqspin)    = selfenergy(:,SOX_,pstate,pqspin)    + selfenergy2(:)
        !$OMP END SINGLE
      enddo  ! astate

      !
      ! C diagrams family
      !

      ! C1   i   a,b,c,d
      do astate=nhomo_G+1,nvirtual_G-1
        do istate=ncore_G+1,nhomo_G
          if( MODULO( istate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)

          !$OMP DO REDUCTION(+:selfenergy1,selfenergy2)
          do bstate=nhomo_G+1,nvirtual_G-1
            eri_paib = eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin)
            eri_pbia = eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)

            do cstate=nhomo_G+1,nvirtual_G-1
              do dstate=nhomo_G+1,nvirtual_G-1
                num2 = eri_eigen(astate,cstate,pqspin,bstate,dstate,pqspin)
                num3 = eri_eigen(qstate,cstate,pqspin,istate,dstate,pqspin)
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
                  denom2 = omega + energy(istate,pqspin) - energy(cstate,pqspin) - energy(dstate,pqspin) + ieta
                  selfenergy1(iomega) = selfenergy1(iomega) + 2.0_dp * eri_paib * num2 * num3 / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) -          eri_pbia * num2 * num3 / ( denom1 * denom2 )
                enddo
              enddo
            enddo
          enddo ! bstate
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,Cd,pstate,pqspin) = selfenergy(:,Cd,pstate,pqspin) + selfenergy1(:)
          selfenergy(:,Cx,pstate,pqspin) = selfenergy(:,Cx,pstate,pqspin) + selfenergy2(:)
          !$OMP END SINGLE

        enddo
      enddo

      ! C2+C3   i,j,k   a,b
      do astate=nhomo_G+1,nvirtual_G-1
        do istate=ncore_G+1,nhomo_G
          if( MODULO( istate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)
          !$OMP DO REDUCTION(+:selfenergy1,selfenergy2)
          do bstate=nhomo_G+1,nvirtual_G-1
            eri_paib = eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin)
            eri_pbia = eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)

            do jstate=ncore_G+1,nhomo_G
              do kstate=ncore_G+1,nhomo_G
                num2 = eri_eigen(astate,jstate,pqspin,bstate,kstate,pqspin)
                num3 = eri_eigen(qstate,jstate,pqspin,istate,kstate,pqspin)
                denom2 = energy(jstate,pqspin) + energy(kstate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
                  selfenergy1(iomega) = selfenergy1(iomega) + 4.0_dp * eri_paib * num2 * num3 / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) - 2.0_dp * eri_pbia * num2 * num3 / ( denom1 * denom2 )
                enddo
              enddo
            enddo ! bstate
          enddo
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,Cd,pstate,pqspin) = selfenergy(:,Cd,pstate,pqspin) + selfenergy1(:)
          selfenergy(:,Cx,pstate,pqspin) = selfenergy(:,Cx,pstate,pqspin) + selfenergy2(:)
          !$OMP END SINGLE
        enddo
      enddo

      ! C4+C5   i,j   a,b,c
      do astate=nhomo_G+1,nvirtual_G-1
        do istate=ncore_G+1,nhomo_G
          if( MODULO( istate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)
          !$OMP DO REDUCTION(+:selfenergy1,selfenergy2)
          do jstate=ncore_G+1,nhomo_G
            eri_pija = eri_eigen(pstate,istate,pqspin,jstate,astate,pqspin)
            eri_pjia = eri_eigen(pstate,jstate,pqspin,istate,astate,pqspin)
            do bstate=nhomo_G+1,nvirtual_G-1
              do cstate=nhomo_G+1,nvirtual_G-1
                num2 = eri_eigen(istate,bstate,pqspin,jstate,cstate,pqspin)
                num3 = eri_eigen(qstate,bstate,pqspin,astate,cstate,pqspin)
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(jstate,pqspin) - ieta
                  denom2 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(bstate,pqspin) - energy(cstate,pqspin)
                  selfenergy1(iomega) = selfenergy1(iomega) + 4.0_dp * eri_pija * num2 * num3 / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) - 2.0_dp * eri_pjia * num2 * num3 / ( denom1 * denom2 )
                enddo
              enddo
            enddo ! jstate
          enddo
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,Cd,pstate,pqspin) = selfenergy(:,Cd,pstate,pqspin) + selfenergy1(:)
          selfenergy(:,Cx,pstate,pqspin) = selfenergy(:,Cx,pstate,pqspin) + selfenergy2(:)
          !$OMP END SINGLE
        enddo
      enddo

      ! C6   i,j,k,l   a
      do astate=nhomo_G+1,nvirtual_G-1
        do kstate=ncore_G+1,nhomo_G
          if( MODULO( kstate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)
          !$OMP DO REDUCTION(+:selfenergy1,selfenergy2)
          do lstate=ncore_G+1,nhomo_G
            eri_pkla = eri_eigen(pstate,kstate,pqspin,lstate,astate,pqspin)
            eri_plka = eri_eigen(pstate,lstate,pqspin,kstate,astate,pqspin)
            do istate=ncore_G+1,nhomo_G
              do jstate=ncore_G+1,nhomo_G
                num2 = eri_eigen(kstate,istate,pqspin,lstate,jstate,pqspin)
                num3 = eri_eigen(qstate,istate,pqspin,astate,jstate,pqspin)
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(jstate,pqspin) - ieta
                  denom2 = omega + energy(astate,pqspin) - energy(kstate,pqspin) - energy(lstate,pqspin) - ieta
                  ! Minus sign from Domcke-Cederbaum book chapter 1977 (forgotten in von niessen review in 1983)
                  selfenergy1(iomega) = selfenergy1(iomega) - 2.0_dp *  eri_pkla * num2 * num3 / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) +           eri_plka * num2 * num3 / ( denom1 * denom2 )
                enddo
              enddo
            enddo ! lstate
          enddo
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,Cd,pstate,pqspin) = selfenergy(:,Cd,pstate,pqspin) + selfenergy1(:)
          selfenergy(:,Cx,pstate,pqspin) = selfenergy(:,Cx,pstate,pqspin) + selfenergy2(:)
          !$OMP END SINGLE
        enddo
      enddo

      !
      ! D diagrams family
      !
      ! D1   i,j   a,b,c
      do astate=nhomo_G+1,nvirtual_G-1
        do istate=ncore_G+1,nhomo_G
          if( MODULO( istate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy0(:) = (0.0_dp,0.0_dp)
          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)
          !$OMP DO REDUCTION(+:selfenergy0,selfenergy1,selfenergy2)
          do cstate=nhomo_G+1,nvirtual_G-1
            do jstate=ncore_G+1,nhomo_G
              eri_iajc = eri_eigen(istate,astate,pqspin,jstate,cstate,pqspin)
              eri_ijac = eri_eigen(istate,jstate,pqspin,astate,cstate,pqspin)
              do bstate=nhomo_G+1,nvirtual_G-1
                eri_paib = eri_eigen(pstate,astate,pqspin,istate,bstate,pqspin)
                eri_pbia = eri_eigen(pstate,bstate,pqspin,istate,astate,pqspin)
                eri_qcjb = eri_eigen(qstate,cstate,pqspin,jstate,bstate,pqspin)
                eri_qbjc = eri_eigen(qstate,bstate,pqspin,jstate,cstate,pqspin)
                num3a = eri_qcjb - 2.0_dp * eri_qbjc
                num3b = eri_qbjc - 2.0_dp * eri_qcjb
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
                  denom2 = omega + energy(jstate,pqspin) - energy(bstate,pqspin) - energy(cstate,pqspin) + ieta
                  selfenergy0(iomega) = selfenergy0(iomega) + ( eri_pbia * eri_iajc * 4.0_dp * eri_qbjc ) / ( denom1 * denom2 )
                  selfenergy1(iomega) = selfenergy1(iomega) &
                             + (  eri_paib * eri_iajc * num3a &
                                + eri_pbia * eri_iajc * (-2.0_dp) * num3a ) / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) &
                             + (  eri_paib * eri_ijac * num3b &
                                + eri_pbia * eri_ijac * num3a )   / ( denom1 * denom2 )
                enddo
              enddo
            enddo
          enddo
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,TWORINGS,pstate,pqspin) = selfenergy(:,TWORINGS,pstate,pqspin) + selfenergy0(:)
          selfenergy(:,DRINGS,pstate,pqspin)   = selfenergy(:,DRINGS,pstate,pqspin)   + selfenergy1(:)
          selfenergy(:,Deh,pstate,pqspin)      = selfenergy(:,Deh,pstate,pqspin)      + selfenergy2(:)
          !$OMP END SINGLE

        enddo
      enddo


      ! D2+D3   i,j   a,b,c
      do astate=nhomo_G+1,nvirtual_G-1
        do istate=ncore_G+1,nhomo_G
          if( MODULO( istate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy0(:) = (0.0_dp,0.0_dp)
          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)
          !$OMP DO REDUCTION(+:selfenergy0,selfenergy1,selfenergy2)
          do bstate=nhomo_G+1,nvirtual_G-1
            do jstate=ncore_G+1,nhomo_G
              num2a = eri_eigen(astate,istate,pqspin,bstate,jstate,pqspin)
              num2b = eri_eigen(astate,jstate,pqspin,bstate,istate,pqspin)
              do cstate=nhomo_G+1,nvirtual_G-1
                num1a = eri_eigen(pstate,cstate,pqspin,istate,astate,pqspin)
                num1b = eri_eigen(pstate,astate,pqspin,istate,cstate,pqspin)
                eri_qcjb = eri_eigen(qstate,cstate,pqspin,jstate,bstate,pqspin)
                eri_qjbc = eri_eigen(qstate,jstate,pqspin,bstate,cstate,pqspin)
                num3a = eri_qjbc - 2.0_dp * eri_qcjb
                num3b = eri_qcjb - 2.0_dp * eri_qjbc
                denom2 = energy(jstate,pqspin) + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(istate,pqspin) - energy(astate,pqspin) - energy(cstate,pqspin) + ieta
                  selfenergy0(iomega) = selfenergy0(iomega) &
                             + 8.0_dp * num1a * num2a * eri_qcjb / ( denom1 * denom2 )
                  selfenergy1(iomega) = selfenergy1(iomega) &
                             + 2.0_dp * (  num1a * num2a * (-2.0_dp) * num3a + num1b * num2a * num3a ) / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) &
                             + 2.0_dp * (  num1a * num2b * num3a + num1b * num2b * num3b )   / ( denom1 * denom2 )
                enddo
              enddo
            enddo
          enddo
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,TWORINGS,pstate,pqspin) = selfenergy(:,TWORINGS,pstate,pqspin) + selfenergy0(:)
          selfenergy(:,DRINGS,pstate,pqspin)   = selfenergy(:,DRINGS,pstate,pqspin)   + selfenergy1(:)
          selfenergy(:,Deh,pstate,pqspin)      = selfenergy(:,Deh,pstate,pqspin)      + selfenergy2(:)
          !$OMP END SINGLE

        enddo
      enddo


      ! D4+D5   i,j,k   a,b
      do astate=nhomo_G+1,nvirtual_G-1
        do istate=ncore_G+1,nhomo_G
          if( MODULO( istate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy0(:) = (0.0_dp,0.0_dp)
          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)
          !$OMP DO REDUCTION(+:selfenergy0,selfenergy1,selfenergy2)
          do bstate=nhomo_G+1,nvirtual_G-1
            do jstate=ncore_G+1,nhomo_G
              num2a = eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin)
              num2b = eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
              do kstate=ncore_G+1,nhomo_G
                num1a = eri_eigen(pstate,kstate,pqspin,astate,jstate,pqspin)
                num1b = eri_eigen(pstate,jstate,pqspin,astate,kstate,pqspin)
                eri_qbik = eri_eigen(qstate,bstate,pqspin,istate,kstate,pqspin)
                eri_qkib = eri_eigen(qstate,kstate,pqspin,istate,bstate,pqspin)
                num3a = eri_qbik - 2.0_dp * eri_qkib
                num3b = eri_qkib - 2.0_dp * eri_qbik
                denom2 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin) + ieta
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(astate,pqspin) - energy(jstate,pqspin) - energy(kstate,pqspin) - ieta
                  selfenergy0(iomega) = selfenergy0(iomega) &
                             + 8.0_dp * num1a * num2a * eri_qkib / ( denom1 * denom2 )
                  selfenergy1(iomega) = selfenergy1(iomega) &
                             + 2.0_dp * (  num1a * num2a * (-2.0_dp) * num3a  &
                                         + num1b * num2a * num3a              ) / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) &
                             + 2.0_dp * (  num1a * num2b * num3a &
                                         + num1b * num2b * num3b )   / ( denom1 * denom2 )
                enddo
              enddo
            enddo
          enddo
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,TWORINGS,pstate,pqspin) = selfenergy(:,TWORINGS,pstate,pqspin) + selfenergy0(:)
          selfenergy(:,DRINGS,pstate,pqspin)   = selfenergy(:,DRINGS,pstate,pqspin)   + selfenergy1(:)
          selfenergy(:,Deh,pstate,pqspin)      = selfenergy(:,Deh,pstate,pqspin)      + selfenergy2(:)
          !$OMP END SINGLE
        enddo
      enddo

      ! D6   i,j,k   a,b
      do astate=nhomo_G+1,nvirtual_G-1
        do istate=ncore_G+1,nhomo_G
          if( MODULO( istate + (astate-1) * (nhomo_G-ncore_G) , ortho%nproc ) /= ortho%rank ) cycle

          selfenergy0(:) = (0.0_dp,0.0_dp)
          selfenergy1(:) = (0.0_dp,0.0_dp)
          selfenergy2(:) = (0.0_dp,0.0_dp)
          !$OMP DO REDUCTION(+:selfenergy0,selfenergy1,selfenergy2)
          do bstate=nhomo_G+1,nvirtual_G-1
            do jstate=ncore_G+1,nhomo_G
              num2a = eri_eigen(istate,astate,pqspin,bstate,jstate,pqspin)
              num2b = eri_eigen(istate,jstate,pqspin,bstate,astate,pqspin)

              do kstate=ncore_G+1,nhomo_G
                num1a = eri_eigen(pstate,kstate,pqspin,astate,istate,pqspin)
                num1b = eri_eigen(pstate,istate,pqspin,astate,kstate,pqspin)
                eri_qjkb = eri_eigen(qstate,jstate,pqspin,kstate,bstate,pqspin)
                eri_qkjb = eri_eigen(qstate,kstate,pqspin,jstate,bstate,pqspin)
                num3a = eri_qjkb - 2.0_dp * eri_qkjb
                num3b = eri_qkjb - 2.0_dp * eri_qjkb
                do iomega=-se%nomega,se%nomega
                  omega = energy(qstate,pqspin) + se%omega(iomega)
                  denom1 = omega + energy(astate,pqspin) - energy(istate,pqspin) - energy(kstate,pqspin) - ieta
                  denom2 = omega + energy(bstate,pqspin) - energy(jstate,pqspin) - energy(kstate,pqspin) - ieta
                  selfenergy0(iomega) = selfenergy0(iomega) - 4.0_dp * num1a * num2a * eri_qkjb / ( denom1 * denom2 )
                  selfenergy1(iomega) = selfenergy1(iomega) &
                             - ( num1a * num2a * (-2.0_dp) * num3a + num1b * num2a * num3a ) / ( denom1 * denom2 )
                  selfenergy2(iomega) = selfenergy2(iomega) &
                             - ( num1a * num2b * num3a + num1b * num2b * num3b ) / ( denom1 * denom2 )
                enddo
              enddo
            enddo
          enddo
          !$OMP END DO
          !$OMP SINGLE
          selfenergy(:,TWORINGS,pstate,pqspin) = selfenergy(:,TWORINGS,pstate,pqspin) + selfenergy0(:)
          selfenergy(:,DRINGS,pstate,pqspin)   = selfenergy(:,DRINGS,pstate,pqspin)   + selfenergy1(:)
          selfenergy(:,Deh,pstate,pqspin)      = selfenergy(:,Deh,pstate,pqspin)      + selfenergy2(:)
          !$OMP END SINGLE
        enddo
      enddo

    endif
    !$OMP END PARALLEL

    call ortho%sum(selfenergy(:,:,pstate,:))

    write(stdout,'(i4,*(2x,f12.4))') pstate, &
                                     SUM(REAL(selfenergy(0,ONERING:SOX_,pstate,pqspin),dp),DIM=1) * Ha_eV, &
                                     REAL(selfenergy(0,ONERING,pstate,pqspin),dp) * Ha_eV,          &
                                     REAL(selfenergy(0,SOX_,pstate,pqspin),dp) * Ha_eV,             &
                                     SUM(REAL(selfenergy(0,Ah:Deh,pstate,:),dp),DIM=1) * Ha_eV,      &
                                     SUM(REAL(selfenergy(0,Ah:Ax,pstate,pqspin),dp),DIM=1) * Ha_eV, &
                                     REAL(selfenergy(0,Ah,pstate,pqspin),dp) * Ha_eV,               &
                                     REAL(selfenergy(0,Ax,pstate,pqspin),dp) * Ha_eV,               &
                                     SUM(REAL(selfenergy(0,Cd:Cx,pstate,pqspin),dp),DIM=1) * Ha_eV, &
                                     REAL(selfenergy(0,Cd,pstate,pqspin),dp) * Ha_eV,               &
                                     REAL(selfenergy(0,Cx,pstate,pqspin),dp) * Ha_eV,               &
                                     SUM(REAL(selfenergy(0,DRINGS:Deh,pstate,pqspin),dp),DIM=1) * Ha_eV, &
                                     REAL(selfenergy(0,TWORINGS,pstate,pqspin),dp) * Ha_eV,         &
                                     REAL(selfenergy(0,DRINGS,pstate,pqspin),dp) * Ha_eV,           &
                                     REAL(selfenergy(0,Deh,pstate,pqspin),dp) * Ha_eV
  enddo


  select case(selfenergy_approx)
  case(PT3)
    se%sigma(:,:,:) = SUM(selfenergy(:,ONERING:Deh,:,:),DIM=2)
  case(GWPT3)
    se%sigma(:,:,:) = SUM(selfenergy(:,SOX_:Deh,:,:),DIM=2) - selfenergy(:,TWORINGS,:,:)
  case(TWO_RINGS)
    se%sigma(:,:,:) = selfenergy(:,ONERING,:,:) + selfenergy(:,TWORINGS,:,:)
  case default
    call die('pt3_selfenergy: invalid choice of diagrams')
  end select

  deallocate(selfenergy)
  deallocate(selfenergy1,selfenergy2,selfenergy0)

  if(has_auxil_basis) then
    call destroy_eri_3center_eigen()
  else
    call destroy_eri_4center_eigen_uks()
  endif

  call stop_clock(timing_pt_self)

end subroutine pt3_selfenergy


!=========================================================================
