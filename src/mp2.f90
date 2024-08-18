!==================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! - MP2 total energy with or without Resolution-of-Identity
! - Single excitation contribution to total energy
!=========================================================================
#include "molgw.h"
subroutine mp2_energy_ri(nstate,basis,occupation,energy,c_matrix,emp2)
  use m_definitions
  use m_mpi
  use m_cart_to_pure
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore,kappa_hybrid
  implicit none

  integer,intent(in)           :: nstate
  type(basis_set),intent(in)   :: basis
  real(dp),intent(in)          :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)          :: c_matrix(basis%nbf,nstate,nspin)
  real(dp),intent(out)         :: emp2
  !====
  integer                    :: astate,bstate,istate,jstate
  integer                    :: iaspin,jbspin
  real(dp)                   :: energy_denom
  real(dp)                   :: tmp_iajb,tmp_ibja
  real(dp)                   :: contrib1,contrib2
  real(dp)                   :: fact
  integer                    :: nocc(nspin)
  integer                    :: ncore,nstate_mp2
  !=====

  call start_clock(timing_mp2_energy)

  write(stdout,'(/,a)') ' RI-MP2 correlation calculation'


  ncore = ncoreg
  if(is_frozencore) then
    if( ncore == 0) ncore = atoms_core_states()
  endif

  call calculate_eri_3center_eigen(c_matrix,ncore+1,nstate,ncore+1,nstate)

  nstate_mp2 = MIN( nvirtualg-1, nstate )



  emp2 = 0.0_dp
  contrib1 = 0.0_dp
  contrib2 = 0.0_dp


  do iaspin=1,nspin
    !
    ! First, set up the list of occupied states
    nocc(iaspin) = ncore
    do istate=ncore+1,nstate
      if( occupation(istate,iaspin) < completely_empty ) cycle
      nocc(iaspin) = istate
    enddo
  enddo


  do iaspin=1,nspin

    do istate=ncore+1,nocc(iaspin)


      write(stdout,'(i4,2x,i4,a,i4)') iaspin,istate-ncore,' / ',nocc(iaspin)-ncore

      do jbspin=1,nspin

        do jstate=ncore+1,nocc(jbspin)

          do astate=ncore+1,nstate_mp2
            if( occupation(astate,iaspin) > spin_fact - completely_empty ) cycle

            do bstate=ncore+1,nstate_mp2
              if( occupation(bstate,jbspin) > spin_fact - completely_empty ) cycle

              fact =  occupation(istate,iaspin) * ( spin_fact - occupation(astate,iaspin) ) &
                    *occupation(jstate,jbspin) * ( spin_fact - occupation(bstate,jbspin) ) / spin_fact**2

              energy_denom = energy(istate,iaspin) + energy(jstate,jbspin) &
                                     - energy(astate,iaspin) - energy(bstate,jbspin)
              ! Avoid the zero denominators
              if( ABS(energy_denom) < 1.d-18) then
                write(stdout,*) 'you skipped something'
                cycle
              endif

              energy_denom =  fact / energy_denom

              tmp_iajb = eri_eigen_ri(istate,astate,iaspin,jstate,bstate,jbspin)

              contrib1 = contrib1 + 0.5_dp * energy_denom * tmp_iajb**2

              if(iaspin==jbspin) then
                tmp_ibja = eri_eigen_ri(istate,bstate,iaspin,jstate,astate,jbspin)
                contrib2 = contrib2 - 0.5_dp * energy_denom * tmp_iajb*tmp_ibja / spin_fact
              endif

            enddo
          enddo
        enddo
      enddo !jbspin

    enddo ! istate

  enddo !iaspin

  if(kappa_hybrid/=zero) then
    write(stdout,'(/,a,f16.10)')       ' MP2 contributions will be scaled by :',kappa_hybrid
    contrib1=kappa_hybrid*contrib1
    contrib2=kappa_hybrid*contrib2
  endif

  emp2 = contrib1 + contrib2
  write(stdout,'(/,a)')       ' MP2 contributions'
  write(stdout,'(a,f16.10)')   ' 2-ring diagram  :',contrib1
  write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib2
  write(stdout,'(a,f16.10,/)') ' MP2 correlation :',emp2

  call destroy_eri_3center_eigen()
  call stop_clock(timing_mp2_energy)

end subroutine mp2_energy_ri


!=========================================================================
subroutine mp2_energy_ri_cmplx(nstate,basis,occupation,energy,c_matrix_cmplx,emp2)
  use m_definitions
  use m_mpi
  use m_cart_to_pure
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore,kappa_hybrid
  implicit none

  integer,intent(in)           :: nstate
  type(basis_set),intent(in)   :: basis
  real(dp),intent(in)          :: occupation(nstate,nspin),energy(nstate,nspin)
  complex(dp),intent(in)       :: c_matrix_cmplx(basis%nbf,nstate,nspin)
  real(dp),intent(out)         :: emp2
  !====
  integer                    :: astate,bstate,istate,jstate
  integer                    :: iaspin,jbspin
  real(dp)                   :: energy_denom
  real(dp)                   :: fact
  complex(dp)                :: tmp_iajb,tmp_ibja
  complex(dp)                :: contrib1,contrib2
  integer                    :: nocc(nspin)
  integer                    :: ncore,nstate_mp2
  !=====

  call start_clock(timing_mp2_energy)

  write(stdout,'(/,a)') ' RI-MP2 correlation calculation'

  ncore = ncoreg
  if(is_frozencore) then
    if( ncore == 0) ncore = atoms_core_states()
  endif

  call calculate_eri_3center_eigen_cmplx(c_matrix_cmplx,ncore+1,nstate,ncore+1,nstate)

  nstate_mp2 = MIN( nvirtualg-1, nstate )

  emp2 = 0.0_dp
  contrib1 = ( 0.0_dp, 0.0_dp )
  contrib2 = ( 0.0_dp, 0.0_dp )

  do iaspin=1,nspin
    !
    ! First, set up the list of occupied states
    nocc(iaspin) = ncore
    do istate=ncore+1,nstate
      if( occupation(istate,iaspin) < completely_empty ) cycle
      nocc(iaspin) = istate
    enddo
  enddo

  do iaspin=1,nspin

    do istate=ncore+1,nocc(iaspin)


      write(stdout,'(i4,2x,i4,a,i4)') iaspin,istate-ncore,' / ',nocc(iaspin)-ncore

      do jbspin=1,nspin

        do jstate=ncore+1,nocc(jbspin)

          do astate=ncore+1,nstate_mp2
            if( occupation(astate,iaspin) > spin_fact - completely_empty ) cycle

            do bstate=ncore+1,nstate_mp2
              if( occupation(bstate,jbspin) > spin_fact - completely_empty ) cycle

              fact =  occupation(istate,iaspin) * ( spin_fact - occupation(astate,iaspin) ) &
                    *occupation(jstate,jbspin) * ( spin_fact - occupation(bstate,jbspin) ) / spin_fact**2

              energy_denom = energy(istate,iaspin) + energy(jstate,jbspin) &
                                     - energy(astate,iaspin) - energy(bstate,jbspin)
              ! Avoid the zero denominators
              if( ABS(energy_denom) < 1.d-18) then
                write(stdout,*) 'you skipped something'
                cycle
              endif

              energy_denom =  fact / energy_denom

              tmp_iajb = eri_eigen_ri_cmplx(istate,astate,iaspin,jstate,bstate,jbspin)

              contrib1 = contrib1 + 0.5_dp * energy_denom * conjg(tmp_iajb)*tmp_iajb

              if(iaspin==jbspin) then
                tmp_ibja = eri_eigen_ri_cmplx(istate,bstate,iaspin,jstate,astate,jbspin)
                contrib2 = contrib2 - 0.5_dp * energy_denom * conjg(tmp_iajb)*tmp_ibja / spin_fact
              endif

            enddo
          enddo
        enddo
      enddo !jbspin

    enddo ! istate

  enddo !iaspin

  if(kappa_hybrid/=zero) then
    write(stdout,'(/,a,f16.10)')       ' MP2 contributions will be scaled by :',kappa_hybrid
    contrib1=kappa_hybrid*contrib1
    contrib2=kappa_hybrid*contrib2
  endif

  emp2 = real(contrib1 + contrib2)
  write(stdout,'(/,a)')       ' MP2 contributions'
  write(stdout,'(a,f16.10)')   ' 2-ring diagram  :',real(contrib1)
  write(stdout,'(a,f16.10)')   ' SOX diagram     :',real(contrib2)
  write(stdout,'(a,f16.10,/)') ' MP2 correlation :',emp2

  call destroy_eri_3center_eigen_cmplx()
  call stop_clock(timing_mp2_energy)

end subroutine mp2_energy_ri_cmplx


!=========================================================================
subroutine mp2_energy_ri_x2c(nstate,nocc,basis,energy,c_matrix_rel,emp2,exx)
  use m_definitions
  use m_mpi
  use m_cart_to_pure
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore,kappa_hybrid
  implicit none

  integer,intent(in)           :: nstate,nocc
  type(basis_set),intent(in)   :: basis
  real(dp),intent(in)          :: energy(basis%nbf,nspin)
  complex(dp),intent(in)       :: c_matrix_rel(nstate,nstate)
  real(dp),intent(out)         :: emp2,exx
  !====
  integer                    :: astate,bstate,istate,jstate
  real(dp)                   :: fact,energy_denom
  complex(dp)                :: tmp_iajb,tmp_ibja
  complex(dp)                :: contrib1
  integer                    :: ncore
  real(dp),allocatable       :: energy_vec(:)
  !=====

  call start_clock(timing_mp2_energy)

  call calculate_eri_x2c(c_matrix_rel,nstate)

  exx = 0.0_dp
  do istate=1,nocc
    do jstate=1,nocc
      exx=exx-real(eri_eigen_ri_x2c(istate,jstate,jstate,istate),dp)
    enddo
  enddo
  exx = 0.5_dp*exx


  write(stdout,'(/,a)') ' X2C RI-MP2 correlation calculation'

  ncore = ncoreg
  if(is_frozencore) then
    if( ncore == 0) ncore = atoms_core_states()
  endif
  ncore = 2*ncore

  emp2 = 0.0_dp
  contrib1 = ( 0.0_dp, 0.0_dp )
  allocate(energy_vec(nstate))

  do istate=1,nstate/2
    energy_vec(2*istate-1)=energy(istate,1)
    energy_vec(2*istate  )=energy(istate,2)
  enddo

  do istate=ncore+1,nocc
        
    write(stdout,'(2x,i4,a,i4)') istate-ncore,' / ',nocc-ncore

    do jstate=ncore+1,nocc

      do astate=nocc+1,nstate

        do bstate=nocc+1,nstate

           energy_denom = energy_vec(istate) + energy_vec(jstate) &
                        - energy_vec(astate) - energy_vec(bstate)
           ! Avoid the zero denominators
           if( ABS(energy_denom) < 1.d-18) then
             write(stdout,*) 'you skipped something'
             cycle
           endif

           tmp_iajb = eri_eigen_ri_x2c(istate,astate,jstate,bstate)
           tmp_ibja = eri_eigen_ri_x2c(istate,bstate,jstate,astate)
      
           contrib1 = contrib1 + conjg(tmp_iajb-tmp_ibja)*(tmp_iajb-tmp_ibja)/energy_denom

        enddo

      enddo

    enddo

  enddo

  if(kappa_hybrid/=zero) then
    write(stdout,'(/,a,f16.10)')       ' MP2 contributions will be scaled by :',kappa_hybrid
    contrib1=kappa_hybrid*contrib1
  endif

  emp2 = 0.25_dp*real(contrib1)

  write(stdout,'(/,a)')       ' MP2 contributions'
  write(stdout,'(a,f16.10,/)') ' MP2 correlation :',emp2

  deallocate(energy_vec)
  call destroy_eri_3center_eigen_x2c()
  call stop_clock(timing_mp2_energy)

end subroutine mp2_energy_ri_x2c


!=========================================================================
subroutine mp3_energy_ri(nstate,basis,occupation,energy,c_matrix,emp3)
  use m_definitions
  use m_mpi
  use m_cart_to_pure
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam,only: nspin,spin_fact,ncoreg,nvirtualg,is_frozencore,kappa_hybrid
  implicit none

  integer,intent(in)         :: nstate
  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
  real(dp),intent(out)       :: emp3
  !====
  integer                    :: astate,bstate,cstate,dstate,istate,jstate,kstate,lstate
  integer                    :: iaspin,jbspin,kcspin
  integer                    :: cspin,dspin,kspin,lspin
  real(dp)                   :: energy_denom
  real(dp)                   :: contrib1,contrib2,contrib3
  real(dp)                   :: denom1,denom2,numer1,numer2
  integer                    :: nocc(nspin)
  integer                    :: ncore,nstate_mp3
  real(dp)                   :: t_ijab_tilde,x_ijab,t_ijcd,t_klab,t_kjac,t_kiac,t_ikac
  !=====

  call start_clock(timing_mp2_energy)

  write(stdout,'(/,a)') ' RI-MP3 correlation calculation'

  if( nspin > 1 ) call die('MP3 not implemented for unrestricted calculations')

  ncore = ncoreg
  if(is_frozencore) then
    if( ncore == 0) ncore = atoms_core_states()
  endif

  call calculate_eri_3center_eigen(c_matrix,ncore+1,nstate,ncore+1,nstate)

  nstate_mp3 = MIN( nvirtualg-1, nstate )



  emp3 = 0.0_dp
  contrib1 = 0.0_dp
  contrib2 = 0.0_dp
  contrib3 = 0.0_dp


  do iaspin=1,nspin
    !
    ! First, set up the list of occupied states
    nocc(iaspin) = ncore
    do istate=ncore+1,nstate
      if( occupation(istate,iaspin) < completely_empty ) cycle
      nocc(iaspin) = istate
    enddo
  enddo


  ! From Helgaker's book
  write(stdout,*) 'From Helgaker book'
  iaspin=1
  jbspin=1

  do istate=ncore+1,nocc(iaspin)
    do astate=nocc(iaspin)+1,nstate_mp3
      do jstate=ncore+1,nocc(iaspin)
        do bstate=nocc(iaspin)+1,nstate_mp3

          t_ijab_tilde = - 2.0_dp * ( 2.0_dp * eri_eigen_ri(istate,astate,iaspin,jstate,bstate,jbspin)  &
                                     - eri_eigen_ri(istate,bstate,iaspin,jstate,astate,jbspin) ) &
                                       / ( energy(astate,iaspin) + energy(bstate,jbspin) &
                                          - energy(istate,iaspin) - energy(jstate,jbspin) )

          !FIXME: Commenting this is inconsistent with Helgaker however yields the correct numerical results
          !if( istate == jstate .AND. astate == bstate ) t_ijab_tilde = t_ijab_tilde * 2.0_dp

          x_ijab = 0.0_dp
          !
          ! Contrib1
          do cstate=nocc(iaspin)+1,nstate_mp3
            do dstate=nocc(iaspin)+1,nstate_mp3

              t_ijcd = - eri_eigen_ri(istate,cstate,iaspin,jstate,dstate,jbspin)  &
                          / ( energy(cstate,iaspin) + energy(dstate,jbspin) - energy(istate,iaspin) - energy(jstate,jbspin) )

              x_ijab = x_ijab + 0.5_dp * eri_eigen_ri(astate,cstate,iaspin,bstate,dstate,jbspin) * t_ijcd

            enddo
          enddo

          !
          ! Contrib2
          do kstate=ncore+1,nocc(iaspin)
            do lstate=ncore+1,nocc(iaspin)

              t_klab = - eri_eigen_ri(kstate,astate,iaspin,lstate,bstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(bstate,jbspin) - energy(kstate,iaspin) - energy(lstate,jbspin) )

              x_ijab = x_ijab + 0.5_dp * eri_eigen_ri(kstate,istate,iaspin,lstate,jstate,jbspin) * t_klab

            enddo
          enddo

          !
          ! Contrib3
          do kstate=ncore+1,nocc(iaspin)
            do cstate=nocc(iaspin)+1,nstate_mp3

              t_ikac = - eri_eigen_ri(istate,astate,iaspin,kstate,cstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(cstate,jbspin) - energy(istate,iaspin) - energy(kstate,jbspin) )

              x_ijab = x_ijab + ( 2.0_dp * eri_eigen_ri(bstate,jstate,iaspin,kstate,cstate,jbspin) &
                                  - eri_eigen_ri(bstate,cstate,iaspin,kstate,jstate,jbspin) ) * t_ikac


              t_kjac = - eri_eigen_ri(kstate,astate,iaspin,jstate,cstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(cstate,jbspin) - energy(kstate,iaspin) - energy(jstate,jbspin) )
              x_ijab = x_ijab - eri_eigen_ri(bstate,cstate,iaspin,kstate,istate,jbspin) * t_kjac

              t_kiac = - eri_eigen_ri(kstate,astate,iaspin,istate,cstate,jbspin)  &
                          / ( energy(astate,iaspin) + energy(cstate,jbspin) - energy(kstate,iaspin) - energy(istate,jbspin) )
              x_ijab = x_ijab - eri_eigen_ri(bstate,jstate,iaspin,kstate,cstate,jbspin) * t_kiac


            enddo
          enddo

          contrib1 = contrib1 + t_ijab_tilde * x_ijab

        enddo
      enddo
    enddo
  enddo

  if(kappa_hybrid/=zero) then
    write(stdout,'(/,a,f16.10)')       ' MP3 contributions will be scaled by :',kappa_hybrid
    contrib1=kappa_hybrid*contrib1
    contrib2=kappa_hybrid*contrib2
    contrib3=kappa_hybrid*contrib3
  endif

  emp3 = contrib1 + contrib2 + contrib3
  ! write(stdout,'(/,a)')       ' MP3 contributions'
  ! write(stdout,'(a,f16.10)')   ' 2-ring diagram  :',contrib1
  ! write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib2
  ! write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib3
  ! write(stdout,'(a,f16.10,/)') ' MP3 correlation :',emp3

  call destroy_eri_3center_eigen()
  call stop_clock(timing_mp2_energy)

end subroutine mp3_energy_ri


!==================================================================
subroutine mp2_energy(nstate,basis,occupation,c_matrix,energy,emp2)
  use m_definitions
  use m_mpi
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam,only: nspin,spin_fact,ncoreg,kappa_hybrid
  implicit none

  integer,intent(in)           :: nstate
  type(basis_set),intent(in)   :: basis
  real(dp),intent(in)          :: occupation(nstate,nspin),energy(nstate,nspin)
  real(dp),intent(in)          :: c_matrix(basis%nbf,nstate,nspin)
  real(dp),intent(out)         :: emp2
  !=====
  integer                    :: astate,bstate,istate,jstate
  integer                    :: ibf,jbf,abf,bbf,iaspin,jbspin
  real(dp)                   :: energy_denom
  real(dp)                   :: tmp_ixjx(basis%nbf,basis%nbf)
  real(dp)                   :: tmp_iajx(basis%nbf),tmp_ixja(basis%nbf)
  real(dp)                   :: tmp_iajb,tmp_ibja
  real(dp)                   :: contrib1,contrib2
  real(dp)                   :: fact
  real(dp),allocatable       :: tmp_ixxx(:,:,:)
  integer                    :: nocc
  !=====

  call start_clock(timing_mp2_energy)

  write(stdout,*) 'starting the MP2 calculation'


  emp2 = 0.0_dp
  contrib1 = 0.0_dp
  contrib2 = 0.0_dp


  allocate(tmp_ixxx(basis%nbf,basis%nbf,basis%nbf))

  do iaspin=1,nspin
    !
    ! First, set up the list of occupied states
    nocc = ncoreg
    do istate=ncoreg+1,nstate
      if( occupation(istate,iaspin) < completely_empty ) cycle
      nocc = istate
    enddo


    do istate=ncoreg+1,nocc


      write(stdout,'(i4,2x,i4,a,i4)') iaspin,istate-ncoreg,' / ',nocc

      tmp_ixxx(:,:,:) = 0.0_dp
      do bbf=1,basis%nbf
        do jbf=1,basis%nbf
          if( negligible_basispair(jbf,bbf) ) cycle
          do abf=1,basis%nbf
            do ibf=1,basis%nbf
              if( negligible_basispair(ibf,abf) ) cycle
              tmp_ixxx(abf,jbf,bbf) = tmp_ixxx(abf,jbf,bbf) &
                + c_matrix(ibf,istate,iaspin) * eri(ibf,abf,jbf,bbf)
            enddo
          enddo
        enddo
      enddo

      do jbspin=1,nspin
        do jstate=ncoreg+1,nstate
          if( occupation(jstate,jbspin) < completely_empty ) cycle

          tmp_ixjx(:,:) = 0.0_dp
          do bbf=1,basis%nbf
            do jbf=1,basis%nbf
              do abf=1,basis%nbf
                tmp_ixjx(abf,bbf) = tmp_ixjx(abf,bbf) + c_matrix(jbf,jstate,jbspin) * tmp_ixxx(abf,jbf,bbf)
              enddo
            enddo
          enddo

          do astate=1,nstate
            if( occupation(astate,iaspin) > spin_fact - completely_empty ) cycle

            tmp_iajx(:) = 0.0_dp
            do bbf=1,basis%nbf
              do abf=1,basis%nbf
                tmp_iajx(bbf) = tmp_iajx(bbf) + c_matrix(abf,astate,iaspin) * tmp_ixjx(abf,bbf)
              enddo
            enddo

            if(iaspin==jbspin) then
              tmp_ixja(:) = 0.0_dp
              do abf=1,basis%nbf
                do bbf=1,basis%nbf
                  tmp_ixja(bbf) = tmp_ixja(bbf) + c_matrix(abf,astate,iaspin) * tmp_ixjx(bbf,abf)
                enddo
              enddo
            endif

            do bstate=1,nstate
              if( occupation(bstate,jbspin) > spin_fact - completely_empty ) cycle

              fact =  occupation(istate,iaspin) * ( spin_fact - occupation(astate,iaspin) ) &
                    *occupation(jstate,jbspin) * ( spin_fact - occupation(bstate,jbspin) ) / spin_fact**2

              energy_denom = energy(istate,iaspin) + energy(jstate,jbspin) &
                                     - energy(astate,iaspin) - energy(bstate,jbspin)
              ! Avoid the zero denominators
              if( ABS(energy_denom) < 1.d-18) then
                write(stdout,*) 'you skipped something'
                cycle
              endif

              energy_denom =  fact / energy_denom

              tmp_iajb = SUM( tmp_iajx(:) * c_matrix(:,bstate,jbspin) )

              contrib1 = contrib1 + 0.5_dp * energy_denom * tmp_iajb**2

              if(iaspin==jbspin) then
                tmp_ibja = SUM( tmp_ixja(:) * c_matrix(:,bstate,jbspin) )
                contrib2 = contrib2 - 0.5_dp * energy_denom * tmp_iajb*tmp_ibja / spin_fact
              endif

            enddo
          enddo
        enddo
      enddo !jbspin

    enddo ! istate

  enddo !iaspin

  deallocate(tmp_ixxx)

  if(kappa_hybrid/=zero) then
    write(stdout,'(/,a,f16.10)')       ' MP2 contributions will be scaled by :',kappa_hybrid
    contrib1=kappa_hybrid*contrib1
    contrib2=kappa_hybrid*contrib2
  endif

  emp2 = contrib1 + contrib2
  write(stdout,'(/,a)')       ' MP2 contributions'
  write(stdout,'(a,f16.10)')   ' 2-ring diagram  :',contrib1
  write(stdout,'(a,f16.10)')   ' SOX diagram     :',contrib2
  write(stdout,'(a,f16.10,/)') ' MP2 correlation :',emp2

  call stop_clock(timing_mp2_energy)

end subroutine mp2_energy


!==================================================================
