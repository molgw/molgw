!=========================================================================
#include "macros.h"
!=========================================================================
subroutine gw_selfenergy(gwmethod,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam
 use m_warning,only: issue_warning,msg
 use m_basis_set
 use m_spectral_function
 use m_eri,only: eri_3center_eigen
 implicit none

 integer,intent(in)                 :: gwmethod
 type(basis_set)                    :: basis,prod_basis
 real(dp),intent(in)                :: occupation(basis%nbf,nspin),energy(basis%nbf,nspin),exchange_m_vxc_diag(basis%nbf,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: selfenergy(basis%nbf,basis%nbf,nspin)
!=====
 logical               :: file_exists=.FALSE.
 logical               :: write_sigma_omega=.FALSE.
 integer               :: nomegai
 integer               :: iomegai
 real(dp),allocatable  :: omegai(:)
 real(dp),allocatable  :: selfenergy_tmp(:,:,:,:)
 integer               :: bbf,ibf,kbf
 integer               :: astate,bstate
 integer               :: istate,ispin,ipole
 real(dp)              :: overlap_tmp
 real(dp)              :: bra(wpol%npole_reso,basis%nbf) ! ,ket(wpol%npole,basis%nbf)
 real(dp)              :: fact_full,fact_empty
 real(dp)              :: zz(nspin)
 real(dp)              :: energy_qp(basis%nbf,nspin)
 real(dp)              :: energy_qp_new(basis%nbf,nspin)
 character(len=3)      :: ctmp
 integer               :: reading_status
!=====

 call start_clock(timing_self)


 WRITE_ME(msg,'(es9.2)') AIMAG(ieta)
 msg='small complex number is '//msg
 call issue_warning(msg)

 WRITE_MASTER(*,*)
 select case(gwmethod)
 case(QS)
   WRITE_MASTER(*,*) 'perform a QP self-consistent GW calculation'
 case(perturbative)
   WRITE_MASTER(*,*) 'perform a one-shot G0W0 calculation'
 case(COHSEX)
   WRITE_MASTER(*,*) 'perform a COHSEX calculation'
 case(QSCOHSEX)
   WRITE_MASTER(*,*) 'perform a self-consistent COHSEX calculation'
 case(GnW0)
   WRITE_MASTER(*,*) 'perform an eigenvalue self-consistent GnW0 calculation'
 case(GnWn)
   WRITE_MASTER(*,*) 'perform an eigenvalue self-consistent GnWn calculation'
 end select

 if(gwmethod==QS .OR. gwmethod==COHSEX .OR. gwmethod==QSCOHSEX .OR.  gwmethod==GnW0 .OR. gwmethod==GnWn) then
   nomegai=1
   allocate(omegai(nomegai))
 else
   ! look for manual omegas
   inquire(file='manual_omega',exist=file_exists)
   if(file_exists) then
     msg='reading frequency file for self-energy evaluation'
     call issue_warning(msg)
     open(12,file='manual_omega',status='old')
     read(12,*) nomegai
     allocate(omegai(nomegai))
     read(12,*) omegai(1)
     read(12,*) omegai(nomegai)
     close(12)
     do iomegai=2,nomegai-1
       omegai(iomegai) = omegai(1) + (omegai(nomegai)-omegai(1)) * (iomegai-1) / DBLE(nomegai-1)
     enddo
     write_sigma_omega=.TRUE.
   else
     nomegai=3
     allocate(omegai(nomegai))
     omegai(1) = -0.01_dp
     omegai(2) =  0.00_dp
     omegai(3) =  0.01_dp
   endif
 endif

 if( gwmethod==GnW0 .OR. gwmethod==GnWn ) then
   call read_energy_qp(nspin,basis%nbf,energy_qp,reading_status)
   if(reading_status/=0) then
     msg='File energy_qp not found: assuming 1st iteration'
     call issue_warning(msg)
     energy_qp(:,:) = energy(:,:)
   endif
 else
   energy_qp(:,:) = energy(:,:)
 endif

 allocate(selfenergy_tmp(nomegai,basis%nbf,basis%nbf,nspin))
 selfenergy_tmp(:,:,:,:) = 0.0_dp

 do ispin=1,nspin
   do istate=1,basis%nbf !INNER LOOP of G

     !
     ! Apply the frozen core and frozen virtual approximation to G
     if(istate <= ncore_G)    cycle
     if(istate >= nvirtual_G) cycle

     !
     ! Prepare the bra and ket with the knowledge of index istate and astate
     if( .NOT. has_auxil_basis) then
       ! Here just grab the precalculated value
       do astate=1,basis%nbf
         kbf = prod_basis%index_prodbasis(istate,astate) + prod_basis%nbf*(ispin-1)
         bra(:,astate) = wpol%residu_left(kbf,:)
!         ket(:,astate) = wpol%residu_right(kbf,:)
       enddo
     else
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,:) = MATMUL( TRANSPOSE(wpol%residu_left(:,:)) , eri_3center_eigen(:,:,istate,ispin) )
!       ket(:,:) = MATMUL( wpol%residu_right(:,:) , eri_3center_eigen(:,:,istate,ispin) )
     endif

     do ipole=1,wpol%npole_reso

       if( wpol%pole(ipole) < 0.0_dp ) then
         fact_empty = (spin_fact - occupation(istate,ispin)) / spin_fact
         fact_full = 0.0_dp
       else
         fact_empty = 0.0_dp
         fact_full = occupation(istate,ispin) / spin_fact
       endif
       if( ABS(fact_empty - fact_full) < 0.0001 ) cycle

       select case(gwmethod)
       case(QS)

         do bstate=1,basis%nbf
           do astate=1,basis%nbf

             selfenergy_tmp(1,astate,bstate,ispin) = selfenergy_tmp(1,astate,bstate,ispin) &
                        - bra(ipole,astate) * bra(ipole,bstate)                                      &    ! * ket(ipole,bstate) &
                          * ( REAL( fact_empty / ( energy_qp(bstate,ispin) + ieta  - energy_qp(istate,ispin) + wpol%pole(ipole) )        &
                                   -fact_full  / ( energy_qp(bstate,ispin) - ieta  - energy_qp(istate,ispin) + wpol%pole(ipole) ) , dp ) &
                             -REAL( fact_empty / ( energy_qp(bstate,ispin) + ieta  - energy_qp(istate,ispin) - wpol%pole(ipole) )        &
                                   -fact_full  / ( energy_qp(bstate,ispin) - ieta  - energy_qp(istate,ispin) - wpol%pole(ipole) ) , dp ) )

           enddo
         enddo

       case(perturbative,GnW0,GnWn)

         do astate=1,basis%nbf
           !
           ! calculate only the diagonal !
           bstate=astate
!           do bstate=1,basis%nbf
             do iomegai=1,nomegai
               selfenergy_tmp(iomegai,astate,bstate,ispin) = selfenergy_tmp(iomegai,astate,bstate,ispin) &
                        - bra(ipole,astate) * bra(ipole,bstate) & ! * ket(ipole,bstate) &
                          * ( REAL(  fact_empty / ( energy_qp(bstate,ispin) + ieta + omegai(iomegai) - energy_qp(istate,ispin) + wpol%pole(ipole)     )          &
                                    -fact_full  / ( energy_qp(bstate,ispin) - ieta + omegai(iomegai) - energy_qp(istate,ispin) + wpol%pole(ipole)     )  , dp )  &
                             -REAL(  fact_empty / ( energy_qp(bstate,ispin) + ieta + omegai(iomegai) - energy_qp(istate,ispin) - wpol%pole(ipole)     )          &
                                    -fact_full  / ( energy_qp(bstate,ispin) - ieta + omegai(iomegai) - energy_qp(istate,ispin) - wpol%pole(ipole)     )  , dp ) )
             enddo
!           enddo
         enddo

       case(COHSEX,QSCOHSEX) 

         !
         ! SEX
         !
         do bstate=1,basis%nbf
           do astate=1,basis%nbf
             selfenergy_tmp(1,astate,bstate,ispin) = selfenergy_tmp(1,astate,bstate,ispin) &
                        - bra(ipole,astate) * bra(ipole,bstate) &   ! * ket(ipole,bstate) &
                          * fact_full  / (-wpol%pole(ipole)) * 2.0_dp
           enddo
         enddo
         !
         ! COH
         !
         do bstate=1,basis%nbf
           do astate=1,basis%nbf
             selfenergy_tmp(1,astate,bstate,ispin) = selfenergy_tmp(1,astate,bstate,ispin) &
                        + bra(ipole,astate) * bra(ipole,bstate) &   !  * ket(ipole,bstate) &
                          * fact_empty  / (-wpol%pole(ipole)) * 2.0_dp
           enddo
         enddo

       case default 
         stop'BUG'
       end select

     enddo !ipole

   enddo !istate
 enddo !ispin

 !
 ! Kotani's Hermitianization trick
 !
 if(gwmethod==QS) then
   do ispin=1,nspin
     selfenergy_tmp(1,:,:,ispin) = 0.5_dp * ( selfenergy_tmp(1,:,:,ispin) + transpose(selfenergy_tmp(1,:,:,ispin)) )
   enddo
 endif

 select case(gwmethod)
 case(QS)
   ! Transform the matrix elements back to the non interacting states
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,ispin) ) , MATMUL( selfenergy_tmp(1,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,ispin)), s_matrix(:,:) ) ) )
   enddo

 case(QSCOHSEX)
   ! Transform the matrix elements back to the non interacting states
   ! do not forget the overlap matrix S
   ! C^T S C = I
   ! the inverse of C is C^T S
   ! the inverse of C^T is S C
   do ispin=1,nspin
     selfenergy(:,:,ispin) = MATMUL( MATMUL( s_matrix(:,:) , c_matrix(:,:,ispin)) , MATMUL( selfenergy_tmp(1,:,:,ispin), &
                             MATMUL( TRANSPOSE(c_matrix(:,:,ispin)), s_matrix(:,:) ) ) )
   enddo


 case(COHSEX) !==========================================================

   selfenergy(:,:,:) = REAL( selfenergy_tmp(1,:,:,:) )
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,*) 'COHSEX Eigenvalues [eV]'
   if(nspin==1) then
     WRITE_MASTER(*,*) '  #          E0        Sigx-Vxc      Sigc          Z         G0W0'
   else
     WRITE_MASTER(*,'(a)') '  #                E0                      Sigx-Vxc                    Sigc                       Z                       G0W0'
   endif
   do astate=1,basis%nbf
     zz(:) = 1.0_dp 
     energy_qp_new(astate,:) = energy_qp(astate,:)+zz(:)*REAL(selfenergy_tmp(1,astate,astate,:) + exchange_m_vxc_diag(astate,:))

     WRITE_MASTER(*,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,exchange_m_vxc_diag(astate,:)*Ha_eV,REAL(selfenergy_tmp(1,astate,astate,:),dp)*Ha_eV,&
           zz(:),energy_qp_new(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp_new)

 case(perturbative) !==========================================================

   if(write_sigma_omega) then

     do astate=1,basis%nbf
       WRITE_ME(ctmp,'(i3.3)') astate
       open(200+astate,file='selfenergy_omega_state'//TRIM(ctmp))
       do iomegai=1,nomegai
         WRITE_MASTER(200+astate,'(20(f12.6,2x))') ( DBLE(omegai(iomegai))+energy_qp(astate,:) )*Ha_eV,&
                                                     ( DBLE(selfenergy_tmp(iomegai,astate,astate,:)) )*Ha_eV,&
                                                     ( DBLE(omegai(iomegai))-exchange_m_vxc_diag(astate,:) )*Ha_eV,&
                                                     ( 1.0_dp/pi/ABS( DBLE(omegai(iomegai))-exchange_m_vxc_diag(astate,:) &
                                                                   - DBLE(selfenergy_tmp(iomegai,astate,astate,:)) ) ) / Ha_eV
       enddo
       WRITE_MASTER(200+astate,*)
     enddo
     close(200+astate)

   else
     selfenergy(:,:,:) = REAL( selfenergy_tmp(2,:,:,:) )
     WRITE_MASTER(*,*)
     WRITE_MASTER(*,*) 'G0W0 Eigenvalues [eV]'
     if(nspin==1) then
       WRITE_MASTER(*,*) '  #          E0        Sigx-Vxc      Sigc          Z         G0W0'
     else
       WRITE_MASTER(*,'(a)') '  #                E0                      Sigx-Vxc                    Sigc                       Z                       G0W0'
     endif
     do astate=1,basis%nbf
       zz(:) = REAL( selfenergy_tmp(3,astate,astate,:) - selfenergy_tmp(1,astate,astate,:) ) / REAL( omegai(3)-omegai(1) )
       zz(:) = 1.0_dp / ( 1.0_dp - zz(:) )
       ! Contrain Z to be in [0:1] to avoid crazy values
       do ispin=1,nspin
         zz(ispin) = MIN( MAX(zz(ispin),0.0_dp) , 1.0_dp )
       enddo
       energy_qp_new(astate,:) = energy_qp(astate,:)+zz(:)*REAL(selfenergy_tmp(2,astate,astate,:) + exchange_m_vxc_diag(astate,:))

       WRITE_MASTER(*,'(i4,x,20(x,f12.6))') astate,energy_qp(astate,:)*Ha_eV,exchange_m_vxc_diag(astate,:)*Ha_eV,REAL(selfenergy_tmp(2,astate,astate,:),dp)*Ha_eV,&
             zz(:),energy_qp_new(astate,:)*Ha_eV
     enddo

     call write_energy_qp(nspin,basis%nbf,energy_qp_new)
   endif

 case(GnW0,GnWn) !==========================================================

     selfenergy(:,:,:) = REAL( selfenergy_tmp(1,:,:,:) )
     WRITE_MASTER(*,*)
     WRITE_MASTER(*,*) 'G0W0 Eigenvalues [eV]'
     if(nspin==1) then
       WRITE_MASTER(*,'(a)') '  #          E0        Sigx-Vxc      Sigc          Z          GW(n-1)       GW(n)'
     else
       WRITE_MASTER(*,'(a)') '  #                E0                      Sigx-Vxc                    Sigc                       Z                       G0W0'
     endif
     do astate=1,basis%nbf
       zz(:) = 1.0_dp 

       energy_qp_new(astate,:) = energy(astate,:)+zz(:)*REAL(selfenergy_tmp(1,astate,astate,:) + exchange_m_vxc_diag(astate,:))

       WRITE_MASTER(*,'(i4,x,20(x,f12.6))') astate,energy(astate,:)*Ha_eV,exchange_m_vxc_diag(astate,:)*Ha_eV,REAL(selfenergy_tmp(2,astate,astate,:),dp)*Ha_eV,&
             zz(:),energy_qp(astate,:)*Ha_eV,energy_qp_new(astate,:)*Ha_eV
     enddo

     call write_energy_qp(nspin,basis%nbf,energy_qp_new)

 end select

 deallocate(omegai)
 deallocate(selfenergy_tmp)

 call stop_clock(timing_self)

end subroutine gw_selfenergy


!=========================================================================
