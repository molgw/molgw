!=========================================================================
#include "macros.h"
!=========================================================================
module m_gw
 use m_definitions
 use m_mpi
 use m_timing 
 use m_inputparam


contains


!=========================================================================
subroutine chi_to_vchiv(nbf,prod_basis,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf
 type(basis_set),intent(in)            :: prod_basis
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 real(dp),intent(in)                   :: c_matrix(nbf,nbf,nspin)
 real(dp),allocatable,intent(inout)    :: eigenvector    (:,:)
 real(dp),allocatable,intent(inout)    :: eigenvector_inv(:,:)
 real(dp),allocatable,intent(inout)    :: eigenvalue     (:)
 type(spectral_function),intent(inout) :: wpol
!=====
 integer                               :: t_kl,klspin,ijspin
 integer                               :: istate,jstate,kstate,lstate,ijstate,ijstate_spin
 integer                               :: klstate_min
 integer                               :: klstate_max
 integer                               :: ipole
 real(dp)                              :: docc_kl
 real(dp)                              :: eri_eigen_klij
 real(dp),allocatable                  :: eri_eigenstate_klmin(:,:,:,:)
!=====

 call start_clock(timing_buildw)

 WRITE_MASTER(*,'(/,a)') ' Build W = v * chi * v'

 if( .NOT. is_auxil_basis ) then
   allocate(eri_eigenstate_klmin(nbf,nbf,nbf,nspin))
   ! Set this to zero and then enforce the calculation of the first array of Coulomb integrals
   eri_eigenstate_klmin(:,:,:,:) = 0.0_dp
 endif

 call allocate_spectral_function(prod_basis%nbf*nspin,wpol)

 wpol%pole(:) = eigenvalue(:)

 wpol%residu_left (:,:) = 0.0_dp
 wpol%residu_right(:,:) = 0.0_dp
 do t_kl=1,wpol%npole
   kstate = wpol%transition_table(1,t_kl)
   lstate = wpol%transition_table(2,t_kl)
   klspin = wpol%transition_table(3,t_kl)
   docc_kl= occupation(kstate,klspin)-occupation(lstate,klspin)

   if( .NOT. is_auxil_basis ) then
     klstate_min = MIN(kstate,lstate)
     klstate_max = MAX(kstate,lstate)
     call transform_eri_basis(nspin,c_matrix,klstate_min,klspin,eri_eigenstate_klmin)
   endif

   do ijspin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(istate,jstate,ijstate_spin,eri_eigen_klij)
     do ijstate=1,prod_basis%nbf
       istate = prod_basis%index_ij(1,ijstate)
       jstate = prod_basis%index_ij(2,ijstate)

       ijstate_spin = ijstate+prod_basis%nbf*(ijspin-1)

       if(is_auxil_basis) then
         eri_eigen_klij = eri_eigen_ri(kstate,lstate,klspin,istate,jstate,ijspin)
       else
         eri_eigen_klij = eri_eigenstate_klmin(klstate_max,istate,jstate,ijspin)
       endif

       wpol%residu_left (:,ijstate_spin)  = wpol%residu_left (:,ijstate_spin) &
                    + eri_eigen_klij * eigenvector(t_kl,:)
       wpol%residu_right(:,ijstate_spin)  = wpol%residu_right(:,ijstate_spin) &
                    + eri_eigen_klij * eigenvector_inv(:,t_kl) * docc_kl

     enddo
!$OMP END DO
!$OMP END PARALLEL
   enddo
 enddo

 deallocate(eigenvector)
 deallocate(eigenvector_inv)
 deallocate(eigenvalue)
 if(allocated(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)

 call stop_clock(timing_buildw)

end subroutine chi_to_vchiv


!=========================================================================
subroutine chi_to_sqrtvchisqrtv_auxil(nbf,nbf_auxil,occupation,c_matrix,eigenvector,eigenvector_inv,eigenvalue,wpol)
 use m_definitions
 use m_warning
 use m_basis_set
 use m_eri
 use m_spectral_function
 implicit none
 
 integer,intent(in)                    :: nbf,nbf_auxil
 real(dp),intent(in)                   :: occupation(nbf,nspin)
 real(dp),intent(in)                   :: c_matrix(nbf,nbf,nspin)
 real(dp),allocatable,intent(inout)    :: eigenvector    (:,:)
 real(dp),allocatable,intent(inout)    :: eigenvector_inv(:,:)
 real(dp),allocatable,intent(inout)    :: eigenvalue     (:)
 type(spectral_function),intent(inout) :: wpol
!=====
 integer                               :: t_ij,t_kl,klspin,ijspin
 integer                               :: kstate,lstate
 real(dp)                              :: docc_kl
 real(dp),allocatable                  :: eri_3center_2index(:,:)
!=====

 call start_clock(timing_buildw)

 WRITE_MASTER(*,'(/,a)') ' Build v^{1/2} * chi * v^{1/2}'

 call allocate_spectral_function(nbf_auxil,wpol)
 wpol%pole(:) = eigenvalue(:)
 deallocate(eigenvalue)

 allocate(eri_3center_2index(nbf_auxil,wpol%npole))
 do t_kl=1,wpol%npole
   kstate = wpol%transition_table(1,t_kl)
   lstate = wpol%transition_table(2,t_kl)
   klspin = wpol%transition_table(3,t_kl)
   eri_3center_2index(:,t_kl) = eri_3center_eigen(:,kstate,lstate,klspin)
 enddo
 wpol%residu_left(:,:)  = TRANSPOSE( MATMUL( eri_3center_2index , eigenvector ) )

 do t_kl=1,wpol%npole
   kstate = wpol%transition_table(1,t_kl)
   lstate = wpol%transition_table(2,t_kl)
   klspin = wpol%transition_table(3,t_kl)
   docc_kl = occupation(kstate,klspin)-occupation(lstate,klspin)
   eri_3center_2index(:,t_kl) = eri_3center_2index(:,t_kl) * docc_kl
 enddo
 wpol%residu_right(:,:) = MATMUL( eigenvector_inv , TRANSPOSE( eri_3center_2index ) )

 deallocate(eri_3center_2index)


 deallocate(eigenvector)
 deallocate(eigenvector_inv)

 call stop_clock(timing_buildw)

end subroutine chi_to_sqrtvchisqrtv_auxil


!=========================================================================
subroutine gw_selfenergy(gwmethod,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,selfenergy)
 use m_definitions
 use m_mpi
 use m_timing 
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
 real(dp)              :: bra(wpol%npole,basis%nbf),ket(wpol%npole,basis%nbf)
 real(dp)              :: fact_full,fact_empty
 real(dp)              :: zz(nspin)
 real(dp)              :: energy_re
 real(dp)              :: energy_qp(basis%nbf,nspin)
 character(len=3)      :: ctmp
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
 end select

 if(gwmethod==QS .OR. gwmethod==COHSEX .OR. gwmethod==QSCOHSEX) then
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
     if( .NOT. is_auxil_basis) then
       ! Here just grab the precalculated value
       do astate=1,basis%nbf
         kbf = prod_basis%index_prodbasis(istate,astate) + prod_basis%nbf*(ispin-1)
         bra(:,astate) = wpol%residu_left (:,kbf)
         ket(:,astate) = wpol%residu_right(:,kbf)
       enddo
     else
       ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
       bra(:,:) = MATMUL( wpol%residu_left (:,:) , eri_3center_eigen(:,:,istate,ispin) )
       ket(:,:) = MATMUL( wpol%residu_right(:,:) , eri_3center_eigen(:,:,istate,ispin) )
     endif

     do ipole=1,wpol%npole

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
                        - bra(ipole,astate) * ket(ipole,bstate) &
                          * REAL(   fact_empty / ( energy(bstate,ispin) + ieta  - energy(istate,ispin) + wpol%pole(ipole) ) &
                                   -fact_full  / ( energy(bstate,ispin) - ieta  - energy(istate,ispin) + wpol%pole(ipole) ) , dp )

           enddo
         enddo

       case(perturbative)

         do astate=1,basis%nbf
           !
           ! calculate only the diagonal !
           bstate=astate
!           do bstate=1,basis%nbf
             do iomegai=1,nomegai
               selfenergy_tmp(iomegai,astate,bstate,ispin) = selfenergy_tmp(iomegai,astate,bstate,ispin) &
                        - bra(ipole,astate) * ket(ipole,bstate) &
                          * REAL(  fact_empty / ( energy(bstate,ispin) + ieta + omegai(iomegai) - energy(istate,ispin) + wpol%pole(ipole)     ) &
                                  -fact_full  / ( energy(bstate,ispin) - ieta + omegai(iomegai) - energy(istate,ispin) + wpol%pole(ipole)     )  , dp )
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
                        - bra(ipole,astate) * ket(ipole,bstate) &
                          * fact_full  / (-wpol%pole(ipole))    
           enddo
         enddo
         !
         ! COH
         !
         do bstate=1,basis%nbf
           do astate=1,basis%nbf
             selfenergy_tmp(1,astate,bstate,ispin) = selfenergy_tmp(1,astate,bstate,ispin) &
                        + bra(ipole,astate) * ket(ipole,bstate) &
                          * fact_empty  / (-wpol%pole(ipole))
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
     energy_qp(astate,:) = energy(astate,:)+zz(:)*REAL(selfenergy_tmp(1,astate,astate,:) + exchange_m_vxc_diag(astate,:))

     WRITE_MASTER(*,'(i4,x,20(x,f12.6))') astate,energy(astate,:)*Ha_eV,exchange_m_vxc_diag(astate,:)*Ha_eV,REAL(selfenergy_tmp(1,astate,astate,:),dp)*Ha_eV,&
           zz(:),energy_qp(astate,:)*Ha_eV
   enddo

   call write_energy_qp(nspin,basis%nbf,energy_qp)

 case(perturbative) !==========================================================

   if(write_sigma_omega) then

     do astate=1,basis%nbf
       WRITE_ME(ctmp,'(i3.3)') astate
       open(200+astate,file='selfenergy_omega_state'//TRIM(ctmp))
       do iomegai=1,nomegai
         WRITE_MASTER(200+astate,'(20(f12.6,2x))') ( DBLE(omegai(iomegai))+energy(astate,:) )*Ha_eV,&
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
       energy_qp(astate,:) = energy(astate,:)+zz(:)*REAL(selfenergy_tmp(2,astate,astate,:) + exchange_m_vxc_diag(astate,:))

       WRITE_MASTER(*,'(i4,x,20(x,f12.6))') astate,energy(astate,:)*Ha_eV,exchange_m_vxc_diag(astate,:)*Ha_eV,REAL(selfenergy_tmp(2,astate,astate,:),dp)*Ha_eV,&
             zz(:),energy_qp(astate,:)*Ha_eV
     enddo

     call write_energy_qp(nspin,basis%nbf,energy_qp)
   endif

 end select

 deallocate(omegai)
 deallocate(selfenergy_tmp)

 call stop_clock(timing_self)

end subroutine gw_selfenergy


!=========================================================================
end module m_gw
!=========================================================================
