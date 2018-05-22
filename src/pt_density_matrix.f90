!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the many-body perturbation theory to obtain the (perturbative) density matrix
!
!=========================================================================
subroutine pt1_density_matrix(nstate,basis,occupation,energy,c_matrix,exchange_m_vxc,p_matrix)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_warning
 use m_timing
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_hamiltonian
 use m_hamiltonian_onebody
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: exchange_m_vxc(nstate,nstate,nspin)
 real(dp),intent(out)       :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer                 :: pstate,istate,astate
 integer                 :: iaspin
 real(dp)                :: denom
 real(dp)                :: p_matrix_pt1(nstate,nstate)
!=====

 write(stdout,'(/,a)') ' Calculate the PT1 density matrix'

 if( nspin /= 1 ) call die('pt1_density_matrix: only implemented for spin restricted calculations')

 iaspin = 1
 p_matrix_pt1(:,:) = 0.0_dp
 do istate=ncore_G+1,nhomo_G
   do astate=nhomo_G+1,nvirtual_G-1

     denom = energy(istate,iaspin) - energy(astate,iaspin)
     p_matrix_pt1(istate,astate) = p_matrix_pt1(istate,astate) + exchange_m_vxc(istate,astate,iaspin) / denom
     p_matrix_pt1(astate,istate) = p_matrix_pt1(astate,istate) + exchange_m_vxc(istate,astate,iaspin) / denom

   enddo
 enddo

! ! Add the SCF density matrix to get to the total density matrix
! do pstate=ncore_G+1,nvirtual_G-1
!   p_matrix_pt1(pstate,pstate) = p_matrix_pt1(pstate,pstate) + occupation(pstate,iaspin)
! enddo

 p_matrix(:,:,iaspin) = MATMUL( c_matrix(:,:,iaspin)  , &
                          MATMUL( p_matrix_pt1(:,:), &
                             TRANSPOSE(c_matrix(:,:,iaspin)) ) )


end subroutine pt1_density_matrix


!=========================================================================
subroutine pt2_density_matrix(nstate,basis,occupation,energy,c_matrix,p_matrix)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_warning
 use m_timing
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_hamiltonian
 use m_hamiltonian_onebody
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)       :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer                 :: pstate,qstate
 integer                 :: istate,jstate,kstate,lstate
 integer                 :: astate,bstate,cstate,dstate
 integer                 :: pqspin
 real(dp)                :: denom1,denom2
 real(dp)                :: num1,num2
 real(dp)                :: p_matrix_pt2(nstate,nstate)
!=====



 write(stdout,'(/,a)') ' Calculate the PT2 density matrix'

 if( nspin /= 1 ) call die('pt2_density_matrix: only implemented for spin restricted calculations')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif


! Full calculation of the MP2 density matrix

 p_matrix_pt2(:,:) = 0.0_dp
 pqspin = 1

 ! A1
 do istate=ncore_G+1,nhomo_G
 do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
       do bstate=nhomo_G+1,nvirtual_G-1

         denom1 = energy(astate,pqspin) + energy(bstate,pqspin) - energy(istate,pqspin) - energy(kstate,pqspin)
         denom2 = energy(astate,pqspin) + energy(bstate,pqspin) - energy(jstate,pqspin) - energy(kstate,pqspin)

         num1 = 2.0_dp * eri_eigen(istate,astate,pqspin,kstate,bstate,pqspin) &
                      - eri_eigen(istate,bstate,pqspin,kstate,astate,pqspin)
         num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,kstate,bstate,pqspin)

         p_matrix_pt2(istate,jstate) = p_matrix_pt2(istate,jstate)  &
                 - num1 * num2 / ( denom1 * denom2 )

       enddo
     enddo
   enddo
 enddo
 enddo

 ! A2
 do astate=nhomo_G+1,nvirtual_G-1
 do bstate=nhomo_G+1,nvirtual_G-1
   do cstate=nhomo_G+1,nvirtual_G-1
     do istate=ncore_G+1,nhomo_G
       do jstate=ncore_G+1,nhomo_G

         denom1 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(astate,pqspin) - energy(cstate,pqspin)
         denom2 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(bstate,pqspin) - energy(cstate,pqspin)

         num1 = 2.0_dp * eri_eigen(astate,istate,pqspin,cstate,jstate,pqspin) &
                      - eri_eigen(astate,jstate,pqspin,cstate,istate,pqspin)
         num2 = 2.0_dp * eri_eigen(bstate,istate,pqspin,cstate,jstate,pqspin)

         p_matrix_pt2(astate,bstate) = p_matrix_pt2(astate,bstate)  &
                 + num1 * num2 / ( denom1 * denom2 )

       enddo
     enddo
   enddo
 enddo
 enddo

 ! A3    P_cj  sum over i,a,b
 ! A4    P_jc  sum over i,a,b
 do cstate=nhomo_G+1,nvirtual_G-1
 do jstate=ncore_G+1,nhomo_G
   do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       do istate=ncore_G+1,nhomo_G
         denom1 = energy(jstate,pqspin) + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         denom2 = energy(jstate,pqspin) - energy(cstate,pqspin)
         num1 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
         num2 = 2.0_dp * eri_eigen(astate,cstate,pqspin,bstate,istate,pqspin)

         p_matrix_pt2(cstate,jstate) = p_matrix_pt2(cstate,jstate)  &
                                           + num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(jstate,cstate) = p_matrix_pt2(jstate,cstate)  &
                                           + num1 * num2 / ( denom1 * denom2 )
       enddo
     enddo
   enddo
 enddo
 enddo

 ! A5   P_bk  sum over i,j,a
 ! A6   P_kb  sum over i,j,a
 do bstate=nhomo_G+1,nvirtual_G-1
 do kstate=ncore_G+1,nhomo_G
   do astate=nhomo_G+1,nvirtual_G-1
     do istate=ncore_G+1,nhomo_G
       do jstate=ncore_G+1,nhomo_G
         denom1 = energy(jstate,pqspin) + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         denom2 = energy(kstate,pqspin) - energy(bstate,pqspin)
         num1 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
         num2 = 2.0_dp * eri_eigen(istate,kstate,pqspin,jstate,astate,pqspin)

         p_matrix_pt2(bstate,kstate) = p_matrix_pt2(bstate,kstate)  &
                                           - num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(kstate,bstate) = p_matrix_pt2(kstate,bstate)  &
                                           - num1 * num2 / ( denom1 * denom2 )
       enddo
     enddo
   enddo
 enddo
 enddo

! open(111,file='p_matrix_pt2.dat',action='write')
! do pstate=1,nstate
!   write(111,'(*(2x,f12.6))') p_matrix_pt2(pstate,:)
! enddo
! close(111)

 ! Add the SCF density matrix to get to the total density matrix
 do pstate=1,nstate
   p_matrix_pt2(pstate,pstate) = p_matrix_pt2(pstate,pstate) + occupation(pstate,pqspin)
 enddo

 p_matrix(:,:,pqspin) = MATMUL( c_matrix(:,:,pqspin)  , &
                          MATMUL( p_matrix_pt2(:,:), &
                             TRANSPOSE(c_matrix(:,:,pqspin)) ) )

!block
! real(dp) :: hh(basis%nbf,basis%nbf)
! real(dp) :: hartree_ii(nstate)
! call calculate_hartree(basis,p_matrix,hh)
! do istate=1,nstate
!   write(stdout,'(1x,a,i5,2x,f12.6)') 'Occ Occ Hartree ii ',istate,DOT_PRODUCT( c_matrix(:,istate,pqspin), MATMUL( hh(:,:), c_matrix(:,istate,pqspin) ) ) * Ha_eV
! enddo
!end block


 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 else
   call destroy_eri_4center_eigen_uks()
 endif


end subroutine pt2_density_matrix


!=========================================================================
subroutine onering_density_matrix(nstate,basis,occupation,energy,c_matrix,p_matrix)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_warning
 use m_timing
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_hamiltonian
 use m_hamiltonian_onebody
 use m_selfenergy_tools
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(out)       :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer                 :: pstate,qstate
 integer                 :: istate,jstate,kstate
 integer                 :: astate,bstate,cstate,dstate
 integer                 :: pqspin
 real(dp)                :: denom1,denom2
 real(dp)                :: num1,num2
 real(dp)                :: p_matrix_pt2(nstate,nstate)
!=====



 write(stdout,'(/,a)') ' Calculate the PT2 density matrix'

 if( nspin /= 1 ) call die('pt2_density_matrix: only implemented for spin restricted calculations')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif


! Full calculation of the MP2 density matrix

 p_matrix_pt2(:,:) = 0.0_dp
 pqspin = 1

 ! A1
 do istate=ncore_G+1,nhomo_G
 do jstate=ncore_G+1,nhomo_G
   do kstate=ncore_G+1,nhomo_G
     do astate=nhomo_G+1,nvirtual_G-1
       do bstate=nhomo_G+1,nvirtual_G-1

         denom1 = energy(astate,pqspin) + energy(bstate,pqspin) - energy(istate,pqspin) - energy(kstate,pqspin)
         denom2 = energy(astate,pqspin) + energy(bstate,pqspin) - energy(jstate,pqspin) - energy(kstate,pqspin)

         num1 = 2.0_dp * eri_eigen(istate,astate,pqspin,kstate,bstate,pqspin)
         num2 = 2.0_dp * eri_eigen(jstate,astate,pqspin,kstate,bstate,pqspin)

         p_matrix_pt2(istate,jstate) = p_matrix_pt2(istate,jstate)  &
                 - num1 * num2 / ( denom1 * denom2 )

       enddo
     enddo
   enddo
 enddo
 enddo

 ! A2
 do astate=nhomo_G+1,nvirtual_G-1
 do bstate=nhomo_G+1,nvirtual_G-1
   do cstate=nhomo_G+1,nvirtual_G-1
     do istate=ncore_G+1,nhomo_G
       do jstate=ncore_G+1,nhomo_G

         denom1 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(astate,pqspin) - energy(cstate,pqspin)
         denom2 = energy(istate,pqspin) + energy(jstate,pqspin) - energy(bstate,pqspin) - energy(cstate,pqspin)

         num1 = 2.0_dp * eri_eigen(astate,istate,pqspin,cstate,jstate,pqspin)
         num2 = 2.0_dp * eri_eigen(bstate,istate,pqspin,cstate,jstate,pqspin)

         p_matrix_pt2(astate,bstate) = p_matrix_pt2(astate,bstate)  &
                 + num1 * num2 / ( denom1 * denom2 )

       enddo
     enddo
   enddo
 enddo
 enddo

 ! A3    P_cj  sum over i,a,b
 ! A4    P_jc  sum over i,a,b
 do cstate=nhomo_G+1,nvirtual_G-1
 do jstate=ncore_G+1,nhomo_G
   do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       do istate=ncore_G+1,nhomo_G
         denom1 = energy(jstate,pqspin) + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         denom2 = energy(jstate,pqspin) - energy(cstate,pqspin)
         num1 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin)
         num2 = 2.0_dp * eri_eigen(astate,cstate,pqspin,bstate,istate,pqspin)

         p_matrix_pt2(cstate,jstate) = p_matrix_pt2(cstate,jstate)  &
                                           + num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(jstate,cstate) = p_matrix_pt2(jstate,cstate)  &
                                           + num1 * num2 / ( denom1 * denom2 )
       enddo
     enddo
   enddo
 enddo
 enddo

 ! A5   P_bk  sum over i,j,a
 ! A6   P_kb  sum over i,j,a
 do bstate=nhomo_G+1,nvirtual_G-1
 do kstate=ncore_G+1,nhomo_G
   do astate=nhomo_G+1,nvirtual_G-1
     do istate=ncore_G+1,nhomo_G
       do jstate=ncore_G+1,nhomo_G
         denom1 = energy(jstate,pqspin) + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         denom2 = energy(kstate,pqspin) - energy(bstate,pqspin)
         num1 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin)
         num2 = 2.0_dp * eri_eigen(istate,kstate,pqspin,jstate,astate,pqspin)

         p_matrix_pt2(bstate,kstate) = p_matrix_pt2(bstate,kstate)  &
                                           - num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(kstate,bstate) = p_matrix_pt2(kstate,bstate)  &
                                           - num1 * num2 / ( denom1 * denom2 )
       enddo
     enddo
   enddo
 enddo
 enddo

! open(111,file='p_matrix_pt2.dat',action='write')
! do pstate=1,nstate
!   write(111,'(*(2x,f12.6))') p_matrix_pt2(pstate,:)
! enddo
! close(111)

 ! Add the SCF density matrix to get to the total density matrix
 do pstate=1,nstate
   p_matrix_pt2(pstate,pstate) = p_matrix_pt2(pstate,pstate) + occupation(pstate,pqspin)
 enddo

 p_matrix(:,:,pqspin) = MATMUL( c_matrix(:,:,pqspin)  , &
                          MATMUL( p_matrix_pt2(:,:), &
                             TRANSPOSE(c_matrix(:,:,pqspin)) ) )

!block
! real(dp) :: hh(basis%nbf,basis%nbf)
! real(dp) :: hartree_ii(nstate)
! call calculate_hartree(basis,p_matrix,hh)
! do istate=1,nstate
!   write(stdout,'(1x,a,i5,2x,f12.6)') 'Occ Occ Hartree ii ',istate,DOT_PRODUCT( c_matrix(:,istate,pqspin), MATMUL( hh(:,:), c_matrix(:,istate,pqspin) ) ) * Ha_eV
! enddo
!end block


 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 else
   call destroy_eri_4center_eigen_uks()
 endif


end subroutine onering_density_matrix


!=========================================================================
subroutine gw_density_matrix(nstate,basis,occupation,energy,c_matrix,wpol,p_matrix)
 use m_definitions
 use m_mpi
 use m_mpi_ortho
 use m_warning
 use m_timing
 use m_basis_set
 use m_eri_ao_mo
 use m_inputparam
 use m_hamiltonian
 use m_hamiltonian_onebody
 use m_selfenergy_tools
 use m_spectral_function
 implicit none

 integer,intent(in)                 :: nstate
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(nstate,nspin),energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer  :: pstate
 integer  :: istate,jstate,kstate,lstate
 integer  :: astate,bstate,cstate,dstate
 integer  :: pqspin
 real(dp) :: denom1,denom2
 real(dp) :: num1,num2
 real(dp) :: p_matrix_gw(nstate,nstate)
 integer  :: ipole
 real(dp),allocatable :: bra_occ(:,:),bra_virt(:,:)
!=====



 write(stdout,'(/,a)') ' Calculate the GW density matrix'

 if( nspin /= 1 ) call die('gw_density_matrix: only implemented for spin restricted calculations')
 if( .NOT. has_auxil_basis)  call die('gw_density_matrix: only implemented without RI')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif


! Full calculation of the MP2 density matrix

 p_matrix_gw(:,:) = 0.0_dp
 pqspin = 1

 allocate(bra_occ(wpol%npole_reso,ncore_G+1:nhomo_G))
 allocate(bra_virt(wpol%npole_reso,nhomo_G+1:nvirtual_G-1))

 ! A1
 do astate=nhomo_G+1,nvirtual_G-1

   bra_occ(:,ncore_G+1:nhomo_G) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,ncore_G+1:nhomo_G,astate,pqspin) )

   do kstate=ncore_G+1,nhomo_G
   do jstate=ncore_G+1,nhomo_G
     do ipole=1,wpol%npole_reso

         denom1 = energy(kstate,pqspin) - energy(astate,pqspin) - wpol%pole(ipole)
         denom2 = energy(jstate,pqspin) - energy(astate,pqspin) - wpol%pole(ipole)

         num1 = 2.0_dp * bra_occ(ipole,jstate)
         num2 = bra_occ(ipole,kstate)

         p_matrix_gw(kstate,jstate) = p_matrix_gw(kstate,jstate) - num1 * num2 / ( denom1 * denom2 )

     enddo
   enddo
   enddo
 enddo

 ! A2
 do istate=ncore_G+1,nhomo_G

   bra_virt(:,nhomo_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,nhomo_G+1:nvirtual_G-1,istate,pqspin) )

   do bstate=nhomo_G+1,nvirtual_G-1
   do cstate=nhomo_G+1,nvirtual_G-1
     do ipole=1,wpol%npole_reso

         denom1 = energy(istate,pqspin) - energy(cstate,pqspin) - wpol%pole(ipole)
         denom2 = energy(istate,pqspin) - energy(bstate,pqspin) - wpol%pole(ipole)

         num1 = 2.0_dp * bra_virt(ipole,cstate)
         num2 = bra_virt(ipole,bstate)

         p_matrix_gw(cstate,bstate) = p_matrix_gw(cstate,bstate) + num1 * num2 / ( denom1 * denom2 )

     enddo
   enddo
   enddo
 enddo

 ! A3    P_cj  sum over i,a,b
 ! A4    P_jc  sum over i,a,b
 do astate=nhomo_G+1,nvirtual_G-1

   bra_occ(:,ncore_G+1:nhomo_G)       = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,ncore_G+1:nhomo_G,astate,pqspin) )
   bra_virt(:,nhomo_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,nhomo_G+1:nvirtual_G-1,astate,pqspin) )

   do cstate=nhomo_G+1,nvirtual_G-1
   do jstate=ncore_G+1,nhomo_G
     do ipole=1,wpol%npole_reso

         denom1 = energy(jstate,pqspin) - energy(astate,pqspin) - wpol%pole(ipole)
         denom2 = energy(jstate,pqspin) - energy(cstate,pqspin)
         num1 = 2.0_dp * bra_virt(ipole,cstate)
         num2 = bra_occ(ipole,jstate)

         p_matrix_gw(cstate,jstate) = p_matrix_gw(cstate,jstate) + num1 * num2 / ( denom1 * denom2 )
         p_matrix_gw(jstate,cstate) = p_matrix_gw(jstate,cstate) + num1 * num2 / ( denom1 * denom2 )

     enddo
   enddo
   enddo
 enddo

 ! A5   P_bk  sum over i,j,a
 ! A6   P_kb  sum over i,j,a
 do istate=ncore_G+1,nhomo_G

   bra_occ(:,ncore_G+1:nhomo_G)       = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,ncore_G+1:nhomo_G,istate,pqspin) )
   bra_virt(:,nhomo_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,nhomo_G+1:nvirtual_G-1,istate,pqspin) )

   do bstate=nhomo_G+1,nvirtual_G-1
   do kstate=ncore_G+1,nhomo_G
     do ipole=1,wpol%npole_reso

         denom1 = energy(istate,pqspin) - energy(bstate,pqspin) - wpol%pole(ipole)
         denom2 = energy(kstate,pqspin) - energy(bstate,pqspin)
         num1 = 2.0_dp * bra_virt(ipole,bstate)
         num2 = bra_occ(ipole,kstate)

         p_matrix_gw(bstate,kstate) = p_matrix_gw(bstate,kstate) - num1 * num2 / ( denom1 * denom2 )
         p_matrix_gw(kstate,bstate) = p_matrix_gw(kstate,bstate) - num1 * num2 / ( denom1 * denom2 )

     enddo
   enddo
   enddo
 enddo

 deallocate(bra_occ)
 deallocate(bra_virt)

! open(111,file='p_matrix_pt2.dat',action='write')
! do pstate=1,nstate
!   write(111,'(*(2x,f12.6))') p_matrix_gw(pstate,:)
! enddo
! close(111)

 ! Add the SCF density matrix to get to the total density matrix
 do pstate=1,nstate
   p_matrix_gw(pstate,pstate) = p_matrix_gw(pstate,pstate) + occupation(pstate,pqspin)
 enddo

 p_matrix(:,:,pqspin) = MATMUL( c_matrix(:,:,pqspin)  , &
                          MATMUL( p_matrix_gw(:,:), &
                             TRANSPOSE(c_matrix(:,:,pqspin)) ) )


 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 else
   call destroy_eri_4center_eigen_uks()
 endif


end subroutine gw_density_matrix


!=========================================================================
