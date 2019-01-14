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
 real(dp),intent(inout)     :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer                 :: pstate,istate,astate
 integer                 :: iaspin
 real(dp)                :: denom
 real(dp)                :: p_matrix_pt1(nstate,nstate)
!=====

 call start_clock(timing_mbpt_dm)

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
! call issue_warning('pt1_density_matrix: this is not correct when starting from something else than HF')
! do pstate=ncore_G+1,nvirtual_G-1
!   p_matrix_pt1(pstate,pstate) = p_matrix_pt1(pstate,pstate) + occupation(pstate,iaspin)
! enddo

 ! Transform from MO to AO
 p_matrix(:,:,iaspin) = MATMUL( c_matrix(:,:,iaspin)  , &
                          MATMUL( p_matrix_pt1(:,:), &
                             TRANSPOSE(c_matrix(:,:,iaspin)) ) )

 call stop_clock(timing_mbpt_dm)

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
 real(dp),intent(inout)     :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====
 integer                 :: file_density_matrix
 integer                 :: pstate,qstate
 integer                 :: istate,jstate,kstate
 integer                 :: astate,bstate,cstate
 integer                 :: pqspin
 real(dp)                :: denom1,denom2
 real(dp)                :: num1,num2
 real(dp)                :: p_matrix_pt2(nstate,nstate,nspin)
!=====


 call start_clock(timing_mbpt_dm)

 write(stdout,'(/,a)') ' Calculate the PT2 density matrix'

 if( nspin /= 1 ) call die('pt2_density_matrix: only implemented for spin restricted calculations')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif


 ! Full calculation of the PT2 density matrix

 p_matrix_pt2(:,:,:) = 0.0_dp
 ! so far, only spin-restricted calculation are possible
 pqspin = 1

 ! A1 P_ij sum over k and a,b
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

         p_matrix_pt2(istate,jstate,pqspin) = p_matrix_pt2(istate,jstate,pqspin)  &
                 - num1 * num2 / ( denom1 * denom2 )

       enddo
     enddo
   enddo
 enddo
 enddo

 ! A2 P_ab sum over  i,j  and c
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

         p_matrix_pt2(astate,bstate,pqspin) = p_matrix_pt2(astate,bstate,pqspin)  &
                 + num1 * num2 / ( denom1 * denom2 )

       enddo
     enddo
   enddo
 enddo
 enddo

 ! A3    P_cj  sum over i, and a,b
 ! A4    P_jc  sum over i, and a,b
 do cstate=nhomo_G+1,nvirtual_G-1
 do jstate=ncore_G+1,nhomo_G
   do astate=nhomo_G+1,nvirtual_G-1
     do bstate=nhomo_G+1,nvirtual_G-1
       do istate=ncore_G+1,nhomo_G
         denom1 = energy(jstate,pqspin) + energy(istate,pqspin) - energy(astate,pqspin) - energy(bstate,pqspin)
         denom2 = energy(jstate,pqspin) - energy(cstate,pqspin)
         num1 = 2.0_dp * eri_eigen(jstate,astate,pqspin,istate,bstate,pqspin) - eri_eigen(jstate,bstate,pqspin,istate,astate,pqspin)
         num2 = 2.0_dp * eri_eigen(astate,cstate,pqspin,bstate,istate,pqspin)

         p_matrix_pt2(cstate,jstate,pqspin) = p_matrix_pt2(cstate,jstate,pqspin)  &
                                           + num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(jstate,cstate,pqspin) = p_matrix_pt2(jstate,cstate,pqspin)  &
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

         p_matrix_pt2(bstate,kstate,pqspin) = p_matrix_pt2(bstate,kstate,pqspin)  &
                                           - num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(kstate,bstate,pqspin) = p_matrix_pt2(kstate,bstate,pqspin)  &
                                           - num1 * num2 / ( denom1 * denom2 )
       enddo
     enddo
   enddo
 enddo
 enddo

 if( ALL( ABS(p_matrix) < 1.0e-6_dp ) ) then
   if( TRIM(calc_type%scf_name) /= 'HF' ) &
               call issue_warning('pt2_density_matrix: this is not correct when starting from something else than HF')
   ! Add the SCF density matrix to get to the total density matrix
   do pstate=1,nstate
     p_matrix_pt2(pstate,pstate,pqspin) = p_matrix_pt2(pstate,pstate,pqspin) + occupation(pstate,pqspin)
   enddo
 endif

 ! Transform from MO to AO
 !p_matrix(:,:,pqspin) = p_matrix(:,:,pqspin) + MATMUL( c_matrix(:,:,pqspin)  , &
 !                                                      MATMUL( p_matrix_pt2(:,:), TRANSPOSE(c_matrix(:,:,pqspin)) ) )
 call matrix_mo_to_ao(c_matrix,p_matrix_pt2,p_matrix)

 if( print_density_matrix_ .AND. is_iomaster ) then
   write(stdout,'(1x,a)') 'Write DENSITY_MATRIX file'
   open(newunit=file_density_matrix,file='DENSITY_MATRIX',form='unformatted',action='write')
   do pqspin=1,nspin
     write(file_density_matrix) p_matrix(:,:,pqspin)
   enddo
   close(file_density_matrix)
 endif


 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 else
   call destroy_eri_4center_eigen_uks()
 endif

 call stop_clock(timing_mbpt_dm)

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
 integer                 :: file_density_matrix
 integer                 :: pstate,qstate
 integer                 :: istate,jstate,kstate
 integer                 :: astate,bstate,cstate
 integer                 :: pqspin
 real(dp)                :: denom1,denom2
 real(dp)                :: num1,num2
 real(dp)                :: p_matrix_pt2(nstate,nstate,nspin)
!=====


 call start_clock(timing_mbpt_dm)

 write(stdout,'(/,a)') ' Calculate the 1-ring density matrix'

 if( nspin /= 1 ) call die('pt2_density_matrix: only implemented for spin restricted calculations')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif


! Full calculation of the 1-ring density matrix

 p_matrix_pt2(:,:,:) = 0.0_dp
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

         p_matrix_pt2(istate,jstate,pqspin) = p_matrix_pt2(istate,jstate,pqspin)  &
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

         p_matrix_pt2(astate,bstate,pqspin) = p_matrix_pt2(astate,bstate,pqspin)  &
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

         p_matrix_pt2(cstate,jstate,pqspin) = p_matrix_pt2(cstate,jstate,pqspin)  &
                                           + num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(jstate,cstate,pqspin) = p_matrix_pt2(jstate,cstate,pqspin)  &
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

         p_matrix_pt2(bstate,kstate,pqspin) = p_matrix_pt2(bstate,kstate,pqspin)  &
                                           - num1 * num2 / ( denom1 * denom2 )
         p_matrix_pt2(kstate,bstate,pqspin) = p_matrix_pt2(kstate,bstate,pqspin)  &
                                           - num1 * num2 / ( denom1 * denom2 )
       enddo
     enddo
   enddo
 enddo
 enddo



 if( ALL( ABS(p_matrix) < 1.0e-6_dp ) ) then
   if( TRIM(calc_type%scf_name) /= 'HF' ) &
               call issue_warning('onering_density_matrix: this is not correct when starting from something else than HF')
   ! Add the SCF density matrix to get to the total density matrix
   do pstate=1,nstate
     p_matrix_pt2(pstate,pstate,pqspin) = p_matrix_pt2(pstate,pstate,pqspin) + occupation(pstate,pqspin)
   enddo
 endif

 ! Transform from MO to AO
 !p_matrix(:,:,pqspin) = p_matrix(:,:,pqspin) + MATMUL( c_matrix(:,:,pqspin)  , &
 !                                                      MATMUL( p_matrix_pt2(:,:), TRANSPOSE(c_matrix(:,:,pqspin)) ) )
 call matrix_mo_to_ao(c_matrix,p_matrix_pt2,p_matrix)

 if( print_density_matrix_ .AND. is_iomaster ) then
   write(stdout,'(1x,a)') 'Write DENSITY_MATRIX file'
   open(newunit=file_density_matrix,file='DENSITY_MATRIX',form='unformatted',action='write')
   do pqspin=1,nspin
     write(file_density_matrix) p_matrix(:,:,pqspin)
   enddo
   close(file_density_matrix)
 endif


 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 else
   call destroy_eri_4center_eigen_uks()
 endif

 call stop_clock(timing_mbpt_dm)

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
 integer  :: pstate,qstate
 integer  :: istate,jstate
 integer  :: astate,bstate
 integer  :: pqspin
 integer  :: ipole
 integer  :: npole_local,ipole_local
 integer  :: nstate_occ,nstate_virt
 integer  :: file_density_matrix
 real(dp) :: p_matrix_gw(nstate,nstate,nspin)
 real(dp),allocatable :: bra_occ(:,:),bra_virt(:,:)
 real(dp),allocatable :: bra_occ_local(:,:),bra_virt_local(:,:)
!=====

 call start_clock(timing_mbpt_dm)

 write(stdout,'(/,a)') ' Calculate the GW density matrix'

 if( nspin /= 1 ) call die('gw_density_matrix: only implemented for spin restricted calculations')
 if( .NOT. has_auxil_basis)  call die('gw_density_matrix: only implemented without RI')

 if(has_auxil_basis) then
   call calculate_eri_3center_eigen(c_matrix,ncore_G+1,nvirtual_G-1,ncore_G+1,nvirtual_G-1)
 else
   call calculate_eri_4center_eigen_uks(c_matrix,ncore_G+1,nvirtual_G-1)
 endif


 ! First order calculation of the GW density matrix

 p_matrix_gw(:,:,:) = 0.0_dp
 pqspin = 1

 nstate_occ  = nhomo_G - ncore_G
 nstate_virt = nvirtual_G - nhomo_G - 1
 allocate(bra_occ(wpol%npole_reso,ncore_G+1:nhomo_G))
 allocate(bra_virt(wpol%npole_reso,nhomo_G+1:nvirtual_G-1))
 npole_local = NUMROC(wpol%npole_reso,1,rank_auxil,0,nproc_auxil)
 allocate(bra_occ_local(npole_local,ncore_G+1:nhomo_G))
 allocate(bra_virt_local(npole_local,nhomo_G+1:nvirtual_G-1))

 do astate=nhomo_G+1,nvirtual_G-1
   if( MODULO( astate - (nhomo_G+1) , nproc_ortho ) /= rank_ortho ) cycle

   ! A1
   !bra_occ(:,ncore_G+1:nhomo_G) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,ncore_G+1:nhomo_G,astate,pqspin) )
   call DGEMM('T','N',wpol%npole_reso,nstate_occ,nauxil_3center, &
                         1.0d0,wpol%residue_left,nauxil_3center, &
                               eri_3center_eigen(1,ncore_G+1,astate,pqspin),nauxil_3center, &
                         0.0_dp,bra_occ(1,ncore_G+1),wpol%npole_reso)
   call xsum_auxil(bra_occ)

   ipole_local = 0
   do ipole=1,wpol%npole_reso
     if( MODULO( ipole-1 , nproc_auxil ) /= rank_auxil ) cycle
     ipole_local = ipole_local + 1
     do jstate=ncore_G+1,nhomo_G
       bra_occ_local(ipole_local,jstate) = bra_occ(ipole,jstate) &
                                            / ( energy(jstate,pqspin) - energy(astate,pqspin) - wpol%pole(ipole) )
     enddo
   enddo

   call DSYRK('U','T',nstate_occ,npole_local,-2.0d0,bra_occ_local,npole_local,1.0d0,p_matrix_gw(ncore_G+1,ncore_G+1,pqspin),nstate)


   ! A3    P_cj  sum over i,a,b
   ! A4    P_jc  sum over i,a,b     ! not actually calculated, but included through the symmetrization step
   !bra_virt(:,nhomo_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,nhomo_G+1:nvirtual_G-1,astate,pqspin) )
   call DGEMM('T','N',wpol%npole_reso,nstate_virt,nauxil_3center, &
                         1.0d0,wpol%residue_left,nauxil_3center,  &
                               eri_3center_eigen(1,nhomo_G+1,astate,pqspin),nauxil_3center, &
                         0.0_dp,bra_virt(1,nhomo_G+1),wpol%npole_reso)
   call xsum_auxil(bra_virt)

   ipole_local = 0
   do ipole=1,wpol%npole_reso
     if( MODULO( ipole-1 , nproc_auxil ) /= rank_auxil ) cycle
     ipole_local = ipole_local + 1
     do bstate=nhomo_G+1,nvirtual_G-1
       bra_virt_local(ipole_local,bstate) = bra_virt(ipole,bstate)
     enddo
   enddo

   call DGEMM('T','N',nstate_occ,nstate_virt,npole_local,2.0d0,bra_occ_local,npole_local, &
                                                               bra_virt_local,npole_local, &
                                                         1.0d0,p_matrix_gw(ncore_G+1,nhomo_G+1,pqspin),nstate)
 enddo

 do istate=ncore_G+1,nhomo_G
   if( MODULO( istate - (ncore_G+1) , nproc_ortho ) /= rank_ortho ) cycle

   ! A2
   !bra_virt(:,nhomo_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,nhomo_G+1:nvirtual_G-1,istate,pqspin) )
   call DGEMM('T','N',wpol%npole_reso,nstate_virt,nauxil_3center, &
                         1.0d0,wpol%residue_left,nauxil_3center,  &
                               eri_3center_eigen(1,nhomo_G+1,istate,pqspin),nauxil_3center, &
                         0.0_dp,bra_virt(1,nhomo_G+1),wpol%npole_reso)
   call xsum_auxil(bra_virt)

   ipole_local = 0
   do ipole=1,wpol%npole_reso
     if( MODULO( ipole-1 , nproc_auxil ) /= rank_auxil ) cycle
     ipole_local = ipole_local + 1
     do bstate=nhomo_G+1,nvirtual_G-1
       bra_virt_local(ipole_local,bstate) = bra_virt(ipole,bstate) &
                                            / ( energy(istate,pqspin) - energy(bstate,pqspin) - wpol%pole(ipole) )
     enddo
   enddo

   call DSYRK('U','T',nstate_virt,npole_local,2.0d0,bra_virt_local,npole_local,1.0d0,p_matrix_gw(nhomo_G+1,nhomo_G+1,pqspin),nstate)

   ! A5   P_bk  sum over i,j,a
   ! A6   P_kb  sum over i,j,a   ! not actually calculated, but included through the symmetrization step
   !bra_occ(:,ncore_G+1:nhomo_G)       = MATMUL( TRANSPOSE(wpol%residue_left(:,:)) , eri_3center_eigen(:,ncore_G+1:nhomo_G,istate,pqspin) )
   call DGEMM('T','N',wpol%npole_reso,nstate_occ,nauxil_3center, &
                         1.0d0,wpol%residue_left,nauxil_3center, &
                               eri_3center_eigen(1,ncore_G+1,istate,pqspin),nauxil_3center, &
                         0.0_dp,bra_occ(1,ncore_G+1),wpol%npole_reso)
   call xsum_auxil(bra_occ)

   ipole_local = 0
   do ipole=1,wpol%npole_reso
     if( MODULO( ipole-1 , nproc_auxil ) /= rank_auxil ) cycle
     ipole_local = ipole_local + 1
     do jstate=ncore_G+1,nhomo_G
       bra_occ_local(ipole_local,jstate) = bra_occ(ipole,jstate)
     enddo
   enddo

   call DGEMM('T','N',nstate_occ,nstate_virt,npole_local,-2.0d0,bra_occ_local,npole_local, &
                                                                bra_virt_local,npole_local, &
                                                          1.0d0,p_matrix_gw(ncore_G+1,nhomo_G+1,pqspin),nstate)
 enddo

 ! A common factor 1/(e_j-e_c) is to be added for the occupied-virtual block (terms A3,A4,A5,A6)
 do bstate=nhomo_G+1,nvirtual_G-1
   do jstate=ncore_G+1,nhomo_G
     p_matrix_gw(jstate,bstate,pqspin) = p_matrix_gw(jstate,bstate,pqspin) / ( energy(jstate,pqspin) - energy(bstate,pqspin) )
   enddo
 enddo

 call xsum_world(p_matrix_gw)

 deallocate(bra_occ)
 deallocate(bra_virt)
 deallocate(bra_occ_local)
 deallocate(bra_virt_local)

 ! Symmetrization of the p_matrix here
 ! Only the upper triangle was set up before
 do pstate=1,nstate
   do qstate=pstate+1,nstate
     p_matrix_gw(qstate,pstate,pqspin) = p_matrix_gw(pstate,qstate,pqspin)
   enddo
 enddo


 if( ALL( ABS(p_matrix) < 1.0e-6_dp ) ) then
   if( TRIM(calc_type%scf_name) /= 'HF' ) &
               call issue_warning('gw_density_matrix: this is not correct when starting from something else than HF')
   ! Add the SCF density matrix to get to the total density matrix
   do pstate=1,nstate
     p_matrix_gw(pstate,pstate,pqspin) = p_matrix_gw(pstate,pstate,pqspin) + occupation(pstate,pqspin)
   enddo
 endif

 !
 ! Transform from MO to AO
 !p_matrix(:,:,pqspin) = p_matrix(:,:,pqspin) + MATMUL( c_matrix(:,:,pqspin)  , &
 !                                                      MATMUL( p_matrix_gw(:,:), TRANSPOSE(c_matrix(:,:,pqspin)) ) )
 call matrix_mo_to_ao(c_matrix,p_matrix_gw,p_matrix)

 if( print_density_matrix_ .AND. is_iomaster ) then
   write(stdout,'(1x,a)') 'Write DENSITY_MATRIX file'
   open(newunit=file_density_matrix,file='DENSITY_MATRIX',form='unformatted',action='write')
   do pqspin=1,nspin
     write(file_density_matrix) p_matrix(:,:,pqspin)
   enddo
   close(file_density_matrix)
 endif

 if(has_auxil_basis) then
   call destroy_eri_3center_eigen()
 else
   call destroy_eri_4center_eigen_uks()
 endif

 call stop_clock(timing_mbpt_dm)

end subroutine gw_density_matrix


!=========================================================================
