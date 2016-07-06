!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the construction the BSE hamiltonian
! or alternatively, the TDDFT "Casida" equations
!
!=========================================================================
subroutine build_amb_apb_common(desc_apb,nmat,nbf,nstate,c_matrix,energy,wpol,alpha_local, &
                                m_apb,n_apb,amb_matrix,apb_matrix,amb_diag_rpa,rpa_correlation)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_tools 
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                 :: desc_apb(ndel),nmat,nbf,nstate
 real(dp),intent(in)                :: energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(in)                :: alpha_local
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(inout)             :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
 real(dp),intent(out)               :: amb_diag_rpa(nmat)
 real(dp),intent(out)               :: rpa_correlation
!=====
 integer              :: t_ia,t_jb,t_ia_global,t_jb_global
 integer              :: istate,astate,jstate,bstate
 integer              :: iaspin,jbspin
 integer              :: jbmin
 real(dp),allocatable :: eri_eigenstate_jbmin(:,:,:,:)
 real(dp)             :: eri_eigen_iajb
 real(dp)             :: eri_eigen_ijab,eri_eigen_ibaj
 logical              :: j_is_jbmin
 integer              :: iprow,ipcol
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: amb_block(:,:)
 real(dp),allocatable :: apb_block(:,:)
!=====

 call start_clock(timing_build_common)

 write(stdout,'(a)') ' Build Common part: Energies + Hartree + possibly Exchange'
 write(stdout,'(a,f8.3)') ' Content of Exchange: ',alpha_local

 if( .NOT. has_auxil_basis) then
   allocate(eri_eigenstate_jbmin(nstate,nstate,nstate,nspin))
   ! Set this to zero and then enforce the calculation of the first series of
   ! Coulomb integrals
   eri_eigenstate_jbmin(:,:,:,:) = 0.0_dp
 endif

 if( nprow_sd * npcol_sd > 1 ) &
    write(stdout,'(a,i4,a,i4)') ' SCALAPACK grid    :',nprow_sd,' x ',npcol_sd

 rpa_correlation = 0.0_dp
 !
 ! Set up energy+hartree+exchange contributions to matrices (A+B) and (A-B) 
 !
 
 ! First loops over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_apb_block = row_block_size(nmat,iprow,nprow_sd)
     n_apb_block = col_block_size(nmat,ipcol,npcol_sd)

     allocate(amb_block(m_apb_block,n_apb_block))
     allocate(apb_block(m_apb_block,n_apb_block))
     amb_block(:,:) = 0.0_dp
     apb_block(:,:) = 0.0_dp

     ! Then loop inside each of the SCALAPACK blocks
     do t_jb=1,n_apb_block
       t_jb_global = colindex_local_to_global(ipcol,npcol_sd,t_jb)
       jstate = wpol%transition_table_apb(1,t_jb_global)
       bstate = wpol%transition_table_apb(2,t_jb_global)
       jbspin = wpol%transition_table_apb(3,t_jb_global)
    
       if( .NOT. has_auxil_basis ) then
         jbmin = MIN(jstate,bstate)
         j_is_jbmin = (jbmin == jstate)
         call calculate_eri_4center_eigen(nbf,nstate,c_matrix,jbmin,jbspin,eri_eigenstate_jbmin)
       endif
    
       do t_ia=1,m_apb_block
         t_ia_global = rowindex_local_to_global(iprow,nprow_sd,t_ia)
         istate = wpol%transition_table_apb(1,t_ia_global)
         astate = wpol%transition_table_apb(2,t_ia_global)
         iaspin = wpol%transition_table_apb(3,t_ia_global)
    
         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ia_global < t_jb_global ) cycle
    
         if(has_auxil_basis) then
           eri_eigen_iajb = eri_eigen_ri_paral(istate,astate,iaspin,jstate,bstate,jbspin)
         else
           if(j_is_jbmin) then ! treating (j,b)
             eri_eigen_iajb = eri_eigenstate_jbmin(bstate,istate,astate,iaspin)
           else                ! treating (b,j)
             eri_eigen_iajb = eri_eigenstate_jbmin(jstate,istate,astate,iaspin)
           endif
         endif
    
         if( .NOT. is_triplet) then
           apb_block(t_ia,t_jb) = 2.0_dp * eri_eigen_iajb * spin_fact
         else
           apb_block(t_ia,t_jb) = 0.0_dp
         endif
         amb_block(t_ia,t_jb) = 0.0_dp
    
         if( alpha_local > 1.0e-6_dp ) then
           if( iaspin == jbspin ) then
             if(has_auxil_basis) then
               eri_eigen_ijab = eri_eigen_ri_paral(istate,jstate,iaspin,astate,bstate,jbspin)
               eri_eigen_ibaj = eri_eigen_ri_paral(istate,bstate,iaspin,astate,jstate,jbspin)
             else
               if(j_is_jbmin) then
                 eri_eigen_ijab = eri_eigenstate_jbmin(istate,astate,bstate,jbspin)
                 eri_eigen_ibaj = eri_eigenstate_jbmin(astate,istate,bstate,jbspin)
               else
                 eri_eigen_ijab = eri_eigenstate_jbmin(astate,istate,jstate,jbspin)
                 eri_eigen_ibaj = eri_eigenstate_jbmin(istate,astate,jstate,jbspin)
               endif
             endif
             apb_block(t_ia,t_jb) = apb_block(t_ia,t_jb) - eri_eigen_ijab * alpha_local - eri_eigen_ibaj * alpha_local
             amb_block(t_ia,t_jb) = amb_block(t_ia,t_jb) - eri_eigen_ijab * alpha_local + eri_eigen_ibaj * alpha_local
           endif
         endif
    
         if( t_ia_global == t_jb_global ) then
           !
           ! Only one proc should add the diagonal
           if( rank_auxil == 0 ) then
             apb_block(t_ia,t_jb) =  apb_block(t_ia,t_jb) + ( energy(bstate,jbspin) - energy(jstate,jbspin) )
             amb_block(t_ia,t_jb) =  amb_block(t_ia,t_jb) + ( energy(bstate,jbspin) - energy(jstate,jbspin) )
           endif
           ! First part of the RPA correlation energy: sum over diagonal terms
           rpa_correlation = rpa_correlation - 0.25_dp * ( apb_block(t_ia,t_jb) + amb_block(t_ia,t_jb) )
         endif
    
       enddo 
    
     enddo 


     call xsum_auxil(amb_block)
     call xsum_auxil(apb_block)

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       amb_matrix(:,:) = amb_block(:,:)
       apb_matrix(:,:) = apb_block(:,:)
     endif


     deallocate(amb_block)
     deallocate(apb_block)
   enddo 
 enddo 

#ifdef HAVE_SCALAPACK
 call xsum_world(rpa_correlation)
#endif

 !
 ! Set up the diagonal of A-B in the RPA approximation
 ! 
 do t_ia_global=1,nmat
   istate = wpol%transition_table_apb(1,t_ia_global)
   astate = wpol%transition_table_apb(2,t_ia_global)
   iaspin = wpol%transition_table_apb(3,t_ia_global)

   amb_diag_rpa(t_ia_global) = energy(astate,iaspin) - energy(istate,iaspin)

 enddo


 if(ALLOCATED(eri_eigenstate_jbmin)) deallocate(eri_eigenstate_jbmin)


 call stop_clock(timing_build_common)

end subroutine build_amb_apb_common


!=========================================================================
subroutine build_amb_apb_diag_auxil(nmat,nstate,energy,wpol,m_apb,n_apb,amb_matrix,apb_matrix,amb_diag_rpa)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_tools 
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                 :: nmat,nstate
 real(dp),intent(in)                :: energy(nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(inout)             :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
 real(dp),intent(out)               :: amb_diag_rpa(nmat)
!=====
 integer              :: t_ia,t_jb,t_jb_global
 integer              :: jstate,bstate
 integer              :: jbspin
!=====

 write(stdout,'(a)') ' Build diagonal part with auxil basis: Energies'

 do t_jb_global=1,nmat
   t_ia = rowindex_global_to_local('S',t_jb_global)
   t_jb = colindex_global_to_local('S',t_jb_global)

   jstate = wpol%transition_table_apb(1,t_jb_global)
   bstate = wpol%transition_table_apb(2,t_jb_global)
   jbspin = wpol%transition_table_apb(3,t_jb_global)
   amb_diag_rpa(t_jb_global) = energy(bstate,jbspin) - energy(jstate,jbspin)

   ! If the diagonal element belongs to this proc, then add it.
   if( t_ia > 0 .AND. t_jb > 0 ) then
     apb_matrix(t_ia,t_jb) =  amb_diag_rpa(t_jb_global)
     amb_matrix(t_ia,t_jb) =  amb_diag_rpa(t_jb_global)
   endif
 enddo


end subroutine build_amb_apb_diag_auxil


!=========================================================================
subroutine get_rpa_correlation(nmat,wpol,m_apb,n_apb,amb_matrix,apb_matrix,rpa_correlation)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_tools 
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                 :: nmat
 type(spectral_function),intent(in) :: wpol
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(in)                :: amb_matrix(m_apb,n_apb)
 real(dp),intent(in)                :: apb_matrix(m_apb,n_apb)
 real(dp),intent(out)               :: rpa_correlation
!=====
 integer                            :: t_ia,t_jb,t_jb_global
!=====

 rpa_correlation = 0.0_dp

 do t_jb_global=1,nmat
   t_ia = rowindex_global_to_local('S',t_jb_global)
   t_jb = colindex_global_to_local('S',t_jb_global)

   ! If the diagonal element belongs to this proc, then add it.
   if( t_ia > 0 .AND. t_jb > 0 ) then
     rpa_correlation = rpa_correlation - 0.25_dp * apb_matrix(t_ia,t_jb)   &
                                       - 0.25_dp * amb_matrix(t_ia,t_jb) 
   endif
 enddo

 call xsum_world(rpa_correlation)


end subroutine get_rpa_correlation


!=========================================================================
subroutine build_apb_hartree_auxil(desc_apb,wpol,m_apb,n_apb,apb_matrix)
 use m_definitions
 use m_timing
 use m_mpi
 use m_spectral_function
 use m_eri_ao_mo
 implicit none

 integer,intent(in)                 :: desc_apb(ndel)
 type(spectral_function),intent(in) :: wpol
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(inout)             :: apb_matrix(m_apb,n_apb)
!=====
 integer              :: nmat
 integer              :: t_ia,t_jb,t_ia_global,t_jb_global
 integer              :: istate,astate,jstate,bstate
 integer              :: iaspin,jbspin
 integer              :: iprow,ipcol
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: apb_block(:,:)
 integer              :: ibf_auxil
 real(dp),allocatable :: eri_3center_left(:),eri_3center_right(:)
!=====

 if( is_triplet) return

 call start_clock(timing_build_common)

 write(stdout,'(a)') ' Build Hartree part with auxil basis'

 nmat = desc_apb(3)


 if( nprow_sd * npcol_sd > 1 ) &
    write(stdout,'(a,i4,a,i4)') ' SCALAPACK grid    :',nprow_sd,' x ',npcol_sd


 
 ! First loops over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_apb_block = row_block_size(nmat,iprow,nprow_sd)
     n_apb_block = col_block_size(nmat,ipcol,npcol_sd)

     ! If block is of size 0, then skip it
     if( m_apb_block == 0 .OR. n_apb_block == 0 ) cycle

     allocate(apb_block(m_apb_block,n_apb_block))

     allocate(eri_3center_left(m_apb_block))
     allocate(eri_3center_right(n_apb_block))

     apb_block(:,:) = 0.0_dp


     do ibf_auxil=1,nauxil_3center


       do t_jb=1,n_apb_block
         t_jb_global = colindex_local_to_global(ipcol,npcol_sd,t_jb)
         jstate = wpol%transition_table_apb(1,t_jb_global)
         bstate = wpol%transition_table_apb(2,t_jb_global)
         jbspin = wpol%transition_table_apb(3,t_jb_global)

         eri_3center_right(t_jb) = eri_3center_eigen(ibf_auxil,jstate,bstate,jbspin)
  
       enddo
    
       do t_ia=1,m_apb_block
         t_ia_global = rowindex_local_to_global(iprow,nprow_sd,t_ia)
         istate = wpol%transition_table_apb(1,t_ia_global)
         astate = wpol%transition_table_apb(2,t_ia_global)
         iaspin = wpol%transition_table_apb(3,t_ia_global)
    
         eri_3center_left(t_ia) = eri_3center_eigen(ibf_auxil,istate,astate,iaspin)

       enddo

       call DGER(m_apb_block,n_apb_block,2.0_dp*spin_fact,eri_3center_left,1,eri_3center_right,1,apb_block,m_apb_block)

     enddo
     
     deallocate(eri_3center_left,eri_3center_right)

     call xsum_auxil(apb_block)

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       apb_matrix(:,:) = apb_matrix(:,:) + apb_block(:,:)
     endif


     deallocate(apb_block)
   enddo 
 enddo 


 call stop_clock(timing_build_common)

end subroutine build_apb_hartree_auxil


!=========================================================================
subroutine build_a_diag_common(nmat,nbf,nstate,c_matrix,energy,wpol,a_diag)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_spectral_function
 use m_eri_ao_mo
 use m_tools 
 implicit none

 integer,intent(in)                 :: nmat,nbf,nstate
 real(dp),intent(in)                :: energy(nstate,nspin)
 real(dp),intent(in)                :: c_matrix(nbf,nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 real(dp),intent(out)               :: a_diag(wpol%npole_reso_spa)
!=====
 integer              :: t_jb,t_jb_global
 integer              :: jstate,bstate
 integer              :: jbspin
 integer              :: jbmin
 real(dp),allocatable :: eri_eigenstate_jbmin(:,:,:,:)
 real(dp)             :: eri_eigen_jbjb
 logical              :: j_is_jbmin
 real(dp),parameter   :: empirical_fact=1.50_dp
 character(len=100)   :: ctmp
!=====

 call start_clock(timing_build_common)

 write(stdout,'(a)') ' Build diagonal of the RPA part: Energies + Hartree'
 if( ABS( empirical_fact - 1.0_dp ) > 1.0e-6_dp ) then
   write(ctmp,'(a,1x,f6.2)') 'Empirical parameter',empirical_fact
   call issue_warning(ctmp)
 endif

 if( .NOT. has_auxil_basis) then
   allocate(eri_eigenstate_jbmin(nstate,nstate,nstate,nspin))
   ! Set this to zero and then enforce the calculation of the first series of
   ! Coulomb integrals
   eri_eigenstate_jbmin(:,:,:,:) = 0.0_dp
 endif

 !
 ! Set up energy+hartree+exchange contributions to matrices (A+B) and (A-B) 
 !


 ! Then loop inside each of the SCALAPACK blocks
 do t_jb=1,wpol%npole_reso_spa
   jstate = wpol%transition_table_spa(1,t_jb)
   bstate = wpol%transition_table_spa(2,t_jb)
   jbspin = wpol%transition_table_spa(3,t_jb)

   if( .NOT. has_auxil_basis ) then
     jbmin = MIN(jstate,bstate)
     j_is_jbmin = (jbmin == jstate)
     call calculate_eri_4center_eigen(nbf,nstate,c_matrix,jbmin,jbspin,eri_eigenstate_jbmin)
   endif

   if(has_auxil_basis) then
     eri_eigen_jbjb = eri_eigen_ri(jstate,bstate,jbspin,jstate,bstate,jbspin)
   else
     if(j_is_jbmin) then ! treating (j,b)
       eri_eigen_jbjb = eri_eigenstate_jbmin(bstate,jstate,bstate,jbspin)
     else                ! treating (b,j)
       eri_eigen_jbjb = eri_eigenstate_jbmin(jstate,jstate,bstate,jbspin)
     endif
   endif

   a_diag(t_jb) = eri_eigen_jbjb * spin_fact + energy(bstate,jbspin) - energy(jstate,jbspin)
   a_diag(t_jb) = a_diag(t_jb) * empirical_fact

 enddo 

 if(ALLOCATED(eri_eigenstate_jbmin)) deallocate(eri_eigenstate_jbmin)


 call stop_clock(timing_build_common)

end subroutine build_a_diag_common


!=========================================================================
subroutine build_apb_tddft(nmat,nstate,basis,c_matrix,occupation,wpol,m_apb,n_apb,apb_matrix)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_spectral_function
 use m_basis_set
 use m_dft_grid
 use m_tddft_fxc
 implicit none

 integer,intent(in)                 :: nmat,nstate
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: occupation(nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(inout)             :: apb_matrix(m_apb,n_apb)
!=====
 integer              :: t_ia,t_jb,t_ia_global,t_jb_global
 integer              :: istate,astate,jstate,bstate
 integer              :: iaspin,jbspin
 real(dp)             :: xctmp
 integer              :: ipcol,iprow
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: apb_block(:,:)
!=====

 call start_clock(timing_build_tddft)

 write(stdout,'(a)') ' Build fxc part'

 !
 ! Prepare TDDFT calculations
 call prepare_tddft(nstate,basis,c_matrix,occupation)


 if( nprow_sd * npcol_sd > 1 ) &
    write(stdout,'(a,i4,a,i4)') ' SCALAPACK grid    :',nprow_sd,' x ',npcol_sd


 ! First loops over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_apb_block = row_block_size(nmat,iprow,nprow_sd)
     n_apb_block = col_block_size(nmat,ipcol,npcol_sd)

     allocate(apb_block(m_apb_block,n_apb_block))
     apb_block(:,:) = 0.0_dp

     !
     ! Set up fxc contributions to matrices (A+B)
     !
     do t_jb=1,n_apb_block
       t_jb_global = colindex_local_to_global(ipcol,npcol_sd,t_jb)
       jstate = wpol%transition_table_apb(1,t_jb_global)
       bstate = wpol%transition_table_apb(2,t_jb_global)
       jbspin = wpol%transition_table_apb(3,t_jb_global)
    
       do t_ia=1,m_apb_block
         t_ia_global = rowindex_local_to_global(iprow,nprow_sd,t_ia)
         istate = wpol%transition_table_apb(1,t_ia_global)
         astate = wpol%transition_table_apb(2,t_ia_global)
         iaspin = wpol%transition_table_apb(3,t_ia_global)
    
         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ia_global < t_jb_global ) cycle
    
         if( nspin == 1 ) then 
           if( .NOT. is_triplet ) then
             xctmp = eval_fxc_rks_singlet(istate,astate,iaspin,jstate,bstate,jbspin)
           else
             xctmp =  eval_fxc_rks_triplet(istate,astate,iaspin,jstate,bstate,jbspin)
           endif
         else
           xctmp =  eval_fxc_uks(istate,astate,iaspin,jstate,bstate,jbspin)
         endif
    
    
         ! The factor two accounts for (A+B), and not A or B.
         apb_block(t_ia,t_jb) = apb_block(t_ia,t_jb) + 2.0_dp * xctmp
    
       enddo
     enddo

     !
     ! real-space integration grid is distributed, one needs to sum contributions here
     call xsum_grid(apb_block)

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       apb_matrix(:,:) = apb_matrix(:,:) + apb_block(:,:)
     endif
     deallocate(apb_block)

   enddo
 enddo


 call destroy_tddft()

 call stop_clock(timing_build_tddft)

end subroutine build_apb_tddft


!=========================================================================
subroutine build_amb_apb_bse(nbf,nstate,wpol,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_spectral_function
 use m_basis_set
 use m_eri_ao_mo
 use m_tools 
 implicit none

 integer,intent(in)                 :: nbf,nstate
 type(spectral_function),intent(in) :: wpol,wpol_static
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(inout)             :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
!=====
 integer              :: t_ia,t_jb,t_ia_global,t_jb_global
 integer              :: istate,astate,jstate,bstate
 integer              :: iaspin,jbspin
 integer              :: nprodbasis
 integer              :: kbf
 real(dp),allocatable :: bra(:),ket(:)
 real(dp)             :: wtmp
!=====

 if( has_auxil_basis ) call die('Have an auxil basis. This should not happen here')

 call start_clock(timing_build_bse)

 write(stdout,'(a)') ' Build W part'

 nprodbasis = index_prodstate(nstate,nstate)

 !
 ! Prepare the bra and ket for BSE
 allocate(bra(wpol_static%npole_reso))
 allocate(ket(wpol_static%npole_reso))

 !
 ! Set up -W contributions to matrices (A+B) and (A-B)
 !
 do t_jb=1,n_apb
   t_jb_global = colindex_local_to_global('C',t_jb)
   jstate = wpol%transition_table_apb(1,t_jb_global)
   bstate = wpol%transition_table_apb(2,t_jb_global)
   jbspin = wpol%transition_table_apb(3,t_jb_global)

   do t_ia=1,m_apb
     t_ia_global = rowindex_local_to_global('C',t_ia)
     istate = wpol%transition_table_apb(1,t_ia_global)
     astate = wpol%transition_table_apb(2,t_ia_global)
     iaspin = wpol%transition_table_apb(3,t_ia_global)

     !
     ! Only calculate the lower triangle
     ! Symmetrization will be performed later (in the diago subroutines)
     if( t_ia_global < t_jb_global ) cycle

     if( iaspin /= jbspin ) cycle

     kbf = index_prodstate(istate,jstate) + nprodbasis * (iaspin-1)
     bra(:) = wpol_static%residu_left(kbf,:)
     kbf = index_prodstate(astate,bstate) + nprodbasis * (jbspin-1)
     ket(:) = wpol_static%residu_left(kbf,:)

     wtmp =  SUM( 2.0_dp * bra(:)*ket(:)/(-wpol_static%pole(:)) )   ! Factor two comes from Resonant and Anti-resonant transitions
     apb_matrix(t_ia,t_jb) =  apb_matrix(t_ia,t_jb) - wtmp
     amb_matrix(t_ia,t_jb) =  amb_matrix(t_ia,t_jb) - wtmp

     kbf = index_prodstate(istate,bstate) + nprodbasis * (iaspin-1)
     bra(:) = wpol_static%residu_left(kbf,:)
     kbf = index_prodstate(astate,jstate) + nprodbasis * (jbspin-1)
     ket(:) = wpol_static%residu_left(kbf,:)  

     wtmp =  SUM( 2.0_dp * bra(:)*ket(:)/(-wpol_static%pole(:)) )   ! Factor two comes from Resonant and Anti-resonant transitions
     apb_matrix(t_ia,t_jb) =  apb_matrix(t_ia,t_jb) - wtmp
     amb_matrix(t_ia,t_jb) =  amb_matrix(t_ia,t_jb) + wtmp


   enddo

 enddo

 deallocate(bra,ket)


 call stop_clock(timing_build_bse)


end subroutine build_amb_apb_bse


!=========================================================================
subroutine build_amb_apb_screened_exchange_auxil(alpha_local,desc_apb,wpol,wpol_static,m_apb,n_apb,amb_matrix,apb_matrix)
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_spectral_function
 use m_basis_set
 use m_eri_ao_mo
 use m_tools 
 implicit none

 real(dp),intent(in)                :: alpha_local
 integer,intent(in)                 :: desc_apb(ndel)
 type(spectral_function),intent(in) :: wpol,wpol_static
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(inout)             :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
!=====
 logical              :: is_bse
 integer              :: nmat
 integer              :: t_ia,t_jb,t_ia_global,t_jb_global
 integer              :: istate,astate,jstate,bstate
 integer              :: jstate_prev
 integer              :: iaspin,jbspin
 integer              :: kbf
 real(dp)             :: wtmp
 integer              :: jstate_min,jstate_max
 integer              :: ipole,ibf_auxil,jbf_auxil,ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable :: vsqrt_chi_vsqrt(:,:)
 real(dp),allocatable :: vsqrt_chi_vsqrt_i(:),residu_i(:),wp0_i(:,:)
 real(dp),allocatable :: wp0(:,:,:,:),w0_local(:)
 integer              :: iprow,ipcol,irank
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: amb_block(:,:)
 real(dp),allocatable :: apb_block(:,:)
!=====

 call start_clock(timing_build_bse)
 if( .NOT. has_auxil_basis ) call die('Does not have auxil basis. This should not happen')

 write(stdout,'(a)')       ' Build W part Auxil' 
 write(stdout,'(a,f8.3)') ' Content of Exchange: ',alpha_local

 nmat = desc_apb(3)
 ! Is it an exchange or a screened exchange calculation
 is_bse = ALLOCATED(wpol_static%w0) .OR. ALLOCATED(wpol_static%residu_left)


 !
 ! Distribution over the "ortho" parallelization direction
 ! 
 jstate_min = ncore_W+1
 jstate_max = MAXVAL( wpol%transition_table_apb(1,1:wpol%npole_reso_apb) )
 do irank=0,rank_ortho
   if( irank > 0 ) jstate_min = jstate_max + 1
   jstate_max = MAXVAL( wpol%transition_table_apb(1,1:wpol%npole_reso_apb) )
   jstate_max = MIN( jstate_min + (jstate_max-jstate_min+1) / (nproc_ortho-irank) - 1 , jstate_max )
 enddo

 call clean_allocate('Temporary array for W',wp0,1,nauxil_3center,ncore_W+1,nvirtual_W-1,jstate_min,jstate_max,1,nspin)
 wp0(:,:,:,:) = 0.0_dp


 if( is_bse ) then
#ifndef HAVE_SCALAPACK
   do iaspin=1,nspin
  
     allocate(vsqrt_chi_vsqrt(nauxil_2center,nauxil_2center))
  
  
     !
     ! Test if w0 is already available or if we need to calculate it first
     if( ALLOCATED(wpol_static%w0) ) then
  
       vsqrt_chi_vsqrt(:,:) = wpol_static%w0(:,:)
  
     else if( ALLOCATED(wpol_static%residu_left) ) then

       vsqrt_chi_vsqrt(:,:) = 0.0_dp
       do ipole=1,wpol_static%npole_reso
         do jbf_auxil=1,nauxil_2center
           vsqrt_chi_vsqrt(:,jbf_auxil) = vsqrt_chi_vsqrt(:,jbf_auxil) &
                    - wpol_static%residu_left(:,ipole)                 &
                      * wpol_static%residu_left(jbf_auxil,ipole) * 2.0_dp / wpol_static%pole(ipole)
         enddo
       enddo
  
     endif
      
     !
     ! The last index of wp0 only runs on occupied states (to save memory and CPU time)
     ! Be careful not to forget it in the following 
     do jstate=jstate_min,jstate_max
       wp0(:,ncore_W+1:nvirtual_W-1,jstate,iaspin) = MATMUL( vsqrt_chi_vsqrt(:,:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,jstate,iaspin) )
     enddo
    
     deallocate(vsqrt_chi_vsqrt)
  
   enddo

#else

   do iaspin=1,nspin
  
     !
     ! Test if w0 is already available or if we need to calculate it first
     if( ALLOCATED(wpol_static%w0) ) then
  
       allocate(wp0_i(ncore_W+1:nvirtual_W-1,jstate_min:jstate_max))
       allocate(w0_local(nauxil_3center))
  
       do ibf_auxil_global=1,nauxil_2center
  
         do jbf_auxil=1,nauxil_3center
           jbf_auxil_global = ibf_auxil_g(jbf_auxil)
           w0_local(jbf_auxil) = wpol_static%w0(ibf_auxil_global,jbf_auxil_global)
         enddo
  
         do jstate=jstate_min,jstate_max
           wp0_i(ncore_W+1:nvirtual_W-1,jstate) = MATMUL( w0_local(:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,jstate,iaspin) )
         enddo
         call xsum_auxil(wp0_i)
  
         if( iproc_ibf_auxil(ibf_auxil_global) == rank_auxil ) then
           wp0(ibf_auxil_l(ibf_auxil_global),:,:,iaspin) = wp0_i(:,:)
         endif
  
       enddo
       deallocate(w0_local)
  
     else if( ALLOCATED(wpol_static%residu_left) ) then
  
       allocate(vsqrt_chi_vsqrt_i(nauxil_3center))
       allocate(residu_i(wpol_static%npole_reso))
       allocate(wp0_i(ncore_W+1:nvirtual_W-1,jstate_min:jstate_max))
      
       do ibf_auxil=1,nauxil_2center
      
         if( iproc_ibf_auxil(ibf_auxil) == rank_auxil ) then
           residu_i(:) = wpol_static%residu_left(ibf_auxil_l(ibf_auxil),:)
         else
           residu_i(:) = 0.0_dp
         endif
         call xsum_auxil(residu_i)
      
         vsqrt_chi_vsqrt_i(:) = 0.0_dp
         do ipole=1,wpol_static%npole_reso
           vsqrt_chi_vsqrt_i(:) = vsqrt_chi_vsqrt_i(:) &
                    - residu_i(ipole) * wpol_static%residu_left(:,ipole) * 2.0_dp / wpol_static%pole(ipole)
         enddo
         !
         ! The last index of wp0 only runs on occupied states (to save memory and CPU time)
         ! Be careful in the following not to forget it
         do jstate=jstate_min,jstate_max
           wp0_i(ncore_W+1:nvirtual_W-1,jstate) = MATMUL( vsqrt_chi_vsqrt_i(:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,jstate,iaspin) )
         enddo
         call xsum_auxil(wp0_i)
      
         if( iproc_ibf_auxil(ibf_auxil) == rank_auxil ) then
           wp0(ibf_auxil_l(ibf_auxil),:,:,iaspin) = wp0_i(:,:)
         endif
      
       enddo
      
       deallocate(vsqrt_chi_vsqrt_i,residu_i,wp0_i)
  
     endif
  
   enddo
#endif
 endif


 ! 
 ! Add the exact exchange here
 if( alpha_local > 1.0e-6_dp ) then

   do iaspin=1,nspin
     do jstate=jstate_min,jstate_max
       wp0(:,ncore_W+1:nvirtual_W-1,jstate,iaspin) = wp0(:,ncore_W+1:nvirtual_W-1,jstate,iaspin) &
                          + alpha_local *  eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,jstate,iaspin)
     enddo
   enddo

 endif
 

 if( nprow_sd * npcol_sd > 1 ) &
    write(stdout,'(a,i4,a,i4)') ' SCALAPACK grid    :',nprow_sd,' x ',npcol_sd
 
 ! First loops over the SCALAPACK grid
 do ipcol=0,npcol_sd-1
   do iprow=0,nprow_sd-1
     m_apb_block = row_block_size(nmat,iprow,nprow_sd)
     n_apb_block = col_block_size(nmat,ipcol,npcol_sd)

     allocate(amb_block(m_apb_block,n_apb_block))
     allocate(apb_block(m_apb_block,n_apb_block))
     apb_block(:,:) = 0.0_dp
     amb_block(:,:) = 0.0_dp


     ! Set up -W contributions to matrices (A+B) and (A-B)
     !
     do t_jb=1,n_apb_block
       t_jb_global = colindex_local_to_global(ipcol,npcol_sd,t_jb)
       jstate = wpol%transition_table_apb(1,t_jb_global)
       bstate = wpol%transition_table_apb(2,t_jb_global)
       jbspin = wpol%transition_table_apb(3,t_jb_global)

!       if( MODULO( jstate, nproc_ortho ) /= rank_ortho ) cycle
       if( jstate < jstate_min .OR. jstate > jstate_max ) cycle

       do t_ia=1,m_apb_block
         t_ia_global = rowindex_local_to_global(iprow,nprow_sd,t_ia)
         istate = wpol%transition_table_apb(1,t_ia_global)
         astate = wpol%transition_table_apb(2,t_ia_global)
         iaspin = wpol%transition_table_apb(3,t_ia_global)

         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ia_global < t_jb_global ) cycle

         if( iaspin /= jbspin ) cycle

         wtmp = DOT_PRODUCT( eri_3center_eigen(:,astate,bstate,iaspin) , wp0(:,istate,jstate,iaspin) )

         apb_block(t_ia,t_jb) = -wtmp
         amb_block(t_ia,t_jb) = -wtmp


         wtmp = DOT_PRODUCT( eri_3center_eigen(:,istate,bstate,iaspin) , wp0(:,astate,jstate,iaspin) )

         apb_block(t_ia,t_jb) =  apb_block(t_ia,t_jb) - wtmp
         amb_block(t_ia,t_jb) =  amb_block(t_ia,t_jb) + wtmp


       enddo

     enddo


     call xsum_world(amb_block)
     call xsum_world(apb_block)

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       amb_matrix(:,:) = amb_matrix(:,:) + amb_block(:,:)
       apb_matrix(:,:) = apb_matrix(:,:) + apb_block(:,:)
     endif
     deallocate(amb_block)
     deallocate(apb_block)


   enddo
 enddo

 call clean_deallocate('Temporary array for W',wp0)


 call stop_clock(timing_build_bse)


end subroutine build_amb_apb_screened_exchange_auxil


!=========================================================================
