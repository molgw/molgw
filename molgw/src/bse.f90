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
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
 integer              :: klmin
 real(dp),allocatable :: eri_eigenstate_klmin(:,:,:,:)
 real(dp)             :: eri_eigen_ijkl
 real(dp)             :: eri_eigen_ikjl,eri_eigen_iljk
 logical              :: k_is_klmin
 integer              :: iprow,ipcol
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: amb_block(:,:)
 real(dp),allocatable :: apb_block(:,:)
!=====

 call start_clock(timing_build_common)

 write(stdout,'(a)') ' Build Common part: Energies + Hartree + possibly Exchange'
 write(stdout,'(a,f8.3)') ' Content of Exchange: ',alpha_local

 if( .NOT. has_auxil_basis) then
   allocate(eri_eigenstate_klmin(nstate,nstate,nstate,nspin))
   ! Set this to zero and then enforce the calculation of the first series of
   ! Coulomb integrals
   eri_eigenstate_klmin(:,:,:,:) = 0.0_dp
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
     do t_kl=1,n_apb_block
       t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)
       kstate = wpol%transition_table_apb(1,t_kl_global)
       lstate = wpol%transition_table_apb(2,t_kl_global)
       klspin = wpol%transition_table_apb(3,t_kl_global)
    
       if( .NOT. has_auxil_basis ) then
         klmin = MIN(kstate,lstate)
         k_is_klmin = (klmin == kstate)
         call calculate_eri_4center_eigen(nbf,nstate,c_matrix,klmin,klspin,eri_eigenstate_klmin)
       endif
    
       do t_ij=1,m_apb_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)
    
         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ij_global < t_kl_global ) cycle
    
         if(has_auxil_basis) then
           eri_eigen_ijkl = eri_eigen_ri_paral(istate,jstate,ijspin,kstate,lstate,klspin)
         else
           if(k_is_klmin) then ! treating (k,l)
             eri_eigen_ijkl = eri_eigenstate_klmin(lstate,istate,jstate,ijspin)
           else                ! treating (l,k)
             eri_eigen_ijkl = eri_eigenstate_klmin(kstate,istate,jstate,ijspin)
           endif
         endif
    
         if( .NOT. is_triplet) then
           apb_block(t_ij,t_kl) = 2.0_dp * eri_eigen_ijkl * spin_fact
         else
           apb_block(t_ij,t_kl) = 0.0_dp
         endif
         amb_block(t_ij,t_kl) = 0.0_dp
    
         if( alpha_local > 1.0e-6_dp ) then
           if(ijspin==klspin) then
             if(has_auxil_basis) then
               eri_eigen_ikjl = eri_eigen_ri_paral(istate,kstate,ijspin,jstate,lstate,klspin)
               eri_eigen_iljk = eri_eigen_ri_paral(istate,lstate,ijspin,jstate,kstate,klspin)
             else
               if(k_is_klmin) then
                 eri_eigen_ikjl = eri_eigenstate_klmin(istate,jstate,lstate,klspin)
                 eri_eigen_iljk = eri_eigenstate_klmin(jstate,istate,lstate,klspin)
               else
                 eri_eigen_ikjl = eri_eigenstate_klmin(jstate,istate,kstate,klspin)
                 eri_eigen_iljk = eri_eigenstate_klmin(istate,jstate,kstate,klspin)
               endif
             endif
             apb_block(t_ij,t_kl) = apb_block(t_ij,t_kl) - eri_eigen_ikjl * alpha_local - eri_eigen_iljk * alpha_local
             amb_block(t_ij,t_kl) = amb_block(t_ij,t_kl) - eri_eigen_ikjl * alpha_local + eri_eigen_iljk * alpha_local
           endif
         endif
    
         if( t_ij_global == t_kl_global ) then
           !
           ! Only one proc should add the diagonal
           if( rank == 0 ) then
             apb_block(t_ij,t_kl) =  apb_block(t_ij,t_kl) + ( energy(lstate,klspin) - energy(kstate,klspin) )
             amb_block(t_ij,t_kl) =  amb_block(t_ij,t_kl) + ( energy(lstate,klspin) - energy(kstate,klspin) )
           endif
           ! First part of the RPA correlation energy: sum over diagonal terms
           rpa_correlation = rpa_correlation - 0.25_dp * ( apb_block(t_ij,t_kl) + amb_block(t_ij,t_kl) )
         endif
    
       enddo 
    
     enddo 


#ifdef HAVE_SCALAPACK
     call DGSUM2D(desc_apb(2),'A',' ',m_apb_block,n_apb_block,amb_block,m_apb_block,iprow,ipcol)
     call DGSUM2D(desc_apb(2),'A',' ',m_apb_block,n_apb_block,apb_block,m_apb_block,iprow,ipcol)
#endif

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       amb_matrix(:,:) = amb_block(:,:)
       apb_matrix(:,:) = apb_block(:,:)
     endif


     deallocate(amb_block)
     deallocate(apb_block)
   enddo 
 enddo 

#ifdef HAVE_SCALAPACK
 call xsum(rpa_correlation)
#endif

 !
 ! Set up the diagonal of A-B in the RPA approximation
 ! 
 do t_ij_global=1,nmat
   istate = wpol%transition_table_apb(1,t_ij_global)
   jstate = wpol%transition_table_apb(2,t_ij_global)
   ijspin = wpol%transition_table_apb(3,t_ij_global)

   amb_diag_rpa(t_ij_global) = energy(jstate,ijspin) - energy(istate,ijspin)

 enddo


 if(ALLOCATED(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)


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
 integer              :: t_ij,t_kl,t_kl_global
 integer              :: kstate,lstate
 integer              :: klspin
!=====

 write(stdout,'(a)') ' Build diagonal part with auxil basis: Energies'

 do t_kl_global=1,nmat
   t_ij = rowindex_global_to_local('S',t_kl_global)
   t_kl = colindex_global_to_local('S',t_kl_global)

   kstate = wpol%transition_table_apb(1,t_kl_global)
   lstate = wpol%transition_table_apb(2,t_kl_global)
   klspin = wpol%transition_table_apb(3,t_kl_global)
   amb_diag_rpa(t_kl_global) = energy(lstate,klspin) - energy(kstate,klspin)

   ! If the diagonal element belongs to this proc, then add it.
   if( t_ij > 0 .AND. t_kl > 0 ) then
     apb_matrix(t_ij,t_kl) =  amb_diag_rpa(t_kl_global)
     amb_matrix(t_ij,t_kl) =  amb_diag_rpa(t_kl_global)
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
 integer              :: t_ij,t_kl,t_kl_global
!=====

 rpa_correlation = 0.0_dp

 do t_kl_global=1,nmat
   t_ij = rowindex_global_to_local('S',t_kl_global)
   t_kl = colindex_global_to_local('S',t_kl_global)

   ! If the diagonal element belongs to this proc, then add it.
   if( t_ij > 0 .AND. t_kl > 0 ) then
     rpa_correlation = rpa_correlation - 0.25_dp * apb_matrix(t_ij,t_kl)   &
                                       - 0.25_dp * amb_matrix(t_ij,t_kl) 
   endif
 enddo

 call xsum(rpa_correlation)


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
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
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


       do t_kl=1,n_apb_block
         t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)
         kstate = wpol%transition_table_apb(1,t_kl_global)
         lstate = wpol%transition_table_apb(2,t_kl_global)
         klspin = wpol%transition_table_apb(3,t_kl_global)

         eri_3center_right(t_kl) = eri_3center_eigen(ibf_auxil,kstate,lstate,klspin)
  
       enddo
    
       do t_ij=1,m_apb_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)
    
         eri_3center_left(t_ij) = eri_3center_eigen(ibf_auxil,istate,jstate,ijspin)

       enddo

       call DGER(m_apb_block,n_apb_block,2.0_dp*spin_fact,eri_3center_left,1,eri_3center_right,1,apb_block,m_apb_block)

     enddo
     
     deallocate(eri_3center_left,eri_3center_right)

#ifdef HAVE_SCALAPACK
     call DGSUM2D(desc_apb(2),'A',' ',m_apb_block,n_apb_block,apb_block,m_apb_block,iprow,ipcol)
#endif

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
 integer              :: t_kl,t_kl_global
 integer              :: kstate,lstate
 integer              :: klspin
 integer              :: klmin
 real(dp),allocatable :: eri_eigenstate_klmin(:,:,:,:)
 real(dp)             :: eri_eigen_klkl
 logical              :: k_is_klmin
 real(dp),parameter   :: empirical_fact=1.5_dp
 character(len=100)   :: ctmp
!=====

 call start_clock(timing_build_common)

 write(stdout,'(a)') ' Build diagonal of the RPA part: Energies + Hartree'
 if( ABS( empirical_fact - 1.0_dp ) > 1.0e-6_dp ) then
   write(ctmp,'(a,1x,f6.2)') 'Empirical parameter',empirical_fact
   call issue_warning(ctmp)
 endif

 if( .NOT. has_auxil_basis) then
   allocate(eri_eigenstate_klmin(nstate,nstate,nstate,nspin))
   ! Set this to zero and then enforce the calculation of the first series of
   ! Coulomb integrals
   eri_eigenstate_klmin(:,:,:,:) = 0.0_dp
 endif

 !
 ! Set up energy+hartree+exchange contributions to matrices (A+B) and (A-B) 
 !


 ! Then loop inside each of the SCALAPACK blocks
 do t_kl=1,wpol%npole_reso_spa
   kstate = wpol%transition_table_spa(1,t_kl)
   lstate = wpol%transition_table_spa(2,t_kl)
   klspin = wpol%transition_table_spa(3,t_kl)

   if( .NOT. has_auxil_basis ) then
     klmin = MIN(kstate,lstate)
     k_is_klmin = (klmin == kstate)
     call calculate_eri_4center_eigen(nbf,nstate,c_matrix,klmin,klspin,eri_eigenstate_klmin)
   endif

   if(has_auxil_basis) then
     eri_eigen_klkl = eri_eigen_ri(kstate,lstate,klspin,kstate,lstate,klspin)
   else
     if(k_is_klmin) then ! treating (k,l)
       eri_eigen_klkl = eri_eigenstate_klmin(lstate,kstate,lstate,klspin)
     else                ! treating (l,k)
       eri_eigen_klkl = eri_eigenstate_klmin(kstate,kstate,lstate,klspin)
     endif
   endif

   a_diag(t_kl) = eri_eigen_klkl * spin_fact + energy(lstate,klspin) - energy(kstate,klspin)
   a_diag(t_kl) = a_diag(t_kl) * empirical_fact

 enddo 

 if(ALLOCATED(eri_eigenstate_klmin)) deallocate(eri_eigenstate_klmin)


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
 implicit none

 integer,intent(in)                 :: nmat,nstate
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)                :: occupation(nstate,nspin)
 type(spectral_function),intent(in) :: wpol
 integer,intent(in)                 :: m_apb,n_apb
 real(dp),intent(inout)             :: apb_matrix(m_apb,n_apb)
!=====
 integer              :: nspin_tddft
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
 logical              :: require_gradient
 real(dp),allocatable :: grad_ij(:,:,:),grad_kl(:,:,:)
 real(dp),allocatable :: dot_ij_kl(:,:),dot_rho_ij(:,:),dot_rho_kl(:,:)
 real(dp),allocatable :: v2rho2(:,:)
 real(dp),allocatable :: vsigma(:,:)
 real(dp),allocatable :: v2rhosigma(:,:)
 real(dp),allocatable :: v2sigma2(:,:)
 real(dp),allocatable :: wf_r(:,:,:)
 real(dp),allocatable :: wf_gradr(:,:,:,:)
 real(dp),allocatable :: rho_gradr(:,:,:)
 real(dp)             :: xctmp
 integer              :: ipcol,iprow
 integer              :: m_apb_block,n_apb_block
 real(dp),allocatable :: apb_block(:,:)
!=====
 interface 
   subroutine prepare_tddft(nspin_tddft,nstate,basis,c_matrix,occupation,v2rho2,vsigma,v2rhosigma,v2sigma2,wf_r,wf_gradr,rho_gradr)
     use,intrinsic ::  iso_c_binding, only: C_INT,C_DOUBLE
     use m_definitions
     use m_timing
     use m_warning
     use m_memory
     use m_inputparam
     use m_mpi
     use m_dft_grid
     use m_basis_set
     use m_hamiltonian,only: setup_density_matrix
#ifdef HAVE_LIBXC
     use libxc_funcs_m
     use xc_f90_lib_m
     use xc_f90_types_m
#endif
     integer,intent(in)               :: nspin_tddft,nstate
     type(basis_set),intent(in)       :: basis
     real(dp),intent(in)              :: c_matrix(basis%nbf,nstate,nspin)
     real(dp),intent(in)              :: occupation(nstate,nspin)
     real(dp),allocatable,intent(out) :: v2rho2(:,:)
     real(dp),allocatable,intent(out) :: vsigma(:,:)
     real(dp),allocatable,intent(out) :: v2rhosigma(:,:)
     real(dp),allocatable,intent(out) :: v2sigma2(:,:)
     real(dp),allocatable,intent(out) :: wf_r(:,:,:)
     real(dp),allocatable,intent(out) :: wf_gradr(:,:,:,:)
     real(dp),allocatable,intent(out) :: rho_gradr(:,:,:)
   end subroutine prepare_tddft
end interface 
!=====

 call start_clock(timing_build_tddft)

 write(stdout,'(a)') ' Build fxc part'

 if( is_triplet ) then
   nspin_tddft = 2
 else
   nspin_tddft = nspin
 endif
 !
 ! Prepare TDDFT calculations
 call prepare_tddft(nspin_tddft,nstate,basis,c_matrix,occupation,v2rho2,vsigma,v2rhosigma,v2sigma2,wf_r,wf_gradr,rho_gradr)
 require_gradient = .FALSE.
 if(ALLOCATED(v2sigma2)) then ! GGA case
   require_gradient = .TRUE.
   allocate(grad_ij(3,ngrid,nspin))
   allocate(grad_kl(3,ngrid,nspin))
   allocate(dot_ij_kl(ngrid,nspin))
   allocate(dot_rho_ij(ngrid,nspin))
   allocate(dot_rho_kl(ngrid,nspin))
 endif


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
     do t_kl=1,n_apb_block
       t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)
       kstate = wpol%transition_table_apb(1,t_kl_global)
       lstate = wpol%transition_table_apb(2,t_kl_global)
       klspin = wpol%transition_table_apb(3,t_kl_global)
    
       do t_ij=1,m_apb_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)
    
         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ij_global < t_kl_global ) cycle
    
    
         if( nspin_tddft == 1 ) then
    
           xctmp = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                              * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                              * v2rho2(:,ijspin) * 2.0_dp )
    
    
           if(require_gradient) then
    
             grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
             dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                                  + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
             dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                                   + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
             dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                                   + rho_gradr(3,:,1) * grad_kl(3,:,klspin)
    
             xctmp = xctmp   &
                   +  SUM( dot_ij_kl(:,1) * 4.0_dp * vsigma(:,1) ) &
                   +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * 8.0_dp * v2sigma2(:,1) ) &
                   +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                           + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                             * 4.0_dp * v2rhosigma(:,1) )
    
           endif
    
         else if( .NOT. is_triplet ) then
    
           xctmp = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                              * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                              * ( v2rho2(:,1) + v2rho2(:,2) ) )
    
    
           if(require_gradient) then
    
             grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
             dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                                  + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
             dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                                   + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
             dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                                   + rho_gradr(3,:,1) * grad_kl(3,:,klspin)
    
             xctmp = xctmp   &
                   +  SUM( dot_ij_kl(:,1) * ( 2.0_dp * vsigma(:,1) + vsigma(:,2) ) ) &
                   +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * ( v2sigma2(:,1) + 0.5_dp * v2sigma2(:,2) +  2.0_dp * v2sigma2(:,3) + v2sigma2(:,4) ) ) &
                   +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                          + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                             * ( v2rhosigma(:,1) + v2rhosigma(:,2) + v2rhosigma(:,3) )  )
           endif
    
    
           
         else ! triplet case
    
           xctmp = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                              * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                              * ( v2rho2(:,1) - v2rho2(:,2) ) )
    
           if(require_gradient) then
    
             grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
             grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
             grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
             dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                                  + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
             dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                                   + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
             dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                                   + rho_gradr(3,:,1) * grad_kl(3,:,klspin)
    
             xctmp = xctmp   &
                   +  SUM( dot_ij_kl(:,1) * ( 2.0_dp * vsigma(:,1) - vsigma(:,2) ) ) &
                   +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * ( v2sigma2(:,1) - v2sigma2(:,3) ) ) &   !FIXME 3 or 5 are working, but only one is correct in principle
                   +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                          + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                             * ( v2rhosigma(:,1) - v2rhosigma(:,4) )  )   !FIXME 3 and 4 are working, but only one is correct in principle
           endif
    
    
    
         endif
    
         ! The factor two accounts for (A+B), and not A or B.
         apb_block(t_ij,t_kl) = apb_block(t_ij,t_kl) + 2.0_dp * xctmp
    
       enddo
     enddo

     call xsum(apb_block)

     if( iprow == iprow_sd .AND. ipcol == ipcol_sd ) then
       apb_matrix(:,:) = apb_matrix(:,:) + apb_block(:,:)
     endif
     deallocate(apb_block)

   enddo
 enddo


 deallocate(wf_r)

 if(require_gradient) then
   deallocate(grad_ij,grad_kl)
   deallocate(dot_ij_kl,dot_rho_ij,dot_rho_kl)
   deallocate(v2rho2)
   deallocate(vsigma)
   deallocate(v2rhosigma)
   deallocate(v2sigma2)
   deallocate(wf_gradr)
   deallocate(rho_gradr)
 endif


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
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: ijspin,klspin
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
 do t_kl=1,n_apb
   t_kl_global = colindex_local_to_global('C',t_kl)
   kstate = wpol%transition_table_apb(1,t_kl_global)
   lstate = wpol%transition_table_apb(2,t_kl_global)
   klspin = wpol%transition_table_apb(3,t_kl_global)

   do t_ij=1,m_apb
     t_ij_global = rowindex_local_to_global('C',t_ij)
     istate = wpol%transition_table_apb(1,t_ij_global)
     jstate = wpol%transition_table_apb(2,t_ij_global)
     ijspin = wpol%transition_table_apb(3,t_ij_global)

     !
     ! Only calculate the lower triangle
     ! Symmetrization will be performed later (in the diago subroutines)
     if( t_ij_global < t_kl_global ) cycle

     if(ijspin/=klspin) cycle

     kbf = index_prodstate(istate,kstate) + nprodbasis * (ijspin-1)
     bra(:) = wpol_static%residu_left(kbf,:)
     kbf = index_prodstate(jstate,lstate) + nprodbasis * (klspin-1)
     ket(:) = wpol_static%residu_left(kbf,:)

     wtmp =  SUM( 2.0_dp * bra(:)*ket(:)/(-wpol_static%pole(:)) )   ! Factor two comes from Resonant and Anti-resonant transitions
     apb_matrix(t_ij,t_kl) =  apb_matrix(t_ij,t_kl) - wtmp
     amb_matrix(t_ij,t_kl) =  amb_matrix(t_ij,t_kl) - wtmp

     kbf = index_prodstate(istate,lstate) + nprodbasis * (ijspin-1)
     bra(:) = wpol_static%residu_left(kbf,:)
     kbf = index_prodstate(jstate,kstate) + nprodbasis * (klspin-1)
     ket(:) = wpol_static%residu_left(kbf,:)  

     wtmp =  SUM( 2.0_dp * bra(:)*ket(:)/(-wpol_static%pole(:)) )   ! Factor two comes from Resonant and Anti-resonant transitions
     apb_matrix(t_ij,t_kl) =  apb_matrix(t_ij,t_kl) - wtmp
     amb_matrix(t_ij,t_kl) =  amb_matrix(t_ij,t_kl) + wtmp


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
 integer              :: t_ij,t_kl,t_ij_global,t_kl_global
 integer              :: istate,jstate,kstate,lstate
 integer              :: kstate_prev
 integer              :: ijspin,klspin
 integer              :: kbf
 real(dp)             :: wtmp
 integer              :: kstate_max
 integer              :: ipole,ibf_auxil,jbf_auxil,ibf_auxil_global,jbf_auxil_global
 real(dp),allocatable :: vsqrt_chi_vsqrt(:,:)
 real(dp),allocatable :: vsqrt_chi_vsqrt_i(:),residu_i(:),wp0_i(:,:)
 real(dp),allocatable :: wp0(:,:,:,:),w0_local(:)
 integer              :: iprow,ipcol
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


 kstate_max = MAXVAL( wpol%transition_table_apb(1,1:wpol%npole_reso_apb) )

 call clean_allocate('Temporary array for W',wp0,1,nauxil_3center,ncore_W+1,nvirtual_W-1,ncore_W+1,kstate_max,1,nspin)
 wp0(:,:,:,:) = 0.0_dp


 if( is_bse ) then
#ifndef HAVE_SCALAPACK
   do ijspin=1,nspin
  
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
     do kstate=ncore_W+1,kstate_max
       wp0(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) = MATMUL( vsqrt_chi_vsqrt(:,:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) )
     enddo
    
     deallocate(vsqrt_chi_vsqrt)
  
   enddo

#else

   do ijspin=1,nspin
  
     !
     ! Test if w0 is already available or if we need to calculate it first
     if( ALLOCATED(wpol_static%w0) ) then
  
       allocate(wp0_i(ncore_W+1:nvirtual_W-1,ncore_W+1:kstate_max))
       allocate(w0_local(nauxil_3center))
  
       do ibf_auxil_global=1,nauxil_2center
  
         do jbf_auxil=1,nauxil_3center
           jbf_auxil_global = ibf_auxil_g(jbf_auxil)
           w0_local(jbf_auxil) = wpol_static%w0(ibf_auxil_global,jbf_auxil_global)
         enddo
  
         do kstate=ncore_W+1,kstate_max
           wp0_i(ncore_W+1:nvirtual_W-1,kstate) = MATMUL( w0_local(:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) )
         enddo
         call xsum(wp0_i)
  
         if( iproc_ibf_auxil(ibf_auxil_global) == rank ) then
           wp0(ibf_auxil_l(ibf_auxil_global),:,:,ijspin) = wp0_i(:,:)
         endif
  
       enddo
       deallocate(w0_local)
  
     else if( ALLOCATED(wpol_static%residu_left) ) then
  
       allocate(vsqrt_chi_vsqrt_i(nauxil_3center))
       allocate(residu_i(wpol_static%npole_reso))
       allocate(wp0_i(ncore_W+1:nvirtual_W-1,ncore_W+1:kstate_max))
      
       do ibf_auxil=1,nauxil_2center
      
         if( iproc_ibf_auxil(ibf_auxil) == rank ) then
           residu_i(:) = wpol_static%residu_left(ibf_auxil_l(ibf_auxil),:)
         else
           residu_i(:) = 0.0_dp
         endif
         call xsum(residu_i)
      
         vsqrt_chi_vsqrt_i(:) = 0.0_dp
         do ipole=1,wpol_static%npole_reso
           vsqrt_chi_vsqrt_i(:) = vsqrt_chi_vsqrt_i(:) &
                    - residu_i(ipole) * wpol_static%residu_left(:,ipole) * 2.0_dp / wpol_static%pole(ipole)
         enddo
         !
         ! The last index of wp0 only runs on occupied states (to save memory and CPU time)
         ! Be careful in the following not to forget it
         do kstate=ncore_W+1,kstate_max
           wp0_i(ncore_W+1:nvirtual_W-1,kstate) = MATMUL( vsqrt_chi_vsqrt_i(:) , eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) )
         enddo
         call xsum(wp0_i)
      
         if( iproc_ibf_auxil(ibf_auxil) == rank ) then
           wp0(ibf_auxil_l(ibf_auxil),:,:,ijspin) = wp0_i(:,:)
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

   do ijspin=1,nspin
     do kstate=ncore_W+1,kstate_max
       wp0(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) = wp0(:,ncore_W+1:nvirtual_W-1,kstate,ijspin) &
                          + alpha_local *  eri_3center_eigen(:,ncore_W+1:nvirtual_W-1,kstate,ijspin)
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


     ! Set up -W contributions to matrices (A+B) and (A-B)
     !
     do t_kl=1,n_apb_block
       t_kl_global = colindex_local_to_global(ipcol,npcol_sd,t_kl)
       kstate = wpol%transition_table_apb(1,t_kl_global)
       lstate = wpol%transition_table_apb(2,t_kl_global)
       klspin = wpol%transition_table_apb(3,t_kl_global)

       do t_ij=1,m_apb_block
         t_ij_global = rowindex_local_to_global(iprow,nprow_sd,t_ij)
         istate = wpol%transition_table_apb(1,t_ij_global)
         jstate = wpol%transition_table_apb(2,t_ij_global)
         ijspin = wpol%transition_table_apb(3,t_ij_global)

         !
         ! Only calculate the lower triangle
         ! Symmetrization will be performed later (in the diago subroutines)
         if( t_ij_global < t_kl_global ) cycle

         if(ijspin/=klspin) cycle

         wtmp = DOT_PRODUCT( eri_3center_eigen(:,jstate,lstate,ijspin) , wp0(:,istate,kstate,ijspin) )

         apb_block(t_ij,t_kl) = -wtmp
         amb_block(t_ij,t_kl) = -wtmp


         wtmp = DOT_PRODUCT( eri_3center_eigen(:,istate,lstate,ijspin) , wp0(:,jstate,kstate,ijspin) )

         apb_block(t_ij,t_kl) =  apb_block(t_ij,t_kl) - wtmp
         amb_block(t_ij,t_kl) =  amb_block(t_ij,t_kl) + wtmp


       enddo

     enddo


#ifdef HAVE_SCALAPACK
     call DGSUM2D(desc_apb(2),'A',' ',m_apb_block,n_apb_block,amb_block,m_apb_block,iprow,ipcol)
     call DGSUM2D(desc_apb(2),'A',' ',m_apb_block,n_apb_block,apb_block,m_apb_block,iprow,ipcol)
#endif

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
subroutine prepare_tddft(nspin_tddft,nstate,basis,c_matrix,occupation,v2rho2,vsigma,v2rhosigma,v2sigma2,wf_r,wf_gradr,rho_gradr)
 use,intrinsic ::  iso_c_binding, only: C_INT,C_DOUBLE
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam
 use m_mpi
 use m_dft_grid
 use m_basis_set
 use m_hamiltonian,only: setup_density_matrix
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 integer,intent(in)               :: nspin_tddft,nstate
 type(basis_set),intent(in)       :: basis
 real(dp),intent(in)              :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)              :: occupation(nstate,nspin)
 real(dp),allocatable,intent(out) :: v2rho2(:,:)
 real(dp),allocatable,intent(out) :: vsigma(:,:)
 real(dp),allocatable,intent(out) :: v2rhosigma(:,:)
 real(dp),allocatable,intent(out) :: v2sigma2(:,:)
 real(dp),allocatable,intent(out) :: wf_r(:,:,:)
 real(dp),allocatable,intent(out) :: wf_gradr(:,:,:,:)
 real(dp),allocatable,intent(out) :: rho_gradr(:,:,:)
!=====
 real(dp),parameter :: kernel_capping=1.0e14_dp
 integer  :: idft_xc,igrid
 integer  :: ispin
 real(dp) :: rr(3)
 real(dp) :: basis_function_r(basis%nbf)
 real(dp) :: basis_function_gradr(3,basis%nbf)
 real(dp) :: rhor_r(nspin)
 real(dp) :: grad_rhor(3,nspin)
 real(dp) :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp) :: max_v2sigma2
 logical  :: require_gradient
 character(len=256) :: string
#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func(ndft_xc),xc_functest
 type(xc_f90_pointer_t) :: xc_info(ndft_xc),xc_infotest
#endif
 real(C_DOUBLE) :: rho_c(nspin_tddft)
 real(C_DOUBLE) :: v2rho2_c(2*nspin_tddft-1)
 real(C_DOUBLE) :: sigma_c(2*nspin_tddft-1)
 real(C_DOUBLE) :: vrho_c(nspin_tddft)
 real(C_DOUBLE) :: vsigma_c(2*nspin_tddft-1)
 real(C_DOUBLE) :: v2rhosigma_c(5*nspin_tddft-4)
 real(C_DOUBLE) :: v2sigma2_c(5*nspin_tddft-4)
!=====

 !
 ! Prepare DFT kernel calculation with Libxc
 !
 require_gradient  =.FALSE.
 do idft_xc=1,ndft_xc
   if(nspin_tddft==1) then
     call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_UNPOLARIZED)
   else
     call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_POLARIZED)
   endif
   call xc_f90_info_name(xc_info(idft_xc),string)
   write(stdout,'(a,i4,a,i6,5x,a)') '   XC functional ',idft_xc,' :  ',xc_f90_info_number(xc_info(idft_xc)),&
         TRIM(string)
   if( MODULO(xc_f90_info_flags( xc_info(idft_xc)),XC_FLAGS_HAVE_FXC*2) < XC_FLAGS_HAVE_FXC ) then
     call die('This functional does not have the kernel implemented in Libxc')
   endif
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_GGA     ) require_gradient  =.TRUE.
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) require_gradient  =.TRUE.
 enddo

 !
 ! calculate rho, grad rho and the kernel
 ! 
 ! Get the density matrix P from C
 call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)

 allocate(v2rho2(ngrid,2*nspin_tddft-1),wf_r(ngrid,basis%nbf,nspin))
 v2rho2(:,:) = 0.0_dp

 if(require_gradient) then
   allocate(vsigma(ngrid,2*nspin_tddft-1))
   allocate(v2rhosigma(ngrid,5*nspin_tddft-4))
   allocate(v2sigma2(ngrid,5*nspin_tddft-4))
   allocate(wf_gradr(3,ngrid,basis%nbf,nspin))
   allocate(rho_gradr(3,ngrid,nspin))
   vsigma(:,:)     = 0.0_dp
   v2rhosigma(:,:) = 0.0_dp
   v2sigma2(:,:)   = 0.0_dp
 endif


 max_v2sigma2 = -1.0_dp
 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)

   if( .NOT. ALLOCATED(bfr) ) call prepare_basis_functions_r(basis)
   if( require_gradient .AND. .NOT. ALLOCATED(bfgr) ) call prepare_basis_functions_gradr(basis)
   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r(basis,igrid,basis_function_r)
   !
   ! store the wavefunction in r
   do ispin=1,nspin
     wf_r(igrid,:,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,:,ispin) )
   enddo

   if( require_gradient ) then
     call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
     !
     ! store the wavefunction in r
     do ispin=1,nspin
       wf_gradr(:,igrid,:,ispin) = MATMUL( basis_function_gradr(:,:) , c_matrix(:,:,ispin) )
     enddo
   endif


   call calc_density_pmatrix(nspin,basis,p_matrix,rr,basis_function_r,rhor_r)
   if( require_gradient ) then
     call calc_density_gradr_pmatrix(nspin,basis%nbf,p_matrix,basis_function_r,basis_function_gradr,grad_rhor)

     rho_gradr(:,igrid,:) = grad_rhor(:,:)

     if( nspin_tddft==1 ) then
       sigma_c(1) = SUM( grad_rhor(:,1)**2 )
     else if( nspin==2 ) then
       sigma_c(2) = SUM( grad_rhor(:,1) * grad_rhor(:,2) )
       sigma_c(3) = SUM( grad_rhor(:,2)**2 )
     else ! triplet excitations from singlet ground-state
       sigma_c(:) = SUM( grad_rhor(:,1)**2 ) * 0.25_dp
     endif

   endif


   if( nspin_tddft==1 ) then
     rho_c(1) = rhor_r(1)
   else if( nspin==2 ) then
     rho_c(:) = rhor_r(:)
   else ! triplet excitations from singlet ground-state
     rho_c(:) = rhor_r(1) * 0.5_dp
   endif

   !
   ! Calculate the kernel
   ! 
   do idft_xc=1,ndft_xc
     select case(xc_f90_info_family(xc_info(idft_xc)))
     case(XC_FAMILY_LDA)
       call xc_f90_lda_fxc(xc_func(idft_xc),1_C_INT,rho_c(1),v2rho2_c(1))
     case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
       call xc_f90_gga_vxc(xc_func(idft_xc),1_C_INT,rho_c(1),sigma_c(1),vrho_c(1),vsigma_c(1))
       call xc_f90_gga_fxc(xc_func(idft_xc),1_C_INT,rho_c(1),sigma_c(1),v2rho2_c(1),v2rhosigma_c(1),v2sigma2_c(1))
     case default
       call die('Other kernels not yet implemented')
     end select
     !
     ! Remove the too large values for stability
     v2rho2_c(:) = MIN( v2rho2_c(:), kernel_capping )
     if(require_gradient) then
       max_v2sigma2 = MAX(ABS(v2sigma2_c(1)),max_v2sigma2)
       vsigma_c(:)     = MIN( vsigma_c(:), kernel_capping )
       v2rhosigma_c(:) = MIN( v2rhosigma_c(:), kernel_capping )
       v2sigma2_c(:)   = MIN( v2sigma2_c(:), kernel_capping )
     endif

     ! Store the result with the weight
     ! Remove too large values for stability
     v2rho2(igrid,:)     = v2rho2(igrid,:) + v2rho2_c(:) * w_grid(igrid) * dft_xc_coef(idft_xc)
     if(require_gradient) then
       vsigma(igrid,:)     = vsigma(igrid,:)     + vsigma_c(:)     * w_grid(igrid) * dft_xc_coef(idft_xc)
       v2rhosigma(igrid,:) = v2rhosigma(igrid,:) + v2rhosigma_c(:) * w_grid(igrid) * dft_xc_coef(idft_xc)
       v2sigma2(igrid,:)   = v2sigma2(igrid,:)   + v2sigma2_c(:)   * w_grid(igrid) * dft_xc_coef(idft_xc)
     endif

   enddo
 enddo
 if(require_gradient) then
   call xmax(max_v2sigma2)
   write(stdout,'(a,e18.6)') ' Maximum numerical value for fxc: ',max_v2sigma2
 endif


end subroutine prepare_tddft


!=========================================================================
