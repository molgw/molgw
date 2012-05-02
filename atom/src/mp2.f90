!=========================================================================
!ERI subroutine mp2_energy(nspin,basis,occupation,c_matrix,energy,eri,emp2)
subroutine mp2_energy(nspin,basis,occupation,c_matrix,energy,emp2)
 use m_definitions
 use m_basis_set
 use m_eri
 implicit none

 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(basis%nbf,nspin),c_matrix(basis%nbf,basis%nbf,nspin),energy(basis%nbf,nspin)
!ERI real(dp),intent(in)        :: eri(basis%nbf,basis%nbf,basis%nbf,basis%nbf)
 real(dp),intent(out)       :: emp2
!=====
 integer                    :: aorbital,borbital,iorbital,jorbital
 integer                    :: ibf,jbf,abf,bbf,ispin,jspin
 real(dp)                   :: energy_denom
 real(dp)                   :: tmp_xaxb(basis%nbf,basis%nbf)
 real(dp)                   :: tmp_iaxb(basis%nbf),tmp_xaib(basis%nbf)
 real(dp)                   :: tmp_iajb,tmp_jaib
 real(dp)                   :: contrib1,contrib2
 real(dp)                   :: fact,spin_fact
 real(dp),allocatable       :: tmp_xaxx(:,:,:)
!=====

 emp2 = 0.0_dp
 contrib1 = 0.0_dp
 contrib2 = 0.0_dp

 spin_fact = REAL(-nspin+3,dp)


 allocate(tmp_xaxx(basis%nbf,basis%nbf,basis%nbf))

 do ispin=1,nspin

   do aorbital=1,basis%nbf
     if(occupation(aorbital,ispin)>spin_fact-completely_empty) cycle
     write(*,'(i4,2x,i4,a,i4)') &
&             ispin,aorbital,' / ',basis%nbf

     tmp_xaxx(:,:,:) = 0.0_dp
     do bbf=1,basis%nbf
       do jbf=1,basis%nbf
         do abf=1,basis%nbf
           do ibf=1,basis%nbf
             tmp_xaxx(ibf,jbf,bbf) = tmp_xaxx(ibf,jbf,bbf) &
&               + c_matrix(abf,aorbital,ispin) * eri(ibf,abf,jbf,bbf)
           enddo
         enddo
       enddo
     enddo

     do jspin=1,nspin
       do borbital=1,basis%nbf
         if(occupation(borbital,jspin)>spin_fact-completely_empty) cycle
       
         tmp_xaxb(:,:) = 0.0_dp
         do bbf=1,basis%nbf
           do jbf=1,basis%nbf
             do ibf=1,basis%nbf
               tmp_xaxb(ibf,jbf) = tmp_xaxb(ibf,jbf) + c_matrix(bbf,borbital,jspin) * tmp_xaxx(ibf,jbf,bbf)
             enddo
           enddo
         enddo

         do iorbital=1,basis%nbf
           if(occupation(iorbital,ispin)<completely_empty) cycle

           tmp_iaxb(:) = 0.0_dp
           tmp_xaib(:) = 0.0_dp
           do jbf=1,basis%nbf
             do ibf=1,basis%nbf
               tmp_iaxb(jbf) = tmp_iaxb(jbf) + c_matrix(ibf,iorbital,ispin) * tmp_xaxb(ibf,jbf)
             enddo
           enddo
           do jbf=1,basis%nbf
             do ibf=1,basis%nbf
               tmp_xaib(jbf) = tmp_xaib(jbf) + c_matrix(ibf,iorbital,ispin) * tmp_xaxb(jbf,ibf)
             enddo
           enddo

           do jorbital=1,basis%nbf
             if(occupation(jorbital,jspin)<completely_empty) cycle
!             fact = ( occupation(iorbital,ispin) - occupation(aorbital,ispin) ) &
!&                     * ( occupation(jorbital,jspin) - occupation(borbital,jspin) )
!             if( ABS(fact) < 1.d-6 ) cycle
             fact =  occupation(iorbital,ispin) * ( 1.0_dp - occupation(aorbital,ispin) ) &
&                   *occupation(jorbital,jspin) * ( 1.0_dp - occupation(borbital,jspin) )

             energy_denom = energy(iorbital,ispin) + energy(jorbital,jspin) &
&                                    - energy(aorbital,ispin) - energy(borbital,jspin)
             ! Avoid the zero denominators
             if( ABS(energy_denom) < 1.d-18) then
               write(*,*) 'write that you skipped something'
               cycle
             endif
             energy_denom =  fact / energy_denom 

             tmp_iajb = SUM( tmp_iaxb(:) * c_matrix(:,jorbital,jspin) )
             tmp_jaib = SUM( tmp_xaib(:) * c_matrix(:,jorbital,jspin) )

             contrib1 = contrib1 + 0.5_dp * energy_denom * tmp_iajb**2 

             if(ispin==jspin) contrib2 = contrib2 - 0.5_dp * energy_denom * tmp_iajb*tmp_jaib / spin_fact


           enddo
         enddo
       enddo
     enddo !jspin

   enddo ! aorbital

 enddo !ispin

 deallocate(tmp_xaxx)


 emp2 = contrib1 + contrib2
 write(*,'(/,a)')       ' MP2 contributions'
 write(*,'(a,f14.8)')   ' 2-ring diagram  :',contrib1
 write(*,'(a,f14.8)')   ' SOX diagram     :',contrib2
 write(*,'(a,f14.8,/)') ' MP2 correlation :',emp2


end subroutine mp2_energy

!==================================================================
subroutine full_ci_2electrons_spin(spinstate,nbf,h_1e,c_matrix)
 use m_definitions
 use m_tools
 use m_eri
! use wavefunction_object
! use basis_object
 implicit none
!
 integer,parameter :: cip=dp
 integer,parameter :: nx=60
!
 integer,intent(in) :: spinstate,nbf
 real(dp),intent(in) :: h_1e(nbf,nbf),c_matrix(nbf,nbf)
!=====
 real(dp) :: h_1e_hf(nbf,nbf)
 integer :: nconf,iconf,jconf
 integer :: ibf,jbf,kbf,lbf
 integer :: istate,jstate,kstate,lstate
 integer :: istate1,istate2,jstate1,jstate2,ispin1,ispin2,jspin1,jspin2
 real(cip),allocatable :: hamiltonian(:,:)
 real(cip),allocatable :: energy(:),eigenvector(:,:)
 real(cip),allocatable :: test1(:),test2(:),test3(:),hphi(:),gr1(:),gr2(:)
 real(cip) :: delta,eigen,hessian,norm_gr1
 integer :: iline,ix
 real(dp) :: rhor(nx),rr(3)
 real(dp) :: rhor_t(nx)
#ifndef LOW_MEMORY2
 real(dp) :: eri_hf(nbf,nbf,nbf,nbf)
#else
 real(dp) :: eri_hf_i(nbf,nbf,nbf,1)
#endif
!
 integer,parameter :: ny=nx,nz=nx
 integer :: iy,iz
 real(dp) :: x(nx),y(ny),z(nz)
 real(dp) :: wx(nx),wy(ny),wz(nz)
 real(dp) :: norm
!
 write(*,*) 
 write(*,*) 'Enter full CI subroutine'
 write(*,*) 

 write(*,*) 'obtain the one-electron Hamiltonian in the HF basis'
 h_1e_hf(:,:) = 0.0_dp
 do jstate=1,nbf
   do istate=1,nbf
     do jbf=1,nbf
       do ibf=1,nbf
         h_1e_hf(istate,jstate) = h_1e_hf(istate,jstate) + h_1e(ibf,jbf) * c_matrix(ibf,istate) * c_matrix(jbf,jstate)
       enddo
     enddo
   enddo
 enddo

#ifndef LOW_MEMORY2
 call transform_eri_basis_fast(nbf,1,c_matrix,eri_hf)
#endif

 select case(spinstate)
 case(0)
   write(*,*) 'calculate spin singlet state'
 case(1)
   write(*,*) 'calculate spin triplet state'
 case default
   stop'BUG: spin state not possible'
 end select

 nconf = ( 2*nbf * (2*nbf -1) ) / 2
 write(*,*)
 write(*,*) 'CI matrix lower than',nconf,' x ',nconf
 allocate(hamiltonian(nconf,nconf))
 hamiltonian(:,:) = 0.0_dp

 iconf=0
 do istate1=1,nbf
#ifdef LOW_MEMORY2
     call transform_eri_basis_lowmem(1,c_matrix,istate1,1,eri_hf_i)
#endif
   do ispin1=1,2

     do istate2=istate1,nbf
       do ispin2=1,2
         if(istate1==istate2 .AND. (ispin2==1 .OR. ispin1==2) ) cycle
         !
         ! S^2 selection
!         if(ABS(ispin1-ispin2)==spinstate) cycle
         !
         ! for two electrons, the two electron wavefunction is even in space
         ! the two one-particle wavefunctions have then to have the parity 
!TODO         if(.NOT.symmetry(istate1,istate2)) cycle
         iconf=iconf+1

         jconf=0
         do jstate1=1,nbf
           do jspin1=1,2
             do jstate2=jstate1,nbf
               do jspin2=1,2
                 if(jstate1==jstate2 .AND. (jspin2==1 .OR. jspin1==2) ) cycle
                 !
                 ! S^2 selection
!                 if(ABS(jspin1-jspin2)==spinstate) cycle
                 !
                 ! for two electrons, the two electron wavefunction is even in space
                 ! the two one-particle wavefunctions have then to have the parity 
!TODO                 if(.NOT.symmetry(jstate1,jstate2)) cycle
                 jconf=jconf+1

!         write(*,'(10(i4,x))') jconf,jstate1,jspin1,jstate2,jspin2


                 if( jspin1==ispin1 .AND. jspin2==ispin2 ) then
                   if(istate1==jstate1) then
                     if(istate2==jstate2) then
                       hamiltonian(iconf,jconf) = h_1e_hf(istate1,jstate1)  + h_1e_hf(istate2,jstate2) 
                     else
                       hamiltonian(iconf,jconf) = h_1e_hf(istate2,jstate2) 
                     endif
                   else
                     if(istate2==jstate2) then
                       hamiltonian(iconf,jconf) = h_1e_hf(istate1,jstate1) 
                     endif
                   endif
                 endif
                 if( jspin1==ispin2 .AND. jspin2==ispin1 ) then
                   if(istate1==jstate2) then
                     if(istate2==jstate1) then
                       hamiltonian(iconf,jconf) = h_1e_hf(istate1,jstate2)  + h_1e_hf(istate2,jstate1) 
                     else
                       hamiltonian(iconf,jconf) = h_1e_hf(istate2,jstate1) 
                     endif
                   else
                     if(istate2==jstate1) then
                       hamiltonian(iconf,jconf) = h_1e_hf(istate1,jstate2) 
                     endif
                   endif
                 endif


                 !
                 ! Painful implementation of the determinant rules as shown in
                 ! p. 70 of "Modern Quantum Chemistry" by A. Szabo and N. S. Ostlung

#ifndef LOW_MEMORY2
                 if( istate1==jstate1 .AND. ispin1==jspin1 .AND.  istate2==jstate2 .AND. ispin2==jspin2 ) then

                   hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate1,jstate1,istate2,jstate2)
                   if(ispin1==ispin2) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) - eri_hf(istate1,jstate2,istate2,jstate1)
                   endif

                 else if ( istate1==jstate1 .AND. ispin1==jspin1 ) then

                   if(ispin2==jspin2) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate2,jstate2,istate1,jstate1)
                   endif
                   if(ispin2==jspin1 .AND. jspin2==ispin1 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate2,jstate1,istate1,jstate2)
                   endif

                 else if ( istate2==jstate2 .AND. ispin2==jspin2 ) then

                   if(ispin1==jspin1) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate1,jstate1,istate2,jstate2)
                   endif
                   if(ispin1==jspin2 .AND. jspin1==ispin2 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate1,jstate2,istate2,jstate1)
                   endif

                 else if ( istate2==jstate1 .AND. ispin2==jspin1 ) then   ! some assymetry due to here

                   if(ispin1==jspin2) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate1,jstate2,istate2,jstate1)
                   endif
                   if(ispin1==jspin1 .AND. jspin2==ispin2 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate1,jstate1,istate2,jstate2)
                   endif

                 else if ( istate1==jstate2 .AND. ispin1==jspin2 ) then  !  corrected here 

                   if(ispin2==jspin1) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate2,jstate1,istate1,jstate2)
                   endif
                   if(ispin2==jspin2 .AND. jspin1==ispin1 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate2,jstate2,istate1,jstate1)
                   endif

                 else ! both states are different

                   if( ispin1==jspin1 .AND. ispin2==jspin2 ) then
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf(istate1,jstate1,istate2,jstate2)
                   endif
                   if( ispin1==jspin2 .AND. ispin2==jspin1 ) then
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) - eri_hf(istate1,jstate2,istate2,jstate1)
                   endif

                 endif

#else

                 if( istate1==jstate1 .AND. ispin1==jspin1 .AND.  istate2==jstate2 .AND. ispin2==jspin2 ) then

                   hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate1,istate2,jstate2,1)
                   if(ispin1==ispin2) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) - eri_hf_i(jstate2,istate2,jstate1,1)
                   endif

                 else if ( istate1==jstate1 .AND. ispin1==jspin1 ) then

                   if(ispin2==jspin2) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate1,istate2,jstate2,1)
                   endif
                   if(ispin2==jspin1 .AND. jspin2==ispin1 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate2,istate2,jstate1,1)
                   endif

                 else if ( istate2==jstate2 .AND. ispin2==jspin2 ) then

                   if(ispin1==jspin1) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate1,istate2,jstate2,1)
                   endif
                   if(ispin1==jspin2 .AND. jspin1==ispin2 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate2,istate2,jstate1,1)
                   endif

                 else if ( istate2==jstate1 .AND. ispin2==jspin1 ) then   ! some assymetry due to here

                   if(ispin1==jspin2) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate2,istate2,jstate1,1)
                   endif
                   if(ispin1==jspin1 .AND. jspin2==ispin2 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate1,istate2,jstate2,1)
                   endif

                 else if ( istate1==jstate2 .AND. ispin1==jspin2 ) then  !  corrected here 

                   if(ispin2==jspin1) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate2,istate2,jstate1,1)
                   endif
                   if(ispin2==jspin2 .AND. jspin1==ispin1 ) then
                     hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate1,istate2,jstate2,1)
                   endif

                 else ! both states are different

                   if( ispin1==jspin1 .AND. ispin2==jspin2 ) then
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) + eri_hf_i(jstate1,istate2,jstate2,1)
                   endif
                   if( ispin1==jspin2 .AND. ispin2==jspin1 ) then
                      hamiltonian(iconf,jconf) =  hamiltonian(iconf,jconf) - eri_hf_i(jstate2,istate2,jstate1,1)
                   endif

                 endif
#endif



               enddo
             enddo
           enddo
         enddo

       enddo
     enddo
   enddo
 enddo

 ! Adjust the real size of the CI hamiltonian
 nconf=iconf
 write(*,*)
 write(*,*) 'CI matrix finally is',nconf,' x ',nconf

 write(*,*) 'Hartree-Fock ground state energy [Ha]',hamiltonian(1,1)
! write(*,*) '=========== H_1e ============== '
! do istate=1,nbf
!   write(*,'(i4,2x,20(x,f12.6))') iconf,h_1e_hf(istate,1:nbf)
! enddo
! write(*,*) '=========== J_ij ============== '
! do istate=1,nbf
!   write(*,'(i4,2x,20(x,f12.6))') iconf,(eri_hf(istate,istate,jstate,jstate),jstate=1,nbf)
! enddo
! write(*,*) '=========== K_ij ============== '
! do istate=1,nbf
!   write(*,'(i4,2x,20(x,f12.6))') iconf,(eri_hf(istate,jstate,jstate,istate),jstate=1,nbf)
! enddo
! write(*,*) '=========== full H ============== '
! do iconf=1,nconf
!   write(*,'(i4,2x,20(x,f12.6))') iconf,hamiltonian(iconf,1:nconf)
! enddo

 allocate(energy(nconf),eigenvector(nconf,nconf))
 !
 ! resize the matrices
 eigenvector(:,:) = hamiltonian(1:nconf,1:nconf)
 deallocate(hamiltonian)
 allocate(hamiltonian(nconf,nconf))
 hamiltonian(:,:) = eigenvector(:,:)

 if(nconf>1000) then ! iterative partial diago or full diago
   ! home made steepest descent algo
   write(*,*)
   write(*,*) 'hamiltonian too big to be diagonalized'
   write(*,*) 'Slow steepest descent algorithm'
   write(*,*) 'home cooked'
   allocate(test1(nconf),test2(nconf),hphi(nconf),gr1(nconf),gr2(nconf))
   test1(:) = 0.001
   test1(1) = 1. 
   test1(:) = test1(:) / SQRT( SUM( test1(:)**2 ) )
   delta = 0.10
  
   do iline=1,500
  
     hphi = matmul(eigenvector,test1)
     eigen = dot_product(test1,hphi)
     gr1(:) = 2. * ( hphi(:) - eigen * test1(:) )
     norm_gr1 = SQRT(dot_product(gr1,gr1))
!     write(*,*) 'x1',norm_gr1
  
!     write(*,*) 'guessed delta',delta
     test2(:) = test1(:) - delta * gr1(:) / norm_gr1
     test2(:) = test2(:) / SQRT( SUM( test2(:)**2 ) )
  
     hphi = matmul(eigenvector,test2)
     eigen = dot_product(test2,hphi)
     gr2(:) = 2. * ( hphi(:) - eigen * test2(:) )
!     write(*,*) 'x2',dot_product(gr2,gr1)/norm_gr1
     hessian = dot_product(gr1,gr2-gr1) / ( delta * norm_gr1 )
     delta = -norm_gr1 / hessian
!     write(*,*) 'optimal delta',delta
     test2(:) = test1(:) - delta * gr1(:) / norm_gr1
     test2(:) = test2(:) / SQRT( SUM( test2(:)**2 ) )
  
     if(modulo(iline,20)==0) write(*,*) 'diff',iline,eigen,SUM(ABS(test2(:)-test1(:)))/DBLE(nconf)
!     write(*,*) '==================='
     test1(:) = test2(:)
  
   enddo

   eigenvector(:,1) =test1(:)
   energy(1) = eigen
   
   deallocate(test1,test2,hphi,gr1,gr2)

 else
   ! full LAPACK diago
   write(*,*) 'starting the diago'
   call diagonalize(nconf,hamiltonian,energy,eigenvector)
   write(*,*) 'diago DONE'

 endif

 write(*,*)
 write(*,*) 'normalization',SUM(eigenvector(:,1)**2)
 write(*,'(i4,2x,20(x,f7.4))') 1,eigenvector(1:min(20,nconf),1)
 write(*,*)
 write(*,'(a30,2x,f12.6)') 'CI ground-state energy [Ha]:',energy(1)
 write(*,'(a30,2x,f12.6)') 'correlation energy [Ha]:',energy(1)-hamiltonian(1,1)
 write(*,*)
  
  ! deallocate(hamiltonian)

#if 0
! Plot the ground state density
 write(*,*)
 write(*,*) 'calculate the density'

 call coeffs_gausslegint(0.0_dp,Lx,x,wx,nx)
 call coeffs_gausslegint(0.0_dp,Ly,y,wy,ny)
 call coeffs_gausslegint(0.0_dp,Lz,z,wz,nz)

! norm=0.0_dp
 rhor(:)=0.0_dp
 iconf=0
 do istate1=1,nbf
   do ispin1=1,2
     do istate2=istate1,nbf
       do ispin2=1,2
         if(istate1==istate2 .AND. (ispin2==1 .OR. ispin1==2) ) cycle
         !
         ! S^2 selection
         if(ABS(ispin1-ispin2)==spinstate) cycle
         !
         ! for two electrons, the two electron wavefunction is even in space
         ! the two one-particle wavefunctions have then to have the parity 
!TODO         if(.NOT.symmetry(istate1,istate2)) cycle
         iconf=iconf+1
         if(modulo(iconf,100)==0) write(*,*) iconf

         jconf=0
         do jstate1=1,nbf
           do jspin1=1,2
             do jstate2=jstate1,nbf
               do jspin2=1,2
                 if(jstate1==jstate2 .AND. (jspin2==1 .OR. jspin1==2) ) cycle
                 !
                 ! S^2 selection
                 if(ABS(jspin1-jspin2)==spinstate) cycle
                 !
                 ! for two electrons, the two electron wavefunction is even in space
                 ! the two one-particle wavefunctions have then to have the parity 
!TODO                 if(.NOT.symmetry(jstate1,jstate2)) cycle
                 jconf=jconf+1

                 do ix=1,nx
!                 do iy=1,ny
!                 do iz=1,nz
!                   rr(1)= x(ix)
!                   rr(2)= y(iy)
!                   rr(3)= z(iz)
                   rr(1)= DBLE(ix-1)/DBLE(nx-1)*Lx
                   rr(2)= DBLE(ix-1)/DBLE(nx-1)*Ly
                   rr(3)= DBLE(ix-1)/DBLE(nx-1)*Lz

                   if( istate2==jstate2 .AND. ispin2==jspin2 .AND. ispin1==jspin1 ) then
                     rhor(ix) = rhor(ix)  &
&                     + eigenvector(iconf,1) * eigenvector(jconf,1) &
&                      * evaluate_wavefunction(basis(istate1),rr) * evaluate_wavefunction(basis(jstate1),rr) !* wy(iy) * wz(iz)
!                     if(istate1==jstate1) norm = norm + eigenvector(iconf,1) * eigenvector(jconf,1) * 0.5
                   endif
                   if( istate1==jstate1 .AND. ispin1==jspin1 .AND. ispin2==jspin2 ) then 
                     rhor(ix) = rhor(ix)  &
&                     + eigenvector(iconf,1) * eigenvector(jconf,1) &
&                      * evaluate_wavefunction(basis(istate2),rr) * evaluate_wavefunction(basis(jstate2),rr) !* wy(iy) * wz(iz)
!                     if(istate2==jstate2) norm = norm + eigenvector(iconf,1) * eigenvector(jconf,1) *0.5
                   endif
                   if( istate2==jstate1 .AND. ispin2==jspin1 .AND. ispin1==jspin2 ) then
                     rhor(ix) = rhor(ix)  &
&                     - eigenvector(iconf,1) * eigenvector(jconf,1) &
&                      * evaluate_wavefunction(basis(istate1),rr) * evaluate_wavefunction(basis(jstate2),rr) !* wy(iy) * wz(iz)
!                     if(istate1==jstate2) norm = norm + eigenvector(iconf,1) * eigenvector(jconf,1) * 0.5
                   endif
                   if( istate1==jstate2 .AND. ispin1==jspin2 .AND. ispin2==jspin1 ) then
                     rhor(ix) = rhor(ix)  &
&                     - eigenvector(iconf,1) * eigenvector(jconf,1) &
&                      * evaluate_wavefunction(basis(istate2),rr) * evaluate_wavefunction(basis(jstate1),rr) !* wy(iy) * wz(iz)
!                     if(istate2==jstate1) norm = norm + eigenvector(iconf,1) * eigenvector(jconf,1) *0.5
                   endif

!                 enddo !iz
!                 enddo !iy
                 enddo !ix


               enddo
             enddo
           enddo
         enddo

       enddo
     enddo
   enddo
 enddo
 rhor(:) = rhor(:) * 0.5_dp

! write(*,*) 'NORM',norm / DBLE(nx*ny*nz)
! write(*,*) 'norm',SUM(rhor(:)*wx(:))

! rhor_t(:)=0.0_dp
! do ix=1,nx
! do iy=1,ny
! do iz=1,nz
!   rr(1)= x(ix)
!   rr(2)= y(iy)
!   rr(3)= z(iz)
!   rhor_t(ix)=rhor_t(ix)+ evaluate_wavefunction(basis(1),rr)**2 * wy(iy) * wz(iz)
! enddo !iz
! enddo !iy
! enddo !ix
! write(*,*) 'norm',SUM(rhor_t(:)*wx(:))
! do ix=1,nx
!   rr(1)= x(ix)
!   write(11,'(5(e12.6,2x))') rr(1),rhor_t(ix)
! enddo

 do ix=1,nx
!   rr(1)= x(ix)
   rr(1)= DBLE(ix-1)/DBLE(nx-1)*Lx
   rr(2)= DBLE(ix-1)/DBLE(nx-1)*Lx
   rr(3)= DBLE(ix-1)/DBLE(nx-1)*Lx
   write(10,'(5(e12.6,2x))') rr(1),(evaluate_wavefunction(basis(1),rr))**2,rhor(ix),&
&         rhor(ix)-(evaluate_wavefunction(basis(1),rr))**2
 enddo

#endif

 ! finalize the CI calculation
 deallocate(energy,eigenvector)

 stop'CI calculation stops here'

end subroutine full_ci_2electrons_spin

