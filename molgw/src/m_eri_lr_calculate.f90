!=========================================================================
! This file is part of MOLGW.
!=========================================================================
module m_eri_lr_calculate
 use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
 use m_definitions
 use m_mpi
 use m_scalapack
 use m_memory
 use m_basis_set
 use m_timing
 use m_inputparam,only: scalapack_block_min
 use m_eri


 integer,protected :: nauxil_2center_lr  ! size of the 2-center matrix
                                         ! 2-center integrals are NOT distributed

 real(prec_eri),private,allocatable :: eri_2center_m1_lr(:,:)


contains


!=========================================================================
subroutine calculate_eri_lr(print_eri_,basis,rcut)
 implicit none
 logical,intent(in)           :: print_eri_
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
!=====

 call start_clock(timing_eri_4center)

 write(stdout,'(/,a,i12)') ' Number of integrals to be stored: ',nsize

 call clean_allocate('4-center LR integrals',eri_4center_lr,nsize)
 eri_4center_lr(:) = 0.0_dp

 if( .NOT. read_eri(rcut) ) call calculate_eri_4center_lr(basis,rcut)


 if( print_eri_ ) then
   call dump_out_eri(rcut)
 endif

 call stop_clock(timing_eri_4center)

end subroutine calculate_eri_lr


!=========================================================================
subroutine calculate_eri_4center_lr(basis,rcut)
 use m_tools,only: boys_function
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
!=====
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: ijshellpair,klshellpair
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ig1,ig2,ig3,ig4
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: iibf
 integer                      :: info
 real(dp)                     :: zeta_12,zeta_34,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
!=====
! variables used to call C
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: rcut_libint
!=====

 write(stdout,'(/,a)') ' Calculate and store all the Electron Repulsion Integrals (ERI)'


 rcut_libint = rcut

 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   !
   ! Order the angular momenta so that libint is pleased
   ! 1) am3+am4 >= am1+am2
   ! 2) am3>=am4
   ! 3) am1>=am2
   amk = shell(kshell)%am
   aml = shell(lshell)%am


   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)

     ami = shell(ishell)%am
     amj = shell(jshell)%am
     if( amk+aml < ami+amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     nk = number_basis_function_am( basis%gaussian_type , amk )
     nl = number_basis_function_am( basis%gaussian_type , aml )


     am1 = shell(ishell)%am
     am2 = shell(jshell)%am
     am3 = shell(kshell)%am
     am4 = shell(lshell)%am
     n1c = number_basis_function_am( 'CART' , ami )
     n2c = number_basis_function_am( 'CART' , amj )
     n3c = number_basis_function_am( 'CART' , amk )
     n4c = number_basis_function_am( 'CART' , aml )
     ng1 = shell(ishell)%ng
     ng2 = shell(jshell)%ng
     ng3 = shell(kshell)%ng
     ng4 = shell(lshell)%ng
     allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
     alpha1(:) = shell(ishell)%alpha(:) 
     alpha2(:) = shell(jshell)%alpha(:)
     alpha3(:) = shell(kshell)%alpha(:)
     alpha4(:) = shell(lshell)%alpha(:)
     x01(:) = shell(ishell)%x0(:)
     x02(:) = shell(jshell)%x0(:)
     x03(:) = shell(kshell)%x0(:)
     x04(:) = shell(lshell)%x0(:)
     allocate(coeff1(shell(ishell)%ng))
     allocate(coeff2(shell(jshell)%ng))
     allocate(coeff3(shell(kshell)%ng))
     allocate(coeff4(shell(lshell)%ng))
     coeff1(:)=shell(ishell)%coeff(:)
     coeff2(:)=shell(jshell)%coeff(:)
     coeff3(:)=shell(kshell)%coeff(:)
     coeff4(:)=shell(lshell)%coeff(:)

     allocate( int_shell( n1c*n2c*n3c*n4c ) )
     allocate( integrals_cart(n1c,n2c,n3c,n4c) )
     allocate( integrals_tmp(n1c,n2c,n3c,n4c) )
     integrals_cart(:,:,:,:) = 0.0_dp


     if(am1+am2+am3+am4==0) then

       do ig4=1,ng4
         do ig3=1,ng3
           do ig2=1,ng2
             do ig1=1,ng1

               zeta_12 = alpha1(ig1) + alpha2(ig2)
               zeta_34 = alpha3(ig3) + alpha4(ig4)
               p(:) = ( alpha1(ig1) * x01(:) + alpha2(ig2) * x02(:) ) / zeta_12 
               q(:) = ( alpha3(ig3) * x03(:) + alpha4(ig4) * x04(:) ) / zeta_34 
               !
               ! Treat carefully the LR only integrals
               rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 + zeta_12*zeta_34*rcut**2 )
               rho1 = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
               tt = rho * SUM( (p(:)-q(:))**2 )
               call boys_function(f0t(0),0,tt)

               integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                     2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                     / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                     / zeta_34 * EXP( -alpha3(ig3)*alpha4(ig4)/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) &
                     * SQRT( rho / rho1 ) &
                     * coeff1(ig1) &
                     * coeff2(ig2) &
                     * coeff3(ig3) &
                     * coeff4(ig4) * cart_to_pure_norm(0)%matrix(1,1)**4

             enddo
           enddo
         enddo
       enddo

     else

        
       info=eval_contr_integral(                &
                               am1,am2,am3,am4, &
                               ng1,ng2,ng3,ng4, &
                               coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
                               alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
                               x01(1),x02(1),x03(1),x04(1),&
                               rcut_libint, &
                               int_shell(1))


       if(info/=0) then
         write(stdout,*) am1,am2,am3,am4
         call die('ERI calculated by libint failed')
       endif

       iibf=0
       do ibf=1,n1c
         do jbf=1,n2c
           do kbf=1,n3c
             do lbf=1,n4c
               iibf=iibf+1
               integrals_cart(ibf,jbf,kbf,lbf) = int_shell(iibf)
             enddo
           enddo
         enddo
       enddo


       do lbf=1,n4c
         do kbf=1,n3c
           do jbf=1,n2c
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1c,jbf,kbf,lbf) * cart_to_pure_norm(shell(ishell)%am)%matrix(1:n1c,ibf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n4c
         do kbf=1,n3c
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2c,kbf,lbf) * cart_to_pure_norm(shell(jshell)%am)%matrix(1:n2c,jbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n4c
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3c,lbf) * cart_to_pure_norm(shell(kshell)%am)%matrix(1:n3c,kbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4c) * cart_to_pure_norm(shell(lshell)%am)%matrix(1:n4c,lbf) )
             enddo
           enddo
         enddo
       enddo

     endif ! is (ss|ss)
     
     do lbf=1,nl
       do kbf=1,nk
         do jbf=1,nj
           do ibf=1,ni
             eri_4center_lr( index_eri(shell(ishell)%istart+ibf-1, &
                                       shell(jshell)%istart+jbf-1, &
                                       shell(kshell)%istart+kbf-1, &
                                       shell(lshell)%istart+lbf-1) ) = integrals_cart(ibf,jbf,kbf,lbf)
           enddo
         enddo
       enddo
     enddo


     deallocate(integrals_cart)
     deallocate(integrals_tmp)
     deallocate(int_shell)
     deallocate(alpha1,alpha2,alpha3,alpha4)
     deallocate(coeff1)
     deallocate(coeff2)
     deallocate(coeff3)
     deallocate(coeff4)


   enddo
 enddo


 write(stdout,'(a,/)') ' All LR ERI have been calculated'


end subroutine calculate_eri_4center_lr


!=========================================================================
subroutine calculate_eri_2center_lr(print_eri_,auxil_basis,rcut)
 use m_tools,only: boys_function
 implicit none
 logical,intent(in)           :: print_eri_
 type(basis_set),intent(in)   :: auxil_basis
 real(dp),intent(in)          :: rcut
!=====
 integer                      :: ishell,kshell
 integer                      :: n1c,n3c
 integer                      :: ig1,ig2,ig3,ig4
 integer                      :: ni,nj,nk
 integer                      :: ami,amk
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: iibf
 integer                      :: info
 real(dp)                     :: zeta_12,zeta_34,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:)
 real(dp),allocatable         :: integrals_cart(:,:)
 real(dp),allocatable         :: eigval(:)
!=====
! variables used to call C
 integer(C_INT)               :: am1,am2,am3,am4
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: rcut_libint
!=====

 call start_clock(timing_eri_2center)


 ! First allocate the 2-center integral array
 !
 ! The 2-center integrals are not distributed because they are small anyway
 ! The size of the 2-center integrals might be reduced later on
 ! by removing the zero eigenvalues
 nauxil_2center_lr = auxil_basis%nbf
 !
 ! 2-CENTER LR INTEGRALS 
 !
 call clean_allocate('2-center LR integrals',eri_2center_m1_lr,auxil_basis%nbf,auxil_basis%nbf)

 eri_2center_m1_lr(:,:) = 0.0_dp

 write(stdout,'(/,a)')    ' Calculate, invert and store the 2-center LR Electron Repulsion Integrals'

 rcut_libint = rcut

 do kshell=1,nshell_auxil

   ! Parallelization over the shell index
   if( MODULO(kshell-1,nproc) /= rank ) cycle

   !
   ! Order the angular momenta so that libint is pleased
   ! 1) am3+am4 >= am1+am2
   ! 2) am3>=am4
   ! 3) am1>=am2
   amk = shell_auxil(kshell)%am
   nk  = number_basis_function_am( auxil_basis%gaussian_type , amk )

   do ishell=1,nshell_auxil
     ami = shell_auxil(ishell)%am
     if( amk < ami ) cycle

     ni = number_basis_function_am( auxil_basis%gaussian_type , ami )

     am1 = shell_auxil(ishell)%am
     am2 = 0
     am3 = shell_auxil(kshell)%am
     am4 = 0
     n1c = number_basis_function_am( 'CART' , ami )
     n3c = number_basis_function_am( 'CART' , amk )
     ng1 = shell_auxil(ishell)%ng
     ng2 = 1
     ng3 = shell_auxil(kshell)%ng
     ng4 = 1
     allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
     alpha1(:) = shell_auxil(ishell)%alpha(:) 
     alpha2(:) = 0.0_dp
     alpha3(:) = shell_auxil(kshell)%alpha(:)
     alpha4(:) = 0.0_dp 
     x01(:) = shell_auxil(ishell)%x0(:)
     x02(:) = shell_auxil(ishell)%x0(:)
     x03(:) = shell_auxil(kshell)%x0(:)
     x04(:) = shell_auxil(kshell)%x0(:)
     allocate(coeff1(shell_auxil(ishell)%ng))
     allocate(coeff2(1))
     allocate(coeff3(shell_auxil(kshell)%ng))
     allocate(coeff4(1))
     coeff1(:)=shell_auxil(ishell)%coeff(:)
     coeff2(:)=1.0_dp
     coeff3(:)=shell_auxil(kshell)%coeff(:)
     coeff4(:)=1.0_dp

     allocate( int_shell( n1c*n3c ) )
     allocate( integrals_cart(n1c,n3c) )
     allocate( integrals_tmp(n1c,n3c) )
     integrals_cart(:,:) = 0.0_dp


     if(am1+am3==0) then

       do ig3=1,ng3
         do ig1=1,ng1

           zeta_12 = alpha1(ig1) 
           zeta_34 = alpha3(ig3) 
           p(:) = x01(:)
           q(:) = x03(:)
           !
           ! Full range or long-range only integrals
           rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 + zeta_12*zeta_34*rcut**2 )
           rho1 = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
           
           tt = rho * SUM( (p(:)-q(:))**2 )
           call boys_function(f0t(0),0,tt)

           integrals_cart(1,1) = integrals_cart(1,1) + &
                 2.0_dp * pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                 / zeta_12 & 
                 / zeta_34 &
                 * SQRT( rho / rho1 ) &
                 * coeff1(ig1)* coeff3(ig3) &
                 * cart_to_pure_norm(0)%matrix(1,1)**4

         enddo
       enddo

     else


       info=eval_contr_integral(                &
                               am1,am2,am3,am4, &
                               ng1,ng2,ng3,ng4, &
                               coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
                               alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
                               x01(1),x02(1),x03(1),x04(1),&
                               rcut_libint, &
                               int_shell(1))


       if(info/=0) then
         write(stdout,*) am1,am2,am3,am4
         call die('ERI calculated by libint failed')
       endif

       iibf=0
       do ibf=1,n1c
         do kbf=1,n3c
           iibf=iibf+1
           integrals_cart(ibf,kbf) = int_shell(iibf)
         enddo
       enddo


       do kbf=1,n3c
         do ibf=1,ni
           integrals_tmp (ibf,kbf) = SUM( integrals_cart(1:n1c,kbf) * cart_to_pure_norm(am1)%matrix(1:n1c,ibf) )
         enddo
       enddo

       do kbf=1,n3c
         do ibf=1,ni
           integrals_cart(ibf,kbf) = integrals_tmp (ibf,kbf) * cart_to_pure_norm(0)%matrix(1,1) 
         enddo
       enddo

       do kbf=1,nk
         do ibf=1,ni
           integrals_tmp (ibf,kbf) = SUM( integrals_cart(ibf,1:n3c) * cart_to_pure_norm(am3)%matrix(1:n3c,kbf) )
         enddo
       enddo

       do kbf=1,nk
         do ibf=1,ni
           integrals_cart(ibf,kbf) = integrals_tmp (ibf,kbf) * cart_to_pure_norm(0)%matrix(1,1) 
         enddo
       enddo

     endif
     

     do kbf=1,nk
       do ibf=1,ni
         eri_2center_m1_lr( shell_auxil(ishell)%istart+ibf-1,    &
                            shell_auxil(kshell)%istart+kbf-1 )    = integrals_cart(ibf,kbf)
         !
         ! And the symmetric too only if it is not already one
         ! When amk == ami  , the symmetric part is already calculated. 
         ! We do not want double counting because of the parallelization
         if( amk > ami ) then
           eri_2center_m1_lr( shell_auxil(kshell)%istart+kbf-1,    &
                              shell_auxil(ishell)%istart+ibf-1 )    = integrals_cart(ibf,kbf)
         endif
       enddo
     enddo

     deallocate(integrals_cart)
     deallocate(integrals_tmp)
     deallocate(int_shell)
     deallocate(alpha1,alpha2,alpha3,alpha4)
     deallocate(coeff1,coeff2,coeff3,coeff4)

   enddo
 enddo

 ! Sum up the contribution from the different procs
 call xsum(eri_2center_m1_lr)


 allocate(eigval(nauxil_2center_lr))
 !
 ! Perform in-place diagonalization here with or without scalapack
 call diagonalize_scalapack(scalapack_block_min,nauxil_2center_lr,eri_2center_m1_lr,eigval)

 !
 ! Skip the too small eigenvalues
 ibf = 0
 do jbf=1,nauxil_2center_lr
   if( eigval(jbf) < TOO_LOW_EIGENVAL ) cycle
   ibf = ibf + 1
   eri_2center_m1_lr(:,ibf) = eri_2center_m1_lr(:,jbf) / SQRT( eigval(jbf) )
 enddo
 deallocate(eigval)

 !
 ! Resize the 2-center LR integral matrix
 ! Now the array is eri_2center_m1_lr(auxil_basis%nbf,nauxil_2center_lr)
 nauxil_2center_lr = ibf

 write(stdout,'(a)')        ' All 2-center LR integrals have been calculated, diagonalized and stored'
 write(stdout,'(a,i6)')     ' Some have been eliminated ',auxil_basis%nbf-nauxil_2center_lr
 write(stdout,'(a,es16.6)') ' because they were lower than:',TOO_LOW_EIGENVAL


 call stop_clock(timing_eri_2center)

end subroutine calculate_eri_2center_lr


!=========================================================================
subroutine calculate_eri_3center_lr(print_eri_,basis,auxil_basis,rcut)
 use m_tools,only: boys_function
 implicit none
 logical,intent(in)           :: print_eri_
 type(basis_set),intent(in)   :: basis
 type(basis_set),intent(in)   :: auxil_basis
 real(dp),intent(in)          :: rcut
!=====
 integer                      :: ishell,kshell,lshell
 integer                      :: klshellpair
 integer                      :: n1,n2,n3,n4
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ig1,ig2,ig3,ig4
 integer                      :: ni,nk,nl
 integer                      :: ami,amk,aml
 integer                      :: ibf,jbf,kbf,lbf
 integer                      :: iibf
 integer                      :: ibf_auxil,jbf_auxil,ipair
 integer                      :: info
 real(dp)                     :: zeta_12,zeta_34,rho,rho1,f0t(0:0),tt
 real(dp)                     :: p(3),q(3)
 real(dp),allocatable         :: integrals_tmp(:,:,:,:)
 real(dp),allocatable         :: integrals_cart(:,:,:,:)
 real(dp),allocatable         :: eri_3tmp(:,:,:)
 real(dp),allocatable         :: eri_2tmp(:,:)
 real(dp),allocatable         :: eri_tmp(:,:,:)
!=====
! variables used to call C
 integer(C_INT)               :: am1,am2,am3,am4
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
 real(C_DOUBLE)               :: rcut_libint
!=====

 call start_clock(timing_eri_3center)

 !
 ! First, copy the part of eri_2center_m1 that is actually needed and deallocate
 ! the rest
 allocate(eri_2tmp(nauxil_3center_lr,auxil_basis%nbf))
 do ibf_auxil=1,nauxil_3center_lr
   jbf_auxil = ibf_auxil_g(ibf_auxil)
   eri_2tmp(ibf_auxil,:) = eri_2center_m1_lr(:,jbf_auxil)
 enddo

 write(stdout,*) 'Now deallocate the 2-center LR integrals: not needed anymore'
 call clean_deallocate('2-center LR integrals',eri_2center_m1_lr)


 ! Second, allocate the LR 3-center integral array
 !
 ! LR 3-CENTER INTEGRALS 
 !
 call clean_allocate('3-center integrals',eri_3center_lr,nauxil_3center_lr,npair)


 write(stdout,'(/,a)')    ' Calculate and store all the 3-center LR Electron Repulsion Integrals'

 rcut_libint = rcut

 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   !
   ! Order the angular momenta so that libint is pleased
   ! 1) am3+am4 >= am1+am2
   ! 2) am3>=am4
   ! 3) am1>=am2
   amk = shell(kshell)%am
   aml = shell(lshell)%am
   nk = number_basis_function_am( basis%gaussian_type , amk )
   nl = number_basis_function_am( basis%gaussian_type , aml )

   allocate(eri_3tmp(auxil_basis%nbf,nk,nl))
   eri_3tmp(:,:,:) = 0.0_dp

   do ishell=1,nshell_auxil

     ! Use the distribution to avoid calculating all the integrals
     ! A summation is performed to propagate eri_3tmp to all processors
     if( MODULO(ishell-1,nproc) /= rank ) cycle

     ami = shell_auxil(ishell)%am
     ni = number_basis_function_am( auxil_basis%gaussian_type , ami )


     if( amk+aml >= ami ) then

       am1 = shell_auxil(ishell)%am
       am2 = 0
       am3 = shell(kshell)%am
       am4 = shell(lshell)%am
       n1c = number_basis_function_am( 'CART' , ami )
       n2c = 1
       n3c = number_basis_function_am( 'CART' , amk )
       n4c = number_basis_function_am( 'CART' , aml )
       n1 = ni
       n2 = 1
       n3 = nk
       n4 = nl
       ng1 = shell_auxil(ishell)%ng
       ng2 = 1
       ng3 = shell(kshell)%ng
       ng4 = shell(lshell)%ng
       allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
       allocate(coeff1(ng1),coeff2(ng2),coeff3(ng3),coeff4(ng4))
       alpha1(:) = shell_auxil(ishell)%alpha(:) 
       alpha2(:) = 0.0_dp
       alpha3(:) = shell(kshell)%alpha(:)
       alpha4(:) = shell(lshell)%alpha(:)
       coeff1(:) = shell_auxil(ishell)%coeff(:)
       coeff2(:) = 1.0_dp
       coeff3(:) = shell(kshell)%coeff(:)
       coeff4(:) = shell(lshell)%coeff(:)
       x01(:) = shell_auxil(ishell)%x0(:)
       x02(:) = shell_auxil(ishell)%x0(:)
       x03(:) = shell(kshell)%x0(:)
       x04(:) = shell(lshell)%x0(:)

     else ! interexchange indexes

       am3 = shell_auxil(ishell)%am
       am4 = 0
       am1 = shell(kshell)%am
       am2 = shell(lshell)%am
       n3c = number_basis_function_am( 'CART' , ami )
       n4c = 1
       n1c = number_basis_function_am( 'CART' , amk )
       n2c = number_basis_function_am( 'CART' , aml )
       n3 = ni
       n4 = 1
       n1 = nk
       n2 = nl
       ng3 = shell_auxil(ishell)%ng
       ng4 = 1
       ng1 = shell(kshell)%ng
       ng2 = shell(lshell)%ng
       allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
       allocate(coeff1(ng1),coeff2(ng2),coeff3(ng3),coeff4(ng4))
       alpha3(:) = shell_auxil(ishell)%alpha(:) 
       alpha4(:) = 0.0_dp 
       alpha1(:) = shell(kshell)%alpha(:)
       alpha2(:) = shell(lshell)%alpha(:)
       coeff3(:) = shell_auxil(ishell)%coeff(:)
       coeff4(:) = 1.0_dp
       coeff1(:) = shell(kshell)%coeff(:)
       coeff2(:) = shell(lshell)%coeff(:)
       x03(:) = shell_auxil(ishell)%x0(:)
       x04(:) = shell_auxil(ishell)%x0(:)
       x01(:) = shell(kshell)%x0(:)
       x02(:) = shell(lshell)%x0(:)

     endif

     allocate( int_shell(n1c*n2c*n3c*n4c) )
     allocate( integrals_cart(n1c,n2c,n3c,n4c) )
     allocate( integrals_tmp (n1c,n2c,n3c,n4c) )
     integrals_cart(:,:,:,:) = 0.0_dp


     if(am1+am2+am3+am4==0) then

       do ig4=1,ng4
         do ig3=1,ng3
           do ig2=1,ng2
             do ig1=1,ng1

               zeta_12 = alpha1(ig1) + alpha2(ig2)
               zeta_34 = alpha3(ig3) + alpha4(ig4)
               p(:) = ( alpha1(ig1) * x01(:) + alpha2(ig2) * x02(:) ) / zeta_12 
               q(:) = ( alpha3(ig3) * x03(:) + alpha4(ig4) * x04(:) ) / zeta_34 
               !
               ! Full range or long-range only integrals
               rho  = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 + zeta_12*zeta_34*rcut**2 )
               rho1 = zeta_12 * zeta_34 / ( zeta_12 + zeta_34 )
               
               tt = rho * SUM( (p(:)-q(:))**2 )
               call boys_function(f0t(0),0,tt)

               integrals_cart(1,1,1,1) = integrals_cart(1,1,1,1) + &
                     2.0_dp*pi**(2.5_dp) / SQRT( zeta_12 + zeta_34 ) * f0t(0) &
                     / zeta_12 * EXP( -alpha1(ig1)*alpha2(ig2)/zeta_12 * SUM( (x01(:)-x02(:))**2 ) ) & 
                     / zeta_34 * EXP( -alpha3(ig3)*alpha4(ig4)/zeta_34 * SUM( (x03(:)-x04(:))**2 ) ) &
                     * SQRT( rho / rho1 ) &
                     * coeff1(ig1) * coeff2(ig2) &
                     * coeff3(ig3) * coeff4(ig4)&
                     * cart_to_pure_norm(0)%matrix(1,1)**4

             enddo
           enddo
         enddo
       enddo

     else

       info=eval_contr_integral(                &
                               am1,am2,am3,am4, &
                               ng1,ng2,ng3,ng4, &
                               coeff1(1),coeff2(1),coeff3(1),coeff4(1),&
                               alpha1(1),alpha2(1),alpha3(1),alpha4(1),&
                               x01(1),x02(1),x03(1),x04(1),&
                               rcut_libint, &
                               int_shell(1))


       if(info/=0) then
         write(stdout,*) am1,am2,am3,am4
         call die('ERI calculated by libint failed')
       endif

       iibf=0
       do ibf=1,n1c
         do jbf=1,n2c
           do kbf=1,n3c
             do lbf=1,n4c
               iibf=iibf+1
               integrals_cart(ibf,jbf,kbf,lbf) = int_shell(iibf)
             enddo
           enddo
         enddo
       enddo


       do lbf=1,n4c
         do kbf=1,n3c
           do jbf=1,n2c
             do ibf=1,n1
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(1:n1c,jbf,kbf,lbf) * cart_to_pure_norm(am1)%matrix(1:n1c,ibf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n4c
         do kbf=1,n3c
           do jbf=1,n2
             do ibf=1,n1
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,1:n2c,kbf,lbf) * cart_to_pure_norm(am2)%matrix(1:n2c,jbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n4c
         do kbf=1,n3
           do jbf=1,n2
             do ibf=1,n1
               integrals_tmp (ibf,jbf,kbf,lbf) = SUM( integrals_cart(ibf,jbf,1:n3c,lbf) * cart_to_pure_norm(am3)%matrix(1:n3c,kbf) )
             enddo
           enddo
         enddo
       enddo

       do lbf=1,n4
         do kbf=1,n3
           do jbf=1,n2
             do ibf=1,n1
               integrals_cart(ibf,jbf,kbf,lbf) = SUM( integrals_tmp (ibf,jbf,kbf,1:n4c) * cart_to_pure_norm(am4)%matrix(1:n4c,lbf) )
             enddo
           enddo
         enddo
       enddo


     endif ! is (ss|ss)
     
     if(amk+aml>=ami) then
       
       do lbf=1,nl
         do kbf=1,nk
           do ibf=1,ni
             ibf_auxil = shell_auxil(ishell)%istart+ibf-1
             eri_3tmp(ibf_auxil,kbf,lbf) = integrals_cart(ibf,1,kbf,lbf)
           enddo
         enddo
       enddo

     else

       do lbf=1,nl
         do kbf=1,nk
           do ibf=1,ni
             ibf_auxil = shell_auxil(ishell)%istart+ibf-1
             eri_3tmp(ibf_auxil,kbf,lbf) = integrals_cart(kbf,lbf,ibf,1)
           enddo
         enddo
       enddo

     endif


     deallocate(integrals_cart)
     deallocate(integrals_tmp)
     deallocate(int_shell)
     deallocate(alpha1,alpha2,alpha3,alpha4)
     deallocate(coeff1,coeff2,coeff3,coeff4)

   enddo

   call barrier()
   call barrier()

   ! Parallelization over the auxiliary shell
   call xsum(eri_3tmp)

   !
   ! Combine the 2-center integral with the 3-center here
   !
   allocate(eri_tmp(nauxil_3center_lr,nk,nl))
   call DGEMM('N','N',nauxil_3center_lr,nk*nl,auxil_basis%nbf,1.0_dp,eri_2tmp,nauxil_3center_lr,eri_3tmp,auxil_basis%nbf,0.0_dp,eri_tmp,nauxil_3center_lr)

   do lbf=1,nl
     do kbf=1,nk
       ipair = index_pair(shell(kshell)%istart+kbf-1,shell(lshell)%istart+lbf-1)
       eri_3center_lr(:,ipair) = eri_tmp(:,kbf,lbf)
     enddo
   enddo
   deallocate(eri_tmp)

   deallocate(eri_3tmp)

 enddo

 write(stdout,'(a)') ' All 3-center LR integrals have been calculated and stored'

 deallocate(eri_2tmp)

 call stop_clock(timing_eri_3center)


end subroutine calculate_eri_3center_lr


!=========================================================================
end module m_eri_lr_calculate
