!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of electronic density derived quantities:
! rho(r), grad rho(r), vxc in specific case, etc..
!
!=========================================================================
subroutine setup_atomic_density(rr,rhor)
 use m_definitions
 use m_atoms
 use m_gaussian
 use m_inputparam
 implicit none

 real(dp),intent(in)  :: rr(3)
 real(dp),intent(out) :: rhor
!=====
 real(dp),parameter   :: bondcharge=1.000_dp
 integer              :: iatom,igau,ngau
 real(dp)             :: dr
 real(dp),allocatable :: alpha(:),coeff(:)
!=====

 rhor = 0.0_dp
 do iatom=1,natom

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),basis_element(iatom),coeff,alpha)

   dr=NORM2( rr(:) - xatom(:,iatom) )

   do igau=1,ngau
     rhor     = rhor     + SQRT(alpha(igau)/pi)**3 * EXP( -alpha(igau)*dr**2) * coeff(igau)
   enddo

   deallocate(alpha,coeff)
 enddo



end subroutine setup_atomic_density


!=========================================================================
subroutine calc_density_pmatrix(nspin,basis,p_matrix,basis_function_r,rhor)
 use m_definitions
 use m_mpi
 use m_basis_set
 implicit none
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)  :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)  :: basis_function_r(basis%nbf)
 real(dp),intent(out) :: rhor(nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 !
 ! Calculate the density rho at point r
 rhor(:)=0.0_dp
 do ispin=1,nspin
   do jbf=1,basis%nbf
!     if( SUM( (basis%bff(jbf)%x0(:) - rr(:))**2 ) > bf_rad2(jbf) ) cycle
     do ibf=1,basis%nbf
       rhor(ispin)=rhor(ispin)+p_matrix(ibf,jbf,ispin)&
                         * basis_function_r(ibf) &
                         * basis_function_r(jbf)
     enddo
   enddo
 enddo


end subroutine calc_density_pmatrix


!=========================================================================
subroutine calc_density_r(nspin,nbf,nstate,occupation,c_matrix,basis_function_r,rhor)
 use m_definitions
 use m_mpi
 use m_basis_set
 implicit none

 integer,intent(in)         :: nspin,nbf,nstate
 real(dp),intent(in)        :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: basis_function_r(nbf)
 real(dp),intent(out)       :: rhor(nspin)
!=====
 integer              :: ispin,istate
 real(dp)             :: phi_ir
 real(dp),allocatable :: phir(:)
 integer              :: nocc
!=====

 !
 ! Calculate the density rho at point r
 rhor(:)=0.0_dp

 do ispin=1,nspin

#if 0
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle

     phi_ir = DOT_PRODUCT( basis_function_r(:) , c_matrix(:,istate,ispin) )
     rhor(ispin) = rhor(ispin) + phi_ir**2 * occupation(istate,ispin)

   enddo
#else
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     nocc = istate
   enddo
   allocate(phir(nocc))
   phir(:) = MATMUL( basis_function_r(:) , c_matrix(:,:nocc,ispin) )

   rhor(ispin) = rhor(ispin) + SUM( phir(:)**2 * occupation(:nocc,ispin) )
   deallocate(phir)

#endif

 enddo


end subroutine calc_density_r

!=========================================================================
subroutine calc_density_r_batch(nspin,nbf,nstate,nr,occupation,c_matrix,basis_function_r,rhor)
 use m_definitions
 use m_mpi
 use m_basis_set
 implicit none

 integer,intent(in)         :: nspin,nbf,nstate,nr
 real(dp),intent(in)        :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: basis_function_r(nbf,nr)
 real(dp),intent(out)       :: rhor(nspin,nr)
!=====
 integer              :: ispin,istate,ir
 real(dp),allocatable :: phir(:,:)
 integer              :: nocc
!=====

 !
 ! Calculate the density rho at points in batch
 rhor(:,:)=0.0_dp

 do ispin=1,nspin

   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     nocc = istate
   enddo

   allocate(phir(nocc,nr))
   phir(:,:) = MATMUL( TRANSPOSE(c_matrix(:,:nocc,ispin)) , basis_function_r(:,:) )

   forall(ir=1:nr)
     rhor(ispin,ir) = rhor(ispin,ir) + SUM( phir(:,ir)**2 * occupation(:nocc,ispin) )
   endforall
   deallocate(phir)


 enddo


end subroutine calc_density_r_batch

!=========================================================================
subroutine calc_density_r_batch_cmplx(nspin,nbf,nstate,nocc,nr,occupation,c_matrix_cmplx,basis_function_r,rhor)
 use m_definitions
 use m_mpi
 use m_basis_set
 implicit none

 integer,intent(in)         :: nspin,nbf,nstate,nr,nocc
 complex(dp),intent(in)     :: c_matrix_cmplx(nbf,nocc,nspin)
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: basis_function_r(nbf,nr)
 real(dp),intent(out)       :: rhor(nspin,nr)
!=====
 integer                 :: ispin,istate,ir
 complex(dp),allocatable :: phir_cmplx(:,:)
!=====

 !
 ! Calculate the density rho at points in batch
 rhor(:,:)=0.0_dp

 do ispin=1,nspin

   allocate(phir_cmplx(nocc,nr))
   phir_cmplx(:,:) = MATMUL( TRANSPOSE(c_matrix_cmplx(:,:nocc,ispin)) , basis_function_r(:,:) )

   forall(ir=1:nr)
     rhor(ispin,ir) = rhor(ispin,ir) + REAL( SUM(phir_cmplx(:,ir) * CONJG(phir_cmplx(:,ir)) * occupation(:nocc,ispin) ) )
   endforall
   deallocate(phir_cmplx)


 enddo


end subroutine calc_density_r_batch_cmplx


!=========================================================================
subroutine calc_density_gradr_pmatrix(nspin,nbf,p_matrix,basis_function_r,basis_function_gradr,grad_rhor)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nspin,nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(in)  :: basis_function_r(nbf)
 real(dp),intent(in)  :: basis_function_gradr(3,nbf)
 real(dp),intent(out) :: grad_rhor(3,nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 grad_rhor(:,:)=0.0_dp
 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       grad_rhor(:,ispin) = grad_rhor(:,ispin) + p_matrix(ibf,jbf,ispin) &
            *( basis_function_gradr(:,ibf) * basis_function_r(jbf) &
             + basis_function_gradr(:,jbf) * basis_function_r(ibf) )
     enddo
   enddo
 enddo

end subroutine calc_density_gradr_pmatrix


!=========================================================================
subroutine calc_density_gradr(nspin,nbf,nstate,occupation,c_matrix,basis_function_r,basis_function_gradr,grad_rhor)
 use m_definitions
 use m_mpi
 use m_basis_set
 implicit none
 integer,intent(in)         :: nspin,nbf,nstate
 real(dp),intent(in)        :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: basis_function_r(nbf)
 real(dp),intent(in)        :: basis_function_gradr(3,nbf)
 real(dp),intent(out)       :: grad_rhor(3,nspin)
!=====
 integer              :: ispin,ibf,istate
 real(dp)             :: phi_ir
 real(dp)             :: grad_phi_ir(3)
!=====

 !
 ! Calculate the density gradient \nabla rho at point r
 grad_rhor(:,:) = 0.0_dp

 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle

     phi_ir         = DOT_PRODUCT( basis_function_r(:) , c_matrix(:,istate,ispin) )
     grad_phi_ir(:) = MATMUL( basis_function_gradr(:,:) , c_matrix(:,istate,ispin) )

     grad_rhor(:,ispin) = grad_rhor(:,ispin) + phi_ir * grad_phi_ir(:) * 2.0_dp * occupation(istate,ispin)

   enddo
 enddo


end subroutine calc_density_gradr

!========================================================================
subroutine calc_density_gradr_batch(nspin,nbf,nstate,nr,occupation,c_matrix,basis_function_r,basis_function_gradr,rhor,grad_rhor)
 use m_definitions
 use m_mpi
 use m_basis_set
 implicit none

 integer,intent(in)         :: nspin,nbf,nstate,nr
 real(dp),intent(in)        :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: basis_function_r(nbf,nr)
 real(dp),intent(in)        :: basis_function_gradr(nbf,nr,3)
 real(dp),intent(out)       :: rhor(nspin,nr)
 real(dp),intent(out)       :: grad_rhor(nspin,nr,3)
!=====
 integer              :: ispin,istate,ir
 real(dp),allocatable :: phir(:,:)
 real(dp),allocatable :: phir_gradx(:,:)
 real(dp),allocatable :: phir_grady(:,:)
 real(dp),allocatable :: phir_gradz(:,:)
 integer              :: nocc
!=====

 !
 ! Calculate rho and grad rho at points in batch
 rhor(:,:)        = 0.0_dp
 grad_rhor(:,:,:) = 0.0_dp

 do ispin=1,nspin

   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     nocc = istate
   enddo

   allocate(phir(nocc,nr))
   allocate(phir_gradx(nocc,nr))
   allocate(phir_grady(nocc,nr))
   allocate(phir_gradz(nocc,nr))
   phir(:,:)       = MATMUL( TRANSPOSE(c_matrix(:,:nocc,ispin)) , basis_function_r(:,:) )
   phir_gradx(:,:) = MATMUL( TRANSPOSE(c_matrix(:,:nocc,ispin)) , basis_function_gradr(:,:,1) )
   phir_grady(:,:) = MATMUL( TRANSPOSE(c_matrix(:,:nocc,ispin)) , basis_function_gradr(:,:,2) )
   phir_gradz(:,:) = MATMUL( TRANSPOSE(c_matrix(:,:nocc,ispin)) , basis_function_gradr(:,:,3) )

   forall(ir=1:nr)
     rhor(ispin,ir)        = rhor(ispin,ir)        + SUM( phir(:,ir)**2 * occupation(:nocc,ispin) )
     grad_rhor(ispin,ir,1) = grad_rhor(ispin,ir,1) + 2.0_dp * SUM(  phir(:,ir) * phir_gradx(:,ir) * occupation(:nocc,ispin) )
     grad_rhor(ispin,ir,2) = grad_rhor(ispin,ir,2) + 2.0_dp * SUM(  phir(:,ir) * phir_grady(:,ir) * occupation(:nocc,ispin) )
     grad_rhor(ispin,ir,3) = grad_rhor(ispin,ir,3) + 2.0_dp * SUM(  phir(:,ir) * phir_gradz(:,ir) * occupation(:nocc,ispin) )
   endforall

   deallocate(phir)
   deallocate(phir_gradx,phir_grady,phir_gradz)


 enddo


end subroutine calc_density_gradr_batch

!========================================================================
subroutine calc_density_gradr_batch_cmplx(nspin,nbf,nstate,nocc,nr,occupation,c_matrix_cmplx,basis_function_r,basis_function_gradr,rhor,grad_rhor)
 use m_definitions
 use m_mpi
 use m_basis_set
 implicit none

 integer,intent(in)         :: nspin,nbf,nstate,nr,nocc
 real(dp),intent(in)        :: c_matrix_cmplx(nbf,nocc,nspin)
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: basis_function_r(nbf,nr)
 real(dp),intent(in)        :: basis_function_gradr(nbf,nr,3)
 real(dp),intent(out)       :: rhor(nspin,nr)
 real(dp),intent(out)       :: grad_rhor(nspin,nr,3)
!=====
 integer              :: ispin,istate,ir
 complex(dp),allocatable :: phir_cmplx(:,:)
 complex(dp),allocatable :: phir_gradx_cmplx(:,:)
 complex(dp),allocatable :: phir_grady_cmplx(:,:)
 complex(dp),allocatable :: phir_gradz_cmplx(:,:)
!=====

 !
 ! Calculate rho and grad rho at points in batch
 rhor(:,:)        = 0.0_dp
 grad_rhor(:,:,:) = 0.0_dp

 do ispin=1,nspin

   allocate(phir_cmplx(nocc,nr))
   allocate(phir_gradx_cmplx(nocc,nr))
   allocate(phir_grady_cmplx(nocc,nr))
   allocate(phir_gradz_cmplx(nocc,nr))
   phir_cmplx(:,:)       = MATMUL( TRANSPOSE(c_matrix_cmplx(:,:nocc,ispin)) , basis_function_r(:,:) )
   phir_gradx_cmplx(:,:) = MATMUL( TRANSPOSE(c_matrix_cmplx(:,:nocc,ispin)) , basis_function_gradr(:,:,1) )
   phir_grady_cmplx(:,:) = MATMUL( TRANSPOSE(c_matrix_cmplx(:,:nocc,ispin)) , basis_function_gradr(:,:,2) )
   phir_gradz_cmplx(:,:) = MATMUL( TRANSPOSE(c_matrix_cmplx(:,:nocc,ispin)) , basis_function_gradr(:,:,3) )

   forall(ir=1:nr)
     rhor(ispin,ir) = rhor(ispin,ir) + REAL( SUM( phir_cmplx(:,ir) * CONJG(phir_cmplx(:,ir)) * occupation(:nocc,ispin) ) )
     grad_rhor(ispin,ir,1) = grad_rhor(ispin,ir,1) + REAL( SUM( (phir_cmplx(:,ir)*CONJG(phir_gradx_cmplx(:,ir)) + &
                                          CONJG(phir_cmplx(:,ir))*phir_gradx_cmplx(:,ir) )*occupation(:nocc,ispin ) ) )
     grad_rhor(ispin,ir,1) = grad_rhor(ispin,ir,1) + REAL( SUM( (phir_cmplx(:,ir)*CONJG(phir_grady_cmplx(:,ir)) + &
                                          CONJG(phir_cmplx(:,ir))*phir_grady_cmplx(:,ir) )*occupation(:nocc,ispin ) ) )
     grad_rhor(ispin,ir,1) = grad_rhor(ispin,ir,1) + REAL( SUM( (phir_cmplx(:,ir)*CONJG(phir_gradz_cmplx(:,ir)) + &
                                          CONJG(phir_cmplx(:,ir))*phir_gradz_cmplx(:,ir) )*occupation(:nocc,ispin ) ) )
   endforall

   deallocate(phir_cmplx)
   deallocate(phir_gradx_cmplx,phir_grady_cmplx,phir_gradz_cmplx)


 enddo


end subroutine calc_density_gradr_batch_cmplx

!=========================================================================
subroutine calc_density_gradr_laplr(nspin,nbf,p_matrix,basis_function_r,basis_function_gradr,basis_function_laplr, &
 grad_rhor,tau,lapl_rhor)
 use m_definitions
 use m_mpi
 implicit none
 integer,intent(in)   :: nspin,nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(in)  :: basis_function_r(nbf)
 real(dp),intent(in)  :: basis_function_gradr(3,nbf)
 real(dp),intent(in)  :: basis_function_laplr(3,nbf)
 real(dp),intent(out) :: grad_rhor(3,nspin)
 real(dp),intent(out) :: tau(nspin)
 real(dp),intent(out) :: lapl_rhor(nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 grad_rhor(:,:)=0.0_dp
 tau(:)        =0.0_dp
 lapl_rhor(:)  =0.0_dp
 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf

       grad_rhor(:,ispin) = grad_rhor(:,ispin) + p_matrix(ibf,jbf,ispin) &
            *( basis_function_gradr(:,ibf) * basis_function_r(jbf) &
             + basis_function_gradr(:,jbf) * basis_function_r(ibf) ) 

       tau(ispin)        = tau(ispin)        + p_matrix(ibf,jbf,ispin) &
            * DOT_PRODUCT( basis_function_gradr(:,ibf) , basis_function_gradr(:,jbf) )

       lapl_rhor(ispin)  = lapl_rhor(ispin)  + p_matrix(ibf,jbf,ispin) &
                          * (  SUM( basis_function_laplr(:,ibf) ) * basis_function_r(jbf)               &
                             + basis_function_r(ibf) * SUM( basis_function_laplr(:,jbf) )               &
                             + 2.0_dp * DOT_PRODUCT( basis_function_gradr(:,ibf) , basis_function_gradr(:,jbf) ) )

     enddo
   enddo
 enddo


end subroutine calc_density_gradr_laplr


!=========================================================================
subroutine teter_lda_vxc_exc(nr,rhor,vxc,exc)
 use m_definitions
 implicit none

 integer,intent(in)   :: nr
 real(dp),intent(in)  :: rhor(nr)
 real(dp),intent(out) :: vxc(nr),exc(nr)
!=====
 integer :: ir
 !
 ! The usual full LDA parameters of Teter
 real(dp),parameter :: a0p=.4581652932831429_dp
 real(dp),parameter :: a1p=2.217058676663745_dp
 real(dp),parameter :: a2p=0.7405551735357053_dp
 real(dp),parameter :: a3p=0.01968227878617998_dp
 real(dp),parameter :: b1p=1.0_dp
 real(dp),parameter :: b2p=4.504130959426697_dp
 real(dp),parameter :: b3p=1.110667363742916_dp
 real(dp),parameter :: b4p=0.02359291751427506_dp

 real(dp)           :: dd1drs(nr)
 real(dp)           :: dn1drs(nr)
 real(dp)           :: dexcdrs(nr)
 real(dp)           :: n1(nr)
 real(dp)           :: d1(nr)
 real(dp)           :: rs(nr)
!=====

 rs(:) = ( 3.0_dp / (4.0_dp*pi*rhor(:)) )**(1.0_dp/3.0_dp) 
 n1(:) = a0p + rs(:) * (a1p + rs(:) * ( a2p + rs(:) * a3p ) )
 d1(:) = rs(:) * ( b1p + rs(:) * ( b2p + rs(:) * ( b3p + rs(:) * b4p ) ) )

 ! Firstly, exchange-correlation energy
 exc(:) = -n1(:) / d1(:)

 ! Secondly, exchange-correlation potential
 dn1drs(:) = a1p + rs(:) * ( 2.0_dp * a2p + rs(:) * ( 3.0_dp * a3p ) )
 dd1drs(:) = b1p + rs(:) * ( 2.0_dp * b2p + rs(:) * ( 3.0_dp * b3p + rs(:) * ( 4.0_dp * b4p ) ) )

 ! dexcdrs is d(exc)/d(rs)
 dexcdrs(:) = -( dn1drs(:) + exc(:) * dd1drs(:) ) / d1(:)
 vxc(:) = exc(:) - rs(:) * dexcdrs(:) / 3.0_dp


end subroutine teter_lda_vxc_exc


!=========================================================================
subroutine my_lda_exc_vxc(nspin,ixc,rhor,exc,vxc)
 use m_definitions
 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nspin,ixc
 real(dp),intent(in) :: rhor(nspin)
!arrays
 real(dp),intent(out) :: exc,vxc(nspin)

 real(dp) :: a0p
 real(dp) :: a1p
 real(dp) :: a2p
 real(dp) :: a3p
 real(dp) :: b1p
 real(dp) :: b2p
 real(dp) :: b3p
 real(dp) :: b4p
 real(dp) :: da0
 real(dp) :: da1
 real(dp) :: da2
 real(dp) :: da3
 real(dp) :: db1
 real(dp) :: db2
 real(dp) :: db3
 real(dp) :: db4

 real(dp),parameter :: alpha_zeta=1.0_dp - 1.0e-6_dp
 real(dp),parameter :: ft=4._dp/3._dp,rsfac=0.6203504908994000_dp
 real(dp),parameter :: rsfacm3=rsfac**(-3)
 real(dp) :: a0,a1,a2,a3,b1,b2,b3,b4,d1,d1m1,d2d1drs2,d2d1drsdf,d2excdf2
 real(dp) :: d2excdrs2,d2excdrsdf,d2excdz2,d2fxcdz2,d2n1drs2,d2n1drsdf,dd1df
 real(dp) :: dd1drs,dexcdf,dexcdrs,dexcdz,dfxcdz,dn1df,dn1drs,dvxcdrs
 real(dp) :: dvxcpdrho,dvxcpdz,fact,fxc,n1
 real(dp) :: rhom1,rs,vxcp,zet,zetm,zetm_third
 real(dp) :: zetp,zetp_third

! *************************************************************************

 select case(ixc)
 case(1200)
   !
   ! The Slater exchange
   a0p=3.0_dp/8.0_dp*(18.0_dp/pi**2)**(1.0_dp/3.0_dp)
   a1p=0.0_dp
   a2p=0.0_dp
   a3p=0.0_dp
   b1p=1.0_dp
   b2p=0.0_dp
   b3p=0.0_dp
   b4p=0.0_dp
   da0=0.0_dp
   da1=0.0_dp
   da2=0.0_dp
   da3=0.0_dp
   db1=0.0_dp
   db2=0.0_dp
   db3=0.0_dp
   db4=0.0_dp
 case(1100)
   !
   ! The usual full LDA parameters of Teter
   a0p=.4581652932831429_dp
   a1p=2.217058676663745_dp
   a2p=0.7405551735357053_dp
   a3p=0.01968227878617998_dp
   b1p=1.0_dp
   b2p=4.504130959426697_dp
   b3p=1.110667363742916_dp
   b4p=0.02359291751427506_dp
   da0=.119086804055547_dp
   da1=0.6157402568883345_dp
   da2=0.1574201515892867_dp
   da3=0.003532336663397157_dp
   db1=0.0_dp
   db2=0.2673612973836267_dp
   db3=0.2052004607777787_dp
   db4=0.004200005045691381_dp
 case(1000)
   !
   ! full range RPA 
   a0p= 0.00033959499_dp
   a1p= 0.1912460_dp
   a2p= 0.8008790_dp
   a3p= 0.092956297_dp
   b1p=1.0_dp
   b2p= 8.710930_dp
   b3p= 3.945600_dp
   b4p= 0.087989897_dp
   da0=-0.00015974799_dp
   da1=-0.082753003_dp
   da2=-0.3675560_dp
   da3=-0.044320997_dp
   db1=0.0_dp
   db2=-1.113190_dp
   db3=-1.221470_dp
   db4=-0.040220797_dp
 case(1020)
   !
   ! Long-range only RPA parameters with rc=2.0
   a0p=-0.000012396600_dp
   a1p= 0.0014478000_dp
   a2p= 0.021771301_dp
   a3p= 0.00025572101_dp
   b1p=1.0_dp
   b2p= 0.3820980_dp
   b3p= 0.042663701_dp
   b4p= 0.00010346600_dp
   da0= 0.0000018310002_dp
   da1=-0.00021740992_dp
   da2=-0.0077045010_dp
   da3=-0.00020484751_dp
   db1=0.0_dp
   db2= 0.021046013_dp
   db3=-0.018320801_dp
   db4=-0.00012402580_dp
 case(1010)
   !
   ! Long-range only RPA parameters with rc=1.0
   a0p=-0.000028384500_dp
   a1p= 0.0037404201_dp
   a2p= 0.063176401_dp
   a3p= 0.0023404600_dp
   b1p=1.0_dp
   b2p= 0.8482450_dp
   b3p= 0.1845470_dp
   b4p= 0.0016536200_dp
   da0= 0.0000059325994_dp
   da1=-0.00076668011_dp
   da2=-0.024234399_dp
   da3=-0.0014384059_dp
   db1=0.0_dp
   db2= 0.025729001_dp
   db3=-0.071891010_dp
   db4=-0.0010981541_dp
 case(1005)
   !
   ! Long-range only RPA parameters with rc=0.5
   a0p=-5.8032401E-05
   a1p= 8.9607602E-03
   a2p= 0.1718570
   a3p= 1.3439300E-02
   b1p=1.0_dp
   b2p= 1.849290
   b3p= 0.7096860
   b4p= 1.1346900E-02
   da0= 1.3963599E-05
   da1= -2.1155602E-03
   da2= -7.3816001E-02
   da3= -7.0218993E-03
   db1=0.0_dp
   db2=-7.2860003E-02
   db3=-0.2463360
   db4=-5.8700801E-03
 end select

!Although fact is parameter value, some compilers are not able to evaluate
!it at compile time.
 fact=1.0_dp/(2.0_dp**(4.0_dp/3.0_dp)-2.0_dp)


 if (nspin==1) then

       rs=( 3.0_dp / (4.0_dp*pi*rhor(1)) )**(1.0_dp/3.0_dp) 
       n1=a0p+rs*(a1p+rs*(a2p+rs*a3p))
       d1=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
       d1m1=1.0_dp/d1

!      Exchange-correlation energy
       exc=-n1*d1m1

!      Exchange-correlation potential
       dn1drs=a1p+rs*(2._dp*a2p+rs*(3._dp*a3p))
       dd1drs=b1p+rs*(2._dp*b2p+rs*(3._dp*b3p+rs*(4._dp*b4p)))

!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+exc*dd1drs)*d1m1
       vxc(1)=exc-rs*dexcdrs/3.0_dp

 else 

!    Allow for spin polarization. This part could be optimized for speed.

       rs=( 3.0_dp / (4.0_dp*pi*SUM(rhor(:))) )**(1.0_dp/3.0_dp) 
       zet= ( rhor(1) - rhor(2) ) / SUM( rhor(:) )
       zetp=1.0_dp+zet*alpha_zeta
       zetm=1.0_dp-zet*alpha_zeta
       zetp_third=zetp**(1.0_dp/3.0_dp)
       zetm_third=zetm**(1.0_dp/3.0_dp)
!      Exchange energy spin interpolation function f(zeta)
       fxc=( zetp*zetp_third + zetm*zetm_third - 2.0_dp ) *fact

       a0=a0p+fxc*da0
       a1=a1p+fxc*da1
       a2=a2p+fxc*da2
       a3=a3p+fxc*da3
       b1=b1p+fxc*db1
       b2=b2p+fxc*db2
       b3=b3p+fxc*db3
       b4=b4p+fxc*db4

       n1= a0+rs*(a1+rs*(a2+rs*a3))
       d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
       d1m1=1.0_dp/d1

!      Exchange-correlation energy
       exc=-n1*d1m1

!      Exchange-correlation potential
       dn1drs=a1+rs*(2._dp*a2+rs*(3._dp*a3))
       dd1drs=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+exc*dd1drs)*d1m1

!      Only vxcp contributes when paramagnetic
       vxcp=exc-rs*dexcdrs/3.0_dp

!      d(fxc)/d(zeta)  (which is 0 at zeta=0)
       dfxcdz=ft*alpha_zeta*(zetp_third-zetm_third)*fact

!      dn1df=d(n1)/d(fxc) and dd1df=d(d1)/d(fxc)
       dn1df=da0+rs*(da1+rs*(da2+rs*da3))
       dd1df=rs*(db1+rs*(db2+rs*(db3+rs*db4)))

!      dexcdz is d(exc)/d(zeta)
       dexcdf=-(dn1df+exc*dd1df)*d1m1
       dexcdz=dfxcdz*dexcdf

!      Compute Vxc for both spin channels

       vxc(1)=vxcp - (zet-1.0_dp)*dexcdz
       vxc(2)=vxcp - (zet+1.0_dp)*dexcdz

 end if

end subroutine my_lda_exc_vxc


!=========================================================================
subroutine my_lda_exc_vxc_mu(mu,rhor,exc,vxc)
 use m_definitions
 implicit none

 real(dp),intent(in) :: mu
 integer,parameter :: npt=1
 integer,parameter :: order=1
 real(dp),intent(in) :: rhor(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
!=====
 real(dp)           :: efac,rs,vfac
 real(dp)           :: kf
 real(dp)           :: omega
 real(dp)           :: rho,aa,f_aa
 real(dp)           :: exp4aa2,rhodaadrho,dfdaa
!=====


 efac=0.75_dp*(1.5_dp/pi)**(2.0_dp/3.0_dp)
 vfac=4.0/3.0 * efac

 omega = mu
 rs= (3.0_dp/(4.0_dp*pi*rhor(1)))**(1.0_dp/3.0_dp)
 rho = 3.0 / ( 4.0 * pi * rs**3 )
 kf  = ( 3.0 * pi**2 * rho )**(1.0/3.0)
 aa  = omega / ( 2.0 * kf )

 exp4aa2 = EXP(-1.0/(4.0*aa**2))

 f_aa = 8.0/3.0 * aa &
       * ( SQRT(pi) * ERF(0.5/aa) &
          + (2.0 * aa - 4.0 * aa**3) * exp4aa2 &
          - 3.0 * aa + 4.0 * aa**3 )

 rhodaadrho = - omega/ ( 6.0 * kf )

 dfdaa = f_aa / aa + 8.0 * aa * ( 4.0 * aa**2 * ( 1.0 - exp4aa2 ) - 1.0 )

 exc(1) = -efac/rs * (1.0 - f_aa)

 vxc(1) = -vfac/rs * (1.0 - f_aa) + efac/rs * dfdaa * rhodaadrho


end subroutine my_lda_exc_vxc_mu


!=========================================================================
subroutine my_gga_exc_vxc_hjs(omega,nn,sigma,exc,vxc,vsigma)
 use m_definitions
 implicit none

 real(dp),intent(in)  :: omega,nn,sigma
 real(dp),intent(out) :: exc,vxc,vsigma
!=====
 real(dp),parameter :: ss0=2.0
 ! HJS parameters
 real(dp),parameter :: aabar= 0.757211
 real(dp),parameter :: bb   =-0.106364
 real(dp),parameter :: cc   =-0.118649
 real(dp),parameter :: dd   = 0.609650
 real(dp),parameter :: ee   =-0.0477963
 ! PBE parameters
 real(dp),parameter :: a2= 0.0159941
 real(dp),parameter :: a3= 0.0852995
 real(dp),parameter :: a4=-0.160368
 real(dp),parameter :: a5= 0.152645
 real(dp),parameter :: a6=-0.0971263
 real(dp),parameter :: a7= 0.0422061
 real(dp),parameter :: b1= 5.33319
 real(dp),parameter :: b2=-12.4780
 real(dp),parameter :: b3=11.0988
 real(dp),parameter :: b4=-5.11013
 real(dp),parameter :: b5= 1.71468
 real(dp),parameter :: b6=-0.610380
 real(dp),parameter :: b7= 0.307555
 real(dp),parameter :: b8=-0.0770547
 real(dp),parameter :: b9= 0.0334840
!=====
 real(dp) :: efac
 real(dp) :: nn_local,sigma_local,rs
 real(dp) :: ss
 real(dp) :: kf
 real(dp) :: nu
 real(dp) :: chi
 real(dp) :: lambda,eta,zeta
 real(dp) :: hh_s
 real(dp) :: ffbar_s
 real(dp) :: ggbar_s
 real(dp) :: factor_w
 real(dp) :: exc_nn,exc_sigma
 real(dp) :: fx,dfxds,dfxdnu
 real(dp) :: dsdsigma,dsdn,dnudn
!=====

 efac=0.75_dp * (1.5_dp/pi)**(2.0_dp/3.0_dp)


!HOME MADE    !
!HOME MADE    ! first calculation
!HOME MADE    nn_local = nn
!HOME MADE    sigma_local = sigma
!HOME MADE    
!HOME MADE    rs = ( 3.0 / (4.0 *pi * nn_local) )**(1./3.)
!HOME MADE    kf = (9.0_dp * pi / 4.0_dp)**(1.0_dp/3.0_dp) / rs
!HOME MADE    nu = omega / kf
!HOME MADE    ss = SQRT(sigma_local) / ( 2.0_dp * kf * nn_local )
!HOME MADE   
!HOME MADE    hh_s = ( a2*ss**2 + a3*ss**3 + a4*ss**4 + a5*ss**5 + a6*ss**6 + a7*ss**7 ) &
!HOME MADE         / ( 1.0_dp + b1*ss + b2*ss**2 + b3*ss**3 + b4*ss**4 + b5*ss**5 + b6*ss**6 + b7*ss**7 + b8*ss**8 + b9*ss**9 )
!HOME MADE   
!HOME MADE    ffbar_s = 1.0_dp - ss**2 / ( 27.0_dp * cc * (1.0_dp + ss**2/ss0**2 ) ) &
!HOME MADE               - ss**2 * hh_s / (2.0_dp * cc )
!HOME MADE   
!HOME MADE   
!HOME MADE    zeta   = ss**2 * hh_s
!HOME MADE    eta    = aabar + ss**2 * hh_s
!HOME MADE    lambda = dd    + ss**2 * hh_s
!HOME MADE    chi = nu / SQRT( lambda + nu**2)
!HOME MADE   
!HOME MADE    ggbar_s = -2./5.*cc*ffbar_s*lambda -4./15.*bb*lambda**2 - 6./5.*aabar*lambda**3 &
!HOME MADE             -4./5.*SQRT(pi)*lambda**(7./2.) &
!HOME MADE             -12./5.*lambda**(7./2.) * ( SQRT(zeta)-SQRT(eta) )
!HOME MADE    ggbar_s = ggbar_s / ee
!HOME MADE   
!HOME MADE   
!HOME MADE    factor_w = aabar - 4./9.*bb/lambda*(1.0-chi) - 4./9.*cc*ffbar_s/lambda**2 * (1.0 - 1.5*chi+0.5*chi**3)  &
!HOME MADE              -8./9.*ee*ggbar_s/lambda**3 * ( 1.0 - 15./8.*chi + 5./4.*chi**3 -3./8.*chi**5 ) &
!HOME MADE              + 2.*nu   * ( SQRT(zeta+nu**2)- SQRT(eta+nu**2) ) &
!HOME MADE              + 2.*zeta * LOG( ( nu + SQRT(zeta+nu**2) ) / ( nu + SQRT(lambda + nu**2) ) ) &
!HOME MADE              - 2.*eta  * LOG( ( nu + SQRT( eta+nu**2) ) / ( nu + SQRT(lambda + nu**2) ) ) 
!HOME MADE   
!HOME MADE   
!HOME MADE    exc = -efac/rs * factor_w
!HOME MADE
!HOME MADE write(stdout,*) 'exc1=',exc

 !
 ! call to the nwchem subroutine
 !
 rs = ( 3.0 / (4.0*pi*nn) )**(1.0/3.0)
 kf = (9.0_dp * pi / 4.0_dp)**(1.0_dp/3.0_dp) / rs
 ss = SQRT(sigma) / ( 2.0_dp * kf * nn )

 call HSE08Fx(omega,1,nn,ss,fx,dfxds,dfxdnu)

 exc = -efac/rs*fx

 dsdsigma= 1.0_dp / ( 4.0_dp * kf * nn * SQRT(sigma) )

 vsigma = -efac/rs * nn * dfxds * dsdsigma

 dsdn  = SQRT(sigma) / (2.0_dp * (3.0*pi**2)**(1./3.) ) * (-4.0/3.0) * nn**(-7.0/3.0)
 dnudn = omega / (3.0*pi**2)**(1./3.) * (-1.0/3.0) * nn**(-4.0/3.0)
 vxc = -efac/rs * nn * ( dfxds * dsdn + dfxdnu * dnudn ) - (4.0/3.0)*efac/rs *fx


end subroutine my_gga_exc_vxc_hjs


!=========================================================================
subroutine HSE08Fx(omega,ipol,rho,s,Fxhse,d10Fxhse,d01Fxhse)
 use m_definitions
 implicit none

! HSE evaluates the Heyd et al. Screened Coulomb
! Exchange Functional
!
! Calculates the enhancement factor
!
 integer  :: ipol
 real(dp) :: omega
 real(dp) :: rho,s,Fxhse,d10Fxhse,d01Fxhse
 real(dp) :: A,B,C,D,E
 real(dp) :: ha2,ha3,ha4,ha5,ha6,ha7
 real(dp) :: hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9
 real(dp) :: zero,one,two,three,four,five,six,seven,eight
 real(dp) :: nine,ten
 real(dp) :: H,hnum,hden 
 real(dp) :: d1H,d1hnum,d1hden 
 real(dp) :: s2,s3,s4,s5,s6,s7,s8,s9
 real(dp) :: Fs, d1Fs
 real(dp) :: zeta, lambda, eta, kf, nu, chi, lambda2
 real(dp) :: d1zeta,d1lambda,d1eta,d1nu,d1chi,d1lambda2
 real(dp) :: EGs,d1EGs
 real(dp) :: nu2,L2,L3,nu3,nu4,nu5,nu6
 real(dp) :: Js,Ks,Ms,Ns
 real(dp) :: d1Js,d1Ks,d1Ms,d1Ns
 real(dp) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
 real(dp) :: tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15
 real(dp) :: Fxhse1,Fxhse2,Fxhse3,Fxhse4,Fxhse5,Fxhse6
 real(dp) :: d1Fxhse1,d1Fxhse2,d1Fxhse3,d1Fxhse4,d1Fxhse5
 real(dp) :: d1Fxhse6,d1Fxhse7
 real(dp) :: r42,r27,r12,r15,r14,r18,r20,r30,r56,r72
 real(dp) :: r16,r32,r24,r48,r11,r64,r35
 real(dp) :: srpi,s02
 real(dp) :: f12,f13,f32,f52,f72,f92

!
!Constants for HJS hole
!
 Data A,B,C,D,E  &
     / 7.57211D-1,-1.06364D-1,-1.18649D-1, &
       6.09650D-1,-4.77963D-2 /
!
!Constants for fit of H(s) (PBE hole)
!Taken from JCTC_5_754 (2009)
!
 Data ha2,ha3,ha4,ha5,ha6,ha7 &
     / 1.59941D-2,8.52995D-2,-1.60368D-1,1.52645D-1, &
      -9.71263D-2,4.22061D-2 /

 Data hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9 &
      / 5.33319D0,-12.4780D0,11.0988D0,-5.11013D0,&
       1.71468D0,-6.10380D-1,3.07555D-1,-7.70547D-2,&
       3.34840D-2 /

!
!Whole numbers used during evaluation
!
 Data zero,one,two,three,four,five,six,seven,eight,nine,ten &
      / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /
  
 Data r11,r12,r14,r15,r16,r18,r20,r24,r27,r30,r32 &
      / 11D0,12D0,14D0,15D0,16D0,18D0,20D0,24D0,27d0,30D0,32D0 /

 Data r35,r42,r48,r56,r64,r72 &
      / 35D0,42D0,48D0,56D0,64D0,72D0 /
!
!Fractions used during evaluation
!
 Data f12     / 0.5D0 /
!
!General constants
!
 f13   = one/three
 f32   = three/two
 f52   = five/two
 f72   = seven/two
 f92   = nine/two
 srpi = dsqrt(pi)
!
!
!Calculate prelim variables
!
 s2 = s*s
 s02 = s2/four
 s3 = s2*s
 s4 = s3*s
 s5 = s4*s
 s6 = s5*s
 s7 = s6*s
 s8 = s7*s
 s9 = s8*s

!
!Calculate H(s) the model exhange hole
!
 hnum = ha2*s2 + ha3*s3 + ha4*s4 + ha5*s5 + ha6*s6 + ha7*s7 
 hden = one + hb1*s + hb2*s2 + hb3*s3 + hb4*s4 + hb5*s5 + &
        hb6*s6 + hb7*s7 + hb8*s8 + hb9*s9
 H = hnum/hden

!
!Calculate helper variables
!
 zeta = s2*H
 eta = A + zeta
 lambda = D + zeta
 if (ipol.eq.1) then
    kf = (three*pi2*rho)**f13 
 else
    kf = (six*pi2*rho)**f13 
 endif
 nu = omega/kf
 chi = nu/dsqrt(lambda+nu**two)
 lambda2 = (one+chi)*(lambda+nu**two)

!
!Calculate F(H(s)) for the model exhange hole
!
 Fs = one-s2/(r27*C*(one+s02))-zeta/(two*C)

!
!Calculate EG(s) 
!
 EGs = -(two/five)*C*Fs*lambda - (four/r15)*B*lambda**two - &
       (six/five)*A*lambda**three - &
       (four/five)*srpi*lambda**(seven/two) -&
       (r12/five)*(lambda**(seven/two))*(dsqrt(zeta)-dsqrt(eta))
 
!
!Calculate the denominators needed
!

 nu2 = nu*nu
 Js = (dsqrt(zeta+nu2)+dsqrt(eta+nu2))*(dsqrt(zeta+nu2)+nu) 
 Ks = (dsqrt(zeta+nu2)+dsqrt(eta+nu2))*(dsqrt(eta+nu2)+nu) 
 Ms = (dsqrt(zeta+nu2)+dsqrt(lambda+nu2))*(dsqrt(lambda+nu2)+nu) 
 Ns = (dsqrt(eta+nu2)+dsqrt(lambda+nu2))*(dsqrt(lambda+nu2)+nu) 

!
!  The final value for the enhancement factor is
!
 tmp1 = one + f12*chi
 tmp2 = one + (nine/eight)*chi + (three/eight)*chi**two 
 Fxhse1  = A*(zeta/Js + eta/Ks) 
 Fxhse2  = -(four/nine)*B/lambda2
 Fxhse3  = -(four/nine)*C*Fs*tmp1/lambda2**two
 Fxhse4  = -(eight/nine)*EGs*tmp2/lambda2**three
 Fxhse5  = two*zeta*dlog(one -D/Ms)
 Fxhse6  = -two*eta*dlog(one -(D-A)/Ns)

 Fxhse = Fxhse1+Fxhse2+Fxhse3+Fxhse4+Fxhse5+Fxhse6
!
!Calculate the first derivative of H with respect to the
!reduced density gradient, s.
!
 d1hnum = two*ha2*s + three*ha3*s2 + four*ha4*s3 + &
           five*ha5*s4 + six*ha6*s5 + seven*ha7*s6

 d1hden  = hb1 + two*hb2*s +three*hb3*s2 + four*hb4*s3 + &
           five*hb5*s4 + six*hb6*s5 + seven*hb7*s6 +&
           eight*hb8*s7 + nine*hb9*s8 
 d1H =   (hden*d1hnum -hnum*d1hden)/hden**two

!
!calculate first derivative of variables needed with respect to s
!
 d1zeta = two*s*H + s2*d1H
 d1eta  = d1zeta
 d1lambda = d1zeta
 d1chi = -f12*nu*d1zeta/(lambda + nu2)**f32
 d1lambda2 = d1chi*(lambda + nu**two) + (one+chi)*d1lambda
 !d1lambda2 = (d1lambda*(one-chi)+lambda*d1chi)/(one-chi)**two

!
!calculate the first derivative of Fs with respect to s
!
 d1Fs = -two*s/(r27*C*(one+s02)**two) - d1zeta/(two*C)

!
!Calculate the first derivate of EGs with respect to s
!
 d1EGs = -(two/five)*C*(d1Fs*lambda + Fs*d1lambda) -&
         (eight/r15)*B*lambda*d1lambda -&
         (r18/five)*A*lambda*lambda*d1lambda -&
         (r14/five)*srpi*d1lambda*lambda**f52 -&
         (r42/five)*(lambda**f52)*&
         d1lambda*(dsqrt(zeta)-dsqrt(eta))-&
         (six/five)*(lambda**(seven/two))*&
         (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))

!
!Calculate the first derivate of denominators needed with respect
!to s
!
 tmp1 = (dsqrt(zeta+nu2)+nu)/(dsqrt(eta+nu2)) 
 tmp2 = (dsqrt(eta+nu2)+nu)/(dsqrt(zeta+nu2))

 d1Js = f12*d1zeta*(two+tmp1+tmp2)
 d1Ks = d1Js

 tmp3 = (dsqrt(zeta+nu2)+nu)/(dsqrt(lambda+nu2))
 tmp4 = (dsqrt(lambda+nu2)+nu)/(dsqrt(zeta+nu2)) 
 d1Ms = f12*d1zeta*(two +tmp3+tmp4)

 tmp5 = (dsqrt(lambda+nu2)+nu)/(dsqrt(eta+nu2))
 tmp6 = (dsqrt(eta+nu2)+nu)/(dsqrt(lambda+nu2))
 d1Ns = f12*d1zeta*(two + tmp5+tmp6)
!
!Calculate the derivative of the 08-Fxhse with respect to s
!
 L2 = lambda2*lambda2
 L3 = lambda2*lambda2*lambda2
 d1Fxhse1  = A*( (Js*d1zeta - zeta*d1Js)/(Js*Js) +&
                 (Ks*d1zeta - eta*d1Ks)/(Ks*Ks) ) 

 d1Fxhse2  = (four/nine)*B*d1lambda2/L2 

 tmp9 = d1lambda2/lambda2
 tmp7 = d1Fs - two*Fs*tmp9
 tmp8 = one + f12*chi
 tmp10 =  f12*Fs*d1chi

 d1Fxhse3 = -(four*C/(nine*L2))*(tmp7*tmp8+tmp10)


   tmp7 = one + (nine/eight)*chi+(three/eight)*chi*chi
   tmp8 = (nine/eight)*d1chi + (six/eight)*chi*d1chi

  d1Fxhse4 = -(eight/(nine*L3))*((d1EGs-three*EGs*tmp9)*tmp7 &
            + EGs*tmp8)
 d1Fxhse5  = two*d1zeta*dlog(one-D/Ms) + &
            two*zeta*D*d1Ms/(Ms*Ms*(one-D/Ms)) 

 d1Fxhse6  = -two*d1eta*dlog(one- (D-A)/Ns) - &
            two*eta*(D-A)*d1Ns/(Ns*Ns*(one-(D-A)/Ns)) 

 d10Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6
!
!Calculate the derivative of 08-Fxhse with respect to nu
!
 nu3 = nu2*nu

 d1Fxhse1 = -((A*(nu + dsqrt(eta + nu2))*zeta)/ &
             (dsqrt(eta + nu2)*dsqrt(nu2 + zeta)*&
             (nu + dsqrt(nu2 + zeta))*&
             (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))))

 d1Fxhse2 = -((A*eta*(nu/dsqrt(eta + nu2) + nu/ &
             dsqrt(nu2 + zeta)))/&
             ((nu + dsqrt(eta + nu2))*&
             (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)) -&
             (A*eta*(one + nu/dsqrt(eta + nu2)))/&
             ((nu + dsqrt(eta + nu2))**two*&
             (dsqrt(eta + nu2) + dsqrt(nu2 + zeta)))

 d1Fxhse3 = (four*B)/(nine*(lambda + nu2)**(f32))

 d1Fxhse4 = (two*C*Fs)/(three*(lambda + nu2)**(f52))

 d1Fxhse5 = (five*EGs*(lambda**two + four*nu3* &
             (nu + dsqrt(lambda + nu2)) +&
             lambda*nu*(five*nu + three*dsqrt(lambda + nu2))))/&
    (three*(lambda + nu2)**four*(nu + dsqrt(lambda + nu2))**three)

 d1Fxhse6 = (two*D*zeta*(nu + dsqrt(nu2 + zeta)))/&
             (dsqrt(lambda + nu2)*dsqrt(nu2 + zeta)*&
             (-D + lambda + (nu + dsqrt(lambda + nu2))*&
             (nu + dsqrt(nu2 + zeta))))

 d1Fxhse7 = (two*(A - D)*eta*(nu + dsqrt(eta + nu2)))/ &
             (dsqrt(eta + nu2)*dsqrt(lambda + nu2)*&
             (A - D + lambda + nu2 + nu*dsqrt(eta + nu2) +&
             nu*dsqrt(lambda + nu2) +&
             dsqrt(eta + nu2)*dsqrt(lambda + nu2)))


 d01Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6+d1Fxhse7
 
end subroutine HSE08Fx

!=========================================================================
