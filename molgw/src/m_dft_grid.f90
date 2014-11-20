!=========================================================================
#include "macros.h"
!=========================================================================
module m_dft_grid
 use m_definitions
 use m_mpi
 use m_inputparam,only: quadrature_name
 
 real(dp),parameter :: shift=1.e-6_dp ! bohr  some shift used
                                      ! to evaluate numerically the divergence of the gradient
 !
 ! Grid definition
 integer,protected    :: ngrid
 integer,protected    :: nradial
 integer,protected    :: nangular

 real(dp),allocatable :: rr_grid(:,:)
 real(dp),allocatable :: w_grid(:)

 !
 ! Function evaluation storage
 integer,parameter    :: ngrid_max_stored=5000
 integer              :: ngrid_stored=0
 real(dp),allocatable :: bfr(:,:)
 real(dp),allocatable :: bfr_x(:,:)
 real(dp),allocatable :: bfr_y(:,:)
 real(dp),allocatable :: bfr_z(:,:)
 real(dp),allocatable :: bfgr(:,:,:)
 real(dp),allocatable :: bfgr_x(:,:,:)
 real(dp),allocatable :: bfgr_y(:,:,:)
 real(dp),allocatable :: bfgr_z(:,:,:)
 real(dp),allocatable :: bflr(:,:,:)


contains


!=========================================================================
subroutine setup_dft_grid()
 use m_atoms
 use m_tools,only: coeffs_gausslegint
 implicit none

 integer              :: iradial,iatom,iangular,ir,igrid
 integer              :: n1
 real(dp)             :: weight
 real(dp),allocatable :: x1(:)
 real(dp),allocatable :: y1(:)
 real(dp),allocatable :: z1(:)
 real(dp),allocatable :: w1(:)
 real(dp),allocatable :: xa(:),wxa(:)
 real(dp)             :: p_becke(natom),s_becke(natom,natom),fact_becke 
 real(dp)             :: mu
 integer              :: jatom,katom
!=====

 select case(TRIM(quadrature_name))
 case('very low')
   nradial  =  6
   nangular =  6
 case('low')
   nradial  = 10
   nangular = 14
 case('medium')
   nradial  = 20
   nangular = 26
 case('high')
   nradial  = 40
   nangular = 50
 case('very high')
   nradial  = 80
   nangular = 86
 case default
   stop'integration quality not recognized'
 end select

 ! Total number of grid points
 ngrid = natom * nradial * nangular

 call init_grid_distribution(ngrid)

 WRITE_MASTER(*,'(/,a)')       ' Setup the DFT quadrature'
 WRITE_MASTER(*,'(a,i4,x,i4)') ' discretization grid per atom [radial points , angular points] ',nradial,nangular
 WRITE_MASTER(*,'(a,i8)')      ' total number of real-space points for this processor',ngrid

 allocate(xa(nradial),wxa(nradial))
 allocate(x1(nangular),y1(nangular),z1(nangular),w1(nangular))

 !
 ! spherical integration
 ! radial part with Gauss-Legendre
 call coeffs_gausslegint(-1.0_dp,1.0_dp,xa,wxa,nradial)
 !
 ! Transformation from [-1;1] to [0;+\infty[
 ! taken from M. Krack JCP 1998
 wxa(:) = wxa(:) * ( 1.0_dp / log(2.0_dp) / ( 1.0_dp - xa(:) ) )
 xa(:)  = log( 2.0_dp / (1.0_dp - xa(:) ) ) / log(2.0_dp)

 n1 = nangular
 ! angular part with Lebedev - Laikov
 ! (x1,y1,z1) on the unit sphere
 select case(nangular)
 case(6)
   call ld0006(x1,y1,z1,w1,n1)
 case(14)
   call ld0014(x1,y1,z1,w1,n1)
 case(26)
   call ld0026(x1,y1,z1,w1,n1)
 case(38)
   call ld0038(x1,y1,z1,w1,n1)
 case(50)
   call ld0050(x1,y1,z1,w1,n1)
 case(74)
   call ld0074(x1,y1,z1,w1,n1)
 case(86)
   call ld0086(x1,y1,z1,w1,n1)
 case default
   WRITE_MASTER(*,*) 'grid points: ',nangular
   stop'Lebedev grid is not available'
 end select
 
 allocate(rr_grid(3,ngrid),w_grid(ngrid))

 igrid = 0
 ir    = 0
 do iatom=1,natom
   do iradial=1,nradial
     do iangular=1,nangular
       igrid = igrid + 1
       if( .NOT. is_my_grid_task(igrid) ) cycle
       ir = ir + 1

       rr_grid(1,ir) = xa(iradial) * x1(iangular) + x(1,iatom)
       rr_grid(2,ir) = xa(iradial) * y1(iangular) + x(2,iatom)
       rr_grid(3,ir) = xa(iradial) * z1(iangular) + x(3,iatom)

       weight   = wxa(iradial) * w1(iangular) * xa(iradial)**2 * 4.0_dp * pi

       !
       ! Partitionning scheme of Axel Becke, J. Chem. Phys. 88, 2547 (1988).
       !
       s_becke(:,:) = 0.0_dp
       do katom=1,natom
         do jatom=1,natom
           if(katom==jatom) cycle
           mu = ( SQRT( SUM( (rr_grid(:,ir)-x(:,katom) )**2 ) ) - SQRT( SUM( (rr_grid(:,ir)-x(:,jatom) )**2 ) ) ) &
                     / SQRT( SUM( (x(:,katom)-x(:,jatom))**2 ) )
           s_becke(katom,jatom) = 0.5_dp * ( 1.0_dp - smooth_step(smooth_step(smooth_step(mu))) )
         enddo
       enddo
       p_becke(:) = 1.0_dp
       do katom=1,natom
         do jatom=1,natom
           if(katom==jatom) cycle
           p_becke(katom) = p_becke(katom) * s_becke(katom,jatom)
         enddo
       enddo
       fact_becke = p_becke(iatom) / SUM( p_becke(:) )

       w_grid(ir) = weight * fact_becke

     enddo
   enddo
 enddo


end subroutine setup_dft_grid


!=========================================================================
subroutine destroy_dft_grid()
 implicit none

 deallocate(rr_grid)
 deallocate(w_grid)

 if( allocated(bfr   ) ) deallocate(bfr)
 if( allocated(bfr_x ) ) deallocate(bfr_x)
 if( allocated(bfr_y ) ) deallocate(bfr_y)
 if( allocated(bfr_z ) ) deallocate(bfr_z)
 if( allocated(bfgr  ) ) deallocate(bfgr)
 if( allocated(bfgr_x) ) deallocate(bfgr_x)
 if( allocated(bfgr_y) ) deallocate(bfgr_y)
 if( allocated(bfgr_z) ) deallocate(bfgr_z)
 if( allocated(bflr  ) ) deallocate(bflr)

end subroutine destroy_dft_grid


!=========================================================================
function smooth_step(mu)
 real(dp) :: smooth_step
 real(dp),intent(in) :: mu
!=====

 smooth_step = 0.5_dp * mu * ( 3.0_dp - mu**2 )

end function smooth_step


!=========================================================================
subroutine prepare_basis_functions_r(basis)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer                    :: igrid
 real(dp)                   :: rr(3)
 real(dp)                   :: basis_function_r(basis%nbf)
!=====

 ngrid_stored = MIN(ngrid,ngrid_max_stored)
 WRITE_MASTER(*,*) 'Precalculate the functions on N grid points',ngrid_stored
 WRITE_MASTER(*,'(a,2x,f14.2)') ' corresponding to [Mb]:',REAL(basis%nbf,dp)*REAL(ngrid_stored,dp)*REAL(dp,dp)/REAL(1024,dp)**2

 allocate(bfr(basis%nbf,ngrid_stored))

 do igrid=1,ngrid_stored
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_r(basis,rr,basis_function_r)
   bfr(:,igrid) = basis_function_r(:)
 enddo

end subroutine prepare_basis_functions_r


!=========================================================================
subroutine prepare_basis_functions_gradr(basis)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer                    :: igrid
 real(dp)                   :: rr(3)
 real(dp)                   :: basis_function_r_shiftx    (basis%nbf)
 real(dp)                   :: basis_function_r_shifty    (basis%nbf)
 real(dp)                   :: basis_function_r_shiftz    (basis%nbf)
 real(dp)                   :: basis_function_gradr       (3,basis%nbf)
 real(dp)                   :: basis_function_gradr_shiftx(3,basis%nbf)
 real(dp)                   :: basis_function_gradr_shifty(3,basis%nbf)
 real(dp)                   :: basis_function_gradr_shiftz(3,basis%nbf)
!=====

 WRITE_MASTER(*,*) 'Precalculate the gradients on N grid points',ngrid_stored
 WRITE_MASTER(*,'(a,2x,f14.2)') ' corresponding to [Mb]:',REAL(15,dp)*REAL(basis%nbf,dp)*REAL(ngrid_stored,dp)*REAL(dp,dp)/REAL(1024,dp)**2

 allocate(bfr_x (basis%nbf,ngrid_stored))
 allocate(bfr_y (basis%nbf,ngrid_stored))
 allocate(bfr_z (basis%nbf,ngrid_stored))
 allocate(bfgr  (3,basis%nbf,ngrid_stored))
 allocate(bfgr_x(3,basis%nbf,ngrid_stored))
 allocate(bfgr_y(3,basis%nbf,ngrid_stored))
 allocate(bfgr_z(3,basis%nbf,ngrid_stored))

 do igrid=1,ngrid_stored
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_gradr(basis,rr,&
                                        basis_function_gradr,basis_function_r_shiftx,basis_function_r_shifty,basis_function_r_shiftz,&
                                        basis_function_gradr_shiftx,basis_function_gradr_shifty,basis_function_gradr_shiftz)
   bfr_x (:,igrid)   = basis_function_r_shiftx    (:)
   bfr_y (:,igrid)   = basis_function_r_shifty    (:)
   bfr_z (:,igrid)   = basis_function_r_shiftz    (:)
   bfgr  (:,:,igrid) = basis_function_gradr       (:,:)
   bfgr_x(:,:,igrid) = basis_function_gradr_shiftx(:,:)
   bfgr_y(:,:,igrid) = basis_function_gradr_shifty(:,:)
   bfgr_z(:,:,igrid) = basis_function_gradr_shiftz(:,:)
 enddo

end subroutine prepare_basis_functions_gradr


!=========================================================================
subroutine prepare_basis_functions_laplr(basis)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer                    :: igrid
 real(dp)                   :: rr(3)
 real(dp)                   :: basis_function_gradr(3,basis%nbf)
 real(dp)                   :: basis_function_laplr(3,basis%nbf)
!=====

 WRITE_MASTER(*,*) 'Precalculate the laplacians on N grid points',ngrid_stored
 WRITE_MASTER(*,'(a,2x,f14.2)') ' corresponding to [Mb]:',REAL(6,dp)*REAL(basis%nbf,dp)*REAL(ngrid_stored,dp)*REAL(dp,dp)/REAL(1024,dp)**2

 allocate(bfgr(3,basis%nbf,ngrid_stored))
 allocate(bflr(3,basis%nbf,ngrid_stored))

 do igrid=1,ngrid_stored
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_laplr(basis,rr,basis_function_gradr,basis_function_laplr)
   bfgr(:,:,igrid) = basis_function_gradr(:,:)
   bflr(:,:,igrid) = basis_function_laplr(:,:)
 enddo

end subroutine prepare_basis_functions_laplr


!=========================================================================
subroutine get_basis_functions_r(basis,igrid,basis_function_r)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: igrid
 real(dp),intent(out)       :: basis_function_r(basis%nbf)
!=====
 real(dp)                   :: rr(3)
!=====

 if( igrid <= ngrid_stored ) then
   basis_function_r(:) = bfr(:,igrid)
 else
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_r(basis,rr,basis_function_r)
 endif

end subroutine get_basis_functions_r


!=========================================================================
subroutine get_basis_functions_gradr(basis,igrid,basis_function_gradr,&
 basis_function_r_shiftx,basis_function_r_shifty,basis_function_r_shiftz,&
 basis_function_gradr_shiftx,basis_function_gradr_shifty,basis_function_gradr_shiftz)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: igrid
 real(dp),intent(out)       :: basis_function_r_shiftx    (basis%nbf)
 real(dp),intent(out)       :: basis_function_r_shifty    (basis%nbf)
 real(dp),intent(out)       :: basis_function_r_shiftz    (basis%nbf)
 real(dp),intent(out)       :: basis_function_gradr       (3,basis%nbf)
 real(dp),intent(out)       :: basis_function_gradr_shiftx(3,basis%nbf)
 real(dp),intent(out)       :: basis_function_gradr_shifty(3,basis%nbf)
 real(dp),intent(out)       :: basis_function_gradr_shiftz(3,basis%nbf)
!=====
 real(dp)                   :: rr(3)
!=====

 if( igrid <= ngrid_stored ) then
   basis_function_r_shiftx    (:)   = bfr_x (:,igrid)   
   basis_function_r_shifty    (:)   = bfr_y (:,igrid)   
   basis_function_r_shiftz    (:)   = bfr_z (:,igrid)   
   basis_function_gradr       (:,:) = bfgr  (:,:,igrid) 
   basis_function_gradr_shiftx(:,:) = bfgr_x(:,:,igrid) 
   basis_function_gradr_shifty(:,:) = bfgr_y(:,:,igrid) 
   basis_function_gradr_shiftz(:,:) = bfgr_z(:,:,igrid) 
 else
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_gradr(basis,rr,&
                                        basis_function_gradr,basis_function_r_shiftx,basis_function_r_shifty,basis_function_r_shiftz,&
                                        basis_function_gradr_shiftx,basis_function_gradr_shifty,basis_function_gradr_shiftz)
 endif

end subroutine get_basis_functions_gradr


!=========================================================================
subroutine get_basis_functions_laplr(basis,igrid,basis_function_gradr,basis_function_laplr)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: igrid
 real(dp),intent(out)       :: basis_function_gradr(3,basis%nbf)
 real(dp),intent(out)       :: basis_function_laplr(3,basis%nbf)
!=====
 real(dp)                   :: rr(3)
!=====

 if( igrid <= ngrid_stored ) then
   basis_function_gradr(:,:) = bfgr(:,:,igrid) 
   basis_function_laplr(:,:) = bflr(:,:,igrid) 
 else
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_laplr(basis,rr,basis_function_gradr,basis_function_laplr)
 endif

end subroutine get_basis_functions_laplr


!=========================================================================
subroutine calculate_basis_functions_r(basis,rr,basis_function_r)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: rr(3)
 real(dp),intent(out)       :: basis_function_r(basis%nbf)
!=====
 integer              :: ibf,ibf_cart,i_cart
 integer              :: ni,ni_cart,li
 real(dp),allocatable :: basis_function_r_cart(:)
!=====


 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am(CARTESIAN,li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   allocate(basis_function_r_cart(ni_cart))

   do i_cart=1,ni_cart
     basis_function_r_cart(i_cart) = eval_basis_function(basis%bf(ibf_cart+i_cart-1),rr)
   enddo
   basis_function_r(ibf:ibf+ni-1) = MATMUL(  basis_function_r_cart(:) , cart_to_pure(li)%matrix(:,:) )
   deallocate(basis_function_r_cart)

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo


end subroutine calculate_basis_functions_r


!=========================================================================
subroutine calculate_basis_functions_gradr(basis,rr,&
 basis_function_gradr,basis_function_r_shiftx,basis_function_r_shifty,basis_function_r_shiftz,&
 basis_function_gradr_shiftx,basis_function_gradr_shifty,basis_function_gradr_shiftz)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: rr(3)
 real(dp),intent(out)       :: basis_function_gradr(3,basis%nbf)
 real(dp),intent(out)       :: basis_function_r_shiftx(basis%nbf)
 real(dp),intent(out)       :: basis_function_r_shifty(basis%nbf)
 real(dp),intent(out)       :: basis_function_r_shiftz(basis%nbf)
 real(dp),intent(out)       :: basis_function_gradr_shiftx(3,basis%nbf)
 real(dp),intent(out)       :: basis_function_gradr_shifty(3,basis%nbf)
 real(dp),intent(out)       :: basis_function_gradr_shiftz(3,basis%nbf)
!=====
 integer              :: ibf,ibf_cart,i_cart
 integer              :: ni,ni_cart,li
 real(dp)             :: rr_shift(3)
 real(dp),allocatable :: basis_function_r_shiftx_cart(:)
 real(dp),allocatable :: basis_function_r_shifty_cart(:)
 real(dp),allocatable :: basis_function_r_shiftz_cart(:)
 real(dp),allocatable :: basis_function_gradr_cart(:,:)
 real(dp),allocatable :: basis_function_gradr_shiftx_cart(:,:)
 real(dp),allocatable :: basis_function_gradr_shifty_cart(:,:)
 real(dp),allocatable :: basis_function_gradr_shiftz_cart(:,:)
!=====


 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am(CARTESIAN,li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   allocate(basis_function_r_shiftx_cart    (ni_cart))
   allocate(basis_function_r_shifty_cart    (ni_cart))
   allocate(basis_function_r_shiftz_cart    (ni_cart))
   allocate(basis_function_gradr_cart       (3,ni_cart))
   allocate(basis_function_gradr_shiftx_cart(3,ni_cart))
   allocate(basis_function_gradr_shifty_cart(3,ni_cart))
   allocate(basis_function_gradr_shiftz_cart(3,ni_cart))

   do i_cart=1,ni_cart

     basis_function_gradr_cart(:,i_cart)        = eval_basis_function_grad(basis%bf(ibf_cart+i_cart-1),rr)

     rr_shift(:) = rr(:) + (/ shift , 0.0_dp , 0.0_dp /)
     basis_function_r_shiftx_cart(i_cart)       = eval_basis_function(basis%bf(ibf_cart+i_cart-1),rr_shift)
     basis_function_gradr_shiftx_cart(:,i_cart) = eval_basis_function_grad(basis%bf(ibf_cart+i_cart-1),rr_shift)

     rr_shift(:) = rr(:) + (/ 0.0_dp , shift , 0.0_dp /)
     basis_function_r_shifty_cart(i_cart)       = eval_basis_function(basis%bf(ibf_cart+i_cart-1),rr_shift)
     basis_function_gradr_shifty_cart(:,i_cart) = eval_basis_function_grad(basis%bf(ibf_cart+i_cart-1),rr_shift)

     rr_shift(:) = rr(:) + (/ 0.0_dp , 0.0_dp , shift /)
     basis_function_r_shiftz_cart(i_cart)       = eval_basis_function(basis%bf(ibf_cart+i_cart-1),rr_shift)
     basis_function_gradr_shiftz_cart(:,i_cart) = eval_basis_function_grad(basis%bf(ibf_cart+i_cart-1),rr_shift)


   enddo

   basis_function_gradr       (:,ibf:ibf+ni-1) = MATMUL(  basis_function_gradr_cart       (:,:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_gradr_shiftx(:,ibf:ibf+ni-1) = MATMUL(  basis_function_gradr_shiftx_cart(:,:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_gradr_shifty(:,ibf:ibf+ni-1) = MATMUL(  basis_function_gradr_shifty_cart(:,:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_gradr_shiftz(:,ibf:ibf+ni-1) = MATMUL(  basis_function_gradr_shiftz_cart(:,:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_r_shiftx(ibf:ibf+ni-1) = MATMUL(  basis_function_r_shiftx_cart(:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_r_shifty(ibf:ibf+ni-1) = MATMUL(  basis_function_r_shifty_cart(:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_r_shiftz(ibf:ibf+ni-1) = MATMUL(  basis_function_r_shiftz_cart(:) , cart_to_pure(li)%matrix(:,:) )

   deallocate(basis_function_gradr_cart)
   deallocate(basis_function_gradr_shiftx_cart,basis_function_gradr_shifty_cart,basis_function_gradr_shiftz_cart)
   deallocate(basis_function_r_shiftx_cart,basis_function_r_shifty_cart,basis_function_r_shiftz_cart)

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo


end subroutine calculate_basis_functions_gradr


!=========================================================================
subroutine calculate_basis_functions_laplr(basis,rr,basis_function_gradr,basis_function_laplr)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: rr(3)
 real(dp),intent(out)       :: basis_function_gradr(3,basis%nbf)
 real(dp),intent(out)       :: basis_function_laplr(3,basis%nbf)
!=====
 integer              :: ibf,ibf_cart,i_cart
 integer              :: ni,ni_cart,li
 real(dp),allocatable :: basis_function_gradr_cart(:,:)
 real(dp),allocatable :: basis_function_laplr_cart(:,:)
!=====


 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am(CARTESIAN,li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   allocate(basis_function_gradr_cart       (3,ni_cart))
   allocate(basis_function_laplr_cart       (3,ni_cart))

   do i_cart=1,ni_cart

     basis_function_gradr_cart(:,i_cart)        = eval_basis_function_grad(basis%bf(ibf_cart+i_cart-1),rr)
     basis_function_laplr_cart(:,i_cart)        = eval_basis_function_lapl(basis%bf(ibf_cart+i_cart-1),rr)

   enddo

   basis_function_gradr       (:,ibf:ibf+ni-1) = MATMUL(  basis_function_gradr_cart       (:,:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_laplr       (:,ibf:ibf+ni-1) = MATMUL(  basis_function_laplr_cart       (:,:) , cart_to_pure(li)%matrix(:,:) )
   deallocate(basis_function_gradr_cart,basis_function_laplr_cart)

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo


end subroutine calculate_basis_functions_laplr


!=========================================================================
end module m_dft_grid
!=========================================================================
