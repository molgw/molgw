!=========================================================================
#include "macros.h"
!=========================================================================
module m_dft_grid
 use m_definitions
 use m_mpi
 use m_memory
 use m_inputparam,only: grid_level
 
 !
 ! Grid definition
 integer,protected    :: ngrid
 integer,protected    :: nradial
 integer,protected    :: nangular_fine
 integer,protected    :: nangular_coarse

 real(dp),protected,allocatable :: rr_grid(:,:)
 real(dp),protected,allocatable :: w_grid(:)

 real(dp),parameter,private :: pruning_limit = 0.75_dp    ! in terms of covalent radius

 !
 ! Function evaluation storage
 integer,parameter    :: ngrid_max_stored=20000
 integer              :: ngrid_stored=0
 real(dp),allocatable :: bfr(:,:)
 real(dp),allocatable :: bfgr(:,:,:)
 real(dp),allocatable :: bflr(:,:,:)


contains


!=========================================================================
subroutine setup_dft_grid()
 use m_elements
 use m_atoms
 use m_tools,only: coeffs_gausslegint
 implicit none

 integer              :: iradial,iatom,iangular,ir,igrid
 integer              :: n1,n2,nangular
 real(dp)             :: weight,radius
 real(dp),allocatable :: x1(:),x2(:)
 real(dp),allocatable :: y1(:),y2(:)
 real(dp),allocatable :: z1(:),z2(:)
 real(dp),allocatable :: w1(:),w2(:)
 real(dp),allocatable :: xa(:),wxa(:)
 real(dp)             :: p_becke(natom),s_becke(natom,natom),fact_becke 
 real(dp)             :: mu
 integer              :: jatom,katom
!=====

 select case(grid_level)
 case(10)       ! accuracy not guaranted, just for quick test runs
   nradial         =  25
   nangular_fine   =  26
   nangular_coarse =   6
 case(20)    ! 10 meV accuracy on potentials
   nradial         =  40
   nangular_fine   =  50
   nangular_coarse =  14
 case(30)      !  1 meV accuracy on potentials
   nradial         =  60
   nangular_fine   = 110
   nangular_coarse =  38
 case(40) ! almost perfect potentials
   nradial         =  70
   nangular_fine   = 230
   nangular_coarse =  50
 case(50)    ! overdoing a lot
   nradial         = 200
   nangular_fine   = 434
   nangular_coarse = 434
 case default
   stop'integration quality not recognized'
 end select

 allocate(xa(nradial),wxa(nradial))
 allocate(x1(nangular_fine),y1(nangular_fine),z1(nangular_fine),w1(nangular_fine))
 allocate(x2(nangular_coarse),y2(nangular_coarse),z2(nangular_coarse),w2(nangular_coarse))

 !
 ! spherical integration
 ! radial part with Gauss-Legendre
 call coeffs_gausslegint(-1.0_dp,1.0_dp,xa,wxa,nradial)
 !
 ! Transformation from [-1;1] to [0;+\infty[
 ! taken from M. Krack JCP 1998
 wxa(:) = wxa(:) * ( 1.0_dp / log(2.0_dp) / ( 1.0_dp - xa(:) ) )
 xa(:)  = log( 2.0_dp / (1.0_dp - xa(:) ) ) / log(2.0_dp)

 n1 = nangular_fine
 ! angular part with Lebedev - Laikov
 ! (x1,y1,z1) on the unit sphere
 select case(nangular_fine)
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
 case(110)
   call ld0110(x1,y1,z1,w1,n1)
 case(146)
   call ld0146(x1,y1,z1,w1,n1)
 case(170)
   call ld0170(x1,y1,z1,w1,n1)
 case(230)
   call ld0230(x1,y1,z1,w1,n1)
 case(434)
   call ld0434(x1,y1,z1,w1,n1)
 case default
   WRITE_MASTER(*,*) 'grid points: ',nangular_fine
   stop'Lebedev grid is not available'
 end select

 n2 = nangular_coarse
 ! angular part with Lebedev - Laikov
 ! (x2,y2,z2) on the unit sphere
 select case(nangular_coarse)
 case(6)
   call ld0006(x2,y2,z2,w2,n2)
 case(14)
   call ld0014(x2,y2,z2,w2,n2)
 case(26)
   call ld0026(x2,y2,z2,w2,n2)
 case(38)
   call ld0038(x2,y2,z2,w2,n2)
 case(50)
   call ld0050(x2,y2,z2,w2,n2)
 case(74)
   call ld0074(x2,y2,z2,w2,n2)
 case(86)
   call ld0086(x2,y2,z2,w2,n2)
 case(110)
   call ld0110(x2,y2,z2,w2,n2)
 case(146)
   call ld0146(x2,y2,z2,w2,n2)
 case(170)
   call ld0170(x2,y2,z2,w2,n2)
 case(230)
   call ld0230(x2,y2,z2,w2,n2)
 case(434)
   call ld0434(x2,y2,z2,w2,n2)
 case default
   WRITE_MASTER(*,*) 'grid points: ',nangular_coarse
   stop'Lebedev grid is not available'
 end select

 ! Calculate the total number of grid points
 ngrid = 0
 do iatom=1,natom
   radius = element_covalent_radius(zatom(iatom))
   do iradial=1,nradial
     if( xa(iradial) < pruning_limit * radius ) then
       ngrid = ngrid + nangular_coarse
     else
       ngrid = ngrid + nangular_fine
     endif
   enddo
 enddo

 call init_grid_distribution(ngrid)

 WRITE_MASTER(*,'(/,a)')            ' Setup the DFT quadrature'
 WRITE_MASTER(*,'(a,i4,x,i4,x,i4)') ' discretization grid per atom [radial , angular max - angular min] ',nradial,nangular_fine,nangular_coarse
 WRITE_MASTER(*,'(a,i8)')           ' total number of real-space points for this processor',ngrid



 
 allocate(rr_grid(3,ngrid),w_grid(ngrid))

 igrid = 0
 ir    = 0
 do iatom=1,natom
   radius = element_covalent_radius(zatom(iatom))

   do iradial=1,nradial
     if( xa(iradial) < pruning_limit * radius ) then
       nangular = nangular_coarse
     else
       nangular = nangular_fine
     endif

     do iangular=1,nangular
       igrid = igrid + 1
       if( .NOT. is_my_grid_task(igrid) ) cycle
       ir = ir + 1

       if( xa(iradial) < pruning_limit * radius ) then
         rr_grid(1,ir) = xa(iradial) * x2(iangular) + x(1,iatom)
         rr_grid(2,ir) = xa(iradial) * y2(iangular) + x(2,iatom)
         rr_grid(3,ir) = xa(iradial) * z2(iangular) + x(3,iatom)
         weight   = wxa(iradial) * w2(iangular) * xa(iradial)**2 * 4.0_dp * pi
       else
         rr_grid(1,ir) = xa(iradial) * x1(iangular) + x(1,iatom)
         rr_grid(2,ir) = xa(iradial) * y1(iangular) + x(2,iatom)
         rr_grid(3,ir) = xa(iradial) * z1(iangular) + x(3,iatom)
         weight   = wxa(iradial) * w1(iangular) * xa(iradial)**2 * 4.0_dp * pi
       endif


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

 deallocate(x1,y1,z1,w1)
 deallocate(x2,y2,z2,w2)

end subroutine setup_dft_grid


!=========================================================================
subroutine destroy_dft_grid()
 implicit none

 deallocate(rr_grid)
 deallocate(w_grid)

 if( allocated(bfr) ) then
   call clean_deallocate('basis ftns on grid',bfr)
 endif
 if( allocated(bfgr) ) then
   call clean_deallocate('basis grad ftns on grid',bfgr)
 endif
 if( allocated(bflr) ) then
   call clean_deallocate('basis lapl ftns on grid',bflr)
 endif
 call destroy_grid_distribution()

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
 call clean_allocate('basis ftns on grid',bfr,basis%nbf,ngrid_stored)


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
 real(dp)                   :: basis_function_gradr(3,basis%nbf)
!=====

 WRITE_MASTER(*,*) 'Precalculate the gradients on N grid points',ngrid_stored
 call clean_allocate('basis grad ftns on grid',bfgr,3,basis%nbf,ngrid_stored)


 do igrid=1,ngrid_stored
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_gradr(basis,rr,basis_function_gradr)
   bfgr(:,:,igrid) = basis_function_gradr(:,:)
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
 call clean_allocate('basis grad ftns on grid',bfgr,3,basis%nbf,ngrid_stored)
 call clean_allocate('basis lapl ftns on grid',bflr,3,basis%nbf,ngrid_stored)


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
subroutine get_basis_functions_gradr(basis,igrid,basis_function_gradr)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: igrid
 real(dp),intent(out)       :: basis_function_gradr(3,basis%nbf)
!=====
 real(dp)                   :: rr(3)
!=====

 if( igrid <= ngrid_stored ) then
   basis_function_gradr(:,:) = bfgr(:,:,igrid) 
 else
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_gradr(basis,rr,basis_function_gradr)
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
   ni_cart = number_basis_function_am('CART',li)
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
subroutine calculate_basis_functions_gradr(basis,rr,basis_function_gradr)
 use m_basis_set
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: rr(3)
 real(dp),intent(out)       :: basis_function_gradr(3,basis%nbf)
!=====
 integer              :: ibf,ibf_cart,i_cart
 integer              :: ni,ni_cart,li
 real(dp),allocatable :: basis_function_gradr_cart(:,:)
!=====


 ibf_cart = 1
 ibf      = 1
 do while(ibf_cart<=basis%nbf_cart)
   li      = basis%bf(ibf_cart)%am
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   allocate(basis_function_gradr_cart(3,ni_cart))

   do i_cart=1,ni_cart
     basis_function_gradr_cart(:,i_cart) = eval_basis_function_grad(basis%bf(ibf_cart+i_cart-1),rr)
   enddo

   basis_function_gradr(:,ibf:ibf+ni-1) = MATMUL( basis_function_gradr_cart(:,:) , cart_to_pure(li)%matrix(:,:) )

   deallocate(basis_function_gradr_cart)

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
   ni_cart = number_basis_function_am('CART',li)
   ni      = number_basis_function_am(basis%gaussian_type,li)

   allocate(basis_function_gradr_cart(3,ni_cart))
   allocate(basis_function_laplr_cart(3,ni_cart))

   do i_cart=1,ni_cart

     basis_function_gradr_cart(:,i_cart)        = eval_basis_function_grad(basis%bf(ibf_cart+i_cart-1),rr)
     basis_function_laplr_cart(:,i_cart)        = eval_basis_function_lapl(basis%bf(ibf_cart+i_cart-1),rr)

   enddo

   basis_function_gradr(:,ibf:ibf+ni-1) = MATMUL(  basis_function_gradr_cart(:,:) , cart_to_pure(li)%matrix(:,:) )
   basis_function_laplr(:,ibf:ibf+ni-1) = MATMUL(  basis_function_laplr_cart(:,:) , cart_to_pure(li)%matrix(:,:) )
   deallocate(basis_function_gradr_cart,basis_function_laplr_cart)

   ibf      = ibf      + ni
   ibf_cart = ibf_cart + ni_cart
 enddo


end subroutine calculate_basis_functions_laplr


!=========================================================================
end module m_dft_grid
!=========================================================================
