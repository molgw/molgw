!=========================================================================
module m_atoms
 use m_definitions
 use m_warning,only: die
 use m_elements

 real(dp),parameter,private     :: tol_geom=1.0e-5_dp

 integer,public                 :: natom_basis
 integer,public                 :: natom
 integer,public                 :: nghost
 integer,protected              :: natom_type
 integer,protected              :: nbond

 real(dp),allocatable,protected :: zatom(:)
 integer,allocatable,protected  :: basis_element(:)

 real(dp),allocatable,protected :: x(:,:)

 logical,protected              :: inversion=.TRUE.
 logical,protected              :: linear=.TRUE.
 logical,protected              :: planar=.TRUE.
 real(dp),protected             :: xcenter(3)
 real(dp),protected             :: xnormal(3)


contains


!=========================================================================
subroutine init_atoms(natom_read,nghost_read,zatom_read,x_read)
 use m_tools,only: cross_product
 implicit none
 integer,intent(in)  :: natom_read,nghost_read
 real(dp),intent(in) :: zatom_read(natom_read+nghost_read),x_read(3,natom_read+nghost_read)
!=====
 integer  :: iatom,jatom
 real(dp) :: xtmp(3),x21(3),x31(3)
 logical  :: found
!=====

 natom  = natom_read
 nghost = nghost_read
 natom_basis = natom + nghost

 allocate(zatom(natom_basis))
 allocate(basis_element(natom_basis))
 allocate(x(3,natom_basis))

 zatom(1:natom)              = zatom_read(1:natom)
 ! Ghost atoms do not have a positive nucleus
 zatom(natom+1:natom+nghost) = 0.0_dp
 ! But ghost atoms have basis functions centered on them.
 basis_element(:)=NINT(zatom_read(:))

 x(:,:) = x_read(:,:)

 !
 ! Check for atoms too close
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( x(:,iatom)-x(:,jatom) ) < 0.2 ) then
       write(stdout,*) 'Atoms',iatom,jatom
       write(stdout,*) 'are closer than 0.2 bohr'
       call die('atoms too close')
     endif
   enddo
 enddo

 !
 ! Find the covalent bond is a simple distance criterium
 nbond = 0
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( x(:,iatom)-x(:,jatom) ) < 4.0 ) then
       nbond = nbond + 1
     endif
   enddo
 enddo

 !
 ! Does the molecule have inversion symmetry?
 call find_inversion()

 !
 ! Is the molecule linear, planar?
 if( natom > 2 ) then
   x21(:) = x(:,2) - x(:,1)
   do iatom=3,natom
     x31(:) = x(:,iatom) - x(:,1)
     call cross_product(x21,x31,xnormal)
     if( NORM2(xnormal(:)) > tol_geom ) then
       xnormal(:) = xnormal(:) / NORM2(xnormal(:))
       linear=.FALSE.
       exit
     endif
   enddo
   if( .NOT. linear) then
     do iatom=1,natom
       if( ABS(DOT_PRODUCT( x(:,iatom) , xnormal(:) )) > tol_geom ) planar=.FALSE.
     enddo
   else
     planar=.FALSE.
   endif
 else
  ! Molecule is linear
  ! Set planar to FALSE for safety
  linear=.TRUE.
  planar=.FALSE.
 endif





end subroutine init_atoms


!=========================================================================
function atoms_core_states()
 implicit none
 integer :: atoms_core_states
!=====
 integer :: iatom
!=====

 atoms_core_states=0
 do iatom=1,natom
   atoms_core_states = atoms_core_states + element_core(zatom(iatom))
 enddo
 
end function atoms_core_states


!=========================================================================
subroutine get_bondcenter(ibond,xbond)
 implicit none
 integer,intent(in)   :: ibond
 real(dp),intent(out) :: xbond(3)
!=====
 integer :: iatom,jatom,jbond
!=====

 jbond = 0
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( x(:,iatom)-x(:,jatom) ) < 4.0 ) then
       jbond = jbond + 1
       if(jbond==ibond) then
         xbond(:) = 0.5_dp * ( x(:,iatom) - x(:,jatom) )
         return
       endif
     endif
   enddo
 enddo


end subroutine get_bondcenter


!=========================================================================
subroutine destroy_atoms()
 implicit none
!=====

 if(ALLOCATED(zatom)) deallocate(zatom)
 if(ALLOCATED(basis_element)) deallocate(basis_element)
 if(ALLOCATED(x)) deallocate(x)

end subroutine destroy_atoms


!=========================================================================
subroutine nucleus_nucleus_energy(energy)
 implicit none
 real(dp),intent(out) :: energy
!=====
 integer              :: iatom,jatom
!=====

 energy = 0.0_dp
 do iatom=1,natom
   do jatom=iatom+1,natom
     energy = energy + zatom(iatom) * zatom(jatom) / SQRT( SUM( (x(:,iatom) - x(:,jatom))**2) )
   enddo
 enddo

end subroutine nucleus_nucleus_energy


!=========================================================================
subroutine find_inversion()
 implicit none
!=====
 integer  :: iatom,jatom
 logical  :: found
 real(dp) :: xtmp(3)
!=====

 xcenter(:) = 0.0_dp
 do iatom=1,natom
   xcenter(:) = xcenter(:) + x(:,iatom) / REAL(natom,dp)
 enddo

 do iatom=1,natom
   xtmp(:) = 2.0_dp * xcenter(:) - x(:,iatom)
   found = .FALSE.
   do jatom=1,natom
     if( NORM2( xtmp(:) - x(:,jatom) ) < tol_geom ) then
       if( ABS(zatom(iatom)-zatom(jatom)) < tol_geom ) found = .TRUE.
       exit
     endif
   enddo
   inversion = inversion .AND. found
 enddo

end subroutine find_inversion


!=========================================================================
subroutine find_rotations()
 implicit none
!=====
 integer  :: iatom,jatom,katom
 logical  :: found
 real(dp) :: xtmpi(3),xtmpj(3),xtmpk(3),rot_axis(3)
 real(dp) :: angle
 integer  :: order
!=====

 ! testing axis from center to 1 atom
 do iatom=1,natom
   rot_axis(:) = x(:,iatom) - xcenter(:)
   if( NORM2(rot_axis) < 1.0e-5_dp ) cycle
   rot_axis(:) = rot_axis(:) / NORM2(rot_axis)

   do jatom=1,natom
     do katom=1,natom
       if( .NOT. same_element(jatom,katom) ) cycle
       xtmpj(:) = x(:,jatom) - xcenter(:)
       xtmpk(:) = x(:,katom) - xcenter(:)
       xtmpj(:) = xtmpj(:) / NORM2(xtmpj)
       xtmpk(:) = xtmpk(:) / NORM2(xtmpk)
       angle = ACOS( DOT_PRODUCT( xtmpj , xtmpk ) )
       order = get_order(angle)
     enddo
   enddo

 enddo




contains

function get_order(angle)
 implicit none
 real(dp),intent(in) :: angle
 integer :: get_order
!=====
 integer :: itested
!=====

 get_order = 0
 do itested=2,12
   if( ABS( angle * itested - 2.0_dp * pi ) < 1.0e-5_dp ) then 
     get_order = itested
     return
   endif
 enddo

end function get_order

end subroutine find_rotations


!=========================================================================
function valid_rotation(order,rot_axis)
 implicit none
 integer,intent(in)  :: order
 real(dp),intent(in) :: rot_axis(3)
 logical             :: valid_rotation
!=====
 real(dp) :: xin(3),xout(3)
 real(dp) :: costheta,sintheta,dr
 integer  :: iatom
!=====

 if( order == 0 ) then
   valid_rotation = .FALSE.
   return
 endif

 do iatom=1,natom
   xin(:) = x(:,iatom) - xcenter(:)
   dr = NORM2(xin)
   if( dr < 1.0e-5_dp ) then 
     valid_rotation=.TRUE.
     return
   endif
   xin(:) = xin(:) / dr
   costheta = DOT_PRODUCT( rot_axis , xin)
   if( ABS(ABS(costheta)-1.0) < 1.0e-5_dp ) then
     valid_rotation=.TRUE.
     return
   endif
   sintheta = SQRT( 1.0 - costheta**2)

 enddo


end function valid_rotation


!=========================================================================
function same_element(iatom,jatom)
 implicit none
 integer,intent(in) :: iatom,jatom
 logical :: same_element
!=====

 same_element = ( ABS( zatom(iatom) - zatom(jatom) ) < 1.0e-5 )

end function same_element



end module m_atoms
!=========================================================================
