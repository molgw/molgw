subroutine linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,B_lr)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  integer,intent(in)            :: ispin,nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  
! Local variables

  double precision              :: delta_dRPA

  integer                       :: i,j,a,b,ia,jb

! Output variables

  double precision,intent(out)  :: B_lr(nS,nS)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build B matrix for singlet manifold

  if(ispin == 1) then

    ia = 0
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        ia = ia + 1
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
 
            B_lr(ia,jb) = 2d0*lambda*ERI(i,j,a,b) - (1d0 - delta_dRPA)*lambda*ERI(i,j,b,a)
 
          end do
        end do
      end do
    end do

  end if

! Build B matrix for triplet manifold

  if(ispin == 2) then

    ia = 0
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        ia = ia + 1
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
 
            B_lr(ia,jb) = - (1d0 - delta_dRPA)*lambda*ERI(i,j,b,a)
 
          end do
        end do
      end do
    end do

  end if

! Build B matrix for spin orbitals 

  if(ispin == 3) then

    ia = 0
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        ia = ia + 1
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
 
            B_lr(ia,jb) = lambda*ERI(i,j,a,b) - (1d0 - delta_dRPA)*lambda*ERI(i,j,b,a)
 
          end do
        end do
      end do
    end do

  end if

end subroutine linear_response_B_matrix
