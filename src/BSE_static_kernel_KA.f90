subroutine BSE_static_kernel_KA(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,Omega,rho,WA)

! Compute the BSE static kernel for the resonant block

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: WA(nS,nS)

! Initialize 

  WA(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          chi = 0d0
          do kc=1,nS
            eps = Omega(kc)**2 + eta**2
            chi = chi + rho(i,j,kc)*rho(a,b,kc)*Omega(kc)/eps
          enddo

          WA(ia,jb) = WA(ia,jb) + lambda*ERI(i,b,j,a) - 4d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine BSE_static_kernel_KA
