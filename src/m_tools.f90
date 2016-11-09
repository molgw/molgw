!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! some unsorted basic numerical algorithms
!
!=========================================================================
module m_tools
 use m_definitions
 use m_warning,only: die

 integer,save :: idum

 interface invert
   module procedure invert_dp
   module procedure invert_inplace_dp
   module procedure invert_cdp
 end interface

 interface diagonalize_wo_vectors
   module procedure diagonalize_wo_vectors_dp
 end interface

 interface diagonalize
   module procedure diagonalize_cdp
   module procedure diagonalize_dp
   module procedure diagonalize_sp
   module procedure diagonalize_inplace_dp
   module procedure diagonalize_inplace_sp
 end interface

 interface append_to_list
   module procedure append_to_list_i
   module procedure append_to_list_r
 end interface


contains


!=========================================================================
function matrix_trace(matrix)
 real(dp),intent(in) :: matrix(:,:)
 real(dp)            :: matrix_trace
!=====
 integer :: n1,i1
!=====

 n1 = SIZE( matrix , DIM=1 )
 if( n1 /= SIZE( matrix , DIM=2 ) ) call die('matrix_trace: non square matrix')

 matrix_trace = 0.0_dp
 do i1=1,n1
   matrix_trace = matrix_trace + matrix(i1,i1)
 enddo

end function matrix_trace

!=========================================================================
subroutine init_seed(iseed)
 implicit none
 integer,intent(in),optional :: iseed

 integer :: i,j

 if(PRESENT(iseed)) then
  idum=iseed
 else
  call system_clock(idum,i,j)
 endif

 write(stdout,'(a,1x,i12)') 'Random seed set to',idum


end subroutine


!=========================================================================
function ran1()
 implicit none
 integer IA,IM,IQ,IR,NTAB,NDIV
 real(dp) ran1,AM,EPS,RNMX
 parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
 parameter (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
 integer j,k,iv(NTAB),iy
 save iv,iy
 data iv /NTAB*0/, iy /0/

 if (idum.le.0.or.iy.eq.0) then
   idum=max(-idum,1)
   do j=NTAB+8,1,-1
     k=idum/IQ
     idum=IA*(idum-k*IQ)-IR*k
     if (idum.lt.0) idum=idum+IM
     if (j.le.NTAB) iv(j)=idum
   enddo
   iy=iv(1)
 endif
 k=idum/IQ
 idum=IA*(idum-k*IQ)-IR*k
 if (idum.lt.0) idum=idum+IM
 j=1+iy/NDIV
 iy=iv(j)
 iv(j)=idum
 ran1=min(AM*iy,RNMX)

end function ran1


!=========================================================================
function random()
 implicit none
 real(dp) :: random

 call random_number(random)

end function random


!=========================================================================
function matrix_is_symmetric(n,matrix)
 implicit none
 logical             :: matrix_is_symmetric
 integer,intent(in)  :: n
 real(dp),intent(in) :: matrix(n,n)
!===== 
 integer :: ii,jj
!===== 

 matrix_is_symmetric = .TRUE.
 do ii=1,n
   do jj=1,ii-1
     if( ABS( matrix(ii,jj) - matrix(jj,ii) ) > 1.0e-5_dp ) then
       matrix_is_symmetric=.FALSE.
       return
     endif
   enddo
 enddo

end function matrix_is_symmetric


!=========================================================================
subroutine invert_dp(n,matrix,matrix_inv)
 implicit none

 integer,intent(in)   :: n
 real(dp),intent(in)  :: matrix(n,n)
 real(dp),intent(out) :: matrix_inv(n,n)
!=====
 real(dp) :: work(n)
 integer  :: ipiv(n),info
!=====

 matrix_inv(:,:) = matrix(:,:)

 call DGETRF(n,n,matrix_inv,n,ipiv,info)
 if(info/=0) call die('FAILURE in DGETRF')

 call DGETRI(n,matrix_inv,n,ipiv,work,n,info)
 if(info/=0) call die('FAILURE in DGETRI')


end subroutine invert_dp


!=========================================================================
subroutine invert_inplace_dp(n,matrix)
 implicit none

 integer,intent(in)     :: n
 real(dp),intent(inout) :: matrix(n,n)
!=====
 real(dp) :: work(n)
 integer  :: ipiv(n),info
!=====

 call DGETRF(n,n,matrix,n,ipiv,info)
 if(info/=0) call die('FAILURE in DGETRF')

 call DGETRI(n,matrix,n,ipiv,work,n,info)
 if(info/=0) call die('FAILURE in DGETRI')


end subroutine invert_inplace_dp


!=========================================================================
subroutine invert_cdp(n,matrix,matrix_inv)
 implicit none

 integer,intent(in) :: n
 complex(dp),intent(in) :: matrix(n,n)
 complex(dp),intent(out) :: matrix_inv(n,n)
!=====
 complex(dp) :: work(n)
 integer     :: ipiv(n),info
!=====

 matrix_inv(:,:) = matrix(:,:)

 call ZGETRF(n,n,matrix_inv,n,ipiv,info)
 if(info/=0) call die('FAILURE in ZGETRF')

 call ZGETRI(n,matrix_inv,n,ipiv,work,n,info)
 if(info/=0) call die('FAILURE in ZGETRI')


end subroutine invert_cdp


!=========================================================================
subroutine diagonalize_wo_vectors_dp(n,matrix,eigval)
 implicit none
 integer,intent(in) :: n
 real(dp),intent(inout) :: matrix(n,n)
 real(dp),intent(out) :: eigval(n)

 real(dp) :: z(1,n)
 integer  :: iwork(5*n),ifail(n)
 integer :: info
 real(dp) :: work(3*n-1)
! real(dp) :: work(8*n)
 
 call DSYEV('N','U',n,matrix,n,eigval,work,3*n-1,info)

! call DSYEVX('N','A','U',n,matrix,n,0.0_dp,0.0_dp,0,0,&
!                         1.0e-20_dp,n,eigval,z,1,work,8*n,iwork,&
!                         ifail,info)

end subroutine diagonalize_wo_vectors_dp


!=========================================================================
subroutine diagonalize_cdp(n,matrix,eigval,eigvec)
 implicit none
 integer,intent(in) :: n
 complex(dp),intent(in) :: matrix(n,n)
 real(dp),intent(out) :: eigval(n)
 complex(dp),intent(out) :: eigvec(n,n)

 complex(dp) :: work(2*n-1)
 real(dp) :: rwork(3*n-2)
 integer :: info

 eigvec(:,:) = matrix(:,:)

 call ZHEEV('V','U',n,eigvec,n,eigval,work,2*n-1,rwork,info)

end subroutine diagonalize_cdp


!=========================================================================
subroutine diagonalize_dp(n,matrix,eigval,eigvec)
 implicit none
 integer,intent(in) :: n
 real(dp),intent(in) :: matrix(n,n)
 real(dp),intent(out) :: eigval(n)
 real(dp),intent(out) :: eigvec(n,n)

 real(dp) :: work(3*n-1)
 integer :: info
 
 eigvec(:,:) = matrix(:,:)

 call DSYEV('V','U',n,eigvec,n,eigval,work,3*n-1,info)

end subroutine diagonalize_dp


!=========================================================================
subroutine diagonalize_sp(n,matrix,eigval,eigvec)
 implicit none
 integer,intent(in) :: n
 real(sp),intent(in) :: matrix(n,n)
 real(sp),intent(out) :: eigval(n)
 real(sp),intent(out) :: eigvec(n,n)

 real(sp) :: work(3*n-1)
 integer :: info

 eigvec(:,:) = matrix(:,:)

 call SSYEV('V','U',n,eigvec,n,eigval,work,3*n-1,info)

end subroutine diagonalize_sp


!=========================================================================
subroutine diagonalize_inplace_dp(n,matrix,eigval)
 implicit none
 integer,intent(in) :: n
 real(dp),intent(inout) :: matrix(n,n)
 real(dp),intent(out) :: eigval(n)
!=====
 real(dp),allocatable :: work(:)
 integer              :: lwork,info
!=====

 lwork = -1
 allocate(work(1))
 call DSYEV('V','U',n,matrix,n,eigval,work,lwork,info)
 lwork = NINT(work(1))
 deallocate(work)

 allocate(work(lwork))
 call DSYEV('V','U',n,matrix,n,eigval,work,lwork,info)
 deallocate(work)

end subroutine diagonalize_inplace_dp


!=========================================================================
subroutine diagonalize_inplace_sp(n,matrix,eigval)
 implicit none
 integer,intent(in) :: n
 real(sp),intent(inout) :: matrix(n,n)
 real(sp),intent(out) :: eigval(n)
!=====
 real(dp),allocatable :: work(:)
 integer              :: lwork,info
!=====

 lwork = -1
 allocate(work(1))
 call SSYEV('V','U',n,matrix,n,eigval,work,lwork,info)
 lwork = NINT(work(1))
 deallocate(work)

 allocate(work(lwork))
 call SSYEV('V','U',n,matrix,n,eigval,work,lwork,info)
 deallocate(work)

end subroutine diagonalize_inplace_sp


!=========================================================================
!
! Generalized eigenvalue problem
!
subroutine diagonalize_generalized_sym(n,matrix,overlap,eigval,eigvec)
 implicit none
 integer,intent(in) :: n
 real(dp),intent(in) :: matrix(n,n),overlap(n,n)
 real(dp),intent(out) :: eigval(n)
 real(dp),intent(out) :: eigvec(n,n)
!=====
 real(dp) :: work(3*n-1),tmp(n,n)
 integer :: info
!=====

 ! write(stdout,*) 'diagonalize_generalized_sym: Enter'
 eigvec(:,:) = matrix(:,:)
 tmp(:,:) = overlap(:,:)

 ! A*x = lambda * B * x
 call DSYGV(1,'V','U',n,eigvec,n,tmp,n,eigval,work,3*n-1,info)
 if(info/=0) call die('ERROR in the symmetric generalized eigenvalue problem')
! write(stdout,*) 'optimal lwork',REAL(work(1))


 ! write(stdout,*) 'diagonalize_generalized_sym: Exit'
end subroutine diagonalize_generalized_sym


!=========================================================================
subroutine coeffs_gausslegint(xmin,xmax,x,weights,n)
!
! Compute the coefficients (supports and weights)
! for Gauss-Legendre integration.
! Inspired by a routine due to G. Rybicki.
!
! Input:
! xmin=lower bound of integration
! xmax=upper bound of integration
! n=order of integration
!
! Output:
! x(n)=array of support points
! weights(n)=array of integration weights
!

 implicit none
     
 integer,intent(in) :: n 
 real(dp),intent(in) :: xmin,xmax
 real(dp),intent(out) :: x(n),weights(n)
!=====
 real(dp) :: tol,xl,xmean,z,p1,p2,p3,pp,z1
 integer  :: i,j
!=====

 tol=1.d-13

 xl=(xmax-xmin)*0.50_dp
 xmean=(xmax+xmin)*0.50_dp

 do i=1,(n+1)/2
  z = COS(pi*(i-0.250_dp)/(n+0.50_dp))
 
  do 

    p1=1.0_dp
    p2=0.0_dp
 
    do j=1,n
     
     p3=p2
     p2=p1
     p1=((2.0_dp*j - 1.0_dp)*z*p2 - (j-1.0_dp)*p3)/j
   
    enddo
  
    pp=n*(p2-z*p1)/(1.0_dp - z**2)
    z1=z
    z=z1-p1/pp
    
    if(abs(z-z1) < tol) exit

  enddo

  x(i)=xmean-xl*z
  x(n+1-i)=xmean+xl*z
  weights(i)=2.0_dp * xl/((1.0_dp-z**2)*pp**2)
  weights(n+1-i)=weights(i)

 enddo

end subroutine coeffs_gausslegint


!=========================================================================
subroutine gequad(p,w)
 implicit none

!Arguments
 real(dp), intent(out) :: p(89),w(89)
!
!       range [10^(-9),1] and accuracy ~10^(-8);
!
  p(1)=4.96142640560223544d19
  p(2)=1.37454269147978052d19
  p(3)=7.58610013441204679d18
  p(4)=4.42040691347806996d18
  p(5)=2.61986077948367892d18
  p(6)=1.56320138155496681d18
  p(7)=9.35645215863028402d17
  p(8)=5.60962910452691703d17
  p(9)=3.3666225119686761d17
  p(10)=2.0218253197947866d17
  p(11)=1.21477756091902017d17
  p(12)=7.3012982513608503d16
  p(13)=4.38951893556421099d16
  p(14)=2.63949482512262325d16
  p(15)=1.58742054072786174d16
  p(16)=9.54806587737665531d15
  p(17)=5.74353712364571709d15
  p(18)=3.455214877389445d15
  p(19)=2.07871658520326804d15
  p(20)=1.25064667315629928d15
  p(21)=7.52469429541933745d14
  p(22)=4.5274603337253175d14
  p(23)=2.72414006900059548d14
  p(24)=1.63912168349216752d14
  p(25)=9.86275802590865738d13
  p(26)=5.93457701624974985d13
  p(27)=3.5709554322296296d13
  p(28)=2.14872890367310454d13
  p(29)=1.29294719957726902d13
  p(30)=7.78003375426361016d12
  p(31)=4.68148199759876704d12
  p(32)=2.8169955024829868d12
  p(33)=1.69507790481958464d12
  p(34)=1.01998486064607581d12
  p(35)=6.13759486539856459d11
  p(36)=3.69320183828682544d11
  p(37)=2.22232783898905102d11
  p(38)=1.33725247623668682d11
  p(39)=8.0467192739036288d10
  p(40)=4.84199582415144143d10
  p(41)=2.91360091170559564d10
  p(42)=1.75321747475309216d10
  p(43)=1.0549735552210995d10
  p(44)=6.34815321079006586d9
  p(45)=3.81991113733594231d9
  p(46)=2.29857747533101109d9
  p(47)=1.38313653595483694d9
  p(48)=8.32282908580025358d8
  p(49)=5.00814519374587467d8
  p(50)=3.01358090773319025d8
  p(51)=1.81337994217503535d8
  p(52)=1.09117589961086823d8
  p(53)=6.56599771718640323d7
  p(54)=3.95099693638497164d7
  p(55)=2.37745694710665991d7
  p(56)=1.43060135285912813d7
  p(57)=8.60844290313506695d6
  p(58)=5.18000974075383424d6
  p(59)=3.116998193057466d6
  p(60)=1.87560993870024029d6
  p(61)=1.12862197183979562d6
  p(62)=679132.441326077231_dp
  p(63)=408658.421279877969_dp
  p(64)=245904.473450669789_dp
  p(65)=147969.568088321005_dp
  p(66)=89038.612357311147_dp
  p(67)=53577.7362552358895_dp
  p(68)=32239.6513926914668_dp
  p(69)=19399.7580852362791_dp
  p(70)=11673.5323603058634_dp
  p(71)=7024.38438577707758_dp
  p(72)=4226.82479307685999_dp
  p(73)=2543.43254175354295_dp
  p(74)=1530.47486269122675_dp
  p(75)=920.941785160749482_dp
  p(76)=554.163803906291646_dp
  p(77)=333.46029740785694_dp
  p(78)=200.6550575335041_dp
  p(79)=120.741366914147284_dp
  p(80)=72.6544243200329916_dp
  p(81)=43.7187810415471025_dp
  p(82)=26.3071631447061043_dp
  p(83)=15.8299486353816329_dp
  p(84)=9.52493152341244004_dp
  p(85)=5.72200417067776041_dp
  p(86)=3.36242234070940928_dp
  p(87)=1.75371394604499472_dp
  p(88)=0.64705932650658966_dp
  p(89)=0.072765905943708247_dp
!
  w(1)=47.67445484528304247d10
  w(2)=11.37485774750442175d9
  w(3)=78.64340976880190239d8
  w(4)=46.27335788759590498d8
  w(5)=24.7380464827152951d8
  w(6)=13.62904116438987719d8
  w(7)=92.79560029045882433d8
  w(8)=52.15931216254660251d8
  w(9)=31.67018011061666244d8
  w(10)=1.29291036801493046d8
  w(11)=1.00139319988015862d8
  w(12)=7.75892350510188341d7
  w(13)=6.01333567950731271d7
  w(14)=4.66141178654796875d7
  w(15)=3.61398903394911448d7
  w(16)=2.80225846672956389d7
  w(17)=2.1730509180930247d7
  w(18)=1.68524482625876965d7
  w(19)=1.30701489345870338d7
  w(20)=1.01371784832269282d7
  w(21)=7.86264116300379329d6
  w(22)=6.09861667912273717d6
  w(23)=4.73045784039455683d6
  w(24)=3.66928949951594161d6
  w(25)=2.8462050836230259d6
  w(26)=2.20777394798527011d6
  w(27)=1.71256191589205524d6
  w(28)=1.32843556197737076d6
  w(29)=1.0304731275955989d6
  w(30)=799345.206572271448_dp
  w(31)=620059.354143595343_dp
  w(32)=480986.704107449333_dp
  w(33)=373107.167700228515_dp
  w(34)=289424.08337412132_dp
  w(35)=224510.248231581788_dp
  w(36)=174155.825690028966_dp
  w(37)=135095.256919654065_dp
  w(38)=104795.442776800312_dp
  w(39)=81291.4458222430418_dp
  w(40)=63059.0493649328682_dp
  w(41)=48915.9040455329689_dp
  w(42)=37944.8484018048756_dp
  w(43)=29434.4290473253969_dp
  w(44)=22832.7622054490044_dp
  w(45)=17711.743950151233_dp
  w(46)=13739.287867104177_dp
  w(47)=10657.7895710752585_dp
  w(48)=8267.42141053961834_dp
  w(49)=6413.17397520136448_dp
  w(50)=4974.80402838654277_dp
  w(51)=3859.03698188553047_dp
  w(52)=2993.51824493299154_dp
  w(53)=2322.1211966811754_dp
  w(54)=1801.30750964719641_dp
  w(55)=1397.30379659817038_dp
  w(56)=1083.91149143250697_dp
  w(57)=840.807939169209188_dp
  w(58)=652.228524366749422_dp
  w(59)=505.944376983506128_dp
  w(60)=392.469362317941064_dp
  w(61)=304.444930257324312_dp
  w(62)=236.162932842453601_dp
  w(63)=183.195466078603525_dp
  w(64)=142.107732186551471_dp
  w(65)=110.23530215723992_dp
  w(66)=85.5113346705382257_dp
  w(67)=66.3325469806696621_dp
  w(68)=51.4552463353841373_dp
  w(69)=39.9146798429449273_dp
  w(70)=30.9624728409162095_dp
  w(71)=24.018098812215013_dp
  w(72)=18.6312338024296588_dp
  w(73)=14.4525541233150501_dp
  w(74)=11.2110836519105938_dp
  w(75)=8.69662175848497178_dp
  w(76)=6.74611236165731961_dp
  w(77)=5.23307018057529994_dp
  w(78)=4.05937850501539556_dp
  w(79)=3.14892659076635714_dp
  w(80)=2.44267408211071604_dp
  w(81)=1.89482240522855261_dp
  w(82)=1.46984505907050079_dp
  w(83)=1.14019261330527007_dp
  w(84)=0.884791217422925293_dp
  w(85)=0.692686387080616483_dp
  w(86)=0.585244576897023282_dp
  w(87)=0.576182522545327589_dp
  w(88)=0.596688817388997178_dp
  w(89)=0.607879901151108771_dp
!
end subroutine gequad


!=========================================================================
subroutine check_unitarity(n,cmat)
 implicit none
 integer,intent(in) :: n
 complex(dp),intent(in) :: cmat(n,n)
!
 real(dp),parameter :: tol=1.0e-9_dp
 integer :: i,j
 complex(dp) :: cmat_tmp(n,n)
!
  cmat_tmp = MATMUL( cmat , TRANSPOSE(CONJG(cmat)) )
  do i=1,n
    do j=1,n
      if(i==j) then
       if(ABS(cmat_tmp(i,j)-1.0_dp)>tol) then
         write(stdout,*) i,j,cmat_tmp(i,j)
         call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
       endif
      else
       if(ABS(cmat_tmp(i,j))>tol) then
         write(stdout,*) i,j,cmat_tmp(i,j)
         call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
       endif
      endif
    enddo
  enddo
  cmat_tmp = MATMUL( TRANSPOSE(CONJG(cmat)) , cmat )
  do i=1,n
    do j=1,n
      if(i==j) then
       if(ABS(cmat_tmp(i,j)-1.0_dp)>tol) then
         write(stdout,*) i,j,cmat_tmp(i,j)
         call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
       endif
      else
       if(ABS(cmat_tmp(i,j))>tol) then
         write(stdout,*) i,j,cmat_tmp(i,j)
         call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
       endif
      endif
    enddo
  enddo
end subroutine check_unitarity


!==========================================
function gamma_function(rin)
 implicit none
 real(dp),intent(in) :: rin
 real(dp) :: gamma_function
!=====
 integer :: nlocal
!=====

 !
 ! just hard coded for some small half-integers
 if( ABS( rin - NINT(rin) ) - 0.5  > 1.d-6 ) call die('GAMMA FUNCTION NOT CODED')

 nlocal = FLOOR(rin-0.499999)

 select case(nlocal)
 case(0)
  gamma_function = 1.0_dp
 case(1)
  gamma_function = 0.5_dp
 case(2)
  gamma_function = 0.75_dp
 case(3)
  gamma_function = 15.0_dp / 8.0_dp
 case(4)
  gamma_function = 105.0_dp / 16.0_dp
 case(5)
  gamma_function = 945.0_dp / 32.0_dp
 case(6)
  gamma_function = 10395.0_dp / 64.0_dp
 case default
  gamma_function = double_factorial(nlocal+1)/ 2**( 0.5_dp * (nlocal+1)/2)
 end select

 gamma_function = gamma_function * SQRT(pi)

end function gamma_function


!==========================================
function double_factorial(intin)
 implicit none
 integer,intent(in) :: intin
 real(dp) :: double_factorial
!=====
 ! just hard coded for some small integers

 select case (intin)
 case(-1) 
   double_factorial = 1.0_dp
 case( 0) 
   double_factorial = 1.0_dp
 case( 1) 
   double_factorial = 1.0_dp
 case( 2) 
   double_factorial = 2.0_dp
 case( 3) 
   double_factorial = 3.0_dp
 case( 4) 
   double_factorial = 8.0_dp
 case( 5) 
   double_factorial = 15.0_dp
 case( 6) 
   double_factorial = 48.0_dp
 case( 7) 
   double_factorial = 105.0_dp
 case( 8) 
   double_factorial = 384.0_dp
 case( 9) 
   double_factorial = 945.0_dp
 case(10) 
   double_factorial = 3840.0_dp
 case(11) 
   double_factorial = 10395.0_dp
 case(12) 
   double_factorial = 46080.0_dp
 case(13) 
   double_factorial = 135135.0_dp
 case(14) 
   double_factorial = 645120.0_dp
 case(15)
   double_factorial = 2027025.0_dp
 case(16) 
   double_factorial = 10321920.0_dp
 case(17) 
   double_factorial = 34459425.0_dp
 case(18) 
   double_factorial = 185794560.0_dp
 case(19) 
   double_factorial = 654729075.0_dp
 case(20) 
   double_factorial = 3715891200.0_dp
 case(21) 
   double_factorial = 13749310575.0_dp
 case(22)
   double_factorial = 81749606400.0_dp
 case(23) 
   double_factorial = 316234143225.0_dp
 case(25) 
   double_factorial = 7905853580625.0_dp
 case(27) 
   double_factorial = 213458046676875.0_dp
 case(29) 
   double_factorial = 6190283353629375.0_dp
 case(31) 
   double_factorial = 191898783962510625.0_dp
 case default
   write(stdout,*) 'integer =',intin
   write(stdout,*) 'double factorial not coded for this integer value'
   double_factorial = 1
 end select

end function double_factorial


!=========================================================================
subroutine cross_product(u1,u2,u3)
 implicit none
 real(dp),intent(in)  :: u1(3),u2(3)
 real(dp),intent(out) :: u3(3)
!=====

 u3(1) = u1(2) * u2(3) - u1(3) * u2(2)
 u3(2) = u1(3) * u2(1) - u1(1) * u2(3)
 u3(3) = u1(1) * u2(2) - u1(2) * u2(1)

end subroutine cross_product


!=========================================================================
function capitalize(str)
 implicit none
 character(*), intent(in) :: str
 character(LEN(str))      :: capitalize
!=====
 character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
 character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
 integer :: ic, ii
!=====

 capitalize = str
 do ii=1,LEN_TRIM(str)
   ic = INDEX(low,str(ii:ii))
   if (ic > 0) capitalize(ii:ii) = cap(ic:ic)
 end do

end function capitalize


!=========================================================================
function orbital_momentum_number(amc)
 character(len=1),intent(in) :: amc
 integer :: orbital_momentum_number
!=====

 select case(capitalize(amc))
 case('S')
   orbital_momentum_number = 0
 case('P')
   orbital_momentum_number = 1
 case('D')
   orbital_momentum_number = 2
 case('F')
   orbital_momentum_number = 3
 case('G')
   orbital_momentum_number = 4
 case('H')
   orbital_momentum_number = 5
 case('I')
   orbital_momentum_number = 6
 case('K')
   orbital_momentum_number = 7
 case default
   write(stdout,*) amc,capitalize(amc)
   call die('orbital_momentum_number: keyword unknown')
 end select


end function orbital_momentum_number


!=========================================================================
function orbital_momentum_name(am)
 integer,intent(in) :: am
 character(len=1) :: orbital_momentum_name
!=====

 select case(am)
 case(0)
   orbital_momentum_name='s'
 case(1)
   orbital_momentum_name='p'
 case(2)
   orbital_momentum_name='d'
 case(3)
   orbital_momentum_name='f'
 case(4)
   orbital_momentum_name='g'
 case(5)
   orbital_momentum_name='h'
 case(6)
   orbital_momentum_name='i'
 case(7)
   orbital_momentum_name='k'
 case(8)
   orbital_momentum_name='l'
 case(9)
   orbital_momentum_name='m'
 case(10)
   orbital_momentum_name='n'
 case(11)
   orbital_momentum_name='o'
 case(12)
   orbital_momentum_name='q'
 case default
   orbital_momentum_name='x'
 end select

end function orbital_momentum_name


!=========================================================================
subroutine append_to_list_i(new_element,list)
 implicit none

 integer,intent(in)                :: new_element
 integer,allocatable,intent(inout) :: list(:)
!=====
 integer :: nsize
 integer,allocatable :: list_old(:)
!=====

 if( ALLOCATED(list) ) then
   nsize = SIZE(list)
 else
   nsize = 0
 endif

 ! Copy old list and free the list
 allocate(list_old(nsize))
 if( nsize > 0 ) then
   list_old(1:nsize) =list(1:nsize)
   deallocate(list)
 endif

 allocate(list(nsize+1))
 if( nsize > 0 ) then
   list(1:nsize) = list_old(1:nsize)
 endif
 list(nsize+1) = new_element


end subroutine append_to_list_i


!=========================================================================
subroutine append_to_list_r(new_element,list)
 implicit none

 real(dp),intent(in)                :: new_element
 real(dp),allocatable,intent(inout) :: list(:)
!=====
 integer :: nsize
 real(dp),allocatable :: list_old(:)
!=====

 if( ALLOCATED(list) ) then
   nsize = SIZE(list)
 else
   nsize = 0
 endif

 ! Copy old list and free the list
 allocate(list_old(nsize))
 if( nsize > 0 ) then
   list_old(1:nsize) =list(1:nsize)
   deallocate(list)
 endif

 allocate(list(nsize+1))
 if( nsize > 0 ) then
   list(1:nsize) = list_old(1:nsize)
 endif
 list(nsize+1) = new_element


end subroutine append_to_list_r


!=========================================================================
end module m_tools
