
module m_bitwise
 use, intrinsic :: ISO_C_BINDING,only: C_SIZEOF

 integer,parameter :: nbit = 8 * C_SIZEOF(1)


contains


!==================================================================
subroutine bitwise_add_onebit(i,key3)
 implicit none
 integer,intent(in)    :: i
 integer,intent(inout) :: key3(:)
!=====
 integer :: j,n
!=====

 n = SIZE(key3)

 j = i / nbit
 key3(j+1) = key3(j+1) + SHIFTL(1,i-nbit*j)

end subroutine bitwise_add_onebit


!==================================================================
subroutine bitwise_sub_onebit(i,key3)
 implicit none
 integer,intent(in)    :: i
 integer,intent(inout) :: key3(:)
!=====
 integer :: j,n
!=====

 n = SIZE(key3)

 j = i / nbit
 key3(j+1) = key3(j+1) - SHIFTL(1,i-nbit*j)

end subroutine bitwise_sub_onebit


!==================================================================
subroutine bitwise_and(key1,key2,key3)
 implicit none

 integer,intent(in)  :: key1(:),key2(:)
 integer,intent(out) :: key3(:)
!=====
 integer :: j,n
!=====

 n = SIZE(key1)

 forall(j=1:n)
   key3(j) = IAND(key1(j),key2(j))
 endforall


end subroutine bitwise_and


!==================================================================
subroutine bitwise_xor(key1,key2,key3)
 implicit none

 integer,intent(in)  :: key1(:),key2(:)
 integer,intent(out) :: key3(:)
!=====
 integer :: j,n
!=====

 n = SIZE(key1)

 forall(j=1:n)
   key3(j) = IEOR(key1(j),key2(j))
 endforall


end subroutine bitwise_xor


!==================================================================
subroutine bitwise_maskr(i,key3)
 implicit none
 integer,intent(in)  :: i
 integer,intent(out) :: key3(:)
!=====
 integer :: j,n
!=====

 n = SIZE(key3)

 do j=1,n
   key3(j) = MASKR( MAX(0 , MIN( nbit , i-nbit*(j-1)) ) )
 enddo

end subroutine bitwise_maskr


!==================================================================
function bitwise_popcnt(key1)
 implicit none
 integer,intent(in) :: key1(:)
 integer :: bitwise_popcnt
!=====
 integer :: j,n
!=====

 n = SIZE(key1)
 bitwise_popcnt = 0
 do j=1,n
   bitwise_popcnt = bitwise_popcnt + POPCNT(key1(j))
 enddo

end function bitwise_popcnt


!==================================================================
function bitwise_poppar(key1)
 implicit none
 integer,intent(in) :: key1(:)
 integer :: bitwise_poppar
!=====
 integer :: j,n
!=====

 n = SIZE(key1)
 bitwise_poppar = 0
 do j=1,n
   bitwise_poppar = bitwise_poppar + POPPAR(key1(j))
 enddo
 bitwise_poppar = MODULO(bitwise_poppar,2)

end function bitwise_poppar


!==================================================================
function bitwise_trailz(key1)
 implicit none
 integer,intent(in) :: key1(:)
 integer :: bitwise_trailz
!=====
 integer :: itmp,j,n
!=====

 n = SIZE(key1)
 bitwise_trailz = 0
 do j=1,n
   itmp = TRAILZ(key1(j))
   bitwise_trailz = bitwise_trailz + itmp
   if( itmp < nbit ) exit
 enddo

end function bitwise_trailz


!==================================================================
function bitwise_leadz(key1)
 implicit none
 integer,intent(in) :: key1(:)
 integer :: bitwise_leadz
!=====
 integer :: itmp,j,n
!=====

 n = SIZE(key1)
 bitwise_leadz = 0
 do j=n,1,-1
   itmp = LEADZ(key1(j))
   bitwise_leadz = bitwise_leadz + itmp
   if( itmp < nbit ) exit
 enddo

end function bitwise_leadz







end module m_bitwise
!==================================================================
