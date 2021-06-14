/*
 * This file is part of MOLGW.
 * C/C++ wrapper to the libint library
 * One-body gradient terms: overlap, kinetic, nuclear attraction
 * with contracted gaussians
 * Author: F. Bruneval
 */
#include<libint2.h>
#include<libint2/util/memory.h>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<assert.h>
#include "libint_molgw.h"


/* Code */

#if defined(LIBINT2_SUPPORT_ONEBODY) && (LIBINT2_DERIV_ONEBODY_ORDER > 0)

/* ==========================================================================
 *                           Overlap
 * ========================================================================== */

extern "C" {
void libint_overlap_grad(int amA, int contrdepthA , double A [] , double alphaA [], double cA [],
                         int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                         double overlapABx [], double overlapABy [], double overlapABz []) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 const unsigned int ni = nint(amA) * nint(amB) ;
 Libint_overlap1_t * inteval = libint2::malloc<Libint_overlap1_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_overlap1);
 assert(amB <= LIBINT2_MAX_AM_overlap1);

#if !defined(LIBINT2_CONTRACTED_INTS)
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = max(amA,amB) ;
 int am = amA + amB ;


 LIBINT2_PREFIXED_NAME(libint2_init_overlap1)(inteval, ammax, 0);

 double alphaP, ksiP ;
 double P[3];
 double AB_x = A[0] - B[0];
 double AB_y = A[1] - B[1];
 double AB_z = A[2] - B[2];

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     Libint_overlap1_t* int12 = &inteval[icontrdepth2];
     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     int12->AB_x[0] = AB_x ;
     int12->AB_y[0] = AB_y ;
     int12->AB_z[0] = AB_z ;
#if LIBINT2_DEFINED(eri, BA_x) && LIBINT2_DEFINED(eri, BA_y) && LIBINT2_DEFINED(eri, BA_z)
     int12->BA_x[0] = -AB_x ;
     int12->BA_y[0] = -AB_y ;
     int12->BA_z[0] = -AB_z ;
#endif
     int12->PA_x[0] = P[0] - A[0] ;
     int12->PA_y[0] = P[1] - A[1] ;
     int12->PA_z[0] = P[2] - A[2] ;
     int12->PB_x[0] = P[0] - B[0] ;
     int12->PB_y[0] = P[1] - B[1] ;
     int12->PB_z[0] = P[2] - B[2] ;


     int12->_0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_x * AB_x ) * cA[icontrdepthA] * cB[icontrdepthB];
     int12->_0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_y * AB_y );
     int12->_0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_z * AB_z );

     int12->oo2z[0] = 0.5 / alphaP ;

     int12->veclen = 1 ;
     int12->contrdepth = contrdepth2 ;

     int12->two_alpha0_bra[0] = 2.0 * alphaA[icontrdepthA];
     int12->two_alpha0_ket[0] = 2.0 * alphaB[icontrdepthB];

     icontrdepth2++ ;
   }
 }





 LIBINT2_PREFIXED_NAME(libint2_build_overlap1)[amA][amB](inteval);
 for( int i12=0; i12 < ni ; i12++ ) {
   overlapABx[i12] = inteval[0].targets[0][i12+0*ni] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   overlapABy[i12] = inteval[0].targets[0][i12+1*ni] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   overlapABz[i12] = inteval[0].targets[0][i12+2*ni] ;
 }


 LIBINT2_PREFIXED_NAME(libint2_cleanup_overlap1)(inteval);
 free(inteval);

}
}


/* ==========================================================================
 *                           Kinetic
 * ========================================================================== */

extern "C" {
void libint_kinetic_grad(int amA, int contrdepthA , double A [] , double alphaA [], double cA [],
                         int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                         double kineticABx [], double kineticABy [], double kineticABz [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 const unsigned int ni = nint(amA) * nint(amB) ;
 Libint_kinetic1_t * inteval = libint2::malloc<Libint_kinetic1_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_kinetic1);
 assert(amB <= LIBINT2_MAX_AM_kinetic1);

#if !defined(LIBINT2_CONTRACTED_INTS)
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = max(amA,amB) ;
 int am = amA + amB ;


 LIBINT2_PREFIXED_NAME(libint2_init_kinetic1)(inteval, ammax, 0);

 double alphaP, ksiP ;
 double P[3];
 double ab2 ;
 double pfac[contrdepth2] ;
 double AB_x = A[0] - B[0];
 double AB_y = A[1] - B[1];
 double AB_z = A[2] - B[2];

 ab2 =  AB_x * AB_x + AB_y * AB_y + AB_z * AB_z ;

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     Libint_kinetic1_t* int12 = &inteval[icontrdepth2];
     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     int12->AB_x[0] = AB_x ;
     int12->AB_y[0] = AB_y ;
     int12->AB_z[0] = AB_z ;
#if LIBINT2_DEFINED(eri, BA_x) && LIBINT2_DEFINED(eri, BA_y) && LIBINT2_DEFINED(eri, BA_z)
     int12->BA_x[0] = -AB_x ;
     int12->BA_y[0] = -AB_y ;
     int12->BA_z[0] = -AB_z ;
#endif
     int12->PA_x[0] = P[0] - A[0] ;
     int12->PA_y[0] = P[1] - A[1] ;
     int12->PA_z[0] = P[2] - A[2] ;
     int12->PB_x[0] = P[0] - B[0] ;
     int12->PB_y[0] = P[1] - B[1] ;
     int12->PB_z[0] = P[2] - B[2] ;


     pfac[icontrdepth2] = ksiP * ( 3.0 - 2.0 * ksiP * ab2 ) ;

     int12->_0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_x * AB_x ) * cA[icontrdepthA] * cB[icontrdepthB];
     int12->_0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_y * AB_y );
     int12->_0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_z * AB_z );

     int12->oo2z[0] = 0.5 / alphaP ;

     int12->veclen = 1 ;
     int12->contrdepth = contrdepth2 ;

     int12->two_alpha0_bra[0] = 2.0 * alphaA[icontrdepthA];
     int12->two_alpha0_ket[0] = 2.0 * alphaB[icontrdepthB];

     icontrdepth2++ ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_build_kinetic1)[amA][amB](inteval);

 for( int i12=0; i12 < ni ; i12++ ) {
   kineticABx[i12] = inteval[0].targets[0][i12+0*ni] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   kineticABy[i12] = inteval[0].targets[0][i12+1*ni] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   kineticABz[i12] = inteval[0].targets[0][i12+2*ni] ;
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_kinetic1)(inteval);
 free(inteval);

}
}


/* ==========================================================================
 *                           ElecPot
 * ========================================================================== */

extern "C" {
void libint_elecpot_grad(int amA, int contrdepthA , double A [] , double alphaA [], double cA [],
                         int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                         double C [],
                         double elecpotAx [], double elecpotAy [], double elecpotAz [],
                         double elecpotBx [], double elecpotBy [], double elecpotBz []) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB ;
 const unsigned int ni = nint(amA) * nint(amB) ;
 Libint_elecpot1_t * inteval = libint2::malloc<Libint_elecpot1_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_elecpot1);
 assert(amB <= LIBINT2_MAX_AM_elecpot1);

#if !defined(LIBINT2_CONTRACTED_INTS)
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = max(amA,amB) ;
 int am = amA + amB ;
 double U ;
 double F[am+2] ;


 LIBINT2_PREFIXED_NAME(libint2_init_elecpot1)(inteval, ammax, 0);

 double alphaP, ksiP ;
 double P[3];
 double ab2, pfac ;
 double AB_x = A[0] - B[0];
 double AB_y = A[1] - B[1];
 double AB_z = A[2] - B[2];

 ab2 =  AB_x * AB_x + AB_y * AB_y + AB_z * AB_z ;

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     Libint_elecpot1_t* int12 = &inteval[icontrdepth2]  ;

     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     int12->AB_x[0] = AB_x ;
     int12->AB_y[0] = AB_y ;
     int12->AB_z[0] = AB_z ;
#if LIBINT2_DEFINED(eri, BA_x) && LIBINT2_DEFINED(eri, BA_y) && LIBINT2_DEFINED(eri, BA_z)
     int12->BA_x[0] = -AB_x ;
     int12->BA_y[0] = -AB_y ;
     int12->BA_z[0] = -AB_z ;
#endif
     int12->PA_x[0] = P[0] - A[0] ;
     int12->PA_y[0] = P[1] - A[1] ;
     int12->PA_z[0] = P[2] - A[2] ;
     int12->PC_x[0] = P[0] - C[0] ;
     int12->PC_y[0] = P[1] - C[1] ;
     int12->PC_z[0] = P[2] - C[2] ;


     int12->_0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_x * AB_x ) * cA[icontrdepthA] * cB[icontrdepthB];
     int12->_0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_y * AB_y );
     int12->_0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_z * AB_z );

     int12->oo2z[0] = 0.5 / alphaP ;

     int12->veclen = 1 ;
     int12->contrdepth = contrdepth2 ;

     int12->two_alpha0_bra[0] = 2.0 * alphaA[icontrdepthA] ;
     int12->two_alpha0_ket[0] = 2.0 * alphaB[icontrdepthB] ;
     int12->rho12_over_alpha1[0] = alphaB[icontrdepthB] / alphaP ;
#if LIBINT2_DEFINED(eri, rho12_over_alpha2)
     int12->rho12_over_alpha2[0] = alphaA[icontrdepthA] / alphaP ;
#endif


     U = alphaP * ( ( C[0] - P[0] ) * ( C[0] - P[0] ) + ( C[1] - P[1] ) * ( C[1] - P[1] ) + ( C[2] - P[2] ) * ( C[2] - P[2] ) ) ;
     boys_function_c(F, am+1, U);

     pfac = 2.0 * ( M_PI / alphaP ) * exp( - ksiP * ab2 ) * cA[icontrdepthA] * cB[icontrdepthB] ;

     double* elecpot = int12->LIBINT_T_S_ELECPOT_S(0);
     for( int l=0; l <= am+1 ; ++l , ++elecpot ) *elecpot = pfac * F[l] ;


     icontrdepth2++ ;
   }
 }


 LIBINT2_PREFIXED_NAME(libint2_build_elecpot1)[amA][amB](inteval);

 for( int i12=0; i12 < ni ; i12++ ) {
   elecpotAx[i12] = inteval[0].targets[0][i12+ni*0] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   elecpotAy[i12] = inteval[0].targets[0][i12+ni*1] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   elecpotAz[i12] = inteval[0].targets[0][i12+ni*2] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   elecpotBx[i12] = inteval[0].targets[0][i12+ni*3] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   elecpotBy[i12] = inteval[0].targets[0][i12+ni*4] ;
 }
 for( int i12=0; i12 < ni ; i12++ ) {
   elecpotBz[i12] = inteval[0].targets[0][i12+ni*5] ;
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_elecpot1)(inteval);
 free(inteval);

}
}

#endif
/* ========================================================================== */
