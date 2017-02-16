/*
 * This file is part of MOLGW.
 * C/C++ wrapper to the libint library
 * One-body gradient terms: overlap, kinetic, nuclear attraction
 * with contracted gaussians
 * Author: F. Bruneval
 */
#include<libint2.h>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<assert.h>
#include "libint_molgw.h"


/* Code */

#ifdef HAVE_LIBINT_ONEBODY

/* ==========================================================================                    
 *                           Overlap
 * ========================================================================== */

extern "C" {
void libint_overlap_grad(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                         int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                         double overlapABx [], double overlapABy [], double overlapABz []) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_overlap1_t * inteval = libint2::malloc<Libint_overlap1_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_overlap1);
 assert(amB <= LIBINT2_MAX_AM_overlap1);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = max(amA,amB) ;
 int am = amA + amB ;


 LIBINT2_PREFIXED_NAME(libint2_init_overlap1)(inteval, ammax, 0);

 double alphaP ;
 double P[3];

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     inteval[icontrdepth2].AB_x[0] = B[0] - A[0] ;
     inteval[icontrdepth2].AB_y[0] = B[1] - A[1] ;
     inteval[icontrdepth2].AB_z[0] = B[2] - A[2] ;
     inteval[icontrdepth2].PA_x[0] = A[0] - P[0] ;
     inteval[icontrdepth2].PA_y[0] = A[1] - P[1] ;
     inteval[icontrdepth2].PA_z[0] = A[2] - P[2] ;
     inteval[icontrdepth2].PB_x[0] = B[0] - P[0] ;
     inteval[icontrdepth2].PB_y[0] = B[1] - P[1] ;
     inteval[icontrdepth2].PB_z[0] = B[2] - P[2] ;


     inteval[icontrdepth2]._0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * cA[icontrdepthA] * cB[icontrdepthB] * pow(-1,am)
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_x[0] * inteval[icontrdepth2].AB_x[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_y[0] * inteval[icontrdepth2].AB_y[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_z[0] * inteval[icontrdepth2].AB_z[0] / alphaP );
 
     inteval[icontrdepth2].oo2z[0] = 0.5 / alphaP ;
    
     inteval[icontrdepth2].veclen = 1 ;
     inteval[icontrdepth2].contrdepth = contrdepth2 ;

     inteval[icontrdepth2].two_alpha0_bra[0] = 2.0 * alphaA[icontrdepthA];
     inteval[icontrdepth2].two_alpha0_ket[0] = 2.0 * alphaB[icontrdepthB];

     icontrdepth2++ ;
   }
 }





 LIBINT2_PREFIXED_NAME(libint2_build_overlap1)[amA][amB](inteval);
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   overlapABx[i12] = inteval[0].targets[0][i12] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   overlapABy[i12] = inteval[0].targets[0][i12+nint(amA) * nint(amB)] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   overlapABz[i12] = inteval[0].targets[0][i12+2*nint(amA)*nint(amB)] ;
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
 Libint_kinetic1_t * inteval = libint2::malloc<Libint_kinetic1_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_kinetic1);
 assert(amB <= LIBINT2_MAX_AM_kinetic1);

#ifndef LIBINT2_CONTRACTED_INTS
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

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     inteval[icontrdepth2].AB_x[0] = B[0] - A[0] ;
     inteval[icontrdepth2].AB_y[0] = B[1] - A[1] ;
     inteval[icontrdepth2].AB_z[0] = B[2] - A[2] ;
     inteval[icontrdepth2].PA_x[0] = A[0] - P[0] ;
     inteval[icontrdepth2].PA_y[0] = A[1] - P[1] ;
     inteval[icontrdepth2].PA_z[0] = A[2] - P[2] ;
     inteval[icontrdepth2].PB_x[0] = B[0] - P[0] ;
     inteval[icontrdepth2].PB_y[0] = B[1] - P[1] ;
     inteval[icontrdepth2].PB_z[0] = B[2] - P[2] ;

     ab2 =  inteval[icontrdepth2].AB_x[0] * inteval[icontrdepth2].AB_x[0] 
          + inteval[icontrdepth2].AB_y[0] * inteval[icontrdepth2].AB_y[0]
          + inteval[icontrdepth2].AB_z[0] * inteval[icontrdepth2].AB_z[0] ;

     pfac[icontrdepth2] = ksiP * ( 3.0 - 2.0 * ksiP * ab2 ) ;

     inteval[icontrdepth2]._0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * cA[icontrdepthA] * cB[icontrdepthB] * pow(-1,am)
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_x[0] * inteval[icontrdepth2].AB_x[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_y[0] * inteval[icontrdepth2].AB_y[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_z[0] * inteval[icontrdepth2].AB_z[0] / alphaP );

     inteval[icontrdepth2].oo2z[0] = 0.5 / alphaP ;
    
     inteval[icontrdepth2].veclen = 1 ;
     inteval[icontrdepth2].contrdepth = contrdepth2 ;

     inteval[icontrdepth2].two_alpha0_bra[0] = 2.0 * alphaA[icontrdepthA];
     inteval[icontrdepth2].two_alpha0_ket[0] = 2.0 * alphaB[icontrdepthB];

     icontrdepth2++ ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_build_kinetic1)[amA][amB](inteval);

 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   kineticABx[i12] = inteval[0].targets[0][i12] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   kineticABy[i12] = inteval[0].targets[0][i12+nint(amA)*nint(amB)] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   kineticABz[i12] = inteval[0].targets[0][i12+2*nint(amA)*nint(amB)] ;
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

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_elecpot1_t * inteval = libint2::malloc<Libint_elecpot1_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_elecpot1);
 assert(amB <= LIBINT2_MAX_AM_elecpot1);

#ifndef LIBINT2_CONTRACTED_INTS
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

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     inteval[icontrdepth2].AB_x[0] = B[0] - A[0] ;
     inteval[icontrdepth2].AB_y[0] = B[1] - A[1] ;
     inteval[icontrdepth2].AB_z[0] = B[2] - A[2] ;
     inteval[icontrdepth2].PA_x[0] = A[0] - P[0] ;
     inteval[icontrdepth2].PA_y[0] = A[1] - P[1] ;
     inteval[icontrdepth2].PA_z[0] = A[2] - P[2] ;
     inteval[icontrdepth2].PC_x[0] = C[0] - P[0] ;
     inteval[icontrdepth2].PC_y[0] = C[1] - P[1] ;
     inteval[icontrdepth2].PC_z[0] = C[2] - P[2] ;

     ab2 =  inteval[icontrdepth2].AB_x[0] * inteval[icontrdepth2].AB_x[0] 
          + inteval[icontrdepth2].AB_y[0] * inteval[icontrdepth2].AB_y[0]
          + inteval[icontrdepth2].AB_z[0] * inteval[icontrdepth2].AB_z[0] ;

     inteval[icontrdepth2]._0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * cA[icontrdepthA] * cB[icontrdepthB] * pow(-1,am) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_x[0] * inteval[icontrdepth2].AB_x[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_y[0] * inteval[icontrdepth2].AB_y[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_z[0] * inteval[icontrdepth2].AB_z[0] / alphaP );
 
     inteval[icontrdepth2].oo2z[0] = 0.5 / alphaP ;
    
     inteval[icontrdepth2].veclen = 1 ;
     inteval[icontrdepth2].contrdepth = contrdepth2 ;

     inteval[icontrdepth2].two_alpha0_bra[0] = 2.0 * alphaA[icontrdepthA] ;
     inteval[icontrdepth2].two_alpha0_ket[0] = 2.0 * alphaB[icontrdepthB] ;
     inteval[icontrdepth2].rho12_over_alpha1[0] = alphaB[icontrdepthA] / alphaP ;
#if LIBINT2_DEFINED(eri, rho12_over_alpha2)
     inteval[icontrdepth2].rho12_over_alpha2[0] = alphaA[icontrdepthB] / alphaP ;
#endif

     U = alphaP * ( ( C[0] - P[0] ) * ( C[0] - P[0] ) + ( C[1] - P[1] ) * ( C[1] - P[1] ) + ( C[2] - P[2] ) * ( C[2] - P[2] ) ) ;
     boys_function_c(F, am+1, U);

     pfac = 2.0 * ( M_PI / alphaP ) * exp( - ksiP * ab2 ) * cA[icontrdepthA] * cB[icontrdepthB] * pow(-1,am) ;

     /* FIXME what happens when F is evaluated beyond its range? */
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_0
     inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_0[0] = pfac * F[0] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_1
     inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_1[0] = pfac * F[1] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_2
     if( am >=1 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_2[0] = pfac * F[2] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_3
     if( am >=2 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_3[0] = pfac * F[3] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_4
     if( am >=3 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_4[0] = pfac * F[4] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_5
     if( am >=4 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_5[0] = pfac * F[5] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_6
     if( am >=5 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_6[0] = pfac * F[6] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_7
     if( am >=6 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_7[0] = pfac * F[7] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_8
     if( am >=7 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_8[0] = pfac * F[8] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_9
     if( am >=8 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_9[0] = pfac * F[9] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_10
     if( am >=9 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_10[0] = pfac * F[10] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_11
     if( am >=10 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_11[0] = pfac * F[11] ;
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_12
     if( am >=11 ) inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_12[0] = pfac * F[12] ;
#endif


     icontrdepth2++ ;
   }
 }


 LIBINT2_PREFIXED_NAME(libint2_build_elecpot1)[amA][amB](inteval);

 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   elecpotAx[i12]+= inteval[0].targets[0][i12] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   elecpotAy[i12]+= inteval[0].targets[0][i12+nint(amA)*nint(amB)*1] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   elecpotAz[i12]+= inteval[0].targets[0][i12+nint(amA)*nint(amB)*2] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   elecpotBx[i12]+= inteval[0].targets[0][i12+nint(amA)*nint(amB)*3] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   elecpotBy[i12]+= inteval[0].targets[0][i12+nint(amA)*nint(amB)*4] ;
 }
 for( int i12=0; i12 < nint(amA) * nint(amB) ; i12++ ) {
   elecpotBz[i12]+= inteval[0].targets[0][i12+nint(amA)*nint(amB)*5] ;
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_elecpot1)(inteval);
 free(inteval);

}
}

#endif

