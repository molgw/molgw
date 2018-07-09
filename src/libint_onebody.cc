/*
 * This file is part of MOLGW.
 * C/C++ wrapper to the libint library
 * One-body terms: overlap, kinetic, nuclear attraction
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



#if defined(HAVE_LIBINT_ONEBODY)

/* ==========================================================================
 *                           Overlap
 * ========================================================================== */

extern "C" {
void libint_overlap(int amA, int contrdepthA , double A [] , double alphaA [], double cA [],
                    int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                    double overlapAB [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_overlap_t * inteval = libint2::malloc<Libint_overlap_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_overlap);
 assert(amB <= LIBINT2_MAX_AM_overlap);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = max(amA,amB) ;
 const int am = amA + amB ;


 LIBINT2_PREFIXED_NAME(libint2_init_overlap)(inteval, ammax, 0);

 double alphaP, ksiP ;
 double P[3];
 double AB_x = A[0] - B[0];
 double AB_y = A[1] - B[1];
 double AB_z = A[2] - B[2];

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     Libint_overlap_t* int12 = &inteval[icontrdepth2] ;
     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     int12->AB_x[0] = AB_x ;
     int12->AB_y[0] = AB_y ;
     int12->AB_z[0] = AB_z ;
#ifdef LIBINT2_DEFINED(eri, BA_x) && LIBINT2_DEFINED(eri, BA_y) && LIBINT2_DEFINED(eri, BA_z)
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

     icontrdepth2++ ;
   }
 }




 if( amA + amB == 0 ) {

   overlapAB[0] = 0.0 ;
   for( int icontrdepth2=0; icontrdepth2 < contrdepth2; icontrdepth2++) {
     overlapAB[0] +=   inteval[icontrdepth2]._0_Overlap_0_x[0]
                     * inteval[icontrdepth2]._0_Overlap_0_y[0]
                     * inteval[icontrdepth2]._0_Overlap_0_z[0] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_overlap)[amA][amB](inteval);
   for( int i12=0; i12 < nint(amA) * nint(amB) ; ++i12 ) {
     overlapAB[i12] = inteval[0].targets[0][i12] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_overlap)(inteval);
 free(inteval);

}
}


/* ==========================================================================
 *                           Kinetic
 * ========================================================================== */

extern "C" {
void libint_kinetic(int amA, int contrdepthA , double A [] , double alphaA [], double cA [],
                    int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                    double kineticAB [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_kinetic_t * inteval = libint2::malloc<Libint_kinetic_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_kinetic);
 assert(amB <= LIBINT2_MAX_AM_kinetic);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = max(amA,amB) ;
 const int am = amA + amB ;


 LIBINT2_PREFIXED_NAME(libint2_init_kinetic)(inteval, ammax, 0);

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

     Libint_kinetic_t* int12 = &inteval[icontrdepth2] ;
     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     int12->AB_x[0] = AB_x ;
     int12->AB_y[0] = AB_y ;
     int12->AB_z[0] = AB_z ;
     int12->BA_x[0] = -AB_x ;
     int12->BA_y[0] = -AB_y ;
     int12->BA_z[0] = -AB_z ;
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



 if( amA + amB == 0 ) {

   kineticAB[0] = 0.0 ;
   for( int icontrdepth2=0; icontrdepth2 < contrdepth2; icontrdepth2++) {
     kineticAB[0] +=   inteval[icontrdepth2]._0_Overlap_0_x[0]
                     * inteval[icontrdepth2]._0_Overlap_0_y[0]
                     * inteval[icontrdepth2]._0_Overlap_0_z[0]
                     * pfac[icontrdepth2] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_kinetic)[amA][amB](inteval);
   for( int i12=0; i12 < nint(amA) * nint(amB) ; ++i12 ) {
     kineticAB[i12] = inteval[0].targets[0][i12] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_kinetic)(inteval);
 free(inteval);

}
}


/* ==========================================================================
 *                           ElecPot
 * ========================================================================== */

extern "C" {
void libint_elecpot(int amA, int contrdepthA , double A [] , double alphaA [], double cA [],
                    int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                    double C [],
                    double elecpotAB [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_elecpot_t * inteval = libint2::malloc<Libint_elecpot_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_elecpot);
 assert(amB <= LIBINT2_MAX_AM_elecpot);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = max(amA,amB) ;
 const int am = amA + amB ;
 double U ;
 double F[am+1] ;


 LIBINT2_PREFIXED_NAME(libint2_init_elecpot)(inteval, ammax, 0);

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

     Libint_elecpot_t* int12 = &inteval[icontrdepth2] ;

     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     ksiP = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     int12->AB_x[0] = AB_x ;
     int12->AB_y[0] = AB_y ;
     int12->AB_z[0] = AB_z ;
     int12->BA_x[0] = -AB_x ;
     int12->BA_y[0] = -AB_y ;
     int12->BA_z[0] = -AB_z ;

     int12->PA_x[0] = P[0] - A[0] ;
     int12->PA_y[0] = P[1] - A[1] ;
     int12->PA_z[0] = P[2] - A[2] ;
     int12->PB_x[0] = P[0] - B[0] ;
     int12->PB_y[0] = P[1] - B[1] ;
     int12->PB_z[0] = P[2] - B[2] ;
     int12->PC_x[0] = P[0] - C[0] ;
     int12->PC_y[0] = P[1] - C[1] ;
     int12->PC_z[0] = P[2] - C[2] ;


     int12->_0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_x * AB_x ) * cA[icontrdepthA] * cB[icontrdepthB];
     int12->_0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_y * AB_y );
     int12->_0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) * exp( -ksiP * AB_z * AB_z );

     int12->oo2z[0] = 0.5 / alphaP ;

     int12->veclen = 1 ;
     int12->contrdepth = contrdepth2 ;

     U = alphaP * ( ( C[0] - P[0] ) * ( C[0] - P[0] ) + ( C[1] - P[1] ) * ( C[1] - P[1] ) + ( C[2] - P[2] ) * ( C[2] - P[2] ) ) ;
     boys_function_c(F, am, U);

     pfac = 2.0 * ( M_PI / alphaP ) * exp( - ksiP * ab2 ) * cA[icontrdepthA] * cB[icontrdepthB] ;

     double* elecpot = int12->LIBINT_T_S_ELECPOT_S(0);
     for( int l=0; l <= am ; ++l , ++elecpot ) *elecpot = pfac * F[l] ;

     icontrdepth2++ ;
   }
 }



 if( am == 0 ) {

   elecpotAB[0] = 0.0 ;
   for( int icontrdepth2=0; icontrdepth2 < contrdepth2; icontrdepth2++) {
     elecpotAB[0] +=   inteval[icontrdepth2].LIBINT_T_S_ELECPOT_S(0)[0] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_elecpot)[amA][amB](inteval);
   for( int i12=0; i12 < nint(amA) * nint(amB) ; ++i12 ) {
     elecpotAB[i12] = inteval[0].targets[0][i12] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_elecpot)(inteval);
 free(inteval);

}
}


#endif
/* ========================================================================== */
