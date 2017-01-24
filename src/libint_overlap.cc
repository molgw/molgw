/* MOLGW */
#include<libint2.h>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<assert.h>
using namespace std;

int nint(int am);

extern "C" {
void libint_overlap(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                    int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                    double overlapAB [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_overlap_t * overlap = libint2::malloc<Libint_overlap_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_eri);
 assert(amB <= LIBINT2_MAX_AM_eri);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = LIBINT2_MAX_AM_eri ;

// LIBINT2_PREFIXED_NAME(libint2_static_init)();

 LIBINT2_PREFIXED_NAME(libint2_init_overlap)(overlap, ammax, 0);



 double alphaP ;
 double P[3];

 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {

     alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
     P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0] ) / alphaP ;
     P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1] ) / alphaP ;
     P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2] ) / alphaP ;


     overlap[icontrdepth2].AB_x[0] = B[0] - A[0] ;
     overlap[icontrdepth2].AB_y[0] = B[1] - A[1] ;
     overlap[icontrdepth2].AB_z[0] = B[2] - B[2] ;
     overlap[icontrdepth2].PA_x[0] = A[0] - P[0] ;
     overlap[icontrdepth2].PA_y[0] = A[1] - P[1] ;
     overlap[icontrdepth2].PA_z[0] = A[2] - P[2] ;
     overlap[icontrdepth2].PB_x[0] = B[0] - P[0] ;
     overlap[icontrdepth2].PB_y[0] = B[1] - P[1] ;
     overlap[icontrdepth2].PB_z[0] = B[2] - P[2] ;

     overlap[icontrdepth2]._0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * cA[icontrdepthA] * cB[icontrdepthB] * pow(-1,amA+amB)
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * overlap[icontrdepth2].AB_x[0] * overlap[icontrdepth2].AB_x[0] / alphaP );
     overlap[icontrdepth2]._0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * overlap[icontrdepth2].AB_y[0] * overlap[icontrdepth2].AB_y[0] / alphaP );
     overlap[icontrdepth2]._0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * overlap[icontrdepth2].AB_z[0] * overlap[icontrdepth2].AB_z[0] / alphaP );
 
     overlap[icontrdepth2].oo2z[0] = 0.5 / alphaP ;
    
     overlap[icontrdepth2].veclen = 1 ;
     overlap[icontrdepth2].contrdepth = contrdepth2 ;

     icontrdepth2++ ;
   }
 }



 if( amA + amB == 0 ) {

   for( int icontrdepth2=0; icontrdepth2 < contrdepth2; icontrdepth2++) {
     overlapAB[0] +=   overlap[icontrdepth2]._0_Overlap_0_x[0]
                     * overlap[icontrdepth2]._0_Overlap_0_y[0]
                     * overlap[icontrdepth2]._0_Overlap_0_z[0] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_overlap)[amA][amB](overlap);
   for( int i12=0; i12 < nint(amA) * nint(amB) ; ++i12 ) {
     overlapAB[i12] = overlap[0].targets[0][i12] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_overlap)(overlap);



// LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();


}
}
