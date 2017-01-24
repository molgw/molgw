/* MOLGW */
#include<libint2.h>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<assert.h>
using namespace std;


int nint(int am);
extern "C" void boys_function_c(double*, int, double);


/*                    
 *                           Overlap
*/

extern "C" {
void libint_overlap(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                    int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                    double overlapAB [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_overlap_t * inteval = libint2::malloc<Libint_overlap_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_eri);
 assert(amB <= LIBINT2_MAX_AM_eri);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = LIBINT2_MAX_AM_eri ;


 LIBINT2_PREFIXED_NAME(libint2_init_overlap)(inteval, ammax, 0);

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

     inteval[icontrdepth2]._0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * cA[icontrdepthA] * cB[icontrdepthB] * pow(-1,amA+amB)
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_x[0] * inteval[icontrdepth2].AB_x[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_y[0] * inteval[icontrdepth2].AB_y[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_z[0] * inteval[icontrdepth2].AB_z[0] / alphaP );
 
     inteval[icontrdepth2].oo2z[0] = 0.5 / alphaP ;
    
     inteval[icontrdepth2].veclen = 1 ;
     inteval[icontrdepth2].contrdepth = contrdepth2 ;

     icontrdepth2++ ;
   }
 }



 if( amA + amB == 0 ) {

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

}
}


/*                    
 *                           ElecPot         
*/

extern "C" {
void libint_elecpot(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                    int amB, int contrdepthB , double B [] , double alphaB [], double cB [],
                    double C [],
                    double elecpotAB [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthB;
 Libint_elecpot_t * inteval = libint2::malloc<Libint_elecpot_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_eri);
 assert(amB <= LIBINT2_MAX_AM_eri);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
#endif

 const unsigned int ammax = LIBINT2_MAX_AM_eri ;
 int am = amA + amB ;
 double U ;
 double F[am+1] ;


 LIBINT2_PREFIXED_NAME(libint2_init_elecpot)(inteval, ammax, 0);

 double alphaP, ksiP ;
 double P[3];
 double ab2 ;

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

     inteval[icontrdepth2]._0_Overlap_0_x[0] = sqrt( M_PI / alphaP ) * cA[icontrdepthA] * cB[icontrdepthB] * pow(-1,amA+amB) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_x[0] * inteval[icontrdepth2].AB_x[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_y[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_y[0] * inteval[icontrdepth2].AB_y[0] / alphaP );
     inteval[icontrdepth2]._0_Overlap_0_z[0] = sqrt( M_PI / alphaP ) 
                                                * exp( - alphaA[icontrdepthA] * alphaB[icontrdepthB] * inteval[icontrdepth2].AB_z[0] * inteval[icontrdepth2].AB_z[0] / alphaP );
 
     inteval[icontrdepth2].oo2z[0] = 0.5 / alphaP ;
    
     inteval[icontrdepth2].veclen = 1 ;
     inteval[icontrdepth2].contrdepth = contrdepth2 ;

     U = alphaP * ( ( C[0] - P[0] ) * ( C[0] - P[0] ) + ( C[1] - P[1] ) * ( C[1] - P[1] ) + ( C[2] - P[2] ) * ( C[2] - P[2] ) ) ;
     boys_function_c(F, am, U);

     inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_0[0] = 2.0 * ( M_PI / alphaP ) * exp( - ksiP * ab2 ) * F[0] 
                                                                    * cA[icontrdepthA] * cB[icontrdepthB]
                                                                    * pow(-1,amA+amB) ;
     inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_1[0] = 2.0 * ( M_PI / alphaP ) * exp( - ksiP * ab2 ) * F[1] 
                                                                    * cA[icontrdepthA] * cB[icontrdepthB]
                                                                    * pow(-1,amA+amB) ;
     inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_2[0] = 2.0 * ( M_PI / alphaP ) * exp( - ksiP * ab2 ) * F[2] 
                                                                    * cA[icontrdepthA] * cB[icontrdepthB]
                                                                    * pow(-1,amA+amB) ;
     inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_3[0] = 2.0 * ( M_PI / alphaP ) * exp( - ksiP * ab2 ) * F[3] 
                                                                    * cA[icontrdepthA] * cB[icontrdepthB]
                                                                    * pow(-1,amA+amB) ;


     icontrdepth2++ ;
   }
 }



 if( amA + amB == 0 ) {

   for( int icontrdepth2=0; icontrdepth2 < contrdepth2; icontrdepth2++) {
     elecpotAB[0] +=   inteval[icontrdepth2]._aB_s___0___ElecPot_s___0___Ab__up_0[0] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_elecpot)[amA][amB](inteval);
   for( int i12=0; i12 < nint(amA) * nint(amB) ; ++i12 ) {
     elecpotAB[i12] += inteval[0].targets[0][i12] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_elecpot)(inteval);

}
}



