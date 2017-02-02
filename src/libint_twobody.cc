/*
 * This file is part of MOLGW.
 * C/C++ wrapper to the libint library
 * two-body terms: 2-center Coulomb
 * with contracted gaussians
 * Author: F. Bruneval
 */
#include<libint2.h>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<assert.h>
using namespace std;


/* Headers */
extern "C" void boys_function_c(double*, int, double);

/* Code */

/* ==========================================================================                    
 *                        LIBINT initialization
 * ========================================================================== */
extern "C" {
void libint_init(int *ammax, bool *has_onebody, bool *has_gradient) {
 LIBINT2_PREFIXED_NAME(libint2_static_init)();
 *has_onebody  = LIBINT2_SUPPORT_ONEBODY ;
 *has_gradient = ( LIBINT2_DERIV_ERI_ORDER > 0 );
 *ammax = LIBINT2_MAX_AM_eri ;
}
}

/* Number of cartesian function for a given angular momentum */
int inline nint(int am) {
  return (am+1)*(am+2)/2;
}


/* ==========================================================================                    
 *                           2-center integrals
 * ========================================================================== */

extern "C" {
void libint_2center(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                    int amC, int contrdepthC , double C [] , double alphaC [], double cC [],
                    double rcut,
                    double eriAC [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthC;
 Libint_2eri_t* inteval = libint2::malloc<Libint_2eri_t>(contrdepth2);

 assert(amA <= LIBINT2_MAX_AM_2eri);
 assert(amC <= LIBINT2_MAX_AM_2eri);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthC == 1 );
#endif

 const unsigned int ammax = max(amA,amC) ;
 const int am = amA + amC ;


 LIBINT2_PREFIXED_NAME(libint2_init_2eri)(inteval, ammax, 0);

 Libint_2eri_t* int12 ;
 double alphaP, alphaQ ;
 double U ;
 double F[am+1] ;
 double gammapq ;
 double gammapq_rc2 ;
 double gammapq_ratio ;
 double PQx ;
 double PQy ;
 double PQz ;
 double PQ2 ;
 double Wx ;
 double Wy ;
 double Wz ;
 double pfac ;
 constexpr double pi_2p5 = pow(M_PI,2.5) ;


 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthC=0; icontrdepthC < contrdepthC; icontrdepthC++)  {

     int12 = &inteval[icontrdepth2] ;
     alphaP = alphaA[icontrdepthA] ;
     alphaQ = alphaC[icontrdepthC] ;

     int12->PA_x[0] = 0 ;
     int12->PA_y[0] = 0 ;
     int12->PA_z[0] = 0 ;
     int12->PB_x[0] = 0 ;
     int12->PB_y[0] = 0 ;
     int12->PB_z[0] = 0 ;
     int12->AB_x[0] = 0 ;
     int12->AB_y[0] = 0 ;
     int12->AB_z[0] = 0 ;
     int12->oo2z[0] = 0.5 / alphaP ;

     gammapq = alphaP * alphaQ / (alphaP + alphaQ);

     /*
       Treat the LR integrals by simply modifying gammapq into gammapq_rc2
       And introducing gammapq_ratio = gammapq_rc2 / gammapq
     */

     gammapq_rc2   =   alphaP * alphaQ   / ( alphaP + alphaQ + alphaP * alphaQ * rcut * rcut );
     gammapq_ratio = ( alphaP + alphaQ ) / ( alphaP + alphaQ + alphaP * alphaQ * rcut * rcut );


     int12->QC_x[0] = 0 ;
     int12->QC_y[0] = 0 ;
     int12->QC_z[0] = 0 ;
     int12->CD_x[0] = 0 ;
     int12->CD_y[0] = 0 ;
     int12->CD_z[0] = 0 ;
     int12->oo2e[0] = 0.5 / alphaQ;


     PQx = A[0] - C[0];
     PQy = A[1] - C[1];
     PQz = A[2] - C[2];
     PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
     Wx = (alphaP * A[0] + alphaQ * C[0]) / (alphaP + alphaQ);
     Wy = (alphaP * A[1] + alphaQ * C[1]) / (alphaP + alphaQ);
     Wz = (alphaP * A[2] + alphaQ * C[2]) / (alphaP + alphaQ);

     int12->WP_x[0] = Wx - A[0];
     int12->WP_y[0] = Wy - A[1];
     int12->WP_z[0] = Wz - A[2];
     int12->WQ_x[0] = Wx - C[0];
     int12->WQ_y[0] = Wy - C[1];
     int12->WQ_z[0] = Wz - C[2];
     int12->oo2ze[0] = 0.5 / ( alphaP + alphaQ ) ;
#if LIBINT2_DEFINED(eri,roz)
     int12->roz[0] = gammapq/alphaP;
#endif
#if LIBINT2_DEFINED(eri,roe)
     int12->roe[0] = gammapq/alphaQ;
#endif

     pfac = 2 * pi_2p5 / (alphaP * alphaQ * sqrt(alphaP + alphaQ)) * cA[icontrdepthA] * cC[icontrdepthC] ;
     U = PQ2 * gammapq_rc2 ;
   
     boys_function_c(F, am, U);

     // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
     int12->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0] * pow(gammapq_ratio,0.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
     if( am > 0 ) int12->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1] * pow(gammapq_ratio,1.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
     if( am > 1 ) int12->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2] * pow(gammapq_ratio,2.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
     if( am > 2 ) int12->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3] * pow(gammapq_ratio,3.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
     if( am > 3 ) int12->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4] * pow(gammapq_ratio,4.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
     if( am > 4 ) int12->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5] * pow(gammapq_ratio,5.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
     if( am > 5 ) int12->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6] * pow(gammapq_ratio,6.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
     if( am > 6 ) int12->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7] * pow(gammapq_ratio,7.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
     if( am > 7 ) int12->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8] * pow(gammapq_ratio,8.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
     if( am > 8 ) int12->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9] * pow(gammapq_ratio,9.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
     if( am > 9 ) int12->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10] * pow(gammapq_ratio,10.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
     if( am > 10) int12->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11] * pow(gammapq_ratio,11.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
     if( am > 11) int12->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12] * pow(gammapq_ratio,12.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
     if( am > 12) int12->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13] * pow(gammapq_ratio,13.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
     if( am > 13) int12->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14] * pow(gammapq_ratio,14.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
     if( am > 14) int12->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15] * pow(gammapq_ratio,15.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
     if( am > 15) int12->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16] * pow(gammapq_ratio,16.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
     if( am > 16) int12->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17] * pow(gammapq_ratio,17.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
     if( am > 17) int12->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18] * pow(gammapq_ratio,18.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
     if( am > 18) int12->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19] * pow(gammapq_ratio,19.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
     if( am > 19) int12->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20] * pow(gammapq_ratio,20.5);
#endif


     int12->_0_Overlap_0_x[0] = 0.0 ;
     int12->_0_Overlap_0_y[0] = 0.0 ;
     int12->_0_Overlap_0_z[0] = 0.0 ;
 
     int12->veclen = 1 ;
     int12->contrdepth = contrdepth2 ;

     icontrdepth2++ ;
   }
 }




 if( am == 0 ) {

   eriAC[0] = 0.0 ;
   for( int icontrdepth2=0; icontrdepth2 < contrdepth2; icontrdepth2++) {
     eriAC[0] +=  inteval[icontrdepth2].LIBINT_T_SS_EREP_SS(0)[0] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_2eri)[amA][amC](inteval);
   for( int i12=0; i12 < nint(amA) * nint(amC) ; ++i12 ) {
     eriAC[i12] = inteval[0].targets[0][i12] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri)(inteval);

}
}



/* ==========================================================================                    
 *                           3-center integrals
 * ========================================================================== */


extern "C" {
void libint_3center(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                    int amC, int contrdepthC , double C [] , double alphaC [], double cC [], 
                    int amD, int contrdepthD , double D [] , double alphaD [], double cD [],
                    double rcut,
                    double eriACD [] ) {

 const unsigned int contrdepth3 = contrdepthA * contrdepthC * contrdepthD ;
 Libint_3eri_t* inteval = libint2::malloc<Libint_3eri_t>(contrdepth3);

 assert(amA <= LIBINT2_MAX_AM_3eri);
 assert(amC <= LIBINT2_MAX_AM_3eri);
 assert(amD <= LIBINT2_MAX_AM_3eri);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthC == 1 );
 assert( contrdepthD == 1 );
#endif

 const unsigned int ammax = max(max(amA,amC),amD) ;
 const int am = amA + amC + amD ;


 LIBINT2_PREFIXED_NAME(libint2_init_3eri)(inteval, ammax, 0);

 Libint_3eri_t* int12 ;
 double alphaP, alphaQ ;
 double P[3], Q[3] ;
 double U ;
 double F[am+1] ;
 double gammapq ;
 double gammapq_rc2 ;
 double gammapq_ratio ;
 double PQx ;
 double PQy ;
 double PQz ;
 double PQ2 ;
 double Wx ;
 double Wy ;
 double Wz ;
 double AB2 ;
 double CD2 ;
 double pfac ;
 constexpr double pi_2p5 = pow(M_PI,2.5) ;


 int icontrdepth3 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthC=0; icontrdepthC < contrdepthC; icontrdepthC++)  {
     for( int icontrdepthD=0; icontrdepthD < contrdepthD; icontrdepthD++)  {


       int12 = &inteval[icontrdepth3] ;
       alphaP = alphaA[icontrdepthA] ;
       alphaQ = alphaC[icontrdepthC] + alphaD[icontrdepthD];
       P[0] =  A[0] ;
       P[1] =  A[1] ;
       P[2] =  A[2] ;
       Q[0] = (alphaC[icontrdepthC] * C[0] + alphaD[icontrdepthD] * D[0]) / alphaQ ;
       Q[1] = (alphaC[icontrdepthC] * C[1] + alphaD[icontrdepthD] * D[1]) / alphaQ ;
       Q[2] = (alphaC[icontrdepthC] * C[2] + alphaD[icontrdepthD] * D[2]) / alphaQ ;

  
       int12->PA_x[0] = 0.0 ;
       int12->PA_y[0] = 0.0 ;
       int12->PA_z[0] = 0.0 ;
       int12->PB_x[0] = 0.0 ;
       int12->PB_y[0] = 0.0 ;
       int12->PB_z[0] = 0.0 ;
       int12->AB_x[0] = 0.0 ;
       int12->AB_y[0] = 0.0 ;
       int12->AB_z[0] = 0.0 ;
       int12->oo2z[0] = 0.5 / alphaP ;

       AB2 = 0.0 ;
  
       gammapq = alphaP * alphaQ / (alphaP + alphaQ);
  
       /*
         Treat the LR integrals by simply modifying gammapq into gammapq_rc2
         And introducing gammapq_ratio = gammapq_rc2 / gammapq
        */
  
       gammapq_rc2   =   alphaP * alphaQ   / ( alphaP + alphaQ + alphaP * alphaQ * rcut * rcut );
       gammapq_ratio = ( alphaP + alphaQ ) / ( alphaP + alphaQ + alphaP * alphaQ * rcut * rcut );
  
  
       int12->QC_x[0] = Q[0] - C[0] ;
       int12->QC_y[0] = Q[1] - C[1] ;
       int12->QC_z[0] = Q[2] - C[2] ;
       int12->CD_x[0] = C[0] - D[0] ;
       int12->CD_y[0] = C[1] - D[1] ;
       int12->CD_z[0] = C[2] - D[2] ;
       CD2 = int12->CD_x[0] * int12->CD_x[0] 
           + int12->CD_y[0] * int12->CD_y[0] 
           + int12->CD_z[0] * int12->CD_z[0] ;
       int12->oo2e[0] = 0.5 / alphaQ;
  

       PQx = P[0] - Q[0];
       PQy = P[1] - Q[1];
       PQz = P[2] - Q[2];
       PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
       Wx = (alphaP * P[0] + alphaQ * Q[0]) / (alphaP + alphaQ);
       Wy = (alphaP * P[1] + alphaQ * Q[1]) / (alphaP + alphaQ);
       Wz = (alphaP * P[2] + alphaQ * Q[2]) / (alphaP + alphaQ);
  
       int12->WP_x[0] = Wx - P[0];
       int12->WP_y[0] = Wy - P[1];
       int12->WP_z[0] = Wz - P[2];
       int12->WQ_x[0] = Wx - Q[0];
       int12->WQ_y[0] = Wy - Q[1];
       int12->WQ_z[0] = Wz - Q[2];
       int12->oo2ze[0] = 0.5 / ( alphaP + alphaQ ) ;
#if LIBINT2_DEFINED(eri,roz)
       int12->roz[0] = gammapq/alphaP;
#endif
#if LIBINT2_DEFINED(eri,roe)
       int12->roe[0] = gammapq/alphaQ;
#endif


       pfac = 2 * pi_2p5 / (alphaP * alphaQ * sqrt(alphaP + alphaQ))
                 * exp(-alphaC[icontrdepthC] * alphaD[icontrdepthD] * CD2 / alphaQ) 
                  * cA[icontrdepthA] * cC[icontrdepthC] * cD[icontrdepthD] ;
       U = PQ2 * gammapq_rc2 ;
   
       boys_function_c(F, am, U);

       // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
       int12->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0] * pow(gammapq_ratio,0.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
       if( am > 0 ) int12->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1] * pow(gammapq_ratio,1.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
       if( am > 1 ) int12->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2] * pow(gammapq_ratio,2.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
       if( am > 2 ) int12->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3] * pow(gammapq_ratio,3.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
       if( am > 3 ) int12->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4] * pow(gammapq_ratio,4.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
       if( am > 4 ) int12->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5] * pow(gammapq_ratio,5.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
       if( am > 5 ) int12->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6] * pow(gammapq_ratio,6.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
       if( am > 6 ) int12->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7] * pow(gammapq_ratio,7.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
       if( am > 7 ) int12->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8] * pow(gammapq_ratio,8.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
       if( am > 8 ) int12->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9] * pow(gammapq_ratio,9.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
       if( am > 9 ) int12->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10] * pow(gammapq_ratio,10.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
       if( am > 10) int12->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11] * pow(gammapq_ratio,11.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
       if( am > 11) int12->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12] * pow(gammapq_ratio,12.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
       if( am > 12) int12->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13] * pow(gammapq_ratio,13.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
       if( am > 13) int12->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14] * pow(gammapq_ratio,14.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
       if( am > 14) int12->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15] * pow(gammapq_ratio,15.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
       if( am > 15) int12->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16] * pow(gammapq_ratio,16.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
       if( am > 16) int12->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17] * pow(gammapq_ratio,17.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
       if( am > 17) int12->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18] * pow(gammapq_ratio,18.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
       if( am > 18) int12->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19] * pow(gammapq_ratio,19.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
       if( am > 19) int12->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20] * pow(gammapq_ratio,20.5);
#endif


       int12->_0_Overlap_0_x[0] = 0.0 ;
       int12->_0_Overlap_0_y[0] = 0.0 ;
       int12->_0_Overlap_0_z[0] = 0.0 ;
 
       int12->veclen = 1 ;
       int12->contrdepth = contrdepth3 ;

       icontrdepth3++ ;
     }
   }
 }



 if( am == 0 ) {

   eriACD[0] = 0.0 ;
   for( int icontrdepth3=0; icontrdepth3 < contrdepth3; icontrdepth3++) {
     eriACD[0] +=  inteval[icontrdepth3].LIBINT_T_SS_EREP_SS(0)[0] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_3eri)[amA][amC][amD](inteval);
   for( int i123=0; i123 < nint(amA) * nint(amC) * nint(amD) ; ++i123 ) {
     eriACD[i123] = inteval[0].targets[0][i123] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_3eri)(inteval);

}
}


/* ==========================================================================                    
 *                           4-center integrals
 * ========================================================================== */


extern "C" {
void libint_4center(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                    int amB, int contrdepthB , double B [] , double alphaB [], double cB [], 
                    int amC, int contrdepthC , double C [] , double alphaC [], double cC [], 
                    int amD, int contrdepthD , double D [] , double alphaD [], double cD [],
                    double rcut,
                    double eriABCD [] ) {

 const unsigned int contrdepth4 = contrdepthA * contrdepthB * contrdepthC * contrdepthD ;
 Libint_eri_t* inteval = libint2::malloc<Libint_eri_t>(contrdepth4);

 assert(amA <= LIBINT2_MAX_AM_eri);
 assert(amB <= LIBINT2_MAX_AM_eri);
 assert(amC <= LIBINT2_MAX_AM_eri);
 assert(amD <= LIBINT2_MAX_AM_eri);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
 assert( contrdepthC == 1 );
 assert( contrdepthD == 1 );
#endif

 const unsigned int ammax = max(max(max(amA,amB),amC),amD) ;
 const int am = amA + amB + amC + amD ;


 LIBINT2_PREFIXED_NAME(libint2_init_eri)(inteval, ammax, 0);

 Libint_eri_t* int12 ;
 double alphaP, alphaQ ;
 double P[3], Q[3] ;
 double U ;
 double F[am+1] ;
 double gammapq ;
 double gammapq_rc2 ;
 double gammapq_ratio ;
 double PQx ;
 double PQy ;
 double PQz ;
 double PQ2 ;
 double Wx ;
 double Wy ;
 double Wz ;
 double AB2 ;
 double CD2 ;
 double pfac ;
 constexpr double pi_2p5 = pow(M_PI,2.5) ;


 int icontrdepth4 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {
     for( int icontrdepthC=0; icontrdepthC < contrdepthC; icontrdepthC++)  {
       for( int icontrdepthD=0; icontrdepthD < contrdepthD; icontrdepthD++)  {


         int12 = &inteval[icontrdepth4] ;
         alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB];
         alphaQ = alphaC[icontrdepthC] + alphaD[icontrdepthD];
         P[0] = (alphaA[icontrdepthA] * A[0] + alphaB[icontrdepthB] * B[0]) / alphaP ;
         P[1] = (alphaA[icontrdepthA] * A[1] + alphaB[icontrdepthB] * B[1]) / alphaP ;
         P[2] = (alphaA[icontrdepthA] * A[2] + alphaB[icontrdepthB] * B[2]) / alphaP ;
         Q[0] = (alphaC[icontrdepthC] * C[0] + alphaD[icontrdepthD] * D[0]) / alphaQ ;
         Q[1] = (alphaC[icontrdepthC] * C[1] + alphaD[icontrdepthD] * D[1]) / alphaQ ;
         Q[2] = (alphaC[icontrdepthC] * C[2] + alphaD[icontrdepthD] * D[2]) / alphaQ ;

    
         int12->PA_x[0] = P[0] - A[0] ;
         int12->PA_y[0] = P[1] - A[1] ;
         int12->PA_z[0] = P[2] - A[2] ;
         int12->PB_x[0] = P[0] - B[0] ;
         int12->PB_y[0] = P[1] - B[1] ;
         int12->PB_z[0] = P[2] - B[2] ;
         int12->AB_x[0] = A[0] - B[0] ;
         int12->AB_y[0] = A[1] - B[1] ;
         int12->AB_z[0] = A[2] - B[2] ;
         int12->oo2z[0] = 0.5 / alphaP ;

         AB2 = int12->AB_x[0] * int12->AB_x[0] 
             + int12->AB_y[0] * int12->AB_y[0] 
             + int12->AB_z[0] * int12->AB_z[0] ;
    
         gammapq = alphaP * alphaQ / (alphaP + alphaQ);
    
         /*
           Treat the LR integrals by simply modifying gammapq into gammapq_rc2
           And introducing gammapq_ratio = gammapq_rc2 / gammapq
          */
    
         gammapq_rc2   =   alphaP * alphaQ   / ( alphaP + alphaQ + alphaP * alphaQ * rcut * rcut );
         gammapq_ratio = ( alphaP + alphaQ ) / ( alphaP + alphaQ + alphaP * alphaQ * rcut * rcut );
    
    
         int12->QC_x[0] = Q[0] - C[0] ;
         int12->QC_y[0] = Q[1] - C[1] ;
         int12->QC_z[0] = Q[2] - C[2] ;
         int12->CD_x[0] = C[0] - D[0] ;
         int12->CD_y[0] = C[1] - D[1] ;
         int12->CD_z[0] = C[2] - D[2] ;
         CD2 = int12->CD_x[0] * int12->CD_x[0] 
             + int12->CD_y[0] * int12->CD_y[0] 
             + int12->CD_z[0] * int12->CD_z[0] ;
         int12->oo2e[0] = 0.5 / alphaQ;
    

         PQx = P[0] - Q[0];
         PQy = P[1] - Q[1];
         PQz = P[2] - Q[2];
         PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
         Wx = (alphaP * P[0] + alphaQ * Q[0]) / (alphaP + alphaQ);
         Wy = (alphaP * P[1] + alphaQ * Q[1]) / (alphaP + alphaQ);
         Wz = (alphaP * P[2] + alphaQ * Q[2]) / (alphaP + alphaQ);
  
         int12->WP_x[0] = Wx - P[0];
         int12->WP_y[0] = Wy - P[1];
         int12->WP_z[0] = Wz - P[2];
         int12->WQ_x[0] = Wx - Q[0];
         int12->WQ_y[0] = Wy - Q[1];
         int12->WQ_z[0] = Wz - Q[2];
         int12->oo2ze[0] = 0.5 / ( alphaP + alphaQ ) ;
#if LIBINT2_DEFINED(eri,roz)
         int12->roz[0] = gammapq/alphaP;
#endif
#if LIBINT2_DEFINED(eri,roe)
         int12->roe[0] = gammapq/alphaQ;
#endif


         pfac = 2 * pi_2p5 / (alphaP * alphaQ * sqrt(alphaP + alphaQ))
                   * exp(-alphaA[icontrdepthA] * alphaB[icontrdepthB] * AB2 / alphaP) 
                   * exp(-alphaC[icontrdepthC] * alphaD[icontrdepthD] * CD2 / alphaQ) 
                    * cA[icontrdepthA] * cB[icontrdepthB] * cC[icontrdepthC] * cD[icontrdepthD] ;
         U = PQ2 * gammapq_rc2 ;
     
         boys_function_c(F, am, U);

  // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
         int12->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0] * pow(gammapq_ratio,0.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
         if( am > 0 ) int12->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1] * pow(gammapq_ratio,1.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
         if( am > 1 ) int12->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2] * pow(gammapq_ratio,2.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
         if( am > 2 ) int12->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3] * pow(gammapq_ratio,3.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
         if( am > 3 ) int12->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4] * pow(gammapq_ratio,4.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
         if( am > 4 ) int12->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5] * pow(gammapq_ratio,5.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
         if( am > 5 ) int12->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6] * pow(gammapq_ratio,6.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
         if( am > 6 ) int12->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7] * pow(gammapq_ratio,7.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
         if( am > 7 ) int12->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8] * pow(gammapq_ratio,8.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
         if( am > 8 ) int12->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9] * pow(gammapq_ratio,9.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
         if( am > 9 ) int12->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10] * pow(gammapq_ratio,10.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
         if( am > 10) int12->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11] * pow(gammapq_ratio,11.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
         if( am > 11) int12->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12] * pow(gammapq_ratio,12.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
         if( am > 12) int12->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13] * pow(gammapq_ratio,13.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
         if( am > 13) int12->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14] * pow(gammapq_ratio,14.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
         if( am > 14) int12->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15] * pow(gammapq_ratio,15.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
         if( am > 15) int12->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16] * pow(gammapq_ratio,16.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
         if( am > 16) int12->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17] * pow(gammapq_ratio,17.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
         if( am > 17) int12->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18] * pow(gammapq_ratio,18.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
         if( am > 18) int12->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19] * pow(gammapq_ratio,19.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
         if( am > 19) int12->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20] * pow(gammapq_ratio,20.5);
#endif


         int12->_0_Overlap_0_x[0] = 0.0 ;
         int12->_0_Overlap_0_y[0] = 0.0 ;
         int12->_0_Overlap_0_z[0] = 0.0 ;
   
         int12->veclen = 1 ;
         int12->contrdepth = contrdepth4 ;

         icontrdepth4++ ;
       }
     }
   }
 }



 if( am == 0 ) {

   eriABCD[0] = 0.0 ;
   for( int icontrdepth4=0; icontrdepth4 < contrdepth4; icontrdepth4++) {
     eriABCD[0] +=  inteval[icontrdepth4].LIBINT_T_SS_EREP_SS(0)[0] ;
   }

 } else {

   LIBINT2_PREFIXED_NAME(libint2_build_eri)[amA][amB][amC][amD](inteval);
   for( int i1234=0; i1234 < nint(amA) * nint(amB) * nint(amC) * nint(amD) ; ++i1234 ) {
     eriABCD[i1234] = inteval[0].targets[0][i1234] ;
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(inteval);

}
}


/* ========================================================================== */
