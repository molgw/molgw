/*
 * This file is part of MOLGW.
 * C/C++ wrapper to the libint library
 * two-body gradient terms: 2- 3- and 4- center Coulomb integrals
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


/* ==========================================================================                    
 *                           2-center integrals
 * ========================================================================== */

extern "C" {
void libint_2center_grad(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                         int amC, int contrdepthC , double C [] , double alphaC [], double cC [],
                         double rcut,
                         double eriAC [] ) {

 const unsigned int contrdepth2 = contrdepthA * contrdepthC;
 Libint_2eri_t* inteval = libint2::malloc<Libint_2eri_t>(contrdepth2);

 assert( false && "libint_2center_grad not implemented yet" ) ;
 assert(amA <= LIBINT2_MAX_AM_2eri1);
 assert(amC <= LIBINT2_MAX_AM_2eri1);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthC == 1 );
#endif

 const unsigned int ammax = max(amA,amC) ;
 const int am = amA + amC ;


 LIBINT2_PREFIXED_NAME(libint2_init_2eri)(inteval, ammax, 0);

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


 int icontrdepth2 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthC=0; icontrdepthC < contrdepthC; icontrdepthC++)  {

     Libint_2eri_t* int12 = &inteval[icontrdepth2] ;
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

     pfac *= sqrt(gammapq_ratio);
     double* erep = int12->LIBINT_T_SS_EREP_SS(0);
     for( int l=0; l <= am+1 ; ++l , ++erep ) {
       *erep = pfac * F[l] ;
       pfac *= gammapq_ratio ;
     }


     int12->_0_Overlap_0_x[0] = 0.0 ;
     int12->_0_Overlap_0_y[0] = 0.0 ;
     int12->_0_Overlap_0_z[0] = 0.0 ;
 
     int12->veclen = 1 ;
     int12->contrdepth = contrdepth2 ;

     icontrdepth2++ ;
   }
 }




 LIBINT2_PREFIXED_NAME(libint2_build_2eri)[amA][amC](inteval);
 for( int i12=0; i12 < nint(amA) * nint(amC) ; ++i12 ) {
   eriAC[i12] = inteval[0].targets[0][i12] ;
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri)(inteval);
 free(inteval);

}
}



/* ==========================================================================                    
 *                           3-center integrals
 * ========================================================================== */


extern "C" {
void libint_3center_grad(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                         int amC, int contrdepthC , double C [] , double alphaC [], double cC [], 
                         int amD, int contrdepthD , double D [] , double alphaD [], double cD [],
                         double rcut,
                         double eriACD [] ) {

 const unsigned int contrdepth3 = contrdepthA * contrdepthC * contrdepthD ;
 Libint_3eri_t* inteval = libint2::malloc<Libint_3eri_t>(contrdepth3);

 assert( false && "libint_3center_grad not implemented yet" ) ;

 assert(amA <= LIBINT2_MAX_AM_3eri1);
 assert(amC <= LIBINT2_MAX_AM_3eri1);
 assert(amD <= LIBINT2_MAX_AM_3eri1);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthC == 1 );
 assert( contrdepthD == 1 );
#endif

 const unsigned int ammax = max(max(amA,amC),amD) ;
 const int am = amA + amC + amD ;


 LIBINT2_PREFIXED_NAME(libint2_init_3eri)(inteval, ammax, 0);

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


 int icontrdepth3 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthC=0; icontrdepthC < contrdepthC; icontrdepthC++)  {
     for( int icontrdepthD=0; icontrdepthD < contrdepthD; icontrdepthD++)  {


       Libint_3eri_t* int12 = &inteval[icontrdepth3] ;
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

       pfac *= sqrt(gammapq_ratio);
       double* erep = int12->LIBINT_T_SS_EREP_SS(0);
       for( int l=0; l <= am+1 ; ++l , ++erep ) {
         *erep = pfac * F[l] ;
         pfac *= gammapq_ratio ;
       }



       int12->_0_Overlap_0_x[0] = 0.0 ;
       int12->_0_Overlap_0_y[0] = 0.0 ;
       int12->_0_Overlap_0_z[0] = 0.0 ;
 
       int12->veclen = 1 ;
       int12->contrdepth = contrdepth3 ;

       icontrdepth3++ ;
     }
   }
 }



 LIBINT2_PREFIXED_NAME(libint2_build_3eri)[amA][amC][amD](inteval);
 for( int i123=0; i123 < nint(amA) * nint(amC) * nint(amD) ; ++i123 ) {
   eriACD[i123] = inteval[0].targets[0][i123] ;
 }



 LIBINT2_PREFIXED_NAME(libint2_cleanup_3eri)(inteval);
 free(inteval);

}
}


/* ==========================================================================                    
 *                           4-center integrals
 * ========================================================================== */


extern "C" {
void libint_4center_grad(int amA, int contrdepthA , double A [] , double alphaA [], double cA [], 
                         int amB, int contrdepthB , double B [] , double alphaB [], double cB [], 
                         int amC, int contrdepthC , double C [] , double alphaC [], double cC [], 
                         int amD, int contrdepthD , double D [] , double alphaD [], double cD [],
                         double rcut,
                         double eriAx [], double eriAy [], double eriAz [], 
                         double eriBx [], double eriBy [], double eriBz [], 
                         double eriCx [], double eriCy [], double eriCz [], 
                         double eriDx [], double eriDy [], double eriDz [] ) {

 const unsigned int contrdepth4 = contrdepthA * contrdepthB * contrdepthC * contrdepthD ;
 Libint_eri1_t* inteval = libint2::malloc<Libint_eri1_t>(contrdepth4);

 assert(amA <= LIBINT2_MAX_AM_eri1);
 assert(amB <= LIBINT2_MAX_AM_eri1);
 assert(amC <= LIBINT2_MAX_AM_eri1);
 assert(amD <= LIBINT2_MAX_AM_eri1);

#ifndef LIBINT2_CONTRACTED_INTS
 assert( contrdepthA == 1 );
 assert( contrdepthB == 1 );
 assert( contrdepthC == 1 );
 assert( contrdepthD == 1 );
#endif

 const unsigned int ammax = max(max(max(amA,amB),amC),amD) ;
 const int am = amA + amB + amC + amD ;
 const int ni = nint(amA) * nint(amB) * nint(amC) * nint(amD) ;


 LIBINT2_PREFIXED_NAME(libint2_init_eri1)(inteval, ammax, 0);

 double alphaP, alphaQ ;
 double rhoP, rhoQ ;
 double P[3], Q[3] ;
 double U ;
 double F[am+2] ;
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


 int icontrdepth4 = 0 ;
 for( int icontrdepthA=0; icontrdepthA < contrdepthA; icontrdepthA++)  {
   for( int icontrdepthB=0; icontrdepthB < contrdepthB; icontrdepthB++)  {
     for( int icontrdepthC=0; icontrdepthC < contrdepthC; icontrdepthC++)  {
       for( int icontrdepthD=0; icontrdepthD < contrdepthD; icontrdepthD++)  {


         Libint_eri1_t* int12 = &inteval[icontrdepth4] ;
         alphaP = alphaA[icontrdepthA] + alphaB[icontrdepthB] ;
         rhoP   = alphaA[icontrdepthA] * alphaB[icontrdepthB] / alphaP ;
         alphaQ = alphaC[icontrdepthC] + alphaD[icontrdepthD];
         rhoQ   = alphaC[icontrdepthC] * alphaD[icontrdepthD] / alphaQ ;
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

/* For derivatives only */
#if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
            int12->alpha1_rho_over_zeta2[0] = alphaA[icontrdepthA] * gammapq / (alphaP * alphaP);
#endif
#if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
            int12->alpha2_rho_over_zeta2[0] = alphaB[icontrdepthB] * gammapq / (alphaP * alphaP);
#endif
#if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
            int12->alpha3_rho_over_eta2[0] = alphaC[icontrdepthC] * gammapq / (alphaQ * alphaQ);
#endif
#if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
            int12->alpha4_rho_over_eta2[0] = alphaD[icontrdepthD] * gammapq / (alphaQ * alphaQ);
#endif
#if LIBINT2_DEFINED(eri,alpha1_over_zetapluseta)
            int12->alpha1_over_zetapluseta[0] = alphaA[icontrdepthA] / (alphaP + alphaQ);
#endif
#if LIBINT2_DEFINED(eri,alpha2_over_zetapluseta)
            int12->alpha2_over_zetapluseta[0] = alphaB[icontrdepthB] / (alphaP + alphaQ);
#endif
#if LIBINT2_DEFINED(eri,alpha3_over_zetapluseta)
            int12->alpha3_over_zetapluseta[0] = alphaC[icontrdepthC] / (alphaP + alphaQ);
#endif
#if LIBINT2_DEFINED(eri,alpha4_over_zetapluseta)
            int12->alpha4_over_zetapluseta[0] = alphaD[icontrdepthD] / (alphaP + alphaQ);
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
            int12->rho12_over_alpha1[0] = rhoP / alphaA[icontrdepthA];
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
            int12->rho12_over_alpha2[0] = rhoP / alphaB[icontrdepthB];
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha3)
            int12->rho34_over_alpha3[0] = rhoQ / alphaC[icontrdepthC];
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha4)
            int12->rho34_over_alpha4[0] = rhoQ / alphaD[icontrdepthD];
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
            int12->two_alpha0_bra[0] = 2.0 * alphaA[icontrdepthA];
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
            int12->two_alpha0_ket[0] = 2.0 * alphaB[icontrdepthB];
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_bra)
            int12->two_alpha1_bra[0] = 2.0 * alphaC[icontrdepthC];
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_ket)
            int12->two_alpha1_ket[0] = 2.0 * alphaD[icontrdepthD];;
#endif


         pfac = 2 * pi_2p5 / (alphaP * alphaQ * sqrt(alphaP + alphaQ))
                   * exp(-alphaA[icontrdepthA] * alphaB[icontrdepthB] * AB2 / alphaP) 
                   * exp(-alphaC[icontrdepthC] * alphaD[icontrdepthD] * CD2 / alphaQ) 
                    * cA[icontrdepthA] * cB[icontrdepthB] * cC[icontrdepthC] * cD[icontrdepthD] ;
         U = PQ2 * gammapq_rc2 ;
     
         boys_function_c(F, am+1, U);

         pfac *= sqrt(gammapq_ratio);
         double* erep = int12->LIBINT_T_SS_EREP_SS(0);
         for( int l=0; l <= am+1 ; ++l , ++erep ) {
           *erep = pfac * F[l] ;
           pfac *= gammapq_ratio ;
         }

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


 LIBINT2_PREFIXED_NAME(libint2_build_eri1)[amA][amB][amC][amD](inteval);


 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriAx[i1234] = inteval[0].targets[0][i1234 + 0*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriAy[i1234] = inteval[0].targets[0][i1234 + 1*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriAz[i1234] = inteval[0].targets[0][i1234 + 2*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriBx[i1234] = inteval[0].targets[0][i1234 + 3*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriBy[i1234] = inteval[0].targets[0][i1234 + 4*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriBz[i1234] = inteval[0].targets[0][i1234 + 5*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriDx[i1234] = inteval[0].targets[0][i1234 + 6*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriDy[i1234] = inteval[0].targets[0][i1234 + 7*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriDz[i1234] = inteval[0].targets[0][i1234 + 8*ni] ;

/*
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriCx[i1234] = inteval[0].targets[0][i1234 + 9*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriCy[i1234] = inteval[0].targets[0][i1234 +10*ni] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriCz[i1234] = inteval[0].targets[0][i1234 +11*ni] ;
*/
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriCx[i1234] = -eriAx[i1234] -eriBx[i1234] -eriDx[i1234] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriCy[i1234] = -eriAy[i1234] -eriBy[i1234] -eriDy[i1234] ;
 for( int i1234 = 0; i1234 < ni ; ++i1234 ) eriCz[i1234] = -eriAz[i1234] -eriBz[i1234] -eriDz[i1234] ;



 LIBINT2_PREFIXED_NAME(libint2_cleanup_eri1)(inteval);
 free(inteval);

}
}


/* ========================================================================== */
