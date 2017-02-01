/*
 * This file is part of MOLGW.
 * C wrapper to the libint library
 * includes contracted gaussian treatment
 */
#include<stdlib.h>
#include<stdio.h>
#include<libint2.h>
#include<math.h>
  

/* First, the interfaces */
int nint(int am);
int max4(int am0, int am1, int am2, int am3);
extern "C" void boys_function_c(double*, int, double);

extern "C"
int eval_contr_integral(
                        int *am0_in, int *am1_in, int *am2_in, int *am3_in,
                        int *ncontr0, int *ncontr1, int *ncontr2, int *ncontr3,
                        double *contr0, double *contr1, double *contr2, double *contr3,
                        double *alpha0, double *alpha1, double *alpha2, double *alpha3,
                        double *x0, double *x1, double *x2, double *x3,
                        double* rcut,          /* rcut is omega^{-1} */
                        double *integrals
                       );

void prep_libint2_contr(Libint_t *erieval,
                        unsigned int am1, double alpha1, double A[3],
                        unsigned int am2, double alpha2, double B[3],
                        unsigned int am3, double alpha3, double C[3],
                        unsigned int am4, double alpha4, double D[3], 
                        double  rcut,          /* rcut is omega^{-1} */
                        int norm_flag,
                        double norm_prefactor
                       );


/* Then the real coding */

extern "C" {
void libint_init(int *ammax, bool *has_onebody, bool *has_gradient) {
 /* this initializes internal Libint data structures -- must happen once in the program */
 LIBINT2_PREFIXED_NAME(libint2_static_init)();
 *has_onebody = LIBINT2_SUPPORT_ONEBODY ;
 *has_gradient = ( LIBINT2_DERIV_ERI_ORDER > 0 );
 *ammax = LIBINT2_MAX_AM_eri ;
}
}

extern "C" {
int eval_contr_integral(
                        int *am0_in, int *am1_in, int *am2_in, int *am3_in,
                        int *ncontr0, int *ncontr1, int *ncontr2, int *ncontr3,
                        double *contr0, double *contr1, double *contr2, double *contr3,
                        double *alpha0, double *alpha1, double *alpha2, double *alpha3,
                        double *x0, double *x1, double *x2, double *x3,
                        double* rcut,          /* rcut is omega^{-1} */
                        double *integrals
                       ) {

 unsigned int am0 = *am0_in;
 unsigned int am1 = *am1_in;
 unsigned int am2 = *am2_in;
 unsigned int am3 = *am3_in;
 unsigned int ncontr0123=*ncontr0 * *ncontr1 * *ncontr2 * *ncontr3;
 Libint_t inteval[ncontr0123];
 int p0, p1, p2, p3, p0123;
 double contrcoef0123;
 int k0,k1,k2,k3;
 int n0,n1,n2,n3;
 unsigned int ammax;


 ammax = max4(am0,am1,am2,am3);

 /* LIBINT2_MAX_AM_eri is a macro defined in libint2.h that specifies the maximum angular momentum
    this Libint library instance can handle */
 if( ammax > LIBINT2_MAX_AM_eri)
   return 1;

 LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[0], ammax, 0);

/* 
 * Prepare the calculation for each contraction
 */

 p0123=0;
 for(p0=0; p0<*ncontr0; ++p0) {
   for(p1=0; p1<*ncontr1; ++p1) {
     for(p2=0; p2<*ncontr2; ++p2) {
       for(p3=0; p3<*ncontr3; ++p3) {
         contrcoef0123 = contr0[p0] * contr1[p1] * contr2[p2] * contr3[p3];

         prep_libint2_contr(&inteval[p0123], /* data for each primitive combination goes into its own evaluator object */
                      am0, alpha0[p0], x0,
                      am1, alpha1[p1], x1,
                      am2, alpha2[p2], x2,
                      am3, alpha3[p3], x3,
                      *rcut,                /* rcut is omega^{-1} */
                      0, contrcoef0123);
         p0123++;
       }
     }
   }
 }
 inteval[0].contrdepth = ncontr0123;


 if( ammax == 0 ) {
   integrals[0] = 0.0 ;
   for(p0123=0; p0123<ncontr0123; ++p0123) 
     integrals[0] += inteval[p0123].LIBINT_T_SS_EREP_SS(0)[0] ;
   
 } else {
   /* this compute the shell set (quartet) of integrals */
   LIBINT2_PREFIXED_NAME(libint2_build_eri)[am0][am1][am2][am3](&inteval[0]);

   for(p0123 = 0 ; p0123< nint(am0)*nint(am1)*nint(am2)*nint(am3) ; ++p0123) {
     integrals[p0123] = inteval[0].targets[0][p0123];
   }
 }

 // this releases all memory that was allocated for this object
 LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);

 return 0;

}
}

/* Number of cartesian function for a given angular momentum */
int nint(int am) {
  return (am+1)*(am+2)/2;
}

/* Maximum of 4 integers */
int max4(int am0, int am1, int am2, int am3) {
 int ammax;

 if(am0 > am1) {
   ammax=am0;
 } else {
   ammax=am1;
 }
 if(am2 > ammax) {
   ammax=am2;
 } 
 if(am3 > ammax) {
   ammax=am3;
 } 

 return ammax;
}



void prep_libint2_contr(Libint_t* erieval,
                        unsigned int am1, double alpha1, double A[3],
                        unsigned int am2, double alpha2, double B[3],
                        unsigned int am3, double alpha3, double C[3],
                        unsigned int am4, double alpha4, double D[3],
                        double rcut,
                        int norm_flag, double norm_prefactor) {

  const unsigned int am = am1 + am2 + am3 + am4;
  double  F[am + 1];
  int i;

  const double gammap = alpha1 + alpha2;
  const double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
  const double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
  const double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
  const double PAx = Px - A[0];
  const double PAy = Py - A[1];
  const double PAz = Pz - A[2];
  const double PBx = Px - B[0];
  const double PBy = Py - B[1];
  const double PBz = Pz - B[2];
  const double AB_x = A[0] - B[0];
  const double AB_y = A[1] - B[1];
  const double AB_z = A[2] - B[2];
  const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

  for(i=0;i<am+1;++i) F[i] = 0. ;

#if LIBINT2_DEFINED(eri,PA_x)
  erieval->PA_x[0] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
  erieval->PA_y[0] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
  erieval->PA_z[0] = PAz;
#endif
#if LIBINT2_DEFINED(eri,PB_x)
  erieval->PB_x[0] = PBx;
#endif
#if LIBINT2_DEFINED(eri,PB_y)
  erieval->PB_y[0] = PBy;
#endif
#if LIBINT2_DEFINED(eri,PB_z)
  erieval->PB_z[0] = PBz;
#endif

#if LIBINT2_DEFINED(eri,AB_x)
  erieval->AB_x[0] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
  erieval->AB_y[0] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
  erieval->AB_z[0] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
  erieval->oo2z[0] = 0.5/gammap;
#endif

  const double gammaq = alpha3 + alpha4;
  const double gammapq = gammap * gammaq / (gammap + gammaq);
  const double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
  const double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
  const double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
  const double QCx = Qx - C[0];
  const double QCy = Qy - C[1];
  const double QCz = Qz - C[2];
  const double QDx = Qx - D[0];
  const double QDy = Qy - D[1];
  const double QDz = Qz - D[2];
  const double CD_x = C[0] - D[0];
  const double CD_y = C[1] - D[1];
  const double CD_z = C[2] - D[2];
  const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;
/*
  Treat the LR integrals by simply modifying gammapq into gammapq_rc2
  And introducing gammapq_ratio = gammapq_rc2 / gammapq
 */
  const double gammapq_rc2 = gammap * gammaq / (gammap + gammaq + gammap * gammaq * pow(rcut,2) );
  const double gammapq_ratio = ( gammap + gammaq ) / ( gammap + gammaq + gammap * gammaq * pow(rcut,2) );


#if LIBINT2_DEFINED(eri,QC_x)
  erieval->QC_x[0] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
  erieval->QC_y[0] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
  erieval->QC_z[0] = QCz;
#endif
#if LIBINT2_DEFINED(eri,QD_x)
  erieval->QD_x[0] = QDx;
#endif
#if LIBINT2_DEFINED(eri,QD_y)
  erieval->QD_y[0] = QDy;
#endif
#if LIBINT2_DEFINED(eri,QD_z)
  erieval->QD_z[0] = QDz;
#endif

#if LIBINT2_DEFINED(eri,CD_x)
  erieval->CD_x[0] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
  erieval->CD_y[0] = CD_y;
#endif
#if LIBINT2_DEFINED(eri,CD_z)
  erieval->CD_z[0] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
  erieval->oo2e[0] = 0.5/gammaq;
#endif

  // Prefactors for interelecttron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_x)
  erieval->TwoPRepITR_pfac0_0_x[0] = - (alpha2*AB_x + alpha4*CD_x)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_y)
  erieval->TwoPRepITR_pfac0_0_y[0] = - (alpha2*AB_y + alpha4*CD_y)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_z)
  erieval->TwoPRepITR_pfac0_0_z[0] = - (alpha2*AB_z + alpha4*CD_z)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_x)
  erieval->TwoPRepITR_pfac0_1_x[0] = - (alpha2*AB_x + alpha4*CD_x)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_y)
  erieval->TwoPRepITR_pfac0_1_y[0] = - (alpha2*AB_y + alpha4*CD_y)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_z)
  erieval->TwoPRepITR_pfac0_1_z[0] = - (alpha2*AB_z + alpha4*CD_z)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_0)
  erieval->TwoPRepITR_pfac1_0[0] = -gammaq/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_1)
  erieval->TwoPRepITR_pfac1_1[0] = -gammap/gammaq;
#endif

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
  const double Wx = (gammap * Px + gammaq * Qx) / (gammap + gammaq);
  const double Wy = (gammap * Py + gammaq * Qy) / (gammap + gammaq);
  const double Wz = (gammap * Pz + gammaq * Qz) / (gammap + gammaq);

#if LIBINT2_DEFINED(eri,WP_x)
  erieval->WP_x[0] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
  erieval->WP_y[0] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
  erieval->WP_z[0] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
  erieval->WQ_x[0] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
  erieval->WQ_y[0] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
  erieval->WQ_z[0] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
  erieval->oo2ze[0] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
  erieval->roz[0] = gammapq/gammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
  erieval->roe[0] = gammapq/gammaq;
#endif

  double K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
  double K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
  double pfac = 2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq * sqrt(gammap + gammaq));
  pfac *= norm_prefactor;

/*
  if (norm_flag > 0) {
     pfac *= norm_const(l1,m1,n1,alpha1,A);
     pfac *= norm_const(l2,m2,n2,alpha2,B);
     pfac *= norm_const(l3,m3,n3,alpha3,C);
     pfac *= norm_const(l4,m4,n4,alpha4,D);
  }
*/

  boys_function_c(F, am, PQ2 * gammapq_rc2);

  // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0] * pow(gammapq_ratio,0.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1] * pow(gammapq_ratio,1.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2] * pow(gammapq_ratio,2.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3] * pow(gammapq_ratio,3.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4] * pow(gammapq_ratio,4.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
  erieval->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5] * pow(gammapq_ratio,5.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
  erieval->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6] * pow(gammapq_ratio,6.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
  erieval->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7] * pow(gammapq_ratio,7.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
  erieval->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8] * pow(gammapq_ratio,8.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
  erieval->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9] * pow(gammapq_ratio,9.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
  erieval->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10] * pow(gammapq_ratio,10.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
  erieval->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11] * pow(gammapq_ratio,11.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
  erieval->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12] * pow(gammapq_ratio,12.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
  erieval->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13] * pow(gammapq_ratio,13.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
  erieval->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14] * pow(gammapq_ratio,14.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
  erieval->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15] * pow(gammapq_ratio,15.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
  erieval->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16] * pow(gammapq_ratio,16.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
  erieval->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17] * pow(gammapq_ratio,17.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
  erieval->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18] * pow(gammapq_ratio,18.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
  erieval->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19] * pow(gammapq_ratio,19.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
  erieval->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20] * pow(gammapq_ratio,20.5);
#endif

}

