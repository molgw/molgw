/*
 * C wrapper to the libint library
 * includes contracted gaussian treatment
 */
#include<stdlib.h>
#include<stdio.h>
#include<libint2.h>
#include<math.h>
  
#define MAXFAC_BOYS  100
#define EPS_BOYS 1.0E-17     /* Absolute precision in computing Fm(t) */


/* First, the interfaces */
int nint(int am);
int max4(int am0, int am1, int am2, int am3);
void calc_boys(double*, int, double);

int libint_init();

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

int libint_init() {
  /* this initializes internal Libint data structures -- must happen once in the program */
  LIBINT2_PREFIXED_NAME(libint2_static_init)();
  return LIBINT2_MAX_AM_ERI;
}

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

 /* LIBINT2_MAX_AM_ERI is a macro defined in libint2.h that specifies the maximum angular momentum
    this Libint library instance can handle */
 if( ammax > LIBINT2_MAX_AM_ERI)
   return 1;

 LIBINT2_PREFIXED_NAME( libint2_init_eri)(&inteval[0], ammax, 0);

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


 /* this compute the shell set (quartet) of integrals */
 LIBINT2_PREFIXED_NAME(libint2_build_eri)[am0][am1][am2][am3](&inteval[0]);

 for(p0123 = 0 ; p0123< nint(am0)*nint(am1)*nint(am2)*nint(am3) ; ++p0123) {
   integrals[p0123] = inteval[0].targets[0][p0123];
 }

 // this releases all memory that was allocated for this object
 LIBINT2_PREFIXED_NAME( libint2_cleanup_eri)(&inteval[0]);

 return 0;

}




