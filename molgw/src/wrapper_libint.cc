//////////////////////////////////////////////////////////////////////////////////
// This program shows how to use Libint for computing electron repulsion integrals
// (c) 2011 Edward F. Valeev
//////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <libint2.h>
#include <cstddef>


#define MAXFAC 100
#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */
#define MIN(a,b) ((a)>(b) ? (b) : (a))

/// this function computes all data required by Libint to evaluate the integrals efficiently
template<typename LibintEval>
void prep_libint2(LibintEval* erieval, unsigned int am1, double alpha1,
                  double A[3], unsigned int am2, double alpha2, double B[3],
                  unsigned int am3, double alpha3, double C[3],
                  unsigned int am4, double alpha4, double D[3], int norm_flag);

template<typename LibintEval>
void prep_libint2_lr(LibintEval* erieval, unsigned int am1, double alpha1,
                     double A[3], unsigned int am2, double alpha2, double B[3],
                     unsigned int am3, double alpha3, double C[3],
                     unsigned int am4, double alpha4, double D[3], int norm_flag,
                     double omega_range);


/*
 \param norm_flag:  tells what kind of normalization to use,
 0 - no normalization, >0 - normalized ERI
 */


extern "C"{
 int calculate_integral(double *,                                // omega_range
                        int*, int*, int*,int*,                   // am's
                        double*, double*, double*, double*,      // alpha's
                        double*, double*, double*,               // x1
                        double*, double*, double*,               // x2
                        double*, double*, double*,               // x3
                        double*, double*, double*,               // x4
                        double*);
 int calculate_integral_(double *,                                // omega_range
                         int*, int*, int*,int*,                   // am's
                         double*, double*, double*, double*,      // alpha's
                         double*, double*, double*,               // x1
                         double*, double*, double*,               // x2
                         double*, double*, double*,               // x3
                         double*, double*, double*,               // x4
                         double*);
}

int calculate_integral_(double* omega_range,
                        int* am1_in, int* am2_in, int* am3_in, int* am4_in,
                        double* alpha1, double* alpha2, double* alpha3, double* alpha4,
                        double* x1x, double* x1y, double* x1z,
                        double* x2x, double* x2y, double* x2z,
                        double* x3x, double* x3y, double* x3z,
                        double* x4x, double* x4y, double* x4z,
                        double* integrals) {
  return calculate_integral(omega_range,
                            am1_in,am2_in,am3_in,am4_in,alpha1,alpha2,alpha3,alpha4,
                            x1x, x1y, x1z,
                            x2x, x2y, x2z,
                            x3x, x3y, x3z,
                            x4x, x4y, x4z,
                            integrals);
}

int calculate_integral(double* omega_range,
                       int* am1_in, int* am2_in, int* am3_in, int* am4_in,
                       double* alpha1, double* alpha2, double* alpha3, double* alpha4,
                       double* x1x, double* x1y, double* x1z,
                       double* x2x, double* x2y, double* x2z,
                       double* x3x, double* x3y, double* x3z,
                       double* x4x, double* x4y, double* x4z,
                       double* integrals) {

  typedef unsigned int uint;

  // this initializes internal Libint data structures -- must happen once in the program
  LIBINT2_PREFIXED_NAME(libint2_static_init)();

  // This example assumes that your library does not support for vectorization
  const uint veclen = 1;

  // These parameters define 4 primitive Gaussian basis functions
  double alpha[4] ;                // orbital exponents for the 4 functions
  double A[3] = { 0.0, 0.0, 0.0 }; // position of function 1
  double B[3] = { 0.0, 0.0, 0.0 }; // etc.
  double C[3] = { 0.0, 0.0, 0.0 };
  double D[3] = { 0.0, 0.0, 0.0 };

  A[0] = *x1x;
  A[1] = *x1y;
  A[2] = *x1z;
  B[0] = *x2x;
  B[1] = *x2y;
  B[2] = *x2z;
  C[0] = *x3x;
  C[1] = *x3y;
  C[2] = *x3z;
  D[0] = *x4x;
  D[1] = *x4y;
  D[2] = *x4z;

  alpha[0] = *alpha1;
  alpha[1] = *alpha2;
  alpha[2] = *alpha3;
  alpha[3] = *alpha4;
  uint am0 = *am1_in;
  uint am1 = *am2_in;
  uint am2 = *am3_in;
  uint am3 = *am4_in;

  // LIBINT2_MAX_AM_ERI is a macro defined in libint2.h that specifies the maximum angular momentum
  // this Libint library instance can handle
  const unsigned int ammax =  LIBINT2_MAX_AM_ERI;
  // Libint_t is the type of a data structure used to pass basis function data to Libint
  Libint_t inteval;
  // Libint_t objects must be initialized prior to use
  // your program may have several such objects to implement computation of integrals in multiple threads
  // each thread would have its own Libint_t object
  LIBINT2_PREFIXED_NAME( libint2_init_eri)(&inteval, ammax, 0);

  // Libint only provides code to computes a symmetry-unique subset of integrals
  // the default convention is to compute integrals (0 1|2 3) with am0 >= am1, am2 >= am3, and am2+am3 >= am0+am1
  // in general you will need to permute shells to satisfy this ordering
  // most production codes do this by precomputing "significant" shell pair lists
  bool can_compute = (am0 >= am1) && (am2 >= am3) &&
                     (am0 + am1 <= am2 + am3);
  can_compute &= (am0 != 0 || am1 != 0 || am2 != 0 || am3 != 0); // skip (ss|ss) integral -- no need to use Libint for that one
  if (can_compute == false)
    return 1;

  // this function fills in all data expected by Libint to evaluate the given integral set
  if ( *omega_range > 0.999999e6 ) {
    prep_libint2(&inteval, am0, alpha[0], A, am1, alpha[1], B, am2,
                 alpha[2], C, am3, alpha[3], D, 0);
  } 
  else {
    prep_libint2_lr(&inteval, am0, alpha[0], A, am1, alpha[1], B, am2,
                    alpha[2], C, am3, alpha[3], D, 0, *omega_range );
  } 
  


  // these are uncontracted basis functions
  // N.B. to compute integrals over contracted functions allocate an array of Libint_t objects (one for each primitive combination)
  //      and fill each with primitive combination data. You only need to call libint2_init_eri using pointer to the first object!
  inteval.contrdepth = 1;
  // this compute the shell set (quartet) of integrals
  LIBINT2_PREFIXED_NAME
      ( libint2_build_eri)[am0][am1][am2][am3](&inteval);

  bool success = true;
  int ijkl = 0;
  // this doubly-nested loop iterates over functions in Libint standard ordering; for example, for f functions is:
  // xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
  // Libint compiler can produce libraries that support different orderings
  for (uint k0 = 0; k0 <= am0; ++k0) {
    for (uint n0 = 0; n0 <= k0; ++n0) {
      const uint m0 = k0 - n0;
      const uint l0 = am0 - k0;
      // l0, m0, n0 are x,y,z exponents for the first basis function

      for (uint k1 = 0; k1 <= am1; ++k1) {
        for (uint n1 = 0; n1 <= k1; ++n1) {
          const uint m1 = k1 - n1;
          const uint l1 = am1 - k1;

          for (uint k2 = 0; k2 <= am2; ++k2) {
            for (uint n2 = 0; n2 <= k2; ++n2) {
              const uint m2 = k2 - n2;
              const uint l2 = am2 - k2;

              for (uint k3 = 0; k3 <= am3; ++k3) {
                for (uint n3 = 0; n3 <= k3; ++n3) {
                  const uint m3 = k3 - n3;
                  const uint l3 = am3 - k3;


                  integrals[ijkl]= inteval.targets[0][ijkl];

                  ++ijkl;

                }
              }
            }
          }
        }
      }
    }
  } // end of loop over basis functions in the shell quartet



  // this releases all memory that was allocated for this object
  LIBINT2_PREFIXED_NAME( libint2_cleanup_eri)(&inteval);

  return 0;
}

static double *df;

extern "C"{
void calc_f(double*, int, double);
void calc_f_(double*, int*, double*);
}

double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                  double alpha1, const double* A);
double* init_array(unsigned long int size);
double** block_matrix(unsigned long int nrow, unsigned long int ncol);
void free_array(double* array);


/*!
 calc_f()

 This function computes infamous integral Fn(t). For its definition
 see Obara and Saika paper, or Shavitt's chapter in the
 Methods in Computational Physics book (see reference below).
 This piece of code is from Dr. Justin Fermann's program CINTS

 \ingroup (QT)
 */

void calc_f_(double *F, int *n, double *t) {
  calc_f(F,*n,*t);
}

void calc_f(double *F, int n, double t) {
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  static double K = 1.0 / M_2_SQRTPI;
  double et;

  if (df == NULL) {
    df = init_array(2 * MAXFAC);
    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (i = 3; i < MAXFAC * 2; i++) {
      df[i] = (i - 1) * df[i - 2];
    }
  }

  if (t > 20.0) { /* For big t's do upward recursion */
    t2 = 2 * t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K * erf(t) / t;
    for (m = 0; m <= n - 1; m++) {
      F[m + 1] = ((2 * m + 1) * F[m] - et) / (t2);
    }
  } else { /* For smaller t's compute F with highest n using
   asymptotic series (see I. Shavitt in
   Methods in Computational Physics, ed. B. Alder eta l,
   vol 2, 1963, page 8) */

    et = exp(-t);
    t2 = 2 * t;
    m2 = 2 * n;
    num = df[m2];
    i = 0;
    sum = 1.0 / (m2 + 1);
    do {
      i++;
      num = num * t2;
      term1 = num / df[m2 + 2 * i + 2];
      sum += term1;
    } while (fabs(term1) > EPS && i < MAXFAC);
    F[n] = sum * et;
    for (m = n - 1; m >= 0; m--) { /* And then do downward recursion */
      F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
    }
  }
/* std::cout << "F[n] " <<  F[n] << std::endl; */
}

/*!
 norm_const()

 \ingroup (QT)
 */
double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                  double alpha1, const double* A) {
  return pow(2 * alpha1 / M_PI, 0.75) * pow(4 * alpha1, 0.5 * (l1 + m1 + n1))
      / sqrt(df[2 * l1] * df[2 * m1] * df[2 * n1]);
}

double* init_array(unsigned long int size) {
  double* result = new double[size];
  for (int i = 0; i < size; i++)
    result[i] = 0.0;
  return result;
}

double** block_matrix(unsigned long int nrow, unsigned long int ncol) {
  double** rows = new double*[nrow];
  rows[0] = new double[nrow * ncol];
  for (int i = 1; i < nrow; i++)
    rows[i] = rows[i - 1] + ncol;

  return rows;
}

void free_array(double* array) {
  delete[] array;
}

template<typename LibintEval>
void prep_libint2(LibintEval* erieval, unsigned int am1, double alpha1,
                  double A[3], unsigned int am2, double alpha2, double B[3],
                  unsigned int am3, double alpha3, double C[3],
                  unsigned int am4, double alpha4, double D[3], int norm_flag) {

  const unsigned int am = am1 + am2 + am3 + am4;
  double* F = init_array(am + 1);

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

  if (norm_flag > 0) {
    /*    pfac *= norm_const(l1,m1,n1,alpha1,A);
     pfac *= norm_const(l2,m2,n2,alpha2,B);
     pfac *= norm_const(l3,m3,n3,alpha3,C);
     pfac *= norm_const(l4,m4,n4,alpha4,D);*/
  }

  calc_f(F, am, PQ2 * gammapq);

  // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
  erieval->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
  erieval->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
  erieval->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
  erieval->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
  erieval->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
  erieval->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
  erieval->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
  erieval->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
  erieval->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
  erieval->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
  erieval->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
  erieval->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
  erieval->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
  erieval->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
  erieval->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
  erieval->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20];
#endif

  free_array(F);
}

template<typename LibintEval>
void prep_libint2_lr(LibintEval* erieval, unsigned int am1, double alpha1,
                  double A[3], unsigned int am2, double alpha2, double B[3],
                  unsigned int am3, double alpha3, double C[3],
                  unsigned int am4, double alpha4, double D[3], int norm_flag,
                  double omega_range) {

  const unsigned int am = am1 + am2 + am3 + am4;
  double* F = init_array(am + 1);

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
/*
  Here is the usual bare Coulomb interaction factor
*/

  const double gammapq = gammap * gammaq / (gammap + gammaq);

/*
  Then comes the long-range only Coulomb interaction factor
*/
  const double omega2 = pow(omega_range,2);
  const double gammapq_omega2 = gammap * gammaq * omega2 / ( gammap*omega2 + gammaq*omega2 + gammap*gammaq );

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
  double pfac = 2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq
      * sqrt(gammap + gammaq));

  if (norm_flag > 0) {
    /*    pfac *= norm_const(l1,m1,n1,alpha1,A);
     pfac *= norm_const(l2,m2,n2,alpha2,B);
     pfac *= norm_const(l3,m3,n3,alpha3,C);
     pfac *= norm_const(l4,m4,n4,alpha4,D);*/
  }

  calc_f(F, am, PQ2 * gammapq_omega2);

  // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0] * pow(gammapq_omega2/gammapq,0.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1] * pow(gammapq_omega2/gammapq,1.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2] * pow(gammapq_omega2/gammapq,2.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3] * pow(gammapq_omega2/gammapq,3.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4] * pow(gammapq_omega2/gammapq,4.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
  erieval->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5] * pow(gammapq_omega2/gammapq,5.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
  erieval->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6] * pow(gammapq_omega2/gammapq,6.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
  erieval->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7] * pow(gammapq_omega2/gammapq,7.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
  erieval->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8] * pow(gammapq_omega2/gammapq,8.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
  erieval->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9] * pow(gammapq_omega2/gammapq,9.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
  erieval->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10] * pow(gammapq_omega2/gammapq,10.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
  erieval->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11] * pow(gammapq_omega2/gammapq,11.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
  erieval->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12] * pow(gammapq_omega2/gammapq,12.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
  erieval->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13] * pow(gammapq_omega2/gammapq,13.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
  erieval->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14] * pow(gammapq_omega2/gammapq,14.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
  erieval->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15] * pow(gammapq_omega2/gammapq,15.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
  erieval->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16] * pow(gammapq_omega2/gammapq,16.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
  erieval->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17] * pow(gammapq_omega2/gammapq,17.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
  erieval->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18] * pow(gammapq_omega2/gammapq,18.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
  erieval->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19] * pow(gammapq_omega2/gammapq,19.5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
  erieval->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20] * pow(gammapq_omega2/gammapq,20.5);
#endif

  free_array(F);
}
