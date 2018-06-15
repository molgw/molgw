/*
 * This file is part of MOLGW.
 * C/C++ wrapper to the libint library
 *
 * Header
 *
 * Author: F. Bruneval
 */


/* Headers */
extern "C" void boys_function_c(double*, int, double);

/* Number of cartesian function for a given angular momentum */
inline int nint(int am) {
  return (am+1)*(am+2)/2;
}

using namespace std;

#ifndef PI_2P5
#define PI_2P5
const double pi_2p5 = pow(M_PI,2.5) ;
#endif
