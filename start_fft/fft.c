/*
 * fft.c -- demo fft procedures
 *
 * Author: Dylan Vanmali
 * How To Write Fast Code 18-642
 *
 * Starter Code
 * https://userpages.umbc.edu/~squire/f14-455/cs455_l18.shtml
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "fft.h"
#define PI 3.14159265358979

/* bitrev(a, k) -- reverse bits 0 thru k-1 in the integer "a"
 */
int bitrev(int a, int k) {
  unsigned int i, b, p, q;
  for (i=b=0, p = 1, q = 1<<(k-1); i<k;  i++, p <<= 1, q >>= 1 )
    if (a & q)
      b |= p;
  return b;
}

/* ilog2(n) -- return an integer log, base 2
 */
int ilog2(int n) {
  int i;
  for (i=8*sizeof(int)-1; i>=0 && ((1<<i) & n)==0; i--);
  return i;
}

/* fft(b, a2, a1, n, sgn) -- do an n-point fft of complex vector b and 
 * return the result in b.  b consists of 2*n FTYPE elements, organized 
 * as n complex pairs real_1, imag_1, real_2, imag_2, ..., real_n, imag_n.
 * the arrays a1 and a2 are used for working storage and each has n FTYPE 
 * elements.  sgn is 1 for an FFT, and -1 for an IFFT.  The procedure is
 * taken from the Cormen, Leiserson, and Rivest Algorithms text, in the 
 * section on efficient FFT implementations.  This implementation wastes 
 * some space; the b array should probably be dropped and the initial and  
 * final staging done in-place in the "a" arrays.
 */
void fft(FTYPE *b, FTYPE *a2, FTYPE *a1, int n, int sgn) {
  int i, j, k, k2, s, m, log2n;
  double wm1, wm2, w1, w2, t1, t2, u1, u2;

  log2n = ilog2(n);

  /* reorder input and split input into real and complex parts */
  for (i=0; i<n; i++) {
    j = bitrev(i,log2n);
    a1[j] = b[2*i];
    a2[j] = b[2*i+1];
  }

  /* loop on FFT stages */
  for (s=1; s<=log2n; s++) {

    m = 1<<s;			/* m = 2^s */
    wm1 = cos(sgn*2*PI/m);	/* wm = exp(q*2*pi*i/m); */			
    wm2 = sin(sgn*2*PI/m);

    w1 = 1.0;
    w2 = 0.0;

    for (j=0; j<m/2; j++) {

      for (k=j; k<n; k+=m) {
        /* t = w*a[k+m/2]; */
        k2 = k+m/2;				
        t1 = w1 * a1[k2]  -  w2 * a2[k2] ;
        t2 = w1 * a2[k2]  +  w2 * a1[k2] ;

        u1 = a1[k];
        u2 = a2[k];

        a1[k] = u1 + t1;
        a2[k] = u2 + t2;

        a1[k2] = u1 - t1;
        a2[k2] = u2 - t2;

      }

      /* w = w * wm; */
      t1 = w1 * wm1  -  w2 * wm2 ;
      w2 = w1 * wm2  +  w2 * wm1 ;
      w1 = t1;
    }
  }
  
  /* flip the final stage */
  for (i=1; i<n/2; i++) {
    t1 = a1[i]; 
    a1[i] = a1[n-i]; 
    a1[n-i] = t1;
    t2 = a2[i]; 
    a2[i] = a2[n-i]; 
    a2[n-i] = t2;
  }

  /* copy out results */
  if (sgn == -1) {
    /* scale by n for ifft */
    for (i=0; i<n; i++) {
      b[2*i]   = a1[i] /= (FTYPE) n;
      b[2*i+1] = a2[i] /= (FTYPE) n;
    }
  }
  else {
    /* just copy the data out */
    for (i=0; i<n; i++) {
      b[2*i]   = a1[i];
      b[2*i+1] = a2[i];
    }
  }
}