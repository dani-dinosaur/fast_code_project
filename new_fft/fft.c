/*
 *  fft.c -- demo fft procedures
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <smmintrin.h>
#include <immintrin.h>

#include "fft.h"
#define PI 3.14159265358979

/* bitrev(a, k) -- reverse bits 0 thru k-1 in the integer "a"
 */
int bitrev(int a, int k) {
  unsigned int i, b, p, q;
  for (i=b=0, p = 1, q = 1<<(k-1); 
       i<k; 
       i++, p <<= 1, q >>= 1 ) if (a & q) b |= p;
  return b;
}

/* ilog2(n) -- return an integer log, base 2
 */
int ilog2(int n) {
  int i;
  for (i=8*sizeof(int)-1; i>=0 && ((1<<i) & n)==0; i--);
  return i;
}

/* SetUp data needs to be called for all new data samples
 * b = original sample
 * a2 = real numbers
 * a1 = imaginary numbers
 * n = sample_size
 */
void fft_data_setup(FTYPE *b, FTYPE *a2, FTYPE *a1, int n) {
  int i, j, log2n;

  log2n = ilog2(n);

  /* reorder input and split input into real and complex parts */
  for (i=0; i<n; i++) {
    j = bitrev(i,log2n);
    a1[j] = b[2*i];
    a2[j] = b[2*i+1];
  }
}

/* SetUp twittle needs to be called only when sample size changes
 *
 * sets up the fft twittle values and constants independent to the fft
 * n = sample_size
 * sgn = 1 for fft, -1 for ifft
 * twittle = an array of size 2 * n
 */
void fft_twittle_setup(int n, int sgn, float* twittle) {
  int j, s, m, log2n;
  float wm1, wm2, w1, w2, t1;

  log2n = ilog2(n);

  /* calculate the twittle factors */
  int index = 0;
  for (s=1; s<=log2n; s++) {
    m = 1<<s;     /* m = 2^s */
    wm1 = cos(sgn*2*PI/m);  /* wm = exp(q*2*pi*i/m); */     
    wm2 = sin(sgn*2*PI/m);

    w1 = 1.0;
    w2 = 0.0;

    for (j=0; j<m/2; j++) {
      twittle[index] = w1;
      twittle[index+1] = w2;
      index+=2;

      /* w = w * wm; */
      t1 = w1 * wm1  -  w2 * wm2 ;
      w2 = w1 * wm2  +  w2 * wm1 ;
      w1 = t1;
    }
  }
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
void fft(FTYPE *b, FTYPE *a2, FTYPE *a1, int n, int sgn, float* twittle) {
  int i, j, k, k1, k2, k3, k4, s, m, log2n;
  float w1, w2, t1, t2, t3, t4, u1, u2, u3, u4;

  log2n = ilog2(n);

  int index = 0;
  /* loop on FFT stages */
  for (s=1; s<=log2n; s++) {

    m = 1<<s;     /* m = 2^s */

    for (j=0; j<m/2; j++) {
      w1 = twittle[index];
      w2 = twittle[index+1];;
      index+=2;


      // Check number of iterations and send to that optimization type
      switch(n/m) {
        case 1: // Loop unravel for one
          k2 = j+m/2;
          t1 = w1 * a1[k2]  -  w2 * a2[k2] ;
          t2 = w1 * a2[k2]  +  w2 * a1[k2] ;

          u1 = a1[j];
          u2 = a2[j];

          a1[j] = u1 + t1;
          a2[j] = u2 + t2;

          a1[k2] = u1 - t1;
          a2[k2] = u2 - t2;
          break;
        case 2:  // Loop unravel for two
          k1 = j;
          k2 = k1+m/2;
          t1 = w1 * a1[k2]  -  w2 * a2[k2] ;
          t2 = w1 * a2[k2]  +  w2 * a1[k2] ;

          k3 = j+m;
          k4 = k3+m/2;
          t3 = w1 * a1[k4]  -  w2 * a2[k4] ;
          t4 = w1 * a2[k4]  +  w2 * a1[k4] ;

          u1 = a1[k1];
          u2 = a2[k1];
          u3 = a1[k3];
          u4 = a2[k3];

          #ifdef _OPENMP
          float arr_u[8] = {0.f,0.f,0.f,0.f, u1, u2, u3, u4};
          float arr_v[8] = {0.f,0.f,0.f,0.f, t1, t2, t3, t4};
          __m256 u = _mm256_load_ps(arr_u);
          __m256 v = _mm256_load_ps(arr_v);
          __m256 sum = _mm256_add_ps(u,v);
          __m256 sub = _mm256_add_ps(u,v);
          float result_sum[8];
          _mm256_store_ps (result_sum, sum);
          float result_sub[8];
          _mm256_store_ps (result_sub, sub);
          a1[k1] = result_sum[4];
          a2[k1] = result_sum[5];
          a1[k3] = result_sum[6];
          a2[k3] = result_sum[7];
          a1[k2] = result_sub[4];
          a2[k2] = result_sub[5];
          a1[k4] = result_sub[6];
          a2[k4] = result_sub[7];
          #else
          a1[k1] = u1 + t1;
          a2[k1] = u2 + t2;
          a1[k2] = u1 - t1;
          a2[k2] = u2 - t2;
          a1[k3] = u3 + t3;
          a2[k3] = u4 + t4;
          a1[k4] = u3 - t3;
          a2[k4] = u4 - t4;
          break;
          #endif

        default: // Large loops don't unravel (probability this occurs is small)
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
          break;
      }

    }
  }
  
  /* transpose 1D Array */
  for (i=1; i<n/2; i++) {
    t1 = a1[i]; // real
    a1[i] = a1[n-i];
    a1[n-i] = t1;

    t2 = a2[i]; // imag
    a2[i] = a2[n-i];
    a2[n-i] = t2;
  }

  /* copy out results */
  if (sgn == -1)
    /* scale by n for ifft */
    for (i=0; i<n; i++) {
      b[2*i] = a1[i] /= (FTYPE) n;
      b[2*i+1] = a2[i] /= (FTYPE) n;
    }
  else
    /* just copy the data out */
    for (i=0; i<n; i++) {
      b[2*i] = a1[i];
      b[2*i+1] = a2[i];
    }
}
