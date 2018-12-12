#include "fft.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926535897932384626434

complex* DFT_naive(complex* x, int N) {
    complex complex_value0, complex_value1, complex_value2, complex_value3;
    complex X_0, X_1, X_2, X_3;
    double radians0, radians1, radians2, radians3;
    complex* X = (complex*) malloc(sizeof(struct complex_t) * N);
    const double multiplier = -2*PI/N;

    int k, n;
    for(k = 0; k < N; k++) {
        X[k].re = 0.0;
        X[k].im = 0.0;
        for(n = 0; n < N; n+=4) { // N shouldn't be expanded since prime
            radians0 = multiplier*(n+0)*k;
            complex_value0.re = cos(radians0);
            complex_value0.im = sin(radians0);
            X_0 = multiply(x[n+0], complex_value0);

            radians1 = multiplier*(n+1)*k;
            complex_value1.re = cos(radians1);
            complex_value1.im = sin(radians1);
            X_1 = multiply(x[n+1], complex_value1);

            radians2 = multiplier*(n+2)*k;
            complex_value2.re = cos(radians2);
            complex_value2.im = sin(radians2);
            X_2 = multiply(x[n+2], complex_value2);

            radians3 = multiplier*(n+3)*k;
            complex_value3.re = cos(radians3);
            complex_value3.im = sin(radians3);
            X_3 = multiply(x[n+3], complex_value3);

            X[k] = add(X[k], add(add(X_0, X_1), add(X_2, X_3)));
        }
    }
    
    return X;
}

static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

/** Implements the Good-Thomas FFT algorithm. 
  *
  * @expects: N1 and N2 must be relatively prime
  * @expects: N1*N2 = N
  */
// complex* FFT_GoodThomas(complex* input, int N, int N1, int N2) {
//     int i, j, k1, k2, z;

//     /* Allocate columnwise matrix */
//     complex** columns = (complex**) malloc(sizeof(struct complex_t*) * N1);
//     for(k1 = 0; k1 < N1; k1++) {
//         columns[k1] = (complex*) malloc(sizeof(struct complex_t) * N2);
//     }
    
//     /* Allocate rowwise matrix */
//     complex** rows = (complex**) malloc(sizeof(struct complex_t*) * N2);
//     for(k2 = 0; k2 < N2; k2++) {
//         rows[k2] = (complex*) malloc(sizeof(struct complex_t) * N1);
//     }
    
//     /* Reshape input into N1 columns (Using Good-Thomas Indexing) */
//     for(z = 0; z < 30; z++) {
//         k1 = z % N1;
//         k2 = z % N2;
//         columns[k1][k2] = input[z];
//     }
    
//     /* Compute N1 DFTs of length N2 using naive method */
//     for (k1 = 0; k1 < N1; k1++) {
//         columns[k1] = DFT_naive(columns[k1], N2);
//     }
    
//     /* Cache Aware Transpose */
//     for (j = 0; j < N1; j+=2) {
//         for (i = 0; i < N2; i+=4) {
//             *(*rows + N1*(i+0) + j) = *(*columns + N2*j + (i+0));
//             *(*rows + N1*(i+1) + j) = *(*columns + N2*j + (i+1));
//             *(*rows + N1*(i+2) + j) = *(*columns + N2*j + (i+2));
//             *(*rows + N1*(i+3) + j) = *(*columns + N2*j + (i+3));

//             *(*rows + N1*(i+0) + (j+1)) = *(*columns + N2*(j+1) + (i+0));
//             *(*rows + N1*(i+1) + (j+1)) = *(*columns + N2*(j+1) + (i+1));
//             *(*rows + N1*(i+2) + (j+1)) = *(*columns + N2*(j+1) + (i+2));
//             *(*rows + N1*(i+3) + (j+1)) = *(*columns + N2*(j+1) + (i+3));
//         }
//     }
   
//     /* Compute N2 DFTs of length N1 using naive method */
//     for (k2 = 0; k2 < N2; k2++) {
//         rows[k2] = DFT_naive(rows[k2], N1);
//     }
    
//     /* Flatten into single output (Using chinese remainder theorem) */
//     complex* output = (complex*) malloc(sizeof(struct complex_t) * N);
    
//     for(k1 = 0; k1 < N1; k1++) {
//         for (k2 = 0; k2 < N2; k2++) {
//             z = N1*k2 + N2*k1;
//             output[z%N] = rows[k2][k1];
//         }
//     }

//     /* Free all alocated memory except output and input arrays */
//     for(k1 = 0; k1 < N1; k1++) {
//         free(columns[k1]);
//     }
//     for(k2 = 0; k2 < N2; k2++) {
//         free(rows[k2]);
//     }
//     free(columns);
//     free(rows);
//     return output;
// }

/** Implements the Cooley-Tukey FFT algorithm. 
  *
  * @expects: N1*N2 = N
  */
complex* FFT_CooleyTukey(complex* input, int N, int N1, int N2) {
    int i, j, k1, k2;

    unsigned long long st, et;
    

    /* Allocate columnwise matrix */
    complex** columns = (complex**) malloc(sizeof(struct complex_t*) * N1);
    for(k1 = 0; k1 < N1; k1++) {
        columns[k1] = (complex*) malloc(sizeof(struct complex_t) * N2);
    }
    
    /* Allocate rowwise matrix */
    complex** rows = (complex**) malloc(sizeof(struct complex_t*) * N2);
    for(k2 = 0; k2 < N2; k2++) {
        rows[k2] = (complex*) malloc(sizeof(struct complex_t) * N1);
    }
    
    /* Reshape input into N1 columns */
    for (k1 = 0; k1 < N1; k1++) {
        for(k2 = 0; k2 < N2; k2++) {
            columns[k1][k2] = input[N1*k2 + k1];
        }
    }

    st = rdtsc();
    /* Compute N1 DFTs of length N2 using naive method */
    // for (k1 = 0; k1 < N1; k1++) {
    //     columns[k1] = DFT_naive(columns[k1], N2);
    // }
    
    /* Multiply by the twiddle factors  ( e^(-2*pi*j/N * k1*k2)) and Cache-Aware Transpose */
    complex complex_value0, complex_value1, complex_value2, complex_value3,
            complex_value4, complex_value5, complex_value6, complex_value7;
    double radians0, radians1, radians2, radians3, radians4, radians5, radians6, radians7;
    for(j = 0; j < N1; j+=2) {
        /* Compute N1 DFTs of length N2 using naive method */
        columns[j]   = DFT_naive(columns[j], N2);
        columns[j+1] = DFT_naive(columns[j+1], N2);
        /* twiddle factors and transpose */
        for (i = 0; i < N2; i+=4) {
            radians0 = -2.0*PI*(j+0)*(i+0)/N;
            radians1 = -2.0*PI*(j+0)*(i+1)/N;
            radians2 = -2.0*PI*(j+0)*(i+2)/N;
            radians3 = -2.0*PI*(j+0)*(i+3)/N;
            radians4 = -2.0*PI*(j+1)*(i+0)/N;
            radians5 = -2.0*PI*(j+1)*(i+1)/N;
            radians6 = -2.0*PI*(j+1)*(i+2)/N;
            radians7 = -2.0*PI*(j+1)*(i+3)/N;

            complex_value0.re = cos(radians0);
            complex_value0.im = sin(radians0);
            complex_value1.re = cos(radians1);
            complex_value1.im = sin(radians1);
            complex_value2.re = cos(radians2);
            complex_value2.im = sin(radians2);
            complex_value3.re = cos(radians3);
            complex_value3.im = sin(radians3);
            complex_value4.re = cos(radians4);
            complex_value4.im = sin(radians4);
            complex_value5.re = cos(radians5);
            complex_value5.im = sin(radians5);
            complex_value6.re = cos(radians6);
            complex_value6.im = sin(radians6);
            complex_value7.re = cos(radians7);
            complex_value7.im = sin(radians7);

            rows[i+0][j+0] = multiply(complex_value0, columns[j+0][i+0]);
            rows[i+1][j+0] = multiply(complex_value1, columns[j+0][i+1]);
            rows[i+2][j+0] = multiply(complex_value2, columns[j+0][i+2]);
            rows[i+3][j+0] = multiply(complex_value3, columns[j+0][i+3]);
            rows[i+0][j+1] = multiply(complex_value4, columns[j+1][i+0]);
            rows[i+1][j+1] = multiply(complex_value5, columns[j+1][i+1]);
            rows[i+2][j+1] = multiply(complex_value6, columns[j+1][i+2]);
            rows[i+3][j+1] = multiply(complex_value7, columns[j+1][i+3]);
        }
    }

    /* Compute N2 DFTs of length N1 using naive method */
    // Can parallelize this
    for (k2 = 0; k2 < N2; k2++) {
        rows[k2] = DFT_naive(rows[k2], N1);
    }

    et = rdtsc();

    /* Flatten into single output */
    complex* output = (complex*) malloc(sizeof(struct complex_t) * N);
    for(k1 = 0; k1 < N1; k1++) {
        for (k2 = 0; k2 < N2; k2++) {
            output[N2*k1 + k2] = rows[k2][k1];
        }
    }

    /* Free all alocated memory except output and input arrays */
    for(k1 = 0; k1 < N1; k1++) {
        free(columns[k1]);
    }
    for(k2 = 0; k2 < N2; k2++) {
        free(rows[k2]);
    }
    printf ("time to take FFT input: %llu\n", (et-st));

    free(columns);
    free(rows);
    return output;
}