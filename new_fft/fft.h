#define FTYPE float

#ifndef FFT_H
#define FFT_H

// #define _OPENMP // Uncomment to turn on openmp code <- right now lower cycles

int bitrev(int a, int k);

int ilog2(int n);

void fft (FTYPE *b, FTYPE *a2, FTYPE *a1, int n, int sgn, float* twittle);

void fft16(float Z[32]);

void fft_data_setup(FTYPE *b, FTYPE *a2, FTYPE *a1, int n);

void fft_twittle_setup(int n, int sgn, float* twittle);

#endif