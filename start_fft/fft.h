#define FTYPE float

#ifndef FFT_H
#define FFT_H

int bitrev(int a, int k);

int ilog2(int n);

void fft (FTYPE *b, FTYPE *a2, FTYPE *a1, int n, int sgn);

#endif