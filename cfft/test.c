#include "fft.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    const int N1 = 8; // Prime Number for good-thomas required
    const int N2 = 8;  // Prime Number for good-thomas required
    const int size = N1 * N2;

    complex * input1 = (complex*) malloc(sizeof(struct complex_t) * size);
    complex * input2 = (complex*) malloc(sizeof(struct complex_t) * size);
    complex * result1, * result2;
    
    /* Init inputs */
    for (int i=0; i < size; i++) {
        input1[i].re = (double) i;
        input1[i].im = 0.0;
        input2[i].re = (double) i;
        input2[i].im = 0.0;
    }
    
    /* Do FFT */
    result1 = FFT_CooleyTukey(input1, size, N1, N2);
    result2 = FFT_GoodThomas(input2, size, N1, N2);
    
    /* Compare results */
    printf("Index \t Cooley-Tukey Output \t \t Good-Thomas Output \n");
    for (int i=0; i < size; i++) {
        printf("%d: \t %f + %fi \n", i, result1[i].re, result1[i].im);
        //printf("%f + %fi \n", result2[i].re, result2[i].im);
    }
  
    return 0;
}