#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "wave.h"
#include "cfft/fft.h"
#include "cfft/complex.h"

// #include <kiss_fftr.h>
// #include <kiss_fft.h>

// #include <immintrin.h>
// #include <x86intrin.h>
// #include <xmmintrin.h>


static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

// void derivative(float* sample, float* dif_sample, int sample_size) {

//     float constant = 48000.0/2;
//     __m256 multiply_constant = _mm256_set1_ps(constant);
//     __m256 front, back;
//     __m256 front2, front3, front4, front5, front6, front7;
//     int j, i = 1;
//     int loop_limit = (sample_size/32)*32;

//     unsigned long long st_1, et_1;

//     st_1 = rdtsc();

//     dif_sample[0] = sample[0];

//     float *fptr = &sample[i+1], *bptr = &sample[i-1], *difptr = &dif_sample[i];

//     //#pragma omp parallel for 
//     for (i = 1; i<loop_limit; i +=64) {
//         front = _mm256_loadu_ps(fptr);
//         back = _mm256_loadu_ps(bptr);
//         front = _mm256_sub_ps(front, back);
//         front = _mm256_mul_ps(front, multiply_constant);
        
//         _mm256_storeu_ps(difptr, front);

//         front = _mm256_loadu_ps(fptr+8);
//         back = _mm256_loadu_ps(bptr+8);
//         front = _mm256_sub_ps(front, back);
//         front = _mm256_mul_ps(front, multiply_constant);

//         _mm256_storeu_ps(difptr+8, front);

//         front2 = _mm256_loadu_ps(fptr+16);
//         back = _mm256_loadu_ps(bptr+16);
//         front2 = _mm256_sub_ps(front2, back);
//         front2 = _mm256_mul_ps(front2, multiply_constant);

//         _mm256_storeu_ps(difptr+16, front2);

//         front3 = _mm256_loadu_ps(fptr+24);
//         back = _mm256_loadu_ps(bptr+24);
//         front3 = _mm256_sub_ps(front3, back);
//         front3 = _mm256_mul_ps(front3, multiply_constant);

//         _mm256_storeu_ps(difptr+24, front3);

//         front4 = _mm256_loadu_ps(fptr+32);
//         back = _mm256_loadu_ps(bptr+32);
//         front4 = _mm256_sub_ps(front4, back);
//         front4 = _mm256_mul_ps(front4, multiply_constant);

//         _mm256_storeu_ps(difptr+32, front4);

//         front5 = _mm256_loadu_ps(fptr+40);
//         back = _mm256_loadu_ps(bptr+40);
//         front5 = _mm256_sub_ps(front5, back);
//         front5 = _mm256_mul_ps(front5, multiply_constant);

//         _mm256_storeu_ps(difptr+40, front5);

//         front6 = _mm256_loadu_ps(fptr+48);
//         back = _mm256_loadu_ps(bptr+48);
//         front6 = _mm256_sub_ps(front6, back);
//         front6 = _mm256_mul_ps(front6, multiply_constant);

//         _mm256_storeu_ps(difptr+48, front6);

//         front7 = _mm256_loadu_ps(fptr+56);
//         back = _mm256_loadu_ps(bptr+56);
//         front7 = _mm256_sub_ps(front7, back);
//         front7 = _mm256_mul_ps(front7, multiply_constant);

//         _mm256_storeu_ps(difptr+56, front7);

//         fptr+=64;
//         bptr+=64;
//         difptr+=64;
//     }
//     //#pragma omp parallel for 
//     for (j=i; j < sample_size - 1; j++) {
//         dif_sample[i] = constant * (sample[i+1] - sample[i-1]);
//     }

//     dif_sample[sample_size - 1] = sample[sample_size-1];

//     et_1 = rdtsc();

//     printf ("time to take derivative: %llu\n", (et_1-st_1));

// }

void readInWavFile(char* file_name, complex* sample_output, int sample_size) {
    //Open wave file in read mode
    FILE * infile = fopen(file_name,"rb");
    // For counting number of frames in wave file.
    int count = 0;                        
    /// short int used for 16 bit as input data format is 16 bit PCM audio
    short int buff16;
    
    if (infile) {
        fseek(infile,44,SEEK_SET); // Skip of head and get only PCM data

        while (count < sample_size) {
            fread(&buff16,sizeof(buff16),1,infile);        // Reading data in chunks of BUFSIZE
            sample_output[count].re = (double)buff16;
            sample_output[count].im = 0.0;
            count++;
        }
    }
}

void derivative(complex* sample, complex* dif_sample, int sample_size) {
    int i=1;

    const int constant = sample_size/2; // This is the divisor for difference
    //__m256 multiply_constant = _mm256_set1_ps(constant);

    dif_sample[0].re = 0.0; // Start differential is 0
    complex *fptr = &sample[i+1], *bptr = &sample[i-1], *difptr = &dif_sample[i];
    //#pragma omp parallel for 
    for (i = 1; i<sample_size-1; i ++) {
        difptr->re = (fptr->re - bptr->re) / constant;

        // Go to the next set of values
        fptr++;
        bptr++;
        difptr++;
    }
    dif_sample[sample_size - 1].re = 0.0; // End differential is 0
}

// http://www.flipcode.com/misc/BeatDetectionAlgorithms.pdf
int combfilter(complex* fft_array, const int sample_size, int start, int fin, int step) {
    double AmpMax = 1.0;
    int energyCount = (fin - start)/step;
    double E[energyCount];
    unsigned long long st, et;
    int i, k;
    int n = sample_size/2+1;

    printf("Give me a moment... calculating all BPMs from %d to %d with steps of %d\n",
        start, fin, step);

    /******************** Iterate through all BPMs calculating FFT *********************/
    complex *out1 = malloc(energyCount*n*sizeof(complex));
    for (i = 0; i < energyCount; i++) {
        int BPM = start + i * step;
        int Ti = 60 * 44100/BPM;

        complex *l = malloc(sample_size*sizeof(complex));

        for (k = 0; k < sample_size; k++) {
            if ((k % Ti) == 0) {
                l[k].re = (double)AmpMax;
            }
            else {
                l[k].re = 0.0;
            }
            l[k].im = 0.0;
        }
        st = rdtsc();
        complex* tmp = FFT_CooleyTukey(l, sample_size, sqrt(sample_size), sqrt(sample_size));
        et = rdtsc();
        printf("time to take fft: %llu\n", (et-st));

        free(tmp);
        free(l);
        memcpy(&out1[i*n], &(*tmp), sizeof(complex *) * n);
        // for (int a=0; a<n; a++)
        //     printf("%d.%d: \t %f + %fi\n", i, a, tmp[a].re, tmp[a].im);
    }

    // Convert to put energies near each other
    complex *out = malloc(energyCount*n*sizeof(complex));
    for (i = 0; i < energyCount; i++) {
        for (k = 0; k < n; k++){
            out[k*energyCount+i]  = out1[i*n + k];
        }
    }

    // Free out1
    free(out1);

    printf("Thank You for waiting... %d BPMs were calculated\n", energyCount);

    /************************ Code for energy dot product ************************/
    float a0, b0, a1, b1, a2, b2, a3, b3, a4, b4;
    st = rdtsc();
    for (k = 0; k < n; k++) {       
        for (i = 0; i < (energyCount); i+=5) {
            a0 = fft_array[k].re * out[k*energyCount + i    ].re - fft_array[k].im * out[k*energyCount + i    ].im;
            b0 = fft_array[k].re * out[k*energyCount + i    ].im + fft_array[k].im * out[k*energyCount + i    ].re;
            a1 = fft_array[k].re * out[k*energyCount + i + 1].re - fft_array[k].im * out[k*energyCount + i + 1].im;
            b1 = fft_array[k].re * out[k*energyCount + i + 1].im + fft_array[k].im * out[k*energyCount + i + 1].re;
            a2 = fft_array[k].re * out[k*energyCount + i + 2].re - fft_array[k].im * out[k*energyCount + i + 2].im;
            b2 = fft_array[k].re * out[k*energyCount + i + 2].im + fft_array[k].im * out[k*energyCount + i + 2].re;
            a3 = fft_array[k].re * out[k*energyCount + i + 3].re - fft_array[k].im * out[k*energyCount + i + 3].im;
            b3 = fft_array[k].re * out[k*energyCount + i + 3].im + fft_array[k].im * out[k*energyCount + i + 3].re;
            a4 = fft_array[k].re * out[k*energyCount + i + 4].re - fft_array[k].im * out[k*energyCount + i + 4].im;
            b4 = fft_array[k].re * out[k*energyCount + i + 4].im + fft_array[k].im * out[k*energyCount + i + 4].re;

            E[i]   += (a0*a0 + b0*b0);
            E[i+1] += (a1*a1 + b1*b1);
            E[i+2] += (a2*a2 + b2*b2);
            E[i+3] += (a3*a3 + b3*b3);
            E[i+4] += (a4*a4 + b4*b4);
        }
    }
    et = rdtsc();
    free(out);
    printf("time to take energy: %llu\n", (et-st));

    /************************ Calculate max of E[k] ************************/
    double max_val = -1;
    int index = -1; 
    int found = 0;

    st = rdtsc();
    for (i = 0; i < energyCount; i++) {
        if (E[i] >= max_val) {
             if (E[i] >= max_val *1.9 && i != 0 && found == 0){
                max_val = E[i];
                index = i;
                found = 1; 
             }
            else if (found ==  1 && E[i] >= max_val *1.1) {
                max_val = E[i];
                index = i;
            }
        }
    }
    et = rdtsc();
    printf ("time to find max energy: %llu\n", (et-st));
    int final = (int)((float)(start + index * step) *(float)(2.0/3.0));
    printf("final val is %i\n", final);

    return final;
}


/* detect_beat
 * Returns the BPM of the given mp3 file
 * @Params: s - the path to the desired mp3
 */
int detect_beat(complex* sample, int sample_size) {
    unsigned long long st, et;
    unsigned long long st_1, et_1;

    st = rdtsc();

    // Step 1: Take Differential
    st_1 = rdtsc();
    complex *dif_sample = malloc(sample_size*sizeof(complex));
    derivative(sample, dif_sample, sample_size);
    et_1 = rdtsc();
    printf ("time to take Derivative: %llu\n", (et_1-st_1));

    // Step 2: Take FFT of input
    //st_1 = rdtsc();
    complex* result = malloc(sample_size*sizeof(complex));
    result = FFT_CooleyTukey(dif_sample, sample_size, sqrt(sample_size), sqrt(sample_size));
    //et_1 = rdtsc();
    //printf ("time to take FFT input: %llu\n", (et_1-st_1));

    // Step 2: Combfilter (FFT is applied on beats from start to fin range with steps)
    st_1 = rdtsc(); 
    int BPM = combfilter(dif_sample, sample_size, 15, 200, 5);
    et_1 = rdtsc();
    printf ("time to take Combfilter: %llu\n", (et_1-st_1));
    et = rdtsc();

    // Prints a complex[sample_size] array
    // for (int i=0; i < sample_size; i++) {
    //     printf("%d: \t %f + %fi\n", i, sample[i].re, sample[i].im);

    // Free Malloc'd Values
    free(dif_sample);
    free(result);

    printf ("total time: %llu\n", (et-st));
    return BPM;
}

int main(int argc, char* argv[]) {
    // Test CPU Version
    // takes in 10 seconds at a sampling rate of 48000 samples/sec

    if ( argc != 2 ) {
        printf( "Usage: %s filename.wav\n", argv[0]);
        return 0;
    }

    // Step 0: Load wav

    // Method 1) Read and interpret full wave file, accurate input, very long don't recommend
    // const int sample_size = 2709800; // Obtained from reading a wav file
    // complex *wav_data = malloc(sample_size*sizeof(complex)); // Obtained from reading a wav file
    // readWAV(argc, argv, wav_data);

    // Method 2) Short, quick and specify your own sample_size
    const int sample_size = 32768; // 2^15 - Smaller size is quicker
    complex *wav_data = malloc(sample_size*sizeof(complex)); // Obtained from reading a wav file
    readInWavFile(argv[1], wav_data, sample_size); // Assumes input must be a WAV file using 16-bit PCM audio

    // Prints a wav_data[sample_size] array
    // for (int i=0; i<sample_size; i++)
    //     printf("%f+%fi\n", wav_data[i].re, wav_data[i].im);

    // Step 1: Obtain BPM and print result
    int BPM = detect_beat(wav_data, sample_size);
    printf("Final BPM: %i\n", BPM);

    free(wav_data);

    return 0;
}

