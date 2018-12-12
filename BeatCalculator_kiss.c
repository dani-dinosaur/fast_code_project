
#include <stdio.h>
#include <math.h>
#include "kiss_fftr.h"
//#include "kiss_fft.h"

#include <immintrin.h>
#include <x86intrin.h>
#include <xmmintrin.h>
// #include "new_fft.h"

static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

void fftrArray(float* sample, int size, kiss_fft_cpx* out) {
    kiss_fftr_cfg cfg;

    if ((cfg = kiss_fftr_alloc(size, 0, NULL, NULL)) == NULL) {
        printf("Not enough memory to allocate fftr!\n");
        exit(-1);
    }

    unsigned long long st, et;
    st = rdtsc(); 

    kiss_fftr(cfg, (kiss_fft_scalar*)sample, out);
    
    et = rdtsc();

    for (int i = 0; i < size; i++)
        printf("out[%d] = %+f + %+fi\n", i, out[i].r, out[i].i);

    printf ("time to take FFT: %llu\n", (et-st));

    kiss_fft_free(cfg);
}

/* detect_beat
 * Returns the BPM of the given mp3 file
 * @Params: s - the path to the desired mp3
 */
int detect_beat(float* sample, int sample_size) {
    // Step 1: Get a 5-second sample of our desired mp3
    // Assume the max frequency is 4096
    unsigned long long st, et;
    st = rdtsc(); 

    kiss_fft_cpx* out = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * (sample_size/2+1));
    fftrArray(sample, sample_size, out);

    et = rdtsc();
    printf("total time: %llu\n", (et-st));

    free(out);
    
    return -1; // BPM;
}

void readInWavFile(char* file_name, float* sample_output, int sample_size) {
    //Open wave file in read mode
    FILE * infile = fopen(file_name,"rb");
    // For counting number of frames in wave file.
    int count = 0;                        
    /// short int used for 16 bit as input data format is 16 bit PCM audio
    short int buff16;
    
    if (infile) {
        fseek(infile,44,SEEK_SET);

        while (count < sample_size) {
            fread(&buff16,sizeof(buff16),1,infile);        // Reading data in chunks of BUFSIZE
            sample_output[count] = buff16;
            count++;                    
        }
    }
}

int main(int argc, char* argv[]) {
    // Test CPU Version
    //takes in 10 seconds at a sampling rate of 48000 samples/sec
    if ( argc != 3 ) {
        printf( "\nUsage: %s filename.wav number_of_samples\n\n", argv[0]);
        return 0;
    }

    int sample_size = atoi(argv[2]);
    float* sample = (float*)malloc(sizeof(float) * sample_size);
    readInWavFile(argv[1], sample, sample_size);

    int BPM = detect_beat(sample, sample_size);
    printf("Final BPM: %i\n", BPM);

    free(sample);

    return 0;
}