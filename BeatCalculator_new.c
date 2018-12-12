#include <stdio.h>
#include <math.h>
#include "fft.h"

#include <immintrin.h>
#include <x86intrin.h>
#include <xmmintrin.h>

static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

/*
 * detect_beat
 *  Returns the BPM of the given mp3 file
 *  @Params: s - the path to the desired mp3
 */
int detect_beat(float* sample, int sample_size) {
    float* a1 = (float*)malloc(sizeof(float) * sample_size);
    float* a2 = (float*)malloc(sizeof(float) * sample_size);
    float* twittle = (float*)malloc(sizeof(float) * sample_size * 2);
    
    // Setup the data and twittle factors beforehand
    fft_data_setup(sample, a1, a2, sample_size);
    fft_twittle_setup(sample_size, 1, twittle);

    // Call FFT function and time how long it takes to run
    unsigned long long st, et;
    st = rdtsc(); 

    fft(sample, a1, a2, sample_size, 1, twittle);

    et = rdtsc();

    // Print out the computed fft value
    printf("Final FFT: \n");
    for (int i=0; i<sample_size * 2; i+=2)
        printf("%f + %fi\n", sample[i], sample[i+1]);
    printf("\n");

    printf("FFT time: %llu\n", (et-st));

    free(a1);
    free(a2);
    free(twittle);

    return -1; // BPM; // TODO: Merge with Derivative and Combfilter
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
            sample_output[count] = (float)buff16;  // real
            sample_output[count+1] = 0.f;          // imaginary
            count += 2;                    
        }
    }
}

int main(int argc, char* argv[]) {
    // Test CPU Version
    //takes in 10 seconds at a sampling rate of 48000 samples/sec
    if ( argc != 3 ) {
        printf( "Usage: %s filename.wav number_of_samples\n", argv[0]);
        return 0;
    }

    const int sample_size = atoi(argv[2]);
    float* sample = (float*)malloc(sizeof(float) * sample_size * 2); // sample[i]=real, sample[i+1]=imaginary
    readInWavFile(argv[1], sample, sample_size * 2);

    int BPM = detect_beat(sample, sample_size);
    printf("Final BPM: %i\n\n", BPM);

    free(sample);

    return 0;
}