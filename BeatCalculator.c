
#include <stdio.h>
#include <math.h>
#include <kiss_fftr.h>
#include <kiss_fft.h>

#include <immintrin.h>
#include <x86intrin.h>
#include <xmmintrin.h>


static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}


void readInWavFile(char* wav_file, float* output, int sample_size){
	
	//Open wave file in read mode
	FILE * infile = fopen(wav_file,"rb");        
   // For counting number of frames in wave file.
    int count = 0;                        
    /// short int used for 16 bit as input data format is 16 bit PCM audio
    short int buff16;

    if (infile){
        fseek(infile,44,SEEK_SET);
        while (count < sample_size){
            fread(&buff16,sizeof(buff16),1,infile);        // Reading data in chunks of BUFSIZE
            output[count] = buff16;
            count++;                    
    	}
	}
}
/* detect_beat
 * Returns the derivative of input sample
 */
void derivative(float* sample, float* dif_sample, int sample_size) {

    float constant = 48000.0/2;
    __m256 multiply_constant = _mm256_set1_ps(constant);
    __m256 front, back;
    __m256 front1, front2, front3, front4, front5, front6, front7;
    int i = 1;
    int loop_limit = (sample_size/32)*32;

    unsigned long long st_1, et_1;

    st_1 = rdtsc();

    dif_sample[0] = sample[0];

    float *fptr = &sample[i+1], *bptr = &sample[i-1], *difptr = &dif_sample[i];

    for (i = 1; i<loop_limit; i +=64) {
        front = _mm256_loadu_ps(fptr);
        back = _mm256_loadu_ps(bptr);
        front = _mm256_sub_ps(front, back);
        front = _mm256_mul_ps(front, multiply_constant);
        
        _mm256_storeu_ps(difptr, front);

        front = _mm256_loadu_ps(fptr+8);
        back = _mm256_loadu_ps(bptr+8);
        front = _mm256_sub_ps(front, back);
        front = _mm256_mul_ps(front, multiply_constant);

        _mm256_storeu_ps(difptr+8, front);

        front2 = _mm256_loadu_ps(fptr+16);
        back = _mm256_loadu_ps(bptr+16);
        front2 = _mm256_sub_ps(front2, back);
        front2 = _mm256_mul_ps(front2, multiply_constant);

        _mm256_storeu_ps(difptr+16, front2);

        front3 = _mm256_loadu_ps(fptr+24);
        back = _mm256_loadu_ps(bptr+24);
        front3 = _mm256_sub_ps(front3, back);
        front3 = _mm256_mul_ps(front3, multiply_constant);

        _mm256_storeu_ps(difptr+24, front3);

        front4 = _mm256_loadu_ps(fptr+32);
        back = _mm256_loadu_ps(bptr+32);
        front4 = _mm256_sub_ps(front4, back);
        front4 = _mm256_mul_ps(front4, multiply_constant);

        _mm256_storeu_ps(difptr+32, front4);

        front5 = _mm256_loadu_ps(fptr+40);
        back = _mm256_loadu_ps(bptr+40);
        front5 = _mm256_sub_ps(front5, back);
        front5 = _mm256_mul_ps(front5, multiply_constant);

        _mm256_storeu_ps(difptr+40, front5);

        front6 = _mm256_loadu_ps(fptr+48);
        back = _mm256_loadu_ps(bptr+48);
        front6 = _mm256_sub_ps(front6, back);
        front6 = _mm256_mul_ps(front6, multiply_constant);

        _mm256_storeu_ps(difptr+48, front6);

        front7 = _mm256_loadu_ps(fptr+56);
        back = _mm256_loadu_ps(bptr+56);
        front7 = _mm256_sub_ps(front7, back);
        front7 = _mm256_mul_ps(front7, multiply_constant);

        _mm256_storeu_ps(difptr+56, front7);

        fptr+=64;
        bptr+=64;
        difptr+=64;
    }

    for (i=i; i < sample_size - 1; i++) {
        dif_sample[i] = constant * (sample[i+1] - sample[i-1]);
    }
    dif_sample[sample_size - 1] = sample[sample_size-1];

    et_1 = rdtsc();

    printf ("time to take derivative: %lu\n", (et_1-st_1));


}

/* fftrArray
 * Returns the fft of a sample (uses KissFFT code)
 */
void fftrArray(float* sample, int size, kiss_fft_cpx* out) {
    kiss_fftr_cfg cfg;

    int i;

    if ((cfg = kiss_fftr_alloc(size, 0, NULL, NULL)) == NULL) {
        printf("Not enough memory to allocate fftr!\n");
        exit(-1);
    }
    kiss_fftr(cfg, (kiss_fft_scalar*)sample, out);
    free(cfg);

}
/* fftArray
 * Returns the fft of a sample (uses KissFFT code)
 */
void fftArray(unsigned int* sample, int size, kiss_fft_cpx* out) {
  kiss_fft_cpx in[size/2];
  kiss_fft_cfg cfg;
  int i;

  if ((cfg = kiss_fft_alloc(size/2, 0, NULL, NULL)) == NULL) {
    printf("Not Enough Memory?!?");
    exit(-1);
  }
  //set real components to one side of stereo input, complex to other
  for(i=0; i < size; i+=2) {
    in[i/2].r = sample[i];
    in[i/2].i = sample[i+1];
  }

  kiss_fft(cfg, in, out);
  free(cfg);

}
/* generate_comb_filters
 * Returns an array of fft's of comb filters with varying BPM
 */
kiss_fft_cpx * generate_comb_filters( int size, int sample_size, int start, int fin, int step){
    float AmpMax =0.0001;
    int energyCount = (fin - start)/step;
    int count = 0;
    int i, k;
    int n = sample_size/2+1;
    
    kiss_fft_cpx *out1= (kiss_fft_cpx*)malloc(energyCount*(n)*sizeof(kiss_fft_cpx));
   
    // Iterate through all possible BPMs
    for (i = 0; i < energyCount; i++) {
        int BPM = start + i * step;
        int Ti = 60 * 44100/BPM;
        float l[sample_size];
        count = 0;

        for (k = 0; k < sample_size; k++) {
            if ((k % Ti) == 0) {
                count++;
                l[k] = (float)AmpMax;    
            }
            else {
                l[k] = -0.0001;
            }
        }

        fftrArray(l, sample_size, &out1[i*n]);
    }

    kiss_fft_cpx *out= (kiss_fft_cpx*)malloc(energyCount*(n)*sizeof(kiss_fft_cpx));
    for (i = 0; i < (energyCount); i++) {
        for (k = 0; k < n; k++){
            out[k*energyCount+i]  = out1[i*n + k];

        }
    }
    free(out1);

    FILE *f = fopen("combfilters.data", "wb");
    fwrite(out, sizeof(kiss_fft_cpx), energyCount*(n)*sizeof(kiss_fft_cpx), f);
    fclose(f);

    return out;
}
/* combfilter_mult
 * Returns the maximum dot product between an array of comb filters and a sample
 */
int combfilter_mult(kiss_fft_cpx fft_array[], int size, int sample_size, int start, int fin, int step) {
    
    float a, b, a1, b1, a2, b2, a3, b3, a4, b4;
    int i, k;
    int energyCount = (fin - start)/step;
    double E[energyCount];
    int n = sample_size/2+1;
    unsigned long long st, et;

    //kiss_fft_cpx *out= generate_comb_filters(size, sample_size, start, fin, step);
    
    kiss_fft_cpx *out= (kiss_fft_cpx*)malloc(energyCount*(n)*sizeof(kiss_fft_cpx));
    FILE *ifp = fopen("combfilters.data", "rb"); 
    fread(out, sizeof(kiss_fft_cpx), energyCount*(n)*sizeof(kiss_fft_cpx), ifp);

    st = rdtsc();

    //take dot product between comb filters and sample
    for (k = 0; k < n; k++){       
        for (i = 0; i < (energyCount); i+=5) {    
    
            //if (k % 5 == 0){
            a = fft_array[k].r * out[k*energyCount + i].r - fft_array[k].i * out[k*energyCount + i].i;
            b = fft_array[k].r * out[k*energyCount + i].i + fft_array[k].i * out[k*energyCount + i].r;
            a1 = fft_array[k].r * out[k*energyCount + i+ 1].r - fft_array[k].i * out[k*energyCount + i + 1].i;
            b1 = fft_array[k].r * out[k*energyCount + i+ 1].i + fft_array[k].i * out[k*energyCount + i + 1].r;
            a2 = fft_array[k].r * out[k*energyCount + i+ 2].r - fft_array[k].i * out[k*energyCount + i + 2].i;
            b2 = fft_array[k].r * out[k*energyCount + i+ 2].i + fft_array[k].i * out[k*energyCount + i + 2].r;
            a3 = fft_array[k].r * out[k*energyCount + i+ 3].r - fft_array[k].i * out[k*energyCount + i + 3].i;
            b3 = fft_array[k].r * out[k*energyCount + i+ 3].i + fft_array[k].i * out[k*energyCount + i + 3].r;
            a4 = fft_array[k].r * out[k*energyCount + i+ 4].r - fft_array[k].i * out[k*energyCount + i + 4].i;
            b4 = fft_array[k].r * out[k*energyCount + i+ 4].i + fft_array[k].i * out[k*energyCount + i + 4].r;

            E[i] += (a*a + b*b);
            E[i+1] += (a1*a1 + b1*b1);
            E[i+2] += (a2*a2 + b2*b2);
            E[i+3] += (a3*a3 + b3*b3);
            E[i+4] += (a4*a4 + b4*b4);
        }
    }

    et = rdtsc();
    printf ("time to take dot product: %lu\n", (et-st));

    free(out);
    //Calculate max of E[k]
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
            else if (found ==  1 && E[i] >= max_val *1.1){
                max_val = E[i];
                index = i;
            }
        }
    }
    et = rdtsc();
    
    printf ("time to find max energy: %lu\n", (et-st));
    int final = (int)((float)(start + index * step) *(float)(2.0/3.0));
    return final;
}

/* detect_beat
 * Returns the BPM of the given wav file
 */
int detect_beat(char* wav_file, int sample_size) {
 
    unsigned long long st, et;
    st = rdtsc(); 

    //step 1: Load wav file
    float* sample = (float*)malloc(sizeof(float) * sample_size);
    readInWavFile(wav_file, sample, sample_size);

    // Step 2: Differentiate
    float* differentiated_sample = (float*)malloc(sizeof(float) * sample_size);
    derivative(sample, differentiated_sample, sample_size);

    // Step 3: Compute the FFT
    kiss_fft_cpx out[sample_size/2+1];
    fftrArray(differentiated_sample, sample_size, out);
    
    //Step 4: Compute the dot product with comb filters
    int BPM = combfilter_mult(out, sample_size / 2 + 1, sample_size, 15, 200, 1);
  
    free(sample);

    et = rdtsc();

    printf ("total time: %lu\n", (et-st));
    return BPM;
}

/* Main
 * calls BPM and prints result out on commandline
 */
int main(int argc, char* argv[]) {

    //takes in 10 seconds at a sampling rate of 48000 samples/sec (rounded to nearest power of 2)
	int sample_size = 524288;
    int BPM = detect_beat(argv[1],sample_size);
	printf("Final BPM: %i\n", BPM);

	return 0;
}

