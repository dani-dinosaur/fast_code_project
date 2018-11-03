
#include <stdio.h>
#include <math.h>
#include <kiss_fftr.h>
#include <kiss_fft.h>


static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}


void readInWavFile(float* output, int sample_size){
	//Open wave file in read mode
	FILE * infile = fopen("100bpm.wav","rb");
    // For counting number of frames in wave file.
    int count = 0;                        
    /// short int used for 16 bit as input data format is 16 bit PCM audio
    short int buff16;
    
    if (infile) {
        fseek(infile,44,SEEK_SET);

        while (count < sample_size) {
            fread(&buff16,sizeof(buff16),1,infile);        // Reading data in chunks of BUFSIZE
            output[count] = buff16;
            count++;                    
    	}
	}
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
    free(cfg);

    et = rdtsc();
    printf ("time to take FFT: %llu\n", (et-st));

}

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

int combfilter(kiss_fft_cpx fft_array[], int size, int sample_size, int start, int fin, int step) {
    float AmpMax =0.0001;
    int energyCount = (fin - start)/step;
    double E[energyCount];
    unsigned long long st, et;
    int count = 0;
    int i, k;
    kiss_fft_cpx *out[energyCount];
    for (int i = 0; i < energyCount; i++){
        out[i] = (kiss_fft_cpx*)malloc((sample_size/2+1)*sizeof(kiss_fft_cpx));
    }

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
                l[k] = -0.1;
            }
        }
        fftrArray(l, sample_size, out[i]);
    }

    st = rdtsc();
    double sum = 0;
    float a, b;

    for (i = 0; i < energyCount; i++) {
        sum = 0;
        for (k = 0; k < sample_size/2+1; k++) {
            a = fft_array[k].r * out[i][k].r - fft_array[k].i * out[i][k].i;
            b = fft_array[k].r * out[i][k].i + fft_array[k].i * out[i][k].r;
            sum += sqrt((double)a*a + b*b);
        }
        E[i] = sum;
    }
    et = rdtsc();
    printf ("time to take energy: %llu\n", (et-st));

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
    printf ("time to find max energy: %llu\n", (et-st));
    return start + index * step;
}

/* detect_beat
 * Returns the BPM of the given mp3 file
 * @Params: s - the path to the desired mp3
 */
int detect_beat(int sample_size) {
    // Step 1: Get a 5-second sample of our desired mp3
    // Assume the max frequency is 4096
    unsigned long long st, et;
    unsigned long long st_1, et_1;

    st = rdtsc(); 

    // Load mp3
    float* sample = (float*)malloc(sizeof(float) * sample_size);
    readInWavFile(sample, sample_size);

    // Step 2: Differentiate
    float* differentiated_sample = (float*)malloc(sizeof(float) * sample_size);
    int Fs = 44100;

    st_1 = rdtsc(); 
    differentiated_sample[0] = sample[0];
    for (int i = 1; i < sample_size - 1; i++) {
        differentiated_sample[i] = Fs * (sample[i+1]-sample[i-1])/2; 
    }
    differentiated_sample[sample_size - 1] = sample[sample_size-1];
    et_1 = rdtsc();

    printf ("time to take derivative: %llu\n", (et_1-st_1));

    // Step 3: Compute the FFT
    kiss_fft_cpx out[sample_size/2+1];
    fftrArray(differentiated_sample, sample_size, out);
    printf("%s\n", "first FFT");

    // Step 4: Apply Comb Filter
    int BPM = combfilter(out, sample_size / 2 + 1, sample_size, 15, 200, 1);


    // Step 5: Free values allocated
    free(sample);
    free(differentiated_sample);

    et = rdtsc();

    printf ("total time: %llu\n", (et-st));
    return BPM;
}

int main(int argc, char* argv[]) {
    // Test CPU Version
    //takes in 10 seconds at a sampling rate of 48000 samples/sec
	int BPM = detect_beat(480000);
	printf("Final BPM: %i\n", BPM);

	return 0;
}

