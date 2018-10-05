
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
   
    //printf("%s\n","start" );

    if (infile)
    {
       
        fseek(infile,44,SEEK_SET);

        while (count < sample_size){
            fread(&buff16,sizeof(buff16),1,infile);        // Reading data in chunks of BUFSIZE
            output[count] = buff16;
            count++;                    
    	}
	/*printf("the first value is %f", output[2000]);   
    printf("the first value is %f", output[3000]);
    printf("the first value is %f", output[4000]);*/
	}
}

void fftrArray(float* sample, int size, kiss_fft_cpx* out) {
    kiss_fftr_cfg cfg;

    int i;

    if ((cfg = kiss_fftr_alloc(size, 0, NULL, NULL)) == NULL) {
        printf("Not enough memory to allocate fftr!\n");
        exit(-1);
    }
    //#pragma omp parallel for
    //for (int i = 0; i < size; i++) {
    //    in[i] = sample[i];
    //}
    unsigned long long st, et;
    st = rdtsc(); 
  kiss_fftr(cfg, (kiss_fft_scalar*)sample, out);
    free(cfg);

    et = rdtsc();
    printf ("time to take FFT: %lu\n", (et-st));

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
    ////#pragma omp parallel for schedule(dynamic)
    for (i = 0; i < energyCount; i++) {
        int BPM = start + i * step;
        int Ti = 60 * 44100/BPM;
        float l[sample_size];
        count = 0;
        //printf("BPM: %d Sample Size: %d Ti: %d\n", BPM, sample_size, Ti);

        for (k = 0; k < sample_size; k++) {
            if ((k % Ti) == 0) {
                count++;
                l[k] = (float)AmpMax;
                /*l[k+1] = (float)AmpMax;*/
            }
            else {
                l[k] = -0.1;
               /* l[k+1] = 0;*/
            }
        }
        //printf("Number of peaks: %d \n", count);

        //kiss_fft_cpx out[sample_size/2+1];

        fftrArray(l, sample_size, out[i]);
    }

    st = rdtsc();
    double sum = 0;
    float a, b;
    double temp;

    for (i = 0; i < energyCount; i++) {
        sum = 0;
        for (k = 0; k < sample_size/2+1; k++) {
            //printf("sample: %f %f \t %d BPM comb: %f %f \t ",fft_array[k].r, fft_array[k].i, BPM, out[k].r, out[k].i);
            a = fft_array[k].r * out[i][k].r - fft_array[k].i * out[i][k].i;
            b = fft_array[k].r * out[i][k].i + fft_array[k].i * out[i][k].r;
            //printf("a: %f b: %f \t ", a ,b);
            sum += sqrt((double)a*a + b*b);
            //printf("Added %f to E[%d], value of %f\n", temp, i, sum);
        }
        //printf("E[%d]: %f\n", i, sum);
        E[i] = sum;
    }
    et = rdtsc();
    printf ("time to take energy: %lu\n", (et-st));

    //Calculate max of E[k]
    double max_val = -1;
    //double max_val_array[];
    int index = -1; 
    int found = 0;
    //int index_array[];

    st = rdtsc();
    //#pragma omp parallel for reduction(max:max_val)
    for (i = 0; i < energyCount; i++) {
       // printf("BPM: %d \t Energy: %f\n", start + i*step, E[i]);
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
    st = rdtsc(); 

    int max_freq = sample_size/4.4;
    // Load mp3
    float* sample = (float*)malloc(sizeof(float) * sample_size);
    //readMP3(s, sample, sample_size);
    readInWavFile(sample, sample_size);

    //for (int i = 0; i < sample_size; i++) {
    //    printf("Element %i: %i\n", i, sample[i]);
    //}

    // Step 2: Differentiate
    float* differentiated_sample = (float*)malloc(sizeof(float) * sample_size);
    int Fs = 44100;

    unsigned long long st_1, et_1;
    st_1 = rdtsc(); 
    differentiated_sample[0] = sample[0];
  
    for (int i = 1; i < sample_size - 1; i++) {
        differentiated_sample[i] = Fs * (sample[i+1]-sample[i-1])/2; 
    }
    differentiated_sample[sample_size - 1] = sample[sample_size-1];
    et_1 = rdtsc();
    printf ("time to take derivative: %lu\n", (et_1-st_1));

    // Step 3: Compute the FFT
    kiss_fft_cpx out[sample_size/2+1];
    fftrArray(differentiated_sample, sample_size, out);
    printf("%s\n", "first FFT");

    //for (int i = 0; i < sample_size / 2; i++)
    //  printf("out[%2zu] = %+f , %+f\n", i, out[i].r, out[i].i);

    //printf("Combfilter performing...\n");

    int BPM = combfilter(out, sample_size / 2 + 1, sample_size, 15, 200, 1);
   //int BPM = combfilter(out, sample_size / 2 + 1, sample_size, BPM-5, BPM+5, 1);

    //printf("Final BPM: %i\n", BPM);

    // Step 4: Generate Sub-band array values
    free(sample);

    et = rdtsc();

    printf ("total time: %lu\n", (et-st));
    return BPM;
}

int main(int argc, char* argv[]) {
  // Test CPU Version
    //takes in 10 seconds at a sampling rate of 48000 samples/sec
	int BPM = detect_beat(480000);
	printf("Final BPM: %i\n", BPM);

	return 0;
}

