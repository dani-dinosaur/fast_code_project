############################To Run Code####################################
To run our BPM calculator with the optimized Derivative and Dot Product:
make bpm
./beat_calc ./pure_bpm/100bpm.wav

To run our optimized FFT code:
make fft_new
./fft_new ./pure_bpm/100bpm.wav 16 

To run the baseline FFT code (for comparision):
make fft_start
./fft_start ./pure_bpm/100bpm.wav 16 

To run the Kiss FFT code (for comparision):
make fft_start
./fft_kiss ./pure_bpm/100bpm.wav 16 

#############################Setup Instructions############################
Source code explaination:
The file BeatCalculator.c contains the optimized derivative and dot product. The fft used for this file is the kissFFT implementation.
The file BeatCalculator_new.c corresponds to the new, optimized FFT implementation. 
The file BeatCalculator_start.c corresponds to the baseline FFT implementation. 
The file BeatCalculator_kiss.c corresponds to the Kiss FFT (start of the art) implementation.

Array of Comb Filters:

Google Drive link sent to Instructor and TAs (combfilter.data).
This file contains the precomputed comb filters used to perform the BPM analysis. 

Setup Instructions:

Pull code from Github
Download data file from Google Drive
Place combfilters.data in file with BeatCalculator.c
README file contains commands to run each implementation
