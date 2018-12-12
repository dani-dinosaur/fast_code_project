To run the BPM calculator with the optimized Derivative and Dot Product:
	make bpm
	./beat_calc ./pure_bpm/100bpm.wav

To run the optimized FFT code:
	make fft_new
	./fft_new ./pure_bpm/100bpm.wav 16 

To run the baseline FFT code:
	make fft_start
	./fft_start ./pure_bpm/100bpm.wav 16 

To run the Kiss FFT code:
	make fft_start
	./fft_kiss ./pure_bpm/100bpm.wav 16 

