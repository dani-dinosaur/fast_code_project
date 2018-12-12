GCC = gcc
FLAGS = -O2 -std=c99 -m64 -mavx -Wall

.SUFFIXES: .o .c

BPM_EXECUTABLE = beat_calc
FFT_KISS_EXECUTABLE = fft_kiss
FFT_START_EXECUTABLE = fft_start
FFT_NEW_EXECUTABLE = fft_new

BPM_SRCS = BeatCalculator.c ./kiss_fft/kiss_fft.c ./kiss_fft/tools/kiss_fftr.c ./cfft/complex.c ./cfft/fft.c 
BPM_LINKERS = -lm -I ./kiss_fft/ -I ./kiss_fft/tools/ -I ./cfft

FFT_KISS_SRCS = BeatCalculator_kiss.c ./kiss_fft/kiss_fft.c ./kiss_fft/tools/kiss_fftr.c 
FFT_KISS_LINKERS = -lm -I ./kiss_fft/ -I ./kiss_fft/tools/

FFT_START_SRCS = BeatCalculator_start.c ./start_fft/fft.c
FFT_START_LINKERS = -lm -I ./start_fft/

FFT_NEW_SRCS = BeatCalculator_new.c ./new_fft/fft.c
FFT_NEW_LINKERS = -lm -I ./new_fft/

all: bpm fft_kiss fft_start fft_new

bpm:
	$(GCC) $(FLAGS) $(BPM_SRCS) $(BPM_LINKERS) -o $(BPM_EXECUTABLE)

fft_kiss:
	$(GCC) $(FLAGS) $(FFT_KISS_SRCS) $(FFT_KISS_LINKERS) -o $(FFT_KISS_EXECUTABLE)

fft_start:
	$(GCC) $(FLAGS) $(FFT_START_SRCS) $(FFT_START_LINKERS) -o $(FFT_START_EXECUTABLE)

fft_new:
	$(GCC) $(FLAGS) $(FFT_NEW_SRCS) $(FFT_NEW_LINKERS) -o $(FFT_NEW_EXECUTABLE)

clean:
	rm *~ $(BPM_EXECUTABLE) $(FFT_KISS_EXECUTABLE) $(FFT_START_EXECUTABLE) $(FFT_NEW_EXECUTABLE)
