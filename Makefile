GCC = gcc
FLAGS = -O2 -fopenmp -std=c99 -m64 -mavx -Wall

.SUFFIXES: .o .c

BPM_SRCS = BeatCalculator.c ./kiss_fft130/kiss_fft.c ./kiss_fft130/tools/kiss_fftr.c 
BPM_EXECUTABLE = beat_calc
BPM_LINKERS = -lm -I ./kiss_fft130/ -I ./kiss_fft130/tools/

all:
	$(GCC) $(FLAGS) $(BPM_SRCS) $(BPM_LINKERS) -o $(BPM_EXECUTABLE)

bpm:
	$(GCC) $(FLAGS) $(BPM_SRCS) $(BPM_LINKERS) -o $(BPM_EXECUTABLE)

clean:
	rm *~ $(BPM_EXECUTABLE)
