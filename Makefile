GCC = gcc
FLAGS = -std=c99 -m64 -mavx -Wall -g

.SUFFIXES: .o .c

BPM_SRCS = BeatCalculator.c ./kiss_fft130/kiss_fft.c ./kiss_fft130/tools/kiss_fftr.c 
BPM_EXECUTABLE = beat_calc
BPM_LINKERS = -lm -I ./kiss_fft130/ -I ./kiss_fft130/tools/

CFFT_SRCS = BeatCalculator.c ./cfft/fft.c ./cfft/complex.c wave.c
CFFT_EXECUTABLE = beat_calc
CFFT_LINKERS = -lm -I ./cfft/

all:
	$(GCC) $(FLAGS) $(CFFT_SRCS) $(CFFT_LINKERS) -o $(CFFT_EXECUTABLE)

cfft:
	$(GCC) $(FLAGS) $(CFFT_SRCS) $(CFFT_LINKERS) -o $(CFFT_EXECUTABLE)

bpm:
	$(GCC) $(FLAGS) $(BPM_SRCS) $(BPM_LINKERS) -o $(BPM_EXECUTABLE)

clean:
	rm *~ $(BPM_EXECUTABLE) $(CFFT_EXECUTABLE)
