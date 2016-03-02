CFLAGS = -Wall -Wextra -O3 -std=c99 -mfpu=neon -march=native -mtune=native 
CFLAGS += -funroll-loops -mneon-for-64bits

all: test_fft rs bch

COMMON = vector.c vector.h common.c bch_code.c

test_fft: test_fft.c $(COMMON)
	$(CC) $(CFLAGS) $< -o $@

fft_bench: fft_bench.o bch_code.o
	$(CC) $(CFLAGS) bch_code.o fft_bench.o  -o $@


rs : rs.c $(COMMON)
	$(CC) $(CFLAGS) $< -o $@

bch : bch.c $(COMMON)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *.dSYM test_fft rs *.s *.aux *.log *.tex *.pdf
