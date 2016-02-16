CFLAGS = -Wall -Wextra -g -O3 -std=c99 -mfpu=neon 
CFLAGS += -funroll-loops

all: test_fft

COMMON = vector.c vector.h common.c

test_fft: test_fft.c $(COMMON)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *.dSYM test_fft *.s *.aux *.log *.tex *.pdf
