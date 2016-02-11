CFLAGS = -Wall -Wextra -g -O3 -std=c99 -mfpu=neon 
CFLAGS += -funroll-loops

all: test_fft

VECTOR = vector.c vector.h

test_fft: test_fft.c $(VECTOR)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *.dSYM test_fft *.s *.aux *.log *.tex *.pdf
