CFLAGS = -mavx2 -Wall -Wextra -g -O3 -std=c99 -march=native
CFLAGS += -funroll-loops

all: test_fft32 fft_bench32 fft_bench test_juggling

fft_bench: fft_bench.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

test_fft: test_fft.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

test_juggling: test_juggling.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

test_juggling16: test_juggling16.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

test_fft32: test_fft32.c vector.h vector32.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@


fft_bench.s: fft_bench.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -S $< -o $@


new_spring: new_spring.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *.dSYM test_fft32 fft_bench32 fft_bench *.s *.aux *.log *.tex *.pdf