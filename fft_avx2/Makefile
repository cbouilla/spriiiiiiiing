CFLAGS = -Wall -Wextra -O3 -std=c99 -mtune=native -mpclmul
CFLAGS += -funroll-loops
OLD = -msse4
NEW = -mavx2 -mtune=native -march=native

all: spring test_fft32 bch_ctr

spring: new_spring.c vector.h vector32.c
	$(CC) $(NEW) $(CFLAGS) $< -o $@

bch_ctr: bch_ctr.c vector.h vector32.c
	$(CC) $(NEW) $(CFLAGS) $< -o $@


fft_bench: fft_bench.c vector.h vector.c
	$(CC) $(OLD) $(CFLAGS) $< -o $@

fft_bench32: fft_bench32.c vector.h vector32.c
	$(CC) $(NEW) $(CFLAGS) $< -o $@

test_fft: test_fft.c vector.h vector.c
	$(CC) $(OLD) $(CFLAGS) $< -o $@

test_fft32: test_fft32.c vector.h vector32.c
	$(CC) $(NEW)  $(CFLAGS) $< -o $@

fft_bench.s: fft_bench.c vector.h vector.c
	$(CC) $(NEW) $(CFLAGS) -S $< -o $@

clean:
	rm -rf *.o *.dSYM test_fft32 fft_bench32 fft_bench *.s *.aux *.log *.tex *.pdf