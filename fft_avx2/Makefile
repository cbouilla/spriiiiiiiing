CFLAGS = -Wall -Wextra -O3 -std=c99 -mtune=native -mpclmul -msse4 
CFLAGS += -funroll-loops
NEW = -mavx2 -mtune=native -march=native

all: new_spring test_fft32 bch_ctr

#common.o: common.c vector.h vector32.c
	#$(CC) -c $(NEW) $(CFLAGS) $< -o $@

new_spring: new_spring.o 
bch_ctr: bch_ctr.o 

new_spring.o: new_spring.c vector.h vector32.c
	$(CC) -c $(NEW) $(CFLAGS) $< -o $@

bch_ctr.o: bch_ctr.c vector.h vector32.c
	$(CC) -c $(NEW) $(CFLAGS) $< -o $@

fft_bench: fft_bench.c vector.h vector.c
	$(CC) $(CFLAGS) $< -o $@

fft_bench32: fft_bench32.c vector.h vector32.c
	$(CC) $(NEW) $(CFLAGS) $< -o $@

test_fft: test_fft.c vector.h vector.c
	$(CC) $(CFLAGS) $< -o $@

test_fft32: test_fft32.c vector.h vector32.c
	$(CC) $(NEW)  $(CFLAGS) $< -o $@

fft_bench.s: fft_bench.c vector.h vector.c
	$(CC) $(NEW) $(CFLAGS) -S $< -o $@

clean:
	rm -rf *.o *.dSYM new_spring bch_ctr test_fft32 fft_bench32 fft_bench *.s *.aux *.log *.tex *.pdf