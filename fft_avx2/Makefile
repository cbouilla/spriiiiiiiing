CFLAGS = -Wall -Wextra -O3 -std=c99 -g -mtune=native -mpclmul -msse4 
CFLAGS += -funroll-loops
NEW = -mavx2 -mtune=native -march=native

all: new_ctr new_ofb bch_ctr bch_ofb test_fft32 fft_bench32

#common.o: common.c vector.h vector32.c
	#$(CC) -c $(NEW) $(CFLAGS) $< -o $@

#new_spring: new_spring.o 
#bch_ctr: bch_ctr.o 

COMMON = vector.h vector32.c common.c

new_ctr.o: new_ctr.c $(COMMON)
	$(CC) -c $(NEW) $(CFLAGS) $< -o $@

new_ofb.o: new_ofb.c $(COMMON)
	$(CC) -c $(NEW) $(CFLAGS) $< -o $@

bch_ctr.o: bch_ctr.c $(COMMON)
	$(CC) -c $(NEW) $(CFLAGS) $< -o $@

bch_ofb.o: bch_ofb.c $(COMMON)
	$(CC) -c $(NEW) $(CFLAGS) $< -o $@

fft_bench: fft_bench.c $(COMMON)
	$(CC) $(CFLAGS) $< -o $@

fft_bench32: fft_bench32.c $(COMMON)
	$(CC) $(NEW) $(CFLAGS) $< -o $@

test_fft: test_fft.c vector.h vector.c
	$(CC) $(CFLAGS) $< -o $@

test_fft32: test_fft32.c $(COMMON)
	$(CC) $(NEW)  $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *.dSYM new_ctr new_ofb bch_ctr bch_ofb test_fft32 fft_bench32 *.s *.aux *.log *.tex *.pdf