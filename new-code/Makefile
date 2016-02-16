CFLAGS = -Wall -Wextra -g -O3 -std=c99 -march=native
CFLAGS += -funroll-loops -fpeel-loops --param max-completely-peeled-insns=850

new_spring: new_spring.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

test_fft: test_fft.c vector.h vector.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@
