CFLAGS = -Wall -Wextra -g -O3 -std=c99 -march=native
CFLAGS += -funroll-loops -fpeel-loops --param max-completely-peeled-insns=850

all: new_spring

new_spring: new_spring.c poly_eval.h vector.h vector.c compat.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@


spring_bch: spring_bch.c spring.h rand_generated_arrays.c vector.h vector.c compat.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

spring_crt: spring_crt.c spring.h rand_generated_arrays.c vector.h vector.c compat.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@

spring_crt.s: spring_crt.c spring.h rand_generated_arrays.c vector.h vector.c compat.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ -S -fverbose-asm