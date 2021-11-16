#!/bin/sh -x
#cc -DHAVE_AVX2 -g -O3 -mavx2 *.c -lm -pthread -o da9j -falign-functions=32
clang -DHAVE_AVX2 -g -O3 -mavx2 *.c -lm -pthread -o da9j
strip da9j

#cc -DHAVE_AVX2 -DDECODE_ONLY -g -O3 -mavx2 *.c -lm -pthread -o da9j_dec
#strip da9j_dec

