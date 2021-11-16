#!/bin/sh -x
CC=${CC:-clang}

$CC -DHAVE_AVX2 -g -O3 -mavx2 *.c -lm -pthread -o dr4j
strip dr4j

$CC -DHAVE_AVX2 -DDECODE_ONLY -g -O3 -mavx2 *.c -lm -pthread -o dr4j_dec
strip dr4j_dec

