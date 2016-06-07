CPPC     = g++
#CPPC = g++-mp-5
#CPPC     = clang++
FLAGS   += -fopenmp
#FLAGS   += -mavx
FAST    += -march=native #-mfpmath=avx
FAST     = -ffast-math -funroll-loops -O3
WARNINGS = -Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith #-Weffc++
CPPSRCS  = $(wildcard *.cpp)
CPPHDRS  = $(wildcard *.h)
PROGRAM  = SPH

NO_CLR    = \x1b[0m
OK_CLR    = \x1b[32;01m
WARN_CLR  = \x1b[33;01m
ERROR_CLR = \x1b[31;01m

.PHONY:	all clean install

disk:	INTERNAL_MBG_M97_B95_MM97 INTERNAL_MBG_vNR_NO_NO INTERNAL_SPH_M97_CD10_CD10 INTERNAL_SPH_vNR_B95_MM97 INTERNAL_MBG_vNR_B95_MM97 INTERNAL_SPH_M97_B95_MM97 INTERNAL_SPH_M97_NO_NO INTERNAL_SPH_vNR_NO_NO
shock:	INTERNAL_MBG_vNR_NO_NO INTERNAL_SPH_vNR_NO_NO INTERNAL_MBG_vNR_NO_MM97 INTERNAL_SPH_vNR_NO_MM97 INTERNAL_MBG_M97_NO_MM97 INTERNAL_SPH_M97_NO_MM97 INTERNAL_SPH_M97_NO_NO
	
all:	test

include list.txt

clean:
	-rm *.o *.s *.out
	-find ./result/ -name '*.txt' -delete
	-rm *.bin
	-rm -r INTERNAL_*

print:
	cat $(CPPSRCS) $(CPPHDRS) > print.out

