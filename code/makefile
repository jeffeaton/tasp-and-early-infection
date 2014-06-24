## Intel compiler
# CC = icc 
# CFLAGS = -lgsl -lgslcblas -openmp -O2 -std=c++11 -D fullmod=0 -I $(INCLUDE)

## GCC compiler
CC = g++-4.8  # (OpenMP support for gcc >= 4.2)
CFLAGS = -lgsl -lgslcblas -fopenmp -O2 -std=c++11 -D fullmod=0

## Note: OpenMP not supported on LLVM 5.0 compiler (default in Mac OS 10.8)
## Note: may need to add -I <path to GSL headers> and -L <path to GSL libraries> if not in default search path

all: fitpost fitfixed

fitpost:
	$(CC) $(CFLAGS) runimis.cpp imis.cpp likelihood.cpp simprev.cpp model.cpp states.cpp -o runimis

fitfixed:
	$(CC) $(CFLAGS) runimis-fixed-relinfect.cpp imis.cpp likelihood-fixed-relinfect.cpp simprev.cpp model.cpp states.cpp -o runimis-fixed-relinfect
