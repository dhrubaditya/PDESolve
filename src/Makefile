cc=g++ -c
link_fft = -lfftw3 -L/Users/dhruba/Lib
include_fft = -I/Users/dhruba/Lib/include
default:
	$(cc) FFT.cc $(include_fft)
	$(cc) test.c $(include_fft)
	g++ test.o FFT.o -lfftw3 -lm -L/Users/dhruba/Lib/lib
code:
	$(cc) hydro.cc
	$(cc) pde.cc
	$(cc) evolve.cc
