OPTICC=-O3 -w -std=c++11 -axSSE4.2 -ipo
OPTGCC   = -std=c++11 -O3 -w -mtune=corei7 -msse4.2 -fprefetch-loop-arrays -funroll-all-loops -ffast-math
OPT0     = -std=c++11 -O0

LIBDIR=$(HOME)/LIBRARY/

all: INFO

INFO:
	@echo
	@echo "Compilation options are:"
	@echo
	@echo "     INT    : SSE4.2 Optimized Intel compiler ipo executable"
	@echo "     GCC    : SSE4.2 Optimized GCC compiler executable"
	@echo
	@echo "     INT0   : Non-Optimized Intel compiler executable "
	@echo "     GCC0   : Non-Optimized GCC compiler executable"
	@echo


INT:	main.cpp
	icpc $(OPTICC) main.cpp -DINTEL=1 -mkl=sequential -o brok-new.run

INT0:	main.cpp
	icpc $(OPT0) main.cpp -DINTEL=1 -mkl=sequential -lm -o brok-new.run

GCC:	main.cpp
	g++ $(OPTGCC) main.cpp -DINTEL=0 -lboost_filesystem -lboost_system -lfftw3 -lm -o brok-new.run

GCC0:	main.cpp
	g++ $(OPT0) main.cpp -DINTEL=0 -lboost_filesystem -lboost_system -lfftw3 -lm -o brok-new.run

CON:	concatenate.cpp
	g++ $(OPTGCC) concatenate.cpp -o DataCon.run

clean:
	@rm -f *.run *~ *backup*
	@rm -rf results
