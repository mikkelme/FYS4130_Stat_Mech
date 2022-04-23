.PHONY: all clean


CPPflags = c++  -O3 -std=c++11
LIB =  #-larmadillo #-llapack -lblas


PROG = main
all: main main.o system.o clean run


main: 						main.o system.o
										${CPPflags} main.o system.o ${LIB} -o main.out


main.o: 					main.cpp
										${CPPflags} -c main.cpp

system.o: 				system.cpp system.h
										${CPPflags} -c system.cpp


run: 
	./main.out


clean:
	rm -f *.o *~




