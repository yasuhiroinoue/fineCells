
#CC=g++
#LIB=-lm
#OPT=-O2 -fopenmp
CC=g++
LIB=-lm
OPT=-O3 -fopenmp
CPPVER=-std=c++11
#C++11以上の規格でコンパイルする必要がある。

all: a.out

a.out: _parameters.o _class.o ODEsolver.o random.o main.o
	$(CC) $(LIB) $(OPT) $(CPPVER) _parameters.o _class.o ODEsolver.o random.o main.o

main.o: main.cpp ODEsolver.h _class.h _parameters.h random.h
	$(CC) $(LIB) $(OPT) $(CPPVER) -c main.cpp

ODEsolver.o: ODEsolver.cpp ODEsolver.h _class.h _parameters.h
	$(CC) $(LIB) $(OPT) $(CPPVER) -c ODEsolver.cpp

_class.o: _class.cpp _class.h _parameters.h
	$(CC) $(LIB) $(OPT) $(CPPVER) -c _class.cpp

_parameters.o: _parameters.cpp _parameters.h
	$(CC) $(LIB) $(OPT) $(CPPVER) -c _parameters.cpp

random.o : random.cpp random.h
	$(CC) $(LIB) $(OPT) $(CPPVER) -c random.cpp

.PHONY: clean rm

clean:
	rm -f *.o

rm:
	rm -f output/* & rm -f data/*
