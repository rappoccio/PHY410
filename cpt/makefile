# makefile for computational physics toolkit (cpt)
CC=g++

libcpt.a : basalg.o diffeq.o graphics.o linalg.o polyfit.o nonlin.o random.o
	ar rs libcpt.a basalg.o diffeq.o graphics.o linalg.o polyfit.o nonlin.o random.o

basalg.o : basalg.cpp basalg.hpp basalg.icc
	$(CC) -c basalg.cpp

diffeq.o : diffeq.cpp diffeq.hpp matrix.hpp
	$(CC) -c diffeq.cpp

graphics.o : graphics.cpp graphics.hpp matrix.hpp
	$(CC) -c graphics.cpp

linalg.o : linalg.cpp linalg.hpp matrix.hpp
	$(CC) -c linalg.cpp

polyfit.o : polyfit.cpp polyfit.hpp matrix.hpp
	$(CC) -c polyfit.cpp


nonlin.o : nonlin.cpp nonlin.hpp matrix.hpp
	$(CC) -c nonlin.cpp

random.o : random.cpp random.hpp
	$(CC) -c random.cpp

clean : 
	rm -rf *.o *.a 