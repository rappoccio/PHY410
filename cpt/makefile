# makefile for computational physics toolkit (cpt)

libcpt.a : basalg.o diffeq.o graphics.o linalg.o nonlin.o random.o
	ar rs libcpt.a basalg.o diffeq.o graphics.o linalg.o nonlin.o random.o

basalg.o : basalg.cpp basalg.hpp
	c++ -c basalg.cpp

diffeq.o : diffeq.cpp diffeq.hpp matrix.hpp
	c++ -c diffeq.cpp

graphics.o : graphics.cpp graphics.hpp matrix.hpp
	c++ -c graphics.cpp

linalg.o : linalg.cpp linalg.hpp matrix.hpp
	c++ -c linalg.cpp

nonlin.o : nonlin.cpp nonlin.hpp matrix.hpp
	c++ -c nonlin.cpp

random.o : random.cpp random.hpp
	c++ -c random.cpp
