# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-I$(INCLUDEPATH) -L$(LD_LIBRARY_PATH) -lcpt


all: ising

ising: ising.cpp
	$(CC) $^  $(LIBS) -o ising



clean:
	rm -rf *o ising


