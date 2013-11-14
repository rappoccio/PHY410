# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-I$(INCLUDEPATH) -L$(LD_LIBRARY_PATH) -lcpt


all: poisson_fft poisson

poisson: poisson.cpp
	$(CC) $^  $(LIBS) -o poisson

poisson_fft: poisson_fft.cpp
	$(CC) $^  $(LIBS) -o poisson_fft


clean:
	rm -rf *o poisson poisson_fft

