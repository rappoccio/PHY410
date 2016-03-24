# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-I$(CPT_PATH) -L$(CPT_PATH) -lcpt


all: poisson_mg

poisson_mg: poisson_mg.cpp
	$(CC) $^  $(LIBS) -o poisson_mg



clean:
	rm -rf *o poisson_mg

