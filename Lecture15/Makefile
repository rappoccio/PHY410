# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-I$(CPT_PATH) -L$(CPT_PATH) -lcpt


all: triatomic matrixmult

triatomic: triatomic.cpp
	$(CC) $^  $(LIBS) -o triatomic

matrixmult: matrixmult.cpp
	$(CC) $^  $(LIBS) -o matrixmult


clean:
	rm -rf *o triatomic matrixmult
