# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-I$(CPT_PATH) -L$(CPT_PATH) -lcpt 


all: random_test quad rwalk

random_test: random_test.cpp
	$(CC) $^  $(LIBS) -o random_test


quad: quad.cpp
	$(CC) $^  $(LIBS) -o quad


rwalk: rwalk.cpp
	$(CC) $^  $(LIBS) -o rwalk



clean:
	rm -rf *o random_test quad rwalk

