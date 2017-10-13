# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-std=c++11 -I$(CPT_PATH) -L$(CPT_PATH) -lcpt


all: dmc pimc

dmc: dmc.cpp
	$(CC) $^  $(LIBS) -o dmc
pimc: pimc.cpp
	$(CC) $^  $(LIBS) -o pimc

clean:
	rm -rf *o dmc pimc *.pyc *.data *.out
