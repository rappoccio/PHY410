# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-I$(INCLUDEPATH) -L$(LD_LIBRARY_PATH) -lcpt


all: co2_polyfit bvp_example

co2_polyfit: co2_polyfit.cpp
	$(CC) $^  $(LIBS) -o co2_polyfit

bvp_example: bvp_example.cpp
	$(CC) $^  $(LIBS) -o bvp_example


clean:
	rm -rf *o co2_polyfit bvp_example
