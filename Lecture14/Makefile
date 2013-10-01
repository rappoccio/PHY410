# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall
LDFLAGS=-I$(INCLUDEPATH) -L$(LD_LIBRARY_PATH) -lcpt

all: co2_polyfit

$(SRC): $(SRC).o
	$(CC) $(LDFLAGS) $(SRC).o -o $(SRC)

$(SRC).o: $(SRC).cpp
	$(CC) $(CFLAGS)  $(LDFLAGS) $(SRC).cpp

clean:
	rm -rf *o co2_polyfit
