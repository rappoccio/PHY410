# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall
LDFLAGS=-I$(CPT_PATH) -L$(CPT_PATH) -lcpt

all: kronig_penney schroedinger

$(SRC): $(SRC).o
	$(CC) $(LDFLAGS) $(SRC).o -o $(SRC)

$(SRC).o: $(SRC).cpp
	$(CC) $(CFLAGS)  $(LDFLAGS) $(SRC).cpp

clean:
	rm -rf *o kronig_penney schroedinger *.pyc *.data
