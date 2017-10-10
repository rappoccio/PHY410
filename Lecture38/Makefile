# I am a comment, and I want to say that the variable CXX will be
# the compiler to use.
CXX=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CXXFLAGS=-Wall -std=c++11
LDFLAGS=-I$(CPT_PATH) -L$(CPT_PATH) -lcpt

all: dmc pimc 

$(SRC): $(SRC).o
	$(CXX) -c $(CXXFLAGS) $(LDFLAGS) $(SRC).o -o $(SRC)

$(SRC).o: $(SRC).cpp
	$(CXX) $(CXXFLAGS)  $(LDFLAGS) $(SRC).cpp

clean:
	rm -rf *o dmc pimc *.pyc *.data
