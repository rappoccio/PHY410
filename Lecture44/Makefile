# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CXX=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CXXFLAGS= -I/usr/X11/include/ -I/usr/local/include
LDFLAGS= -L/usr/X11/lib/ -L/usr/local/lib -lGL -lGLU -lglut -lm


all: hopfield

vmc: hopfield.cpp
	$(CXX) $^  $(CXXFLAGS) $(LDFLAGS) -o hopfield

clean:
	rm -rf *o hopfield 
