# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
LIBS=-I$(INCLUDEPATH) -L$(LD_LIBRARY_PATH) -lcpt


all: wavepacket_fft wavepacket

wavepacket_fft: wavepacket_fft.cpp
	$(CC) $^  $(LIBS) -o wavepacket_fft

wavepacket: wavepacket.cpp
	$(CC) $^  $(LIBS) -o wavepacket


clean:
	rm -rf *o wavepacket_fft wavepacket

