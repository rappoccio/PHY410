ALL:    build

lib: 
	@echo "Building library..."
	cd cpt; make


build: lib
	@echo "Building system..."
	cd Lecture14; make
	cd Lecture15; make
	cd Lecture17; make
	cd Lecture18; make
	cd Lecture20; make
	cd Lecture21; make
	cd Lecture22; make
	cd Lecture24; make
	cd Lecture25; make
	cd Lecture26; make
	cd Lecture27; make
	cd Lecture28; make
	cd Lecture29; make
	cd Lecture30; make
	cd Lecture31; make
	cd Lecture32; make

clean: 
	@echo "Cleaing system..."
	cd cpt; make clean
	cd Lecture14; make clean
	cd Lecture15; make clean
	cd Lecture17; make clean
	cd Lecture18; make clean
	cd Lecture20; make clean
	cd Lecture21; make clean
	cd Lecture22; make clean
	cd Lecture24; make clean
	cd Lecture25; make clean
	cd Lecture26; make clean
	cd Lecture27; make clean
	cd Lecture28; make clean
	cd Lecture29; make clean
	cd Lecture30; make clean
	cd Lecture31; make clean
	cd Lecture32; make clean
