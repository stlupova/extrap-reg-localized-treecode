
# -Wno-c++11-extensions —— turns off warnings about recent c++ features

CPP = clang++
# First line below for serial runs; second line for parallel runs on a Mac
#CPPFLAGS = -c -O2 -Wno-c++11-extensions
CPPFLAGS = -c -O2 -Wno-c++11-extensions -Xclang -fopenmp

main: main.o stokeslets.o KITC.o surface.o utilities.o 
	$(CPP) main.o stokeslets.o KITC.o surface.o utilities.o -lomp -o main.out
	rm *.o

main.o: main.cpp
	$(CPP) $(CPPFLAGS) main.cpp


stokeslets.o: stokeslets.cpp stokeslets.h
	$(CPP) $(CPPFLAGS) stokeslets.cpp


KITC.o: KITC.cpp KITC.h
	$(CPP) $(CPPFLAGS) KITC.cpp


surface.o: surface.cpp surface.h
	$(CPP) $(CPPFLAGS) surface.cpp

utilities.o: utilities.cpp utilities.h
	$(CPP) $(CPPFLAGS) utilities.cpp

clean:
	rm *.o *~

#when name is different from makefile:
#make -f makefile main
