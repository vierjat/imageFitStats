CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS += 
OBJECTS = imageFitStats.o 
HEADERS = globalConstants.h

ALL : imageFitStats.exe
	echo "Listo!"

imageFitStats.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o imageFitStats.exe $(LIBS) $(GLIBS) $(CFLAGS)

imageFitStats.o : imageFitStats.cc $(HEADERS)
	$(CPP) -c imageFitStats.cc -o imageFitStats.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
