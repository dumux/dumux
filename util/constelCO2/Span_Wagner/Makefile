# Makefile f"ur C-Script Pointer-Beispiele

SHELL   = /bin/sh
CC      = gcc
CXX     = g++
G77     = gfortran
LDFLAGS = -lc -lm -lgfortran
OBJECTS = extractproperties.o \
	  carb_dio.o \
	  solubilities.o

extractproperties: $(OBJECTS)
	$(CXX)	$(CFLAGS) -o extractproperties $(OBJECTS) $(LDFLAGS) 


all: $(OBJECTS)

.cc.o:
	$(CXX) -O2 -c  $<

.c.o:
	$(CC) -O2 -c  $<

.f.o:
	$(G77) -O2 -c  $<

clean:
	rm -f *.o core a.out *.*~ extractproperties
