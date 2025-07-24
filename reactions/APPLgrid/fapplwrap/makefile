
# ROOTC=$(shell root-config --cflags)
# ROOTL=$(shell root-config --libs)

APPLC=$(shell applgrid-config --cxxflags)
APPLL=$(shell applgrid-config --ldflags)
APPLFL=$(shell applgrid-config --ldfflags)

LHAC=$(shell lhapdf-config --cflags)
LHAL=$(shell lhapdf-config --libs)

XTRA=-I$(HOME)/nlo/appl-conv

UNAME=$(shell uname)

ifeq ($(UNAME),Linux)
    CSTDLIB=-lstdc++
else ifeq ($(UNAME),Darwin)
    CSTDLIB=-lc++
else
    # this will cause the compilation to fail with unresolved tokens
    CSTDLIB=
endif

install : libs
	cp -v libfapplwrap.a libapplwrap.a $(HOME)/lib
	cp -v *.h  $(HOME)/include


libs : libfapplwrap.a libapplwrap.a

%.o : %.cxx 
	g++ -O2 -c $< -o $@ $(APPLC) $(LHAC) $(XTRA)

%.o : %.f90
	gfortran -O2 -c $< -o $@ $(APPLC) $(LHAC) $(XTRA)

libfapplwrap.a : applwrap.o smooth.o fapplwrap.o
	ar -r $@  $^

libapplwrap.a : applwrap.o smooth.o 
	ar -r $@  $^


fpdf : fpdf.f90 libs
	gfortran -O2 $< -o $@ $(APPLFL) $(LHAL) -L$(HOME)/lib -lfapplwrap $(CSTDLIB)

dist:
	( cd .. ; tar -czf fapplwrap/fapplwrap.tgz fapplwrap/applwrap.cxx fapplwrap/applwrap.h fapplwrap/smooth.h fapplwrap/smooth.cxx fapplwrap/fapplwrap.cxx fapplwrap/fpdf.f90 fapplwrap/makefile ; pwd )

clean :
	rm *.o *.a

