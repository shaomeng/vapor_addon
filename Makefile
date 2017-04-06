ARCH=$(shell uname)
USER=$(shell whoami)

CXXFLAGS=-O2 -std=c++11 -Wall -DMODELS -m64 -pthread
LDFLAGS=


ifeq ($(ARCH), Linux)
ifeq ($(USER), shaomeng)  #NCAR
VAPOR_INSTALL=/glade/u/home/shaomeng/lib_gnu/vapor-git
else											#Alaska
VAPOR_INSTALL=/home/users/samuelli/Install/vapor-git
endif
CC=gcc
CXX=g++
VAPOR_DEP_LIB=/glade/p/DASG/VAPOR/third-party/apps-2014/Linux_x86_64
CXXFLAGS+=-DLINUX -D__USE_LARGEFILE64 -D__USE_LARGEFILE64 -pthread -DLinux -DOPENMP -fopenmp
LDFLAGS+=-lrt -pthread -fopenmp -Wl,-rpath,$(VAPOR_INSTALL)/lib -Wl,-rpath,$(VAPOR_DEP_LIB)/lib

else ifeq ($(ARCH), Darwin)
CC=clang
CXX=clang++
VAPOR_DEP_LIB=/glade/p/DASG/VAPOR/third-party/apps-2017/Darwin_x86_64
VAPOR_INSTALL=/Users/shaomeng/Install/vapor-git
endif

john_time_comp: john_time_comp.cpp
	$(CXX) -c john_time_comp.cpp -o bin/john_time_comp.o $(CXXFLAGS) -I${VAPOR_INSTALL}/include -I. 
	$(CXX) bin/john_time_comp.o -o bin/john_time_comp $(LDFLAGS) -L$(VAPOR_INSTALL)/lib -lvdf -lcommon -lexpat -L$(VAPOR_DEP_LIB)/lib -ludunits2 -lproj -lnetcdf -lvdf

cube3d.o: cube3d.cpp cube3d.h
	$(CXX) -c cube3d.cpp -o cube3d.o $(CXXFLAGS) -I${VAPOR_INC} -I. 
#	$(CXX) cube3d.cpp -o bin/cube3d $(CXXFLAGS) $(LDFLAGS) -I${VAPOR_INC} -L${VAPOR_BIN} -L$(VAPOR_LIB) -lvdf -lcommon -I.

slicegroup.o: slicegroup.cpp slicegroup.h
	$(CXX) -c slicegroup.cpp  -o slicegroup.o $(CXXFLAGS) -I${VAPOR_INC}  -I.

temporal: temporal.cpp bin/cube3d.o bin/slicegroup.o 
	$(CXX) -c temporal.cpp -o bin/temporal.o $(CXXFLAGS) -I${VAPOR_INC} -I. 
	$(CXX) bin/temporal.o bin/cube3d.o bin/slicegroup.o -o bin/temporal $(LDFLAGS) -L${VAPOR_BIN} -L$(VAPOR_LIB) -lvdf -lcommon 

wavelet4d.o: wavelet4d.h wavelet4d.cpp
	$(CXX) -c wavelet4d.cpp -o wavelet4d.o $(CXXFLAGS) -I. -I${VAPOR_INC} -fopenmp 
#	$(CXX) wavelet4d.o cube3d.o slicegroup.o -o wavelet4d $(LDFLAGS) -L$(VAPOR_BIN) -L$(VAPOR_LIB) -lvdf -lcommon -fopenmp 

libwavelet4d.a: wavelet4d.o cube3d.o slicegroup.o
	ar -rsv libwavelet4d.a cube3d.o slicegroup.o wavelet4d.o

testwavelet4d: testwavelet4d.cpp libwavelet4d.a
	$(CXX) $(CXXFLAGS) -I. -I${VAPOR_INC} testwavelet4d.cpp -o testwave -L. -lwavelet4d $(LDFLAGS) -L$(VAPOR_BIN) -L$(VAPOR_LIB) -lvdf -lcommon -fopenmp



clean:
	rm bin/*.o bin/john_time_comp bin/temporal
