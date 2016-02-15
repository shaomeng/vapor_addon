CC=gcc
CXX=g++
CXXFLAGS=-O3 -std=c++0x -Wall -DMODELS -m64 -g
LDFLAGS=-m64

ARCH=$(shell uname)

ifeq ($(ARCH), Linux)
VAPOR_INSTALL=/home/users/samuelli/Tools/vapor-git/Install
CXXFLAGS+=-DLINUX -D__USE_LARGEFILE64 -D__USE_LARGEFILE64 -pthread -DLinux
LDFLAGS+=-lrt -pthread -Wl,-rpath,$(VAPOR_INSTALL)/bin -Wl,-rpath,$(VAPOR_INSTALL)/lib
endif

ifeq ($(ARCH), Darwin)
CC=gcc-5
CXX=g++-5
VAPOR_INSTALL=/Users/samuel/Tools/vapor-git/install
CXXFLAGS+=-DDARWIN -DDarwin
LDFLAGS+=-headerpad_max_install_names  -framework CoreFoundation
endif


#ifeq ($(ARCH), Darwin)
#CC=clang
#CXX=clang++
#VAPOR_INSTALL=/Users/samuel/Tools/vapor-git/install
#CXXFLAGS+=-DDARWIN -DDarwin -stdlib=libc++
#LDFLAGS+=-stdlib=libc++ -headerpad_max_install_names  -framework CoreFoundation
#endif



VAPOR_INC=${VAPOR_INSTALL}/include
VAPOR_LIB=${VAPOR_INSTALL}/lib
VAPOR_BIN=${VAPOR_INSTALL}/bin

john_time_comp: john_time_comp.cpp
	$(CXX) john_time_comp.cpp -o bin/john_time_comp.o $(CXXFLAGS) -I${VAPOR_INC} -I. 
	$(CXX) bin/john_time_comp.o -o bin/john_time_comp $(LDFLAGS) -L${VAPOR_BIN} -L$(VAPOR_LIB) -lwasp -lcommon 

cube3d.o: cube3d.cpp cube3d.h
	$(CXX) -c cube3d.cpp -o cube3d.o $(CXXFLAGS) -I${VAPOR_INC} -I. 
#	$(CXX) cube3d.cpp -o bin/cube3d $(CXXFLAGS) $(LDFLAGS) -I${VAPOR_INC} -L${VAPOR_BIN} -L$(VAPOR_LIB) -lwasp -lcommon -I.

slicegroup.o: slicegroup.cpp slicegroup.h
	$(CXX) -c slicegroup.cpp  -o slicegroup.o $(CXXFLAGS) -I${VAPOR_INC}  -I.

temporal: temporal.cpp bin/cube3d.o bin/slicegroup.o 
	$(CXX) -c temporal.cpp -o bin/temporal.o $(CXXFLAGS) -I${VAPOR_INC} -I. 
	$(CXX) bin/temporal.o bin/cube3d.o bin/slicegroup.o -o bin/temporal $(LDFLAGS) -L${VAPOR_BIN} -L$(VAPOR_LIB) -lwasp -lcommon 

wavelet4d.o: wavelet4d.h wavelet4d.cpp
	$(CXX) -c wavelet4d.cpp -o wavelet4d.o $(CXXFLAGS) -I. -I${VAPOR_INC} -fopenmp 
#	$(CXX) wavelet4d.o cube3d.o slicegroup.o -o wavelet4d $(LDFLAGS) -L$(VAPOR_BIN) -L$(VAPOR_LIB) -lwasp -lcommon -fopenmp 

libwavelet4d.a: wavelet4d.o cube3d.o slicegroup.o
	ar -rsv libwavelet4d.a cube3d.o slicegroup.o wavelet4d.o

testwavelet4d: testwavelet4d.cpp libwavelet4d.a
	$(CXX) $(CXXFLAGS) -I. -I${VAPOR_INC} testwavelet4d.cpp -o testwave -L. -lwavelet4d $(LDFLAGS) -L$(VAPOR_BIN) -L$(VAPOR_LIB) -lwasp -lcommon -fopenmp



clean:
	rm bin/*.o bin/john_time_comp bin/temporal
