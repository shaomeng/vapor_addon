CC=gcc
CXX=g++
CXXFLAGS=-O2 -c -std=c++0x -Wall -DMODELS -fPIC -m64 -g
LDFLAGS=-m64 -lwasp -lcommon 

ARCH=$(shell uname)

ifeq ($(ARCH), Linux)
VAPOR_INSTALL=/home/users/samuelli/Tools/vapor-git/install
CXXFLAGS+=-DLINUX -D__USE_LARGEFILE64 -D__USE_LARGEFILE64 -pthread -DLinux
LDFLAGS+=-m64 -lrt -pthread -Wl,-rpath
endif

ifeq ($(ARCH), Darwin)
CC=clang
CXX=clang++
VAPOR_INSTALL=/Users/samuel/Tools/vapor-git/install
CXXFLAGS+=-DDARWIN -DDarwin -stdlib=libc++
LDFLAGS+=-stdlib=libc++ -headerpad_max_install_names  -framework CoreFoundation
endif


VAPOR_INC=${VAPOR_INSTALL}/include
VAPOR_LIB=${VAPOR_INSTALL}/lib
VAPOR_BIN=${VAPOR_INSTALL}/bin

john_time_comp: john_time_comp.cpp
	$(CXX) $(CXXFLAGS) -I${VAPOR_INC} -I. john_time_comp.cpp -o bin/john_time_comp.o
	$(CXX) $(LDFLAGS) -L${VAPOR_LIB} bin/john_time_comp.o -o bin/john_time_comp
#	g++ -O2 -o bin/john_time_comp.o -c -std=c++0x -DLINUX -Wall -Wno-sign-compare  -D__USE_LARGEFILE64 -pthread -fPIC -m64 -I${VAPOR_INC} -I. -DMODELS john_time_comp.cpp
#	g++ bin/john_time_comp.o -o bin/john_time_comp -m64 -lrt -pthread   -Wl,-rpath, -L${VAPOR_LIB} -lwasp -lcommon 

cube3d.o: cube3d.cpp cube3d.h
	g++ -o bin/cube3d.o -c -std=c++0x -DLINUX -Wall -Wno-sign-compare  -D__USE_LARGEFILE64 -pthread -fPIC -m64 -g -I${VAPOR_INC} -I. -DMODELS cube3d.cpp

slicegroup.o: slicegroup.cpp slicegroup.h
	g++ -o bin/slicegroup.o -c -std=c++0x -DLINUX -Wall -Wno-sign-compare  -D__USE_LARGEFILE64 -pthread -fPIC -m64 -g -I${VAPOR_INC}  -I. -DMODELS slicegroup.cpp

temporal: temporal.cpp bin/cube3d.o bin/slicegroup.o
	g++ -o bin/temporal.o -c -std=c++0x -DLINUX -Wall -Wno-sign-compare  -D__USE_LARGEFILE64 -pthread -fPIC -m64 -g -I${VAPOR_INC} -I. -DMODELS temporal.cpp
	g++ bin/temporal.o bin/cube3d.o bin/slicegroup.o -o bin/temporal -m64 -lrt -pthread   -Wl,-rpath, -L${VAPOR_LIB} -lvdc -lwasp -lcommon 

clean:
	rm john_time_comp.o john_time_comp cube3d.o slicegroup.o temporal.o temporal
