
john_time_comp: john_time_comp.cpp
	g++ -o john_time_comp.o -c -std=c++0x -DLINUX -Wall -Wno-sign-compare  -D__USE_LARGEFILE64 -pthread -fPIC -m64 -g -I/home/samuel/Tools/vapor-git/install/include -I. -DMODELS john_time_comp.cpp
	g++ john_time_comp.o -o john_time_comp -m64 -lrt -pthread   -Wl,-rpath, -L/home/samuel/Tools/vapor-git/install/lib -L/home/samuel/Tools/vapor-git/install/bin  -lvdc -lwasp -lcommon 

cube3d: cube3d.cpp cube3d.h
	g++ -o cube3d.o -c -std=c++0x -DLINUX -Wall -Wno-sign-compare  -D__USE_LARGEFILE64 -pthread -fPIC -m64 -g -I/home/samuel/Tools/vapor-git/install/include -I. -DMODELS cube3d.cpp

slicegroup: slicegroup.cpp slicegroup.h
	g++ -o slicegroup.o -c -std=c++0x -DLINUX -Wall -Wno-sign-compare  -D__USE_LARGEFILE64 -pthread -fPIC -m64 -g -I/home/samuel/Tools/vapor-git/install/include -I. -DMODELS slicegroup.cpp

clean:
	rm john_time_comp.o john_time_comp
