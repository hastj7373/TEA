BUILDING INSTRUCTIONS

STEP 1: build mybamtools
	install CMake version >= 2.6.4
	extract bzipped tarball
		tar xjvf mybamtools.tbz
		cd mybamtools
	follow instructions in README:
		mkdir build
		cd build
		cmake ..
		make
	this will populate the lib folder with the requisite shared libraries
	the src folder IS still necessary--needed header files are stored here

STEP 2: build dre and bamreader
	must have mybamtools compiled
	extract bzipped tarball
		tar xjvf bamreader.tbz
		cd bamreader
	edit Makefile and set BTROOT to the path to which mybamtools was extracted
		vi Makefile
		BTROOT = /path/to/mybamtools
	compile the code
		make
	add path to mybamtools' compiled libraries to LD_LIBRARY_PATH
		in bash: export LD_LIBRARY_PATH=/path/to/mybamtools/lib
		in csh: setenv LD_LIBRARY_PATH /path/to/mybamtools/lib

clus:
	extract bzipped tarball
		tar xjvf sclus.tbz
		cd sclus
	no need to edit makefile
	compile the code
		make


