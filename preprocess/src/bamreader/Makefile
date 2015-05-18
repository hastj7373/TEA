CXX = g++

BTROOT = /data/home/jluquette/lixing/from-alice/mybamtools/
BTINC = $(BTROOT)/src
BTLIB = $(BTROOT)/lib

PROFILE = #-pg

# Use this for debugging
#CXXFLAGS = -c -g $(PROFILE) -I $(BTINC)
CXXFLAGS = -c -O3 -I $(BTINC)

all: bamreader

bamreader: bamreader.o gzstream.o ReadGroup.o Histogram.o
	$(CXX) -L $(BTLIB) $(PROFILE) -o bamreader Histogram.o ReadGroup.o bamreader.o gzstream.o -lbamtools -lbamtools-utils

bamreader.o: bamreader.cpp
	$(CXX) $(CXXFLAGS) bamreader.cpp
gzstream.o: gzstream.cpp
	$(CXX) $(CXXFLAGS) gzstream.cpp
Histogram.o: Histogram.cpp Histogram.h
	$(CXX) $(CXXFLAGS) Histogram.cpp
ReadGroup.o: ReadGroup.cpp ReadGroup.h
	$(CXX) $(CXXFLAGS) ReadGroup.cpp

clean:
	rm *.o bamreader
