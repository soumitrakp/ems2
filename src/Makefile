#CC=gcc
CXX=g++
RM=rm -f
#CXXFLAGS=-Wall -std=c++11 -g
CXXFLAGS=-Wall -std=c++11 -O3
LDFLAGS=
LDLIBS=

SRCS=ems.hpp \
		 utils.h \
		 ems1.hpp \
		 motif_tree.hpp \
		 ems2.hpp

.PHONY: all clean dist-clean

all: ems

ems.o: ems.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) -c $<

ems: ems.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

clean:
	$(RM) ems.o

dist-clean: clean
	$(RM) ems 
