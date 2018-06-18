CC = g++
CXXFLAGS = -std=c++0x -Wall -O3
LDFLAGS = -I /home/jordangumm/Inbox/uofm/devel/autofizkin/dependencies/miniconda/pkgs/jellyfish-2.2.7-h14c3975_0/include/jellyfish-2.2.7
#LDFLAGS = -Wl,--rpath=$(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L//g')
#LDLIBS = $(shell pkg-config --libs jellyfish-2.0)

all: query_per_sequence
query_per_sequence: query_per_sequence.cc sequence_mers.hpp
clean:
	rm -f *.o query_per_sequence


