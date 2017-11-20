CPP=g++ 
CC=gcc
OPTIMISER = -O3
CPPFLAGS =  -Wall -Wshadow -fmessage-length=0 -std=c++11 -msse2 -mfpmath=sse
CFLAGS = 

INC = -I/usr/local/include
PROGRAM = seqfilter

# Headers
HDR = hmm.h SeqFilter.h Sequence.h ZorroInterface.h

# Source
CPPS = SeqFilter.cpp Options.cpp Sequence.cpp ZorroInterface.cpp
CPPO = SeqFilter.o Options.o Sequence.o ZorroInterface.o
SOURCE = hmm.c
COBJS = hmm.o



all : $(PROGRAM)

$(COBJS) : $(HDR) $(SOURCE)
	$(CC) $(OPTIMISER) $(CFLAGS) $(INC) -c $(SOURCE)

$(CPPO) : $(HDR) $(CPPS)
	$(CPP) $(OPTIMISER) $(CPPFLAGS) $(INC) -c $(CPPS)

seqfilter : $(CPPO) $(COBJS)
	$(CPP) $(CPPFLAGS) $(OPTIMISER) $(INC) $(LIB) $(COBJS) $(CPPO) -o $(PROGRAM)


clean:
	rm -f $(COBJS) $(CPPO)
	rm -f $(PROGRAM)
	rm -f core

