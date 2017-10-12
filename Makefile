CPP=g++ 
CC=gcc
CPPFLAGS = -O3 -Wall -Wmissing-prototypes -Wshadow -fmessage-length=0 -std=c++11 -msse2 -mfpmath=sse
CFLAGS = -O3

INC = -I/usr/local/include
PROGRAM = seqfilter

# Headers
HDR = utils.h hmm.h SeqFilter.h Sequence.h

# Source
CPPS = SeqFilter.cpp Options.cpp Sequence.cpp
CPPO = SeqFilter.o Options.o Sequence.o
SOURCE = utils.c hmm.c
COBJS = utils.o hmm.o



all : $(PROGRAM)

$(COBJS) : $(HDR) $(SOURCE)
	$(CC) $(CFLAGS) $(INC) -c $(SOURCE)

$(CPPO) : $(HDR) $(CPPS)
	$(CPP) $(CPPFLAGS) $(INC) -c $(CPPS)

seqfilter : $(CPPO) $(COBJS)
	$(CPP) $(CPPFLAGS) $(INC) $(LIB) $(COBJS) $(CPPO) -o $(PROGRAM)


clean:
	rm -f $(COBJS) $(CPPO)
	rm -f $(PROGRAM)
	rm -f core

