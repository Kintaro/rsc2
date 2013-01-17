OUT = build/rsc2
CC = mpic++
IDIR = .
ODIR = obj
SDIR = src
INC = 
LIBS =
CFLAGS = -std=c++0x -g -Wall -Isrc/ -Wno-c++11-extensions
LFLAGS = -lboost_serialization -lboost_mpi

_OBJS = Daemon.o Options.o FileUtil.o VecDataBlock.o rsc2.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

all: $(OBJS) $(OUT)


$(ODIR)/%.o: $(SDIR)/%.cpp 
	$(CC) -c $(INC) -o $@ $< $(CFLAGS)

$(OUT): $(OBJS)
	$(MAKE) -C src/index_structures/sash
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS) $(LIBS)
	strip $(OUT)

clean:
	rm -rf $(ODIR)/*
	rm -rf build/*

.PHONY: clean

install:
	sudo cp $(OUT) /usr/local/bin
