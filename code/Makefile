OUT = build/rsc2
CC = mpic++
IDIR = .
ODIR = obj
SDIR = src
INC = 
LIBS =
CFLAGS = -std=c++0x -g -Wall
LFLAGS = -lboost_serialization-mt -lboost_mpi-mt

_OBJS = Daemon.o Options.o FileUtil.o VecDataBlock.o rsc2.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

all: $(OBJS) $(OUT)


$(ODIR)/%.o: $(SDIR)/%.cpp 
	$(CC) -c $(INC) -o $@ $< $(CFLAGS)

$(OUT): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS) $(LIBS)
	strip $(OUT)

clean:
	rm -rf $(ODIR)/*
	rm -rf build/*

.PHONY: clean

install:
	sudo cp $(OUT) /usr/local/bin
