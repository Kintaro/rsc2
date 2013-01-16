OUT = ../../../build/libsash.so
CC = mpic++
IDIR = .
ODIR = ../../../obj
SDIR = .
INC = 
LIBS =
CFLAGS = -std=c++0x -g -Wall -shared -fPIC
LFLAGS = -lboost_serialization -lboost_mpi

_OBJS = Sash.o mtrand.o

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
