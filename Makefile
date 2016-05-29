#
#
# 'make'            make executable file 'mspar'
# 'make clean'      removes all .o and executable files
#

# Compiler
CC=mpicc

# Compilation flags
CFLAGS=-O3 -std=gnu99 -I.

# define any libraries to link into executable:
LIBS=-lm

# Dependencies
DEPS=ms.h mspar.h

# Folder to put the generated binaries
BIN=./bin

# Object files
OBJ=$(BIN)/mspar.o $(BIN)/ms.o $(BIN)/streec.o

# Random functions using drand48()
RND_48=rand1.c

# Random functions using rand()
RND=rand2.c

.PHONY: clean

$(BIN)/%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

default: $(BIN)/mspar

# download: packages
#	wget http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.2.tar.gz
#	tar -xf openmpi-1.8.2.tar.gz -C $(CURDIR)/packages

#packages:
#	mkdir packages

clean:
	rm -f $(BIN)/*
	@echo ""
	@echo "*** All resources were cleaned-up ***"
	@echo ""

$(BIN)/mspar: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(RND_48) $(LIBS)
	@echo ""
	@echo "*** make complete: generated executable 'mspar' ***"
