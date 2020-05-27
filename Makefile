MF=	Makefile

CC=	mpicc -cc=icc
CFLAGS=	-O3

LFLAGS= $(CFLAGS)

EXE=	percolate

INC= \
	percolate.h

SRC= \
	percolate.c \
	uni.c \
	arralloc.c \
	initmap.c \
	deliversmallmap.c \
	percprocess.c \
	percwritedynamic.c \
	setparameters.c 


#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ):	$(MF)

clean:
	rm -f $(EXE) $(OBJ) core *.o map.pgm
