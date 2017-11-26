#
# linux
LIBS = -llapack -lblas -lm


#CFLAGS = -O2
CFLAGS  = -g
CC = gcc -c
LD = gcc


HEADERS = gkGlobal.h  gk.h
OBJ = gkDisplay.o gkSetup.o input.o gkGlobal.o expan.o moments.o numQuad.o fmm.o kernel.o solver.o direct.o

%.o: %.c $(HEADERS) Makefile
	$(CC) $(CFLAGS) -o $@ $<

coulomb: coulomb.o $(OBJ)
	$(LD) -o coulomb coulomb.o $(OBJ) $(LIBS) $(LIBSLA)

clean:
	\rm -f *.o
