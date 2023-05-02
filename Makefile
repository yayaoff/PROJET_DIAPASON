GMSH_DIR := ../gmsh-sdk
CC := gcc
LIB := $(GMSH_DIR)/lib/libgmsh.so

PROG := project.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c
OBJS := project.o elasticity.o lu.o matrix.o eigen.o design.o
LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib -lm

fec:
	make clean
	$(CC) -c $(PROG)
	$(CC) -o project $(OBJS) $(LIB) $(LDFLAGS)
	rm -f *.o

clean:
	rm -f project
	rm -f *.o