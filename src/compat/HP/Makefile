CFLAGS=-g -I. -Aa

OBJS=dlfcn.o vsystem.o

LIB= ../../lib/libHPcompat.a

$(LIB): ${OBJS}
	ar cru ${LIB} ${OBJS}
	cp dlfcn.h ../../include

clean:
	rm -f *.o
