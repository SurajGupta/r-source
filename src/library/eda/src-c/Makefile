LIBNAME=eda
LD=ld

OBJS=	line.o \
	smooth.o

eda.so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)
	@mv $(LIBNAME).so ../../../../lib

clean:
	@rm -f *.o *.so

realclean:
	@rm -f Makefile *.o *.so
