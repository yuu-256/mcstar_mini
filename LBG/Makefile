include Mkinclude

OBJS = $(GOBJS)

.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FFLAGS) $<

all: $(OBJS)
	cp $(OBJS) *.mod ../

clean: 
	rm -f *.o *.mod a.out *~ *__genmod.*

