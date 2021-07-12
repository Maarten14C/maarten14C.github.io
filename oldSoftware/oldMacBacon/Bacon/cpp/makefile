
OFILES =  input.o  Matrix.o  ranfun.o vector.o kernel.o

HEAD = bacon.h input.h cal.h   ranfun.h Matrix.h twalk.h
OPTI =  -O -D_FILE_OFFSET_BITS=64 # to use 64 bit file access  -DHAVE_INLINE  # -DHIDE_INLINE_STATIC
#-g debugger, -O to optimize

#include to save accepted iterations only
##DEFS = -DSAVEACCONLY #-DFLUSHEVERYIT to flush at every iteration
## No longer used, now is a parameter in the twalk object 

LIBS = -lgsl -lgslcblas -lm 

COMPILER = g++ 

CFLAGS = $(DEFS) $(OPTI) -L/usr/local/lib  -Wno-write-strings
#-Wno-deprecated 
###### Targets
all :  bacon hist events
	echo "Done all."


bacon : $(HEAD) bacon.o $(OFILES) 
	$(COMPILER)  $(CFLAGS) -o bacon bacon.o $(OFILES) $(LIBS)
	mv bacon ../bin/.

hist : hist.o
	$(COMPILER)  $(CFLAGS) -o hist hist.o $(OFILES) $(LIBS)
	mv hist ../bin/.

events : events.o
	$(COMPILER)  $(CFLAGS) -o events events.o $(OFILES) $(LIBS)
	mv events ../bin/.


%.o: %.c $(HEAD)
	$(COMPILER) $(CFLAGS) -c -o $@ $< 

