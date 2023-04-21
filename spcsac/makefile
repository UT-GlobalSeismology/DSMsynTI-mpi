# makefile for spcsac

CC = /opt/intel/cc/9.0/bin/icc
CFLAGS = -L/usr/lib -g
OBJS = spcsac.o

all:
	make spcsac

spcsac: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o spcsac;
.c.o: 
	$(CC) -c $(FFLAGS) $<

clean:
	rm -f $(OBJS) spcsac

