# makefile for spcsac

CC = cc
#CFLAGS = -L/usr/lib -g
LDLIBS = -lm
OBJS = spcsac.o

all:
	make spcsac

spcsac: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDLIBS) -o spcsac

.c.o: 
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f $(OBJS) spcsac

