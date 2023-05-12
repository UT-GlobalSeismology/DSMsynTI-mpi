# makefile for spcsac

CC = icx
#CFLAGS = -L/usr/lib -g
OBJS = spcsac.o

all:
	make spcsac

spcsac: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o spcsac;

.c.o: 
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f $(OBJS) spcsac

