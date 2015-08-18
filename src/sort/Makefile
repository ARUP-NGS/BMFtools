CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		-DHAVE_CONFIG_H
OBJS=		error.o lh3sort.o
PROG=		lh3sort
LIBS=		

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all:$(PROG)

sort:$(OBJS)
		$(CC) $(CXXFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)


clean:
		rm -f *.o a.out *~ *.a $(PROG)
