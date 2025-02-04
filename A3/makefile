CC = gcc
LD = gcc
CFLAGS = -O3 -Wall #-Werror
LDFLAGS = -L/usr/X11R6/lib -lX11 -lm
RM = /bin/rm -f
OBJS = galsim.o
EXECUTABLE = galsim

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(EXECUTABLE)

galsim.o: galsim.c 
	$(CC) $(CFLAGS) -c galsim.c 

clean:
	$(RM) $(EXECUTABLE) $(OBJS)