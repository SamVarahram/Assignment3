CC = gcc
LD = gcc
CFLAGS = -O3 -g -Wall -Werror
LDFLAGS = 
RM = /bin/rm -f
OBJS = galsim.o
EXECUTABLE = galsim

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE) -lm

galsim.o: galsim.c 
	$(CC) $(CFLAGS) -c galsim.c 

clean:
	$(RM) $(EXECUTABLE) $(OBJS)