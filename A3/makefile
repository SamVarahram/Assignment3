CC = gcc
LD = gcc
CFLAGS = -O3 -Wall #-Werror
INCLUDES = -I../graphics
LDFLAGS = -L../graphics -lX11 -lm
RM = /bin/rm -f
OBJS = galsim.o ../graphics/graphics.o
EXECUTABLE = galsim

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE) -lm

galsim.o: galsim.c 
	$(CC) $(CFLAGS) $(INCLUDES) -c galsim.c 

clean:
	$(RM) $(EXECUTABLE) $(OBJS)