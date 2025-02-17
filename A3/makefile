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

run:
	./$(EXECUTABLE) 3000 /home/sava3149/Assignment3/input_data/ellipse_N_03000.gal 100 1e-5 0