#define rappture src directory
RP_BASE	= /apps/rappture/current

#defining the compiling environment
CC = gcc
CFLAGS = -I$(RP_BASE)/include

LIB_RAPPTURE = -L$(RP_BASE)/lib -lrappture -lexpat -lz

OBJECTS = archimedes-read-gui.o

archimedes-read-gui: $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LIB_RAPPTURE) -lm

install: archimedes-read-gui
	mv $< ../bin

clean: 
	rm -f *.o archimedes-read-gui

