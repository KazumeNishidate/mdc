CC= gcc
CFLAGS= -Wall -O2

###  where is the X include file "Xlib.h" and "Xutil.h"?
INC=  -I/usr/include/X11R5 -I/usr/include/X11 -I/usr/X11R6/include/X11 -I/opt/X11/include

###  where is the X library "libX11.a"?
LIBS= -lm -L/usr/X11R6/lib -lX11

## to display in X Window
RUN_OPTIONS= -DXD

SOURCES= main.c ext.c pt.c control.c init.c\
	 files.c real.c rcprcl.c rv.c msd.c\
	 ./xsrc/xdmain.c ./xsrc/xdinit.c ./xsrc/xdicon.c\
	 ./xsrc/xhist.c ./xsrc/xmsd.c ./xsrc/xnet.c\
	 ./xsrc/network.c ./xsrc/xunit.c 
OBJECTS= main.o ext.o pt.o control.o init.o\
	 files.o real.o rcprcl.o rv.o msd.o\
	 ./xsrc/xdmain.o ./xsrc/xdinit.o ./xsrc/xdicon.o\
	 ./xsrc/xhist.o ./xsrc/xmsd.o ./xsrc/xnet.o\
	 ./xsrc/network.o ./xsrc/xunit.o 

TARGET= md 

.c.o	:$(SOURCES)
	$(CC) -c $(CFLAGS) $(INC) $(RUN_OPTIONS) $< -o $@ 

$(TARGET)	:$(OBJECTS)
	$(CC) $(CFLAGS) $(RUN_OPTIONS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean	:
	rm -f $(TARGET) $(OBJECTS) *~ "#"* *.ln ./xsrc/*~ ./xsrc/"#"* ./xsrc/*.ln

