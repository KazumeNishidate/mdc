CC= gcc
CFLAGS=  -Wall -O 
LIBS= -lm
SOURCES= fft1.c control.c files.c main.c 
OBJECTS= fft1.o control.o files.o main.o 
TARGET= skw
.c.o	:$(HEADERS) $(SOURCES)
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

$(TARGET)	:$(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean	:
	rm -f $(OBJECTS) *~ "#"* *.ln out
