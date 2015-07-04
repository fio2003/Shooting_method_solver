#This should bild everything
CC = gcc
CFLAGS = -O3 -flto -Wall -fmessage-length=0  -std=gnu11
LIBS = -lm

all: shooting.o 
	$(CC) $(CFLAGS) shooting.o -o shooting.run $(LIBS)
green: shooting.c
	$(CC) $(CFLAGS) -c shooting.c
clean:
	rm -f *.o *.run