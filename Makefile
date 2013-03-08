CC = gcc
CFLAGS = -Wall -O3

gless: gless.o main.c
	$(CC) $^ -o $@

gless.o: gless.c gless.h
	$(CC) gless.c -c

clean:
	rm -f *.o

