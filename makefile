CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -ggdb -fprofile-arcs -ftest-coverage
SOURCES = functions.c main.c
EXECUTABLE = main
LIBS = -lm

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(SOURCES) -o $(EXECUTABLE) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(EXECUTABLE)