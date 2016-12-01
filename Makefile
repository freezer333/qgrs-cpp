CC = g++
CFLAGS = --std=c++11
SRC = src/
BIN = bin/

all: qgrs

qgrs: default.o qgrs.o
	$(CC) $(BIN)default.o $(BIN)qgrs.o -o $(BIN)qgrs $(CFLAGS)

default.o: $(SRC)default.cpp
	$(CC) -c $(SRC)default.cpp -o $(BIN)default.o $(CFLAGS)

qgrs.o: $(SRC)qgrs.cpp
	$(CC) -c $(SRC)qgrs.cpp  -o $(BIN)qgrs.o $(CFLAGS)

clean:
	rm $(BIN)*
