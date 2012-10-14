CC=g++
SOURCES=main.cpp test/test.cpp
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))
FLAGS=-Wl
LDFLAGS=-ldl
CFLAGS=-g -Wall -m32 -std=c++11
OUTPUT=../bin/svmasm
all:	$(OUTPUT)

$(OUTPUT) : $(OBJECTS)
	$(CC) -g -m32 -Wall $(LDFLAGS) -o $(OUTPUT) $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS);
