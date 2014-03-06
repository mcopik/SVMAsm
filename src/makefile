CC=g++
SOURCES=main.cpp
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))
FLAGS=-Wl
LDFLAGS=-ldl -pthread
CFLAGS=-g -m32 -Wall -std=c++11
OUTPUT=../bin/svmasm
all: ../bin $(OUTPUT)
	make -C libasm/ all

$(OUTPUT) : $(OBJECTS)
	$(CC) -g -m32 -Wall -o $(OUTPUT) $(OBJECTS) $(LDFLAGS)

../bin:
	mkdir $@

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS);
