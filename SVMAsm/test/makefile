CC=g++
SOURCES=AllTests.cpp
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))
FLAGS=-Wl
LDFLAGS=-ldl -pthread #../contrib/gtest.a #../contrib/gtest_main.a
CFLAGS=-g -Wall -m32 -std=c++11  -I../contrib -I..
OUTPUT=tests
all:	$(OUTPUT)

$(OUTPUT) : $(OBJECTS)
	$(CC) -g -m32 -Wall $(LDFLAGS) -o $(OUTPUT) $(OBJECTS) ../contrib/gtest.a

.cpp.o:
	$(CC) $(CFLAGS) ../contrib/gtest.a -c $< -o $@

clean:
	rm -rf $(OBJECTS);
