CC=g++
SOURCES=foo.cpp
OBJECTS=foo.o
FLAGS=-Wl
LDFLAGS=-shared
CFLAGS=-g -fPIC -Wall
OUTPUT=../bin/libcpp.so
TARGET=libcpp.so
all:	$(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) -shared -Wl,-soname,$(TARGET).1 \
	-o $(OUTPUT) $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@
