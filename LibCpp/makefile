CC=g++
SOURCES=foo.cpp
OBJECTS=foo.o
FLAGS=-Wl
LDFLAGS=-shared
CFLAGS=-g -fPIC -Wall -m32
OUTPUT=../bin/libcpp.so
TARGET=libcpp.so
all:	$(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) -shared -m32 -Wl,-soname,$(TARGET).1 \
	-o $(OUTPUT) $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@
